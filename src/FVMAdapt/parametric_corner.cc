// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/parametric_cell.h>
#include "parametric_mixed.h"
#include <FVMAdapt/parametric/face_plan.h>
#include <FVMAdapt/parametric/surface_map.h>
#include "library/parametric/record_io.h"
#include "library/parametric/exact_geometry.h"
#include <functional>
#include <set>

namespace Loci {
namespace FVMAdapt {
  namespace {
    constexpr auto absent = noBoundary ;
    void require(bool ok, const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    RecordLimits historyRecordLimits(HistoryLimits h) {
      RecordLimits r ; r.maxBytes=h.exact.maxBytes ; r.maxNodes=h.maxNodes ; r.maxSplits=h.maxSplits ; r.maxDepth=h.maxDepth ;
      r.geometry=h.geometry ; r.geometry.exact=h.exact ; r.jacobian=h.jacobian ; return r ;
    }
    void bytes(detail::RecordWriter &w,const ParametricRecord &r) { w.u64(r.size()) ; w.append(r) ; }
    ParametricRecord bytes(detail::RecordReader &r,std::size_t limit) {
      ParametricRecord out(r.count(limit)) ; for(auto &b:out) b=r.byte() ; return out ;
    }

    bool nearlyStraightCorner(const std::vector<Point> &points) {
      if(points.size()<4) return false ;
      long double scale=0 ;
      for(const auto &p:points) for(double x:p) scale=std::max(scale,std::abs(static_cast<long double>(x))) ;
      const auto tolerance=64*std::numeric_limits<double>::epsilon()*scale ;
      for(unsigned i=0;i<points.size();++i) {
        std::array<long double,3> a,b,cross ; long double dot=0 ;
        for(unsigned d=0;d<3;++d) {
          a[d]=static_cast<long double>(points[i][d])-points[(i+points.size()-1)%points.size()][d] ;
          b[d]=static_cast<long double>(points[(i+1)%points.size()][d])-points[i][d] ; dot+=a[d]*b[d] ;
        }
        if(dot<=0) continue ;
        for(unsigned d=0;d<3;++d) cross[d]=a[(d+1)%3]*b[(d+2)%3]-a[(d+2)%3]*b[(d+1)%3] ;
        const auto first=std::hypot(a[0],a[1],a[2]),second=std::hypot(b[0],b[1],b[2]) ;
        const auto area=std::hypot(cross[0],cross[1],cross[2]) ;
        // Coordinate roundoff alone must not classify a resolved turn as straight.
        if(area<=tolerance*(first+second) && area<=std::sqrt(std::numeric_limits<double>::epsilon())*first*second) return true ;
      }
      return false ;
    }
    bool needsFaceFan(Shape shape,const std::vector<Point> &points) {
      if(!nearlyStraightCorner(points)) return false ;
      // Keep a planar quad's original coordinates when its interior is certified regular.
      return shape!=Shape::quad || !detail::hasStraightQuadCorner(points) ||
        geometryValidity(CornerMap(shape,points))!=GeometryValidity::valid ;
    }
    CornerCell sourceBoundary(CornerBoundary &state,const CellTopology &t,const std::vector<Point> &corners,bool patches=false) {
      CornerCell root ;
      require(corners.size()==t.vertexCount(),"source corner count differs from its shape") ;
      for(const auto &p:corners) root.vertices.push_back(state.addPoint(p)) ;
      std::map<std::pair<BoundaryId,BoundaryId>,BoundaryId> edges ;
      for(unsigned i=0;i<t.edgeCount();++i) {
        const auto ends=t.edgeVertices(i) ; const auto a=ends[0],b=ends[1] ;
        edges.emplace(state.endpoints(a,b),state.addEdge(a,b,std::make_shared<CornerMap>(Shape::segment,std::vector<Point>{corners[a],corners[b]}))) ;
      }
      for(unsigned f=0;f<t.faceCount();++f) {
        const auto loop=t.faceVertices(f) ; const auto shape=t.faceShape(f) ;
        std::vector<BoundaryId> vertices(loop.size()),faceEdges ; std::vector<Point> points(loop.size()) ;
        for(unsigned v=0;v<loop.size();++v) {
          const auto j=orientedVertex(unsigned(loop.size()),t.faceOrientation(f),v) ; vertices[j]=loop[v] ; points[j]=corners[loop[v]] ;
        }
        for(unsigned i=0;i<vertices.size();++i) faceEdges.push_back(edges.at(state.endpoints(vertices[i],vertices[(i+1)%vertices.size()]))) ;
        std::shared_ptr<const GeometryMap> geometry=patches ? sourceFaceGeometry(shape,points):shape==Shape::polygon
          ? std::shared_ptr<const GeometryMap>(std::make_shared<PolygonMap>(points)):std::make_shared<CornerMap>(shape,points) ;
        root.faces.push_back({state.addSurface(std::move(vertices),std::move(faceEdges),std::move(geometry)),t.faceOrientation(f)}) ;
      }
      return root ;
    }
    CornerSubdivision boundarySubdivision(const CellTopology &t,const std::vector<Point> &corners,HistoryLimits limits) {
      CornerBoundary boundary ; auto root=sourceBoundary(boundary,t,corners,true) ;
      Point interior=corners.at(0) ; std::array<long double,3> sum{} ;
      for(const auto &p:corners) for(unsigned d=0;d<3;++d) sum[d]+=static_cast<long double>(p[d])-corners[0][d] ;
      for(unsigned d=0;d<3;++d) interior[d]+=sum[d]/corners.size() ;
      return CornerSubdivision(t.id(),std::move(boundary),std::move(root),interior,limits) ;
    }

  }

  std::shared_ptr<const GeometryMap> sourceFaceGeometry(Shape shape,const std::vector<Point> &points) {
    if(needsFaceFan(shape,points)) {
      Point center=points.front() ; std::array<long double,3> sum{} ;
      for(const auto &p:points) for(unsigned d=0;d<3;++d) sum[d]+=static_cast<long double>(p[d])-points.front()[d] ;
      for(unsigned d=0;d<3;++d) center[d]+=sum[d]/points.size() ;
      return std::make_shared<FanMap>(shape,points,center) ;
    }
    return shape==Shape::polygon ? std::shared_ptr<const GeometryMap>(std::make_shared<PolygonMap>(points)):
      std::make_shared<CornerMap>(shape,points) ;
  }
  bool needsFacePatches(const CellTopology &topology,const std::vector<Point> &corners) {
    require(corners.size()==topology.vertexCount(),"source corner count differs from its topology") ;
    for(unsigned f=0;f<topology.faceCount();++f) {
      const auto loop=topology.faceVertices(f) ; std::vector<Point> points(loop.size()) ;
      for(unsigned v=0;v<loop.size();++v)
        points[orientedVertex(unsigned(loop.size()),topology.faceOrientation(f),v)]=corners[loop[v]] ;
      if(needsFaceFan(topology.faceShape(f),points)) return true ;
    }
    return false ;
  }

  ParametricCorner::ParametricCorner(const CellTopology &t,const std::vector<Point> &corners,HistoryLimits limits)
    : topology(t),subdivision(std::make_shared<CornerMap>(t.shape(),corners),t.id(),0,limits) {
    require(t.shape()==Shape::hex || t.shape()==Shape::prism || t.shape()==Shape::tet,"corner splitting requires hex, prism or tet source cells") ;
    require(validParametricSource(CornerMap(t.shape(),corners),limits.jacobian),"invalid or unresolved corner source geometry") ;
    options.maxNodes=limits.maxNodes ; options.maxSnapshotBytes=limits.exact.maxBytes ; options.jacobianLimits=limits.jacobian ;
    cells.push_back(sourceBoundary(*this,t,corners)) ;
  }
  BoundaryId ParametricCorner::faceId(EntityId face) const {
    for(unsigned f=0;f<cells[0].faces.size();++f) if(topology.faceId(f)==face) return cells[0].faces[f].surface ;
    throw std::invalid_argument("corner face is not incident to the original cell") ;
  }
  void ParametricCorner::refine(TreeNodeId leaf,CenterMethod method) {
    require(subdivision.history().active(leaf) && subdivision.node(leaf).isLeaf(),"corner refinement requires an active leaf") ;
    // Ordinary refinement reuses an accepted split, including after complete coarsening.
    for(const auto &split:subdivision.history().splits()) if(split.parent==leaf) { subdivision.reactivate(leaf,*split.construction) ; return ; }
    std::vector<std::size_t> selected(boundaries.faces.size(),absent) ;
    for(unsigned f=0;f<boundaries.faces.size();++f) if(!boundaries.faces[f].splits.empty()) selected[f]=0 ;
    splitCell(subdivision,cells,leaf,method,selected) ;
  }
  ParametricRecord ParametricCorner::faceRecord(EntityId face) const {
    return encodeFacePlan(boundaries,faceId(face),boundaryLimits()) ;
  }
  bool ParametricCorner::faceChanged(EntityId face) const {
    const auto &f=boundaries.faces[faceId(face)] ;
    if(!f.splits.empty()) return true ;
    for(auto e:f.edges) if(!boundaries.edges[e].cuts.empty()) return true ;
    return false ;
  }
  bool ParametricCorner::mergeFace(EntityId face,const ParametricRecord &record) {
    return applyFacePlan(*this,faceId(face),record) ;
  }
  ParametricRecord ParametricCorner::record() const {
    const auto limits=historyRecordLimits(subdivision.history().limits()) ; detail::RecordWriter w(limits.maxBytes,limits.geometry.exact) ;
    w.u32(0x434f524e) ; w.u32(1) ; const auto t=encodeParametricTopology(topology) ; w.u64(t.size()) ; for(int i:t) w.u32(i) ;
    std::vector<TreeNodeId> encodedIds ;
    bytes(w,encodeRefinement(subdivision.history(),encodedIds,limits)) ; bytes(w,encodeBoundaries(boundaries,boundaryLimits())) ;
    std::vector<TreeNodeId> order(encodedIds.size()) ;
    for(TreeNodeId id=0;id<encodedIds.size();++id) order[encodedIds[id]]=id ;
    w.u64(order.size()) ;
    for(TreeNodeId id:order) {
      const auto &cell=cells.at(id) ;
      w.u64(cell.vertices.size()) ; for(auto v:cell.vertices) w.u64(v) ;
      w.u64(cell.faces.size()) ; for(auto f:cell.faces) { w.u64(f.surface) ; w.u32(f.orientation) ; }
    }
    return w.take() ;
  }
  ParametricCorner ParametricCorner::restore(const ParametricRecord &record,HistoryLimits limits) {
    const auto bounds=historyRecordLimits(limits) ; detail::RecordReader r(record,limits.exact.maxBytes,limits.exact) ;
    require(r.u32()==0x434f524e && r.u32()==1,"unsupported corner history version") ;
    std::vector<int> t(r.count(128,4)) ; for(auto &i:t) i=r.u32() ;
    auto topology=decodeParametricTopology(t) ; auto history=decodeRefinement(bytes(r,bounds.maxBytes),bounds) ;
    require(history.rootKey()==topology.id() && history.geometryRevision()==0,"corner history source differs") ;
    ParametricCorner cell(std::move(topology),MappedVolumeSubdivision(std::move(history))) ;
    cell.options.maxNodes=limits.maxNodes ; cell.options.maxSnapshotBytes=limits.exact.maxBytes ; cell.options.jacobianLimits=limits.jacobian ;
    cell.boundaries=decodeBoundaries(bytes(r,bounds.maxBytes),cell.boundaryLimits()) ;
    const auto count=r.count(limits.maxNodes) ; require(count==cell.subdivision.nodes().size(),"corner connectivity count differs from history") ;
    for(TreeNodeId id=0;id<count;++id) {
      CornerCell c ; c.vertices.resize(r.count(8,8)) ; for(auto &v:c.vertices) { v=r.u64() ; require(v<cell.boundaries.points.size(),"invalid corner vertex") ; }
      c.faces.resize(r.count(6,12)) ; const auto &e=referenceElement(cell.subdivision.geometry(id)->shape()) ;
      require(c.vertices.size()==e.vertexCount && c.faces.size()==e.boundaryCount,"corner connectivity shape differs") ;
      for(unsigned f=0;f<c.faces.size();++f) {
        auto &use=c.faces[f] ; use.surface=r.u64() ; use.orientation=r.u32() ;
        require(use.surface<cell.boundaries.faces.size() && use.orientation<orientationCount(e.boundaries[f].shape),"invalid corner face use") ;
        const auto &face=cell.boundaries.faces[use.surface] ; require(face.vertices.size()==e.boundaries[f].vertexCount,"corner face shape differs") ;
        for(unsigned v=0;v<face.vertices.size();++v)
          require(c.vertices[e.boundaries[f].vertices[v]]==face.vertices[orientedVertex(e.boundaries[f].shape,use.orientation,v)],"corner face incidence differs") ;
      }
      cell.cells.push_back(std::move(c)) ;
    }
    r.finish() ; cell.validateCells(cell.subdivision,cell.cells) ;
    const auto &h=cell.subdivision.history() ; const auto *root=std::get_if<CornerRecord>(&h.geometry().record(h.node(0).map)) ;
    require(root && validParametricSource(CornerMap(root->shape,root->corners),limits.jacobian),"invalid or unresolved corner source geometry") ;
    require(root && root->shape==cell.topology.shape(),"corner source requires its original corner map") ;
    const auto &e=referenceElement(root->shape) ;
    for(unsigned f=0;f<e.boundaryCount;++f) {
      const auto use=cell.cells[0].faces[f] ; const auto &face=cell.boundaries.faces[use.surface] ;
      const auto *map=std::get_if<CornerRecord>(&cell.boundaries.geometry.record(face.geometry)) ;
      require(face.parent==absent && map && map->shape==e.boundaries[f].shape,"corner source face has no original map") ;
      for(unsigned v=0;v<e.boundaries[f].vertexCount;++v)
        require(map->corners[orientedVertex(map->shape,use.orientation,v)]==root->corners[e.boundaries[f].vertices[v]],"corner source face map differs from its cell") ;
    }
    return cell ;
  }

  ParametricBoundary::ParametricBoundary(const CellTopology &t,const std::vector<Point> &corners,HistoryLimits limits)
    : CornerSubdivision(boundarySubdivision(t,corners,limits)),topology(t) {}
  std::vector<Point> ParametricBoundary::sourceCorners() const {
    std::vector<Point> result ; for(auto v:node(0).cell.vertices) result.push_back(boundaries.points[v]) ; return result ;
  }
  ParametricRecord ParametricBoundary::splitRecord(TreeNodeId parent) const {
    const auto &split=splits().at(node(parent).split) ; const auto &cell=node(parent).cell ;
    detail::RecordWriter out(limits().exact.maxBytes,limits().exact) ;
    out.u32(unsigned(split.centerMethod)) ; out.u64(split.children.size()) ;
    for(auto x:boundaries.points[split.centerVertex]) out.real(x) ;
    out.u64(cell.faces.size()) ;
    for(unsigned f=0;f<cell.faces.size();++f)
      for(auto x:boundaries.faces[cell.faces[f].surface].splits[split.faceSplits[f]].referenceCenter) out.real(x) ;
    return out.take() ;
  }
  ParametricRecord ParametricBoundary::record() const {
    detail::RecordWriter out(limits().exact.maxBytes,limits().exact) ; out.u32(0x42434f52) ; out.u32(1) ;
    const auto t=encodeParametricTopology(topology) ; out.u64(t.size()) ; for(int i:t) out.u32(i) ;
    bytes(out,encode()) ; return out.take() ;
  }
  ParametricBoundary ParametricBoundary::restore(const ParametricRecord &record,HistoryLimits limits) {
    detail::RecordReader in(record,limits.exact.maxBytes,limits.exact) ;
    require(in.u32()==0x42434f52 && in.u32()==1,"unsupported boundary cell version") ;
    std::vector<int> t(in.count(limits.exact.maxBytes/4,4)) ; for(auto &i:t) i=in.u32() ;
    auto topology=decodeParametricTopology(t) ;
    auto cells=CornerSubdivision::decode(bytes(in,limits.exact.maxBytes),limits) ; in.finish() ;
    ParametricBoundary result(std::move(topology),std::move(cells)) ;
    require(result.rootKey()==result.topology.id(),"boundary history source differs") ;
    const auto &root=result.node(0).cell ; const auto &source=result.topology ;
    require(root.vertices.size()==source.vertexCount() && root.faces.size()==source.faceCount(),"boundary source shape differs") ;
    for(unsigned f=0;f<root.faces.size();++f) {
      const auto use=root.faces[f] ; const auto &face=result.boundaries.faces[use.surface] ; const auto expected=source.faceVertices(f) ; const auto shape=source.faceShape(f) ;
      require(use.orientation==result.topology.faceOrientation(f) && face.parent==absent,"boundary source face orientation differs") ;
      const auto &record=result.boundaries.geometry.record(face.geometry) ;
      const auto *map=std::get_if<CornerRecord>(&record) ; const auto *fan=std::get_if<FanRecord>(&record) ;
      require((map && map->shape==shape) || (fan && fan->shape==shape),"boundary source face map differs") ;
      const auto &corners=map ? map->corners:fan->corners() ;
      require(face.vertices.size()==expected.size() && corners.size()==expected.size(),"boundary source face size differs") ;
      for(unsigned v=0;v<expected.size();++v) {
        const auto canonical=orientedVertex(unsigned(expected.size()),use.orientation,v) ; const auto vertex=root.vertices[expected[v]] ;
        require(face.vertices[canonical]==vertex && corners[canonical]==result.boundaries.points[vertex],"boundary source incidence differs") ;
      }
      const auto sourceGeometry=sourceFaceGeometry(shape,corners) ;
      const auto *sourceFan=dynamic_cast<const FanMap *>(sourceGeometry.get()) ;
      require(bool(fan)==bool(sourceFan) && (!fan || fan->center()==sourceFan->center()),"saved boundary source geometry differs") ;
    }
    return result ;
  }

  namespace {
    void addOriginalEdgePoints(ParametricCell &cell,const std::vector<ParametricNodeKey> &points) {
      const auto &topology=cell.cornerTopology() ; const auto &root=cell.cornerCell(0) ;
      const auto edges=cell.cornerGeometry().cellEdges(root) ; CornerBoundary::Transaction transaction(cell.cornerGeometry()) ;
      for(const auto &key:points) if(key.reference.kind==EntityKind::edge) {
        for(unsigned e=0;e<topology.edgeCount();++e) if(cell.cornerTopology().edgeId(e)==key.reference.entity) {
          const auto ends=topology.edgeVertices(e) ; const auto a=ends[0],b=ends[1] ; const bool forward=cell.cornerTopology().vertexId(a)<cell.cornerTopology().vertexId(b) ;
          auto edge=edges.at(cell.cornerGeometry().endpoints(root.vertices[a],root.vertices[b])),start=root.vertices[forward ? a:b] ;
          Fraction lo=0,hi=1 ; const auto target=key.reference.point[0] ;
          require(target>lo && target<hi,"original edge point is outside its edge") ;
          for(unsigned depth=0;;++depth) {
            require(depth<31,"corner edge point exceeds supported depth") ;
            const auto midpoint=cell.cornerGeometry().splitEdge(edge,transaction) ; const auto mid=(lo+hi)/Fraction(2) ; if(target==mid) break ;
            const auto &record=cell.cornerGeometry().boundaries.edges[edge] ; const bool same=record.vertices[0]==start ;
            if(target<mid) { edge=record.children[same ? 0:1].edge ; hi=mid ; }
            else { edge=record.children[same ? 1:0].edge ; start=midpoint ; lo=mid ; }
          }
        }
      }
      transaction.committed=true ;
    }
    struct CornerView {
      const ParametricCell &cell ; const BoundaryRecords &b ;
      std::vector<ParametricNodeKey> keys ;
      struct Cover { int leaf ; bool right,reverse ; } ;
      std::vector<std::vector<Cover>> covers ; std::vector<bool> active,edgeSplit ;
      explicit CornerView(const ParametricCell &c) : cell(c),b(c.cornerGeometry().boundaries),covers(b.faces.size()),active(b.faces.size()),edgeSplit(b.edges.size()) {
        for(unsigned v=0;v<b.points.size();++v) keys.push_back({{EntityKind::cell,c.cornerTopology().id(),{{0,0,0}}},int(v)}) ;
        const auto &topology=c.cornerTopology() ; const auto &root=c.cornerCell(0) ;
        for(unsigned v=0;v<root.vertices.size();++v) keys[root.vertices[v]]={{EntityKind::vertex,c.cornerTopology().vertexId(v),{{0,0,0}}},-1} ;
        const auto edges=c.cornerGeometry().cellEdges(root) ;
        for(unsigned i=0;i<topology.edgeCount();++i) {
          const auto ends=topology.edgeVertices(i) ; const auto a=ends[0],d=ends[1] ; const bool forward=c.cornerTopology().vertexId(a)<c.cornerTopology().vertexId(d) ;
          edgeKeys(edges.at(c.cornerGeometry().endpoints(root.vertices[a],root.vertices[d])),root.vertices[forward ? a:d],c.cornerTopology().edgeId(i),0,1) ;
        }
        for(unsigned f=0;f<root.faces.size();++f) { int next=0 ; faceKeys(root.faces[f].surface,c.cornerTopology().faceId(f),next) ; }
        const auto leaves=c.leafIds() ;
        for(unsigned i=0;i<leaves.size();++i) for(auto use:c.cornerCell(leaves[i]).faces)
          addCover(use.surface,{int(i),false,use.orientation>=b.faces[use.surface].vertices.size()}) ;
      }
      void edgeKeys(BoundaryId edge,BoundaryId start,EntityId owner,Fraction lo,Fraction hi) {
        const auto &e=b.edges[edge] ; if(e.cuts.empty()) return ;
        const auto mid=(lo+hi)/Fraction(2) ; keys[e.cuts[0].vertex]={{EntityKind::edge,owner,{{mid,0,0}}},-1} ;
        const bool forward=e.vertices[0]==start ;
        edgeKeys(e.children[forward ? 0:1].edge,start,owner,lo,mid) ;
        edgeKeys(e.children[forward ? 1:0].edge,e.cuts[0].vertex,owner,mid,hi) ;
      }
      void facePoint(BoundaryId vertex,EntityId owner,int &next) {
        if(keys[vertex].reference.kind==EntityKind::cell) keys[vertex]={{EntityKind::face,owner,{{0,0,0}}},next++} ;
      }
      void faceEdgeKeys(BoundaryId edge,BoundaryId start,EntityId owner,int &next) {
        const auto &e=b.edges[edge] ; if(e.cuts.empty()) return ;
        facePoint(e.cuts[0].vertex,owner,next) ; const bool forward=e.vertices[0]==start ;
        faceEdgeKeys(e.children[forward ? 0:1].edge,start,owner,next) ;
        faceEdgeKeys(e.children[forward ? 1:0].edge,e.cuts[0].vertex,owner,next) ;
      }
      void faceKeys(BoundaryId face,EntityId owner,int &next) {
        const auto &f=b.faces[face] ;
        for(unsigned e=0;e<f.edges.size();++e) faceEdgeKeys(f.edges[e],f.vertices[e],owner,next) ;
        for(const auto &split:f.splits) {
          facePoint(split.centerVertex,owner,next) ;
          for(auto child:split.children) faceKeys(child,owner,next) ;
        }
      }
      void addCover(BoundaryId face,Cover cover) {
        covers.at(face).push_back(cover) ;
        while(b.faces[face].parent!=absent) { face=b.faces[face].parent ; active[face]=true ; }
      }
      void addRight(EntityId face,const ParametricCell &right) {
        std::map<BoundaryId,BoundaryId> ids ;
        std::function<void(BoundaryId,BoundaryId)> match=[&](BoundaryId a,BoundaryId r) {
          ids.emplace(r,a) ; const auto &rf=right.cornerGeometry().boundaries.faces[r] ;
          if(!rf.splits.empty()) for(unsigned i=0;i<rf.splits[0].children.size();++i)
            match(b.faces[a].splits[0].children[i],rf.splits[0].children[i]) ;
        } ;
        match(cell.cornerFaceId(face),right.cornerFaceId(face)) ; const auto leaves=right.leafIds() ;
        for(unsigned i=0;i<leaves.size();++i) for(auto use:right.cornerCell(leaves[i]).faces) {
          const auto found=ids.find(use.surface) ;
          if(found!=ids.end()) addCover(found->second,{int(i),true,use.orientation>=right.cornerGeometry().boundaries.faces[use.surface].vertices.size()}) ;
        }
      }
      void splitParents(BoundaryId edge) {
        while(b.edges[edge].parent!=absent) { edge=b.edges[edge].parent ; edgeSplit[edge]=true ; }
      }
      std::map<BoundaryId,std::vector<long long>> surfaceNames() const {
        std::map<BoundaryId,std::vector<long long>> names ;
        for(unsigned f=0;f<cell.cornerCell(0).faces.size();++f) names[cell.cornerCell(0).faces[f].surface]={2,static_cast<long long>(cell.cornerTopology().faceId(f))} ;
        for(std::size_t i=0;i<cell.splitCount();++i) {
          const auto path=cell.path(cell.splitParent(i)) ;
          std::vector<long long> name{4,static_cast<long long>(cell.cornerTopology().id()),static_cast<long long>(path.size())} ;
          for(auto child:path) name.push_back(child) ;
          unsigned ordinal=0 ;
          for(auto child:cell.splitChildren(i)) for(auto use:cell.cornerCell(child).faces) {
            if(b.faces[use.surface].parent!=absent || names.count(use.surface)) continue ;
            auto key=name ; key.push_back(ordinal++) ; names.emplace(use.surface,std::move(key)) ;
          }
        }
        return names ;
      }
      FacePatch facePatch(BoundaryId id,bool reverse,const std::map<BoundaryId,std::vector<long long>> &names) const {
        std::vector<unsigned> path ; auto root=id ;
        while(b.faces[root].parent!=absent) {
          const auto &face=b.faces[root] ; const auto &siblings=b.faces[face.parent].splits.at(face.parentSplit).children ;
          path.push_back(std::find(siblings.begin(),siblings.end(),root)-siblings.begin()) ; root=face.parent ;
        }
        using Pair=std::array<long long,2> ; constexpr long long extent=1LL<<30 ;
        std::array<Pair,4> corners{{{{0,0}},{{extent,0}},{{extent,extent}},{{0,extent}}}} ;
        FacePatch patch ; patch.surface=names.at(root) ; patch.sector=b.faces[root].vertices.size()==4 ? 0:-int(b.faces[root].vertices.size()) ;
        const auto midpoint=[](Pair a,Pair d) { return Pair{{a[0]+(d[0]-a[0])/2,a[1]+(d[1]-a[1])/2}} ; } ;
        // Integer coordinates encode the accepted tree, independently of physical center placement.
        for(auto it=path.rbegin();it!=path.rend();++it) {
          if(patch.sector<0) { patch.sector=*it ; continue ; }
          const auto old=corners ; const unsigned c=*it ;
          require(c<4 && std::abs(old[0][0]-old[2][0])>=2 && std::abs(old[0][1]-old[2][1])>=2,"corner patch exceeds coordinate depth") ;
          corners={{old[c],midpoint(old[c],old[(c+1)%4]),midpoint(old[0],old[2]),midpoint(old[c],old[(c+3)%4])}} ;
        }
        patch.bounds={{std::min(corners[0][0],corners[2][0]),std::min(corners[0][1],corners[2][1]),
                       std::max(corners[0][0],corners[2][0]),std::max(corners[0][1],corners[2][1])}} ;
        patch.orientation=reverse ? -1:1 ; patch.mappedArea=true ;
        const auto measure=geometryMeasure(*b.geometry.map(b.faces[id].geometry),cell.cornerGeometry().options.assessment.measure) ;
        require(measure.converged && measure.measure>0,"corner face area integration did not converge") ; patch.area=measure.measure ;
        return patch ;
      }
      std::vector<ParametricOutputFace> faces(EntityId sourceFace,bool interior,const std::vector<ParametricNodeKey> &required={},bool neighbor=false,
                                               bool recordPatches=false) {
        struct Patch { BoundaryId id ; Cover left ; int right ; } ; std::vector<Patch> patches ;
        std::vector<bool> selected(b.faces.size(),interior) ;
        if(interior) for(auto use:cell.cornerCell(0).faces) selected[use.surface]=false ;
        else selected[cell.cornerFaceId(sourceFace)]=true ;
        // Records are ordered parent before children; inherited cell covers need one forward pass.
        for(BoundaryId f=0;f<b.faces.size();++f) {
          const auto &face=b.faces[f] ; if(face.parent!=absent) selected[f]=selected[face.parent] ;
          if(!selected[f]) continue ;
          auto &incident=covers[f] ;
          if(active[f]) {
            require(face.splits.size()==1,"active corner face has no unique construction") ;
            for(auto child:face.splits[0].children) covers[child].insert(covers[child].end(),incident.begin(),incident.end()) ;
          } else if(!incident.empty()) {
            require(incident.size()==(interior || neighbor ? 2u:1u),"corner face coverage is incomplete") ;
            if(incident.front().right) std::reverse(incident.begin(),incident.end()) ;
            if(incident.size()==2) require(incident[0].reverse!=incident[1].reverse,"corner face orientations disagree") ;
            patches.push_back({f,incident[0],incident.size()==2 ? incident[1].leaf:-1}) ;
            for(auto edge:face.edges) splitParents(edge) ;
          }
        }
        const std::set<ParametricNodeKey> needed(required.begin(),required.end()) ;
        for(BoundaryId e=0;e<b.edges.size();++e) for(const auto &cut:b.edges[e].cuts) if(needed.count(keys[cut.vertex])) { edgeSplit[e]=true ; splitParents(e) ; }
        const auto names=recordPatches ? surfaceNames():std::map<BoundaryId,std::vector<long long>>{} ;
        std::vector<ParametricOutputFace> result ;
        for(const auto &p:patches) {
          const auto &f=b.faces[p.id] ; std::vector<BoundaryId> vertices ;
          for(unsigned e=0;e<f.edges.size();++e) cell.cornerGeometry().edgeVertices(f.edges[e],f.vertices[e],edgeSplit,vertices) ;
          if(p.left.reverse) std::reverse(vertices.begin(),vertices.end()) ;
          ParametricOutputFace face{{},p.left.leaf,p.right} ; for(auto v:vertices) face.vertices.push_back(keys[v]) ;
          if(recordPatches) face.patch=facePatch(p.id,p.left.reverse,names) ;
          result.push_back(std::move(face)) ;
        }
        return result ;
      }
    } ;
  }
  std::vector<ParametricOutputFace> cornerCellFaces(const ParametricCell &cell,const std::vector<ParametricNodeKey> &required,bool patches) {
    auto copy=cell ; addOriginalEdgePoints(copy,required) ; return CornerView(copy).faces(0,true,required,false,patches) ;
  }
  std::vector<ParametricOutputFace> cornerFaceFaces(EntityId face,const ParametricCell &left,
      const std::vector<ParametricNodeKey> &required,const ParametricCell *right,bool patches) {
    auto merged=left ; if(right) merged.mergeFace(face,parametricFacePlan(face,*right)) ;
    addOriginalEdgePoints(merged,required) ; CornerView view(merged) ; if(right) view.addRight(face,*right) ; return view.faces(face,false,required,bool(right),patches) ;
  }
  std::vector<std::vector<ParametricNodeKey>> cornerLeafVertices(const ParametricCell &cell) {
    const CornerView view(cell) ; std::vector<std::vector<ParametricNodeKey>> result ;
    for(auto leaf:cell.leafIds()) {
      result.emplace_back() ;
      for(auto vertex:cell.cornerCell(leaf).vertices) result.back().push_back(view.keys.at(vertex)) ;
    }
    return result ;
  }
  std::vector<ParametricNodeKey> cornerNodes(const ParametricCell &cell,EntityKind kind,EntityId owner,const ParametricCell *right) {
    const auto faces=kind==EntityKind::cell ? cornerCellFaces(cell,{}) : cornerFaceFaces(owner,cell,{},right) ;
    std::set<ParametricNodeKey> keys ;
    for(const auto &f:faces) for(const auto &key:f.vertices) if(key.reference.kind==kind && key.reference.entity==owner) keys.insert(key) ;
    return {keys.begin(),keys.end()} ;
  }
  std::vector<Point> cornerPositions(const ParametricCell &cell,const std::vector<ParametricNodeKey> &keys) {
    const CornerView view(cell) ; std::map<ParametricNodeKey,Point> positions ;
    for(unsigned i=0;i<view.keys.size();++i) positions.emplace(view.keys[i],view.b.points[i]) ;
    std::vector<Point> result ; for(const auto &key:keys) result.push_back(positions.at(key)) ; return result ;
  }
  ParametricTrace cornerEdgeTrace(const ParametricCell &c) {
    const CornerView view(c) ;
    const auto &t=c.cornerTopology() ;
    const auto leaves=c.leafIds() ;
    ParametricTrace out ;
    for(unsigned i=0;i<leaves.size();++i) for(unsigned edge=0;edge<t.edgeCount();++edge) {
      const auto vertices=t.edgeVertices(edge) ; const auto a=t.vertexId(vertices[0]),b=t.vertexId(vertices[1]) ;
      const auto ends=std::minmax(a,b) ; std::vector<Fraction> cuts ;
      for(auto v:c.cornerCell(leaves[i]).vertices) {
        const auto &p=view.keys[v].reference ;
        if(p.kind==EntityKind::edge && p.entity==t.edgeId(edge)) cuts.push_back(p.point[0]) ;
        if(p.kind==EntityKind::vertex && (p.entity==ends.first || p.entity==ends.second)) cuts.push_back(p.entity==ends.first ? 0:1) ;
      }
      if(cuts.empty()) continue ;
      const auto range=std::minmax_element(cuts.begin(),cuts.end()) ;
      out.push_back({int(t.edgeId(edge)),int(t.id()),int(i),int(c.path(leaves[i]).size()),
        checkedExactInteger<int>(range.first->numerator()),checkedExactInteger<int>(range.first->denominator()),
        checkedExactInteger<int>(range.second->numerator()),checkedExactInteger<int>(range.second->denominator())}) ;
    }
    return out ;
  }
  std::vector<int> cornerVertexDepths(const ParametricCell &c) {
    const auto &root=c.cornerCell(0).vertices ; std::vector<int> out(root.size()) ;
    for(auto leaf:c.leafIds()) {
      const auto &vertices=c.cornerCell(leaf).vertices ; const int depth=c.path(leaf).size() ;
      for(unsigned v=0;v<root.size();++v)
        if(std::find(vertices.begin(),vertices.end(),root[v])!=vertices.end()) out[v]=std::max(out[v],depth) ;
    }
    return out ;
  }
  std::vector<int> cornerVertexRequests(const ParametricCell &c,const std::map<EntityId,int> &depth,unsigned difference) {
    std::vector<int> out ; const auto leaves=c.leafIds() ; const auto &root=c.cornerCell(0).vertices ;
    for(unsigned i=0;i<leaves.size();++i) {
      const auto &vertices=c.cornerCell(leaves[i]).vertices ; const int own=c.path(leaves[i]).size() ;
      for(unsigned v=0;v<root.size();++v) {
        const int other=depth.at(c.cornerTopology().vertexId(v)) ;
        if(other>own && unsigned(other-own)>difference && std::find(vertices.begin(),vertices.end(),root[v])!=vertices.end()) {
          out.push_back(i) ; break ;
        }
      }
    }
    return out ;
  }
  std::map<ParametricNodeKey,ReferencePoint> cornerFaceCoordinates(EntityId face,const ParametricCell &cell) {
    const CornerView view(cell) ; const auto &b=view.b ; const auto root=cell.cornerFaceId(face) ;
    const unsigned count=b.faces[root].vertices.size() ;
    require(count==3 || count==4,"mixed faces require three or four source vertices") ;
    std::map<BoundaryId,ReferencePoint> points ;
    const auto assign=[&](BoundaryId id,const ReferencePoint &p) {
      const auto old=points.emplace(id,p) ; require(old.second || old.first->second==p,"corner face coordinates disagree") ;
    } ;
    for(unsigned v=0;v<count;++v) assign(b.faces[root].vertices[v],referenceVertex(count==3 ? Shape::triangle:Shape::quad,v)) ;
    std::function<void(BoundaryId)> edgePoints=[&](BoundaryId id) {
      const auto &edge=b.edges[id] ; const auto a=points.at(edge.vertices[0]),c=points.at(edge.vertices[1]) ;
      for(const auto &cut:edge.cuts) {
        ReferencePoint p ; for(unsigned d=0;d<3;++d) p[d]=a[d]+cut.position*(c[d]-a[d]) ; assign(cut.vertex,p) ;
      }
      for(const auto &child:edge.children) edgePoints(child.edge) ;
    } ;
    std::function<void(BoundaryId)> facePoints=[&](BoundaryId id) {
      const auto &f=b.faces[id] ; for(auto edge:f.edges) edgePoints(edge) ;
      for(const auto &split:f.splits) {
        const Point expected=f.vertices.size()==3 ? Point{{1./3.,1./3.,0}}:Point{{.5,.5,0}} ;
        require(split.referenceCenter==expected,"mixed faces currently require midpoint face centers") ;
        ReferencePoint center{} ;
        for(auto vertex:f.vertices) for(unsigned d=0;d<3;++d) center[d]=center[d]+points.at(vertex)[d]/Fraction(f.vertices.size()) ;
        assign(split.centerVertex,center) ;
        for(auto child:split.children) facePoints(child) ;
      }
    } ;
    facePoints(root) ; std::map<ParametricNodeKey,ReferencePoint> result ;
    for(const auto &entry:points) result.emplace(view.keys.at(entry.first),entry.second) ;
    return result ;
  }

}
}
