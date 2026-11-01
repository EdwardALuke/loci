// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include "parametric_mixed.h"

#include <FVMAdapt/parametric_cell.h>
#include "parametric_hex.h"
#include <FVMAdapt/parametric/face_plan.h>
#include <FVMAdapt/parametric/geometry_quality.h>
#include <FVMAdapt/parametric/geometry_measure.h>
#include "library/parametric/record_io.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>

namespace Loci {
namespace FVMAdapt {
  namespace {
    void require(bool ok, const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    RecordLimits recordLimits(const HistoryLimits &h) {
      RecordLimits r ; r.maxBytes=h.exact.maxBytes ; r.maxNodes=h.maxNodes ; r.maxSplits=h.maxSplits ;
      r.maxDepth=h.maxDepth ; r.geometry=h.geometry ; r.geometry.exact=h.exact ; r.jacobian=h.jacobian ; return r ;
    }
    constexpr std::uint32_t extrusionMagic=0x45585452 ;
    constexpr std::uint32_t sourceMagic=0x53524345 ;

    struct ExtrusionLayout {
      std::array<unsigned,2> caps ;
      std::array<double,2> planes ;
      std::vector<unsigned> lower, upper, sides ;
      std::vector<Point> polygon ;
      ExtrusionLayout(Shape shape, const std::vector<Point> &corners, int axis) {
        require(axis>=0 && axis<3,"center-based extrusion requires axis 0, 1 or 2") ;
        parametricGradingSplit(CornerMap(shape,corners),axis) ;
        const auto &e=referenceElement(shape) ;
        planes={{corners.front()[axis],corners.front()[axis]}} ;
        for(const auto &p:corners) { planes[0]=std::min(planes[0],p[axis]) ; planes[1]=std::max(planes[1],p[axis]) ; }
        require(planes[1]>planes[0],"extrusion has no thickness") ;
        for(unsigned layer=0;layer<2;++layer) {
          bool found=false ;
          for(unsigned f=0;f<e.boundaryCount;++f) {
            const auto &face=e.boundaries[f] ; bool match=true ;
            for(unsigned j=0;j<face.vertexCount;++j) match=match && corners[face.vertices[j]][axis]==planes[layer] ;
            if(match) { require(!found,"extrusion has duplicate end faces") ; caps[layer]=f ; found=true ; }
          }
          require(found,"extrusion requires fixed end planes") ;
        }
        const auto &bottom=e.boundaries[caps[0]] ;
        lower.assign(bottom.vertices.begin(),bottom.vertices.begin()+bottom.vertexCount) ; std::reverse(lower.begin(),lower.end()) ;
        for(unsigned v:lower) {
          unsigned partner=e.vertexCount ;
          for(unsigned w=0;w<e.vertexCount;++w)
            if(corners[w][axis]==planes[1] && corners[v][(axis+1)%3]==corners[w][(axis+1)%3] &&
               corners[v][(axis+2)%3]==corners[w][(axis+2)%3]) partner=w ;
          require(partner<e.vertexCount,"extrusion end vertices do not align") ; upper.push_back(partner) ;
          polygon.push_back({{corners[v][(axis+1)%3],corners[v][(axis+2)%3],0}}) ;
        }
        for(unsigned i=0;i<lower.size();++i) {
          unsigned side=e.boundaryCount ;
          for(unsigned f=0;f<e.boundaryCount;++f) if(f!=caps[0] && f!=caps[1]) {
            const auto &face=e.boundaries[f] ; std::set<unsigned> vertices(face.vertices.begin(),face.vertices.begin()+face.vertexCount) ;
            if(vertices.count(lower[i]) && vertices.count(lower[(i+1)%lower.size()])) side=f ;
          }
          require(side<e.boundaryCount,"extrusion side face is missing") ; sides.push_back(side) ;
        }
      }
    } ;
    std::vector<ParametricNodeKey> outputKeys(const std::vector<VertexKey> &keys) {
      std::vector<ParametricNodeKey> out ; for(const auto &key:keys) out.push_back({key,-1}) ; return out ;
    }
    std::vector<VertexKey> referenceKeys(const std::vector<ParametricNodeKey> &keys) {
      std::vector<VertexKey> out ; for(const auto &key:keys) { require(key.local<0,"native output received a center-based node") ; out.push_back(key.reference) ; } return out ;
    }
    std::vector<ParametricOutputFace> outputFaces(const std::vector<ParametricFace> &faces) {
      std::vector<ParametricOutputFace> out ; for(const auto &f:faces) out.push_back({outputKeys(f.vertices),f.left,f.right}) ; return out ;
    }
    // Native cuts lie in exact reference planes, even when their physical surfaces curve.
    void nativeFacePatches(std::vector<ParametricOutputFace> &faces,const CellTopology &topology,
                           const AffineSubdivision &cell,const EntityId *original = nullptr) {
      const auto geometry=cell.history().geometry().map(cell.history().node(0).map) ;
      unsigned boundary=0 ;
      if(original) {
        while(boundary<topology.faceCount() && topology.faceId(boundary)!=*original) ++boundary ;
        require(boundary<topology.faceCount(),"face is absent from native root") ;
      }
      for(auto &face:faces) {
        std::vector<ReferencePoint> points ;
        for(const auto &key:face.vertices) {
          ReferencePoint point ; require(topology.locateVertex(key.reference,point),"native face point is outside its root") ;
          points.push_back(point) ;
        }
        ReferencePoint normal{} ;
        for(unsigned i=1;i+1<points.size();++i) for(unsigned d=0;d<3;++d) {
          const unsigned a=(d+1)%3,b=(d+2)%3 ;
          normal[d]=normal[d]+(points[i][a]-points[0][a])*(points[i+1][b]-points[0][b])-
                     (points[i][b]-points[0][b])*(points[i+1][a]-points[0][a]) ;
        }
        unsigned axis=0 ; while(axis<3 && normal[axis]==0) ++axis ;
        require(axis<3,"native face has zero reference area") ;
        auto &patch=face.patch ;
        std::vector<ReferencePoint> chart ;
        if(original) {
          patch.surface={2,static_cast<long long>(*original)} ;
          for(const auto &point:points) chart.push_back(orientReferencePoint(topology.faceShape(boundary),topology.faceOrientation(boundary),
            boundaryCoordinates(topology.shape(),boundary,point))) ;
        } else {
          patch.surface={5,static_cast<long long>(topology.id())} ;
          const auto scale=normal[axis] ; Fraction offset ;
          for(unsigned d=0;d<3;++d) { normal[d]=normal[d]/scale ; offset=offset+normal[d]*points[0][d] ; }
          for(const auto &value:{normal[0],normal[1],normal[2],offset}) {
            patch.surface.push_back(checkedExactInteger<long long>(value.numerator())) ;
            patch.surface.push_back(checkedExactInteger<long long>(value.denominator())) ;
          }
          for(const auto &point:points) chart.push_back({{point[(axis+1)%3],point[(axis+2)%3],0}}) ;
        }
        Fraction area ;
        for(unsigned i=0;i<chart.size();++i) {
          const auto &a=chart[i],&b=chart[(i+1)%chart.size()] ; area=area+a[0]*b[1]-a[1]*b[0] ;
          patch.polygon.push_back({{a[0].value(),a[1].value()}}) ;
        }
        require(area!=0,"native face has zero chart area") ;
        patch.orientation=area>0 ? 1:-1 ;
        if(area<0) std::reverse(patch.polygon.begin(),patch.polygon.end()) ;
        constexpr long long extent=1LL<<30 ; patch.bounds={{extent,extent,0,0}} ;
        for(const auto &p:patch.polygon) for(unsigned d=0;d<2;++d) {
          patch.bounds[d]=std::min(patch.bounds[d],static_cast<long long>(std::floor(p[d]*extent))) ;
          patch.bounds[d+2]=std::max(patch.bounds[d+2],static_cast<long long>(std::ceil(p[d]*extent))) ;
        }
        // Integrate the saved root map over a triangle fan; edge-only vertices contribute no area.
        patch.mappedArea=true ;
        for(unsigned i=1;i+1<points.size();++i) {
          const auto &a=chart[0],&b=chart[i],&c=chart[i+1] ;
          if((b[0]-a[0])*(c[1]-a[1])==(b[1]-a[1])*(c[0]-a[0])) continue ;
          auto placement=std::make_shared<CornerMap>(Shape::triangle,
            std::vector<Point>{approximate(points[0]),approximate(points[i]),approximate(points[i+1])}) ;
          const auto measure=geometryMeasure(RestrictedMap(geometry,placement)) ;
          require(measure.converged && measure.measure>0,"native face area integration did not converge") ;
          patch.area+=measure.measure ;
        }
      }
    }
    struct EdgeUse { int leaf, depth ; PlanarEdgeUse use ; } ;
    struct PlanarEdgeLess { bool operator()(const PlanarEdgeKey &a,const PlanarEdgeKey &b) const { return a<b ; } } ;
    constexpr long long faceExtent=1LL<<30 ;
    void appendConstruction(std::vector<long long>& key,const SplitConstruction& construction) {
      require(bool(construction.planar),"extruded face requires a planar construction") ;
      const auto& split=*construction.planar ;
      key.insert(key.end(),{static_cast<long long>(construction.version),int(split.shape),int(split.centerMethod),static_cast<long long>(split.points.size())}) ;
      for(const auto& p:split.points) for(double x:p) {
        if(x==0) x=0 ;
        long long bits ; std::memcpy(&bits,&x,sizeof(double)) ; key.push_back(bits) ;
      }
    }

    // A local view of accepted polygons and their exact edge intervals.
    struct ExtrudedView {
      const CellTopology &topology ; const ExtrudedSubdivision &cell ; const PlanarSubdivision &tree ;
      ExtrusionLayout layout ;
      std::map<PlanarEdgeKey,std::vector<EdgeUse>,PlanarEdgeLess> uses ;
      std::map<PlanarEdgeKey,std::set<Fraction>,PlanarEdgeLess> cuts ;
      std::vector<PlanarVertexKey> interior ;
      std::vector<TreeNodeId> leaves ;
      ExtrudedView(const CellTopology &t,const ExtrudedSubdivision &c) : topology(t),cell(c),tree(c.planar),layout(c.shape,c.corners,c.axis),leaves(tree.leafIds()) {
        std::set<PlanarVertexKey> vertices ;
        for(unsigned i=0;i<leaves.size();++i) {
          const auto &node=tree.node(leaves[i]) ;
          for(const auto &use:node.edges) {
            uses[use.edge].push_back({int(i),int(tree.history().path(leaves[i]).size()),use}) ;
            cuts[use.edge].insert(use.begin) ; cuts[use.edge].insert(use.end) ;
          }
          for(const auto &p:node.corners)
            if(p.kind==PlanarVertexKey::Kind::center || (p.kind==PlanarVertexKey::Kind::edge && !p.edge.original)) vertices.insert(p) ;
        }
        interior.assign(vertices.begin(),vertices.end()) ;
      }
      VertexKey edgePoint(unsigned edge, Fraction p, unsigned layer) const {
        const auto &vertices=layer ? layout.upper:layout.lower ; const auto shape=topology.shape() ;
        const auto a=referenceVertex(shape,vertices.at(edge)), b=referenceVertex(shape,vertices.at((edge+1)%vertices.size())) ;
        ReferencePoint q ; for(unsigned d=0;d<3;++d) q[d]=a[d]+p*(b[d]-a[d]) ; return topology.vertexKey(q) ;
      }
      ParametricNodeKey key(const PlanarVertexKey &p,unsigned layer) const {
        if(p.kind==PlanarVertexKey::Kind::corner) return {edgePoint(p.corner,0,layer),-1} ;
        if(p.kind==PlanarVertexKey::Kind::edge && p.edge.original) return {edgePoint(p.edge.index,p.position,layer),-1} ;
        const auto it=std::lower_bound(interior.begin(),interior.end(),p) ;
        require(it!=interior.end() && *it==p,"accepted cap vertex is missing") ;
        return {{EntityKind::face,topology.faceId(layout.caps[layer]),{{0,0,0}}},int(it-interior.begin())} ;
      }
      ParametricNodeKey key(const PlanarEdgeKey &edge,Fraction p,unsigned layer) const { return key(tree.vertexAt(edge,p),layer) ; }
      Point position(unsigned index,unsigned layer) const {
        const auto p=tree.evaluate(interior.at(index)) ; Point out ; out[cell.axis]=layout.planes[layer] ;
        out[(cell.axis+1)%3]=p[0] ; out[(cell.axis+2)%3]=p[1] ; return out ;
      }
      int cap(EntityId face) const {
        for(unsigned i=0;i<2;++i) if(topology.faceId(layout.caps[i])==face) return i ;
        return -1 ;
      }
      unsigned side(EntityId face) const {
        for(unsigned i=0;i<layout.sides.size();++i) if(topology.faceId(layout.sides[i])==face) return i ;
        throw std::invalid_argument("face is not an extrusion side") ;
      }
      bool reversed(unsigned edge) const { return topology.vertexId(layout.lower[edge])>topology.vertexId(layout.lower[(edge+1)%layout.lower.size()]) ; }
      Fraction coordinate(unsigned edge,Fraction p) const { return reversed(edge) ? Fraction(1)-p:p ; }
      std::vector<EdgeUse> sideUses(unsigned edge) const {
        const auto found=uses.find(PlanarEdgeKey{true,{},edge}) ; require(found!=uses.end(),"extrusion side has no leaves") ; return found->second ;
      }
      FacePatch sidePatch(std::vector<long long> surface,Fraction lo,Fraction hi,int orientation,double length) const {
        FacePatch p ; p.surface=std::move(surface) ;
        p.bounds={{checkedExactInteger<long long>((lo*Fraction(faceExtent)).numerator()),0,
                   checkedExactInteger<long long>((hi*Fraction(faceExtent)).numerator()),faceExtent}} ;
        require((lo*Fraction(faceExtent)).denominator()==1 && (hi*Fraction(faceExtent)).denominator()==1,"face interval exceeds coordinate depth") ;
        p.orientation=orientation ; p.area=length*(hi-lo).value()*(layout.planes[1]-layout.planes[0]) ; return p ;
      }
      FacePatch spokePatch(const PlanarEdgeKey& edge,Fraction lo,Fraction hi,int orientation) const {
        std::vector<long long> surface{3,static_cast<long long>(tree.history().rootKey()),static_cast<long long>(tree.history().geometryRevision()),static_cast<long long>(edge.owner.parent.lineage.size())} ;
        for(const auto& step:edge.owner.parent.lineage) { appendConstruction(surface,*step.construction) ; surface.push_back(step.child) ; }
        appendConstruction(surface,*edge.owner.construction) ; surface.push_back(edge.index) ;
        const auto a=tree.evaluate(tree.vertexAt(edge,0)),b=tree.evaluate(tree.vertexAt(edge,1)) ;
        return sidePatch(std::move(surface),lo,hi,orientation,std::hypot(b[0]-a[0],b[1]-a[1])) ;
      }
      void addBoundaryPoints(const std::vector<ParametricNodeKey> &points) {
        for(unsigned edge=0;edge<layout.sides.size();++edge) {
          for(unsigned layer=0;layer<2;++layer) {
            const auto midpoint=edgePoint(edge,Fraction(1,2),layer) ; const auto &vertices=layer ? layout.upper:layout.lower ;
            for(const auto &p:points) if(p.local<0 && p.reference.kind==EntityKind::edge && p.reference.entity==midpoint.entity) {
              Fraction local=p.reference.point[0] ;
              if(topology.vertexId(vertices[edge])>topology.vertexId(vertices[(edge+1)%vertices.size()])) local=Fraction(1)-local ;
              cuts[PlanarEdgeKey{true,{},edge}].insert(local) ;
            }
          }
        }
      }
      std::vector<ParametricOutputFace> caps(unsigned layer) const {
        std::vector<ParametricOutputFace> result ;
        Point low=layout.polygon.front(),high=low ;
        for(const auto& p:layout.polygon) for(unsigned d=0;d<2;++d) { low[d]=std::min(low[d],p[d]) ; high[d]=std::max(high[d],p[d]) ; }
        for(unsigned i=0;i<leaves.size();++i) {
          ParametricOutputFace f{{},int(i),-1} ;
          for(const auto &use:tree.node(leaves[i]).edges) {
            const auto &all=cuts.at(use.edge) ;
            std::vector<Fraction> part(all.lower_bound(std::min(use.begin,use.end)),all.upper_bound(std::max(use.begin,use.end))) ;
            if(use.end<use.begin) std::reverse(part.begin(),part.end()) ;
            part.pop_back() ;
            for(const auto &p:part) f.vertices.push_back(key(use.edge,p,layer)) ;
          }
          if(!layer) std::reverse(f.vertices.begin(),f.vertices.end()) ;
          f.patch.surface={2,static_cast<long long>(topology.faceId(layout.caps[layer]))} ; f.patch.orientation=layer ? 1:-1 ;
          f.patch.bounds={{faceExtent,faceExtent,0,0}} ;
          const auto& node=tree.node(leaves[i]) ; const auto origin=tree.evaluate(node.corners[0]) ;
          long double area=0 ;
          for(unsigned j=0;j<node.corners.size();++j) {
            const auto p=tree.evaluate(node.corners[j]),q=tree.evaluate(node.corners[(j+1)%node.corners.size()]) ;
            f.patch.polygon.push_back({{(p[0]-low[0])/(high[0]-low[0]),(p[1]-low[1])/(high[1]-low[1])}}) ;
            // Outward bounds prune the correspondence search; polygons determine overlap.
            for(unsigned d=0;d<2;++d) {
              const double scaled=f.patch.polygon.back()[d]*faceExtent ;
              f.patch.bounds[d]=std::min(f.patch.bounds[d],static_cast<long long>(std::floor(scaled))) ;
              f.patch.bounds[d+2]=std::max(f.patch.bounds[d+2],static_cast<long long>(std::ceil(scaled))) ;
            }
            area+=(static_cast<long double>(p[0])-origin[0])*(static_cast<long double>(q[1])-origin[1])-(static_cast<long double>(p[1])-origin[1])*(static_cast<long double>(q[0])-origin[0]) ;
          }
          f.patch.area=double(area/2) ;
          result.push_back(std::move(f)) ;
        }
        return result ;
      }
      std::vector<ParametricOutputFace> internal() const {
        std::vector<ParametricOutputFace> result ;
        for(const auto &entry:uses) if(!entry.first.original) {
          const auto &all=cuts.at(entry.first) ;
          for(auto a=all.begin(),b=std::next(a);b!=all.end();++a,++b) {
            std::vector<const EdgeUse *> incident ;
            for(const auto &use:entry.second)
              if(std::min(use.use.begin,use.use.end)<=*a && std::max(use.use.begin,use.use.end)>=*b) incident.push_back(&use) ;
            require(incident.size()==2,"internal extrusion edge must have two incident cells") ;
            require((incident[0]->use.begin<incident[0]->use.end)!=(incident[1]->use.begin<incident[1]->use.end),"internal extrusion edges have equal orientations") ;
            const bool forward=incident[0]->use.begin<incident[0]->use.end ; const auto x=forward ? *a:*b,y=forward ? *b:*a ;
            result.push_back({{key(entry.first,x,0),key(entry.first,y,0),key(entry.first,y,1),key(entry.first,x,1)},incident[0]->leaf,incident[1]->leaf}) ;
            result.back().patch=spokePatch(entry.first,*a,*b,forward ? 1:-1) ;
          }
        }
        return result ;
      }
    } ;
    std::vector<int> depths(const ParametricCell &cell) {
      if(cell.source()) return {0} ;
      std::vector<int> out ; for(auto id:cell.leafIds()) out.push_back(cell.path(id).size()) ; return out ;
    }
    void gradePair(std::set<int> &out,int a,int b,const std::vector<int> &depth,unsigned difference) {
      if(depth[a]>depth[b] && unsigned(depth[a]-depth[b])>difference) out.insert(b) ;
      if(depth[b]>depth[a] && unsigned(depth[b]-depth[a])>difference) out.insert(a) ;
    }
    ParametricOutputFace sourceFace(EntityId face,const ParametricSource &source,
        const std::vector<ParametricNodeKey> &points,bool neighbor,bool patches) {
      const auto &t=source.topology ;
      unsigned f=0 ; while(f<t.faceCount() && t.faceId(f)!=face) ++f ;
      require(f<t.faceCount(),"face is not incident to the original cell") ;
      const auto loop=t.faceVertices(f) ; const auto shape=t.faceShape(f) ; const bool reverse=t.faceOrientation(f)>=loop.size() ;
      std::vector<unsigned> vertices(loop.size()) ; std::vector<Point> positions(loop.size()) ;
      for(unsigned i=0;i<loop.size();++i) {
        const auto j=orientedVertex(unsigned(loop.size()),t.faceOrientation(f),i) ; vertices[j]=loop[i] ; positions[j]=source.corners[loop[i]] ;
      }
      ParametricOutputFace out{{},0,neighbor ? 0:-1} ;
      for(unsigned i=0;i<vertices.size();++i) {
        const auto a=vertices[i],d=vertices[(i+1)%vertices.size()] ;
        out.vertices.push_back({{EntityKind::vertex,t.vertexId(a),{{0,0,0}}},-1}) ;
        const auto ends=std::minmax(a,d) ; unsigned edge=0 ;
        while(edge<t.edgeCount()) {
          const auto vertices=t.edgeVertices(edge) ; if(std::minmax(vertices[0],vertices[1])==ends) break ; ++edge ;
        }
        require(edge<t.edgeCount(),"source face edge is missing") ;
        std::map<Fraction,ParametricNodeKey> cuts ;
        for(const auto &key:points) if(key.reference.kind==EntityKind::edge && key.reference.entity==t.edgeId(edge)) {
          const auto q=key.reference.point[0] ; require(key.local<0 && q>0 && q<1,"invalid source edge point") ;
          cuts.emplace(t.vertexId(a)<t.vertexId(d) ? q:Fraction(1)-q,key) ;
        }
        for(const auto &cut:cuts) out.vertices.push_back(cut.second) ;
      }
      if(reverse) std::reverse(out.vertices.begin(),out.vertices.end()) ;
      if(patches) {
        auto &p=out.patch ; p.surface={2,static_cast<long long>(face)} ; p.sector=loop.size()==4 ? 0:-int(loop.size()) ;
        p.bounds={{0,0,faceExtent,faceExtent}} ; p.orientation=reverse ? -1:1 ; p.mappedArea=true ;
        const auto measure=shape==Shape::polygon ? geometryMeasure(PolygonMap(positions),MeasureOptions{}):
          geometryMeasure(CornerMap(shape,positions),MeasureOptions{}) ;
        require(measure.converged && measure.measure>0,"corner face area integration did not converge") ; p.area=measure.measure ;
      }
      return out ;
    }
  }

  ParametricCell ParametricCell::original(const CellTopology &topology,const std::vector<Point> &corners,HistoryLimits limits) {
    require(limits.maxNodes>0,"parametric source requires a positive node limit") ;
    if(topology.isStandard()) {
      const auto shape=topology.shape() ;
      require(shape==Shape::hex || shape==Shape::prism || shape==Shape::tet || shape==Shape::pyramid,"unsupported corner source cell") ;
      const CornerMap geometry(shape,corners) ;
      require(geometryQuality(geometry).status==GeometryQualityStatus::sampled,"invalid or unresolved parametric root geometry") ;
    } else {
      const ParametricBoundary checked(topology,corners,limits) ;
    }
    return ParametricCell(std::make_shared<const ParametricSource>(ParametricSource{topology,corners,std::move(limits)})) ;
  }
  ParametricCell::ParametricCell(const ParametricCell &other) : source_(other.source_),
    cell_(other.cell_ ? std::make_unique<Subdivision>(*other.cell_):nullptr) {}
  ParametricCell &ParametricCell::operator=(const ParametricCell &other) {
    if(this!=&other) { ParametricCell copy(other) ; source_.swap(copy.source_) ; cell_.swap(copy.cell_) ; }
    return *this ;
  }
  void ParametricCell::expand() {
    if(const auto *p=source()) {
      auto cell=!p->topology.isStandard() || p->topology.shape()==Shape::pyramid
        ? std::make_unique<Subdivision>(ParametricBoundary(p->topology,p->corners,p->limits))
        : std::make_unique<Subdivision>(ParametricCorner(p->topology,p->corners,p->limits)) ;
      cell_=std::move(cell) ; source_.reset() ;
    }
  }
  std::size_t ParametricCell::leafCount() const { return source() ? 1:leafIds().size() ; }
  std::vector<TreeNodeId> ParametricCell::leafIds() const { return source() ? std::vector<TreeNodeId>{0}:boundary() ? boundary()->leafIds():history().leafIds() ; }
  const HistoryLimits &ParametricCell::limits() const { return source() ? source()->limits:boundary() ? boundary()->limits():history().limits() ; }
  std::uint64_t ParametricCell::rootKey() const { return source() ? source()->topology.id():boundary() ? boundary()->rootKey():history().rootKey() ; }
  CellPath ParametricCell::path(TreeNodeId id) const {
    if(source()) { require(id==0,"original cell has no child") ; return {} ; }
    return boundary() ? boundary()->path(id):history().path(id) ;
  }
  TreeNodeId ParametricCell::parent(TreeNodeId id) const {
    if(source()) { require(id==0,"original cell has no child") ; return noTreeNode ; }
    return boundary() ? boundary()->node(id).parent:history().node(id).parent ;
  }
  NodeRange ParametricCell::children(TreeNodeId id) const {
    if(source()) { require(id==0,"original cell has no child") ; return {} ; }
    return boundary() ? (boundary()->node(id).refined ? boundary()->children(id):NodeRange()):history().children(id) ;
  }
  bool ParametricCell::active(TreeNodeId id) const { return source() ? id==0:boundary() ? boundary()->active(id):history().active(id) ; }
  bool ParametricCell::leaf(TreeNodeId id) const {
    return active(id) && (source() || (boundary() ? !boundary()->node(id).refined:history().node(id).isLeaf())) ;
  }
  std::size_t ParametricCell::splitCount() const { return source() ? 0:boundary() ? boundary()->splits().size():history().splits().size() ; }
  TreeNodeId ParametricCell::splitParent(std::size_t i) const { return boundary() ? boundary()->splits().at(i).parent:history().splits().at(i).parent ; }
  NodeRange ParametricCell::splitChildren(std::size_t i) const { return boundary() ? boundary()->splits().at(i).children:history().splits().at(i).children ; }
  const CellTopology &ParametricCell::cornerTopology() const {
    if(source()) return source()->topology ;
    require(hasCornerGeometry(),"cell has no corner topology") ; return corner() ? corner()->topology:boundary()->topology ;
  }
  const CornerCell &ParametricCell::cornerCell(TreeNodeId id) const {
    require(hasCornerGeometry(),"cell has no shared corner geometry") ; return corner() ? corner()->cells.at(id):boundary()->node(id).cell ;
  }
  const CornerBoundary &ParametricCell::cornerGeometry() const {
    if(corner()) return *corner() ;
    require(boundary(),"cell has no shared corner geometry") ; return *boundary() ;
  }
  CornerBoundary &ParametricCell::cornerGeometry() {
    if(auto *p=std::get_if<ParametricCorner>(cell_.get())) return *p ;
    auto *p=std::get_if<ParametricBoundary>(cell_.get()) ; require(p,"cell has no shared corner geometry") ; return *p ;
  }
  BoundaryId ParametricCell::cornerFaceId(EntityId face) const {
    const auto &root=cornerCell(0) ;
    for(unsigned f=0;f<root.faces.size();++f) if(cornerTopology().faceId(f)==face) return root.faces[f].surface ;
    throw std::invalid_argument("face is not incident to the original corner cell") ;
  }
  const RefinementHistory &ParametricCell::history() const {
    require(!source() && !boundary(),"cell has no standard-element history") ;
    if(const auto *p=native()) return p->history() ;
    if(const auto *p=corner()) return p->subdivision.history() ;
    return extruded()->planar.history() ;
  }
  CenterMethod parametricCenterMethod(const std::string &name) {
    if(name=="reference") return CenterMethod::reference ;
    if(name=="weighted") return CenterMethod::weighted ;
    if(name=="integrated") return CenterMethod::integrated ;
    throw std::invalid_argument("unknown parametric center method") ;
  }
  void ParametricCell::refine(TreeNodeId leaf,const SplitOperation &request,CenterMethod method) {
    if(auto *p=std::get_if<AffineSubdivision>(cell_.get())) { require(method==CenterMethod::reference,"native cuts support only reference centers") ; p->refine(leaf,request) ; }
    else throw std::invalid_argument("corner and extruded splitting require an isotropic request") ;
  }

  void ParametricCell::refineIsotropic(TreeNodeId leaf,CenterMethod method,int axis) {
    if(source()) {
      auto candidate=*this ; candidate.expand() ; candidate.refineIsotropic(leaf,method,axis) ;
      *this=std::move(candidate) ; return ;
    }
    if(const auto *p=native()) refine(leaf,axis<0 ? isotropicSplit(p->node(leaf).cell.shape):parametricGradingSplit(*p,axis),method) ;
    else if(auto *p=std::get_if<ParametricCorner>(cell_.get())) { require(axis<0,"corner splitting requires full 3D refinement") ; p->refine(leaf,method) ; }
    else if(auto *p=std::get_if<ParametricBoundary>(cell_.get())) { require(axis<0,"corner splitting requires full 3D refinement") ; p->refine(leaf,method) ; }
    else std::get<ExtrudedSubdivision>(*cell_).planar.splitLeaf(leaf,method) ;
  }
  void ParametricCell::refineForBalance(TreeNodeId leaf,CenterMethod method,int axis,unsigned directions) {
    require(directions<8,"invalid balancing directions") ;
    const auto *p=native() ;
    if(!p || directions==0) { refineIsotropic(leaf,method,axis) ; return ; }
    auto request=axis<0 ? isotropicSplit(p->node(leaf).cell.shape):parametricGradingSplit(*p,axis) ;
    if(auto *tensor=std::get_if<TensorSplit>(&request)) {
      unsigned allowed=0 ; for(unsigned d=0;d<3;++d) if(!tensor->cuts[d].empty()) allowed|=1u<<d ;
      if(directions&allowed) for(unsigned d=0;d<3;++d) if(!(directions&(1u<<d))) tensor->cuts[d].clear() ;
    } else if(auto *prism=std::get_if<PrismSplit>(&request)) {
      const unsigned allowed=(prism->base.kind!=TriangleSplit::Kind::none ? 3:0) | (!prism->axialCuts.empty() ? 4:0) ;
      if(directions&allowed) {
        if(!(directions&3)) prism->base=TriangleSplit() ;
        if(!(directions&4)) prism->axialCuts.clear() ;
      }
    }
    refine(leaf,request,method) ;
  }
  void ParametricCell::coarsen(TreeNodeId parent) {
    require(!source(),"original cell has no children to coarsen") ;
    if(auto *p=std::get_if<AffineSubdivision>(cell_.get())) p->coarsen(parent) ;
    else if(auto *p=std::get_if<ParametricCorner>(cell_.get())) p->subdivision.coarsen(parent) ;
    else if(auto *p=std::get_if<ParametricBoundary>(cell_.get())) p->coarsen(parent) ;
    else std::get<ExtrudedSubdivision>(*cell_).planar.coarsen(parent) ;
  }
  bool ParametricCell::restoreSplit(TreeNodeId parent,const RefinementHistory &previous) {
    TreeNodeId old ;
    try { old=previous.nodeId(history().key(parent)) ; } catch(const std::out_of_range &) { return false ; }
    if(!previous.active(old) || previous.node(old).isLeaf()) return false ;
    const auto &split=*previous.splits().at(previous.node(old).activeSplit).construction ;
    if(auto *p=std::get_if<AffineSubdivision>(cell_.get())) p->reactivate(parent,split) ;
    else if(auto *p=std::get_if<ParametricCorner>(cell_.get())) p->subdivision.reactivate(parent,split) ;
    else std::get<ExtrudedSubdivision>(*cell_).planar.reactivate(parent,split) ;
    return true ;
  }
  bool ParametricCell::restoreSplit(TreeNodeId parent,const ParametricCell &previous) {
    if(previous.source()) return false ;
    if(!boundary()) { require(!previous.boundary(),"previous split family differs") ; return restoreSplit(parent,previous.history()) ; }
    require(previous.boundary(),"previous split family differs") ;
    TreeNodeId old=0 ;
    for(auto child:path(parent)) {
      const auto children=previous.boundary()->children(old) ; if(child>=children.size()) return false ; old=children[child] ;
    }
    if(!previous.active(old) || previous.leaf(old)) return false ;
    require(boundary()->splitRecord(parent)==previous.boundary()->splitRecord(old),"previous accepted split differs") ;
    std::get<ParametricBoundary>(*cell_).refine(parent) ; return true ;
  }
  void ParametricCell::validate(int axis) const {
    if(source()) { require(axis<0,"corner splitting requires full 3D refinement") ; return ; }
    if(boundary()) { require(axis<0,"corner splitting requires full 3D refinement") ; return ; }
    if(const auto *p=native()) parametricGradingSplit(*p,axis) ;
    else if(corner()) require(axis<0,"corner splitting requires full 3D refinement") ;
    else require(axis==extruded()->axis,"saved extrusion axis differs") ;
    for(auto id:history().leafIds())
      require(geometryQuality(*history().geometry().map(history().node(id).map)).status==GeometryQualityStatus::sampled,"invalid or unresolved parametric leaf geometry") ;
  }
  bool ParametricNodeKey::operator<(const ParametricNodeKey &o) const { return reference==o.reference ? local<o.local:reference<o.reference ; }
  bool ParametricNodeKey::operator==(const ParametricNodeKey &o) const { return reference==o.reference && local==o.local ; }

  ParametricRecord encodeParametricCell(const ParametricCell &cell) {
    if(const auto *p=cell.source()) {
      detail::RecordWriter out(p->limits.exact.maxBytes,p->limits.exact) ; out.u32(sourceMagic) ; out.u32(1) ;
      const auto topology=encodeParametricTopology(p->topology) ; out.u64(topology.size()) ; for(int i:topology) out.u32(i) ;
      for(const auto &point:p->corners) for(double x:point) out.real(x) ;
      return out.take() ;
    }
    if(const auto *p=cell.native()) return encodeParametricHistory(*p) ;
    if(const auto *p=cell.boundary()) return p->record() ;
    if(const auto *p=cell.corner()) return p->topology.shape()==Shape::hex ? ParametricHex(*p).record() : p->record() ;
    const auto &p=*cell.extruded() ; const auto limits=recordLimits(cell.history().limits()) ;
    detail::RecordWriter out(limits.maxBytes,limits.geometry.exact) ; out.u32(extrusionMagic) ; out.u32(1) ;
    out.u32(unsigned(p.shape)) ; out.u32(p.axis) ; out.u64(p.corners.size()) ;
    for(const auto &point:p.corners) for(double x:point) out.real(x) ;
    const auto history=encodeRefinement(cell.history(),limits) ; out.u64(history.size()) ; out.append(history) ; return out.take() ;
  }
  ParametricCell decodeParametricCell(const ParametricRecord &record,HistoryLimits limits) {
    if(ParametricHex::saved(record)) return ParametricHex::decode(record,limits).expand() ;
    detail::RecordReader in(record,limits.exact.maxBytes,limits.exact) ;
    const auto magic=in.u32() ;
    if(magic==sourceMagic) {
      require(in.u32()==1,"unknown source cell version") ;
      std::vector<int> t(in.count(limits.exact.maxBytes/4,4)) ; for(auto &i:t) i=in.u32() ; const auto topology=decodeParametricTopology(t) ;
      std::vector<Point> corners(topology.vertexCount()) ; for(auto &p:corners) for(auto &x:p) x=in.real() ;
      in.finish() ; return ParametricCell::original(topology,corners,limits) ;
    }
    if(magic==0x434f524e) return ParametricCorner::restore(record,limits) ;
    if(magic==0x42434f52) return ParametricBoundary::restore(record,limits) ;
    if(magic!=extrusionMagic) return decodeParametricHistory(record,limits) ;
    require(in.u32()==1,"unknown extruded history version") ; const auto shape=Shape(in.u32()) ; const int axis=in.u32() ;
    require(shape==Shape::hex || shape==Shape::prism,"extruded history needs a hex or prism") ;
    std::vector<Point> corners(in.count(8,24)) ; require(corners.size()==referenceElement(shape).vertexCount,"invalid extruded corner count") ;
    for(auto &p:corners) for(double &x:p) x=in.real() ;
    ParametricRecord bytes(in.count(limits.exact.maxBytes)) ; for(auto &b:bytes) b=in.byte() ; in.finish() ;
    PlanarSubdivision tree(decodeRefinement(bytes,recordLimits(limits))) ; const ExtrusionLayout layout(shape,corners,axis) ;
    const auto &root=tree.history() ; const auto *geometry=std::get_if<CornerRecord>(&root.geometry().record(root.node(0).map)) ;
    require(geometry && geometry->corners==layout.polygon,"saved planar geometry differs from its extrusion") ;
    return ExtrudedSubdivision{shape,std::move(corners),axis,std::move(tree)} ;
  }
  ParametricCell makeParametricCell(const CellTopology &topology,const std::vector<Point> &corners,
      const ParametricRecord &record,HistoryLimits limits,const std::string &requested,int axis) {
    const auto family=parametricCellFamily(topology,requested) ;
    if(requested=="mixed" && axis>=0) throw std::invalid_argument("mixed splitting requires full 3D refinement") ;
    require(family=="native" || family=="extruded" || family=="corner","unknown parametric split family") ;
    if(family=="native") return makeParametricHistory(topology,corners,record,limits) ;
    if(family=="corner") {
      require(axis<0,"corner splitting requires full 3D refinement") ;
      if(record.empty()) {
        if(!topology.isStandard() || topology.shape()==Shape::pyramid) return ParametricBoundary(topology,corners,limits) ;
        return ParametricCorner(topology,corners,limits) ;
      }
      auto saved=decodeParametricCell(record,limits) ; require(saved.hasCornerGeometry(),"saved history has a different split family") ;
      require(saved.cornerTopology()==topology,"saved corner source topology differs") ;
      if(saved.boundary()) {
        require(saved.boundary()->sourceCorners()==corners,"saved boundary source geometry differs") ; return saved ;
      }
      const auto &p=*saved.corner() ; const auto &h=p.subdivision.history() ;
      const auto *root=std::get_if<CornerRecord>(&h.geometry().record(h.node(0).map)) ;
      require(root && root->corners==corners,"saved corner source geometry differs") ;
      return saved ;
    }
    const ExtrusionLayout layout(topology.shape(),corners,axis) ;
    if(!record.empty()) {
      auto saved=decodeParametricCell(record,limits) ; const auto *p=saved.extruded() ;
      require(p && p->shape==topology.shape() && p->corners==corners && p->axis==axis && saved.history().rootKey()==topology.id() && saved.history().geometryRevision()==0,
              "saved extrusion differs from original geometry") ; return saved ;
    }
    return ExtrudedSubdivision{topology.shape(),corners,axis,PlanarSubdivision(layout.polygon,topology.id(),0,limits)} ;
  }

  ParametricRecord initialParametricCell(const CellTopology &topology,const std::vector<Point> &corners,
      const ParametricRecord &record,HistoryLimits limits,const std::string &family,int axis) {
    const auto cell=makeParametricCell(topology,corners,record,limits,family,axis) ;
    return record.empty() ? encodeParametricCell(cell):record ;
  }

  bool GeometricRefinement::refine(ParametricCell &cell,CenterMethod method,int axis,bool automatic,double aspectRatio) const {
    if(!select) return false ;
    const double cutoff=minEdgeLength ;
    require(std::isfinite(cutoff) && cutoff>=0,"geometric refinement requires a nonnegative finite minEdgeLength") ;
    if(automatic) require(axis==-1 && std::isfinite(aspectRatio) && aspectRatio>1,
      "automatic geometric refinement requires 3D cells and an aspect ratio greater than one") ;
    auto pending=cell.leafIds() ; bool changed=false ;
    while(!pending.empty()) {
      const auto leaf=pending.back() ; pending.pop_back() ;
      std::vector<Point> points ; double shortest=std::numeric_limits<double>::max() ;
      const auto edge=[&](const Point &a,const Point &b) {
        double full=0,plane=0 ;
        for(int d=0;d<3;++d) { const double square=(b[d]-a[d])*(b[d]-a[d]) ; full+=square ; if(d!=axis) plane+=square ; }
        if(axis>=0 && plane<=1e-24*full) return ;
        shortest=std::min(shortest,std::sqrt(plane)) ;
      } ;
      if(const auto *source=cell.source()) {
        points=source->corners ;
        for(unsigned e=0;e<source->topology.edgeCount();++e) {
          const auto ends=source->topology.edgeVertices(e) ; edge(points[ends[0]],points[ends[1]]) ;
        }
      } else if(cell.hasCornerGeometry()) {
        const auto &geometry=cell.cornerGeometry() ; const auto &part=cell.cornerCell(leaf) ;
        for(auto vertex:part.vertices) points.push_back(geometry.boundaries.points[vertex]) ;
        for(const auto &ends:geometry.cellEdges(part)) edge(geometry.boundaries.points[ends.first.first],geometry.boundaries.points[ends.first.second]) ;
      } else if(const auto *extruded=cell.extruded()) {
        const ExtrusionLayout layout(extruded->shape,extruded->corners,extruded->axis) ;
        const auto &corners=extruded->planar.node(leaf).corners ;
        for(unsigned layer=0;layer<2;++layer) for(const auto &corner:corners) {
          const auto p=extruded->planar.evaluate(corner) ; Point q{} ; q[extruded->axis]=layout.planes[layer] ;
          q[(extruded->axis+1)%3]=p[0] ; q[(extruded->axis+2)%3]=p[1] ; points.push_back(q) ;
        }
        for(unsigned v=0;v<corners.size();++v) edge(points[v],points[(v+1)%corners.size()]) ;
      } else {
        const auto &history=cell.history() ; const auto map=history.geometry().map(history.node(leaf).map) ;
        const auto &element=referenceElement(map->shape()) ;
        for(unsigned v=0;v<element.vertexCount;++v) points.push_back(map->evaluate(element.vertices[v]).position) ;
        for(unsigned e=0;e<element.edgeCount;++e) edge(points[element.edges[e][0]],points[element.edges[e][1]]) ;
      }
      if(!select(points,shortest)) continue ;
      if(automatic && cell.native()) {
        const auto &history=cell.history() ; const auto map=history.geometry().map(history.node(leaf).map) ;
        const auto split=automaticSplit(*map,aspectRatio,cutoff) ; if(!split) continue ; cell.refine(leaf,*split,method) ;
      } else {
        if(shortest<2*cutoff) continue ;
        cell.refineIsotropic(leaf,method,axis) ;
      }
      changed=true ; for(auto child:cell.children(leaf)) pending.push_back(child) ;
    }
    return changed ;
  }

  GeometricRefinement GeometricRefinement::region(const std::function<bool(const Point &)> &selected,double cutoff) {
    require(std::isfinite(cutoff) && cutoff>0,"XML refinement requires a positive finite minEdgeLength") ;
    return {[selected](const std::vector<Point> &points,double) { return std::any_of(points.begin(),points.end(),selected) ; },cutoff} ;
  }
  GeometricRefinement GeometricRefinement::sources(const std::vector<SpacingSource> &sources,double cutoff) {
    require(std::isfinite(cutoff) && cutoff>=0,"geometric refinement requires a nonnegative finite minEdgeLength") ;
    return {[&sources](const std::vector<Point> &points,double shortest) { return refineForSpacing(points,sources,shortest) ; },cutoff} ;
  }
  bool refineParametricRegion(ParametricCell &cell,const std::function<bool(const Point &)> &selected,double cutoff,
      CenterMethod method,int axis,bool automatic,double aspectRatio) {
    return GeometricRefinement::region(selected,cutoff).refine(cell,method,axis,automatic,aspectRatio) ;
  }
  bool refineParametricSources(ParametricCell &cell,const std::vector<SpacingSource> &sources,double cutoff,
      CenterMethod method,int axis,bool automatic,double aspectRatio) {
    return GeometricRefinement::sources(sources,cutoff).refine(cell,method,axis,automatic,aspectRatio) ;
  }

  std::vector<std::vector<ParametricNodeKey>> parametricLeafVertices(const CellTopology &t,const ParametricCell &c) {
    if(c.source()) {
      std::vector<ParametricNodeKey> vertices ;
      for(unsigned v=0;v<t.vertexCount();++v) vertices.push_back({{EntityKind::vertex,t.vertexId(v),{{0,0,0}}},-1}) ;
      return {vertices} ;
    }
    if(c.hasCornerGeometry()) return cornerLeafVertices(c) ;
    std::vector<std::vector<ParametricNodeKey>> result ;
    if(const auto *native=c.native()) {
      for(auto leaf:c.leafIds()) {
        const auto &child=native->node(leaf).cell ; result.emplace_back() ;
        for(unsigned v=0;v<referenceElement(child.shape).vertexCount;++v)
          result.back().push_back({t.vertexKey(child.placement.apply(referenceVertex(child.shape,v))),-1}) ;
      }
    } else {
      const ExtrudedView view(t,*c.extruded()) ;
      for(auto leaf:view.leaves) {
        result.emplace_back() ;
        for(unsigned layer=0;layer<2;++layer) for(const auto &vertex:view.tree.node(leaf).corners)
          result.back().push_back(view.key(vertex,layer)) ;
      }
    }
    return result ;
  }
  std::vector<ParametricNodeKey> parametricEdgeNodes(EntityId edge,const ParametricTrace &trace) {
    std::set<Fraction> coordinates ;
    for(const auto &interval:trace) {
      require(interval.size()==8 && interval[0]==int(edge),"invalid parametric edge interval") ;
      coordinates.insert(Fraction(interval[4],interval[5])) ; coordinates.insert(Fraction(interval[6],interval[7])) ;
    }
    std::vector<ParametricNodeKey> keys ;
    for(const auto &p:coordinates) if(p>Fraction(0) && p<Fraction(1))
      keys.push_back({{EntityKind::edge,edge,{{p,0,0}}},-1}) ;
    return keys ;
  }
  std::vector<ParametricNodeKey> parametricCellNodes(const CellTopology &t,const ParametricCell &c) {
    if(c.source()) return {} ;
    if(const auto *p=c.native()) return outputKeys(parametricCellNodes(t,*p)) ;
    if(c.hasCornerGeometry()) return cornerNodes(c,EntityKind::cell,t.id()) ;
    return {} ;
  }
  std::vector<ParametricNodeKey> parametricFaceNodes(EntityId face,const CellTopology &t,const ParametricCell &c,const CellTopology *r,const ParametricCell *rc) {
    if(c.source() && (!rc || rc->source())) return {} ;
    if(c.source() || (rc && rc->source())) {
      auto left=c ; left.expand() ; if(!rc) return parametricFaceNodes(face,t,left) ;
      auto right=*rc ; right.expand() ; return parametricFaceNodes(face,t,left,r,&right) ;
    }
    if(mixedFace(c,rc)) return mixedFaceNodes(face,t,c,*r,*rc) ;
    if(const auto *p=c.native()) { require(!rc || rc->native(),"mixed native and center-based boundaries are unsupported") ; return outputKeys(parametricFaceNodes(face,t,*p,r,rc ? rc->native():nullptr)) ; }
    if(c.hasCornerGeometry()) {
      require(!rc || rc->hasCornerGeometry(),"mixed corner and other split families are unsupported") ;
      return cornerNodes(c,EntityKind::face,face,rc) ;
    }
    require(!rc || rc->extruded(),"mixed native and center-based boundaries are unsupported") ;
    const ExtrudedView v(t,*c.extruded()) ; const int layer=v.cap(face) ;
    if(layer<0) return {} ;
    require(!rc,"center-based extrusion end faces must be boundary faces") ;
    std::vector<ParametricNodeKey> keys ; for(const auto &p:v.interior) keys.push_back(v.key(p,layer)) ; return keys ;
  }
  std::vector<Point> parametricFacePositions(EntityId face,const CellTopology &t,const ParametricCell &c,const CellTopology *r,const ParametricCell *rc) {
    if(c.source() && (!rc || rc->source())) return {} ;
    if(c.source() || (rc && rc->source())) {
      auto left=c ; left.expand() ; if(!rc) return parametricFacePositions(face,t,left) ;
      auto right=*rc ; right.expand() ; return parametricFacePositions(face,t,left,r,&right) ;
    }
    if(mixedFace(c,rc)) return mixedFacePositions(face,t,c,*r,*rc) ;
    const auto keys=parametricFaceNodes(face,t,c,r,rc) ; std::vector<Point> points ;
    if(c.hasCornerGeometry()) {
      auto merged=c ; if(rc) merged.mergeFace(face,parametricFacePlan(face,*rc)) ; return cornerPositions(merged,keys) ;
    }
    if(const auto *p=c.extruded()) {
      const ExtrudedView v(t,*p) ; const int layer=v.cap(face) ; for(const auto &key:keys) points.push_back(v.position(key.local,layer)) ;
    } else {
      const auto geometry=c.history().geometry().map(c.history().node(0).map) ;
      for(const auto &key:keys) { ReferencePoint q ; require(t.locateVertex(key.reference,q),"face point is outside source topology") ; points.push_back(geometry->evaluate(approximate(q)).position) ; }
    }
    return points ;
  }
  std::vector<ParametricOutputFace> parametricCellFaces(const CellTopology &t,const ParametricCell &c,const std::vector<ParametricNodeKey> &points,bool patches) {
    if(c.source()) return {} ;
    if(const auto *p=c.native()) {
      auto faces=outputFaces(parametricCellFaces(t,*p,referenceKeys(points))) ;
      if(patches) nativeFacePatches(faces,t,*p) ;
      return faces ;
    }
    if(c.hasCornerGeometry()) {
      auto faces=cornerCellFaces(c,points,patches) ; mixedCornerFaces(c,points,faces) ; return faces ;
    }
    return ExtrudedView(t,*c.extruded()).internal() ;
  }
  std::vector<ParametricOutputFace> parametricFaceFaces(EntityId face,const CellTopology &t,const ParametricCell &c,const std::vector<ParametricNodeKey> &points,
      const CellTopology *r,const ParametricCell *rc,bool patches) {
    require(bool(r)==bool(rc),"parametric neighbor needs topology and history") ;
    if(c.source() && (!rc || rc->source())) {
      auto out=sourceFace(face,*c.source(),points,bool(rc),patches) ;
      if(rc) {
        const auto other=sourceFace(face,*rc->source(),points,true,false) ;
        auto a=out.vertices,b=other.vertices ; std::reverse(b.begin(),b.end()) ;
        const auto first=std::find(b.begin(),b.end(),a.front()) ;
        require(first!=b.end(),"shared source face vertices differ") ; std::rotate(b.begin(),first,b.end()) ;
        require(a==b,"shared source face connectivity or orientation differs") ;
      }
      return {std::move(out)} ;
    }
    if(c.source() || (rc && rc->source())) {
      auto left=c ; left.expand() ; if(!rc) return parametricFaceFaces(face,t,left,points,nullptr,nullptr,patches) ;
      auto right=*rc ; right.expand() ; return parametricFaceFaces(face,t,left,points,r,&right,patches) ;
    }
    if(mixedFace(c,rc)) return mixedFaceFaces(face,t,c,*r,*rc,points,patches) ;
    if(const auto *p=c.native()) {
      require(!rc || rc->native(),"mixed native and center-based boundaries are unsupported") ;
      auto faces=outputFaces(parametricFaceFaces(face,t,*p,referenceKeys(points),r,rc ? rc->native():nullptr)) ;
      if(patches) nativeFacePatches(faces,t,*p,&face) ;
      return faces ;
    }
    if(c.hasCornerGeometry()) {
      require(!rc || rc->hasCornerGeometry(),"mixed corner and other split families are unsupported") ;
      return cornerFaceFaces(face,c,points,rc,patches) ;
    }
    require(!rc || rc->extruded(),"mixed native and center-based boundaries are unsupported") ;
    ExtrudedView left(t,*c.extruded()) ; const int layer=left.cap(face) ;
    if(layer>=0) { require(!rc,"center-based extrusion end faces must be boundary faces") ; left.addBoundaryPoints(points) ; return left.caps(layer) ; }
    const unsigned edge=left.side(face) ; const auto lu=left.sideUses(edge) ; std::vector<EdgeUse> ru ; std::set<Fraction> cuts ;
    for(const auto &u:lu) { cuts.insert(left.coordinate(edge,u.use.begin)) ; cuts.insert(left.coordinate(edge,u.use.end)) ; }
    if(rc) {
      const ExtrudedView right(*r,*rc->extruded()) ; const unsigned re=right.side(face) ; ru=right.sideUses(re) ;
      require(left.layout.planes==right.layout.planes && c.extruded()->axis==rc->extruded()->axis,"neighbor extrusion planes differ") ;
      for(auto &u:ru) { u.use.begin=right.coordinate(re,u.use.begin) ; u.use.end=right.coordinate(re,u.use.end) ; cuts.insert(u.use.begin) ; cuts.insert(u.use.end) ; }
    }
    std::vector<ParametricOutputFace> result ;
    for(auto a=cuts.begin(),b=std::next(a);b!=cuts.end();++a,++b) {
      int li=-1,ri=-1 ;
      for(const auto &u:lu) { const auto x=left.coordinate(edge,u.use.begin),y=left.coordinate(edge,u.use.end) ; if(std::min(x,y)<=*a && std::max(x,y)>=*b) { require(li<0,"overlapping extrusion side cells") ; li=u.leaf ; } }
      for(const auto &u:ru) if(std::min(u.use.begin,u.use.end)<=*a && std::max(u.use.begin,u.use.end)>=*b) { require(ri<0,"overlapping neighbor extrusion cells") ; ri=u.leaf ; }
      require(li>=0 && (!rc || ri>=0),"extrusion side coverage is incomplete") ;
      const auto x=left.coordinate(edge,*a),y=left.coordinate(edge,*b) ; const auto lo=std::min(x,y),hi=std::max(x,y) ;
      result.push_back({{{left.edgePoint(edge,lo,0),-1},{left.edgePoint(edge,hi,0),-1},{left.edgePoint(edge,hi,1),-1},{left.edgePoint(edge,lo,1),-1}},li,ri}) ;
      const auto& pa=left.layout.polygon[edge] ; const auto& pb=left.layout.polygon[(edge+1)%left.layout.polygon.size()] ;
      result.back().patch=left.sidePatch({2,static_cast<long long>(face)},*a,*b,left.reversed(edge) ? -1:1,std::hypot(pb[0]-pa[0],pb[1]-pa[1])) ;
    }
    return result ;
  }

  bool ParametricCell::mergeFace(EntityId face,const ParametricRecord &record) {
    if(record.empty()) return false ;
    if(source()) {
      auto candidate=*this ; candidate.expand() ; const bool changed=candidate.mergeFace(face,record) ;
      if(changed) *this=std::move(candidate) ;
      return changed ;
    }
    if(native()) { validateMixedFacePlan(record) ; return false ; }
    if(hasCornerGeometry()) return applyFacePlan(cornerGeometry(),cornerFaceId(face),record) ;
    else require(record.empty(),"shared corner construction supplied to another split family") ;
    return false ;
  }
  ParametricRecord parametricFacePlan(EntityId face,const ParametricCell &cell) {
    if(!cell.hasCornerGeometry()) return {} ;
    const auto &boundary=cell.cornerGeometry() ;
    return encodeFacePlan(boundary.boundaries,cell.cornerFaceId(face),boundary.boundaryLimits()) ;
  }
  std::vector<Point> parametricCellPositions(const ParametricCell &cell,const std::vector<ParametricNodeKey> &keys) {
    if(cell.source()) { require(keys.empty(),"original cell has no interior nodes") ; return {} ; }
    if(cell.hasCornerGeometry()) return cornerPositions(cell,keys) ;
    const auto geometry=cell.history().geometry().map(cell.history().node(0).map) ; std::vector<Point> out ;
    for(const auto &key:keys) out.push_back(geometry->evaluate(approximate(key.reference.point)).position) ;
    return out ;
  }
  namespace {
    std::set<int> cornerContacts(const std::vector<ParametricOutputFace> &faces,const std::vector<int> &depth,
                                unsigned difference,NeighborGrading grading,int rightOffset=0) {
      std::set<int> out ; std::map<std::vector<ParametricNodeKey>,std::set<int>> contacts ;
      for(const auto &face:faces) {
        if(face.right>=0) gradePair(out,face.left,rightOffset+face.right,depth,difference) ;
        if(grading==NeighborGrading::faces) continue ;
        for(unsigned i=0;i<face.vertices.size();++i) {
          std::vector<ParametricNodeKey> key{face.vertices[i]} ;
          if(grading==NeighborGrading::edges) { key.push_back(face.vertices[(i+1)%face.vertices.size()]) ; std::sort(key.begin(),key.end()) ; }
          auto &incident=contacts[key] ; incident.insert(face.left) ; if(face.right>=0) incident.insert(rightOffset+face.right) ;
        }
      }
      for(const auto &entry:contacts) for(int a:entry.second) for(int b:entry.second) gradePair(out,a,b,depth,difference) ;
      return out ;
    }
  }

  FaceBalance parametricFaceBalance(const std::string &name) {
    if(name=="none") return FaceBalance::none ;
    if(name=="majority") return FaceBalance::majority ;
    if(name=="opposite") return FaceBalance::opposite ;
    throw std::invalid_argument("parametricFaceBalance must be none, majority or opposite") ;
  }

  namespace {
    struct CornerFaceCover { int leaf, face ; std::vector<unsigned> path ; } ;
    std::vector<CornerFaceCover> cornerFaceCover(EntityId face, const ParametricCell &cell) {
      std::vector<CornerFaceCover> out ;
      if(cell.source()) {
        for(unsigned f=0;f<cell.source()->topology.faceCount();++f)
          if(cell.source()->topology.faceId(f)==face) out.push_back({0,int(f),{}}) ;
        return out ;
      }
      const auto root=cell.cornerFaceId(face) ; const auto &b=cell.cornerGeometry().boundaries ; const auto leaves=cell.leafIds() ;
      for(unsigned i=0;i<leaves.size();++i) {
        const auto &uses=cell.cornerCell(leaves[i]).faces ;
        for(unsigned f=0;f<uses.size();++f) {
          auto id=uses[f].surface ; std::vector<unsigned> path ;
          while(id!=root && b.faces[id].parent!=noBoundary) {
            const auto &child=b.faces[id] ; const auto &siblings=b.faces[child.parent].splits.at(child.parentSplit).children ;
            path.push_back(std::find(siblings.begin(),siblings.end(),id)-siblings.begin()) ; id=child.parent ;
          }
          if(id==root) { std::reverse(path.begin(),path.end()) ; out.push_back({int(i),int(f),std::move(path)}) ; }
        }
      }
      return out ;
    }
    std::vector<FaceBalanceRequest> localDividedFaces(const CellTopology &t,const ParametricCell &cell,bool directions) {
      if(cell.source()) return {} ;
      if(cell.native()) return parametricDividedFaces(t,*cell.native(),directions) ;
      const auto leaves=cell.leafIds() ; std::vector<FaceBalanceRequest> out ;
      if(cell.hasCornerGeometry()) {
        const auto &b=cell.cornerGeometry().boundaries ; std::vector<bool> divided(b.faces.size(),false) ;
        for(auto leaf:leaves) for(const auto &use:cell.cornerCell(leaf).faces) {
          auto parent=b.faces[use.surface].parent ;
          while(parent!=noBoundary && !divided[parent]) { divided[parent]=true ; parent=b.faces[parent].parent ; }
        }
        for(unsigned i=0;i<leaves.size();++i) {
          const auto &faces=cell.cornerCell(leaves[i]).faces ;
          for(unsigned f=0;f<faces.size();++f) if(divided[faces[f].surface]) out.push_back({int(t.id()),int(i),int(f)}) ;
        }
      } else {
        const ExtrudedView v(t,*cell.extruded()) ;
        for(unsigned i=0;i<leaves.size();++i) {
          const auto &edges=v.tree.node(leaves[i]).edges ;
          for(unsigned e=0;e<edges.size();++e) {
            const auto &use=edges[e] ; const auto &cuts=v.cuts.at(use.edge) ;
            const auto next=cuts.upper_bound(std::min(use.begin,use.end)) ;
            if(next!=cuts.end() && *next<std::max(use.begin,use.end)) out.push_back({int(t.id()),int(i),int(e)}) ;
          }
        }
      }
      return out ;
    }
    std::vector<std::vector<unsigned>> leafFaceVertices(const CellTopology &t,const ParametricCell &cell,TreeNodeId leaf) {
      std::vector<std::vector<unsigned>> faces ;
      if(cell.source()) for(unsigned f=0;f<t.faceCount();++f) faces.push_back(t.faceVertices(f)) ;
      else if(cell.native()) {
        const auto &e=referenceElement(cell.native()->node(leaf).cell.shape) ;
        for(unsigned f=0;f<e.boundaryCount;++f) {
          const auto &b=e.boundaries[f] ; faces.emplace_back(b.vertices.begin(),b.vertices.begin()+b.vertexCount) ;
        }
      } else if(cell.hasCornerGeometry()) {
        const auto &c=cell.cornerCell(leaf) ; const auto &b=cell.cornerGeometry().boundaries ;
        for(auto use:c.faces) {
          std::vector<unsigned> vertices ;
          for(auto v:b.faces[use.surface].vertices) vertices.push_back(std::find(c.vertices.begin(),c.vertices.end(),v)-c.vertices.begin()) ;
          faces.push_back(std::move(vertices)) ;
        }
      } else {
        const unsigned n=cell.extruded()->planar.node(leaf).corners.size() ;
        for(unsigned v=0;v<n;++v) faces.push_back({v,(v+1)%n}) ;
      }
      return faces ;
    }
  }

  std::vector<FaceBalanceRequest> parametricDividedFaces(EntityId face,const CellTopology &t,const ParametricCell &cell,
      const CellTopology &right,const ParametricCell &other,bool directions) {
    if(cell.native() && other.native()) return parametricDividedFaces(face,t,*cell.native(),right,*other.native(),directions) ;
    if((cell.native() && other.source()) || (cell.source() && other.native())) {
      auto left=cell,next=other ; left.expand() ; next.expand() ;
      return mixedDividedFaces(face,t,left,right,next) ;
    }
    if(mixedFace(cell,&other)) return mixedDividedFaces(face,t,cell,right,other) ;
    std::vector<FaceBalanceRequest> out ;
    if((cell.source() || cell.hasCornerGeometry()) && (other.source() || other.hasCornerGeometry())) {
      const auto a=cornerFaceCover(face,cell), b=cornerFaceCover(face,other) ; std::set<std::vector<unsigned>> divided ;
      for(const auto *cover:{&a,&b}) for(const auto &use:*cover) {
        auto path=use.path ; while(!path.empty()) { path.pop_back() ; divided.insert(path) ; }
      }
      for(const auto &use:a) if(divided.count(use.path)) out.push_back({int(t.id()),use.leaf,use.face}) ;
      for(const auto &use:b) if(divided.count(use.path)) out.push_back({int(right.id()),use.leaf,use.face}) ;
    } else {
      require(cell.extruded() && other.extruded(),"mixed split families do not support face balancing") ;
      const auto faces=parametricFaceFaces(face,t,cell,{},&right,&other) ; std::map<int,unsigned> a,b ;
      for(const auto &piece:faces) { ++a[piece.left] ; ++b[piece.right] ; }
      const auto collect=[&](const CellTopology &topology,const ParametricCell &c,const std::map<int,unsigned> &counts) {
        const ExtrudedView v(topology,*c.extruded()) ; const auto side=v.side(face) ;
        for(const auto &entry:counts) if(entry.second>1) {
          const auto &edges=v.tree.node(v.leaves.at(entry.first)).edges ;
          for(unsigned e=0;e<edges.size();++e)
            if(edges[e].edge.original && edges[e].edge.index==side) out.push_back({int(topology.id()),entry.first,int(e)}) ;
        }
      } ;
      collect(t,cell,a) ; collect(right,other,b) ;
    }
    return out ;
  }
  std::vector<int> parametricFaceBalanceRequests(const CellTopology &t,const ParametricCell &cell,
      const std::vector<FaceBalanceRequest> &shared,FaceBalance policy,std::map<int,unsigned> *directions) {
    if(policy==FaceBalance::none && !directions) return {} ;
    if(directions) directions->clear() ;
    if(!cell.native()) directions=nullptr ;
    std::map<int,std::vector<unsigned>> divided ; std::set<int> fullSplit ;
    const auto collect=[&](const FaceBalanceRequest &entry) {
      require(entry.face>=0,"face balancing received a level-grading request") ; divided[entry.leaf].push_back(entry.face) ;
      if(directions) {
        require(entry.directions<8,"invalid balancing directions") ;
        if(entry.directions==0) fullSplit.insert(entry.leaf) ;
        (*directions)[entry.leaf]|=entry.directions ;
      }
    } ;
    for(const auto &entry:localDividedFaces(t,cell,directions!=nullptr)) collect(entry) ;
    for(const auto &entry:shared) if(entry.root==int(t.id())) collect(entry) ;
    if(directions) for(int leaf:fullSplit) (*directions)[leaf]=0 ;
    const auto leaves=cell.leafIds() ; std::vector<int> out ;
    for(const auto &entry:divided)
      if(faceBalanceRequired(leafFaceVertices(t,cell,leaves.at(entry.first)),entry.second,policy)) out.push_back(entry.first) ;
    return out ;
  }

  std::vector<int> parametricLocalRequests(const CellTopology &t,const ParametricCell &c,unsigned difference,NeighborGrading grading) {
    if(c.source()) return {} ;
    if(const auto *p=c.native()) return parametricLocalRequests(t,*p,difference,grading) ;
    if(c.hasCornerGeometry()) {
      auto faces=cornerCellFaces(c,{}) ;
      for(unsigned f=0;f<c.cornerCell(0).faces.size();++f) {
        const auto boundary=cornerFaceFaces(t.faceId(f),c,{}) ; faces.insert(faces.end(),boundary.begin(),boundary.end()) ;
      }
      const auto needed=cornerContacts(faces,depths(c),difference,grading) ; return {needed.begin(),needed.end()} ;
    }
    const ExtrudedView v(t,*c.extruded()) ; const auto d=depths(c) ; std::set<int> result ;
    for(const auto &f:v.internal()) gradePair(result,f.left,f.right,d,difference) ;
    if(grading!=NeighborGrading::faces) {
      std::map<ParametricNodeKey,std::set<int>> contacts ;
      for(const auto &f:v.caps(0)) for(const auto &p:f.vertices) contacts[p].insert(f.left) ;
      for(const auto &entry:contacts) for(int a:entry.second) for(int b:entry.second) gradePair(result,a,b,d,difference) ;
    }
    return {result.begin(),result.end()} ;
  }
  ParametricTrace parametricFaceRequests(EntityId face,const CellTopology &t,const ParametricCell &c,const CellTopology &r,const ParametricCell &rc,unsigned difference,NeighborGrading grading) {
    if(c.source() && rc.source()) return {} ;
    if(c.source() || rc.source()) {
      auto left=c,right=rc ; left.expand() ; right.expand() ; return parametricFaceRequests(face,t,left,r,right,difference,grading) ;
    }
    if(mixedFace(c,&rc)) return mixedFaceRequests(face,t,c,r,rc,difference,grading) ;
    if(const auto *p=c.native()) { require(rc.native(),"mixed native and center-based boundaries are unsupported") ; return parametricFaceRequests(face,t,*p,r,*rc.native(),difference,grading) ; }
    if(c.hasCornerGeometry()) {
      require(rc.hasCornerGeometry(),"mixed corner and other split families are unsupported") ;
      auto lc=c,rr=rc ;
      const auto shared=mergeFacePlans(parametricFacePlan(face,c),parametricFacePlan(face,rc),c.cornerGeometry().boundaryLimits()) ;
      lc.mergeFace(face,shared) ; rr.mergeFace(face,shared) ;
      const auto faces=parametricFaceFaces(face,t,lc,{},&r,&rr) ; auto d=depths(c) ; const int offset=d.size() ; const auto rd=depths(rc) ; d.insert(d.end(),rd.begin(),rd.end()) ;
      const auto needed=cornerContacts(faces,d,difference,grading,offset) ; ParametricTrace out ;
      for(int i:needed) out.push_back({int(i<offset ? t.id():r.id()),i<offset ? i:i-offset}) ;
      return out ;
    }
    require(rc.extruded(),"mixed native and center-based boundaries are unsupported") ;
    const auto faces=parametricFaceFaces(face,t,c,{},&r,&rc) ; auto d=depths(c) ; const int offset=d.size() ; const auto rd=depths(rc) ; d.insert(d.end(),rd.begin(),rd.end()) ;
    std::set<int> needed ; std::map<ParametricNodeKey,std::set<int>> contacts ;
    for(const auto &f:faces) {
      gradePair(needed,f.left,offset+f.right,d,difference) ;
      if(grading!=NeighborGrading::faces) for(const auto &p:f.vertices) { contacts[p].insert(f.left) ; contacts[p].insert(offset+f.right) ; }
    }
    for(const auto &entry:contacts) for(int a:entry.second) for(int b:entry.second) gradePair(needed,a,b,d,difference) ;
    ParametricTrace result ; for(int i:needed) result.push_back({int(i<offset ? t.id():r.id()),i<offset ? i:i-offset}) ; return result ;
  }
  ParametricTrace parametricEdgeTrace(const CellTopology &t,const ParametricCell &c) {
    if(c.source()) {
      ParametricTrace out ; for(unsigned edge=0;edge<t.edgeCount();++edge)
        out.push_back({int(t.edgeId(edge)),int(t.id()),0,0,0,1,1,1}) ;
      return out ;
    }
    if(const auto *p=c.native()) return parametricEdgeTrace(t,*p) ;
    if(c.hasCornerGeometry()) return cornerEdgeTrace(c) ;
    const ExtrudedView v(t,*c.extruded()) ; ParametricTrace out ; const auto &e=referenceElement(t.shape()) ;
    for(unsigned i=0;i<v.leaves.size();++i) {
      const int depth=c.history().path(v.leaves[i]).size() ; std::vector<VertexKey> points ;
      for(const auto &p:v.tree.node(v.leaves[i]).corners) for(unsigned layer=0;layer<2;++layer) { const auto key=v.key(p,layer) ; if(key.local<0) points.push_back(key.reference) ; }
      for(unsigned edge=0;edge<e.edgeCount;++edge) {
        const auto va=t.vertexId(e.edges[edge][0]),vb=t.vertexId(e.edges[edge][1]) ; const auto ends=std::minmax(va,vb) ; std::vector<Fraction> cuts ;
        for(const auto &p:points) {
          if(p.kind==EntityKind::edge && p.entity==t.edgeId(edge)) cuts.push_back(p.point[0]) ;
          if(p.kind==EntityKind::vertex && (p.entity==ends.first || p.entity==ends.second)) cuts.push_back(p.entity==ends.first ? 0:1) ;
        }
        if(cuts.empty()) continue ;
        const auto range=std::minmax_element(cuts.begin(),cuts.end()) ;
        out.push_back({int(t.edgeId(edge)),int(t.id()),int(i),depth,checkedExactInteger<int>(range.first->numerator()),checkedExactInteger<int>(range.first->denominator()),
                       checkedExactInteger<int>(range.second->numerator()),checkedExactInteger<int>(range.second->denominator())}) ;
      }
    }
    return out ;
  }
  std::vector<int> parametricVertexDepths(const CellTopology &t,const ParametricCell &c) {
    if(c.source()) return std::vector<int>(t.vertexCount(),0) ;
    if(const auto *p=c.native()) return parametricVertexDepths(t,*p) ;
    if(c.hasCornerGeometry()) return cornerVertexDepths(c) ;
    const ExtrudedView v(t,*c.extruded()) ; std::vector<int> out(t.vertexCount(),0) ;
    for(auto id:v.leaves) for(const auto &p:v.tree.node(id).corners) if(p.kind==PlanarVertexKey::Kind::corner)
      for(unsigned layer=0;layer<2;++layer) { const unsigned i=(layer ? v.layout.upper:v.layout.lower).at(p.corner) ; out[i]=std::max(out[i],int(c.history().path(id).size())) ; }
    return out ;
  }
  std::vector<int> parametricVertexRequests(const CellTopology &t,const ParametricCell &c,const std::map<EntityId,int> &depth,unsigned difference) {
    if(c.source()) {
      for(unsigned v=0;v<t.vertexCount();++v)
        if(depth.at(t.vertexId(v))>int(difference)) return {0} ;
      return {} ;
    }
    if(const auto *p=c.native()) return parametricVertexRequests(t,*p,depth,difference) ;
    if(c.hasCornerGeometry()) return cornerVertexRequests(c,depth,difference) ;
    const ExtrudedView v(t,*c.extruded()) ; std::set<int> out ;
    for(unsigned i=0;i<v.leaves.size();++i) {
      const int own=c.history().path(v.leaves[i]).size() ;
      for(const auto &p:v.tree.node(v.leaves[i]).corners) if(p.kind==PlanarVertexKey::Kind::corner)
        for(unsigned layer=0;layer<2;++layer) { const int other=depth.at(v.key(p,layer).reference.entity) ; if(other>own && unsigned(other-own)>difference) out.insert(i) ; }
    }
    return {out.begin(),out.end()} ;
  }

}
}
