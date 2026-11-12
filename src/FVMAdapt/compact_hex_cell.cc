// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
// Cell refinement and cell-interior output from compact hex plans.
#include "compact_hex.h"
#include <FVMAdapt/parametric_cell.h>
#include "library/parametric/face_plan_io.h"
#include "library/parametric/geometry_validity.h"
#include <set>

namespace Loci {
namespace FVMAdapt {
  namespace {
    constexpr unsigned absent = std::numeric_limits<unsigned>::max() ;
    using Quad = std::array<unsigned,4> ;
    void require(bool ok,const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    Point midpoint(const Point &a,const Point &b) {
      Point p ; for(unsigned d=0;d<3;++d) p[d]=(a[d]+b[d])/2 ; return p ;
    }
    const VertexVolumeSplit &pattern() {
      static const auto split=splitVertexVolume(Shape::hex,std::vector<Point>(6,{{.5,.5,0}}),{{.5,.5,.5}}) ; return split ;
    }
    std::shared_ptr<const GeometryMap> placedChild(std::shared_ptr<const GeometryMap> parent,
                                                  const VertexVolumeSplit &split,unsigned child) {
      std::vector<Point> points ; for(auto v:split.children.at(child).vertices) points.push_back(split.points[v]) ;
      return std::make_shared<RestrictedMap>(std::move(parent),std::make_shared<CornerMap>(Shape::hex,std::move(points))) ;
    }
    // Expand only for a boundary-singularity proof; ordinary cells use compact evaluators.
    std::shared_ptr<const GeometryMap> placedCell(std::shared_ptr<const GeometryMap> source,const HexPlan &plan,unsigned id) {
      std::vector<unsigned> path ;
      for(auto cell=id;cell!=HexPlan::root;cell=plan.parent(cell)) path.push_back(cell) ;
      std::shared_ptr<const GeometryMap> map=std::move(source) ;
      for(auto i=path.rbegin();i!=path.rend();++i) {
        const auto parent=plan.parent(*i) ; const auto siblings=plan.children(parent) ;
        map=placedChild(std::move(map),plan.construction(parent),std::find(siblings.begin(),siblings.end(),*i)-siblings.begin()) ;
      }
      return map ;
    }
    unsigned orientation(const Quad &local,const Quad &canonical) {
      for(unsigned o=0;o<8;++o) {
        bool match=true ; for(unsigned v=0;v<4;++v) match=match && local[v]==canonical[orientedVertex(Shape::quad,o,v)] ;
        if(match) return o ;
      }
      throw std::invalid_argument("hex output face incidence differs") ;
    }

    /// Temporary cell/face/edge incidence shared by plan edits and mesh output.
    /// Geometry is borrowed from the source and accepted placements.
    struct HexBoundary {
      /// Coordinates on a face, or in a cell when face is absent.
      struct Location { unsigned face,cell ; Point point ; } ;
      /// Segment endpoints in the owning face's coordinates; children retain that map.
      struct Edge {
        std::array<unsigned,2> vertices,children{{absent,absent}} ;
        unsigned parent=absent,middle=absent,face ; Point a,b ; bool active=false ;
      } ;
      /// Current leaf index and orientation relative to a face.
      struct Cover { int leaf ; bool reverse ; } ;
      /// Accepted quad subdivision and the active leaves using it.
      struct Face {
        Quad vertices,edges{},children{{absent,absent,absent,absent}} ;
        unsigned parent=absent,corner=0,centerVertex=absent,cell=0,boundary=0 ;
        bool exterior=false,active=false ; Point center{} ; CenterMethod method=CenterMethod::reference ; std::vector<Cover> covers ;
      } ;
      /// Face identity and native-to-face orientation.
      struct Use { unsigned face,orientation ; } ;
      /// Hex incidence in native vertex and face order.
      struct Cell { std::array<unsigned,8> vertices ; std::array<Use,6> faces ; } ;
      const CompactHexCell &saved ;
      const std::shared_ptr<const GeometryMap> source ;
      std::array<std::shared_ptr<const GeometryMap>,6> sourceFaces ;
      std::vector<Location> locations ; std::vector<ParametricNodeKey> keys ;
      std::vector<Edge> edges ; std::vector<Face> faces ; std::vector<Cell> cells ;
      std::map<std::pair<unsigned,unsigned>,unsigned> edgeIds ;

      void checkEntity() const {
        if(locations.size()+edges.size()+faces.size()>=BoundaryRecordLimits().maxEntities)
          throw std::length_error("hex output entity limit exceeded") ;
      }
      unsigned point(Location location) {
        checkEntity() ; const unsigned id=locations.size() ; locations.push_back(location) ;
        keys.push_back({{EntityKind::cell,saved.topology.id(),{{0,0,0}}},int(id)}) ; return id ;
      }
      unsigned edge(unsigned a,unsigned b,unsigned face,Point p,Point q,unsigned parent=absent) {
        const std::pair<unsigned,unsigned> ends=std::minmax(a,b) ; const auto found=edgeIds.find(ends) ; if(found!=edgeIds.end()) return found->second ;
        if(b<a) { std::swap(a,b) ; std::swap(p,q) ; }
        checkEntity() ; const unsigned id=edges.size() ; Edge e ; e.vertices={{a,b}} ; e.parent=parent ; e.face=face ; e.a=p ; e.b=q ;
        edges.push_back(e) ; edgeIds.emplace(ends,id) ; return id ;
      }
      unsigned addFace(Face f) {
        checkEntity() ; const unsigned id=faces.size() ; faces.push_back(f) ; const auto &quad=referenceElement(Shape::quad) ;
        for(unsigned e=0;e<4;++e) faces[id].edges[e]=edge(f.vertices[e],f.vertices[(e+1)%4],id,quad.vertices[e],quad.vertices[(e+1)%4]) ;
        return id ;
      }
      unsigned splitEdge(unsigned id) {
        if(edges[id].middle!=absent) return edges[id].middle ;
        const auto e=edges[id] ; const auto p=midpoint(e.a,e.b) ; const auto mid=point({e.face,absent,p}) ; edges[id].middle=mid ;
        const auto first=edge(e.vertices[0],mid,e.face,e.a,p,id),second=edge(mid,e.vertices[1],e.face,p,e.b,id) ;
        edges[id].children={{first,second}} ; return mid ;
      }
      void readEdge(detail::FacePlanReader &in,unsigned id,unsigned start,unsigned depth) {
        if(!in.edge(depth)) return ;
        const auto mid=splitEdge(id) ; const auto e=edges[id] ; const bool reverse=e.vertices[0]!=start ;
        readEdge(in,e.children[reverse ? 1:0],start,depth+1) ; readEdge(in,e.children[reverse ? 0:1],mid,depth+1) ;
      }
      void divideFace(unsigned id,const Point &center,CenterMethod method) {
        const auto f=faces[id] ;
        if(f.centerVertex!=absent) { require(f.center==center,"hex face centers disagree") ; return ; }
        faces[id].center=center ; faces[id].method=method ;
        const auto vertex=point({id,absent,center}) ; faces[id].centerVertex=vertex ;
        Quad mid ; for(unsigned e=0;e<4;++e) mid[e]=splitEdge(f.edges[e]) ;
        const auto &q=referenceElement(Shape::quad) ;
        for(unsigned e=0;e<4;++e) edge(mid[e],vertex,id,midpoint(q.vertices[e],q.vertices[(e+1)%4]),center) ;
        for(unsigned v=0;v<4;++v) {
          Face child ; child.parent=id ; child.corner=v ; child.exterior=f.exterior ;
          child.vertices={{f.vertices[v],mid[v],vertex,mid[(v+3)%4]}} ;
          const auto next=addFace(child) ; faces[id].children[v]=next ;
        }
      }
      void readFace(detail::FacePlanReader &in,unsigned id,unsigned depth) {
        const auto division=in.face(4,depth) ; const auto f=faces[id] ;
        for(unsigned e=0;e<4;++e) readEdge(in,f.edges[e],f.vertices[e],depth) ;
        if(!division.split) return ;
        divideFace(id,division.center,division.method) ;
        const auto children=faces[id].children ; for(auto child:children) readFace(in,child,depth+1) ;
      }
      void read(unsigned face,const FacePlan &plan) {
        if(plan.empty()) return ;
        BoundaryRecordLimits limits ; limits.maxBytes=saved.limits.exact.maxBytes ;
        detail::FacePlanReader in(plan,limits) ; require(in.vertices==4,"hex output requires quad faces") ; readFace(in,face,0) ; in.input.finish() ;
      }
      Quad faceVertices(const Cell &cell,unsigned f) const {
        Quad vertices ; const auto &boundary=referenceElement(Shape::hex).boundaries[f] ;
        for(unsigned v=0;v<4;++v) vertices[v]=cell.vertices[boundary.vertices[v]] ;
        return vertices ;
      }
      void construct(unsigned split) {
        const auto parent=saved.plan.splitParent(split) ; const auto c=cells.at(parent) ;
        const auto &hex=referenceElement(Shape::hex) ; const auto &p=pattern() ; std::array<unsigned,27> vertices ;
        std::copy(c.vertices.begin(),c.vertices.end(),vertices.begin()) ;
        for(unsigned e=0;e<12;++e) vertices[8+e]=splitEdge(edgeIds.at(std::minmax(c.vertices[hex.edges[e][0]],c.vertices[hex.edges[e][1]]))) ;
        for(unsigned f=0;f<6;++f) {
          vertices[20+f]=faces[c.faces[f].face].centerVertex ; require(vertices[20+f]!=absent,"hex output is missing an accepted face split") ;
        }
        vertices[26]=point({absent,parent,saved.plan.centers(parent).cell}) ; const auto ids=saved.plan.children(parent) ;
        require(ids[0]==cells.size(),"hex output children are not contiguous") ; cells.resize(cells.size()+8) ;
        for(unsigned i=0;i<8;++i) for(unsigned v=0;v<8;++v) cells[ids[i]].vertices[v]=vertices[p.children[i].vertices[v]] ;
        const auto assign=[&](FaceUse use,unsigned face) {
          auto &cell=cells[ids[use.child]] ; cell.faces[use.boundary]={face,orientation(faceVertices(cell,use.boundary),faces[face].vertices)} ;
        } ;
        for(const auto &piece:p.boundaryPieces) {
          const auto use=c.faces[piece.boundary] ; const auto v=orientedVertex(Shape::quad,use.orientation,piece.corner) ;
          assign(piece.face,faces[use.face].children[v]) ;
        }
        for(const auto &pair:p.internalFaces) {
          Face f ; f.cell=ids[pair.first.child] ; f.boundary=pair.first.boundary ; f.vertices=faceVertices(cells[f.cell],f.boundary) ;
          const auto id=addFace(f) ; assign(pair.first,id) ; assign(pair.second,id) ;
        }
      }
      void edgeKeys(unsigned id,unsigned start,EntityId owner,Fraction lo,Fraction hi) {
        const auto &e=edges[id] ; if(e.middle==absent) return ; const auto mid=(lo+hi)/Fraction(2) ;
        keys[e.middle]={{EntityKind::edge,owner,{{mid,0,0}}},-1} ; const bool reverse=e.vertices[0]!=start ;
        edgeKeys(e.children[reverse ? 1:0],start,owner,lo,mid) ; edgeKeys(e.children[reverse ? 0:1],e.middle,owner,mid,hi) ;
      }
      void facePoint(unsigned vertex,EntityId owner,int &next) {
        if(keys[vertex].reference.kind==EntityKind::cell) keys[vertex]={{EntityKind::face,owner,{{0,0,0}}},next++} ;
      }
      void faceEdgeKeys(unsigned id,unsigned start,EntityId owner,int &next) {
        const auto &e=edges[id] ; if(e.middle==absent) return ; facePoint(e.middle,owner,next) ; const bool reverse=e.vertices[0]!=start ;
        faceEdgeKeys(e.children[reverse ? 1:0],start,owner,next) ; faceEdgeKeys(e.children[reverse ? 0:1],e.middle,owner,next) ;
      }
      void faceKeys(unsigned id,EntityId owner,int &next) {
        const auto &f=faces[id] ; for(unsigned e=0;e<4;++e) faceEdgeKeys(f.edges[e],f.vertices[e],owner,next) ;
        if(f.centerVertex==absent) return ;
        facePoint(f.centerVertex,owner,next) ; for(auto child:f.children) faceKeys(child,owner,next) ;
      }
      explicit HexBoundary(const CompactHexCell &state,const std::vector<ParametricNodeKey> &required={},bool forOutput=true)
        : saved(state),source(state.sourceMap()) {
        const auto &hex=referenceElement(Shape::hex) ; const auto &t=saved.topology ; Cell root ;
        for(unsigned v=0;v<8;++v) { root.vertices[v]=point({absent,0,hex.vertices[v]}) ; keys[v]={{EntityKind::vertex,t.vertexId(v),{}},-1} ; }
        for(unsigned f=0;f<6;++f) {
          sourceFaces[f]=saved.faceMap(f) ;
          Face face ; face.exterior=true ; face.boundary=f ; const auto &boundary=hex.boundaries[f] ;
          for(unsigned v=0;v<4;++v) face.vertices[orientedVertex(Shape::quad,t.faceOrientation(f),v)]=boundary.vertices[v] ;
          root.faces[f]={addFace(face),t.faceOrientation(f)} ;
        }
        cells.push_back(root) ; for(unsigned f=0;f<6;++f) read(f,saved.faces[f]) ;
        std::size_t next=0 ;
        for(unsigned i=0;i<saved.plan.splitCount();++i) {
          construct(i) ;
          if(next<saved.interior.size() && saved.interior[next].split==i) {
            std::vector<unsigned> internal ;
            for(auto child:saved.plan.children(saved.plan.splitParent(i))) for(auto use:cells[child].faces)
              if(faces[use.face].parent==absent && !faces[use.face].exterior &&
                 std::find(internal.begin(),internal.end(),use.face)==internal.end()) internal.push_back(use.face) ;
            while(next<saved.interior.size() && saved.interior[next].split==i) {
              read(internal.at(saved.interior[next].face),saved.interior[next].plan) ; ++next ;
            }
          }
        }
        require(next==saved.interior.size(),"unused hex output face plan") ;
        if(!forOutput) return ;
        // Extra original-edge points can be required by cells meeting only along that edge.
        for(unsigned e=0;e<12;++e) {
          const unsigned a=hex.edges[e][0],b=hex.edges[e][1] ; const bool reverse=t.vertexId(a)>t.vertexId(b) ;
          const auto edgeOwner=t.edgeId(e) ; const auto id=edgeIds.at(std::minmax(a,b)) ;
          for(const auto &key:required) if(key.reference.kind==EntityKind::edge && key.reference.entity==edgeOwner) {
            auto current=id,start=reverse ? b:a ; Fraction lo=0,hi=1 ; const auto target=key.reference.point[0] ;
            require(target>lo && target<hi,"original edge point is outside its edge") ;
            for(unsigned depth=0;;++depth) {
              require(depth<31,"corner edge point exceeds supported depth") ; const auto mid=splitEdge(current) ; const auto value=(lo+hi)/Fraction(2) ;
              if(value==target) break ;
              const auto edge=edges[current] ; const bool flipped=edge.vertices[0]!=start ;
              if(target<value) { current=edge.children[flipped ? 1:0] ; hi=value ; }
              else { current=edge.children[flipped ? 0:1] ; start=mid ; lo=value ; }
            }
          }
          edgeKeys(id,reverse ? b:a,edgeOwner,0,1) ;
        }
        for(unsigned f=0;f<6;++f) { int ordinal=0 ; faceKeys(f,t.faceId(f),ordinal) ; }
        const auto leaves=saved.plan.leaves() ;
        for(unsigned i=0;i<leaves.size();++i) for(auto use:cells[leaves[i]].faces) {
          faces[use.face].covers.push_back({int(i),use.orientation>=4}) ; auto f=use.face ;
          while(faces[f].parent!=absent) { f=faces[f].parent ; faces[f].active=true ; }
        }
      }
      std::shared_ptr<const GeometryMap> faceGeometry(unsigned id) const {
        const auto &f=faces[id] ; const auto &q=referenceElement(Shape::quad) ;
        if(f.parent!=absent) {
          const auto v=f.corner ;
          auto placement=std::make_shared<CornerMap>(quadChildMap(v,faces[f.parent].center)) ;
          return std::make_shared<RestrictedMap>(faceGeometry(f.parent),std::move(placement)) ;
        }
        std::vector<Point> corners ;
        if(f.exterior) return sourceFaces[f.boundary] ;
        for(unsigned v=0;v<4;++v) corners.push_back(boundaryPoint(Shape::hex,f.boundary,q.vertices[v])) ;
        return std::make_shared<RestrictedMap>(std::make_shared<HexCellMap>(*source,saved.plan,f.cell),
          std::make_shared<CornerMap>(Shape::quad,std::move(corners))) ;
      }
      void writeEdge(detail::RecordWriter &out,unsigned id,unsigned start,detail::FacePlanBudget &budget,unsigned depth) const {
        budget.visit(depth) ; const auto &e=edges[id] ; const bool divided=e.middle!=absent ;
        budget.division(divided) ; out.byte(divided) ; if(!divided) return ;
        const bool reverse=e.vertices[0]!=start ;
        writeEdge(out,e.children[reverse ? 1:0],start,budget,depth+1) ;
        writeEdge(out,e.children[reverse ? 0:1],e.middle,budget,depth+1) ;
      }
      void writeFace(detail::RecordWriter &out,unsigned id,detail::FacePlanBudget &budget,unsigned depth) const {
        budget.visit(depth) ; const auto &f=faces[id] ; const bool divided=f.centerVertex!=absent ; budget.division(divided) ;
        detail::writeDivision(out,{divided,f.method,f.center},4) ;
        for(unsigned e=0;e<4;++e) writeEdge(out,f.edges[e],f.vertices[e],budget,depth) ;
        if(divided) for(auto child:f.children) writeFace(out,child,budget,depth+1) ;
      }
      FacePlan record(unsigned id) const {
        const auto &f=faces[id] ; bool divided=f.centerVertex!=absent ;
        for(auto e:f.edges) divided=divided || edges[e].middle!=absent ;
        if(!divided) return {} ;
        BoundaryRecordLimits limits ; limits.maxBytes=saved.limits.exact.maxBytes ; limits.maxSplits=limits.maxEntities ;
        detail::RecordWriter out(limits.maxBytes,limits.geometry.exact) ; detail::header(out,4) ;
        detail::FacePlanBudget budget{limits} ; writeFace(out,id,budget,0) ; return out.take() ;
      }
      void save(CompactHexCell &target) const {
        for(unsigned f=0;f<6;++f) target.faces[f]=record(f) ;
        target.interior.clear() ;
        for(unsigned i=0;i<target.plan.splitCount();++i) {
          std::vector<unsigned> roots ;
          for(auto child:target.plan.children(target.plan.splitParent(i))) for(auto use:cells[child].faces)
            if(faces[use.face].parent==absent && !faces[use.face].exterior &&
               std::find(roots.begin(),roots.end(),use.face)==roots.end()) roots.push_back(use.face) ;
          require(roots.size()==12,"hex construction requires twelve interior faces") ;
          for(unsigned f=0;f<roots.size();++f) { auto bytes=record(roots[f]) ; if(!bytes.empty()) target.interior.push_back({i,f,std::move(bytes)}) ; }
        }
      }
      void checkCell(unsigned id) const {
        const HexCellMap map(*source,saved.plan,id) ; const auto &hex=referenceElement(Shape::hex) ;
        std::array<Point,8> corners ;
        for(unsigned v=0;v<8;++v) corners[v]=map.position(hex.vertices[v]) ;
        const double extent=detail::pointExtent(8,[&](std::size_t v) { return corners[v] ; }) ;
        const auto close=[extent](const Point &a,const Point &b) {
          require(detail::mappedPointAgrees(a,b,extent),"saved shared geometry disagrees with cell geometry") ;
        } ;
        for(unsigned v=0;v<8;++v) close(corners[v],position(cells[id].vertices[v])) ;
        const std::array<Point,5> probes{{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{.5,.5,0}}}} ;
        for(unsigned f=0;f<6;++f) {
          const auto use=cells[id].faces[f] ; const auto geometry=faceGeometry(use.face) ;
          for(const auto &p:probes) close(map.position(boundaryPoint(Shape::hex,f,p)),
            geometry->position(orientPoint(Shape::quad,use.orientation,p))) ;
        }
      }
      void refine(CompactHexCell &target,unsigned id,CenterMethod method) {
        require(target.plan.leaf(id),"hex refinement requires an active leaf") ;
        if(target.plan.accepted(id)) { target.plan.refine(id) ; return ; }
        const auto &limits=target.limits ;
        if(target.plan.depth(id)>=limits.maxDepth) throw std::length_error("refinement exceeds its depth limit") ;
        if(limits.maxNodes<9 || target.plan.splitCount()>=(limits.maxNodes-1)/8 || target.plan.splitCount()>=limits.maxSplits)
          throw std::length_error("refinement exceeds its node or split limit") ;
        HexCenters centers ; std::vector<Point> faceCenters ;
        for(unsigned f=0;f<6;++f) {
          const auto use=cells[id].faces[f] ;
          if(faces[use.face].centerVertex==absent) divideFace(use.face,faceCenterCoordinates(*faceGeometry(use.face),method),method) ;
          const auto p=orientPoint(Shape::quad,inverseOrientation(Shape::quad,use.orientation),faces[use.face].center) ;
          centers.faces[f]={{p[0],p[1]}} ; faceCenters.push_back(p) ;
        }
        const HexCellMap parent(*source,target.plan,id) ; centers.cell=cellCenterCoordinates(parent,faceCenters,method) ;
        const auto split=target.plan.splitCount() ; target.plan.refine(id,centers,limits.jacobian) ;
        std::vector<std::shared_ptr<const GeometryMap>> children ;
        for(auto child:target.plan.children(id)) children.push_back(std::make_shared<HexCellMap>(*source,target.plan,child)) ;
        auto assessment=splitQuality(parent,children) ;
        bool finite=!assessment.parent.quality.invalidEvaluation ;
        for(const auto &child:assessment.children) finite=finite && !child.quality.invalidEvaluation ;
        if(assessment.status!=GeometryQualityStatus::sampled && finite && geometryValidity(*source,limits.jacobian)==GeometryValidity::valid) {
          const auto geometry=placedCell(source,target.plan,id) ; const auto construction=target.plan.construction(id) ;
          children.clear() ; for(unsigned child=0;child<8;++child) children.push_back(placedChild(geometry,construction,child)) ;
          assessment=splitQuality(*geometry,children) ;
        }
        requireSplitQuality(assessment) ;
        construct(split) ; for(auto child:target.plan.children(id)) checkCell(child) ;
      }
      Point position(unsigned vertex) const {
        const auto &location=locations.at(vertex) ;
        if(location.face==absent) return HexCellMap(*source,saved.plan,location.cell).position(location.point) ;
        auto id=location.face ; auto p=location.point ;
        while(faces[id].parent!=absent) {
          const auto &f=faces[id] ; const auto v=f.corner ;
          p=quadChildMap(v,faces[f.parent].center).position(p) ; id=f.parent ;
        }
        const auto &f=faces[id] ;
        if(f.exterior) return sourceFaces[f.boundary]->position(p) ;
        return HexCellMap(*source,saved.plan,f.cell).position(boundaryPoint(Shape::hex,f.boundary,p)) ;
      }
      void splitParents(unsigned id) {
        while(edges[id].parent!=absent) { id=edges[id].parent ; edges[id].active=true ; }
      }
      template<class Visit> void visitEdgeVertices(unsigned id,unsigned start,const Visit &visit) const {
        const auto &e=edges[id] ; if(!e.active) { visit(start) ; return ; }
        require(e.middle!=absent,"active hex edge has no midpoint") ; const bool reverse=e.vertices[0]!=start ;
        visitEdgeVertices(e.children[reverse ? 1:0],start,visit) ; visitEdgeVertices(e.children[reverse ? 0:1],e.middle,visit) ;
      }
      /// Propagate coverage and activate edges once, before this temporary's output query.
      std::vector<unsigned> prepareFaces() {
        std::vector<unsigned> selected ;
        for(unsigned id=0;id<faces.size();++id) {
          const auto &f=faces[id] ; if(f.exterior) continue ;
          if(f.active) for(auto child:f.children) {
            require(child!=absent,"active hex face has no children") ; faces[child].covers.insert(faces[child].covers.end(),f.covers.begin(),f.covers.end()) ;
          } else if(!f.covers.empty()) {
            require(f.covers.size()==2 && f.covers[0].reverse!=f.covers[1].reverse,"hex interior face coverage differs") ;
            selected.push_back(id) ; for(auto edge:f.edges) splitParents(edge) ;
          }
        }
        return selected ;
      }
      std::vector<ParametricNodeKey> nodes() {
        const auto selected=prepareFaces() ; std::vector<unsigned char> used(keys.size(),0) ; std::size_t count=0 ;
        const auto mark=[&](unsigned vertex) {
          if(!used[vertex] && keys[vertex].reference.kind==EntityKind::cell) { used[vertex]=1 ; ++count ; }
        } ;
        for(auto id:selected) {
          const auto &f=faces[id] ; for(unsigned e=0;e<4;++e) visitEdgeVertices(f.edges[e],f.vertices[e],mark) ;
        }
        // Cell-owned keys use the vertex index as their local ID.
        std::vector<ParametricNodeKey> result ; result.reserve(count) ;
        for(unsigned vertex=0;vertex<keys.size();++vertex) if(used[vertex]) result.push_back(keys[vertex]) ;
        return result ;
      }
      std::vector<ParametricOutputFace> output(const std::vector<ParametricNodeKey> &required) {
        const auto selected=prepareFaces() ;
        const std::set<ParametricNodeKey> needed(required.begin(),required.end()) ;
        for(unsigned e=0;e<edges.size();++e) if(edges[e].middle!=absent && needed.count(keys[edges[e].middle])) { edges[e].active=true ; splitParents(e) ; }
        std::vector<ParametricOutputFace> result ; result.reserve(selected.size()) ;
        for(auto id:selected) {
          const auto &f=faces[id] ; ParametricOutputFace face{{},f.covers[0].leaf,f.covers[1].leaf} ;
          const auto append=[&](unsigned vertex) { face.vertices.push_back(keys[vertex]) ; } ;
          for(unsigned e=0;e<4;++e) visitEdgeVertices(f.edges[e],f.vertices[e],append) ;
          if(f.covers[0].reverse) std::reverse(face.vertices.begin(),face.vertices.end()) ;
          result.push_back(std::move(face)) ;
        }
        return result ;
      }
    } ;
  }
  bool CompactHexCell::edit(const std::vector<ParametricLeafRequest> &requests,int levels,CenterMethod method,
      const std::array<const FacePlan *,6> &incoming,const GeometricRefinement *selection) {
    std::unique_ptr<HexBoundary> boundary ; bool changed=!requests.empty() ;
    for(unsigned f=0;f<6;++f) if(incoming[f] && !incoming[f]->empty() && *incoming[f]!=faces[f]) {
      if(!boundary) boundary=std::make_unique<HexBoundary>(*this,std::vector<ParametricNodeKey>{},false) ;
      const auto points=boundary->locations.size() ; boundary->read(f,*incoming[f]) ;
      for(unsigned v=points;v<boundary->locations.size();++v) for(double x:boundary->position(v))
        require(std::isfinite(x),"nonfinite shared boundary point") ;
      changed=changed || boundary->locations.size()!=points ;
    }
    const auto refine=[&](HexPlan::Cell leaf) {
      if(plan.accepted(leaf)) plan.refine(leaf) ;
      else {
        if(!boundary) boundary=std::make_unique<HexBoundary>(*this,std::vector<ParametricNodeKey>{},false) ;
        boundary->refine(*this,leaf,method) ;
      }
    } ;
    for(const auto &request:requests) {
      if(request.action==ParametricAction::coarsen) { plan.coarsen(request.cell) ; continue ; }
      std::vector<HexPlan::Cell> leaves{HexPlan::Cell(request.cell)} ;
      for(int level=0;level<levels;++level) {
        std::vector<HexPlan::Cell> next ;
        for(auto leaf:leaves) {
          refine(leaf) ;
          if(level+1<levels) { const auto children=plan.children(leaf) ; next.insert(next.end(),children.begin(),children.end()) ; }
        }
        leaves.swap(next) ;
      }
    }
    if(selection && selection->select) {
      auto pending=plan.leaves() ; std::vector<Point> points ; points.reserve(8) ;
      const auto &hex=referenceElement(Shape::hex) ;
      while(!pending.empty()) {
        const auto leaf=pending.back() ; pending.pop_back() ; points.clear() ;
        if(leaf==0) points=corners ;
        else {
          if(!boundary) boundary=std::make_unique<HexBoundary>(*this,std::vector<ParametricNodeKey>{},false) ;
          for(auto vertex:boundary->cells[leaf].vertices) points.push_back(boundary->position(vertex)) ;
        }
        double shortest=std::numeric_limits<double>::max() ;
        for(unsigned e=0;e<12;++e) {
          const auto &a=points[hex.edges[e][0]],&b=points[hex.edges[e][1]] ; double square=0 ;
          for(unsigned d=0;d<3;++d) square+=(b[d]-a[d])*(b[d]-a[d]) ;
          shortest=std::min(shortest,std::sqrt(square)) ;
        }
        if(!selection->select(points,shortest) || shortest<2*selection->minEdgeLength) continue ;
        refine(leaf) ; changed=true ; const auto children=plan.children(leaf) ;
        pending.insert(pending.end(),children.begin(),children.end()) ;
      }
    }
    if(boundary && changed) boundary->save(*this) ;
    return changed ;
  }
  std::vector<ParametricNodeKey> CompactHexCell::cellNodes(std::vector<Point> *positions) const {
    if(positions) positions->clear() ;
    if(plan.leaf(0)) return {} ;
    HexBoundary view(*this) ; auto nodes=view.nodes() ;
    if(positions) {
      positions->reserve(nodes.size()) ;
      for(const auto &key:nodes) positions->push_back(view.position(key.local)) ;
    }
    return nodes ;
  }
  std::vector<Point> CompactHexCell::cellPositions(const std::vector<ParametricNodeKey> &keys) const {
    if(keys.empty()) return {} ;
    const HexBoundary view(*this) ; std::vector<Point> result ; result.reserve(keys.size()) ;
    for(const auto &key:keys) {
      require(key.local>=0 && unsigned(key.local)<view.keys.size() && view.keys[key.local]==key,"invalid hex interior node key") ;
      result.push_back(view.position(key.local)) ;
    }
    return result ;
  }
  std::vector<ParametricOutputFace> CompactHexCell::cellFaces(const std::vector<ParametricNodeKey> &required) const {
    if(plan.leaf(0)) return {} ;
    return HexBoundary(*this,required).output(required) ;
  }
}
}
