// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include "parametric_hex.h"
#include <FVMAdapt/parametric_cell.h>
#include "library/parametric/face_plan_io.h"
#include <set>

namespace Loci {
namespace FVMAdapt {
  namespace {
    constexpr unsigned absent = std::numeric_limits<unsigned>::max() ;
    using IndexPoint = std::array<std::uint32_t,2> ;
    void require(bool ok,const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    template<class T> T midpoint(const T &a,const T &b) {
      T p ; for(unsigned d = 0; d < p.size(); ++d) p[d] = (a[d]+b[d])/2 ; return p ;
    }
    unsigned faceIndex(const CellTopology &topology,EntityId face) {
      for(unsigned f = 0; f < 6; ++f) if(topology.faceId(f) == face) return f ;
      throw std::invalid_argument("face is not incident to the hex") ;
    }

    // Temporary connectivity for one original face. Geometry is evaluated along its saved face path.
    struct FaceOutput {
      struct Edge {
        std::array<unsigned,2> vertices,children{{absent,absent}} ;
        unsigned parent = absent, middle = absent ; bool active = false ;
      } ;
      struct Cover { int leaf ; bool right,reverse ; } ;
      struct Face {
        std::array<unsigned,4> vertices,edges{},children{{absent,absent,absent,absent}} ;
        // Midpoint boxes identify adjoining leaves; accepted centers determine physical positions.
        std::array<IndexPoint,4> box ;
        unsigned parent = absent,corner = 0,centerVertex = absent ;
        Point center{{0,0,0}} ; bool active = false ; std::vector<Cover> covers ;
      } ;
      EntityId owner ;
      struct Location { unsigned face ; Point point ; } ;
      std::vector<Point> source ; std::vector<Location> locations ;
      std::vector<ParametricNodeKey> keys ;
      std::vector<Edge> edges ; std::vector<Face> faces ;
      std::map<std::pair<unsigned,unsigned>,unsigned> edgeIds ;
      std::map<std::array<std::uint32_t,3>,unsigned> faceIds ;

      unsigned edge(unsigned a,unsigned b,unsigned parent = absent) {
        const auto key = std::minmax(a,b) ; const auto found = edgeIds.find(key) ;
        if(found != edgeIds.end()) return found->second ;
        const unsigned id = edges.size() ; Edge e ; e.vertices = {{a,b}} ; e.parent = parent ;
        edges.push_back(e) ; edgeIds.emplace(key,id) ; return id ;
      }
      unsigned addFace(Face face) {
        for(unsigned e = 0; e < 4; ++e) face.edges[e] = edge(face.vertices[e],face.vertices[(e+1)%4]) ;
        const auto &a = face.box[0], &b = face.box[2] ;
        const unsigned id = faces.size() ;
        faceIds.emplace(std::array<std::uint32_t,3>{{std::min(a[0],b[0]),std::min(a[1],b[1]),std::max(a[0],b[0])-std::min(a[0],b[0])}},id) ;
        faces.push_back(std::move(face)) ; return id ;
      }
      Point position(unsigned face,Point p) const {
        while(faces[face].parent != absent) {
          const auto &f = faces[face] ; const unsigned v = f.corner ;
          p = quadChildMap(v,faces[f.parent].center).evaluate(p).position ;
          face = f.parent ;
        }
        return CornerMap(Shape::quad,source).evaluate(p).position ;
      }
      unsigned point(unsigned face,const Point &p) {
        const unsigned id = locations.size() ; locations.push_back({face,p}) ;
        keys.push_back({{EntityKind::face,owner,{{0,0,0}}},-1}) ; return id ;
      }
      unsigned splitEdge(unsigned id,unsigned face,const Point &a,const Point &b) {
        if(edges[id].middle != absent) return edges[id].middle ;
        const auto ends = edges[id].vertices ; const unsigned mid = point(face,midpoint(a,b)) ;
        edges[id].middle = mid ;
        const auto first = edge(ends[0],mid,id),second = edge(mid,ends[1],id) ;
        edges[id].children = {{first,second}} ; return mid ;
      }
      void readEdge(detail::FacePlanReader &in,unsigned id,unsigned start,unsigned face,Point a,Point b,unsigned depth) {
        if(!in.edge(depth)) return ;
        const unsigned mid = splitEdge(id,face,a,b) ; const auto e = edges[id] ; const bool reverse = e.vertices[0] != start ;
        const auto p = midpoint(a,b) ;
        readEdge(in,e.children[reverse ? 1:0],start,face,a,p,depth+1) ;
        readEdge(in,e.children[reverse ? 0:1],mid,face,p,b,depth+1) ;
      }
      void readFace(detail::FacePlanReader &in,unsigned id,unsigned depth) {
        const auto division = in.face(4,depth) ; const auto face = faces[id] ; const auto &q = referenceElement(Shape::quad) ;
        for(unsigned e = 0; e < 4; ++e) readEdge(in,face.edges[e],face.vertices[e],id,q.vertices[e],q.vertices[(e+1)%4],depth) ;
        if(!division.split) return ;
        faces[id].center = division.center ; const unsigned center = point(id,division.center) ; faces[id].centerVertex = center ;
        std::array<unsigned,4> mid ;
        for(unsigned e = 0; e < 4; ++e) mid[e] = splitEdge(face.edges[e],id,q.vertices[e],q.vertices[(e+1)%4]) ;
        for(unsigned v = 0; v < 4; ++v) {
          Face child ; child.parent = id ; child.corner = v ;
          child.vertices = {{face.vertices[v],mid[v],center,mid[(v+3)%4]}} ;
          child.box = {{face.box[v],midpoint(face.box[v],face.box[(v+1)%4]),midpoint(face.box[0],face.box[2]),midpoint(face.box[v],face.box[(v+3)%4])}} ;
          const auto next = addFace(std::move(child)) ; faces[id].children[v] = next ;
        }
        const auto children = faces[id].children ; for(auto child : children) readFace(in,child,depth+1) ;
      }
      void edgeKeys(unsigned id,unsigned start,EntityId edgeOwner,Fraction lo,Fraction hi) {
        const auto &e = edges[id] ; if(e.middle == absent) return ;
        const auto mid = (lo+hi)/Fraction(2) ; keys[e.middle] = {{EntityKind::edge,edgeOwner,{{mid,0,0}}},-1} ;
        const bool reverse = e.vertices[0] != start ;
        edgeKeys(e.children[reverse ? 1:0],start,edgeOwner,lo,mid) ; edgeKeys(e.children[reverse ? 0:1],e.middle,edgeOwner,mid,hi) ;
      }
      void facePoint(unsigned vertex,int &next) { if(keys[vertex].reference.kind == EntityKind::face && keys[vertex].local < 0) keys[vertex].local = next++ ; }
      void faceEdgeKeys(unsigned id,unsigned start,int &next) {
        const auto &e = edges[id] ; if(e.middle == absent) return ;
        facePoint(e.middle,next) ; const bool reverse = e.vertices[0] != start ;
        faceEdgeKeys(e.children[reverse ? 1:0],start,next) ; faceEdgeKeys(e.children[reverse ? 0:1],e.middle,next) ;
      }
      void faceKeys(unsigned id,int &next) {
        const auto &face = faces[id] ; for(unsigned e = 0; e < 4; ++e) faceEdgeKeys(face.edges[e],face.vertices[e],next) ;
        if(face.centerVertex == absent) return ;
        facePoint(face.centerVertex,next) ; for(auto child : face.children) faceKeys(child,next) ;
      }
      void cover(const ParametricHex &cell,bool right) {
        const bool reverse = cell.topology.faceOrientation(faceIndex(cell.topology,owner)) >= 4 ;
        for(const auto &leaf : cell.faceLeaves(owner)) {
          auto id = faceIds.at({{leaf.lower[0],leaf.lower[1],leaf.width}}) ; faces[id].covers.push_back({leaf.leaf,right,reverse}) ;
          while(faces[id].parent != absent) { id = faces[id].parent ; faces[id].active = true ; }
        }
      }
      FaceOutput(EntityId face,const ParametricHex &left,const ParametricHex *right,const std::vector<ParametricNodeKey> &required = {}) : owner(face) {
        const auto f = faceIndex(left.topology,face),orientation = left.topology.faceOrientation(f) ;
        const auto &hex = referenceElement(Shape::hex) ; const auto &boundary = hex.boundaries[f] ;
        source.resize(4) ; keys.resize(4) ; std::array<unsigned,4> native ;
        for(unsigned v = 0; v < 4; ++v) {
          const unsigned i = orientedVertex(Shape::quad,orientation,v),vertex = boundary.vertices[v] ; native[i] = vertex ;
          source[i] = left.corners[vertex] ; keys[i] = {{EntityKind::vertex,left.topology.vertexId(vertex),{{0,0,0}}},-1} ;
        }
        const auto &q = referenceElement(Shape::quad) ;
        for(unsigned v = 0; v < 4; ++v) locations.push_back({0,q.vertices[v]}) ;
        Face root ; root.vertices = {{0,1,2,3}} ; constexpr auto n = HexPlan::extent ;
        root.box = {{{{0,0}},{{n,0}},{{n,n}},{{0,n}}}} ; addFace(root) ;
        const auto &a = left.faces[f] ;
        const auto &b = right ? right->faces[faceIndex(right->topology,face)] : a ;
        const auto record = a == b ? a : mergeFacePlans(a,b) ;
        if(!record.empty()) {
          BoundaryRecordLimits limits ; limits.maxBytes = left.limits.exact.maxBytes ; limits.maxSplits = limits.maxEntities ;
          detail::FacePlanReader in(record,limits) ; require(in.vertices == 4,"hex face requires four corners") ;
          readFace(in,0,0) ; in.input.finish() ;
        }
        for(unsigned e = 0; e < 4; ++e) {
          unsigned local = 0 ; const auto ends = std::minmax(native[e],native[(e+1)%4]) ;
          while(local < 12 && std::minmax(unsigned(hex.edges[local][0]),unsigned(hex.edges[local][1])) != ends) ++local ;
          require(local < 12,"invalid original hex edge") ; const auto edgeOwner = left.topology.edgeId(local) ;
          const bool reverse = keys[e].reference.entity > keys[(e+1)%4].reference.entity ;
          // A third incident cell can require original-edge points absent from this face plan.
          for(const auto &key : required) if(key.reference.kind == EntityKind::edge && key.reference.entity == edgeOwner) {
            auto id = faces[0].edges[e],start = e ; Point lo = q.vertices[e],hi = q.vertices[(e+1)%4] ;
            Fraction low = 0,high = 1,target = reverse ? Fraction(1)-key.reference.point[0] : key.reference.point[0] ;
            require(target > low && target < high,"original edge point is outside its edge") ;
            for(unsigned depth = 0;; ++depth) {
              require(depth < 31,"corner edge point exceeds supported depth") ; const auto mid = splitEdge(id,0,lo,hi) ;
              const auto value = (low+high)/Fraction(2) ; if(value == target) break ;
              const bool reversed = edges[id].vertices[0] != start ; const auto p = midpoint(lo,hi) ;
              if(target < value) { id = edges[id].children[reversed ? 1:0] ; high = value ; hi = p ; }
              else { id = edges[id].children[reversed ? 0:1] ; start = mid ; low = value ; lo = p ; }
            }
          }
          edgeKeys(faces[0].edges[e],reverse ? (e+1)%4:e,edgeOwner,0,1) ;
        }
        int next = 0 ; faceKeys(0,next) ; cover(left,false) ; if(right) cover(*right,true) ;
      }
      void splitParents(unsigned id) {
        while(edges[id].parent != absent) { id = edges[id].parent ; edges[id].active = true ; }
      }
      void edgeVertices(unsigned id,unsigned start,std::vector<ParametricNodeKey> &vertices) const {
        const auto &e = edges[id] ;
        if(!e.active) { vertices.push_back(keys[start]) ; return ; }
        require(e.middle != absent,"active face edge has no midpoint") ; const bool reverse = e.vertices[0] != start ;
        edgeVertices(e.children[reverse ? 1:0],start,vertices) ; edgeVertices(e.children[reverse ? 0:1],e.middle,vertices) ;
      }
      std::vector<ParametricOutputFace> output(bool neighbor,const std::vector<ParametricNodeKey> &required = {}) {
        std::vector<unsigned> selected ;
        for(unsigned id = 0; id < faces.size(); ++id) {
          const auto &face = faces[id] ;
          if(face.active) for(auto child : face.children) {
            require(child != absent,"active face has no children") ;
            faces[child].covers.insert(faces[child].covers.end(),face.covers.begin(),face.covers.end()) ;
          } else if(!face.covers.empty()) {
            require(face.covers.size() == (neighbor ? 2u:1u),"hex face coverage is incomplete") ; selected.push_back(id) ;
            for(auto edge : face.edges) splitParents(edge) ;
          }
        }
        const std::set<ParametricNodeKey> needed(required.begin(),required.end()) ;
        for(unsigned e = 0; e < edges.size(); ++e) if(edges[e].middle != absent && needed.count(keys[edges[e].middle])) {
          edges[e].active = true ; splitParents(e) ;
        }
        std::vector<ParametricOutputFace> result ;
        for(auto id : selected) {
          auto &face = faces[id] ; if(face.covers.front().right) std::reverse(face.covers.begin(),face.covers.end()) ;
          if(neighbor) require(face.covers[0].reverse != face.covers[1].reverse,"hex face orientations disagree") ;
          ParametricOutputFace output{{},face.covers[0].leaf,neighbor ? face.covers[1].leaf:-1} ;
          for(unsigned e = 0; e < 4; ++e) edgeVertices(face.edges[e],face.vertices[e],output.vertices) ;
          if(face.covers[0].reverse) std::reverse(output.vertices.begin(),output.vertices.end()) ;
          result.push_back(std::move(output)) ;
        }
        return result ;
      }
      std::vector<ParametricNodeKey> nodes(bool neighbor) {
        std::set<ParametricNodeKey> result ;
        for(const auto &face : output(neighbor)) for(const auto &key : face.vertices)
          if(key.reference.kind == EntityKind::face) result.insert(key) ;
        return {result.begin(),result.end()} ;
      }
    } ;
  }
  std::vector<ParametricNodeKey> ParametricHex::faceNodes(EntityId face,const ParametricHex *right) const {
    return FaceOutput(face,*this,right).nodes(right != nullptr) ;
  }
  std::vector<Point> ParametricHex::facePositions(EntityId face,const ParametricHex *right) const {
    FaceOutput output(face,*this,right) ; const auto nodes = output.nodes(right != nullptr) ;
    std::map<ParametricNodeKey,unsigned> points ;
    for(unsigned i = 0; i < output.keys.size(); ++i) points.emplace(output.keys[i],i) ;
    // Numbering and connectivity need only local coordinates; project the requested nodes here.
    std::vector<Point> result ;
    for(const auto &key : nodes) {
      const auto &p = output.locations[points.at(key)] ; result.push_back(output.position(p.face,p.point)) ;
    }
    return result ;
  }
  std::vector<ParametricOutputFace> ParametricHex::faceFaces(EntityId face,const std::vector<ParametricNodeKey> &required,const ParametricHex *right) const {
    return FaceOutput(face,*this,right,required).output(right != nullptr,required) ;
  }
}
}
