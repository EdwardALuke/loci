// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
// Shared-face nodes and polygons from compact hex plans.
#include "compact_hex.h"
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
      std::shared_ptr<const GeometryMap> source ; std::vector<Location> locations ;
      std::vector<ParametricNodeKey> keys ;
      // Face-local node number to temporary vertex index.
      std::vector<unsigned> localNodes ;
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
      Point position(unsigned face,Point p,const GeometryMap &geometry) const {
        while(faces[face].parent != absent) {
          const auto &f = faces[face] ; const unsigned v = f.corner ;
          p = quadChildMap(v,faces[f.parent].center).position(p) ;
          face = f.parent ;
        }
        return geometry.position(p) ;
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
      void facePoint(unsigned vertex) {
        if(keys[vertex].reference.kind == EntityKind::face && keys[vertex].local < 0) {
          keys[vertex].local = localNodes.size() ; localNodes.push_back(vertex) ;
        }
      }
      void faceEdgeKeys(unsigned id,unsigned start) {
        const auto &e = edges[id] ; if(e.middle == absent) return ;
        facePoint(e.middle) ; const bool reverse = e.vertices[0] != start ;
        faceEdgeKeys(e.children[reverse ? 1:0],start) ; faceEdgeKeys(e.children[reverse ? 0:1],e.middle) ;
      }
      void faceKeys(unsigned id) {
        const auto &face = faces[id] ; for(unsigned e = 0; e < 4; ++e) faceEdgeKeys(face.edges[e],face.vertices[e]) ;
        if(face.centerVertex == absent) return ;
        facePoint(face.centerVertex) ; for(auto child : face.children) faceKeys(child) ;
      }
      void cover(const CompactHexCell &cell,bool right) {
        const bool reverse = cell.topology.faceOrientation(faceIndex(cell.topology,owner)) >= 4 ;
        for(const auto &leaf : cell.faceLeaves(owner)) {
          auto id = faceIds.at({{leaf.lower[0],leaf.lower[1],leaf.width}}) ; faces[id].covers.push_back({leaf.leaf,right,reverse}) ;
          while(faces[id].parent != absent) { id = faces[id].parent ; faces[id].active = true ; }
        }
      }
      FaceOutput(EntityId face,const CompactHexCell &left,const CompactHexCell *right,const std::vector<ParametricNodeKey> &required = {}) : owner(face) {
        const auto f = faceIndex(left.topology,face),orientation = left.topology.faceOrientation(f) ;
        const auto &hex = referenceElement(Shape::hex) ; const auto &boundary = hex.boundaries[f] ;
        source=left.faceMap(f) ; keys.resize(4) ; std::array<unsigned,4> native ;
        for(unsigned v = 0; v < 4; ++v) {
          const unsigned i = orientedVertex(Shape::quad,orientation,v),vertex = boundary.vertices[v] ; native[i] = vertex ;
          keys[i] = {{EntityKind::vertex,left.topology.vertexId(vertex),{{0,0,0}}},-1} ;
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
        faceKeys(0) ; cover(left,false) ; if(right) cover(*right,true) ;
      }
      void splitParents(unsigned id) {
        while(edges[id].parent != absent) { id = edges[id].parent ; edges[id].active = true ; }
      }
      template<class Visit> void visitEdgeVertices(unsigned id,unsigned start,const Visit &visit) const {
        const auto &e = edges[id] ;
        if(!e.active) { visit(start) ; return ; }
        require(e.middle != absent,"active face edge has no midpoint") ; const bool reverse = e.vertices[0] != start ;
        visitEdgeVertices(e.children[reverse ? 1:0],start,visit) ; visitEdgeVertices(e.children[reverse ? 0:1],e.middle,visit) ;
      }
      /// Propagate coverage and activate edges once, before this temporary's output query.
      std::vector<unsigned> prepareFaces(bool neighbor) {
        std::vector<unsigned> selected ;
        for(unsigned id = 0; id < faces.size(); ++id) {
          const auto &face = faces[id] ;
          if(face.active) for(auto child : face.children) {
            require(child != absent,"active face has no children") ;
            faces[child].covers.insert(faces[child].covers.end(),face.covers.begin(),face.covers.end()) ;
          } else if(!face.covers.empty()) {
            require(face.covers.size() == (neighbor ? 2u:1u),"hex face coverage is incomplete") ;
            if(neighbor) require(face.covers[0].reverse != face.covers[1].reverse,"hex face orientations disagree") ;
            selected.push_back(id) ;
            for(auto edge : face.edges) splitParents(edge) ;
          }
        }
        return selected ;
      }
      std::vector<ParametricOutputFace> output(bool neighbor,const std::vector<ParametricNodeKey> &required = {}) {
        const auto selected = prepareFaces(neighbor) ;
        const std::set<ParametricNodeKey> needed(required.begin(),required.end()) ;
        for(unsigned e = 0; e < edges.size(); ++e) if(edges[e].middle != absent && needed.count(keys[edges[e].middle])) {
          edges[e].active = true ; splitParents(e) ;
        }
        std::vector<ParametricOutputFace> result ; result.reserve(selected.size()) ;
        for(auto id : selected) {
          auto &face = faces[id] ; if(face.covers.front().right) std::reverse(face.covers.begin(),face.covers.end()) ;
          ParametricOutputFace output{{},face.covers[0].leaf,neighbor ? face.covers[1].leaf:-1} ;
          const auto append = [&](unsigned vertex) { output.vertices.push_back(keys[vertex]) ; } ;
          for(unsigned e = 0; e < 4; ++e) visitEdgeVertices(face.edges[e],face.vertices[e],append) ;
          if(face.covers[0].reverse) std::reverse(output.vertices.begin(),output.vertices.end()) ;
          result.push_back(std::move(output)) ;
        }
        return result ;
      }
      std::vector<ParametricNodeKey> nodes(bool neighbor) {
        const auto selected = prepareFaces(neighbor) ; std::vector<unsigned char> used(keys.size(),0) ; std::size_t count = 0 ;
        const auto mark = [&](unsigned vertex) {
          if(!used[vertex] && keys[vertex].reference.kind == EntityKind::face) { used[vertex] = 1 ; ++count ; }
        } ;
        for(auto id : selected) {
          const auto &face = faces[id] ; for(unsigned e = 0; e < 4; ++e) visitEdgeVertices(face.edges[e],face.vertices[e],mark) ;
        }
        std::vector<ParametricNodeKey> result ; result.reserve(count) ;
        for(auto vertex : localNodes) if(used[vertex]) result.push_back(keys[vertex]) ;
        return result ;
      }
      std::vector<Point> positions(const std::vector<ParametricNodeKey> &nodes) const {
        if(nodes.empty()) return {} ;
        std::vector<Point> result ; result.reserve(nodes.size()) ;
        for(const auto &key:nodes) {
          const auto vertex=localNodes.at(key.local) ;
          require(keys[vertex]==key,"face node does not belong to the output face") ;
          const auto &p=locations[vertex] ; result.push_back(position(p.face,p.point,*source)) ;
        }
        return result ;
      }
    } ;
  }
  std::vector<ParametricNodeKey> CompactHexCell::faceNodes(EntityId face,const CompactHexCell *right,std::vector<Point> *positions) const {
    FaceOutput output(face,*this,right) ; auto nodes=output.nodes(right!=nullptr) ;
    if(positions) *positions=output.positions(nodes) ;
    return nodes ;
  }
  std::vector<Point> CompactHexCell::facePositions(EntityId face,const CompactHexCell *right) const {
    std::vector<Point> positions ; faceNodes(face,right,&positions) ; return positions ;
  }
  std::vector<ParametricOutputFace> CompactHexCell::faceFaces(EntityId face,const std::vector<ParametricNodeKey> &required,const CompactHexCell *right) const {
    return FaceOutput(face,*this,right,required).output(right != nullptr,required) ;
  }
}
}
