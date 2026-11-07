// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include "parametric_hex.h"
#include <FVMAdapt/parametric/geometry_quality.h>
#include "library/parametric/face_plan_io.h"
#include "library/parametric/geometry_validity.h"

namespace Loci {
namespace FVMAdapt {
  namespace {
    constexpr unsigned absent = std::numeric_limits<unsigned>::max() ;
    void require(bool ok, const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    Point midpoint(const Point &a, const Point &b) {
      Point p ; for(unsigned d = 0; d < 3; ++d) p[d] = .5*a[d]+.5*b[d] ; return p ;
    }
    void finite(const Point &p) {
      for(double x : p) require(std::isfinite(x),"nonfinite shared boundary point") ;
    }

    struct FaceNode {
      Point center{{0,0,0}} ;
      std::array<unsigned,4> children{{absent,absent,absent,absent}} ;
      unsigned parent = absent, corner = 0 ;
      bool split() const { return children[0] != absent ; }
    } ;
    // A received face tree borrows its root map; child maps are evaluated along the path.
    struct FaceTree {
      const CornerMap &source ;
      const HexPlan &plan ;
      HexPlan::Cell owner ;
      unsigned face, orientation ;
      bool original ;
      std::vector<FaceNode> nodes{{}} ;

      Point position(unsigned id, Point p) const {
        while(nodes[id].parent != absent) {
          const auto &node = nodes[id] ; const auto &parent = nodes[node.parent] ; const auto v = node.corner ;
          p = quadChildMap(v,parent.center).position(p) ;
          id = node.parent ;
        }
        if(original) {
          const auto &boundary = referenceElement(Shape::hex).boundaries[face] ; std::vector<Point> corners(4) ;
          for(unsigned v = 0; v < 4; ++v) corners[orientedVertex(Shape::quad,orientation,v)] = source.corners()[boundary.vertices[v]] ;
          return CornerMap(Shape::quad,std::move(corners)).position(p) ;
        }
        return HexCellMap(source,plan,owner).position(boundaryPoint(Shape::hex,face,p)) ;
      }
      void readEdge(detail::FacePlanReader &in, unsigned id, Point a, Point b, unsigned depth) {
        if(!in.edge(depth)) return ;
        const auto p = midpoint(a,b) ; finite(position(id,p)) ;
        readEdge(in,id,a,p,depth+1) ; readEdge(in,id,p,b,depth+1) ;
      }
      void readFace(detail::FacePlanReader &in, unsigned id, unsigned depth) {
        const auto division = in.face(4,depth) ; const auto &quad = referenceElement(Shape::quad) ;
        for(unsigned e = 0; e < 4; ++e) readEdge(in,id,quad.vertices[e],quad.vertices[(e+1)%4],depth) ;
        if(!division.split) return ;
        nodes[id].center = division.center ; finite(position(id,division.center)) ;
        for(unsigned v = 0; v < 4; ++v) {
          finite(position(id,midpoint(quad.vertices[v],quad.vertices[(v+1)%4]))) ;
          const unsigned child = nodes.size() ; nodes[id].children[v] = child ;
          FaceNode node ; node.parent = id ; node.corner = v ; nodes.push_back(node) ;
          readFace(in,child,depth+1) ;
        }
      }
      void read(const FacePlan &record, const BoundaryRecordLimits &limits) {
        if(record.empty()) return ;
        detail::FacePlanReader in(record,limits) ; require(in.vertices == 4,"shared face plan shapes differ") ;
        readFace(in,0,0) ; in.input.finish() ;
      }
    } ;
    struct PlanFace {
      unsigned tree = absent, node = 0, orientation = 0 ;
    } ;
    using CellFaces = std::array<PlanFace,6> ;
    using FaceCorners = std::array<unsigned,4> ;

    const VertexVolumeSplit &pattern() {
      static const auto split = splitVertexVolume(Shape::hex,std::vector<Point>(6,{{.5,.5,0}}),{{.5,.5,.5}}) ; return split ;
    }
    FaceCorners childFace(unsigned child, unsigned face) {
      const auto &boundary = referenceElement(Shape::hex).boundaries[face] ; FaceCorners vertices ;
      for(unsigned v = 0; v < 4; ++v) vertices[v] = pattern().children[child].vertices[boundary.vertices[v]] ;
      return vertices ;
    }
    unsigned orientation(const FaceCorners &local, const FaceCorners &canonical) {
      for(unsigned o = 0; o < 8; ++o) {
        bool match = true ;
        for(unsigned v = 0; v < 4; ++v) match = match && local[v] == canonical[orientedVertex(Shape::quad,o,v)] ;
        if(match) return o ;
      }
      throw std::logic_error("invalid hex split face orientation") ;
    }
    unsigned edgePoint(int a, int b) {
      const auto &hex = referenceElement(Shape::hex) ;
      for(unsigned e = 0; e < 12; ++e)
        if((hex.edges[e][0] == a && hex.edges[e][1] == b) || (hex.edges[e][0] == b && hex.edges[e][1] == a)) return 8+e ;
      throw std::logic_error("invalid hex split edge") ;
    }
    FaceCorners boundaryChild(unsigned face, unsigned o, unsigned corner) {
      FaceCorners vertices ; const auto &boundary = referenceElement(Shape::hex).boundaries[face] ;
      for(unsigned v = 0; v < 4; ++v) vertices[orientedVertex(Shape::quad,o,v)] = boundary.vertices[v] ;
      const auto a = vertices[corner] ;
      return {{a,edgePoint(a,vertices[(corner+1)%4]),20+face,edgePoint(a,vertices[(corner+3)%4])}} ;
    }
    const std::array<InternalFace,12> &interiorOrder() {
      // Records number interior roots by first appearance in child-face order.
      static const auto faces = [] {
        std::array<InternalFace,12> result ; unsigned next = 0 ;
        for(unsigned c = 0; c < 8; ++c) for(unsigned f = 0; f < 6; ++f)
          for(const auto &pair : pattern().internalFaces)
            if(pair.first.child == c && pair.first.boundary == f) result[next++] = pair ;
        require(next == 12,"hex construction requires twelve interior faces") ; return result ;
      }() ; return faces ;
    }
  }

  void CompactHexCell::validate() const {
    const CornerMap source(Shape::hex,corners) ;
    BoundaryRecordLimits boundaryLimits ; boundaryLimits.maxSplits = boundaryLimits.maxEntities ;
    boundaryLimits.maxBytes = limits.exact.maxBytes ;
    std::vector<FaceTree> trees ; trees.reserve(6+12*plan.splitCount()) ;
    std::vector<CellFaces> cells(1+8*plan.splitCount()) ;
    for(unsigned f = 0; f < 6; ++f) {
      trees.push_back({source,plan,0,f,topology.faceOrientation(f),true}) ; trees.back().read(faces[f],boundaryLimits) ;
      cells[0][f] = {f,0,topology.faceOrientation(f)} ;
    }
    std::size_t next = 0 ;
    for(std::size_t i = 0; i < plan.splitCount(); ++i) {
      const auto parent = plan.splitParent(i) ; const auto centers = plan.centers(parent) ; const auto children = plan.children(parent) ;
      for(unsigned f = 0; f < 6; ++f) {
        const auto use = cells[parent][f] ; const auto &node = trees.at(use.tree).nodes.at(use.node) ;
        require(node.split(),"accepted hex face choice is missing") ;
        const auto p = orientPoint(Shape::quad,inverseOrientation(Shape::quad,use.orientation),node.center) ;
        require(centers.faces[f] == std::array<double,2>{{p[0],p[1]}},"hex and face plan centers disagree") ;
      }
      for(const auto &piece : pattern().boundaryPieces) {
        const auto use = cells[parent][piece.boundary] ;
        const auto corner = orientedVertex(Shape::quad,use.orientation,piece.corner) ;
        const auto node = trees[use.tree].nodes[use.node].children[corner] ;
        const auto o = orientation(childFace(piece.face.child,piece.face.boundary),boundaryChild(piece.boundary,use.orientation,corner)) ;
        cells[children[piece.face.child]][piece.face.boundary] = {use.tree,node,o} ;
      }
      for(unsigned f = 0; f < 12; ++f) {
        const auto &pair = interiorOrder()[f] ; const unsigned tree = trees.size() ;
        trees.push_back({source,plan,children[pair.first.child],pair.first.boundary,0,false}) ;
        if(next < interior.size() && interior[next].split == i && interior[next].face == f)
          trees.back().read(interior[next++].plan,boundaryLimits) ;
        cells[children[pair.first.child]][pair.first.boundary] = {tree,0,0} ;
        cells[children[pair.second.child]][pair.second.boundary] = {tree,0,inverseOrientation(Shape::quad,pair.orientation)} ;
      }
    }
    require(next == interior.size(),"unused hex interior face plan") ;
    // Check retained (including dormant) boundaries and the selected physical leaves.
    const std::array<Point,5> probes{{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{.5,.5,0}}}} ;
    bool sourceProved=false ;
    for(HexPlan::Cell id = 0; id < cells.size(); ++id) {
      const HexCellMap map(source,plan,id) ;
      const auto &hex = referenceElement(Shape::hex) ;
      const double extent = detail::pointExtent(8,[&](std::size_t v) { return map.position(hex.vertices[v]) ; }) ;
      for(unsigned f = 0; f < 6; ++f) {
        const auto use = cells[id][f] ;
        for(const auto &p : probes)
          require(detail::mappedPointAgrees(map.position(boundaryPoint(Shape::hex,f,p)),
            trees.at(use.tree).position(use.node,orientPoint(Shape::quad,use.orientation,p)),extent),
            "saved shared geometry disagrees with cell geometry") ;
      }
      if(plan.active(id) && plan.leaf(id)) {
        const auto quality=geometryQuality(map) ;
        require(!quality.invalidEvaluation,"invalid compact hex leaf evaluation") ;
        if(quality.status!=GeometryQualityStatus::sampled) {
          if(!sourceProved) sourceProved=geometryValidity(source,limits.jacobian)==GeometryValidity::valid ;
          require(sourceProved,"invalid compact hex leaf geometry") ;
        }
      }
    }
  }
}
}
