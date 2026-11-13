// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
// Neighbor queries and cell correspondence from compact hex plans.
#include "compact_hex.h"
#include <algorithm>
#include <numeric>
#include <set>
#include <unordered_map>

namespace Loci {
namespace FVMAdapt {
  namespace {
    constexpr auto extent = HexPlan::extent ;
    using IndexPoint = std::array<std::uint32_t,3> ;
    unsigned contactAxes(NeighborGrading grading) {
      return grading == NeighborGrading::faces ? 1 : grading == NeighborGrading::edges ? 2 : 3 ;
    }
    // Translate retained plan-node IDs to the current output leaf order.
    struct LeafIndex {
      const HexPlan &plan ;
      std::vector<HexPlan::LeafBox> boxes ;
      std::unordered_map<HexPlan::Cell,unsigned> indices ;
      explicit LeafIndex(const HexPlan &p) : plan(p), boxes(p.leafBoxes()) {
        for(unsigned i = 0; i < boxes.size(); ++i) indices.emplace(boxes[i].cell,i) ;
      }
      unsigned at(const IndexPoint &point) const { return indices.at(plan.findLeaf(point)) ; }
      void grade(std::set<int> &required, const IndexPoint &point, unsigned depth, unsigned difference) const {
        const auto i = at(point) ; const auto own = boxes[i].depth ;
        if(depth > own && depth-own > difference) required.insert(i) ;
      }
    } ;
    // Integer point just beyond either side, or inside the interval when side is zero.
    bool probe(std::uint32_t low, std::uint32_t width, int side, std::uint32_t &point) {
      if(side < 0) { if(low == 0) return false ; point = low-1 ; }
      else if(side > 0) { point = low+width ; if(point == extent) return false ; }
      else point = low+width/2 ;
      return true ;
    }
    IndexPoint vertexPoint(unsigned vertex) {
      IndexPoint p ; const auto &v = referenceElement(Shape::hex).vertices[vertex] ;
      for(unsigned d = 0; d < 3; ++d) p[d] = v[d] ? extent-1 : 0 ;
      return p ;
    }
    // Canonical face coordinates use the same vertex orientation as shared face plans.
    struct FaceAxes {
      unsigned normal, u, v ; bool high, reverseU, reverseV ;
      FaceAxes(const CellTopology &t, EntityId face) {
        const auto &e = referenceElement(Shape::hex) ; unsigned f = 0 ;
        while(f < e.boundaryCount && t.faceId(f) != face) ++f ;
        if(f == e.boundaryCount) throw std::invalid_argument("face is not incident to the hex") ;
        std::array<Point,4> corners ;
        for(unsigned i = 0; i < 4; ++i) corners[orientedVertex(Shape::quad,t.faceOrientation(f),i)] = e.vertices[e.boundaries[f].vertices[i]] ;
        for(unsigned d = 0; d < 3; ++d) {
          if(corners[0][d] != corners[1][d]) { u = d ; reverseU = corners[0][d] > corners[1][d] ; }
          else if(corners[0][d] != corners[3][d]) { v = d ; reverseV = corners[0][d] > corners[3][d] ; }
          else { normal = d ; high = corners[0][d] != 0 ; }
        }
      }
      bool incident(const HexPlan::LeafBox &box) const {
        return high ? box.lower[normal]+box.width() == extent : box.lower[normal] == 0 ;
      }
      std::array<std::uint32_t,2> lower(const HexPlan::LeafBox &box) const {
        return {{reverseU ? extent-box.lower[u]-box.width() : box.lower[u],
                 reverseV ? extent-box.lower[v]-box.width() : box.lower[v]}} ;
      }
      IndexPoint point(std::uint32_t x, std::uint32_t y) const {
        IndexPoint result ; result[normal] = high ? extent-1 : 0 ;
        result[u] = reverseU ? extent-1-x : x ; result[v] = reverseV ? extent-1-y : y ; return result ;
      }
    } ;
  }
  std::vector<int> CompactHexCell::localRequests(unsigned difference, NeighborGrading grading) const {
    const LeafIndex index(plan) ; std::set<int> required ;
    // A fine box's contact with a coarser box fits on one face, edge or corner.
    // Probe from the fine side; the boxes nest regardless of the physical centers.
    for(const auto &box : index.boxes) for(int x = -1; x <= 1; ++x) for(int y = -1; y <= 1; ++y) for(int z = -1; z <= 1; ++z) {
      const unsigned axes = (x != 0)+(y != 0)+(z != 0) ;
      if(!axes || axes > contactAxes(grading)) continue ;
      IndexPoint point ;
      if(probe(box.lower[0],box.width(),x,point[0]) && probe(box.lower[1],box.width(),y,point[1]) && probe(box.lower[2],box.width(),z,point[2]))
        index.grade(required,point,box.depth,difference) ;
    }
    return {required.begin(),required.end()} ;
  }
  std::vector<FaceBalanceRequest> hexDividedFaces(const CellTopology &topology,const HexPlan &plan) {
    const LeafIndex index(plan) ; const auto &element=referenceElement(Shape::hex) ; std::set<std::pair<int,int>> divided ;
    for(unsigned f=0;f<6;++f) {
      const FaceAxes axes(topology,topology.faceId(f)) ;
      for(const auto &box:index.boxes) {
        IndexPoint point ;
        for(unsigned d=0;d<3;++d) point[d]=box.lower[d]+box.width()/2 ;
        if(!probe(box.lower[axes.normal],box.width(),axes.high ? 1:-1,point[axes.normal])) continue ;
        const auto other=index.at(point) ; if(index.boxes[other].depth>=box.depth) continue ;
        for(unsigned opposite=0;opposite<6;++opposite) {
          const auto &vertices=element.boundaries[opposite] ;
          bool match=true ;
          for(unsigned v=0;v<vertices.vertexCount;++v)
            match=match && (element.vertices[vertices.vertices[v]][axes.normal]==(axes.high ? 0:1)) ;
          if(match) divided.emplace(other,opposite) ;
        }
      }
    }
    std::vector<FaceBalanceRequest> out ;
    for(const auto &use:divided) out.push_back({int(topology.id()),use.first,use.second}) ;
    return out ;
  }
  std::vector<FaceBalanceRequest> hexDividedFaces(EntityId face,const CellTopology &left,const HexPlan &lp,
      const CellTopology &right,const HexPlan &rp) {
    const LeafIndex a(lp),b(rp) ; const FaceAxes af(left,face),bf(right,face) ; std::vector<FaceBalanceRequest> out ;
    const auto visit=[&](const LeafIndex &from,const FaceAxes &ff,const LeafIndex &to,const FaceAxes &tf,const CellTopology &t) {
      unsigned f=0 ; while(t.faceId(f)!=face) ++f ; std::set<int> divided ;
      for(const auto &box:from.boxes) if(ff.incident(box)) {
        const auto low=ff.lower(box) ; const auto other=to.at(tf.point(low[0]+box.width()/2,low[1]+box.width()/2)) ;
        if(to.boxes[other].depth<box.depth) divided.insert(other) ;
      }
      for(int leaf:divided) out.push_back({int(t.id()),leaf,int(f)}) ;
    } ;
    visit(a,af,b,bf,right) ; visit(b,bf,a,af,left) ; return out ;
  }
  std::vector<CompactHexCell::FaceLeaf> CompactHexCell::faceLeaves(EntityId face) const {
    const auto boxes = plan.leafBoxes() ; const FaceAxes axes(topology,face) ; std::vector<FaceLeaf> result ;
    for(unsigned i = 0; i < boxes.size(); ++i)
      if(axes.incident(boxes[i])) result.push_back({axes.lower(boxes[i]),boxes[i].width(),int(i)}) ;
    return result ;
  }
  ParametricTrace CompactHexCell::edgeTrace() const {
    const auto boxes = plan.leafBoxes() ; const auto &e = referenceElement(Shape::hex) ; ParametricTrace result ;
    for(unsigned i = 0; i < boxes.size(); ++i) for(unsigned edge = 0; edge < e.edgeCount; ++edge) {
      const auto &box = boxes[i] ; const auto &a = e.vertices[e.edges[edge][0]], &b = e.vertices[e.edges[edge][1]] ;
      bool touches = true ; unsigned axis = 0 ;
      for(unsigned d = 0; d < 3; ++d) {
        if(a[d] != b[d]) axis = d ;
        else touches = touches && (a[d] ? box.lower[d]+box.width() == extent : box.lower[d] == 0) ;
      }
      if(!touches) continue ;
      auto lo = box.lower[axis], hi = lo+box.width() ;
      const bool forward = (topology.vertexId(e.edges[edge][0]) < topology.vertexId(e.edges[edge][1])) == (a[axis] < b[axis]) ;
      if(!forward) { const auto old = lo ; lo = extent-hi ; hi = extent-old ; }
      const auto lg = std::gcd(lo,extent), hg = std::gcd(hi,extent) ;
      result.push_back({int(topology.edgeId(edge)),int(topology.id()),int(i),int(box.depth),int(lo/lg),int(extent/lg),int(hi/hg),int(extent/hg)}) ;
    }
    return result ;
  }
  std::vector<int> CompactHexCell::vertexDepths() const {
    std::vector<int> result ;
    for(unsigned v = 0; v < 8; ++v) result.push_back(plan.depth(plan.findLeaf(vertexPoint(v)))) ;
    return result ;
  }
  std::vector<int> CompactHexCell::vertexRequests(const std::map<EntityId,int> &depth, unsigned difference) const {
    const LeafIndex index(plan) ; std::set<int> result ;
    for(unsigned v = 0; v < 8; ++v) {
      const int other = depth.at(topology.vertexId(v)) ;
      if(other > 0) index.grade(result,vertexPoint(v),other,difference) ;
    }
    return {result.begin(),result.end()} ;
  }
  ParametricTrace hexFaceRequests(EntityId face, const CellTopology &left, const HexPlan &lp,
      const CellTopology &right, const HexPlan &rp, unsigned difference, NeighborGrading grading) {
    const LeafIndex a(lp), b(rp) ; const FaceAxes af(left,face), bf(right,face) ; std::set<int> ar, br ;
    const auto visit = [&](const LeafIndex &from, const FaceAxes &ff, std::set<int> &own,
                           const LeafIndex &to, const FaceAxes &tf, std::set<int> &other) {
      for(const auto &box : from.boxes) if(ff.incident(box)) {
        const auto low = ff.lower(box) ;
        for(int x = -1; x <= 1; ++x) for(int y = -1; y <= 1; ++y) {
          if(unsigned((x != 0)+(y != 0)) >= contactAxes(grading)) continue ;
          std::uint32_t u, v ; if(!probe(low[0],box.width(),x,u) || !probe(low[1],box.width(),y,v)) continue ;
          to.grade(other,tf.point(u,v),box.depth,difference) ;
          if(grading != NeighborGrading::faces) from.grade(own,ff.point(u,v),box.depth,difference) ;
        }
      }
    } ;
    visit(a,af,ar,b,bf,br) ; visit(b,bf,br,a,af,ar) ; ParametricTrace result ;
    for(int i : ar) result.push_back({int(left.id()),i}) ;
    for(int i : br) result.push_back({int(right.id()),i}) ;
    return result ;
  }
  std::vector<std::vector<int>> hexCellDonors(const HexPlan &oldPlan, const HexPlan &newPlan) {
    std::vector<std::vector<int>> result ; int oldIndex = 0 ;
    const auto collectOld = [&](const auto &self, HexPlan::Cell cell, std::vector<int> &donors) -> void {
      if(oldPlan.leaf(cell)) donors.push_back(oldIndex++) ;
      else for(auto child : oldPlan.children(cell)) self(self,child,donors) ;
    } ;
    const auto collectNew = [&](const auto &self, HexPlan::Cell cell) -> void {
      if(newPlan.leaf(cell)) result.push_back({oldIndex}) ;
      else for(auto child : newPlan.children(cell)) self(self,child) ;
    } ;
    const auto visit = [&](const auto &self, HexPlan::Cell oldCell, HexPlan::Cell newCell) -> void {
      if(newPlan.leaf(newCell)) { result.emplace_back() ; collectOld(collectOld,oldCell,result.back()) ; }
      else if(oldPlan.leaf(oldCell)) { collectNew(collectNew,newCell) ; ++oldIndex ; }
      else {
        if(!(oldPlan.centers(oldCell) == newPlan.centers(newCell))) throw std::invalid_argument("new cell has no ancestral donors") ;
        const auto a = oldPlan.children(oldCell), b = newPlan.children(newCell) ;
        for(unsigned i = 0; i < 8; ++i) self(self,a[i],b[i]) ;
      }
    } ;
    visit(visit,0,0) ; return result ;
  }
}
}
