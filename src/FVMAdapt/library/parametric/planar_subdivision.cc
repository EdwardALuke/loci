//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#include <FVMAdapt/parametric/planar_subdivision.h>

#include <algorithm>
#include <stdexcept>
#include <type_traits>

namespace Loci {
namespace FVMAdapt {
  namespace {
    Point interpolate(const Point &a, const Point &b, const Fraction &position) {
      Point result ;
      const long double t=position.rational().convert_to<long double>() ;
      const long double s=(Fraction(1)-position).rational().convert_to<long double>() ;
      for(unsigned d=0; d<3; ++d) result[d]=static_cast<double>(s*a[d]+t*b[d]) ;
      return result ;
    }
    PlanarVertexKey cornerPoint(unsigned corner) {
      PlanarVertexKey result ; result.corner=corner ; return result ;
    }
    PlanarVertexKey centerPoint(const PlanarSplitKey &owner) {
      PlanarVertexKey result ; result.kind=PlanarVertexKey::Kind::center ; result.owner=owner ; return result ;
    }
    RefinementHistory initialHistory(std::vector<Point> corners, std::uint64_t root,
                                     std::uint64_t revision, HistoryLimits limits) {
      validatePlanarPolygon(corners) ;
      const auto shape=corners.size() == 3 ? Shape::triangle : Shape::quad ;
      GeometryStore geometry(limits.geometry) ;
      const auto map=geometry.append(CornerRecord{shape,std::move(corners)}) ;
      return RefinementHistory(root,std::move(geometry),map,revision,limits) ;
    }
    HistoryLimits nodeLimit(std::size_t maxNodes) {
      HistoryLimits limits ; limits.maxNodes=maxNodes ; return limits ;
    }
  }

  bool operator==(const PlanarSplitKey &a, const PlanarSplitKey &b) {
    return a.parent == b.parent && bool(a.construction) == bool(b.construction) &&
      (!a.construction || *a.construction == *b.construction) ;
  }
  bool operator<(const PlanarSplitKey &a, const PlanarSplitKey &b) {
    if(a.parent != b.parent) return a.parent < b.parent ;
    if(bool(a.construction) != bool(b.construction)) return bool(a.construction) < bool(b.construction) ;
    return a.construction && *a.construction < *b.construction ;
  }
  bool operator==(const PlanarEdgeKey &a, const PlanarEdgeKey &b) {
    return a.original == b.original && a.index == b.index && (a.original || a.owner == b.owner) ;
  }
  bool operator<(const PlanarEdgeKey &a, const PlanarEdgeKey &b) {
    if(a.original != b.original) return a.original < b.original ;
    if(!a.original && !(a.owner == b.owner)) return a.owner < b.owner ;
    return a.index < b.index ;
  }
  bool operator==(const PlanarVertexKey &a, const PlanarVertexKey &b) {
    if(a.kind != b.kind) return false ;
    switch(a.kind) {
    case PlanarVertexKey::Kind::corner: return a.corner == b.corner ;
    case PlanarVertexKey::Kind::edge: return a.edge == b.edge && a.position == b.position ;
    case PlanarVertexKey::Kind::center: return a.owner == b.owner ;
    }
    return false ;
  }
  bool operator<(const PlanarVertexKey &a, const PlanarVertexKey &b) {
    if(a.kind != b.kind) return a.kind < b.kind ;
    switch(a.kind) {
    case PlanarVertexKey::Kind::corner: return a.corner < b.corner ;
    case PlanarVertexKey::Kind::edge: return a.edge < b.edge || (a.edge == b.edge && a.position < b.position) ;
    case PlanarVertexKey::Kind::center: return a.owner < b.owner ;
    }
    return false ;
  }

  PlanarSubdivision::PlanarSubdivision(std::vector<Point> corners, std::size_t maxNodes)
    : PlanarSubdivision(std::move(corners),0,0,nodeLimit(maxNodes)) {}

  PlanarSubdivision::PlanarSubdivision(std::vector<Point> corners, std::uint64_t root,
                                       std::uint64_t revision, HistoryLimits limits)
    : PlanarSubdivision(initialHistory(std::move(corners),root,revision,limits)) {}

  PlanarSubdivision::PlanarSubdivision(RefinementHistory saved) : edit_(std::move(saved)) {
    const auto &history=this->history() ;
    const auto *rootMap=dynamic_cast<const CornerMap*>(history.geometry().map(history.node(0).map).get()) ;
    if(!rootMap || (rootMap->shape() != Shape::triangle && rootMap->shape() != Shape::quad))
      throw std::invalid_argument("planar subdivision requires a physical root triangle or quad") ;
    validatePlanarPolygon(rootMap->corners()) ;
    std::vector<std::size_t> owners(history.nodes().size(),noTreeNode) ;
    for(std::size_t i=0; i<history.splits().size(); ++i) {
      const auto &split=history.splits()[i] ;
      if(!split.construction->planar) throw std::invalid_argument("planar subdivision requires planar constructions") ;
      for(TreeNodeId child : split.children) owners[child]=i ;
    }
    nodes_.reserve(history.nodes().size()) ;
    PlanarNode root ; root.geometry=history.geometry().map(history.node(0).map) ;
    for(unsigned i=0; i<rootMap->corners().size(); ++i) {
      const auto vertex=cornerPoint(i) ; root.corners.push_back(vertex) ;
      root.edges.push_back({{true,{},i},Fraction(0),Fraction(1)}) ;
      points_.emplace(vertex,rootMap->corners()[i]) ;
    }
    const auto roots=history.children(0) ; root.children.assign(roots.begin(),roots.end()) ;
    nodes_.push_back(std::move(root)) ;
    while(nodes_.size()<history.nodes().size()) {
      const auto &split=history.splits().at(owners[nodes_.size()]) ;
      PlanarSplitKey owner{history.key(split.parent),split.construction} ;
      const auto &accepted=*split.construction->planar ;
      points_.emplace(centerPoint(owner),accepted.points.back()) ;
      auto children=childNodes(split.parent,owner) ;
      for(auto &child : children) {
        const auto id=nodes_.size() ; child.geometry=history.geometry().map(history.node(id).map) ;
        const auto selected=history.children(id) ; child.children.assign(selected.begin(),selected.end()) ;
        nodes_.push_back(std::move(child)) ;
      }
      const auto &parent=nodes_[split.parent] ;
      for(unsigned edge=0; edge<parent.edges.size(); ++edge) {
        const auto &use=parent.edges[edge] ;
        if(evaluate(vertexAt(use.edge,(use.begin+use.end)/Fraction(2))) != accepted.points[parent.corners.size()+edge])
          throw std::invalid_argument("accepted planar midpoint differs from its retained edge") ;
      }
    }
  }

  TreeNodeId PlanarSubdivision::nodeId(const CellPath &path) const {
    TreeNodeId id=0 ;
    for(unsigned ordinal : path) {
      const auto &children=nodes_[id].children ;
      if(ordinal >= children.size()) throw std::out_of_range("unknown planar subdivision path") ;
      id=children[ordinal] ;
    }
    return id ;
  }

  std::vector<TreeNodeId> PlanarSubdivision::leafIds() const { return history().leafIds() ; }

  TreeNodeId PlanarSubdivision::ownerNode(const PlanarSplitKey &owner) const {
    if(!owner.construction || !owner.construction->planar)
      throw std::invalid_argument("planar edge has no accepted construction") ;
    CellKey child=owner.parent ; child.lineage.push_back({owner.construction,0}) ;
    return history().node(history().nodeId(child)).parent ;
  }

  PlanarVertexKey PlanarSubdivision::vertexAt(const PlanarEdgeKey &edge, Fraction position) const {
    ExactBudget budget(history().limits().exact) ; return vertexAt(edge,std::move(position),budget) ;
  }

  PlanarVertexKey PlanarSubdivision::vertexAt(const PlanarEdgeKey &edge, Fraction position, ExactBudget &budget) const {
    budget.operations(1+edge.owner.parent.lineage.size()) ; budget.check(position) ;
    if(position < Fraction(0) || position > Fraction(1))
      throw std::invalid_argument("planar vertex position is outside its edge") ;
    if(edge.original) {
      if(edge.owner.construction || !edge.owner.parent.lineage.empty() || edge.owner.parent.root != 0 ||
         edge.index >= nodes_[0].corners.size())
        throw std::invalid_argument("invalid original planar subdivision edge") ;
      if(position == Fraction(0)) return nodes_[0].corners[edge.index] ;
      if(position == Fraction(1)) return nodes_[0].corners[(edge.index+1)%nodes_[0].corners.size()] ;
    } else {
      const auto &owner=nodes_.at(ownerNode(edge.owner)) ;
      if(edge.index >= owner.edges.size()) throw std::invalid_argument("planar subdivision spoke does not exist") ;
      if(position == Fraction(0)) return centerPoint(edge.owner) ;
      if(position == Fraction(1)) {
        const auto &use=owner.edges[edge.index] ;
        return vertexAt(use.edge,(use.begin+use.end)/Fraction(2),budget) ;
      }
    }
    PlanarVertexKey result ; result.kind=PlanarVertexKey::Kind::edge ; result.edge=edge ; result.position=position ;
    return result ;
  }

  Point PlanarSubdivision::evaluate(const PlanarVertexKey &vertex) const {
    ExactBudget budget(history().limits().exact) ; return evaluate(vertex,budget) ;
  }

  Point PlanarSubdivision::evaluate(const PlanarVertexKey &vertex, ExactBudget &budget) const {
    budget.operations() ;
    const auto stored=points_.find(vertex) ;
    if(stored != points_.end()) return stored->second ;
    if(vertex.kind != PlanarVertexKey::Kind::edge) throw std::out_of_range("unknown planar subdivision point") ;
    const auto canonical=vertexAt(vertex.edge,vertex.position,budget) ;
    if(!(canonical == vertex)) return evaluate(canonical,budget) ;
    return interpolate(evaluate(vertexAt(vertex.edge,Fraction(0),budget),budget),
                       evaluate(vertexAt(vertex.edge,Fraction(1),budget),budget),vertex.position) ;
  }

  PlanarSplit PlanarSubdivision::proposedSplit(TreeNodeId parent, CenterMethod method, MeasureOptions options) const {
    std::vector<Point> corners, midpoints ;
    for(const auto &vertex : node(parent).corners) corners.push_back(evaluate(vertex)) ;
    for(const auto &edge : node(parent).edges)
      midpoints.push_back(evaluate(vertexAt(edge.edge,(edge.begin+edge.end)/Fraction(2)))) ;
    return splitPlanarPolygon(corners,midpoints,method,options) ;
  }

  std::vector<PlanarNode> PlanarSubdivision::childNodes(TreeNodeId parentId, const PlanarSplitKey &owner) const {
    const auto &parent=node(parentId) ; const unsigned count=parent.corners.size() ;
    std::vector<PlanarNode> children ; children.reserve(count) ;
    std::vector<PlanarVertexKey> midpoints ; std::vector<Fraction> cuts ;
    for(const auto &edge : parent.edges) {
      const auto cut=(edge.begin+edge.end)/Fraction(2) ;
      cuts.push_back(cut) ; midpoints.push_back(vertexAt(edge.edge,cut)) ;
    }
    const auto center=centerPoint(owner) ;
    for(unsigned i=0; i<count; ++i) {
      const unsigned previous=(i+count-1)%count ; PlanarNode child ;
      child.parent=parentId ; child.path=parent.path ; child.path.push_back(i) ;
      child.corners={parent.corners[i],midpoints[i],center,midpoints[previous]} ;
      child.edges={
        {parent.edges[i].edge,parent.edges[i].begin,cuts[i]},
        {{false,owner,i},Fraction(1),Fraction(0)},
        {{false,owner,previous},Fraction(0),Fraction(1)},
        {parent.edges[previous].edge,cuts[previous],parent.edges[previous].end}
      } ;
      children.push_back(std::move(child)) ;
    }
    return children ;
  }

  std::vector<TreeNodeId> PlanarSubdivision::splitLeaf(TreeNodeId leaf, CenterMethod method, MeasureOptions options) {
    static_assert(std::is_nothrow_move_constructible<PlanarNode>::value,"planar publication requires nonthrowing moves") ;
    if(!node(leaf).isLeaf() || !history().active(leaf))
      throw std::invalid_argument("planar subdivision requires an active leaf") ;
    auto split=std::make_shared<const PlanarSplit>(proposedSplit(leaf,method,options)) ;
    auto construction=std::make_shared<const SplitConstruction>(SplitConstruction{split->shape,SplitFamily::planar,1,{},nullptr,split}) ;
    PlanarSplitKey owner{history().key(leaf),construction} ;
    auto additions=childNodes(leaf,owner) ;
    std::map<PlanarVertexKey,Point> points ; points.emplace(centerPoint(owner),split->points.back()) ;
    std::vector<TreeNodeId> ids(additions.size()), result(additions.size()) ;
    const std::size_t maximum=std::min(maxNodes(),nodes_.max_size()) ;
    const std::size_t needed=nodes_.size()+std::min(additions.size(),maximum-nodes_.size()) ;
    if(nodes_.capacity()<needed) {
      const auto grown=nodes_.capacity()>maximum/2 ? maximum : 2*nodes_.capacity() ;
      nodes_.reserve(std::max(needed,grown)) ;
    }
    edit_.refine(leaf,*split) ;
    const auto accepted=history().children(leaf) ;
    for(unsigned ordinal=0; ordinal<accepted.size(); ++ordinal) {
      ids[ordinal]=result[ordinal]=accepted[ordinal] ;
      if(accepted[ordinal] >= nodes_.size()) {
        additions[ordinal].geometry=history().geometry().map(history().node(accepted[ordinal]).map) ;
        nodes_.push_back(std::move(additions[ordinal])) ;
      }
    }
    points_.merge(points) ; nodes_[leaf].children=std::move(ids) ; return result ;
  }

  CellChange PlanarSubdivision::coarsen(TreeNodeId parent) {
    auto result=edit_.coarsen(parent) ; nodes_[parent].children.clear() ; return result ;
  }

  CellChange PlanarSubdivision::replace(TreeNodeId parent, CenterMethod method, MeasureOptions options) {
    HistoryEdit candidate(history()) ; auto change=candidate.replace(parent,proposedSplit(parent,method,options)) ;
    PlanarSubdivision accepted(std::move(candidate).finish()) ; swap(accepted) ; return change ;
  }

  CellChange PlanarSubdivision::reactivate(TreeNodeId parent, const SplitConstruction &construction) {
    const auto retained = history().retainedChildren(parent,construction) ;
    std::vector<TreeNodeId> children(retained.begin(),retained.end()) ;
    auto change = edit_.reactivate(parent,construction) ;
    nodes_[parent].children = std::move(children) ; return change ;
  }

  void PlanarSubdivision::compact() { PlanarSubdivision accepted(compactHistory(history())) ; swap(accepted) ; }

  void PlanarSubdivision::restoreHistory(RefinementHistory saved) {
    if(saved.rootKey() != history().rootKey() || saved.geometryRevision() != history().geometryRevision())
      throw std::invalid_argument("planar history belongs to a different root or geometry revision") ;
    auto accepted=RefinementHistory::restore(saved.rootKey(),saved.geometryRevision(),saved.geometry(),
      saved.nodes(),saved.splits(),history().limits()) ;
    const auto *source=dynamic_cast<const CornerMap*>(node(0).geometry.get()) ;
    const auto *target=dynamic_cast<const CornerMap*>(accepted.geometry().map(accepted.node(0).map).get()) ;
    if(!source || !target || source->shape() != target->shape() || source->corners() != target->corners())
      throw std::invalid_argument("planar history has different root geometry") ;
    PlanarSubdivision result(std::move(accepted)) ; swap(result) ;
  }

  void PlanarSubdivision::swap(PlanarSubdivision &other) noexcept {
    edit_.swap(other.edit_) ; nodes_.swap(other.nodes_) ; points_.swap(other.points_) ;
  }

}
}
