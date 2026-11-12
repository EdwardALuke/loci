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

#include <FVMAdapt/parametric/affine_subdivision.h>
#include "subdivision_internal.h"

#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace Loci {
namespace FVMAdapt {
  namespace {

    template<class T> void reserveGrowth(std::vector<T> &values, std::size_t needed,
                                        std::size_t limit) {
      if(values.capacity() >= needed) return ;
      limit = std::min(limit,values.max_size()) ;
      const std::size_t grown = values.capacity() > limit/2 ? limit : 2*values.capacity() ;
      values.reserve(std::max(needed,grown)) ;
    }

  }

  namespace {
    RefinementHistory referenceHistory(Shape shape, std::size_t maxNodes) {
      const auto &element = referenceElement(shape) ;
      GeometryStore geometry ;
      const MapId root = geometry.append(CornerRecord{shape,
        std::vector<Point>(element.vertices.begin(),element.vertices.begin()+element.vertexCount)}) ;
      HistoryLimits limits ; limits.maxNodes = maxNodes ; limits.maxSplits = maxNodes ;
      return RefinementHistory(0,std::move(geometry),root,0,limits) ;
    }
  }

  AffineSubdivision::AffineSubdivision(Shape shape, std::size_t maxNodes)
    : AffineSubdivision(referenceHistory(shape,maxNodes)) {}

  AffineSubdivision::AffineSubdivision(RefinementHistory history) : edit_(std::move(history)) {
    nodes_.reserve(this->history().nodes().size()) ;
    for(TreeNodeId id = 0; id < this->history().nodes().size(); ++id) {
      const CellNode &cell = this->history().node(id) ;
      nodes_.push_back(AffineNode(cell.parent,{this->history().localCell(id).shape,this->history().rootPlacement(id)})) ;
    }
  }

  void AffineSubdivision::checkLeaf(TreeNodeId leaf) const {
    if(!isLeaf(leaf) || !history().active(leaf))
      throw std::invalid_argument("subdivision node is not an active leaf") ;
  }

  std::vector<TreeNodeId> AffineSubdivision::leafIds() const {
    auto result = history().leafIds() ; std::sort(result.begin(),result.end()) ; return result ;
  }

  CellChange AffineSubdivision::append(TreeNodeId leaf, std::vector<Child> children, SplitOperation operation) {
    static_assert(std::is_nothrow_move_constructible<AffineNode>::value,
                  "root-coordinate publication requires nonthrowing node moves") ;
    std::vector<AffineNode> additions ; additions.reserve(children.size()) ;
    for(const Child &child : children) {
      Child placed = {child.shape,nodes_[leaf].cell.placement.compose(child.placement)} ;
      for(const Fraction &value : placed.placement.origin) checkExact(value,history().limits().exact) ;
      for(const auto &column : placed.placement.columns)
        for(const Fraction &value : column) checkExact(value,history().limits().exact) ;
      additions.push_back(AffineNode(leaf,std::move(placed))) ;
    }
    // A retained construction needs no new cache entries, even at the node limit.
    if(children.size() <= maxNodes()-nodes_.size())
      reserveGrowth(nodes_,nodes_.size()+children.size(),maxNodes()) ;
    CellChange change = edit_.append(leaf,std::move(operation),std::move(children)) ;
    if(history().nodes().size() != nodes_.size())
      for(auto &node : additions) nodes_.push_back(std::move(node)) ;
    return change ;
  }

  namespace {
    std::vector<Child> construct(Shape shape, const SplitOperation &request) {
      return std::visit([&](const auto &operation) {
        using T = std::decay_t<decltype(operation)> ;
        if constexpr(std::is_same_v<T,TensorSplit>) return detail::splitChildren(shape,operation) ;
        else return detail::splitChildren(operation) ;
      },request) ;
    }
  }

  CellChange AffineSubdivision::refine(TreeNodeId leaf, const SplitOperation &request) {
    checkLeaf(leaf) ;
    if(edit_.preflight(leaf,request) == 1) return edit_.refine(leaf,request) ;
    return append(leaf,construct(node(leaf).cell.shape,request),request) ;
  }

  CellChange AffineSubdivision::coarsen(TreeNodeId parent) {
    return edit_.coarsen(parent) ;
  }

  CellChange AffineSubdivision::replace(TreeNodeId parent, const SplitOperation &request) {
    if(edit_.preflight(parent,request) == 1)
      throw std::invalid_argument("replacement requires a nonempty split") ;
    return append(parent,construct(node(parent).cell.shape,request),request) ;
  }

  CellChange AffineSubdivision::reactivate(TreeNodeId parent, const SplitConstruction &construction) {
    return edit_.reactivate(parent,construction) ;
  }

  void AffineSubdivision::compact() {
    AffineSubdivision next(compactHistory(history())) ; swap(next) ;
  }

  void AffineSubdivision::swap(AffineSubdivision &other) noexcept {
    edit_.swap(other.edit_) ; nodes_.swap(other.nodes_) ;
  }

  NodeRange AffineSubdivision::splitLeaf(TreeNodeId leaf, const TensorSplit &request) {
    return refine(leaf,request).changed() ? children(leaf) : NodeRange() ;
  }

  NodeRange AffineSubdivision::splitLeaf(TreeNodeId leaf, const TriangleSplit &request) {
    return refine(leaf,request).changed() ? children(leaf) : NodeRange() ;
  }

  NodeRange AffineSubdivision::splitLeaf(TreeNodeId leaf, const PrismSplit &request) {
    return refine(leaf,request).changed() ? children(leaf) : NodeRange() ;
  }

  NodeRange AffineSubdivision::splitLeaf(TreeNodeId leaf, const TetSplit &request) {
    return refine(leaf,request).changed() ? children(leaf) : NodeRange() ;
  }

  NodeRange AffineSubdivision::splitLeaf(TreeNodeId leaf, const PyramidSplit &request) {
    return refine(leaf,request).changed() ? children(leaf) : NodeRange() ;
  }

}
}
