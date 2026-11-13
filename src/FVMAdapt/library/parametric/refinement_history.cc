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

#include <FVMAdapt/parametric/refinement_history.h>
#include "subdivision_internal.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace Loci {
namespace FVMAdapt {
  namespace {
    void checkShape(Shape shape) {
      if(referenceElement(shape).dimension == 0)
        throw std::invalid_argument("refinement requires a non-point root") ;
    }
    bool samePlacement(const Placement &a, const Placement &b) {
      return a.origin == b.origin && a.columns == b.columns ;
    }
    void checkPlacement(const Placement &placement, ExactBudget &budget) {
      for(const auto &value : placement.origin) budget.check(value) ;
      for(const auto &column : placement.columns) for(const auto &value : column) budget.check(value) ;
    }
    std::size_t product(std::size_t a, std::size_t b) {
      if(b && a > std::numeric_limits<std::size_t>::max()/b)
        throw std::length_error("refinement count exceeds addressable size") ;
      return a*b ;
    }
    void checkCuts(const std::vector<Fraction> &cuts, ExactBudget &budget) {
      budget.operations(cuts.size()) ;
      Fraction previous = 0 ;
      for(const Fraction &cut : cuts) {
        budget.check(cut) ;
        if(cut <= previous || cut >= Fraction(1))
          throw std::invalid_argument("refinement cuts must increase strictly inside (0,1)") ;
        previous = cut ;
      }
    }
    unsigned triangleCount(const TriangleSplit &request, ExactBudget &budget) {
      switch(request.kind) {
      case TriangleSplit::Kind::none: return 1 ;
      case TriangleSplit::Kind::uniform: return 4 ;
      case TriangleSplit::Kind::bisect:
        budget.check(request.position) ;
        if(request.edge >= 3 || request.position <= Fraction(0) || request.position >= Fraction(1))
          throw std::invalid_argument("invalid triangle bisection") ;
        return 2 ;
      }
      throw std::invalid_argument("unknown triangle split") ;
    }
    void normalizeRequest(TriangleSplit &request) {
      if(request.kind == TriangleSplit::Kind::none) request = TriangleSplit() ;
      else if(request.kind == TriangleSplit::Kind::uniform) request = TriangleSplit::uniform4() ;
    }
    void normalizeRequest(SplitOperation &operation) {
      std::visit([](auto &request) {
        using Request = std::decay_t<decltype(request)> ;
        if constexpr(std::is_same_v<Request,TriangleSplit>) normalizeRequest(request) ;
        else if constexpr(std::is_same_v<Request,PrismSplit>) normalizeRequest(request.base) ;
        else if constexpr(std::is_same_v<Request,TetSplit>) {
          if(request.kind == TetSplit::Kind::none) request = TetSplit() ;
          else if(request.kind == TetSplit::Kind::uniform) request = TetSplit::uniform8(request.diagonal) ;
          else request.diagonal = 0 ;
        }
      },operation) ;
    }
    std::size_t requestCount(Shape shape, const SplitOperation &operation, ExactBudget &budget) {
      return std::visit([&](const auto &request)->std::size_t {
        using Request = std::decay_t<decltype(request)> ;
        if constexpr(std::is_same_v<Request,TensorSplit>) {
          if(shape != Shape::segment && shape != Shape::quad && shape != Shape::hex)
            throw std::invalid_argument("tensor split requires a segment, quad or hex") ;
          const unsigned dimension = referenceElement(shape).dimension ;
          std::size_t count = 1 ;
          for(unsigned d = 0; d < 3; ++d) {
            if(d >= dimension && !request.cuts[d].empty()) throw std::invalid_argument("cuts in an unused direction") ;
            checkCuts(request.cuts[d],budget) ;
            if(request.cuts[d].size() == std::numeric_limits<std::size_t>::max())
              throw std::length_error("refinement cut count exceeds addressable size") ;
            count = product(count,request.cuts[d].size()+1) ;
          }
          return count ;
        } else if constexpr(std::is_same_v<Request,TriangleSplit>) {
          if(shape != Shape::triangle) throw std::invalid_argument("triangle split requires a triangle") ;
          return triangleCount(request,budget) ;
        } else if constexpr(std::is_same_v<Request,PrismSplit>) {
          if(shape != Shape::prism) throw std::invalid_argument("prism split requires a prism") ;
          checkCuts(request.axialCuts,budget) ;
          if(request.axialCuts.size() == std::numeric_limits<std::size_t>::max())
            throw std::length_error("refinement cut count exceeds addressable size") ;
          return product(triangleCount(request.base,budget),request.axialCuts.size()+1) ;
        } else if constexpr(std::is_same_v<Request,TetSplit>) {
          if(shape != Shape::tet) throw std::invalid_argument("tet split requires a tet") ;
          switch(request.kind) {
          case TetSplit::Kind::none: return 1 ;
          case TetSplit::Kind::uniform:
            if(request.diagonal > 2) throw std::invalid_argument("invalid tet diagonal") ;
            return 8 ;
          case TetSplit::Kind::bisect:
            budget.check(request.position) ;
            if(request.edge >= 6 || request.position <= Fraction(0) || request.position >= Fraction(1))
              throw std::invalid_argument("invalid tet bisection") ;
            return 2 ;
          }
          throw std::invalid_argument("unknown tet split") ;
        } else {
          if(shape != Shape::pyramid) throw std::invalid_argument("pyramid split requires a pyramid") ;
          switch(request.kind) {
          case PyramidSplit::Kind::none: return 1 ;
          case PyramidSplit::Kind::uniform: return 10 ;
          }
          throw std::invalid_argument("unknown pyramid split") ;
        }
      },operation) ;
    }
    SplitFamily family(const SplitOperation &operation) {
      return std::visit([](const auto &request) {
        using Request = std::decay_t<decltype(request)> ;
        if constexpr(std::is_same_v<Request,TensorSplit>) return SplitFamily::tensor ;
        else if constexpr(std::is_same_v<Request,TriangleSplit>) return SplitFamily::triangle ;
        else if constexpr(std::is_same_v<Request,PrismSplit>) return SplitFamily::prism ;
        else if constexpr(std::is_same_v<Request,TetSplit>) return SplitFamily::tet ;
        else return SplitFamily::pyramid ;
      },operation) ;
    }
    std::vector<Child> construct(Shape shape, const SplitOperation &operation) {
      return std::visit([&](const auto &request) {
        if constexpr(std::is_same_v<std::decay_t<decltype(request)>,TensorSplit>)
          return detail::splitChildren(shape,request) ;
        else return detail::splitChildren(request) ;
      },operation) ;
    }
    template<class T> void reserveGrowth(std::vector<T> &values, std::size_t needed, std::size_t limit) {
      if(values.capacity() >= needed) return ;
      limit = std::min(limit,values.max_size()) ;
      if(needed > limit) throw std::length_error("refinement table exceeds its limit") ;
      const auto grown = values.capacity() > limit/2 ? limit : 2*values.capacity() ;
      values.reserve(std::max(needed,grown)) ;
    }
    std::size_t integerBytes(Fraction::Integer value) {
      if(value < 0) value = -value ;
      return value == 0 ? 1 : boost::multiprecision::msb(value)/8+1 ;
    }
    std::size_t placementBytes(const Placement &placement) {
      std::size_t bytes = 0 ;
      const auto count = [&](const Fraction &value) {
        return integerBytes(value.numerator())+integerBytes(value.denominator()) ;
      } ;
      for(const auto &value : placement.origin) bytes += count(value) ;
      for(const auto &column : placement.columns) for(const auto &value : column) bytes += count(value) ;
      return bytes ;
    }
    std::size_t operationBytes(const SplitOperation &operation) {
      return std::visit([](const auto &request) {
        using Request = std::decay_t<decltype(request)> ;
        std::size_t bytes = 0 ;
        const auto count = [&](const Fraction &value) {
          return integerBytes(value.numerator())+integerBytes(value.denominator()) ;
        } ;
        if constexpr(std::is_same_v<Request,TensorSplit>) {
          for(const auto &cuts : request.cuts) for(const auto &cut : cuts) bytes += sizeof(Fraction)+count(cut) ;
        } else if constexpr(std::is_same_v<Request,PrismSplit>) {
          for(const auto &cut : request.axialCuts) bytes += sizeof(Fraction)+count(cut) ;
          bytes += count(request.base.position) ;
        } else if constexpr(std::is_same_v<Request,TriangleSplit> || std::is_same_v<Request,TetSplit>) {
          bytes += count(request.position) ;
        }
        return bytes ;
      },operation) ;
    }
    std::size_t mappedBytes(const VertexVolumeSplit &split) {
      return sizeof(VertexVolumeSplit)+split.points.size()*sizeof(Point)+
        split.children.size()*(sizeof(VertexVolumeChild)+8*sizeof(Point))+
        split.internalFaces.size()*sizeof(InternalFace)+split.boundaryPieces.size()*sizeof(VertexVolumeBoundary) ;
    }
    std::shared_ptr<const VertexVolumeSplit> validatedMapped(const VertexVolumeSplit &split,
                                                           JacobianCheckLimits limits) {
      if(split.shape != Shape::hex && split.shape != Shape::prism && split.shape != Shape::tet)
        throw std::invalid_argument("vertex volume subdivision requires a hex, prism or tet") ;
      const auto &element=referenceElement(split.shape) ;
      if(split.children.size() != element.vertexCount ||
         split.points.size() != element.vertexCount+element.edgeCount+element.boundaryCount+1)
        throw std::invalid_argument("invalid accepted vertex volume counts") ;
      std::vector<std::array<unsigned,8>> children ;
      for(const auto &child : split.children) children.push_back(child.vertices) ;
      auto result = restoreVertexVolumeSplit(split.shape,split.points,children,limits) ;
      for(std::size_t i=0; i<children.size(); ++i) {
        const auto *actual = dynamic_cast<const CornerMap*>(split.children[i].placement.get()) ;
        const auto &expected = static_cast<const CornerMap&>(*result.children[i].placement) ;
        if(!actual || actual->shape() != Shape::hex || actual->corners() != expected.corners())
          throw std::invalid_argument("accepted vertex volume chart differs from its corners") ;
      }
      if(split.internalFaces.size() != result.internalFaces.size() || split.boundaryPieces.size() != result.boundaryPieces.size())
        throw std::invalid_argument("accepted vertex volume has inconsistent face counts") ;
      for(std::size_t i=0; i<split.internalFaces.size(); ++i) {
        const auto &a=split.internalFaces[i], &b=result.internalFaces[i] ;
        if(a.first.child != b.first.child || a.first.boundary != b.first.boundary ||
           a.second.child != b.second.child || a.second.boundary != b.second.boundary || a.orientation != b.orientation)
          throw std::invalid_argument("accepted vertex volume has inconsistent interior faces") ;
      }
      for(std::size_t i=0; i<split.boundaryPieces.size(); ++i) {
        const auto &a=split.boundaryPieces[i], &b=result.boundaryPieces[i] ;
        if(a.face.child != b.face.child || a.face.boundary != b.face.boundary || a.boundary != b.boundary ||
           a.corner != b.corner || a.orientation != b.orientation)
          throw std::invalid_argument("accepted vertex volume has inconsistent boundary faces") ;
      }
      return std::make_shared<const VertexVolumeSplit>(std::move(result)) ;
    }
    std::size_t planarBytes(const PlanarSplit &split) {
      return sizeof(PlanarSplit)+split.points.size()*sizeof(Point)+split.children.size()*sizeof(std::array<unsigned,4>) ;
    }
    std::shared_ptr<const PlanarSplit> validatedPlanar(const PlanarSplit &split, const GeometryMap &parent) {
      const auto *map=dynamic_cast<const CornerMap*>(&parent) ;
      if(!map || (map->shape() != Shape::triangle && map->shape() != Shape::quad))
        throw std::invalid_argument("planar subdivision requires a physical triangle or quad map") ;
      validatePlanarSplit(split,map->corners()) ;
      return std::make_shared<const PlanarSplit>(split) ;
    }
    std::size_t addBytes(std::size_t bytes, std::size_t addition, const ExactLimits &limits) {
      if(bytes > limits.maxBytes || addition > limits.maxBytes-bytes)
        throw std::length_error("refinement exceeds retained payload limit") ;
      return bytes+addition ;
    }

  }

  RefinementHistory::RefinementHistory(std::uint64_t root, GeometryStore geometry, MapId rootMap,
                                       std::uint64_t revision, HistoryLimits limits)
    : rootKey_(root), geometryRevision_(revision), limits_(limits),
      rootCell_(Child{geometry.map(rootMap)->shape(),Placement()}), geometry_(limits.geometry) {
    if(limits.maxNodes == 0) throw std::invalid_argument("refinement node limit must include the root") ;
    checkShape(rootCell_.shape) ;
    ExactBudget budget(limits.exact) ; checkPlacement(rootCell_.placement,budget) ;
    retainedBytes_ = addBytes(0,sizeof(Child)+placementBytes(rootCell_.placement)+
                             sizeof(CellNode)+2*sizeof(std::size_t),limits.exact) ;
    geometry.setLimits(limits.geometry,limits.exact.maxBytes-retainedBytes_) ;
    geometry_.swap(geometry) ;
    nodes_.push_back({noTreeNode,0,rootMap,noTreeNode}) ;
    owners_.push_back(noTreeNode) ; firstSplit_.push_back(noTreeNode) ;
  }

  const Child &RefinementHistory::localCell(TreeNodeId id) const {
    const CellNode &cell = node(id) ;
    if(id != 0 && (splits_.at(owners_.at(id)).construction->vertexVolume ||
                   splits_.at(owners_.at(id)).construction->planar))
      throw std::invalid_argument("mapped cell has no affine placement") ;
    return id == 0 ? rootCell_ : splits_.at(owners_.at(id)).construction->children.at(cell.ordinal) ;
  }

  std::size_t RefinementHistory::splitId(TreeNodeId parent, const SplitConstruction &construction) const {
    node(parent) ;
    for(std::size_t id = firstSplit_[parent]; id != noTreeNode; id = nextSplit_[id])
      if(*splits_[id].construction == construction) return id ;
    return noTreeNode ;
  }

  NodeRange RefinementHistory::retainedChildren(TreeNodeId parent, const SplitConstruction &construction) const {
    const auto id = splitId(parent,construction) ;
    return id == noTreeNode ? NodeRange() : splits_[id].children ;
  }

  NodeRange RefinementHistory::children(TreeNodeId id) const {
    const CellNode &cell = node(id) ;
    return cell.isLeaf() ? NodeRange() : splits_.at(cell.activeSplit).children ;
  }

  const SplitOperation *RefinementHistory::operation(TreeNodeId id) const {
    const CellNode &cell = node(id) ;
    if(cell.isLeaf()) return nullptr ;
    const auto &operation = splits_.at(cell.activeSplit).operation ;
    return operation ? &*operation : nullptr ;
  }

  bool RefinementHistory::active(TreeNodeId id) const {
    node(id) ;
    while(id != 0) {
      const CellNode &current = nodes_[id] ;
      const auto selected = children(current.parent) ;
      if(current.ordinal >= selected.size() || selected[current.ordinal] != id) return false ;
      id = current.parent ;
    }
    return true ;
  }

  std::vector<TreeNodeId> RefinementHistory::leafIds() const {
    std::vector<TreeNodeId> result, pending = {0} ;
    while(!pending.empty()) {
      const TreeNodeId id = pending.back() ; pending.pop_back() ;
      const auto selected = children(id) ;
      if(selected.empty()) result.push_back(id) ;
      else for(std::size_t c = selected.size(); c > 0; --c) pending.push_back(selected[c-1]) ;
    }
    return result ;
  }

  CellPath RefinementHistory::path(TreeNodeId id) const {
    node(id) ; CellPath result ;
    while(id != 0) { result.push_back(nodes_[id].ordinal) ; id = nodes_[id].parent ; }
    std::reverse(result.begin(),result.end()) ; return result ;
  }

  CellKey RefinementHistory::key(TreeNodeId id) const {
    node(id) ; CellKey result = {rootKey_,{}} ;
    while(id != 0) {
      const CellNode &current = nodes_[id] ;
      result.lineage.push_back({splits_[owners_[id]].construction,current.ordinal}) ; id = current.parent ;
    }
    std::reverse(result.lineage.begin(),result.lineage.end()) ; return result ;
  }

  TreeNodeId RefinementHistory::nodeId(const CellKey &key) const {
    if(key.root != rootKey_) throw std::out_of_range("cell key belongs to another root") ;
    TreeNodeId node = 0 ;
    for(const CellStep &step : key.lineage) {
      if(!step.construction) throw std::out_of_range("cell key has no accepted construction") ;
      const auto split = splitId(node,*step.construction) ;
      if(split == noTreeNode || step.child >= splits_[split].children.size())
        throw std::out_of_range("cell key is not retained in this history") ;
      node = splits_[split].children[step.child] ;
    }
    return node ;
  }

  Placement RefinementHistory::rootPlacement(TreeNodeId id) const {
    node(id) ;
    std::vector<TreeNodeId> ancestors ;
    while(id != 0) { ancestors.push_back(id) ; id = nodes_[id].parent ; }
    Placement result ; ExactBudget budget(limits_.exact) ;
    for(auto it = ancestors.rbegin(); it != ancestors.rend(); ++it) {
      budget.operations(84) ; result = result.compose(localCell(*it).placement) ; checkPlacement(result,budget) ;
    }
    return result ;
  }

  void RefinementHistory::swap(RefinementHistory &other) noexcept {
    static_assert(std::is_nothrow_swappable<Child>::value,"root placement must support a no-throw swap") ;
    using std::swap ;
    swap(rootKey_,other.rootKey_) ; swap(geometryRevision_,other.geometryRevision_) ;
    swap(limits_,other.limits_) ; swap(rootCell_,other.rootCell_) ; geometry_.swap(other.geometry_) ;
    nodes_.swap(other.nodes_) ; splits_.swap(other.splits_) ; owners_.swap(other.owners_) ;
    firstSplit_.swap(other.firstSplit_) ; nextSplit_.swap(other.nextSplit_) ;
    swap(retainedBytes_,other.retainedBytes_) ;
  }

  RefinementHistory RefinementHistory::restore(std::uint64_t root, std::uint64_t revision,
                                               GeometryStore geometry, std::vector<CellNode> nodes,
                                               std::vector<AcceptedSplit> splits, HistoryLimits limits) {
    return restoreRecords(root,revision,std::move(geometry),std::move(nodes),std::move(splits),limits,false) ;
  }

  RefinementHistory RefinementHistory::restoreRecords(std::uint64_t root, std::uint64_t revision,
      GeometryStore geometry, std::vector<CellNode> nodes, std::vector<AcceptedSplit> splits,
      HistoryLimits limits, bool mappedValidated) {
    if(nodes.empty()) throw std::invalid_argument("refinement record has no root") ;
    if(nodes.size() > limits.maxNodes || splits.size() > limits.maxSplits)
      throw std::length_error("refinement record exceeds table limits") ;
    if(nodes[0].parent != noTreeNode || nodes[0].ordinal != 0)
      throw std::invalid_argument("invalid refinement root") ;
    RefinementHistory result(root,std::move(geometry),nodes[0].map,revision,limits) ;
    std::size_t retained = addBytes(sizeof(Child)+placementBytes(result.rootCell_.placement),
      product(nodes.size(),sizeof(CellNode)+2*sizeof(std::size_t)),limits.exact) ;
    addBytes(retained,result.geometry_.retainedBytes(),limits.exact) ;
    std::vector<unsigned> depths(nodes.size(),0) ;
    std::vector<std::size_t> owners(nodes.size(),noTreeNode), first(nodes.size(),noTreeNode), next(splits.size(),noTreeNode) ;
    for(TreeNodeId id = 0; id < nodes.size(); ++id) {
      const CellNode &cell = nodes[id] ;
      checkShape(result.geometry_.map(cell.map)->shape()) ;
      if(id != 0) {
        if(cell.parent >= id) throw std::invalid_argument("refinement parent must precede child") ;
        if(depths[cell.parent] >= limits.maxDepth) throw std::length_error("refinement record exceeds depth limit") ;
        depths[id] = depths[cell.parent]+1 ;
      }
      if(cell.activeSplit != noTreeNode &&
         (cell.activeSplit >= splits.size() || splits[cell.activeSplit].parent != id))
        throw std::invalid_argument("refinement active split has the wrong parent") ;
    }
    for(std::size_t index = 0; index < splits.size(); ++index) {
      AcceptedSplit &split = splits[index] ;
      if(split.operation) normalizeRequest(*split.operation) ;
      if(split.parent >= nodes.size() || split.children.size() < 2 ||
         split.children.size() > std::numeric_limits<unsigned>::max() ||
         split.children.front() <= split.parent || split.children.front() > nodes.size() ||
         split.children.size() > nodes.size()-split.children.front())
        throw std::invalid_argument("invalid refinement child range") ;
      const auto parentShape=result.geometry_.map(nodes[split.parent].map)->shape() ;
      ExactBudget budget(limits.exact) ;
      if(!split.construction || split.construction->version != 1 ||
         split.construction->parentShape != parentShape || split.construction->childCount() != split.children.size())
        throw std::invalid_argument("invalid accepted refinement construction") ;
      SplitConstruction construction = *split.construction ;
      retained = addBytes(retained,sizeof(AcceptedSplit)+sizeof(SplitConstruction)+sizeof(std::size_t),limits.exact) ;
      if(construction.family == SplitFamily::vertexVolume) {
        if(split.operation || !construction.children.empty() || construction.planar || !construction.vertexVolume ||
           construction.vertexVolume->shape != parentShape)
          throw std::invalid_argument("invalid accepted mapped construction") ;
        if(!mappedValidated) construction.vertexVolume = validatedMapped(*construction.vertexVolume,limits.jacobian) ;
        retained = addBytes(retained,mappedBytes(*construction.vertexVolume),limits.exact) ;
      } else if(construction.family == SplitFamily::planar) {
        if(split.operation || !construction.children.empty() || construction.vertexVolume || !construction.planar)
          throw std::invalid_argument("invalid accepted planar construction") ;
        construction.planar=validatedPlanar(*construction.planar,*result.geometry_.map(nodes[split.parent].map)) ;
        retained=addBytes(retained,planarBytes(*construction.planar),limits.exact) ;
      } else {
        if(!split.operation || construction.vertexVolume || construction.planar || construction.family != family(*split.operation) ||
           requestCount(parentShape,*split.operation,budget) != split.children.size())
          throw std::invalid_argument("invalid accepted native construction") ;
        retained = addBytes(retained,operationBytes(*split.operation),limits.exact) ;
        for(const Child &child : construction.children) {
          checkPlacement(child.placement,budget) ;
          retained = addBytes(retained,sizeof(Child)+placementBytes(child.placement),limits.exact) ;
        }
        detail::validateSubdivision(parentShape,construction.children,*split.operation,budget) ;
      }
      for(std::size_t ordinal = 0; ordinal < split.children.size(); ++ordinal) {
        const TreeNodeId id = split.children[ordinal] ;
        if(owners[id] != noTreeNode || nodes[id].parent != split.parent || nodes[id].ordinal != ordinal)
          throw std::invalid_argument("refinement child has inconsistent ownership") ;
        owners[id] = index ;
      }
      split.construction = std::make_shared<const SplitConstruction>(std::move(construction)) ;
      for(auto previous = first[split.parent]; previous != noTreeNode; previous = next[previous]) {
        if(*splits[previous].construction == *split.construction)
          throw std::invalid_argument("duplicate accepted refinement construction") ;
      }
      next[index] = first[split.parent] ; first[split.parent] = index ;
    }
    if(std::find(owners.begin()+1,owners.end(),noTreeNode) != owners.end())
      throw std::invalid_argument("refinement record contains orphan cells") ;
    addBytes(retained,result.geometry_.retainedBytes(),limits.exact) ;
    result.nodes_ = std::move(nodes) ; result.splits_ = std::move(splits) ; result.owners_ = std::move(owners) ;
    result.firstSplit_ = std::move(first) ; result.nextSplit_ = std::move(next) ; result.retainedBytes_ = retained ;
    for(TreeNodeId id = 1; id < result.nodes_.size(); ++id) {
      const auto &cell = result.nodes_[id] ;
      const auto &construction = *result.splits_[result.owners_[id]].construction ;
      if(construction.planar) {
        const auto &planar=*construction.planar ;
        const auto *record=std::get_if<CornerRecord>(&result.geometry_.record(cell.map)) ;
        std::vector<Point> corners ;
        for(unsigned vertex : planar.children[cell.ordinal]) corners.push_back(planar.points[vertex]) ;
        if(!record || record->shape != Shape::quad || record->corners != corners)
          throw std::invalid_argument("refinement geometry differs from accepted planar polygon") ;
        if(!result.active(id) && !cell.isLeaf())
          throw std::invalid_argument("inactive refinement branch has active descendants") ;
        continue ;
      }
      const auto *composed = std::get_if<ComposeRecord>(&result.geometry_.record(cell.map)) ;
      if(!composed || composed->parent != result.nodes_[cell.parent].map)
        throw std::invalid_argument("refinement child does not retain parent geometry") ;
      if(construction.vertexVolume) {
        const auto &mapped = *construction.vertexVolume ;
        const auto *placement = std::get_if<CornerRecord>(&result.geometry_.record(composed->placement)) ;
        std::vector<Point> corners ;
        for(unsigned vertex : mapped.children[cell.ordinal].vertices) corners.push_back(mapped.points[vertex]) ;
        if(!placement || placement->shape != Shape::hex || placement->corners != corners)
          throw std::invalid_argument("refinement geometry differs from accepted mapped chart") ;
      } else {
        const auto &local = result.localCell(id) ;
        const auto *placement = std::get_if<AffineRecord>(&result.geometry_.record(composed->placement)) ;
        if(!placement || placement->shape != local.shape || !samePlacement(placement->placement(),local.placement))
          throw std::invalid_argument("refinement geometry differs from accepted placement") ;
      }
      if(!result.active(id) && !cell.isLeaf())
        throw std::invalid_argument("inactive refinement branch has active descendants") ;
    }
    return result ;
  }

  HistoryEdit::MappedRefinement HistoryEdit::beforeMappedRefinement(TreeNodeId leaf) const {
    const auto &node = history_.node(leaf) ;
    if(!node.isLeaf() || !history_.active(leaf))
      throw std::invalid_argument("mapped refinement requires an active mapped leaf") ;
    return {leaf,history_.nodes_.size(),history_.splits_.size(),history_.geometry_.size(),
      history_.firstSplit_[leaf],history_.retainedBytes_} ;
  }

  void HistoryEdit::undoMappedRefinement(const MappedRefinement &before) noexcept {
    history_.nodes_[before.parent].activeSplit = noTreeNode ;
    history_.firstSplit_[before.parent] = before.firstSplit ;
    history_.nodes_.resize(before.nodes) ; history_.owners_.resize(before.nodes) ; history_.firstSplit_.resize(before.nodes) ;
    history_.splits_.resize(before.splits) ; history_.nextSplit_.resize(before.splits) ;
    history_.geometry_.truncate(before.maps) ; history_.retainedBytes_ = before.retainedBytes ;
  }

  void HistoryEdit::checkFamily(TreeNodeId parent) const {
    if(!history_.active(parent)) throw std::invalid_argument("refinement requires an active cell") ;
    for(TreeNodeId child : history_.children(parent))
      if(!history_.node(child).isLeaf())
        throw std::invalid_argument("refinement requires a complete active leaf family") ;
  }

  CellChange HistoryEdit::change(TreeNodeId parent,
                                const std::shared_ptr<const SplitConstruction> &construction) const {
    CellChange result ; result.parent = history_.key(parent) ;
    if(history_.node(parent).isLeaf()) result.oldLeaves.push_back(result.parent) ;
    else for(TreeNodeId child : history_.children(parent)) result.oldLeaves.push_back(history_.key(child)) ;
    if(!construction) result.newLeaves.push_back(result.parent) ;
    else for(std::size_t ordinal = 0; ordinal < construction->childCount(); ++ordinal) {
      CellKey key = result.parent ; key.lineage.push_back({construction,static_cast<unsigned>(ordinal)}) ;
      result.newLeaves.push_back(std::move(key)) ;
    }
    return result ;
  }

  CellChange HistoryEdit::select(TreeNodeId parent, std::size_t split) {
    CellChange result = change(parent,split == noTreeNode ? nullptr : history_.splits_[split].construction) ;
    history_.nodes_[parent].activeSplit = split ; return result ;
  }

  std::size_t HistoryEdit::preflight(TreeNodeId parent, const SplitOperation &request) const {
    checkFamily(parent) ;
    ExactBudget budget(history_.limits_.exact) ;
    const auto count = requestCount(history_.geometry_.map(history_.node(parent).map)->shape(),request,budget) ;
    if(count == 1) return count ;
    if(history_.path(parent).size() >= history_.limits_.maxDepth)
      throw std::length_error("refinement exceeds its depth limit") ;
    budget.operations(product(count,256)) ;
    if(count > std::numeric_limits<unsigned>::max() || count >= history_.limits_.maxNodes)
      throw std::length_error("refinement child count exceeds its range") ;
    addBytes(0,product(count,sizeof(Child)),history_.limits_.exact) ;
    return count ;
  }

  CellChange HistoryEdit::refine(TreeNodeId leaf, const SplitOperation &request) {
    if(!history_.node(leaf).isLeaf()) throw std::invalid_argument("refinement requires an active leaf") ;
    if(preflight(leaf,request) == 1) return change(leaf,nullptr) ;
    return append(leaf,request,construct(history_.geometry_.map(history_.node(leaf).map)->shape(),request)) ;
  }

  CellChange HistoryEdit::coarsen(TreeNodeId parent) {
    checkFamily(parent) ;
    if(history_.node(parent).isLeaf()) throw std::invalid_argument("coarsening requires an active leaf family") ;
    return select(parent,noTreeNode) ;
  }

  CellChange HistoryEdit::replace(TreeNodeId parent, const SplitOperation &request) {
    if(preflight(parent,request) == 1) throw std::invalid_argument("replacement requires a nonempty split") ;
    return append(parent,request,construct(history_.geometry_.map(history_.node(parent).map)->shape(),request)) ;
  }

  CellChange HistoryEdit::reactivate(TreeNodeId parent, const SplitConstruction &construction) {
    checkFamily(parent) ;
    const auto split = history_.splitId(parent,construction) ;
    if(split == noTreeNode) throw std::invalid_argument("refinement construction is not retained for this parent") ;
    return select(parent,split) ;
  }

  CellChange HistoryEdit::append(TreeNodeId parentId, SplitOperation request, std::vector<Child> children) {
    static_assert(std::is_nothrow_move_constructible<CellNode>::value &&
                  std::is_nothrow_move_constructible<AcceptedSplit>::value &&
                  std::is_nothrow_move_constructible<CellChange>::value,
                  "refinement publication requires nonthrowing moves") ;
    checkFamily(parentId) ;
    const auto &parent = history_.node(parentId) ;
    if(children.size() == 1)
      return change(parentId,parent.isLeaf() ? nullptr : history_.splits_[parent.activeSplit].construction) ;
    const Shape shape = history_.geometry_.map(history_.node(parentId).map)->shape() ;
    ExactBudget budget(history_.limits_.exact) ;
    if(requestCount(shape,request,budget) != children.size())
      throw std::invalid_argument("refinement child count differs from request") ;
    normalizeRequest(request) ; detail::validateSubdivision(shape,children,request,budget) ;
    const auto tag = family(request) ;
    return publish(parentId,std::move(request),{shape,tag,1,std::move(children),nullptr,nullptr}) ;
  }

  CellChange HistoryEdit::refine(TreeNodeId leaf, const VertexVolumeSplit &split, JacobianCheckLimits limits) {
    if(!history_.node(leaf).isLeaf()) throw std::invalid_argument("refinement requires an active leaf") ;
    return replace(leaf,split,limits) ;
  }

  CellChange HistoryEdit::replace(TreeNodeId parent, const VertexVolumeSplit &split, JacobianCheckLimits limits) {
    checkFamily(parent) ;
    if(split.shape != history_.geometry_.map(history_.node(parent).map)->shape())
      throw std::invalid_argument("mapped construction has the wrong parent shape") ;
    limits.maxDepth=std::min(limits.maxDepth,history_.limits_.jacobian.maxDepth) ;
    limits.maxBoxes=std::min(limits.maxBoxes,history_.limits_.jacobian.maxBoxes) ;
    return publish(parent,std::nullopt,{split.shape,SplitFamily::vertexVolume,1,{},validatedMapped(split,limits),nullptr}) ;
  }

  CellChange HistoryEdit::refine(TreeNodeId leaf, const PlanarSplit &split) {
    if(!history_.node(leaf).isLeaf()) throw std::invalid_argument("refinement requires an active leaf") ;
    return replace(leaf,split) ;
  }

  CellChange HistoryEdit::replace(TreeNodeId parent, const PlanarSplit &split) {
    checkFamily(parent) ;
    auto accepted=validatedPlanar(split,*history_.geometry_.map(history_.node(parent).map)) ;
    return publish(parent,std::nullopt,{split.shape,SplitFamily::planar,1,{},nullptr,std::move(accepted)}) ;
  }

  CellChange HistoryEdit::publish(TreeNodeId parentId, std::optional<SplitOperation> request,
                                 SplitConstruction proposed) {
    const auto saved = history_.splitId(parentId,proposed) ;
    if(saved != noTreeNode) return select(parentId,saved) ;
    const auto count = proposed.childCount() ;
    auto &nodes = history_.nodes_ ; auto &splits = history_.splits_ ;
    if(count > history_.limits_.maxNodes-nodes.size())
      throw std::length_error("refinement node limit exceeded") ;
    if(splits.size() >= history_.limits_.maxSplits)
      throw std::length_error("refinement split limit exceeded") ;
    if(history_.path(parentId).size() >= history_.limits_.maxDepth || count > std::numeric_limits<unsigned>::max())
      throw std::length_error("refinement exceeds its depth or ordinal limit") ;
    std::size_t retained = addBytes(history_.retainedBytes_,sizeof(AcceptedSplit)+sizeof(SplitConstruction)+
      sizeof(std::size_t)+(request ? operationBytes(*request) : 0),history_.limits_.exact) ;
    retained = addBytes(retained,product(count,sizeof(CellNode)+2*sizeof(std::size_t)),history_.limits_.exact) ;
    for(const Child &child : proposed.children)
      retained = addBytes(retained,sizeof(Child)+placementBytes(child.placement),history_.limits_.exact) ;
    if(proposed.vertexVolume) retained = addBytes(retained,mappedBytes(*proposed.vertexVolume),history_.limits_.exact) ;
    if(proposed.planar) retained = addBytes(retained,planarBytes(*proposed.planar),history_.limits_.exact) ;
    addBytes(retained,history_.geometry_.retainedBytes(),history_.limits_.exact) ;
    const auto construction = std::make_shared<const SplitConstruction>(std::move(proposed)) ;
    CellChange result = change(parentId,construction) ;
    const NodeRange range(nodes.size(),count) ;
    AcceptedSplit accepted = {parentId,range,std::move(request),construction} ;
    std::vector<GeometryRecord> maps ; std::vector<CellNode> additions ;
    maps.reserve(product(count,2)) ; additions.reserve(count) ;
    const MapId firstMap = history_.geometry_.size() ;
    for(std::size_t ordinal = 0; ordinal < count; ++ordinal) {
      const MapId parentMap=nodes[parentId].map ;
      const MapId placement = firstMap+maps.size() ;
      if(construction->planar) {
        const auto &planar=*construction->planar ;
        std::vector<Point> corners ;
        for(unsigned vertex : planar.children[ordinal]) corners.push_back(planar.points[vertex]) ;
        maps.push_back(CornerRecord{Shape::quad,std::move(corners)}) ;
        additions.push_back({parentId,static_cast<unsigned>(ordinal),placement,noTreeNode}) ;
        continue ;
      }
      if(construction->vertexVolume) {
        const auto &mapped = *construction->vertexVolume ;
        std::vector<Point> corners ;
        for(unsigned vertex : mapped.children[ordinal].vertices) corners.push_back(mapped.points[vertex]) ;
        maps.push_back(CornerRecord{Shape::hex,std::move(corners)}) ;
      } else {
        const Child &child = construction->children[ordinal] ;
        maps.push_back(AffineRecord{child.shape,child.placement}) ;
      }
      maps.push_back(ComposeRecord{parentMap,placement}) ;
      additions.push_back({parentId,static_cast<unsigned>(ordinal),placement+1,noTreeNode}) ;
    }
    // All allocating work precedes publication of geometry and the selected split.
    reserveGrowth(nodes,nodes.size()+additions.size(),history_.limits_.maxNodes) ;
    reserveGrowth(splits,splits.size()+1,history_.limits_.maxSplits) ;
    reserveGrowth(history_.owners_,nodes.size()+additions.size(),history_.limits_.maxNodes) ;
    reserveGrowth(history_.firstSplit_,nodes.size()+additions.size(),history_.limits_.maxNodes) ;
    reserveGrowth(history_.nextSplit_,splits.size()+1,history_.limits_.maxSplits) ;
    history_.geometry_.append(std::move(maps),history_.limits_.exact.maxBytes-retained) ;
    const auto selected = splits.size() ; splits.push_back(std::move(accepted)) ;
    history_.nextSplit_.push_back(history_.firstSplit_[parentId]) ; history_.firstSplit_[parentId] = selected ;
    for(auto &cell : additions) {
      nodes.push_back(std::move(cell)) ; history_.owners_.push_back(selected) ; history_.firstSplit_.push_back(noTreeNode) ;
    }
    history_.retainedBytes_ = retained ; nodes[parentId].activeSplit = selected ;
    return result ;
  }

}
}
