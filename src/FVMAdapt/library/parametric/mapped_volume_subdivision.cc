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

#include <FVMAdapt/parametric/mapped_volume_subdivision.h>
#include <FVMAdapt/parametric/refinement_record.h>

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace Loci {
namespace FVMAdapt {

  namespace {
    /// Retain the original volume map as an unsplit cell history.
    RefinementHistory initialHistory(std::shared_ptr<const GeometryMap> map, std::uint64_t root,
                                     std::uint64_t revision, HistoryLimits limits) {
      if(!map || (map->shape() != Shape::hex && map->shape() != Shape::prism && map->shape() != Shape::tet))
        throw std::invalid_argument("mapped-volume subdivision requires a hex, prism or tet map") ;
      GeometryStore geometry(limits.geometry) ; const auto rootMap = geometry.capture(std::move(map)) ;
      return RefinementHistory(root,std::move(geometry),rootMap,revision,limits) ;
    }

    HistoryLimits nodeLimit(std::size_t maxNodes) {
      HistoryLimits limits ; limits.maxNodes=maxNodes ; return limits ;
    }
  }

  MappedVolumeSubdivision::MappedVolumeSubdivision(std::shared_ptr<const GeometryMap> geometry, std::size_t maxNodes)
    : MappedVolumeSubdivision(std::move(geometry),0,0,nodeLimit(maxNodes)) {}

  MappedVolumeSubdivision::MappedVolumeSubdivision(std::shared_ptr<const GeometryMap> geometry, std::uint64_t root,
                                                   std::uint64_t revision, HistoryLimits limits)
    : MappedVolumeSubdivision(initialHistory(std::move(geometry),root,revision,limits)) {}

  /// Require retained vertex-volume constructions before adopting the history.
  MappedVolumeSubdivision::MappedVolumeSubdivision(RefinementHistory history) : edit_(std::move(history)) {
    const auto &saved = this->history() ;
    const Shape shape = saved.geometry().map(saved.node(0).map)->shape() ;
    if(shape != Shape::hex && shape != Shape::prism && shape != Shape::tet)
      throw std::invalid_argument("mapped-volume subdivision requires a hex, prism or tet history") ;
    for(const auto &split : saved.splits())
      if(!split.construction->vertexVolume)
        throw std::invalid_argument("mapped-volume subdivision requires vertex-volume constructions") ;
  }

  /// Follow child ordinals through the selected constructions.
  TreeNodeId MappedVolumeSubdivision::nodeId(const CellPath &path) const {
    TreeNodeId id = 0 ;
    for(unsigned ordinal : path) {
      const auto selected = children(id) ;
      if(ordinal >= selected.size()) throw std::out_of_range("unknown mapped-volume subdivision path") ;
      id = selected[ordinal] ;
    }
    return id ;
  }

  /// List active leaf IDs in ascending node order.
  std::vector<TreeNodeId> MappedVolumeSubdivision::leafIdsInCreationOrder() const {
    auto result = history().leafIds() ; std::sort(result.begin(),result.end()) ; return result ;
  }

  /// Construct vertex-centered children using the supplied face and cell centers.
  std::vector<TreeNodeId> MappedVolumeSubdivision::splitLeaf(TreeNodeId leaf, const std::vector<Point> &faceCenters,
                                                 const Point &cellCenter, JacobianCheckLimits limits) {
    const auto &parent = node(leaf) ;
    if(!parent.isLeaf() || !history().active(leaf)) throw std::invalid_argument("mapped-volume subdivision requires an active leaf") ;
    return append(leaf,std::make_shared<const VertexVolumeSplit>(
      splitVertexVolume(geometry(leaf)->shape(),faceCenters,cellCenter,limits))) ;
  }

  /// Save history extents so an incomplete mapped append can be undone.
  MappedVolumeSubdivision::AppendGuard::AppendGuard(MappedVolumeSubdivision &s, TreeNodeId leaf)
    : state(s), before(s.edit_.beforeMappedRefinement(leaf)) {}

  /// Undo newly appended history unless publication completed.
  MappedVolumeSubdivision::AppendGuard::~AppendGuard() {
    if(committed) return ;
    state.edit_.undoMappedRefinement(before) ;
  }

  /// Publish a validated volume split and return its child IDs.
  std::vector<TreeNodeId> MappedVolumeSubdivision::append(TreeNodeId leaf, std::shared_ptr<const VertexVolumeSplit> split) {
    if(!node(leaf).isLeaf() || !history().active(leaf))
      throw std::invalid_argument("mapped-volume subdivision requires an active leaf") ;
    if(!split || split->shape != geometry(leaf)->shape())
      throw std::invalid_argument("mapped-volume construction has the wrong shape") ;
    std::vector<TreeNodeId> result(split->children.size()) ;
    edit_.refine(leaf,*split) ;
    const auto accepted = children(leaf) ;
    std::copy(accepted.begin(),accepted.end(),result.begin()) ;
    return result ;
  }

  /// Select the parent of a complete active leaf family.
  CellChange MappedVolumeSubdivision::coarsen(TreeNodeId parent) {
    return edit_.coarsen(parent) ;
  }

  /// Replace a leaf family using the supplied face and cell centers.
  CellChange MappedVolumeSubdivision::replace(TreeNodeId parent, const std::vector<Point> &faceCenters,
                                              const Point &cellCenter, JacobianCheckLimits limits) {
    HistoryEdit candidate(history()) ;
    auto split=splitVertexVolume(geometry(parent)->shape(),faceCenters,cellCenter,limits) ;
    auto change=candidate.replace(parent,split,limits) ;
    MappedVolumeSubdivision accepted(std::move(candidate).finish()) ; swap(accepted) ; return change ;
  }

  /// Select an already retained volume construction.
  CellChange MappedVolumeSubdivision::reactivate(TreeNodeId parent, const SplitConstruction &construction) {
    return edit_.reactivate(parent,construction) ;
  }

  /// Discard inactive history while preserving the selected cell geometry.
  void MappedVolumeSubdivision::compact() {
    MappedVolumeSubdivision result(compactHistory(history())) ; swap(result) ;
  }

  /// Validate a replacement history before installing it.
  void MappedVolumeSubdivision::restoreHistory(RefinementHistory history) {
    if(history.rootKey() != this->history().rootKey() || history.geometryRevision() != this->history().geometryRevision())
      throw std::invalid_argument("mapped-volume history belongs to a different root or geometry revision") ;
    auto accepted=RefinementHistory::restore(history.rootKey(),history.geometryRevision(),history.geometry(),
      history.nodes(),history.splits(),this->history().limits()) ;
    const auto rootGeometry=[](const RefinementHistory &source) {
      auto root=RefinementHistory(source.rootKey(),source.geometry(),source.node(0).map,source.geometryRevision(),source.limits()) ;
      return encodeGeometry(compactHistory(root).geometry()) ;
    } ;
    if(rootGeometry(accepted) != rootGeometry(this->history()))
      throw std::invalid_argument("mapped-volume history has different root geometry") ;
    MappedVolumeSubdivision result(std::move(accepted)) ; swap(result) ;
  }

  void MappedVolumeSubdivision::swap(MappedVolumeSubdivision &other) noexcept {
    edit_.swap(other.edit_) ;
  }

}
}
