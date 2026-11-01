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

#ifndef LOCI_FVMADAPT_REFINEMENT_HISTORY_H
#define LOCI_FVMADAPT_REFINEMENT_HISTORY_H

#include "history_types.h"
#include "geometry_record.h"
#include <optional>
#include <set>
#include <map>
#include <string>

namespace Loci {
namespace FVMAdapt {

  /// Bounds include retained ancestors, accepted geometry and exact local coordinates.
  /// exact.maxOperations applies to each construction or root-coordinate query;
  /// maxNodes and exact.maxBytes bound the complete snapshot.
  struct HistoryLimits {
    std::size_t maxNodes = 1000000 ;
    std::size_t maxSplits = 1000000 ;
    unsigned maxDepth = 256 ;
    ExactLimits exact ;
    GeometryLimits geometry ;
    JacobianCheckLimits jacobian ;
  } ;

  /// Cell and its explicitly selected subdivision.
  struct CellNode {
    TreeNodeId parent = noTreeNode ;
    unsigned ordinal = 0 ;
    MapId map = 0 ;
    std::size_t activeSplit = noTreeNode ;
    bool isLeaf() const { return activeSplit == noTreeNode ; }
  } ;

  /// One accepted construction; selection belongs to its parent node.
  struct AcceptedSplit {
    TreeNodeId parent ;
    NodeRange children ;
    /// Original native request; mapped constructions retain their accepted charts instead.
    std::optional<SplitOperation> operation ;
    std::shared_ptr<const SplitConstruction> construction ;
  } ;

  /// Immutable refinement snapshot. Geometry and placements survive request replay.
  class RefinementHistory {
    std::uint64_t rootKey_, geometryRevision_ ;
    HistoryLimits limits_ ;
    Child rootCell_ ;
    GeometryStore geometry_ ;
    std::vector<CellNode> nodes_ ;
    std::vector<AcceptedSplit> splits_ ;
    std::vector<std::size_t> owners_ ;
    std::vector<std::size_t> firstSplit_, nextSplit_ ;
    std::size_t retainedBytes_ = 0 ;
    std::size_t splitId(TreeNodeId parent, const SplitConstruction &construction) const ;
    friend class HistoryEdit ;
  public:
    RefinementHistory(std::uint64_t root, GeometryStore geometry, MapId rootMap,
                      std::uint64_t revision = 0, HistoryLimits limits = HistoryLimits()) ;
    std::uint64_t rootKey() const { return rootKey_ ; }
    std::uint64_t geometryRevision() const { return geometryRevision_ ; }
    const HistoryLimits &limits() const { return limits_ ; }
    const GeometryStore &geometry() const { return geometry_ ; }
    const std::vector<CellNode> &nodes() const { return nodes_ ; }
    const std::vector<AcceptedSplit> &splits() const { return splits_ ; }
    /// Accounted record payload, excluding allocator capacity and evaluator overhead.
    std::size_t retainedBytes() const { return retainedBytes_+geometry_.retainedBytes() ; }
    const CellNode &node(TreeNodeId id) const { return nodes_.at(id) ; }
    /// Accepted affine placement; mapped children throw invalid_argument. The root is the identity.
    const Child &localCell(TreeNodeId id) const ;
    NodeRange children(TreeNodeId id) const ;
    /// Children of a retained construction, selected or dormant; empty if it is not saved.
    NodeRange retainedChildren(TreeNodeId parent, const SplitConstruction &construction) const ;
    const SplitOperation *operation(TreeNodeId id) const ;
    /// Active root-reachable leaves in child order; dormant records are excluded.
    std::vector<TreeNodeId> leafIds() const ;
    bool active(TreeNodeId id) const ;
    CellPath path(TreeNodeId id) const ;
    CellKey key(TreeNodeId id) const ;
    /// Find an active or dormant cell by identity; missing keys throw out_of_range.
    TreeNodeId nodeId(const CellKey &key) const ;
    Placement rootPlacement(TreeNodeId id) const ;
    void swap(RefinementHistory &other) noexcept ;
    /// Import accepted records; validate topology, local geometry and limits before publication.
    /// Child construction is not rerun. Invalid data leaves the caller's snapshot unchanged.
    static RefinementHistory restore(std::uint64_t root, std::uint64_t revision,
                                     GeometryStore geometry, std::vector<CellNode> nodes,
                                     std::vector<AcceptedSplit> splits,
                                     HistoryLimits limits = HistoryLimits()) ;
  } ;

  /// Old and new active leaf families; geometric transfer weights are separate.
  struct CellChange {
    CellKey parent ;
    std::vector<CellKey> oldLeaves, newLeaves ;
    bool changed() const { return oldLeaves != newLeaves ; }
  } ;

  /// Keep active ancestors and their geometry in a new snapshot; all local IDs may change.
  RefinementHistory compactHistory(const RefinementHistory &history) ;

  /// Candidate edit; failures preserve its previous contents and the input snapshot.
  class HistoryEdit {
    RefinementHistory history_ ;
    void checkFamily(TreeNodeId parent) const ;
    CellChange change(TreeNodeId parent, const std::shared_ptr<const SplitConstruction> &construction) const ;
    CellChange select(TreeNodeId parent, std::size_t split) ;
    std::size_t preflight(TreeNodeId leaf, const SplitOperation &request) const ;
    CellChange append(TreeNodeId parent, SplitOperation request, std::vector<Child> children) ;
    CellChange publish(TreeNodeId parent, std::optional<SplitOperation> request, SplitConstruction construction) ;
    friend class AffineSubdivision ;
    /// One mapped leaf edit changes its parent selection and appends to these tables.
    struct MappedRefinement {
      TreeNodeId parent ;
      std::size_t nodes, splits, maps, firstSplit, retainedBytes ;
    } ;
    MappedRefinement beforeMappedRefinement(TreeNodeId leaf) const ;
    void undoMappedRefinement(const MappedRefinement &) noexcept ;
    friend class MappedVolumeSubdivision ;
  public:
    explicit HistoryEdit(RefinementHistory history) : history_(std::move(history)) {}
    /// View of the candidate; copy it to retain a snapshot across later edits.
    const RefinementHistory &snapshot() const { return history_ ; }
    void swap(HistoryEdit &other) noexcept { history_.swap(other.history_) ; }
    /// Refine an active leaf; an empty request leaves it intact.
    CellChange refine(TreeNodeId leaf, const SplitOperation &request) ;
    /// Refine with accepted nonlinear reference charts; physical geometry is composed unchanged.
    CellChange refine(TreeNodeId leaf, const VertexVolumeSplit &split, JacobianCheckLimits limits = {}) ;
    /// Refine with accepted physical planar polygons and independent child corner charts.
    CellChange refine(TreeNodeId leaf, const PlanarSplit &split) ;
    /// Collapse a complete active leaf family; accepted constructions remain available.
    CellChange coarsen(TreeNodeId parent) ;
    /// Select a nonempty split of an active leaf or complete active leaf family.
    CellChange replace(TreeNodeId parent, const SplitOperation &request) ;
    CellChange replace(TreeNodeId parent, const VertexVolumeSplit &split, JacobianCheckLimits limits = {}) ;
    CellChange replace(TreeNodeId parent, const PlanarSplit &split) ;
    /// Select retained children as leaves; deeper selections are not restored.
    CellChange reactivate(TreeNodeId parent, const SplitConstruction &construction) ;
    /// Publish the candidate by move. References into the edit must then be discarded.
    RefinementHistory finish() && { return std::move(history_) ; }
  } ;

}
}

#endif
