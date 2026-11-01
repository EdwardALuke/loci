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

#ifndef LOCI_FVMADAPT_MAPPED_VOLUME_SUBDIVISION_H
#define LOCI_FVMADAPT_MAPPED_VOLUME_SUBDIVISION_H

#include "refinement_history.h"
#include "vertex_volume_split.h"

namespace Loci {
namespace FVMAdapt {

  /// Vertex-centered hex/prism/tet refinement with immutable composed geometry.
  /// Nodes and selected children belong to history(); edits can invalidate node references.
  class MappedVolumeSubdivision {
    HistoryEdit edit_ ;
    std::vector<TreeNodeId> append(TreeNodeId leaf, std::shared_ptr<const VertexVolumeSplit> split) ;
    friend struct CornerBoundary ;
    /// Keep one leaf append reversible until its shared boundary construction succeeds.
    struct AppendGuard {
      MappedVolumeSubdivision &state ;
      HistoryEdit::MappedRefinement before ;
      bool committed = false ;
      AppendGuard(MappedVolumeSubdivision &state, TreeNodeId leaf) ;
      AppendGuard(const AppendGuard &) = delete ;
      AppendGuard &operator=(const AppendGuard &) = delete ;
      ~AppendGuard() ;
    } ;
  public:
    /// A hex, prism or tet map and a positive node limit are required. Physical map validity is the caller's responsibility.
    explicit MappedVolumeSubdivision(std::shared_ptr<const GeometryMap> geometry, std::size_t maxNodes = 1000000) ;
    MappedVolumeSubdivision(std::shared_ptr<const GeometryMap> geometry, std::uint64_t root,
                            std::uint64_t revision, HistoryLimits limits) ;
    /// Restore accepted geometry and selections without choosing new centers.
    explicit MappedVolumeSubdivision(RefinementHistory history) ;
    std::size_t maxNodes() const { return history().limits().maxNodes ; }
    const RefinementHistory &history() const { return edit_.snapshot() ; }
    const std::vector<CellNode> &nodes() const { return history().nodes() ; }
    const CellNode &node(TreeNodeId id) const { return history().node(id) ; }
    NodeRange children(TreeNodeId id) const { return history().children(id) ; }
    /// Retained physical map, shared with the history and boundary constructions.
    std::shared_ptr<const GeometryMap> geometry(TreeNodeId id) const {
      return history().geometry().map(node(id).map) ;
    }
    /// Resolve the selected path; dormant alternatives require history().nodeId(CellKey).
    TreeNodeId nodeId(const CellPath &path) const ;
    std::vector<TreeNodeId> leafIds() const ;
    CellChange coarsen(TreeNodeId parent) ;
    CellChange replace(TreeNodeId parent, const std::vector<Point> &faceCenters,
                       const Point &cellCenter, JacobianCheckLimits limits = {}) ;
    CellChange reactivate(TreeNodeId parent, const SplitConstruction &construction) ;
    /// Discard dormant constructions; local indices may change, persistent cell keys do not.
    void compact() ;
    void restoreHistory(RefinementHistory history) ;
    void swap(MappedVolumeSubdivision &other) noexcept ;
    /// Split one leaf using boundary-local face centers and a parent-local cell center.
    /// Invalid requests, certification failures and limits preserve the hierarchy.
    /// Child ordinal is the original vertex index; all children are hexes.
    std::vector<TreeNodeId> splitLeaf(TreeNodeId leaf, const std::vector<Point> &faceCenters,
                                const Point &cellCenter, JacobianCheckLimits limits = {}) ;
  } ;

}
}

#endif
