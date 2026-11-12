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

#ifndef LOCI_FVMADAPT_AFFINE_SUBDIVISION_H
#define LOCI_FVMADAPT_AFFINE_SUBDIVISION_H

#include "refinement_history.h"

namespace Loci {
namespace FVMAdapt {

  /// Native cell in root reference coordinates; local accepted geometry lives in the history.
  struct AffineNode {
    TreeNodeId parent ;
    Child cell ;
  private:
    AffineNode(TreeNodeId parent, Child cell)
      : parent(parent), cell(std::move(cell)) {}
    friend class AffineSubdivision ;
  } ;

  /// Root-coordinate view of a native refinement history for boundary assembly.
  /// Node references may be invalidated by splitting; IDs remain valid while this tree grows.
  class AffineSubdivision {
    HistoryEdit edit_ ;
    std::vector<AffineNode> nodes_ ;
    void checkLeaf(TreeNodeId leaf) const ;
    CellChange append(TreeNodeId parent, std::vector<Child> children, SplitOperation operation) ;
  public:
    /// Root geometry is the reference shape; maxNodes includes retained ancestors.
    explicit AffineSubdivision(Shape shape, std::size_t maxNodes = 1000000) ;
    /// Assemble a root-coordinate view without regenerating the accepted children.
    explicit AffineSubdivision(RefinementHistory history) ;
    Shape shape() const { return nodes_[0].cell.shape ; }
    std::size_t maxNodes() const { return history().limits().maxNodes ; }
    const RefinementHistory &history() const { return edit_.snapshot() ; }
    const std::vector<AffineNode> &nodes() const { return nodes_ ; }
    const AffineNode &node(TreeNodeId id) const { return nodes_.at(id) ; }
    /// No selected split; use history().active(id) to exclude dormant cells.
    bool isLeaf(TreeNodeId id) const { return history().node(id).isLeaf() ; }
    NodeRange children(TreeNodeId id) const { return history().children(id) ; }
    const SplitOperation *operation(TreeNodeId id) const { return history().operation(id) ; }
    /// Active leaf IDs in creation order; placements remain in node(id).cell.
    std::vector<TreeNodeId> leafIds() const ;
    CellChange refine(TreeNodeId leaf, const SplitOperation &request) ;
    CellChange coarsen(TreeNodeId parent) ;
    CellChange replace(TreeNodeId parent, const SplitOperation &request) ;
    CellChange reactivate(TreeNodeId parent, const SplitConstruction &construction) ;
    /// Discard dormant records; persistent cell keys remain unchanged.
    void compact() ;
    void swap(AffineSubdivision &other) noexcept ;
    /// Refine an active leaf; an empty request leaves it intact. Invalid requests or
    /// resource exhaustion preserve the tree's contents.
    NodeRange splitLeaf(TreeNodeId leaf, const TensorSplit &request) ;
    NodeRange splitLeaf(TreeNodeId leaf, const TriangleSplit &request) ;
    NodeRange splitLeaf(TreeNodeId leaf, const PrismSplit &request) ;
    NodeRange splitLeaf(TreeNodeId leaf, const TetSplit &request) ;
    NodeRange splitLeaf(TreeNodeId leaf, const PyramidSplit &request) ;
  } ;

}
}

#endif
