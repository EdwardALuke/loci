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

#ifndef LOCI_FVMADAPT_PLANAR_SUBDIVISION_H
#define LOCI_FVMADAPT_PLANAR_SUBDIVISION_H

#include "refinement_history.h"
#include "geometry_map.h"
#include "reference_coordinate.h"

#include <map>

namespace Loci {
namespace FVMAdapt {

  /// Parent cell and accepted construction; independent of history table indices.
  struct PlanarSplitKey {
    CellKey parent = {0,{}} ;
    std::shared_ptr<const SplitConstruction> construction ;
  } ;
  bool operator==(const PlanarSplitKey &a, const PlanarSplitKey &b) ;
  bool operator<(const PlanarSplitKey &a, const PlanarSplitKey &b) ;

  /// Original polygon edge, or a spoke created by splitting owner toward edge index.
  /// Original edges follow the root loop; spokes run from the owner center to its edge midpoint.
  struct PlanarEdgeKey {
    bool original = true ;
    PlanarSplitKey owner ;
    unsigned index = 0 ;
  } ;
  bool operator==(const PlanarEdgeKey &a, const PlanarEdgeKey &b) ;
  bool operator<(const PlanarEdgeKey &a, const PlanarEdgeKey &b) ;

  /// Exact construction identity, independent of physical coordinates.
  struct PlanarVertexKey {
    enum class Kind { corner, edge, center } ;
    Kind kind = Kind::corner ;
    unsigned corner = 0 ;
    PlanarEdgeKey edge ;
    Fraction position ;
    PlanarSplitKey owner ;
  } ;
  bool operator==(const PlanarVertexKey &a, const PlanarVertexKey &b) ;
  bool operator<(const PlanarVertexKey &a, const PlanarVertexKey &b) ;

  /// Directed subinterval of one original edge or introduced spoke.
  struct PlanarEdgeUse {
    PlanarEdgeKey edge ;
    Fraction begin, end ;
  } ;

  /// One accepted planar polygon; corner and edge arrays follow its oriented loop.
  struct PlanarNode {
    TreeNodeId parent = noTreeNode ;
    CellPath path ;
    std::vector<PlanarVertexKey> corners ;
    std::vector<PlanarEdgeUse> edges ;
    std::shared_ptr<const GeometryMap> geometry ;
    std::vector<TreeNodeId> children ;
    bool isLeaf() const { return children.empty() ; }
    Shape shape() const { return corners.size() == 3 ? Shape::triangle : Shape::quad ; }
  } ;

  /// Vertex-centered refinement of a strictly convex planar triangle or quad.
  /// Each accepted polygon has its own corner chart. Children preserve its physical
  /// polygon, but their charts are not restrictions of a fixed root bilinear map.
  /// Edits invalidate node references. Compaction may change IDs; history().key(id) persists.
  class PlanarSubdivision {
    HistoryEdit edit_ ;
    std::vector<PlanarNode> nodes_ ;
    std::map<PlanarVertexKey,Point> points_ ;
    PlanarSplit proposedSplit(TreeNodeId parent, CenterMethod method, MeasureOptions options) const ;
    std::vector<PlanarNode> childNodes(TreeNodeId parent, const PlanarSplitKey &owner) const ;
    TreeNodeId ownerNode(const PlanarSplitKey &owner) const ;
    PlanarVertexKey vertexAt(const PlanarEdgeKey &edge, Fraction position, ExactBudget &budget) const ;
    Point evaluate(const PlanarVertexKey &vertex, ExactBudget &budget) const ;
  public:
    /// Root corners follow the polygon loop. Invalid geometry or a zero node limit throws invalid_argument.
    explicit PlanarSubdivision(std::vector<Point> corners, std::size_t maxNodes = 1000000) ;
    PlanarSubdivision(std::vector<Point> corners, std::uint64_t root, std::uint64_t revision, HistoryLimits limits) ;
    /// Restore accepted polygons without choosing new centers.
    explicit PlanarSubdivision(RefinementHistory history) ;
    std::size_t maxNodes() const { return history().limits().maxNodes ; }
    const RefinementHistory &history() const { return edit_.snapshot() ; }
    const std::vector<PlanarNode> &nodes() const { return nodes_ ; }
    const PlanarNode &node(TreeNodeId id) const { return nodes_.at(id) ; }
    /// Resolve the selected path; dormant alternatives require history().nodeId(CellKey).
    TreeNodeId nodeId(const CellPath &path) const ;
    std::vector<TreeNodeId> leafIds() const ;
    /// Original edge endpoints and spoke endpoints alias their existing vertex identities.
    /// Invalid edges or positions outside [0,1] throw invalid_argument or out_of_range.
    PlanarVertexKey vertexAt(const PlanarEdgeKey &edge, Fraction position) const ;
    /// Evaluate a retained corner/center or any exact point on a retained edge support.
    Point evaluate(const PlanarVertexKey &vertex) const ;
    /// One quad per corner; accepted edge midpoints and center define physical child charts.
    /// Invalid requests, limits or invalid children preserve the hierarchy.
    std::vector<TreeNodeId> splitLeaf(TreeNodeId leaf, CenterMethod method = CenterMethod::weighted,
                                     MeasureOptions options = {}) ;
    CellChange coarsen(TreeNodeId parent) ;
    CellChange replace(TreeNodeId parent, CenterMethod method, MeasureOptions options = {}) ;
    CellChange reactivate(TreeNodeId parent, const SplitConstruction &construction) ;
    /// Discard dormant constructions; local indices may change, persistent keys do not.
    void compact() ;
    /// Require matching source geometry and root identity; apply this object's limits.
    void restoreHistory(RefinementHistory history) ;
    void swap(PlanarSubdivision &other) noexcept ;
  } ;

}
}

#endif
