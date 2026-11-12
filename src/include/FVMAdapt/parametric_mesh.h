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

#ifndef LOCI_FVMADAPT_PARAMETRIC_MESH_H
#define LOCI_FVMADAPT_PARAMETRIC_MESH_H

#include "cell_boundary.h"
#include "parametric/refinement_record.h"
#include <map>
#include <string>

namespace Loci {
namespace FVMAdapt {

  /// Source-cell admission, preserving nonzero edges and faces when boundary derivatives are singular.
  bool validParametricSource(const CornerMap &geometry, JacobianCheckLimits limits = {}) ;

  /// Serialized accepted history or persistent refinement requests.
  typedef std::vector<unsigned char> ParametricRecord ;

  /// Two interior cubic controls per source edge, ordered by increasing vertex ID; empty rows are unspecified.
  using ParametricEdgeControls = std::vector<std::vector<Point>> ;
  /// Borrowed cubic source edges in reference-element order; the map owns their lifetime.
  struct CubicEdgeView {
    const CubicEdgeUse *data = nullptr ;
    std::size_t size = 0 ;
    const CubicEdgeUse &operator[](std::size_t i) const { return data[i] ; }
  } ;
  /// Empty for maps other than cubic tet, prism, pyramid and hex source geometry.
  CubicEdgeView parametricCubicEdges(const GeometryMap &) ;
  /// Construct a cubic standard volume; return null when every edge is straight.
  std::shared_ptr<const GeometryMap> parametricCubicGeometry(const CellTopology &,const std::vector<Point> &,
      const ParametricEdgeControls &) ;
  /// Accepted canonical controls; exactly linear edges have empty rows.
  ParametricEdgeControls parametricSourceEdgeControls(const CellTopology &,const GeometryMap &) ;
  /// Check original corners, supplied controls and source validity; omitted controls retain saved curves.
  void checkParametricSource(const CellTopology &,const std::vector<Point> &,const GeometryMap &,
      const ParametricEdgeControls &controls = {},JacobianCheckLimits limits = {}) ;
  /// Compare complete accepted corner or cubic source geometry in reference order.
  bool sameParametricSource(const GeometryMap &,const GeometryMap &) ;

  /// Neighbor contacts that require a bounded difference in refinement depth.
  enum class NeighborGrading { faces, edges, vertices } ;
  NeighborGrading parametricGrading(const std::string &name) ;
  /// Module bounds: depth in [0,30], positive retained-node count, 64 MiB per root.
  HistoryLimits parametricHistoryLimits(int maxDepth, int maxNodes) ;
  /// Serialize accepted geometry and retained native split constructions.
  std::vector<unsigned char> encodeParametricHistory(const AffineSubdivision &history) ;
  /// Restore the accepted history without reconstructing its split geometry.
  AffineSubdivision decodeParametricHistory(const std::vector<unsigned char> &record, HistoryLimits limits = {}) ;
  /// Build or load a native cell, verifying its original geometry and identity.
  AffineSubdivision makeParametricHistory(const CellTopology &topology, const std::vector<Point> &corners,
      const ParametricRecord &record, HistoryLimits limits = {},const ParametricEdgeControls &controls = {}) ;
  /// Initialize an empty record, or verify a saved history against the original native geometry.
  std::vector<unsigned char> initialParametricHistory(const CellTopology &topology, const std::vector<Point> &corners,
      const std::vector<unsigned char> &record, HistoryLimits limits = {},const ParametricEdgeControls &controls = {}) ;

  /// Explicit native cut, family-specific isotropic refinement, or complete-family merge.
  enum class ParametricAction { split, isotropic, coarsen } ;

  /// An edit identified independently of storage indices; split is used only for action == split.
  struct ParametricRequest {
    CellKey cell ;
    SplitOperation split ;
    ParametricAction action = ParametricAction::split ;
  } ;
  /// A request resolved against one accepted history snapshot.
  struct ParametricLeafRequest {
    TreeNodeId cell ;
    SplitOperation split ;
    ParametricAction action = ParametricAction::split ;
  } ;
  /// Empty bytes mean no requests. Duplicate or inactive cell keys reject; split validity is checked when applied.
  std::vector<unsigned char> encodeParametricRequests(const RefinementHistory &history, const std::vector<ParametricRequest> &requests) ;
  /// Convenience overload for the native isotropic pattern of each named cell.
  std::vector<unsigned char> encodeParametricRequests(const RefinementHistory &history, const std::vector<CellKey> &cells) ;
  std::vector<ParametricLeafRequest> decodeParametricRequests(const std::vector<unsigned char> &record, const RefinementHistory &history) ;

  /// Native midpoint grading split. axis == -1 refines in 3D; 0/1/2 preserves a straight X/Y/Z extrusion.
  /// Extrusion mode requires hexes or prisms with matching end faces and no transverse shear.
  SplitOperation parametricGradingSplit(const GeometryMap &root, int axis) ;
  /// Also require every active child to span the original extrusion thickness.
  SplitOperation parametricGradingSplit(const AffineSubdivision &history, int axis) ;

  /// Compact module records exchanged on original faces and edges.
  typedef std::vector<std::vector<int>> ParametricTrace ;
  /// A leaf refinement request (face == -1), or one actively divided cell-local face.
  /// Leaf indices follow leafIds(); root is the original cell file number.
  struct FaceBalanceRequest {
    int root, leaf, face ;
    unsigned directions = 0 ; ///< Optional leaf-local X/Y/Z bits; prism X/Y select the whole base.
  } ;
  /// Active leaf faces covered by multiple neighboring pieces, locally or across an original face.
  /// When requested, also identify the directions divided on each face.
  std::vector<FaceBalanceRequest> parametricDividedFaces(const CellTopology &, const AffineSubdivision &, bool directions = false) ;
  std::vector<FaceBalanceRequest> parametricDividedFaces(EntityId, const CellTopology &, const AffineSubdivision &,
      const CellTopology &, const AffineSubdivision &, bool directions = false) ;

  /// An emitted polygon, oriented from left to right; cell indices use history().leafIds().
  /// right == -1 denotes an original boundary; the caller supplies its boundary tag.
  struct ParametricFace {
    std::vector<VertexKey> vertices ;
    int left, right ;
  } ;

  /// Recover a standard element or a closed polygon shell from outward face loops.
  CellTopology parametricCellTopology(EntityId cell, const std::vector<EntityId> &faces,
      const std::vector<std::vector<EntityId>> &outwardFaces,
      const std::vector<EntityId> &edges, const std::vector<std::array<EntityId,2>> &edgeEnds) ;
  /// Original IDs and either a standard shape or explicit face/edge connectivity. IDs fit nonnegative int.
  std::vector<int> encodeParametricTopology(const CellTopology &topology) ;
  CellTopology decodeParametricTopology(const std::vector<int> &record) ;

  /// Indices in history().leafIds() requiring a split within one original cell.
  std::vector<int> parametricLocalRequests(const CellTopology &topology, const AffineSubdivision &history,
      unsigned difference, NeighborGrading grading) ;
  /// [original cell ID, index in history().leafIds()] for unequal levels across a shared face.
  ParametricTrace parametricFaceRequests(EntityId face, const CellTopology &left, const AffineSubdivision &leftHistory,
      const CellTopology &right, const AffineSubdivision &rightHistory, unsigned difference, NeighborGrading grading) ;
  /// Rows [edge ID, cell ID, leaf index, depth, lower numerator/denominator, upper numerator/denominator].
  ParametricTrace parametricEdgeTrace(const CellTopology &topology, const AffineSubdivision &history) ;
  std::vector<int> parametricEdgeRequests(EntityId root, const ParametricTrace &trace,
      unsigned difference, NeighborGrading grading) ;
  /// Maximum touching leaf depth at each original vertex, in native vertex order.
  std::vector<int> parametricVertexDepths(const CellTopology &topology, const AffineSubdivision &history) ;
  std::vector<int> parametricVertexRequests(const CellTopology &topology, const AffineSubdivision &history,
      const std::map<EntityId,int> &depths, unsigned difference) ;

  /// Sorted unique nodes owned by the original cell or face, including matching-cut intersections.
  std::vector<VertexKey> parametricCellNodes(const CellTopology &topology, const AffineSubdivision &history) ;
  std::vector<VertexKey> parametricFaceNodes(EntityId face, const CellTopology &left,
      const AffineSubdivision &leftHistory, const CellTopology *right = nullptr,
      const AffineSubdivision *rightHistory = nullptr) ;
  /// Interior polygons; supply agreed original face and edge nodes.
  std::vector<ParametricFace> parametricCellFaces(const CellTopology &topology,
      const AffineSubdivision &history, const std::vector<VertexKey> &boundaryNodes) ;
  /// Common original-face polygons; supply agreed original edge nodes.
  std::vector<ParametricFace> parametricFaceFaces(EntityId face, const CellTopology &left,
      const AffineSubdivision &leftHistory, const std::vector<VertexKey> &edgeNodes,
      const CellTopology *right = nullptr, const AffineSubdivision *rightHistory = nullptr) ;

}
}

#endif
