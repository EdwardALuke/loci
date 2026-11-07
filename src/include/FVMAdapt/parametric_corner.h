// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_PARAMETRIC_CORNER_H
#define LOCI_FVMADAPT_PARAMETRIC_CORNER_H

#include "parametric_mesh.h"
#include "parametric/corner_boundary.h"
#include "parametric/corner_subdivision.h"

namespace Loci {
namespace FVMAdapt {
  /// Shared source-face geometry, selected from its canonical vertex order.
  std::shared_ptr<const GeometryMap> sourceFaceGeometry(Shape,const std::vector<Point> &) ;
  /// Extra nearly straight face corners require a piecewise surface representation.
  bool needsFacePatches(const CellTopology &,const std::vector<Point> &) ;

  /// Vertex-based splits of a standard mapped cell and its accepted boundary constructions.
  /// Shared source faces are reconciled by the module before acceptance.
  struct ParametricCorner : CornerBoundary {
    CellTopology topology ;
    MappedVolumeSubdivision subdivision ;
    std::vector<CornerCell> cells ;
    ParametricCorner(const CellTopology &, const std::vector<Point> &, HistoryLimits) ;
    ParametricCorner(CellTopology t, MappedVolumeSubdivision s) : topology(std::move(t)), subdivision(std::move(s)) {}
    void refine(TreeNodeId leaf, CenterMethod method) ;
    bool mergeFace(EntityId face, const ParametricRecord &record) ;
    bool faceChanged(EntityId face) const ;
    ParametricRecord faceRecord(EntityId face) const ;
    ParametricRecord record() const ;
    static ParametricCorner restore(const ParametricRecord &, HistoryLimits) ;
    BoundaryId faceId(EntityId face) const ;
  } ;

  /// Cell topology and vertex-based splits defined by shared face maps.
  struct ParametricBoundary : CornerSubdivision {
    CellTopology topology ;
    ParametricBoundary(const CellTopology &, const std::vector<Point> &, HistoryLimits) ;
    ParametricBoundary(CellTopology t, CornerSubdivision s) : CornerSubdivision(std::move(s)), topology(std::move(t)) {}
    std::vector<Point> sourceCorners() const ;
    /// Accepted centers in parent face order, independent of local boundary IDs.
    ParametricRecord splitRecord(TreeNodeId parent) const ;
    ParametricRecord record() const ;
    static ParametricBoundary restore(const ParametricRecord &, HistoryLimits) ;
  } ;

  class ParametricCell ;
  struct ParametricNodeKey ;
  struct ParametricOutputFace ;
  /// Active local or shared face polygons, with agreed hanging edge points.
  std::vector<ParametricOutputFace> cornerCellFaces(const ParametricCell &, const std::vector<ParametricNodeKey> &, bool patches = false) ;
  std::vector<ParametricOutputFace> cornerFaceFaces(EntityId, const ParametricCell &, const std::vector<ParametricNodeKey> &,
                                                  const ParametricCell *right = nullptr, bool patches = false) ;
  /// Own vertex keys in leaf order.
  std::vector<std::vector<ParametricNodeKey>> cornerLeafVertices(const ParametricCell &) ;
  std::vector<ParametricNodeKey> cornerNodes(const ParametricCell &, EntityKind, EntityId,
                                            const ParametricCell *right = nullptr) ;
  std::vector<Point> cornerPositions(const ParametricCell &, const std::vector<ParametricNodeKey> &) ;
  ParametricTrace cornerEdgeTrace(const ParametricCell &) ;
  std::vector<int> cornerVertexDepths(const ParametricCell &) ;
  std::vector<int> cornerVertexRequests(const ParametricCell &, const std::map<EntityId,int> &, unsigned) ;
}
}
#endif
