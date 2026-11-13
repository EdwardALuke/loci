// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_PARAMETRIC_MIXED_H
#define LOCI_FVMADAPT_PARAMETRIC_MIXED_H
#include <FVMAdapt/parametric_cell.h>

namespace Loci { namespace FVMAdapt {
  /// Working geometry selected for a source cell; compact storage is chosen separately.
  enum class CellRepresentation { native, mappedCorner, boundary, extruded } ;
  struct CellSelection {
    CellRepresentation representation ;
    bool facePatches ; ///< Source faces need piecewise maps and cannot use deferred corner geometry.
    bool corner() const {
      return representation==CellRepresentation::mappedCorner || representation==CellRepresentation::boundary ;
    }
  } ;
  /// Resolve split family, source-face representation and extrusion compatibility.
  CellSelection selectParametricCell(const CellTopology &topology,const std::vector<Point> &corners,const std::string &family,int axis) ;

  /// Construct using the selection already resolved for these source coordinates.
  /// Saved-history and physical-geometry checks still apply.
  ParametricCell makeParametricCell(const CellTopology &topology,const std::vector<Point> &corners,
      const ParametricRecord &record,HistoryLimits limits,const CellSelection &selection,int axis,
      const ParametricEdgeControls &controls = {}) ;

  /// Original-face coordinates of corner points; requires midpoint face constructions.
  std::map<ParametricNodeKey,ReferencePoint> cornerFaceCoordinates(EntityId face,const ParametricCell &cell) ;
  void validateMixedFacePlan(const ParametricRecord &record) ;
  bool mixedFace(const ParametricCell &leftCell,const ParametricCell *rightCell) ;
  /// Shared-face queries require matching source-edge geometry, checked by the module before balancing.
  std::vector<ParametricOutputFace> mixedFaceFaces(EntityId face,const CellTopology &leftTopology,const ParametricCell &leftCell,
      const CellTopology &rightTopology,const ParametricCell &rightCell,const std::vector<ParametricNodeKey> &requiredNodes,bool patches = false) ;
  std::vector<ParametricNodeKey> mixedFaceNodes(EntityId face,const CellTopology &leftTopology,const ParametricCell &leftCell,
      const CellTopology &rightTopology,const ParametricCell &rightCell,std::vector<Point> *positions = nullptr) ;
  ParametricTrace mixedFaceRequests(EntityId face,const CellTopology &leftTopology,const ParametricCell &leftCell,
      const CellTopology &rightTopology,const ParametricCell &rightCell,unsigned difference,NeighborGrading grading) ;
  std::vector<FaceBalanceRequest> mixedDividedFaces(EntityId face,const CellTopology &leftTopology,const ParametricCell &leftCell,
      const CellTopology &rightTopology,const ParametricCell &rightCell) ;
  /// Include agreed mixed-face points in adjoining corner-cell faces and node tags.
  void mixedCornerFaces(const ParametricCell &cell,const std::vector<ParametricNodeKey> &requiredNodes,std::vector<ParametricOutputFace> &faces) ;
  void mixedCornerVertices(const ParametricCell &cell,const std::vector<ParametricNodeKey> &requiredNodes,std::vector<std::vector<ParametricNodeKey>> &vertices) ;
} }
#endif
