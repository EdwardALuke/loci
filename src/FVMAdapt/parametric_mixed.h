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
  CellSelection selectParametricCell(const CellTopology &,const std::vector<Point> &,const std::string &family,int axis) ;

  /// Construct using the selection already resolved for these source coordinates.
  /// Saved-history and physical-geometry checks still apply.
  ParametricCell makeParametricCell(const CellTopology &topology,const std::vector<Point> &corners,
      const ParametricRecord &record,HistoryLimits limits,const CellSelection &,int axis) ;

  /// Original-face coordinates of corner points; requires midpoint face constructions.
  std::map<ParametricNodeKey,ReferencePoint> cornerFaceCoordinates(EntityId,const ParametricCell &) ;
  void validateMixedFacePlan(const ParametricRecord &) ;
  bool mixedFace(const ParametricCell &,const ParametricCell *) ;
  std::vector<ParametricOutputFace> mixedFaceFaces(EntityId,const CellTopology &,const ParametricCell &,
      const CellTopology &,const ParametricCell &,const std::vector<ParametricNodeKey> &,bool patches = false) ;
  std::vector<ParametricNodeKey> mixedFaceNodes(EntityId,const CellTopology &,const ParametricCell &,
      const CellTopology &,const ParametricCell &,std::vector<Point> *positions = nullptr) ;
  ParametricTrace mixedFaceRequests(EntityId,const CellTopology &,const ParametricCell &,
      const CellTopology &,const ParametricCell &,unsigned,NeighborGrading) ;
  std::vector<FaceBalanceRequest> mixedDividedFaces(EntityId,const CellTopology &,const ParametricCell &,
      const CellTopology &,const ParametricCell &) ;
  /// Include agreed mixed-face points in adjoining corner-cell faces and node tags.
  void mixedCornerFaces(const ParametricCell &,const std::vector<ParametricNodeKey> &,std::vector<ParametricOutputFace> &) ;
  void mixedCornerVertices(const ParametricCell &,const std::vector<ParametricNodeKey> &,std::vector<std::vector<ParametricNodeKey>> &) ;
} }
#endif
