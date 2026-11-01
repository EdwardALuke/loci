// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_PARAMETRIC_MIXED_H
#define LOCI_FVMADAPT_PARAMETRIC_MIXED_H
#include <FVMAdapt/parametric_cell.h>

namespace Loci { namespace FVMAdapt {
  /// Select directional standard-cell or isotropic general-cell splitting.
  std::string parametricCellFamily(const CellTopology &,const std::string &) ;
  /// Original-face coordinates of corner points; requires midpoint face constructions.
  std::map<ParametricNodeKey,ReferencePoint> cornerFaceCoordinates(EntityId,const ParametricCell &) ;
  void validateMixedFacePlan(const ParametricRecord &) ;
  bool mixedFace(const ParametricCell &,const ParametricCell *) ;
  std::vector<ParametricOutputFace> mixedFaceFaces(EntityId,const CellTopology &,const ParametricCell &,
      const CellTopology &,const ParametricCell &,const std::vector<ParametricNodeKey> &,bool patches = false) ;
  std::vector<ParametricNodeKey> mixedFaceNodes(EntityId,const CellTopology &,const ParametricCell &,
      const CellTopology &,const ParametricCell &) ;
  std::vector<Point> mixedFacePositions(EntityId,const CellTopology &,const ParametricCell &,
      const CellTopology &,const ParametricCell &) ;
  ParametricTrace mixedFaceRequests(EntityId,const CellTopology &,const ParametricCell &,
      const CellTopology &,const ParametricCell &,unsigned,NeighborGrading) ;
  std::vector<FaceBalanceRequest> mixedDividedFaces(EntityId,const CellTopology &,const ParametricCell &,
      const CellTopology &,const ParametricCell &) ;
  /// Include agreed mixed-face points in adjoining corner-cell faces and node tags.
  void mixedCornerFaces(const ParametricCell &,const std::vector<ParametricNodeKey> &,std::vector<ParametricOutputFace> &) ;
  void mixedCornerVertices(const ParametricCell &,const std::vector<ParametricNodeKey> &,std::vector<std::vector<ParametricNodeKey>> &) ;
} }
#endif
