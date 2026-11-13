// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_FACE_PLAN_H
#define LOCI_FVMADAPT_FACE_PLAN_H

#include "boundary_record.h"

namespace Loci {
namespace FVMAdapt {
  struct CornerBoundary ;

  /// Accepted polygon divisions in canonical face order. Empty means no divisions.
  /// Each face has at most one accepted construction; conflicting centers reject.
  /// Contains centers and midpoint edge trees, without physical coordinates or geometry maps.
  using FacePlan = std::vector<unsigned char> ;

  FacePlan encodeFacePlan(const BoundaryRecords &boundaries, BoundaryId face, const BoundaryRecordLimits &limits = {}) ;
  /// Union accepted divisions; conflicting centers reject. Both inputs use the same root face order.
  FacePlan mergeFacePlans(const FacePlan &first, const FacePlan &second, const BoundaryRecordLimits &limits = {}) ;
  /// Apply accepted choices to existing source geometry. Failure leaves the boundary unchanged.
  /// Return true when retained divisions are added; false when all are already present.
  bool applyFacePlan(CornerBoundary &boundary, BoundaryId face, const FacePlan &plan) ;
}
}
#endif
