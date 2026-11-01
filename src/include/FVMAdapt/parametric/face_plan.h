// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_FACE_PLAN_H
#define LOCI_FVMADAPT_FACE_PLAN_H

#include "boundary_record.h"

namespace Loci {
namespace FVMAdapt {
  struct CornerBoundary ;

  /// Accepted polygon divisions in canonical face order. Empty means no divisions.
  /// Contains centers and midpoint edge trees, without physical coordinates or geometry maps.
  using FacePlan = std::vector<unsigned char> ;

  FacePlan encodeFacePlan(const BoundaryRecords &, BoundaryId face, const BoundaryRecordLimits &limits = {}) ;
  /// Union accepted divisions; conflicting centers reject. Both inputs use the same root face order.
  FacePlan mergeFacePlans(const FacePlan &, const FacePlan &, const BoundaryRecordLimits &limits = {}) ;
  /// Apply accepted choices to existing source geometry. Failure leaves the boundary unchanged.
  bool applyFacePlan(CornerBoundary &, BoundaryId face, const FacePlan &) ;
}
}
#endif
