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

#ifndef LOCI_FVMADAPT_GEOMETRY_MEASURE_H
#define LOCI_FVMADAPT_GEOMETRY_MEASURE_H

#include "geometry_map.h"

#include <cstddef>
#include <limits>

namespace Loci {
namespace FVMAdapt {

  /// Bounds for successive Gaussian rules of orders 2, 4, ... up to maxOrder.
  struct MeasureOptions {
    /// Relative measure accuracy and centroid accuracy relative to geometry size, plus coordinate roundoff.
    double relativeTolerance = 1e-10 ;
    /// Inclusive bound in [2,256]; only even orders are evaluated.
    unsigned maxOrder = 16 ;
    /// Total map evaluations allowed for one measurement request.
    std::size_t maxEvaluations = 100000 ;
  } ;

  /// Physical length, area or volume and its first-moment centroid.
  struct GeometryMeasure {
    double measure = 0 ;
    Point centroid{{0,0,0}} ;
    /// Successive rules agree within tolerance and coordinate roundoff; this is not an error bound or validity proof.
    bool converged = false ;
    /// Highest completed rule order and cumulative map evaluations.
    unsigned order = 0 ;
    std::size_t evaluations = 0 ;
    /// Successive-rule or summed absolute estimates, not certified error bounds.
    double estimatedMeasureError = std::numeric_limits<double>::infinity() ;
    double estimatedCentroidError = std::numeric_limits<double>::infinity() ;
  } ;

  /// Integrate a retained map using its Jacobian; point maps are unsupported.
  /// A curved face's physical centroid need not lie on the face.
  /// Exhausted bounds return the last complete estimate with converged=false; no first rule throws length_error.
  /// Invalid options or sampled geometry throw invalid_argument; an unrepresentable result throws overflow_error.
  /// Sampled singularities and negative volume Jacobians throw; unsampled defects may remain.
  GeometryMeasure geometryMeasure(const GeometryMap &geometry, MeasureOptions options = {}) ;

}
}

#endif
