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

#ifndef LOCI_FVMADAPT_GEOMETRY_CENTER_H
#define LOCI_FVMADAPT_GEOMETRY_CENTER_H

#include <FVMAdapt/parametric/geometry_measure.h>

namespace Loci {
namespace FVMAdapt {

  /// Reference center, weighted target, or integrated physical centroid on the retained map.
  enum class CenterMethod { reference, weighted, integrated } ;

  /// Return surface reference coordinates without changing the retained map.
  /// weighted projects the images of reference edge midpoints, weighted by
  /// endpoint chord lengths; it does not use arc-length midpoints or quadrature.
  /// integrated projects geometryMeasure()'s area centroid; options bound that integration.
  /// It finds an interior stationary point, not a certified global nearest point.
  /// Invalid input, singular derivatives or failed convergence throw invalid_argument.
  /// Unresolved integration throws length_error.
  /// This computation does not certify the face or its proposed children.
  Point faceCenterCoordinates(const GeometryMap &geometry,
                              CenterMethod method = CenterMethod::weighted, MeasureOptions options = {}) ;

  /// Return volume reference coordinates for the selected center construction.
  /// faceCenters is required only for weighted: boundary-local coordinates in native boundary order.
  /// weighted uses their physical positions, weighted by magnitudes of corner-polygon area vectors.
  /// integrated inverts geometryMeasure()'s volume centroid; options bound that integration.
  /// The target must have an interior inverse with positive Jacobian; otherwise throws
  /// invalid_argument. Unresolved integration throws length_error. No global injectivity or child validity is certified.
  Point cellCenterCoordinates(const GeometryMap &geometry, const std::vector<Point> &faceCenters,
                              CenterMethod method = CenterMethod::weighted, MeasureOptions options = {}) ;

}
}

#endif
