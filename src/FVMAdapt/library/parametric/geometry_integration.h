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

#ifndef LOCI_FVMADAPT_GEOMETRY_INTEGRATION_H
#define LOCI_FVMADAPT_GEOMETRY_INTEGRATION_H

#include <FVMAdapt/parametric/geometry_measure.h>
#include <functional>

namespace Loci { namespace FVMAdapt { namespace detail {

  /// Physical position and positive measure per unit reference volume, area or length.
  struct MeasureSample {
    Point position ;
    long double density ;
  } ;

  /// One Gaussian rule, with first moments measured from the supplied origin.
  struct MeasureIntegral {
    long double measure = 0 ;
    std::array<long double,3> moment{{0,0,0}} ;
    long double radius = 0 ;
  } ;

  MeasureIntegral integrateGeometryRule(Shape shape,
    const std::function<MeasureSample(const Point&)> &sample, unsigned order, const Point &origin) ;

  /// Shared Gaussian driver for maps and factored reference transformations.
  GeometryMeasure integrateGeometry(Shape shape,
    const std::function<MeasureSample(const Point&)> &sample, MeasureOptions options) ;

} } }
#endif
