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

  /// Physical position and positive measure density per unit reference volume, area or length.
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

  /// Origin and physical scales for a whole-domain centroid error check.
  struct MeasureScale { Point origin ; long double radius=0,magnitude=0 ; } ;

  /// Sum positive pieces under one error tolerance and evaluation budget.
  /// A null scale selects the first centroid and the measured domain size.
  /// Initial-estimate budget failure throws; unresolved continuation returns converged=false.
  GeometryMeasure measureGeometryPieces(const std::vector<std::shared_ptr<const GeometryMap>> &pieces,
    MeasureOptions options,const MeasureScale *scale=nullptr) ;

  MeasureIntegral integrateGeometryRule(Shape shape,
    const std::function<MeasureSample(const Point&)> &sample, unsigned order, const Point &origin) ;

  /// Factor a radial triangle's cone into at most two curve integrals; false for other maps.
  bool radialConeSamples(const ConeMap &, std::vector<std::function<MeasureSample(const Point&)>> &samples) ;

  /// Shared Gaussian driver for maps and factored reference transformations.
  GeometryMeasure integrateGeometry(Shape shape,
    const std::function<MeasureSample(const Point&)> &sample, MeasureOptions options) ;

} } }
#endif
