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

#ifndef LOCI_FVMADAPT_PLANAR_SPLIT_H
#define LOCI_FVMADAPT_PLANAR_SPLIT_H

#include "geometry_center.h"

namespace Loci {
namespace FVMAdapt {

  /// Accepted physical polygon fan: corners, edge midpoints, then the center.
  /// Each child is a quad with its own physical corner chart.
  struct PlanarSplit {
    Shape shape ;
    CenterMethod centerMethod = CenterMethod::weighted ;
    std::vector<Point> points ;
    std::vector<std::array<unsigned,4>> children ;
  } ;

  /// Require a finite, strictly convex planar triangle or quad.
  void validatePlanarPolygon(const std::vector<Point> &corners) ;
  /// Construct quads from accepted boundary midpoints and the selected physical center.
  /// Integration limits apply only to the integrated center; unresolved integration throws.
  PlanarSplit splitPlanarPolygon(const std::vector<Point> &corners,
                                 const std::vector<Point> &midpoints,
                                 CenterMethod method = CenterMethod::weighted,
                                 MeasureOptions options = {}) ;
  /// Validate accepted points and child connectivity without choosing a new center.
  void validatePlanarSplit(const PlanarSplit &split, const std::vector<Point> &parentCorners) ;

}
}

#endif
