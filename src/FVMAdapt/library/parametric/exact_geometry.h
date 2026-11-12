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

#ifndef LOCI_FVMADAPT_EXACT_GEOMETRY_H
#define LOCI_FVMADAPT_EXACT_GEOMETRY_H

#include <FVMAdapt/parametric/reference_coordinate.h>
#include <vector>

namespace Loci {
namespace FVMAdapt {
namespace detail {

  using Exact = Fraction::Rational ;
  using ExactPoint = std::array<Exact,3> ;

  inline const Exact &exact(const Fraction &value) {
    return value.rational() ;
  }

  /// Whether a quad has two consecutive edges continuing along exactly one line.
  inline bool hasStraightQuadCorner(const std::vector<Point> &points) {
    if(points.size()!=4) return false ;
    for(unsigned i=0;i<4;++i) {
      ExactPoint before,after ; Exact dot=0 ;
      for(unsigned d=0;d<3;++d) {
        before[d]=Exact(points[i][d])-Exact(points[(i+3)%4][d]) ;
        after[d]=Exact(points[(i+1)%4][d])-Exact(points[i][d]) ;
        dot+=before[d]*after[d] ;
      }
      if(dot>0 && before[1]*after[2]==before[2]*after[1] &&
         before[2]*after[0]==before[0]*after[2] && before[0]*after[1]==before[1]*after[0]) return true ;
    }
    return false ;
  }

}
}
}

#endif
