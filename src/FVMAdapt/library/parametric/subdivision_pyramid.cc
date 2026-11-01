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

#include <FVMAdapt/parametric/subdivision.h>
#include "subdivision_internal.h"
#include "exact_geometry.h"

#include <stdexcept>

namespace Loci {
namespace FVMAdapt {
  namespace {

    Child pyramid(const ReferencePoint &a, const ReferencePoint &b,
                  const ReferencePoint &d, const ReferencePoint &apex) {
      Child child ; child.shape = Shape::pyramid ; child.placement.origin = a ;
      const std::array<ReferencePoint,3> ends = {{b,d,apex}} ;
      for(unsigned i = 0; i < 3; ++i)
        for(unsigned j = 0; j < 3; ++j)
          child.placement.columns[i][j] = Fraction(detail::exact(ends[i][j])-detail::exact(a[j])) ;
      return child ;
    }

  }

  PyramidSplit PyramidSplit::uniform10() {
    PyramidSplit request ; request.kind = Kind::uniform ; return request ;
  }

  std::vector<Child> detail::splitChildren(const PyramidSplit &request) {
    switch(request.kind) {
    case PyramidSplit::Kind::none: return {{Shape::pyramid,Placement()}} ;
    case PyramidSplit::Kind::uniform: break ;
    default: throw std::invalid_argument("unknown pyramid split pattern") ;
    }
    const Fraction half(1,2) ;
    const ReferencePoint a = {{0,0,0}}, b = {{1,0,0}}, d = {{0,1,0}}, apex = {{0,0,1}} ;
    const ReferencePoint ab = {{half,0,0}}, bc = {{1,half,0}}, cd = {{half,1,0}}, da = {{0,half,0}} ;
    const ReferencePoint ae = {{0,0,half}}, be = {{half,0,half}}, ce = {{half,half,half}}, de = {{0,half,half}} ;
    const ReferencePoint center = {{half,half,0}} ;
    // Edge midpoints and the base center form six pyramids with four tetrahedral gaps.
    return {pyramid(a,ab,da,ae), pyramid(ab,b,center,be),
            pyramid(center,bc,cd,ce), pyramid(da,center,d,de),
            pyramid(ae,be,de,apex), pyramid(ae,de,be,center),
            tetChild({{ab,center,ae,be}}), tetChild({{bc,center,be,ce}}),
            tetChild({{cd,center,ce,de}}), tetChild({{da,center,de,ae}})} ;
  }

  SplitResult split(const PyramidSplit &request) {
    SplitResult result ; result.shape = Shape::pyramid ;
    result.children = detail::splitChildren(request) ;
    detail::completeSubdivision(result) ; return result ;
  }

  SplitResult split(Shape shape, const PyramidSplit &request) {
    if(shape != Shape::pyramid) throw std::invalid_argument("pyramid split requires a pyramid") ;
    return split(request) ;
  }

}
}
