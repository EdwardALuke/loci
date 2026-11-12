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

#include <limits>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {

  std::vector<Child> detail::splitChildren(const PrismSplit &request) {
    const std::vector<Child> base = splitChildren(request.base) ;
    if(request.axialCuts.size() >= std::numeric_limits<unsigned>::max()/base.size())
      throw std::invalid_argument("prism subdivision has too many children") ;
    Fraction previous = 0 ;
    for(const Fraction &cut : request.axialCuts) {
      if(cut <= previous || cut >= Fraction(1))
        throw std::invalid_argument("prism axial cuts must increase strictly inside (0,1)") ;
      previous = cut ;
    }
    std::vector<Child> children ;
    const std::size_t count = (request.axialCuts.size()+1)*base.size() ;
    if(count > children.max_size()) throw std::invalid_argument("prism subdivision exceeds container capacity") ;
    children.reserve(count) ;
    Fraction lower = 0 ;
    for(std::size_t layer = 0; layer <= request.axialCuts.size(); ++layer) {
      const Fraction upper = layer < request.axialCuts.size() ? request.axialCuts[layer] : Fraction(1) ;
      for(const Child &triangle : base) {
        Child prism = triangle ; prism.shape = Shape::prism ;
        prism.placement.origin[2] = lower ;
        prism.placement.columns[2][2] = Fraction(detail::exact(upper)-detail::exact(lower)) ;
        children.push_back(prism) ;
      }
      lower = upper ;
    }
    return children ;
  }

  SplitResult split(const PrismSplit &request) {
    SplitResult result ; result.shape = Shape::prism ;
    result.children = detail::splitChildren(request) ;
    detail::completeSubdivision(result) ; return result ;
  }

  SplitResult split(Shape shape, const PrismSplit &request) {
    if(shape != Shape::prism) throw std::invalid_argument("prism split requires a prism") ;
    return split(request) ;
  }

}
}
