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

    Child triangle(const ReferencePoint &a, const ReferencePoint &b, const ReferencePoint &c) {
      Child child ; child.shape = Shape::triangle ; child.placement.origin = a ;
      for(unsigned d = 0; d < 3; ++d) {
        child.placement.columns[0][d] = Fraction(detail::exact(b[d])-detail::exact(a[d])) ;
        child.placement.columns[1][d] = Fraction(detail::exact(c[d])-detail::exact(a[d])) ;
      }
      return child ;
    }

  }

  TriangleSplit TriangleSplit::bisectEdge(unsigned edge, Fraction position) {
    TriangleSplit request ; request.kind = Kind::bisect ; request.edge = edge ; request.position = position ;
    return request ;
  }

  TriangleSplit TriangleSplit::uniform4() {
    TriangleSplit request ; request.kind = Kind::uniform ; return request ;
  }

  std::vector<Child> detail::splitChildren(const TriangleSplit &request) {
    std::vector<Child> children ;
    switch(request.kind) {
    case TriangleSplit::Kind::none:
      children.push_back({Shape::triangle,Placement()}) ; break ;
    case TriangleSplit::Kind::bisect: {
      if(request.edge >= 3 || request.position <= Fraction(0) || request.position >= Fraction(1))
        throw std::invalid_argument("triangle bisection requires an edge in [0,3) and position in (0,1)") ;
      const ReferencePoint a = referenceVertex(Shape::triangle,request.edge) ;
      const ReferencePoint b = referenceVertex(Shape::triangle,(request.edge+1)%3) ;
      const ReferencePoint c = referenceVertex(Shape::triangle,(request.edge+2)%3) ;
      ReferencePoint middle ;
      for(unsigned d = 0; d < 3; ++d)
        middle[d] = Fraction(detail::exact(a[d])+detail::exact(request.position)*
                                  (detail::exact(b[d])-detail::exact(a[d]))) ;
      children.push_back(triangle(a,middle,c)) ;
      children.push_back(triangle(middle,b,c)) ;
      break ;
    }
    case TriangleSplit::Kind::uniform: {
      const ReferencePoint a = {{0,0,0}}, b = {{1,0,0}}, c = {{0,1,0}} ;
      const ReferencePoint ab = {{Fraction(1,2),0,0}}, bc = {{Fraction(1,2),Fraction(1,2),0}}, ca = {{0,Fraction(1,2),0}} ;
      children.push_back(triangle(a,ab,ca)) ;
      children.push_back(triangle(ab,b,bc)) ;
      children.push_back(triangle(ca,bc,c)) ;
      children.push_back(triangle(ab,bc,ca)) ;
      break ;
    }
    default: throw std::invalid_argument("unknown triangle split pattern") ;
    }
    return children ;
  }

  SplitResult split(const TriangleSplit &request) {
    SplitResult result ; result.shape = Shape::triangle ;
    result.children = detail::splitChildren(request) ;
    detail::completeSubdivision(result) ; return result ;
  }

  SplitResult split(Shape shape, const TriangleSplit &request) {
    if(shape != Shape::triangle) throw std::invalid_argument("triangle split requires a triangle") ;
    return split(request) ;
  }

}
}
