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

#include <FVMAdapt/parametric/reference_coordinate.h>
#include "exact_geometry.h"

#include <stdexcept>
#include <string>

namespace Loci {
namespace FVMAdapt {
  namespace {

    std::size_t bitCount(const Fraction::Integer &value) {
      if(value == 0) return 0 ;
      const Fraction::Integer magnitude = value < 0 ? -value : value ;
      return boost::multiprecision::msb(magnitude)+1 ;
    }

    void charge(std::size_t &used, std::size_t count, std::size_t limit, const char *resource) {
      if(count > limit-used)
        throw std::length_error(std::string("exact construction exceeds ")+resource+" limit: used="+
          std::to_string(used)+", requested="+std::to_string(count)+", limit="+std::to_string(limit)) ;
      used += count ;
    }

    ReferencePoint difference(const ReferencePoint &a, const ReferencePoint &b) {
      return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}} ;
    }

    // Every native boundary and boundary orientation is affine in reference space.
    ReferencePoint affinePoint(const ReferencePoint &origin,
                               const ReferencePoint &u, const ReferencePoint &v,
                               const ReferencePoint &point) {
      ReferencePoint result ;
      for(unsigned d = 0; d < 3; ++d)
        result[d] = Fraction(detail::exact(origin[d])+detail::exact(point[0])*detail::exact(u[d])
                                  +detail::exact(point[1])*detail::exact(v[d])) ;
      return result ;
    }

    const Boundary &checkedBoundary(const ReferenceElement &element, unsigned boundary) {
      if(boundary >= element.boundaryCount)
        throw std::out_of_range("reference boundary index") ;
      return element.boundaries[boundary] ;
    }

    void boundaryBasis(Shape shape, const Boundary &boundary, ReferencePoint &origin,
                       ReferencePoint &u, ReferencePoint &v) {
      origin = referenceVertex(shape,boundary.vertices[0]) ;
      u = ReferencePoint() ;
      v = ReferencePoint() ;
      if(boundary.vertexCount > 1)
        u = difference(referenceVertex(shape,boundary.vertices[1]),origin) ;
      if(boundary.vertexCount > 2)
        v = difference(referenceVertex(shape,boundary.vertices[boundary.vertexCount-1]),origin) ;
    }

  }

  Fraction::Fraction(std::int64_t numerator, std::int64_t denominator)
    : Fraction(Integer(numerator),Integer(denominator)) {}

  Fraction::Fraction(const Integer &numerator, const Integer &denominator) {
    if(denominator == 0) throw std::invalid_argument("fraction denominator is zero") ;
    value_ = Rational(numerator,denominator) ;
  }

  Fraction::Fraction(const Rational &value)
    : Fraction(boost::multiprecision::numerator(value),boost::multiprecision::denominator(value)) {}

  void checkExact(const Fraction &value, const ExactLimits &limits) {
    if(bitCount(value.numerator()) > limits.maxBits || bitCount(value.denominator()) > limits.maxBits)
      throw std::length_error("exact coordinate exceeds its bit limit") ;
  }

  void ExactBudget::operations(std::size_t count) {
    charge(operations_,count,limits_.maxOperations,"operation") ;
  }

  void ExactBudget::retain(std::size_t bytes) {
    charge(bytes_,bytes,limits_.maxBytes,"byte") ;
  }

  Fraction operator+(const Fraction &a, const Fraction &b) {
    return Fraction(Fraction::Rational(a.rational()+b.rational())) ;
  }

  Fraction operator-(const Fraction &a, const Fraction &b) {
    return Fraction(Fraction::Rational(a.rational()-b.rational())) ;
  }

  Fraction operator-(const Fraction &a) {
    return Fraction(Fraction::Rational(-a.rational())) ;
  }

  Fraction operator*(const Fraction &a, const Fraction &b) {
    return Fraction(Fraction::Rational(a.rational()*b.rational())) ;
  }

  Fraction operator/(const Fraction &a, const Fraction &b) {
    if(b.rational() == 0) throw std::domain_error("division by zero fraction") ;
    return Fraction(Fraction::Rational(a.rational()/b.rational())) ;
  }

  bool operator==(const Fraction &a, const Fraction &b) {
    return a.rational() == b.rational() ;
  }

  bool operator!=(const Fraction &a, const Fraction &b) { return !(a == b) ; }

  bool operator<(const Fraction &a, const Fraction &b) {
    return a.rational() < b.rational() ;
  }

  bool operator>(const Fraction &a, const Fraction &b) { return b < a ; }
  bool operator<=(const Fraction &a, const Fraction &b) { return !(b < a) ; }
  bool operator>=(const Fraction &a, const Fraction &b) { return !(a < b) ; }

  Point approximate(const ReferencePoint &point) {
    return {{point[0].value(),point[1].value(),point[2].value()}} ;
  }

  ReferencePoint referenceVertex(Shape shape, unsigned vertex) {
    const ReferenceElement &element = referenceElement(shape) ;
    if(vertex >= element.vertexCount) throw std::out_of_range("reference vertex index") ;
    const Point &point = element.vertices[vertex] ;
    return {{Fraction(std::int64_t(point[0])),Fraction(std::int64_t(point[1])),
             Fraction(std::int64_t(point[2]))}} ;
  }

  bool containsExact(Shape shape, const ReferencePoint &point) {
    const ReferenceElement &element = referenceElement(shape) ;
    for(unsigned d = 0; d < element.dimension; ++d)
      if(point[d] < 0 || point[d] > 1) return false ;
    for(unsigned d = element.dimension; d < 3; ++d)
      if(point[d] != 0) return false ;
    const Fraction &x = point[0], &y = point[1], &z = point[2] ;
    switch(shape) {
    case Shape::point:
    case Shape::segment:
    case Shape::quad:
    case Shape::hex: return true ;
    case Shape::triangle:
    case Shape::prism: return x <= 1-y ;
    case Shape::tet: return x+y+z <= 1 ;
    case Shape::pyramid: return x <= 1-z && y <= 1-z ;
    case Shape::polygon: break ;
    }
    return false ;
  }

  bool onBoundary(Shape shape, unsigned boundary, const ReferencePoint &point) {
    const ReferenceElement &element = referenceElement(shape) ;
    const Boundary &face = checkedBoundary(element,boundary) ;
    if(!containsExact(shape,point)) return false ;
    ReferencePoint origin, u, v ;
    boundaryBasis(shape,face,origin,u,v) ;
    if(element.dimension == 1) return point == origin ;

    ReferencePoint normal = {{u[1],-u[0],0}} ;
    if(element.dimension == 3)
      for(unsigned d = 0; d < 3; ++d) {
        const unsigned a = (d+1)%3, b = (d+2)%3 ;
        normal[d] = u[a]*v[b]-u[b]*v[a] ;
      }
    Fraction residual ;
    for(unsigned d = 0; d < 3; ++d) residual = residual+normal[d]*(point[d]-origin[d]) ;
    return residual == 0 ;
  }

  ReferencePoint boundaryReferencePoint(Shape shape, unsigned boundary,
                                       const ReferencePoint &point) {
    const ReferenceElement &element = referenceElement(shape) ;
    const Boundary &face = checkedBoundary(element,boundary) ;
    if(!containsExact(face.shape,point))
      throw std::invalid_argument("point is outside the reference boundary") ;
    ReferencePoint origin, u, v ;
    boundaryBasis(shape,face,origin,u,v) ;
    return affinePoint(origin,u,v,point) ;
  }

  ReferencePoint boundaryCoordinates(Shape shape, unsigned boundary,
                                     const ReferencePoint &point) {
    const ReferenceElement &element = referenceElement(shape) ;
    const Boundary &face = checkedBoundary(element,boundary) ;
    if(!onBoundary(shape,boundary,point))
      throw std::invalid_argument("point is not on the chosen reference boundary") ;
    ReferencePoint origin, u, v, result = {{0,0,0}} ;
    boundaryBasis(shape,face,origin,u,v) ;
    detail::ExactPoint offset, eu, ev ;
    for(unsigned d = 0; d < 3; ++d) {
      offset[d] = detail::exact(point[d])-detail::exact(origin[d]) ;
      eu[d] = detail::exact(u[d]) ; ev[d] = detail::exact(v[d]) ;
    }
    const unsigned dimension = referenceElement(face.shape).dimension ;
    if(dimension == 1) {
      for(unsigned d = 0; d < 3; ++d)
        if(u[d] != 0) {
          result[0] = Fraction(offset[d]/eu[d]) ;
          break ;
        }
    } else if(dimension == 2) {
      bool found = false ;
      for(unsigned a = 0; a < 3 && !found; ++a)
        for(unsigned b = a+1; b < 3 && !found; ++b) {
          const detail::Exact determinant = eu[a]*ev[b]-eu[b]*ev[a] ;
          if(determinant == 0) continue ;
          result[0] = Fraction((offset[a]*ev[b]-offset[b]*ev[a])/determinant) ;
          result[1] = Fraction((eu[a]*offset[b]-eu[b]*offset[a])/determinant) ;
          found = true ;
        }
    }
    if(!containsExact(face.shape,result) || affinePoint(origin,u,v,result) != point)
      throw std::invalid_argument("point is not on the chosen reference boundary") ;
    return result ;
  }

  ReferencePoint orientReferencePoint(Shape shape, unsigned orientation,
                                     const ReferencePoint &point) {
    if(orientation >= orientationCount(shape))
      throw std::out_of_range("reference orientation index") ;
    if(!containsExact(shape,point))
      throw std::invalid_argument("point is outside the oriented reference shape") ;
    const ReferenceElement &element = referenceElement(shape) ;
    const ReferencePoint origin = referenceVertex(shape,orientedVertex(shape,orientation,0)) ;
    ReferencePoint u = {{0,0,0}}, v = {{0,0,0}} ;
    if(element.vertexCount > 1)
      u = difference(referenceVertex(shape,orientedVertex(shape,orientation,1)),origin) ;
    if(element.vertexCount > 2)
      v = difference(referenceVertex(shape,orientedVertex(shape,orientation,element.vertexCount-1)),origin) ;
    return affinePoint(origin,u,v,point) ;
  }

}
}
