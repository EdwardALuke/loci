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

#include <FVMAdapt/parametric/reference_element.h>

#include <cmath>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {
  namespace {

    const ReferenceElement pointElement = {
      Shape::point, 0, 1, {{{{0,0,0}}}}, 0, {{}}, 0, {{}}
    } ;

    const ReferenceElement segmentElement = {
      Shape::segment, 1, 2, {{{{0,0,0}}, {{1,0,0}}}},
      1, {{{{0,1}}}},
      2, {{{Shape::point,1,{{0}}}, {Shape::point,1,{{1}}}}}
    } ;

    const ReferenceElement triangleElement = {
      Shape::triangle, 2, 3, {{{{0,0,0}}, {{1,0,0}}, {{0,1,0}}}},
      3, {{{{0,1}}, {{1,2}}, {{2,0}}}},
      3, {{{Shape::segment,2,{{0,1}}}, {Shape::segment,2,{{1,2}}},
           {Shape::segment,2,{{2,0}}}}}
    } ;

    const ReferenceElement quadElement = {
      Shape::quad, 2, 4, {{{{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}}}},
      4, {{{{0,1}}, {{1,2}}, {{2,3}}, {{3,0}}}},
      4, {{{Shape::segment,2,{{0,1}}}, {Shape::segment,2,{{1,2}}},
           {Shape::segment,2,{{2,3}}}, {Shape::segment,2,{{3,0}}}}}
    } ;

    const ReferenceElement tetElement = {
      Shape::tet, 3, 4, {{{{0,0,0}}, {{1,0,0}}, {{0,1,0}}, {{0,0,1}}}},
      6, {{{{0,1}}, {{1,2}}, {{2,0}}, {{0,3}}, {{1,3}}, {{2,3}}}},
      4, {{{Shape::triangle,3,{{0,2,1}}}, {Shape::triangle,3,{{0,1,3}}},
           {Shape::triangle,3,{{1,2,3}}}, {Shape::triangle,3,{{2,0,3}}}}}
    } ;

    const ReferenceElement hexElement = {
      Shape::hex, 3, 8,
      {{{{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}},
        {{0,0,1}}, {{1,0,1}}, {{1,1,1}}, {{0,1,1}}}},
      12, {{{{0,1}}, {{1,2}}, {{2,3}}, {{3,0}}, {{4,5}}, {{5,6}},
            {{6,7}}, {{7,4}}, {{0,4}}, {{1,5}}, {{2,6}}, {{3,7}}}},
      6, {{{Shape::quad,4,{{0,3,2,1}}}, {Shape::quad,4,{{4,5,6,7}}},
           {Shape::quad,4,{{0,1,5,4}}}, {Shape::quad,4,{{1,2,6,5}}},
           {Shape::quad,4,{{2,3,7,6}}}, {Shape::quad,4,{{3,0,4,7}}}}}
    } ;

    const ReferenceElement prismElement = {
      Shape::prism, 3, 6,
      {{{{0,0,0}}, {{1,0,0}}, {{0,1,0}}, {{0,0,1}}, {{1,0,1}}, {{0,1,1}}}},
      9, {{{{0,1}}, {{1,2}}, {{2,0}}, {{3,4}}, {{4,5}}, {{5,3}},
           {{0,3}}, {{1,4}}, {{2,5}}}},
      5, {{{Shape::triangle,3,{{0,2,1}}}, {Shape::triangle,3,{{3,4,5}}},
           {Shape::quad,4,{{0,1,4,3}}}, {Shape::quad,4,{{1,2,5,4}}},
           {Shape::quad,4,{{2,0,3,5}}}}}
    } ;

    const ReferenceElement pyramidElement = {
      Shape::pyramid, 3, 5,
      {{{{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}}, {{0,0,1}}}},
      8, {{{{0,1}}, {{1,2}}, {{2,3}}, {{3,0}}, {{0,4}}, {{1,4}}, {{2,4}}, {{3,4}}}},
      5, {{{Shape::quad,4,{{0,3,2,1}}}, {Shape::triangle,3,{{0,1,4}}},
           {Shape::triangle,3,{{1,2,4}}}, {Shape::triangle,3,{{2,3,4}}},
           {Shape::triangle,3,{{3,0,4}}}}}
    } ;

    bool finite(const Point &point) {
      return std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]) ;
    }

    std::array<double,4> boundaryWeights(Shape shape, const Point &point) {
      if(!finite(point))
        throw std::invalid_argument("reference coordinates must be finite") ;
      const double x = point[0], y = point[1] ;
      switch(shape) {
      case Shape::point: return {{1,0,0,0}} ;
      case Shape::segment: return {{1-x,x,0,0}} ;
      case Shape::triangle: return {{1-x-y,x,y,0}} ;
      case Shape::quad: return {{(1-x)*(1-y),x*(1-y),x*y,(1-x)*y}} ;
      default: throw std::invalid_argument("expected a point, segment, triangle or quad") ;
      }
    }

  }

  const ReferenceElement &referenceElement(Shape shape) {
    switch(shape) {
    case Shape::point: return pointElement ;
    case Shape::segment: return segmentElement ;
    case Shape::triangle: return triangleElement ;
    case Shape::quad: return quadElement ;
    case Shape::tet: return tetElement ;
    case Shape::hex: return hexElement ;
    case Shape::prism: return prismElement ;
    case Shape::pyramid: return pyramidElement ;
    case Shape::polygon: throw std::invalid_argument("polygon reference topology requires its vertex count") ;
    }
    throw std::invalid_argument("unknown reference shape") ;
  }

  bool contains(Shape shape, const Point &point, double tolerance) {
    const ReferenceElement &element = referenceElement(shape) ;
    if(!std::isfinite(tolerance) || tolerance < 0)
      throw std::invalid_argument("reference tolerance must be finite and nonnegative") ;
    if(!finite(point)) return false ;
    for(unsigned i = 0; i < element.dimension; ++i)
      if(point[i] < -tolerance) return false ;
    for(unsigned i = element.dimension; i < 3; ++i)
      if(std::abs(point[i]) > tolerance) return false ;
    const double x = point[0], y = point[1], z = point[2], upper = 1+tolerance ;
    switch(shape) {
    case Shape::point: return true ;
    case Shape::segment: return x <= upper ;
    case Shape::triangle: return x+y <= upper ;
    case Shape::quad: return x <= upper && y <= upper ;
    case Shape::tet: return x+y+z <= upper ;
    case Shape::hex: return x <= upper && y <= upper && z <= upper ;
    case Shape::prism: return x+y <= upper && z <= upper ;
    case Shape::pyramid: return z <= upper && x <= upper-z && y <= upper-z ;
    case Shape::polygon: break ;
    }
    return false ;
  }

  Point boundaryPoint(Shape shape, unsigned boundary, const Point &point) {
    const ReferenceElement &element = referenceElement(shape) ;
    if(boundary >= element.boundaryCount)
      throw std::out_of_range("reference boundary index") ;
    const Boundary &face = element.boundaries[boundary] ;
    const std::array<double,4> weights = boundaryWeights(face.shape,point) ;
    Point result = {{0,0,0}} ;
    for(unsigned i = 0; i < face.vertexCount; ++i)
      for(unsigned d = 0; d < 3; ++d)
        result[d] += weights[i]*element.vertices[face.vertices[i]][d] ;
    return result ;
  }

  unsigned orientationCount(Shape shape) {
    switch(shape) {
    case Shape::point: return 1 ;
    case Shape::segment: return 2 ;
    case Shape::triangle: return 6 ;
    case Shape::quad: return 8 ;
    default: throw std::invalid_argument("orientations require a point, segment, triangle or quad") ;
    }
  }

  unsigned orientedVertex(unsigned vertices, unsigned orientation, unsigned vertex) {
    if(vertices<3 || orientation/vertices>=2 || vertex>=vertices)
      throw std::out_of_range("polygon orientation or vertex index") ;
    return orientation<vertices ? (vertex+orientation)%vertices : (orientation-vertex)%vertices ;
  }

  unsigned orientedVertex(Shape shape, unsigned orientation, unsigned vertex) {
    if(orientation >= orientationCount(shape))
      throw std::out_of_range("reference orientation index") ;
    const unsigned count = referenceElement(shape).vertexCount ;
    if(vertex >= count) throw std::out_of_range("reference vertex index") ;
    if(orientation < count) return (vertex+orientation)%count ;
    return (orientation-vertex)%count ;
  }

  Point orientPoint(Shape shape, unsigned orientation, const Point &point) {
    if(orientation >= orientationCount(shape))
      throw std::out_of_range("reference orientation index") ;
    const ReferenceElement &element = referenceElement(shape) ;
    const std::array<double,4> weights = boundaryWeights(shape,point) ;
    Point result = {{0,0,0}} ;
    for(unsigned i = 0; i < element.vertexCount; ++i) {
      const Point &vertex = element.vertices[orientedVertex(shape,orientation,i)] ;
      for(unsigned d = 0; d < 3; ++d) result[d] += weights[i]*vertex[d] ;
    }
    return result ;
  }

  unsigned inverseOrientation(Shape shape, unsigned orientation) {
    if(orientation >= orientationCount(shape))
      throw std::out_of_range("reference orientation index") ;
    const unsigned count = referenceElement(shape).vertexCount ;
    return orientation < count ? (count-orientation)%count : orientation ;
  }

}
}
