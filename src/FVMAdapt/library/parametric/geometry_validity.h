// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_GEOMETRY_VALIDITY_INTERNAL_H
#define LOCI_FVMADAPT_GEOMETRY_VALIDITY_INTERNAL_H
#include "exact_geometry.h"
#include <FVMAdapt/parametric/geometry_map.h>

namespace Loci {
namespace FVMAdapt {
namespace detail {
  /// Algebra retained by temporary piece evaluators; unknown callbacks have no view.
  struct ExactMapView {
    enum class Kind { unknown, affine, delegate, coons } ;
    Kind kind = Kind::unknown ;
    ExactPoint origin{} ;
    std::array<ExactPoint,3> columns{} ;
    const GeometryMap *delegate = nullptr ;
    const std::array<std::shared_ptr<const GeometryMap>,4> *edges = nullptr ;
    const std::array<Point,4> *corners = nullptr ;
  } ;
  ExactMapView exactMapView(const GeometryMap &) ;
  /// Transform a segment only when the retained mapping is exactly affine along it.
  bool affineTrace(const GeometryMap &, std::array<ExactPoint,2> &, std::size_t &remaining) ;

  using BernsteinCoefficients = std::array<Exact,27> ;
  template<class Number> std::array<Number,3> jacobianDifference(const std::array<Number,3> &a,const std::array<Number,3> &b) {
    return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}} ;
  }
  template<class Number> std::array<Number,3> jacobianMidpoint(const std::array<Number,3> &a,const std::array<Number,3> &b) {
    return {{(a[0]+b[0])/2,(a[1]+b[1])/2,(a[2]+b[2])/2}} ;
  }
  template<class Number> Number jacobianDeterminant(const std::array<Number,3> &a,const std::array<Number,3> &b,const std::array<Number,3> &c) {
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]) ;
  }
  inline unsigned jacobianIndex(unsigned x, unsigned y, unsigned z) {
    return x+3*y+9*z ;
  }

  template<class Number> std::array<Number,27> hexJacobianCoefficients(const std::vector<Point> &corners) {
    const auto &element = referenceElement(Shape::hex) ;
    // Tensor order: x+2*y+4*z.
    std::array<std::array<Number,3>,8> points ;
    for(unsigned v = 0; v < 8; ++v) {
      const auto &q = element.vertices[v] ;
      points[unsigned(q[0])+2*unsigned(q[1])+4*unsigned(q[2])] = {{Number(corners[v][0]),Number(corners[v][1]),Number(corners[v][2])}} ;
    }

    // Each derivative is bilinear in the other two coordinates. Its corner
    // values are edge differences; midpoints fill the remaining 3x3 samples.
    std::array<std::array<std::array<Number,3>,9>,3> derivatives ;
    for(unsigned axis = 0; axis < 3; ++axis) {
      const unsigned u = (axis+1)%3, v = (axis+2)%3 ;
      auto &samples = derivatives[axis] ;
      for(unsigned j = 0; j < 2; ++j) for(unsigned i = 0; i < 2; ++i) {
        const unsigned start = (i<<u)+(j<<v) ;
        samples[2*i+6*j] = jacobianDifference(points[start+(1u<<axis)],points[start]) ;
      }
      for(unsigned row : {0u,6u}) samples[row+1] = jacobianMidpoint(samples[row],samples[row+2]) ;
      for(unsigned i = 0; i < 3; ++i) samples[i+3] = jacobianMidpoint(samples[i],samples[i+6]) ;
    }

    std::array<Number,27> coefficients{} ;
    // A trilinear map's determinant has degree at most two in each coordinate.
    // Exact values on the 3x3x3 tensor grid determine that polynomial completely.
    for(unsigned i = 0; i < 27; ++i) {
      const unsigned x = i%3, y = i/3%3, z = i/9 ;
      coefficients[i] = jacobianDeterminant(derivatives[0][y+3*z],derivatives[1][z+3*x],derivatives[2][x+3*y]) ;
    }
    // Convert the tensor quadratic to Bernstein coefficients, which bound its full range.
    for(unsigned step : {1u,3u,9u})
      for(unsigned i = 0; i < 27; ++i) if(i/step%3 == 0)
        coefficients[i+step] = 2*coefficients[i+step]-(coefficients[i]+coefficients[i+2*step])/2 ;
    return coefficients ;
  }

  inline std::array<BernsteinCoefficients,2> bisectJacobianCoefficients(const BernsteinCoefficients &coefficients, unsigned axis) {
    std::array<BernsteinCoefficients,2> result{} ;
    const unsigned step = axis == 0 ? 1 : axis == 1 ? 3 : 9 ;
    for(unsigned i = 0; i < 27; ++i) if(i/step%3 == 0) {
      const Exact &a = coefficients[i], &b = coefficients[i+step], &c = coefficients[i+2*step] ;
      const Exact center = (a+2*b+c)/4 ;
      result[0][i] = a ; result[0][i+step] = (a+b)/2 ; result[0][i+2*step] = center ;
      result[1][i] = center ; result[1][i+step] = (b+c)/2 ; result[1][i+2*step] = c ;
    }
    return result ;
  }

}
}
}
#endif
