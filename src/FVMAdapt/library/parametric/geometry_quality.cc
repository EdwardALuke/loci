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
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>.
//#
//#############################################################################

#include <FVMAdapt/parametric/geometry_quality.h>
#include "geometry_bounds.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {
  bool positiveHexJacobian(const std::vector<Point> &corners) {
    if(corners.size() != 8) throw std::invalid_argument("hex Jacobian bounds require eight corners") ;
    using detail::Bounds ;
    // Tensor order x+2*y+4*z. Each derivative is a convex blend of four edge differences.
    constexpr unsigned vertex[8] = {0,1,3,2,4,5,7,6} ;
    std::array<std::array<Bounds,3>,3> columns ;
    for(unsigned axis = 0; axis < 3; ++axis) {
      const unsigned u = (axis+1)%3, v = (axis+2)%3 ;
      for(unsigned j = 0; j < 2; ++j) for(unsigned i = 0; i < 2; ++i) {
        const unsigned start = (i<<u)+(j<<v) ;
        for(unsigned d = 0; d < 3; ++d) {
          const double a = corners[vertex[start]][d], b = corners[vertex[start+(1u<<axis)]][d] ;
          if(!std::isfinite(a) || !std::isfinite(b)) return false ;
          const auto edge = Bounds(b)-Bounds(a) ; auto &range = columns[axis][d] ;
          if(i == 0 && j == 0) range = edge ;
          else { range.low = std::min(range.low,edge.low) ; range.high = std::max(range.high,edge.high) ; }
        }
      }
    }
    const auto &a = columns[0], &b = columns[1], &c = columns[2] ;
    const auto determinant = a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]) ;
    return determinant.low > 0 && std::isfinite(determinant.high) ;
  }

  namespace {

    using WidePoint = std::array<long double,3> ;

    std::vector<Point> defaultSamples(Shape shape) {
      std::vector<Point> samples ;
      if(shape == Shape::triangle || shape == Shape::tet || shape == Shape::prism) {
        const unsigned layers = shape == Shape::triangle ? 1 : 3 ;
        for(unsigned k = 0; k < layers; ++k)
          for(unsigned j = 0; j <= 2; ++j)
            for(unsigned i = 0; i+j <= 2; ++i)
              if(shape != Shape::tet || i+j+k <= 2) samples.push_back({{0.5*i,0.5*j,0.5*k}}) ;
        if(shape == Shape::triangle) samples.push_back({{1./3,1./3,0}}) ;
        else if(shape == Shape::tet) samples.push_back({{0.25,0.25,0.25}}) ;
        else for(double z : {0.,0.5,1.}) samples.push_back({{1./3,1./3,z}}) ;
      } else if(shape == Shape::pyramid) {
        for(double z : {0.,0.5,0.9375})
          for(double y : {0.,0.5,1.})
            for(double x : {0.,0.5,1.}) samples.push_back({{x*(1-z),y*(1-z),z}}) ;
      } else {
        const unsigned dimension = referenceElement(shape).dimension ;
        for(unsigned k = 0; k < (dimension == 3 ? 3u : 1u); ++k)
          for(unsigned j = 0; j < (dimension >= 2 ? 3u : 1u); ++j)
            for(unsigned i = 0; i < 3; ++i) samples.push_back({{0.5*i,0.5*j,0.5*k}}) ;
      }
      return samples ;
    }

    bool finite(const Point &point) {
      for(double value : point) if(!std::isfinite(value)) return false ;
      return true ;
    }

    long double dot(const WidePoint &a, const WidePoint &b) {
      return a[0]*b[0]+a[1]*b[1]+a[2]*b[2] ;
    }

    WidePoint cross(const WidePoint &a, const WidePoint &b) {
      return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}} ;
    }

    void recordProblem(GeometryQuality &quality, const Point &point, bool invalid) {
      if(quality.invalidCount == 0 && quality.unresolvedCount == 0) quality.problemPoint = point ;
      if(invalid) { ++quality.invalidCount ; quality.status = GeometryQualityStatus::invalid ; }
      else {
        ++quality.unresolvedCount ;
        if(quality.status != GeometryQualityStatus::invalid) quality.status = GeometryQualityStatus::unresolved ;
      }
    }

    void measure(GeometryQuality &quality, const MapValue &value, const Point &point, unsigned dimension) {
      if(!finite(value.position)) { recordProblem(quality,point,true) ; return ; }
      if(!value.hasJacobian) { recordProblem(quality,point,false) ; return ; }
      std::array<WidePoint,3> tangents{} ;
      std::array<long double,3> norms{}, lengths{} ;
      for(unsigned d = 0; d < dimension; ++d) {
        if(!finite(value.jacobian[d])) { recordProblem(quality,point,true) ; return ; }
        long double scale = 0 ;
        for(double component : value.jacobian[d]) scale = std::max(scale,std::abs(static_cast<long double>(component))) ;
        if(scale == 0) { quality.collapsed = true ; recordProblem(quality,point,true) ; return ; }
        // Normalize each tangent separately to preserve very small and very large length scales.
        for(unsigned i = 0; i < 3; ++i) tangents[d][i] = value.jacobian[d][i]/scale ;
        norms[d] = std::sqrt(dot(tangents[d],tangents[d])) ;
        lengths[d] = scale*norms[d] ;
      }

      long double jacobian = 1, alignment = 0 ;
      for(unsigned d = 0; d < dimension; ++d)
        for(unsigned e = d+1; e < dimension; ++e)
          alignment = std::max(alignment,std::abs(dot(tangents[d],tangents[e]))/(norms[d]*norms[e])) ;
      if(dimension >= 2) {
        const auto normal = cross(tangents[0],tangents[1]) ;
        jacobian = dimension == 2 ? std::sqrt(dot(normal,normal))/(norms[0]*norms[1]) :
          dot(normal,tangents[2])/(norms[0]*norms[1]*norms[2]) ;
      }
      if(jacobian <= 0) {
        quality.inverted = quality.inverted || jacobian < 0 ;
        quality.collapsed = quality.collapsed || jacobian == 0 ;
        recordProblem(quality,point,true) ;
      }
      const auto limits = std::minmax_element(lengths.begin(),lengths.begin()+dimension) ;
      const long double ratio = *limits.second / *limits.first ;
      if(!std::isfinite(ratio) || ratio > std::numeric_limits<double>::max()) {
        recordProblem(quality,point,false) ; return ;
      }
      const double measuredJacobian = static_cast<double>(std::max(-1.L,std::min(1.L,jacobian))) ;
      if(jacobian != 0 && measuredJacobian == 0) { recordProblem(quality,point,false) ; return ; }
      const double measuredAlignment = static_cast<double>(std::min(1.L,alignment)) ;
      if(quality.measuredCount == 0 || measuredJacobian < quality.minScaledJacobian) {
        quality.minScaledJacobian = measuredJacobian ; quality.scaledJacobianPoint = point ;
      }
      if(quality.measuredCount == 0 || measuredAlignment > quality.maxAlignment) {
        quality.maxAlignment = measuredAlignment ; quality.alignmentPoint = point ;
      }
      if(quality.measuredCount == 0 || ratio > quality.maxDerivativeRatio) {
        quality.maxDerivativeRatio = static_cast<double>(ratio) ; quality.derivativeRatioPoint = point ;
      }
      ++quality.measuredCount ;
    }

  }

  GeometryQuality geometryQuality(const GeometryMap &geometry, const GeometryQualityOptions &options) {
    const Shape shape = geometry.shape() ;
    const unsigned dimension = geometry.dimension() ;
    if(dimension == 0) throw std::invalid_argument("geometry quality requires a positive-dimensional shape") ;
    if(options.maxSamples == 0 || options.samples.size() > options.maxSamples)
      throw std::invalid_argument("geometry quality sample count exceeds its limit") ;
    std::vector<Point> defaults ;
    if(options.samples.empty()) {
      if(shape!=Shape::polygon) defaults=defaultSamples(shape) ;
      else {
        const auto corners=geometry.referenceCorners() ; defaults.assign(corners.begin(),corners.end()) ; const auto n=defaults.size() ;
        for(unsigned i=0;i<n;++i) defaults.push_back({{.25+.5*defaults[i][0],.25+.5*defaults[i][1],0}}) ;
        defaults.push_back({{.5,.5,0}}) ;
      }
    }
    const auto &samples = options.samples.empty() ? defaults : options.samples ;
    if(samples.size() > options.maxSamples) throw std::invalid_argument("geometry quality sample count exceeds its limit") ;
    for(const auto &point : samples)
      if(!geometry.containsReference(point,shape==Shape::polygon ? 64*std::numeric_limits<double>::epsilon():0))
        throw std::invalid_argument("geometry quality sample must lie in the reference element") ;

    GeometryQuality quality ; quality.sampleCount = samples.size() ;
    for(const auto &point : samples) {
      try { measure(quality,geometry.evaluate(point),point,dimension) ; }
      catch(const std::bad_alloc &) { throw ; }
      catch(const std::exception &) { recordProblem(quality,point,false) ; }
    }
    return quality ;
  }

}
}
