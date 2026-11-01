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

#include "geometry_integration.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {
  namespace {

    using WidePoint = std::array<long double,3> ;

    struct Sum {
      long double value = 0, correction = 0 ;
      void add(long double term) {
        const long double next = value+term ;
        correction += std::abs(value) >= std::abs(term) ? (value-next)+term : (term-next)+value ;
        value = next ;
      }
      long double total() const { return value+correction ; }
    } ;

    struct GaussPoint { long double coordinate, weight ; } ;

    std::vector<GaussPoint> gaussRule(unsigned order) {
      std::vector<GaussPoint> rule(order) ;
      const long double pi = std::acos(-1.L) ;
      for(unsigned i = 0; i < (order+1)/2; ++i) {
        long double root = std::cos(pi*(i+.75L)/(order+.5L)), derivative = 0 ;
        bool resolved = false ;
        for(unsigned iteration = 0; iteration < 64; ++iteration) {
          long double previous = 1, polynomial = root ;
          for(unsigned degree = 2; degree <= order; ++degree) {
            const long double next = ((2*degree-1)*root*polynomial-(degree-1)*previous)/degree ;
            previous = polynomial ; polynomial = next ;
          }
          derivative = order*(root*polynomial-previous)/(root*root-1) ;
          const long double step = polynomial/derivative ;
          root -= step ;
          if(std::abs(step) <= 8*std::numeric_limits<long double>::epsilon()) { resolved = true ; break ; }
        }
        if(!resolved) throw std::length_error("Gaussian quadrature nodes did not converge") ;
        const long double weight = 1/((1-root*root)*derivative*derivative) ;
        rule[i] = {(1-root)/2,weight} ;
        rule[order-1-i] = {(1+root)/2,weight} ;
      }
      return rule ;
    }

    WidePoint cross(const WidePoint &a, const WidePoint &b) {
      return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}} ;
    }

    long double norm(const WidePoint &a) { return std::hypot(a[0],a[1],a[2]) ; }

    // Product Gaussian rules reach simplices and the pyramid through collapsed coordinates.
    Point referencePoint(Shape shape, long double u, long double v, long double w, long double &weight) {
      switch(shape) {
      case Shape::triangle: case Shape::prism:
        weight *= 1-u ; v *= 1-u ; break ;
      case Shape::tet:
        weight *= (1-u)*(1-u)*(1-v) ; w *= (1-u)*(1-v) ; v *= 1-u ; break ;
      case Shape::pyramid:
        weight *= (1-w)*(1-w) ; u *= 1-w ; v *= 1-w ; break ;
      default: break ;
      }
      return {{static_cast<double>(u),static_cast<double>(v),static_cast<double>(w)}} ;
    }

    struct Integral {
      long double measure = 0 ;
      WidePoint centroid{{0,0,0}} ;
      WidePoint moment{{0,0,0}} ;
    } ;

    Integral integrate(Shape shape, const std::function<detail::MeasureSample(const Point&)> &sample,
                       const std::vector<GaussPoint> &rule,
                       unsigned dimension, Point &origin, bool &hasOrigin, long double &size) {
      Sum measure ;
      std::array<Sum,3> moment ;
      const unsigned count = rule.size(), secondCount = dimension > 1 ? count : 1, thirdCount = dimension > 2 ? count : 1 ;
      for(unsigned i = 0; i < count; ++i) for(unsigned j = 0; j < secondCount; ++j) for(unsigned k = 0; k < thirdCount; ++k) {
        long double weight = rule[i].weight ;
        if(dimension > 1) weight *= rule[j].weight ;
        if(dimension > 2) weight *= rule[k].weight ;
        const Point q = referencePoint(shape,rule[i].coordinate,dimension > 1 ? rule[j].coordinate : 0,
                                        dimension > 2 ? rule[k].coordinate : 0,weight) ;
        const auto value = sample(q) ;
        for(double x: value.position) if(!std::isfinite(x))
          throw std::invalid_argument("geometry measure requires finite positions") ;
        if(!(value.density > 0) || !std::isfinite(value.density))
          throw std::invalid_argument("geometry measure requires positive nonsingular sampled Jacobians") ;
        if(!hasOrigin) { origin = value.position ; hasOrigin = true ; }
        WidePoint offset ;
        for(unsigned d = 0; d < 3; ++d) offset[d] = static_cast<long double>(value.position[d])-origin[d] ;
        size = std::max(size,norm(offset)) ;
        weight *= value.density ; measure.add(weight) ;
        for(unsigned d = 0; d < 3; ++d) moment[d].add(weight*offset[d]) ;
      }
      Integral result ; result.measure = measure.total() ;
      if(!(result.measure > 0) || !std::isfinite(result.measure))
        throw std::overflow_error("geometry measure integral is unrepresentable") ;
      for(unsigned d = 0; d < 3; ++d) {
        result.moment[d] = moment[d].total() ;
        result.centroid[d] = result.moment[d]/result.measure ;
        if(!std::isfinite(result.centroid[d])) throw std::overflow_error("geometry centroid is unrepresentable") ;
      }
      return result ;
    }

    double checkedDouble(long double value) {
      if(!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
        throw std::overflow_error("geometry measure result is unrepresentable") ;
      return static_cast<double>(value) ;
    }

  }

  detail::MeasureIntegral detail::integrateGeometryRule(Shape shape,
    const std::function<MeasureSample(const Point&)> &sample, unsigned order, const Point &origin) {
    const unsigned dimension = referenceElement(shape).dimension ;
    if(dimension == 0 || order < 2 || order > 256)
      throw std::invalid_argument("invalid geometry integration rule") ;
    Point fixedOrigin = origin ; bool hasOrigin = true ; long double radius = 0 ;
    const auto value = integrate(shape,sample,gaussRule(order),dimension,fixedOrigin,hasOrigin,radius) ;
    return {value.measure,value.moment,radius} ;
  }

  GeometryMeasure detail::integrateGeometry(Shape shape,
    const std::function<MeasureSample(const Point&)> &sample, MeasureOptions options) {
    const unsigned dimension = referenceElement(shape).dimension ;
    if(dimension == 0) throw std::invalid_argument("geometry measure requires a curve, surface or volume") ;
    if(!std::isfinite(options.relativeTolerance) || options.relativeTolerance <= 0 ||
       options.maxOrder < 2 || options.maxOrder > 256 || options.maxEvaluations == 0)
      throw std::invalid_argument("invalid geometry measure bounds or tolerance") ;

    GeometryMeasure result ;
    Integral previous ;
    Point origin{{0,0,0}} ;
    bool hasOrigin = false ;
    long double size = 0 ;
    for(unsigned order = 2; order <= options.maxOrder; order += 2) {
      std::size_t evaluations = order ;
      for(unsigned d = 1; d < dimension; ++d) evaluations *= order ;
      if(evaluations > options.maxEvaluations-result.evaluations) break ;
      const auto current = integrate(shape,sample,gaussRule(order),dimension,origin,hasOrigin,size) ;
      if(result.order != 0) {
        WidePoint change ;
        for(unsigned d = 0; d < 3; ++d) change[d] = current.centroid[d]-previous.centroid[d] ;
        long double magnitude = size ;
        for(double x: origin) magnitude = std::max(magnitude,std::abs(static_cast<long double>(x))) ;
        const long double roundoff = 8*std::numeric_limits<double>::epsilon()*magnitude ;
        const long double measureError = std::abs(current.measure-previous.measure), centerError = norm(change) ;
        result.estimatedMeasureError = checkedDouble(measureError) ;
        result.estimatedCentroidError = checkedDouble(centerError) ;
        result.converged = measureError <= options.relativeTolerance*std::max(current.measure,previous.measure) &&
          centerError <= options.relativeTolerance*size+roundoff ;
      }
      result.order = order ; result.evaluations += evaluations ; previous = current ;
      result.measure = checkedDouble(current.measure) ;
      if(result.measure == 0) throw std::overflow_error("geometry measure underflows its result type") ;
      for(unsigned d = 0; d < 3; ++d) result.centroid[d] = checkedDouble(static_cast<long double>(origin[d])+current.centroid[d]) ;
      if(result.converged) return result ;
    }
    if(result.order == 0) throw std::length_error("geometry measure budget cannot evaluate its first rule") ;
    return result ;
  }

  GeometryMeasure geometryMeasure(const GeometryMap &geometry, MeasureOptions options) {
    if(const auto *cone=dynamic_cast<const ConeMap *>(&geometry)) {
      // Integrate the linear cone direction exactly: volume weight 1/3, first moment 1/4.
      return detail::integrateGeometry(cone->base()->shape(),[&](const Point &q) {
        const auto base=cone->base()->evaluate(q) ;
        if(!base.hasJacobian) throw std::invalid_argument("cone measure requires surface derivatives") ;
        WidePoint a,b,height ; Point centroid ;
        for(unsigned d=0;d<3;++d) {
          a[d]=base.jacobian[0][d] ; b[d]=base.jacobian[1][d] ;
          height[d]=static_cast<long double>(cone->apex()[d])-base.position[d] ;
          centroid[d]=cone->apex()[d]-.75L*height[d] ;
        }
        const auto normal=cross(a,b) ;
        return detail::MeasureSample{centroid,(normal[0]*height[0]+normal[1]*height[1]+normal[2]*height[2])/3} ;
      },options) ;
    }
    if(geometry.shape()==Shape::polygon) {
      const auto corners=geometry.referenceCorners() ;
      const auto borrowed=std::shared_ptr<const GeometryMap>(&geometry,[](const GeometryMap *) {}) ;
      GeometryMeasure result ; result.converged=true ; result.estimatedMeasureError=result.estimatedCentroidError=0 ;
      const auto origin=geometry.evaluate({{.5,.5,0}}).position ;
      long double total=0,centroidError=0 ; std::array<long double,3> moment{} ;
      for(unsigned i=0;i<corners.size();++i) {
        auto placement=std::make_shared<CornerMap>(Shape::triangle,std::vector<Point>{{{.5,.5,0}},corners[i],corners[(i+1)%corners.size()]}) ;
        auto bounds=options ; bounds.maxEvaluations=options.maxEvaluations-result.evaluations ;
        const auto value=geometryMeasure(RestrictedMap(borrowed,placement),bounds) ;
        total+=value.measure ;
        for(unsigned d=0;d<3;++d) moment[d]+=value.measure*(static_cast<long double>(value.centroid[d])-origin[d]) ;
        centroidError+=value.measure*value.estimatedCentroidError+value.estimatedMeasureError*
          std::hypot(value.centroid[0]-origin[0],value.centroid[1]-origin[1],value.centroid[2]-origin[2]) ;
        result.converged=result.converged && value.converged ; result.order=std::max(result.order,value.order) ;
        result.evaluations+=value.evaluations ; result.estimatedMeasureError+=value.estimatedMeasureError ;
      }
      result.measure=checkedDouble(total) ;
      if(result.measure==0) throw std::overflow_error("polygon measure underflows its result type") ;
      for(unsigned d=0;d<3;++d) result.centroid[d]=checkedDouble(origin[d]+moment[d]/total) ;
      result.estimatedCentroidError=(centroidError+result.estimatedMeasureError*
        std::hypot(moment[0]/total,moment[1]/total,moment[2]/total))/total ;
      return result ;
    }
    const unsigned dimension = geometry.dimension() ;
    return detail::integrateGeometry(geometry.shape(),[&](const Point &q) {
      const auto value = geometry.evaluate(q) ;
      if(!value.hasJacobian) throw std::invalid_argument("geometry measure requires derivatives") ;
      std::array<WidePoint,3> columns{} ;
      for(unsigned d = 0; d < dimension; ++d) {
        for(unsigned component = 0; component < 3; ++component) {
          const double entry = value.jacobian[d][component] ;
          if(!std::isfinite(entry)) throw std::invalid_argument("geometry measure requires finite derivatives") ;
          columns[d][component] = entry ;
        }
      }
      long double jacobian = norm(columns[0]) ;
      if(dimension == 2) jacobian = norm(cross(columns[0],columns[1])) ;
      if(dimension == 3) {
        const auto product = cross(columns[1],columns[2]) ;
        jacobian = columns[0][0]*product[0]+columns[0][1]*product[1]+columns[0][2]*product[2] ;
      }
      return detail::MeasureSample{value.position,jacobian} ;
    },options) ;
  }

}
}
