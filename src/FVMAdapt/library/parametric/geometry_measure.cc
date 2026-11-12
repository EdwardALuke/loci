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
#include <FVMAdapt/parametric/surface_map.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
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

    std::vector<GaussPoint> makeGaussRule(unsigned order) {
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

    const std::vector<GaussPoint> &gaussRule(unsigned order,std::vector<GaussPoint> &work) {
      // Common orders share immutable tables; other orders use local storage.
      static const auto tables=[]() {
        std::array<std::vector<GaussPoint>,8> result ;
        for(unsigned i=0;i<result.size();++i) result[i]=makeGaussRule(2*(i+1)) ;
        return result ;
      }() ;
      if(order<=16 && order%2==0) return tables[order/2-1] ;
      work=makeGaussRule(order) ; return work ;
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
        if(!std::isfinite(value.density) || !(value.density > 0))
          throw std::invalid_argument("geometry measure requires positive nonsingular sampled Jacobians") ;
        if(!hasOrigin) { origin = value.position ; hasOrigin = true ; }
        WidePoint offset ;
        for(unsigned d = 0; d < 3; ++d) offset[d] = static_cast<long double>(value.position[d])-origin[d] ;
        size = std::max(size,norm(offset)) ;
        weight *= value.density ; measure.add(weight) ;
        for(unsigned d = 0; d < 3; ++d) moment[d].add(weight*offset[d]) ;
      }
      Integral result ; result.measure = measure.total() ;
      if(!std::isfinite(result.measure) || !(result.measure > 0))
        throw std::overflow_error("geometry measure integral is unrepresentable") ;
      for(unsigned d = 0; d < 3; ++d) {
        result.moment[d] = moment[d].total() ;
        if(!std::isfinite(result.moment[d])) throw std::overflow_error("geometry moment is unrepresentable") ;
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

    detail::MeasureSample coneSample(const ConeMap &cone,const Point &q) {
      const auto base=cone.base()->evaluate(q) ;
      if(!base.hasJacobian) throw std::invalid_argument("cone measure requires surface derivatives") ;
      WidePoint a,b,height ; Point centroid ;
      for(unsigned d=0;d<3;++d) {
        a[d]=base.jacobian[0][d] ; b[d]=base.jacobian[1][d] ;
        height[d]=static_cast<long double>(cone.apex()[d])-base.position[d] ;
        centroid[d]=cone.apex()[d]-.75L*height[d] ;
      }
      const auto normal=cross(a,b) ;
      return {centroid,(normal[0]*height[0]+normal[1]*height[1]+normal[2]*height[2])/3} ;
    }

    detail::MeasureSample mapSample(const GeometryMap &geometry,const Point &q,unsigned dimension) {
      const auto value=geometry.evaluate(q) ;
      if(!value.hasJacobian) throw std::invalid_argument("geometry measure requires derivatives") ;
      std::array<WidePoint,3> columns{} ;
      for(unsigned d=0;d<dimension;++d) for(unsigned c=0;c<3;++c) {
        const double entry=value.jacobian[d][c] ;
        if(!std::isfinite(entry)) throw std::invalid_argument("geometry measure requires finite derivatives") ;
        columns[d][c]=entry ;
      }
      long double density=norm(columns[0]) ;
      if(dimension==2) density=norm(cross(columns[0],columns[1])) ;
      if(dimension==3) {
        const auto normal=cross(columns[1],columns[2]) ;
        density=columns[0][0]*normal[0]+columns[0][1]*normal[1]+columns[0][2]*normal[2] ;
      }
      return {value.position,density} ;
    }

    struct MeasureDomain {
      Shape shape=Shape::point ;
      std::vector<std::function<detail::MeasureSample(const Point&)>> samples ;
      bool prepare(std::shared_ptr<const GeometryMap> map) {
        if(const auto *cone=dynamic_cast<const ConeMap *>(map.get())) {
          if(hasSurfacePieces(*cone->base())) return false ;
          if(detail::radialConeSamples(*cone,samples)) shape=Shape::segment ;
          else { shape=cone->base()->shape() ; samples.push_back([map,cone](const Point &q) { return coneSample(*cone,q) ; }) ; }
        } else {
          shape=map->shape() ; const auto dimension=map->dimension() ;
          if(dimension>2 || hasSurfacePieces(*map)) return false ;
          samples.push_back([map,dimension](const Point &q) { return mapSample(*map,q,dimension) ; }) ;
        }
        return shape==Shape::segment || shape==Shape::triangle || shape==Shape::quad ;
      }
    } ;

    struct Estimate {
      long double measure=0,measureError=0,momentError=0 ;
      WidePoint moment{} ;
    } ;
    struct MeasureTotals {
      Sum measure,measureError,momentError ; std::array<Sum,3> moment ;
      unsigned unknownMeasure=0,unknownMoment=0 ;
      void add(const Estimate &value,int sign=1) {
        measure.add(sign*value.measure) ;
        for(unsigned d=0;d<3;++d) moment[d].add(sign*value.moment[d]) ;
        if(std::isfinite(value.measureError)) measureError.add(sign*value.measureError) ; else unknownMeasure+=sign ;
        if(std::isfinite(value.momentError)) momentError.add(sign*value.momentError) ; else unknownMoment+=sign ;
      }
      void finish(GeometryMeasure &result,const detail::MeasureScale &scale,MeasureOptions options) const {
        const auto volume=measure.total() ;
        result.measure=checkedDouble(volume) ;
        if(!(result.measure>0)) throw std::overflow_error("piecewise measure is not positive and representable") ;
        WidePoint offset ;
        for(unsigned d=0;d<3;++d) {
          offset[d]=moment[d].total()/volume ; result.centroid[d]=checkedDouble(scale.origin[d]+offset[d]) ;
        }
        result.estimatedMeasureError=unknownMeasure ? std::numeric_limits<double>::infinity():checkedDouble(std::max(0.L,measureError.total())) ;
        result.estimatedCentroidError=unknownMeasure || unknownMoment ? std::numeric_limits<double>::infinity():
          checkedDouble((std::max(0.L,momentError.total())+result.estimatedMeasureError*norm(offset))/volume) ;
        result.converged=result.estimatedMeasureError<=options.relativeTolerance*volume &&
          result.estimatedCentroidError<=options.relativeTolerance*scale.radius+
            8*std::numeric_limits<double>::epsilon()*std::max(scale.radius,scale.magnitude) ;
      }
    } ;

    struct MeasurePanel {
      Estimate value ; std::size_t source=0,serial=0 ; unsigned interval=0 ; long double priority=0 ; bool expanded=false ;
      std::array<Point,3> vertices{{{{0,0,0}},{{1,0,0}},{{0,1,0}}}} ;
      bool operator<(const MeasurePanel &other) const {
        return priority<other.priority || (priority==other.priority && serial>other.serial) ;
      }
    } ;

    bool bisectPanel(const MeasurePanel &parent,Shape shape,std::array<MeasurePanel,2> &children) {
      children={{parent,parent}} ; const auto &p=parent.vertices ;
      unsigned edge=0 ; long double longest=-1 ;
      const unsigned count=shape==Shape::triangle ? 3:shape==Shape::quad ? 2:1 ;
      for(unsigned e=0;e<count;++e) {
        const auto &a=p[shape==Shape::triangle ? e:0],&b=p[shape==Shape::triangle ? (e+1)%3:e+1] ;
        long double length=0 ; for(unsigned d=0;d<2;++d) { const long double x=static_cast<long double>(b[d])-a[d] ; length+=x*x ; }
        if(length>longest) { longest=length ; edge=e ; }
      }
      const auto &a=p[shape==Shape::triangle ? edge:0],&b=p[shape==Shape::triangle ? (edge+1)%3:edge+1] ;
      Point midpoint ; for(unsigned d=0;d<3;++d) midpoint[d]=(static_cast<long double>(a[d])+b[d])/2 ;
      if(midpoint==a || midpoint==b) return false ;
      if(shape==Shape::triangle) {
        const auto &c=p[(edge+2)%3] ; children[0].vertices={{a,midpoint,c}} ; children[1].vertices={{midpoint,b,c}} ;
      } else {
        const unsigned axis=edge+1 ; children[0].vertices[axis]=midpoint ; children[1].vertices[0]=midpoint ;
        if(shape==Shape::quad) for(unsigned d=0;d<3;++d) children[1].vertices[3-axis][d]+=midpoint[d]-a[d] ;
      }
      return true ;
    }

    Estimate measurePanel(const MeasurePanel &panel,const MeasureDomain &domain,const Point &origin,
                           const std::vector<GaussPoint> &lowRule,const std::vector<GaussPoint> &highRule) {
      const unsigned dimension=domain.shape==Shape::segment ? 1:2 ;
      const auto &a=panel.vertices[0],&b=panel.vertices[1],&c=panel.vertices[2] ;
      const long double jacobian=dimension==1 ? std::abs(static_cast<long double>(b[0])-a[0]):
        std::abs((static_cast<long double>(b[0])-a[0])*(static_cast<long double>(c[1])-a[1])-
                 (static_cast<long double>(b[1])-a[1])*(static_cast<long double>(c[0])-a[0])) ;
      if(!(jacobian>0)) throw std::length_error("integration subdivision has collapsed") ;
      const auto sample=[&](const Point &q) {
        Point p ; for(unsigned d=0;d<3;++d) p[d]=a[d]+static_cast<long double>(q[0])*(static_cast<long double>(b[d])-a[d])+
          (dimension==2 ? static_cast<long double>(q[1])*(static_cast<long double>(c[d])-a[d]):0) ;
        auto value=domain.samples[panel.interval](p) ; value.density*=jacobian ; return value ;
      } ;
      Point fixedOrigin=origin ; bool hasOrigin=true ; long double size=0 ;
      const auto low=integrate(domain.shape,sample,lowRule,dimension,fixedOrigin,hasOrigin,size) ;
      const auto high=integrate(domain.shape,sample,highRule,dimension,fixedOrigin,hasOrigin,size) ;
      WidePoint change ; for(unsigned d=0;d<3;++d) change[d]=high.moment[d]-low.moment[d] ;
      return {high.measure,std::abs(high.measure-low.measure),norm(change),high.moment} ;
    }

  }

  GeometryMeasure detail::measureGeometryPieces(
    const std::vector<std::shared_ptr<const GeometryMap>> &pieces,MeasureOptions options,const MeasureScale *scale) {
    if(pieces.empty()) throw std::invalid_argument("piecewise measure requires at least one piece") ;
    if(!std::isfinite(options.relativeTolerance) || options.relativeTolerance<=0 ||
       options.maxOrder<2 || options.maxOrder>256 || options.maxEvaluations==0)
      throw std::invalid_argument("invalid geometry measure bounds or tolerance") ;
    MeasureScale frame=scale ? *scale:MeasureScale{} ;
    for(double x:frame.origin) if(!std::isfinite(x)) throw std::invalid_argument("invalid integration origin") ;
    if(!std::isfinite(frame.radius) || !std::isfinite(frame.magnitude) || frame.radius<0 || frame.magnitude<0)
      throw std::invalid_argument("invalid integration scale") ;
    GeometryMeasure result ; MeasureTotals totals ; std::vector<Estimate> estimates ;
    estimates.reserve(pieces.size()) ;
    for(const auto &piece:pieces) {
      if(!piece || piece->dimension()==0) throw std::invalid_argument("invalid integration piece") ;
      if(result.evaluations>=options.maxEvaluations)
        throw std::length_error("piecewise measure exhausted its initial evaluation budget") ;
      auto bounds=options ; bounds.maxEvaluations-=result.evaluations ;
      const auto value=geometryMeasure(*piece,bounds) ;
      if(!scale && estimates.empty()) frame.origin=value.centroid ;
      Estimate estimate ; estimate.measure=value.measure ; estimate.measureError=value.estimatedMeasureError ;
      WidePoint offset ;
      for(unsigned d=0;d<3;++d) {
        offset[d]=static_cast<long double>(value.centroid[d])-frame.origin[d] ;
        estimate.moment[d]=value.measure*offset[d] ;
      }
      estimate.momentError=std::isfinite(value.estimatedMeasureError) && std::isfinite(value.estimatedCentroidError) ?
        value.measure*value.estimatedCentroidError+value.estimatedMeasureError*norm(offset):
        std::numeric_limits<long double>::infinity() ;
      if(!scale) frame.radius=std::max(frame.radius,norm(offset)) ;
      totals.add(estimate) ; estimates.push_back(estimate) ;
      result.order=std::max(result.order,value.order) ; result.evaluations+=value.evaluations ;
    }
    if(!scale) {
      frame.radius=std::max(frame.radius,std::pow(totals.measure.total(),1.L/pieces.front()->dimension())) ;
      frame.magnitude=frame.radius ;
      for(unsigned d=0;d<3;++d) frame.magnitude=std::max(frame.magnitude,
        std::abs(frame.origin[d]+totals.moment[d].total()/totals.measure.total())) ;
    }
    totals.finish(result,frame,options) ;
    if(result.converged || options.maxOrder<8) return result ;

    // Refine the largest absolute error contributions under the remaining shared budget.
    const long double measureTarget=options.relativeTolerance*totals.measure.total() ;
    const long double momentTarget=(options.relativeTolerance*frame.radius+
      8*std::numeric_limits<double>::epsilon()*std::max(frame.radius,frame.magnitude))*totals.measure.total() ;
    WidePoint offset ; for(unsigned d=0;d<3;++d) offset[d]=totals.moment[d].total()/totals.measure.total() ;
    const long double centerOffset=norm(offset) ;
    const auto priority=[&](const Estimate &value) {
      if(!std::isfinite(value.measureError) || !std::isfinite(value.momentError))
        return std::numeric_limits<long double>::infinity() ;
      const long double momentError=value.momentError+value.measureError*centerOffset ;
      return std::max(value.measureError/measureTarget,momentTarget>0 ? momentError/momentTarget:
        momentError>0 ? std::numeric_limits<long double>::infinity():0.L) ;
    } ;
    std::priority_queue<MeasurePanel> pending ; std::size_t serial=0 ;
    for(std::size_t i=0;i<estimates.size();++i) {
      MeasurePanel panel ; panel.value=estimates[i] ; panel.source=i ; panel.serial=serial++ ;
      panel.priority=priority(panel.value) ; pending.push(panel) ;
    }
    std::vector<MeasureDomain> domains(pieces.size()) ;
    std::vector<GaussPoint> lowWork,highWork ;
    const auto &lowRule=gaussRule(4,lowWork),&highRule=gaussRule(8,highWork) ;
    while(!pending.empty()) {
      auto parent=pending.top() ; pending.pop() ; auto &domain=domains[parent.source] ;
      if(!parent.expanded && !domain.prepare(pieces[parent.source])) continue ;
      const std::size_t cost=domain.shape==Shape::segment ? 12:80 ;
      std::array<MeasurePanel,2> children ;
      const unsigned count=parent.expanded ? 2:domain.samples.size() ;
      if(count==0 || count>children.size()) throw std::logic_error("invalid integration piece count") ;
      if(parent.expanded) { if(!bisectPanel(parent,domain.shape,children)) continue ; }
      else for(unsigned i=0;i<count;++i) { children[i]=parent ; children[i].interval=i ; }
      if(count*cost>options.maxEvaluations-result.evaluations) break ;
      for(unsigned i=0;i<count;++i) {
        auto &child=children[i] ; child.expanded=true ; child.serial=serial++ ;
        child.value=measurePanel(child,domain,frame.origin,lowRule,highRule) ; result.evaluations+=cost ;
        child.priority=priority(child.value) ;
      }
      totals.add(parent.value,-1) ;
      for(unsigned i=0;i<count;++i) { totals.add(children[i].value) ; pending.push(children[i]) ; }
      result.order=std::max(result.order,8U) ; totals.finish(result,frame,options) ;
      if(result.converged) break ;
    }
    return result ;
  }

  detail::MeasureIntegral detail::integrateGeometryRule(Shape shape,
    const std::function<MeasureSample(const Point&)> &sample, unsigned order, const Point &origin) {
    const unsigned dimension = referenceElement(shape).dimension ;
    if(dimension == 0 || order < 2 || order > 256)
      throw std::invalid_argument("invalid geometry integration rule") ;
    Point fixedOrigin = origin ; bool hasOrigin = true ; long double radius = 0 ;
    std::vector<GaussPoint> work ;
    const auto value = integrate(shape,sample,gaussRule(order,work),dimension,fixedOrigin,hasOrigin,radius) ;
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
    std::vector<GaussPoint> work ;
    for(unsigned order = 2; order <= options.maxOrder; order += 2) {
      std::size_t evaluations = order ;
      for(unsigned d = 1; d < dimension; ++d) evaluations *= order ;
      if(evaluations > options.maxEvaluations-result.evaluations) break ;
      const auto current = integrate(shape,sample,gaussRule(order,work),dimension,origin,hasOrigin,size) ;
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
    if(!std::isfinite(options.relativeTolerance) || options.relativeTolerance<=0 ||
       options.maxOrder<2 || options.maxOrder>256 || options.maxEvaluations==0)
      throw std::invalid_argument("invalid geometry measure bounds or tolerance") ;
    if(hasSurfacePieces(geometry)) {
      const auto borrowed=std::shared_ptr<const GeometryMap>(&geometry,[](const GeometryMap *) {}) ;
      if(geometry.dimension()==1) return detail::measureGeometryPieces(smoothCurvePieces(borrowed),options) ;
      if(geometry.dimension()==2) {
        std::vector<std::shared_ptr<const GeometryMap>> pieces ;
        for(const auto &piece:smoothSurfacePieces(borrowed).pieces) pieces.push_back(piece.geometry) ;
        return detail::measureGeometryPieces(pieces,options) ;
      }
    }
    if(const auto *cone=dynamic_cast<const ConeMap *>(&geometry)) {
      if(hasSurfacePieces(*cone->base())) {
        std::vector<std::shared_ptr<const GeometryMap>> pieces ;
        for(const auto &piece:smoothSurfacePieces(cone->base()).pieces) pieces.push_back(std::make_shared<ConeMap>(piece.geometry,cone->apex())) ;
        return detail::measureGeometryPieces(pieces,options) ;
      }
      // Integrate the linear cone direction exactly: volume weight 1/3, first moment 1/4.
      return detail::integrateGeometry(cone->base()->shape(),[&](const Point &q) { return coneSample(*cone,q) ; },options) ;
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
    const unsigned dimension=geometry.dimension() ;
    return detail::integrateGeometry(geometry.shape(),[&](const Point &q) { return mapSample(geometry,q,dimension) ; },options) ;
  }

}
}
