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

#include <FVMAdapt/parametric/geometry_map.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Loci {
namespace FVMAdapt {
  namespace {

    void validateControls(const std::vector<Point> &controls) {
      for(const Point &point : controls)
        for(double value : point)
          if(!std::isfinite(value))
            throw std::invalid_argument("geometry controls must be finite") ;
    }

    void validatePoint(const Point &point) {
      for(double coordinate : point)
        if(!std::isfinite(coordinate))
          throw std::invalid_argument("geometry coordinates must be finite") ;
    }

    MapValue emptyValue() {
      return {{{0,0,0}}, {{{{0,0,0}},{{0,0,0}},{{0,0,0}}}}, true} ;
    }

    // The nearest point to zero in the control bounds keeps every offset finite.
    template<class Controls> Point controlOrigin(const Controls &controls) {
      Point origin = {{0,0,0}} ;
      if(controls.empty()) return origin ;
      for(unsigned d = 0; d < 3; ++d) {
        double low = controls[0][d], high = low ;
        for(const Point &point : controls) {
          if(point[d] < low) low = point[d] ;
          if(point[d] > high) high = point[d] ;
        }
        origin[d] = low > 0 ? low : (high < 0 ? high : 0) ;
      }
      return origin ;
    }

    // Partition of unity permits accumulation relative to a fixed origin.
    struct MapSum {
      const Point &origin ;
      Point position = {{0,0,0}} ;
      std::array<Point,3> derivative = {} ;

      void add(const Point &control, double weight, const Point &gradient) {
        for(unsigned component = 0; component < 3; ++component) {
          const double offset = control[component]-origin[component] ;
          position[component] += weight*offset ;
          for(unsigned direction = 0; direction < 3; ++direction)
            derivative[direction][component] += gradient[direction]*offset ;
        }
      }

      MapValue value() const {
        MapValue result = emptyValue() ;
        for(unsigned component = 0; component < 3; ++component) {
          result.position[component] = origin[component]+position[component] ;
          for(unsigned direction = 0; direction < 3; ++direction)
            result.jacobian[direction][component] = derivative[direction][component] ;
        }
        return result ;
      }
    } ;

    MapValue hexValue(const Point *corners, const Point &origin, const Point &point) {
      MapSum sum{origin} ;
      const auto &element = referenceElement(Shape::hex) ;
      for(unsigned i = 0; i < 8; ++i) {
        Point factors, signs, gradient ;
        for(unsigned d = 0; d < 3; ++d) {
          const bool upper = element.vertices[i][d] == 1 ;
          factors[d] = upper ? point[d] : 1-point[d] ; signs[d] = upper ? 1 : -1 ;
        }
        for(unsigned d = 0; d < 3; ++d) gradient[d] = signs[d]*factors[(d+1)%3]*factors[(d+2)%3] ;
        sum.add(corners[i],factors[0]*factors[1]*factors[2],gradient) ;
      }
      return sum.value() ;
    }

    std::vector<double> bernstein(unsigned degree, double coordinate) {
      std::vector<double> values(static_cast<std::size_t>(degree)+1,0) ;
      values[0] = 1 ;
      for(std::size_t level = 1; level <= degree; ++level) {
        for(std::size_t i = level; i > 0; --i)
          values[i] = (1-coordinate)*values[i]+coordinate*values[i-1] ;
        values[0] *= 1-coordinate ;
      }
      return values ;
    }

  }

  MapValue evaluateHexCorners(const std::array<Point,8> &corners, const Point &point) {
    validatePoint(point) ;
    return hexValue(corners.data(),controlOrigin(corners),point) ;
  }

  ReferenceCorners GeometryMap::referenceCorners() const {
    const auto &element=referenceElement(shape()) ;
    return {element.vertices.data(),element.vertexCount} ;
  }
  bool GeometryMap::containsReference(const Point &q,double tolerance) const { return contains(shape(),q,tolerance) ; }

  std::vector<Point> polygonReferenceCorners(unsigned vertices) {
    if(vertices<3) throw std::invalid_argument("polygon requires at least three vertices") ;
    std::vector<Point> result ; result.reserve(vertices) ;
    const double angle=2*std::acos(-1.)/vertices ;
    for(unsigned i=0;i<vertices;++i) result.push_back({{.5+.5*std::cos(i*angle),.5+.5*std::sin(i*angle),0}}) ;
    return result ;
  }
  bool containsReferencePolygon(unsigned vertices,const Point &q,double tolerance) {
    if(!std::isfinite(tolerance) || tolerance<0) throw std::invalid_argument("invalid polygon containment tolerance") ;
    for(auto x:q) if(!std::isfinite(x)) return false ;
    if(std::abs(q[2])>tolerance) return false ;
    const auto corners=polygonReferenceCorners(vertices) ;
    for(unsigned i=0;i<vertices;++i) {
      const auto &a=corners[i],&b=corners[(i+1)%vertices] ;
      if((b[0]-a[0])*(q[1]-a[1])-(b[1]-a[1])*(q[0]-a[0]) < -tolerance*std::hypot(b[0]-a[0],b[1]-a[1])) return false ;
    }
    return true ;
  }
  bool insideReferencePolygon(unsigned vertices,const Point &q) {
    for(auto x:q) if(!std::isfinite(x)) return false ;
    if(q[2]!=0) return false ;
    const auto corners=polygonReferenceCorners(vertices) ;
    for(unsigned i=0;i<vertices;++i) {
      const auto &a=corners[i],&b=corners[(i+1)%vertices] ;
      if((b[0]-a[0])*(q[1]-a[1])-(b[1]-a[1])*(q[0]-a[0])<=0) return false ;
    }
    return true ;
  }
  PolygonMap::PolygonMap(std::vector<Point> corners)
    : corners_(std::move(corners)),reference_(polygonReferenceCorners(corners_.size())) {
    if(corners_.size()<5) throw std::invalid_argument("polygon map requires at least five corners") ;
    validateControls(corners_) ;
  }
  bool PolygonMap::containsReference(const Point &q,double tolerance) const {
    return containsReferencePolygon(corners_.size(),q,tolerance) ;
  }
  MapValue PolygonMap::evaluate(const Point &q) const {
    validatePoint(q) ; const auto n=corners_.size() ;
    std::vector<std::array<long double,3>> areas(n) ; std::vector<long double> scales(n) ;
    for(unsigned e=0;e<n;++e) {
      const auto &a=reference_[e],&b=reference_[(e+1)%n] ;
      const long double dx=static_cast<long double>(b[0])-a[0],dy=static_cast<long double>(b[1])-a[1] ;
      scales[e]=dx*(.5L-a[1])-dy*(.5L-a[0]) ;
      areas[e]={{(dx*(static_cast<long double>(q[1])-a[1])-dy*(static_cast<long double>(q[0])-a[0]))/scales[e],
                  -dy/scales[e],dx/scales[e]}} ;
    }
    // Product weights avoid division by a vanishing edge distance at edges and vertices.
    std::array<long double,3> total{} ; std::array<std::array<long double,3>,3> moments{} ;
    for(unsigned i=0;i<n;++i) {
      const auto prev=(i+n-1)%n,next=(i+1)%n ;
      const long double ax=static_cast<long double>(reference_[i][0])-reference_[prev][0],
                        ay=static_cast<long double>(reference_[i][1])-reference_[prev][1],
                        bx=static_cast<long double>(reference_[next][0])-reference_[i][0],
                        by=static_cast<long double>(reference_[next][1])-reference_[i][1] ;
      std::array<long double,3> w{{(ax*by-ay*bx)/(scales[prev]*scales[i]),0,0}} ;
      for(unsigned e=0;e<n;++e) if(e!=i && e!=prev) {
        w[1]=w[1]*areas[e][0]+w[0]*areas[e][1] ; w[2]=w[2]*areas[e][0]+w[0]*areas[e][2] ; w[0]*=areas[e][0] ;
      }
      for(unsigned k=0;k<3;++k) {
        total[k]+=w[k] ;
        for(unsigned d=0;d<3;++d) moments[k][d]+=w[k]*(static_cast<long double>(corners_[i][d])-corners_[0][d]) ;
      }
    }
    if(total[0]==0 || !std::isfinite(total[0])) throw std::invalid_argument("polygon map has singular coordinates") ;
    auto result=emptyValue() ;
    for(unsigned d=0;d<3;++d) {
      const auto position=moments[0][d]/total[0] ; result.position[d]=corners_[0][d]+position ;
      for(unsigned k=0;k<2;++k) result.jacobian[k][d]=(moments[k+1][d]-position*total[k+1])/total[0] ;
    }
    return result ;
  }

  CornerMap::CornerMap(Shape shape, std::vector<Point> corners)
    : shape_(shape), corners_(std::move(corners)), origin_(controlOrigin(corners_)) {
    if(corners_.size() != referenceElement(shape_).vertexCount)
      throw std::invalid_argument("geometry corner count does not match the reference shape") ;
    validateControls(corners_) ;
  }

  MapValue CornerMap::evaluate(const Point &point) const {
    validatePoint(point) ;
    MapSum sum{origin_} ;
    const double x = point[0], y = point[1], z = point[2] ;
    switch(shape_) {
    case Shape::polygon: throw std::invalid_argument("polygon requires a PolygonMap") ;
    case Shape::point:
      break ;
    case Shape::segment:
      sum.add(corners_[0],1-x,{{-1,0,0}}) ;
      sum.add(corners_[1],x,{{1,0,0}}) ;
      break ;
    case Shape::triangle:
      sum.add(corners_[0],1-x-y,{{-1,-1,0}}) ;
      sum.add(corners_[1],x,{{1,0,0}}) ;
      sum.add(corners_[2],y,{{0,1,0}}) ;
      break ;
    case Shape::tet:
      sum.add(corners_[0],1-x-y-z,{{-1,-1,-1}}) ;
      sum.add(corners_[1],x,{{1,0,0}}) ;
      sum.add(corners_[2],y,{{0,1,0}}) ;
      sum.add(corners_[3],z,{{0,0,1}}) ;
      break ;
    case Shape::hex:
      return hexValue(corners_.data(),origin_,point) ;
    case Shape::quad: {
      const ReferenceElement &element = referenceElement(shape_) ;
      for(unsigned i = 0; i < element.vertexCount; ++i) {
        Point factors = {{1,1,1}}, signs = {{0,0,0}}, gradient = {{0,0,0}} ;
        for(unsigned d = 0; d < element.dimension; ++d) {
          const bool upper = element.vertices[i][d] == 1 ;
          factors[d] = upper ? point[d] : 1-point[d] ;
          signs[d] = upper ? 1 : -1 ;
        }
        for(unsigned d = 0; d < element.dimension; ++d)
          gradient[d] = signs[d]*factors[(d+1)%3]*factors[(d+2)%3] ;
        sum.add(corners_[i],factors[0]*factors[1]*factors[2],gradient) ;
      }
      break ;
    }
    case Shape::prism: {
      const double triangle[3] = {1-x-y,x,y} ;
      const Point gradients[3] = {{{-1,-1,0}},{{1,0,0}},{{0,1,0}}} ;
      for(unsigned layer = 0; layer < 2; ++layer) {
        const double height = layer ? z : 1-z, sign = layer ? 1 : -1 ;
        for(unsigned i = 0; i < 3; ++i)
          sum.add(corners_[3*layer+i],triangle[i]*height,
                  {{gradients[i][0]*height,gradients[i][1]*height,sign*triangle[i]}}) ;
      }
      break ;
    }
    case Shape::pyramid: {
      if(z == 1) {
        if(x != 0 || y != 0)
          throw std::invalid_argument("pyramid coordinates are singular at z = 1 outside the apex") ;
        MapValue apex = emptyValue() ;
        apex.position = corners_[4] ;
        apex.hasJacobian = false ;
        return apex ;
      }
      // Within the reference pyramid, collapsed square coordinates keep the basis bounded near the apex.
      const double height = 1-z, u = x/height, v = y/height ;
      const ReferenceElement &element = referenceElement(Shape::quad) ;
      for(unsigned i = 0; i < 4; ++i) {
        const double a = element.vertices[i][0] == 1 ? u : 1-u ;
        const double b = element.vertices[i][1] == 1 ? v : 1-v ;
        const double du = (element.vertices[i][0] == 1 ? 1 : -1)*b ;
        const double dv = (element.vertices[i][1] == 1 ? 1 : -1)*a ;
        sum.add(corners_[i],height*a*b,{{du,dv,-a*b+u*du+v*dv}}) ;
      }
      sum.add(corners_[4],z,{{0,0,1}}) ;
      break ;
    }
    }
    return sum.value() ;
  }

  BezierMap::BezierMap(Shape shape, std::array<unsigned,3> degree, std::vector<Point> controls)
    : shape_(shape), degree_(degree), controls_(std::move(controls)), origin_(controlOrigin(controls_)) {
    if(shape_ != Shape::segment && shape_ != Shape::quad && shape_ != Shape::hex)
      throw std::invalid_argument("Bezier maps require a segment, quad or hex") ;
    const unsigned dimension = referenceElement(shape_).dimension ;
    std::size_t count = 1 ;
    for(unsigned d = 0; d < 3; ++d) {
      if((d < dimension && degree_[d] == 0) || (d >= dimension && degree_[d] != 0))
        throw std::invalid_argument("Bezier degrees do not match the reference dimension") ;
      const std::size_t degree = degree_[d] ;
      if(degree == std::numeric_limits<std::size_t>::max() ||
         count > std::numeric_limits<std::size_t>::max()/(degree+1))
        throw std::invalid_argument("Bezier control count exceeds the addressable size") ;
      count *= degree+1 ;
    }
    if(controls_.size() != count)
      throw std::invalid_argument("Bezier control count does not match its degrees") ;
    validateControls(controls_) ;
  }

  MapValue BezierMap::evaluate(const Point &point) const {
    validatePoint(point) ;
    std::array<std::vector<double>,3> basis, derivative ;
    for(unsigned d = 0; d < 3; ++d) {
      basis[d] = bernstein(degree_[d],point[d]) ;
      derivative[d].assign(static_cast<std::size_t>(degree_[d])+1,0) ;
      if(degree_[d] == 0) continue ;
      const std::vector<double> lower = bernstein(degree_[d]-1,point[d]) ;
      for(std::size_t i = 0; i <= degree_[d]; ++i)
        derivative[d][i] = degree_[d]*((i > 0 ? lower[i-1] : 0)-(i < degree_[d] ? lower[i] : 0)) ;
    }
    MapSum sum{origin_} ;
    std::size_t index = 0 ;
    for(std::size_t k = 0; k <= degree_[2]; ++k)
      for(std::size_t j = 0; j <= degree_[1]; ++j)
        for(std::size_t i = 0; i <= degree_[0]; ++i, ++index)
          sum.add(controls_[index],basis[0][i]*basis[1][j]*basis[2][k],
                  {{derivative[0][i]*basis[1][j]*basis[2][k],
                    basis[0][i]*derivative[1][j]*basis[2][k],
                    basis[0][i]*basis[1][j]*derivative[2][k]}}) ;
    return sum.value() ;
  }

  RestrictedMap::RestrictedMap(std::shared_ptr<const GeometryMap> parent,
                               std::shared_ptr<const GeometryMap> referenceMap)
    : parent_(std::move(parent)), referenceMap_(std::move(referenceMap)) {
    if(!parent_ || !referenceMap_)
      throw std::invalid_argument("geometry restriction requires two maps") ;
    if(referenceMap_->dimension() > parent_->dimension())
      throw std::invalid_argument("geometry restriction dimension exceeds the parent dimension") ;
  }

  MapValue RestrictedMap::evaluate(const Point &point) const {
    const MapValue reference = referenceMap_->evaluate(point) ;
    const MapValue physical = parent_->evaluate(reference.position) ;
    MapValue value = emptyValue() ;
    value.position = physical.position ;
    value.hasJacobian = reference.hasJacobian && physical.hasJacobian ;
    if(value.hasJacobian)
      for(unsigned direction = 0; direction < 3; ++direction)
        for(unsigned component = 0; component < 3; ++component)
          for(unsigned inner = 0; inner < 3; ++inner)
            value.jacobian[direction][component] += physical.jacobian[inner][component]*reference.jacobian[direction][inner] ;
    return value ;
  }

  BoundarySurfaceMap::BoundarySurfaceMap(std::array<std::shared_ptr<const GeometryMap>,4> edges)
    : edges_(std::move(edges)) {
    for(unsigned e=0; e<4; ++e) {
      if(!edges_[e] || edges_[e]->shape()!=Shape::segment)
        throw std::invalid_argument("boundary surface requires four segment maps") ;
      corners_[e]=edges_[e]->evaluate({{0,0,0}}).position ;
      validatePoint(corners_[e]) ;
    }
    for(unsigned e=0; e<4; ++e) {
      const auto end=edges_[e]->evaluate({{1,0,0}}).position ; validatePoint(end) ;
      for(unsigned d=0; d<3; ++d)
        if(std::abs(end[d]-corners_[(e+1)%4][d]) > 64*std::numeric_limits<double>::epsilon()*
           std::max({1.,std::abs(end[d]),std::abs(corners_[(e+1)%4][d])}))
          throw std::invalid_argument("boundary surface edges do not meet") ;
    }
  }

  MapValue BoundarySurfaceMap::evaluate(const Point &point) const {
    validatePoint(point) ; const double u=point[0], v=point[1] ;
    const std::array<MapValue,4> edges{{edges_[0]->evaluate({{u,0,0}}),edges_[1]->evaluate({{v,0,0}}),
      edges_[2]->evaluate({{1-u,0,0}}),edges_[3]->evaluate({{1-v,0,0}})}} ;
    // Coons blending preserves each edge curve, subtracting the corner blend counted twice.
    MapValue result=emptyValue() ;
    for(const auto &edge:edges) result.hasJacobian=result.hasJacobian && edge.hasJacobian ;
    for(unsigned d=0; d<3; ++d) {
      const double a=edges[0].position[d]-corners_[0][d], b=edges[1].position[d]-corners_[0][d],
                   c=edges[2].position[d]-corners_[0][d], f=edges[3].position[d]-corners_[0][d] ;
      const double p=corners_[1][d]-corners_[0][d], q=corners_[2][d]-corners_[0][d], r=corners_[3][d]-corners_[0][d] ;
      const double blend=u*(1-v)*p+u*v*q+(1-u)*v*r ;
      result.position[d]=corners_[0][d]+(1-v)*a+u*b+v*c+(1-u)*f-blend ;
      if(result.hasJacobian) {
        result.jacobian[0][d]=(1-v)*edges[0].jacobian[0][d]+b-v*edges[2].jacobian[0][d]-f-((1-v)*p+v*(q-r)) ;
        result.jacobian[1][d]=-a+u*edges[1].jacobian[0][d]+c-(1-u)*edges[3].jacobian[0][d]-(u*(q-p)+(1-u)*r) ;
      }
    }
    return result ;
  }

  ConeMap::ConeMap(std::shared_ptr<const GeometryMap> base, Point apex)
    : base_(std::move(base)), apex_(apex) {
    if(!base_ || (base_->shape()!=Shape::triangle && base_->shape()!=Shape::quad))
      throw std::invalid_argument("cone base must be a triangle or quad") ;
    validatePoint(apex_) ;
  }

  MapValue ConeMap::evaluate(const Point &point) const {
    validatePoint(point) ; MapValue result=emptyValue() ;
    const double height=1-point[2] ;
    if(height==0) {
      if(point[0]!=0 || point[1]!=0)
        throw std::invalid_argument("cone coordinates are singular outside the apex") ;
      result.position=apex_ ; result.hasJacobian=false ; return result ;
    }
    const double u=point[0]/height, v=point[1]/height ;
    const auto base=base_->evaluate({{u,v,0}}) ; result.hasJacobian=base.hasJacobian ;
    for(unsigned d=0; d<3; ++d) {
      result.position[d]=base.position[d]+point[2]*(apex_[d]-base.position[d]) ;
      if(result.hasJacobian) {
        result.jacobian[0][d]=base.jacobian[0][d] ; result.jacobian[1][d]=base.jacobian[1][d] ;
        result.jacobian[2][d]=apex_[d]-base.position[d]+u*base.jacobian[0][d]+v*base.jacobian[1][d] ;
      }
    }
    return result ;
  }

}
}
