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

    template<class Controls> void validateControls(const Controls &controls) {
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

    const Point &cubicEdgeControl(const CubicEdgeUse &use, unsigned index) {
      return use.edge->controls()[use.reversed ? 3-index:index] ;
    }

    template<std::size_t N> void validateCubicEdges(const std::array<CubicEdgeUse,N> &edges) {
      for(const auto &use:edges) {
        if(!use.edge || use.edge->shape()!=Shape::segment ||
           use.edge->degree()!=std::array<unsigned,3>{{3,0,0}})
          throw std::invalid_argument("cubic geometry requires cubic segment edges") ;
        validateControls(use.edge->controls()) ;
      }
    }

    template<std::size_t N> void validateCubicFaceEdges(const std::array<CubicEdgeUse,N> &edges) {
      validateCubicEdges(edges) ;
      for(std::size_t e=0; e<N; ++e)
        if(cubicEdgeControl(edges[e],3)!=cubicEdgeControl(edges[(e+1)%N],0))
          throw std::invalid_argument("cubic face edges require identical stored endpoints") ;
    }

    template<std::size_t N> void validateCubicCellEdges(Shape shape, const std::array<CubicEdgeUse,N> &edges) {
      validateCubicEdges(edges) ;
      const auto &element=referenceElement(shape) ;
      std::array<const Point *,8> corners{} ;
      for(unsigned e=0; e<N; ++e)
        for(unsigned endpoint=0; endpoint<2; ++endpoint) {
          const unsigned vertex=element.edges[e][endpoint] ;
          const Point &point=cubicEdgeControl(edges[e],3*endpoint) ;
          if(corners[vertex] && *corners[vertex]!=point)
            throw std::invalid_argument("cubic cell edges require identical stored endpoints") ;
          corners[vertex]=&point ;
        }
    }

    template<std::size_t N> Point triangleInterior(const std::array<CubicEdgeUse,N> &edges,
                                                  const ReferenceElement &element, const Boundary &face) {
      std::array<CubicEdgeUse,3> boundary ;
      for(unsigned i=0; i<3; ++i) {
        const int a=face.vertices[i], b=face.vertices[(i+1)%3] ;
        for(unsigned e=0; e<N; ++e) {
          const auto ends=element.edges[e] ;
          if((ends[0]==a && ends[1]==b) || (ends[0]==b && ends[1]==a)) {
            boundary[i]=edges[e] ;
            if(ends[0]!=a) boundary[i].reversed=!boundary[i].reversed ;
            break ;
          }
        }
      }
      return cubicTriangleFromEdges(boundary).controls()[9] ;
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

    template<std::size_t N, std::size_t F>
    Point cubicCellOrigin(const std::array<CubicEdgeUse,N> &edges, const std::array<Point,F> &faceControls) {
      std::array<Point,4*N+F> controls ;
      for(unsigned e=0; e<N; ++e)
        for(unsigned i=0; i<4; ++i) controls[4*e+i]=cubicEdgeControl(edges[e],i) ;
      std::copy(faceControls.begin(),faceControls.end(),controls.begin()+4*N) ;
      return controlOrigin(controls) ;
    }

    // Partition of unity permits accumulation relative to a fixed origin.
    template<bool withJacobian> struct MapSum {
      const Point &origin ;
      Point position = {{0,0,0}} ;
      std::array<Point,3> derivative = {} ;

      void add(const Point &control, double weight, const Point &gradient) {
        for(unsigned component = 0; component < 3; ++component) {
          const double offset = control[component]-origin[component] ;
          position[component] += weight*offset ;
          if constexpr(withJacobian)
            for(unsigned direction = 0; direction < 3; ++direction)
              derivative[direction][component] += gradient[direction]*offset ;
        }
      }

      auto value() const {
        Point result ;
        for(unsigned component = 0; component < 3; ++component)
          result[component] = origin[component]+position[component] ;
        if constexpr(withJacobian) return MapValue{result,derivative,true} ;
        else return result ;
      }
    } ;

    template<bool withJacobian> auto hexValue(const Point *corners, const Point &origin, const Point &point) {
      MapSum<withJacobian> sum{origin} ;
      const auto &element = referenceElement(Shape::hex) ;
      for(unsigned i = 0; i < 8; ++i) {
        Point factors, signs, gradient{} ;
        for(unsigned d = 0; d < 3; ++d) {
          const bool upper = element.vertices[i][d] == 1 ;
          factors[d] = upper ? point[d] : 1-point[d] ; signs[d] = upper ? 1 : -1 ;
        }
        if constexpr(withJacobian)
          for(unsigned d = 0; d < 3; ++d) gradient[d] = signs[d]*factors[(d+1)%3]*factors[(d+2)%3] ;
        sum.add(corners[i],factors[0]*factors[1]*factors[2],gradient) ;
      }
      return sum.value() ;
    }

    template<bool withJacobian> auto cubicHexValue(const std::array<CubicEdgeUse,12> &edges,
                                                  const Point &origin, const Point &point) {
      validatePoint(point) ;
      std::array<std::array<double,4>,3> basis, derivative{} ;
      for(unsigned d=0; d<3; ++d) {
        const double t=point[d], s=1-t, tt=t*t, ss=s*s ;
        basis[d]={{ss*s,3*ss*t,3*s*tt,tt*t}} ;
        if constexpr(withJacobian) derivative[d]={{-3*ss,3*ss-6*s*t,6*s*t-3*tt,3*tt}} ;
      }
      MapSum<withJacobian> sum{origin} ;
      const auto &element=referenceElement(Shape::hex) ;
      // The six Coons faces reduce the volume blend to three edge blends minus twice the corner blend.
      for(unsigned e=0; e<12; ++e) {
        const Point &start=element.vertices[element.edges[e][0]], &end=element.vertices[element.edges[e][1]] ;
        unsigned axis=0 ;
        while(axis<2 && start[axis]==end[axis]) ++axis ;
        const unsigned a=(axis+1)%3, b=(axis+2)%3 ;
        const double wa=start[a] ? point[a]:1-point[a], wb=start[b] ? point[b]:1-point[b] ;
        for(unsigned i=0; i<4; ++i) {
          Point gradient{} ;
          if constexpr(withJacobian) {
            gradient[axis]=derivative[axis][i]*wa*wb ;
            gradient[a]=basis[axis][i]*(start[a] ? 1:-1)*wb ;
            gradient[b]=basis[axis][i]*wa*(start[b] ? 1:-1) ;
          }
          sum.add(cubicEdgeControl(edges[e],start[axis] ? 3-i:i),basis[axis][i]*wa*wb,gradient) ;
        }
      }
      for(unsigned i=0; i<8; ++i) {
        Point weight, gradient{} ;
        for(unsigned d=0; d<3; ++d) weight[d]=element.vertices[i][d] ? point[d]:1-point[d] ;
        if constexpr(withJacobian)
          for(unsigned d=0; d<3; ++d)
            gradient[d]=-2*(element.vertices[i][d] ? 1:-1)*weight[(d+1)%3]*weight[(d+2)%3] ;
        // The first eight native edges start at corners 0 through 7.
        sum.add(cubicEdgeControl(edges[i],0),-2*weight[0]*weight[1]*weight[2],gradient) ;
      }
      return sum.value() ;
    }

    template<bool withJacobian> auto cornerValue(Shape shape, const Point *corners, const Point &origin, const Point &point) {
      validatePoint(point) ;
      MapSum<withJacobian> sum{origin} ;
      const double x = point[0], y = point[1], z = point[2] ;
      switch(shape) {
      case Shape::polygon: throw std::invalid_argument("polygon requires a PolygonMap") ;
      case Shape::point:
        break ;
      case Shape::segment:
        sum.add(corners[0],1-x,{{-1,0,0}}) ;
        sum.add(corners[1],x,{{1,0,0}}) ;
        break ;
      case Shape::triangle:
        sum.add(corners[0],1-x-y,{{-1,-1,0}}) ;
        sum.add(corners[1],x,{{1,0,0}}) ;
        sum.add(corners[2],y,{{0,1,0}}) ;
        break ;
      case Shape::tet:
        sum.add(corners[0],1-x-y-z,{{-1,-1,-1}}) ;
        sum.add(corners[1],x,{{1,0,0}}) ;
        sum.add(corners[2],y,{{0,1,0}}) ;
        sum.add(corners[3],z,{{0,0,1}}) ;
        break ;
      case Shape::hex:
        return hexValue<withJacobian>(corners,origin,point) ;
      case Shape::quad: {
        const ReferenceElement &element = referenceElement(shape) ;
        for(unsigned i = 0; i < element.vertexCount; ++i) {
          Point factors = {{1,1,1}}, signs = {{0,0,0}}, gradient = {{0,0,0}} ;
          for(unsigned d = 0; d < element.dimension; ++d) {
            const bool upper = element.vertices[i][d] == 1 ;
            factors[d] = upper ? point[d] : 1-point[d] ;
            signs[d] = upper ? 1 : -1 ;
          }
          if constexpr(withJacobian)
            for(unsigned d = 0; d < element.dimension; ++d)
              gradient[d] = signs[d]*factors[(d+1)%3]*factors[(d+2)%3] ;
          sum.add(corners[i],factors[0]*factors[1]*factors[2],gradient) ;
        }
        break ;
      }
      case Shape::prism: {
        const double triangle[3] = {1-x-y,x,y} ;
        const Point gradients[3] = {{{-1,-1,0}},{{1,0,0}},{{0,1,0}}} ;
        for(unsigned layer = 0; layer < 2; ++layer) {
          const double height = layer ? z : 1-z, sign = layer ? 1 : -1 ;
          for(unsigned i = 0; i < 3; ++i) {
            Point gradient{} ;
            if constexpr(withJacobian) gradient = {{gradients[i][0]*height,gradients[i][1]*height,sign*triangle[i]}} ;
            sum.add(corners[3*layer+i],triangle[i]*height,gradient) ;
          }
        }
        break ;
      }
      case Shape::pyramid: {
        if(z == 1) {
          if(x != 0 || y != 0)
            throw std::invalid_argument("pyramid coordinates are singular at z = 1 outside the apex") ;
          if constexpr(withJacobian) return MapValue{corners[4],{},false} ;
          else return corners[4] ;
        }
        // Within the reference pyramid, collapsed square coordinates keep the basis bounded near the apex.
        const double height = 1-z, u = x/height, v = y/height ;
        const ReferenceElement &element = referenceElement(Shape::quad) ;
        for(unsigned i = 0; i < 4; ++i) {
          const double a = element.vertices[i][0] == 1 ? u : 1-u ;
          const double b = element.vertices[i][1] == 1 ? v : 1-v ;
          Point gradient{} ;
          if constexpr(withJacobian) {
            const double du = (element.vertices[i][0] == 1 ? 1 : -1)*b ;
            const double dv = (element.vertices[i][1] == 1 ? 1 : -1)*a ;
            gradient = {{du,dv,-a*b+u*du+v*dv}} ;
          }
          sum.add(corners[i],height*a*b,gradient) ;
        }
        sum.add(corners[4],z,{{0,0,1}}) ;
        break ;
      }
      }
      return sum.value() ;
    }

    // Degree elevation builds the Bernstein basis in zero-initialized storage.
    void bernstein(unsigned degree, double coordinate, double *values) {
      values[0] = 1 ;
      for(std::size_t level = 1; level <= degree; ++level) {
        for(std::size_t i = level; i > 0; --i)
          values[i] = (1-coordinate)*values[i]+coordinate*values[i-1] ;
        values[0] *= 1-coordinate ;
      }
    }

    std::vector<double> bernstein(unsigned degree, double coordinate) {
      std::vector<double> values(static_cast<std::size_t>(degree)+1) ;
      bernstein(degree,coordinate,values.data()) ;
      return values ;
    }

    template<bool withJacobian>
    auto cubicSegmentValue(const std::vector<Point> &controls, const Point &origin, const Point &point) {
      validatePoint(point) ;
      std::array<double,4> basis{}, derivative{} ;
      bernstein(3,point[0],basis.data()) ;
      if constexpr(withJacobian) {
        std::array<double,3> lower{} ;
        bernstein(2,point[0],lower.data()) ;
        for(unsigned i=0; i<4; ++i)
          derivative[i]=3*((i>0 ? lower[i-1]:0)-(i<3 ? lower[i]:0)) ;
      }
      MapSum<withJacobian> sum{origin} ;
      for(unsigned i=0; i<4; ++i) sum.add(controls[i],basis[i],{{derivative[i],0,0}}) ;
      return sum.value() ;
    }

    template<bool withJacobian>
    void cubicTriangleBasis(const Point &point, std::array<double,10> &weights, std::array<Point,10> &gradients) {
      const double u = point[0], v = point[1], w = 1-u-v ;
      const double uu = u*u, vv = v*v, ww = w*w ;
      weights={{ww*w,uu*u,vv*v,3*ww*u,3*w*uu,3*uu*v,3*u*vv,3*w*vv,3*ww*v,6*w*u*v}} ;
      if constexpr(withJacobian) {
        const double wu = w*u, wv = w*v, uv = u*v ;
        gradients={{
          {{-3*ww,-3*ww,0}}, {{3*uu,0,0}}, {{0,3*vv,0}},
          {{3*ww-6*wu,-6*wu,0}}, {{6*wu-3*uu,-3*uu,0}},
          {{6*uv,3*uu,0}}, {{3*vv,6*uv,0}},
          {{-3*vv,6*wv-3*vv,0}}, {{-6*wv,3*ww-6*wv,0}},
          {{6*v*(w-u),6*u*(w-v),0}}
        }} ;
      }
    }

    template<bool withJacobian> auto cubicTriangleValue(const std::array<Point,10> &controls,
                                                       const Point &origin, const Point &point) {
      validatePoint(point) ;
      std::array<double,10> weights ; std::array<Point,10> gradients{} ;
      cubicTriangleBasis<withJacobian>(point,weights,gradients) ;
      MapSum<withJacobian> sum{origin} ;
      for(unsigned i=0; i<controls.size(); ++i) sum.add(controls[i],weights[i],gradients[i]) ;
      return sum.value() ;
    }

    template<bool withJacobian>
    auto cubicTetValue(const std::array<CubicEdgeUse,6> &edges, const std::array<Point,4> &faceControls,
                       const Point &origin, const Point &point) {
      validatePoint(point) ;
      const std::array<double,4> weights{{1-point[0]-point[1]-point[2],point[0],point[1],point[2]}} ;
      const std::array<Point,4> gradients{{{{-1,-1,-1}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}} ;
      const auto &element=referenceElement(Shape::tet) ;
      MapSum<withJacobian> sum{origin} ;
      // A degree-three tetrahedron has vertex, edge and face controls, but no volume-interior control.
      for(unsigned i=0; i<4; ++i) {
        const double w=weights[i] ; Point gradient{} ;
        if constexpr(withJacobian)
          for(unsigned d=0; d<3; ++d) gradient[d]=3*w*w*gradients[i][d] ;
        const Point &corner=cubicEdgeControl(edges[i],i==3 ? 3:0) ;
        sum.add(corner,w*w*w,gradient) ;
      }
      for(unsigned e=0; e<6; ++e)
        for(unsigned endpoint=0; endpoint<2; ++endpoint) {
          const unsigned a=element.edges[e][endpoint], b=element.edges[e][1-endpoint] ;
          const double u=weights[a], v=weights[b] ; Point gradient{} ;
          if constexpr(withJacobian)
            for(unsigned d=0; d<3; ++d) gradient[d]=3*u*(2*v*gradients[a][d]+u*gradients[b][d]) ;
          sum.add(cubicEdgeControl(edges[e],1+endpoint),3*u*u*v,gradient) ;
        }
      for(unsigned f=0; f<4; ++f) {
        const auto &vertices=element.boundaries[f].vertices ;
        const unsigned a=vertices[0], b=vertices[1], c=vertices[2] ;
        const double u=weights[a], v=weights[b], w=weights[c] ; Point gradient{} ;
        if constexpr(withJacobian)
          for(unsigned d=0; d<3; ++d)
            gradient[d]=6*(v*w*gradients[a][d]+u*w*gradients[b][d]+u*v*gradients[c][d]) ;
        sum.add(faceControls[f],6*u*v*w,gradient) ;
      }
      return sum.value() ;
    }

    template<bool withJacobian>
    auto cubicPrismValue(const std::array<CubicEdgeUse,9> &edges, const std::array<Point,2> &faceControls,
                         const Point &origin, const Point &point) {
      validatePoint(point) ;
      std::array<double,10> weights ; std::array<Point,10> gradients{} ;
      cubicTriangleBasis<withJacobian>(point,weights,gradients) ;
      MapSum<withJacobian> sum{origin} ;
      const double z=point[2], s=1-z ;
      // Blend the two triangular ends, then add each axial edge's departure from its chord.
      for(unsigned layer=0; layer<2; ++layer) {
        const double height=layer ? z:s, sign=layer ? 1:-1 ;
        auto add=[&](const Point &control, unsigned i) {
          Point gradient{} ;
          if constexpr(withJacobian) gradient={{height*gradients[i][0],height*gradients[i][1],sign*weights[i]}} ;
          sum.add(control,height*weights[i],gradient) ;
        } ;
        for(unsigned e=0; e<3; ++e) {
          add(cubicEdgeControl(edges[3*layer+e],0),e) ;
          add(cubicEdgeControl(edges[3*layer+e],1),3+2*e) ;
          add(cubicEdgeControl(edges[3*layer+e],2),4+2*e) ;
        }
        add(faceControls[layer],9) ;
      }
      const std::array<double,3> triangle{{1-point[0]-point[1],point[0],point[1]}} ;
      const std::array<Point,3> triangleGradient{{{{-1,-1,0}},{{1,0,0}},{{0,1,0}}}} ;
      // Subtracting the linear endpoints in the basis avoids separate large-coordinate sums.
      const std::array<double,4> deviation{{s*s*s-s,3*s*s*z,3*s*z*z,z*z*z-z}} ;
      const std::array<double,4> derivative{{1-3*s*s,3*s*s-6*s*z,6*s*z-3*z*z,3*z*z-1}} ;
      for(unsigned e=0; e<3; ++e)
        for(unsigned i=0; i<4; ++i) {
          Point gradient{} ;
          if constexpr(withJacobian)
            gradient={{triangleGradient[e][0]*deviation[i],triangleGradient[e][1]*deviation[i],triangle[e]*derivative[i]}} ;
          sum.add(cubicEdgeControl(edges[6+e],i),triangle[e]*deviation[i],gradient) ;
        }
      return sum.value() ;
    }

    template<bool withJacobian>
    auto cubicPyramidValue(const std::array<CubicEdgeUse,8> &edges, const std::array<Point,4> &faceControls,
                           const Point &origin, const Point &point) {
      validatePoint(point) ;
      const double x=point[0], y=point[1], z=point[2], s=1-z ;
      const Point &apex=cubicEdgeControl(edges[4],3) ;
      if(s==0) {
        if(x!=0 || y!=0) throw std::invalid_argument("pyramid coordinates are singular at z = 1 outside the apex") ;
        if constexpr(withJacobian) return MapValue{apex,{},false} ;
        else return apex ;
      }
      // u,v cover the unit square as the pyramid section shrinks toward its apex.
      const double u=x/s, v=y/s, ss=s*s, zz=z*z ;
      MapSum<withJacobian> sum{origin} ;
      // Coons sections blend the triangular sides. Cancel height factors before
      // differentiating collapsed coordinates to keep derivatives bounded near the apex.
      const auto add=[&](const Point &control, double scale, double overHeight,
                         double heightDerivative, double weight, double du, double dv) {
        Point gradient{} ;
        if constexpr(withJacobian)
          gradient={{overHeight*du,overHeight*dv,heightDerivative*weight+overHeight*(u*du+v*dv)}} ;
        sum.add(control,scale*weight,gradient) ;
      } ;
      // Base corners and the controls on the four edges leading to the apex.
      const auto &quad=referenceElement(Shape::quad) ;
      for(unsigned i=0; i<4; ++i) {
        const double a=quad.vertices[i][0] ? u:1-u, b=quad.vertices[i][1] ? v:1-v ;
        const double da=quad.vertices[i][0] ? 1:-1, db=quad.vertices[i][1] ? 1:-1 ;
        const double n=a*b, f=a*a+b*b-1, g=a+b-1 ;
        add(cubicEdgeControl(edges[i],0),ss*s,ss,-3*ss,n*f,da*b*(f+2*a*a),db*a*(f+2*b*b)) ;
        add(cubicEdgeControl(edges[4+i],1),3*ss*z,3*s*z,3*ss-6*s*z,n*g,da*b*(g+a),db*a*(g+b)) ;
        add(cubicEdgeControl(edges[4+i],2),3*s*zz,3*zz,6*s*z-3*zz,n,da*b,db*a) ;
      }
      // Base-edge controls and the interior controls of the triangular sides.
      const double t[4]={u,v,1-u,1-v}, w[4]={1-v,u,v,1-u} ;
      const double tu[4]={1,0,-1,0}, tv[4]={0,1,0,-1}, wu[4]={0,1,0,-1}, wv[4]={-1,0,1,0} ;
      for(unsigned e=0; e<4; ++e) {
        const double a=1-t[e], b=t[e] ;
        const auto blend=[&](const Point &control,double scale,double overHeight,double heightDerivative,
                            double weight,double derivative) {
          add(control,scale,overHeight,heightDerivative,w[e]*weight,
              wu[e]*weight+w[e]*derivative*tu[e],wv[e]*weight+w[e]*derivative*tv[e]) ;
        } ;
        blend(cubicEdgeControl(edges[e],1),3*ss*s,3*ss,-9*ss,a*a*b,a*a-2*a*b) ;
        blend(cubicEdgeControl(edges[e],2),3*ss*s,3*ss,-9*ss,a*b*b,2*a*b-b*b) ;
        blend(faceControls[e],6*ss*z,6*s*z,6*ss-12*s*z,a*b,a-b) ;
      }
      // The common apex contributes once to the completed blend.
      sum.add(apex,zz*z,{{0,0,3*zz}}) ;
      return sum.value() ;
    }

  }

  MapValue evaluateHexCorners(const std::array<Point,8> &corners, const Point &point) {
    validatePoint(point) ;
    return hexValue<true>(corners.data(),controlOrigin(corners),point) ;
  }
  Point hexCornerPosition(const std::array<Point,8> &corners, const Point &point) {
    validatePoint(point) ;
    return hexValue<false>(corners.data(),controlOrigin(corners),point) ;
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
    return cornerValue<true>(shape_,corners_.data(),origin_,point) ;
  }
  Point CornerMap::position(const Point &point) const {
    return cornerValue<false>(shape_,corners_.data(),origin_,point) ;
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
    if(shape_==Shape::segment && degree_[0]==3)
      return cubicSegmentValue<true>(controls_,origin_,point) ;
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
    MapSum<true> sum{origin_} ;
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

  Point BezierMap::position(const Point &point) const {
    if(shape_==Shape::segment && degree_[0]==3)
      return cubicSegmentValue<false>(controls_,origin_,point) ;
    return evaluate(point).position ;
  }

  CubicTriangleMap::CubicTriangleMap(std::array<Point,10> controls)
    : controls_(std::move(controls)), origin_(controlOrigin(controls_)) {
    validateControls(controls_) ;
  }

  MapValue CubicTriangleMap::evaluate(const Point &point) const {
    return cubicTriangleValue<true>(controls_,origin_,point) ;
  }
  Point CubicTriangleMap::position(const Point &point) const {
    return cubicTriangleValue<false>(controls_,origin_,point) ;
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
  Point RestrictedMap::position(const Point &point) const {
    return parent_->position(referenceMap_->position(point)) ;
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
    return detail::evaluateBoundarySurface(edges_,corners_,point) ;
  }

  CubicTriangleMap cubicTriangleFromEdges(const std::array<CubicEdgeUse,3> &edges) {
    validateCubicFaceEdges(edges) ;
    std::array<Point,10> controls ;
    for(unsigned e=0; e<3; ++e) {
      controls[e]=cubicEdgeControl(edges[e],0) ;
      controls[3+2*e]=cubicEdgeControl(edges[e],1) ;
      controls[4+2*e]=cubicEdgeControl(edges[e],2) ;
    }
    // C=1.5*mean(edge controls)-0.5*mean(vertices) is symmetric and reproduces quadratics.
    for(unsigned d=0; d<3; ++d) {
      const long double origin=controls[0][d] ;
      long double edgeSum=0, vertexSum=0 ;
      for(unsigned i=3; i<9; ++i) edgeSum+=static_cast<long double>(controls[i][d])-origin ;
      for(unsigned i=0; i<3; ++i) vertexSum+=static_cast<long double>(controls[i][d])-origin ;
      const long double interior=origin+edgeSum/4-vertexSum/6 ;
      if(!std::isfinite(interior) || std::abs(interior)>std::numeric_limits<double>::max())
        throw std::invalid_argument("cubic triangle interior control must be finite") ;
      controls[9][d]=static_cast<double>(interior) ;
    }
    return CubicTriangleMap(controls) ;
  }

  BoundarySurfaceMap cubicQuadFromEdges(const std::array<CubicEdgeUse,4> &edges) {
    validateCubicFaceEdges(edges) ;
    std::array<std::shared_ptr<const GeometryMap>,4> curves ;
    for(unsigned e=0; e<4; ++e) {
      curves[e]=edges[e].edge ;
      if(edges[e].reversed)
        curves[e]=std::make_shared<RestrictedMap>(edges[e].edge,
          std::make_shared<CornerMap>(Shape::segment,std::vector<Point>{{{1,0,0}},{{0,0,0}}})) ;
    }
    return BoundarySurfaceMap(std::move(curves)) ;
  }

  CubicTetMap::CubicTetMap(std::array<CubicEdgeUse,6> edges) : edges_(std::move(edges)) {
    validateCubicCellEdges(Shape::tet,edges_) ;
    const auto &element=referenceElement(Shape::tet) ;
    for(unsigned f=0; f<4; ++f) faceControls_[f]=triangleInterior(edges_,element,element.boundaries[f]) ;
    origin_=cubicCellOrigin(edges_,faceControls_) ;
  }

  MapValue CubicTetMap::evaluate(const Point &point) const {
    return cubicTetValue<true>(edges_,faceControls_,origin_,point) ;
  }
  Point CubicTetMap::position(const Point &point) const {
    return cubicTetValue<false>(edges_,faceControls_,origin_,point) ;
  }

  CubicPrismMap::CubicPrismMap(std::array<CubicEdgeUse,9> edges) : edges_(std::move(edges)) {
    validateCubicCellEdges(Shape::prism,edges_) ;
    const auto &element=referenceElement(Shape::prism) ;
    for(unsigned f=0; f<2; ++f) faceControls_[f]=triangleInterior(edges_,element,element.boundaries[f]) ;
    origin_=cubicCellOrigin(edges_,faceControls_) ;
  }

  MapValue CubicPrismMap::evaluate(const Point &point) const {
    return cubicPrismValue<true>(edges_,faceControls_,origin_,point) ;
  }
  Point CubicPrismMap::position(const Point &point) const {
    return cubicPrismValue<false>(edges_,faceControls_,origin_,point) ;
  }

  CubicPyramidMap::CubicPyramidMap(std::array<CubicEdgeUse,8> edges) : edges_(std::move(edges)) {
    validateCubicCellEdges(Shape::pyramid,edges_) ;
    const auto &element=referenceElement(Shape::pyramid) ;
    for(unsigned f=0; f<4; ++f) faceControls_[f]=triangleInterior(edges_,element,element.boundaries[f+1]) ;
    origin_=cubicCellOrigin(edges_,faceControls_) ;
  }

  MapValue CubicPyramidMap::evaluate(const Point &point) const {
    return cubicPyramidValue<true>(edges_,faceControls_,origin_,point) ;
  }
  Point CubicPyramidMap::position(const Point &point) const {
    return cubicPyramidValue<false>(edges_,faceControls_,origin_,point) ;
  }

  CubicHexMap::CubicHexMap(std::array<CubicEdgeUse,12> edges) : edges_(std::move(edges)) {
    validateCubicCellEdges(Shape::hex,edges_) ;
    origin_=cubicCellOrigin(edges_,std::array<Point,0>{}) ;
  }

  MapValue CubicHexMap::evaluate(const Point &point) const {
    return cubicHexValue<true>(edges_,origin_,point) ;
  }
  Point CubicHexMap::position(const Point &point) const {
    return cubicHexValue<false>(edges_,origin_,point) ;
  }

  MapValue detail::evaluateBoundarySurface(const std::array<std::shared_ptr<const GeometryMap>,4> &curves,
                                          const std::array<Point,4> &corners,const Point &point) {
    validatePoint(point) ; const double u=point[0], v=point[1] ;
    const std::array<MapValue,4> edges{{curves[0]->evaluate({{u,0,0}}),curves[1]->evaluate({{v,0,0}}),
      curves[2]->evaluate({{1-u,0,0}}),curves[3]->evaluate({{1-v,0,0}})}} ;
    // Coons blending preserves each edge curve, subtracting the corner blend counted twice.
    MapValue result=emptyValue() ;
    for(const auto &edge:edges) result.hasJacobian=result.hasJacobian && edge.hasJacobian ;
    for(unsigned d=0; d<3; ++d) {
      const double a=edges[0].position[d]-corners[0][d], b=edges[1].position[d]-corners[0][d],
                   c=edges[2].position[d]-corners[0][d], f=edges[3].position[d]-corners[0][d] ;
      const double p=corners[1][d]-corners[0][d], q=corners[2][d]-corners[0][d], r=corners[3][d]-corners[0][d] ;
      const double blend=u*(1-v)*p+u*v*q+(1-u)*v*r ;
      result.position[d]=corners[0][d]+(1-v)*a+u*b+v*c+(1-u)*f-blend ;
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
