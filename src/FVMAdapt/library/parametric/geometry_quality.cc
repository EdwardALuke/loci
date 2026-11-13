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
#include "geometry_validity.h"

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
      if(!finite(value.position)) { quality.invalidEvaluation=true ; recordProblem(quality,point,true) ; return ; }
      if(!value.hasJacobian) { recordProblem(quality,point,false) ; return ; }
      for(unsigned d = 0; d < dimension; ++d)
        if(!finite(value.jacobian[d])) { quality.invalidEvaluation=true ; recordProblem(quality,point,true) ; return ; }
      std::array<WidePoint,3> tangents{} ;
      std::array<long double,3> norms{}, lengths{} ;
      for(unsigned d = 0; d < dimension; ++d) {
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

  namespace {
    using detail::Exact ;
    using detail::ExactPoint ;
    using detail::BernsteinCoefficients ;
    using detail::jacobianDifference ;
    using detail::jacobianDeterminant ;
    using detail::jacobianIndex ;

    ExactPoint exactPoint(const Point &point) { return {{Exact(point[0]),Exact(point[1]),Exact(point[2])}} ; }
    ExactPoint exactCross(const ExactPoint &a,const ExactPoint &b) {
      return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}} ;
    }
    Exact exactSquaredNorm(const ExactPoint &a) { return a[0]*a[0]+a[1]*a[1]+a[2]*a[2] ; }
    std::vector<ExactPoint> exactCorners(const CornerMap &map) {
      std::vector<ExactPoint> result ; for(const auto &p:map.corners()) result.push_back(exactPoint(p)) ; return result ;
    }

    // Convert endpoint/midpoint samples to quadratic Bernstein coefficients along each axis.
    template<class Function> BernsteinCoefficients coefficients(Function value) {
      BernsteinCoefficients result ;
      for(unsigned k=0;k<3;++k) for(unsigned j=0;j<3;++j) for(unsigned i=0;i<3;++i)
        result[jacobianIndex(i,j,k)]=value(Exact(i)/2,Exact(j)/2,Exact(k)/2) ;
      for(unsigned step:{1u,3u,9u}) for(unsigned i=0;i<27;++i) if(i/step%3==0)
        result[i+step]=2*result[i+step]-(result[i]+result[i+2*step])/2 ;
      return result ;
    }

    GeometryValidity certify(BernsteinCoefficients values,unsigned dimension,JacobianCheckLimits limits,std::size_t &remaining) {
      struct Box {
        BernsteinCoefficients values ;
        unsigned depth ;
        /// Bits 2*d and 2*d+1 mark sides on the original domain's lower and upper boundaries.
        unsigned boundarySides ;
      } ;
      std::vector<Box> pending{{std::move(values),0,(1u<<(2*dimension))-1}} ;
      while(!pending.empty()) {
        if(remaining==0) return GeometryValidity::unresolved ;
        --remaining ;
        auto box=std::move(pending.back()) ;
        pending.pop_back() ;
        bool nonnegative=true ;
        for(const auto &v:box.values)
          if(v<0) { nonnegative=false ; break ; }

        // Corner coefficients are exact values. Only original-domain boundaries may vanish.
        for(unsigned k=0;k<(dimension==3 ? 3u:1u);k+=2)
          for(unsigned j=0;j<(dimension>=2 ? 3u:1u);j+=2)
            for(unsigned i:{0u,2u}) {
              const auto &v=box.values[jacobianIndex(i,j,k)] ;
              if(v<0) return GeometryValidity::invalid ;
              bool exterior=false ;
              const unsigned q[3]={i,j,k} ;
              for(unsigned d=0;d<dimension;++d)
                exterior=exterior || (box.boundarySides&(1u<<(2*d+q[d]/2))) ;
              if(v==0 && !exterior) return GeometryValidity::invalid ;
            }
        if(nonnegative) {
          // Nonnegative coefficients bound the polynomial below by zero. Require a positive
          // contribution on each interior part, including sides introduced by proof subdivision.
          const unsigned count=dimension==3 ? 27:dimension==2 ? 9:3 ;
          for(unsigned part=0;part<count;++part) {
            unsigned code=part ;
            std::array<unsigned,3> side{} ;
            bool exterior=false ;
            for(unsigned d=0;d<dimension;++d) {
              // 0 spans the axis; 1 and 2 select its lower and upper ends.
              side[d]=code%3 ;
              code/=3 ;
              if(side[d]) exterior=exterior || (box.boundarySides&(1u<<(2*d+side[d]-1))) ;
            }
            if(exterior) continue ;
            bool positive=false ;
            for(unsigned n=0;n<27 && !positive;++n) {
              unsigned index=n ;
              bool member=true ;
              for(unsigned d=0;d<dimension;++d) {
                const unsigned q=index%3 ;
                index/=3 ;
                if(side[d] && q!=2*(side[d]-1)) member=false ;
              }
              if(member && box.values[n]>0) positive=true ;
            }
            if(!positive) return GeometryValidity::invalid ;
          }
          continue ;
        }
        if(box.depth>=limits.maxDepth) return GeometryValidity::unresolved ;
        const unsigned axis=box.depth%dimension ;
        auto children=detail::bisectJacobianCoefficients(box.values,axis) ;
        pending.push_back({std::move(children[1]),box.depth+1,box.boundarySides&~(1u<<(2*axis))}) ;
        pending.push_back({std::move(children[0]),box.depth+1,box.boundarySides&~(1u<<(2*axis+1))}) ;
      }
      return GeometryValidity::valid ;
    }

    struct SurfaceValue { ExactPoint point,u,v ; } ;
    SurfaceValue surfaceValue(const std::vector<ExactPoint> &p,const Exact &u,const Exact &v) {
      SurfaceValue result ; const bool quad=p.size()==4 ;
      for(unsigned d=0;d<3;++d) {
        const Exact a=p[1][d]-p[0][d],b=p[quad ? 3:2][d]-p[0][d] ;
        const Exact twist=quad ? Exact(p[2][d]-p[1][d]-p[3][d]+p[0][d]):Exact(0) ;
        result.point[d]=p[0][d]+u*a+v*b+u*v*twist ;
        result.u[d]=a+v*twist ; result.v[d]=b+u*twist ;
      }
      return result ;
    }
    GeometryValidity cornerValidity(const CornerMap &map,JacobianCheckLimits limits,std::size_t &remaining) {
      const auto shape=map.shape() ;
      if(shape==Shape::hex) {
        if(positiveHexJacobian(map.corners())) return GeometryValidity::valid ;
        return certify(detail::hexJacobianCoefficients<Exact>(map.corners()),3,limits,remaining) ;
      }
      const auto p=exactCorners(map) ;
      switch(shape) {
      case Shape::segment: return exactSquaredNorm(jacobianDifference(p[1],p[0]))>0 ? GeometryValidity::valid:GeometryValidity::invalid ;
      case Shape::triangle: return exactSquaredNorm(exactCross(jacobianDifference(p[1],p[0]),jacobianDifference(p[2],p[0])))>0
          ? GeometryValidity::valid:GeometryValidity::invalid ;
      case Shape::tet: return jacobianDeterminant(jacobianDifference(p[1],p[0]),jacobianDifference(p[2],p[0]),jacobianDifference(p[3],p[0]))>0
          ? GeometryValidity::valid:GeometryValidity::invalid ;
      case Shape::quad:
        return certify(coefficients([&](const Exact &u,const Exact &v,const Exact &) {
          const auto s=surfaceValue(p,u,v) ; return exactSquaredNorm(exactCross(s.u,s.v)) ;
        }),2,limits,remaining) ;
      case Shape::pyramid: {
        const std::vector<ExactPoint> base(p.begin(),p.begin()+4) ;
        return certify(coefficients([&](const Exact &u,const Exact &v,const Exact &) {
          const auto s=surfaceValue(base,u,v) ; return jacobianDeterminant(s.u,s.v,jacobianDifference(p[4],s.point)) ;
        }),2,limits,remaining) ;
      }
      case Shape::prism: {
        const auto a=jacobianDifference(p[3],p[0]),b=jacobianDifference(p[1],p[0]),c=jacobianDifference(p[2],p[0]) ;
        const auto db=jacobianDifference(jacobianDifference(p[4],p[3]),b),dc=jacobianDifference(jacobianDifference(p[5],p[3]),c) ;
        return certify(coefficients([&](const Exact &u,const Exact &q,const Exact &z) {
          const Exact v=(1-u)*q ; ExactPoint x,y,w ;
          for(unsigned d=0;d<3;++d) { x[d]=b[d]+z*db[d] ; y[d]=c[d]+z*dc[d] ; w[d]=a[d]+u*db[d]+v*dc[d] ; }
          return jacobianDeterminant(x,y,w) ;
        }),3,limits,remaining) ;
      }
      default: return GeometryValidity::unresolved ;
      }
    }

    bool inside(Shape shape,const ExactPoint &p) {
      const auto &x=p[0],&y=p[1],&z=p[2] ;
      if(x<0 || y<0 || z<0) return false ;
      switch(shape) {
      case Shape::segment: return x<=1 && y==0 && z==0 ;
      case Shape::triangle: return x+y<=1 && z==0 ;
      case Shape::quad: return x<=1 && y<=1 && z==0 ;
      case Shape::tet: return x+y+z<=1 ;
      case Shape::prism: return x+y<=1 && z<=1 ;
      case Shape::hex: return x<=1 && y<=1 && z<=1 ;
      case Shape::pyramid: return x+z<=1 && y+z<=1 && z<=1 ;
      default: return false ;
      }
    }
    bool contained(const CornerMap &placement,Shape parent) {
      for(const auto &p:placement.corners()) if(!inside(parent,exactPoint(p))) return false ;
      return true ;
    }
    ExactPoint affineValue(const detail::ExactMapView &map,const ExactPoint &q) {
      auto p=map.origin ;
      for(unsigned d=0;d<3;++d) for(unsigned c=0;c<3;++c) p[d]+=map.columns[c][d]*q[c] ;
      return p ;
    }
    using detail::affineTrace ;
    bool bilinearSurfaceCorners(const std::array<std::shared_ptr<const GeometryMap>,4> &curves,
                                const std::array<Point,4> &blend,std::vector<ExactPoint> &p,std::size_t &remaining) {
      std::array<std::array<ExactPoint,2>,4> edges ;
      for(unsigned edge=0;edge<4;++edge) {
        edges[edge]={{exactPoint({{0,0,0}}),exactPoint({{1,0,0}})}} ;
        if(!affineTrace(*curves[edge],edges[edge],remaining)) return false ;
      }
      // Affine edges make the Coons expression bilinear, including any endpoint mismatch.
      // Derive its coefficients from that expression rather than fitting evaluated positions.
      for(unsigned edge=0;edge<4;++edge) {
        ExactPoint corner ;
        for(unsigned d=0;d<3;++d)
          corner[d]=edges[edge][0][d]+edges[(edge+3)%4][1][d]-Exact(blend[edge][d]) ;
        p.push_back(std::move(corner)) ;
      }
      return true ;
    }
    bool surfaceCorners(const GeometryMap &geometry,std::vector<ExactPoint> &p,std::size_t &remaining) {
      const GeometryMap *map=&geometry ;
      while(remaining) {
        --remaining ;
        if(const auto *corners=dynamic_cast<const CornerMap *>(map)) {
          if(map->shape()!=Shape::triangle && map->shape()!=Shape::quad) return false ;
          p=exactCorners(*corners) ; return true ;
        }
        if(const auto *surface=dynamic_cast<const BoundarySurfaceMap *>(map))
          return bilinearSurfaceCorners(surface->edges(),surface->corners(),p,remaining) ;
        const auto view=detail::exactMapView(*map) ;
        if(view.kind==detail::ExactMapView::Kind::delegate) { map=view.delegate ; continue ; }
        if(view.kind==detail::ExactMapView::Kind::coons)
          return bilinearSurfaceCorners(*view.edges,*view.corners,p,remaining) ;
        if(view.kind!=detail::ExactMapView::Kind::affine || map->dimension()!=2) return false ;
        // A frozen affine sector extends over the parameter plane; its determinant is constant.
        const auto &element=referenceElement(map->shape()==Shape::quad ? Shape::quad:Shape::triangle) ;
        for(unsigned v=0;v<element.vertexCount;++v) p.push_back(affineValue(view,exactPoint(element.vertices[v]))) ;
        return true ;
      }
      return false ;
    }
    int surfaceOrientation(const std::vector<ExactPoint> &p,JacobianCheckLimits limits,std::size_t &remaining) {
      auto values=coefficients([&](const Exact &u,const Exact &v,const Exact &) {
        const auto s=surfaceValue(p,u,v) ; return Exact(s.u[0]*s.v[1]-s.u[1]*s.v[0]) ;
      }) ;
      if(certify(values,2,limits,remaining)==GeometryValidity::valid) return 1 ;
      for(auto &v:values) v=-v ;
      return certify(std::move(values),2,limits,remaining)==GeometryValidity::valid ? -1:0 ;
    }
    GeometryValidity coneValidity(const ConeMap &cone,JacobianCheckLimits limits,std::size_t &remaining) {
      const GeometryMap *base=cone.base().get() ; int orientation=1 ; bool restricted=false ; std::size_t maps=limits.maxBoxes ;
      while(true) {
        if(maps==0) return GeometryValidity::unresolved ;
        --maps ;
        const auto view=detail::exactMapView(*base) ;
        if(view.kind==detail::ExactMapView::Kind::delegate) { base=view.delegate ; continue ; }
        const auto *map=dynamic_cast<const RestrictedMap *>(base) ; if(!map) break ;
        std::vector<ExactPoint> placement ;
        if(map->parent()->dimension()!=2 || !surfaceCorners(*map->referenceMap(),placement,maps)) return GeometryValidity::unresolved ;
        if(detail::exactMapView(*map->parent()).kind!=detail::ExactMapView::Kind::affine)
          for(const auto &q:placement) if(!inside(map->parent()->shape(),q)) return GeometryValidity::unresolved ;
        const int sign=surfaceOrientation(placement,limits,remaining) ; if(sign==0) return GeometryValidity::unresolved ;
        orientation*=sign ; base=map->parent().get() ; restricted=true ;
      }
      std::vector<ExactPoint> p ;
      if(!surfaceCorners(*base,p,maps)) return GeometryValidity::unresolved ;
      const auto apex=exactPoint(cone.apex()) ;
      const auto result=certify(coefficients([&](const Exact &u,const Exact &v,const Exact &) {
        const auto s=surfaceValue(p,u,v) ; return Exact(Exact(orientation)*jacobianDeterminant(s.u,s.v,jacobianDifference(apex,s.point))) ;
      }),2,limits,remaining) ;
      // Failure on the whole parent surface does not establish failure on a restriction.
      return restricted && result==GeometryValidity::invalid ? GeometryValidity::unresolved:result ;
    }
  }

  GeometryValidity geometryValidity(const GeometryMap &geometry,JacobianCheckLimits limits) {
    if(limits.maxBoxes==0) throw std::invalid_argument("geometry validity requires a positive box limit") ;
    if(geometry.dimension()==0) throw std::invalid_argument("geometry validity requires a positive-dimensional shape") ;
    const GeometryMap *map=&geometry ; std::size_t remaining=limits.maxBoxes,maps=limits.maxBoxes ; bool restricted=false ;
    while(maps) {
      --maps ; GeometryValidity result=GeometryValidity::unresolved ;
      if(const auto *corners=dynamic_cast<const CornerMap *>(map)) result=cornerValidity(*corners,limits,remaining) ;
      else if(const auto *cone=dynamic_cast<const ConeMap *>(map)) result=coneValidity(*cone,limits,remaining) ;
      else if(const auto *child=dynamic_cast<const RestrictedMap *>(map)) {
        const auto *placement=dynamic_cast<const CornerMap *>(child->referenceMap().get()) ;
        if(!placement || child->dimension()!=child->parent()->dimension() || !contained(*placement,child->parent()->shape()) ||
           cornerValidity(*placement,limits,remaining)!=GeometryValidity::valid) return GeometryValidity::unresolved ;
        map=child->parent().get() ; restricted=true ; continue ;
      }
      return restricted && result==GeometryValidity::invalid ? GeometryValidity::unresolved:result ;
    }
    return GeometryValidity::unresolved ;
  }

  bool validGeometry(const GeometryMap &geometry,const GeometryQualityOptions &options,JacobianCheckLimits limits) {
    if(limits.maxBoxes==0) throw std::invalid_argument("geometry validity requires a positive box limit") ;
    const auto quality=geometryQuality(geometry,options) ;
    return quality.status==GeometryQualityStatus::sampled || (!quality.invalidEvaluation && geometryValidity(geometry,limits)==GeometryValidity::valid) ;
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
      catch(const std::exception &) { quality.invalidEvaluation=true ; recordProblem(quality,point,false) ; }
    }
    return quality ;
  }

}
}
