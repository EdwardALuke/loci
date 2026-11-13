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

#include <FVMAdapt/parametric/geometry_center.h>
#include <FVMAdapt/parametric/surface_map.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {
  namespace {

    using WidePoint = std::array<long double,3> ;

    // Bound local searches; failure leaves the proposed split unaccepted.
    constexpr unsigned maxFaceIterations=128, maxCellIterations=64, maxBacktracks=48 ;

    WidePoint difference(const Point &a, const Point &b) {
      return {{static_cast<long double>(a[0])-b[0],static_cast<long double>(a[1])-b[1],
               static_cast<long double>(a[2])-b[2]}} ;
    }

    WidePoint cross(const WidePoint &a, const WidePoint &b) {
      return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}} ;
    }

    long double accuracy(const Point &origin, long double scale) {
      long double magnitude = 0 ;
      for(double x: origin) magnitude = std::max(magnitude,std::abs(static_cast<long double>(x))) ;
      return 8*std::numeric_limits<double>::epsilon()*(1+magnitude/scale) ;
    }

    long double dot(const WidePoint &a, const WidePoint &b) {
      return a[0]*b[0]+a[1]*b[1]+a[2]*b[2] ;
    }

    WidePoint integratedTarget(const GeometryMap &geometry, const Point &origin, long double scale,
                               MeasureOptions options) {
      const auto measure = geometryMeasure(geometry,options) ;
      if(!measure.converged) throw std::length_error("center integration did not converge within its bounds") ;
      auto target = difference(measure.centroid,origin) ;
      for(auto &coordinate: target) coordinate /= scale ;
      return target ;
    }

    bool interior(Shape shape, const Point &q) {
      if(!(q[0] > 0 && q[1] > 0)) return false ;
      switch(shape) {
      case Shape::triangle: return q[0]+q[1] < 1 ;
      case Shape::quad: return q[0] < 1 && q[1] < 1 ;
      case Shape::tet: return q[2] > 0 && q[0]+q[1]+q[2] < 1 ;
      case Shape::hex: return q[2] > 0 && q[0] < 1 && q[1] < 1 && q[2] < 1 ;
      case Shape::prism: return q[2] > 0 && q[0]+q[1] < 1 && q[2] < 1 ;
      case Shape::pyramid: return q[2] > 0 && q[0]+q[2] < 1 && q[1]+q[2] < 1 ;
      default: return false ;
      }
    }

    MapValue checkedEvaluation(const GeometryMap &map, const Point &q) {
      auto value = map.evaluate(q) ;
      for(double x: value.position) if(!std::isfinite(x))
        throw std::invalid_argument("center evaluation requires finite geometry") ;
      return value ;
    }

    void scaledColumns(const MapValue &value, long double scale,
                       WidePoint &a, WidePoint &b) {
      if(!value.hasJacobian)
        throw std::invalid_argument("center evaluation requires a Jacobian") ;
      for(unsigned d = 0; d < 3; ++d) {
        if(!std::isfinite(value.jacobian[0][d]) || !std::isfinite(value.jacobian[1][d]))
          throw std::invalid_argument("center evaluation requires finite Jacobian entries") ;
        a[d] = value.jacobian[0][d]/scale ;
        b[d] = value.jacobian[1][d]/scale ;
      }
    }

    // Targets use (physical position - first corner)/scale. The face target also supplies a starting coordinate.
    WidePoint weightedFaceTarget(const GeometryMap &geometry, const std::vector<Point> &corners,
                                 long double scale, Point &initialCoordinates) {
      const auto referenceCorners=geometry.referenceCorners() ;
      WidePoint target{{0,0,0}} ;
      long double total = 0 ;
      WidePoint reference{{0,0,0}} ;
      for(unsigned e = 0; e < referenceCorners.size(); ++e) {
        const unsigned first = e, second = (e+1)%referenceCorners.size() ;
        auto delta = difference(corners[second],corners[first]) ;
        const long double weight = std::hypot(delta[0],delta[1],delta[2])/scale ;
        Point midpoint ;
        for(unsigned d = 0; d < 3; ++d)
          midpoint[d] = (referenceCorners[first][d]+referenceCorners[second][d])/2 ;
        delta = difference(checkedEvaluation(geometry,midpoint).position,corners[0]) ;
        for(unsigned d = 0; d < 3; ++d) {
          reference[d] += weight*midpoint[d] ;
          target[d] += weight*delta[d]/scale ;
        }
        total += weight ;
      }
      for(unsigned d = 0; d < 3; ++d) {
        initialCoordinates[d] = static_cast<double>(reference[d]/total) ;
        target[d] /= total ;
      }
      return target ;
    }

    // Minimize physical distance to the target with Gauss-Newton steps kept inside the face.
    Point projectFaceTarget(const GeometryMap &geometry, const Point &origin, long double scale,
                            const WidePoint &target, Point q, CenterMethod method) {
      const Shape shape = geometry.shape() ;
      const auto inside=[&](const Point &p) { return shape==Shape::polygon ? insideReferencePolygon(geometry.referenceCorners().size(),p):interior(shape,p) ; } ;
      if(!inside(q)) throw std::invalid_argument("face center has no interior starting point") ;
      const long double tolerance = 256*std::numeric_limits<double>::epsilon() ;
      for(unsigned iteration = 0; iteration < maxFaceIterations; ++iteration) {
        const auto value = checkedEvaluation(geometry,q) ;
        WidePoint a, b ;
        scaledColumns(value,scale,a,b) ;
        const long double aa = dot(a,a), ab = dot(a,b), bb = dot(b,b) ;
        const WidePoint normal{{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}} ;
        const long double determinant = dot(normal,normal) ;
        if(!(determinant > 0) || !std::isfinite(determinant))
          throw std::invalid_argument("face center has singular derivatives") ;
        if(method == CenterMethod::reference) return q ;
        auto residual = difference(value.position,origin) ;
        for(unsigned d = 0; d < 3; ++d) residual[d] = residual[d]/scale-target[d] ;
        const long double gradient0 = dot(residual,a), gradient1 = dot(residual,b) ;
        const long double step0 = (bb*gradient0-ab*gradient1)/determinant ;
        const long double step1 = (aa*gradient1-ab*gradient0)/determinant ;
        if(!std::isfinite(step0) || !std::isfinite(step1))
          throw std::invalid_argument("face center has singular derivatives") ;
        const long double dataAccuracy = accuracy(origin,scale) ;
        const long double referenceAccuracy = dataAccuracy*std::sqrt((aa+bb)/determinant) ;
        // Stop when the reference step is below the accuracy of the physical input.
        if(std::hypot(step0,step1) <= std::max(tolerance,referenceAccuracy)) return q ;
        const long double objective = dot(residual,residual) ;
        const long double descent = gradient0*step0+gradient1*step1 ;
        bool accepted = false ;
        for(unsigned backtrack = 0; backtrack < maxBacktracks; ++backtrack) {
          const long double fraction = std::ldexp(1.L,-static_cast<int>(backtrack)) ;
          const Point trial{{static_cast<double>(q[0]-fraction*step0),
                             static_cast<double>(q[1]-fraction*step1),0}} ;
          if(!inside(trial)) continue ;
          auto trialResidual = difference(checkedEvaluation(geometry,trial).position,origin) ;
          for(unsigned d = 0; d < 3; ++d) trialResidual[d] = trialResidual[d]/scale-target[d] ;
          const long double roundoff = 4*dataAccuracy*std::sqrt(objective) ;
          if(dot(trialResidual,trialResidual) <= objective-1e-4L*fraction*descent+roundoff) {
            q = trial ;
            accepted = true ;
            break ;
          }
        }
        if(!accepted) break ;
      }
      throw std::invalid_argument("face center projection did not converge in the face interior") ;
    }

    Point closestTriangle(const std::array<WidePoint,3> &p,const WidePoint &target) {
      WidePoint a,b,r ;
      for(unsigned d=0;d<3;++d) { a[d]=p[1][d]-p[0][d] ; b[d]=p[2][d]-p[0][d] ; r[d]=target[d]-p[0][d] ; }
      const long double aa=dot(a,a),ab=dot(a,b),bb=dot(b,b),det=dot(cross(a,b),cross(a,b)) ;
      if(!(det>0)) throw std::invalid_argument("surface piece projection has singular derivatives") ;
      const long double u=(bb*dot(a,r)-ab*dot(b,r))/det,v=(aa*dot(b,r)-ab*dot(a,r))/det ;
      if(u>=0 && v>=0 && u+v<=1) return {{double(u),double(v),0}} ;
      const std::array<Point,3> q{{{{0,0,0}},{{1,0,0}},{{0,1,0}}}} ;
      long double best=std::numeric_limits<long double>::infinity() ; Point result{} ;
      for(unsigned e=0;e<3;++e) {
        const unsigned next=(e+1)%3 ; WidePoint edge,residual ;
        for(unsigned d=0;d<3;++d) { edge[d]=p[next][d]-p[e][d] ; residual[d]=target[d]-p[e][d] ; }
        const long double t=std::max(0.L,std::min(1.L,dot(residual,edge)/dot(edge,edge))) ;
        for(unsigned d=0;d<3;++d) residual[d]-=t*edge[d] ;
        const auto value=dot(residual,residual) ;
        if(value<best) { best=value ; for(unsigned d=0;d<3;++d) result[d]=q[e][d]+t*(q[next][d]-q[e][d]) ; }
      }
      return result ;
    }

    // Piece boundaries are allowed here; an outer face boundary is rejected after piece selection.
    Point projectTriangle(const GeometryMap &map,const Point &origin,long double scale,const WidePoint &target) {
      Point q{{1./3,1./3,0}} ;
      for(unsigned iteration=0;iteration<maxFaceIterations;++iteration) {
        const auto value=checkedEvaluation(map,q) ; WidePoint a,b ; scaledColumns(value,scale,a,b) ;
        auto residual=difference(value.position,origin) ; std::array<WidePoint,3> linear ;
        for(unsigned d=0;d<3;++d) {
          residual[d]/=scale ; linear[0][d]=residual[d]-q[0]*a[d]-q[1]*b[d] ;
          linear[1][d]=linear[0][d]+a[d] ; linear[2][d]=linear[0][d]+b[d] ; residual[d]-=target[d] ;
        }
        const auto next=closestTriangle(linear,target) ;
        const long double step=std::hypot(static_cast<long double>(next[0])-q[0],static_cast<long double>(next[1])-q[1]) ;
        if(step<=256*std::numeric_limits<double>::epsilon()) return q ;
        const long double objective=dot(residual,residual),roundoff=16*accuracy(origin,scale)*std::sqrt(objective) ; bool accepted=false ;
        for(unsigned backtrack=0;backtrack<maxBacktracks;++backtrack) {
          const long double t=std::ldexp(1.L,-int(backtrack)) ; const Point trial{{double(q[0]+t*(next[0]-q[0])),double(q[1]+t*(next[1]-q[1])),0}} ;
          auto r=difference(checkedEvaluation(map,trial).position,origin) ; for(unsigned d=0;d<3;++d) r[d]=r[d]/scale-target[d] ;
          if(dot(r,r)<=objective+roundoff) { q=trial ; accepted=true ; break ; }
        }
        if(!accepted) break ;
      }
      throw std::invalid_argument("surface piece projection did not converge") ;
    }

    Point projectSurfacePieces(const GeometryMap &geometry,const Point &origin,long double scale,const WidePoint &target) {
      const auto borrowed=std::shared_ptr<const GeometryMap>(&geometry,[](const GeometryMap *) {}) ;
      const auto partition=smoothSurfacePieces(borrowed) ; long double best=std::numeric_limits<long double>::infinity() ; Point chart{} ;
      for(const auto &piece:partition.pieces) {
        const auto q=projectTriangle(*piece.geometry,origin,scale,target) ;
        auto residual=difference(checkedEvaluation(*piece.geometry,q).position,origin) ; for(unsigned d=0;d<3;++d) residual[d]=residual[d]/scale-target[d] ;
        const auto objective=dot(residual,residual) ;
        if(objective<best) { best=objective ; chart=piece.coordinates->position(q) ; }
      }
      const auto local=inverseSurfaceCoordinates(*partition.placement,chart) ;
      const bool inside=geometry.shape()==Shape::polygon ? insideReferencePolygon(geometry.referenceCorners().size(),local):interior(geometry.shape(),local) ;
      if(!inside) throw std::invalid_argument("surface center projects onto its outer boundary") ;
      return local ;
    }

    WidePoint weightedCellTarget(const GeometryMap &geometry, const std::vector<Point> &faceCenters,
                                 const std::vector<Point> &corners, long double scale) {
      const Shape shape = geometry.shape() ;
      const auto &element = referenceElement(shape) ;
      WidePoint target{{0,0,0}} ;
      long double total = 0 ;
      for(unsigned f = 0; f < element.boundaryCount; ++f) {
        const auto &face = element.boundaries[f] ;
        WidePoint area{{0,0,0}} ;
        for(unsigned i = 0; i < face.vertexCount; ++i) {
          auto a = difference(corners[face.vertices[i]],corners[face.vertices[0]]) ;
          auto b = difference(corners[face.vertices[(i+1)%face.vertexCount]],corners[face.vertices[0]]) ;
          for(unsigned d = 0; d < 3; ++d) { a[d] /= scale ; b[d] /= scale ; }
          const auto product = cross(a,b) ;
          for(unsigned d = 0; d < 3; ++d) area[d] += product[d]/2 ;
        }
        const long double weight = std::hypot(area[0],area[1],area[2]) ;
        const Point point = boundaryPoint(shape,f,faceCenters[f]) ;
        const auto delta = difference(checkedEvaluation(geometry,point).position,corners[0]) ;
        for(unsigned d = 0; d < 3; ++d) target[d] += weight*delta[d]/scale ;
        total += weight ;
      }
      if(!(total > 0)) throw std::invalid_argument("volume center requires nonzero boundary areas") ;
      for(auto &value: target) value /= total ;
      return target ;
    }

    // Newton inversion with decreasing physical residual and reference-interior trial points.
    Point invertCellTarget(const GeometryMap &geometry, const Point &origin, long double scale,
                           const WidePoint &target, Point q, CenterMethod method) {
      const Shape shape = geometry.shape() ;
      const long double dataAccuracy = accuracy(origin,scale) ;
      for(unsigned iteration = 0; iteration < maxCellIterations; ++iteration) {
        const auto value = checkedEvaluation(geometry,q) ;
        WidePoint a, b, c ;
        scaledColumns(value,scale,a,b) ;
        for(unsigned d = 0; d < 3; ++d) {
          if(!std::isfinite(value.jacobian[2][d]))
            throw std::invalid_argument("volume center requires finite derivatives") ;
          c[d] = value.jacobian[2][d]/scale ;
        }
        const auto bc = cross(b,c), ca = cross(c,a), ab = cross(a,b) ;
        const long double determinant = dot(a,bc) ;
        if(!(determinant > 0) || !std::isfinite(determinant))
          throw std::invalid_argument("volume center requires a positive Jacobian") ;
        if(method == CenterMethod::reference) return q ;
        auto residual = difference(value.position,origin) ;
        for(unsigned d = 0; d < 3; ++d) residual[d] = residual[d]/scale-target[d] ;
        const WidePoint step{{dot(residual,bc)/determinant,dot(residual,ca)/determinant,
                              dot(residual,ab)/determinant}} ;
        for(auto coordinate: step) if(!std::isfinite(coordinate))
          throw std::invalid_argument("volume center has singular derivatives") ;
        const long double objective = dot(residual,residual) ;
        if(std::sqrt(objective) <= dataAccuracy) return q ;
        bool accepted = false ;
        for(unsigned backtrack = 0; backtrack < maxBacktracks; ++backtrack) {
          const long double fraction = std::ldexp(1.L,-static_cast<int>(backtrack)) ;
          Point trial ;
          for(unsigned d = 0; d < 3; ++d) trial[d] = static_cast<double>(q[d]-fraction*step[d]) ;
          if(!interior(shape,trial)) continue ;
          auto trialResidual = difference(checkedEvaluation(geometry,trial).position,origin) ;
          for(unsigned d = 0; d < 3; ++d) trialResidual[d] = trialResidual[d]/scale-target[d] ;
          const long double roundoff = 4*dataAccuracy*std::sqrt(objective) ;
          if(dot(trialResidual,trialResidual) <= objective*(1-1e-4L*fraction)+roundoff) {
            q = trial ;
            accepted = true ;
            break ;
          }
        }
        if(!accepted) break ;
      }
      throw std::invalid_argument("volume center target has no resolved interior inverse") ;
    }

  }

  Point faceCenterCoordinates(const GeometryMap &geometry, CenterMethod method, MeasureOptions options) {
    const Shape shape = geometry.shape() ;
    if(geometry.dimension()!=2)
      throw std::invalid_argument("face center requires a surface") ;
    if(method != CenterMethod::reference && method != CenterMethod::weighted && method != CenterMethod::integrated)
      throw std::invalid_argument("unknown face center method") ;
    const auto referenceCorners=geometry.referenceCorners() ;
    Point q = shape == Shape::triangle ? Point{{1./3,1./3,0}} : Point{{.5,.5,0}} ;
    std::vector<Point> corners ;
    for(const auto &p:referenceCorners) corners.push_back(checkedEvaluation(geometry,p).position) ;
    long double scale = 0 ;
    for(const auto &p: corners) {
      const auto delta = difference(p,corners[0]) ;
      scale = std::max(scale,std::hypot(delta[0],delta[1],delta[2])) ;
    }
    if(!(scale > 0)) throw std::invalid_argument("face center requires distinct corners") ;
    WidePoint target{{0,0,0}} ;
    if(method == CenterMethod::weighted) target = weightedFaceTarget(geometry,corners,scale,q) ;
    if(method == CenterMethod::integrated) target = integratedTarget(geometry,corners[0],scale,options) ;
    if(hasSurfacePieces(geometry)) return method==CenterMethod::reference ? q:projectSurfacePieces(geometry,corners[0],scale,target) ;
    return projectFaceTarget(geometry,corners[0],scale,target,q,method) ;
  }

  Point cellCenterCoordinates(const GeometryMap &geometry, const std::vector<Point> &faceCenters,
                              CenterMethod method, MeasureOptions options) {
    const Shape shape = geometry.shape() ;
    const auto &element = referenceElement(shape) ;
    if(element.dimension != 3)
      throw std::invalid_argument("volume center requires a volume element") ;
    if(method != CenterMethod::reference && method != CenterMethod::weighted && method != CenterMethod::integrated)
      throw std::invalid_argument("unknown volume center method") ;
    if(method == CenterMethod::weighted) {
      if(faceCenters.size() != element.boundaryCount)
        throw std::invalid_argument("volume center requires one center per boundary") ;
      for(unsigned f = 0; f < element.boundaryCount; ++f)
        if(!contains(element.boundaries[f].shape,faceCenters[f]))
          throw std::invalid_argument("volume face center is outside its reference boundary") ;
    }
    std::vector<Point> corners ;
    Point q{{0,0,0}} ;
    for(unsigned v = 0; v < element.vertexCount; ++v) {
      corners.push_back(checkedEvaluation(geometry,element.vertices[v]).position) ;
      for(unsigned d = 0; d < 3; ++d) q[d] += element.vertices[v][d]/element.vertexCount ;
    }
    long double scale = 0 ;
    for(const auto &p: corners) {
      const auto delta = difference(p,corners[0]) ;
      scale = std::max(scale,std::hypot(delta[0],delta[1],delta[2])) ;
    }
    if(!(scale > 0)) throw std::invalid_argument("volume center requires distinct corners") ;
    WidePoint target{{0,0,0}} ;
    if(method == CenterMethod::weighted) target = weightedCellTarget(geometry,faceCenters,corners,scale) ;
    if(method == CenterMethod::integrated) target = integratedTarget(geometry,corners[0],scale,options) ;
    return invertCellTarget(geometry,corners[0],scale,target,q,method) ;
  }

}
}
