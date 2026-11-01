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

#include <FVMAdapt/parametric/planar_split.h>
#include <FVMAdapt/parametric/reference_coordinate.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {
  namespace {
    using WidePoint = std::array<long double,3> ;

    WidePoint difference(const Point &a, const Point &b) {
      return {{static_cast<long double>(a[0])-b[0],static_cast<long double>(a[1])-b[1],static_cast<long double>(a[2])-b[2]}} ;
    }

    WidePoint cross(const WidePoint &a, const WidePoint &b) {
      return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}} ;
    }

    long double dot(const WidePoint &a, const WidePoint &b) {
      return a[0]*b[0]+a[1]*b[1]+a[2]*b[2] ;
    }

    long double length(const WidePoint &a) {
      return std::hypot(a[0],a[1],a[2]) ;
    }

    void validatePolygon(const std::vector<Point> &points) {
      if(points.size() != 3 && points.size() != 4)
        throw std::invalid_argument("planar subdivision requires a triangle or quad") ;
      long double scale = 0 ;
      for(const auto &point: points) {
        for(double value: point) if(!std::isfinite(value))
          throw std::invalid_argument("planar subdivision requires finite corners") ;
        scale = std::max(scale,length(difference(point,points[0]))) ;
      }
      if(!(scale > 0)) throw std::invalid_argument("planar subdivision requires distinct corners") ;
      std::vector<WidePoint> local ;
      for(const auto &point: points) {
        auto value = difference(point,points[0]) ;
        for(auto &coordinate: value) coordinate /= scale ;
        local.push_back(value) ;
      }
      WidePoint normal{{0,0,0}} ;
      long double magnitude = 0 ;
      for(std::size_t i = 1; i < local.size(); ++i)
        for(std::size_t j = i+1; j < local.size(); ++j) {
          const auto candidate = cross(local[i],local[j]) ;
          const long double candidateLength = length(candidate) ;
          if(candidateLength > magnitude) { normal = candidate ; magnitude = candidateLength ; }
        }
      if(!(magnitude > 0)) throw std::invalid_argument("planar subdivision has collinear corners") ;
      const long double tolerance = 128*std::numeric_limits<double>::epsilon()*magnitude ;
      for(const auto &point: local)
        if(std::abs(dot(normal,point)) > tolerance)
          throw std::invalid_argument("planar subdivision requires planar corners") ;
      for(std::size_t i = 0; i < points.size(); ++i) {
        WidePoint first, second ;
        for(unsigned d = 0; d < 3; ++d) {
          first[d] = local[(i+1)%points.size()][d]-local[i][d] ;
          second[d] = local[(i+2)%points.size()][d]-local[(i+1)%points.size()][d] ;
        }
        if(!(dot(cross(first,second),normal) > 0))
          throw std::invalid_argument("planar subdivision requires a strictly convex polygon") ;
      }
    }

    Point interpolate(const Point &a, const Point &b, const Fraction &position) {
      Point result ;
      const long double t = position.rational().convert_to<long double>() ;
      const long double s = (Fraction(1)-position).rational().convert_to<long double>() ;
      for(unsigned d = 0; d < 3; ++d)
        result[d] = static_cast<double>(s*a[d]+t*b[d]) ;
      return result ;
    }

    Point weightedCenter(const std::vector<Point> &corners) {
      std::vector<long double> weights ;
      long double maximum = 0 ;
      for(std::size_t i = 0; i < corners.size(); ++i) {
        const long double weight = length(difference(corners[(i+1)%corners.size()],corners[i])) ;
        weights.push_back(weight) ;
        maximum = std::max(maximum,weight) ;
      }
      WidePoint sum{{0,0,0}} ;
      long double total = 0 ;
      for(std::size_t i = 0; i < corners.size(); ++i) {
        const long double weight = weights[i]/maximum ;
        const auto first = difference(corners[i],corners[0]) ;
        const auto second = difference(corners[(i+1)%corners.size()],corners[0]) ;
        for(unsigned d = 0; d < 3; ++d) sum[d] += weight*(first[d]+second[d])/2 ;
        total += weight ;
      }
      Point center ;
      for(unsigned d = 0; d < 3; ++d)
        center[d] = static_cast<double>(static_cast<long double>(corners[0][d])+sum[d]/total) ;
      return center ;
    }

  }

  void validatePlanarPolygon(const std::vector<Point> &corners) { validatePolygon(corners) ; }

  void validatePlanarSplit(const PlanarSplit &split, const std::vector<Point> &parentCorners) {
    validatePolygon(parentCorners) ;
    const unsigned n=parentCorners.size() ;
    if(split.shape != (n == 3 ? Shape::triangle : Shape::quad) ||
       split.points.size() != 2*n+1 || split.children.size() != n)
      throw std::invalid_argument("invalid accepted planar construction counts") ;
    if(split.centerMethod != CenterMethod::reference && split.centerMethod != CenterMethod::weighted &&
       split.centerMethod != CenterMethod::integrated)
      throw std::invalid_argument("unknown accepted planar center method") ;
    for(const auto &point : split.points) for(double value : point)
      if(!std::isfinite(value)) throw std::invalid_argument("accepted planar construction is not finite") ;
    for(unsigned i=0; i<n; ++i) {
      if(split.points[i] != parentCorners[i])
        throw std::invalid_argument("accepted planar corners differ from the parent") ;
      const auto expected=interpolate(parentCorners[i],parentCorners[(i+1)%n],Fraction(1,2)) ;
      for(unsigned d=0; d<3; ++d) {
        const double scale=std::max({std::abs(parentCorners[i][d]),std::abs(parentCorners[(i+1)%n][d]),
                                    std::abs(expected[d]),std::numeric_limits<double>::min()}) ;
        if(std::abs(split.points[n+i][d]-expected[d]) > 8*std::numeric_limits<double>::epsilon()*scale)
          throw std::invalid_argument("accepted planar midpoint differs from its parent edge") ;
      }
      const std::array<unsigned,4> indices{{i,n+i,2*n,n+(i+n-1)%n}} ;
      if(split.children[i] != indices)
        throw std::invalid_argument("accepted planar child has inconsistent connectivity") ;
      std::vector<Point> child ;
      for(unsigned vertex : indices) child.push_back(split.points[vertex]) ;
      validatePolygon(child) ;
      const auto parentNormal=cross(difference(parentCorners[1],parentCorners[0]),
                                    difference(parentCorners[2],parentCorners[0])) ;
      const auto childNormal=cross(difference(child[1],child[0]),difference(child[2],child[0])) ;
      if(!(dot(parentNormal,childNormal)>0))
        throw std::invalid_argument("accepted planar child reverses the parent orientation") ;
    }
  }

  PlanarSplit splitPlanarPolygon(const std::vector<Point> &corners, const std::vector<Point> &midpoints,
                                 CenterMethod method, MeasureOptions options) {
    validatePolygon(corners) ;
    const unsigned n=corners.size() ;
    if(midpoints.size() != n) throw std::invalid_argument("planar split requires one midpoint per edge") ;
    PlanarSplit split ; split.shape=n == 3 ? Shape::triangle : Shape::quad ; split.centerMethod=method ;
    Point center ;
    switch(method) {
    case CenterMethod::weighted: center=weightedCenter(corners) ; break ;
    case CenterMethod::reference: {
      CornerMap map(split.shape,corners) ;
      center=map.evaluate(n == 3 ? Point{{1./3,1./3,0}} : Point{{.5,.5,0}}).position ; break ;
    }
    case CenterMethod::integrated: {
      const auto measure=geometryMeasure(CornerMap(split.shape,corners),options) ;
      if(!measure.converged) throw std::length_error("planar center integration did not converge") ;
      center=measure.centroid ; break ;
    }
    default: throw std::invalid_argument("unknown planar center method") ;
    }
    split.points=corners ; split.points.insert(split.points.end(),midpoints.begin(),midpoints.end()) ;
    split.points.push_back(center) ;
    for(unsigned i=0; i<n; ++i) split.children.push_back({{i,n+i,2*n,n+(i+n-1)%n}}) ;
    validatePlanarSplit(split,corners) ; return split ;
  }

}
}
