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

#include <FVMAdapt/parametric/boundary_partition.h>
#include "exact_geometry.h"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace Loci {
namespace FVMAdapt {
  namespace {
    using detail::Exact ;
    using detail::ExactPoint ;
    using detail::exact ;
    using Polygon = std::vector<ExactPoint> ;

    struct Piece {
      Polygon vertices ;
      BoundaryUse use ;
      Exact measure ; ///< Length in 1D; twice the positive polygon area in 2D.
      std::array<Exact,2> lower, upper ;
    } ;

    /// Visit piece pairs whose bounds overlap, within one cover or across two covers.
    template<class Visit>
    void overlappingPieces(const std::vector<Piece> &first, const std::vector<Piece> *second,
                           unsigned dimension, Visit visit) {
      struct Entry { const Piece *piece ; unsigned side ; } ;
      std::vector<Entry> entries ;
      for(const Piece &piece : first) entries.push_back({&piece,0}) ;
      if(second)
        for(const Piece &piece : *second) entries.push_back({&piece,1}) ;
      std::array<Exact,2> length = {{0,0}} ;
      for(const Entry &entry : entries)
        for(unsigned d = 0; d < dimension; ++d) length[d] += entry.piece->upper[d]-entry.piece->lower[d] ;
      const unsigned axis = dimension == 2 && length[1] < length[0] ? 1 : 0 ;
      std::stable_sort(entries.begin(),entries.end(),[&](const Entry &a, const Entry &b) {
        return a.piece->lower[axis] < b.piece->lower[axis] ;
      }) ;
      std::array<std::set<const Piece *>,2> active ;
      std::multimap<Exact,Entry> ends ;
      for(const Entry &entry : entries) {
        while(!ends.empty() && ends.begin()->first <= entry.piece->lower[axis]) {
          const Entry old = ends.begin()->second ;
          active[old.side].erase(old.piece) ; ends.erase(ends.begin()) ;
        }
        for(const Piece *other : active[second ? 1-entry.side : 0]) {
          bool overlap = true ;
          for(unsigned d = 0; d < dimension; ++d)
            if(entry.piece->lower[d] >= other->upper[d] || other->lower[d] >= entry.piece->upper[d]) {
              overlap = false ; break ;
            }
          if(overlap) {
            if(entry.side == 0) visit(*entry.piece,*other) ;
            else visit(*other,*entry.piece) ;
          }
        }
        active[entry.side].insert(entry.piece) ;
        ends.emplace(entry.piece->upper[axis],entry) ;
      }
    }

    Exact cross(const ExactPoint &a, const ExactPoint &b, const ExactPoint &c) {
      return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]) ;
    }

    /// Return twice the signed area in shared boundary coordinates.
    Exact areaTwice(const Polygon &polygon) {
      Exact area = 0 ;
      for(std::size_t i = 0; i < polygon.size(); ++i) {
        const ExactPoint &a = polygon[i], &b = polygon[(i+1)%polygon.size()] ;
        area += a[0]*b[1]-a[1]*b[0] ;
      }
      return area ;
    }

    /// Remove redundant polygon vertices and choose a consistent traversal.
    void canonicalize(Polygon &polygon) {
      polygon.erase(std::unique(polygon.begin(),polygon.end()),polygon.end()) ;
      if(polygon.size() > 1 && polygon.front() == polygon.back()) polygon.pop_back() ;
      const Exact area = areaTwice(polygon) ;
      if(area == 0) {
        polygon.clear() ;
        return ;
      }
      if(area < 0) std::reverse(polygon.begin(),polygon.end()) ;
      for(std::size_t i = 0; i < polygon.size();) {
        const std::size_t previous = (i+polygon.size()-1)%polygon.size() ;
        const std::size_t next = (i+1)%polygon.size() ;
        if(cross(polygon[previous],polygon[i],polygon[next]) == 0) {
          polygon.erase(polygon.begin()+i) ;
          i = 0 ;
        } else ++i ;
      }
      std::rotate(polygon.begin(),std::min_element(polygon.begin(),polygon.end()),polygon.end()) ;
    }

    /// Intersect exact boundary segments or convex polygons.
    Polygon intersect(const Polygon &first, const Polygon &second, unsigned dimension) {
      if(dimension == 1) {
        const ExactPoint lower = std::max(first[0],second[0]) ;
        const ExactPoint upper = std::min(first[1],second[1]) ;
        return lower[0] < upper[0] ? Polygon{lower,upper} : Polygon{} ;
      }
      Polygon result = first ;
      for(std::size_t edge = 0; edge < second.size() && !result.empty(); ++edge) {
        const ExactPoint &a = second[edge], &b = second[(edge+1)%second.size()] ;
        Polygon clipped ;
        ExactPoint previous = result.back() ;
        Exact previousSide = cross(a,b,previous) ;
        for(const ExactPoint &current : result) {
          const Exact currentSide = cross(a,b,current) ;
          if((previousSide < 0) != (currentSide < 0)) {
            const Exact t = previousSide/(previousSide-currentSide) ;
            ExactPoint crossing ;
            for(unsigned d = 0; d < 3; ++d)
              crossing[d] = previous[d]+t*(current[d]-previous[d]) ;
            clipped.push_back(std::move(crossing)) ;
          }
          if(currentSide >= 0) clipped.push_back(current) ;
          previous = current ; previousSide = currentSide ;
        }
        result = std::move(clipped) ;
      }
      canonicalize(result) ;
      return result ;
    }

    /// Return exact segment length or twice polygon area in boundary coordinates.
    Exact measure(const Polygon &vertices, unsigned dimension) {
      if(vertices.empty()) return 0 ;
      return dimension == 1 ? Exact(vertices[1][0]-vertices[0][0]) : areaTwice(vertices) ;
    }

    /// Test whether exact coordinates lie in the original boundary.
    bool inFrame(Shape shape, const ExactPoint &point) {
      if(point[2] != 0 || point[0] < 0 || point[0] > 1) return false ;
      if(shape == Shape::segment) return point[1] == 0 ;
      if(point[1] < 0 || point[1] > 1) return false ;
      return shape != Shape::triangle || point[0]+point[1] <= 1 ;
    }

    /// Express one child boundary in its original boundary coordinates.
    Piece placedPiece(const BoundaryUse &use, Shape frame, unsigned dimension) {
      if((dimension == 1 && use.shape != Shape::segment) ||
         (dimension == 2 && use.shape != Shape::triangle && use.shape != Shape::quad))
        throw std::invalid_argument("boundary piece shape has the wrong dimension") ;
      const Placement &placement = use.placement ;
      ExactPoint origin, a, b ;
      for(unsigned d = 0; d < 3; ++d) {
        origin[d] = exact(placement.origin[d]) ;
        a[d] = exact(placement.columns[0][d]) ;
        b[d] = exact(placement.columns[1][d]) ;
      }
      for(unsigned d = dimension; d < 3; ++d) {
        if(origin[d] != 0 || a[d] != 0 || (dimension == 2 && b[d] != 0))
          throw std::invalid_argument("boundary placement has a nonzero unused coordinate") ;
      }
      const Exact determinant = dimension == 1 ? a[0] : a[0]*b[1]-a[1]*b[0] ;
      if(determinant == 0)
        throw std::invalid_argument("boundary placement is degenerate") ;
      Piece result ; result.use = use ;
      const unsigned count = referenceElement(use.shape).vertexCount ;
      result.vertices.reserve(count) ;
      for(unsigned i = 0; i < count; ++i) {
        const ReferencePoint vertex = referenceVertex(use.shape,i) ;
        ExactPoint point = origin ;
        for(unsigned local = 0; local < dimension; ++local)
          for(unsigned d = 0; d < dimension; ++d)
            point[d] += exact(placement.columns[local][d])*exact(vertex[local]) ;
        if(!inFrame(frame,point))
          throw std::invalid_argument("boundary placement is outside the reference frame") ;
        result.vertices.push_back(std::move(point)) ;
      }
      if(dimension == 1)
        std::sort(result.vertices.begin(),result.vertices.end()) ;
      else
        canonicalize(result.vertices) ;
      result.measure = measure(result.vertices,dimension) ;
      for(unsigned d = 0; d < 2; ++d) {
        result.lower[d] = result.upper[d] = result.vertices.front()[d] ;
        for(const ExactPoint &point : result.vertices) {
          result.lower[d] = std::min(result.lower[d],point[d]) ;
          result.upper[d] = std::max(result.upper[d],point[d]) ;
        }
      }
      return result ;
    }

    /// Validate a boundary cover and construct its nonoverlapping pieces.
    std::vector<Piece> pieces(const std::vector<BoundaryUse> &cover, Shape frame, unsigned dimension) {
      std::vector<Piece> result ;
      result.reserve(cover.size()) ;
      std::set<std::pair<unsigned,unsigned>> faces ;
      for(const BoundaryUse &use : cover) {
        if(!faces.insert({use.face.child,use.face.boundary}).second)
          throw std::invalid_argument("a boundary occurs more than once in a cover") ;
        result.push_back(placedPiece(use,frame,dimension)) ;
      }
      overlappingPieces(result,nullptr,dimension,[&](const Piece &a, const Piece &b) {
        if(!intersect(a.vertices,b.vertices,dimension).empty())
          throw std::invalid_argument("boundary cover pieces overlap") ;
      }) ;
      return result ;
    }

  }

  /// Convert an original-boundary point to the local coordinates of this piece.
  ReferencePoint BoundaryUse::coordinates(const ReferencePoint &point) const {
    if(shape != Shape::segment && shape != Shape::triangle && shape != Shape::quad)
      throw std::invalid_argument("boundary coordinates require a segment, triangle or quad") ;
    const unsigned dimension = shape == Shape::segment ? 1 : 2 ;
    for(unsigned d = dimension; d < 3; ++d) {
      if(point[d] != 0 || placement.origin[d] != 0)
        throw std::invalid_argument("boundary coordinates have a nonzero unused component") ;
      for(unsigned local = 0; local < dimension; ++local)
        if(placement.columns[local][d] != 0)
          throw std::invalid_argument("boundary placement has a nonzero unused component") ;
    }
    const Exact a = exact(placement.columns[0][0]), b = exact(placement.columns[1][0]) ;
    const Exact c = exact(placement.columns[0][1]), d = exact(placement.columns[1][1]) ;
    const Exact x = exact(point[0])-exact(placement.origin[0]) ;
    const Exact y = exact(point[1])-exact(placement.origin[1]) ;
    const Exact determinant = dimension == 1 ? a : a*d-b*c ;
    if(determinant == 0) throw std::invalid_argument("boundary placement is degenerate") ;
    ExactPoint local = {{0,0,0}} ;
    if(dimension == 1) local[0] = x/determinant ;
    else {
      local[0] = (d*x-b*y)/determinant ;
      local[1] = (a*y-c*x)/determinant ;
    }
    if(!inFrame(shape,local)) throw std::invalid_argument("point is outside the boundary piece") ;
    return {{Fraction(local[0]),Fraction(local[1]),Fraction(local[2])}} ;
  }

  /// Intersect two complete boundary covers and retain both owners of each piece.
  std::vector<BoundaryMatch> reconcileBoundary(Shape shape,
      const std::vector<BoundaryUse> &first, const std::vector<BoundaryUse> &second,
      std::size_t maxPieces) {
    if(shape != Shape::segment && shape != Shape::triangle && shape != Shape::quad)
      throw std::invalid_argument("boundary reconciliation requires a segment, triangle or quad") ;
    const unsigned dimension = shape == Shape::segment ? 1 : 2 ;
    const std::vector<Piece> a = pieces(first,shape,dimension), b = pieces(second,shape,dimension) ;
    Exact firstMeasure = 0, secondMeasure = 0, commonMeasure = 0 ;
    for(const Piece &piece : a) firstMeasure += piece.measure ;
    for(const Piece &piece : b) secondMeasure += piece.measure ;
    std::vector<BoundaryMatch> result ;
    overlappingPieces(a,&b,dimension,[&](const Piece &left, const Piece &right) {
        const Polygon intersection = intersect(left.vertices,right.vertices,dimension) ;
        if(intersection.empty()) return ;
        commonMeasure += measure(intersection,dimension) ;
        if(result.size() >= maxPieces || result.size() == result.max_size())
          throw std::length_error("boundary reconciliation exceeds the piece limit") ;
        BoundaryMatch match ;
        match.first = left.use ; match.second = right.use ;
        for(const ExactPoint &point : intersection)
          match.vertices.push_back({{Fraction(point[0]),Fraction(point[1]),Fraction(point[2])}}) ;
        result.push_back(std::move(match)) ;
    }) ;
    // Nonoverlapping covers agree when their entire measure belongs to the intersection.
    if(commonMeasure != firstMeasure || commonMeasure != secondMeasure)
      throw std::invalid_argument("boundary covers occupy different regions") ;
    std::sort(result.begin(),result.end(),[](const BoundaryMatch &a, const BoundaryMatch &b) {
      return a.vertices < b.vertices ;
    }) ;
    return result ;
  }

  /// Match polygon pieces from two complete covers of an original face.
  std::vector<PolygonMatch> reconcilePolygons(Shape frame,const std::vector<BoundaryPolygon> &first,
      const std::vector<BoundaryPolygon> &second,std::size_t maxPieces) {
    if(frame!=Shape::triangle && frame!=Shape::quad)
      throw std::invalid_argument("polygon reconciliation requires a triangle or quad frame") ;

    // Normalize each cover and check its pieces before intersecting the two covers.
    const auto prepareCover=[&](const std::vector<BoundaryPolygon> &input) {
      if(input.size()>maxPieces) throw std::length_error("boundary polygon limit exceeded") ;
      std::vector<Piece> result ;
      for(const auto &polygon:input) {
        Piece piece ; piece.use.face=polygon.face ;
        for(const auto &vertex:polygon.vertices) {
          ExactPoint point{{exact(vertex[0]),exact(vertex[1]),exact(vertex[2])}} ;
          if(!inFrame(frame,point)) throw std::invalid_argument("boundary polygon is outside its frame") ;
          piece.vertices.push_back(std::move(point)) ;
        }
        canonicalize(piece.vertices) ;
        if(piece.vertices.size()<3) throw std::invalid_argument("boundary polygon has zero area") ;
        for(unsigned i=0;i<piece.vertices.size();++i)
          if(cross(piece.vertices[i],piece.vertices[(i+1)%piece.vertices.size()],piece.vertices[(i+2)%piece.vertices.size()])<=0)
            throw std::invalid_argument("boundary polygon is not convex") ;
        piece.measure=measure(piece.vertices,2) ;
        for(unsigned d=0;d<2;++d) {
          piece.lower[d]=piece.upper[d]=piece.vertices.front()[d] ;
          for(const auto &point:piece.vertices) {
            piece.lower[d]=std::min(piece.lower[d],point[d]) ;
            piece.upper[d]=std::max(piece.upper[d],point[d]) ;
          }
        }
        result.push_back(std::move(piece)) ;
      }
      overlappingPieces(result,nullptr,2,[&](const Piece &a,const Piece &b) {
        if(!intersect(a.vertices,b.vertices,2).empty()) throw std::invalid_argument("boundary cover pieces overlap") ;
      }) ;
      return result ;
    } ;

    const auto firstPieces=prepareCover(first),secondPieces=prepareCover(second) ;
    Exact firstMeasure=0,secondMeasure=0,commonMeasure=0 ;
    for(const auto &piece:firstPieces) firstMeasure+=piece.measure ;
    for(const auto &piece:secondPieces) secondMeasure+=piece.measure ;
    std::vector<PolygonMatch> result ;
    overlappingPieces(firstPieces,&secondPieces,2,[&](const Piece &left,const Piece &right) {
      const auto polygon=intersect(left.vertices,right.vertices,2) ;
      if(polygon.empty()) return ;
      if(result.size()>=maxPieces) throw std::length_error("boundary polygon limit exceeded") ;
      PolygonMatch match ; match.first=left.use.face ; match.second=right.use.face ;
      for(const auto &point:polygon)
        match.vertices.push_back({{Fraction(point[0]),Fraction(point[1]),Fraction(point[2])}}) ;
      commonMeasure+=measure(polygon,2) ;
      result.push_back(std::move(match)) ;
    }) ;

    // All input area must be represented by the common pieces.
    if(commonMeasure!=firstMeasure || commonMeasure!=secondMeasure)
      throw std::invalid_argument("boundary covers have different areas") ;
    std::sort(result.begin(),result.end(),[](const PolygonMatch &a,const PolygonMatch &b) { return a.vertices<b.vertices ; }) ;
    return result ;
  }

}
}
