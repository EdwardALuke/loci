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

#include <FVMAdapt/parametric/vertex_volume_split.h>
#include "exact_geometry.h"
#include "geometry_bounds.h"
#include "geometry_validity.h"
#include <FVMAdapt/parametric/geometry_quality.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Loci {
namespace FVMAdapt {
  namespace {

    using detail::Exact ;
    using detail::ExactPoint ;
    using detail::BernsteinCoefficients ;
    using detail::hexJacobianCoefficients ;
    using detail::bisectJacobianCoefficients ;
    using detail::jacobianIndex ;
    using FaceCorners = std::array<unsigned,4> ;

    using detail::Bounds ;

    ExactPoint exactPoint(const Point &point) {
      for(double value : point)
        if(!std::isfinite(value))
          throw std::invalid_argument("vertex volume centers must be finite") ;
      return {{Exact(point[0]),Exact(point[1]),Exact(point[2])}} ;
    }

    Point approximate(const ExactPoint &point) {
      return {{point[0].convert_to<double>(),point[1].convert_to<double>(),point[2].convert_to<double>()}} ;
    }

    template<class Number> std::array<Number,3> difference(const std::array<Number,3> &a, const std::array<Number,3> &b) {
      return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}} ;
    }

    template<class Number> std::array<Number,3> midpoint(const std::array<Number,3> &a, const std::array<Number,3> &b) {
      return {{(a[0]+b[0])/2,(a[1]+b[1])/2,(a[2]+b[2])/2}} ;
    }

    template<class Number> Number determinant(const std::array<Number,3> &a, const std::array<Number,3> &b, const std::array<Number,3> &c) {
      return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+
             a[2]*(b[0]*c[1]-b[1]*c[0]) ;
    }

    void certifyChildJacobian(const std::vector<Point> &corners, JacobianCheckLimits limits) {
      if(positiveHexJacobian(corners)) return ;
      const auto bounds = hexJacobianCoefficients<Bounds>(corners) ;
      if(std::all_of(bounds.begin(),bounds.end(),[](Bounds b) { return b.low > 0 && std::isfinite(b.high) ; })) return ;
      struct Box { BernsteinCoefficients coefficients ; unsigned depth ; } ;
      std::vector<Box> pending{{hexJacobianCoefficients<Exact>(corners),0}} ;
      std::size_t examined = 0 ;
      while(!pending.empty()) {
        if(++examined > limits.maxBoxes)
          throw std::length_error("vertex volume Jacobian certification exceeded its box limit") ;
        auto box = std::move(pending.back()) ; pending.pop_back() ;
        bool positive = true ;
        for(const auto &coefficient : box.coefficients) if(coefficient <= 0) { positive = false ; break ; }
        if(positive) continue ;
        for(unsigned i : {0u,2u}) for(unsigned j : {0u,2u}) for(unsigned k : {0u,2u})
          if(box.coefficients[jacobianIndex(i,j,k)] <= 0)
            throw std::invalid_argument("vertex volume child has a nonpositive Jacobian") ;
        if(box.depth >= limits.maxDepth)
          throw std::length_error("vertex volume Jacobian certification exceeded its depth limit") ;
        auto children = bisectJacobianCoefficients(box.coefficients,box.depth%3) ;
        pending.push_back({std::move(children[1]),box.depth+1}) ;
        pending.push_back({std::move(children[0]),box.depth+1}) ;
      }
    }

    ExactPoint faceCenterCoordinates(Shape shape, unsigned boundary, const Point &point) {
      const auto local = exactPoint(point) ;
      const auto &element = referenceElement(shape) ;
      const auto &face = element.boundaries[boundary] ;
      const auto &reference = referenceElement(face.shape) ;
      if(local[2] != 0 || local[0] <= 0 || local[1] <= 0 ||
         (face.shape == Shape::triangle ? local[0]+local[1] >= 1 : local[0] >= 1 || local[1] >= 1))
        throw std::invalid_argument("vertex volume face center must lie strictly inside its boundary") ;
      // Corner patches are convex when the center lies inside the midpoint polygon.
      const Bounds x(point[0]), y(point[1]), sum = x+y, delta = x-y ;
      const bool convex = face.shape == Shape::triangle ? point[0] < .5 && point[1] < .5 && sum.low > .5
        : sum.low > .5 && sum.high < 1.5 && delta.low > -.5 && delta.high < .5 ;
      for(unsigned v = 0; !convex && v < face.vertexCount; ++v) {
        const auto corner = exactPoint(reference.vertices[v]) ;
        const std::array<ExactPoint,4> patch{{corner,
          midpoint(corner,exactPoint(reference.vertices[(v+1)%face.vertexCount])), local,
          midpoint(corner,exactPoint(reference.vertices[(v+face.vertexCount-1)%face.vertexCount]))}} ;
        for(unsigned p = 0; p < 4; ++p) {
          const auto first = difference(patch[(p+1)%4],patch[p]) ;
          const auto second = difference(patch[(p+2)%4],patch[(p+1)%4]) ;
          if(first[0]*second[1]-first[1]*second[0] <= 0)
            throw std::invalid_argument("vertex volume boundary patch must be strictly convex") ;
        }
      }
      ExactPoint result{} ;
      // Native triangle and quad boundaries have affine reference embeddings.
      const auto &origin = element.vertices[face.vertices[0]] ;
      const auto &u = element.vertices[face.vertices[1]], &v = element.vertices[face.vertices[face.vertexCount-1]] ;
      for(unsigned d = 0; d < 3; ++d)
        result[d] = Exact(origin[d])+local[0]*Exact(u[d]-origin[d])+local[1]*Exact(v[d]-origin[d]) ;
      return result ;
    }

    unsigned edgeIndex(const ReferenceElement &element, unsigned a, unsigned b) {
      for(unsigned edge = 0; edge < element.edgeCount; ++edge) {
        const auto &ends = element.edges[edge] ;
        if((ends[0] == static_cast<int>(a) && ends[1] == static_cast<int>(b)) || (ends[0] == static_cast<int>(b) && ends[1] == static_cast<int>(a))) return edge ;
      }
      throw std::logic_error("vertex volume edge incidence") ;
    }

    unsigned faceIndex(const ReferenceElement &element, unsigned a, unsigned b, unsigned c) {
      for(unsigned face = 0; face < element.boundaryCount; ++face) {
        const auto &vertices = element.boundaries[face].vertices ;
        const auto end = vertices.begin()+element.boundaries[face].vertexCount ;
        if(std::find(vertices.begin(),end,a) != end && std::find(vertices.begin(),end,b) != end &&
           std::find(vertices.begin(),end,c) != end) return face ;
      }
      throw std::logic_error("vertex volume face incidence") ;
    }

    std::array<unsigned,8> childVertices(const ReferenceElement &element, unsigned vertex) {
      std::vector<unsigned> neighbors ;
      for(unsigned edge = 0; edge < element.edgeCount; ++edge)
        for(unsigned end = 0; end < 2; ++end)
          if(element.edges[edge][end] == static_cast<int>(vertex)) neighbors.push_back(element.edges[edge][1-end]) ;
      std::sort(neighbors.begin(),neighbors.end()) ;
      if(neighbors.size() != 3) throw std::logic_error("vertex volume requires three incident edges") ;
      const auto origin = exactPoint(element.vertices[vertex]) ;
      if(determinant(difference(exactPoint(element.vertices[neighbors[0]]),origin),
                     difference(exactPoint(element.vertices[neighbors[1]]),origin),
                     difference(exactPoint(element.vertices[neighbors[2]]),origin)) < 0)
        std::swap(neighbors[1],neighbors[2]) ;
      const unsigned a = neighbors[0], b = neighbors[1], c = neighbors[2] ;
      return {{vertex,element.vertexCount+edgeIndex(element,vertex,a),
        element.vertexCount+element.edgeCount+faceIndex(element,vertex,a,b),element.vertexCount+edgeIndex(element,vertex,b),
        element.vertexCount+edgeIndex(element,vertex,c),element.vertexCount+element.edgeCount+faceIndex(element,vertex,a,c),
        element.vertexCount+element.edgeCount+element.boundaryCount,
        element.vertexCount+element.edgeCount+faceIndex(element,vertex,b,c)}} ;
    }

    // Connectivity depends only on the parent shape, not on the accepted centers.
    struct VertexPattern {
      std::vector<std::array<unsigned,8>> children ;
      std::vector<InternalFace> internalFaces ;
      std::vector<VertexVolumeBoundary> boundaryPieces ;
    } ;

    FaceCorners sorted(FaceCorners corners) {
      std::sort(corners.begin(),corners.end()) ; return corners ;
    }

    unsigned orientation(const FaceCorners &first, const FaceCorners &second) {
      for(unsigned o = 0; o < 8; ++o) {
        bool matches = true ;
        for(unsigned v = 0; v < 4; ++v) if(first[v] != second[orientedVertex(Shape::quad,o,v)]) matches = false ;
        if(matches) return o ;
      }
      throw std::logic_error("vertex volume face orientation") ;
    }

    void connect(Shape shape, VertexPattern &result) {
      const auto &parent = referenceElement(shape) ;
      const auto &hex = referenceElement(Shape::hex) ;
      struct OuterFace { unsigned boundary, corner ; FaceCorners vertices ; } ;
      struct InnerFace { FaceUse face ; FaceCorners vertices ; } ;
      std::map<FaceCorners,OuterFace> boundary ;
      for(unsigned b = 0; b < parent.boundaryCount; ++b) {
        const auto &face = parent.boundaries[b] ;
        for(unsigned v = 0; v < face.vertexCount; ++v) {
          const unsigned corner = face.vertices[v] ;
          const FaceCorners points{{corner,
            parent.vertexCount+edgeIndex(parent,corner,face.vertices[(v+1)%face.vertexCount]),
            parent.vertexCount+parent.edgeCount+b,
            parent.vertexCount+edgeIndex(parent,corner,face.vertices[(v+face.vertexCount-1)%face.vertexCount])}} ;
          boundary.emplace(sorted(points),OuterFace{b,v,points}) ;
        }
      }
      std::map<FaceCorners,InnerFace> pending ;
      for(unsigned c = 0; c < result.children.size(); ++c)
        for(unsigned b = 0; b < hex.boundaryCount; ++b) {
          FaceCorners corners ;
          for(unsigned v = 0; v < 4; ++v) corners[v] = result.children[c][hex.boundaries[b].vertices[v]] ;
          const auto key = sorted(corners) ;
          const auto outer = boundary.find(key) ;
          if(outer != boundary.end()) {
            const unsigned o = orientation(corners,outer->second.vertices) ;
            if(o >= 4) throw std::logic_error("vertex volume exterior face orientation") ;
            result.boundaryPieces.push_back({{c,b},outer->second.boundary,outer->second.corner,o}) ;
            boundary.erase(outer) ;
          } else {
            const auto inner = pending.find(key) ;
            if(inner == pending.end()) pending.emplace(key,InnerFace{{c,b},corners}) ;
            else {
              const unsigned o = orientation(inner->second.vertices,corners) ;
              if(o < 4) throw std::logic_error("vertex volume interior face orientation") ;
              result.internalFaces.push_back({inner->second.face,{c,b},o}) ;
              pending.erase(inner) ;
            }
          }
        }
      if(!pending.empty() || !boundary.empty()) throw std::logic_error("vertex volume incomplete boundary partition") ;
    }

    const VertexPattern &childPattern(Shape shape) {
      const auto build = [](Shape shape) {
        const auto &element = referenceElement(shape) ; VertexPattern result ;
        for(unsigned v = 0; v < element.vertexCount; ++v) result.children.push_back(childVertices(element,v)) ;
        connect(shape,result) ;
        return result ;
      } ;
      switch(shape) {
      case Shape::hex: { static const auto value = build(Shape::hex) ; return value ; }
      case Shape::prism: { static const auto value = build(Shape::prism) ; return value ; }
      case Shape::tet: { static const auto value = build(Shape::tet) ; return value ; }
      default: throw std::invalid_argument("vertex volume subdivision requires a hex, prism or tet") ;
      }
    }

    std::vector<Point> splitPoints(Shape shape, const std::vector<Point> &faceCenters,
                                   const Point &cellCenter, JacobianCheckLimits limits) {
      if(shape != Shape::hex && shape != Shape::prism && shape != Shape::tet)
        throw std::invalid_argument("vertex volume subdivision requires a hex, prism or tet") ;
      if(limits.maxBoxes == 0) throw std::invalid_argument("vertex volume certification requires a positive box limit") ;
      const auto &element = referenceElement(shape) ;
      if(faceCenters.size() != element.boundaryCount)
        throw std::invalid_argument("vertex volume subdivision requires one center per face") ;
      const auto center = exactPoint(cellCenter) ;
      if(center[0] <= 0 || center[1] <= 0 || center[2] <= 0 || center[2] >= 1 ||
         (shape == Shape::hex ? center[0] >= 1 || center[1] >= 1 :
          shape == Shape::prism ? center[0]+center[1] >= 1 : center[0]+center[1]+center[2] >= 1))
        throw std::invalid_argument("vertex volume cell center must lie strictly inside the parent") ;
      std::vector<ExactPoint> points ;
      for(unsigned v = 0; v < element.vertexCount; ++v) points.push_back(exactPoint(element.vertices[v])) ;
      for(unsigned e = 0; e < element.edgeCount; ++e)
        points.push_back(midpoint(points[element.edges[e][0]],points[element.edges[e][1]])) ;
      for(unsigned b = 0; b < element.boundaryCount; ++b) points.push_back(faceCenterCoordinates(shape,b,faceCenters[b])) ;
      points.push_back(center) ;
      std::vector<Point> result ; for(const auto &point: points) result.push_back(approximate(point)) ; return result ;
    }

  }

  VertexVolumeSplit splitVertexVolume(Shape shape, const std::vector<Point> &faceCenters,
                                     const Point &cellCenter, JacobianCheckLimits limits) {
    VertexVolumeSplit result ; result.shape = shape ;
    result.points = splitPoints(shape,faceCenters,cellCenter,limits) ;
    const auto &pattern=childPattern(shape) ;
    for(const auto &vertices: pattern.children) {
      std::vector<Point> corners ; for(unsigned vertex: vertices) corners.push_back(result.points[vertex]) ;
      certifyChildJacobian(corners,limits) ;
      result.children.push_back({vertices,std::make_shared<CornerMap>(Shape::hex,std::move(corners))}) ;
    }
    result.internalFaces=pattern.internalFaces ; result.boundaryPieces=pattern.boundaryPieces ;
    return result ;
  }

  void validateVertexVolumeSplit(Shape shape, const std::vector<Point> &faceCenters,
                                 const Point &cellCenter, JacobianCheckLimits limits) {
    const auto points = splitPoints(shape,faceCenters,cellCenter,limits) ;
    std::vector<Point> corners(8) ;
    for(const auto &vertices: childPattern(shape).children) {
      for(unsigned v = 0; v < 8; ++v) corners[v] = points[vertices[v]] ;
      certifyChildJacobian(corners,limits) ;
    }
  }

  VertexVolumeSplit restoreVertexVolumeSplit(Shape shape, std::vector<Point> points,
      const std::vector<std::array<unsigned,8>> &children, JacobianCheckLimits limits) {
    if(shape != Shape::hex && shape != Shape::prism && shape != Shape::tet)
      throw std::invalid_argument("vertex volume subdivision requires a hex, prism or tet") ;
    if(limits.maxBoxes == 0) throw std::invalid_argument("vertex volume certification requires a positive box limit") ;
    const auto &element = referenceElement(shape) ;
    const unsigned faceOffset = element.vertexCount+element.edgeCount ;
    if(points.size() != faceOffset+element.boundaryCount+1 || children.size() != element.vertexCount)
      throw std::invalid_argument("invalid accepted vertex volume counts") ;
    for(const auto &point : points) exactPoint(point) ;
    for(unsigned vertex = 0; vertex < element.vertexCount; ++vertex)
      if(points[vertex] != element.vertices[vertex])
        throw std::invalid_argument("accepted vertex volume moved an original corner") ;
    for(unsigned edge = 0; edge < element.edgeCount; ++edge)
      if(points[element.vertexCount+edge] != approximate(midpoint(
          exactPoint(points[element.edges[edge][0]]),exactPoint(points[element.edges[edge][1]]))))
        throw std::invalid_argument("accepted vertex volume moved an edge midpoint") ;
    for(unsigned boundary = 0; boundary < element.boundaryCount; ++boundary) {
      const auto &face = element.boundaries[boundary] ;
      const Point &origin = element.vertices[face.vertices[0]] ;
      const Point &u = element.vertices[face.vertices[1]] ;
      const Point &v = element.vertices[face.vertices[face.vertexCount-1]] ;
      Point a{}, b{}, q{} ; double aa=0, ab=0, bb=0, aq=0, bq=0 ;
      for(unsigned d=0; d<3; ++d) {
        a[d]=u[d]-origin[d] ; b[d]=v[d]-origin[d] ; q[d]=points[faceOffset+boundary][d]-origin[d] ;
        aa+=a[d]*a[d] ; ab+=a[d]*b[d] ; bb+=b[d]*b[d] ; aq+=a[d]*q[d] ; bq+=b[d]*q[d] ;
      }
      const double determinant = aa*bb-ab*ab ;
      const Point local{{(aq*bb-bq*ab)/determinant,(bq*aa-aq*ab)/determinant,0}} ;
      const Point embedded = approximate(faceCenterCoordinates(shape,boundary,local)) ;
      for(unsigned d=0; d<3; ++d)
        if(std::abs(embedded[d]-points[faceOffset+boundary][d]) > 8*std::numeric_limits<double>::epsilon())
          throw std::invalid_argument("accepted vertex volume face center is outside its boundary") ;
    }
    const auto &center = points.back() ;
    if(center[0] <= 0 || center[1] <= 0 || center[2] <= 0 || center[2] >= 1 ||
       (shape == Shape::hex ? center[0] >= 1 || center[1] >= 1 :
        shape == Shape::prism ? center[0]+center[1] >= 1 : center[0]+center[1]+center[2] >= 1))
      throw std::invalid_argument("accepted vertex volume center is outside its cell") ;
    VertexVolumeSplit result ; result.shape=shape ; result.points=std::move(points) ;
    const auto &pattern=childPattern(shape) ;
    for(unsigned ordinal=0; ordinal<children.size(); ++ordinal) {
      if(children[ordinal] != pattern.children[ordinal])
        throw std::invalid_argument("accepted vertex volume child has invalid corner incidence") ;
      std::vector<Point> corners ;
      for(unsigned vertex : children[ordinal]) corners.push_back(result.points[vertex]) ;
      certifyChildJacobian(corners,limits) ;
      result.children.push_back({children[ordinal],std::make_shared<CornerMap>(Shape::hex,std::move(corners))}) ;
    }
    result.internalFaces=pattern.internalFaces ; result.boundaryPieces=pattern.boundaryPieces ;
    return result ;
  }

  CornerMap quadChildMap(unsigned corner,const Point &center) {
    if(corner>=4) throw std::out_of_range("quad child corner") ;
    const auto &q=referenceElement(Shape::quad).vertices ;
    const auto midpoint=[](const Point &a,const Point &b) {
      Point p ; for(unsigned d=0;d<3;++d) p[d]=.5*a[d]+.5*b[d] ; return p ;
    } ;
    return CornerMap(Shape::quad,{q[corner],midpoint(q[corner],q[(corner+1)%4]),center,midpoint(q[corner],q[(corner+3)%4])}) ;
  }

}
}
