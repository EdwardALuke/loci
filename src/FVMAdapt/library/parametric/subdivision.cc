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

#include <FVMAdapt/parametric/subdivision.h>
#include "exact_geometry.h"
#include "subdivision_internal.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Loci {
namespace FVMAdapt {

  SplitOperation isotropicSplit(Shape shape) {
    switch(shape) {
    case Shape::hex: {
      TensorSplit request ;
      for(auto &cuts: request.cuts) cuts = {Fraction(1,2)} ;
      return request ;
    }
    case Shape::prism: return PrismSplit{TriangleSplit::uniform4(),{Fraction(1,2)}} ;
    case Shape::tet: return TetSplit::uniform8(0) ;
    case Shape::pyramid: return PyramidSplit::uniform10() ;
    default: throw std::invalid_argument("isotropic refinement requires a native volume cell") ;
    }
  }

  std::optional<SplitOperation> automaticSplit(const GeometryMap &geometry, double aspectRatio, double minEdgeLength) {
    if(!std::isfinite(aspectRatio) || aspectRatio<=1 || !std::isfinite(minEdgeLength) || minEdgeLength<0)
      throw std::invalid_argument("automatic splitting requires a finite aspect ratio greater than one and a nonnegative edge cutoff") ;
    const auto shape=geometry.shape() ;
    if(shape!=Shape::hex && shape!=Shape::prism && shape!=Shape::tet && shape!=Shape::pyramid)
      throw std::invalid_argument("automatic splitting requires a native volume cell") ;
    const auto &element=referenceElement(shape) ;
    std::array<Point,8> points ;
    for(unsigned v=0;v<element.vertexCount;++v) points[v]=geometry.evaluate(element.vertices[v]).position ;
    std::array<double,3> mean{}, shortest ; shortest.fill(std::numeric_limits<double>::infinity()) ;
    std::array<unsigned,3> count{} ;
    for(unsigned e=0;e<element.edgeCount;++e) {
      const auto a=element.edges[e][0], b=element.edges[e][1] ; unsigned group=0 ;
      if(shape==Shape::hex) {
        while(element.vertices[a][group]==element.vertices[b][group]) ++group ;
      } else if(shape==Shape::prism) group=element.vertices[a][2]!=element.vertices[b][2] ;
      const auto &p=points[a], &q=points[b] ;
      const double length=std::hypot(p[0]-q[0],p[1]-q[1],p[2]-q[2]) ;
      if(!std::isfinite(length) || length<=0) throw std::invalid_argument("automatic splitting requires positive finite edge lengths") ;
      ++count[group] ; mean[group]+=(length-mean[group])/count[group] ;
      shortest[group]=std::min(shortest[group],length) ;
    }
    const unsigned groups=shape==Shape::hex ? 3:shape==Shape::prism ? 2:1 ;
    const double minimum=*std::min_element(mean.begin(),mean.begin()+groups) ;
    std::array<bool,3> selected{} ; bool stretched=false, any=false ;
    for(unsigned g=0;g<groups;++g) stretched=stretched || mean[g]/minimum>aspectRatio ;
    for(unsigned g=0;g<groups;++g) {
      selected[g]=(!stretched || mean[g]/minimum>aspectRatio) && shortest[g]/2>=minEdgeLength ;
      any=any || selected[g] ;
    }
    if(!any) return std::nullopt ;
    if(shape==Shape::hex) {
      TensorSplit request ;
      for(unsigned g=0;g<3;++g) if(selected[g]) request.cuts[g]={Fraction(1,2)} ;
      return request ;
    }
    if(shape==Shape::prism) {
      PrismSplit request ;
      if(selected[0]) request.base=TriangleSplit::uniform4() ;
      if(selected[1]) request.axialCuts={Fraction(1,2)} ;
      return request ;
    }
    return isotropicSplit(shape) ;
  }

  Placement::Placement() : origin{{0,0,0}},
    columns{{{{1,0,0}}, {{0,1,0}}, {{0,0,1}}}} {}

  ReferencePoint Placement::apply(const ReferencePoint &point) const {
    ReferencePoint result ;
    for(unsigned d = 0; d < 3; ++d) {
      detail::Exact value = detail::exact(origin[d]) ;
      for(unsigned i = 0; i < 3; ++i)
        value += detail::exact(columns[i][d])*detail::exact(point[i]) ;
      result[d] = Fraction(value) ;
    }
    return result ;
  }

  Placement Placement::compose(const Placement &child) const {
    Placement result ;
    result.origin = apply(child.origin) ;
    for(unsigned i = 0; i < 3; ++i)
      for(unsigned d = 0; d < 3; ++d) {
        detail::Exact value = 0 ;
        for(unsigned j = 0; j < 3; ++j)
          value += detail::exact(columns[j][d])*detail::exact(child.columns[i][j]) ;
        result.columns[i][d] = Fraction(value) ;
      }
    return result ;
  }

  namespace {

    class AffineMap final : public GeometryMap {
      const Shape shape_ ;
      const unsigned dimension_ ;
      const Point origin_ ;
      std::array<Point,3> jacobian_ ;
    public:
      AffineMap(Shape shape, const Placement &placement) : shape_(shape),
        dimension_(referenceElement(shape).dimension), origin_(approximate(placement.origin)), jacobian_{} {
        for(unsigned d = 0; d < dimension_; ++d)
          jacobian_[d] = approximate(placement.columns[d]) ;
      }
      Shape shape() const override { return shape_ ; }
      MapValue evaluate(const Point &point) const override {
        for(double component : point)
          if(!std::isfinite(component))
            throw std::invalid_argument("reference coordinates must be finite") ;
        MapValue result = {origin_,jacobian_,true} ;
        for(unsigned i = 0; i < dimension_; ++i)
          for(unsigned d = 0; d < 3; ++d)
            result.position[d] += jacobian_[i][d]*point[i] ;
        return result ;
      }
    } ;

    struct BoundaryDirection {
      unsigned axis ;
      unsigned side ;
    } ;

    std::size_t checkedProduct(std::size_t a, std::size_t b) {
      if(b && a > std::numeric_limits<std::size_t>::max()/b)
        throw std::invalid_argument("subdivision count exceeds addressable size") ;
      return a*b ;
    }

    std::size_t checkedSum(std::size_t a, std::size_t b) {
      if(a > std::numeric_limits<std::size_t>::max()-b)
        throw std::invalid_argument("subdivision count exceeds addressable size") ;
      return a+b ;
    }

    template<class T> void checkCapacity(const std::vector<T> &values, std::size_t count) {
      if(count > values.max_size())
        throw std::invalid_argument("subdivision count exceeds container capacity") ;
    }

    BoundaryDirection boundaryDirection(const ReferenceElement &element, unsigned boundary) {
      const Boundary &face = element.boundaries[boundary] ;
      for(unsigned axis = 0; axis < element.dimension; ++axis) {
        const double side = element.vertices[face.vertices[0]][axis] ;
        bool fixed = true ;
        for(unsigned v = 1; v < face.vertexCount; ++v)
          fixed = fixed && element.vertices[face.vertices[v]][axis] == side ;
        if(fixed) return {axis,unsigned(side)} ;
      }
      throw std::logic_error("tensor boundary has no fixed coordinate") ;
    }

    unsigned edgeAxis(const ReferenceElement &element, unsigned edge) {
      const Point &first = element.vertices[element.edges[edge][0]] ;
      const Point &second = element.vertices[element.edges[edge][1]] ;
      for(unsigned d = 0; d < element.dimension; ++d)
        if(first[d] != second[d]) return d ;
      throw std::logic_error("tensor edge has no varying coordinate") ;
    }

    unsigned interfaceOrientation(const ReferenceElement &element, unsigned first,
                                  unsigned second, unsigned axis) {
      const Boundary &a = element.boundaries[first], &b = element.boundaries[second] ;
      for(unsigned orientation = 0; orientation < orientationCount(a.shape); ++orientation) {
        bool matches = true ;
        for(unsigned v = 0; v < a.vertexCount; ++v) {
          const Point &p = element.vertices[a.vertices[v]] ;
          const Point &q = element.vertices[b.vertices[orientedVertex(b.shape,orientation,v)]] ;
          for(unsigned d = 0; d < element.dimension; ++d)
            if(d != axis && p[d] != q[d]) matches = false ;
        }
        if(matches) return orientation ;
      }
      throw std::logic_error("tensor boundaries have incompatible orientations") ;
    }

    struct TensorGrid {
      std::array<std::size_t,3> count = {{1,1,1}} ;
      std::array<std::vector<Fraction>,3> coordinates ;
      std::size_t childCount = 1 ;
      TensorGrid(Shape shape, const TensorSplit &request) {
        if(shape != Shape::segment && shape != Shape::quad && shape != Shape::hex)
          throw std::invalid_argument("tensor subdivision requires a segment, quad or hex") ;
        const ReferenceElement &element = referenceElement(shape) ;
        const unsigned dimension = element.dimension ;
        for(unsigned d = 0; d < 3; ++d) {
          const std::vector<Fraction> &cuts = request.cuts[d] ;
          if(d >= dimension && !cuts.empty())
            throw std::invalid_argument("unused subdivision direction has cuts") ;
          if(cuts.size() >= std::numeric_limits<unsigned>::max())
            throw std::invalid_argument("subdivision has too many cuts") ;
          Fraction previous = 0 ;
          for(const Fraction &cut : cuts) {
            if(cut <= previous || cut >= Fraction(1))
              throw std::invalid_argument("subdivision cuts must increase strictly inside (0,1)") ;
            previous = cut ;
          }
          count[d] = cuts.size()+1 ; childCount = checkedProduct(childCount,count[d]) ;
        }
        if(childCount > std::numeric_limits<unsigned>::max())
          throw std::invalid_argument("subdivision has too many children") ;
        for(unsigned d = 0; d < 3; ++d) {
          checkCapacity(coordinates[d],checkedSum(count[d],1)) ;
          coordinates[d].reserve(count[d]+1) ; coordinates[d].push_back(0) ;
          coordinates[d].insert(coordinates[d].end(),request.cuts[d].begin(),request.cuts[d].end()) ;
          coordinates[d].push_back(1) ;
        }
      }
    } ;

    std::vector<Child> tensorChildren(Shape shape, const TensorGrid &grid) {
      const unsigned dimension = referenceElement(shape).dimension ;
      std::vector<Child> children ;
      checkCapacity(children,grid.childCount) ; children.reserve(grid.childCount) ;
      for(std::size_t k = 0; k < grid.count[2]; ++k)
        for(std::size_t j = 0; j < grid.count[1]; ++j)
          for(std::size_t i = 0; i < grid.count[0]; ++i) {
            const std::array<std::size_t,3> index = {{i,j,k}} ;
            Child child ; child.shape = shape ;
            for(unsigned d = 0; d < dimension; ++d) {
              child.placement.origin[d] = grid.coordinates[d][index[d]] ;
              child.placement.columns[d][d] = Fraction(detail::exact(grid.coordinates[d][index[d]+1])-
                                                           detail::exact(grid.coordinates[d][index[d]])) ;
            }
            children.push_back(child) ;
          }
      return children ;
    }

  }

  std::shared_ptr<const GeometryMap> Placement::map(Shape shape) const {
    return std::make_shared<AffineMap>(shape,*this) ;
  }

  std::vector<Child> detail::splitChildren(Shape shape, const TensorSplit &request) {
    return tensorChildren(shape,TensorGrid(shape,request)) ;
  }

  SplitResult split(Shape shape, const TensorSplit &request) {
    const TensorGrid grid(shape,request) ;
    const ReferenceElement &element = referenceElement(shape) ;
    const unsigned dimension = element.dimension ;
    const auto &count = grid.count ;
    const std::size_t childCount = grid.childCount ;

    std::array<BoundaryDirection,6> directions ;
    std::array<std::array<unsigned,2>,3> axisBoundaries = {} ;
    for(unsigned b = 0; b < element.boundaryCount; ++b) {
      directions[b] = boundaryDirection(element,b) ;
      axisBoundaries[directions[b].axis][directions[b].side] = b ;
    }
    std::array<unsigned,3> orientations = {{0,0,0}} ;
    std::size_t internalCount = 0, boundaryCount = 0, edgeCount = 0 ;
    for(unsigned d = 0; d < dimension; ++d) {
      const std::size_t crossSection = childCount/count[d] ;
      internalCount = checkedSum(internalCount,checkedProduct(count[d]-1,crossSection)) ;
      boundaryCount = checkedSum(boundaryCount,checkedProduct(2,crossSection)) ;
      orientations[d] = interfaceOrientation(element,axisBoundaries[d][1],axisBoundaries[d][0],d) ;
    }
    for(unsigned e = 0; e < element.edgeCount; ++e)
      edgeCount = checkedSum(edgeCount,count[edgeAxis(element,e)]) ;

    SplitResult result ; result.shape = shape ;
    checkCapacity(result.children,childCount) ; checkCapacity(result.internalFaces,internalCount) ;
    checkCapacity(result.boundaryPieces,boundaryCount) ; checkCapacity(result.edgePieces,edgeCount) ;
    result.children = tensorChildren(shape,grid) ; result.internalFaces.reserve(internalCount) ;
    result.boundaryPieces.reserve(boundaryCount) ; result.edgePieces.reserve(edgeCount) ;
    const std::array<std::size_t,3> stride = {{1,count[0],checkedProduct(count[0],count[1])}} ;

    for(std::size_t k = 0; k < count[2]; ++k)
      for(std::size_t j = 0; j < count[1]; ++j)
        for(std::size_t i = 0; i < count[0]; ++i) {
          const std::array<std::size_t,3> index = {{i,j,k}} ;
          const unsigned child = unsigned(i+count[0]*(j+count[1]*k)) ;
          const Child &current = result.children[child] ;

          for(unsigned b = 0; b < element.boundaryCount; ++b) {
            const unsigned axis = directions[b].axis, side = directions[b].side ;
            if(side == 0 && index[axis] != 0) continue ;
            if(side == 1 && index[axis]+1 < count[axis]) {
              const FaceUse first = {child,b} ;
              const FaceUse second = {unsigned(child+stride[axis]),axisBoundaries[axis][0]} ;
              result.internalFaces.push_back({first,second,orientations[axis]}) ;
            } else {
              BoundaryPiece piece ; piece.face = {child,b} ; piece.boundary = b ;
              const Boundary &face = element.boundaries[b] ;
              piece.vertices.reserve(face.vertexCount) ;
              for(unsigned v = 0; v < face.vertexCount; ++v) {
                const ReferencePoint point = current.placement.apply(referenceVertex(shape,face.vertices[v])) ;
                piece.vertices.push_back(boundaryCoordinates(shape,b,point)) ;
              }
              result.boundaryPieces.push_back(std::move(piece)) ;
            }
          }

          for(unsigned e = 0; e < element.edgeCount; ++e) {
            const unsigned axis = edgeAxis(element,e) ;
            const ReferencePoint first = referenceVertex(shape,element.edges[e][0]) ;
            bool onParentEdge = true ;
            for(unsigned d = 0; d < dimension; ++d) {
              if(d == axis) continue ;
              if((first[d] == Fraction(0) && index[d] != 0) ||
                 (first[d] == Fraction(1) && index[d]+1 != count[d])) onParentEdge = false ;
            }
            if(!onParentEdge) continue ;
            const ReferencePoint second = referenceVertex(shape,element.edges[e][1]) ;
            const ReferencePoint a = current.placement.apply(first), b = current.placement.apply(second) ;
            const Fraction direction = second[axis]-first[axis] ;
            result.edgePieces.push_back({child,e,e,{{(a[axis]-first[axis])/direction,
                                                    (b[axis]-first[axis])/direction}}}) ;
          }
        }
    return result ;
  }

}
}
