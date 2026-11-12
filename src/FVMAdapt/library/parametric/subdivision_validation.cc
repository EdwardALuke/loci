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

#include "subdivision_internal.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace Loci {
namespace FVMAdapt {
namespace detail {
  namespace {
    using Vertices = std::vector<ReferencePoint> ;
    using Pattern = std::vector<Vertices> ;

    ReferencePoint midpoint(const ReferencePoint &a, const ReferencePoint &b) {
      return {{(a[0]+b[0])/Fraction(2),(a[1]+b[1])/Fraction(2),(a[2]+b[2])/Fraction(2)}} ;
    }

    std::size_t product(std::size_t a, std::size_t b) {
      if(b && a > std::numeric_limits<std::size_t>::max()/b)
        throw std::length_error("refinement count exceeds addressable size") ;
      return a*b ;
    }

    void count(const std::vector<Child> &children, std::size_t expected) {
      if(children.size() != expected)
        throw std::invalid_argument("accepted refinement count differs from request") ;
    }

    void cuts(const std::vector<Fraction> &values, ExactBudget &budget) {
      budget.operations(values.size()) ;
      Fraction previous = 0 ;
      for(const Fraction &value : values) {
        budget.check(value) ;
        if(value <= previous || value >= Fraction(1))
          throw std::invalid_argument("refinement cuts must increase strictly inside (0,1)") ;
        previous = value ;
      }
    }

    std::array<Fraction,2> interval(const std::vector<Fraction> &cuts, std::size_t index) {
      return {{index ? cuts[index-1] : Fraction(0), index < cuts.size() ? cuts[index] : Fraction(1)}} ;
    }

    void childMatches(const Child &child, Shape shape, Vertices expected, ExactBudget &budget) {
      if(child.shape != shape)
        throw std::invalid_argument("accepted refinement child has the wrong shape") ;
      const auto &element = referenceElement(shape) ;
      budget.operations(32+24*element.vertexCount) ;
      const auto &columns = child.placement.columns ;
      for(const auto &value : child.placement.origin) budget.check(value) ;
      for(const auto &column : columns) for(const auto &value : column) budget.check(value) ;
      Fraction determinant = columns[0][0] ;
      if(element.dimension == 2)
        determinant = columns[0][0]*columns[1][1]-columns[0][1]*columns[1][0] ;
      if(element.dimension == 3)
        determinant = columns[0][0]*(columns[1][1]*columns[2][2]-columns[1][2]*columns[2][1])
                    - columns[0][1]*(columns[1][0]*columns[2][2]-columns[1][2]*columns[2][0])
                    + columns[0][2]*(columns[1][0]*columns[2][1]-columns[1][1]*columns[2][0]) ;
      if(determinant <= Fraction(0))
        throw std::invalid_argument("accepted refinement child has nonpositive orientation") ;
      Vertices actual ; actual.reserve(element.vertexCount) ;
      for(unsigned v = 0; v < element.vertexCount; ++v) {
        auto point = child.placement.apply(referenceVertex(shape,v)) ;
        for(const auto &value : point) budget.check(value) ;
        actual.push_back(std::move(point)) ;
      }
      std::sort(actual.begin(),actual.end()) ; std::sort(expected.begin(),expected.end()) ;
      if(actual != expected)
        throw std::invalid_argument("accepted refinement vertices differ from request") ;
    }

    Pattern triangles(const TriangleSplit &request, ExactBudget &budget) {
      const ReferencePoint a = {{0,0,0}}, b = {{1,0,0}}, c = {{0,1,0}} ;
      switch(request.kind) {
      case TriangleSplit::Kind::none: return {{a,b,c}} ;
      case TriangleSplit::Kind::uniform: {
        const auto ab = midpoint(a,b), bc = midpoint(b,c), ca = midpoint(c,a) ;
        return {{a,ab,ca},{ab,b,bc},{ca,bc,c},{ab,bc,ca}} ;
      }
      case TriangleSplit::Kind::bisect: {
        budget.check(request.position) ; budget.operations(12) ;
        if(request.edge >= 3 || request.position <= Fraction(0) || request.position >= Fraction(1))
          throw std::invalid_argument("invalid triangle bisection") ;
        const std::array<ReferencePoint,3> vertices = {{a,b,c}} ;
        const auto &first = vertices[request.edge], &second = vertices[(request.edge+1)%3] ;
        const auto &opposite = vertices[(request.edge+2)%3] ;
        ReferencePoint middle ;
        for(unsigned d = 0; d < 3; ++d)
          middle[d] = first[d]+request.position*(second[d]-first[d]) ;
        return {{first,middle,opposite},{middle,second,opposite}} ;
      }
      }
      throw std::invalid_argument("unknown triangle split") ;
    }

    Pattern tetrahedra(const TetSplit &request, ExactBudget &budget) {
      const auto &element = referenceElement(Shape::tet) ;
      Vertices vertices ;
      for(unsigned v = 0; v < 4; ++v) vertices.push_back(referenceVertex(Shape::tet,v)) ;
      switch(request.kind) {
      case TetSplit::Kind::none: return {vertices} ;
      case TetSplit::Kind::bisect: {
        budget.check(request.position) ; budget.operations(12) ;
        if(request.edge >= 6 || request.position <= Fraction(0) || request.position >= Fraction(1))
          throw std::invalid_argument("invalid tetrahedron bisection") ;
        const auto &edge = element.edges[request.edge] ;
        ReferencePoint middle ;
        for(unsigned d = 0; d < 3; ++d)
          middle[d] = vertices[edge[0]][d]+request.position*(vertices[edge[1]][d]-vertices[edge[0]][d]) ;
        auto first = vertices, second = vertices ;
        first[edge[1]] = middle ; second[edge[0]] = middle ;
        return {first,second} ;
      }
      case TetSplit::Kind::uniform: {
        if(request.diagonal >= 3) throw std::invalid_argument("invalid tetrahedron diagonal") ;
        Pattern result ;
        for(unsigned corner = 0; corner < 4; ++corner) {
          Vertices local ;
          for(const auto &vertex : vertices) local.push_back(midpoint(vertices[corner],vertex)) ;
          result.push_back(std::move(local)) ;
        }
        std::array<ReferencePoint,6> middle ;
        for(unsigned e = 0; e < 6; ++e)
          middle[e] = midpoint(vertices[element.edges[e][0]],vertices[element.edges[e][1]]) ;
        const unsigned diagonal[3][2] = {{0,5},{1,3},{2,4}} ;
        const unsigned ring[3][4] = {{2,3,4,1},{0,4,5,2},{0,1,5,3}} ;
        const unsigned choice = request.diagonal ;
        for(unsigned i = 0; i < 4; ++i)
          result.push_back({middle[diagonal[choice][0]],middle[diagonal[choice][1]],
                            middle[ring[choice][i]],middle[ring[choice][(i+1)%4]]}) ;
        return result ;
      }
      }
      throw std::invalid_argument("unknown tetrahedron split") ;
    }

    void tensor(Shape parent, const std::vector<Child> &children,
                const TensorSplit &request, ExactBudget &budget) {
      if(parent != Shape::segment && parent != Shape::quad && parent != Shape::hex)
        throw std::invalid_argument("tensor split requires a segment, quad or hex") ;
      const auto &element = referenceElement(parent) ;
      std::array<std::size_t,3> counts = {{1,1,1}} ; std::size_t expected = 1 ;
      for(unsigned d = 0; d < 3; ++d) {
        if(d >= element.dimension && !request.cuts[d].empty())
          throw std::invalid_argument("cuts in an unused direction") ;
        cuts(request.cuts[d],budget) ;
        if(request.cuts[d].size() == std::numeric_limits<std::size_t>::max())
          throw std::length_error("refinement interval count exceeds addressable size") ;
        counts[d] = request.cuts[d].size()+1 ; expected = product(expected,counts[d]) ;
      }
      count(children,expected) ;
      for(std::size_t ordinal = 0; ordinal < children.size(); ++ordinal) {
        std::size_t index = ordinal ; std::array<std::array<Fraction,2>,3> intervals ;
        for(unsigned d = 0; d < 3; ++d) {
          intervals[d] = interval(request.cuts[d],index%counts[d]) ; index /= counts[d] ;
        }
        Vertices vertices ;
        for(unsigned v = 0; v < element.vertexCount; ++v) {
          ReferencePoint point = {{0,0,0}} ;
          for(unsigned d = 0; d < element.dimension; ++d)
            point[d] = intervals[d][element.vertices[v][d] == 0 ? 0 : 1] ;
          vertices.push_back(std::move(point)) ;
        }
        childMatches(children[ordinal],parent,std::move(vertices),budget) ;
      }
    }
  }

  void validateSubdivision(Shape parent, const std::vector<Child> &children,
                           const SplitOperation &operation, ExactBudget &budget) {
    if(children.size() < 2) throw std::invalid_argument("refinement requires at least two children") ;
    // Matching every ordinal to the requested native partition certifies coverage
    // without pairwise intersection tests. The accepted positive chart is retained.
    std::visit([&](const auto &request) {
      using Request = std::decay_t<decltype(request)> ;
      if constexpr(std::is_same_v<Request,TensorSplit>) {
        tensor(parent,children,request,budget) ;
      } else if constexpr(std::is_same_v<Request,TriangleSplit>) {
        if(parent != Shape::triangle) throw std::invalid_argument("triangle split requires a triangle") ;
        auto pattern = triangles(request,budget) ; count(children,pattern.size()) ;
        for(std::size_t i = 0; i < children.size(); ++i)
          childMatches(children[i],Shape::triangle,std::move(pattern[i]),budget) ;
      } else if constexpr(std::is_same_v<Request,PrismSplit>) {
        if(parent != Shape::prism) throw std::invalid_argument("prism split requires a prism") ;
        const auto base = triangles(request.base,budget) ; cuts(request.axialCuts,budget) ;
        if(request.axialCuts.size() == std::numeric_limits<std::size_t>::max())
          throw std::length_error("refinement interval count exceeds addressable size") ;
        count(children,product(base.size(),request.axialCuts.size()+1)) ;
        for(std::size_t i = 0; i < children.size(); ++i) {
          const auto axial = interval(request.axialCuts,i/base.size()) ; Vertices vertices ;
          for(const auto &height : axial) for(auto vertex : base[i%base.size()]) {
            vertex[2] = height ; vertices.push_back(std::move(vertex)) ;
          }
          childMatches(children[i],Shape::prism,std::move(vertices),budget) ;
        }
      } else if constexpr(std::is_same_v<Request,TetSplit>) {
        if(parent != Shape::tet) throw std::invalid_argument("tet split requires a tet") ;
        auto pattern = tetrahedra(request,budget) ; count(children,pattern.size()) ;
        for(std::size_t i = 0; i < children.size(); ++i)
          childMatches(children[i],Shape::tet,std::move(pattern[i]),budget) ;
      } else {
        if(parent != Shape::pyramid || request.kind != PyramidSplit::Kind::uniform)
          throw std::invalid_argument("accepted pyramid refinement requires uniform10") ;
        const Fraction h(1,2) ;
        const ReferencePoint a={{0,0,0}}, b={{1,0,0}}, c={{1,1,0}}, d={{0,1,0}}, top={{0,0,1}} ;
        const ReferencePoint ab={{h,0,0}}, bc={{1,h,0}}, cd={{h,1,0}}, da={{0,h,0}} ;
        const ReferencePoint ae={{0,0,h}}, be={{h,0,h}}, ce={{h,h,h}}, de={{0,h,h}}, center={{h,h,0}} ;
        Pattern pattern = {{a,ab,center,da,ae},{ab,b,bc,center,be},{center,bc,c,cd,ce},
                           {da,center,cd,d,de},{ae,be,ce,de,top},{ae,de,ce,be,center},
                           {ab,center,ae,be},{bc,center,be,ce},{cd,center,ce,de},{da,center,de,ae}} ;
        count(children,pattern.size()) ;
        for(std::size_t i = 0; i < children.size(); ++i)
          childMatches(children[i],i < 6 ? Shape::pyramid : Shape::tet,std::move(pattern[i]),budget) ;
      }
    },operation) ;
  }

}
}
}
