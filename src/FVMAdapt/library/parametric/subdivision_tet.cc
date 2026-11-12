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
#include "subdivision_internal.h"
#include "exact_geometry.h"

#include <algorithm>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {

  Child detail::tetChild(std::array<ReferencePoint,4> vertices) {
    std::array<ExactPoint,3> columns ;
    for(unsigned i = 0; i < 3; ++i)
      for(unsigned d = 0; d < 3; ++d)
        columns[i][d] = exact(vertices[i+1][d])-exact(vertices[0][d]) ;
    const auto &a = columns[0], &b = columns[1], &c = columns[2] ;
    const Exact determinant = a[0]*(b[1]*c[2]-b[2]*c[1])
                            -a[1]*(b[0]*c[2]-b[2]*c[0])
                            +a[2]*(b[0]*c[1]-b[1]*c[0]) ;
    if(determinant == 0) throw std::logic_error("degenerate tetrahedron split template") ;
    if(determinant < 0) std::swap(columns[1],columns[2]) ;
    Child child ; child.shape = Shape::tet ; child.placement.origin = vertices[0] ;
    for(unsigned i = 0; i < 3; ++i)
      for(unsigned d = 0; d < 3; ++d) child.placement.columns[i][d] = Fraction(columns[i][d]) ;
    return child ;
  }

  TetSplit TetSplit::bisectEdge(unsigned edge, Fraction position) {
    TetSplit request ; request.kind = Kind::bisect ; request.edge = edge ; request.position = position ;
    return request ;
  }

  TetSplit TetSplit::uniform8(unsigned diagonal) {
    TetSplit request ; request.kind = Kind::uniform ; request.diagonal = diagonal ; return request ;
  }

  std::vector<Child> detail::splitChildren(const TetSplit &request) {
    std::vector<Child> children ;
    const ReferenceElement &element = referenceElement(Shape::tet) ;
    std::array<ReferencePoint,4> vertices ;
    for(unsigned v = 0; v < 4; ++v) vertices[v] = referenceVertex(Shape::tet,v) ;
    switch(request.kind) {
    case TetSplit::Kind::none:
      children.push_back({Shape::tet,Placement()}) ; break ;
    case TetSplit::Kind::bisect: {
      if(request.edge >= 6 || request.position <= Fraction(0) || request.position >= Fraction(1))
        throw std::invalid_argument("tetrahedron bisection requires an edge in [0,6) and position in (0,1)") ;
      const auto &edge = element.edges[request.edge] ;
      ReferencePoint middle ;
      for(unsigned d = 0; d < 3; ++d)
        middle[d] = Fraction(exact(vertices[edge[0]][d])+exact(request.position)*
                           (exact(vertices[edge[1]][d])-exact(vertices[edge[0]][d]))) ;
      auto first = vertices, second = vertices ;
      first[edge[1]] = middle ; second[edge[0]] = middle ;
      children.push_back(tetChild(first)) ; children.push_back(tetChild(second)) ;
      break ;
    }
    case TetSplit::Kind::uniform: {
      if(request.diagonal >= 3)
        throw std::invalid_argument("uniform tetrahedron split requires a diagonal in [0,3)") ;
      std::array<ReferencePoint,6> midpoints ;
      for(unsigned e = 0; e < 6; ++e)
        for(unsigned d = 0; d < 3; ++d)
          midpoints[e][d] = Fraction((exact(vertices[element.edges[e][0]][d])+
                                    exact(vertices[element.edges[e][1]][d]))/2) ;
      for(unsigned corner = 0; corner < 4; ++corner) {
        auto local = vertices ;
        for(unsigned v = 0; v < 4; ++v)
          for(unsigned d = 0; d < 3; ++d)
            local[v][d] = Fraction((exact(vertices[corner][d])+exact(vertices[v][d]))/2) ;
        children.push_back(tetChild(local)) ;
      }
      // Opposite midpoint pairs and the four edge-connected vertices around each diagonal.
      const unsigned diagonal[3][2] = {{0,5},{1,3},{2,4}} ;
      const unsigned ring[3][4] = {{2,3,4,1},{0,4,5,2},{0,1,5,3}} ;
      const unsigned choice = request.diagonal ;
      const unsigned order[3][4] = {{0,1,2,3},{2,0,1,3},{1,2,0,3}} ;
      for(unsigned i = 0; i < 4; ++i) {
        const Child central = tetChild({{midpoints[diagonal[choice][0]],midpoints[diagonal[choice][1]],
                                          midpoints[ring[choice][i]],midpoints[ring[choice][(i+1)%4]]}}) ;
        // Repeated use of this local diagonal retains the central children's refinement pattern.
        std::array<ReferencePoint,4> local ;
        for(unsigned v = 0; v < 4; ++v)
          local[v] = central.placement.apply(referenceVertex(Shape::tet,order[choice][v])) ;
        children.push_back(tetChild(local)) ;
      }
      break ;
    }
    default: throw std::invalid_argument("unknown tetrahedron split pattern") ;
    }
    return children ;
  }

  SplitResult split(const TetSplit &request) {
    SplitResult result ; result.shape = Shape::tet ;
    result.children = detail::splitChildren(request) ;
    detail::completeSubdivision(result) ; return result ;
  }

  SplitResult split(Shape shape, const TetSplit &request) {
    if(shape != Shape::tet) throw std::invalid_argument("tetrahedron split requires a tetrahedron") ;
    return split(request) ;
  }

}
}
