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

#ifndef LOCI_FVMADAPT_REFERENCE_ELEMENT_H
#define LOCI_FVMADAPT_REFERENCE_ELEMENT_H

#include <array>

namespace Loci {
namespace FVMAdapt {

  /// Three coordinates or vector components; unused reference components are zero.
  typedef std::array<double,3> Point ;

  /// Reference shapes; a polygon supplies its vertex count through its geometry map.
  enum class Shape { point, segment, triangle, quad, tet, hex, prism, pyramid, polygon } ;

  /// A codimension-one boundary, ordered as the vertices of its own shape.
  struct Boundary {
    const Shape shape ;
    const unsigned vertexCount ;
    const std::array<int,4> vertices ;
  } ;

  /// Fixed reference topology. Only entries below each count are used.
  /// Polygon boundaries run counterclockwise; volume boundary faces point outward.
  /// Edges have a fixed local direction, independent of boundary orientation.
  struct ReferenceElement {
    const Shape shape ;
    const unsigned dimension ;
    const unsigned vertexCount ;
    const std::array<Point,8> vertices ;
    const unsigned edgeCount ;
    const std::array<std::array<int,2>,12> edges ;
    const unsigned boundaryCount ;
    const std::array<Boundary,6> boundaries ;
  } ;

  /// Unit segment, square and cube; simplex coordinates are nonnegative with
  /// sum <= 1. The prism is triangle x [0,1]. The pyramid has x,y,z >= 0,
  /// x+z <= 1 and y+z <= 1, with apex (0,0,1). Throws invalid_argument for an unknown shape.
  const ReferenceElement &referenceElement(Shape shape) ;

  /// Test the reference-domain inequalities with an absolute coordinate slack.
  /// Nonfinite points return false; invalid shapes or negative/nonfinite tolerances throw invalid_argument.
  bool contains(Shape shape, const Point &point, double tolerance = 0) ;

  /// Embed boundary coordinates in its parent reference element using its corner basis.
  /// Finite coordinates may extrapolate outside the boundary; unused components are ignored.
  /// Invalid shape/coordinates throw invalid_argument; an invalid boundary index throws out_of_range.
  Point boundaryPoint(Shape shape, unsigned boundary, const Point &point) ;

  /// Number of vertex-preserving orientations: 1, 2, 6, 8 for point, segment, triangle, quad.
  /// Volume shapes and unknown shapes throw invalid_argument.
  unsigned orientationCount(Shape shape) ;

  /// Map local vertex i to its vertex index in the oriented reference shape.
  /// For polygons, codes 0..n-1 map i to (i+code)%n; codes n..2*n-1 map i to (code-i)%n.
  /// Segment codes are identity and reversal. Invalid indices throw out_of_range.
  unsigned orientedVertex(Shape shape, unsigned orientation, unsigned vertex) ;
  /// Cyclic face orientation for any vertex count >= 3.
  unsigned orientedVertex(unsigned vertices, unsigned orientation, unsigned vertex) ;

  /// Apply the orientation to coordinates; local vertex i maps to orientedVertex(...,i).
  /// Finite coordinates may extrapolate; unused components are ignored. Nonfinite inputs throw invalid_argument.
  Point orientPoint(Shape shape, unsigned orientation, const Point &point) ;

  /// Return the orientation that undoes the supplied orientation.
  unsigned inverseOrientation(Shape shape, unsigned orientation) ;

}
}

#endif
