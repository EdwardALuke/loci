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

#ifndef LOCI_FVMADAPT_SUBDIVISION_H
#define LOCI_FVMADAPT_SUBDIVISION_H

#include "reference_coordinate.h"
#include "geometry_map.h"

#include <variant>
#include <optional>

namespace Loci {
namespace FVMAdapt {

  /// Exact affine placement: origin + sum(columns[i]*q[i]); default is the identity.
  struct Placement {
    ReferencePoint origin ;
    std::array<ReferencePoint,3> columns ;
    Placement() ;
    /// Evaluate any exact point.
    ReferencePoint apply(const ReferencePoint &point) const ;
    /// Compose parent and child placements as parent(child(q)).
    Placement compose(const Placement &child) const ;
    /// Numerical affine evaluator on the supplied shape, with constant derivatives.
    /// Unused coordinates are ignored; unused derivative columns are zero.
    /// Containment and injectivity are the caller's responsibility.
    std::shared_ptr<const GeometryMap> map(Shape shape) const ;
  } ;

  /// A child shape and its placement in the parent's reference coordinates.
  struct Child {
    Shape shape ;
    Placement placement ;
  } ;

  /// Interior cuts in each reference direction. Each list must be strictly increasing.
  /// Empty active lists leave that axis intact; unused directions must have no cuts.
  struct TensorSplit {
    std::array<std::vector<Fraction>,3> cuts ;
  } ;

  /// Triangle refinement in its directed native edge ordering. Default leaves the triangle intact.
  struct TriangleSplit {
    enum class Kind { none, bisect, uniform } ;
    Kind kind ;
    unsigned edge ;
    Fraction position ;
    TriangleSplit() : kind(Kind::none), edge(0), position(1,2) {}
    /// Divide the selected edge at position in (0,1) and join it to the opposite vertex.
    static TriangleSplit bisectEdge(unsigned edge, Fraction position = Fraction(1,2)) ;
    /// Four triangles defined by the three edge midpoints.
    static TriangleSplit uniform4() ;
  } ;

  /// Triangle refinement times strictly increasing interior cuts in the axial direction.
  struct PrismSplit {
    TriangleSplit base ;
    std::vector<Fraction> axialCuts ;
  } ;

  /// Tetrahedral refinement in its directed native edge ordering. Default leaves the tet intact.
  struct TetSplit {
    enum class Kind { none, bisect, uniform } ;
    Kind kind ;
    unsigned edge ;
    Fraction position ;
    unsigned diagonal ;
    TetSplit() : kind(Kind::none), edge(0), position(1,2), diagonal(0) {}
    /// Divide the selected edge at position in (0,1) and join it to the opposite edge.
    static TetSplit bisectEdge(unsigned edge, Fraction position = Fraction(1,2)) ;
    /// Eight midpoint tets. The central diagonal joins midpoint pairs (01,23),
    /// (12,03) or (20,13) for diagonal 0, 1 or 2. No quality policy is implied.
    static TetSplit uniform8(unsigned diagonal) ;
  } ;

  /// Pyramid refinement; uniform refinement gives six pyramids and four tetrahedra.
  struct PyramidSplit {
    enum class Kind { none, uniform } ;
    Kind kind ;
    PyramidSplit() : kind(Kind::none) {}
    /// Four base-corner pyramids, top and inverted middle pyramids, then four side tets.
    static PyramidSplit uniform10() ;
  } ;

  /// Native refinement request; accepted child placements are stored separately.
  typedef std::variant<TensorSplit,TriangleSplit,PrismSplit,TetSplit,PyramidSplit> SplitOperation ;

  /// Midpoint hex/prism cuts, tet uniform8 (diagonal zero), or pyramid uniform10.
  /// Other shapes reject; native patterns have fixed reference centers.
  SplitOperation isotropicSplit(Shape shape) ;

  /// Choose midpoint cuts along long hex directions or across/along a prism.
  /// Use all directions when their mean edge lengths differ by at most aspectRatio.
  /// Suppress cuts whose shortest current edge would halve below minEdgeLength.
  /// Tets and pyramids use isotropic cuts, subject to the same edge cutoff.
  /// Empty means no eligible cut. Uses mapped edge endpoints, not hanging nodes.
  std::optional<SplitOperation> automaticSplit(const GeometryMap &, double aspectRatio = 2, double minEdgeLength = 0) ;

  /// A codimension-one boundary of a child, indexed in referenceElement(child.shape).
  struct FaceUse {
    unsigned child ;
    unsigned boundary ;
  } ;

  /// Paired child boundaries; orientation maps first's local coordinates to second's.
  struct InternalFace {
    FaceUse first ;
    FaceUse second ;
    unsigned orientation ;
  } ;

  /// A child boundary on a parent boundary; vertices use parent-boundary coordinates
  /// and follow the child boundary's vertex order.
  struct BoundaryPiece {
    FaceUse face ;
    unsigned boundary ;
    std::vector<ReferencePoint> vertices ;
  } ;

  /// A child edge on a parent edge; endpoints use the parent's directed edge coordinates.
  struct EdgePiece {
    unsigned child ;
    unsigned edge ;
    unsigned parentEdge ;
    std::array<Fraction,2> endpoints ;
  } ;

  /// Complete local subdivision, with every child boundary paired or assigned to its parent boundary.
  struct SplitResult {
    Shape shape ;
    std::vector<Child> children ;
    std::vector<InternalFace> internalFaces ;
    std::vector<BoundaryPiece> boundaryPieces ;
    std::vector<EdgePiece> edgePieces ;
  } ;

  /// Split a segment, quad or hex at exact reference cuts; children are ordered with x varying fastest.
  /// An empty request returns one identity child. Invalid shapes, cuts or unrepresentable counts
  /// throw invalid_argument. Physical quality is not tested.
  SplitResult split(Shape shape, const TensorSplit &request) ;

  /// Split a triangle; children retain the native triangle shape and positive orientation.
  /// Invalid edge indices, positions or patterns throw invalid_argument.
  SplitResult split(const TriangleSplit &request) ;
  /// Split a prism with base children varying fastest, then increasing axial coordinate.
  /// Invalid cuts or unrepresentable child counts throw invalid_argument.
  SplitResult split(const PrismSplit &request) ;
  /// Shape-checked overloads; mismatched shapes throw invalid_argument.
  SplitResult split(Shape shape, const TriangleSplit &request) ;
  SplitResult split(Shape shape, const PrismSplit &request) ;

  /// Split a tet; uniform children start with the four original corners, then the central ring.
  /// Invalid edges, positions, diagonals or patterns throw invalid_argument.
  SplitResult split(const TetSplit &request) ;
  /// Split a pyramid; all children have positive affine placements in the parent reference cell.
  SplitResult split(const PyramidSplit &request) ;
  SplitResult split(Shape shape, const TetSplit &request) ;
  SplitResult split(Shape shape, const PyramidSplit &request) ;

}
}

#endif
