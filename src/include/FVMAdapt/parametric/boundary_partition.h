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

#ifndef LOCI_FVMADAPT_BOUNDARY_PARTITION_H
#define LOCI_FVMADAPT_BOUNDARY_PARTITION_H

#include "subdivision.h"
#include <cstddef>

namespace Loci {
namespace FVMAdapt {

  /// An affine segment, triangle or quad placed in a common reference frame.
  struct BoundaryUse {
    FaceUse face ;
    Shape shape ;
    Placement placement ;
    /// Convert a point in the placed piece to its local coordinates. Outside points or
    /// degenerate placements throw invalid_argument.
    ReferencePoint coordinates(const ReferencePoint &point) const ;
  } ;

  /// One common piece, with canonical ordered vertices and coordinates on both leaf boundaries.
  struct BoundaryMatch {
    std::vector<ReferencePoint> vertices ;
    BoundaryUse first, second ;
  } ;

  /// One convex polygon in common triangle/quad face coordinates.
  struct BoundaryPolygon {
    FaceUse face ;
    std::vector<ReferencePoint> vertices ;
  } ;
  /// Common polygon and its incident leaf faces.
  struct PolygonMatch {
    std::vector<ReferencePoint> vertices ;
    FaceUse first, second ;
  } ;
  /// Reconcile convex covers without requiring an affine map for each polygon.
  /// Both covers must cover the same part of the face; empty intersections are omitted.
  std::vector<PolygonMatch> reconcilePolygons(Shape frame,
      const std::vector<BoundaryPolygon> &,const std::vector<BoundaryPolygon> &,std::size_t maxPieces) ;

  /// Intersect two nonoverlapping affine covers of a segment, triangle or quad frame.
  /// Covers may occupy a subset of the unit domain. Empty covers match only each other.
  /// Each face occurs once per input cover. Polygon loops are counterclockwise, rotated
  /// to their smallest vertex; results are sorted by these canonical vertex lists.
  /// Invalid covers/placements throw invalid_argument; excessive output throws length_error.
  /// Inputs are never changed.
  std::vector<BoundaryMatch> reconcileBoundary(Shape shape,
      const std::vector<BoundaryUse> &first, const std::vector<BoundaryUse> &second,
      std::size_t maxPieces) ;

}
}

#endif
