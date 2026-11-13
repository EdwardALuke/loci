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

#ifndef LOCI_FVMADAPT_VERTEX_VOLUME_SPLIT_H
#define LOCI_FVMADAPT_VERTEX_VOLUME_SPLIT_H

#include "subdivision.h"
#include "geometry_quality.h"

#include <cstddef>

namespace Loci {
namespace FVMAdapt {

  /// Quad child at a parent corner, bounded by two edge midpoints and the accepted center.
  CornerMap quadChildMap(unsigned corner, const Point &center) ;

  /// One hex chart in parent reference coordinates; vertices address VertexVolumeSplit::points.
  struct VertexVolumeChild {
    std::array<unsigned,8> vertices ;
    std::shared_ptr<const GeometryMap> placement ;
  } ;

  /// A child face on the quad patch around one corner of a parent boundary.
  /// orientation maps child-face coordinates into that patch's coordinates.
  struct VertexVolumeBoundary {
    FaceUse face ; ///< Child and its local face index.
    unsigned boundary ; ///< Parent boundary index.
    unsigned corner ; ///< Corner index within that parent boundary.
    unsigned orientation ;
  } ;

  /// Vertex-centered hex, prism or tet subdivision, with one hex child per original vertex.
  struct VertexVolumeSplit {
    Shape shape ;
    /// Original corners, edge midpoints, face centers, then the cell center.
    std::vector<Point> points ;
    /// Child ordinal is its original vertex index.
    std::vector<VertexVolumeChild> children ;
    std::vector<InternalFace> internalFaces ;
    std::vector<VertexVolumeBoundary> boundaryPieces ;
  } ;

  /// Construct one hex child per vertex of a hex, prism or tet in reference space.
  /// faceCenters use each native boundary's local coordinates; cellCenter uses the parent coordinates.
  /// Every boundary patch must be strictly convex. Bounds with exact fallback certify positive child
  /// Jacobians; unresolved bounds throw length_error, invalid geometry throws invalid_argument.
  /// Maps are numerical trilinear evaluators. Physical geometry is retained by composing them
  /// with the parent GeometryMap; this function does not test physical mesh quality.
  VertexVolumeSplit splitVertexVolume(Shape shape, const std::vector<Point> &faceCenters,
                                     const Point &cellCenter, JacobianCheckLimits limits = {}) ;

  /// Check center placement and positive child Jacobians without constructing geometry maps.
  void validateVertexVolumeSplit(Shape shape, const std::vector<Point> &faceCenters,
                                 const Point &cellCenter, JacobianCheckLimits limits = {}) ;

  /// Restore accepted points and child corners, validating topology and positive child Jacobians.
  /// Face connections and trilinear evaluators are rebuilt without choosing new centers.
  VertexVolumeSplit restoreVertexVolumeSplit(Shape shape, std::vector<Point> points,
      const std::vector<std::array<unsigned,8>> &children, JacobianCheckLimits limits = {}) ;

}
}

#endif
