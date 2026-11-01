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

#ifndef LOCI_FVMADAPT_COMPOSITE_CELL_H
#define LOCI_FVMADAPT_COMPOSITE_CELL_H

#include "entity_identity.h"
#include "geometry_map.h"

#include <limits>
#include <optional>

namespace Loci {
namespace FVMAdapt {

  constexpr unsigned noCompositeIndex = std::numeric_limits<unsigned>::max() ;

  /// Original outward face loop; vertex indices address the input points.
  struct CompositeFace {
    EntityId id ;
    std::vector<unsigned> vertices ;
  } ;

  /// A piece boundary belongs either to an original face patch or to another piece.
  struct PieceBoundary {
    unsigned sourceFace = noCompositeIndex ;
    unsigned sourcePatch = noCompositeIndex ;
    unsigned neighborPiece = noCompositeIndex ;
    unsigned neighborBoundary = noCompositeIndex ;
  } ;

  /// Native geometry piece; its internal boundaries are not automatically solver faces.
  struct GeometryPiece {
    Shape shape ;
    /// Indices into CompositeCell::vertices, in native vertex order.
    std::vector<unsigned> vertices ;
    std::shared_ptr<const GeometryMap> geometry ;
    /// Correspondence in native boundary order.
    std::vector<PieceBoundary> boundaries ;
  } ;

  /// One source cell represented by connected native maps.
  struct CompositeCell {
    /// Input points followed by one interior point.
    std::vector<Point> vertices ;
    /// Outward patches of each input face; triangles and quads remain intact.
    std::vector<std::vector<std::vector<unsigned> > > facePatches ;
    std::vector<GeometryPiece> pieces ;
  } ;

  /// Triangulate a face whose area-normal projection is convex, retaining boundary vertices.
  /// Result indices address the input loop. Vertex IDs select the same diagonals under
  /// loop rotation or reversal. For warped polygons the triangles define its surface.
  /// Degenerate or nonconvex projected loops throw invalid_argument.
  std::vector<std::array<unsigned,3> > triangulateFace(
    const std::vector<EntityId> &vertexIds, const std::vector<Point> &points) ;

  /// Cone an outward closed spherical shell to a certified interior point.
  /// Triangles produce tets; quads produce pyramids retaining their bilinear surface.
  /// Larger faces use the canonical triangular surface from triangulateFace().
  /// All input points must occur in the shell. IDs must be unique within each entity kind.
  /// The optional point must see every face patch with a positive Jacobian. Otherwise the
  /// rounded vertex average is tried; failure requires another supplied point or decomposition.
  /// Exact predicates reject overlapping radial sectors and invalid incidence. They do not
  /// impose a mesh-quality threshold. Invalid input throws invalid_argument.
  CompositeCell makeCompositeCell(const std::vector<Point> &points,
                                 const std::vector<EntityId> &vertexIds,
                                 const std::vector<CompositeFace> &faces,
                                 std::optional<Point> interiorPoint = std::nullopt) ;

}
}

#endif
