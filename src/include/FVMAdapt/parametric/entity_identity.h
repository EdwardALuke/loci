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

#ifndef LOCI_FVMADAPT_ENTITY_IDENTITY_H
#define LOCI_FVMADAPT_ENTITY_IDENTITY_H

#include "reference_coordinate.h"

#include <vector>
#include <memory>
#include <optional>

namespace Loci {
namespace FVMAdapt {

  /// Original entity ID supplied by the caller, stable across renumbering and migration.
  typedef std::uint64_t EntityId ;
  enum class EntityKind { vertex, edge, face, cell } ;

  /// Extra refinement for divided faces: none, more than half, or also a disjoint pair.
  enum class FaceBalance { none, majority, opposite } ;
  /// Face loops use cell-local vertex numbers. In 2D, supply the polygon's edges.
  /// Only active subdivisions count; saved but unused cuts must be excluded.
  bool faceBalanceRequired(const std::vector<std::vector<unsigned>> &faces,
      const std::vector<unsigned> &divided, FaceBalance policy) ;

  /// A vertex's smallest containing original entity and canonical coordinates within it.
  /// Original vertices have zero coordinates; unused coordinates are zero.
  struct VertexKey {
    EntityKind kind ;
    EntityId entity ;
    ReferencePoint point ;
  } ;

  bool operator==(const VertexKey &a, const VertexKey &b) ;
  bool operator!=(const VertexKey &a, const VertexKey &b) ;
  bool operator<(const VertexKey &a, const VertexKey &b) ;

  /// Local incidence and orientations of one original cell; contains no physical geometry.
  /// Shared IDs must name the same entities and vertex connectivity in every incident cell.
  class CellTopology {
    std::optional<Shape> shape_ ;
    EntityId id_ ;
    std::vector<EntityId> vertices_, edges_, faces_ ;
    std::vector<unsigned> edgeOrientations_, faceOrientations_ ;
    struct Incidence ;
    std::shared_ptr<const Incidence> incidence_ ;
  public:
    /// IDs follow referenceElement(shape) order. Supply edges for dimension >= 2 and
    /// faces for dimension 3; otherwise those lists are empty. Counts/duplicates throw invalid_argument.
    CellTopology(Shape shape, EntityId cell, std::vector<EntityId> vertices,
                 std::vector<EntityId> edges = {}, std::vector<EntityId> faces = {}) ;
    /// A closed polygon shell. Edge endpoints and outward face loops index vertices.
    CellTopology(EntityId cell, std::vector<EntityId> vertices, std::vector<EntityId> edges,
                 std::vector<EntityId> faces, std::vector<std::array<unsigned,2>> edgeVertices,
                 std::vector<std::vector<unsigned>> faceVertices) ;
    /// Same original IDs, local ordering and connectivity.
    bool operator==(const CellTopology &) const ;
    bool operator!=(const CellTopology &other) const { return !(*this==other) ; }
    bool isStandard() const { return shape_.has_value() ; }
    /// Reference shape; throws for cells defined only by their bounding faces.
    Shape shape() const ;
    EntityId id() const { return id_ ; }
    std::size_t vertexCount() const { return vertices_.size() ; }
    std::size_t edgeCount() const { return edges_.size() ; }
    std::size_t faceCount() const { return faces_.size() ; }
    std::array<unsigned,2> edgeVertices(unsigned edge) const ;
    std::vector<unsigned> faceVertices(unsigned face) const ;
    Shape faceShape(unsigned face) const ;
    /// Original IDs in local topology order; invalid indices throw out_of_range.
    EntityId vertexId(unsigned vertex) const { return vertices_.at(vertex) ; }
    EntityId edgeId(unsigned edge) const { return edges_.at(edge) ; }
    EntityId faceId(unsigned face) const { return faces_.at(face) ; }
    /// Map local edge coordinates to increasing endpoint-ID order; invalid indices throw out_of_range.
    unsigned edgeOrientation(unsigned edge) const ;
    /// Map local face coordinates to the lexicographically least vertex-ID loop among its orientations.
    /// Invalid indices throw out_of_range.
    unsigned faceOrientation(unsigned face) const ;
    /// Reduce to an original vertex, edge or face when possible; otherwise use cell coordinates.
    /// Requires a standard cell. Outside points throw invalid_argument.
    VertexKey vertexKey(const ReferencePoint &point) const ;
    /// Locate a canonical key in this cell; return false unchanged if its entity is not incident.
    /// Requires a standard cell. Invalid/noncanonical coordinates throw invalid_argument.
    bool locateVertex(const VertexKey &key, ReferencePoint &point) const ;
  } ;

}
}

#endif
