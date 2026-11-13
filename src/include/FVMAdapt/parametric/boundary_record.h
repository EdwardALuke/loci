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

#ifndef LOCI_FVMADAPT_BOUNDARY_RECORD_H
#define LOCI_FVMADAPT_BOUNDARY_RECORD_H

#include "geometry_center.h"
#include "geometry_record.h"

namespace Loci {
namespace FVMAdapt {

  /// Local index in a boundary point, edge, or face table.
  using BoundaryId = std::size_t ;
  constexpr BoundaryId noBoundary = std::numeric_limits<BoundaryId>::max() ;

  /// An accepted interior point in its parent edge coordinates.
  struct EdgeCut {
    Fraction position ;
    BoundaryId vertex ;
  } ;

  /// A child edge covers [begin,end] in its parent's coordinates.
  struct EdgeInterval {
    Fraction begin, end ;
    BoundaryId edge ;
  } ;

  /// Accepted edge cuts; the endpoint order matches the stored geometry.
  struct EdgeRecord {
    std::array<BoundaryId,2> vertices{{noBoundary,noBoundary}} ;
    MapId geometry = noBoundary ;
    BoundaryId parent = noBoundary ;
    std::vector<EdgeCut> cuts ;
    std::vector<EdgeInterval> children ;
  } ;

  /// One accepted corner construction, with one quad child per parent vertex.
  /// Children use [corner, next midpoint, center, previous midpoint] order.
  struct FaceSplitRecord {
    CenterMethod centerMethod = CenterMethod::weighted ;
    Point referenceCenter{{0,0,0}} ;
    BoundaryId centerVertex = noBoundary ;
    std::vector<BoundaryId> children ;
  } ;

  /// Faces retain alternatives; current cell uses select the active one.
  struct FaceRecord {
    BoundaryId parent = noBoundary ;
    /// Index in the parent face's splits; noBoundary for a root face.
    std::size_t parentSplit = noBoundary ;
    MapId geometry = noBoundary ;
    std::vector<BoundaryId> vertices, edges ;
    std::vector<FaceSplitRecord> splits ;
    /// Caller-supplied boundary label, inherited by child faces.
    int boundary = 0 ;
  } ;

  /// Shared local geometry and connectivity, independent of cell ownership.
  struct BoundaryRecords {
    GeometryStore geometry ;
    std::vector<Point> points ;
    std::vector<EdgeRecord> edges ;
    std::vector<FaceRecord> faces ;
  } ;

  /// Bounds on all points/edges/faces, accepted constructions, and ancestry depth.
  /// maxBytes bounds the encoded snapshot; geometry.exact.maxBytes bounds retained payloads.
  struct BoundaryRecordLimits {
    std::size_t maxEntities = 10000000 ;
    std::size_t maxSplits = 1000000 ;
    std::size_t maxDepth = 1024 ;
    std::size_t maxBytes = 64*1024*1024 ;
    GeometryLimits geometry ;
  } ;

  /// Check IDs, parent/child links, edge partitions, face loops, and caller limits.
  /// Parent IDs precede their children; unused constructions remain valid records.
  void validateBoundaries(const BoundaryRecords &records, const BoundaryRecordLimits &limits = {}) ;

  /// Preserve point/edge/face IDs and accepted maps; local map IDs may change on decode.
  std::vector<unsigned char> encodeBoundaries(const BoundaryRecords &records,
                                             const BoundaryRecordLimits &limits = {}) ;
  /// Also return the geometry-ID remapping used by the encoded tables.
  std::vector<unsigned char> encodeBoundaries(const BoundaryRecords &records, std::vector<MapId> &mapIds,
                                             const BoundaryRecordLimits &limits = {}) ;
  /// Decode into temporary state; malformed records throw invalid_argument, limits throw length_error.
  BoundaryRecords decodeBoundaries(const std::vector<unsigned char> &record,
                                   const BoundaryRecordLimits &limits = {}) ;

}
}

#endif
