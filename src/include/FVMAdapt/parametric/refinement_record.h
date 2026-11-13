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

#ifndef LOCI_FVMADAPT_REFINEMENT_RECORD_H
#define LOCI_FVMADAPT_REFINEMENT_RECORD_H

#include "geometry_record.h"
#include "refinement_history.h"

namespace Loci {
namespace FVMAdapt {

  /// Bounds applied before allocation and throughout snapshot encoding/decoding.
  struct RecordLimits {
    std::size_t maxBytes = 64*1024*1024 ;
    std::size_t maxNodes = 1000000 ;
    std::size_t maxSplits = 1000000 ;
    std::size_t maxDepth = 1024 ;
    GeometryLimits geometry ;
    JacobianCheckLimits jacobian ;
  } ;

  /// Versioned geometry table; equivalent descriptions have identical bytes.
  /// Map IDs are local handles and may change when records are decoded.
  std::vector<unsigned char> encodeGeometry(const GeometryStore &geometry,
                                           const RecordLimits &limits = {}) ;
  /// Also return the original-to-decoded map ID correspondence.
  std::vector<unsigned char> encodeGeometry(const GeometryStore &geometry,
                                           std::vector<MapId> &encodedIds, const RecordLimits &limits = {}) ;
  GeometryStore decodeGeometry(const std::vector<unsigned char> &record,
                               const RecordLimits &limits = {}) ;

  /// Versioned accepted geometry, constructions and current selections.
  /// Local node/map IDs may change on decode; persistent cell keys do not.
  std::vector<unsigned char> encodeRefinement(const RefinementHistory &history,
                                             const RecordLimits &limits = {}) ;
  /// Also return the original-to-decoded node ID correspondence for accompanying cell data.
  std::vector<unsigned char> encodeRefinement(const RefinementHistory &history,
                                             std::vector<TreeNodeId> &encodedIds, const RecordLimits &limits = {}) ;
  /// Decode into temporary state. Invalid data throws invalid_argument; limits throw length_error.
  RefinementHistory decodeRefinement(const std::vector<unsigned char> &record,
                                     const RecordLimits &limits = {}) ;

}
}

#endif
