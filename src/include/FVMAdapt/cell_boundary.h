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

#ifndef LOCI_FVMADAPT_CELL_BOUNDARY_H
#define LOCI_FVMADAPT_CELL_BOUNDARY_H

#include "parametric/boundary_partition.h"
#include "parametric/entity_identity.h"
#include "parametric/affine_subdivision.h"

namespace Loci {
namespace FVMAdapt {

  /// Original topology and a subdivision built through validated local split operations.
  struct CellPartition {
    CellTopology topology ;
    AffineSubdivision subdivision ;
  } ;

  /// A native boundary of one leaf; root indexes the inputs and leaf is its subdivision node ID.
  struct LeafUse {
    unsigned root ;
    TreeNodeId leaf ;
    unsigned boundary ;
  } ;

  /// Face vertices in an incident root's coordinates, in the shared face's vertex order.
  /// reversed requests the opposite loop direction for this leaf's outward boundary.
  struct FaceSide {
    LeafUse cell ;
    std::vector<ReferencePoint> points ;
    bool reversed ;
  } ;

  /// One emitted face piece; one or two incident leaves use the same vertex keys.
  /// A single side means the boundary of the supplied collection, not a physical boundary label.
  struct CellFace {
    std::vector<VertexKey> vertices ;
    std::vector<FaceSide> sides ;
  } ;

  /// A face index and its direction in a leaf's outward boundary.
  struct FaceReference {
    std::size_t face ;
    bool reversed ;
  } ;

  /// Assembled faces for one leaf; indices refer to CellBoundaries::faces.
  struct LeafBoundary {
    unsigned root ;
    TreeNodeId leaf ;
    std::vector<FaceReference> faces ;
  } ;

  /// Limits on output faces and the sum of their vertex counts.
  struct BoundaryLimits {
    std::size_t faces, vertices ;
  } ;

  /// Local connectivity only; geometry maps and original subdivisions remain caller-owned.
  struct CellBoundaries {
    std::vector<CellFace> faces ;
    std::vector<LeafBoundary> leaves ;
  } ;

  /// Reconcile triangle, quad, tet, pyramid, prism and hex subdivisions, including all incident edge points.
  /// Shared entity IDs/connectivity and opposite face orientations must agree. Geometry compatibility
  /// is the caller's responsibility. Unpaired original faces delimit this local collection.
  /// Invalid partitions throw invalid_argument, output limits length_error, exact overflow overflow_error.
  /// Input order changes incidence indices, but not canonical face connectivity. Inputs are unchanged.
  CellBoundaries assembleBoundaries(const std::vector<CellPartition> &cells, BoundaryLimits limits) ;

}
}

#endif
