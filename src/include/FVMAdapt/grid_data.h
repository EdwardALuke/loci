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

#ifndef LOCI_FVMADAPT_GRID_DATA_H
#define LOCI_FVMADAPT_GRID_DATA_H

#include <Loci>
#include "parametric/entity_identity.h"
#include "parametric/history_types.h"
#include <string>
#include <utility>
#include <vector>

namespace Loci {
namespace FVMAdapt {

  /// A retained source geometry piece and a region in its refinement history.
  struct ParametricRegion {
    /// Stable topology ID, not an index into the imported geometry pieces.
    EntityId pieceId ;
    /// Adapter-specific leaf or fragment handle; see the emitting grid's interface.
    TreeNodeId leaf ;
  } ;

  /// One emitted solver cell and its source geometry regions.
  /// The nodal cell approximates these regions; it need not have their exact volume.
  struct ParametricCell {
    int source ;
    /// Loci cell entity; subtract ParametricGridData::cells.Min() for the file cell number.
    int output ;
    std::vector<ParametricRegion> regions ;
  } ;

  /// Serial nodal mesh ready for writeVOG; volume tags use zero-based file cell numbers.
  struct ParametricGridData {
    store<vector3d<double>> pos ;
    Map cl, cr ;
    multiMap face2node ;
    entitySet cells ;
    std::vector<std::pair<int,std::string>> boundary_ids ;
    std::vector<std::pair<std::string,entitySet>> volTags ;
    std::vector<ParametricCell> correspondence ;
    ParametricGridData() = default ;
    ParametricGridData(const ParametricGridData &) = delete ;
    ParametricGridData &operator=(const ParametricGridData &) = delete ;
    ParametricGridData(ParametricGridData &&other) { *this = std::move(other) ; }
    ParametricGridData &operator=(ParametricGridData &&other) {
      if(this != &other) {
        pos.setRep(other.pos.Rep()) ; cl.setRep(other.cl.Rep()) ;
        cr.setRep(other.cr.Rep()) ; face2node.setRep(other.face2node.Rep()) ;
        cells = std::move(other.cells) ; boundary_ids = std::move(other.boundary_ids) ;
        volTags = std::move(other.volTags) ; correspondence = std::move(other.correspondence) ;
      }
      return *this ;
    }
  } ;

}
}

#endif
