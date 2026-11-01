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

#ifndef LOCI_FVMADAPT_NATIVE_CELL_H
#define LOCI_FVMADAPT_NATIVE_CELL_H

#include "entity_identity.h"

namespace Loci {
namespace FVMAdapt {

  /// A native cell ordered from its oriented boundary loops.
  struct NativeCell {
    Shape shape ;
    /// Original vertex IDs in referenceElement(shape) order.
    std::vector<EntityId> vertices ;
    /// Input face index for each native boundary.
    std::vector<unsigned> faces ;
    /// Native boundary coordinates to input loop coordinates, using orientPoint.
    std::vector<unsigned> faceOrientations ;
  } ;

  /// Recognize a tet, pyramid, prism or hex from complete outward face loops.
  /// Selects the lexicographically least native vertex order, preserving loop handedness.
  /// For native boundary f and vertex i, input[faces[f]][orientedVertex(...,i)]
  /// equals vertices[referenceElement(shape).boundaries[f].vertices[i]].
  /// Unsupported shapes or invalid connectivity throw invalid_argument.
  /// Coordinates and geometric validity are not checked.
  NativeCell recognizeNativeCell(const std::vector<std::vector<EntityId> > &outwardFaces) ;

}
}

#endif
