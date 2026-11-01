// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_FACE_PATCH_H
#define LOCI_FVMADAPT_FACE_PATCH_H
#include <array>
#include <vector>

namespace Loci {
  /// Identifies a face region so subdivisions can be matched between meshes.
  /// Coordinates describe logical refinement, not physical positions.
  struct FacePatch {
    std::vector<long long> surface ; // Opaque surface identity within one base mesh.
    std::array<long long,4> bounds{} ; // xmin, ymin, xmax, ymax
    int sector = 0 ; // Negative values denote an unsplit polygon's sector count.
    int orientation = 1 ;
    /// Convex polygon in a common planar chart; empty for integer rectangles/sectors.
    std::vector<std::array<double,2>> polygon ;
    /// Physical area; zero when unavailable.
    double area = 0 ;
    /// Area integrated on a retained map. Logical patches must be nested.
    bool mappedArea = false ;
  } ;

}
#endif
