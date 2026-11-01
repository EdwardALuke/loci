// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_SPACING_SOURCE_H
#define LOCI_FVMADAPT_SPACING_SOURCE_H

#include "parametric/reference_element.h"
#include <string>
#include <vector>

namespace Loci {
namespace FVMAdapt {
  /// Spacing around a segment; lengths use the grid's physical units.
  struct SpacingSource {
    Point start, end ; ///< Segment endpoints; coincident endpoints define a point source.
    double innerRadius, innerSpacing ; ///< Inner cylinder; excludes the segment's end planes.
    double outerRadius, outerSpacing ; ///< Outer capsule, including the endpoints.
    double exponent ; ///< Nonnegative power-law growth beyond outerRadius.
  } ;
  /// Read a positive count followed by endpoints, r0, s0, r1, s1 and exponent per source.
  /// Throws on missing, incomplete or invalid input.
  std::vector<SpacingSource> readSpacingSources(const std::string &filename) ;
  /// Evaluate a validated source at a physical point.
  double sourceSpacing(const Point &, const SpacingSource &) ;
  /// Minimum source spacing at own vertices and their average position.
  double minimumSourceSpacing(const std::vector<Point> &, const std::vector<SpacingSource> &) ;
  /// Sample spacing and resolve nearby source features using the shortest own edge.
  bool refineForSpacing(const std::vector<Point> &, const std::vector<SpacingSource> &, double shortestEdge) ;
}
}
#endif
