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
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>.
//#
//#############################################################################

#ifndef LOCI_FVMADAPT_GEOMETRY_QUALITY_H
#define LOCI_FVMADAPT_GEOMETRY_QUALITY_H

#include "geometry_map.h"

#include <cstddef>

namespace Loci {
namespace FVMAdapt {

  /// Work limits for exact Bernstein bounds on a map's Jacobian.
  struct JacobianCheckLimits {
    unsigned maxDepth = 18 ;
    std::size_t maxBoxes = 4096 ;
  } ;

  /// Sampled success is not a validity certificate; invalid takes precedence over unresolved.
  enum class GeometryQualityStatus { sampled, invalid, unresolved } ;

  /// An empty sample list selects fixed native-shape samples, omitting the pyramid apex.
  struct GeometryQualityOptions {
    std::vector<Point> samples ;
    std::size_t maxSamples = 4096 ;
  } ;

  /// Dimensionless derivative measurements in the map's reference coordinates.
  /// Metrics and their locations describe complete measurements only; inspect measuredCount first.
  /// Values depend on the native reference shape; no ideal-element normalization or thresholds are imposed.
  struct GeometryQuality {
    GeometryQualityStatus status = GeometryQualityStatus::sampled ;
    /// Counts may overlap: an inverted sample can also have complete finite measurements.
    std::size_t sampleCount = 0, measuredCount = 0, invalidCount = 0, unresolvedCount = 0 ;
    /// Signed normalized determinant in 3D, unsigned tangent-area ratio in 2D, one in 1D.
    double minScaledJacobian = 1 ;
    /// Largest absolute cosine between distinct reference tangents; zero in 1D.
    double maxAlignment = 0 ;
    /// Largest ratio of tangent lengths, separate from angular distortion; one in 1D.
    double maxDerivativeRatio = 1 ;
    Point scaledJacobianPoint = {{0,0,0}}, alignmentPoint = {{0,0,0}}, derivativeRatioPoint = {{0,0,0}} ;
    /// First invalid or unresolved sample; meaningful when either count is nonzero.
    Point problemPoint = {{0,0,0}} ;
    bool inverted = false, collapsed = false ;
    /// Nonfinite evaluated positions/derivatives or an exception from the map evaluator.
    bool invalidEvaluation = false ;
  } ;

  /// Evaluate finite samples inside a segment, triangle, quad, tet, prism, pyramid or hex.
  /// Nonfinite geometry, zero tangents, zero tangent area or nonpositive 3D determinants are invalid.
  /// Missing derivatives, evaluation exceptions and unrepresentable metrics are unresolved.
  /// Surface orientation is not tested. Explicit pyramid-apex samples may be unresolved.
  /// Invalid shapes, sample coordinates or work limits throw invalid_argument; allocation failures propagate.
  GeometryQuality geometryQuality(const GeometryMap &geometry, const GeometryQualityOptions &options = {}) ;

  /// Interior rank/orientation proof; boundary singularities are permitted.
  enum class GeometryValidity { valid, invalid, unresolved } ;

  /// Certify positive interior volume Jacobians or full interior surface rank for supported maps.
  /// Exact input doubles define the geometry. Unsupported maps and exhausted bounds are unresolved.
  /// This does not certify boundary incidence, global injectivity or mesh quality.
  GeometryValidity geometryValidity(const GeometryMap &geometry, JacobianCheckLimits limits = {}) ;

  /// Accept sampled regular geometry or a certified positive interior with boundary singularities.
  /// This retains the ordinary sampled check; it is not a global injectivity certificate.
  bool validGeometry(const GeometryMap &geometry, const GeometryQualityOptions &options = {}, JacobianCheckLimits limits = {}) ;

  /// Prove a positive Jacobian throughout a trilinear hex from bounds on its edge derivatives.
  /// False means inconclusive. Eight finite corners in native order are required.
  bool positiveHexJacobian(const std::vector<Point> &corners) ;

}
}

#endif
