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

#ifndef LOCI_FVMADAPT_SPLIT_QUALITY_H
#define LOCI_FVMADAPT_SPLIT_QUALITY_H

#include "geometry_measure.h"
#include "geometry_quality.h"

#include <limits>
#include <string>

namespace Loci {
namespace FVMAdapt {

  constexpr std::size_t noSplitChild = std::numeric_limits<std::size_t>::max() ;

  /// Per-map integration and sampling bounds, with a bound on proposed children.
  struct SplitQualityOptions {
    MeasureOptions measure ;
    GeometryQualityOptions quality ;
    std::size_t maxChildren = 1024 ;
  } ;

  /// Physical measure and sampled quality of one retained map.
  struct ChildQuality {
    GeometryMeasure measure ;
    GeometryQuality quality ;
  } ;

  /// Local proposal measurements; no hierarchy or geometry is changed.
  struct SplitQuality {
    GeometryQualityStatus status = GeometryQualityStatus::sampled ;
    ChildQuality parent ;
    std::vector<ChildQuality> children ;
    /// Consistency fields are available only when every measure converged and the fields are representable.
    bool measuresResolved = false ;
    /// Absolute difference between summed child measures and parent measure, divided by parent measure.
    double relativeMeasureError = 0 ;
    /// Measure-weighted child centroid minus the parent centroid, in physical coordinates.
    Point centroidDifference{{0,0,0}} ;
    /// Centered first-moment residual divided by parent measure times parent measure^(1/dimension).
    double relativeMomentError = 0 ;
    /// Indices identify the child's corresponding worst reference point in GeometryQuality.
    std::size_t minScaledJacobianChild = noSplitChild, maxAlignmentChild = noSplitChild, maxDerivativeRatioChild = noSplitChild ;
    /// First highest-severity problem; parentProblem distinguishes a parent failure from a child failure.
    std::size_t problemChild = noSplitChild ;
    bool parentProblem = false ;
    std::string reason ;
  } ;

  /// Assess physical child maps with the parent's dimension. Null maps, invalid options or dimensions throw invalid_argument.
  /// Invalid geometry and unresolved measurements are reported; allocation failures propagate.
  /// Agreement of aggregate measures and moments does not prove coverage, disjointness or valid geometry between samples.
  SplitQuality splitQuality(const GeometryMap &parent, const std::vector<std::shared_ptr<const GeometryMap>> &children,
                            const SplitQualityOptions &options = {}) ;

  /// Optional child-quality and aggregate-consistency limits; defaults impose no anisotropy or consistency bound.
  struct SplitQualityLimits {
    double minScaledJacobian = 0 ;
    double maxAlignment = 1 ;
    double maxDerivativeRatio = std::numeric_limits<double>::infinity() ;
    double maxRelativeMeasureError = std::numeric_limits<double>::infinity() ;
    double maxRelativeMomentError = std::numeric_limits<double>::infinity() ;
  } ;

  /// Require complete sampled quality and converged measures, then apply the supplied child and consistency limits.
  /// Invalid limits throw invalid_argument, rejected geometry/limits domain_error, and unresolved evidence length_error.
  void requireSplitQuality(const SplitQuality &report, const SplitQualityLimits &limits = {}) ;

}
}

#endif
