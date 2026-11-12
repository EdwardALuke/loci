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

#include <FVMAdapt/parametric/split_quality.h>

#include <algorithm>
#include <cmath>
#include <new>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {
  namespace {

    struct Sum {
      long double value = 0, correction = 0 ;
      void add(long double term) {
        const long double next = value+term ;
        correction += std::abs(value) >= std::abs(term) ? (value-next)+term : (term-next)+value ;
        value = next ;
      }
      long double total() const { return value+correction ; }
    } ;

    // Distinguish an unavailable evaluator from invalid geometry observed by integration.
    class MeasuredMap final : public GeometryMap {
      const GeometryMap &map_ ;
    public:
      mutable bool unavailable = false ;
      explicit MeasuredMap(const GeometryMap &map) : map_(map) {}
      Shape shape() const override { return map_.shape() ; }
      MapValue evaluate(const Point &point) const override {
        try {
          const auto value = map_.evaluate(point) ;
          unavailable = unavailable || !value.hasJacobian ;
          return value ;
        } catch(const std::bad_alloc &) { throw ; }
        catch(...) { unavailable = true ; throw ; }
      }
    } ;

    std::string label(bool parent, std::size_t child) {
      return parent ? "parent" : "child "+std::to_string(child) ;
    }

    void problem(SplitQuality &report, GeometryQualityStatus status, bool parent, std::size_t child,
                 const std::string &reason) {
      if(report.status == GeometryQualityStatus::sampled ||
         (report.status == GeometryQualityStatus::unresolved && status == GeometryQualityStatus::invalid)) {
        report.status = status ; report.parentProblem = parent ; report.problemChild = parent ? noSplitChild : child ;
        report.reason = label(parent,child)+": "+reason ;
      }
    }

    ChildQuality assess(const GeometryMap &map, const SplitQualityOptions &options, SplitQuality &report,
                        bool parent, std::size_t child) {
      ChildQuality result ;
      result.quality = geometryQuality(map,options.quality) ;
      if(result.quality.status != GeometryQualityStatus::sampled) {
        if(!result.quality.invalidEvaluation) result.validity = geometryValidity(map) ;
        if(result.validity != GeometryValidity::valid) {
          problem(report,result.quality.status,parent,child,result.quality.status == GeometryQualityStatus::invalid ?
            "invalid sampled geometry" : "unresolved geometry quality") ;
          return result ;
        }
      }
      MeasuredMap measured(map) ;
      try {
        result.measure = geometryMeasure(measured,options.measure) ;
        if(!result.measure.converged)
          problem(report,GeometryQualityStatus::unresolved,parent,child,"integration did not converge within its bounds") ;
      } catch(const std::bad_alloc &) { throw ; }
      catch(const std::invalid_argument &error) {
        problem(report,measured.unavailable ? GeometryQualityStatus::unresolved : GeometryQualityStatus::invalid,parent,child,error.what()) ;
      } catch(const std::exception &error) {
        problem(report,GeometryQualityStatus::unresolved,parent,child,error.what()) ;
      }
      return result ;
    }

    double finiteDouble(long double value) {
      if(!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
        throw std::overflow_error("split consistency diagnostic is unrepresentable") ;
      return static_cast<double>(value) ;
    }

    void consistency(SplitQuality &report, unsigned dimension) {
      if(!report.parent.measure.converged) return ;
      Sum weight ;
      std::array<Sum,3> moment ;
      const auto &parent = report.parent.measure ;
      for(const auto &child: report.children) {
        if(!child.measure.converged) return ;
        const long double ratio = static_cast<long double>(child.measure.measure)/parent.measure ;
        weight.add(ratio) ;
        for(unsigned d = 0; d < 3; ++d)
          moment[d].add(ratio*(static_cast<long double>(child.measure.centroid[d])-parent.centroid[d])) ;
      }
      const long double total = weight.total(), scale = std::pow(static_cast<long double>(parent.measure),1.L/dimension) ;
      report.relativeMeasureError = finiteDouble(std::abs(total-1)) ;
      for(unsigned d = 0; d < 3; ++d) report.centroidDifference[d] = finiteDouble(moment[d].total()/total) ;
      report.relativeMomentError = finiteDouble(std::hypot(moment[0].total(),moment[1].total(),moment[2].total())/scale) ;
      report.measuresResolved = true ;
    }

    void requireComplete(const ChildQuality &cell, const std::string &name) {
      const auto &quality = cell.quality ;
      if(cell.validity != GeometryValidity::valid) {
        if(quality.status == GeometryQualityStatus::invalid || quality.invalidCount != 0)
          throw std::domain_error(name+": invalid sampled geometry") ;
        if(quality.status != GeometryQualityStatus::sampled || quality.sampleCount == 0 ||
           quality.measuredCount != quality.sampleCount || quality.unresolvedCount != 0)
          throw std::length_error(name+": incomplete geometry assessment") ;
      }
      if(quality.invalidEvaluation || !std::isfinite(quality.minScaledJacobian) ||
         !std::isfinite(quality.maxAlignment) || !std::isfinite(quality.maxDerivativeRatio))
        throw std::length_error(name+": unrepresentable geometry assessment") ;
      if(!cell.measure.converged)
        throw std::length_error(name+": incomplete geometry assessment") ;
      if(!std::isfinite(cell.measure.measure) || cell.measure.measure <= 0)
        throw std::length_error(name+": unrepresentable geometry assessment") ;
    }

  }

  SplitQuality splitQuality(const GeometryMap &parent, const std::vector<std::shared_ptr<const GeometryMap>> &children,
                            const SplitQualityOptions &options) {
    const unsigned dimension = referenceElement(parent.shape()).dimension ;
    if(dimension == 0 || children.empty() || options.maxChildren == 0 || children.size() > options.maxChildren)
      throw std::invalid_argument("split quality requires a positive-dimensional parent and a bounded nonempty child list") ;
    if(!std::isfinite(options.measure.relativeTolerance) || options.measure.relativeTolerance <= 0 ||
       options.measure.maxOrder < 2 || options.measure.maxOrder > 256 || options.measure.maxEvaluations == 0)
      throw std::invalid_argument("invalid split measure bounds or tolerance") ;
    for(const auto &child: children)
      if(!child || referenceElement(child->shape()).dimension != dimension)
        throw std::invalid_argument("split quality requires child maps with the parent's dimension") ;

    SplitQuality report ;
    report.parent = assess(parent,options,report,true,noSplitChild) ;
    report.children.reserve(children.size()) ;
    for(std::size_t i = 0; i < children.size(); ++i) {
      report.children.push_back(assess(*children[i],options,report,false,i)) ;
      const auto &quality = report.children.back().quality ;
      if(quality.measuredCount == 0) continue ;
      if(report.minScaledJacobianChild == noSplitChild || quality.minScaledJacobian < report.children[report.minScaledJacobianChild].quality.minScaledJacobian)
        report.minScaledJacobianChild = i ;
      if(report.maxAlignmentChild == noSplitChild || quality.maxAlignment > report.children[report.maxAlignmentChild].quality.maxAlignment)
        report.maxAlignmentChild = i ;
      if(report.maxDerivativeRatioChild == noSplitChild || quality.maxDerivativeRatio > report.children[report.maxDerivativeRatioChild].quality.maxDerivativeRatio)
        report.maxDerivativeRatioChild = i ;
    }
    try { consistency(report,dimension) ; }
    catch(const std::overflow_error &error) { problem(report,GeometryQualityStatus::unresolved,true,noSplitChild,error.what()) ; }
    return report ;
  }

  void requireSplitQuality(const SplitQuality &report, const SplitQualityLimits &limits) {
    if(!std::isfinite(limits.minScaledJacobian) || limits.minScaledJacobian < 0 || limits.minScaledJacobian > 1 ||
       !std::isfinite(limits.maxAlignment) || limits.maxAlignment < 0 || limits.maxAlignment > 1 ||
       std::isnan(limits.maxDerivativeRatio) || limits.maxDerivativeRatio < 1 ||
       std::isnan(limits.maxRelativeMeasureError) || limits.maxRelativeMeasureError < 0 ||
       std::isnan(limits.maxRelativeMomentError) || limits.maxRelativeMomentError < 0)
      throw std::invalid_argument("invalid split quality limits") ;
    if(report.status == GeometryQualityStatus::invalid) throw std::domain_error(report.reason) ;
    if(report.status == GeometryQualityStatus::unresolved) throw std::length_error(report.reason) ;
    if(report.children.empty() || !report.measuresResolved)
      throw std::length_error("split quality requires complete child and consistency measurements") ;
    requireComplete(report.parent,"parent") ;
    for(std::size_t i = 0; i < report.children.size(); ++i) {
      const auto &child = report.children[i] ; const auto name = label(false,i) ;
      requireComplete(child,name) ;
      if((limits.minScaledJacobian > 0 || limits.maxAlignment < 1 || std::isfinite(limits.maxDerivativeRatio)) &&
         (child.quality.sampleCount == 0 || child.quality.measuredCount != child.quality.sampleCount || child.quality.unresolvedCount != 0))
        throw std::length_error(name+": quality limits require complete boundary measurements") ;
      if(limits.minScaledJacobian > 0 && child.quality.minScaledJacobian < limits.minScaledJacobian)
        throw std::domain_error(name+": scaled Jacobian is below its limit") ;
      if(child.quality.maxAlignment > limits.maxAlignment)
        throw std::domain_error(name+": tangent alignment exceeds its limit") ;
      if(child.quality.maxDerivativeRatio > limits.maxDerivativeRatio)
        throw std::domain_error(name+": derivative ratio exceeds its limit") ;
    }
    if(!std::isfinite(report.relativeMeasureError) || !std::isfinite(report.relativeMomentError))
      throw std::length_error("split consistency measurements are unresolved") ;
    if(report.relativeMeasureError > limits.maxRelativeMeasureError)
      throw std::domain_error("split measure mismatch exceeds its limit") ;
    if(report.relativeMomentError > limits.maxRelativeMomentError)
      throw std::domain_error("split first-moment mismatch exceeds its limit") ;
  }

}
}
