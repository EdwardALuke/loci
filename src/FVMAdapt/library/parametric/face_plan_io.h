// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_FACE_PLAN_IO_H
#define LOCI_FVMADAPT_FACE_PLAN_IO_H
#include <FVMAdapt/parametric/face_plan.h>
#include "record_io.h"

namespace Loci {
namespace FVMAdapt {
namespace detail {
  constexpr unsigned facePlanMagic = 0x46504c4e ;
  constexpr unsigned facePlanVersion = 2 ;

  inline void requireFacePlan(bool ok, const char *message) {
    if(!ok) throw std::invalid_argument(message) ;
  }

  // Every face and edge has an explicit leaf marker. Children follow corner/endpoint order.
  struct FacePlanBudget {
    const BoundaryRecordLimits &limits ;
    std::size_t nodes = 0, splits = 0 ;

    /// Count one visited face or edge and enforce the nesting limit.
    void visit(unsigned depth) {
      if(depth > limits.maxDepth || depth > 1024) throw std::length_error("face plan depth limit exceeded") ;
      if(++nodes > limits.maxEntities) throw std::length_error("face plan node limit exceeded") ;
    }

    /// Count accepted divisions against the split limit.
    void division(bool split) {
      if(split && ++splits > limits.maxSplits) throw std::length_error("face plan split limit exceeded") ;
    }
  } ;

  /// Return the default split point in reference-face coordinates.
  inline Point midpoint(unsigned vertices) { return vertices == 3 ? Point{{1./3.,1./3.,0}} : Point{{.5,.5,0}} ; }

  /// A face leaf marker or an accepted split point and center method.
  struct FaceDivision {
    bool split = false ;
    CenterMethod method = CenterMethod::reference ;
    Point center{{0,0,0}} ;
  } ;

  /// Require an accepted center strictly inside its reference face.
  inline void checkFaceDivision(const FaceDivision &division, unsigned vertices) {
    if(!division.split) return ;
    requireFacePlan(division.method == CenterMethod::reference || division.method == CenterMethod::weighted ||
                    division.method == CenterMethod::integrated,"unknown face plan center method") ;
    const auto &center = division.center ;
    requireFacePlan(std::isfinite(center[0]) && std::isfinite(center[1]) && center[2] == 0 && center[0] > 0 && center[1] > 0 &&
                    (vertices == 3 ? center[0]+center[1] < 1 :
                     vertices == 4 ? center[0] < 1 && center[1] < 1 : insideReferencePolygon(vertices,center)),
                    "face plan center is outside its face") ;
  }

  /// Write a leaf marker or the accepted face-center construction.
  inline void writeDivision(detail::RecordWriter &out, const FaceDivision &division, unsigned vertices) {
    checkFaceDivision(division,vertices) ;
    if(!division.split) {
      out.byte(0) ;
      return ;
    }
    const bool fixed = division.center == midpoint(vertices) ;
    out.byte(fixed ? 1 : 2) ;
    out.byte(unsigned(division.method)) ;
    if(!fixed) {
      out.real(division.center[0]) ;
      out.real(division.center[1]) ;
    }
  }

  /// Write the face-plan format and original vertex count.
  inline void header(detail::RecordWriter &out, unsigned vertices) {
    out.u32(detail::facePlanMagic) ;
    out.u32(facePlanVersion) ;
    out.u32(vertices) ;
  }

  /// Read face and edge divisions with shared size and depth limits.
  struct FacePlanReader {
    detail::RecordReader input ;
    FacePlanBudget budget ;
    unsigned vertices ;

    FacePlanReader(const FacePlan &plan, const BoundaryRecordLimits &limits)
      : input(plan,limits.maxBytes,limits.geometry.exact), budget{limits} {
      requireFacePlan(input.u32() == facePlanMagic && input.u32() == facePlanVersion,"unsupported face plan version") ;
      vertices = input.u32() ;
      requireFacePlan(vertices >= 3,"invalid face plan vertex count") ;
      if(vertices > limits.maxEntities) throw std::length_error("face plan vertex limit exceeded") ;
      requireFacePlan(vertices < input.remaining(),"truncated face plan edges") ;
    }

    /// Read a face marker and its optional reference-space split point.
    FaceDivision face(unsigned count, unsigned depth) {
      budget.visit(depth) ;
      const auto flag = input.byte() ;
      requireFacePlan(flag <= 2,"invalid face plan division") ;
      budget.division(flag != 0) ;
      if(!flag) return {} ;

      FaceDivision division{true,CenterMethod(input.byte()),midpoint(count)} ;
      if(flag == 2) {
        division.center[0] = input.real() ;
        division.center[1] = input.real() ;
      }
      checkFaceDivision(division,count) ;
      return division ;
    }

    /// Read whether the next edge is bisected.
    bool edge(unsigned depth) {
      budget.visit(depth) ;
      const auto split = input.byte() ;
      requireFacePlan(split <= 1,"invalid face plan edge division") ;
      budget.division(split) ;
      return split ;
    }
  } ;
}
}
}
#endif
