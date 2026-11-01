// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_FACE_PLAN_IO_H
#define LOCI_FVMADAPT_FACE_PLAN_IO_H
#include <FVMAdapt/parametric/face_plan.h>
#include "record_io.h"

namespace Loci { namespace FVMAdapt { namespace detail {
    constexpr unsigned facePlanMagic = 0x46504c4e ;
    inline void requireFacePlan(bool ok, const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    // Every face and edge has an explicit leaf marker. Children follow corner/endpoint order.
    struct FacePlanBudget {
      const BoundaryRecordLimits &limits ;
      std::size_t nodes = 0, splits = 0 ;
      void visit(unsigned depth) {
        if(depth > limits.maxDepth || depth > 1024) throw std::length_error("face plan depth limit exceeded") ;
        if(++nodes > limits.maxEntities) throw std::length_error("face plan node limit exceeded") ;
      }
      void division(bool split) {
        if(split && ++splits > limits.maxSplits) throw std::length_error("face plan split limit exceeded") ;
      }
    } ;
    inline Point midpoint(unsigned vertices) { return vertices == 3 ? Point{{1./3.,1./3.,0}} : Point{{.5,.5,0}} ; }
    struct FaceDivision {
      bool split = false ;
      CenterMethod method = CenterMethod::reference ;
      Point center{{0,0,0}} ;
    } ;
    inline void checkFaceDivision(const FaceDivision &division, unsigned vertices) {
      if(!division.split) return ;
      requireFacePlan(division.method == CenterMethod::reference || division.method == CenterMethod::weighted ||
              division.method == CenterMethod::integrated,"unknown face plan center method") ;
      const auto &p = division.center ;
      requireFacePlan(std::isfinite(p[0]) && std::isfinite(p[1]) && p[2] == 0 && p[0] > 0 && p[1] > 0 &&
              (vertices == 3 ? p[0]+p[1] < 1 : vertices == 4 ? p[0] < 1 && p[1] < 1:insideReferencePolygon(vertices,p)),"face plan center is outside its face") ;
    }
    inline void writeDivision(detail::RecordWriter &out, const FaceDivision &division, unsigned vertices) {
      checkFaceDivision(division,vertices) ;
      if(!division.split) { out.byte(0) ; return ; }
      const bool fixed = division.center == midpoint(vertices) ;
      out.byte(fixed ? 1 : 2) ; out.byte(unsigned(division.method)) ;
      if(!fixed) { out.real(division.center[0]) ; out.real(division.center[1]) ; }
    }
    inline void header(detail::RecordWriter &out, unsigned vertices) { out.u32(detail::facePlanMagic) ; out.u32(2) ; out.u32(vertices) ; }
    struct FacePlanReader {
      detail::RecordReader input ;
      FacePlanBudget budget ;
      unsigned vertices ;
      FacePlanReader(const FacePlan &plan, const BoundaryRecordLimits &limits) : input(plan,limits.maxBytes,limits.geometry.exact), budget{limits} {
        requireFacePlan(input.u32() == facePlanMagic && input.u32() == 2,"unsupported face plan version") ;
        vertices = input.u32() ; requireFacePlan(vertices>=3,"invalid face plan vertex count") ;
        if(vertices>limits.maxEntities) throw std::length_error("face plan vertex limit exceeded") ;
        requireFacePlan(vertices<input.remaining(),"truncated face plan edges") ;
      }
      FaceDivision face(unsigned count, unsigned depth) {
        budget.visit(depth) ; const auto flag = input.byte() ;
        requireFacePlan(flag <= 2,"invalid face plan division") ;
        budget.division(flag != 0) ;
        if(!flag) return {} ;
        FaceDivision division{true,CenterMethod(input.byte()),midpoint(count)} ;
        if(flag == 2) { division.center[0] = input.real() ; division.center[1] = input.real() ; }
        checkFaceDivision(division,count) ; return division ;
      }
      bool edge(unsigned depth) {
        budget.visit(depth) ; const auto split = input.byte() ; requireFacePlan(split <= 1,"invalid face plan edge division") ;
        budget.division(split) ; return split ;
      }
    } ;
}}}
#endif
