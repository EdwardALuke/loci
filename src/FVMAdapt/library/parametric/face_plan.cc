// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/parametric/face_plan.h>
#include <FVMAdapt/parametric/corner_boundary.h>
#include "face_plan_io.h"
#include <optional>

namespace Loci {
namespace FVMAdapt {
  namespace {
    void require(bool ok, const char *message) { if(!ok) throw std::invalid_argument(message) ; }

    using Budget = detail::FacePlanBudget ;
    using Division = detail::FaceDivision ;
    using Reader = detail::FacePlanReader ;
    using detail::midpoint ;
    using detail::checkFaceDivision ;
    using detail::writeDivision ;
    using detail::header ;

    void writeEdge(detail::RecordWriter &out, const BoundaryRecords &b, BoundaryId id, BoundaryId start, Budget &budget, unsigned depth) {
      budget.visit(depth) ; const auto &edge = b.edges.at(id) ;
      require(start == edge.vertices[0] || start == edge.vertices[1],"face plan edge endpoint differs") ;
      budget.division(!edge.cuts.empty()) ; out.byte(!edge.cuts.empty()) ; if(edge.cuts.empty()) return ;
      require(edge.cuts.size() == 1 && edge.children.size() == 2 && edge.cuts[0].position == Fraction(1,2),
              "face plans require midpoint edges") ;
      const bool reverse = start == edge.vertices[1] ;
      writeEdge(out,b,edge.children[reverse ? 1 : 0].edge,start,budget,depth+1) ;
      writeEdge(out,b,edge.children[reverse ? 0 : 1].edge,edge.cuts[0].vertex,budget,depth+1) ;
    }
    void writeFace(detail::RecordWriter &out, const BoundaryRecords &b, BoundaryId id, Budget &budget, unsigned depth) {
      budget.visit(depth) ; const auto &face = b.faces.at(id) ; const auto n = face.vertices.size() ;
      require(n >= 3 && face.edges.size() == n && face.splits.size() <= 1,"face plan requires one polygon construction") ;
      Division division ;
      if(!face.splits.empty()) division = {true,face.splits[0].centerMethod,face.splits[0].referenceCenter} ;
      budget.division(division.split) ;
      writeDivision(out,division,n) ;
      for(unsigned e = 0; e < n; ++e) writeEdge(out,b,face.edges[e],face.vertices[e],budget,depth) ;
      if(division.split) {
        require(face.splits[0].children.size() == n,"face plan child count differs") ;
        for(auto child : face.splits[0].children) writeFace(out,b,child,budget,depth+1) ;
      }
    }
    bool mergeEdge(Reader *a, Reader *b, detail::RecordWriter &out, Budget &budget, unsigned depth) {
      budget.visit(depth) ; const bool left = a && a->edge(depth), right = b && b->edge(depth) ;
      budget.division(left || right) ;
      out.byte(left || right) ;
      if(left || right) for(unsigned child = 0; child < 2; ++child)
        mergeEdge(left ? a : nullptr,right ? b : nullptr,out,budget,depth+1) ;
      return left || right ;
    }
    bool mergeFace(Reader *a, Reader *b, detail::RecordWriter &out, Budget &budget, unsigned n, unsigned depth) {
      budget.visit(depth) ; const auto left = a ? a->face(n,depth) : Division{}, right = b ? b->face(n,depth) : Division{} ;
      auto division = left.split ? left : right ;
      if(left.split && right.split) {
        require(left.center == right.center,"shared face plan centers disagree") ;
        // The accepted point defines the construction; retain a deterministic method label.
        division.method = CenterMethod(std::min(unsigned(left.method),unsigned(right.method))) ;
      }
      budget.division(division.split) ;
      writeDivision(out,division,n) ; bool changed = division.split ;
      for(unsigned e = 0; e < n; ++e) changed = mergeEdge(a,b,out,budget,depth) || changed ;
      if(division.split) for(unsigned child = 0; child < n; ++child)
        mergeFace(left.split ? a : nullptr,right.split ? b : nullptr,out,budget,4,depth+1) ;
      return changed ;
    }
    void applyEdge(Reader &in, CornerBoundary &b, BoundaryId id, BoundaryId start, CornerBoundary::Transaction &edit, unsigned depth) {
      if(!in.edge(depth)) return ;
      const auto middle = b.splitEdge(id,edit) ; const auto edge = b.boundaries.edges[id] ;
      const bool reverse = start == edge.vertices[1] ;
      applyEdge(in,b,edge.children[reverse ? 1 : 0].edge,start,edit,depth+1) ;
      applyEdge(in,b,edge.children[reverse ? 0 : 1].edge,middle,edit,depth+1) ;
    }
    void applyFace(Reader &in, CornerBoundary &b, BoundaryId id, CornerBoundary::Transaction &edit, unsigned depth) {
      const auto face = b.boundaries.faces.at(id) ; const auto n = face.vertices.size() ;
      const auto division = in.face(n,depth) ;
      require(face.splits.size() <= 1,"face plan requires one retained construction") ;
      if(division.split && !face.splits.empty())
        require(face.splits[0].referenceCenter == division.center,"shared face plan centers disagree") ;
      for(unsigned e = 0; e < n; ++e) applyEdge(in,b,face.edges[e],face.vertices[e],edit,depth) ;
      if(!division.split) return ;
      const FaceSplitRecord accepted{division.method,division.center,noBoundary,{}} ;
      const auto selected = b.splitSurface(id,division.method,{},edit,&accepted) ;
      const auto children = b.boundaries.faces[id].splits[selected].children ;
      for(auto child : children) applyFace(in,b,child,edit,depth+1) ;
    }
  }

  FacePlan encodeFacePlan(const BoundaryRecords &b, BoundaryId id, const BoundaryRecordLimits &limits) {
    const auto &face = b.faces.at(id) ; bool divided = !face.splits.empty() ;
    for(auto edge : face.edges) divided = divided || !b.edges.at(edge).cuts.empty() ;
    if(!divided) return {} ;
    detail::RecordWriter out(limits.maxBytes,limits.geometry.exact) ; header(out,face.vertices.size()) ;
    Budget budget{limits} ; writeFace(out,b,id,budget,0) ; return out.take() ;
  }
  FacePlan mergeFacePlans(const FacePlan &a, const FacePlan &b, const BoundaryRecordLimits &limits) {
    if(a.empty() && b.empty()) return {} ;
    std::optional<Reader> left, right ;
    if(!a.empty()) left.emplace(a,limits) ;
    if(!b.empty()) right.emplace(b,limits) ;
    const auto n = left ? left->vertices : right->vertices ;
    require(!right || right->vertices == n,"shared face plan shapes differ") ;
    detail::RecordWriter out(limits.maxBytes,limits.geometry.exact) ; header(out,n) ; Budget budget{limits} ;
    const bool divided = mergeFace(left ? &*left : nullptr,right ? &*right : nullptr,out,budget,n,0) ;
    if(left) left->input.finish() ;
    if(right) right->input.finish() ;
    return divided ? out.take() : FacePlan{} ;
  }
  bool applyFacePlan(CornerBoundary &b, BoundaryId id, const FacePlan &plan) {
    if(plan.empty()) return false ;
    const auto limits = b.boundaryLimits() ; Reader in(plan,limits) ;
    require(b.boundaries.faces.at(id).vertices.size() == in.vertices,"shared face plan shapes differ") ;
    const auto points = b.boundaries.points.size() ; CornerBoundary::Transaction edit(b) ;
    applyFace(in,b,id,edit,0) ; in.input.finish() ; edit.committed = true ;
    return b.boundaries.points.size() != points ;
  }
}
}
