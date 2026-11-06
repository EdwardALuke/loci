// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include "parametric_hex.h"
#include <FVMAdapt/parametric/geometry_quality.h>
#include "library/parametric/record_io.h"

namespace Loci {
namespace FVMAdapt {
  namespace {
    constexpr std::uint32_t magic = 0x48455843 ;
    void require(bool ok, const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    void bytes(detail::RecordWriter &out, const ParametricRecord &record) { out.u64(record.size()) ; out.append(record) ; }
    ParametricRecord bytes(detail::RecordReader &in, std::size_t limit) {
      ParametricRecord record(in.count(limit)) ; for(auto &b : record) b = in.byte() ; return record ;
    }
    std::vector<BoundaryId> interiorFaces(const ParametricCorner &cell, NodeRange children) {
      std::vector<BoundaryId> result ;
      for(auto child : children) for(auto use : cell.cells.at(child).faces)
        if(cell.boundaries.faces.at(use.surface).parent == noBoundary &&
           std::find(result.begin(),result.end(),use.surface) == result.end()) result.push_back(use.surface) ;
      require(result.size() == 12,"hex construction requires twelve interior faces") ; return result ;
    }
    HexCenters centers(const ParametricCorner &cell, const AcceptedSplit &split) {
      require(split.construction->vertexVolume && split.construction->vertexVolume->shape == Shape::hex,
              "compact hex history requires vertex-based splits") ;
      HexCenters result ; const auto &parent = cell.cells.at(split.parent) ;
      for(unsigned f = 0; f < 6; ++f) {
        const auto use = parent.faces.at(f) ; const auto &face = cell.boundaries.faces.at(use.surface) ;
        require(face.splits.size() == 1,"compact hex history requires one accepted face construction") ;
        const auto p = orientPoint(Shape::quad,inverseOrientation(Shape::quad,use.orientation),face.splits[0].referenceCenter) ;
        result.faces[f] = {{p[0],p[1]}} ;
      }
      result.cell = split.construction->vertexVolume->points.back() ; return result ;
    }
  }

  ParametricHex::ParametricHex(const ParametricCorner &cell) : topology(cell.topology), limits(cell.subdivision.history().limits()) {
    require(topology.shape() == Shape::hex,"compact hex history requires a hex source") ;
    const auto &history = cell.subdivision.history() ;
    const auto *source = std::get_if<CornerRecord>(&history.geometry().record(history.node(0).map)) ;
    require(source && source->shape == Shape::hex && history.geometryRevision() == 0,"compact hex requires original corner geometry") ;
    corners = source->corners ;
    for(unsigned f = 0; f < 6; ++f) faces[f] = cell.faceRecord(topology.faceId(f)) ;
    for(unsigned i = 0; i < history.splits().size(); ++i) {
      const auto &split = history.splits()[i] ;
      require(split.children.size() == 8 && split.children.front() == 1+8*std::size_t(i),"hex history has noncontiguous children") ;
      require(split.parent <= 8*std::size_t(i) && !plan.accepted(split.parent),"hex history has multiple accepted constructions") ;
      plan.refine(split.parent,centers(cell,split),limits.jacobian) ;
      const auto construction = plan.construction(split.parent) ;
      require(construction.points == split.construction->vertexVolume->points,"hex plan does not reproduce accepted split points") ;
      for(unsigned c = 0; c < 8; ++c)
        require(construction.children[c].vertices == split.construction->vertexVolume->children[c].vertices,"hex child order differs") ;
      const auto shared = interiorFaces(cell,split.children) ;
      for(unsigned f = 0; f < shared.size(); ++f) {
        auto face = encodeFacePlan(cell.boundaries,shared[f],cell.boundaryLimits()) ;
        if(!face.empty()) interior.push_back({i,f,std::move(face)}) ;
      }
    }
    for(std::size_t i = plan.splitCount(); i > 0; --i) {
      const auto parent = plan.splitParent(i-1) ;
      if(history.node(parent).isLeaf()) plan.coarsen(parent) ;
    }
    const auto leaves = history.leafIds() ;
    require(plan.leaves() == std::vector<HexPlan::Cell>(leaves.begin(),leaves.end()),"hex selections differ") ;
  }
  ParametricCorner ParametricHex::expand() const {
    ParametricCorner cell(topology,corners,limits) ;
    for(unsigned f = 0; f < 6; ++f) cell.mergeFace(topology.faceId(f),faces[f]) ;
    std::size_t next = 0 ;
    for(std::size_t i = 0; i < plan.splitCount(); ++i) {
      const auto parent = plan.splitParent(i) ;
      const auto expected = plan.centers(parent) ;
      for(unsigned f = 0; f < 6; ++f) {
        const auto use = cell.cells.at(parent).faces[f] ; const auto &face = cell.boundaries.faces.at(use.surface) ;
        require(face.splits.size() == 1,"accepted hex face choice is missing") ;
        const auto p = orientPoint(Shape::quad,inverseOrientation(Shape::quad,use.orientation),face.splits[0].referenceCenter) ;
        require(expected.faces[f] == std::array<double,2>{{p[0],p[1]}},"hex and face plan centers disagree") ;
      }
      cell.restoreHexCell(cell.subdivision,cell.cells,parent,plan) ;
      if(next < interior.size() && interior[next].split == i) {
        const auto shared = interiorFaces(cell,cell.subdivision.children(parent)) ;
        while(next < interior.size() && interior[next].split == i) {
          applyFacePlan(cell,shared.at(interior[next].face),interior[next].plan) ; ++next ;
        }
      }
    }
    require(next == interior.size(),"unused hex interior face plan") ;
    for(std::size_t i = plan.splitCount(); i > 0; --i) {
      const auto parent = plan.splitParent(i-1) ;
      if(!plan.active(parent) || plan.leaf(parent)) cell.subdivision.coarsen(parent) ;
    }
    return cell ;
  }
  ParametricRecord ParametricHex::record() const {
    detail::RecordWriter out(limits.exact.maxBytes,limits.exact) ; out.u32(magic) ; out.u32(1) ;
    const auto t = encodeParametricTopology(topology) ; out.u64(t.size()) ; for(int v : t) out.u32(v) ;
    for(const auto &p : corners) for(double x : p) out.real(x) ;
    bytes(out,plan.encode()) ; for(const auto &face : faces) bytes(out,face) ;
    out.u64(interior.size()) ;
    for(const auto &face : interior) { out.u32(face.split) ; out.byte(face.face) ; bytes(out,face.plan) ; }
    return out.take() ;
  }
  bool ParametricHex::saved(const ParametricRecord &record) {
    if(record.size() < 4) return false ;
    unsigned tag = 0 ; for(unsigned b = 0; b < 4; ++b) tag |= unsigned(record[b])<<(8*b) ; return tag == magic ;
  }
  ParametricHex ParametricHex::decode(const ParametricRecord &record, HistoryLimits limits) {
    detail::RecordReader in(record,limits.exact.maxBytes,limits.exact) ;
    require(in.u32() == magic && in.u32() == 1,"unsupported compact hex record") ;
    std::vector<int> t(in.count(128,4)) ; for(auto &v : t) v = in.u32() ;
    const auto topology = decodeParametricTopology(t) ; require(topology.shape() == Shape::hex,"compact source is not a hex") ;
    std::vector<Point> corners(8) ; for(auto &p : corners) for(auto &x : p) x = in.real() ;
    require(validParametricSource(CornerMap(Shape::hex,corners),limits.jacobian),"invalid compact hex source geometry") ;
    ParametricHex result(topology,std::move(corners),limits) ;
    require(limits.maxNodes > 0,"compact hex requires a positive node limit") ;
    result.plan = HexPlan::decode(bytes(in,limits.exact.maxBytes),std::min(limits.maxSplits,(limits.maxNodes-1)/8),limits.jacobian) ;
    for(std::size_t i = 0; i < result.plan.splitCount(); ++i)
      if(result.plan.depth(result.plan.splitParent(i)) >= limits.maxDepth) throw std::length_error("compact hex depth limit exceeded") ;
    for(auto &face : result.faces) face = bytes(in,limits.exact.maxBytes) ;
    const auto count = in.count(12*result.plan.splitCount(),13) ;
    for(std::size_t i = 0; i < count; ++i) {
      const auto split = in.u32() ; const unsigned face = in.byte() ;
      require(split < result.plan.splitCount() && face < 12,"invalid hex interior face") ;
      if(!result.interior.empty()) {
        const auto &last = result.interior.back() ;
        require(std::make_pair(last.split,last.face) < std::make_pair(split,face),"hex interior faces are not ordered") ;
      }
      auto plan = bytes(in,limits.exact.maxBytes) ; require(!plan.empty(),"empty hex interior face plan") ;
      result.interior.push_back({split,face,std::move(plan)}) ;
    }
    in.finish() ;
    result.validate() ;
    return result ;
  }
}
}
