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

#include <FVMAdapt/parametric/refinement_record.h>

#include "record_io.h"

#include <algorithm>
#include <map>
#include <type_traits>

namespace Loci {
namespace FVMAdapt {
namespace {

  using detail::RecordReader ;
  using detail::RecordWriter ;
  // Little-endian headers: "FVAG" for geometry and "FVAR" for refinement.
  constexpr std::uint32_t geometryMagic=0x47415646, historyMagic=0x52415646 ;
  constexpr std::uint32_t geometryVersion=1, historyVersion=7 ;
  constexpr std::size_t minFractionBytes=19 ; // Zero numerator and unit denominator.
  constexpr std::size_t minChildBytes=1+12*minFractionBytes ;
  constexpr std::size_t minNodeBytes=28, minSplitBytes=30 ;
  constexpr std::uint64_t absent=std::numeric_limits<std::uint64_t>::max() ;

  unsigned char shapeTag(Shape shape) {
    switch(shape) {
    case Shape::point: return 1 ;
    case Shape::segment: return 2 ;
    case Shape::triangle: return 3 ;
    case Shape::quad: return 4 ;
    case Shape::tet: return 5 ;
    case Shape::hex: return 6 ;
    case Shape::prism: return 7 ;
    case Shape::pyramid: return 8 ;
    case Shape::polygon: return 9 ;
    }
    throw std::invalid_argument("unknown snapshot shape") ;
  }

  Shape readShape(RecordReader &in) {
    switch(in.byte()) {
    case 1: return Shape::point ;
    case 2: return Shape::segment ;
    case 3: return Shape::triangle ;
    case 4: return Shape::quad ;
    case 5: return Shape::tet ;
    case 6: return Shape::hex ;
    case 7: return Shape::prism ;
    case 8: return Shape::pyramid ;
    case 9: return Shape::polygon ;
    }
    throw std::invalid_argument("unknown snapshot shape tag") ;
  }

  void writePlacement(RecordWriter &out, const Placement &placement) {
    for(const auto &value:placement.origin) out.fraction(value) ;
    for(const auto &column:placement.columns)
      for(const auto &value:column) out.fraction(value) ;
  }

  Placement readPlacement(RecordReader &in) {
    Placement result ;
    for(auto &value:result.origin) value=in.fraction() ;
    for(auto &column:result.columns)
      for(auto &value:column) value=in.fraction() ;
    return result ;
  }

  void writePoints(RecordWriter &out, const std::vector<Point> &points) {
    out.u64(points.size()) ;
    for(const auto &point:points) for(double value:point) out.real(value) ;
  }

  std::vector<Point> readPoints(RecordReader &in, std::size_t limit) {
    std::vector<Point> result(in.count(limit,24)) ;
    for(auto &point:result) for(double &value:point) value=in.real() ;
    return result ;
  }

  void writeMap(RecordWriter &out, const GeometryRecord &record,
                const std::vector<MapId> &ids) {
    std::visit([&](const auto &value) {
      using T=std::decay_t<decltype(value)> ;
      if constexpr(std::is_same_v<T,CornerRecord>) {
        out.byte(1) ; out.byte(shapeTag(value.shape)) ; writePoints(out,value.corners) ;
      } else if constexpr(std::is_same_v<T,TensorBezierRecord>) {
        out.byte(2) ; out.byte(shapeTag(value.shape)) ;
        for(unsigned degree:value.degree) out.u32(degree) ;
        writePoints(out,value.controls) ;
      } else if constexpr(std::is_same_v<T,AffineRecord>) {
        out.byte(3) ; out.byte(shapeTag(value.shape)) ; writePlacement(out,value.placement()) ;
      } else if constexpr(std::is_same_v<T,ComposeRecord>) {
        out.byte(4) ; out.u64(ids.at(value.parent)) ; out.u64(ids.at(value.placement)) ;
      } else if constexpr(std::is_same_v<T,BoundarySurfaceRecord>) {
        out.byte(5) ; for(auto edge:value.edges) out.u64(ids.at(edge)) ;
      } else if constexpr(std::is_same_v<T,ConeRecord>) {
        out.byte(6) ; out.u64(ids.at(value.base)) ; for(double coordinate:value.apex) out.real(coordinate) ;
      } else if constexpr(std::is_same_v<T,FanRecord>) {
        out.byte(7) ; out.byte(shapeTag(value.shape)) ; writePoints(out,value.corners()) ;
        for(double coordinate:value.center()) out.real(coordinate) ;
      } else if constexpr(std::is_same_v<T,SurfacePatchRecord>) {
        out.byte(8) ; out.u64(ids.at(value.source)) ; writePoints(out,value.chartCorners) ;
      } else if constexpr(std::is_same_v<T,CubicTriangleRecord>) {
        out.byte(10) ; writePoints(out,value.controls) ;
      } else if constexpr(isCubicVolumeRecord<T>) {
        static_assert(std::is_same_v<T,CubicHexRecord> || std::is_same_v<T,CubicTetRecord> ||
                      std::is_same_v<T,CubicPrismRecord> || std::is_same_v<T,CubicPyramidRecord>,
                      "cubic volume needs an explicit wire tag") ;
        constexpr unsigned tag=std::is_same_v<T,CubicHexRecord> ? 11 :
          (std::is_same_v<T,CubicTetRecord> ? 12 : (std::is_same_v<T,CubicPrismRecord> ? 13 : 14)) ;
        out.byte(tag) ; for(auto edge:value.edges()) out.u64(ids.at(edge)) ;
        out.byte(value.reversed() & 255) ; out.byte(value.reversed()>>8) ;
      } else {
        static_assert(std::is_same_v<T,BoundaryFanRecord>,"geometry record needs an explicit wire tag") ;
        out.byte(9) ; for(auto edge:value.edges()) out.u64(ids.at(edge)) ;
        for(double coordinate:value.center()) out.real(coordinate) ;
      }
    },record) ;
  }

  template<std::size_t EdgeCount>
  GeometryRecord readCubicVolume(RecordReader &in, const GeometryLimits &limits, std::size_t id) {
    if(limits.maxDegree < 3) throw std::length_error("snapshot degree limit exceeded") ;
    std::array<MapId,EdgeCount> edges ;
    for(auto &edge:edges) {
      const auto value=in.u64() ;
      if(value>=id) throw std::invalid_argument("snapshot geometry reference must precede its use") ;
      edge=static_cast<MapId>(value) ;
    }
    const unsigned low=in.byte(), high=in.byte() ;
    const unsigned reversed=low | (high<<8) ;
    if(reversed>=(1u<<EdgeCount)) throw std::invalid_argument("snapshot cubic volume has an invalid edge reversal mask") ;
    return CubicVolumeRecord<EdgeCount>{edges,static_cast<std::uint16_t>(reversed)} ;
  }

  GeometryRecord readMap(RecordReader &in, const GeometryLimits &limits, std::size_t id) {
    switch(in.byte()) {
    case 1: {
      const auto shape=readShape(in) ;
      return CornerRecord{shape,readPoints(in,shape==Shape::polygon ? limits.maxControls:referenceElement(shape).vertexCount)} ;
    }
    case 2: {
      const auto shape=readShape(in) ;
      std::array<unsigned,3> degree ;
      for(auto &value:degree) {
        value=in.u32() ;
        if(value>limits.maxDegree) throw std::length_error("snapshot degree limit exceeded") ;
      }
      return TensorBezierRecord{shape,degree,readPoints(in,limits.maxControls)} ;
    }
    case 3: {
      const auto shape=readShape(in) ;
      return AffineRecord{shape,readPlacement(in)} ;
    }
    case 4: {
      const auto parent=in.u64(), placement=in.u64() ;
      if(parent>=id || placement>=id)
        throw std::invalid_argument("snapshot geometry reference must precede its use") ;
      return ComposeRecord{static_cast<MapId>(parent),static_cast<MapId>(placement)} ;
    }
    case 5: {
      BoundarySurfaceRecord result ;
      for(auto &edge:result.edges) {
        const auto value=in.u64() ;
        if(value>=id) throw std::invalid_argument("snapshot geometry reference must precede its use") ;
        edge=static_cast<MapId>(value) ;
      }
      return result ;
    }
    case 6: {
      const auto base=in.u64() ;
      if(base>=id) throw std::invalid_argument("snapshot geometry reference must precede its use") ;
      Point apex ; for(double &coordinate:apex) coordinate=in.real() ;
      return ConeRecord{static_cast<MapId>(base),apex} ;
    }
    case 7: {
      const auto shape=readShape(in) ;
      auto corners=readPoints(in,limits.maxControls) ;
      Point center ; for(double &coordinate:center) coordinate=in.real() ;
      return FanRecord{shape,std::move(corners),center} ;
    }
    case 8: {
      const auto source=in.u64() ;
      if(source>=id) throw std::invalid_argument("snapshot geometry reference must precede its use") ;
      return SurfacePatchRecord{static_cast<MapId>(source),readPoints(in,limits.maxControls)} ;
    }
    case 9: {
      std::array<MapId,4> edges ;
      for(auto &edge:edges) {
        const auto value=in.u64() ;
        if(value>=id) throw std::invalid_argument("snapshot geometry reference must precede its use") ;
        edge=static_cast<MapId>(value) ;
      }
      Point center ; for(auto &coordinate:center) coordinate=in.real() ;
      return BoundaryFanRecord{edges,center} ;
    }
    case 10: {
      if(limits.maxDegree < 3) throw std::length_error("snapshot degree limit exceeded") ;
      return CubicTriangleRecord{readPoints(in,std::min(std::size_t(10),limits.maxControls))} ;
    }
    case 11: return readCubicVolume<12>(in,limits,id) ;
    case 12: return readCubicVolume<6>(in,limits,id) ;
    case 13: return readCubicVolume<9>(in,limits,id) ;
    case 14: return readCubicVolume<8>(in,limits,id) ;
    }
    throw std::invalid_argument("unknown snapshot geometry tag") ;
  }

  struct GeometryOrder {
    std::vector<MapId> ids ;
    std::vector<std::vector<unsigned char>> records ;
    std::size_t fractions=0 ;
  } ;

  GeometryOrder geometryOrder(const GeometryStore &geometry, const RecordLimits &limits) {
    if(geometry.size()>limits.geometry.maxRecords)
      throw std::length_error("snapshot map limit exceeded") ;
    if(geometry.retainedBytes()>limits.geometry.exact.maxBytes)
      throw std::length_error("snapshot geometry storage limit exceeded") ;
    std::size_t controls=0 ;
    for(const auto &record:geometry.records()) std::visit([&](const auto &value) {
      using T=std::decay_t<decltype(value)> ;
      std::size_t count=0 ;
      if constexpr(std::is_same_v<T,CornerRecord>) count=value.corners.size() ;
      else if constexpr(std::is_same_v<T,TensorBezierRecord>) {
        count=value.controls.size() ;
        for(unsigned degree:value.degree)
          if(degree>limits.geometry.maxDegree) throw std::length_error("snapshot degree limit exceeded") ;
      }
      else if constexpr(std::is_same_v<T,CubicTriangleRecord>) {
        count=value.controls.size() ;
        if(limits.geometry.maxDegree < 3) throw std::length_error("snapshot degree limit exceeded") ;
      }
      else if constexpr(isCubicVolumeRecord<T>) {
        if(limits.geometry.maxDegree < 3) throw std::length_error("snapshot degree limit exceeded") ;
      }
      else if constexpr(std::is_same_v<T,ConeRecord> || std::is_same_v<T,BoundaryFanRecord>) count=1 ;
      else if constexpr(std::is_same_v<T,FanRecord>) {
        if(value.corners().size()>=limits.geometry.maxControls-controls)
          throw std::length_error("snapshot control limit exceeded") ;
        count=value.corners().size()+1 ;
      } else if constexpr(std::is_same_v<T,SurfacePatchRecord>) count=value.chartCorners.size() ;
      if(count>limits.geometry.maxControls-controls)
        throw std::length_error("snapshot control limit exceeded") ;
      controls+=count ;
    },record) ;
    std::vector<std::size_t> depths(geometry.size(),0) ;
    std::map<std::size_t,std::vector<MapId>> levels ;
    for(MapId id=0; id<geometry.size(); ++id) {
      if(geometry.evaluationCost(id)>limits.geometry.exact.maxOperations)
        throw std::length_error("snapshot geometry evaluation limit exceeded") ;
      for(auto dependency:geometryDependencies(geometry.record(id))) {
        if(dependency>=id) throw std::invalid_argument("snapshot geometry reference must precede its use") ;
        depths[id]=std::max(depths[id],depths[dependency]+1) ;
      }
      if(depths[id]>limits.geometry.maxCompositionDepth)
        throw std::length_error("snapshot composition depth limit exceeded") ;
      levels[depths[id]].push_back(id) ;
    }
    GeometryOrder result ;
    result.ids.resize(geometry.size()) ;
    std::size_t retained=0 ;
    for(const auto &level:levels) {
      std::vector<std::pair<std::vector<unsigned char>,MapId>> ordered ;
      ordered.reserve(level.second.size()) ;
      for(MapId id:level.second) {
        RecordWriter out(limits.maxBytes,limits.geometry.exact) ;
        writeMap(out,geometry.record(id),result.ids) ;
        if(out.bytes().size()>limits.maxBytes-retained)
          throw std::length_error("snapshot byte limit exceeded") ;
        retained+=out.bytes().size() ;
        ordered.emplace_back(out.take(),id) ;
      }
      std::sort(ordered.begin(),ordered.end()) ;
      for(auto &item:ordered) {
        if(result.records.empty() || result.records.back()!=item.first) {
          if(item.first.front()==3) result.fractions+=12 ;
          result.records.push_back(std::move(item.first)) ;
        }
        result.ids[item.second]=result.records.size()-1 ;
      }
    }
    return result ;
  }

  void writeGeometry(RecordWriter &out, const GeometryOrder &order) {
    out.operations(order.fractions) ;
    out.u64(order.records.size()) ;
    for(const auto &record:order.records) out.append(record) ;
  }

  GeometryStore readGeometry(RecordReader &in, const RecordLimits &limits) {
    GeometryStore result(limits.geometry) ;
    const auto count=in.count(limits.geometry.maxRecords,2) ;
    for(std::size_t id=0; id<count; ++id) result.append(readMap(in,limits.geometry,id)) ;
    return result ;
  }

  unsigned char familyTag(SplitFamily family) {
    switch(family) {
    case SplitFamily::tensor: return 1 ;
    case SplitFamily::triangle: return 2 ;
    case SplitFamily::prism: return 3 ;
    case SplitFamily::tet: return 4 ;
    case SplitFamily::pyramid: return 5 ;
    case SplitFamily::vertexVolume: return 6 ;
    case SplitFamily::planar: return 7 ;
    }
    throw std::invalid_argument("unknown snapshot split family") ;
  }

  SplitFamily readFamily(RecordReader &in) {
    switch(in.byte()) {
    case 1: return SplitFamily::tensor ;
    case 2: return SplitFamily::triangle ;
    case 3: return SplitFamily::prism ;
    case 4: return SplitFamily::tet ;
    case 5: return SplitFamily::pyramid ;
    case 6: return SplitFamily::vertexVolume ;
    case 7: return SplitFamily::planar ;
    }
    throw std::invalid_argument("unknown snapshot split family") ;
  }

  void writeCuts(RecordWriter &out, const std::vector<Fraction> &cuts) {
    out.u64(cuts.size()) ;
    for(const auto &cut:cuts) out.fraction(cut) ;
  }

  std::vector<Fraction> readCuts(RecordReader &in, std::size_t limit) {
    std::vector<Fraction> result(in.count(limit,minFractionBytes)) ;
    Fraction previous=0 ;
    for(auto &cut:result) {
      cut=in.fraction() ;
      if(cut<=previous || cut>=Fraction(1))
        throw std::invalid_argument("invalid snapshot cuts") ;
      previous=cut ;
    }
    return result ;
  }

  void writeTriangle(RecordWriter &out, const TriangleSplit &request) {
    switch(request.kind) {
    case TriangleSplit::Kind::none: out.byte(0) ; break ;
    case TriangleSplit::Kind::bisect:
      out.byte(1) ; out.u32(request.edge) ; out.fraction(request.position) ; break ;
    case TriangleSplit::Kind::uniform: out.byte(2) ; break ;
    default: throw std::invalid_argument("unknown snapshot triangle split") ;
    }
  }

  TriangleSplit readTriangle(RecordReader &in) {
    switch(in.byte()) {
    case 0: return TriangleSplit() ;
    case 1: {
      const auto edge=in.u32() ;
      const auto position=in.fraction() ;
      if(edge>=3 || position<=Fraction(0) || position>=Fraction(1))
        throw std::invalid_argument("invalid snapshot triangle cut") ;
      return TriangleSplit::bisectEdge(edge,position) ;
    }
    case 2: return TriangleSplit::uniform4() ;
    }
    throw std::invalid_argument("unknown snapshot triangle split tag") ;
  }

}
namespace detail {
  void writeSplitOperation(RecordWriter &out, const SplitOperation &operation) {
    std::visit([&](const auto &request) {
      using T=std::decay_t<decltype(request)> ;
      if constexpr(std::is_same_v<T,TensorSplit>) {
        out.byte(1) ;
        for(const auto &cuts:request.cuts) writeCuts(out,cuts) ;
      } else if constexpr(std::is_same_v<T,TriangleSplit>) {
        out.byte(2) ; writeTriangle(out,request) ;
      } else if constexpr(std::is_same_v<T,PrismSplit>) {
        out.byte(3) ; writeTriangle(out,request.base) ; writeCuts(out,request.axialCuts) ;
      } else if constexpr(std::is_same_v<T,TetSplit>) {
        out.byte(4) ;
        switch(request.kind) {
        case TetSplit::Kind::none: out.byte(0) ; break ;
        case TetSplit::Kind::bisect:
          out.byte(1) ; out.u32(request.edge) ; out.fraction(request.position) ; break ;
        case TetSplit::Kind::uniform: out.byte(2) ; out.u32(request.diagonal) ; break ;
        default: throw std::invalid_argument("unknown snapshot tetrahedron split") ;
        }
      } else {
        static_assert(std::is_same_v<T,PyramidSplit>,"split operation needs an explicit wire tag") ;
        out.byte(5) ;
        switch(request.kind) {
        case PyramidSplit::Kind::none: out.byte(0) ; break ;
        case PyramidSplit::Kind::uniform: out.byte(1) ; break ;
        default: throw std::invalid_argument("unknown snapshot pyramid split") ;
        }
      }
    },operation) ;
  }

  SplitOperation readSplitOperation(RecordReader &in, std::size_t limit) {
    switch(in.byte()) {
    case 1: {
      TensorSplit result ;
      for(auto &cuts:result.cuts) cuts=readCuts(in,limit) ;
      return result ;
    }
    case 2: return readTriangle(in) ;
    case 3: {
      PrismSplit result ;
      result.base=readTriangle(in) ; result.axialCuts=readCuts(in,limit) ;
      return result ;
    }
    case 4:
      switch(in.byte()) {
      case 0: return TetSplit() ;
      case 1: {
        const auto edge=in.u32() ;
        const auto position=in.fraction() ;
        if(edge>=6 || position<=Fraction(0) || position>=Fraction(1))
          throw std::invalid_argument("invalid snapshot tetrahedron cut") ;
        return TetSplit::bisectEdge(edge,position) ;
      }
      case 2: {
        const auto diagonal=in.u32() ;
        if(diagonal>=3) throw std::invalid_argument("invalid snapshot tetrahedron diagonal") ;
        return TetSplit::uniform8(diagonal) ;
      }
      }
      break ;
    case 5:
      switch(in.byte()) {
      case 0: return PyramidSplit() ;
      case 1: return PyramidSplit::uniform10() ;
      }
      break ;
    }
    throw std::invalid_argument("unknown snapshot operation tag") ;
  }

}
namespace {
  struct HistoryOrder {
    std::vector<TreeNodeId> nodes, nodeIds ;
    std::vector<std::size_t> splits, splitIds ;
  } ;

  HistoryOrder historyOrder(const RefinementHistory &history, const RecordLimits &limits) {
    const auto &nodes=history.nodes() ;
    const auto &splits=history.splits() ;
    if(nodes.empty()) throw std::invalid_argument("snapshot has no root") ;
    if(nodes.size()>limits.maxNodes || splits.size()>limits.maxSplits)
      throw std::length_error("snapshot history count limit exceeded") ;
    HistoryOrder order ;
    order.nodeIds.resize(nodes.size(),noTreeNode) ;
    order.splitIds.resize(splits.size(),noTreeNode) ;
    std::map<TreeNodeId,std::vector<std::size_t>> children ;
    for(std::size_t id=0; id<splits.size(); ++id) {
      if(splits[id].parent>=nodes.size() || !splits[id].construction)
        throw std::invalid_argument("invalid snapshot split parent or construction") ;
      children[splits[id].parent].push_back(id) ;
    }
    for(auto &item:children) {
      auto &list=item.second ;
      std::sort(list.begin(),list.end(),[&](std::size_t a, std::size_t b) {
        return *splits[a].construction<*splits[b].construction ;
      }) ;
      for(std::size_t i=1; i<list.size(); ++i)
        if(*splits[list[i-1]].construction==*splits[list[i]].construction)
          throw std::invalid_argument("duplicate snapshot construction") ;
    }
    std::vector<std::size_t> depths(nodes.size(),0) ;
    order.nodes.push_back(0) ; order.nodeIds[0]=0 ;
    for(std::size_t cursor=0; cursor<order.nodes.size(); ++cursor) {
      const TreeNodeId parent=order.nodes[cursor] ;
      const auto found=children.find(parent) ;
      if(found==children.end()) continue ;
      if(depths[parent]>=limits.maxDepth) throw std::length_error("snapshot depth limit exceeded") ;
      for(std::size_t splitId:found->second) {
        order.splitIds[splitId]=order.splits.size() ; order.splits.push_back(splitId) ;
        for(TreeNodeId child:splits[splitId].children) {
          if(child>=nodes.size() || nodes[child].parent!=parent || order.nodeIds[child]!=noTreeNode)
            throw std::invalid_argument("invalid snapshot child range") ;
          order.nodeIds[child]=order.nodes.size() ; order.nodes.push_back(child) ;
          depths[child]=depths[parent]+1 ;
        }
      }
    }
    if(order.nodes.size()!=nodes.size() || order.splits.size()!=splits.size())
      throw std::invalid_argument("unreachable snapshot records") ;
    return order ;
  }

  std::size_t readIndex(RecordReader &in, std::size_t count, bool allowAbsent=false) {
    const auto index=in.u64() ;
    if(allowAbsent && index==absent) return noTreeNode ;
    if(index>=count) throw std::invalid_argument("snapshot index out of range") ;
    return static_cast<std::size_t>(index) ;
  }

  void writeHeader(RecordWriter &out, std::uint32_t magic, std::uint32_t version) {
    out.u32(magic) ; out.u32(version) ;
  }

  void readHeader(RecordReader &in, std::uint32_t magic, std::uint32_t version) {
    if(in.u32()!=magic || in.u32()!=version)
      throw std::invalid_argument("unknown snapshot header or version") ;
  }

}

  std::vector<unsigned char> encodeGeometry(const GeometryStore &geometry, const RecordLimits &limits) {
    std::vector<MapId> ids ; return encodeGeometry(geometry,ids,limits) ;
  }

  std::vector<unsigned char> encodeGeometry(const GeometryStore &geometry, std::vector<MapId> &encodedIds,
                                            const RecordLimits &limits) {
    const auto order=geometryOrder(geometry,limits) ;
    RecordWriter out(limits.maxBytes,limits.geometry.exact) ;
    writeHeader(out,geometryMagic,geometryVersion) ;
    writeGeometry(out,order) ;
    encodedIds=order.ids ;
    return out.take() ;
  }

  GeometryStore decodeGeometry(const std::vector<unsigned char> &record, const RecordLimits &limits) {
    RecordReader in(record,limits.maxBytes,limits.geometry.exact) ;
    readHeader(in,geometryMagic,geometryVersion) ;
    auto result=readGeometry(in,limits) ;
    in.finish() ;
    return result ;
  }

  std::vector<unsigned char> encodeRefinement(const RefinementHistory &history, const RecordLimits &limits) {
    std::vector<TreeNodeId> ids ; return encodeRefinement(history,ids,limits) ;
  }

  std::vector<unsigned char> encodeRefinement(const RefinementHistory &history, std::vector<TreeNodeId> &encodedIds,
                                              const RecordLimits &limits) {
    if(history.retainedBytes()>limits.geometry.exact.maxBytes)
      throw std::length_error("snapshot history storage limit exceeded") ;
    const auto geometry=geometryOrder(history.geometry(),limits) ;
    const auto order=historyOrder(history,limits) ;
    RecordWriter out(limits.maxBytes,limits.geometry.exact) ;
    writeHeader(out,historyMagic,historyVersion) ;
    out.u64(history.rootKey()) ; out.u64(history.geometryRevision()) ;
    writeGeometry(out,geometry) ;
    out.u64(order.nodes.size()) ; out.u64(order.splits.size()) ;
    for(TreeNodeId id:order.nodes) {
      const auto &node=history.node(id) ;
      out.u64(node.parent==noTreeNode ? absent : order.nodeIds.at(node.parent)) ;
      out.u32(node.ordinal) ;
      out.u64(geometry.ids.at(node.map)) ;
      out.u64(node.activeSplit==noTreeNode ? absent : order.splitIds.at(node.activeSplit)) ;
    }
    for(std::size_t id:order.splits) {
      const auto &split=history.splits()[id] ;
      out.u64(order.nodeIds.at(split.parent)) ;
      out.byte(familyTag(split.construction->family)) ; out.u32(split.construction->version) ;
      if(split.operation) detail::writeSplitOperation(out,*split.operation) ;
      else if(split.construction->family != SplitFamily::vertexVolume && split.construction->family != SplitFamily::planar)
        throw std::invalid_argument("native snapshot construction has no request") ;
      if(split.children.empty()) throw std::invalid_argument("empty snapshot construction") ;
      out.u64(order.nodeIds.at(split.children.front())) ; out.u64(split.children.size()) ;
      if(split.construction->childCount()!=split.children.size())
        throw std::invalid_argument("snapshot construction child count differs") ;
      if(split.construction->vertexVolume) {
        const auto &mapped=*split.construction->vertexVolume ;
        writePoints(out,mapped.points) ;
        for(const auto &child:mapped.children) for(unsigned vertex:child.vertices) out.u32(vertex) ;
      } else if(split.construction->planar) {
        const auto &planar=*split.construction->planar ;
        switch(planar.centerMethod) {
        case CenterMethod::reference: out.byte(1) ; break ;
        case CenterMethod::weighted: out.byte(2) ; break ;
        case CenterMethod::integrated: out.byte(3) ; break ;
        default: throw std::invalid_argument("unknown snapshot planar center method") ;
        }
        writePoints(out,planar.points) ;
        for(const auto &child:planar.children) for(unsigned vertex:child) out.u32(vertex) ;
      } else for(const auto &child:split.construction->children) {
        out.byte(shapeTag(child.shape)) ; writePlacement(out,child.placement) ;
      }
    }
    encodedIds=order.nodeIds ;
    return out.take() ;
  }

  RefinementHistory decodeRefinement(const std::vector<unsigned char> &record, const RecordLimits &limits) {
    RecordReader in(record,limits.maxBytes,limits.geometry.exact) ;
    readHeader(in,historyMagic,historyVersion) ;
    const auto root=in.u64(), revision=in.u64() ;
    auto geometry=readGeometry(in,limits) ;
    const auto nodeCount=in.count(limits.maxNodes,minNodeBytes) ;
    const auto splitCount=in.count(limits.maxSplits,minSplitBytes) ;
    if(nodeCount==0) throw std::invalid_argument("snapshot has no root") ;
    if(nodeCount>in.remaining()/minNodeBytes || splitCount>(in.remaining()-nodeCount*minNodeBytes)/minSplitBytes)
      throw std::invalid_argument("truncated snapshot tables") ;
    ExactBudget retained(limits.geometry.exact) ;
    retained.retain(geometry.retainedBytes()) ;
    if(nodeCount>limits.geometry.exact.maxBytes/sizeof(CellNode) ||
       splitCount>limits.geometry.exact.maxBytes/sizeof(AcceptedSplit))
      throw std::length_error("snapshot history storage limit exceeded") ;
    retained.retain(nodeCount*sizeof(CellNode)) ;
    retained.retain(splitCount*sizeof(AcceptedSplit)) ;
    std::vector<CellNode> nodes(nodeCount) ;
    std::vector<std::size_t> depths(nodeCount,0) ;
    for(TreeNodeId id=0; id<nodeCount; ++id) {
      auto &node=nodes[id] ;
      node.parent=readIndex(in,id,true) ;
      if((node.parent==noTreeNode)!=(id==0))
        throw std::invalid_argument("snapshot nodes must follow their parent") ;
      if(id) {
        if(depths[node.parent]>=limits.maxDepth) throw std::length_error("snapshot depth limit exceeded") ;
        depths[id]=depths[node.parent]+1 ;
      }
      node.ordinal=in.u32() ;
      node.map=readIndex(in,geometry.size()) ;
      node.activeSplit=readIndex(in,splitCount,true) ;
    }
    std::vector<AcceptedSplit> splits ;
    splits.reserve(splitCount) ;
    std::vector<bool> claimed(nodeCount,false) ;
    claimed[0]=true ;
    for(std::size_t id=0; id<splitCount; ++id) {
      const auto parent=readIndex(in,nodeCount) ;
      auto construction=std::make_shared<SplitConstruction>() ;
      construction->parentShape=geometry.map(nodes[parent].map)->shape() ;
      construction->family=readFamily(in) ; construction->version=in.u32() ;
      if(construction->version!=1) throw std::invalid_argument("unknown snapshot construction version") ;
      std::optional<SplitOperation> operation ;
      if(construction->family != SplitFamily::vertexVolume && construction->family != SplitFamily::planar)
        operation=detail::readSplitOperation(in,limits.maxNodes) ;
      const auto first=readIndex(in,nodeCount) ;
      const bool mapped=construction->family == SplitFamily::vertexVolume ;
      const bool planar=construction->family == SplitFamily::planar ;
      const auto count=in.count(nodeCount-first,mapped ? 32 : (planar ? 16 : minChildBytes)) ;
      if(count==0 || first<=parent) throw std::invalid_argument("invalid snapshot child range") ;
      for(std::size_t child=first; child<first+count; ++child) {
        if(claimed[child] || nodes[child].parent!=parent || nodes[child].ordinal!=child-first)
          throw std::invalid_argument("conflicting snapshot child range") ;
        claimed[child]=true ;
      }
      if(mapped) {
        const auto &element=referenceElement(construction->parentShape) ;
        const auto expected=element.vertexCount+element.edgeCount+element.boundaryCount+1 ;
        auto points=readPoints(in,expected) ;
        if(count != element.vertexCount || points.size() != expected)
          throw std::invalid_argument("invalid snapshot mapped construction counts") ;
        retained.retain(points.size()*sizeof(Point)+count*(sizeof(VertexVolumeChild)+8*sizeof(Point))) ;
        std::vector<std::array<unsigned,8>> corners(count) ;
        for(auto &child:corners) for(unsigned &vertex:child) vertex=in.u32() ;
        construction->vertexVolume=std::make_shared<const VertexVolumeSplit>(
          restoreVertexVolumeSplit(construction->parentShape,std::move(points),corners,limits.jacobian)) ;
      } else if(planar) {
        const auto shape=construction->parentShape ;
        if(shape != Shape::triangle && shape != Shape::quad)
          throw std::invalid_argument("snapshot planar construction requires a triangle or quad") ;
        const auto n=referenceElement(shape).vertexCount ;
        if(count != n) throw std::invalid_argument("invalid snapshot planar child count") ;
        PlanarSplit split ; split.shape=shape ;
        switch(in.byte()) {
        case 1: split.centerMethod=CenterMethod::reference ; break ;
        case 2: split.centerMethod=CenterMethod::weighted ; break ;
        case 3: split.centerMethod=CenterMethod::integrated ; break ;
        default: throw std::invalid_argument("unknown snapshot planar center method") ;
        }
        split.points=readPoints(in,2*n+1) ;
        if(split.points.size() != 2*n+1) throw std::invalid_argument("invalid snapshot planar point count") ;
        retained.retain(sizeof(PlanarSplit)+split.points.size()*sizeof(Point)+count*sizeof(std::array<unsigned,4>)) ;
        split.children.resize(count) ;
        for(auto &child:split.children) for(unsigned &vertex:child) vertex=in.u32() ;
        construction->planar=std::make_shared<const PlanarSplit>(std::move(split)) ;
      } else {
        if(count>limits.geometry.exact.maxBytes/sizeof(Child))
          throw std::length_error("snapshot construction storage limit exceeded") ;
        retained.retain(count*sizeof(Child)) ;
        construction->children.reserve(count) ;
        for(std::size_t ordinal=0; ordinal<count; ++ordinal) {
          const auto shape=readShape(in) ;
          construction->children.push_back({shape,readPlacement(in)}) ;
        }
      }
      splits.push_back({parent,NodeRange(first,count),std::move(operation),std::move(construction)}) ;
    }
    in.finish() ;
    if(std::find(claimed.begin(),claimed.end(),false)!=claimed.end())
      throw std::invalid_argument("snapshot has unclaimed child nodes") ;
    HistoryLimits historyLimits ;
    historyLimits.maxNodes=limits.maxNodes ; historyLimits.maxSplits=limits.maxSplits ;
    historyLimits.maxDepth=static_cast<unsigned>(std::min<std::size_t>(limits.maxDepth,std::numeric_limits<unsigned>::max())) ;
    historyLimits.exact=limits.geometry.exact ; historyLimits.geometry=limits.geometry ;
    historyLimits.jacobian=limits.jacobian ;
    return RefinementHistory::restore(root,revision,std::move(geometry),std::move(nodes),
                                      std::move(splits),historyLimits) ;
  }

}
}
