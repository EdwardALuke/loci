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

#include <FVMAdapt/parametric/boundary_record.h>
#include <FVMAdapt/parametric/refinement_record.h>
#include <FVMAdapt/parametric/surface_map.h>

#include "record_io.h"
#include "geometry_validity.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <type_traits>

namespace Loci {
namespace FVMAdapt {
namespace {

  using detail::RecordReader ;
  using detail::RecordWriter ;
  constexpr std::uint32_t boundaryMagic=0x42415646, boundaryVersion=1 ; // "FVAB"
  constexpr std::uint64_t absent=std::numeric_limits<std::uint64_t>::max() ;

  void require(bool condition, const char *message) {
    if(!condition) throw std::invalid_argument(message) ;
  }

  void count(std::size_t amount, std::size_t &used, std::size_t limit) {
    if(used>limit || amount>limit-used)
      throw std::length_error("boundary record count limit exceeded") ;
    used+=amount ;
  }

  RecordLimits geometryLimits(const BoundaryRecordLimits &limits) {
    RecordLimits result ;
    result.maxBytes=limits.maxBytes ; result.geometry=limits.geometry ;
    return result ;
  }

  unsigned char centerTag(CenterMethod method) {
    switch(method) {
    case CenterMethod::reference: return 1 ;
    case CenterMethod::weighted: return 2 ;
    case CenterMethod::integrated: return 3 ;
    }
    throw std::invalid_argument("unknown boundary center method") ;
  }

  CenterMethod readCenter(RecordReader &in) {
    switch(in.byte()) {
    case 1: return CenterMethod::reference ;
    case 2: return CenterMethod::weighted ;
    case 3: return CenterMethod::integrated ;
    }
    throw std::invalid_argument("unknown boundary center method tag") ;
  }

  void writeId(RecordWriter &out, BoundaryId id) {
    out.u64(id==noBoundary ? absent : id) ;
  }

  BoundaryId readId(RecordReader &in) {
    const auto id=in.u64() ;
    if(id==absent) return noBoundary ;
    if(id>=noBoundary) throw std::length_error("boundary ID exceeds local range") ;
    return static_cast<BoundaryId>(id) ;
  }

  void writeIds(RecordWriter &out, const std::vector<BoundaryId> &ids) {
    out.u64(ids.size()) ;
    for(const auto id:ids) writeId(out,id) ;
  }


  bool connects(const EdgeRecord &edge, BoundaryId first, BoundaryId second) {
    return (edge.vertices[0]==first && edge.vertices[1]==second) ||
           (edge.vertices[1]==first && edge.vertices[0]==second) ;
  }

  void checkGeometryLimits(const GeometryStore &geometry, const GeometryLimits &limits) {
    if(geometry.size()>limits.maxRecords || geometry.retainedBytes()>limits.exact.maxBytes)
      throw std::length_error("boundary geometry limit exceeded") ;
    std::vector<std::size_t> depths(geometry.size(),0) ;
    std::size_t controls=0 ;
    for(MapId id=0; id<geometry.size(); ++id) {
      if(geometry.evaluationCost(id)>limits.exact.maxOperations)
        throw std::length_error("boundary geometry evaluation limit exceeded") ;
      std::visit([&](const auto &record) {
        using T=std::decay_t<decltype(record)> ;
        if constexpr(std::is_same_v<T,CornerRecord>) {
          count(record.corners.size(),controls,limits.maxControls) ;
        } else if constexpr(std::is_same_v<T,FanRecord>) {
          count(record.corners().size(),controls,limits.maxControls) ;
          count(1,controls,limits.maxControls) ;
        } else if constexpr(std::is_same_v<T,SurfacePatchRecord>) {
          count(record.chartCorners.size(),controls,limits.maxControls) ;
        } else if constexpr(std::is_same_v<T,BoundaryFanRecord>) {
          count(1,controls,limits.maxControls) ;
        } else if constexpr(std::is_same_v<T,TensorBezierRecord>) {
          count(record.controls.size(),controls,limits.maxControls) ;
          for(unsigned degree:record.degree)
            if(degree>limits.maxDegree) throw std::length_error("boundary geometry degree limit exceeded") ;
        } else if constexpr(std::is_same_v<T,AffineRecord>) {
          for(const auto &coordinate:record.placement().origin) checkExact(coordinate,limits.exact) ;
          for(const auto &column:record.placement().columns)
            for(const auto &coordinate:column) checkExact(coordinate,limits.exact) ;
        } else if constexpr(std::is_same_v<T,ConeRecord>) count(1,controls,limits.maxControls) ;
      },geometry.record(id)) ;
      for(auto dependency:geometryDependencies(geometry.record(id))) {
        require(dependency<id,"boundary map references must precede their use") ;
        depths[id]=std::max(depths[id],depths[dependency]+1) ;
      }
      if(depths[id]>limits.maxCompositionDepth)
        throw std::length_error("boundary geometry depth limit exceeded") ;
    }
  }

  std::size_t integerBytes(Fraction::Integer value) {
    if(value<0) value=-value ;
    if(value==0) return 0 ;
    const auto bits=std::size_t(boost::multiprecision::msb(value))+1 ;
    return bits/8+(bits%8!=0) ;
  }

}

  void validateBoundaries(const BoundaryRecords &records, const BoundaryRecordLimits &limits) {
    std::size_t entities=0, splits=0 ;
    count(records.points.size(),entities,limits.maxEntities) ;
    count(records.edges.size(),entities,limits.maxEntities) ;
    count(records.faces.size(),entities,limits.maxEntities) ;
    checkGeometryLimits(records.geometry,limits.geometry) ;
    ExactBudget budget(limits.geometry.exact) ;
    const auto payload=[&](std::size_t amount, std::size_t size) {
      if(amount>limits.geometry.exact.maxBytes/size)
        throw std::length_error("boundary record payload limit exceeded") ;
      budget.retain(amount*size) ;
    } ;
    budget.retain(records.geometry.retainedBytes()) ;
    payload(records.points.size(),sizeof(Point)) ;
    payload(records.edges.size(),sizeof(EdgeRecord)) ;
    payload(records.faces.size(),sizeof(FaceRecord)) ;
    const auto fraction=[&](const Fraction &value) {
      budget.check(value) ; budget.operations() ;
      budget.retain(integerBytes(value.numerator())) ;
      budget.retain(integerBytes(value.denominator())) ;
    } ;
    for(const auto &point:records.points)
      for(double coordinate:point)
        require(std::isfinite(coordinate),"nonfinite boundary point") ;
    const auto position=[&](MapId map, const Point &coordinate) {
      return records.geometry.map(map)->evaluate(coordinate).position ;
    } ;
    const auto samePoint=[](const Point &first, const Point &second, double extent) {
      require(detail::storedPointAgrees(first,second,extent),"boundary points disagree with their accepted geometry") ;
    } ;

    std::vector<unsigned char> edgeLinks(records.edges.size(),0), faceLinks(records.faces.size(),0) ;
    std::vector<std::size_t> edgeDepth(records.edges.size(),0), faceDepth(records.faces.size(),0) ;
    for(BoundaryId id=0; id<records.edges.size(); ++id) {
      const auto &edge=records.edges[id] ;
      require(edge.vertices[0]<records.points.size() && edge.vertices[1]<records.points.size() &&
              edge.vertices[0]!=edge.vertices[1],"invalid boundary edge endpoints") ;
      require(edge.geometry<records.geometry.size(),"invalid boundary edge map") ;
      require(records.geometry.map(edge.geometry)->shape()==Shape::segment,"boundary edge map must be a segment") ;
      const double extent=detail::pointExtent(2,[&](std::size_t v) { return records.points[edge.vertices[v]] ; }) ;
      for(unsigned v=0; v<2; ++v)
        samePoint(records.points[edge.vertices[v]],position(edge.geometry,{{double(v),0,0}}),extent) ;
      if(edge.parent!=noBoundary) {
        require(edge.parent<id,"boundary edge parent must precede its child") ;
        edgeDepth[id]=edgeDepth[edge.parent]+1 ;
        if(edgeDepth[id]>limits.maxDepth) throw std::length_error("boundary edge depth limit exceeded") ;
      }
      count(edge.cuts.size(),splits,limits.maxSplits) ;
      payload(edge.cuts.size(),sizeof(EdgeCut)) ; payload(edge.children.size(),sizeof(EdgeInterval)) ;
      require(edge.cuts.empty() ? edge.children.empty() : edge.children.size()==edge.cuts.size()+1,
              "boundary edge cuts must have a complete child partition") ;
      Fraction previous=0 ;
      std::set<BoundaryId> vertices{edge.vertices[0],edge.vertices[1]} ;
      for(const auto &cut:edge.cuts) {
        fraction(cut.position) ;
        require(cut.position>previous && cut.position<Fraction(1),"boundary edge cuts must be sorted and interior") ;
        require(cut.vertex<records.points.size() && vertices.insert(cut.vertex).second,"invalid boundary edge cut vertex") ;
        samePoint(records.points[cut.vertex],position(edge.geometry,{{cut.position.value(),0,0}}),extent) ;
        previous=cut.position ;
      }
      for(std::size_t c=0; c<edge.children.size(); ++c) {
        const auto &child=edge.children[c] ;
        fraction(child.begin) ; fraction(child.end) ;
        const Fraction begin=c ? edge.cuts[c-1].position : Fraction(0) ;
        const Fraction end=c<edge.cuts.size() ? edge.cuts[c].position : Fraction(1) ;
        require(child.begin==begin && child.end==end,"boundary edge children must partition their parent") ;
        require(child.edge>id && child.edge<records.edges.size(),"invalid boundary child edge") ;
        require(records.edges[child.edge].parent==id && !edgeLinks[child.edge],"inconsistent boundary edge ancestry") ;
        edgeLinks[child.edge]=1 ;
        const auto first=c ? edge.cuts[c-1].vertex : edge.vertices[0] ;
        const auto second=c<edge.cuts.size() ? edge.cuts[c].vertex : edge.vertices[1] ;
        require(connects(records.edges[child.edge],first,second),"boundary child edge endpoints do not match its interval") ;
        // Child edge maps may include a final reversal into canonical endpoint order.
        const auto &childEdge=records.edges[child.edge] ;
        MapId map=childEdge.geometry ;
        std::array<Point,2> endpoints{{{{0,0,0}},{{1,0,0}}}} ;
        while(map!=edge.geometry) {
          require(map<records.geometry.size(),"invalid boundary child edge map") ;
          const auto *composition=std::get_if<ComposeRecord>(&records.geometry.record(map)) ;
          require(composition!=nullptr,"boundary child edge must restrict its parent map") ;
          const auto *placement=std::get_if<CornerRecord>(&records.geometry.record(composition->placement)) ;
          require(placement && placement->shape==Shape::segment &&
                  placement->corners[0][1]==0 && placement->corners[0][2]==0 &&
                  placement->corners[1][1]==0 && placement->corners[1][2]==0,
                  "boundary child edge placement must follow its parent coordinates") ;
          for(auto &endpoint:endpoints) endpoint=position(composition->placement,endpoint) ;
          map=composition->parent ;
        }
        const bool forward=childEdge.vertices[0]==first ;
        samePoint(endpoints[0],{{(forward ? begin : end).value(),0,0}},1) ;
        samePoint(endpoints[1],{{(forward ? end : begin).value(),0,0}},1) ;
      }
    }
    for(BoundaryId id=0; id<records.edges.size(); ++id)
      require(bool(edgeLinks[id])==(records.edges[id].parent!=noBoundary),"unlinked boundary child edge") ;

    for(BoundaryId id=0; id<records.faces.size(); ++id) {
      const auto &face=records.faces[id] ;
      require(face.geometry<records.geometry.size(),"invalid boundary face map") ;
      const auto map=records.geometry.map(face.geometry) ; const auto shape=map->shape() ;
      require(map->dimension()==2,"boundary face map must be a surface") ;
      const auto coordinates=map->referenceCorners() ; const auto vertexCount=coordinates.size() ;
      const auto along=[&](unsigned e,double t) {
        Point q ; for(unsigned d=0;d<3;++d) q[d]=(1-t)*coordinates[e][d]+t*coordinates[(e+1)%vertexCount][d] ; return q ;
      } ;
      require(face.vertices.size()==vertexCount && face.edges.size()==vertexCount,"invalid boundary face loop size") ;
      payload(face.vertices.size()+face.edges.size(),sizeof(BoundaryId)) ;
      std::set<BoundaryId> vertices ;
      for(auto vertex:face.vertices)
        require(vertex<records.points.size() && vertices.insert(vertex).second,"invalid boundary face vertex") ;
      const double extent=detail::pointExtent(vertexCount,[&](std::size_t v) { return records.points[face.vertices[v]] ; }) ;
      for(std::size_t v=0; v<face.vertices.size(); ++v) {
        require(face.edges[v]<records.edges.size(),"invalid boundary face edge") ;
        const auto &edge=records.edges[face.edges[v]] ;
        require(connects(edge,face.vertices[v],face.vertices[(v+1)%vertexCount]),
                "boundary face edges must follow its vertex loop") ;
        samePoint(records.points[face.vertices[v]],position(face.geometry,coordinates[v]),extent) ;
        for(double coordinate:{.25,.5,.75}) {
          const double alongEdge=edge.vertices[0]==face.vertices[v] ? coordinate : 1-coordinate ;
          samePoint(position(edge.geometry,{{alongEdge,0,0}}),
                    position(face.geometry,along(v,coordinate)),extent) ;
        }
      }
      if(face.parent==noBoundary) {
        require(face.parentSplit==noBoundary,"root boundary face has a parent construction") ;
      } else {
        require(face.parent<id,"boundary face parent must precede its child") ;
        require(face.parentSplit<records.faces[face.parent].splits.size(),"invalid boundary parent construction") ;
        require(face.boundary==records.faces[face.parent].boundary,"boundary child face changed its label") ;
        faceDepth[id]=faceDepth[face.parent]+1 ;
        if(faceDepth[id]>limits.maxDepth) throw std::length_error("boundary face depth limit exceeded") ;
      }
      count(face.splits.size(),splits,limits.maxSplits) ;
      payload(face.splits.size(),sizeof(FaceSplitRecord)) ;
      for(std::size_t s=0; s<face.splits.size(); ++s) {
        const auto &split=face.splits[s] ;
        centerTag(split.centerMethod) ;
        const auto &center=split.referenceCenter ;
        require((shape==Shape::polygon ? insideReferencePolygon(vertexCount,center):map->containsReference(center)) && center[0]>0 && center[1]>0 &&
                (shape==Shape::triangle ? center[0]+center[1]<1 : center[0]<1 && center[1]<1),
                "boundary face center must lie inside its reference element") ;
        require(split.centerVertex<records.points.size() && !vertices.count(split.centerVertex),"invalid boundary face center vertex") ;
        samePoint(records.points[split.centerVertex],position(face.geometry,split.referenceCenter),extent) ;
        require(split.children.size()==vertexCount,"boundary corner construction needs one child per vertex") ;
        payload(split.children.size(),sizeof(BoundaryId)) ;
        for(std::size_t c=0; c<split.children.size(); ++c) {
          const auto child=split.children[c] ;
          require(child>id && child<records.faces.size(),"invalid boundary child face") ;
          const auto &childFace=records.faces[child] ;
          require(childFace.parent==id && childFace.parentSplit==s && !faceLinks[child],"inconsistent boundary face ancestry") ;
          const auto previous=(c+vertexCount-1)%vertexCount ;
          const auto midpoint=[&](std::size_t e) {
            const auto &cuts=records.edges[face.edges[e]].cuts ;
            const auto found=std::find_if(cuts.begin(),cuts.end(),[](const EdgeCut &cut) { return cut.position==Fraction(1,2) ; }) ;
            require(found!=cuts.end(),"boundary corner construction needs edge midpoints") ;
            return found->vertex ;
          } ;
          require(childFace.vertices==std::vector<BoundaryId>{face.vertices[c],midpoint(c),split.centerVertex,midpoint(previous)},
                  "boundary corner child has an inconsistent vertex order") ;
          require(childFace.geometry<records.geometry.size(),"invalid boundary child face map") ;
          std::vector<Point> corners={coordinates[c],along(c,.5),split.referenceCenter,along(previous,.5)} ;
          const auto &childRecord=records.geometry.record(childFace.geometry) ;
          if(const auto *childPatch=std::get_if<SurfacePatchRecord>(&childRecord)) {
            MapId source=face.geometry ;
            if(const auto *parentPatch=std::get_if<SurfacePatchRecord>(&records.geometry.record(face.geometry))) {
              source=parentPatch->source ;
              const auto *patch=dynamic_cast<const SurfacePatchMap *>(map.get()) ;
              require(patch,"boundary patch has an inconsistent evaluator") ;
              for(auto &point:corners) point=patch->chartPosition(point) ;
            }
            require(childPatch->source==source && childPatch->chartCorners==corners,
                    "boundary child patch does not match its accepted construction") ;
          } else {
            require(!hasSurfacePieces(*map),"piecewise boundary children require a retained source patch") ;
            const auto *composition=std::get_if<ComposeRecord>(&childRecord) ;
            require(composition && composition->parent==face.geometry,"boundary child face must restrict its parent map") ;
            const auto *placement=std::get_if<CornerRecord>(&records.geometry.record(composition->placement)) ;
            require(placement && placement->shape==Shape::quad && placement->corners==corners,
                    "boundary child face map does not match its accepted construction") ;
          }
          faceLinks[child]=1 ;
        }
      }
    }
    for(BoundaryId id=0; id<records.faces.size(); ++id)
      require(bool(faceLinks[id])==(records.faces[id].parent!=noBoundary),"unlinked boundary child face") ;
  }

  std::vector<unsigned char> encodeBoundaries(const BoundaryRecords &records, const BoundaryRecordLimits &limits) {
    validateBoundaries(records,limits) ;
    std::vector<MapId> mapIds ;
    const auto geometry=encodeGeometry(records.geometry,mapIds,geometryLimits(limits)) ;
    RecordWriter out(limits.maxBytes,limits.geometry.exact) ;
    out.u32(boundaryMagic) ; out.u32(boundaryVersion) ;
    out.u64(geometry.size()) ; out.append(geometry) ;
    out.u64(records.points.size()) ;
    for(const auto &point:records.points) for(double coordinate:point) out.real(coordinate) ;
    out.u64(records.edges.size()) ;
    for(const auto &edge:records.edges) {
      for(const auto vertex:edge.vertices) writeId(out,vertex) ;
      out.u64(mapIds.at(edge.geometry)) ; writeId(out,edge.parent) ;
      out.u64(edge.cuts.size()) ;
      for(const auto &cut:edge.cuts) { out.fraction(cut.position) ; writeId(out,cut.vertex) ; }
      out.u64(edge.children.size()) ;
      for(const auto &child:edge.children) {
        out.fraction(child.begin) ; out.fraction(child.end) ; writeId(out,child.edge) ;
      }
    }
    out.u64(records.faces.size()) ;
    for(const auto &face:records.faces) {
      writeId(out,face.parent) ; writeId(out,face.parentSplit) ; out.u64(mapIds.at(face.geometry)) ;
      static_assert(sizeof(int)<=4,"boundary labels must fit the signed 32-bit wire format") ;
      out.u32(static_cast<std::uint32_t>(face.boundary)) ;
      writeIds(out,face.vertices) ; writeIds(out,face.edges) ;
      out.u64(face.splits.size()) ;
      for(const auto &split:face.splits) {
        out.byte(centerTag(split.centerMethod)) ;
        for(double coordinate:split.referenceCenter) out.real(coordinate) ;
        writeId(out,split.centerVertex) ; writeIds(out,split.children) ;
      }
    }
    return out.take() ;
  }

  BoundaryRecords decodeBoundaries(const std::vector<unsigned char> &record, const BoundaryRecordLimits &limits) {
    RecordReader in(record,limits.maxBytes,limits.geometry.exact) ;
    require(in.u32()==boundaryMagic,"unknown boundary record header") ;
    require(in.u32()==boundaryVersion,"unsupported boundary record version") ;
    std::vector<unsigned char> geometry(in.count(limits.maxBytes)) ;
    for(auto &byte:geometry) byte=in.byte() ;
    BoundaryRecords result ;
    result.geometry=decodeGeometry(geometry,geometryLimits(limits)) ;
    std::size_t entities=0, splits=0, childEdges=0, childFaces=0 ;
    ExactBudget payload(limits.geometry.exact) ;
    payload.retain(result.geometry.retainedBytes()) ;
    const auto allocation=[&](std::size_t amount, std::size_t size) {
      if(amount>limits.geometry.exact.maxBytes/size)
        throw std::length_error("boundary record payload limit exceeded") ;
      payload.retain(amount*size) ;
      return amount ;
    } ;
    const auto readFraction=[&]() {
      auto value=in.fraction() ;
      payload.retain(integerBytes(value.numerator())) ;
      payload.retain(integerBytes(value.denominator())) ;
      return value ;
    } ;
    const auto ids=[&](std::size_t limit) {
      const auto size=allocation(in.count(limit,8),sizeof(BoundaryId)) ;
      std::vector<BoundaryId> values(size) ;
      for(auto &value:values) value=readId(in) ;
      return values ;
    } ;
    result.points.resize(allocation(in.count(limits.maxEntities,24),sizeof(Point))) ;
    count(result.points.size(),entities,limits.maxEntities) ;
    for(auto &point:result.points) for(double &coordinate:point) coordinate=in.real() ;
    result.edges.resize(allocation(in.count(limits.maxEntities-entities,48),sizeof(EdgeRecord))) ;
    count(result.edges.size(),entities,limits.maxEntities) ;
    for(auto &edge:result.edges) {
      for(auto &vertex:edge.vertices) vertex=readId(in) ;
      edge.geometry=readId(in) ; edge.parent=readId(in) ;
      edge.cuts.resize(allocation(in.count(limits.maxSplits-splits,27),sizeof(EdgeCut))) ;
      count(edge.cuts.size(),splits,limits.maxSplits) ;
      for(auto &cut:edge.cuts) { cut.position=readFraction() ; cut.vertex=readId(in) ; }
      edge.children.resize(allocation(in.count(result.edges.size()-childEdges,46),sizeof(EdgeInterval))) ;
      count(edge.children.size(),childEdges,result.edges.size()) ;
      for(auto &child:edge.children) { child.begin=readFraction() ; child.end=readFraction() ; child.edge=readId(in) ; }
    }
    result.faces.resize(allocation(in.count(limits.maxEntities-entities,52),sizeof(FaceRecord))) ;
    for(auto &face:result.faces) {
      face.parent=readId(in) ; face.parentSplit=readId(in) ; face.geometry=readId(in) ;
      const auto label=in.u32() ;
      face.boundary=static_cast<int>(label<0x80000000u ? std::int64_t(label) : std::int64_t(label)-0x100000000LL) ;
      face.vertices=ids(limits.maxEntities) ; face.edges=ids(limits.maxEntities) ;
      face.splits.resize(allocation(in.count(limits.maxSplits-splits,57),sizeof(FaceSplitRecord))) ;
      count(face.splits.size(),splits,limits.maxSplits) ;
      for(auto &split:face.splits) {
        split.centerMethod=readCenter(in) ;
        for(double &coordinate:split.referenceCenter) coordinate=in.real() ;
        split.centerVertex=readId(in) ; split.children=ids(result.faces.size()-childFaces) ;
        count(split.children.size(),childFaces,result.faces.size()) ;
      }
    }
    in.finish() ; validateBoundaries(result,limits) ;
    return result ;
  }

}
}
