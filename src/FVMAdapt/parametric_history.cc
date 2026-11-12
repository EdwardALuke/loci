// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.

#include <FVMAdapt/parametric_cell.h>
#include "parametric_hex.h"
#include <FVMAdapt/parametric/geometry_quality.h>
#include <FVMAdapt/parametric/planar_subdivision.h>
#include "library/parametric/record_io.h"
#include <algorithm>
#include <cmath>
#include <set>

namespace Loci {
namespace FVMAdapt {
  namespace {
    RecordLimits recordLimits(const HistoryLimits &history) {
      RecordLimits limits ; limits.maxBytes = history.exact.maxBytes ;
      limits.maxNodes = history.maxNodes ; limits.maxSplits = history.maxSplits ;
      limits.maxDepth = history.maxDepth ; limits.geometry = history.geometry ;
      limits.geometry.exact = history.exact ; limits.jacobian = history.jacobian ;
      return limits ;
    }
    TreeNodeId requestedLeaf(const RefinementHistory &history, const CellKey &key) {
      const auto id = history.nodeId(key) ;
      if(!history.active(id) || !history.node(id).isLeaf())
        throw std::invalid_argument("refinement request is not an active leaf") ;
      return id ;
    }
    void writeConstruction(detail::RecordWriter &out, const SplitConstruction &split) {
      out.u32(unsigned(split.parentShape)) ; out.u32(unsigned(split.family)) ; out.u32(split.version) ;
      out.byte(split.vertexVolume ? 2 : (split.planar ? 1 : 0)) ;
      if(split.vertexVolume) {
        const auto &p=*split.vertexVolume ; out.u64(p.points.size()) ;
        for(const auto &point:p.points) for(double x:point) out.real(x) ;
        out.u64(p.children.size()) ; for(const auto &child:p.children) for(unsigned v:child.vertices) out.u32(v) ; return ;
      }
      if(split.planar) {
        const auto &p=*split.planar ; out.u32(unsigned(p.centerMethod)) ; out.u64(p.points.size()) ;
        for(const auto &point:p.points) for(double x:point) out.real(x) ;
        out.u64(p.children.size()) ; for(const auto &child:p.children) for(unsigned v:child) out.u32(v) ;
        return ;
      }
      out.u64(split.children.size()) ;
      for(const auto &child: split.children) {
        out.u32(unsigned(child.shape)) ;
        for(const auto &q: child.placement.origin) out.fraction(q) ;
        for(const auto &column: child.placement.columns) for(const auto &q: column) out.fraction(q) ;
      }
    }
    // Encode validated requests in persistent-key order, sharing repeated constructions.
    ParametricRecord requestRecord(EntityId root,std::uint64_t revision,const HistoryLimits &limits,
        const std::vector<ParametricRequest> &requests) {
      if(requests.empty()) return {} ;
      std::map<CellKey,const ParametricRequest *> keys ; std::map<SplitConstruction,unsigned> constructions ;
      for(const auto &request:requests) {
        if(!keys.emplace(request.cell,&request).second) throw std::invalid_argument("duplicate refinement request") ;
        for(const auto &step:request.cell.lineage) constructions.emplace(*step.construction,0) ;
      }
      detail::RecordWriter out(limits.exact.maxBytes,limits.exact) ;
      out.u32(3) ; out.u64(root) ; out.u64(revision) ;
      out.u64(constructions.size()) ; unsigned index=0 ;
      for(auto &entry:constructions) { entry.second=index++ ; writeConstruction(out,entry.first) ; }
      out.u64(keys.size()) ;
      for(const auto &entry:keys) {
        const auto &key=entry.first ; out.u64(key.lineage.size()) ;
        for(const auto &step:key.lineage) { out.u32(constructions.at(*step.construction)) ; out.u32(step.child) ; }
        out.byte(unsigned(entry.second->action)) ;
        if(entry.second->action==ParametricAction::split) detail::writeSplitOperation(out,entry.second->split) ;
      }
      return out.take() ;
    }
    SplitConstruction readConstruction(detail::RecordReader &in, const HistoryLimits &limits) {
      SplitConstruction split ; split.parentShape = Shape(in.u32()) ;
      split.family = SplitFamily(in.u32()) ; split.version = in.u32() ;
      const unsigned kind=in.byte() ;
      if(kind>2) throw std::invalid_argument("unknown request construction kind") ;
      if(kind==2) {
        std::vector<Point> points(in.count(27,24)) ; for(auto &p:points) for(double &x:p) x=in.real() ;
        std::vector<std::array<unsigned,8>> children(in.count(8,32)) ; for(auto &c:children) for(auto &v:c) v=in.u32() ;
        split.vertexVolume=std::make_shared<const VertexVolumeSplit>(restoreVertexVolumeSplit(split.parentShape,std::move(points),children,limits.jacobian)) ;
        return split ;
      }
      if(kind==1) {
        PlanarSplit p ; p.shape=split.parentShape ; p.centerMethod=CenterMethod(in.u32()) ;
        if(p.shape!=Shape::triangle && p.shape!=Shape::quad) throw std::invalid_argument("invalid planar request shape") ;
        if(p.centerMethod!=CenterMethod::reference && p.centerMethod!=CenterMethod::weighted && p.centerMethod!=CenterMethod::integrated)
          throw std::invalid_argument("invalid planar request center") ;
        p.points.resize(in.count(9,24)) ; for(auto &point:p.points) for(double &x:point) x=in.real() ;
        p.children.resize(in.count(4,16)) ; for(auto &child:p.children) for(unsigned &v:child) v=in.u32() ;
        split.planar=std::make_shared<const PlanarSplit>(std::move(p)) ; return split ;
      }
      split.children.resize(in.count(limits.maxNodes,4)) ;
      for(auto &child: split.children) {
        child.shape = Shape(in.u32()) ;
        for(auto &q: child.placement.origin) q = in.fraction() ;
        for(auto &column: child.placement.columns) for(auto &q: column) q = in.fraction() ;
      }
      return split ;
    }
    // Decode the shared request format; each representation resolves its own cell keys.
    template<class Resolve,class Parent>
    std::vector<ParametricLeafRequest> readRequests(const ParametricRecord &record,EntityId root,std::uint64_t revision,
        const HistoryLimits &limits,Resolve resolve,Parent parent) {
      if(record.empty()) return {} ;
      detail::RecordReader in(record,limits.exact.maxBytes,limits.exact) ;
      if(in.u32() != 3 || in.u64() != root || in.u64() != revision)
        throw std::invalid_argument("parametric request version, root or geometry revision differs") ;
      std::vector<std::shared_ptr<const SplitConstruction>> constructions ;
      const auto count = in.count(limits.maxSplits,20) ;
      for(std::size_t i = 0; i < count; ++i)
        constructions.push_back(std::make_shared<SplitConstruction>(readConstruction(in,limits))) ;
      std::vector<ParametricLeafRequest> result ; std::set<TreeNodeId> unique ;
      const auto requests = in.count(limits.maxNodes,8) ;
      for(std::size_t i = 0; i < requests; ++i) {
        CellKey key{root,{}} ; const auto depth = in.count(limits.maxDepth,8) ;
        for(std::size_t d = 0; d < depth; ++d) {
          const auto construction = in.u32(), child = in.u32() ;
          if(construction >= constructions.size() || child >= constructions[construction]->childCount())
            throw std::invalid_argument("invalid parametric request lineage") ;
          key.lineage.push_back({constructions[construction],child}) ;
        }
        const unsigned action=in.byte() ;
        if(action>2) throw std::invalid_argument("unknown parametric request action") ;
        const auto id=resolve(key,ParametricAction(action)) ;
        if(!unique.insert(id).second) throw std::invalid_argument("duplicate refinement request") ;
        result.push_back({id,action==0 ? detail::readSplitOperation(in,limits.maxNodes):SplitOperation(TensorSplit()),ParametricAction(action)}) ;
      }
      in.finish() ;
      for(auto id:unique) for(auto ancestor=parent(id);ancestor!=noTreeNode;ancestor=parent(ancestor))
        if(unique.count(ancestor)) throw std::invalid_argument("overlapping ancestor and descendant requests") ;
      return result ;
    }
    std::shared_ptr<const SplitConstruction> hexConstruction(const HexPlan &plan,HexPlan::Cell parent) {
      SplitConstruction split ; split.parentShape=Shape::hex ; split.family=SplitFamily::vertexVolume ;
      split.vertexVolume=std::make_shared<const VertexVolumeSplit>(plan.construction(parent)) ;
      return std::make_shared<const SplitConstruction>(std::move(split)) ;
    }
  }

  bool validParametricSource(const CornerMap &geometry,JacobianCheckLimits limits) {
    if(geometry.dimension()!=3 || limits.maxBoxes==0) throw std::invalid_argument("invalid source cell shape or Jacobian bounds") ;
    const auto quality=geometryQuality(geometry) ;
    if(quality.status==GeometryQualityStatus::sampled) return true ;
    if(quality.invalidEvaluation || geometryValidity(geometry,limits)!=GeometryValidity::valid) return false ;
    const auto &element=referenceElement(geometry.shape()) ; const auto &points=geometry.corners() ;
    for(unsigned e=0;e<element.edgeCount;++e)
      if(points[element.edges[e][0]]==points[element.edges[e][1]]) return false ;
    for(unsigned f=0;f<element.boundaryCount;++f) {
      const auto &face=element.boundaries[f] ; std::vector<Point> corners ;
      for(unsigned v=0;v<face.vertexCount;++v) corners.push_back(points[face.vertices[v]]) ;
      if(geometryValidity(CornerMap(face.shape,std::move(corners)),limits)!=GeometryValidity::valid) return false ;
    }
    return true ;
  }

  namespace {
    std::array<unsigned,2> canonicalEdge(const CellTopology &topology,unsigned edge) {
      auto ends=topology.edgeVertices(edge) ;
      if(topology.vertexId(ends[1])<topology.vertexId(ends[0])) std::swap(ends[0],ends[1]) ;
      return ends ;
    }
    std::array<Point,2> linearEdgeControls(const Point &first,const Point &last) {
      std::array<Point,2> result ;
      for(unsigned d=0;d<3;++d) {
        const long double delta=static_cast<long double>(last[d])-first[d] ;
        result[0][d]=first[d]+delta/3 ; result[1][d]=first[d]+2*delta/3 ;
      }
      return result ;
    }
    bool suppliedCurves(const CellTopology &topology,const std::vector<Point> &corners,const ParametricEdgeControls &controls) {
      if(corners.size()!=topology.vertexCount() || (!controls.empty() && controls.size()!=topology.edgeCount()))
        throw std::invalid_argument("source corner or edge-control count differs from topology") ;
      bool curved=false ;
      for(unsigned e=0;e<controls.size();++e) if(!controls[e].empty()) {
        if(controls[e].size()!=2) throw std::invalid_argument("a cubic edge requires two interior controls") ;
        for(const auto &point:controls[e]) for(double x:point)
          if(!std::isfinite(x)) throw std::invalid_argument("cubic edge controls must be finite") ;
        const auto ends=canonicalEdge(topology,e) ; const auto linear=linearEdgeControls(corners[ends[0]],corners[ends[1]]) ;
        curved=curved || controls[e][0]!=linear[0] || controls[e][1]!=linear[1] ;
      }
      return curved ;
    }
    const Point &sourceEdgeControl(CubicEdgeView edges,unsigned edge,unsigned control) {
      const auto &use=edges[edge] ; return use.edge->controls()[use.reversed ? 3-control:control] ;
    }
    template<class Map,std::size_t N>
    std::shared_ptr<const GeometryMap> cubicSource(const CellTopology &topology,const std::vector<Point> &corners,
        const ParametricEdgeControls &controls) {
      std::array<CubicEdgeUse,N> edges ;
      for(unsigned e=0;e<N;++e) {
        const auto ends=canonicalEdge(topology,e) ; const auto linear=linearEdgeControls(corners[ends[0]],corners[ends[1]]) ;
        const auto &row=controls[e] ;
        std::vector<Point> points{corners[ends[0]],row.empty() ? linear[0]:row[0],row.empty() ? linear[1]:row[1],corners[ends[1]]} ;
        edges[e]={std::make_shared<BezierMap>(Shape::segment,std::array<unsigned,3>{3,0,0},std::move(points)),topology.edgeOrientation(e)!=0} ;
      }
      return std::make_shared<Map>(std::move(edges)) ;
    }
  }

  CubicEdgeView parametricCubicEdges(const GeometryMap &geometry) {
    if(const auto *tet=dynamic_cast<const CubicTetMap *>(&geometry)) return {tet->edges().data(),tet->edges().size()} ;
    if(const auto *prism=dynamic_cast<const CubicPrismMap *>(&geometry)) return {prism->edges().data(),prism->edges().size()} ;
    if(const auto *pyramid=dynamic_cast<const CubicPyramidMap *>(&geometry)) return {pyramid->edges().data(),pyramid->edges().size()} ;
    if(const auto *hex=dynamic_cast<const CubicHexMap *>(&geometry)) return {hex->edges().data(),hex->edges().size()} ;
    return {} ;
  }

  std::shared_ptr<const GeometryMap> parametricCubicGeometry(const CellTopology &topology,const std::vector<Point> &corners,
      const ParametricEdgeControls &controls) {
    if(!suppliedCurves(topology,corners,controls)) return {} ;
    if(topology.isStandard()) switch(topology.shape()) {
      case Shape::tet: return cubicSource<CubicTetMap,6>(topology,corners,controls) ;
      case Shape::prism: return cubicSource<CubicPrismMap,9>(topology,corners,controls) ;
      case Shape::pyramid: return cubicSource<CubicPyramidMap,8>(topology,corners,controls) ;
      case Shape::hex: return cubicSource<CubicHexMap,12>(topology,corners,controls) ;
      default: break ;
    }
    throw std::invalid_argument("curved source edges require a standard tet, prism, pyramid or hex") ;
  }

  ParametricEdgeControls parametricSourceEdgeControls(const CellTopology &topology,const GeometryMap &geometry) {
    if(!topology.isStandard() || topology.shape()!=geometry.shape())
      throw std::invalid_argument("source geometry shape differs from topology") ;
    ParametricEdgeControls controls(topology.edgeCount()) ;
    if(dynamic_cast<const CornerMap *>(&geometry)) return controls ;
    const auto edges=parametricCubicEdges(geometry) ;
    if(!edges.size) throw std::invalid_argument("unsupported parametric source geometry") ;
    for(unsigned e=0;e<edges.size;++e) {
      const bool reverse=topology.edgeOrientation(e)!=0 ;
      const auto control=[&](unsigned i)->const Point & { return sourceEdgeControl(edges,e,reverse ? 3-i:i) ; } ;
      const auto linear=linearEdgeControls(control(0),control(3)) ;
      if(control(1)!=linear[0] || control(2)!=linear[1]) controls[e]={control(1),control(2)} ;
    }
    return controls ;
  }

  void checkParametricSource(const CellTopology &topology,const std::vector<Point> &corners,const GeometryMap &geometry,
      const ParametricEdgeControls &controls,JacobianCheckLimits limits) {
    suppliedCurves(topology,corners,controls) ;
    if(!topology.isStandard() || topology.shape()!=geometry.shape())
      throw std::invalid_argument("source geometry shape differs from topology") ;
    const auto *corner=dynamic_cast<const CornerMap *>(&geometry) ;
    const auto edges=parametricCubicEdges(geometry) ;
    if(corner) {
      if(corner->corners()!=corners) throw std::invalid_argument("saved source corners differ") ;
    } else if(edges.size) {
      for(unsigned e=0;e<edges.size;++e) {
        const auto ends=topology.edgeVertices(e) ;
        if(sourceEdgeControl(edges,e,0)!=corners[ends[0]] || sourceEdgeControl(edges,e,3)!=corners[ends[1]])
          throw std::invalid_argument("saved cubic edge endpoints differ from source corners") ;
      }
    } else throw std::invalid_argument("unsupported parametric source geometry") ;
    if(!controls.empty()) {
      const auto accepted=parametricSourceEdgeControls(topology,geometry) ;
      for(unsigned e=0;e<controls.size();++e) if(!controls[e].empty()) {
        const auto ends=canonicalEdge(topology,e) ; const auto linear=linearEdgeControls(corners[ends[0]],corners[ends[1]]) ;
        if(accepted[e].empty() ? controls[e][0]!=linear[0] || controls[e][1]!=linear[1] : controls[e]!=accepted[e])
          throw std::invalid_argument("supplied cubic controls differ from saved source geometry") ;
      }
    }
    if(corner ? !validParametricSource(*corner,limits):!validGeometry(geometry,{},limits))
      throw std::invalid_argument("invalid or unresolved parametric root geometry") ;
  }

  bool sameParametricSource(const GeometryMap &first,const GeometryMap &second) {
    if(const auto *a=dynamic_cast<const CornerMap *>(&first)) {
      const auto *b=dynamic_cast<const CornerMap *>(&second) ;
      return b && a->shape()==b->shape() && a->corners()==b->corners() ;
    }
    const auto a=parametricCubicEdges(first),b=parametricCubicEdges(second) ;
    if(!a.size || a.size!=b.size || first.shape()!=second.shape()) return false ;
    for(unsigned e=0;e<a.size;++e) for(unsigned i=0;i<4;++i)
      if(sourceEdgeControl(a,e,i)!=sourceEdgeControl(b,e,i)) return false ;
    return true ;
  }

  HistoryLimits parametricHistoryLimits(int maxDepth, int maxNodes) {
    if(maxDepth < 0 || maxDepth > 30 || maxNodes <= 0)
      throw std::invalid_argument("parametric limits require depth in [0,30] and positive node count") ;
    HistoryLimits limits ; limits.maxDepth = maxDepth ; limits.maxNodes = maxNodes ;
    limits.maxSplits = maxNodes ; return limits ;
  }
  std::vector<unsigned char> encodeParametricHistory(const AffineSubdivision &history) {
    return encodeRefinement(history.history(),recordLimits(history.history().limits())) ;
  }
  AffineSubdivision decodeParametricHistory(const std::vector<unsigned char> &record, HistoryLimits limits) {
    return AffineSubdivision(decodeRefinement(record,recordLimits(limits))) ;
  }
  AffineSubdivision makeParametricHistory(const CellTopology &topology, const std::vector<Point> &corners,
      const std::vector<unsigned char> &record, HistoryLimits limits,const ParametricEdgeControls &controls) {
    if(!record.empty()) {
      auto saved = decodeParametricHistory(record,limits) ; const auto &history = saved.history() ;
      if(history.rootKey() != topology.id() || history.geometryRevision() != 0)
        throw std::invalid_argument("parametric history does not match the original cell geometry") ;
      checkParametricSource(topology,corners,*history.geometry().map(history.node(0).map),controls,limits.jacobian) ;
      return saved ;
    }
    GeometryStore geometry(limits.geometry) ;
    const auto cubic=parametricCubicGeometry(topology,corners,controls) ;
    const auto root = cubic ? geometry.capture(cubic):geometry.append(CornerRecord{topology.shape(),corners}) ;
    checkParametricSource(topology,corners,*geometry.map(root),{},limits.jacobian) ;
    return AffineSubdivision(RefinementHistory(topology.id(),std::move(geometry),root,0,limits)) ;
  }

  std::vector<unsigned char> initialParametricHistory(const CellTopology &topology, const std::vector<Point> &corners,
      const std::vector<unsigned char> &record, HistoryLimits limits,const ParametricEdgeControls &controls) {
    const auto cell=makeParametricHistory(topology,corners,record,limits,controls) ;
    return record.empty() ? encodeParametricHistory(cell):record ;
  }

  std::vector<unsigned char> encodeParametricRequests(const RefinementHistory &history, const std::vector<ParametricRequest> &requests) {
    if(requests.empty()) return {} ;
    for(const auto &request:requests) {
      if(request.action==ParametricAction::coarsen) {
        const auto id=history.nodeId(request.cell) ;
        if(!history.active(id) || history.node(id).isLeaf()) throw std::invalid_argument("coarsening request is not an active parent") ;
      } else requestedLeaf(history,request.cell) ;
    }
    return requestRecord(history.rootKey(),history.geometryRevision(),history.limits(),requests) ;
  }

  std::vector<unsigned char> encodeParametricRequests(const RefinementHistory &history, const std::vector<CellKey> &cells) {
    std::vector<ParametricRequest> requests ; requests.reserve(cells.size()) ;
    for(const auto &key: cells) {
      const auto id = requestedLeaf(history,key) ;
      requests.push_back({key,isotropicSplit(history.localCell(id).shape)}) ;
    }
    return encodeParametricRequests(history,requests) ;
  }
  std::vector<ParametricLeafRequest> decodeParametricRequests(const std::vector<unsigned char> &record, const RefinementHistory &history) {
    return readRequests(record,history.rootKey(),history.geometryRevision(),history.limits(),
      [&](const CellKey &key,ParametricAction action) {
        const bool coarsen=action==ParametricAction::coarsen ;
        const auto id=coarsen ? history.nodeId(key):requestedLeaf(history,key) ;
        if(!history.active(id) || (coarsen && history.node(id).isLeaf())) throw std::invalid_argument("coarsening request is not an active parent") ;
        return id ;
      },[&](TreeNodeId id) { return history.node(id).parent ; }) ;
  }

  namespace {
    void requestBytes(detail::RecordWriter &out,const ParametricRecord &bytes) { out.u64(bytes.size()) ; out.append(bytes) ; }
    ParametricRecord requestBytes(detail::RecordReader &in,std::size_t limit) {
      ParametricRecord result(in.count(limit)) ; for(auto &b:result) b=in.byte() ; return result ;
    }
    ParametricRecord boundarySource(const ParametricBoundary &cell) {
      detail::RecordWriter out(cell.limits().exact.maxBytes,cell.limits().exact) ;
      const auto topology=encodeParametricTopology(cell.topology) ; out.u64(topology.size()) ; for(int x:topology) out.u32(x) ;
      for(const auto &point:cell.sourceCorners()) for(auto x:point) out.real(x) ;
      if(const auto source=cell.sourceGeometry()) {
        GeometryStore geometry(cell.limits().geometry) ; geometry.capture(source) ;
        requestBytes(out,encodeGeometry(geometry,recordLimits(cell.limits()))) ;
      }
      return out.take() ;
    }
    void checkRequests(const ParametricCell &cell,const std::vector<ParametricLeafRequest> &requests) {
      std::set<TreeNodeId> unique ;
      for(const auto &r:requests) {
        if(r.action!=ParametricAction::split && r.action!=ParametricAction::isotropic && r.action!=ParametricAction::coarsen)
          throw std::invalid_argument("unknown parametric request action") ;
        const bool merge=r.action==ParametricAction::coarsen ;
        if(!cell.active(r.cell) || merge==cell.leaf(r.cell)) throw std::invalid_argument("request is not an active leaf or coarsening parent") ;
        if(!unique.insert(r.cell).second) throw std::invalid_argument("duplicate refinement request") ;
      }
      for(auto id:unique) for(auto parent=cell.parent(id);parent!=noTreeNode;parent=cell.parent(parent))
        if(unique.count(parent)) throw std::invalid_argument("overlapping ancestor and descendant requests") ;
    }
  }
  ParametricRecord encodeParametricRequests(const ParametricCell &cell,const std::vector<ParametricLeafRequest> &requests) {
    if(requests.empty()) return {} ;
    if(cell.source()) { auto expanded=cell ; expanded.expand() ; return encodeParametricRequests(expanded,requests) ; }
    checkRequests(cell,requests) ;
    if(!cell.boundary()) {
      std::vector<ParametricRequest> keys ;
      for(const auto &r:requests) keys.push_back({cell.history().key(r.cell),r.split,r.action}) ;
      return encodeParametricRequests(cell.history(),keys) ;
    }
    const auto &tree=*cell.boundary() ; const auto &limits=tree.limits() ;
    detail::RecordWriter out(limits.exact.maxBytes,limits.exact) ; out.u32(4) ; out.u64(cell.rootKey()) ;
    requestBytes(out,boundarySource(tree)) ;
    std::map<CellPath,const ParametricLeafRequest*> ordered ; for(const auto &r:requests) ordered.emplace(cell.path(r.cell),&r) ;
    out.u64(ordered.size()) ;
    for(const auto &entry:ordered) {
      out.u64(entry.first.size()) ; TreeNodeId parent=0 ;
      for(auto child:entry.first) {
        out.u32(child) ; requestBytes(out,tree.splitRecord(parent)) ; parent=tree.children(parent)[child] ;
      }
      const auto &r=*entry.second ; out.byte(unsigned(r.action)) ;
      if(r.action==ParametricAction::split) detail::writeSplitOperation(out,r.split) ;
    }
    return out.take() ;
  }
  std::vector<ParametricLeafRequest> decodeParametricRequests(const ParametricRecord &record,const ParametricCell &cell) {
    if(record.empty()) return {} ;
    if(cell.source()) { auto expanded=cell ; expanded.expand() ; return decodeParametricRequests(record,expanded) ; }
    if(!cell.boundary()) return decodeParametricRequests(record,cell.history()) ;
    const auto &tree=*cell.boundary() ; const auto &limits=tree.limits() ;
    detail::RecordReader in(record,limits.exact.maxBytes,limits.exact) ;
    if(in.u32()!=4 || in.u64()!=cell.rootKey() || requestBytes(in,limits.exact.maxBytes)!=boundarySource(tree))
      throw std::invalid_argument("parametric request version or source geometry differs") ;
    std::vector<ParametricLeafRequest> result ; const auto count=in.count(limits.maxNodes,9) ;
    for(std::size_t i=0;i<count;++i) {
      TreeNodeId id=0 ; const auto depth=in.count(limits.maxDepth,12) ;
      for(std::size_t d=0;d<depth;++d) {
        const auto child=in.u32() ; const auto split=requestBytes(in,limits.exact.maxBytes) ; const auto children=tree.children(id) ;
        if(child>=children.size() || tree.splitRecord(id)!=split) throw std::invalid_argument("parametric request accepted split differs") ;
        id=children[child] ;
      }
      const unsigned action=in.byte() ; if(action>2) throw std::invalid_argument("unknown parametric request action") ;
      result.push_back({id,action==0 ? detail::readSplitOperation(in,limits.maxNodes):SplitOperation(TensorSplit()),ParametricAction(action)}) ;
    }
    in.finish() ; checkRequests(cell,result) ; return result ;
  }
  std::vector<ParametricLeafRequest> hexRequests(const ParametricRecord &record,const HexPlan &plan,EntityId root,const HistoryLimits &limits) {
    std::map<HexPlan::Cell,std::shared_ptr<const SplitConstruction>> constructions ;
    return readRequests(record,root,0,limits,[&](const CellKey &key,ParametricAction action) {
      auto id=HexPlan::root ;
      for(const auto &step:key.lineage) {
        if(!plan.accepted(id)) throw std::out_of_range("cell key is not retained in this history") ;
        auto &construction=constructions[id] ; if(!construction) construction=hexConstruction(plan,id) ;
        if(!(*construction==*step.construction)) throw std::out_of_range("cell key is not retained in this history") ;
        id=plan.children(id).at(step.child) ;
      }
      if(!plan.active(id) || (action==ParametricAction::coarsen ? plan.leaf(id):!plan.leaf(id)))
        throw std::invalid_argument("request is not an active leaf or coarsening parent") ;
      return TreeNodeId(id) ;
    },[&](TreeNodeId id) { return id==HexPlan::root ? noTreeNode:TreeNodeId(plan.parent(id)) ; }) ;
  }
  ParametricRecord hexTagRequests(const HexPlan &plan,EntityId root,const HistoryLimits &limits,const std::vector<char> &values) {
    const auto leaves=plan.leaves() ;
    if(values.size()!=leaves.size()) throw std::invalid_argument("parametric tag count differs from accepted leaves") ;
    std::vector<std::pair<HexPlan::Cell,ParametricAction>> selected ;
    std::map<HexPlan::Cell,unsigned> coarsen ;
    for(unsigned i=0;i<leaves.size();++i) {
      const auto tag=values[i] ; if(tag<0 || tag>2) throw std::invalid_argument("unknown parametric refinement tag") ;
      if(tag==1) selected.emplace_back(leaves[i],ParametricAction::isotropic) ;
      if(tag==2 && leaves[i]!=HexPlan::root) ++coarsen[plan.parent(leaves[i])] ;
    }
    // Eight marked active leaves are exactly one complete sibling family.
    for(const auto &entry:coarsen) if(entry.second==8) selected.emplace_back(entry.first,ParametricAction::coarsen) ;
    std::map<HexPlan::Cell,std::shared_ptr<const SplitConstruction>> constructions ;
    std::vector<ParametricRequest> requests ; requests.reserve(selected.size()) ;
    for(const auto &entry:selected) {
      CellKey key{root,{}} ;
      for(auto cell=entry.first;cell!=HexPlan::root;) {
        const auto parent=plan.parent(cell) ; auto &construction=constructions[parent] ;
        if(!construction) construction=hexConstruction(plan,parent) ;
        const auto children=plan.children(parent) ;
        key.lineage.push_back({construction,unsigned(cell-children.front())}) ; cell=parent ;
      }
      std::reverse(key.lineage.begin(),key.lineage.end()) ;
      requests.push_back({std::move(key),TensorSplit(),entry.second}) ;
    }
    return requestRecord(root,0,limits,requests) ;
  }

  ParametricRecord encodeParametricTags(const ParametricCell &cell,const std::vector<char> &values) {
    if(values.size()!=cell.leafCount()) throw std::invalid_argument("parametric tag count differs from accepted leaves") ;
    const auto leaves=cell.leafIds() ; std::map<TreeNodeId,int> tags ; std::set<TreeNodeId> parents ;
    std::vector<ParametricLeafRequest> requests ;
    for(unsigned i=0;i<leaves.size();++i) {
      const auto tag=values[i] ; if(tag<0 || tag>2) throw std::invalid_argument("unknown parametric refinement tag") ;
      tags.emplace(leaves[i],tag) ;
      if(tag==1) requests.push_back({leaves[i],TensorSplit(),ParametricAction::isotropic}) ;
      if(tag==2 && cell.parent(leaves[i])!=noTreeNode) parents.insert(cell.parent(leaves[i])) ;
    }
    for(auto parent:parents) {
      bool complete=true ;
      for(auto child:cell.children(parent)) if(!tags.count(child) || tags.at(child)!=2) complete=false ;
      if(complete) requests.push_back({parent,TensorSplit(),ParametricAction::coarsen}) ;
    }
    return encodeParametricRequests(cell,requests) ;
  }

  namespace {
    void require(bool condition,const char *message) { if(!condition) throw std::invalid_argument(message) ; }
  }
  std::vector<std::vector<int>> parametricCellDonors(const RefinementHistory &oldHistory,const RefinementHistory &newHistory) {
    require(oldHistory.rootKey()==newHistory.rootKey() && oldHistory.geometryRevision()==newHistory.geometryRevision(),
            "cell correspondence requires the same original geometry") ;
    require(sameParametricSource(*oldHistory.geometry().map(oldHistory.node(0).map),*newHistory.geometry().map(newHistory.node(0).map)),
            "cell correspondence requires the same source geometry") ;
    std::vector<CellKey> old ; for(auto id:oldHistory.leafIds()) old.push_back(oldHistory.key(id)) ; std::vector<std::vector<int>> result ;
    for(auto id:newHistory.leafIds()) {
      const auto next=newHistory.key(id) ; std::vector<int> donors ;
      for(unsigned i=0;i<old.size();++i) {
        const auto &a=old[i].lineage,&b=next.lineage ; const auto n=std::min(a.size(),b.size()) ;
        if(std::equal(a.begin(),a.begin()+n,b.begin(),[](const CellStep &x,const CellStep &y) {
          return x.child==y.child && *x.construction==*y.construction ;
        })) donors.push_back(i) ;
      }
      require(!donors.empty(),"new cell has no ancestral donors") ; result.push_back(std::move(donors)) ;
    }
    return result ;
  }
  std::vector<std::vector<int>> parametricCellDonors(const ParametricCell &oldCell,const ParametricCell &newCell) {
    if(oldCell.boundary() || newCell.boundary()) {
      require(newCell.boundary(),"accepted boundary history cannot change split families") ;
      const auto &next=*newCell.boundary() ;
      if(oldCell.source()) {
        require(!next.sourceGeometry() && oldCell.source()->topology==next.topology &&
                oldCell.source()->corners==next.sourceCorners(),"cell correspondence requires the same original geometry") ;
        return std::vector<std::vector<int>>(newCell.leafCount(),std::vector<int>{0}) ;
      }
      require(oldCell.boundary(),"cell correspondence split families differ") ; const auto &old=*oldCell.boundary() ;
      require(old.topology==next.topology && old.sourceCorners()==next.sourceCorners(),
              "cell correspondence requires the same original geometry") ;
      const auto oldSource=old.sourceGeometry(),newSource=next.sourceGeometry() ;
      require(bool(oldSource)==bool(newSource) && (!oldSource || sameParametricSource(*oldSource,*newSource)),
              "cell correspondence requires the same source geometry") ;
      for(const auto &split:old.splits()) {
        const auto parent=next.nodeId(old.path(split.parent)) ;
        require(next.node(parent).split!=noTreeNode && old.splitRecord(split.parent)==next.splitRecord(parent),"accepted boundary split differs") ;
      }
      std::map<CellPath,int> donors ; int index=0 ; for(auto id:old.leafIds()) donors.emplace(old.path(id),index++) ;
      std::vector<std::vector<int>> result ;
      for(auto id:next.leafIds()) {
        const auto path=next.path(id) ; auto ancestor=path ; std::vector<int> indices ;
        for(;;) {
          const auto found=donors.find(ancestor) ;
          if(found!=donors.end()) { indices.push_back(found->second) ; break ; }
          if(ancestor.empty()) break ;
          ancestor.pop_back() ;
        }
        if(indices.empty()) for(auto it=donors.lower_bound(path);it!=donors.end();++it) {
          if(it->first.size()<path.size() || !std::equal(path.begin(),path.end(),it->first.begin())) break ;
          indices.push_back(it->second) ;
        }
        require(!indices.empty(),"new boundary cell has no ancestral donors") ; result.push_back(std::move(indices)) ;
      }
      return result ;
    }
    if(!oldCell.source() && !newCell.source()) return parametricCellDonors(oldCell.history(),newCell.history()) ;
    const auto root=[](const ParametricCell &c) { return c.source() ? c.source()->topology.id():c.history().rootKey() ; } ;
    require(root(oldCell)==root(newCell),"cell correspondence requires the same original geometry") ;
    if(oldCell.source()) {
      const auto &source=*oldCell.source() ;
      if(newCell.source()) require(source.topology==newCell.source()->topology && source.corners==newCell.source()->corners,
                                  "cell correspondence requires the same original geometry") ;
      else require(sameParametricSource(CornerMap(source.topology.shape(),source.corners),
                    *newCell.history().geometry().map(newCell.history().node(0).map)),
                   "cell correspondence requires the same source geometry") ;
      return std::vector<std::vector<int>>(newCell.leafCount(),std::vector<int>{0}) ;
    }
    throw std::invalid_argument("accepted refinement history cannot be replaced by an untouched source") ;
  }
}
}
