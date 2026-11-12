// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
// Compact hex state, record validation and temporary working geometry.
#include "compact_hex.h"
#include <FVMAdapt/parametric/geometry_quality.h>
#include "library/parametric/record_io.h"
#include "library/parametric/face_plan_io.h"
#include "library/parametric/geometry_validity.h"

namespace Loci {
namespace FVMAdapt {
  namespace {
    constexpr std::uint32_t magic = 0x48455843 ;
    void require(bool ok, const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    void bytes(detail::RecordWriter &out, const ParametricRecord &record) { out.u64(record.size()) ; out.append(record) ; }
    ParametricRecord bytes(detail::RecordReader &in, std::size_t limit) {
      ParametricRecord record(in.count(limit)) ; for(auto &b : record) b = in.byte() ; return record ;
    }
    RecordLimits sourceLimits(const HistoryLimits &limits) {
      RecordLimits result ; result.maxBytes=limits.exact.maxBytes ; result.geometry=limits.geometry ;
      result.geometry.exact=limits.exact ; result.jacobian=limits.jacobian ; return result ;
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

    constexpr unsigned absent = std::numeric_limits<unsigned>::max() ;
    Point midpoint(const Point &a, const Point &b) {
      Point p ; for(unsigned d = 0; d < 3; ++d) p[d] = .5*a[d]+.5*b[d] ; return p ;
    }
    void finite(const Point &p) {
      for(double x : p) require(std::isfinite(x),"nonfinite shared boundary point") ;
    }

    struct FaceNode {
      Point center{{0,0,0}} ;
      std::array<unsigned,4> children{{absent,absent,absent,absent}} ;
      unsigned parent = absent, corner = 0 ;
      bool split() const { return children[0] != absent ; }
    } ;
    // A received face tree borrows its root map; child maps are evaluated along the path.
    struct FaceTree {
      const GeometryMap &source ;
      const HexPlan &plan ;
      HexPlan::Cell owner ;
      unsigned face ;
      std::shared_ptr<const GeometryMap> original ;
      std::vector<FaceNode> nodes{{}} ;

      Point position(unsigned id, Point p) const {
        while(nodes[id].parent != absent) {
          const auto &node = nodes[id] ; const auto &parent = nodes[node.parent] ; const auto v = node.corner ;
          p = quadChildMap(v,parent.center).position(p) ;
          id = node.parent ;
        }
        if(original) return original->position(p) ;
        return HexCellMap(source,plan,owner).position(boundaryPoint(Shape::hex,face,p)) ;
      }
      void readEdge(detail::FacePlanReader &in, unsigned id, Point a, Point b, unsigned depth) {
        if(!in.edge(depth)) return ;
        const auto p = midpoint(a,b) ; finite(position(id,p)) ;
        readEdge(in,id,a,p,depth+1) ; readEdge(in,id,p,b,depth+1) ;
      }
      void readFace(detail::FacePlanReader &in, unsigned id, unsigned depth) {
        const auto division = in.face(4,depth) ; const auto &quad = referenceElement(Shape::quad) ;
        for(unsigned e = 0; e < 4; ++e) readEdge(in,id,quad.vertices[e],quad.vertices[(e+1)%4],depth) ;
        if(!division.split) return ;
        nodes[id].center = division.center ; finite(position(id,division.center)) ;
        for(unsigned v = 0; v < 4; ++v) {
          finite(position(id,midpoint(quad.vertices[v],quad.vertices[(v+1)%4]))) ;
          const unsigned child = nodes.size() ; nodes[id].children[v] = child ;
          FaceNode node ; node.parent = id ; node.corner = v ; nodes.push_back(node) ;
          readFace(in,child,depth+1) ;
        }
      }
      void read(const FacePlan &record, const BoundaryRecordLimits &limits) {
        if(record.empty()) return ;
        detail::FacePlanReader in(record,limits) ; require(in.vertices == 4,"shared face plan shapes differ") ;
        readFace(in,0,0) ; in.input.finish() ;
      }
    } ;
    struct PlanFace {
      unsigned tree = absent, node = 0, orientation = 0 ;
    } ;
    using CellFaces = std::array<PlanFace,6> ;
    using FaceCorners = std::array<unsigned,4> ;

    const VertexVolumeSplit &pattern() {
      static const auto split = splitVertexVolume(Shape::hex,std::vector<Point>(6,{{.5,.5,0}}),{{.5,.5,.5}}) ; return split ;
    }
    FaceCorners childFace(unsigned child, unsigned face) {
      const auto &boundary = referenceElement(Shape::hex).boundaries[face] ; FaceCorners vertices ;
      for(unsigned v = 0; v < 4; ++v) vertices[v] = pattern().children[child].vertices[boundary.vertices[v]] ;
      return vertices ;
    }
    unsigned orientation(const FaceCorners &local, const FaceCorners &canonical) {
      for(unsigned o = 0; o < 8; ++o) {
        bool match = true ;
        for(unsigned v = 0; v < 4; ++v) match = match && local[v] == canonical[orientedVertex(Shape::quad,o,v)] ;
        if(match) return o ;
      }
      throw std::logic_error("invalid hex split face orientation") ;
    }
    unsigned edgePoint(int a, int b) {
      const auto &hex = referenceElement(Shape::hex) ;
      for(unsigned e = 0; e < 12; ++e)
        if((hex.edges[e][0] == a && hex.edges[e][1] == b) || (hex.edges[e][0] == b && hex.edges[e][1] == a)) return 8+e ;
      throw std::logic_error("invalid hex split edge") ;
    }
    FaceCorners boundaryChild(unsigned face, unsigned o, unsigned corner) {
      FaceCorners vertices ; const auto &boundary = referenceElement(Shape::hex).boundaries[face] ;
      for(unsigned v = 0; v < 4; ++v) vertices[orientedVertex(Shape::quad,o,v)] = boundary.vertices[v] ;
      const auto a = vertices[corner] ;
      return {{a,edgePoint(a,vertices[(corner+1)%4]),20+face,edgePoint(a,vertices[(corner+3)%4])}} ;
    }
    const std::array<InternalFace,12> &interiorOrder() {
      // Records number interior roots by first appearance in child-face order.
      static const auto faces = [] {
        std::array<InternalFace,12> result ; unsigned next = 0 ;
        for(unsigned c = 0; c < 8; ++c) for(unsigned f = 0; f < 6; ++f)
          for(const auto &pair : pattern().internalFaces)
            if(pair.first.child == c && pair.first.boundary == f) result[next++] = pair ;
        require(next == 12,"hex construction requires twelve interior faces") ; return result ;
      }() ; return faces ;
    }
  }

  CompactHexCell::CompactHexCell(CellTopology t,std::vector<Point> p,HistoryLimits h,std::shared_ptr<const CubicHexMap> geometry)
    : topology(std::move(t)),corners(std::move(p)),cubic(std::move(geometry)),limits(std::move(h)) {
    require(topology.shape()==Shape::hex && corners.size()==8,"compact source is not a hex") ;
    if(cubic) checkParametricSource(topology,corners,*cubic,{},limits.jacobian) ;
  }
  std::shared_ptr<const GeometryMap> CompactHexCell::sourceMap() const {
    return cubic ? std::shared_ptr<const GeometryMap>(cubic):std::make_shared<CornerMap>(Shape::hex,corners) ;
  }
  std::shared_ptr<const GeometryMap> CompactHexCell::faceMap(unsigned face) const {
    return sourceFaceGeometry(topology,face,corners,cubic.get()) ;
  }

  CompactHexCell::CompactHexCell(const ParametricCorner &cell) : topology(cell.topology), limits(cell.subdivision.history().limits()) {
    require(topology.shape() == Shape::hex,"compact hex history requires a hex source") ;
    const auto &history = cell.subdivision.history() ;
    const auto source=history.geometry().map(history.node(0).map) ;
    const auto *linear=dynamic_cast<const CornerMap *>(source.get()) ; cubic=std::dynamic_pointer_cast<const CubicHexMap>(source) ;
    require((linear || cubic) && source->shape()==Shape::hex && history.geometryRevision()==0,"compact hex requires original source geometry") ;
    if(linear) corners=linear->corners() ;
    else for(unsigned v=0;v<8;++v) {
      const auto &edge=cubic->edges()[v] ; corners.push_back(edge.edge->controls()[edge.reversed ? 3:0]) ;
    }
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
  ParametricCorner CompactHexCell::expand() const {
    ParametricCorner cell(topology,corners,limits,cubic) ;
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
  ParametricRecord CompactHexCell::record() const {
    detail::RecordWriter out(limits.exact.maxBytes,limits.exact) ; out.u32(magic) ; out.u32(cubic ? 2:1) ;
    const auto t = encodeParametricTopology(topology) ; out.u64(t.size()) ; for(int v : t) out.u32(v) ;
    for(const auto &p : corners) for(double x : p) out.real(x) ;
    if(cubic) {
      const auto bounds=sourceLimits(limits) ; GeometryStore geometry(bounds.geometry) ; const auto root=geometry.capture(cubic) ;
      std::vector<MapId> ids ; const auto record=encodeGeometry(geometry,ids,bounds) ;
      out.u64(ids.at(root)) ; bytes(out,record) ;
    }
    bytes(out,plan.encode()) ; for(const auto &face : faces) bytes(out,face) ;
    out.u64(interior.size()) ;
    for(const auto &face : interior) { out.u32(face.split) ; out.byte(face.face) ; bytes(out,face.plan) ; }
    return out.take() ;
  }
  bool CompactHexCell::saved(const ParametricRecord &record) {
    if(record.size() < 4) return false ;
    unsigned tag = 0 ; for(unsigned b = 0; b < 4; ++b) tag |= unsigned(record[b])<<(8*b) ; return tag == magic ;
  }
  CompactHexCell CompactHexCell::decode(const ParametricRecord &record, HistoryLimits limits) {
    detail::RecordReader in(record,limits.exact.maxBytes,limits.exact) ;
    require(in.u32()==magic,"unsupported compact hex record") ; const auto version=in.u32() ;
    require(version==1 || version==2,"unsupported compact hex record") ;
    std::vector<int> t(in.count(128,4)) ; for(auto &v : t) v = in.u32() ;
    const auto topology = decodeParametricTopology(t) ; require(topology.shape() == Shape::hex,"compact source is not a hex") ;
    std::vector<Point> corners(8) ; for(auto &p : corners) for(auto &x : p) x = in.real() ;
    std::shared_ptr<const CubicHexMap> cubic ;
    if(version==2) {
      const auto root=in.u64() ; const auto geometry=decodeGeometry(bytes(in,limits.exact.maxBytes),sourceLimits(limits)) ;
      require(root<geometry.size(),"invalid compact hex source map") ;
      cubic=std::dynamic_pointer_cast<const CubicHexMap>(geometry.map(root)) ; require(bool(cubic),"compact hex source is not cubic") ;
    } else require(validParametricSource(CornerMap(Shape::hex,corners),limits.jacobian),"invalid compact hex source geometry") ;
    CompactHexCell result(topology,std::move(corners),limits,std::move(cubic)) ;
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

  void CompactHexCell::validate() const {
    const auto source=sourceMap() ;
    BoundaryRecordLimits boundaryLimits ; boundaryLimits.maxSplits = boundaryLimits.maxEntities ;
    boundaryLimits.maxBytes = limits.exact.maxBytes ;
    std::vector<FaceTree> trees ; trees.reserve(6+12*plan.splitCount()) ;
    std::vector<CellFaces> cells(1+8*plan.splitCount()) ;
    for(unsigned f = 0; f < 6; ++f) {
      trees.push_back({*source,plan,0,f,faceMap(f)}) ; trees.back().read(faces[f],boundaryLimits) ;
      cells[0][f] = {f,0,topology.faceOrientation(f)} ;
    }
    std::size_t next = 0 ;
    for(std::size_t i = 0; i < plan.splitCount(); ++i) {
      const auto parent = plan.splitParent(i) ; const auto centers = plan.centers(parent) ; const auto children = plan.children(parent) ;
      for(unsigned f = 0; f < 6; ++f) {
        const auto use = cells[parent][f] ; const auto &node = trees.at(use.tree).nodes.at(use.node) ;
        require(node.split(),"accepted hex face choice is missing") ;
        const auto p = orientPoint(Shape::quad,inverseOrientation(Shape::quad,use.orientation),node.center) ;
        require(centers.faces[f] == std::array<double,2>{{p[0],p[1]}},"hex and face plan centers disagree") ;
      }
      for(const auto &piece : pattern().boundaryPieces) {
        const auto use = cells[parent][piece.boundary] ;
        const auto corner = orientedVertex(Shape::quad,use.orientation,piece.corner) ;
        const auto node = trees[use.tree].nodes[use.node].children[corner] ;
        const auto o = orientation(childFace(piece.face.child,piece.face.boundary),boundaryChild(piece.boundary,use.orientation,corner)) ;
        cells[children[piece.face.child]][piece.face.boundary] = {use.tree,node,o} ;
      }
      for(unsigned f = 0; f < 12; ++f) {
        const auto &pair = interiorOrder()[f] ; const unsigned tree = trees.size() ;
        trees.push_back({*source,plan,children[pair.first.child],pair.first.boundary,{}}) ;
        if(next < interior.size() && interior[next].split == i && interior[next].face == f)
          trees.back().read(interior[next++].plan,boundaryLimits) ;
        cells[children[pair.first.child]][pair.first.boundary] = {tree,0,0} ;
        cells[children[pair.second.child]][pair.second.boundary] = {tree,0,inverseOrientation(Shape::quad,pair.orientation)} ;
      }
    }
    require(next == interior.size(),"unused hex interior face plan") ;
    // Check retained (including dormant) boundaries and the selected physical leaves.
    const std::array<Point,5> probes{{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{.5,.5,0}}}} ;
    bool sourceProved=false ;
    for(HexPlan::Cell id = 0; id < cells.size(); ++id) {
      const HexCellMap map(*source,plan,id) ;
      const auto &hex = referenceElement(Shape::hex) ;
      const double extent = detail::pointExtent(8,[&](std::size_t v) { return map.position(hex.vertices[v]) ; }) ;
      for(unsigned f = 0; f < 6; ++f) {
        const auto use = cells[id][f] ;
        for(const auto &p : probes)
          require(detail::mappedPointAgrees(map.position(boundaryPoint(Shape::hex,f,p)),
            trees.at(use.tree).position(use.node,orientPoint(Shape::quad,use.orientation,p)),extent),
            "saved shared geometry disagrees with cell geometry") ;
      }
      if(plan.active(id) && plan.leaf(id)) {
        const auto quality=geometryQuality(map) ;
        require(!quality.invalidEvaluation,"invalid compact hex leaf evaluation") ;
        if(quality.status!=GeometryQualityStatus::sampled) {
          if(!sourceProved) sourceProved=geometryValidity(*source,limits.jacobian)==GeometryValidity::valid ;
          require(sourceProved,"invalid compact hex leaf geometry") ;
        }
      }
    }
  }
}
}
