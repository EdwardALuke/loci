// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/parametric/hex_plan.h>
#include <FVMAdapt/parametric/geometry_quality.h>
#include "record_io.h"
#include <algorithm>
#include <cmath>

namespace Loci {
namespace FVMAdapt {
  namespace {
    void require(bool ok, const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    std::vector<Point> facePoints(const HexCenters &centers) {
      std::vector<Point> result ; result.reserve(6) ;
      for(const auto &p: centers.faces) result.push_back({{p[0],p[1],0}}) ;
      return result ;
    }
    const VertexVolumeSplit &pattern() {
      static const auto value = splitVertexVolume(Shape::hex,facePoints(HexCenters()),{{.5,.5,.5}}) ;
      return value ;
    }
    MapValue midpointValue(unsigned ordinal, const Point &point) {
      static const auto placements = [] {
        std::array<MapValue,8> result ;
        for(unsigned i = 0; i < 8; ++i) result[i] = pattern().children[i].placement->evaluate({{0,0,0}}) ;
        return result ;
      }() ;
      auto result = placements[ordinal] ;
      for(unsigned d = 0; d < 3; ++d) for(unsigned j = 0; j < 3; ++j)
        result.position[d] += result.jacobian[j][d]*point[j] ;
      return result ;
    }
    // Child corner order is fixed even when the accepted centers are off midpoint.
    const std::array<std::array<unsigned,8>,8> &childOrientations() {
      static const auto result = [] {
        std::array<std::array<unsigned,8>,8> order{} ;
        const auto &e = referenceElement(Shape::hex) ;
        for(unsigned c = 0; c < 8; ++c) for(unsigned v = 0; v < 8; ++v) {
          auto p = midpointValue(c,e.vertices[v]).position ;
          for(unsigned d = 0; d < 3; ++d) p[d] = 2*p[d]-e.vertices[c][d] ;
          const auto found = std::find(e.vertices.begin(),e.vertices.begin()+8,p) ;
          require(found != e.vertices.begin()+8,"invalid hex child orientation") ;
          order[c][v] = found-e.vertices.begin() ;
        }
        return order ;
      }() ;
      return result ;
    }
    struct LogicalBox {
      std::array<std::uint32_t,3> lower{} ;
      std::uint32_t width = HexPlan::extent ;
      std::array<unsigned,8> orientation{{0,1,2,3,4,5,6,7}} ;
      LogicalBox child(unsigned ordinal) const {
        LogicalBox result ; result.width = width/2 ;
        const auto &corner = referenceElement(Shape::hex).vertices[orientation[ordinal]] ;
        for(unsigned d = 0; d < 3; ++d) result.lower[d] = lower[d]+unsigned(corner[d])*result.width ;
        for(unsigned v = 0; v < 8; ++v) result.orientation[v] = orientation[childOrientations()[ordinal][v]] ;
        return result ;
      }
      unsigned childAt(const std::array<std::uint32_t,3> &point) const {
        const auto &e = referenceElement(Shape::hex) ;
        for(unsigned c = 0; c < 8; ++c) {
          bool match = true ;
          for(unsigned d = 0; d < 3; ++d)
            match = match && (point[d] >= lower[d]+width/2) == bool(e.vertices[orientation[c]][d]) ;
          if(match) return c ;
        }
        throw std::logic_error("invalid hex adjacency point") ;
      }
    } ;
    bool exactSum(double a, double b, double &sum) {
      sum = a+b ; const double part = sum-a ;
      if(!std::isfinite(sum) || !std::isfinite(part)) return false ;
      return (a-(sum-part))+(b-part) == 0 ;
    }
    bool affineCorners(const GeometryMap &map) {
      const auto *corner = dynamic_cast<const CornerMap*>(&map) ;
      if(!corner || corner->shape() != Shape::hex) return false ;
      const auto &p = corner->corners() ; const auto &e = referenceElement(Shape::hex) ;
      std::array<Point,3> edges ; const unsigned ends[3] = {1,3,4} ;
      // Require exact edge differences and sums; uncertain cases use the full geometry checks.
      for(unsigned j = 0; j < 3; ++j) for(unsigned d = 0; d < 3; ++d)
        if(!exactSum(p[ends[j]][d],-p[0][d],edges[j][d])) return false ;
      for(unsigned v = 0; v < 8; ++v) for(unsigned d = 0; d < 3; ++d) {
        double value = p[0][d] ;
        for(unsigned j = 0; j < 3; ++j)
          if(e.vertices[v][j] == 1 && !exactSum(value,edges[j][d],value)) return false ;
        if(value != p[v][d]) return false ;
      }
      return true ;
    }
    bool affineHex(const GeometryMap &map) {
      if(const auto *cell = dynamic_cast<const HexCellMap*>(&map)) return cell->affine() ;
      return affineCorners(map) ;
    }
    Point facePoint(unsigned face, const std::array<double,2> &point) {
      const auto &e = referenceElement(Shape::hex) ; const auto &f = e.boundaries[face] ;
      Point result = e.vertices[f.vertices[0]] ;
      for(unsigned d = 0; d < 3; ++d) for(unsigned j = 0; j < 2; ++j) {
        const auto direction = e.vertices[f.vertices[j == 0 ? 1 : 3]][d]-e.vertices[f.vertices[0]][d] ;
        if(direction) result[d] = direction > 0 ? point[j] : 1-point[j] ;
      }
      return result ;
    }
    std::array<Point,8> childCorners(const HexCenters &centers, unsigned ordinal) {
      const auto &fixed = pattern() ; std::array<Point,8> corners ;
      for(unsigned v = 0; v < 8; ++v) {
        const auto i = fixed.children[ordinal].vertices[v] ;
        if(i < 20) corners[v] = fixed.points[i] ;
        else if(i < 26) {
          corners[v] = facePoint(i-20,centers.faces[i-20]) ;
        } else corners[v] = centers.cell ;
      }
      return corners ;
    }
    void certify(const HexCenters &centers, JacobianCheckLimits limits) {
      require(limits.maxBoxes > 0,"hex certification requires a positive box limit") ;
      if(centers == HexCenters()) { pattern() ; return ; }
      validateVertexVolumeSplit(Shape::hex,facePoints(centers),centers.cell,limits) ;
    }
    MapValue compose(const MapValue &outer, const MapValue &inner) {
      MapValue result{outer.position,{},outer.hasJacobian && inner.hasJacobian} ;
      if(result.hasJacobian)
        for(unsigned d = 0; d < 3; ++d) for(unsigned c = 0; c < 3; ++c) for(unsigned j = 0; j < 3; ++j)
          result.jacobian[d][c] += outer.jacobian[j][c]*inner.jacobian[d][j] ;
      return result ;
    }
    class HexFaceMap final : public GeometryMap {
      const GeometryMap &cell_ ;
      unsigned face_ ;
    public:
      HexFaceMap(const GeometryMap &cell, unsigned face) : cell_(cell), face_(face) {}
      Shape shape() const override { return Shape::quad ; }
      MapValue evaluate(const Point &point) const override {
        const auto &element = referenceElement(Shape::hex) ; const auto &face = element.boundaries[face_] ;
        MapValue local{boundaryPoint(Shape::hex,face_,point),{},true} ;
        for(unsigned d = 0; d < 3; ++d) {
          local.jacobian[0][d] = element.vertices[face.vertices[1]][d]-element.vertices[face.vertices[0]][d] ;
          local.jacobian[1][d] = element.vertices[face.vertices[3]][d]-element.vertices[face.vertices[0]][d] ;
        }
        return compose(cell_.evaluate(local.position),local) ;
      }
    } ;
    class HexChildMap final : public GeometryMap {
      const GeometryMap &parent_ ;
      std::array<Point,8> corners_ ;
      unsigned ordinal_ ;
      bool midpoint_ ;
    public:
      HexChildMap(const GeometryMap &parent, const HexCenters &centers, unsigned ordinal)
        : parent_(parent), corners_(childCorners(centers,ordinal)), ordinal_(ordinal), midpoint_(centers == HexCenters()) {}
      Shape shape() const override { return Shape::hex ; }
      MapValue evaluate(const Point &point) const override {
        const auto local = midpoint_ ? midpointValue(ordinal_,point) : evaluateHexCorners(corners_,point) ;
        return compose(parent_.evaluate(local.position),local) ;
      }
    } ;
    void requireQuality(const GeometryMap &map, bool affine = false) {
      static const GeometryQualityOptions constantJacobian{{{{.5,.5,.5}}},1} ;
      const auto quality = affine ? geometryQuality(map,constantJacobian) : geometryQuality(map) ;
      if(quality.status == GeometryQualityStatus::invalid) throw std::domain_error("hex split has invalid physical geometry") ;
      if(quality.status != GeometryQualityStatus::sampled) throw std::length_error("hex split physical quality is unresolved") ;
    }
  }

  std::uint32_t HexPlan::splitId(Cell cell) const {
    if(cell == root) return empty() ? absent : 0 ;
    if(!data_ || (cell-1)/8 >= data_->splits.size()) throw std::out_of_range("hex cell does not exist") ;
    return data_->splits[(cell-1)/8].children[(cell-1)%8] ;
  }
  HexPlan::Cell HexPlan::parent(Cell cell) const {
    require(cell != root,"root hex has no parent") ;
    splitId(cell) ; return data_->splits[(cell-1)/8].parent ;
  }
  unsigned HexPlan::depth(Cell cell) const {
    unsigned result = 0 ;
    while(cell != root) { cell = parent(cell) ; ++result ; }
    return result ;
  }
  bool HexPlan::active(Cell cell) const {
    splitId(cell) ;
    while(cell != root) {
      const auto &owner = data_->splits[(cell-1)/8] ;
      if(!owner.active) return false ;
      cell = owner.parent ;
    }
    return true ;
  }
  bool HexPlan::leaf(Cell cell) const {
    const auto split = splitId(cell) ;
    return active(cell) && (split == absent || !data_->splits[split].active) ;
  }
  std::array<HexPlan::Cell,8> HexPlan::children(Cell cell) const {
    const auto split = splitId(cell) ; require(split != absent,"hex has no accepted children") ;
    std::array<Cell,8> result ; for(unsigned i = 0; i < 8; ++i) result[i] = 1+8*split+i ; return result ;
  }
  HexCenters HexPlan::centers(Cell cell) const {
    const auto split = splitId(cell) ; require(split != absent,"hex has no accepted centers") ;
    const auto id = data_->splits[split].centers ; return id == absent ? HexCenters() : data_->centers[id] ;
  }
  std::vector<HexPlan::Cell> HexPlan::leaves() const {
    std::vector<Cell> result ;
    const auto visit = [&](const auto &self, Cell cell) -> void {
      const auto split = splitId(cell) ;
      if(split == absent || !data_->splits[split].active) result.push_back(cell) ;
      else for(unsigned i = 0; i < 8; ++i) self(self,1+8*split+i) ;
    } ;
    visit(visit,root) ; return result ;
  }
  std::vector<HexPlan::LeafBox> HexPlan::leafBoxes() const {
    std::vector<LeafBox> result ;
    const auto visit = [&](const auto &self, Cell cell, const LogicalBox &box, unsigned depth) -> void {
      const auto split = splitId(cell) ;
      if(split == absent || !data_->splits[split].active) result.push_back({cell,box.lower,depth}) ;
      else for(unsigned c = 0; c < 8; ++c) self(self,1+8*split+c,box.child(c),depth+1) ;
    } ;
    visit(visit,root,{},0) ; return result ;
  }
  HexPlan::Cell HexPlan::findLeaf(const std::array<std::uint32_t,3> &point) const {
    for(auto x : point) require(x < extent,"hex adjacency point is outside the root") ;
    Cell cell = root ; LogicalBox box ;
    for(;;) {
      const auto split = splitId(cell) ;
      if(split == absent || !data_->splits[split].active) return cell ;
      const auto c = box.childAt(point) ; cell = 1+8*split+c ; box = box.child(c) ;
    }
  }
  std::size_t HexPlan::leafCount() const {
    const auto visit = [&](const auto &self, Cell cell) -> std::size_t {
      const auto split = splitId(cell) ;
      if(split == absent || !data_->splits[split].active) return 1 ;
      std::size_t count = 0 ; for(unsigned i = 0; i < 8; ++i) count += self(self,1+8*split+i) ; return count ;
    } ;
    return visit(visit,root) ;
  }
  void HexPlan::refine(Cell cell, const HexCenters &centers, JacobianCheckLimits limits) {
    require(leaf(cell),"hex refinement requires an active leaf") ;
    auto id = splitId(cell) ;
    if(id != absent) { data_->splits[id].active = true ; return ; }
    if(depth(cell) >= maxDepth) throw std::length_error("hex refinement depth limit exceeded") ;
    if(splitCount() >= (absent-1)/8) throw std::length_error("hex plan exceeds its cell index range") ;
    certify(centers,limits) ;
    auto created = data_ ? nullptr : std::make_unique<Data>() ;
    auto &data = data_ ? *data_ : *created ;
    Split split{cell,{},absent,true} ; split.children.fill(absent) ;
    const bool custom = !(centers == HexCenters()) ;
    if(custom) { split.centers = data.centers.size() ; data.centers.push_back(centers) ; }
    id = data.splits.size() ;
    try { data.splits.push_back(split) ; }
    catch(...) { if(custom) data.centers.pop_back() ; throw ; }
    if(cell != root) data.splits[(cell-1)/8].children[(cell-1)%8] = id ;
    if(created) data_ = std::move(created) ;
  }

  void HexPlan::coarsen(Cell cell) {
    const auto id = splitId(cell) ;
    require(active(cell) && id != absent && data_->splits[id].active,"hex coarsening requires an active parent") ;
    for(auto child: children(cell)) require(leaf(child),"hex coarsening requires eight leaf children") ;
    data_->splits[id].active = false ;
  }
  HexPlan::Cell HexPlan::splitParent(std::size_t index) const {
    if(!data_) throw std::out_of_range("hex has no accepted splits") ;
    return data_->splits.at(index).parent ;
  }
  VertexVolumeSplit HexPlan::construction(Cell cell) const {
    const auto accepted = centers(cell) ; auto result = pattern() ;
    if(accepted == HexCenters()) return result ;
    for(unsigned f = 0; f < 6; ++f) result.points[20+f] = facePoint(f,accepted.faces[f]) ;
    result.points.back() = accepted.cell ;
    for(unsigned c = 0; c < 8; ++c) {
      const auto corners = childCorners(accepted,c) ;
      result.children[c].placement = std::make_shared<CornerMap>(Shape::hex,std::vector<Point>(corners.begin(),corners.end())) ;
    }
    return result ;
  }
  std::size_t HexPlan::storageBytes() const {
    return sizeof(*this)+(data_ ? sizeof(Data)+data_->splits.capacity()*sizeof(Split)+data_->centers.capacity()*sizeof(HexCenters) : 0) ;
  }

  std::vector<unsigned char> HexPlan::encode() const {
    if(empty()) return {} ;
    detail::RecordWriter out(std::numeric_limits<std::size_t>::max(),{}) ;
    out.u32(0x4858504c) ; out.u32(1) ; out.u64(data_->splits.size()) ; out.u64(data_->centers.size()) ;
    for(const auto &split: data_->splits) {
      out.u32(split.parent) ; for(auto child: split.children) out.u32(child) ;
      out.u32(split.centers) ; out.byte(split.active) ;
    }
    for(const auto &center: data_->centers) {
      for(const auto &face: center.faces) for(double x: face) out.real(x) ;
      for(double x: center.cell) out.real(x) ;
    }
    return out.take() ;
  }
  HexPlan HexPlan::decode(const std::vector<unsigned char> &bytes, std::size_t maxSplits, JacobianCheckLimits limits) {
    HexPlan result ; if(bytes.empty()) return result ;
    detail::RecordReader in(bytes,bytes.size(),{}) ;
    require(in.u32() == 0x4858504c && in.u32() == 1,"unsupported hex plan format") ;
    const auto count = in.count(std::min(maxSplits,std::size_t((absent-1)/8)),41) ;
    const auto centers = in.count(count,120) ;
    require(count > 0,"noncanonical empty hex plan") ;
    require(in.remaining() == 41*count+120*centers,"invalid hex plan length") ;
    result.data_ = std::make_unique<Data>() ;
    result.data_->splits.resize(count) ; result.data_->centers.resize(centers) ;
    std::size_t nextCenter = 0 ;
    for(std::size_t i = 0; i < count; ++i) {
      auto &split = result.data_->splits[i] ; split.parent = in.u32() ;
      for(auto &child: split.children) {
        child = in.u32() ; require(child == absent || (child > i && child < count),"invalid hex child link") ;
      }
      split.centers = in.u32() ; const auto active = in.byte() ;
      require(active <= 1,"invalid hex selection") ; split.active = active ;
      if(split.centers != absent) require(split.centers == nextCenter++ && split.centers < centers,"invalid hex center index") ;
      if(i == 0) require(split.parent == root,"invalid hex root") ;
      else {
        require(split.parent > 0 && (split.parent-1)/8 < i,"invalid hex parent link") ;
        const auto &owner = result.data_->splits[(split.parent-1)/8] ;
        require(owner.children[(split.parent-1)%8] == i,"hex parent and child links disagree") ;
        require(!split.active || owner.active,"active hex below a coarsened parent") ;
      }
      if(result.depth(split.parent) >= maxDepth) throw std::length_error("hex refinement depth limit exceeded") ;
    }
    require(nextCenter == centers,"unused hex center records") ;
    for(std::size_t i = 0; i < count; ++i) for(unsigned c = 0; c < 8; ++c) {
      const auto child = result.data_->splits[i].children[c] ;
      if(child != absent) require(result.data_->splits[child].parent == 1+8*i+c,"hex child has multiple parents") ;
    }
    for(auto &center: result.data_->centers) {
      for(auto &face: center.faces) for(double &x: face) x = in.real() ;
      for(double &x: center.cell) x = in.real() ;
      require(!(center == HexCenters()),"redundant midpoint center record") ; certify(center,limits) ;
    }
    in.finish() ; return result ;
  }

  HexCellMap::HexCellMap(const GeometryMap &source, const HexPlan &plan, HexPlan::Cell cell)
    : source_(source), plan_(plan), cell_(cell) {
    require(source.shape() == Shape::hex,"hex plan requires hex source geometry") ; plan.splitId(cell) ;
  }
  bool HexCellMap::affine() const {
    if(!affineCorners(source_)) return false ;
    auto cell = cell_ ;
    while(cell != HexPlan::root) {
      const auto &split = plan_.data_->splits[(cell-1)/8] ;
      if(split.centers != HexPlan::absent) return false ;
      cell = split.parent ;
    }
    return true ;
  }
  MapValue HexCellMap::evaluate(const Point &point) const {
    MapValue local{point,{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}},true} ;
    auto cell = cell_ ;
    while(cell != HexPlan::root) {
      const auto &split = plan_.data_->splits[(cell-1)/8] ; const unsigned ordinal = (cell-1)%8 ;
      const auto placement = split.centers == HexPlan::absent ? midpointValue(ordinal,local.position)
        : evaluateHexCorners(childCorners(plan_.data_->centers[split.centers],ordinal),local.position) ;
      local = compose(placement,local) ; cell = split.parent ;
    }
    return compose(source_.evaluate(local.position),local) ;
  }
  HexCenters hexCenters(const GeometryMap &geometry, CenterMethod method, MeasureOptions options) {
    require(geometry.shape() == Shape::hex,"hex centers require hex geometry") ; HexCenters result ;
    require(method == CenterMethod::reference || method == CenterMethod::weighted || method == CenterMethod::integrated,"unknown hex center method") ;
    if(method == CenterMethod::reference || affineHex(geometry)) return result ;
    for(unsigned f = 0; f < 6; ++f) {
      const auto p = faceCenterCoordinates(HexFaceMap(geometry,f),method,options) ; result.faces[f] = {{p[0],p[1]}} ;
    }
    result.cell = cellCenterCoordinates(geometry,facePoints(result),method,options) ; return result ;
  }
  void refineHex(HexPlan &plan, const GeometryMap &source, HexPlan::Cell cell, const HexCenters &centers, JacobianCheckLimits limits) {
    require(plan.leaf(cell),"hex refinement requires an active leaf") ;
    require(source.shape() == Shape::hex,"hex plan requires hex source geometry") ;
    if(plan.accepted(cell)) { plan.refine(cell) ; return ; }
    const auto *corners = dynamic_cast<const CornerMap*>(&source) ;
    // Positive source and reference-map Jacobians stay positive under composition.
    if(!corners || !positiveHexJacobian(corners->corners())) {
      const HexCellMap parent(source,plan,cell) ; requireQuality(parent,parent.affine()) ;
      for(unsigned child = 0; child < 8; ++child) requireQuality(HexChildMap(parent,centers,child)) ;
    }
    plan.refine(cell,centers,limits) ;
  }
  void refineHex(HexPlan &plan, const GeometryMap &source, HexPlan::Cell cell, CenterMethod method, MeasureOptions options,
                 JacobianCheckLimits limits) {
    const HexCellMap parent(source,plan,cell) ;
    const auto centers = plan.accepted(cell) ? plan.centers(cell) : hexCenters(parent,method,options) ;
    refineHex(plan,source,cell,centers,limits) ;
  }
}
}
