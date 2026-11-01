// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include "parametric_mixed.h"
#include "parametric_cell_data.h"
#include "parametric_hex.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace Loci {
namespace FVMAdapt {
  std::vector<ParametricNodeKey> decodeParametricNodeKeys(const Loci::FineFaces &record) {
    std::vector<ParametricNodeKey> keys ; keys.reserve(record.size()) ;
    for(const auto &row : record) {
      if(row.size() != 9 || row[0] < 0 || row[0] > 3 || row[1] < 0)
        throw std::invalid_argument("invalid parametric node record") ;
      keys.push_back({{EntityKind(row[0]),EntityId(row[1]),
        {{Fraction(row[2],row[3]),Fraction(row[4],row[5]),Fraction(row[6],row[7])} } },row[8]}) ;
    }
    return keys ;
  }

  Loci::FineFaces encodeParametricNodeKeys(const std::vector<ParametricNodeKey> &keys) {
    Loci::FineFaces record ; record.reserve(keys.size()) ;
    for(const auto &key : keys) {
      std::vector<int> row = {int(key.reference.kind),int(key.reference.entity)} ;
      for(const auto &coordinate : key.reference.point) {
        row.push_back(checkedExactInteger<int>(coordinate.numerator())) ;
        row.push_back(checkedExactInteger<int>(coordinate.denominator())) ;
      }
      row.push_back(key.local) ; record.push_back(std::move(row)) ;
    }
    return record ;
  }

  namespace {
    void includeEdge(std::array<double,2> &lengths, const Point &a, const Point &b) {
      const double xy=std::hypot(b[0]-a[0],b[1]-a[1]), xyz=std::hypot(xy,b[2]-a[2]) ;
      if(!std::isfinite(xyz) || xyz<=0) throw std::invalid_argument("invalid leaf edge length") ;
      lengths[0]=std::max(lengths[0],xyz) ; lengths[1]=std::max(lengths[1],xy) ;
    }
    std::array<double,2> mapEdgeLengths(const GeometryMap &map) {
      const auto &element=referenceElement(map.shape()) ;
      std::array<Point,8> points ; std::array<double,2> lengths{} ;
      for(unsigned i=0;i<element.vertexCount;++i) points[i]=map.evaluate(element.vertices[i]).position ;
      for(unsigned i=0;i<element.edgeCount;++i)
        includeEdge(lengths,points[element.edges[i][0]],points[element.edges[i][1]]) ;
      return lengths ;
    }
  }
  struct ParametricCellData::Data {
    ParametricRecord record ;
    std::optional<ParametricCell> cell ;
    std::unique_ptr<ParametricHex> hex ;
    std::string error ;
    int maxDepth=10,maxNodes=1000000 ;
    void retain(ParametricCell value) {
      const auto *corner = value.corner() ;
      if(corner && corner->topology.shape() == Shape::hex) {
        value.validate(-1) ;
        hex = std::make_unique<ParametricHex>(*corner) ; ParametricRecord().swap(record) ;
      } else {
        if(!value.source()) record = encodeParametricCell(value) ;
        else ParametricRecord().swap(record) ;
        cell = std::move(value) ;
      }
    }
  } ;

  ParametricCellData::ParametricCellData(ParametricRecord record,int maxDepth,int maxNodes) {
    auto data=std::make_shared<Data>() ; data->maxDepth=maxDepth ; data->maxNodes=maxNodes ;
    try {
      const auto limits=parametricHistoryLimits(maxDepth,maxNodes) ;
      if(ParametricHex::saved(record)) data->hex=std::make_unique<ParametricHex>(ParametricHex::decode(record,limits)) ;
      else data->retain(decodeParametricCell(record,limits)) ;
    } catch(const std::exception &error) { data->error=error.what() ; }
    data_=std::move(data) ;
  }
  ParametricCellData::ParametricCellData(ParametricCell cell) {
    auto data=std::make_shared<Data>() ; const auto &limits=cell.limits() ;
    if(limits.maxDepth>30 || limits.maxNodes>std::size_t(std::numeric_limits<int>::max()))
      throw std::length_error("parametric cell limits exceed module range") ;
    data->maxDepth=limits.maxDepth ; data->maxNodes=limits.maxNodes ; data->retain(std::move(cell)) ; data_=std::move(data) ;
  }
  ParametricCellData ParametricCellData::initial(const CellTopology &topology,const std::vector<Point> &corners,
      const ParametricRecord &record,int maxDepth,int maxNodes,const std::string &requested,int axis) {
    const auto family=parametricCellFamily(topology,requested) ;
    if(requested=="mixed" && axis>=0) throw std::invalid_argument("mixed splitting requires full 3D refinement") ;
    const auto limits=parametricHistoryLimits(maxDepth,maxNodes) ;
    if(ParametricHex::saved(record)) {
      ParametricCellData result(record,maxDepth,maxNodes) ;
      if(!result.data_->hex) throw std::invalid_argument(result.data_->error) ;
      const auto &hex=*result.data_->hex ;
      if(family!="corner" || axis>=0 || hex.topology!=topology || hex.corners!=corners)
        throw std::invalid_argument("saved compact hex source or split family differs") ;
      return result ;
    }
    if(record.empty() && family=="corner") {
      if(axis>=0) throw std::invalid_argument("corner splitting requires full 3D refinement") ;
      return ParametricCellData(ParametricCell::original(topology,corners,limits)) ;
    }
    return ParametricCellData(makeParametricCell(topology,corners,record,limits,family,axis)) ;
  }
  std::shared_ptr<const ParametricCell> ParametricCellData::workingCell() const {
    if(!data_) throw std::invalid_argument("parametric cell is not initialized") ;
    if(data_->hex) return std::make_shared<const ParametricCell>(data_->hex->expand()) ;
    if(!data_->cell) throw std::invalid_argument(data_->error) ;
    return std::shared_ptr<const ParametricCell>(data_,&*data_->cell) ;
  }
  ParametricCell ParametricCellData::workingCopy() const {
    if(data_ && data_->hex) return data_->hex->expand() ;
    return *workingCell() ;
  }
  std::vector<ParametricLeafRequest> ParametricCellData::decodeRequests(const ParametricRecord &record) const {
    if(const auto plan=hexPlan()) {
      const auto root=compactHex() ? data_->hex->topology.id() : data_->cell->source()->topology.id() ;
      const auto &limits=compactHex() ? data_->hex->limits : data_->cell->limits() ;
      return hexRequests(record,*plan,root,limits) ;
    }
    return decodeParametricRequests(record,*workingCell()) ;
  }
  std::optional<ParametricCellData> ParametricCellData::editPlan(const std::vector<ParametricLeafRequest> &requests,int levels,CenterMethod method,
      const std::vector<std::pair<EntityId,const ParametricRecord *>> &incoming,bool *changed,const GeometricRefinement *selection) const {
    if(changed) *changed=false ;
    if(!hexPlan()) { workingCell() ; return std::nullopt ; }
    if(levels<1 || levels>30) throw std::invalid_argument("parametricRefinementLevels must be in [1,30]") ;
    if(method!=CenterMethod::reference && method!=CenterMethod::weighted && method!=CenterMethod::integrated)
      throw std::invalid_argument("unknown parametric center method") ;
    for(const auto &request:requests)
      if(request.action!=ParametricAction::isotropic && request.action!=ParametricAction::coarsen) return std::nullopt ;
    if(selection && selection->select && (!std::isfinite(selection->minEdgeLength) || selection->minEdgeLength<0))
      throw std::invalid_argument("geometric refinement requires a nonnegative finite minEdgeLength") ;
    if(requests.empty() && incoming.empty() && !(selection && selection->select)) return *this ;
    auto data=std::make_shared<Data>() ; data->maxDepth=data_->maxDepth ; data->maxNodes=data_->maxNodes ;
    if(compactHex()) data->hex=std::make_unique<ParametricHex>(*data_->hex) ;
    else {
      const auto &source=*data_->cell->source() ;
      data->hex=std::make_unique<ParametricHex>(source.topology,source.corners,source.limits) ;
    }
    std::array<const FacePlan *,6> faces{} ;
    for(const auto &input:incoming) {
      unsigned f=0 ; while(f<6 && data->hex->topology.faceId(f)!=input.first) ++f ;
      if(f==6 || !input.second || faces[f]) throw std::invalid_argument("invalid or repeated incoming hex face") ;
      faces[f]=input.second ;
    }
    if(!data->hex->edit(requests,levels,method,faces,selection)) return *this ;
    ParametricCellData result ; result.data_=std::move(data) ;
    if(changed) *changed=true ;
    return result ;
  }
  ParametricCellData ParametricCellData::applyRequests(const ParametricRecord &record,const CellRefinementOptions &options,
      const IncomingFacePlans &incoming,const GeometricRefinement &selection) const {
    if(!compactHex() && workingCell()->native() && options.center!=CenterMethod::reference)
      throw std::invalid_argument("native cuts support only reference centers") ;
    if(options.levels<1 || options.levels>30) throw std::invalid_argument("parametricRefinementLevels must be in [1,30]") ;
    if(record.empty() && incoming.empty() && !selection.select) return *this ;
    const bool compact=hexPlan()!=nullptr ;
    auto requests=compact ? decodeRequests(record) : std::vector<ParametricLeafRequest>{} ;
    if(compact && options.extrusionAxis<0) {
      const auto edited=editPlan(requests,options.levels,options.center,incoming,nullptr,&selection) ;
      if(edited) return *edited ;
    }
    auto cell=workingCopy() ;
    for(const auto &face:incoming) cell.mergeFace(face.first,*face.second) ;
    if(!record.empty()) cell.expand() ;
    if(!compact) requests=decodeParametricRequests(record,cell) ;
    for(const auto &request:requests) {
      if(request.action==ParametricAction::coarsen) cell.coarsen(request.cell) ;
      else if(request.action==ParametricAction::isotropic) {
        std::vector<TreeNodeId> leaves{request.cell} ;
        for(int level=0;level<options.levels;++level) {
          std::vector<TreeNodeId> next ;
          for(auto leaf:leaves) {
            if(options.automatic && cell.native()) {
              const auto &history=cell.history() ;
              const auto split=automaticSplit(*history.geometry().map(history.node(leaf).map),options.aspectRatio,options.minEdgeLength) ;
              if(!split) continue ;
              cell.refine(leaf,*split,options.center) ;
            } else cell.refineIsotropic(leaf,options.center,options.extrusionAxis) ;
            if(level+1<options.levels) for(auto child:cell.children(leaf)) next.push_back(child) ;
          }
          leaves.swap(next) ;
        }
      } else cell.refine(request.cell,request.split,options.center) ;
    }
    selection.refine(cell,options.center,options.extrusionAxis,options.automatic,options.aspectRatio) ;
    return !requests.empty() || !cell.source() ? ParametricCellData(std::move(cell)) : *this ;
  }

  ParametricCellData ParametricCellData::balance(const std::set<int> &required,const std::map<int,unsigned> &directions,
      const ParametricCellData &previous,CenterMethod method,int axis,const IncomingFacePlans &incoming,bool &changed) const {
    changed=false ;
    if(required.empty() && incoming.empty()) return *this ;
    if(const auto plan=hexPlan()) {
      const auto leaves=plan->leaves() ; std::vector<ParametricLeafRequest> requests ;
      for(int index:required) requests.push_back({leaves.at(index),TensorSplit(),ParametricAction::isotropic}) ;
      const auto edited=editPlan(requests,1,method,incoming,&changed) ;
      if(edited) return *edited ;
    }
    auto cell=workingCopy() ; bool faceChanged=false ;
    for(const auto &face:incoming) faceChanged=cell.mergeFace(face.first,*face.second) || faceChanged ;
    if(!required.empty()) {
      cell.expand() ; const auto leaves=cell.leafIds() ; const auto old=previous.workingCell() ;
      for(int index:required) if(old->source() || !cell.restoreSplit(leaves.at(index),*old)) {
        const auto selected=directions.find(index) ;
        cell.refineForBalance(leaves.at(index),method,axis,selected==directions.end() ? 0:selected->second) ;
      }
    }
    if(required.empty() && !faceChanged) return *this ;
    ParametricCellData result(std::move(cell)) ; changed=true ; return result ;
  }

  ParametricRecord ParametricCellData::tagRequests(const std::vector<char> &tags) const {
    if(const auto plan=hexPlan()) {
      const auto root=compactHex() ? data_->hex->topology.id() : data_->cell->source()->topology.id() ;
      const auto &limits=compactHex() ? data_->hex->limits : data_->cell->limits() ;
      return hexTagRequests(*plan,root,limits,tags) ;
    }
    return encodeParametricTags(*workingCell(),tags) ;
  }
  std::vector<ParametricNodeKey> ParametricCellData::cellNodes(const CellTopology &topology) const {
    return compactHex() ? data_->hex->cellNodes() : parametricCellNodes(topology,*workingCell()) ;
  }
  std::vector<Point> ParametricCellData::cellPositions(const std::vector<ParametricNodeKey> &keys) const {
    return compactHex() ? data_->hex->cellPositions(keys) : parametricCellPositions(*workingCell(),keys) ;
  }
  std::vector<ParametricOutputFace> ParametricCellData::cellFaces(const CellTopology &topology,const std::vector<ParametricNodeKey> &required) const {
    return compactHex() ? data_->hex->cellFaces(required) : parametricCellFaces(topology,*workingCell(),required) ;
  }
  std::vector<ParametricNodeKey> ParametricCellData::faceNodes(EntityId face,const CellTopology &topology,
      const ParametricCellData *right,const CellTopology *rightTopology) const {
    if(compactHex() && (!right || right->compactHex())) return data_->hex->faceNodes(face,right ? right->data_->hex.get():nullptr) ;
    const auto leftCell = workingCell(),rightCell = right ? right->workingCell():nullptr ;
    return parametricFaceNodes(face,topology,*leftCell,rightTopology,rightCell.get()) ;
  }
  std::vector<Point> ParametricCellData::facePositions(EntityId face,const CellTopology &topology,
      const ParametricCellData *right,const CellTopology *rightTopology) const {
    if(compactHex() && (!right || right->compactHex())) return data_->hex->facePositions(face,right ? right->data_->hex.get():nullptr) ;
    const auto leftCell = workingCell(),rightCell = right ? right->workingCell():nullptr ;
    return parametricFacePositions(face,topology,*leftCell,rightTopology,rightCell.get()) ;
  }
  std::vector<ParametricOutputFace> ParametricCellData::faceFaces(EntityId face,const CellTopology &topology,
      const std::vector<ParametricNodeKey> &required,const ParametricCellData *right,const CellTopology *rightTopology) const {
    if(compactHex() && (!right || right->compactHex())) return data_->hex->faceFaces(face,required,right ? right->data_->hex.get():nullptr) ;
    const auto leftCell = workingCell(),rightCell = right ? right->workingCell():nullptr ;
    return parametricFaceFaces(face,topology,*leftCell,required,rightTopology,rightCell.get()) ;
  }
  const HexPlan *ParametricCellData::hexPlan() const {
    if(data_ && data_->hex) return &data_->hex->plan ;
    if(data_ && data_->cell && data_->cell->source() && data_->cell->source()->topology.isStandard() && data_->cell->source()->topology.shape() == Shape::hex) {
      static const HexPlan empty ; return &empty ;
    }
    return nullptr ;
  }
  std::vector<int> ParametricCellData::localRequests(const CellTopology &topology, unsigned difference, NeighborGrading grading, int axis) const {
    if(compactHex()) {
      if(axis >= 0) throw std::invalid_argument("corner splitting requires full 3D refinement") ;
      return data_->hex->localRequests(difference,grading) ;
    }
    const auto cell = workingCell() ; cell->validate(axis) ; return parametricLocalRequests(topology,*cell,difference,grading) ;
  }
  ParametricTrace ParametricCellData::edgeTrace(const CellTopology &topology) const {
    return compactHex() ? data_->hex->edgeTrace() : parametricEdgeTrace(topology,*workingCell()) ;
  }
  std::vector<int> ParametricCellData::vertexDepths(const CellTopology &topology) const {
    return compactHex() ? data_->hex->vertexDepths() : parametricVertexDepths(topology,*workingCell()) ;
  }
  std::vector<int> ParametricCellData::vertexRequests(const CellTopology &topology, const std::map<EntityId,int> &depth, unsigned difference) const {
    return compactHex() ? data_->hex->vertexRequests(depth,difference) : parametricVertexRequests(topology,*workingCell(),depth,difference) ;
  }
  ParametricTrace ParametricCellData::faceRequests(EntityId face, const CellTopology &left,
      const ParametricCellData &rightData, const CellTopology &right, unsigned difference, NeighborGrading grading) const {
    const auto a = hexPlan(), b = rightData.hexPlan() ;
    if(a && b) {
      mergeFacePlans(facePlan(face),rightData.facePlan(face)) ;
      return hexFaceRequests(face,left,*a,right,*b,difference,grading) ;
    }
    const auto l = workingCell(), r = rightData.workingCell() ;
    return parametricFaceRequests(face,left,*l,right,*r,difference,grading) ;
  }
  std::vector<FaceBalanceRequest> ParametricCellData::dividedFaces(EntityId face,const CellTopology &left,
      const ParametricCellData &rightData,const CellTopology &right,bool directions) const {
    const auto a=hexPlan(), b=rightData.hexPlan() ;
    if(a && b && !directions) return hexDividedFaces(face,left,*a,right,*b) ;
    return parametricDividedFaces(face,left,*workingCell(),right,*rightData.workingCell(),directions) ;
  }
  std::vector<int> ParametricCellData::faceBalanceRequests(const CellTopology &topology,
      const std::vector<FaceBalanceRequest> &shared,FaceBalance policy,std::map<int,unsigned> *directions) const {
    if(policy==FaceBalance::none && !directions) return {} ;
    const auto plan=hexPlan() ;
    if(!plan || directions) return parametricFaceBalanceRequests(topology,*workingCell(),shared,policy,directions) ;
    std::map<int,std::vector<unsigned>> divided ;
    for(const auto &use:hexDividedFaces(topology,*plan)) divided[use.leaf].push_back(use.face) ;
    for(const auto &use:shared) if(use.root==int(topology.id())) {
      if(use.face<0) throw std::invalid_argument("face balancing received a level-grading request") ;
      divided[use.leaf].push_back(use.face) ;
    }
    const auto &element=referenceElement(Shape::hex) ; std::vector<std::vector<unsigned>> faces ;
    for(unsigned f=0;f<6;++f) {
      const auto &b=element.boundaries[f] ; faces.emplace_back(b.vertices.begin(),b.vertices.begin()+b.vertexCount) ;
    }
    const auto count=plan->leafCount() ; std::vector<int> out ;
    for(const auto &use:divided) {
      if(use.first<0 || std::size_t(use.first)>=count) throw std::invalid_argument("divided face has invalid leaf index") ;
      if(faceBalanceRequired(faces,use.second,policy)) out.push_back(use.first) ;
    }
    return out ;
  }
  std::vector<std::vector<int>> ParametricCellData::donorsFrom(const ParametricCellData &previous) const {
    const auto old = previous.hexPlan(), next = hexPlan() ;
    if(old && next) {
      const auto topology = [](const Data &d) -> const CellTopology & { return d.hex ? d.hex->topology : d.cell->source()->topology ; } ;
      const auto corners = [](const Data &d) -> const std::vector<Point> & { return d.hex ? d.hex->corners : d.cell->source()->corners ; } ;
      if(topology(*data_)!=topology(*previous.data_) || corners(*data_)!=corners(*previous.data_))
        throw std::invalid_argument("cell correspondence requires the same original geometry") ;
      if(!data_->hex && previous.data_->hex) throw std::invalid_argument("accepted refinement history cannot be replaced by an untouched source") ;
      return hexCellDonors(*old,*next) ;
    }
    const auto oldCell = previous.workingCell(), newCell = workingCell() ; return parametricCellDonors(*oldCell,*newCell) ;
  }
  bool ParametricCellData::compactHex() const { return data_ && bool(data_->hex) ; }
  std::size_t ParametricCellData::leafCount() const { return compactHex() ? data_->hex->plan.leafCount() : workingCell()->leafCount() ; }
  std::vector<int> ParametricCellData::leafDepths() const {
    std::vector<int> depths ;
    if(compactHex()) {
      for(auto leaf:data_->hex->plan.leaves()) depths.push_back(data_->hex->plan.depth(leaf)) ;
    } else {
      const auto cell=workingCell() ;
      if(cell->source()) return {0} ;
      for(auto leaf:cell->leafIds()) depths.push_back(cell->path(leaf).size()) ;
    }
    return depths ;
  }
  ParametricRecord ParametricCellData::facePlan(EntityId face) const {
    if(compactHex()) {
      for(unsigned f=0;f<6;++f) if(data_->hex->topology.faceId(f)==face) return data_->hex->faces[f] ;
      throw std::invalid_argument("face is not incident to the compact hex") ;
    }
    return parametricFacePlan(face,*workingCell()) ;
  }
  bool ParametricCellData::needsFaceMerge(EntityId face, const ParametricRecord &record) const {
    if(record.empty()) return false ;
    if(compactHex()) {
      for(unsigned f=0;f<6;++f) if(data_->hex->topology.faceId(f)==face) return data_->hex->faces[f]!=record ;
      throw std::invalid_argument("face is not incident to the compact hex") ;
    }
    if(workingCell()->native()) { validateMixedFacePlan(record) ; return false ; }
    return facePlan(face)!=record ;
  }
  std::vector<std::array<double,2>> ParametricCellData::edgeLengths() const {
    std::vector<std::array<double,2>> lengths ;
    if(compactHex()) {
      const auto &hex=*data_->hex ; CornerMap source(Shape::hex,hex.corners) ;
      lengths.reserve(hex.plan.leafCount()) ;
      for(auto leaf:hex.plan.leaves()) lengths.push_back(mapEdgeLengths(HexCellMap(source,hex.plan,leaf))) ;
      return lengths ;
    }
    const auto cell=workingCell() ;
    if(const auto *source=cell->source()) {
      std::array<double,2> value{} ;
      for(unsigned edge=0;edge<source->topology.edgeCount();++edge) {
        const auto vertices=source->topology.edgeVertices(edge) ; includeEdge(value,source->corners[vertices[0]],source->corners[vertices[1]]) ;
      }
      return {value} ;
    }
    lengths.reserve(cell->leafCount()) ;
    if(const auto *extruded=cell->extruded()) {
      const int axis=extruded->axis ;
      double low=extruded->corners.front()[axis],high=low ;
      for(const auto &point:extruded->corners) { low=std::min(low,point[axis]) ; high=std::max(high,point[axis]) ; }
      Point height{} ; height[axis]=high-low ;
      for(auto leaf:extruded->planar.leafIds()) {
        std::array<double,2> value{} ; includeEdge(value,Point{},height) ;
        const auto &corners=extruded->planar.node(leaf).corners ;
        const auto position=[&](const PlanarVertexKey &key) {
          const auto p=extruded->planar.evaluate(key) ; Point q{} ;
          q[(axis+1)%3]=p[0] ; q[(axis+2)%3]=p[1] ; return q ;
        } ;
        for(std::size_t i=0;i<corners.size();++i)
          includeEdge(value,position(corners[i]),position(corners[(i+1)%corners.size()])) ;
        lengths.push_back(value) ;
      }
    } else if(cell->boundary()) {
      const auto &state=cell->cornerGeometry() ;
      for(auto leaf:cell->leafIds()) {
        std::array<double,2> value{} ;
        for(const auto &edge:state.cellEdges(cell->cornerCell(leaf)))
          includeEdge(value,state.boundaries.points[edge.first.first],state.boundaries.points[edge.first.second]) ;
        lengths.push_back(value) ;
      }
    } else {
      const auto &history=cell->history() ;
      for(auto leaf:cell->leafIds())
        lengths.push_back(mapEdgeLengths(*history.geometry().map(history.node(leaf).map))) ;
    }
    return lengths ;
  }
  ParametricRecord ParametricCellData::record() const {
    return compactHex() ? data_->hex->record() : (data_ ? data_->record : ParametricRecord{}) ;
  }
  int ParametricCellData::maxDepth() const { return data_ ? data_->maxDepth:10 ; }
  int ParametricCellData::maxNodes() const { return data_ ? data_->maxNodes:1000000 ; }

  void ParametricCellDataConverter::encode() {
    record_=value_.record() ;
    if(value_.data_ && value_.data_->cell && value_.data_->cell->source()) record_=encodeParametricCell(*value_.data_->cell) ;
  }
  ParametricCellDataConverter::ParametricCellDataConverter(ParametricCellData &value) : value_(value) { encode() ; }
  int ParametricCellDataConverter::getSize() const {
    if(record_.size()>64*1024*1024) throw std::length_error("parametric cell record exceeds module limit") ;
    return 8+record_.size() ;
  }
  void ParametricCellDataConverter::getState(unsigned char *buffer,int &size) const {
    size=getSize() ; const unsigned limits[2]={unsigned(value_.maxDepth()),unsigned(value_.maxNodes())} ;
    for(unsigned i=0;i<2;++i) for(unsigned b=0;b<4;++b) buffer[4*i+b]=(limits[i]>>(8*b))&255 ;
    std::copy(record_.begin(),record_.end(),buffer+8) ;
  }
  void ParametricCellDataConverter::setState(const unsigned char *buffer,int size) {
    record_.clear() ;
    if(size<8 || size>8+64*1024*1024) { value_=ParametricCellData({},0,1) ; return ; }
    unsigned limits[2]={0,0} ;
    for(unsigned i=0;i<2;++i) for(unsigned b=0;b<4;++b) limits[i]|=unsigned(buffer[4*i+b])<<(8*b) ;
    const int depth=limits[0]<=30 ? int(limits[0]):-1 ;
    const int nodes=limits[1]<=unsigned(std::numeric_limits<int>::max()) ? int(limits[1]):0 ;
    value_=ParametricCellData(ParametricRecord(buffer+8,buffer+size),depth,nodes) ; encode() ;
  }
  std::ostream &operator<<(std::ostream &out,const ParametricCellData &value) {
    return out << "parametric cell (" << value.record().size() << " saved bytes)" ;
  }
  std::istream &operator>>(std::istream &in,ParametricCellData &) { in.setstate(std::ios::failbit) ; return in ; }
}
}
