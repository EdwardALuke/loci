// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include "parametric_mixed.h"
#include "parametric_cell_data.h"
#include "compact_hex.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
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
      std::vector<int> row ; row.reserve(9) ;
      row.push_back(int(key.reference.kind)) ; row.push_back(int(key.reference.entity)) ;
      for(const auto &coordinate : key.reference.point) {
        row.push_back(checkedExactInteger<int>(coordinate.numerator())) ;
        row.push_back(checkedExactInteger<int>(coordinate.denominator())) ;
      }
      row.push_back(key.local) ; record.push_back(std::move(row)) ;
    }
    return record ;
  }

  namespace {
    struct NodeRecordLess {
      template<class A,class B> bool operator()(const A &a,const B &b) const {
        for(unsigned i=0;i<2;++i) if(a[i]!=b[i]) return a[i]<b[i] ;
        // Products of signed 32-bit numerators and positive denominators fit in 64 bits.
        static_assert(std::numeric_limits<int>::digits<=31,"node records require at most 32-bit integers") ;
        for(unsigned i=2;i<8;i+=2) {
          const auto left=std::int64_t(a[i])*b[i+1],right=std::int64_t(b[i])*a[i+1] ;
          if(left!=right) return left<right ;
        }
        return a[8]<b[8] ;
      }
    } ;
    bool nodeOwnerLess(const std::vector<int> &row,EntityKind kind,EntityId entity) {
      return row[0]!=int(kind) ? row[0]<int(kind) : EntityId(row[1])<entity ;
    }
  }

  void ParametricNodeNumbers::add(const Loci::FineFaces &record,int offset) {
    if(record.empty()) return ;
    if(offset<0 || record.size()-1>std::size_t(std::numeric_limits<int>::max()-offset))
      throw std::overflow_error("parametric node numbers exceed output range") ;
    const NodeRecordLess less ; bool local=true ;
    for(std::size_t i=0;i<record.size();++i) {
      const auto &row=record[i] ;
      if(row.size()!=9 || row[0]<0 || row[0]>3 || row[1]<0 || row[3]<=0 || row[5]<=0 || row[7]<=0)
        throw std::invalid_argument("invalid parametric node record") ;
      if(row[0]!=record[0][0] || row[1]!=record[0][1] || (i && !less(record[i-1],row)))
        throw std::invalid_argument("parametric owner nodes are not strictly ordered") ;
      local=local && row[8]>=0 && row[2]==0 && row[4]==0 && row[6]==0 ;
    }
    const auto kind=EntityKind(record[0][0]) ; const auto entity=EntityId(record[0][1]) ;
    const auto owner=std::lower_bound(owners_.begin(),owners_.end(),entity,[&](const Owner &item,EntityId id) {
      return nodeOwnerLess(item.record->front(),kind,id) ;
    }) ;
    if(owner!=owners_.end() && owner->record->front()[0]==int(kind) && EntityId(owner->record->front()[1])==entity) {
      if(owner->offset!=offset || owner->record->size()!=record.size())
        throw std::logic_error("parametric node owner has conflicting output numbers") ;
      if(owner->record!=&record) for(std::size_t i=0;i<record.size();++i)
        if(less((*owner->record)[i],record[i]) || less(record[i],(*owner->record)[i]))
          throw std::logic_error("parametric node owner has conflicting output numbers") ;
      return ;
    }
    owners_.insert(owner,{&record,offset,local}) ;
  }

  int ParametricNodeNumbers::at(const ParametricNodeKey &key) const {
    const auto kind=key.reference.kind ; const auto entity=key.reference.entity ;
    const auto owner=std::lower_bound(owners_.begin(),owners_.end(),entity,[&](const Owner &item,EntityId id) {
      return nodeOwnerLess(item.record->front(),kind,id) ;
    }) ;
    if(owner==owners_.end() || owner->record->front()[0]!=int(kind) || EntityId(owner->record->front()[1])!=entity)
      throw std::out_of_range("parametric output node owner is absent") ;
    const auto &record=*owner->record ; auto node=record.end() ;
    if(owner->local && key.local>=0 && key.reference.point[0].rational()==0 &&
        key.reference.point[1].rational()==0 && key.reference.point[2].rational()==0) {
      node=std::lower_bound(record.begin(),record.end(),key.local,[](const auto &row,int local) { return row[8]<local ; }) ;
      if(node==record.end() || (*node)[8]!=key.local) throw std::out_of_range("parametric output node is absent") ;
    } else {
      std::array<int,9> row{{int(kind),int(entity),0,1,0,1,0,1,key.local}} ;
      for(unsigned d=0;d<3;++d) {
        row[2+2*d]=checkedExactInteger<int>(key.reference.point[d].numerator()) ;
        row[3+2*d]=checkedExactInteger<int>(key.reference.point[d].denominator()) ;
      }
      const NodeRecordLess less ; node=std::lower_bound(record.begin(),record.end(),row,less) ;
      if(node==record.end() || less(row,*node)) throw std::out_of_range("parametric output node is absent") ;
    }
    // Persistent local IDs may have gaps; output numbers follow the owned row order.
    return owner->offset+int(node-record.begin()) ;
  }

  namespace {
    void checkIncomingFaces(const IncomingFacePlans &incoming) {
      if(incoming.size()<2) return ;
      std::vector<EntityId> faces ; faces.reserve(incoming.size()) ;
      for(const auto &input:incoming) faces.push_back(input.face) ;
      std::sort(faces.begin(),faces.end()) ;
      if(std::adjacent_find(faces.begin(),faces.end())!=faces.end())
        throw std::invalid_argument("repeated incoming face") ;
    }
    void includeEdge(std::array<double,2> &lengths, const Point &a, const Point &b) {
      const double xy=std::hypot(b[0]-a[0],b[1]-a[1]), xyz=std::hypot(xy,b[2]-a[2]) ;
      if(!std::isfinite(xyz) || xyz<=0) throw std::invalid_argument("invalid leaf edge length") ;
      lengths[0]=std::max(lengths[0],xyz) ; lengths[1]=std::max(lengths[1],xy) ;
    }
    std::array<double,2> mapEdgeLengths(const GeometryMap &map) {
      const auto &element=referenceElement(map.shape()) ;
      std::array<Point,8> points ; std::array<double,2> lengths{} ;
      for(unsigned i=0;i<element.vertexCount;++i) points[i]=map.position(element.vertices[i]) ;
      for(unsigned i=0;i<element.edgeCount;++i)
        includeEdge(lengths,points[element.edges[i][0]],points[element.edges[i][1]]) ;
      return lengths ;
    }
  }
  struct ParametricCellData::Data {
    // Accepted state owns either a working cell or a compact hex; failed receives own only an error.
    ParametricRecord savedRecord ; // Cached encoding for expanded, non-source cells.
    std::optional<ParametricCell> cell ;
    std::unique_ptr<CompactHexCell> hex ;
    std::string error ;
    int maxDepth=10,maxNodes=1000000 ;
    void retain(ParametricCell value) {
      const auto *corner = value.corner() ;
      if(corner && corner->topology.shape() == Shape::hex) {
        value.validate(-1) ;
        hex = std::make_unique<CompactHexCell>(*corner) ;
        ParametricRecord().swap(savedRecord) ;
      } else {
        if(!value.source()) savedRecord = encodeParametricCell(value) ;
        else ParametricRecord().swap(savedRecord) ;
        cell = std::move(value) ;
      }
    }
  } ;

  ParametricCellData::ParametricCellData(ParametricRecord record,int maxDepth,int maxNodes) {
    auto data=std::make_shared<Data>() ; data->maxDepth=maxDepth ; data->maxNodes=maxNodes ;
    try {
      const auto limits=parametricHistoryLimits(maxDepth,maxNodes) ;
      if(CompactHexCell::saved(record)) {
        auto hex=std::make_unique<CompactHexCell>(CompactHexCell::decode(record,limits)) ;
        if(needsFacePatches(hex->topology,hex->corners))
          throw std::invalid_argument("saved compact hex geometry differs from the source face representation") ;
        data->hex=std::move(hex) ;
      }
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
      const ParametricRecord &record,int maxDepth,int maxNodes,const std::string &requested,int axis,const ParametricEdgeControls &controls) {
    const auto selection=selectParametricCell(topology,corners,requested,axis) ;
    const auto limits=parametricHistoryLimits(maxDepth,maxNodes) ;
    const bool unspecified=controls.empty() || (controls.size()==topology.edgeCount() &&
      std::all_of(controls.begin(),controls.end(),[](const auto &row) { return row.empty() ; })) ;
    if(CompactHexCell::saved(record)) {
      if(selection.facePatches)
        throw std::invalid_argument("saved compact hex geometry differs from the source face representation") ;
      ParametricCellData result(record,maxDepth,maxNodes) ;
      if(!result.data_->hex) throw std::invalid_argument(result.data_->error) ;
      const auto &hex=*result.data_->hex ;
      if(!selection.corner() || axis>=0 || hex.topology!=topology || hex.corners!=corners)
        throw std::invalid_argument("saved compact hex source or split family differs") ;
      if(!unspecified) checkParametricSource(topology,corners,*hex.sourceMap(),controls,limits.jacobian) ;
      return result ;
    }
    if(record.empty() && selection.corner() && unspecified) {
      return ParametricCellData(ParametricCell::original(topology,corners,limits,selection)) ;
    }
    return ParametricCellData(makeParametricCell(topology,corners,record,limits,selection,axis,controls)) ;
  }
  ParametricEdgeControls ParametricCellData::sourceEdgeControls(const CellTopology &topology) const {
    if(compactHex()) return data_->hex->cubic ? parametricSourceEdgeControls(topology,*data_->hex->cubic):ParametricEdgeControls(topology.edgeCount()) ;
    const auto cell=workingCell() ;
    if(const auto *boundary=cell->boundary()) {
      const auto source=boundary->sourceGeometry() ;
      return source ? parametricSourceEdgeControls(topology,*source):ParametricEdgeControls(topology.edgeCount()) ;
    }
    if(cell->source() || cell->extruded()) return ParametricEdgeControls(topology.edgeCount()) ;
    const auto &history=cell->history() ;
    return parametricSourceEdgeControls(topology,*history.geometry().map(history.node(0).map)) ;
  }
  std::shared_ptr<const ParametricCell> ParametricCellData::workingCell() const {
    if(!data_) throw std::invalid_argument("parametric cell is not initialized") ;
    if(data_->hex) return std::make_shared<const ParametricCell>(data_->hex->expand()) ;
    if(!data_->cell) throw std::invalid_argument(data_->error) ;
    // Alias the accepted cell while keeping its owning Data alive.
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
      const IncomingFacePlans &incoming,bool *changed,const GeometricRefinement *selection) const {
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
    if(compactHex()) data->hex=std::make_unique<CompactHexCell>(*data_->hex) ;
    else {
      const auto &source=*data_->cell->source() ;
      data->hex=std::make_unique<CompactHexCell>(source.topology,source.corners,source.limits) ;
    }
    std::array<const FacePlan *,6> faces{} ;
    for(const auto &input:incoming) {
      unsigned f=0 ; while(f<6 && data->hex->topology.faceId(f)!=input.face) ++f ;
      if(f==6 || faces[f]) throw std::invalid_argument("invalid or repeated incoming hex face") ;
      faces[f]=&input.record.get() ;
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
    checkIncomingFaces(incoming) ; auto cell=workingCopy() ;
    for(const auto &face:incoming) cell.mergeFace(face.face,face.record.get()) ;
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
    checkIncomingFaces(incoming) ; auto cell=workingCopy() ; bool faceChanged=false ;
    for(const auto &face:incoming) faceChanged=cell.mergeFace(face.face,face.record.get()) || faceChanged ;
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
  ParametricRecord ParametricCellData::nodeTagRequests(const std::vector<int> &topology,const std::map<ParametricNodeKey,char> &marking) const {
    const auto cell=workingCell() ; auto vertices=parametricLeafVertices(decodeParametricTopology(topology),*cell) ;
    if(cell->hasCornerGeometry()) {
      std::vector<ParametricNodeKey> keys ; keys.reserve(marking.size()) ;
      for(const auto &entry:marking) keys.push_back(entry.first) ;
      mixedCornerVertices(*cell,keys,vertices) ;
    }
    std::vector<char> tags ; tags.reserve(vertices.size()) ;
    for(const auto &leaf:vertices) {
      bool refine=false,coarsen=true ;
      for(const auto &vertex:leaf) {
        const int tag=marking.at(vertex) ;
        if(tag<0 || tag>2) throw std::invalid_argument("unknown parametric node tag") ;
        refine=refine || tag==1 ; coarsen=coarsen && tag==2 ;
      }
      tags.push_back(refine ? 1:coarsen ? 2:0) ;
    }
    return tagRequests(tags) ;
  }
  std::vector<ParametricNodeKey> ParametricCellData::cellNodes(const std::vector<int> &topology,std::vector<Point> *positions) const {
    if(compactHex()) return data_->hex->cellNodes(positions) ;
    const auto cell=workingCell() ; auto keys=parametricCellNodes(decodeParametricTopology(topology),*cell) ;
    if(positions) *positions=parametricCellPositions(*cell,keys) ;
    return keys ;
  }
  std::vector<CellPath> ParametricCellData::leafPaths() const {
    std::vector<CellPath> result ;
    if(const auto *plan=hexPlan()) {
      for(auto leaf:plan->leaves()) {
        CellPath path ;
        while(leaf!=0) { const auto parent=plan->parent(leaf) ; const auto children=plan->children(parent) ;
          path.push_back(std::find(children.begin(),children.end(),leaf)-children.begin()) ; leaf=parent ; }
        std::reverse(path.begin(),path.end()) ; result.push_back(std::move(path)) ;
      }
    } else { const auto cell=workingCell() ; for(auto leaf:cell->leafIds()) result.push_back(cell->path(leaf)) ; }
    return result ;
  }
  std::vector<Point> ParametricCellData::cellPositions(const std::vector<ParametricNodeKey> &keys) const {
    return compactHex() ? data_->hex->cellPositions(keys) : parametricCellPositions(*workingCell(),keys) ;
  }
  std::vector<ParametricOutputFace> ParametricCellData::cellFaces(const std::vector<int> &topology,const std::vector<ParametricNodeKey> &required) const {
    return compactHex() ? data_->hex->cellFaces(required) : parametricCellFaces(decodeParametricTopology(topology),*workingCell(),required) ;
  }
  std::vector<ParametricNodeKey> ParametricCellData::faceNodes(EntityId face,const std::vector<int> &topology,
      const ParametricCellData *right,const std::vector<int> *rightTopology,std::vector<Point> *positions) const {
    if(compactHex() && (!right || right->compactHex())) return data_->hex->faceNodes(face,right ? right->data_->hex.get():nullptr,positions) ;
    const auto left=decodeParametricTopology(topology) ;
    const auto other=rightTopology ? std::make_optional(decodeParametricTopology(*rightTopology)) : std::nullopt ;
    const auto leftCell = workingCell(),rightCell = right ? right->workingCell():nullptr ;
    return parametricFaceNodes(face,left,*leftCell,other ? &*other:nullptr,rightCell.get(),positions) ;
  }
  std::vector<Point> ParametricCellData::facePositions(EntityId face,const std::vector<int> &topology,
      const ParametricCellData *right,const std::vector<int> *rightTopology) const {
    if(compactHex() && (!right || right->compactHex())) return data_->hex->facePositions(face,right ? right->data_->hex.get():nullptr) ;
    const auto left=decodeParametricTopology(topology) ;
    const auto other=rightTopology ? std::make_optional(decodeParametricTopology(*rightTopology)) : std::nullopt ;
    const auto leftCell=workingCell(),rightCell=right ? right->workingCell():nullptr ;
    return parametricFacePositions(face,left,*leftCell,other ? &*other:nullptr,rightCell.get()) ;
  }
  std::vector<ParametricOutputFace> ParametricCellData::faceFaces(EntityId face,const std::vector<int> &topology,
      const std::vector<ParametricNodeKey> &required,const ParametricCellData *right,const std::vector<int> *rightTopology) const {
    if(compactHex() && (!right || right->compactHex())) return data_->hex->faceFaces(face,required,right ? right->data_->hex.get():nullptr) ;
    const auto left=decodeParametricTopology(topology) ;
    const auto other=rightTopology ? std::make_optional(decodeParametricTopology(*rightTopology)) : std::nullopt ;
    const auto leftCell = workingCell(),rightCell = right ? right->workingCell():nullptr ;
    return parametricFaceFaces(face,left,*leftCell,required,other ? &*other:nullptr,rightCell.get()) ;
  }
  std::vector<FacePatch> ParametricCellData::cellPatches(const std::vector<int> &topology) const {
    const auto cell=workingCell() ;
    const auto faces=parametricCellFaces(decodeParametricTopology(topology),*cell,{},true) ;
    std::vector<FacePatch> patches ; patches.reserve(faces.size()) ;
    for(const auto &face:faces) patches.push_back(face.patch) ;
    return patches ;
  }
  std::vector<FacePatch> ParametricCellData::facePatches(EntityId face,const std::vector<int> &topology,
      const ParametricCellData *right,const std::vector<int> *rightTopology) const {
    if(bool(right)!=bool(rightTopology)) throw std::invalid_argument("shared face requires both neighbor state and topology") ;
    const auto left=decodeParametricTopology(topology) ;
    const auto other=rightTopology ? std::make_optional(decodeParametricTopology(*rightTopology)) : std::nullopt ;
    const auto leftCell=workingCell(),rightCell=right ? right->workingCell():nullptr ;
    const auto faces=parametricFaceFaces(face,left,*leftCell,{},other ? &*other:nullptr,rightCell.get(),true) ;
    std::vector<FacePatch> patches ; patches.reserve(faces.size()) ;
    for(const auto &piece:faces) patches.push_back(piece.patch) ;
    return patches ;
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
      const auto leftPlan=facePlan(face),rightPlan=rightData.facePlan(face) ;
      if(leftPlan!=rightPlan) mergeFacePlans(leftPlan,rightPlan) ;
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
      const auto *first=data_->hex ? data_->hex->cubic.get():nullptr ;
      const auto *second=previous.data_->hex ? previous.data_->hex->cubic.get():nullptr ;
      if(bool(first)!=bool(second) || (first && !sameParametricSource(*first,*second)))
        throw std::invalid_argument("cell correspondence requires the same source curves") ;
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
      const GeometryMap &geometry=hex.cubic ? static_cast<const GeometryMap &>(*hex.cubic):source ;
      lengths.reserve(hex.plan.leafCount()) ;
      for(auto leaf:hex.plan.leaves()) lengths.push_back(mapEdgeLengths(HexCellMap(geometry,hex.plan,leaf))) ;
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
          includeEdge(value,cell->cornerPosition(edge.first.first),cell->cornerPosition(edge.first.second)) ;
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
    return compactHex() ? data_->hex->record() : (data_ ? data_->savedRecord : ParametricRecord{}) ;
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
