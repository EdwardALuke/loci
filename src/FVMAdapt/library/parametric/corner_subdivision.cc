// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/parametric/corner_subdivision.h>
#include "geometry_integration.h"
#include "record_io.h"
#include <set>

namespace Loci {
namespace FVMAdapt {
  namespace {
    template<class T> void reserveRecords(std::vector<T> &records,std::size_t count,std::size_t limit) {
      if(count>records.capacity()) {
        const auto growth=std::min(records.capacity(),limit-records.capacity()) ;
        records.reserve(std::max(count,records.capacity()+growth)) ;
      }
    }
    void checkMethod(CenterMethod method) {
      CornerBoundary::require(method==CenterMethod::reference || method==CenterMethod::weighted || method==CenterMethod::integrated,
                              "unknown corner center method") ;
    }
    unsigned methodCode(CenterMethod method) {
      checkMethod(method) ; return method==CenterMethod::reference ? 0 : method==CenterMethod::weighted ? 1 : 2 ;
    }
    CenterMethod readMethod(unsigned code) {
      CornerBoundary::require(code<3,"unknown saved corner center method") ;
      return code==0 ? CenterMethod::reference : code==1 ? CenterMethod::weighted : CenterMethod::integrated ;
    }
  }

  CornerSubdivision::CornerSubdivision(std::uint64_t root,HistoryLimits limits) : root_(root),limits_(limits) {
    require(limits.maxNodes>0,"corner history requires a positive node limit") ;
    options.maxNodes=limits.maxNodes ; options.maxSnapshotBytes=limits.exact.maxBytes ; options.jacobianLimits=limits.jacobian ;
  }
  CornerSubdivision::CornerSubdivision(std::uint64_t root,CornerBoundary boundary,CornerCell cell,Point interior,HistoryLimits limits)
    : CornerSubdivision(root,{},std::move(boundary),std::move(cell),interior,limits) {}
  CornerSubdivision::CornerSubdivision(std::uint64_t root,std::shared_ptr<const GeometryMap> source,
      CornerBoundary boundary,CornerCell cell,Point interior,HistoryLimits limits)
    : CornerSubdivision(root,limits) {
    options.maxEntities=boundary.options.maxEntities ; options.assessment=boundary.options.assessment ;
    boundaries=std::move(boundary.boundaries) ;
    boundaries.geometry.setLimits(savedBoundaryLimits().geometry) ;
    if(source) {
      require(source->dimension()==3,"corner source requires a volume map") ;
      source_=boundaries.geometry.capture(std::move(source)) ;
    }
    nodes_.push_back({noTreeNode,0,noTreeNode,false,std::move(cell),interior}) ; validate() ;
  }
  CornerSubdivision::CornerSubdivision(std::uint64_t root,std::shared_ptr<const GeometryMap> source,HistoryLimits limits)
    : CornerSubdivision(root,limits) {
    require(source && source->dimension()==3,"corner source requires a volume map") ;
    boundaries.geometry.setLimits(savedBoundaryLimits().geometry) ; source_=boundaries.geometry.capture(source) ;
    const auto &element=referenceElement(source->shape()) ; CornerCell cell ; Point interior{} ;
    for(unsigned v=0;v<element.vertexCount;++v) {
      cell.vertices.push_back(addPoint(element.vertices[v])) ;
      for(unsigned d=0;d<3;++d) interior[d]+=element.vertices[v][d]/element.vertexCount ;
    }
    std::map<std::pair<BoundaryId,BoundaryId>,BoundaryId> edges ;
    for(unsigned f=0;f<element.boundaryCount;++f) {
      const auto &face=element.boundaries[f] ; std::vector<BoundaryId> vertices,faceEdges ; std::vector<Point> corners ;
      for(unsigned v=0;v<face.vertexCount;++v) {
        const auto a=face.vertices[v],b=face.vertices[(v+1)%face.vertexCount] ; vertices.push_back(a) ; corners.push_back(element.vertices[a]) ;
        const auto key=endpoints(a,b) ; auto found=edges.find(key) ;
        if(found==edges.end()) found=edges.emplace(key,addEdge(key.first,key.second,std::make_shared<CornerMap>(Shape::segment,
          std::vector<Point>{element.vertices[key.first],element.vertices[key.second]}))).first ;
        faceEdges.push_back(found->second) ;
      }
      cell.faces.push_back({addSurface(std::move(vertices),std::move(faceEdges),
        std::make_shared<CornerMap>(face.shape,std::move(corners))),0}) ;
    }
    nodes_.push_back({noTreeNode,0,noTreeNode,false,std::move(cell),interior}) ; validate() ;
  }
  std::shared_ptr<const GeometryMap> CornerSubdivision::sourceGeometry() const {
    return source_==noBoundary ? nullptr : boundaries.geometry.map(source_) ;
  }
  std::shared_ptr<const GeometryMap> CornerSubdivision::project(std::shared_ptr<const GeometryMap> map) const {
    return source_==noBoundary ? map : std::make_shared<RestrictedMap>(sourceGeometry(),std::move(map)) ;
  }
  Point CornerSubdivision::position(BoundaryId vertex) const {
    const auto &point=boundaries.points.at(vertex) ; return source_==noBoundary ? point : sourceGeometry()->position(point) ;
  }
  std::shared_ptr<const GeometryMap> CornerSubdivision::edgeGeometry(BoundaryId edge) const {
    return project(boundaries.geometry.map(boundaries.edges.at(edge).geometry)) ;
  }
  std::shared_ptr<const GeometryMap> CornerSubdivision::faceGeometry(BoundaryId face) const {
    return project(boundaries.geometry.map(boundaries.faces.at(face).geometry)) ;
  }
  std::vector<std::shared_ptr<const GeometryMap>> CornerSubdivision::cellPieces(TreeNodeId id) const {
    const auto &n=node(id) ; auto pieces=CornerBoundary::cellPieces(n.cell,n.interior) ;
    if(source_!=noBoundary) for(auto &piece:pieces) piece=project(std::move(piece)) ;
    return pieces ;
  }
  GeometryMeasure CornerSubdivision::measureCellGeometry(const CornerCell &cell,const Point &interior) const {
    if(source_==noBoundary) return measureCell(cell,interior) ;
    auto pieces=CornerBoundary::cellPieces(cell,interior) ;
    for(auto &piece:pieces) {
      piece=project(std::move(piece)) ; auto samples=options.assessment.quality ;
      if(samples.samples.empty()) {
        // Interior cone samples avoid the pyramid apex's nonunique derivative.
        for(double z:{.125,.5,.875}) for(double v:{.125,.5,.875}) for(double u:{.125,.5,.875}) {
          if(piece->shape()==Shape::tet && u+v>=1) continue ;
          samples.samples.push_back({{u*(1-z),v*(1-z),z}}) ;
        }
      }
      require(validGeometry(*piece,samples,limits_.jacobian),"mapped corner piece has invalid or unresolved sampled geometry") ;
    }
    const auto result=detail::measureGeometryPieces(pieces,options.assessment.measure) ;
    if(!result.converged) throw std::length_error("mapped corner volume integration did not converge") ;
    return result ;
  }
  BoundaryRecordLimits CornerSubdivision::savedBoundaryLimits() const {
    auto result=boundaryLimits() ; result.geometry=limits_.geometry ; result.geometry.exact=limits_.exact ;
    result.maxBytes=std::min(result.maxBytes,limits_.exact.maxBytes) ; return result ;
  }
  bool CornerSubdivision::active(TreeNodeId id) const {
    node(id) ;
    while(id!=0) { id=node(id).parent ; if(!node(id).refined) return false ; }
    return true ;
  }
  unsigned CornerSubdivision::depth(TreeNodeId id) const {
    unsigned result=0 ; node(id) ;
    while(id!=0) { id=node(id).parent ; ++result ; }
    return result ;
  }
  CellPath CornerSubdivision::path(TreeNodeId id) const {
    CellPath result ; node(id) ;
    while(id!=0) { result.push_back(node(id).ordinal) ; id=node(id).parent ; }
    std::reverse(result.begin(),result.end()) ; return result ;
  }
  TreeNodeId CornerSubdivision::nodeId(const CellPath &path) const {
    TreeNodeId result=0 ;
    for(auto ordinal:path) {
      const auto range=children(result) ; require(ordinal<range.size(),"corner path has no accepted child") ; result=range[ordinal] ;
    }
    return result ;
  }
  NodeRange CornerSubdivision::children(TreeNodeId id) const {
    const auto split=node(id).split ; return split==noTreeNode ? NodeRange() : splits_.at(split).children ;
  }
  std::vector<TreeNodeId> CornerSubdivision::leafIds() const {
    std::vector<TreeNodeId> result,pending{0} ;
    while(!pending.empty()) {
      const auto id=pending.back() ; pending.pop_back() ;
      if(!node(id).refined) result.push_back(id) ;
      else { const auto range=children(id) ; for(auto i=range.size();i>0;--i) pending.push_back(range[i-1]) ; }
    }
    return result ;
  }
  NodeRange CornerSubdivision::refine(TreeNodeId id,CenterMethod method) {
    checkMethod(method) ; require(leaf(id),"corner refinement requires an active leaf") ;
    if(node(id).split!=noTreeNode) {
      const auto range=children(id) ; for(auto child:range) nodes_[child].refined=false ;
      nodes_[id].refined=true ; return range ;
    }
    require(source_==noBoundary || method==CenterMethod::reference,"mapped corner refinement currently requires reference centers") ;
    const auto count=node(id).cell.vertices.size() ;
    if(depth(id)>=limits_.maxDepth || splits_.size()>=limits_.maxSplits || count>limits_.maxNodes-nodes_.size())
      throw std::length_error("corner history limit exceeded") ;
    reserveRecords(nodes_,nodes_.size()+count,limits_.maxNodes) ; reserveRecords(splits_,splits_.size()+1,limits_.maxSplits) ;
    std::vector<std::size_t> selected(boundaries.faces.size(),noBoundary) ;
    for(std::size_t f=0;f<boundaries.faces.size();++f) if(!boundaries.faces[f].splits.empty()) {
      require(boundaries.faces[f].splits.size()==1,"corner history requires one accepted construction per face") ; selected[f]=0 ;
    }
    Transaction transaction(*this) ;
    auto accepted=constructCell(node(id).cell,node(id).interior,method,selected,transaction) ;
    if(source_!=noBoundary) for(unsigned c=0;c<accepted.children.size();++c)
      measureCellGeometry(accepted.children[c],accepted.interiors[c]) ;
    const NodeRange range(nodes_.size(),count) ;
    CornerConstruction construction{id,range,method,accepted.centerVertex,std::move(accepted.faceSplits)} ;
    // Capacity is reserved before geometry changes; publishing only moves accepted records.
    for(unsigned c=0;c<count;++c)
      nodes_.push_back({id,c,noTreeNode,false,std::move(accepted.children[c]),accepted.interiors[c]}) ;
    nodes_[id].split=splits_.size() ; nodes_[id].refined=true ; splits_.push_back(std::move(construction)) ;
    transaction.committed=true ; return range ;
  }
  void CornerSubdivision::coarsen(TreeNodeId id) {
    require(active(id) && node(id).refined,"corner coarsening requires an active parent") ;
    for(auto child:children(id)) require(leaf(child),"corner coarsening requires a complete leaf family") ;
    nodes_[id].refined=false ;
  }

  void CornerSubdivision::validateSplit(const CornerConstruction &split) const {
    checkMethod(split.centerMethod) ; const auto &parent=node(split.parent).cell ;
    require(source_==noBoundary || split.centerMethod==CenterMethod::reference,"mapped corner history requires reference centers") ;
    require(split.children.size()==parent.vertices.size() && split.faceSplits.size()==parent.faces.size(),"corner construction size differs from its parent") ;
    require(split.centerVertex<boundaries.points.size(),"corner construction has no center") ;
    std::map<BoundaryId,unsigned> ordinal ;
    for(unsigned c=0;c<parent.vertices.size();++c) ordinal.emplace(parent.vertices[c],c) ;
    std::vector<std::map<BoundaryId,unsigned>> exterior(parent.vertices.size()) ;
    std::map<BoundaryId,std::vector<BoundaryId>> edgeCenters ;
    for(unsigned f=0;f<parent.faces.size();++f) {
      const auto use=parent.faces[f] ; const auto &face=boundaries.faces.at(use.surface) ;
      require(split.faceSplits[f]<face.splits.size(),"corner construction has no accepted face split") ;
      const auto &cut=face.splits[split.faceSplits[f]] ;
      for(unsigned c=0;c<face.vertices.size();++c)
        exterior[ordinal.at(face.vertices[c])].emplace(cut.children[c],use.orientation>=face.vertices.size() ? 4u:0u) ;
      for(auto edge:face.edges) edgeCenters[edge].push_back(cut.centerVertex) ;
    }
    std::map<BoundaryId,std::vector<std::pair<unsigned,unsigned>>> interior ;
    for(unsigned c=0;c<split.children.size();++c) {
      const auto &child=node(split.children[c]).cell ; auto expected=exterior[c] ;
      for(auto use:child.faces) {
        const auto found=expected.find(use.surface) ;
        if(found==expected.end()) interior[use.surface].push_back({c,use.orientation}) ;
        else { require(use.orientation==found->second,"corner child boundary orientation differs") ; expected.erase(found) ; }
      }
      require(expected.empty(),"corner child omits a parent face patch") ;
    }
    const auto edges=cellEdges(parent) ; std::set<BoundaryId> used ;
    for(const auto &entry:interior) {
      const auto &uses=entry.second ; const auto &face=boundaries.faces.at(entry.first) ;
      require(uses.size()==2 && (uses[0].second<4)!=(uses[1].second<4),"corner interior face requires opposite child uses") ;
      const auto edge=edges.find(endpoints(parent.vertices[uses[0].first],parent.vertices[uses[1].first])) ;
      require(edge!=edges.end() && used.insert(edge->second).second,"corner interior face does not match a parent edge") ;
      const auto &record=boundaries.edges.at(edge->second) ; const auto &centers=edgeCenters.at(edge->second) ;
      require(record.cuts.size()==1 && record.cuts[0].position==Fraction(1,2) && centers.size()==2,"corner parent edge has no midpoint") ;
      const std::set<BoundaryId> vertices{record.cuts[0].vertex,centers[0],split.centerVertex,centers[1]} ;
      require(face.parent==noBoundary && face.vertices.size()==4 && vertices.size()==4 &&
              std::set<BoundaryId>(face.vertices.begin(),face.vertices.end())==vertices,"corner interior face has incorrect vertices") ;
    }
    require(used.size()==edges.size(),"corner construction omits an interior face") ;
  }
  void CornerSubdivision::validate() const {
    require(!nodes_.empty(),"corner history has no root") ;
    if(nodes_.size()>limits_.maxNodes || splits_.size()>limits_.maxSplits) throw std::length_error("corner history limit exceeded") ;
    validateBoundaries(boundaries,savedBoundaryLimits()) ;
    for(const auto &face:boundaries.faces) {
      require(face.splits.size()<=1,"corner history requires one accepted construction per face") ;
      if(source_!=noBoundary) for(const auto &split:face.splits)
        require(split.centerMethod==CenterMethod::reference,"mapped corner history requires reference face centers") ;
    }
    require(node(0).parent==noTreeNode && node(0).ordinal==0,"invalid corner root") ;
    if(source_!=noBoundary) {
      require(source_<boundaries.geometry.size(),"saved corner source map is missing") ;
      const auto source=sourceGeometry() ; require(source->dimension()==3,"corner source requires a volume map") ;
      require(validGeometry(*source,options.assessment.quality,limits_.jacobian),"corner source has invalid or unresolved sampled geometry") ;
      const auto &element=referenceElement(source->shape()) ; const auto &root=node(0).cell ;
      require(root.vertices.size()==element.vertexCount && root.faces.size()==element.boundaryCount,"saved corner source topology differs") ;
      for(unsigned v=0;v<element.vertexCount;++v)
        require(boundaries.points.at(root.vertices[v])==element.vertices[v],"saved corner source coordinates differ") ;
      for(unsigned f=0;f<element.boundaryCount;++f) {
        const auto &reference=element.boundaries[f] ; const auto use=root.faces[f] ; const auto &face=boundaries.faces.at(use.surface) ;
        const auto geometry=boundaries.geometry.map(face.geometry) ; const auto *map=dynamic_cast<const CornerMap *>(geometry.get()) ;
        require(face.parent==noBoundary && face.vertices.size()==reference.vertexCount &&
          use.orientation<orientationCount(reference.shape) && map && map->shape()==reference.shape,
          "saved corner source face differs") ;
        for(unsigned v=0;v<reference.vertexCount;++v) {
          require(face.vertices[orientedVertex(reference.shape,use.orientation,v)]==root.vertices[reference.vertices[v]],
            "saved corner source face incidence differs") ;
          require(map->corners()[v]==boundaries.points.at(face.vertices[v]),"saved corner source face coordinates differ") ;
        }
      }
      for(const auto &entry:cellEdges(root)) {
        const auto &edge=boundaries.edges[entry.second] ; const auto geometry=boundaries.geometry.map(edge.geometry) ;
        const auto *map=dynamic_cast<const CornerMap *>(geometry.get()) ;
        require(edge.parent==noBoundary && map && map->shape()==Shape::segment &&
          map->corners()[0]==boundaries.points.at(edge.vertices[0]) && map->corners()[1]==boundaries.points.at(edge.vertices[1]),
          "saved corner source edge differs") ;
      }
    }
    std::vector<unsigned> depths(nodes_.size()) ; std::vector<bool> owned(nodes_.size()),used(splits_.size()) ; owned[0]=true ;
    for(TreeNodeId id=0;id<nodes_.size();++id) {
      const auto &n=node(id) ;
      if(id!=0) {
        require(n.parent<id,"corner parent must precede its child") ; depths[id]=depths[n.parent]+1 ;
        if(depths[id]>limits_.maxDepth) throw std::length_error("corner history depth exceeded") ;
      }
      measureCell(n.cell,n.interior) ;
      if(source_!=noBoundary) measureCellGeometry(n.cell,n.interior) ;
      if(n.split==noTreeNode) require(!n.refined,"corner leaf has no accepted subdivision") ;
      else {
        require(n.split<splits_.size() && !used[n.split],"invalid corner subdivision link") ; used[n.split]=true ;
        const auto &split=splits_[n.split] ; const auto range=split.children ;
        require(split.parent==id && range.front()>id && range.front()<=nodes_.size() && range.size()<=nodes_.size()-range.front(),
                "invalid corner child range") ;
        for(unsigned c=0;c<range.size();++c) {
          const auto child=range[c] ; require(!owned[child] && node(child).parent==id && node(child).ordinal==c,"invalid corner child ownership") ;
          owned[child]=true ;
          require(n.refined || !node(child).refined,"dormant corner child has an active subdivision") ;
        }
        validateSplit(split) ;
      }
    }
    require(std::find(owned.begin(),owned.end(),false)==owned.end() && std::find(used.begin(),used.end(),false)==used.end(),
            "corner history has unowned records") ;
  }

  std::vector<unsigned char> CornerSubdivision::encode() const {
    validate() ; detail::RecordWriter out(limits_.exact.maxBytes,limits_.exact) ;
    out.u32(0x43535542) ; out.u32(source_==noBoundary ? 1:2) ; out.u64(root_) ;
    std::vector<MapId> mapIds ; const auto geometry=encodeBoundaries(boundaries,mapIds,savedBoundaryLimits()) ;
    out.u64(geometry.size()) ; out.append(geometry) ;
    if(source_!=noBoundary) out.u64(mapIds.at(source_)) ;
    out.u64(nodes_.size()) ;
    for(const auto &n:nodes_) {
      out.u64(n.parent) ; out.u32(n.ordinal) ; out.u64(n.split) ; out.byte(n.refined) ; for(auto x:n.interior) out.real(x) ;
      out.u64(n.cell.vertices.size()) ; for(auto v:n.cell.vertices) out.u64(v) ;
      out.u64(n.cell.faces.size()) ; for(auto use:n.cell.faces) { out.u64(use.surface) ; out.u32(use.orientation) ; }
    }
    out.u64(splits_.size()) ;
    for(const auto &split:splits_) {
      out.u64(split.parent) ; out.u64(split.children.front()) ; out.u64(split.children.size()) ;
      out.u32(methodCode(split.centerMethod)) ; out.u64(split.centerVertex) ; out.u64(split.faceSplits.size()) ; for(auto f:split.faceSplits) out.u64(f) ;
    }
    return out.take() ;
  }
  CornerSubdivision CornerSubdivision::decode(const std::vector<unsigned char> &bytes,HistoryLimits limits) {
    detail::RecordReader in(bytes,limits.exact.maxBytes,limits.exact) ;
    require(in.u32()==0x43535542,"unsupported corner subdivision header") ; const auto version=in.u32() ;
    require(version==1 || version==2,"unsupported corner subdivision version") ;
    CornerSubdivision result(in.u64(),limits) ;
    std::vector<unsigned char> geometry(in.count(limits.exact.maxBytes)) ; for(auto &b:geometry) b=in.byte() ;
    result.boundaries=decodeBoundaries(geometry,result.savedBoundaryLimits()) ;
    if(version==2) {
      const auto source=in.u64() ; require(source<result.boundaries.geometry.size(),"saved corner source map is missing") ;
      result.source_=static_cast<MapId>(source) ;
    }
    const auto count=in.count(limits.maxNodes,61) ; result.nodes_.reserve(count) ;
    for(std::size_t i=0;i<count;++i) {
      CornerNode n ; n.parent=in.u64() ; n.ordinal=in.u32() ; n.split=in.u64() ; const auto active=in.byte() ;
      require(active<=1,"invalid saved corner selection") ; n.refined=active ; for(auto &x:n.interior) x=in.real() ;
      n.cell.vertices.resize(in.count(result.boundaries.points.size(),8)) ; for(auto &v:n.cell.vertices) v=in.u64() ;
      n.cell.faces.resize(in.count(result.boundaries.faces.size(),12)) ; for(auto &use:n.cell.faces) { use.surface=in.u64() ; use.orientation=in.u32() ; }
      result.nodes_.push_back(std::move(n)) ;
    }
    const auto splits=in.count(limits.maxSplits,44) ; result.splits_.reserve(splits) ;
    for(std::size_t i=0;i<splits;++i) {
      const auto parent=in.u64(),first=in.u64(),size=in.u64() ; const auto method=readMethod(in.u32()) ; const auto center=in.u64() ;
      std::vector<std::size_t> faces(in.count(result.boundaries.faces.size(),8)) ; for(auto &f:faces) f=in.u64() ;
      result.splits_.push_back({parent,NodeRange(first,size),method,center,std::move(faces)}) ;
    }
    in.finish() ; result.validate() ; return result ;
  }
}
}
