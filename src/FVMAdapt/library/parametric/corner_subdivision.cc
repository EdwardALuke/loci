// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/parametric/corner_subdivision.h>
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
    : CornerSubdivision(root,limits) {
    options.maxEntities=boundary.options.maxEntities ; options.assessment=boundary.options.assessment ;
    boundaries=std::move(boundary.boundaries) ;
    boundaries.geometry.setLimits(savedBoundaryLimits().geometry) ;
    nodes_.push_back({noTreeNode,0,noTreeNode,false,std::move(cell),interior}) ; validate() ;
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
    for(const auto &face:boundaries.faces) require(face.splits.size()<=1,"corner history requires one accepted construction per face") ;
    require(node(0).parent==noTreeNode && node(0).ordinal==0,"invalid corner root") ;
    std::vector<unsigned> depths(nodes_.size()) ; std::vector<bool> owned(nodes_.size()),used(splits_.size()) ; owned[0]=true ;
    for(TreeNodeId id=0;id<nodes_.size();++id) {
      const auto &n=node(id) ;
      if(id!=0) {
        require(n.parent<id,"corner parent must precede its child") ; depths[id]=depths[n.parent]+1 ;
        if(depths[id]>limits_.maxDepth) throw std::length_error("corner history depth exceeded") ;
      }
      measureCell(n.cell,n.interior) ;
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
    out.u32(0x43535542) ; out.u32(1) ; out.u64(root_) ;
    const auto geometry=encodeBoundaries(boundaries,savedBoundaryLimits()) ; out.u64(geometry.size()) ; out.append(geometry) ;
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
    require(in.u32()==0x43535542 && in.u32()==1,"unsupported corner subdivision version") ;
    CornerSubdivision result(in.u64(),limits) ;
    std::vector<unsigned char> geometry(in.count(limits.exact.maxBytes)) ; for(auto &b:geometry) b=in.byte() ;
    result.boundaries=decodeBoundaries(geometry,result.savedBoundaryLimits()) ;
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
