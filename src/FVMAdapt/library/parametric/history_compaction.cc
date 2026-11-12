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

#include <FVMAdapt/parametric/refinement_history.h>

#include <algorithm>
#include <map>
#include <tuple>
#include <type_traits>

namespace Loci {
namespace FVMAdapt {
namespace {

  bool mapLess(const GeometryRecord &a, const GeometryRecord &b,
               const std::vector<MapId> &ids) {
    if(a.index()!=b.index()) return a.index()<b.index() ;
    return std::visit([&](const auto &first) {
      using T=std::decay_t<decltype(first)> ;
      const auto &second=std::get<T>(b) ;
      if constexpr(std::is_same_v<T,CornerRecord>) {
        return std::tie(first.shape,first.corners)<std::tie(second.shape,second.corners) ;
      } else if constexpr(std::is_same_v<T,TensorBezierRecord>) {
        return std::tie(first.shape,first.degree,first.controls)<
               std::tie(second.shape,second.degree,second.controls) ;
      } else if constexpr(std::is_same_v<T,CubicTriangleRecord>) {
        return first.controls < second.controls ;
      } else if constexpr(std::is_same_v<T,AffineRecord>) {
        return std::tie(first.shape,first.placement().origin,first.placement().columns)<
               std::tie(second.shape,second.placement().origin,second.placement().columns) ;
      } else if constexpr(std::is_same_v<T,ConeRecord>) {
        return std::tie(ids[first.base],first.apex)<std::tie(ids[second.base],second.apex) ;
      } else if constexpr(std::is_same_v<T,FanRecord>) {
        return std::tie(first.shape,first.corners(),first.center())<std::tie(second.shape,second.corners(),second.center()) ;
      } else if constexpr(std::is_same_v<T,SurfacePatchRecord>) {
        return std::tie(ids[first.source],first.chartCorners)<std::tie(ids[second.source],second.chartCorners) ;
      } else if constexpr(std::is_same_v<T,BoundaryFanRecord>) {
        auto a=first.edges(),b=second.edges() ;
        for(auto &id:a) id=ids[id] ;
        for(auto &id:b) id=ids[id] ;
        return std::tie(a,first.center())<std::tie(b,second.center()) ;
      } else if constexpr(isCubicVolumeRecord<T>) {
        auto a=first.edges(),b=second.edges() ;
        for(auto &id:a) id=ids[id] ;
        for(auto &id:b) id=ids[id] ;
        if(a!=b) return a<b ;
        return first.reversed()<second.reversed() ;
      } else {
        auto firstIds=geometryDependencies(a), secondIds=geometryDependencies(b) ;
        for(auto &id:firstIds) id=ids[id] ;
        for(auto &id:secondIds) id=ids[id] ;
        return firstIds<secondIds ;
      }
    },a) ;
  }

  GeometryRecord remapRecord(GeometryRecord record, const std::vector<MapId> &ids) {
    std::visit([&](auto &value) {
      using T=std::decay_t<decltype(value)> ;
      if constexpr(std::is_same_v<T,ComposeRecord>) {
        value.parent=ids[value.parent] ; value.placement=ids[value.placement] ;
      } else if constexpr(std::is_same_v<T,BoundarySurfaceRecord>) {
        for(auto &edge:value.edges) edge=ids[edge] ;
      } else if constexpr(std::is_same_v<T,BoundaryFanRecord>) {
        auto edges=value.edges() ; auto center=value.center() ;
        for(auto &edge:edges) edge=ids[edge] ;
        for(auto &coordinate:center) if(coordinate==0) coordinate=0 ;
        value=BoundaryFanRecord(edges,center) ;
      } else if constexpr(isCubicVolumeRecord<T>) {
        auto edges=value.edges() ;
        for(auto &edge:edges) edge=ids[edge] ;
        value=T(edges,value.reversed()) ;
      } else if constexpr(std::is_same_v<T,ConeRecord>) {
        value.base=ids[value.base] ;
        for(double &coordinate:value.apex) if(coordinate==0) coordinate=0 ;
      } else if constexpr(std::is_same_v<T,SurfacePatchRecord>) {
        value.source=ids[value.source] ;
        for(auto &point:value.chartCorners) for(double &coordinate:point)
          if(coordinate==0) coordinate=0 ;
      } else if constexpr(std::is_same_v<T,FanRecord>) {
        auto corners=value.corners() ; auto center=value.center() ;
        for(auto &point:corners) for(double &coordinate:point) if(coordinate==0) coordinate=0 ;
        for(double &coordinate:center) if(coordinate==0) coordinate=0 ;
        value=FanRecord(value.shape,std::move(corners),center) ;
      } else if constexpr(std::is_same_v<T,CornerRecord> || std::is_same_v<T,TensorBezierRecord> ||
                          std::is_same_v<T,CubicTriangleRecord>) {
        auto &points=[&]() -> std::vector<Point>& {
          if constexpr(std::is_same_v<T,CornerRecord>) return value.corners ;
          else return value.controls ;
        }() ;
        for(auto &point:points) for(double &coordinate:point)
          if(coordinate==0) coordinate=0 ;
      }
    },record) ;
    return record ;
  }

  GeometryStore compactGeometry(const GeometryStore &geometry, std::vector<CellNode> &nodes,
                                const GeometryLimits &limits) {
    std::vector<bool> retained(geometry.size(),false) ;
    for(const auto &node:nodes) retained[node.map]=true ;
    // Dependencies precede their users, so one reverse pass closes the retained set.
    for(MapId id=geometry.size(); id-- >0;)
      if(retained[id]) for(auto dependency:geometryDependencies(geometry.record(id))) retained[dependency]=true ;

    std::vector<std::size_t> depths(geometry.size(),0) ;
    std::map<std::size_t,std::vector<MapId>> levels ;
    for(MapId id=0; id<geometry.size(); ++id) if(retained[id]) {
      for(auto dependency:geometryDependencies(geometry.record(id))) depths[id]=std::max(depths[id],depths[dependency]+1) ;
      levels[depths[id]].push_back(id) ;
    }
    GeometryStore result(limits) ;
    std::vector<MapId> ids(geometry.size(),noTreeNode) ;
    const auto less=[&](MapId a, MapId b) {
      return mapLess(geometry.record(a),geometry.record(b),ids) ;
    } ;
    for(auto &level:levels) {
      auto &records=level.second ;
      std::sort(records.begin(),records.end(),less) ;
      MapId previous=noTreeNode ;
      for(MapId id:records) {
        if(previous==noTreeNode || less(previous,id))
          result.append(remapRecord(geometry.record(id),ids)) ;
        ids[id]=result.size()-1 ;
        previous=id ;
      }
    }
    for(auto &node:nodes) node.map=ids[node.map] ;
    return result ;
  }

}

  RefinementHistory compactHistory(const RefinementHistory &history) {
    std::vector<TreeNodeId> oldIds(1,0) ;
    std::vector<CellNode> nodes(1,history.node(0)) ;
    std::vector<AcceptedSplit> splits ;
    nodes[0].activeSplit=noTreeNode ;
    // Breadth-first child order gives each selected construction a contiguous range.
    for(TreeNodeId parent=0; parent<nodes.size(); ++parent) {
      const auto &oldNode=history.node(oldIds[parent]) ;
      if(oldNode.isLeaf()) continue ;
      const auto &oldSplit=history.splits()[oldNode.activeSplit] ;
      const TreeNodeId first=nodes.size() ;
      nodes[parent].activeSplit=splits.size() ;
      splits.push_back({parent,NodeRange(first,oldSplit.children.size()),
                        oldSplit.operation,oldSplit.construction}) ;
      for(TreeNodeId oldChild:oldSplit.children) {
        CellNode child=history.node(oldChild) ;
        child.parent=parent ; child.activeSplit=noTreeNode ;
        oldIds.push_back(oldChild) ; nodes.push_back(std::move(child)) ;
      }
    }
    auto geometry=compactGeometry(history.geometry(),nodes,history.limits().geometry) ;
    return RefinementHistory::restore(history.rootKey(),history.geometryRevision(),
                                      std::move(geometry),std::move(nodes),std::move(splits),
                                      history.limits()) ;
  }

}
}
