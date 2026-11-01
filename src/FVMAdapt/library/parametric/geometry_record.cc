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

#include <FVMAdapt/parametric/geometry_record.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace Loci {
namespace FVMAdapt {
  namespace {

    template<class T> void reserve(std::vector<T> &values, std::size_t count) {
      if(count <= values.capacity()) return ;
      if(count > values.max_size())
        throw std::length_error("geometry table exceeds container capacity") ;
      const std::size_t growth = values.capacity() > values.max_size()/2 ?
        values.max_size() : 2*values.capacity() ;
      values.reserve(std::max(count,growth)) ;
    }

    std::size_t integerBytes(Fraction::Integer value) {
      if(value < 0) value = -value ;
      if(value == 0) return 0 ;
      const std::size_t bits = std::size_t(boost::multiprecision::msb(value))+1 ;
      return bits/8+(bits%8 != 0) ;
    }

    template<class Retain> void checkPlacement(const Placement &placement,
                                              const ExactLimits &limits, Retain retain) {
      const auto check = [&](const ReferencePoint &point) {
        for(const Fraction &value : point) {
          checkExact(value,limits) ;
          if(!std::isfinite(value.value()))
            throw std::invalid_argument("affine geometry must have finite numerical coordinates") ;
          retain(integerBytes(value.numerator())) ;
          retain(integerBytes(value.denominator())) ;
        }
      } ;
      check(placement.origin) ;
      for(const ReferencePoint &column : placement.columns) check(column) ;
    }

  }

  std::vector<MapId> geometryDependencies(const GeometryRecord &record) {
    return std::visit([](const auto &value) -> std::vector<MapId> {
      using T=std::decay_t<decltype(value)> ;
      if constexpr(std::is_same_v<T,ComposeRecord>) return {value.parent,value.placement} ;
      else if constexpr(std::is_same_v<T,BoundarySurfaceRecord>) return {value.edges.begin(),value.edges.end()} ;
      else if constexpr(std::is_same_v<T,ConeRecord>) return {value.base} ;
      else return {} ;
    },record) ;
  }

  MapId GeometryStore::append(GeometryRecord record) {
    std::vector<GeometryRecord> records ;
    records.push_back(std::move(record)) ;
    return append(std::move(records)) ;
  }

  void GeometryStore::setLimits(GeometryLimits limits, std::size_t maxBytes) {
    if(size() > limits.maxRecords || controls_ > limits.maxControls ||
       bytes_ > std::min(maxBytes,limits.exact.maxBytes))
      throw std::length_error("geometry table exceeds its destination limits") ;
    for(std::size_t id=0; id<records_.size(); ++id) {
      if(maps_[id].depth > limits.maxCompositionDepth || maps_[id].evaluationCost > limits.exact.maxOperations)
        throw std::length_error("geometry evaluation exceeds its destination limits") ;
      if(const auto *tensor=std::get_if<TensorBezierRecord>(&records_[id])) {
        for(unsigned degree : tensor->degree) if(degree > limits.maxDegree)
          throw std::length_error("geometry degree limit exceeded") ;
      } else if(const auto *affine=std::get_if<AffineRecord>(&records_[id])) {
        checkPlacement(affine->placement(),limits.exact,[](std::size_t) {}) ;
      }
    }
    limits_=limits ;
  }

  MapId GeometryStore::capture(std::shared_ptr<const GeometryMap> map) {
    const auto first = size() ;
    std::vector<GeometryRecord> records ;
    std::vector<std::shared_ptr<const GeometryMap>> sources ;
    std::map<const GeometryMap*,MapId> ids ;
    const auto retain = [&](GeometryRecord record, std::shared_ptr<const GeometryMap> source) {
      if(records.size() >= limits_.maxRecords-size())
        throw std::length_error("geometry record limit exceeded") ;
      const MapId id = size()+records.size() ;
      records.push_back(std::move(record)) ; ids.emplace(source.get(),id) ; sources.push_back(std::move(source)) ; return id ;
    } ;
    std::function<MapId(const std::shared_ptr<const GeometryMap>&,std::size_t)> visit ;
    visit = [&](const std::shared_ptr<const GeometryMap> &value, std::size_t depth) -> MapId {
      if(!value) throw std::invalid_argument("cannot retain a null geometry map") ;
      const auto saved = captured_.find(value.get()) ;
      if(saved != captured_.end()) return saved->second ;
      const auto found = ids.find(value.get()) ;
      if(found != ids.end()) return found->second ;
      if(depth > limits_.maxCompositionDepth)
        throw std::length_error("geometry composition depth limit exceeded") ;
      if(const auto *polygon = dynamic_cast<const PolygonMap*>(value.get()))
        return retain(CornerRecord{Shape::polygon,polygon->corners()},value) ;
      if(const auto *corner = dynamic_cast<const CornerMap*>(value.get()))
        return retain(CornerRecord{corner->shape(),corner->corners()},value) ;
      if(const auto *bezier = dynamic_cast<const BezierMap*>(value.get()))
      {
        if(bezier->controls().size() > limits_.maxControls-controls_)
          throw std::length_error("geometry control limit exceeded") ;
        for(unsigned degree : bezier->degree()) if(degree > limits_.maxDegree)
          throw std::length_error("geometry degree limit exceeded") ;
        return retain(TensorBezierRecord{bezier->shape(),bezier->degree(),bezier->controls()},value) ;
      }
      if(const auto *composed = dynamic_cast<const RestrictedMap*>(value.get())) {
        const MapId parent = visit(composed->parent(),depth+1) ;
        const MapId placement = visit(composed->referenceMap(),depth+1) ;
        return retain(ComposeRecord{parent,placement},value) ;
      }
      if(const auto *surface=dynamic_cast<const BoundarySurfaceMap*>(value.get())) {
        BoundarySurfaceRecord record ;
        for(unsigned e=0; e<4; ++e) record.edges[e]=visit(surface->edges()[e],depth+1) ;
        return retain(record,value) ;
      }
      if(const auto *cone=dynamic_cast<const ConeMap*>(value.get()))
        return retain(ConeRecord{visit(cone->base(),depth+1),cone->apex()},value) ;
      throw std::invalid_argument("geometry map has no retained description") ;
    } ;
    const MapId result = visit(map,0) ; append(std::move(records)) ;
    try {
      for(std::size_t i=0; i<sources.size(); ++i) {
        captured_.erase(maps_[first+i].value.get()) ;
        maps_[first+i].value=std::move(sources[i]) ;
        captured_.emplace(maps_[first+i].value.get(),first+i) ;
      }
    } catch(...) { truncate(first) ; throw ; }
    return result ;
  }

  MapId GeometryStore::append(std::vector<GeometryRecord> records, std::size_t maxBytes) {
    maxBytes = std::min(maxBytes,limits_.exact.maxBytes) ;
    if(bytes_ > maxBytes) throw std::length_error("geometry payload limit exceeded") ;
    const MapId first = size() ;
    if(records.size() > limits_.maxRecords-first)
      throw std::length_error("geometry record limit exceeded") ;
    std::size_t controls = controls_, bytes = bytes_ ;
    const auto retain = [&](std::size_t count) {
      if(count > maxBytes-bytes)
        throw std::length_error("geometry payload limit exceeded") ;
      bytes += count ;
    } ;
    const auto countControls = [&](const std::vector<Point> &points) {
      if(points.size() > limits_.maxControls-controls)
        throw std::length_error("geometry control limit exceeded") ;
      // Evaluators keep a second immutable copy of their physical controls.
      if(points.size() > (maxBytes-bytes)/(2*sizeof(Point)))
        throw std::length_error("geometry payload limit exceeded") ;
      retain(2*sizeof(Point)*points.size()) ;
      controls += points.size() ;
      for(const Point &point : points)
        for(double value : point)
          if(!std::isfinite(value))
            throw std::invalid_argument("geometry controls must be finite") ;
    } ;
    std::vector<std::pair<std::size_t,std::size_t>> retained ;
    retained.reserve(records.size()) ;
    for(const GeometryRecord &record : records) {
      retain(sizeof(GeometryRecord)+sizeof(CachedMap)+sizeof(std::pair<const GeometryMap* const,MapId>)) ;
      std::visit([&](const auto &description) {
        using Record = std::decay_t<decltype(description)> ;
        if constexpr(std::is_same<Record,CornerRecord>::value) {
          countControls(description.corners) ;
          if(description.shape==Shape::polygon) retain(description.corners.size()*sizeof(Point)) ;
        } else if constexpr(std::is_same<Record,TensorBezierRecord>::value) {
          for(unsigned degree : description.degree)
            if(degree > limits_.maxDegree)
              throw std::length_error("geometry degree limit exceeded") ;
          countControls(description.controls) ;
        } else if constexpr(std::is_same<Record,AffineRecord>::value) {
          retain(sizeof(Placement)) ;
          checkPlacement(description.placement(),limits_.exact,retain) ;
        } else if constexpr(std::is_same_v<Record,ConeRecord>) {
          countControls({description.apex}) ;
        }
      },record) ;
      retained.emplace_back(controls,bytes) ;
    }

    std::vector<CachedMap> maps ;
    maps.reserve(records.size()) ;
    for(const GeometryRecord &record : records) {
      const auto lookup = [&](MapId id) -> const CachedMap & {
        if(id < first) return maps_[id] ;
        if(id-first < maps.size()) return maps[id-first] ;
        throw std::invalid_argument("geometry references must precede their composition") ;
      } ;
      maps.push_back(std::visit([&](const auto &description) -> CachedMap {
        using Record = std::decay_t<decltype(description)> ;
        std::size_t cost = 0 ;
        const auto work = [&](std::size_t amount) {
          if(amount > limits_.exact.maxOperations-cost)
            throw std::length_error("geometry evaluation work limit exceeded") ;
          cost += amount ;
        } ;
        work(1) ;
        if constexpr(std::is_same<Record,CornerRecord>::value) {
          work(description.corners.size()) ;
          if(description.shape==Shape::polygon) {
            const auto n=description.corners.size() ;
            if(n && n>(limits_.exact.maxOperations-cost)/n) throw std::length_error("polygon evaluation work limit exceeded") ;
            work(n*n) ; return {std::make_shared<PolygonMap>(description.corners),0,cost} ;
          }
          return {std::make_shared<CornerMap>(description.shape,description.corners),0,cost} ;
        } else if constexpr(std::is_same<Record,TensorBezierRecord>::value) {
          work(description.controls.size()) ;
          for(std::size_t degree : description.degree) {
            if(degree && degree > (limits_.exact.maxOperations-cost)/degree)
              throw std::length_error("geometry evaluation work limit exceeded") ;
            work(degree*degree) ;
          }
          return {std::make_shared<BezierMap>(description.shape,description.degree,description.controls),0,cost} ;
        } else if constexpr(std::is_same<Record,AffineRecord>::value) {
          return {description.placement().map(description.shape),0,cost} ;
        } else {
          std::size_t depth=0 ;
          for(auto id:geometryDependencies(record)) {
            const auto &dependency=lookup(id) ; depth=std::max(depth,dependency.depth) ; work(dependency.evaluationCost) ;
          }
          if(depth >= limits_.maxCompositionDepth)
            throw std::length_error("geometry composition depth limit exceeded") ;
          if constexpr(std::is_same_v<Record,ComposeRecord>)
            return {std::make_shared<RestrictedMap>(lookup(description.parent).value,lookup(description.placement).value),depth+1,cost} ;
          else if constexpr(std::is_same_v<Record,BoundarySurfaceRecord>) {
            std::array<std::shared_ptr<const GeometryMap>,4> edges ;
            for(unsigned e=0; e<4; ++e) edges[e]=lookup(description.edges[e]).value ;
            work(4) ; return {std::make_shared<BoundarySurfaceMap>(std::move(edges)),depth+1,cost} ;
          } else
            return {std::make_shared<ConeMap>(lookup(description.base).value,description.apex),depth+1,cost} ;
        }
      },record)) ;
      maps.back().controls = retained[maps.size()-1].first ;
      maps.back().bytes = retained[maps.size()-1].second ;
    }

    // All allocation and validation precede the no-throw append.
    static_assert(std::is_nothrow_move_constructible<GeometryRecord>::value,
                  "geometry records must support a no-throw commit") ;
    reserve(maps_,first+maps.size()) ;
    reserve(records_,first+records.size()) ;
    try {
      for(std::size_t i=0; i<maps.size(); ++i) captured_.emplace(maps[i].value.get(),first+i) ;
    } catch(...) {
      for(const auto &map : maps) captured_.erase(map.value.get()) ;
      throw ;
    }
    for(GeometryRecord &record : records) records_.push_back(std::move(record)) ;
    for(CachedMap &map : maps) maps_.push_back(std::move(map)) ;
    controls_ = controls ;
    bytes_ = bytes ;
    return first ;
  }

  void GeometryStore::truncate(std::size_t count) {
    if(count > size()) throw std::out_of_range("geometry truncation exceeds the table") ;
    controls_ = count ? maps_[count-1].controls : 0 ;
    bytes_ = count ? maps_[count-1].bytes : 0 ;
    for(std::size_t id=count; id<maps_.size(); ++id) {
      captured_.erase(maps_[id].value.get()) ;
    }
    maps_.resize(count) ; records_.resize(count) ;
  }

  void GeometryStore::swap(GeometryStore &other) noexcept {
    using std::swap ;
    swap(limits_,other.limits_) ;
    records_.swap(other.records_) ;
    maps_.swap(other.maps_) ;
    captured_.swap(other.captured_) ;
    swap(controls_,other.controls_) ;
    swap(bytes_,other.bytes_) ;
  }

}
}
