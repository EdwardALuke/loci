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

#ifndef LOCI_FVMADAPT_GEOMETRY_RECORD_H
#define LOCI_FVMADAPT_GEOMETRY_RECORD_H

#include "subdivision.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <variant>
#include <unordered_map>
#include <vector>

namespace Loci {
namespace FVMAdapt {

  /// Index in one geometry store; it is not a persistent mesh identity.
  using MapId = std::size_t ;

  /// Accepted physical corners in standard reference order or canonical polygon order.
  struct CornerRecord {
    Shape shape ;
    std::vector<Point> corners ;
  } ;

  /// Accepted tensor Bernstein controls in BezierMap order.
  struct TensorBezierRecord {
    Shape shape ;
    std::array<unsigned,3> degree ;
    std::vector<Point> controls ;
  } ;

  /// Ten cubic triangular controls in CubicTriangleMap order.
  struct CubicTriangleRecord { std::vector<Point> controls ; } ;

  /// Shared cubic edges in reference-element order; bit e reverses edge e.
  template<std::size_t EdgeCount> class CubicVolumeRecord {
    struct Data { std::array<MapId,EdgeCount> edges ; std::uint16_t reversed ; } ;
    std::shared_ptr<const Data> data_ ;
    const Data &data() const {
      if(!data_) throw std::invalid_argument("cubic volume has no edges") ;
      return *data_ ;
    }
  public:
    CubicVolumeRecord(std::array<MapId,EdgeCount> edges,std::uint16_t reversed=0)
      : data_(std::make_shared<const Data>(Data{edges,reversed})) {}
    const std::array<MapId,EdgeCount> &edges() const { return data().edges ; }
    std::uint16_t reversed() const { return data().reversed ; }
    std::size_t payloadBytes() const { return sizeof(Data) ; }
  } ;
  using CubicTetRecord = CubicVolumeRecord<6> ;
  using CubicPrismRecord = CubicVolumeRecord<9> ;
  using CubicPyramidRecord = CubicVolumeRecord<8> ;
  using CubicHexRecord = CubicVolumeRecord<12> ;

  template<class T> inline constexpr bool isCubicVolumeRecord = false ;
  template<std::size_t N> inline constexpr bool isCubicVolumeRecord<CubicVolumeRecord<N>> = true ;

  /// Exact affine map with immutable placement data shared by record copies.
  struct AffineRecord {
    Shape shape ;
  private:
    std::shared_ptr<const Placement> placement_ ;
  public:
    AffineRecord(Shape shape, Placement placement)
      : shape(shape), placement_(std::make_shared<const Placement>(std::move(placement))) {}
    const Placement &placement() const {
      if(!placement_) throw std::invalid_argument("affine geometry has no placement") ;
      return *placement_ ;
    }
  } ;

  /// Composition parent(placement(q)); both IDs precede this record.
  struct ComposeRecord {
    MapId parent ;
    MapId placement ;
  } ;

  /// Four oriented edge maps defining a quad surface.
  struct BoundarySurfaceRecord { std::array<MapId,4> edges ; } ;
  /// Cone from a retained surface to an interior point.
  struct ConeRecord { MapId base ; Point apex ; } ;

  /// Piecewise planar face with an accepted physical center.
  struct FanRecord {
    Shape shape ;
  private:
    struct Data { std::vector<Point> corners ; Point center ; } ;
    std::shared_ptr<const Data> data_ ;
    const Data &data() const {
      if(!data_) throw std::invalid_argument("face fan has no controls") ;
      return *data_ ;
    }
  public:
    FanRecord(Shape shape, std::vector<Point> corners, Point center)
      : shape(shape), data_(std::make_shared<const Data>(Data{std::move(corners),center})) {}
    const std::vector<Point> &corners() const { return data().corners ; }
    const Point &center() const { return data().center ; }
  } ;
  /// Face domain in the coordinates of an earlier retained source surface.
  struct SurfacePatchRecord { MapId source ; std::vector<Point> chartCorners ; } ;
  /// Four retained boundary curves joined to one accepted physical center.
  class BoundaryFanRecord {
    struct Data { std::array<MapId,4> edges ; Point center ; } ;
    std::shared_ptr<const Data> data_ ;
    const Data &data() const {
      if(!data_) throw std::invalid_argument("boundary fan has no controls") ;
      return *data_ ;
    }
  public:
    BoundaryFanRecord(std::array<MapId,4> edges,Point center)
      : data_(std::make_shared<const Data>(Data{edges,center})) {}
    const std::array<MapId,4> &edges() const { return data().edges ; }
    const Point &center() const { return data().center ; }
  } ;

  using GeometryRecord = std::variant<CornerRecord,TensorBezierRecord,AffineRecord,ComposeRecord,
                                      BoundarySurfaceRecord,ConeRecord,FanRecord,SurfacePatchRecord,BoundaryFanRecord,
                                      CubicTriangleRecord,CubicHexRecord,CubicTetRecord,CubicPrismRecord,CubicPyramidRecord> ;
  /// Dependencies in construction order; each must precede its user in the store.
  std::vector<MapId> geometryDependencies(const GeometryRecord &record) ;

  /// Bounds on one retained table; maxControls counts physical and chart control points.
  /// exact.maxBytes bounds record/cache payloads, excluding allocator and evaluator-object overhead.
  /// exact.maxOperations bounds estimated work for one map evaluation, including repeated dependencies.
  struct GeometryLimits {
    std::size_t maxRecords = 1000000 ;
    std::size_t maxControls = 1000000 ;
    unsigned maxDegree = 64 ;
    std::size_t maxCompositionDepth = 1024 ;
    ExactLimits exact ;
  } ;

  /// Append-only accepted descriptions and their shared immutable evaluators.
  /// Arbitrary callback maps cannot enter this store by corner sampling.
  class GeometryStore {
    struct CachedMap {
      std::shared_ptr<const GeometryMap> value ;
      std::size_t depth, evaluationCost ;
      std::size_t controls = 0, bytes = 0 ;
    } ;
    GeometryLimits limits_ ;
    std::vector<GeometryRecord> records_ ;
    std::vector<CachedMap> maps_ ;
    std::unordered_map<const GeometryMap*,MapId> captured_ ;
    std::size_t controls_ = 0, bytes_ = 0 ;
  public:
    explicit GeometryStore(GeometryLimits limits = {}) : limits_(limits) {}
    const GeometryLimits &limits() const { return limits_ ; }
    std::size_t size() const { return records_.size() ; }
    std::size_t retainedBytes() const { return bytes_ ; }
    std::size_t evaluationCost(MapId id) const { return maps_.at(id).evaluationCost ; }
    const std::vector<GeometryRecord> &records() const { return records_ ; }
    const GeometryRecord &record(MapId id) const { return records_.at(id) ; }
    std::shared_ptr<const GeometryMap> map(MapId id) const { return maps_.at(id).value ; }
    /// Validate before append; invalid descriptions or exhausted limits leave the store unchanged.
    MapId append(GeometryRecord record) ;
    /// Retain supported geometry descriptions exactly. Unknown map types throw.
    /// Invalid inputs and exhausted limits leave the store unchanged.
    MapId capture(std::shared_ptr<const GeometryMap> map) ;
    /// Validate existing records against new limits, preserving shared evaluators and map IDs.
    /// Rejection leaves the store and its limits unchanged.
    void setLimits(GeometryLimits limits,
                   std::size_t maxBytes = std::numeric_limits<std::size_t>::max()) ;
    /// Remove a suffix, preserving all remaining IDs and maps. Used to roll back an edit.
    void truncate(std::size_t count) ;
    /// Atomically append a dependency-ordered batch, returning its first ID.
    /// References may use earlier records in this batch. An empty batch returns size().
    /// maxBytes may impose a tighter payload ceiling for this append.
    MapId append(std::vector<GeometryRecord> records,
                 std::size_t maxBytes = std::numeric_limits<std::size_t>::max()) ;
    void swap(GeometryStore &other) noexcept ;
  } ;

}
}

#endif
