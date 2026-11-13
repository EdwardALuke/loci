// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_SURFACE_MAP_H
#define LOCI_FVMADAPT_SURFACE_MAP_H
#include "geometry_map.h"

namespace Loci {
namespace FVMAdapt {

  /// A retained surface made of triangles joining each boundary edge to one center.
  /// Reference corners use the native triangle, quad or regular polygon coordinates.
  class FanMap final : public GeometryMap {
    Shape shape_ ;
    std::vector<Point> corners_, reference_ ;
    Point center_ ;
  public:
    FanMap(Shape shape, std::vector<Point> corners, Point center) ;
    Shape shape() const override { return shape_ ; }
    const auto &corners() const { return corners_ ; }
    const Point &center() const { return center_ ; }
    ReferenceCorners referenceCorners() const override { return {reference_.data(),reference_.size()} ; }
    bool containsReference(const Point &, double tolerance = 0) const override ;
    /// Seams use the first incident triangle's one-sided derivative.
    MapValue evaluate(const Point &) const override ;
  } ;

  /// A quadrilateral surface joining four retained edge curves to one center.
  class BoundaryFanMap final : public GeometryMap {
    std::array<std::shared_ptr<const GeometryMap>,4> edges_ ;
    Point center_ ;
  public:
    BoundaryFanMap(std::array<std::shared_ptr<const GeometryMap>,4> edges, Point center) ;
    Shape shape() const override { return Shape::quad ; }
    const auto &edges() const { return edges_ ; }
    const Point &center() const { return center_ ; }
    /// Seams use the first incident sector; the center uses its middle radial direction.
    MapValue evaluate(const Point &) const override ;
  } ;

  /// A face polygon in one retained surface's coordinates. Descendants refer to that same source.
  class SurfacePatchMap final : public GeometryMap {
    std::shared_ptr<const GeometryMap> source_, placement_ ;
    std::vector<Point> corners_ ;
  public:
    /// Source must be a surface, not another patch; chart corners form a convex polygon.
    SurfacePatchMap(std::shared_ptr<const GeometryMap> source, std::vector<Point> chartCorners) ;
    Shape shape() const override { return placement_->shape() ; }
    const auto &source() const { return source_ ; }
    const auto &chartCorners() const { return corners_ ; }
    const auto &placement() const { return placement_ ; }
    ReferenceCorners referenceCorners() const override { return placement_->referenceCorners() ; }
    bool containsReference(const Point &, double tolerance = 0) const override ;
    Point chartPosition(const Point &) const ;
    MapValue evaluate(const Point &) const override ;
    Point position(const Point &) const override ;
  } ;

  /// source(coordinates(q)) and geometry(q) represent the same physical point.
  struct SurfacePiece {
    /// Smooth physical evaluator over a native triangle or quad.
    std::shared_ptr<const GeometryMap> geometry ;
    /// The same piece's coordinates in SurfacePieces::source.
    std::shared_ptr<const GeometryMap> coordinates ;
  } ;
  struct SurfacePieces {
    std::shared_ptr<const GeometryMap> source ; ///< Retained surface shared by all pieces.
    /// Map from the requested surface's local coordinates to the retained source coordinates.
    std::shared_ptr<const GeometryMap> placement ;
    std::vector<SurfacePiece> pieces ;
  } ;

  /// Whether this map contains retained seams requiring separate integration.
  bool hasSurfacePieces(const GeometryMap &) ;
  /// Partition a surface at retained seams without changing its geometry.
  /// Nonlinear generic restrictions crossing seams are unsupported; use SurfacePatchMap.
  SurfacePieces smoothSurfacePieces(std::shared_ptr<const GeometryMap>, std::size_t maxPieces = 4096) ;
  /// Reverse a smooth piece while preserving its stable integration coordinates.
  std::shared_ptr<const GeometryMap> reverseSurfacePiece(std::shared_ptr<const GeometryMap>) ;
  /// Split a retained edge curve at surface seams, preserving its evaluator on each interval.
  std::vector<std::shared_ptr<const GeometryMap>> smoothCurvePieces(std::shared_ptr<const GeometryMap>, std::size_t maxPieces = 4096) ;
  /// Invert a planar chart placement; failure to find a bounded inverse throws invalid_argument.
  Point inverseSurfaceCoordinates(const GeometryMap &placement, const Point &point) ;

}
}
#endif
