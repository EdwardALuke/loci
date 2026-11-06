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

#ifndef LOCI_FVMADAPT_GEOMETRY_MAP_H
#define LOCI_FVMADAPT_GEOMETRY_MAP_H

#include "reference_element.h"

#include <memory>
#include <vector>

namespace Loci {
namespace FVMAdapt {

  /// Position and Jacobian columns with respect to reference coordinates.
  /// Unused columns are zero. All columns are zero and unavailable when hasJacobian is false.
  struct MapValue {
    Point position ;
    std::array<Point,3> jacobian ;
    bool hasJacobian ;
  } ;

  /// Trilinear interpolation in native hex corner order, without allocating a map object.
  /// Corners must be finite. Reference coordinates may extrapolate.
  MapValue evaluateHexCorners(const std::array<Point,8> &corners, const Point &point) ;
  /// Position-only form of evaluateHexCorners, with the same coordinate contract.
  Point hexCornerPosition(const std::array<Point,8> &corners, const Point &point) ;

  /// Read-only corners owned by a map or a standard reference table.
  struct ReferenceCorners {
    const Point *data ; std::size_t count ;
    const Point *begin() const { return data ; }
    const Point *end() const { return data+count ; }
    std::size_t size() const { return count ; }
    const Point &operator[](std::size_t i) const { return data[i] ; }
  } ;

  /// Immutable local geometry evaluator; topology and subdivision state are separate.
  class GeometryMap {
  public:
    virtual ~GeometryMap() = default ;
    virtual Shape shape() const = 0 ;
    unsigned dimension() const { return shape()==Shape::polygon ? 2:referenceElement(shape()).dimension ; }
    /// Ordered reference corners; polygon maps supply their own reference domain.
    virtual ReferenceCorners referenceCorners() const ;
    virtual bool containsReference(const Point &, double tolerance = 0) const ;
    /// Evaluate finite reference coordinates, allowing extrapolation; unused components are ignored.
    /// Nonfinite coordinates throw invalid_argument. Use containsReference() separately for domain membership.
    /// Evaluation does not certify injectivity or positive orientation.
    virtual MapValue evaluate(const Point &point) const = 0 ;
    /// Evaluate only the position, with the same coordinate contract as evaluate().
    virtual Point position(const Point &point) const { return evaluate(point).position ; }
  } ;

  /// Regular reference polygon, centered at (1/2,1/2), with radius 1/2.
  std::vector<Point> polygonReferenceCorners(unsigned vertices) ;
  bool containsReferencePolygon(unsigned vertices, const Point &, double tolerance = 0) ;
  /// Strict interior test for accepted face split points.
  bool insideReferencePolygon(unsigned vertices, const Point &) ;

  /// Rational polygon surface using Wachspress coordinates on a regular reference polygon.
  /// Interpolates its corners and straight edges; warped corners define a retained curved surface.
  class PolygonMap final : public GeometryMap {
    const std::vector<Point> corners_, reference_ ;
  public:
    explicit PolygonMap(std::vector<Point> corners) ;
    Shape shape() const override { return Shape::polygon ; }
    const auto &corners() const { return corners_ ; }
    ReferenceCorners referenceCorners() const override { return {reference_.data(),reference_.size()} ; }
    bool containsReference(const Point &, double tolerance = 0) const override ;
    MapValue evaluate(const Point &) const override ;
  } ;

  /// Corner interpolation: affine simplices, bilinear quad, trilinear hex,
  /// triangle-times-linear prism and rational pyramid. Corner order is referenceElement(shape).vertices.
  class CornerMap final : public GeometryMap {
    const Shape shape_ ;
    const std::vector<Point> corners_ ;
    const Point origin_ ;
  public:
    /// Invalid shapes, corner counts or nonfinite coordinates throw invalid_argument.
    CornerMap(Shape shape, std::vector<Point> corners) ;
    Shape shape() const override { return shape_ ; }
    const std::vector<Point> &corners() const { return corners_ ; }
    /// At the pyramid apex the position is defined, but the rational basis has no unique derivative.
    /// Other points with z == 1 are singular and throw invalid_argument.
    MapValue evaluate(const Point &point) const override ;
    Point position(const Point &point) const override ;
  } ;

  /// Tensor-product Bernstein map for a segment, quad or hex. Controls are Bezier controls, not samples.
  /// Control (i,j,k) has index i+(degree[0]+1)*(j+(degree[1]+1)*k).
  class BezierMap final : public GeometryMap {
    const Shape shape_ ;
    const std::array<unsigned,3> degree_ ;
    const std::vector<Point> controls_ ;
    const Point origin_ ;
  public:
    /// Active degrees must be >= 1; unused degrees must be zero. Invalid shape, degrees,
    /// control count or nonfinite control coordinates throw invalid_argument.
    BezierMap(Shape shape, std::array<unsigned,3> degree, std::vector<Point> controls) ;
    Shape shape() const override { return shape_ ; }
    const std::array<unsigned,3> &degree() const { return degree_ ; }
    const std::vector<Point> &controls() const { return controls_ ; }
    MapValue evaluate(const Point &point) const override ;
  } ;

  /// Composition X(R(q)): referenceMap places a child or boundary in the parent's reference coordinates.
  /// Both maps remain alive through shared const ownership; the parent geometry is never refitted.
  class RestrictedMap final : public GeometryMap {
    const std::shared_ptr<const GeometryMap> parent_ ;
    const std::shared_ptr<const GeometryMap> referenceMap_ ;
  public:
    /// Null maps or a child dimension larger than the parent throw invalid_argument.
    /// The caller must ensure a valid placement; composition does not certify containment or injectivity.
    RestrictedMap(std::shared_ptr<const GeometryMap> parent,
                  std::shared_ptr<const GeometryMap> referenceMap) ;
    Shape shape() const override { return referenceMap_->shape() ; }
    ReferenceCorners referenceCorners() const override { return referenceMap_->referenceCorners() ; }
    bool containsReference(const Point &q, double tolerance = 0) const override { return referenceMap_->containsReference(q,tolerance) ; }
    std::shared_ptr<const GeometryMap> parent() const { return parent_ ; }
    std::shared_ptr<const GeometryMap> referenceMap() const { return referenceMap_ ; }
    /// Compose Jacobians by the chain rule; unavailable derivatives remain unavailable.
    MapValue evaluate(const Point &point) const override ;
    Point position(const Point &point) const override ;
  } ;

  /// Quad surface interpolating four edge curves in counterclockwise boundary order.
  class BoundarySurfaceMap final : public GeometryMap {
    const std::array<std::shared_ptr<const GeometryMap>,4> edges_ ;
    std::array<Point,4> corners_ ;
  public:
    /// Edges are segment maps with matching consecutive endpoints.
    explicit BoundarySurfaceMap(std::array<std::shared_ptr<const GeometryMap>,4> edges) ;
    Shape shape() const override { return Shape::quad ; }
    const auto &edges() const { return edges_ ; }
    const auto &corners() const { return corners_ ; }
    MapValue evaluate(const Point &point) const override ;
  } ;

  namespace detail {
    /// Shared Coons evaluator; retained pieces use their own one-sided edge derivatives.
    MapValue evaluateBoundarySurface(const std::array<std::shared_ptr<const GeometryMap>,4> &edges,
                                    const std::array<Point,4> &corners,const Point &point) ;
  }

  /// Linear cone from a triangle or quad surface to an apex at reference z=1.
  /// The base orientation points into the cone for a positive volume Jacobian.
  class ConeMap final : public GeometryMap {
    const std::shared_ptr<const GeometryMap> base_ ;
    const Point apex_ ;
  public:
    ConeMap(std::shared_ptr<const GeometryMap> base, Point apex) ;
    Shape shape() const override { return base_->shape()==Shape::triangle ? Shape::tet : Shape::pyramid ; }
    const auto &base() const { return base_ ; }
    const Point &apex() const { return apex_ ; }
    /// At the apex the position is defined, but its derivative is unavailable.
    MapValue evaluate(const Point &point) const override ;
  } ;

}
}

#endif
