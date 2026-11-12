// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_HEX_PLAN_H
#define LOCI_FVMADAPT_HEX_PLAN_H

#include "geometry_center.h"
#include "vertex_volume_split.h"
#include <cstdint>
#include <limits>

namespace Loci {
namespace FVMAdapt {

  /// Accepted centers in face-local and cell-local coordinates. Defaults give midpoint splits.
  struct HexCenters {
    std::array<std::array<double,2>,6> faces{{{{.5,.5}},{{.5,.5}},{{.5,.5}},{{.5,.5}},{{.5,.5}},{{.5,.5}}}} ;
    Point cell{{.5,.5,.5}} ;
    bool operator==(const HexCenters &other) const { return faces == other.faces && cell == other.cell ; }
  } ;

  /// One root's eight-child subdivisions. Source geometry and mesh ownership belong to the caller.
  /// Coarsening retains accepted centers and child identities. An empty plan allocates no storage.
  class HexPlan {
  public:
    using Cell = std::uint32_t ;
    static constexpr Cell root = 0 ;
    static constexpr unsigned maxDepth = 30 ;
    static constexpr std::uint32_t extent = 1u << maxDepth ;
    /// Integer boxes describe adjacency, not physical positions or selected centers.
    struct LeafBox {
      Cell cell ;
      std::array<std::uint32_t,3> lower ;
      unsigned depth ;
      std::uint32_t width() const { return extent >> depth ; }
    } ;
  private:
    static constexpr std::uint32_t absent = std::numeric_limits<std::uint32_t>::max() ;
    struct Split {
      Cell parent ;
      /// Accepted split at each child; absent denotes an unsplit child.
      std::array<std::uint32_t,8> children ;
      /// Optional coordinates; midpoint centers need no record.
      std::uint32_t centers ;
      bool active ;
    } ;
    struct Data {
      std::vector<Split> splits ;
      std::vector<HexCenters> centers ;
    } ;
    std::unique_ptr<Data> data_ ;
    std::uint32_t splitId(Cell cell) const ;
    friend class HexCellMap ;
  public:
    HexPlan() = default ;
    HexPlan(const HexPlan &other) : data_(other.data_ ? std::make_unique<Data>(*other.data_) : nullptr) {}
    HexPlan &operator=(const HexPlan &other) { HexPlan copy(other) ; data_.swap(copy.data_) ; return *this ; }
    HexPlan(HexPlan &&) noexcept = default ;
    HexPlan &operator=(HexPlan &&) noexcept = default ;
    bool empty() const { return !data_ ; }
    bool active(Cell cell) const ;
    bool leaf(Cell cell) const ;
    bool accepted(Cell cell) const { return splitId(cell) != absent ; }
    Cell parent(Cell cell) const ;
    unsigned depth(Cell cell) const ;
    /// Active leaves in depth-first child order. IDs remain stable across coarsening and reload.
    std::vector<Cell> leaves() const ;
    /// Boxes in leaf order, within [0,extent) on each axis.
    std::vector<LeafBox> leafBoxes() const ;
    /// Active leaf containing an integer point in the adjacency boxes.
    Cell findLeaf(const std::array<std::uint32_t,3> &point) const ;
    std::size_t leafCount() const ;
    /// Accepted children, including dormant ones; throws for a cell never split.
    std::array<Cell,8> children(Cell cell) const ;
    HexCenters centers(Cell cell) const ;
    /// Certify the reference subdivision before appending. Existing subdivisions reuse their centers.
    /// The caller checks physical geometry and supplies agreed shared-face centers before this operation.
    /// Reactivation selects eight leaves; it does not restore their deeper selections.
    void refine(Cell cell, const HexCenters &centers = {}, JacobianCheckLimits limits = {}) ;
    /// Collapse eight active leaves. Reject incomplete families without changing the plan.
    void coarsen(Cell cell) ;
    /// Parent of an accepted split, in creation order.
    Cell splitParent(std::size_t index) const ;
    /// Temporary child placements from the accepted centers; no geometry is retained in the plan.
    VertexVolumeSplit construction(Cell cell) const ;
    std::size_t splitCount() const { return data_ ? data_->splits.size() : 0 ; }
    std::size_t centerCount() const { return data_ ? data_->centers.size() : 0 ; }
    /// Live object and vector capacity, excluding allocator metadata and caller-owned source geometry.
    std::size_t storageBytes() const ;
    /// Explicit links and accepted centers only; the empty plan has an empty encoding.
    std::vector<unsigned char> encode() const ;
    static HexPlan decode(const std::vector<unsigned char> &bytes, std::size_t maxSplits = 1000000,
                          JacobianCheckLimits limits = {}) ;
  } ;

  /// Temporary evaluator borrowing an unchanged source map and plan; neither may be destroyed or edited during use.
  /// A descendant composes its accepted child placements with the original map, without refitting corners.
  class HexCellMap final : public GeometryMap {
    const GeometryMap &source_ ;
    const HexPlan &plan_ ;
    HexPlan::Cell cell_ ;
  public:
    HexCellMap(const GeometryMap &source, const HexPlan &plan, HexPlan::Cell cell) ;
    Shape shape() const override { return Shape::hex ; }
    /// An affine corner source and midpoint-only ancestry give a constant Jacobian.
    bool affine() const ;
    MapValue evaluate(const Point &point) const override ;
    Point position(const Point &point) const override ;
  } ;

  /// Choose centers on one hex. Callers reconcile shared face choices before accepting a mesh-wide split.
  HexCenters hexCenters(const GeometryMap &geometry, CenterMethod method, MeasureOptions options = {}) ;
  /// Certify reference children; bound the physical Jacobian where possible, otherwise sample it.
  /// This does not perform neighbor agreement, balancing, or certify global injectivity.
  void refineHex(HexPlan &plan, const GeometryMap &source, HexPlan::Cell cell,
                 const HexCenters &centers, JacobianCheckLimits limits = {}) ;
  /// Choose centers only for new splits; reactivation keeps the accepted construction.
  void refineHex(HexPlan &plan, const GeometryMap &source, HexPlan::Cell cell,
                 CenterMethod method, MeasureOptions options = {}, JacobianCheckLimits limits = {}) ;

}
}
#endif
