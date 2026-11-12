// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.

#ifndef LOCI_FVMADAPT_PARAMETRIC_CELL_H
#define LOCI_FVMADAPT_PARAMETRIC_CELL_H

#include "parametric_mesh.h"
#include "parametric_corner.h"
#include "face_patch.h"
#include "parametric/planar_subdivision.h"
#include <variant>
#include <functional>
#include "spacing_source.h"

namespace Loci {
namespace FVMAdapt {

  /// A planar subdivision extended between two fixed planes along a physical axis.
  struct ExtrudedSubdivision {
    Shape shape ;
    std::vector<Point> corners ;
    int axis ;
    PlanarSubdivision planar ;
  } ;

  struct CellSelection ;

  /// Original corner-cell geometry, with no subdivision or accepted split choices.
  struct ParametricSource {
    CellTopology topology ;
    std::vector<Point> corners ;
    HistoryLimits limits ;
  } ;

  /// Local cell operations used by the distributed rules; contains no MPI state.
  class ParametricCell {
    using Subdivision = std::variant<AffineSubdivision,ExtrudedSubdivision,ParametricCorner,ParametricBoundary> ;
    std::shared_ptr<const ParametricSource> source_ ;
    std::unique_ptr<Subdivision> cell_ ;
    explicit ParametricCell(std::shared_ptr<const ParametricSource> source) : source_(std::move(source)) {}
    static ParametricCell original(const CellTopology &,const std::vector<Point> &,HistoryLimits,const CellSelection &) ;
    friend class ParametricCellData ;
  public:
    ParametricCell(AffineSubdivision cell) : cell_(std::make_unique<Subdivision>(std::move(cell))) {}
    ParametricCell(ExtrudedSubdivision cell) : cell_(std::make_unique<Subdivision>(std::move(cell))) {}
    ParametricCell(ParametricCorner cell) : cell_(std::make_unique<Subdivision>(std::move(cell))) {}
    ParametricCell(ParametricBoundary cell) : cell_(std::make_unique<Subdivision>(std::move(cell))) {}
    ParametricCell(const ParametricCell &) ;
    ParametricCell &operator=(const ParametricCell &) ;
    ParametricCell(ParametricCell &&) noexcept = default ;
    ParametricCell &operator=(ParametricCell &&) noexcept = default ;
    static ParametricCell original(const CellTopology &,const std::vector<Point> &,HistoryLimits) ;
    const ParametricSource *source() const { return source_.get() ; }
    /// Construct a working subdivision before editing an original cell.
    void expand() ;
    std::size_t leafCount() const ;
    std::vector<TreeNodeId> leafIds() const ;
    const HistoryLimits &limits() const ;
    const ParametricCorner *corner() const { return std::get_if<ParametricCorner>(cell_.get()) ; }
    const ParametricBoundary *boundary() const { return std::get_if<ParametricBoundary>(cell_.get()) ; }
    bool hasCornerGeometry() const { return corner() || boundary() ; }
    /// Shape-independent cell ancestry; untouched sources contain only node zero.
    std::uint64_t rootKey() const ;
    CellPath path(TreeNodeId) const ;
    TreeNodeId parent(TreeNodeId) const ;
    NodeRange children(TreeNodeId) const ;
    bool active(TreeNodeId) const ;
    bool leaf(TreeNodeId) const ;
    /// Retained constructions, including dormant subdivisions.
    std::size_t splitCount() const ;
    TreeNodeId splitParent(std::size_t) const ;
    NodeRange splitChildren(std::size_t) const ;
    /// Shared-boundary operations used by corner output and grading.
    const CellTopology &cornerTopology() const ;
    const CornerCell &cornerCell(TreeNodeId) const ;
    /// Construction tables; mapped boundary cells store source-reference coordinates.
    const CornerBoundary &cornerGeometry() const ;
    CornerBoundary &cornerGeometry() ;
    /// Physical geometry, including the retained source map of a boundary-defined child.
    Point cornerPosition(BoundaryId vertex) const ;
    std::shared_ptr<const GeometryMap> cornerEdgeGeometry(BoundaryId edge) const ;
    std::shared_ptr<const GeometryMap> cornerFaceGeometry(BoundaryId face) const ;
    BoundaryId cornerFaceId(EntityId) const ;
    /// Accept shared boundary divisions; return whether the saved construction changed.
    bool mergeFace(EntityId face, const ParametricRecord &record) ;
    const AffineSubdivision *native() const { return std::get_if<AffineSubdivision>(cell_.get()) ; }
    const ExtrudedSubdivision *extruded() const { return std::get_if<ExtrudedSubdivision>(cell_.get()) ; }
    /// Standard-element history; unavailable for untouched or boundary-defined cells.
    const RefinementHistory &history() const ;
    void refine(TreeNodeId leaf, const SplitOperation &request, CenterMethod method) ;
    void refineIsotropic(TreeNodeId leaf, CenterMethod method, int axis) ;
    /// Restrict native hex/prism balancing to X/Y/Z bits, preserving the extrusion axis.
    /// Zero or no eligible direction uses the full balancing split; other shapes stay isotropic.
    void refineForBalance(TreeNodeId leaf, CenterMethod method, int axis, unsigned directions) ;
    void coarsen(TreeNodeId parent) ;
    /// Undo a merge refused by balancing by selecting its previous accepted children.
    bool restoreSplit(TreeNodeId parent, const RefinementHistory &previous) ;
    bool restoreSplit(TreeNodeId parent, const ParametricCell &previous) ;
    void validate(int axis) const ;
  } ;

  /// Resolve edits without requiring a standard-element map for every child.
  ParametricRecord encodeParametricRequests(const ParametricCell &, const std::vector<ParametricLeafRequest> &) ;
  std::vector<ParametricLeafRequest> decodeParametricRequests(const ParametricRecord &, const ParametricCell &) ;
  /// Tags: 0 keep, 1 refine, 2 coarsen a complete immediate family.
  ParametricRecord encodeParametricTags(const ParametricCell &, const std::vector<char> &) ;

  CenterMethod parametricCenterMethod(const std::string &name) ;
  FaceBalance parametricFaceBalance(const std::string &name) ;
  /// Active divided faces on both sides of an original face; unused saved cuts do not count.
  std::vector<FaceBalanceRequest> parametricDividedFaces(EntityId, const CellTopology &, const ParametricCell &,
      const CellTopology &, const ParametricCell &, bool directions = false) ;
  /// Combine local face subdivisions with votes from original-cell neighbors.
  /// Optional output collects native leaf directions even when the face policy is none.
  std::vector<int> parametricFaceBalanceRequests(const CellTopology &, const ParametricCell &,
      const std::vector<FaceBalanceRequest> &, FaceBalance, std::map<int,unsigned> *directions = nullptr) ;
  ParametricCell decodeParametricCell(const ParametricRecord &record, HistoryLimits limits = {}) ;
  ParametricRecord encodeParametricCell(const ParametricCell &cell) ;
  /// Build or load a cell, verifying its original geometry and requested split family.
  ParametricCell makeParametricCell(const CellTopology &topology, const std::vector<Point> &corners,
      const ParametricRecord &record, HistoryLimits limits, const std::string &family, int axis,
      const ParametricEdgeControls &controls = {}) ;
  ParametricRecord initialParametricCell(const CellTopology &topology, const std::vector<Point> &corners,
      const ParametricRecord &record, HistoryLimits limits, const std::string &family, int axis,
      const ParametricEdgeControls &controls = {}) ;

  /// Output key within an original entity. local >= 0 indexes accepted corner or extruded points.
  /// Local output indices are rebuilt per mesh; persistent geometry remains in the history.
  struct ParametricNodeKey {
    VertexKey reference ;
    int local = -1 ;
    bool operator<(const ParametricNodeKey &other) const ;
    bool operator==(const ParametricNodeKey &other) const ;
  } ;
  struct ParametricOutputFace {
    std::vector<ParametricNodeKey> vertices ;
    int left, right ;
    FacePatch patch{} ;
  } ;

  std::vector<int> parametricLocalRequests(const CellTopology &, const ParametricCell &, unsigned, NeighborGrading) ;
  ParametricTrace parametricFaceRequests(EntityId face, const CellTopology &, const ParametricCell &,
      const CellTopology &, const ParametricCell &, unsigned, NeighborGrading) ;
  ParametricTrace parametricEdgeTrace(const CellTopology &, const ParametricCell &) ;
  std::vector<int> parametricVertexDepths(const CellTopology &, const ParametricCell &) ;
  std::vector<int> parametricVertexRequests(const CellTopology &, const ParametricCell &, const std::map<EntityId,int> &, unsigned) ;
  /// Recursive selection from own physical vertices and the shortest own edge; an empty selector makes no requests.
  struct GeometricRefinement {
    std::function<bool(const std::vector<Point> &,double)> select ;
    double minEdgeLength = 0 ; ///< Physical edge cutoff; isotropic splits require shortestEdge >= 2*minEdgeLength.
    static GeometricRefinement region(const std::function<bool(const Point &)> &,double minEdgeLength) ;
    /// Borrow a validated source list, which must outlive this selection.
    static GeometricRefinement sources(const std::vector<SpacingSource> &,double minEdgeLength) ;
    bool refine(ParametricCell &,CenterMethod,int axis,bool automatic = false,double aspectRatio = 2) const ;
  } ;
  /// Refine selected leaves and descendants to the edge cutoff; neighbor grading is separate.
  bool refineParametricRegion(ParametricCell &, const std::function<bool(const Point &)> &, double minEdgeLength,
      CenterMethod, int axis, bool automatic = false, double aspectRatio = 2) ;
  /// Refine descendants using source spacing; cutoff may be zero to use only the source limits.
  bool refineParametricSources(ParametricCell &, const std::vector<SpacingSource> &, double minEdgeLength,
      CenterMethod, int axis, bool automatic = false, double aspectRatio = 2) ;
  /// Own vertices in leaf order; excludes hanging vertices introduced by neighbors.
  std::vector<std::vector<ParametricNodeKey>> parametricLeafVertices(const CellTopology &, const ParametricCell &) ;
  /// Sorted interior nodes of an original edge from its accepted intervals.
  std::vector<ParametricNodeKey> parametricEdgeNodes(EntityId, const ParametricTrace &) ;
  std::vector<ParametricNodeKey> parametricCellNodes(const CellTopology &, const ParametricCell &) ;
  /// Original-face nodes, optionally with physical positions in the same order.
  std::vector<ParametricNodeKey> parametricFaceNodes(EntityId face, const CellTopology &, const ParametricCell &,
      const CellTopology *right = nullptr, const ParametricCell *rightCell = nullptr, std::vector<Point> *positions = nullptr) ;
  /// Interior polygons, optionally with retained-surface identity and physical area.
  std::vector<ParametricOutputFace> parametricCellFaces(const CellTopology &, const ParametricCell &,
      const std::vector<ParametricNodeKey> &boundaryNodes, bool patches = false) ;
  std::vector<ParametricOutputFace> parametricFaceFaces(EntityId face, const CellTopology &, const ParametricCell &,
      const std::vector<ParametricNodeKey> &edgeNodes, const CellTopology *right = nullptr, const ParametricCell *rightCell = nullptr,
      bool patches = false) ;
  /// Physical positions in the supplied original-cell-owned output-key order.
  std::vector<Point> parametricCellPositions(const ParametricCell &, const std::vector<ParametricNodeKey> &) ;
  /// Accepted source-face divisions; no geometry maps. Empty for other split families.
  ParametricRecord parametricFacePlan(EntityId face, const ParametricCell &) ;
  /// Positions for original-face-owned nodes, in parametricFaceNodes order.
  std::vector<Point> parametricFacePositions(EntityId face, const CellTopology &, const ParametricCell &,
      const CellTopology *right = nullptr, const ParametricCell *rightCell = nullptr) ;
  /// Immediate old leaf indices contributing to each new leaf; complete sibling merges are supported.
  std::vector<std::vector<int>> parametricCellDonors(const RefinementHistory &oldHistory, const RefinementHistory &newHistory) ;
  std::vector<std::vector<int>> parametricCellDonors(const ParametricCell &oldCell, const ParametricCell &newCell) ;

}
}
#endif
