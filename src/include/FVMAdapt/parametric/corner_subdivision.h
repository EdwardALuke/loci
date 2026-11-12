// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_CORNER_SUBDIVISION_H
#define LOCI_FVMADAPT_CORNER_SUBDIVISION_H

#include "corner_boundary.h"

namespace Loci {
namespace FVMAdapt {

  /// One retained cell and its selected subdivision. The root has no parent.
  struct CornerNode {
    TreeNodeId parent = noTreeNode ;
    unsigned ordinal = 0 ;
    std::size_t split = noTreeNode ;
    bool refined = false ;
    CornerCell cell ;
    Point interior ;
  } ;

  /// Accepted centers and child range; geometry is shared through the boundary tables.
  struct CornerConstruction {
    TreeNodeId parent ;
    NodeRange children ;
    CenterMethod centerMethod ;
    BoundaryId centerVertex ;
    std::vector<std::size_t> faceSplits ;
  } ;

  /// Isotropic subdivision of pyramids and general cells bounded by polygon faces.
  /// Boundary tables use source-reference coordinates when mapped, physical coordinates otherwise.
  /// Each construction face must be visible from the supplied interior point.
  /// One child per vertex; accepted faces and cell identities survive coarsening and reload.
  /// Mesh ownership, neighbor communication and output numbering belong to the module.
  class CornerSubdivision : public CornerBoundary {
    std::uint64_t root_ ;
    HistoryLimits limits_ ;
    MapId source_ = noBoundary ;
    std::vector<CornerNode> nodes_ ;
    std::vector<CornerConstruction> splits_ ;
    CornerSubdivision(std::uint64_t root, HistoryLimits limits) ;
    BoundaryRecordLimits savedBoundaryLimits() const ;
    void validate() const ;
    void validateSplit(const CornerConstruction &) const ;
    GeometryMeasure measureCellGeometry(const CornerCell &, const Point &) const ;
    std::shared_ptr<const GeometryMap> project(std::shared_ptr<const GeometryMap>) const ;
  public:
    CornerSubdivision(std::uint64_t root, CornerBoundary boundary, CornerCell cell, Point interior, HistoryLimits limits = {}) ;
    /// Split a standard source volume in its reference coordinates; retain the source map once.
    /// Inherited boundary tables hold construction coordinates. Use the accessors below for physical geometry.
    /// New mapped splits currently require CenterMethod::reference.
    CornerSubdivision(std::uint64_t root, std::shared_ptr<const GeometryMap> source, HistoryLimits limits = {}) ;
    /// Use caller-oriented reference faces of a standard source volume.
    /// A null source leaves the supplied boundary coordinates physical.
    CornerSubdivision(std::uint64_t root, std::shared_ptr<const GeometryMap> source,
        CornerBoundary boundary, CornerCell cell, Point interior, HistoryLimits limits = {}) ;
    /// Null for a cell supplied directly through its physical boundary.
    std::shared_ptr<const GeometryMap> sourceGeometry() const ;
    Point position(BoundaryId vertex) const ;
    std::shared_ptr<const GeometryMap> edgeGeometry(BoundaryId edge) const ;
    std::shared_ptr<const GeometryMap> faceGeometry(BoundaryId face) const ;
    /// Physical volume pieces; mapped descendants all compose with the original source.
    std::vector<std::shared_ptr<const GeometryMap>> cellPieces(TreeNodeId id) const ;
    std::uint64_t rootKey() const { return root_ ; }
    const HistoryLimits &limits() const { return limits_ ; }
    const std::vector<CornerNode> &nodes() const { return nodes_ ; }
    const std::vector<CornerConstruction> &splits() const { return splits_ ; }
    /// Refinement can invalidate node references.
    const CornerNode &node(TreeNodeId id) const { return nodes_.at(id) ; }
    bool active(TreeNodeId id) const ;
    bool leaf(TreeNodeId id) const { return active(id) && !node(id).refined ; }
    unsigned depth(TreeNodeId id) const ;
    CellPath path(TreeNodeId id) const ;
    /// Resolve an active or dormant path. Child order is the parent vertex order.
    TreeNodeId nodeId(const CellPath &) const ;
    /// Retained children, including dormant ones; empty for a cell never split.
    NodeRange children(TreeNodeId id) const ;
    std::vector<TreeNodeId> leafIds() const ;
    /// New splits reuse accepted shared faces. Reactivation uses the entire accepted construction.
    /// Failed edits preserve both the tree and its boundary geometry.
    NodeRange refine(TreeNodeId leaf, CenterMethod method) ;
    /// Reference centers for mapped sources, weighted centers for physical boundary cells.
    NodeRange refine(TreeNodeId leaf) { return refine(leaf,source_==noBoundary ? CenterMethod::weighted:CenterMethod::reference) ; }
    /// Collapse a complete active leaf family; neighboring cuts remain in the shared boundary tables.
    void coarsen(TreeNodeId parent) ;
    /// Physical volume and centroid; sampled Jacobian checks are not a global validity certificate.
    GeometryMeasure measure(TreeNodeId id) const { const auto &n=node(id) ; return measureCellGeometry(n.cell,n.interior) ; }
    /// Save explicit cell links, accepted constructions and shared geometry; IDs keep their order.
    std::vector<unsigned char> encode() const ;
    static CornerSubdivision decode(const std::vector<unsigned char> &, HistoryLimits limits = {}) ;
  } ;
}
}
#endif
