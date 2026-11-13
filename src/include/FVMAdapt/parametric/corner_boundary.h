// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_CORNER_BOUNDARY_H
#define LOCI_FVMADAPT_CORNER_BOUNDARY_H

#include "boundary_record.h"
#include "hex_plan.h"
#include "mapped_volume_subdivision.h"
#include "split_quality.h"
#include "refinement_record.h"
#include <algorithm>
#include <cmath>
#include <map>

namespace Loci {
namespace FVMAdapt {
  /// Construction and assessment limits for corner subdivision.
  struct CornerBoundaryOptions {
    /// Retained history nodes per source cell.
    std::size_t maxNodes = 1000000 ;
    /// Combined boundary point, edge and face count.
    std::size_t maxEntities = 10000000 ;
    /// Maximum encoded snapshot size, including shared geometry.
    std::size_t maxSnapshotBytes = 64*1024*1024 ;
    JacobianCheckLimits jacobianLimits ;
    SplitQualityOptions assessment ;
  } ;
  struct CornerFaceUse {
    BoundaryId surface ;
    /// Maps the cell's outward vertex order to the stored face order.
    unsigned orientation ;
  } ;
  /// One cell's vertices and outward uses of shared faces.
  struct CornerCell { std::vector<BoundaryId> vertices ; std::vector<CornerFaceUse> faces ; } ;

  /// Accepted corner split of a boundary-defined cell; children follow parent vertex order.
  struct CornerSplit {
    CenterMethod centerMethod = CenterMethod::weighted ;
    BoundaryId centerVertex = noBoundary ;
    /// FaceRecord::splits indices in parent-cell face order.
    std::vector<std::size_t> faceSplits ;
    std::vector<CornerCell> children ;
    /// Interior points used to measure each child's retained surface geometry.
    std::vector<Point> interiors ;
  } ;

  /// Local shared geometry and corner construction; contains no mesh ownership or MPI.
  struct CornerBoundary {
    using VertexId = BoundaryId ;
    using EdgeId = BoundaryId ;
    using SurfaceId = BoundaryId ;
    static constexpr auto absent = noBoundary ;
    CornerBoundaryOptions options ;
    BoundaryRecords boundaries ;
    static void require(bool ok, const char *message) { if(!ok) throw std::invalid_argument(message) ; }
    static auto endpoints(VertexId a, VertexId b) { return std::make_pair(std::min(a,b),std::max(a,b)) ; }
    BoundaryRecordLimits boundaryLimits() const ;
    RecordLimits recordLimits() const ;
    /// Shared geometry is appended; edits to existing records are journaled.
    struct Transaction {
      CornerBoundary &state ;
      std::size_t points, edges, faces, maps ;
      std::vector<std::pair<EdgeId,EdgeRecord>> oldEdges ;
      std::vector<std::pair<SurfaceId,FaceRecord>> oldFaces ;
      bool committed = false ;
      explicit Transaction(CornerBoundary &s) : state(s), points(s.boundaries.points.size()), edges(s.boundaries.edges.size()),
        faces(s.boundaries.faces.size()), maps(s.boundaries.geometry.size()) {}
      void saveEdge(EdgeId id) { if(id < edges) oldEdges.emplace_back(id,state.boundaries.edges.at(id)) ; }
      void saveSurface(SurfaceId id) { if(id < faces) oldFaces.emplace_back(id,state.boundaries.faces.at(id)) ; }
      ~Transaction() {
        if(committed) return ;
        auto &b = state.boundaries ;
        for(auto i = oldEdges.rbegin(); i != oldEdges.rend(); ++i) b.edges[i->first] = std::move(i->second) ;
        for(auto i = oldFaces.rbegin(); i != oldFaces.rend(); ++i) b.faces[i->first] = std::move(i->second) ;
        b.points.resize(points) ; b.edges.resize(edges) ; b.faces.resize(faces) ; b.geometry.truncate(maps) ;
      }
    } ;
    void checkEntity() const ;
    VertexId addPoint(const Point &point) ;
    EdgeId addEdge(VertexId a, VertexId b, std::shared_ptr<const GeometryMap> geometry, EdgeId parent = absent) ;
    SurfaceId addSurface(std::vector<VertexId> vertices, std::vector<EdgeId> edges,
                         std::shared_ptr<const GeometryMap> geometry, int boundary = 0,
                         SurfaceId parent = absent, std::size_t parentSplit = absent) ;
    VertexId splitEdge(EdgeId id, Transaction &transaction) ;
    EdgeId edgeSegment(EdgeId parent, VertexId endpoint) const ;
    /// active[surface] selects a retained split, overriding method; absent or missing leaves it unspecified.
    std::size_t splitSurface(SurfaceId id, CenterMethod method, const std::vector<std::size_t> &active,
      Transaction &transaction, const FaceSplitRecord *accepted = nullptr) ;
    /// Include only cuts required by emitted neighboring edges.
    void edgeVertices(EdgeId id, VertexId start, const std::vector<bool> &split, std::vector<VertexId> &vertices) const ;
    std::map<std::pair<VertexId,VertexId>,EdgeId> cellEdges(const CornerCell &cell) const ;

    /// Check a connected, closed polygon boundary with consistent face orientations.
    void validateCell(const CornerCell &cell) const ;
    /// Cone pieces joining the oriented surfaces to an interior point; one logical cell.
    std::vector<std::shared_ptr<const GeometryMap>> cellPieces(const CornerCell &cell, const Point &interior) const ;
    /// Require positive sampled cone Jacobians and converged volume integration.
    /// This is not a global injectivity certificate.
    GeometryMeasure measureCell(const CornerCell &cell, const Point &interior) const ;
    /// Split arbitrary vertex valence using the shared face constructions.
    /// Reference centers use interior; weighted uses area-vector weights; integrated uses the volume centroid.
    /// Failure preserves all boundary records. The caller retains the accepted split for reuse after coarsening.
    CornerSplit splitCell(const CornerCell &cell, Point interior, CenterMethod method,
                         const std::vector<std::size_t> &active = {}) ;

    /// Check cell incidence and sampled agreement with retained face geometry.
    void validateCells(const MappedVolumeSubdivision &subdivision, const std::vector<CornerCell> &cells) const ;

    /// Split a cell using the selected shared face constructions. Failure preserves its state.
    std::vector<TreeNodeId> splitCell(MappedVolumeSubdivision &subdivision, std::vector<CornerCell> &cells,
        TreeNodeId leaf, CenterMethod method, const std::vector<std::size_t> &active,
        const SplitQualityLimits &limits = {}, SplitQuality *report = nullptr) ;
    /// Rebuild an accepted hex split. Its shared face choices must already be present.
    void restoreHexCell(MappedVolumeSubdivision &subdivision, std::vector<CornerCell> &cells,
                        HexPlan::Cell leaf, const HexPlan &plan) ;
  private:
    friend class CornerSubdivision ;
    /// Selected split indices follow the cell's face order; absent leaves a face unspecified.
    CornerSplit constructCell(const CornerCell &cell, Point interior, CenterMethod method,
                             const std::vector<std::size_t> &selected, Transaction &transaction) ;
    std::vector<TreeNodeId> constructCell(MappedVolumeSubdivision &subdivision, std::vector<CornerCell> &cells,
        TreeNodeId leaf, CenterMethod method, const std::vector<std::size_t> &active,
        const SplitQualityLimits &limits, SplitQuality *report, const VertexVolumeSplit *accepted) ;
  } ;
}
}
#endif
