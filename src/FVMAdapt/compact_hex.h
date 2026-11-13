// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_COMPACT_HEX_H
#define LOCI_FVMADAPT_COMPACT_HEX_H

#include <FVMAdapt/parametric_corner.h>
#include <FVMAdapt/parametric/face_plan.h>
#include <FVMAdapt/parametric/hex_plan.h>

namespace Loci {
namespace FVMAdapt {
  struct GeometricRefinement ;
  /// Module-owned hex source, compact cell plan and accepted shared-face choices.
  /// Working geometry is constructed locally and never retained here.
  struct CompactHexCell {
    struct InteriorFace {
      std::uint32_t split ;
      unsigned face ;
      FacePlan plan ;
    } ;
    CellTopology topology ;
    std::vector<Point> corners ;
    /// Curved source geometry; null retains the inexpensive corner-only source.
    std::shared_ptr<const CubicHexMap> cubic ;
    HistoryLimits limits ;
    HexPlan plan ;
    std::array<FacePlan,6> faces ;
    /// Only divided interior faces need records; ordinals follow child-face order.
    std::vector<InteriorFace> interior ;

    explicit CompactHexCell(const ParametricCorner &cell) ;
    /// Unsplit hex from validated source topology and coordinates.
    CompactHexCell(CellTopology topology, std::vector<Point> corners, HistoryLimits limits, std::shared_ptr<const CubicHexMap> cubic = {}) ;
    /// Temporary source evaluator; curved sources retain their shared edge maps.
    std::shared_ptr<const GeometryMap> sourceMap() const ;
    /// Original face geometry in canonical face coordinates.
    std::shared_ptr<const GeometryMap> faceMap(unsigned face) const ;
    /// Merge borrowed face plans, then apply decoded requests to a private candidate; return whether it changed.
    /// Geometric selection follows explicit requests; failure may leave the candidate partially edited.
    bool edit(const std::vector<ParametricLeafRequest> &requests,int levels,CenterMethod method,const std::array<const FacePlan *,6> &incoming,
        const GeometricRefinement *selection) ;
    /// Build temporary geometry within the working-history allocation limits.
    ParametricCorner expand() const ;
    ParametricRecord record() const ;
    static bool saved(const ParametricRecord &record) ;
    /// Check bounded plans, shared choices and physical geometry without expanding the cell.
    static CompactHexCell decode(const ParametricRecord &record, HistoryLimits limits) ;
    /// Grading queries use only the split topology; no physical geometry is built.
    std::vector<int> localRequests(unsigned difference, NeighborGrading grading) const ;
    ParametricTrace edgeTrace() const ;
    std::vector<int> vertexDepths() const ;
    std::vector<int> vertexRequests(const std::map<EntityId,int> &depths, unsigned difference) const ;
    /// Active leaf coverage in canonical face coordinates; boxes describe adjacency only.
    struct FaceLeaf {
      std::array<std::uint32_t,2> lower ;
      std::uint32_t width ;
      int leaf ; ///< Current output ordinal, not a retained plan-node ID.
    } ;
    std::vector<FaceLeaf> faceLeaves(EntityId face) const ;
    /// Cell-interior nodes, optionally with positions in the same order, from accepted connectivity.
    std::vector<ParametricNodeKey> cellNodes(std::vector<Point> *positions = nullptr) const ;
    std::vector<Point> cellPositions(const std::vector<ParametricNodeKey> &nodes) const ;
    std::vector<ParametricOutputFace> cellFaces(const std::vector<ParametricNodeKey> &required) const ;
    /// Original-face nodes, optionally with positions in the same order, without expanding cell geometry.
    std::vector<ParametricNodeKey> faceNodes(EntityId face,const CompactHexCell *right = nullptr,
        std::vector<Point> *positions = nullptr) const ;
    std::vector<Point> facePositions(EntityId face,const CompactHexCell *right = nullptr) const ;
    std::vector<ParametricOutputFace> faceFaces(EntityId face,const std::vector<ParametricNodeKey> &required,
        const CompactHexCell *right = nullptr) const ;
  private:
    void validate() const ;
  } ;
  /// Resolve persistent requests against accepted hex constructions without expanding geometry.
  std::vector<ParametricLeafRequest> hexRequests(const ParametricRecord &record,const HexPlan &plan,EntityId root,const HistoryLimits &limits) ;
  /// Cell tags to persistent requests using active leaves and accepted reference constructions.
  ParametricRecord hexTagRequests(const HexPlan &plan,EntityId root,const HistoryLimits &limits,const std::vector<char> &tags) ;
  ParametricTrace hexFaceRequests(EntityId face, const CellTopology &leftTopology, const HexPlan &leftPlan,
      const CellTopology &rightTopology, const HexPlan &rightPlan, unsigned difference, NeighborGrading grading) ;
  std::vector<FaceBalanceRequest> hexDividedFaces(const CellTopology &topology, const HexPlan &plan) ;
  std::vector<FaceBalanceRequest> hexDividedFaces(EntityId face, const CellTopology &leftTopology, const HexPlan &leftPlan,
      const CellTopology &rightTopology, const HexPlan &rightPlan) ;
  /// Old leaf indices for each new leaf; compare accepted centers along shared ancestry.
  std::vector<std::vector<int>> hexCellDonors(const HexPlan &oldPlan, const HexPlan &newPlan) ;
}
}
#endif
