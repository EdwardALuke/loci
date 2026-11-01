// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_PARAMETRIC_CELL_DATA_H
#define LOCI_FVMADAPT_PARAMETRIC_CELL_DATA_H

#include <Loci.h>
#include <FVMAdapt/defines.h>
#include <FVMAdapt/parametric_cell.h>
#include <memory>
#include <optional>

namespace Loci {
namespace FVMAdapt {
  inline std::ostream &operator<<(std::ostream &out,const FaceBalanceRequest &r) { return out<<r.root<<' '<<r.leaf<<' '<<r.face<<' '<<r.directions ; }
  inline std::istream &operator>>(std::istream &in,FaceBalanceRequest &r) { return in>>r.root>>r.leaf>>r.face>>r.directions ; }
  class HexPlan ;

  /// Node keys in original-entity ownership order, for Loci stores.
  Loci::FineFaces encodeParametricNodeKeys(const std::vector<ParametricNodeKey> &) ;
  std::vector<ParametricNodeKey> decodeParametricNodeKeys(const Loci::FineFaces &) ;

  /// Cell-local controls; neighbor policies and stopping conditions belong to the rules.
  struct CellRefinementOptions {
    CenterMethod center = CenterMethod::reference ;
    int levels = 1, extrusionAxis = -1 ;
    bool automatic = false ;
    double aspectRatio = 2, minEdgeLength = 0 ;
  } ;
  /// Face records borrowed for the duration of a cell edit.
  using IncomingFacePlans = std::vector<std::pair<EntityId,const ParametricRecord *>> ;

  /// Immutable cell state shared by rules; hexes retain plans, not expanded geometry.
  class ParametricCellData {
    struct Data ;
    std::shared_ptr<const Data> data_ ;
    friend class ParametricCellDataConverter ;
  public:
    ParametricCellData() = default ;
    /// Validate received state with the current limits; workingCell() reports a failed receive.
    ParametricCellData(ParametricRecord record,int maxDepth,int maxNodes) ;
    /// Retain accepted state; hex working geometry is reduced to compact plans.
    explicit ParametricCellData(ParametricCell cell) ;
    /// Original corner cells retain source geometry and an empty saved history.
    static ParametricCellData initial(const CellTopology &,const std::vector<Point> &,const ParametricRecord &,
        int maxDepth,int maxNodes,const std::string &family,int axis) ;
    /// Read-only working geometry; compact cells expand for the lifetime of this handle.
    std::shared_ptr<const ParametricCell> workingCell() const ;
    /// An independent working cell for an edit.
    ParametricCell workingCopy() const ;
    /// Resolve persistent requests against this cell's accepted history.
    std::vector<ParametricLeafRequest> decodeRequests(const ParametricRecord &) const ;
    /// Merge shared faces and apply requests atomically; failures preserve this state.
    ParametricCellData applyRequests(const ParametricRecord &, const CellRefinementOptions &,
        const IncomingFacePlans & = {}, const GeometricRefinement & = {}) const ;
    /// Refine required leaf indices, restoring previous splits when coarsening was refused.
    /// changed reports a published edit; neighbor votes are supplied by the rules.
    ParametricCellData balance(const std::set<int> &required, const std::map<int,unsigned> &directions,
        const ParametricCellData &previous, CenterMethod, int axis, const IncomingFacePlans &, bool &changed) const ;
    /// Apply isotropic hex splits and coarsening on an independent plan; nullopt for other representations or split families.
    /// Requests name this state's cells; incoming face plans are borrowed for this call and merged before cell edits.
    /// Geometric selection then refines descendants to its cutoff. changed is true only after publishing a real edit.
    std::optional<ParametricCellData> editPlan(const std::vector<ParametricLeafRequest> &,int refinementLevels = 1,
        CenterMethod method = CenterMethod::reference,const IncomingFacePlans &incoming = {},
        bool *changed = nullptr,const GeometricRefinement *selection = nullptr) const ;
    /// Tags 0/1/2 keep, refine or request complete-family coarsening in active leaf order.
    ParametricRecord tagRequests(const std::vector<char> &) const ;
    std::size_t leafCount() const ;
    std::vector<int> leafDepths() const ;
    /// Longest leaf edge chord in 3D and XY; ignores neighboring edge subdivisions.
    std::vector<std::array<double,2>> edgeLengths() const ;
    ParametricRecord facePlan(EntityId face) const ;
    /// False when incoming divisions are empty or identical to the retained face plan.
    /// Other records still need merging and validation against the shared boundary.
    bool needsFaceMerge(EntityId face, const ParametricRecord &) const ;
    /// Read accepted topology directly for hexes; other split families use their existing queries.
    std::vector<int> localRequests(const CellTopology &, unsigned difference, NeighborGrading, int axis) const ;
    ParametricTrace edgeTrace(const CellTopology &) const ;
    std::vector<int> vertexDepths(const CellTopology &) const ;
    std::vector<int> vertexRequests(const CellTopology &, const std::map<EntityId,int> &, unsigned difference) const ;
    ParametricTrace faceRequests(EntityId, const CellTopology &, const ParametricCellData &, const CellTopology &,
        unsigned difference, NeighborGrading) const ;
    /// Active divided-face queries; compact hexes read split topology without physical geometry.
    std::vector<FaceBalanceRequest> dividedFaces(EntityId,const CellTopology &,const ParametricCellData &,const CellTopology &,bool directions = false) const ;
    std::vector<int> faceBalanceRequests(const CellTopology &,const std::vector<FaceBalanceRequest> &,FaceBalance,
        std::map<int,unsigned> *directions = nullptr) const ;
    std::vector<std::vector<int>> donorsFrom(const ParametricCellData &previous) const ;
    /// Cell-interior output; compact hexes read accepted plans without expanding geometry.
    std::vector<ParametricNodeKey> cellNodes(const CellTopology &) const ;
    std::vector<Point> cellPositions(const std::vector<ParametricNodeKey> &) const ;
    std::vector<ParametricOutputFace> cellFaces(const CellTopology &,const std::vector<ParametricNodeKey> &) const ;
    /// Original-face output; compact hex pairs read accepted plans directly.
    std::vector<ParametricNodeKey> faceNodes(EntityId,const CellTopology &,const ParametricCellData *right = nullptr,
        const CellTopology *rightTopology = nullptr) const ;
    std::vector<Point> facePositions(EntityId,const CellTopology &,const ParametricCellData *right = nullptr,
        const CellTopology *rightTopology = nullptr) const ;
    std::vector<ParametricOutputFace> faceFaces(EntityId,const CellTopology &,const std::vector<ParametricNodeKey> &,
        const ParametricCellData *right = nullptr,const CellTopology *rightTopology = nullptr) const ;
    /// Compact plan, or an empty plan for an original hex; null for other representations.
    const HexPlan *hexPlan() const ;
    bool compactHex() const ;
    /// Accepted history bytes; empty for an untouched corner cell.
    ParametricRecord record() const ;
    int maxDepth() const ;
    int maxNodes() const ;
  } ;

  class ParametricCellDataConverter {
    ParametricCellData &value_ ;
    ParametricRecord record_ ;
    void encode() ;
  public:
    explicit ParametricCellDataConverter(ParametricCellData &value) ;
    int getSize() const ;
    void getState(unsigned char *buffer,int &size) const ;
    void setState(const unsigned char *buffer,int size) ;
  } ;

  std::ostream &operator<<(std::ostream &,const ParametricCellData &) ;
  std::istream &operator>>(std::istream &,ParametricCellData &) ;
}

  template<> struct data_schema_traits<FVMAdapt::FaceBalanceRequest> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct=CompoundFactory(FVMAdapt::FaceBalanceRequest()) ;
      LOCI_INSERT_TYPE(ct,FVMAdapt::FaceBalanceRequest,root) ;
      LOCI_INSERT_TYPE(ct,FVMAdapt::FaceBalanceRequest,leaf) ;
      LOCI_INSERT_TYPE(ct,FVMAdapt::FaceBalanceRequest,face) ;
      LOCI_INSERT_TYPE(ct,FVMAdapt::FaceBalanceRequest,directions) ;
      return DatatypeP(ct) ;
    }
  } ;
  template<> struct data_schema_traits<FVMAdapt::ParametricCellData> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef unsigned char Converter_Base_Type ;
    typedef FVMAdapt::ParametricCellDataConverter Converter_Type ;
  } ;
}
#endif
