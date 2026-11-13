// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_PARAMETRIC_CELL_DATA_H
#define LOCI_FVMADAPT_PARAMETRIC_CELL_DATA_H

#include <Loci.h>
#include <FVMAdapt/defines.h>
#include <FVMAdapt/parametric_cell.h>
#include <FVMAdapt/quality.h>
#include <functional>
#include <memory>
#include <optional>

namespace Loci {
namespace FVMAdapt {
  inline std::ostream &operator<<(std::ostream &out,const FaceBalanceRequest &r) { return out<<r.root<<' '<<r.leaf<<' '<<r.face<<' '<<r.directions ; }
  inline std::istream &operator>>(std::istream &in,FaceBalanceRequest &r) { return in>>r.root>>r.leaf>>r.face>>r.directions ; }
  class HexPlan ;

  /// Node keys in original-entity ownership order, for Loci stores.
  /// Each row: entity kind, original ID, three numerator/denominator pairs, retained local ID.
  Loci::FineFaces encodeParametricNodeKeys(const std::vector<ParametricNodeKey> &keys) ;
  std::vector<ParametricNodeKey> decodeParametricNodeKeys(const Loci::FineFaces &record) ;

  /// Output numbers from sorted owner records, borrowed unchanged for this object's lifetime.
  class ParametricNodeNumbers {
    struct Owner {
      const Loci::FineFaces *record ;
      int offset ; ///< Output number of the first owned row.
      bool local ; ///< Every row uses a retained local ID rather than reference coordinates.
    } ;
    std::vector<Owner> owners_ ;
  public:
    void add(const Loci::FineFaces &record,int offset) ;
    int at(const ParametricNodeKey &key) const ;
  } ;

  /// Cell-local controls; neighbor policies and stopping conditions belong to the rules.
  struct CellRefinementOptions {
    CenterMethod center = CenterMethod::reference ;
    int levels = 1, extrusionAxis = -1 ;
    bool automatic = false ;
    double aspectRatio = 2, minEdgeLength = 0 ;
  } ;
  /// Parsed controls shared by every cell in one adaptation query.
  struct AdaptationControls {
    CellRefinementOptions refinement ;
    NeighborGrading grading = NeighborGrading::vertices ;
    FaceBalance faceBalance = FaceBalance::none ;
    int levelDifference = 1, maxRounds = 32 ;
    bool directionalBalance = false ;
    MeshQualityOptions quality ;
  } ;
  /// Face records borrowed for the duration of a cell edit.
  struct IncomingFacePlan {
    EntityId face ;
    std::reference_wrapper<const ParametricRecord> record ;
  } ;
  using IncomingFacePlans = std::vector<IncomingFacePlan> ;

  /// Immutable cell state shared by rules; vertex-based hexes retain compact plans.
  class ParametricCellData {
    struct Data ;
    std::shared_ptr<const Data> data_ ;
    friend class ParametricCellDataConverter ;
  public:
    // Accepted state and temporary working geometry.
    ParametricCellData() = default ;
    /// Validate received state with the current limits; workingCell() reports a failed receive.
    ParametricCellData(ParametricRecord record,int maxDepth,int maxNodes) ;
    /// Retain accepted state; hex working geometry is reduced to compact plans.
    explicit ParametricCellData(ParametricCell cell) ;
    /// Original corner cells retain source geometry and an empty saved history.
    static ParametricCellData initial(const CellTopology &topology,const std::vector<Point> &corners,const ParametricRecord &record,
        int maxDepth,int maxNodes,const std::string &family,int axis,const ParametricEdgeControls &controls = {}) ;
    /// Accepted source curves in canonical edge order; straight edges have empty rows.
    ParametricEdgeControls sourceEdgeControls(const CellTopology &topology) const ;
    /// Read-only working geometry; each compact-cell call allocates a temporary expansion.
    /// Other cells share accepted storage. The returned handle owns either lifetime.
    std::shared_ptr<const ParametricCell> workingCell() const ;
    /// An independent working cell for an edit.
    ParametricCell workingCopy() const ;

    // Candidate edits preserve the accepted state on failure.
    /// Resolve persistent requests against this cell's accepted history.
    std::vector<ParametricLeafRequest> decodeRequests(const ParametricRecord &record) const ;
    /// Merge shared faces and apply requests atomically; failures preserve this state.
    ParametricCellData applyRequests(const ParametricRecord &record, const CellRefinementOptions &options,
        const IncomingFacePlans &incoming = {}, const GeometricRefinement &selection = {}) const ;
    /// Refine required leaf indices, restoring previous splits when coarsening was refused.
    /// changed reports a published edit; neighbor votes are supplied by the rules.
    ParametricCellData balance(const std::set<int> &required, const std::map<int,unsigned> &directions,
        const ParametricCellData &previous, CenterMethod method, int axis, const IncomingFacePlans &incoming, bool &changed) const ;
    /// Apply isotropic hex splits and coarsening on an independent plan; nullopt for other representations or split families.
    /// Requests name this state's cells; incoming face plans are borrowed for this call and merged before cell edits.
    /// Geometric selection then refines descendants to its cutoff. changed is true only after publishing a real edit.
    std::optional<ParametricCellData> editPlan(const std::vector<ParametricLeafRequest> &requests,int refinementLevels = 1,
        CenterMethod method = CenterMethod::reference,const IncomingFacePlans &incoming = {},
        bool *changed = nullptr,const GeometricRefinement *selection = nullptr) const ;
    /// Tags 0/1/2 keep, refine or request complete-family coarsening in active leaf order.
    ParametricRecord tagRequests(const std::vector<char> &tags) const ;
    /// Refine if any own vertex is tagged 1; coarsen only if all are tagged 2.
    ParametricRecord nodeTagRequests(const std::vector<int> &topology,const std::map<ParametricNodeKey,char> &tags) const ;

    // Read accepted state without expanding compact geometry where possible.
    std::size_t leafCount() const ;
    std::vector<int> leafDepths() const ;
    /// Active child-ordinal paths in output leaf order, without expanding compact hex geometry.
    std::vector<CellPath> leafPaths() const ;
    /// Longest leaf edge chord in 3D and XY; ignores neighboring edge subdivisions.
    std::vector<std::array<double,2>> edgeLengths() const ;
    ParametricRecord facePlan(EntityId face) const ;
    /// False when incoming divisions are empty or identical to the retained face plan.
    /// Other records still need merging and validation against the shared boundary.
    bool needsFaceMerge(EntityId face, const ParametricRecord &incoming) const ;
    /// Read accepted topology directly for hexes; other split families use their existing queries.
    std::vector<int> localRequests(const CellTopology &topology, unsigned difference, NeighborGrading grading, int axis) const ;
    ParametricTrace edgeTrace(const CellTopology &topology) const ;
    std::vector<int> vertexDepths(const CellTopology &topology) const ;
    std::vector<int> vertexRequests(const CellTopology &topology, const std::map<EntityId,int> &depths, unsigned difference) const ;
    ParametricTrace faceRequests(EntityId face, const CellTopology &leftTopology, const ParametricCellData &rightCell,
        const CellTopology &rightTopology, unsigned difference, NeighborGrading grading) const ;
    /// Active divided-face queries; compact hexes read split topology without physical geometry.
    std::vector<FaceBalanceRequest> dividedFaces(EntityId face,const CellTopology &leftTopology,const ParametricCellData &rightCell,
        const CellTopology &rightTopology,bool directions = false) const ;
    std::vector<int> faceBalanceRequests(const CellTopology &topology,const std::vector<FaceBalanceRequest> &shared,FaceBalance policy,
        std::map<int,unsigned> *directions = nullptr) const ;
    std::vector<std::vector<int>> donorsFrom(const ParametricCellData &previous) const ;
    /// Output uses encoded source topology; compact hexes use their retained topology.
    /// Cell-interior nodes, optionally with positions in the same order.
    std::vector<ParametricNodeKey> cellNodes(const std::vector<int> &topology,std::vector<Point> *positions = nullptr) const ;
    std::vector<Point> cellPositions(const std::vector<ParametricNodeKey> &nodes) const ;
    std::vector<ParametricOutputFace> cellFaces(const std::vector<int> &topology,const std::vector<ParametricNodeKey> &requiredNodes) const ;
    /// Original-face nodes, optionally with positions in the same order; compact hex pairs read accepted plans directly.
    std::vector<ParametricNodeKey> faceNodes(EntityId face,const std::vector<int> &topology,const ParametricCellData *right = nullptr,
        const std::vector<int> *rightTopology = nullptr,std::vector<Point> *positions = nullptr) const ;
    std::vector<Point> facePositions(EntityId face,const std::vector<int> &topology,const ParametricCellData *right = nullptr,
        const std::vector<int> *rightTopology = nullptr) const ;
    std::vector<ParametricOutputFace> faceFaces(EntityId face,const std::vector<int> &topology,const std::vector<ParametricNodeKey> &requiredNodes,
        const ParametricCellData *right = nullptr,const std::vector<int> *rightTopology = nullptr) const ;
    /// Retained face patches in the same order as cellFaces() and faceFaces().
    std::vector<FacePatch> cellPatches(const std::vector<int> &topology) const ;
    std::vector<FacePatch> facePatches(EntityId face,const std::vector<int> &topology,const ParametricCellData *right = nullptr,
        const std::vector<int> *rightTopology = nullptr) const ;
    /// Compact plan, or an empty plan for an original hex; null for other representations.
    const HexPlan *hexPlan() const ;
    bool compactHex() const ;
    // Serialization retains accepted constructions, not temporary expanded views.
    /// Accepted history bytes; empty for an untouched corner cell.
    ParametricRecord record() const ;
    int maxDepth() const ;
    int maxNodes() const ;
  } ;

  /// Loci transport of accepted history; invalid receives remain errors until collective checking.
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
