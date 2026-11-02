// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef FVMADAPT_FACE_HISTORY_DATA_H
#define FVMADAPT_FACE_HISTORY_DATA_H

// Implementation records; solver interfaces are in face_history.h.
#include <FVMAdapt/face_history.h>
#include <array>
#include <mpi.h>

namespace Loci {
  /// Identifies a face region so subdivisions can be matched between meshes.
  /// Coordinates describe logical refinement, not physical positions.
  struct FacePatch {
    std::vector<long long> surface ; // Opaque surface identity within one base mesh.
    std::array<long long,4> bounds{} ; // xmin, ymin, xmax, ymax
    int sector = 0 ; // Negative values denote an unsplit polygon's sector count.
    int orientation = 1 ;
  } ;

  typedef std::vector<FacePatch> FacePatches ;

  class FacePatchConverter {
    FacePatch& value ;
  public:
    explicit FacePatchConverter(FacePatch& v): value(v) {}
    int getSize() const { return 6+value.surface.size() ; }
    void getState(long long* buffer, int& size) const ;
    void setState(const long long* buffer, int size) ;
  } ;

  template<> struct data_schema_traits<FacePatch> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef long long Converter_Base_Type ;
    typedef FacePatchConverter Converter_Type ;
  } ;

  class FacePatchesConverter {
    FacePatches& value ;
  public:
    explicit FacePatchesConverter(FacePatches& v): value(v) {}
    int getSize() const ;
    void getState(long long* buffer, int& size) const ;
    void setState(const long long* buffer, int size) ;
  } ;

  template<> struct data_schema_traits<FacePatches> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef long long Converter_Base_Type ;
    typedef FacePatchesConverter Converter_Type ;
  } ;

  std::ostream& operator<<(std::ostream&, const FacePatch&) ;
  std::istream& operator>>(std::istream&, FacePatch&) ;

  /// One oriented face; right is a cell ID or a negative boundary ID.
  struct HistoryFace {
    FaceId face ;
    CellId left, right ;
    FacePatch patch ;
  } ;

  class HistoryFaceConverter {
    HistoryFace& value ;
  public:
    explicit HistoryFaceConverter(HistoryFace& v): value(v) {}
    int getSize() const { return 3+FacePatchConverter(value.patch).getSize() ; }
    void getState(long long* buffer, int& size) const ;
    void setState(const long long* buffer, int size) ;
  } ;

  template<> struct data_schema_traits<HistoryFace> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef long long Converter_Base_Type ;
    typedef HistoryFaceConverter Converter_Type ;
  } ;

  std::ostream& operator<<(std::ostream&, const HistoryFace&) ;
  std::istream& operator>>(std::istream&, HistoryFace&) ;

  /// Sources, change, and origin follow the local new-face order.
  struct FaceMatches {
    /// Contributors in each row are sorted by old-face ID.
    std::vector<std::vector<FaceSource> > sources ;
    std::vector<FaceChange> change ;
    /// Previous cell for FACE_CREATED; -1 for all other faces.
    std::vector<CellId> origin ;
    std::vector<FaceId> removed ; // Old internal faces removed by cell coarsening.
    std::vector<CellId> removedInto ; // New cell, in removed-face order.
  } ;

  struct FaceExchangeStats {
    size_t localFaces = 0, matchingFaces = 0, contributions = 0 ;
    size_t sentWords = 0, receivedWords = 0 ;
  } ;

  /// Collective matching of locally owned faces from two partitions of one base mesh.
  /// Cell pairs may be distributed independently; each (new, old) pair occurs once.
  /// Rows follow local new faces; removed faces belong to their old owner.
  /// Faces must be uniquely owned. Cell IDs are nonnegative; boundary IDs are negative.
  /// This matches topology only: no periodic transforms, mesh motion, or field reconstruction.
  bool buildFaceHistory(const std::vector<HistoryFace>& previous,
                        const std::vector<HistoryFace>& next,
                        const std::vector<std::pair<CellId,CellId> >& cell2parent,
                        FaceMatches& history, std::string& error, MPI_Comm comm,
                        FaceExchangeStats* stats = nullptr) ;

  /// Reusable requests by face ID; owned IDs follow the input store's domain order.
  class FaceDataExchange {
    struct Impl ;
    std::shared_ptr<Impl> impl ;
  public:
    bool setup(const std::vector<FaceId>& owned, const std::vector<FaceId>& needed,
               std::string& error, MPI_Comm comm) ;
    const std::vector<FaceId>& sources() const ;
    /// Output follows sources(), with no duplicate donors.
    bool exchange(storeRepP input, storeRepP output, const entitySet& owned,
                  int keyspace, int width, std::string& error) const ;
  } ;


  class refinedGridData ;
  store<FacePatch>& recordedFacePatches(refinedGridData& grid) ;
  const store<FacePatch>& recordedFacePatches(const refinedGridData& grid) ;
}
#endif
