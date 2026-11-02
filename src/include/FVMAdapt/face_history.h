//#############################################################################
//#
//# Copyright 2026, Mississippi State University
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
#ifndef FVMADAPT_FACE_HISTORY_H
#define FVMADAPT_FACE_HISTORY_H

#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <store.h>
#include <storeVec.h>
#include <multiStore.h>

namespace Loci {
  typedef long long FaceId ;
  typedef long long CellId ;

  enum FaceChange {
    FACE_RETAINED, // Covers the same patch as one old face.
    FACE_SPLIT, // Covers part of one old face.
    FACE_MERGED, // Combines whole old faces.
    FACE_CREATED, // New face inside one old cell; no old-face source.
    FACE_RESHAPED // Old and new subdivisions overlap many-to-many.
  } ;

  /// One old face contributing to a new face.
  struct FaceSource {
    FaceId source ; // Old mesh face file number.
    int orientation ; // +1 for matching vertex-loop normals, -1 for opposing.
    double referenceFraction ; // Overlap / old reference-patch measure; not physical area.
  } ;

  std::ostream& operator<<(std::ostream&, const FaceSource&) ;
  std::istream& operator>>(std::istream&, FaceSource&) ;

  template<> struct data_schema_traits<FaceSource> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() ;
  } ;

  class fact_db ;
  struct FaceHistoryAccess ;

  /// One mesh's face information retained for the next adaptation; no solver fields.
  class FaceState {
    struct Impl ;
    std::shared_ptr<const Impl> impl ;
    friend struct FaceHistoryAccess ;
    friend class FaceHistory ;
    FaceState() = default ;
  } ;

  /// Relates new faces to the preceding mesh after refinement or coarsening.
  /// Stores use installed new-face entities.
  /// Face and cell IDs are signed file numbers in their respective meshes.
  /// Replacing existing split directions, periodic transforms, and mesh-motion weights are unsupported.
  class FaceHistory {
    struct Impl ;
    std::shared_ptr<const Impl> impl ;
    friend struct FaceHistoryAccess ;
    FaceHistory() = default ;
    bool gather(const FaceState& previous, storeRepP input, storeRepP output,
                int width, std::string& error) const ;
  public:
    const_multiStore<FaceSource> sources() const ;
    const_store<int> change() const ;
    /// Old cell file number for FACE_CREATED; -1 otherwise. Do not test its sign.
    const_store<CellId> origin() const ;
    /// Removed old faces and their new cell file numbers, on the old owning rank.
    const std::vector<std::pair<FaceId,CellId> >& removed() const ;
    /// Rank-local index in a gathered field; never store this index as face metadata.
    int donorIndex(FaceId source) const ;

    /// Retrieve only locally needed donors, without weighting or reorientation.
    /// Collective outside queries; replaces output in donor-index order.
    /// Output contents are unspecified on failure; input and output must not alias.
    /// Input uses the retained mesh's installed numbering on MPI_COMM_WORLD.
    template<class T>
    bool gather(const FaceState& previous, const_store<T> input,
                store<T>& donors, std::string& error) const {
      return gather(previous,input.Rep(),donors.Rep(),1,error) ;
    }
    template<class T>
    bool gather(const FaceState& previous, const_storeVec<T> input,
                storeVec<T>& donors, std::string& error) const {
      return gather(previous,input.Rep(),donors.Rep(),input.vecSize(),error) ;
    }
  } ;

  /// Face state to retain for the next adaptation; null when recording was not requested.
  std::shared_ptr<const FaceState> getFaceState(fact_db& facts) ;
  /// History of the latest adaptation; null when correspondence was not requested.
  std::shared_ptr<const FaceHistory> getFaceHistory(fact_db& facts) ;
  /// Record an unrefined VOG mesh's face state in memory before its first adaptation.
  /// Pass as previousFaceState with faceHistory enabled; later use getFaceState() from the new mesh.
  /// Collective outside queries; does not write a file.
  std::shared_ptr<const FaceState> recordBaseFaces(fact_db& facts, const std::string& casename) ;
  typedef std::shared_ptr<const FaceState> FaceStateP ;
  typedef std::shared_ptr<const FaceHistory> FaceHistoryP ;

}
#endif
