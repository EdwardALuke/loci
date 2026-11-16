// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef FVMADAPT_FACE_HISTORY_HELPERS_H
#define FVMADAPT_FACE_HISTORY_HELPERS_H

#include <FVMAdapt/face_history_internal.h>
#include <functional>

namespace Loci {
  class Map ;
  class refinedGridData ;

  /// Prepare face patches and cell donors for installation of the accepted mesh.
  void prepareParametricFaceState(refinedGridData& grid, const std::string& base, FaceStateP previous,
                                  storeRepP oldHistory, storeRepP newHistory, storeRepP cellParents = 0) ;

  /// Record original face patches in the installed mesh's numbering and orientation.
  FaceStateP recordBaseParametricFaces(fact_db& facts,const refinedGridData& grid,const std::string& name,storeRepP history) ;

  /// Redistribute recorded faces and install their state and donor correspondence.
  bool installFaceHistory(fact_db& facts, const refinedGridData& grid,
                          const Map& faceNumbers, std::string& error) ;

  namespace face_history_detail {
    /// Check patch coordinates, orientation and available physical area.
    bool validPatch(const FacePatch& patch) ;

    /// Test complete patch coverage without gaps or overlapping interiors.
    bool covers(const FacePatch& face, const std::vector<FacePatch>& parts,
                std::string* reason = nullptr) ;

    /// Area in reference coordinates, including all covered sectors.
    long double patchMeasure(const FacePatch& patch) ;

    /// Reference area shared by two patches on the same surface.
    long double commonMeasure(const FacePatch& a, const FacePatch& b) ;

    /// Test whether an overlap covers the full reference patch.
    bool wholePatch(long double common,const FacePatch& patch) ;

    /// Old-to-new orientation and old-face overlap fractions; unavailable physical fractions are -1.
    FaceSource contribution(const HistoryFace& old, const HistoryFace& next) ;

    /// Call match with the indices of each overlapping old and new face pair.
    void matchPatches(const std::vector<HistoryFace>& previous,
      const std::vector<HistoryFace>& next,
      const std::function<void(size_t,size_t)>& match) ;
  }
}
#endif
