// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef FVMADAPT_FACE_HISTORY_INTERNAL_H
#define FVMADAPT_FACE_HISTORY_INTERNAL_H

#include <FVMAdapt/face_history_internal.h>
#include <functional>

namespace Loci {
  class Map ;
  class refinedGridData ;
  struct RefinementOptions ;
  void checkFaceState(const RefinementOptions& options, const std::string& base, int level) ;
  void prepareFaceState(refinedGridData& grid, const std::string& base, int level,
                        const RefinementOptions& options, bool hasParents) ;
  bool installFaceHistory(fact_db& facts, const refinedGridData& grid,
                          const Map& faceNumbers, std::string& error) ;
  namespace face_history_detail {
    bool validPatch(const FacePatch& patch) ;
    bool covers(const FacePatch& face, const std::vector<FacePatch>& parts,
                std::string* reason = nullptr) ;
    long double patchMeasure(const FacePatch& patch) ;
    long double commonMeasure(const FacePatch& a, const FacePatch& b) ;
    FaceSource contribution(const HistoryFace& old, const HistoryFace& next) ;
    void matchPatches(const std::vector<HistoryFace>& previous,
      const std::vector<HistoryFace>& next,
      const std::function<void(size_t,size_t)>& match) ;
  }
}
#endif
