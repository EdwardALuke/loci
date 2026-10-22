// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.

#ifndef FVMADAPT_FACE_PATCH_H
#define FVMADAPT_FACE_PATCH_H

#include <map>
#include <FVMAdapt/face_history_internal.h>

class Face ;
class QuadFace ;
class Cell ;
class DiamondCell ;
class HexCell ;
class Prism ;

namespace Loci {

  /// Temporary lookup from tree faces to their refinement patches.
  /// Does not own the tree faces.
  class FacePatchMap {
    std::map<const Face*,FacePatch> faces ;
    std::map<const QuadFace*,FacePatch> quads ;
  public:
    void add(const Face* face, int faceNumber) ;
    void add(const QuadFace* face, int faceNumber) ;
    void add(const Cell* cell, int root) ;
    void add(const HexCell* cell, int root) ;
    void add(const Prism* cell, int root) ;
    const FacePatch& at(const Face* face) const { return faces.at(face) ; }
    const FacePatch& at(const QuadFace* face) const { return quads.at(face) ; }
  private:
    struct PrismRegion ;
    void add(const Face* face, const std::vector<long long>& surface) ;
    void polygon(const Face* face, const FacePatch& patch,
                 const std::array<std::array<long long,2>,4>& corners) ;
    void quad(const QuadFace* face, const FacePatch& patch) ;
    void general(const std::vector<const DiamondCell*>& children,
                 const std::vector<long long>& parent) ;
    void hex(const HexCell* cell, int root, const std::array<long long,6>& box) ;
    void prism(const Prism* cell, int root, const PrismRegion& region) ;
  } ;

}
#endif
