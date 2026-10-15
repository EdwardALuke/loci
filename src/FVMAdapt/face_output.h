// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef FVMADAPT_FACE_OUTPUT_H
#define FVMADAPT_FACE_OUTPUT_H

#include <list>
#include <optional>
#include <FVMAdapt/defines.h>
#include "face_patch.h"

class Node ;
class Edge ;

namespace Loci {
  /// Append connectivity and optional patches in the same face order.
  class FaceOutput {
    FineFaces& faces ;
    FacePatches* patches ;
    std::optional<FacePatchMap> patchMap ;
    int root ;
  public:
    explicit FaceOutput(FineFaces& faces) ;
    FaceOutput(FineFaces& faces, FacePatches& patches, int root) ;
    FaceOutput(const FaceOutput&) = delete ;
    FaceOutput& operator=(const FaceOutput&) = delete ;
    void record(const Cell* cell) ;
    void record(const HexCell* cell) ;
    void record(const Prism* cell) ;
    void reserve(size_t count) ;
    void append(Face* face, int left, int right) ;
    void append(QuadFace* face, int left, int right) ;
  private:
    void append(const std::list<int32>& vertices, int left, int right,
                const FacePatch* patch) ;
  } ;

  /// Split a local cell tree and collect its internal faces. The caller owns the tree and lists.
  void make_hex_cell_faces(HexCell* cell, const std::vector<char>& plan,
    int nodeOffset, int cellOffset, std::list<Node*>& nodes,
    std::list<Edge*>& edges, std::list<QuadFace*>& faces, FaceOutput& output) ;
  void make_prism_cell_faces(Prism* cell, const std::vector<char>& plan,
    int nodeOffset, int cellOffset, std::list<Node*>& nodes,
    std::list<Edge*>& edges, std::list<QuadFace*>& quadFaces,
    std::list<Face*>& generalFaces, FaceOutput& output) ;
  void make_general_cell_faces(Cell* cell, const std::vector<char>& plan,
    int nodeOffset, int cellOffset, std::list<Node*>& nodes,
    std::list<Edge*>& edges, std::list<Face*>& faces, FaceOutput& output) ;
}
#endif
