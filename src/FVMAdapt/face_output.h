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
  /// Construction clears the supplied output containers.
  class FaceOutput {
    FineFaces& faces ;
    FacePatches* patches ;
    std::optional<FacePatchMap> patchMap ;
    int root ; // Base cell or face file number.
  public:
    explicit FaceOutput(FineFaces& faces) ;
    FaceOutput(FineFaces& faces, FacePatches& patches, int root) ;
    FaceOutput(const FaceOutput&) = delete ;
    FaceOutput& operator=(const FaceOutput&) = delete ;
    void record(const Cell* cell) ;
    void record(const HexCell* cell) ;
    void record(const Prism* cell) ;
    void record(const Face* face) ;
    void record(const QuadFace* face) ;
    /// Split a face and number its new nodes in coordinate-output order.
    std::vector<Face*> split(Face* face, const std::vector<char>& plan,
      int nodeOffset, std::list<Node*>& nodes, std::list<Edge*>& edges) ;
    std::vector<QuadFace*> split(QuadFace* face, const std::vector<char>& plan,
      int nodeOffset, std::list<Node*>& nodes, std::list<Edge*>& edges) ;
    void reserve(size_t count) ;
    void append(Face* face, int left, int right) ;
    void append(QuadFace* face, int left, int right) ;
    /// Convert adjacent fine-cell indices to file numbers; empty plans select the root.
    void append(const std::vector<Face*>& leaves,
      const std::vector<int32>& left, const std::vector<char>& leftPlan, int leftOffset,
      const std::vector<int32>& right, const std::vector<char>& rightPlan, int rightOffset) ;
    void append(const std::vector<QuadFace*>& leaves,
      const std::vector<int32>& left, const std::vector<char>& leftPlan, int leftOffset,
      const std::vector<int32>& right, const std::vector<char>& rightPlan, int rightOffset) ;
    /// Boundary IDs are already negative in the connectivity output.
    void append(const std::vector<Face*>& leaves,
      const std::vector<int32>& left, const std::vector<char>& leftPlan, int leftOffset,
      int boundary) ;
    void append(const std::vector<QuadFace*>& leaves,
      const std::vector<int32>& left, const std::vector<char>& leftPlan, int leftOffset,
      int boundary) ;
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
