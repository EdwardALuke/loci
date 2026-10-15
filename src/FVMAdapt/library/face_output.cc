// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.

#include <algorithm>
#include "../face_output.h"
#include "diamondcell.h"
#include "hexcell.h"
#include "prism.h"

namespace Loci {
  FaceOutput::FaceOutput(FineFaces& faces): faces(faces), patches(0), root(0) {
    FineFaces().swap(faces) ;
  }

  FaceOutput::FaceOutput(FineFaces& faces, FacePatches& patches, int root):
    faces(faces), patches(&patches), patchMap(std::in_place), root(root) {
    FineFaces().swap(faces) ;
    patches.clear() ;
  }

  void FaceOutput::record(const Cell* cell) {
    if(patches != 0) patchMap->add(cell,root) ;
  }

  void FaceOutput::record(const HexCell* cell) {
    if(patches != 0) patchMap->add(cell,root) ;
  }

  void FaceOutput::record(const Prism* cell) {
    if(patches != 0) patchMap->add(cell,root) ;
  }

  void FaceOutput::reserve(size_t count) {
    faces.reserve(count) ;
    if(patches != 0) patches->reserve(count) ;
  }

  void FaceOutput::append(Face* face, int left, int right) {
    std::list<int32> vertices ;
    face->set_f2n(vertices) ;
    append(vertices,left,right,patches != 0 ? &patchMap->at(face) : 0) ;
  }

  void FaceOutput::append(QuadFace* face, int left, int right) {
    std::list<int32> vertices ;
    face->set_f2n(vertices) ;
    append(vertices,left,right,patches != 0 ? &patchMap->at(face) : 0) ;
  }

  void FaceOutput::append(const std::list<int32>& vertices, int left, int right,
                          const FacePatch* patch) {
    faces.emplace_back(vertices.size()+2) ;
    auto& row = faces.back() ;
    row[0] = left ; row[1] = right ;
    std::copy(vertices.begin(),vertices.end(),row.begin()+2) ;
    if(patch != 0) patches->push_back(*patch) ;
  }

  void make_hex_cell_faces(HexCell* cell, const std::vector<char>& plan,
    int nodeOffset, int cellOffset, std::list<Node*>& nodes,
    std::list<Edge*>& edges, std::list<QuadFace*>& faces, FaceOutput& output) {
    std::vector<HexCell*> cells ;
    cell->resplit(plan,nodes,edges,faces,cells) ;
    output.record(cell) ;

    // Number new nodes in the order used by the coordinate output.
    int index = 0 ;
    for(Node* node : nodes) node->index = nodeOffset+index++ ;

    std::map<QuadFace*,NeibIndex> innerFaces ;
    set_hex_faces(cells,innerFaces) ;
    output.reserve(innerFaces.size()) ;
    for(const auto& face : innerFaces)
      output.append(face.first,cellOffset+face.second.c1,cellOffset+face.second.c2) ;
  }

  void make_prism_cell_faces(Prism* cell, const std::vector<char>& plan,
    int nodeOffset, int cellOffset, std::list<Node*>& nodes,
    std::list<Edge*>& edges, std::list<QuadFace*>& quadFaces,
    std::list<Face*>& generalFaces, FaceOutput& output) {
    std::vector<Prism*> cells ;
    cell->resplit(plan,nodes,edges,quadFaces,generalFaces,cells) ;
    output.record(cell) ;

    int index = 0 ;
    for(Node* node : nodes) node->index = nodeOffset+index++ ;

    std::map<QuadFace*,NeibIndex> innerQuads ;
    std::map<Face*,NeibIndex> innerFaces ;
    set_prism_faces(cells,innerQuads,innerFaces) ;
    output.reserve(innerQuads.size()+innerFaces.size()) ;
    for(const auto& face : innerQuads)
      output.append(face.first,cellOffset+face.second.c1,cellOffset+face.second.c2) ;
    for(const auto& face : innerFaces)
      output.append(face.first,cellOffset+face.second.c1,cellOffset+face.second.c2) ;
  }

  void make_general_cell_faces(Cell* cell, const std::vector<char>& plan,
    int nodeOffset, int cellOffset, std::list<Node*>& nodes,
    std::list<Edge*>& edges, std::list<Face*>& faces, FaceOutput& output) {
    std::vector<DiamondCell*> cells ;
    cell->resplit(plan,nodes,edges,faces,cells) ;
    output.record(cell) ;

    int index = 0 ;
    for(Node* node : nodes) node->index = nodeOffset+index++ ;

    auto nodeEdges = cell->set_n2e() ;
    std::list<std::pair<Face*,NeibIndex> > innerFaces ;
    set_general_faces(cell,cells,nodeEdges,innerFaces) ;
    output.reserve(innerFaces.size()) ;
    for(const auto& face : innerFaces)
      output.append(face.first,cellOffset+face.second.c1,cellOffset+face.second.c2) ;
  }
}
