//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
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

#include "fvmadapt_core.h"

// These are small test-owned cell topologies.
// They keep the nodes, edges, faces, and root cell together for cleanup.

namespace {

/// Allocate a test node and store it for cleanup.
template<class Fixture>
Node* add_node(Fixture& fixture, double x, double y, double z) {
  vect3d p(x, y, z) ;
  Node* node = new Node(p, int32(fixture.nodes.size() + 1)) ;
  fixture.nodes.push_back(node) ;
  return node ;
}

/// Allocate a test edge and store it for cleanup.
template<class Fixture>
Edge* add_edge(Fixture& fixture, Node* head, Node* tail) {
  Edge* edge = new Edge(head, tail) ;
  fixture.edges.push_back(edge) ;
  return edge ;
}
} // namespace

HexFixture::HexFixture(): root(0) {}

HexFixture::~HexFixture() {
  delete root ;
  cleanup_list(nodes, edges, faces) ;
}

PrismFixture::PrismFixture(): root(0) {}

PrismFixture::~PrismFixture() {
  delete root ;
  cleanup_list(nodes, edges) ;
  cleanup_list(qfaces) ;
  cleanup_list(gfaces) ;
}

GeneralFixture::GeneralFixture(): root(0) {}

GeneralFixture::~GeneralFixture() {
  delete root ;
  cleanup_list(nodes, edges, faces) ;
}

/// Build a unit hexahedron with the face ordering expected by HexCell.
void build_unit_hex(HexFixture& fixture) {
  Node* node[8] ;
  node[0] = add_node(fixture, 0.0, 0.0, 0.0) ;
  node[1] = add_node(fixture, 0.0, 0.0, 1.0) ;
  node[2] = add_node(fixture, 0.0, 1.0, 0.0) ;
  node[3] = add_node(fixture, 0.0, 1.0, 1.0) ;
  node[4] = add_node(fixture, 1.0, 0.0, 0.0) ;
  node[5] = add_node(fixture, 1.0, 0.0, 1.0) ;
  node[6] = add_node(fixture, 1.0, 1.0, 0.0) ;
  node[7] = add_node(fixture, 1.0, 1.0, 1.0) ;

  const int edge_head[12] = {0, 1, 2, 3, 0, 1, 4, 5, 0, 2, 4, 6} ;
  const int edge_tail[12] = {4, 5, 6, 7, 2, 3, 6, 7, 1, 3, 5, 7} ;

  Edge* edge[12] ;
  for(int i = 0; i < 12; ++i) {
    edge[i] = add_edge(fixture, node[edge_head[i]], node[edge_tail[i]]) ;
  }

  // Wire the 12 unit-hex edges into the six QuadFace objects.
  const int face_to_edge[6][4] = {
    {6, 11, 7, 10},
    {4, 9, 5, 8},
    {2, 11, 3, 9},
    {0, 10, 1, 8},
    {1, 7, 3, 5},
    {0, 6, 2, 4}
  } ;

  QuadFace** face = new QuadFace*[6] ;
  for(int i = 0; i < 6; ++i) {
    face[i] = new QuadFace(4) ;
    fixture.faces.push_back(face[i]) ;
    for(int j = 0; j < 4; ++j) {
      face[i]->edge[j] = edge[face_to_edge[i][j]] ;
    }
  }

  fixture.root = new HexCell(face) ;
}

/// Build a unit triangular prism with two general faces and three quad faces.
void build_unit_prism(PrismFixture& fixture) {
  Node* node[6] ;
  node[0] = add_node(fixture, 0.0, 0.0, 0.0) ;
  node[1] = add_node(fixture, 1.0, 0.0, 0.0) ;
  node[2] = add_node(fixture, 0.0, 1.0, 0.0) ;
  node[3] = add_node(fixture, 0.0, 0.0, 1.0) ;
  node[4] = add_node(fixture, 1.0, 0.0, 1.0) ;
  node[5] = add_node(fixture, 0.0, 1.0, 1.0) ;

  const int edge_head[9] = {0, 1, 2, 3, 4, 5, 0, 1, 2} ;
  const int edge_tail[9] = {1, 2, 0, 4, 5, 3, 3, 4, 5} ;

  Edge* edge[9] ;
  for(int i = 0; i < 9; ++i) {
    edge[i] = add_edge(fixture, node[edge_head[i]], node[edge_tail[i]]) ;
  }

  const int gface_to_edge[2][3] = {{0, 1, 2}, {3, 4, 5}} ;
  const int qface_to_edge[3][4] = {{0, 7, 3, 6}, {1, 8, 4, 7}, {2, 6, 5, 8}} ;

  Prism* prism = new Prism(3) ;
  for(int i = 0; i < 2; ++i) {
    Face* face = new Face(3) ;
    fixture.gfaces.push_back(face) ;
    for(int j = 0; j < 3; ++j) {
      face->edge[j] = edge[gface_to_edge[i][j]] ;
      face->needReverse[j] = false ;
    }
    prism->setFace(i, face) ;
  }

  for(int i = 0; i < 3; ++i) {
    QuadFace* face = new QuadFace(4) ;
    fixture.qfaces.push_back(face) ;
    for(int j = 0; j < 4; ++j) {
      face->edge[j] = edge[qface_to_edge[i][j]] ;
    }
    prism->setFace(i, face) ;
  }

  fixture.root = prism ;
}

/// Add an owned general Face to a general-cell fixture.
Face* add_general_face(GeneralFixture& fixture,
                       const std::vector<Edge*>& edges,
                       const std::vector<bool>& need_reverse) {
  Face* face = new Face(int(edges.size())) ;
  fixture.faces.push_back(face) ;
  for(size_t i = 0; i < edges.size(); ++i) {
    face->edge[i] = edges[i] ;
    face->needReverse[i] = need_reverse[i] ;
  }
  return face ;
}

/// Build a tetrahedron represented by the general Cell path.
void build_tetra_cell(GeneralFixture& fixture) {
  Node* node[4] ;
  node[0] = add_node(fixture, 0.0, 0.0, 0.0) ;
  node[1] = add_node(fixture, 1.0, 0.0, 0.0) ;
  node[2] = add_node(fixture, 0.0, 1.0, 0.0) ;
  node[3] = add_node(fixture, 0.0, 0.0, 1.0) ;

  Edge* edge[6] ;
  edge[0] = add_edge(fixture, node[0], node[1]) ;
  edge[1] = add_edge(fixture, node[1], node[2]) ;
  edge[2] = add_edge(fixture, node[2], node[0]) ;
  edge[3] = add_edge(fixture, node[0], node[3]) ;
  edge[4] = add_edge(fixture, node[1], node[3]) ;
  edge[5] = add_edge(fixture, node[2], node[3]) ;

  Face** face = new Face*[4] ;
  face[0] = add_general_face(fixture, {edge[0], edge[1], edge[2]},
                             {false, false, false}) ;
  face[1] = add_general_face(fixture, {edge[0], edge[4], edge[3]},
                             {false, false, true}) ;
  face[2] = add_general_face(fixture, {edge[1], edge[5], edge[4]},
                             {false, false, true}) ;
  face[3] = add_general_face(fixture, {edge[2], edge[3], edge[5]},
                             {false, false, true}) ;

  Node** cell_nodes = new Node*[4] ;
  for(int i = 0; i < 4; ++i) {
    cell_nodes[i] = node[i] ;
  }

  Edge** cell_edges = new Edge*[6] ;
  for(int i = 0; i < 6; ++i) {
    cell_edges[i] = edge[i] ;
  }

  char* face_orient = new char[4] ;
  for(int i = 0; i < 4; ++i) {
    face_orient[i] = 0 ;
  }

  fixture.root = new Cell(4, 6, 4, cell_nodes, cell_edges, face, face_orient) ;
}

/// Build a square pyramid represented by the general Cell path.
void build_pyramid_cell(GeneralFixture& fixture) {
  Node* node[5] ;
  node[0] = add_node(fixture, 0.0, 0.0, 0.0) ;
  node[1] = add_node(fixture, 1.0, 0.0, 0.0) ;
  node[2] = add_node(fixture, 1.0, 1.0, 0.0) ;
  node[3] = add_node(fixture, 0.0, 1.0, 0.0) ;
  node[4] = add_node(fixture, 0.5, 0.5, 1.0) ;

  Edge* edge[8] ;
  edge[0] = add_edge(fixture, node[0], node[1]) ;
  edge[1] = add_edge(fixture, node[1], node[2]) ;
  edge[2] = add_edge(fixture, node[2], node[3]) ;
  edge[3] = add_edge(fixture, node[3], node[0]) ;
  edge[4] = add_edge(fixture, node[0], node[4]) ;
  edge[5] = add_edge(fixture, node[1], node[4]) ;
  edge[6] = add_edge(fixture, node[2], node[4]) ;
  edge[7] = add_edge(fixture, node[3], node[4]) ;

  Face** face = new Face*[5] ;
  face[0] = add_general_face(fixture, {edge[0], edge[1], edge[2], edge[3]},
                             {false, false, false, false}) ;
  face[1] = add_general_face(fixture, {edge[0], edge[5], edge[4]},
                             {false, false, true}) ;
  face[2] = add_general_face(fixture, {edge[1], edge[6], edge[5]},
                             {false, false, true}) ;
  face[3] = add_general_face(fixture, {edge[2], edge[7], edge[6]},
                             {false, false, true}) ;
  face[4] = add_general_face(fixture, {edge[3], edge[4], edge[7]},
                             {false, false, true}) ;

  Node** cell_nodes = new Node*[5] ;
  for(int i = 0; i < 5; ++i) {
    cell_nodes[i] = node[i] ;
  }

  Edge** cell_edges = new Edge*[8] ;
  for(int i = 0; i < 8; ++i) {
    cell_edges[i] = edge[i] ;
  }

  char* face_orient = new char[5] ;
  for(int i = 0; i < 5; ++i) {
    face_orient[i] = 0 ;
  }

  fixture.root = new Cell(5, 8, 5, cell_nodes, cell_edges, face, face_orient) ;
}
