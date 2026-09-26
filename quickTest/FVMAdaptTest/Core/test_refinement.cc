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

#include <list>
#include <memory>

#include <Loci.h>
#include <FVMAdapt/hexcell.h>
#include <FVMAdapt/prism.h>
#include <doctest.h>

namespace {

// A cell owns its child cells, but shares its geometry. Keep the root geometry
// and anything created by resplit() alive until after the cell is destroyed.
struct CellGeometry {
  std::list<Node*> nodes;
  std::list<Edge*> edges;
  std::list<QuadFace*> quad_faces;
  std::list<Face*> faces;

  ~CellGeometry() {
    cleanup_list(nodes, edges);
    cleanup_list(quad_faces);
    cleanup_list(faces);
  }

  Node* add_node(double x, double y, double z) {
    nodes.push_back(new Node(vect3d(x, y, z)));
    return nodes.back();
  }

  Edge* add_edge(Node* head, Node* tail) {
    edges.push_back(new Edge(head, tail));
    return edges.back();
  }
};

// Unit cube, using the local node/edge/face numbering expected by HexCell.
std::unique_ptr<HexCell> unit_hex(CellGeometry& geometry) {
  Node* nodes[8] = {
    geometry.add_node(0, 0, 0), geometry.add_node(0, 0, 1),
    geometry.add_node(0, 1, 0), geometry.add_node(0, 1, 1),
    geometry.add_node(1, 0, 0), geometry.add_node(1, 0, 1),
    geometry.add_node(1, 1, 0), geometry.add_node(1, 1, 1)
  };
  const int edge_nodes[12][2] = {
    {0, 4}, {1, 5}, {2, 6}, {3, 7}, {0, 2}, {1, 3},
    {4, 6}, {5, 7}, {0, 1}, {2, 3}, {4, 5}, {6, 7}
  };
  Edge* edges[12];
  for(int i = 0; i < 12; ++i)
    edges[i] = geometry.add_edge(nodes[edge_nodes[i][0]], nodes[edge_nodes[i][1]]);

  const int face_edges[6][4] = {
    {6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9},
    {0, 10, 1, 8}, {1, 7, 3, 5}, {0, 6, 2, 4}
  };
  QuadFace** faces = new QuadFace*[6];
  for(int i = 0; i < 6; ++i) {
    faces[i] = new QuadFace(4);
    geometry.quad_faces.push_back(faces[i]);
    for(int j = 0; j < 4; ++j)
      faces[i]->edge[j] = edges[face_edges[i][j]];
  }
  return std::make_unique<HexCell>(faces);
}

// Right triangular prism, with two triangular ends and three quad faces.
std::unique_ptr<Prism> unit_prism(CellGeometry& geometry) {
  Node* nodes[6] = {
    geometry.add_node(0, 0, 0), geometry.add_node(1, 0, 0),
    geometry.add_node(0, 1, 0), geometry.add_node(0, 0, 1),
    geometry.add_node(1, 0, 1), geometry.add_node(0, 1, 1)
  };
  const int edge_nodes[9][2] = {
    {0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3},
    {0, 3}, {1, 4}, {2, 5}
  };
  Edge* edges[9];
  for(int i = 0; i < 9; ++i)
    edges[i] = geometry.add_edge(nodes[edge_nodes[i][0]], nodes[edge_nodes[i][1]]);

  auto cell = std::make_unique<Prism>(3);
  const int end_edges[2][3] = {{0, 1, 2}, {3, 4, 5}};
  for(int i = 0; i < 2; ++i) {
    Face* face = new Face(3);
    geometry.faces.push_back(face);
    for(int j = 0; j < 3; ++j) {
      face->edge[j] = edges[end_edges[i][j]];
      face->needReverse[j] = false;
    }
    cell->setFace(i, face);
  }
  const int side_edges[3][4] = {{0, 7, 3, 6}, {1, 8, 4, 7}, {2, 6, 5, 8}};
  for(int i = 0; i < 3; ++i) {
    QuadFace* face = new QuadFace(4);
    geometry.quad_faces.push_back(face);
    for(int j = 0; j < 4; ++j)
      face->edge[j] = edges[side_edges[i][j]];
    cell->setFace(i, face);
  }
  return cell;
}

} // namespace

/// Two levels must refine all eight hex children again, producing 64 leaves,
/// rather than spending the requested levels on only the first few cells.
TEST_CASE("two refinement levels reach every hex child") {
  CellGeometry geometry;
  auto cell = unit_hex(geometry);
  cell->resplit(2, geometry.nodes, geometry.edges, geometry.quad_faces);

  std::list<HexCell*> leaves;
  cell->sort_leaves(leaves);
  CHECK(leaves.size() == 64);
  REQUIRE(cell->numChildren() == 8);
  for(int child = 0; child < 8; ++child) {
    CAPTURE(child);
    CHECK(cell->getChildCell(child)->numChildren() == 8);
  }
}

/// A triangular prism first makes six children with quadrilateral end faces;
/// refining every child again must produce 48 leaves, not just refine one branch.
TEST_CASE("two refinement levels reach every prism child") {
  CellGeometry geometry;
  auto cell = unit_prism(geometry);
  cell->resplit(2, geometry.nodes, geometry.edges, geometry.quad_faces, geometry.faces);

  std::list<Prism*> leaves;
  cell->sort_leaves(leaves);
  CHECK(leaves.size() == 48);
  REQUIRE(cell->numChildren() == 6);
  for(int child = 0; child < 6; ++child) {
    CAPTURE(child);
    CHECK(cell->getChildCell(child)->numChildren() == 8);
  }
}
