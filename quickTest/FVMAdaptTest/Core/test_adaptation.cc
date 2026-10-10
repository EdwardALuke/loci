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
#include <vector>

#include <Loci.h>
#include <FVMAdapt/hexcell.h>
#include <FVMAdapt/prism.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

// This library helper is defined in build_general_cell.cc, without a header
// declaration. Keep the regression here without changing the library API.
std::vector<Entity> reorder_edges(const const_store<int>& node_remap,
                                 const const_MapVec<2>& edge2node,
                                 const entitySet& edges);

/// Splitting a hex in one, two, or three local directions creates 2, 4, or 8
/// fine cells. Saving and rebuilding its plan must preserve that count.
TEST_CASE("hex split plans preserve fine-cell counts") {
  const int fine_cells[] = {1, 2, 2, 4, 2, 4, 4, 8};
  for(char split_code = 0; split_code < 8; ++split_code) {
    CAPTURE(int(split_code));
    HexCell cell;
    CHECK(cell.empty_resplit({split_code}) == fine_cells[int(split_code)]);

    std::list<HexCell*> leaves;
    cell.sort_leaves(leaves);
    CHECK(leaves.size() == fine_cells[int(split_code)]);

    HexCell restored;
    CHECK(restored.empty_resplit(cell.make_cellplan()) == leaves.size());
  }
}

/// A triangular prism splits axially into 2 cells, across its end faces into
/// 3, or both ways into 6. Its saved plan must reproduce the fine-cell count.
TEST_CASE("prism split plans preserve fine-cell counts") {
  const int fine_cells[] = {1, 2, 3, 6};
  for(char split_code = 0; split_code < 4; ++split_code) {
    CAPTURE(int(split_code));
    Prism cell;
    CHECK(cell.empty_resplit({split_code}) == fine_cells[int(split_code)]);

    std::list<Prism*> leaves;
    cell.sort_leaves(leaves);
    CHECK(leaves.size() == fine_cells[int(split_code)]);

    Prism restored;
    CHECK(restored.empty_resplit(cell.make_cellplan()) == leaves.size());
  }
}

/// Refine one of eight hex children, then derefine that child and the root.
/// Check the resulting leaf counts, not the bytes used to save each plan.
TEST_CASE("nested hex refinement can be replayed and derefined") {
  HexCell cell;
  // Split the root in all directions, then split its first child in one.
  REQUIRE(cell.empty_resplit({7, 1}) == 9);

  HexCell restored;
  CHECK(restored.empty_resplit(cell.make_cellplan()) == 9);
  REQUIRE(restored.numChildren() == 8);
  CHECK(restored.getChildCell(0)->numChildren() == 2);

  restored.getChildCell(0)->derefine();
  std::list<HexCell*> leaves;
  restored.sort_leaves(leaves);
  CHECK(leaves.size() == 8);

  restored.derefine();
  leaves.clear();
  restored.sort_leaves(leaves);
  REQUIRE(leaves.size() == 1);
  CHECK(leaves.front() == &restored);
  HexCell coarse;
  CHECK(coarse.empty_resplit(restored.make_cellplan()) == 1);
}

/// After splitting a triangular prism, its children have quadrilateral end
/// faces: refining one fully creates 8 grandchildren. Derefinement reverses it.
TEST_CASE("nested prism refinement can be replayed and derefined") {
  Prism cell;
  // Six children, with the first replaced by eight: thirteen leaves.
  REQUIRE(cell.empty_resplit({3, 3}) == 13);

  Prism restored;
  CHECK(restored.empty_resplit(cell.make_cellplan()) == 13);
  REQUIRE(restored.numChildren() == 6);
  CHECK(restored.getChildCell(0)->numChildren() == 8);

  restored.getChildCell(0)->derefine();
  std::list<Prism*> leaves;
  restored.sort_leaves(leaves);
  CHECK(leaves.size() == 6);

  restored.derefine();
  leaves.clear();
  restored.sort_leaves(leaves);
  REQUIRE(leaves.size() == 1);
  CHECK(leaves.front() == &restored);
  Prism coarse;
  CHECK(coarse.empty_resplit(restored.make_cellplan()) == 1);
}

/// A polygon splits into one quadrilateral per edge; refining one of those
/// children replaces it with four faces, regardless of the original polygon.
TEST_CASE("polygon face plans refine quadrilateral children") {
  for(int edges : {3, 4, 5}) {
    CAPTURE(edges);
    Face face(edges);
    CHECK(face.empty_resplit({1}) == edges);
    CHECK(face.empty_resplit({1, 1}) == edges + 3);
    CHECK(face.get_num_leaves() == edges + 3);
  }
}

/// Neighbors requesting perpendicular splits need four shared subfaces and
/// every boundary edge bisected, regardless of which neighbor is supplied first.
TEST_CASE("quad-face merging retains both neighbors' refinement") {
  for(bool reverse_neighbors : {false, true}) {
    CAPTURE(reverse_neighbors);
    std::vector<char> first = {1}; // Split in the face's local y direction.
    std::vector<char> second = {2}; // Split in the face's local x direction.
    const auto merged = reverse_neighbors
      ? merge_quad_face(second, 0, first, 0)
      : merge_quad_face(first, 0, second, 0);

    QuadFace face;
    std::vector<QuadFace*> fine_faces;
    face.empty_resplit(merged, 0, fine_faces);
    CHECK(fine_faces.size() == 4);

    for(unsigned int edge_id = 0; edge_id < 4; ++edge_id) {
      CAPTURE(edge_id);
      std::vector<char> edge_plan;
      extract_quad_edge(merged, edge_plan, edge_id);
      Node head(vect3d(0.0, 0.0, 0.0));
      Node tail(vect3d(1.0, 0.0, 0.0));
      std::list<Node*> new_nodes;
      {
        Edge edge(&head, &tail);
        edge.resplit(edge_plan, new_nodes);
        std::list<Edge*> leaves;
        edge.sort_leaves(leaves);
        CHECK(leaves.size() == 2);
        for(Edge* leaf : leaves)
          CHECK(leaf->get_length() == doctest::Approx(0.5));
      }
      cleanup_list(new_nodes);
    }
  }
}

/// When edges share their first endpoint, the other endpoint's remapped node
/// number must determine their order, rather than the edges' mesh entity IDs.
TEST_CASE("general-cell edge ordering uses both endpoint node numbers") {
  const Entity lower_id_edge = 4, higher_id_edge = 9;
  const Entity shared_node = 20, later_node = 21, earlier_node = 22;
  entitySet edges;
  edges += lower_id_edge;
  edges += higher_id_edge;
  entitySet nodes;
  nodes += shared_node;
  nodes += later_node;
  nodes += earlier_node;

  store<int> node_remap;
  node_remap.allocate(nodes);
  node_remap[shared_node] = 0;
  node_remap[later_node] = 2;
  node_remap[earlier_node] = 1;

  MapVec<2> edge2node;
  edge2node.allocate(edges);
  edge2node[lower_id_edge][0] = shared_node;
  edge2node[lower_id_edge][1] = later_node;
  edge2node[higher_id_edge][0] = shared_node;
  edge2node[higher_id_edge][1] = earlier_node;

  const_store<int> remap_view(node_remap.Rep());
  const_MapVec<2> edge_view(edge2node.Rep());
  const auto ordered = reorder_edges(remap_view, edge_view, edges);
  REQUIRE(ordered.size() == 2);
  CHECK(ordered[0] == higher_id_edge);
  CHECK(ordered[1] == lower_id_edge);
}

int main(int argc, char** argv) {
  Loci::Init(&argc, &argv);
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();
  Loci::Finalize();
  return result;
}
