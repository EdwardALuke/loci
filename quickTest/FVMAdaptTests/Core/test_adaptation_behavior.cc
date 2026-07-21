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

#include <doctest.h>

#include <cmath>
#include <initializer_list>
#include <set>

namespace {

/// Keep split plans readable while the library continues to store codes as char.
std::vector<char> plan(std::initializer_list<int> codes) {
  std::vector<char> result ;
  for(std::initializer_list<int>::const_iterator code = codes.begin() ;
      code != codes.end(); ++code) {
    result.push_back(char(*code)) ;
  }
  return result ;
}

} // namespace


/// Edges sharing one endpoint must still be ordered by the remapped number of
/// the other endpoint, rather than falling back to their mesh entity numbers.
TEST_CASE("general-cell edge ordering uses both endpoint node numbers") {
  const Entity lower_id_edge = 4 ;
  const Entity higher_id_edge = 9 ;
  const Entity shared_node = 20 ;
  const Entity later_node = 21 ;
  const Entity earlier_node = 22 ;

  entitySet edges ;
  edges += lower_id_edge ;
  edges += higher_id_edge ;

  entitySet nodes ;
  nodes += shared_node ;
  nodes += later_node ;
  nodes += earlier_node ;

  store<int> node_remap ;
  node_remap.allocate(nodes) ;
  node_remap[shared_node] = 0 ;
  node_remap[later_node] = 2 ;
  node_remap[earlier_node] = 1 ;

  MapVec<2> edge2node ;
  edge2node.allocate(edges) ;
  edge2node[lower_id_edge][0] = shared_node ;
  edge2node[lower_id_edge][1] = later_node ;
  edge2node[higher_id_edge][0] = shared_node ;
  edge2node[higher_id_edge][1] = earlier_node ;

  const_store<int> remap_view(node_remap.Rep()) ;
  const_MapVec<2> edge_view(edge2node.Rep()) ;
  const std::vector<Entity> ordered =
    reorder_edges(remap_view, edge_view, edges) ;

  REQUIRE(ordered.size() == 2) ;
  CHECK(ordered[0] == higher_id_edge) ;
  CHECK(ordered[1] == lower_id_edge) ;
}


/// A cell is selected for refinement when any node is tagged `1`, selected for
/// derefinement only when every node is tagged `2`, and otherwise is unchanged.
/// This is the common node-tag contract for hex, prism, and general cells.
TEST_CASE("node tags classify each supported cell type consistently") {
  SUBCASE("HexCell") {
    HexFixture hex ;
    build_unit_hex(hex) ;

    CHECK(hex.root->get_tagged() == 0) ;

    hex.nodes.front()->tag = 2 ;
    CHECK(hex.root->get_tagged() == 0) ;

    hex.nodes.front()->tag = 1 ;
    CHECK(hex.root->get_tagged() == 1) ;

    for(std::list<Node*>::iterator node = hex.nodes.begin() ;
        node != hex.nodes.end(); ++node) {
      (*node)->tag = 2 ;
    }
    CHECK(hex.root->get_tagged() == 2) ;
  }

  SUBCASE("Prism") {
    PrismFixture prism ;
    build_unit_prism(prism) ;

    CHECK(prism.root->get_tagged() == 0) ;

    prism.nodes.front()->tag = 2 ;
    CHECK(prism.root->get_tagged() == 0) ;

    prism.nodes.front()->tag = 1 ;
    CHECK(prism.root->get_tagged() == 1) ;

    for(std::list<Node*>::iterator node = prism.nodes.begin() ;
        node != prism.nodes.end(); ++node) {
      (*node)->tag = 2 ;
    }
    CHECK(prism.root->get_tagged() == 2) ;
  }

  SUBCASE("general Cell") {
    GeneralFixture tetra ;
    build_tetra_cell(tetra) ;

    CHECK(tetra.root->get_tagged() == 0) ;

    tetra.nodes.front()->tag = 2 ;
    CHECK(tetra.root->get_tagged() == 0) ;

    tetra.nodes.front()->tag = 1 ;
    CHECK(tetra.root->get_tagged() == 1) ;

    for(std::list<Node*>::iterator node = tetra.nodes.begin() ;
        node != tetra.nodes.end(); ++node) {
      (*node)->tag = 2 ;
    }
    CHECK(tetra.root->get_tagged() == 2) ;
  }
}


/// Quad-face orientation changes local edge ids at cell/face boundaries. The
/// two conversion functions must remain exact inverses for every orientation.
TEST_CASE("quad-face edge orientation round-trips between cell and face ids") {
  for(char orientation = 0; orientation < 8; ++orientation) {
    for(char edge = 0; edge < 4; ++edge) {
      CAPTURE(int(orientation)) ;
      CAPTURE(int(edge)) ;
      CHECK(orient_edgeID_f2c(orient_edgeID_c2f(edge, orientation),
                              orientation) == edge) ;
      CHECK(orient_edgeID_c2f(orient_edgeID_f2c(edge, orientation),
                              orientation) == edge) ;
    }
  }
}


/// Merging two neighboring quad-face requests must retain every requested
/// boundary refinement. Reversing which neighbor is supplied first must give
/// the same merged behavior without assuming a particular merged split code.
TEST_CASE("quad-face merging preserves both requests regardless of side order") {
  std::vector<char> first_request = plan({1}) ;
  std::vector<char> second_request = plan({2}) ;

  std::vector<char> first_copy = first_request ;
  std::vector<char> second_copy = second_request ;
  const std::vector<char> merged =
    merge_quad_face(first_copy, char(0), second_copy, char(0)) ;

  first_copy = first_request ;
  second_copy = second_request ;
  const std::vector<char> reverse_merged =
    merge_quad_face(second_copy, char(0), first_copy, char(0)) ;

  // Commutativity is checked against the other merge result, not a hard-coded
  // representation of the combined face plan.
  CHECK(merged == reverse_merged) ;

  for(unsigned int edge = 0; edge < 4; ++edge) {
    std::vector<char> first_edge ;
    std::vector<char> second_edge ;
    std::vector<char> merged_edge ;
    std::vector<char> reverse_merged_edge ;
    extract_quad_edge(first_request, first_edge, edge) ;
    extract_quad_edge(second_request, second_edge, edge) ;
    extract_quad_edge(merged, merged_edge, edge) ;
    extract_quad_edge(reverse_merged, reverse_merged_edge, edge) ;

    CAPTURE(edge) ;
    if(!first_edge.empty() || !second_edge.empty()) {
      CHECK_FALSE(merged_edge.empty()) ;
    }
    CHECK(merged_edge == reverse_merged_edge) ;
  }
}


/// With matching one-level cell and face plans, each fine face must receive a
/// valid, distinct leaf index. The check intentionally ignores index ordering.
TEST_CASE("matching fine faces receive valid hex and prism leaf indices") {
  SUBCASE("HexCell boundary") {
    const std::vector<char> cell_plan = plan({7}) ;
    const std::vector<char> face_plan = plan({3}) ;
    const std::vector<int32> owners =
      get_c1_hex(cell_plan, face_plan, char(0), char(RIGHT)) ;

    HexCell counter ;
    const int fine_cell_count = counter.num_fine_cells(cell_plan) ;
    QuadFace face ;
    std::vector<QuadFace*> fine_faces ;
    face.empty_resplit(face_plan, char(0), fine_faces) ;

    REQUIRE(owners.size() == fine_faces.size()) ;
    std::set<int32> distinct_owners(owners.begin(), owners.end()) ;
    CHECK(distinct_owners.size() == owners.size()) ;
    for(size_t i = 0; i < owners.size(); ++i) {
      CAPTURE(i) ;
      CAPTURE(owners[i]) ;
      CHECK(owners[i] >= 1) ;
      CHECK(owners[i] <= fine_cell_count) ;
    }
  }

  SUBCASE("Prism quad boundary") {
    const std::vector<char> cell_plan = plan({3}) ;
    const std::vector<char> face_plan = plan({3}) ;
    const std::vector<int32> owners =
      get_c1_prism(cell_plan, face_plan, char(0), 2) ;

    Prism counter ;
    const int fine_cell_count = counter.empty_resplit(cell_plan) ;
    QuadFace face ;
    std::vector<QuadFace*> fine_faces ;
    face.empty_resplit(face_plan, char(0), fine_faces) ;

    REQUIRE(owners.size() == fine_faces.size()) ;
    std::set<int32> distinct_owners(owners.begin(), owners.end()) ;
    CHECK(distinct_owners.size() == owners.size()) ;
    for(size_t i = 0; i < owners.size(); ++i) {
      CAPTURE(i) ;
      CAPTURE(owners[i]) ;
      CHECK(owners[i] >= 1) ;
      CHECK(owners[i] <= fine_cell_count) ;
    }
  }
}


/// When the current mesh is refined but its parent plan is unsplit, traverse()
/// must map every current leaf back to the one coarse parent cell.
TEST_CASE("traverse maps refined leaves to an unsplit parent cell") {
  const std::vector<char> coarse_plan ;

  SUBCASE("HexCell") {
    HexFixture hex ;
    build_unit_hex(hex) ;
    hex.root->empty_resplit(coarse_plan) ;
    const int fine_cell_count = hex.root->empty_resplit(plan({7})) ;

    std::vector<std::pair<int32, int32> > index_map ;
    CHECK(hex.root->traverse(coarse_plan, index_map) == 1) ;
    REQUIRE(index_map.size() == size_t(fine_cell_count)) ;

    std::set<int32> fine_indices ;
    for(size_t i = 0; i < index_map.size(); ++i) {
      fine_indices.insert(index_map[i].first) ;
      CHECK(index_map[i].second == 1) ;
    }
    CHECK(fine_indices.size() == index_map.size()) ;
    for(int index = 1; index <= fine_cell_count; ++index) {
      CHECK(fine_indices.count(index) == 1) ;
    }
  }

  SUBCASE("Prism") {
    PrismFixture prism ;
    build_unit_prism(prism) ;
    prism.root->empty_resplit(coarse_plan) ;
    const int fine_cell_count = prism.root->empty_resplit(plan({3})) ;

    std::vector<std::pair<int32, int32> > index_map ;
    CHECK(prism.root->traverse(coarse_plan, index_map) == 1) ;
    REQUIRE(index_map.size() == size_t(fine_cell_count)) ;

    std::set<int32> fine_indices ;
    for(size_t i = 0; i < index_map.size(); ++i) {
      fine_indices.insert(index_map[i].first) ;
      CHECK(index_map[i].second == 1) ;
    }
    CHECK(fine_indices.size() == index_map.size()) ;
    for(int index = 1; index <= fine_cell_count; ++index) {
      CHECK(fine_indices.count(index) == 1) ;
    }
  }

  SUBCASE("general Cell") {
    GeneralFixture tetra ;
    build_tetra_cell(tetra) ;
    tetra.root->empty_resplit(coarse_plan) ;
    const int fine_cell_count = tetra.root->empty_resplit(plan({1})) ;

    std::vector<std::pair<int32, int32> > index_map ;
    CHECK(tetra.root->traverse(coarse_plan, index_map) == 1) ;
    REQUIRE(index_map.size() == size_t(fine_cell_count)) ;

    std::set<int32> fine_indices ;
    for(size_t i = 0; i < index_map.size(); ++i) {
      fine_indices.insert(index_map[i].first) ;
      CHECK(index_map[i].second == 1) ;
    }
    CHECK(fine_indices.size() == index_map.size()) ;
    for(int index = 1; index <= fine_cell_count; ++index) {
      CHECK(fine_indices.count(index) == 1) ;
    }
  }
}


/// A successful split must leave every resulting child non-null with a finite,
/// positive minimum edge length.
TEST_CASE("refined hex, prism, and general-cell children have nondegenerate edges") {
  SUBCASE("HexCell children") {
    HexFixture hex ;
    build_unit_hex(hex) ;
    std::vector<HexCell*> leaves ;
    hex.root->resplit(plan({7}), hex.nodes, hex.edges, hex.faces, leaves) ;

    REQUIRE_FALSE(leaves.empty()) ;
    for(size_t i = 0; i < leaves.size(); ++i) {
      CAPTURE(i) ;
      REQUIRE(leaves[i] != 0) ;
      const double minimum_edge = leaves[i]->get_min_edge_length() ;
      CAPTURE(minimum_edge) ;
      CHECK(std::isfinite(minimum_edge)) ;
      CHECK(minimum_edge > 0.0) ;
    }
  }

  SUBCASE("Prism children") {
    PrismFixture prism ;
    build_unit_prism(prism) ;
    std::vector<Prism*> leaves ;
    prism.root->resplit(plan({3}), prism.nodes, prism.edges, prism.qfaces,
                        prism.gfaces, leaves) ;

    REQUIRE_FALSE(leaves.empty()) ;
    for(size_t i = 0; i < leaves.size(); ++i) {
      CAPTURE(i) ;
      REQUIRE(leaves[i] != 0) ;
      const double minimum_edge = leaves[i]->get_min_edge_length() ;
      CAPTURE(minimum_edge) ;
      CHECK(std::isfinite(minimum_edge)) ;
      CHECK(minimum_edge > 0.0) ;
    }
  }

  SUBCASE("general Cell children") {
    GeneralFixture tetra ;
    build_tetra_cell(tetra) ;
    std::vector<DiamondCell*> leaves ;
    tetra.root->resplit(plan({1}), tetra.nodes, tetra.edges, tetra.faces,
                        leaves) ;

    REQUIRE_FALSE(leaves.empty()) ;
    for(size_t i = 0; i < leaves.size(); ++i) {
      CAPTURE(i) ;
      REQUIRE(leaves[i] != 0) ;
      const double minimum_edge = leaves[i]->get_min_edge_length() ;
      CAPTURE(minimum_edge) ;
      CHECK(std::isfinite(minimum_edge)) ;
      CHECK(minimum_edge > 0.0) ;
    }
  }
}


/// Requesting two refinement levels must refine every branch to the same depth.
/// This is checked for both HexCell and Prism without assuming child order.
TEST_CASE("level refinement reaches every hex and prism branch") {
  SUBCASE("HexCell") {
    HexFixture hex ;
    build_unit_hex(hex) ;
    const std::vector<char> expected_plan = hex.root->make_cellplan(2) ;

    hex.root->resplit(2, hex.nodes, hex.edges, hex.faces) ;
    CHECK(hex.root->make_cellplan() == expected_plan) ;
  }

  SUBCASE("Prism") {
    PrismFixture prism ;
    build_unit_prism(prism) ;
    const std::vector<char> expected_plan = prism.root->make_cellplan(2) ;

    prism.root->resplit(2, prism.nodes, prism.edges, prism.qfaces,
                        prism.gfaces) ;
    CHECK(prism.root->make_cellplan() == expected_plan) ;
  }
}
