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

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <cmath>
#include <initializer_list>
#include <iostream>

namespace {

// Keeps test cases readable while plan storage remains char-based. For tests
// writing the split codes as integers is more readable, so this function
// just bridges that gap between readability and the actual library codes.
std::vector<char> plan(std::initializer_list<int> codes) {
  std::vector<char> result ;
  for(std::initializer_list<int>::const_iterator i = codes.begin() ;
      i != codes.end(); ++i) {
    result.push_back(char(*i)) ;
  }
  return result ;
}

/// Print the custom inventory used by the Makefile list target.
void print_test_inventory(std::ostream& out) {
  out << "FVMAdapt core checks\n" ;
  out << "  root split counts: HexCell, Prism, and general Cell\n" ;
  out << "  plan replay: omitted children, trailing zeros, and canonical plans\n" ;
  out << "  plan-transfer: transfer_plan_g2q/transfer_plan_q2g round trips\n" ;
  out << "  plan projection: cell-to-face and face-to-edge extraction\n" ;
  out << "  hex-cell plans: resplit, sort_leaves, num_fine_cells, make_cellplan\n" ;
  out << "  prism-cell plans: resplit, sort_leaves, empty_resplit, make_cellplan\n" ;
  out << "  general-cell plans: Cell/DiamondCell resplit, sort_leaves, "
      << "empty_resplit, make_cellplan\n" ;
  out << "  tree state: child tags, derefine checks, and derefine cleanup\n" ;
  out << "  replay traces: independent breadth-first leaf counts for all cases\n" ;
  out << "\n" ;
  out << "Example artifacts written by --write-dir\n" ;
  out << "  VTK/wireframe files for selected plans\n" ;
  out << "  core_split_summary.dat\n" ;
  out << "  core_plan_audit.dat\n" ;
  out << "  core_plan_replay.dat\n" ;
}

/// Run the checks and then write the optional example artifacts.
void write_examples(const std::string& output_dir) {
  make_directory(output_dir) ;

  const std::vector<CaseResult> results =
    run_cell_plan_tests(true, output_dir) ;

  check_plan_replay_counts(results) ;
  write_example_artifacts(output_dir, results) ;
}

/// Check the independent replay trace against the cell-plan results.
void check_replay_traces() {
  const std::vector<CaseResult> results = run_cell_plan_tests(false, "") ;
  check_plan_replay_counts(results) ;
}

} // namespace

/// Checks that a four-sided original mesh face keeps the same refined pattern
/// when its plan is translated between general Face and QuadFace formats.
/// Mixed hex/general interfaces use this so both cells split their copy of the
/// same original face consistently.
TEST_CASE("FVMAdapt face plans round-trip between general and quad encodings") {
  CHECK_NOTHROW(check_plan_transfer_round_trips()) ;
}

/// A one-entry cell plan says how the parent cell itself is split. With no
/// child entries following it, every child is a leaf, so the library should
/// report exactly the child count implied by that root split code.
TEST_CASE("root split codes produce the expected fine-cell counts") {
  // Hex root split codes produce the fixed HexCell leaf-count table.
  const int hex_counts[8] = {1, 2, 2, 4, 2, 4, 4, 8} ;
  for(int code = 0; code < 8; ++code) {
    CAPTURE(code) ;
    HexCell counter ;
    CHECK(counter.num_fine_cells(plan({code})) == hex_counts[code]) ;
  }

  // Prism root split codes produce the fixed Prism leaf-count table.
  const int prism_counts[4] = {1, 2, 3, 6} ;
  for(int code = 0; code < 4; ++code) {
    CAPTURE(code) ;
    Prism counter ;
    CHECK(counter.empty_resplit(plan({code})) == prism_counts[code]) ;
  }

  // A tetrahedron root split creates one DiamondCell child per original node.
  GeneralFixture tetra ;
  build_tetra_cell(tetra) ;
  CHECK(tetra.root->empty_resplit(plan({1})) == tetra.root->numNode) ;

  // A pyramid root split follows the same general-cell child-count rule.
  GeneralFixture pyramid ;
  build_pyramid_cell(pyramid) ;
  CHECK(pyramid.root->empty_resplit(plan({1})) == pyramid.root->numNode) ;
}

/// Checks the compact cell-plan convention for omitted child entries. When a
/// plan ends before every child has an explicit code, the missing children are
/// unsplit leaves, and make_cellplan() should preserve that compact form.
TEST_CASE("short cell plans treat omitted child entries as leaves") {
  // Only one hex child is split explicitly; omitted siblings remain leaves.
  HexFixture hex ;
  build_unit_hex(hex) ;
  std::vector<HexCell*> hex_leaves ;
  hex.root->resplit(plan({7, 1}), hex.nodes, hex.edges, hex.faces, hex_leaves) ;
  CHECK(hex_leaves.size() == 9) ;
  CHECK(hex.root->make_cellplan() == plan({7, 1})) ;

  // Prism plans use the same breadth-first convention for omitted children.
  PrismFixture prism ;
  build_unit_prism(prism) ;
  std::vector<Prism*> prism_leaves ;
  prism.root->resplit(plan({3, 1}), prism.nodes, prism.edges, prism.qfaces,
                      prism.gfaces, prism_leaves) ;
  CHECK(prism_leaves.size() == 7) ;
  CHECK(prism.root->make_cellplan() == plan({3, 1})) ;

  // General-cell plans also allow the trailing child entries to be implicit.
  GeneralFixture tetra ;
  build_tetra_cell(tetra) ;
  std::vector<DiamondCell*> tetra_leaves ;
  tetra.root->resplit(plan({1, 1}), tetra.nodes, tetra.edges, tetra.faces,
                      tetra_leaves) ;
  CHECK(tetra_leaves.size() > size_t(tetra.root->numNode)) ;
  CHECK(tetra.root->make_cellplan() == plan({1, 1})) ;
}

/// Checks that explicit trailing leaf markers do not become part of the saved
/// cell plan. A producer may spell out every child of a split parent as `0`,
/// but make_cellplan() should serialize the same compact tree as if those leaf
/// children had been omitted. This keeps equivalent plans stable across
/// HexCell, Prism, and general Cell replay paths.
TEST_CASE("trailing zero child plans canonicalize to the split tree") {
  // A fully split hex with all eight children marked as leaves saves as `{7}`.
  HexFixture hex ;
  build_unit_hex(hex) ;
  std::vector<HexCell*> hex_leaves ;
  hex.root->resplit(plan({7, 0, 0, 0, 0, 0, 0, 0, 0}),
                    hex.nodes,
                    hex.edges,
                    hex.faces,
                    hex_leaves) ;
  CHECK(hex.root->make_cellplan() == plan({7})) ;

  // Prism full-split plans trim the same explicit child-leaf markers.
  PrismFixture prism ;
  build_unit_prism(prism) ;
  std::vector<Prism*> prism_leaves ;
  prism.root->resplit(plan({3, 0, 0, 0, 0, 0, 0}),
                      prism.nodes,
                      prism.edges,
                      prism.qfaces,
                      prism.gfaces,
                      prism_leaves) ;
  CHECK(prism.root->make_cellplan() == plan({3})) ;

  // General-cell plans also save only the parent split once children are leaves.
  GeneralFixture pyramid ;
  build_pyramid_cell(pyramid) ;
  std::vector<DiamondCell*> pyramid_leaves ;
  pyramid.root->resplit(plan({1, 0, 0, 0, 0, 0}),
                        pyramid.nodes,
                        pyramid.edges,
                        pyramid.faces,
                        pyramid_leaves) ;
  CHECK(pyramid.root->make_cellplan() == plan({1})) ;
}

TEST_CASE("general cells keep the empty plan as the unsplit-cell convention") {
  // resplit(empty) leaves the output list empty for an unsplit general cell.
  GeneralFixture tetra_for_resplit ;
  build_tetra_cell(tetra_for_resplit) ;
  std::vector<DiamondCell*> leaves ;
  tetra_for_resplit.root->resplit(std::vector<char>(),
                                  tetra_for_resplit.nodes,
                                  tetra_for_resplit.edges,
                                  tetra_for_resplit.faces,
                                  leaves) ;
  CHECK(leaves.empty()) ;

  // empty_resplit(empty) still counts the parent as one unsplit cell.
  GeneralFixture tetra_for_count ;
  build_tetra_cell(tetra_for_count) ;
  CHECK(tetra_for_count.root->empty_resplit(std::vector<char>()) == 1) ;
  // make_cellplan() preserves the empty-plan representation.
  CHECK(tetra_for_count.root->make_cellplan().empty()) ;
}

/// Checks the helpers that derive lower-dimensional plans from a cell or face
/// plan. Face-balancing code uses these projections to decide how much of a
/// cell split is visible on one shared parent face, and edge-balancing code
/// does the same for one edge of a quad face. These checks anchor a few small
/// mappings for hex faces, prism triangular/quad faces, and quad-face edges.
TEST_CASE("cell and face plans project to representative boundary plans") {
  // An unsplit hex cell induces no split on any selected face.
  CHECK(extract_hex_face(std::vector<char>(), RIGHT).empty()) ;

  // A full hex split induces a full quad-face split on the right face.
  CHECK(extract_hex_face(plan({7}), RIGHT) == plan({3})) ;

  // A one-direction hex split is visible on a side face it cuts.
  CHECK(extract_hex_face(plan({1}), RIGHT) == plan({1})) ;

  // The same one-direction split does not cut the top face.
  CHECK(extract_hex_face(plan({1}), UP).empty()) ;

  // An unsplit prism cell induces no split on the selected triangular face.
  CHECK(extract_prism_face(std::vector<char>(), 0).empty()) ;

  // A full prism split induces a one-level split on a triangular end face.
  CHECK(extract_prism_face(plan({3}), 0) == plan({1})) ;

  // The same full prism split induces a full split on a quad side face.
  CHECK(extract_prism_face(plan({3}), 2) == plan({3})) ;

  // Prism split code 1 leaves the triangular end face unsplit.
  CHECK(extract_prism_face(plan({1}), 0).empty()) ;

  // Prism split code 1 is visible on the selected quad side face.
  CHECK(extract_prism_face(plan({1}), 2) == plan({1})) ;

  // A full quad-face split cuts the selected edge.
  std::vector<char> edge_plan ;
  extract_quad_edge(plan({3}), edge_plan, 0) ;
  CHECK(edge_plan == plan({1})) ;

  // Quad-face split code 1 does not cut edge 0.
  edge_plan.clear() ;
  extract_quad_edge(plan({1}), edge_plan, 0) ;
  CHECK(edge_plan.empty()) ;

  // The same quad-face split code does cut edge 1.
  edge_plan.clear() ;
  extract_quad_edge(plan({1}), edge_plan, 1) ;
  CHECK(edge_plan == plan({1})) ;
}

/// Checks the smallest fixed mapping between the general Face plan encoding
/// and the QuadFace plan encoding. Mixed general/hex and general/prism paths
/// convert a four-sided general face into QuadFace form before merging face
/// plans. In that conversion, a one-level general Face split uses code `1`,
/// while the equivalent QuadFace split uses code `3`.
TEST_CASE("small face-plan transfers have fixed general and quad encodings") {
  std::vector<char> general_split = plan({1}) ;

  // General Face code 1 becomes QuadFace code 3 for the same four-way split.
  CHECK(transfer_plan_g2q(general_split) == plan({3})) ;

  // QuadFace code 3 converts back to general Face code 1.
  CHECK(transfer_plan_q2g(plan({3})) == plan({1})) ;
}

/// Checks the direct-tag derefinement predicate for split cells. A parent is
/// ready for `needDerefine_ctag()` only when its immediate children are tagged
/// with `2` and those children are leaves. If that test passes, `derefine()`
/// should collapse the parent back to an unsplit cell plan.
TEST_CASE("derefine checks require tagged immediate leaf children") {
  // A split hex with every immediate child tagged 2 is ready to derefine.
  HexFixture hex ;
  build_unit_hex(hex) ;
  std::vector<HexCell*> hex_leaves ;
  hex.root->resplit(plan({7}), hex.nodes, hex.edges, hex.faces, hex_leaves) ;
  for(int child = 0; child < hex.root->numChildren(); ++child) {
    hex.root->getChildCell(child)->setTag(char(2)) ;
  }
  CHECK(hex.root->needDerefine_ctag()) ;

  // derefine() removes the root split and returns the plan to unsplit.
  hex.root->derefine() ;
  CHECK(hex.root->getMySplitCode() == 0) ;
  CHECK(hex.root->make_cellplan().empty()) ;

  // Tagged immediate children are not enough when one child is still split.
  HexFixture nested_hex ;
  build_unit_hex(nested_hex) ;
  std::vector<HexCell*> nested_hex_leaves ;
  nested_hex.root->resplit(plan({7, 1}),
                           nested_hex.nodes,
                           nested_hex.edges,
                           nested_hex.faces,
                           nested_hex_leaves) ;
  for(int child = 0; child < nested_hex.root->numChildren(); ++child) {
    nested_hex.root->getChildCell(child)->setTag(char(2)) ;
  }
  CHECK_FALSE(nested_hex.root->needDerefine_ctag()) ;

  // Prism derefine readiness uses the same direct-tag and leaf-child rule.
  PrismFixture prism ;
  build_unit_prism(prism) ;
  std::vector<Prism*> prism_leaves ;
  prism.root->resplit(plan({3}),
                      prism.nodes,
                      prism.edges,
                      prism.qfaces,
                      prism.gfaces,
                      prism_leaves) ;
  for(int child = 0; child < prism.root->numChildren(); ++child) {
    prism.root->getChildCell(child)->setTag(char(2)) ;
  }
  CHECK(prism.root->needDerefine_ctag()) ;
}

/// Smoke-tests that representative split operations create usable child cells,
/// not just the right plan structure. The checks verify that the split returns
/// the expected number of non-null leaf pointers and that the resulting parent
/// still reports a finite, positive minimum edge length.
TEST_CASE("representative splits leave finite positive geometry") {
  // A full hex split produces eight non-null leaves.
  HexFixture hex ;
  build_unit_hex(hex) ;
  std::vector<HexCell*> hex_leaves ;
  hex.root->resplit(plan({7}), hex.nodes, hex.edges, hex.faces, hex_leaves) ;
  REQUIRE(hex_leaves.size() == 8) ;
  for(size_t i = 0; i < hex_leaves.size(); ++i) {
    CAPTURE(i) ;
    CHECK(hex_leaves[i] != 0) ;
  }
  // Split hex geometry still has a finite, positive minimum edge length.
  const double hex_min_edge = hex.root->get_min_edge_length() ;
  CHECK(std::isfinite(hex_min_edge)) ;
  CHECK(hex_min_edge > 0.0) ;

  // A general-cell split produces one non-null child per original node.
  GeneralFixture tetra ;
  build_tetra_cell(tetra) ;
  std::vector<DiamondCell*> tetra_leaves ;
  tetra.root->resplit(plan({1}),
                      tetra.nodes,
                      tetra.edges,
                      tetra.faces,
                      tetra_leaves) ;
  REQUIRE(tetra_leaves.size() == size_t(tetra.root->numNode)) ;
  for(size_t i = 0; i < tetra_leaves.size(); ++i) {
    CAPTURE(i) ;
    CHECK(tetra_leaves[i] != 0) ;
  }
  // Split general-cell geometry also keeps a finite, positive minimum edge.
  const double tetra_min_edge = tetra.root->get_min_edge_length() ;
  CHECK(std::isfinite(tetra_min_edge)) ;
  CHECK(tetra_min_edge > 0.0) ;
}

/// Runs the representative HexCell plan table through the real HexCell replay
/// paths. Each case must give the same leaf count from resplit(), sort_leaves(),
/// and num_fine_cells(), and make_cellplan() must return the expected compact
/// canonical plan.
TEST_CASE("HexCell plans preserve resplit leaves, sorted leaves, and canonical plans") {
  CHECK_NOTHROW(run_hex_cell_plan_tests(false, "")) ;
}

/// Runs the representative Prism plan table through the real Prism replay
/// paths. Each case must give the same leaf count from resplit(), sort_leaves(),
/// and empty_resplit(), and make_cellplan() must return the expected compact
/// canonical plan.
TEST_CASE("Prism plans preserve resplit leaves, sorted leaves, and canonical plans") {
  CHECK_NOTHROW(run_prism_cell_plan_tests(false, "")) ;
}

/// Runs the representative general-cell plan table on a tetrahedron. This
/// checks that Cell/DiamondCell resplit(), sort_leaves(), empty_resplit(), and
/// make_cellplan() agree on the leaf count and canonical plan.
TEST_CASE("general tetrahedron plans preserve resplit leaves, sorted leaves, and canonical plans") {
  CHECK_NOTHROW(run_tetra_cell_plan_tests(false, "")) ;
}

/// Runs the same general-cell plan table on a pyramid. This gives the generic
/// Cell/DiamondCell replay checks a second parent-cell topology.
TEST_CASE("general pyramid plans preserve resplit leaves, sorted leaves, and canonical plans") {
  CHECK_NOTHROW(run_pyramid_cell_plan_tests(false, "")) ;
}

/// Checks the independent breadth-first replay tracer against the real cell
/// plan results. The tracer does not call HexCell, Prism, or Cell replay; it
/// reads the plan vector with the same child-count rules and verifies that its
/// leaf count matches the library result for every representative plan case.
TEST_CASE("plan replay traces reproduce the counted leaves for every cell case") {
  CHECK_NOTHROW(check_replay_traces()) ;
}

int main(int argc, char** argv) {
  bool write_example_files = false ;
  bool list_tests = false ;
  std::string output_dir = "output" ;
  std::vector<char*> doctest_args ;
  doctest_args.push_back(argv[0]) ;

  for(int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]) ;
    if(arg == "--check-only") {
      continue ;
    } else if(arg == "--list") {
      list_tests = true ;
    } else if(arg == "--write-dir") {
      if(i + 1 >= argc) {
        std::cerr << "--write-dir requires a path\n" ;
        std::cout << "FAILURE!\n" ;
        return 1 ;
      }
      write_example_files = true ;
      output_dir = argv[++i] ;
    } else {
      doctest_args.push_back(argv[i]) ;
    }
  }

  if(list_tests) {
    print_test_inventory(std::cout) ;
    return 0 ;
  }

  doctest::Context context ;
  context.applyCommandLine(int(doctest_args.size()), doctest_args.data()) ;
  const int result = context.run() ;

  if(context.shouldExit() || result != 0) {
    return result ;
  }

  if(!write_example_files) {
    return result ;
  }

  try {
    write_examples(output_dir) ;
    std::cout << "SUCCESS!\n" ;
    return result ;
  } catch(const std::exception& error) {
    std::cerr << error.what() << "\n" ;
    std::cout << "FAILURE!\n" ;
    return 1 ;
  }
}
