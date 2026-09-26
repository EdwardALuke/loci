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

#include <set>
#include <utility>
#include <vector>

#include <Loci.h>
#include <LociGridReaders.h>
#include <FVMAdapt/defines.h>
#include <FVMAdapt/gridInterface.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

namespace {
  using CellPair = std::pair<int32, int32>; // (new cell, old cell), local indices

  // Use the real mesh setup and mapping rules, not a reconstructed rule body.
  void check_mapping(const std::vector<char>& old_plan,
                     const std::vector<char>& new_plan,
                     int old_count, int new_count,
                     const std::set<CellPair>& expected) {
    fact_db facts;
    REQUIRE(Loci::setupFVMGrid(facts, "tet.vog"));
    Loci::createLowerUpper(facts);
    Loci::createEdgesPar(facts);
    Loci::parallelClassifyCell(facts);

    constraint cells = facts.get_variable("gnrlcells");
    const int local_roots = (*cells).size();
    int roots = 0;
    MPI_Allreduce(&local_roots, &roots, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    REQUIRE(roots == 1);

    store<std::vector<char>> parent_plan, balanced_plan;
    parent_plan.allocate(*cells);
    balanced_plan.allocate(*cells);
    for(Entity cell : *cells) {
      parent_plan[cell] = old_plan;
      balanced_plan[cell] = new_plan;
    }
    facts.create_fact("parentPlan", parent_plan);
    facts.create_fact("balancedCellPlan", balanced_plan);

    rule_db rules;
    Loci::load_module("fvmadapt", rules);
    REQUIRE(Loci::makeQuery(rules, facts,
      "indexMap,parent_num_fine_cells,balanced_num_fine_cells"));

    const_store<int> old_cells = facts.get_variable("parent_num_fine_cells");
    const_store<int> new_cells = facts.get_variable("balanced_num_fine_cells");
    const_store<std::vector<CellPair>> mapping = facts.get_variable("indexMap");
    for(Entity cell : *cells) {
      CHECK(old_cells[cell] == old_count);
      CHECK(new_cells[cell] == new_count);
      CHECK(mapping[cell].size() == expected.size());
      const std::set<CellPair> actual(mapping[cell].begin(), mapping[cell].end());
      for(const CellPair& relation : expected)
        CHECK_MESSAGE(actual.count(relation) == 1,
                      "Missing mapping: new cell ", relation.first,
                      ", old cell ", relation.second);
    }
  }
}

/// Collapsing a refined tetrahedron must map every old cell to the one new
/// root, including when one of its four children was refined again.
TEST_CASE("general-cell collapse maps every old cell to the root") {
  for(bool nested : {false, true}) {
    CAPTURE(nested);
    // One root split gives four cells; splitting one child gives eleven.
    const std::vector<char> old_plan = nested ? std::vector<char>{1, 1}
                                             : std::vector<char>{1};
    const int old_count = nested ? 11 : 4;
    std::set<CellPair> expected;
    for(int old_cell = 1; old_cell <= old_count; ++old_cell)
      expected.emplace(1, old_cell);
    check_mapping(old_plan, {}, old_count, 1, expected);
  }
}

/// Retention, refinement, and partial derefinement must keep their existing
/// relations; only a complete collapse sends all old cells to the root.
TEST_CASE("general-cell mappings retain refinement and partial derefinement") {
  SUBCASE("unchanged root") {
    check_mapping({}, {}, 1, 1, {{1, 1}});
  }
  SUBCASE("refine the root") {
    check_mapping({}, {1}, 1, 4, {{1, 1}, {2, 1}, {3, 1}, {4, 1}});
  }
  SUBCASE("unchanged refined cells") {
    check_mapping({1}, {1}, 4, 4, {{1, 1}, {2, 2}, {3, 3}, {4, 4}});
  }
  SUBCASE("collapse only the refined child") {
    // Breadth-first numbering puts the three unchanged old cells first.
    std::set<CellPair> expected = {{2, 1}, {3, 2}, {4, 3}};
    for(int old_cell = 4; old_cell <= 11; ++old_cell)
      expected.emplace(1, old_cell);
    check_mapping({1, 1}, {1}, 11, 4, expected);
  }
}

int main(int argc, char** argv) {
  Loci::Init(&argc, &argv);
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  context.setOption("abort-after", 0);
  const int local_failure = context.run() == 0 ? 0 : 1;
  int failure = 0;
  MPI_Allreduce(&local_failure, &failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  Loci::Finalize();
  return failure;
}
