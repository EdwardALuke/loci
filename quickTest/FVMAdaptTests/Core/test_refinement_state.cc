//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# This program is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//#############################################################################

#include "fvmadapt_core.h"

#include "refinement_state_internal.h"

#include <doctest.h>

#include <algorithm>
#include <initializer_list>

namespace {

std::vector<char> state_plan(std::initializer_list<int> codes) {
  std::vector<char> result ;
  for(std::initializer_list<int>::const_iterator code = codes.begin();
      code != codes.end(); ++code)
    result.push_back(char(*code)) ;
  return result ;
}

void check_depth_range(const std::vector<int>& depths,
                       int expectedLeaves,
                       int minimumDepth,
                       int maximumDepth) {
  REQUIRE(depths.size() == size_t(expectedLeaves)) ;
  REQUIRE_FALSE(depths.empty()) ;
  CHECK(*std::min_element(depths.begin(), depths.end()) == minimumDepth) ;
  CHECK(*std::max_element(depths.begin(), depths.end()) == maximumDepth) ;
}

} // namespace


/// Accepted plans expose one depth value for every final leaf, independent of
/// whether the leaf came from a general, hexahedral, or prism cell.
TEST_CASE("accepted cell plans expose final-leaf refinement depth") {
  SUBCASE("unsplit roots have depth zero") {
    HexCell hex ;
    Prism prism ;
    GeneralFixture general ;
    build_tetra_cell(general) ;
    std::vector<int> depths ;

    CHECK(hex.empty_resplit(std::vector<char>()) == 1) ;
    REQUIRE(Loci::detail::getLeafRefinementDepths(&hex, depths)) ;
    CHECK(depths == std::vector<int>(1, 0)) ;

    CHECK(prism.empty_resplit(std::vector<char>()) == 1) ;
    REQUIRE(Loci::detail::getLeafRefinementDepths(&prism, depths)) ;
    CHECK(depths == std::vector<int>(1, 0)) ;

    CHECK(general.root->empty_resplit(std::vector<char>()) == 1) ;
    REQUIRE(Loci::detail::getLeafRefinementDepths(general.root, depths)) ;
    CHECK(depths == std::vector<int>(1, 0)) ;
  }

  SUBCASE("one split puts every accepted leaf at depth one") {
    HexCell hex ;
    Prism prism ;
    GeneralFixture general ;
    build_tetra_cell(general) ;
    std::vector<int> depths ;

    const int hexLeaves = hex.empty_resplit(state_plan({7})) ;
    REQUIRE(Loci::detail::getLeafRefinementDepths(&hex, depths)) ;
    check_depth_range(depths, hexLeaves, 1, 1) ;

    const int prismLeaves = prism.empty_resplit(state_plan({3})) ;
    REQUIRE(Loci::detail::getLeafRefinementDepths(&prism, depths)) ;
    check_depth_range(depths, prismLeaves, 1, 1) ;

    const int generalLeaves = general.root->empty_resplit(state_plan({1})) ;
    REQUIRE(Loci::detail::getLeafRefinementDepths(general.root, depths)) ;
    check_depth_range(depths, generalLeaves, 1, 1) ;
  }

  SUBCASE("uneven plans preserve shallow and deeper accepted leaves") {
    HexCell hex ;
    Prism prism ;
    GeneralFixture general ;
    build_tetra_cell(general) ;
    std::vector<int> depths ;

    const int hexLeaves = hex.empty_resplit(state_plan({1, 1, 0})) ;
    REQUIRE(Loci::detail::getLeafRefinementDepths(&hex, depths)) ;
    check_depth_range(depths, hexLeaves, 1, 2) ;
    CHECK(depths == std::vector<int>({1, 2, 2})) ;

    const int prismLeaves = prism.empty_resplit(state_plan({1, 1, 0})) ;
    REQUIRE(Loci::detail::getLeafRefinementDepths(&prism, depths)) ;
    check_depth_range(depths, prismLeaves, 1, 2) ;
    CHECK(depths == std::vector<int>({1, 2, 2})) ;

    const int generalLeaves =
      general.root->empty_resplit(state_plan({1, 1, 0, 0, 0})) ;
    REQUIRE(Loci::detail::getLeafRefinementDepths(general.root, depths)) ;
    check_depth_range(depths, generalLeaves, 1, 2) ;
  }
}


/// Relation cardinality distinguishes retained, refined, and derefined output
/// cells without depending on topology-specific plan encodings.
TEST_CASE("old/new cell relations classify the latest adaptation result") {
  std::vector<int> result ;

  REQUIRE(Loci::detail::classifyAdaptResult(
    std::vector<std::pair<int32, int32> >({{1, 1}, {2, 2}}),
    1, 2, result)) ;
  CHECK(result == std::vector<int>({Loci::adapt_result::retained,
                                    Loci::adapt_result::retained})) ;

  REQUIRE(Loci::detail::classifyAdaptResult(
    std::vector<std::pair<int32, int32> >({{1, 1}, {2, 1}}),
    1, 2, result)) ;
  CHECK(result == std::vector<int>({Loci::adapt_result::refined,
                                    Loci::adapt_result::refined})) ;

  REQUIRE(Loci::detail::classifyAdaptResult(
    std::vector<std::pair<int32, int32> >({{1, 1}, {1, 2}}),
    1, 1, result)) ;
  CHECK(result == std::vector<int>(1, Loci::adapt_result::derefined)) ;

  REQUIRE(Loci::detail::classifyAdaptResult(
    std::vector<std::pair<int32, int32> >({{5, 8}, {6, 9}, {7, 9}}),
    5, 3, result)) ;
  CHECK(result == std::vector<int>({Loci::adapt_result::retained,
                                    Loci::adapt_result::refined,
                                    Loci::adapt_result::refined})) ;
}


/// Incomplete, duplicate, or many-to-many relations cannot describe one
/// unambiguous result for every generated cell and must be rejected.
TEST_CASE("invalid old/new cell relations are rejected") {
  std::vector<int> result ;

  CHECK_FALSE(Loci::detail::classifyAdaptResult(
    std::vector<std::pair<int32, int32> >({{1, 1}}), 1, 2, result)) ;
  CHECK_FALSE(Loci::detail::classifyAdaptResult(
    std::vector<std::pair<int32, int32> >({{1, 1}, {1, 1}}),
    1, 1, result)) ;
  CHECK_FALSE(Loci::detail::classifyAdaptResult(
    std::vector<std::pair<int32, int32> >(
      {{1, 1}, {1, 2}, {2, 1}, {2, 2}}),
    1, 2, result)) ;
}
