//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
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

#include <Loci.h>
#include <RStore.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <vector>
#include <sstream>

using namespace Loci ;
using namespace std ;

// =========================================================================
// rstore tests
// =========================================================================

TEST_CASE("rstore: bulk_insert and const access") {
  rstore<double> rs ;
  vector<pair<Entity,double>> entries ;
  entries.push_back({10, 1.0}) ;
  entries.push_back({5, 2.0}) ;
  entries.push_back({20, 3.0}) ;
  entries.push_back({15, 4.0}) ;
  rs.bulk_insert(entries) ;

  CHECK(rs.domain().size() == 4) ;
  CHECK(rs[5] == 2.0) ;
  CHECK(rs[10] == 1.0) ;
  CHECK(rs[15] == 4.0) ;
  CHECK(rs[20] == 3.0) ;
  CHECK(rs.elem(5) == 2.0) ;
}

TEST_CASE("rstore: bulk_insert overwrites existing entries") {
  rstore<int> rs ;
  vector<pair<Entity,int>> entries = {{1, 100}, {2, 200}, {3, 300}} ;
  rs.bulk_insert(entries) ;
  CHECK(rs[2] == 200) ;

  vector<pair<Entity,int>> updates = {{2, 999}, {4, 400}} ;
  rs.bulk_insert(updates) ;

  CHECK(rs.domain().size() == 4) ;
  CHECK(rs[1] == 100) ;
  CHECK(rs[2] == 999) ;
  CHECK(rs[3] == 300) ;
  CHECK(rs[4] == 400) ;
}

TEST_CASE("rstore: bulk_erase") {
  rstore<double> rs ;
  vector<pair<Entity,double>> entries = {
    {1, 10.0}, {2, 20.0}, {3, 30.0}, {4, 40.0}, {5, 50.0}
  } ;
  rs.bulk_insert(entries) ;

  entitySet to_erase ;
  to_erase += 2 ;
  to_erase += 4 ;
  rs.bulk_erase(to_erase) ;

  CHECK(rs.domain().size() == 3) ;
  CHECK(rs[1] == 10.0) ;
  CHECK(rs[3] == 30.0) ;
  CHECK(rs[5] == 50.0) ;
}

TEST_CASE("rstore: allocate sets domain") {
  rstore<int> rs ;
  entitySet ptn ;
  ptn += interval(10, 14) ;
  rs.allocate(ptn) ;

  CHECK(rs.domain() == ptn) ;
  CHECK(rs[10] == 0) ;
  CHECK(rs[14] == 0) ;

  // Reallocate to a different domain
  entitySet ptn2 ;
  ptn2 += interval(12, 16) ;
  rs.allocate(ptn2) ;

  CHECK(rs.domain() == ptn2) ;
}

TEST_CASE("rstore: const_rstore provides read-only access") {
  rstore<double> rs ;
  vector<pair<Entity,double>> entries = {{1, 1.5}, {2, 2.5}, {3, 3.5}} ;
  rs.bulk_insert(entries) ;

  storeRepP rep = rs.Rep() ;
  const_rstore<double> crs(rep) ;

  CHECK(crs[1] == 1.5) ;
  CHECK(crs[2] == 2.5) ;
  CHECK(crs.elem(3) == 3.5) ;
  CHECK(crs.domain().size() == 3) ;
}

TEST_CASE("rstore: Print and Input round-trip") {
  rstore<int> rs ;
  vector<pair<Entity,int>> entries = {{5, 50}, {10, 100}, {15, 150}} ;
  rs.bulk_insert(entries) ;

  stringstream ss ;
  rs.Print(ss) ;

  rstore<int> rs2 ;
  rs2.Input(ss) ;

  CHECK(rs2.domain() == rs.domain()) ;
  CHECK(rs2[5] == 50) ;
  CHECK(rs2[10] == 100) ;
  CHECK(rs2[15] == 150) ;
}

TEST_CASE("rstore: empty store") {
  rstore<double> rs ;
  CHECK(rs.domain() == EMPTY) ;

  entitySet empty_erase ;
  rs.bulk_erase(empty_erase) ;
  CHECK(rs.domain() == EMPTY) ;

  vector<pair<Entity,double>> empty_insert ;
  rs.bulk_insert(empty_insert) ;
  CHECK(rs.domain() == EMPTY) ;
}

// =========================================================================
// rmultiStore tests
// =========================================================================

TEST_CASE("rmultiStore: basic bulk_insert and access") {
  rmultiStore<int> rms ;
  vector<pair<Entity,int>> entries = {
    {1, 10}, {1, 20}, {1, 30},
    {2, 40},
    {3, 50}, {3, 60}
  } ;
  rms.bulk_insert(entries) ;

  CHECK(rms.domain().size() == 3) ;
  CHECK(rms.num_elems(1) == 3) ;
  CHECK(rms.num_elems(2) == 1) ;
  CHECK(rms.num_elems(3) == 2) ;

  rmultiStore_ref<int> v1 = rms[1] ;
  CHECK(v1.size() == 3) ;
  CHECK(v1[0] == 10) ;
  CHECK(v1[1] == 20) ;
  CHECK(v1[2] == 30) ;

  rmultiStore_ref<int> v2 = rms[2] ;
  CHECK(v2.size() == 1) ;
  CHECK(v2[0] == 40) ;
}

TEST_CASE("rmultiStore: stable sort preserves insertion order") {
  rmultiStore<int> rms ;
  // Insert values for entity 5 in a specific order across two batches
  vector<pair<Entity,int>> batch1 = {
    {5, 100}, {5, 200}, {5, 300}
  } ;
  rms.bulk_insert(batch1) ;

  vector<pair<Entity,int>> batch2 = {
    {5, 400}, {5, 500}
  } ;
  rms.bulk_insert(batch2) ;

  CHECK(rms.num_elems(5) == 5) ;
  rmultiStore_ref<int> vals = rms[5] ;
  // batch1 entries come first, then batch2 entries
  CHECK(vals[0] == 100) ;
  CHECK(vals[1] == 200) ;
  CHECK(vals[2] == 300) ;
  CHECK(vals[3] == 400) ;
  CHECK(vals[4] == 500) ;
}

TEST_CASE("rmultiStore: allocate with zero values") {
  rmultiStore<double> rms ;
  entitySet ptn ;
  ptn += 1 ;
  ptn += 2 ;
  ptn += 3 ;
  rms.allocate(ptn) ;

  // All three entities are in the domain with zero values
  CHECK(rms.domain() == ptn) ;
  CHECK(rms.num_elems(1) == 0) ;
  CHECK(rms.num_elems(2) == 0) ;
  CHECK(rms.num_elems(3) == 0) ;

  // Now insert some values for entity 2 only
  vector<pair<Entity,double>> entries = {{2, 1.5}, {2, 2.5}} ;
  rms.bulk_insert(entries) ;

  // Entity 1 and 3 still have zero values but remain in domain
  CHECK(rms.domain() == ptn) ;
  CHECK(rms.num_elems(1) == 0) ;
  CHECK(rms.num_elems(2) == 2) ;
  CHECK(rms.num_elems(3) == 0) ;
}

TEST_CASE("rmultiStore: bulk_erase removes entities and values") {
  rmultiStore<int> rms ;
  vector<pair<Entity,int>> entries = {
    {1, 10}, {1, 20},
    {2, 30},
    {3, 40}, {3, 50}, {3, 60}
  } ;
  rms.bulk_insert(entries) ;

  entitySet to_erase ;
  to_erase += 1 ;
  to_erase += 3 ;
  rms.bulk_erase(to_erase) ;

  CHECK(rms.domain().size() == 1) ;
  CHECK(rms.num_elems(2) == 1) ;
  CHECK(rms[2][0] == 30) ;
  CHECK(rms.num_elems(1) == 0) ;
  CHECK(rms.num_elems(3) == 0) ;
}

TEST_CASE("rmultiStore: allocate shrinks domain") {
  rmultiStore<int> rms ;
  vector<pair<Entity,int>> entries = {
    {1, 10}, {2, 20}, {3, 30}
  } ;
  rms.bulk_insert(entries) ;

  entitySet new_domain ;
  new_domain += 2 ;
  rms.allocate(new_domain) ;

  CHECK(rms.domain() == new_domain) ;
  CHECK(rms.num_elems(2) == 1) ;
  CHECK(rms[2][0] == 20) ;
}

TEST_CASE("rmultiStore: const_rmultiStore read-only access") {
  rmultiStore<double> rms ;
  vector<pair<Entity,double>> entries = {
    {1, 1.5}, {1, 2.5},
    {2, 3.5}
  } ;
  rms.bulk_insert(entries) ;

  storeRepP rep = rms.Rep() ;
  const_rmultiStore<double> crms(rep) ;

  CHECK(crms.num_elems(1) == 2) ;
  CHECK(crms[1][0] == 1.5) ;
  CHECK(crms[1][1] == 2.5) ;
  CHECK(crms[2][0] == 3.5) ;
}

TEST_CASE("rmultiStore: range-based for over ref") {
  rmultiStore<int> rms ;
  vector<pair<Entity,int>> entries = {
    {1, 10}, {1, 20}, {1, 30}
  } ;
  rms.bulk_insert(entries) ;

  vector<int> collected ;
  for(const auto &val : rms[1])
    collected.push_back(val) ;

  REQUIRE(collected.size() == 3) ;
  CHECK(collected[0] == 10) ;
  CHECK(collected[1] == 20) ;
  CHECK(collected[2] == 30) ;
}

TEST_CASE("rmultiStore: Print and Input round-trip") {
  rmultiStore<int> rms ;
  entitySet ptn ;
  ptn += 5 ;
  ptn += 10 ;
  ptn += 15 ;
  rms.allocate(ptn) ;

  vector<pair<Entity,int>> entries = {
    {5, 50}, {5, 55},
    {15, 150}
  } ;
  rms.bulk_insert(entries) ;

  stringstream ss ;
  rms.Print(ss) ;

  rmultiStore<int> rms2 ;
  rms2.Input(ss) ;

  CHECK(rms2.domain() == rms.domain()) ;
  CHECK(rms2.num_elems(5) == 2) ;
  CHECK(rms2[5][0] == 50) ;
  CHECK(rms2[5][1] == 55) ;
  CHECK(rms2.num_elems(10) == 0) ;
  CHECK(rms2.num_elems(15) == 1) ;
  CHECK(rms2[15][0] == 150) ;
}

TEST_CASE("rmultiStore: empty store") {
  rmultiStore<double> rms ;
  CHECK(rms.domain() == EMPTY) ;
  CHECK(rms.num_elems(1) == 0) ;

  rmultiStore_ref<double> vals = rms[1] ;
  CHECK(vals.size() == 0) ;
}

TEST_CASE("rmultiStore: mixed entities with varying value counts") {
  rmultiStore<int> rms ;
  entitySet ptn ;
  ptn += interval(1, 5) ;
  rms.allocate(ptn) ;

  vector<pair<Entity,int>> entries ;
  // Entity 1: 0 values (already in domain via allocate)
  // Entity 2: 1 value
  entries.push_back({2, 20}) ;
  // Entity 3: 3 values
  entries.push_back({3, 30}) ;
  entries.push_back({3, 31}) ;
  entries.push_back({3, 32}) ;
  // Entity 4: 0 values
  // Entity 5: 2 values
  entries.push_back({5, 50}) ;
  entries.push_back({5, 51}) ;

  rms.bulk_insert(entries) ;

  CHECK(rms.domain() == ptn) ;
  CHECK(rms.num_elems(1) == 0) ;
  CHECK(rms.num_elems(2) == 1) ;
  CHECK(rms.num_elems(3) == 3) ;
  CHECK(rms.num_elems(4) == 0) ;
  CHECK(rms.num_elems(5) == 2) ;
  CHECK(rms[3][0] == 30) ;
  CHECK(rms[3][1] == 31) ;
  CHECK(rms[3][2] == 32) ;
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv) ;

  doctest::Context context ;
  context.applyCommandLine(argc, argv) ;
  const int result = context.run() ;

  Loci::Finalize() ;
  return result ;
}
