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
#include <RMap.h>
#include <RMultiMap.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <vector>
#include <sstream>

using namespace Loci ;
using namespace std ;

// =========================================================================
// rMap tests
// =========================================================================

TEST_CASE("rMap: bulk_insert and const access") {
  rMap rm ;
  vector<pair<Entity,Entity>> entries ;
  entries.push_back({10, 100}) ;
  entries.push_back({5, 50}) ;
  entries.push_back({20, 200}) ;
  entries.push_back({15, 150}) ;
  rm.bulk_insert(entries) ;

  CHECK(rm.domain().size() == 4) ;
  CHECK(rm[5] == 50) ;
  CHECK(rm[10] == 100) ;
  CHECK(rm[15] == 150) ;
  CHECK(rm[20] == 200) ;
}

TEST_CASE("rMap: bulk_insert overwrites") {
  rMap rm ;
  vector<pair<Entity,Entity>> entries = {{1, 10}, {2, 20}, {3, 30}} ;
  rm.bulk_insert(entries) ;

  vector<pair<Entity,Entity>> updates = {{2, 99}, {4, 40}} ;
  rm.bulk_insert(updates) ;

  CHECK(rm.domain().size() == 4) ;
  CHECK(rm[1] == 10) ;
  CHECK(rm[2] == 99) ;
  CHECK(rm[3] == 30) ;
  CHECK(rm[4] == 40) ;
}

TEST_CASE("rMap: bulk_erase") {
  rMap rm ;
  vector<pair<Entity,Entity>> entries = {
    {1, 10}, {2, 20}, {3, 30}, {4, 40}, {5, 50}
  } ;
  rm.bulk_insert(entries) ;

  entitySet to_erase ;
  to_erase += 2 ;
  to_erase += 4 ;
  rm.bulk_erase(to_erase) ;

  CHECK(rm.domain().size() == 3) ;
  CHECK(rm[1] == 10) ;
  CHECK(rm[3] == 30) ;
  CHECK(rm[5] == 50) ;
}

TEST_CASE("rMap: image") {
  rMap rm ;
  vector<pair<Entity,Entity>> entries = {
    {1, 100}, {2, 200}, {3, 100}, {4, 300}
  } ;
  rm.bulk_insert(entries) ;

  entitySet query ;
  query += 1 ;
  query += 3 ;
  entitySet img = rm.image(query) ;

  CHECK(img.size() == 1) ;
  CHECK(img.inSet(100)) ;
}

TEST_CASE("rMap: preimage") {
  rMap rm ;
  vector<pair<Entity,Entity>> entries = {
    {1, 100}, {2, 200}, {3, 100}, {4, 300}
  } ;
  rm.bulk_insert(entries) ;

  entitySet codomain ;
  codomain += 100 ;
  auto result = rm.preimage(codomain) ;

  CHECK(result.first.size() == 2) ;
  CHECK(result.first.inSet(1)) ;
  CHECK(result.first.inSet(3)) ;
}

TEST_CASE("rMap: const_rMap read-only access") {
  rMap rm ;
  vector<pair<Entity,Entity>> entries = {{1, 10}, {2, 20}} ;
  rm.bulk_insert(entries) ;

  storeRepP rep = rm.Rep() ;
  const_rMap crm(rep) ;

  CHECK(crm[1] == 10) ;
  CHECK(crm[2] == 20) ;
  CHECK(crm.domain().size() == 2) ;
}

TEST_CASE("rMap: allocate sets domain") {
  rMap rm ;
  entitySet ptn ;
  ptn += interval(10, 14) ;
  rm.allocate(ptn) ;

  CHECK(rm.domain() == ptn) ;
  CHECK(rm[10] == 0) ;
  CHECK(rm[14] == 0) ;
}

TEST_CASE("rMap: Print and Input round-trip") {
  rMap rm ;
  vector<pair<Entity,Entity>> entries = {{5, 50}, {10, 100}, {15, 150}} ;
  rm.bulk_insert(entries) ;

  stringstream ss ;
  rm.Print(ss) ;

  rMap rm2 ;
  rm2.Input(ss) ;

  CHECK(rm2.domain() == rm.domain()) ;
  CHECK(rm2[5] == 50) ;
  CHECK(rm2[10] == 100) ;
  CHECK(rm2[15] == 150) ;
}

TEST_CASE("rMap: image of full domain") {
  rMap rm ;
  vector<pair<Entity,Entity>> entries = {
    {1, 100}, {2, 200}, {3, 300}
  } ;
  rm.bulk_insert(entries) ;

  entitySet img = rm.image(rm.domain()) ;
  CHECK(img.size() == 3) ;
  CHECK(img.inSet(100)) ;
  CHECK(img.inSet(200)) ;
  CHECK(img.inSet(300)) ;
}

// =========================================================================
// rmultiMap tests
// =========================================================================

TEST_CASE("rmultiMap: basic bulk_insert and access") {
  rmultiMap rmm ;
  vector<pair<Entity,Entity>> entries = {
    {1, 10}, {1, 20}, {1, 30},
    {2, 40},
    {3, 50}, {3, 60}
  } ;
  rmm.bulk_insert(entries) ;

  CHECK(rmm.domain().size() == 3) ;
  CHECK(rmm.num_elems(1) == 3) ;
  CHECK(rmm.num_elems(2) == 1) ;
  CHECK(rmm.num_elems(3) == 2) ;

  rmultiMap_ref v1 = rmm[1] ;
  CHECK(v1.size() == 3) ;
  CHECK(v1[0] == 10) ;
  CHECK(v1[1] == 20) ;
  CHECK(v1[2] == 30) ;
}

TEST_CASE("rmultiMap: stable sort preserves insertion order") {
  rmultiMap rmm ;
  vector<pair<Entity,Entity>> batch1 = {
    {5, 100}, {5, 200}, {5, 300}
  } ;
  rmm.bulk_insert(batch1) ;

  vector<pair<Entity,Entity>> batch2 = {
    {5, 400}, {5, 500}
  } ;
  rmm.bulk_insert(batch2) ;

  CHECK(rmm.num_elems(5) == 5) ;
  rmultiMap_ref vals = rmm[5] ;
  CHECK(vals[0] == 100) ;
  CHECK(vals[1] == 200) ;
  CHECK(vals[2] == 300) ;
  CHECK(vals[3] == 400) ;
  CHECK(vals[4] == 500) ;
}

TEST_CASE("rmultiMap: allocate with zero mappings") {
  rmultiMap rmm ;
  entitySet ptn ;
  ptn += 1 ;
  ptn += 2 ;
  ptn += 3 ;
  rmm.allocate(ptn) ;

  CHECK(rmm.domain() == ptn) ;
  CHECK(rmm.num_elems(1) == 0) ;
  CHECK(rmm.num_elems(2) == 0) ;
  CHECK(rmm.num_elems(3) == 0) ;

  vector<pair<Entity,Entity>> entries = {{2, 100}, {2, 200}} ;
  rmm.bulk_insert(entries) ;

  CHECK(rmm.domain() == ptn) ;
  CHECK(rmm.num_elems(1) == 0) ;
  CHECK(rmm.num_elems(2) == 2) ;
  CHECK(rmm.num_elems(3) == 0) ;
}

TEST_CASE("rmultiMap: bulk_erase") {
  rmultiMap rmm ;
  vector<pair<Entity,Entity>> entries = {
    {1, 10}, {1, 20},
    {2, 30},
    {3, 40}, {3, 50}, {3, 60}
  } ;
  rmm.bulk_insert(entries) ;

  entitySet to_erase ;
  to_erase += 1 ;
  to_erase += 3 ;
  rmm.bulk_erase(to_erase) ;

  CHECK(rmm.domain().size() == 1) ;
  CHECK(rmm.num_elems(2) == 1) ;
  CHECK(rmm[2][0] == 30) ;
}

TEST_CASE("rmultiMap: image") {
  rmultiMap rmm ;
  vector<pair<Entity,Entity>> entries = {
    {1, 100}, {1, 200},
    {2, 200}, {2, 300},
    {3, 400}
  } ;
  rmm.bulk_insert(entries) ;

  entitySet query ;
  query += 1 ;
  query += 2 ;
  entitySet img = rmm.image(query) ;

  CHECK(img.size() == 3) ;
  CHECK(img.inSet(100)) ;
  CHECK(img.inSet(200)) ;
  CHECK(img.inSet(300)) ;
}

TEST_CASE("rmultiMap: preimage") {
  rmultiMap rmm ;
  vector<pair<Entity,Entity>> entries = {
    {1, 100}, {1, 200},
    {2, 300},
    {3, 100}, {3, 400}
  } ;
  rmm.bulk_insert(entries) ;

  entitySet codomain ;
  codomain += 100 ;
  auto result = rmm.preimage(codomain) ;

  CHECK(result.first.size() == 2) ;
  CHECK(result.first.inSet(1)) ;
  CHECK(result.first.inSet(3)) ;
}

TEST_CASE("rmultiMap: const_rmultiMap read-only access") {
  rmultiMap rmm ;
  vector<pair<Entity,Entity>> entries = {
    {1, 10}, {1, 20},
    {2, 30}
  } ;
  rmm.bulk_insert(entries) ;

  storeRepP rep = rmm.Rep() ;
  const_rmultiMap crmm(rep) ;

  CHECK(crmm.num_elems(1) == 2) ;
  CHECK(crmm[1][0] == 10) ;
  CHECK(crmm[1][1] == 20) ;
  CHECK(crmm[2][0] == 30) ;
}

TEST_CASE("rmultiMap: range-based for over ref") {
  rmultiMap rmm ;
  vector<pair<Entity,Entity>> entries = {
    {1, 10}, {1, 20}, {1, 30}
  } ;
  rmm.bulk_insert(entries) ;

  vector<Entity> collected ;
  for(const auto &val : rmm[1])
    collected.push_back(val) ;

  REQUIRE(collected.size() == 3) ;
  CHECK(collected[0] == 10) ;
  CHECK(collected[1] == 20) ;
  CHECK(collected[2] == 30) ;
}

TEST_CASE("rmultiMap: Print and Input round-trip") {
  rmultiMap rmm ;
  entitySet ptn ;
  ptn += 5 ;
  ptn += 10 ;
  ptn += 15 ;
  rmm.allocate(ptn) ;

  vector<pair<Entity,Entity>> entries = {
    {5, 50}, {5, 55},
    {15, 150}
  } ;
  rmm.bulk_insert(entries) ;

  stringstream ss ;
  rmm.Print(ss) ;

  rmultiMap rmm2 ;
  rmm2.Input(ss) ;

  CHECK(rmm2.domain() == rmm.domain()) ;
  CHECK(rmm2.num_elems(5) == 2) ;
  CHECK(rmm2[5][0] == 50) ;
  CHECK(rmm2[5][1] == 55) ;
  CHECK(rmm2.num_elems(10) == 0) ;
  CHECK(rmm2.num_elems(15) == 1) ;
  CHECK(rmm2[15][0] == 150) ;
}

TEST_CASE("rmultiMap: empty map") {
  rmultiMap rmm ;
  CHECK(rmm.domain() == EMPTY) ;
  CHECK(rmm.num_elems(1) == 0) ;

  rmultiMap_ref vals = rmm[1] ;
  CHECK(vals.size() == 0) ;
}

TEST_CASE("rmultiMap: mixed entities with varying target counts") {
  rmultiMap rmm ;
  entitySet ptn ;
  ptn += interval(1, 5) ;
  rmm.allocate(ptn) ;

  vector<pair<Entity,Entity>> entries ;
  entries.push_back({2, 20}) ;
  entries.push_back({3, 30}) ;
  entries.push_back({3, 31}) ;
  entries.push_back({3, 32}) ;
  entries.push_back({5, 50}) ;
  entries.push_back({5, 51}) ;

  rmm.bulk_insert(entries) ;

  CHECK(rmm.domain() == ptn) ;
  CHECK(rmm.num_elems(1) == 0) ;
  CHECK(rmm.num_elems(2) == 1) ;
  CHECK(rmm.num_elems(3) == 3) ;
  CHECK(rmm.num_elems(4) == 0) ;
  CHECK(rmm.num_elems(5) == 2) ;
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv) ;

  doctest::Context context ;
  context.applyCommandLine(argc, argv) ;
  const int result = context.run() ;

  Loci::Finalize() ;
  return result ;
}
