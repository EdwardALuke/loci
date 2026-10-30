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
#include <RMapVec.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <vector>
#include <sstream>

using namespace Loci ;
using namespace std ;

typedef Array<Entity,2> VEC2 ;
typedef Array<Entity,3> VEC3 ;

// =========================================================================
// rMapVec<2> tests (e.g. edge-to-node mapping)
// =========================================================================

TEST_CASE("rMapVec<2>: bulk_insert and const access") {
  rMapVec<2> rm ;
  vector<pair<Entity,VEC2>> entries ;

  VEC2 v1 ; v1[0] = 100 ; v1[1] = 200 ;
  VEC2 v2 ; v2[0] = 300 ; v2[1] = 400 ;
  VEC2 v3 ; v3[0] = 500 ; v3[1] = 600 ;
  entries.push_back({10, v1}) ;
  entries.push_back({5, v2}) ;
  entries.push_back({20, v3}) ;
  rm.bulk_insert(entries) ;

  CHECK(rm.domain().size() == 3) ;
  CHECK(rm[5][0] == 300) ;
  CHECK(rm[5][1] == 400) ;
  CHECK(rm[10][0] == 100) ;
  CHECK(rm[10][1] == 200) ;
  CHECK(rm[20][0] == 500) ;
  CHECK(rm[20][1] == 600) ;
}

TEST_CASE("rMapVec<2>: bulk_insert overwrites") {
  rMapVec<2> rm ;
  VEC2 v1 ; v1[0] = 10 ; v1[1] = 20 ;
  VEC2 v2 ; v2[0] = 30 ; v2[1] = 40 ;
  vector<pair<Entity,VEC2>> entries = {{1, v1}, {2, v2}} ;
  rm.bulk_insert(entries) ;

  VEC2 v3 ; v3[0] = 99 ; v3[1] = 88 ;
  vector<pair<Entity,VEC2>> updates = {{2, v3}} ;
  rm.bulk_insert(updates) ;

  CHECK(rm[1][0] == 10) ;
  CHECK(rm[1][1] == 20) ;
  CHECK(rm[2][0] == 99) ;
  CHECK(rm[2][1] == 88) ;
}

TEST_CASE("rMapVec<2>: bulk_erase") {
  rMapVec<2> rm ;
  VEC2 v1 ; v1[0] = 10 ; v1[1] = 20 ;
  VEC2 v2 ; v2[0] = 30 ; v2[1] = 40 ;
  VEC2 v3 ; v3[0] = 50 ; v3[1] = 60 ;
  vector<pair<Entity,VEC2>> entries = {{1, v1}, {2, v2}, {3, v3}} ;
  rm.bulk_insert(entries) ;

  entitySet to_erase ;
  to_erase += 2 ;
  rm.bulk_erase(to_erase) ;

  CHECK(rm.domain().size() == 2) ;
  CHECK(rm[1][0] == 10) ;
  CHECK(rm[3][0] == 50) ;
}

TEST_CASE("rMapVec<2>: allocate") {
  rMapVec<2> rm ;
  entitySet ptn ;
  ptn += interval(1, 3) ;
  rm.allocate(ptn) ;

  CHECK(rm.domain() == ptn) ;
  CHECK(rm[1][0] == 0) ;
  CHECK(rm[1][1] == 0) ;
  CHECK(rm[3][0] == 0) ;
}

TEST_CASE("rMapVec<2>: image") {
  rMapVec<2> rm ;
  VEC2 v1 ; v1[0] = 100 ; v1[1] = 200 ;
  VEC2 v2 ; v2[0] = 200 ; v2[1] = 300 ;
  VEC2 v3 ; v3[0] = 400 ; v3[1] = 100 ;
  vector<pair<Entity,VEC2>> entries = {{1, v1}, {2, v2}, {3, v3}} ;
  rm.bulk_insert(entries) ;

  entitySet query ;
  query += 1 ;
  query += 2 ;
  entitySet img = rm.image(query) ;

  CHECK(img.size() == 3) ;
  CHECK(img.inSet(100)) ;
  CHECK(img.inSet(200)) ;
  CHECK(img.inSet(300)) ;
}

TEST_CASE("rMapVec<2>: preimage") {
  rMapVec<2> rm ;
  VEC2 v1 ; v1[0] = 100 ; v1[1] = 200 ;
  VEC2 v2 ; v2[0] = 300 ; v2[1] = 400 ;
  VEC2 v3 ; v3[0] = 100 ; v3[1] = 400 ;
  vector<pair<Entity,VEC2>> entries = {{1, v1}, {2, v2}, {3, v3}} ;
  rm.bulk_insert(entries) ;

  entitySet codomain ;
  codomain += 100 ;
  codomain += 200 ;
  auto result = rm.preimage(codomain) ;

  // domaini: entities where ALL targets are in codomain -> entity 1
  CHECK(result.first.size() == 1) ;
  CHECK(result.first.inSet(1)) ;

  // domainu: entities where ANY target is in codomain -> entities 1, 3
  CHECK(result.second.size() == 2) ;
  CHECK(result.second.inSet(1)) ;
  CHECK(result.second.inSet(3)) ;
}

TEST_CASE("rMapVec<2>: const_rMapVec read-only access") {
  rMapVec<2> rm ;
  VEC2 v1 ; v1[0] = 10 ; v1[1] = 20 ;
  VEC2 v2 ; v2[0] = 30 ; v2[1] = 40 ;
  vector<pair<Entity,VEC2>> entries = {{1, v1}, {2, v2}} ;
  rm.bulk_insert(entries) ;

  storeRepP rep = rm.Rep() ;
  const_rMapVec<2> crm(rep) ;

  CHECK(crm[1][0] == 10) ;
  CHECK(crm[1][1] == 20) ;
  CHECK(crm[2][0] == 30) ;
  CHECK(crm[2][1] == 40) ;
  CHECK(crm.domain().size() == 2) ;
}

TEST_CASE("rMapVec<2>: Print and Input round-trip") {
  rMapVec<2> rm ;
  VEC2 v1 ; v1[0] = 50 ; v1[1] = 55 ;
  VEC2 v2 ; v2[0] = 100 ; v2[1] = 105 ;
  vector<pair<Entity,VEC2>> entries = {{5, v1}, {10, v2}} ;
  rm.bulk_insert(entries) ;

  stringstream ss ;
  rm.Print(ss) ;

  rMapVec<2> rm2 ;
  rm2.Input(ss) ;

  CHECK(rm2.domain() == rm.domain()) ;
  CHECK(rm2[5][0] == 50) ;
  CHECK(rm2[5][1] == 55) ;
  CHECK(rm2[10][0] == 100) ;
  CHECK(rm2[10][1] == 105) ;
}

// =========================================================================
// rMapVec<3> tests (e.g. face-to-node for triangles)
// =========================================================================

TEST_CASE("rMapVec<3>: basic operations") {
  rMapVec<3> rm ;
  VEC3 v1 ; v1[0] = 10 ; v1[1] = 20 ; v1[2] = 30 ;
  VEC3 v2 ; v2[0] = 40 ; v2[1] = 50 ; v2[2] = 60 ;
  vector<pair<Entity,VEC3>> entries = {{1, v1}, {2, v2}} ;
  rm.bulk_insert(entries) ;

  CHECK(rm.domain().size() == 2) ;
  CHECK(rm[1][0] == 10) ;
  CHECK(rm[1][1] == 20) ;
  CHECK(rm[1][2] == 30) ;
  CHECK(rm[2][0] == 40) ;
  CHECK(rm[2][1] == 50) ;
  CHECK(rm[2][2] == 60) ;
}

TEST_CASE("rMapVec<3>: image spans all M targets") {
  rMapVec<3> rm ;
  VEC3 v1 ; v1[0] = 100 ; v1[1] = 200 ; v1[2] = 300 ;
  VEC3 v2 ; v2[0] = 200 ; v2[1] = 400 ; v2[2] = 500 ;
  vector<pair<Entity,VEC3>> entries = {{1, v1}, {2, v2}} ;
  rm.bulk_insert(entries) ;

  entitySet img = rm.image(rm.domain()) ;

  CHECK(img.size() == 5) ;
  CHECK(img.inSet(100)) ;
  CHECK(img.inSet(200)) ;
  CHECK(img.inSet(300)) ;
  CHECK(img.inSet(400)) ;
  CHECK(img.inSet(500)) ;
}

TEST_CASE("rMapVec<3>: Print and Input round-trip") {
  rMapVec<3> rm ;
  VEC3 v1 ; v1[0] = 1 ; v1[1] = 2 ; v1[2] = 3 ;
  VEC3 v2 ; v2[0] = 4 ; v2[1] = 5 ; v2[2] = 6 ;
  vector<pair<Entity,VEC3>> entries = {{10, v1}, {20, v2}} ;
  rm.bulk_insert(entries) ;

  stringstream ss ;
  rm.Print(ss) ;

  rMapVec<3> rm2 ;
  rm2.Input(ss) ;

  CHECK(rm2.domain() == rm.domain()) ;
  CHECK(rm2[10][0] == 1) ;
  CHECK(rm2[10][1] == 2) ;
  CHECK(rm2[10][2] == 3) ;
  CHECK(rm2[20][0] == 4) ;
  CHECK(rm2[20][1] == 5) ;
  CHECK(rm2[20][2] == 6) ;
}

TEST_CASE("rMapVec<2>: empty map") {
  rMapVec<2> rm ;
  CHECK(rm.domain() == EMPTY) ;

  entitySet empty_erase ;
  rm.bulk_erase(empty_erase) ;
  CHECK(rm.domain() == EMPTY) ;

  vector<pair<Entity,VEC2>> empty_insert ;
  rm.bulk_insert(empty_insert) ;
  CHECK(rm.domain() == EMPTY) ;
}

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv) ;

  doctest::Context context ;
  context.applyCommandLine(argc, argv) ;
  const int result = context.run() ;

  Loci::Finalize() ;
  return result ;
}
