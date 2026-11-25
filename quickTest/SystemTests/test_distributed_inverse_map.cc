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
#include <distribute_container.h>
#include <partition.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

using namespace Loci;
using std::vector;
using std::pair;
using std::make_pair;
using std::set;

namespace {

  // Number of entities per processor in tests.
  const int N = 20;

  // Build a simple contiguous partition: proc r owns [r*N, (r+1)*N).
  struct TestPartition {
    int n, p, r;
    entitySet local_entities;
    dataPartitionP dp;

    TestPartition(int entities_per_proc = N,
                  MPI_Comm comm = MPI_COMM_WORLD) {
      n = entities_per_proc;
      MPI_Comm_size(comm, &p);
      MPI_Comm_rank(comm, &r);
      local_entities = entitySet(interval(r * n, (r + 1) * n - 1));
      vector<entitySet> ptn(p);
      for (int i = 0; i < p; ++i)
        ptn[i] = entitySet(interval(i * n, (i + 1) * n - 1));
      dp = std::make_shared<dataPartitionGeneral>(ptn, comm);
    }
  };

  // Collect a multiMap into a set of (key, value) pairs for comparison.
  set<pair<int,int>> multimap_to_set(const multiMap &mm) {
    set<pair<int,int>> s;
    entitySet dom = mm.domain();
    FORALL(dom, i) {
      int sz = mm.end(i) - mm.begin(i);
      for (int k = 0; k < sz; ++k)
        s.insert(make_pair(i, mm[i][k]));
    } ENDFORALL;
    return s;
  }

  // Collect a dmultiMap into a set of (key, value) pairs.
  set<pair<int,int>> dmultimap_to_set(const dmultiMap &dm) {
    set<pair<int,int>> s;
    entitySet dom = dm.domain();
    FORALL(dom, i) {
      for (size_t k = 0; k < dm[i].size(); ++k)
        s.insert(make_pair(i, dm[i][k]));
    } ENDFORALL;
    return s;
  }

  // Gather a set of pairs from all ranks into a globally consistent set.
  set<pair<int,int>> gather_global_set(const set<pair<int,int>> &local) {
    int local_sz = local.size();
    vector<int> buf(local_sz * 2);
    int j = 0;
    for (set<pair<int,int>>::const_iterator it = local.begin();
         it != local.end(); ++it) {
      buf[j++] = it->first;
      buf[j++] = it->second;
    }

    int p = MPI_processes;
    vector<int> counts(p), displs(p);
    int send_count = local_sz * 2;
    MPI_Allgather(&send_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    displs[0] = 0;
    for (int i = 1; i < p; ++i)
      displs[i] = displs[i-1] + counts[i-1];
    int total = displs[p-1] + counts[p-1];
    vector<int> rbuf(total);
    MPI_Allgatherv(buf.data(), send_count, MPI_INT,
                   rbuf.data(), counts.data(), displs.data(), MPI_INT,
                   MPI_COMM_WORLD);

    set<pair<int,int>> global;
    for (int i = 0; i < total; i += 2)
      global.insert(make_pair(rbuf[i], rbuf[i+1]));
    return global;
  }

  // Build the global pair set of a Map across all ranks.
  set<pair<int,int>> map_to_global_set(const Map &m, entitySet dom) {
    set<pair<int,int>> local;
    FORALL(dom, i) {
      local.insert(make_pair(i, m[i]));
    } ENDFORALL;
    return gather_global_set(local);
  }

} // namespace


// ===========================================================================
// Idempotent test for Map (bijection): inverse(inverse(M)) == M
//
// For a bijective Map f, let inv = inverse(f).  Then inverse(inv) should
// recover exactly the pairs of f.  We test this with the dataPartitionP API.
// ===========================================================================

TEST_CASE("idempotent: inverse(inverse(Map)) == Map for reversed bijection") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  // f(i) = total - 1 - i  (bijection)
  Map fwd;
  fwd.allocate(tp.local_entities);
  FORALL(tp.local_entities, i) {
    fwd[i] = total - 1 - i;
  } ENDFORALL;

  // inv = inverse(fwd)
  multiMap inv;
  distributed_inverseMap(inv, fwd, all, tp.local_entities, tp.dp);

  // inv2 = inverse(inv) — should recover fwd
  multiMap inv2;
  distributed_inverseMap(inv2, inv, all, tp.local_entities, tp.dp);

  set<pair<int,int>> original = gather_global_set(multimap_to_set(inv2));
  set<pair<int,int>> expected = map_to_global_set(fwd, tp.local_entities);
  CHECK(original == expected);
  CHECK(expected.size() == (size_t)total);
}


TEST_CASE("idempotent: inverse(inverse(Map)) == Map for shift-by-1") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  // f(i) = (i + 1) % total  (bijection, exercises cross-processor routing)
  Map fwd;
  fwd.allocate(tp.local_entities);
  FORALL(tp.local_entities, i) {
    fwd[i] = (i + 1) % total;
  } ENDFORALL;

  multiMap inv;
  distributed_inverseMap(inv, fwd, all, tp.local_entities, tp.dp);

  multiMap inv2;
  distributed_inverseMap(inv2, inv, all, tp.local_entities, tp.dp);

  set<pair<int,int>> recovered = gather_global_set(multimap_to_set(inv2));
  set<pair<int,int>> expected = map_to_global_set(fwd, tp.local_entities);
  CHECK(recovered == expected);
}


TEST_CASE("idempotent: inverse(inverse(Map)) == Map for identity") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  Map fwd;
  fwd.allocate(tp.local_entities);
  FORALL(tp.local_entities, i) {
    fwd[i] = i;
  } ENDFORALL;

  multiMap inv;
  distributed_inverseMap(inv, fwd, all, tp.local_entities, tp.dp);

  multiMap inv2;
  distributed_inverseMap(inv2, inv, all, tp.local_entities, tp.dp);

  set<pair<int,int>> recovered = gather_global_set(multimap_to_set(inv2));
  set<pair<int,int>> expected = map_to_global_set(fwd, tp.local_entities);
  CHECK(recovered == expected);
}


// ===========================================================================
// Idempotent test for many-to-one Map: inverse(inverse(M)) == M
//
// When f is many-to-one, inverse(f) is one-to-many, and
// inverse(inverse(f)) should recover the original pairs of f.
// ===========================================================================

TEST_CASE("idempotent: inverse(inverse(Map)) == Map for many-to-one") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  // f(i) = i / 2 * 2 (maps pairs of consecutive entities to the same value)
  // e.g. 0->0, 1->0, 2->2, 3->2, 4->4, 5->4, ...
  Map fwd;
  fwd.allocate(tp.local_entities);
  FORALL(tp.local_entities, i) {
    fwd[i] = (i / 2) * 2;
  } ENDFORALL;

  multiMap inv;
  distributed_inverseMap(inv, fwd, all, tp.local_entities, tp.dp);

  multiMap inv2;
  distributed_inverseMap(inv2, inv, all, tp.local_entities, tp.dp);

  set<pair<int,int>> recovered = gather_global_set(multimap_to_set(inv2));
  set<pair<int,int>> expected = map_to_global_set(fwd, tp.local_entities);
  CHECK(recovered == expected);
}


// ===========================================================================
// Idempotent test for multiMap: inverse(inverse(M)) == M
//
// For a general multiMap, inverse(inverse(M)) should recover the original
// set of (key, value) pairs.
// ===========================================================================

TEST_CASE("idempotent: inverse(inverse(multiMap)) == multiMap") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  // Each entity i maps to {i % total, (i+1) % total}
  store<int> sizes;
  sizes.allocate(tp.local_entities);
  FORALL(tp.local_entities, i) {
    sizes[i] = 2;
  } ENDFORALL;

  multiMap fwd;
  fwd.allocate(sizes);
  FORALL(tp.local_entities, i) {
    fwd[i][0] = i % total;
    fwd[i][1] = (i + 1) % total;
  } ENDFORALL;

  // inv = inverse(fwd)
  multiMap inv;
  distributed_inverseMap(inv, fwd, all, tp.local_entities, tp.dp);

  // inv2 = inverse(inv) — should recover fwd
  multiMap inv2;
  distributed_inverseMap(inv2, inv, all, tp.local_entities, tp.dp);

  set<pair<int,int>> recovered = gather_global_set(multimap_to_set(inv2));
  set<pair<int,int>> expected =
      gather_global_set(multimap_to_set(fwd));
  CHECK(recovered == expected);
  CHECK(expected.size() > 0);
}


// ===========================================================================
// Idempotent test for dmultiMap input -> dmultiMap result
// ===========================================================================

TEST_CASE("idempotent: inverse(inverse(dmultiMap)) == dmultiMap") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  // Each entity i maps to {(i*3) % total, (i*3+1) % total}
  dmultiMap fwd;
  FORALL(tp.local_entities, i) {
    vector<int, malloc_alloc<int>> v(2);
    v[0] = (i * 3) % total;
    v[1] = (i * 3 + 1) % total;
    v.swap(fwd[i]);
  } ENDFORALL;

  dmultiMap inv;
  distributed_inverseMap(inv, fwd, all, tp.local_entities, tp.dp);

  dmultiMap inv2;
  distributed_inverseMap(inv2, inv, all, tp.local_entities, tp.dp);

  set<pair<int,int>> recovered = gather_global_set(dmultimap_to_set(inv2));
  set<pair<int,int>> expected =
      gather_global_set(dmultimap_to_set(fwd));
  CHECK(recovered == expected);
  CHECK(expected.size() > 0);
}


// ===========================================================================
// Idempotent test for dMap -> dmultiMap
// ===========================================================================

TEST_CASE("idempotent: inverse(inverse(dMap)) == dMap") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  // dMap: f(i) = (i + N) % total (shift by N — bijection)
  dMap fwd;
  FORALL(tp.local_entities, i) {
    fwd[i] = (i + tp.n) % total;
  } ENDFORALL;

  // inv = inverse(fwd) as dmultiMap
  dmultiMap inv;
  distributed_inverseMap(inv, fwd, all, tp.local_entities, tp.dp);

  // inv2 = inverse(inv) — should recover fwd
  dmultiMap inv2;
  distributed_inverseMap(inv2, inv, all, tp.local_entities, tp.dp);

  // Build expected set from original dMap
  set<pair<int,int>> local_expected;
  FORALL(tp.local_entities, i) {
    local_expected.insert(make_pair(i, fwd[i]));
  } ENDFORALL;

  set<pair<int,int>> recovered = gather_global_set(dmultimap_to_set(inv2));
  set<pair<int,int>> expected = gather_global_set(local_expected);
  CHECK(recovered == expected);
}


// ===========================================================================
// Correctness: inverse of identity is identity
// ===========================================================================

TEST_CASE("inverse of identity Map is identity") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  Map fwd;
  fwd.allocate(tp.local_entities);
  FORALL(tp.local_entities, i) {
    fwd[i] = i;
  } ENDFORALL;

  multiMap inv;
  distributed_inverseMap(inv, fwd, all, tp.local_entities, tp.dp);

  // For identity, inverse[i] == {i} for every local entity
  FORALL(tp.local_entities, i) {
    int sz = inv.end(i) - inv.begin(i);
    CAPTURE(i);
    CHECK(sz == 1);
    if (sz == 1)
      CHECK(inv[i][0] == i);
  } ENDFORALL;
}


// ===========================================================================
// Correctness: cross-processor shift
// ===========================================================================

TEST_CASE("inverse of shift-by-1 routes correctly across processors") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  Map fwd;
  fwd.allocate(tp.local_entities);
  FORALL(tp.local_entities, i) {
    fwd[i] = (i + 1) % total;
  } ENDFORALL;

  multiMap inv;
  distributed_inverseMap(inv, fwd, all, tp.local_entities, tp.dp);

  // inverse[j] should be {(j-1+total) % total}
  FORALL(tp.local_entities, j) {
    int expected_source = (j - 1 + total) % total;
    int sz = inv.end(j) - inv.begin(j);
    CAPTURE(j);
    CHECK(sz == 1);
    if (sz == 1)
      CHECK(inv[j][0] == expected_source);
  } ENDFORALL;
}


// ===========================================================================
// Many-to-one: inverse collects all sources
// ===========================================================================

TEST_CASE("inverse of many-to-one Map collects all sources") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  // Every entity maps to entity 0
  Map fwd;
  fwd.allocate(tp.local_entities);
  FORALL(tp.local_entities, i) {
    fwd[i] = 0;
  } ENDFORALL;

  multiMap inv;
  distributed_inverseMap(inv, fwd, all, tp.local_entities, tp.dp);

  // Gather total number of inverse entries across all ranks
  int local_entries = 0;
  entitySet res_dom = inv.domain();
  FORALL(res_dom, i) {
    local_entries += inv.end(i) - inv.begin(i);
  } ENDFORALL;

  int global_entries = 0;
  MPI_Allreduce(&local_entries, &global_entries, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  CHECK(global_entries == total);

  // On rank 0, inv[0] should have exactly 'total' entries
  if (tp.r == 0) {
    CHECK(res_dom.inSet(0));
    if (res_dom.inSet(0)) {
      int sz = inv.end(0) - inv.begin(0);
      CHECK(sz == total);
    }
  }
}


// ===========================================================================
// Restricted image: inverse only on a subset
// ===========================================================================

TEST_CASE("inverse with restricted image returns only requested subset") {
  TestPartition tp;

  Map fwd;
  fwd.allocate(tp.local_entities);
  FORALL(tp.local_entities, i) {
    fwd[i] = i;
  } ENDFORALL;

  // Only collect inverse for the first half of local entities
  int half = tp.n / 2;
  entitySet image_subset = entitySet(interval(tp.r * tp.n,
                                              tp.r * tp.n + half - 1));

  multiMap inv;
  distributed_inverseMap(inv, fwd, image_subset, tp.local_entities, tp.dp);

  entitySet res_dom = inv.domain();
  CHECK((res_dom - image_subset) == EMPTY);

  FORALL(res_dom, i) {
    int sz = inv.end(i) - inv.begin(i);
    CHECK(sz == 1);
    if (sz == 1)
      CHECK(inv[i][0] == i);
  } ENDFORALL;
}


// ===========================================================================
// Empty preimage: no input produces zero entries
// ===========================================================================

TEST_CASE("inverse with empty preimage produces zero entries") {
  TestPartition tp;
  int total = tp.n * tp.p;
  entitySet all = entitySet(interval(0, total - 1));

  Map fwd;
  fwd.allocate(EMPTY);

  multiMap inv;
  distributed_inverseMap(inv, fwd, all, EMPTY, tp.dp);

  int local_entries = 0;
  entitySet res_dom = inv.domain();
  FORALL(res_dom, i) {
    local_entries += inv.end(i) - inv.begin(i);
  } ENDFORALL;

  int global_entries = 0;
  MPI_Allreduce(&local_entries, &global_entries, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  CHECK(global_entries == 0);
}


// ===========================================================================
// dataPartitionSplits: idempotent test with split-based partition
// ===========================================================================

TEST_CASE("idempotent inverse works with dataPartitionSplits") {
  int p = MPI_processes;
  int r = MPI_rank;
  int total = N * p;

  entitySet local = entitySet(interval(r * N, (r + 1) * N - 1));
  entitySet all = entitySet(interval(0, total - 1));

  vector<int> splits(p);
  for (int i = 0; i < p; ++i)
    splits[i] = i * N;

  dataPartitionP dp = createPartition(splits, MPI_COMM_WORLD);

  // Reversed bijection: f(i) = total - 1 - i
  Map fwd;
  fwd.allocate(local);
  FORALL(local, i) {
    fwd[i] = total - 1 - i;
  } ENDFORALL;

  multiMap inv;
  distributed_inverseMap(inv, fwd, all, local, dp);

  multiMap inv2;
  distributed_inverseMap(inv2, inv, all, local, dp);

  set<pair<int,int>> recovered = gather_global_set(multimap_to_set(inv2));
  set<pair<int,int>> expected = map_to_global_set(fwd, local);
  CHECK(recovered == expected);
}


// ===========================================================================
// dataPartitionComputed: idempotent test with algorithmic partition
// ===========================================================================

TEST_CASE("idempotent inverse works with dataPartitionComputed") {
  int p = MPI_processes;
  int r = MPI_rank;
  int total = N * p;

  entitySet local = entitySet(interval(r * N, (r + 1) * N - 1));
  entitySet all = entitySet(interval(0, total - 1));

  dataPartitionP dp = createPartition(0, N, MPI_COMM_WORLD);

  // Shift-by-N: f(i) = (i + N) % total
  Map fwd;
  fwd.allocate(local);
  FORALL(local, i) {
    fwd[i] = (i + N) % total;
  } ENDFORALL;

  multiMap inv;
  distributed_inverseMap(inv, fwd, all, local, dp);

  multiMap inv2;
  distributed_inverseMap(inv2, inv, all, local, dp);

  set<pair<int,int>> recovered = gather_global_set(multimap_to_set(inv2));
  set<pair<int,int>> expected = map_to_global_set(fwd, local);
  CHECK(recovered == expected);
}


int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  doctest::Context context;
  context.applyCommandLine(argc, argv);
  context.setOption("abort-after", 0);
  const int local_failure = context.run() == 0 ? 0 : 1;
  int failure = 0;
  MPI_Allreduce(&local_failure, &failure, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);

  std::cout.flush();
  std::cerr.flush();
  MPI_Barrier(MPI_COMM_WORLD);
  if (failure == 0 && MPI_rank == 0)
    std::cout << "SUCCESS!" << std::endl;

  Loci::Finalize();
  return failure;
}
