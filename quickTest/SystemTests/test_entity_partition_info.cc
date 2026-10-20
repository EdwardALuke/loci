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
#include <distribute.h>
#include <distribute_io.h>
#include <mpi_containerIO.h>
#include <partition.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

using namespace Loci;
using std::vector;
using std::pair;
using std::make_pair;

namespace {

  // Number of entities per processor for the default test configuration.
  const int ENTITIES_PER_PROC = 20;

  // Build an entityPartitionInfoP from a simple contiguous assignment
  // of global entities to processors.  Each processor owns
  // ENTITIES_PER_PROC entities in global numbering, and the local
  // numbering starts from 0 on each processor.
  //
  // Global numbering: proc r owns [r*N, (r+1)*N)
  // File numbering:   identical to global numbering (identity mapping)
  // Local numbering:  each processor renumbers its entities starting at 0
  //
  // All entities belong to key domain 0.
  struct TestDistribution {
    entitySet local_entities;
    entitySet global_entities;
    Map l2g;
    Map l2f;
    store<unsigned char> key_domain;
    fact_db::distribute_infoP dist;
    entityPartitionInfoP ptn_info;
    vector<dataPartitionP> ptn_list;

    TestDistribution(int entities_per_proc = ENTITIES_PER_PROC,
                     MPI_Comm comm = MPI_COMM_WORLD) {
      int p, r;
      MPI_Comm_size(comm, &p);
      MPI_Comm_rank(comm, &r);

      int n = entities_per_proc;
      int global_start = r * n;
      global_entities = entitySet(interval(global_start,
                                           global_start + n - 1));
      local_entities = entitySet(interval(0, n - 1));

      l2g.allocate(local_entities);
      l2f.allocate(local_entities);
      key_domain.allocate(local_entities);

      FORALL(local_entities, i) {
        l2g[i] = global_start + i;
        l2f[i] = global_start + i;
        key_domain[i] = 0;
      } ENDFORALL;

      dist = new fact_db::distribute_info;
      dist->myid = r;
      dist->isDistributed = 1;
      dist->my_entities = local_entities;
      dist->comp_entities = local_entities;
      dist->copy_total_size = 0;
      dist->xmit_total_size = 0;
      dist->l2g = l2g.Rep();
      dist->l2f = l2f.Rep();
      dist->key_domain = key_domain.Rep();

      vector<entitySet> global_ptn(p);
      for (int i = 0; i < p; ++i)
        global_ptn[i] = entitySet(interval(i * n, (i + 1) * n - 1));

      dataPartitionP dp = std::make_shared<dataPartitionGeneral>(
          global_ptn, comm);
      ptn_list.push_back(dp);

      ptn_info = std::make_shared<entityPartitionInfo>(
          local_entities, l2g, key_domain, l2f, ptn_list, comm);
    }
  };

  // Build a distribution with non-trivial file numbering.
  // File numbers are reversed within each processor's block.
  // For processor r owning global [r*N, (r+1)*N - 1], file numbers
  // are [(r+1)*N - 1, ..., r*N], i.e. reversed.
  // This exercises a non-identity local-to-file mapping even on a
  // single processor.
  struct ReversedFileDistribution {
    entitySet local_entities;
    Map l2g;
    Map l2f;
    store<unsigned char> key_domain;
    entityPartitionInfoP ptn_info;
    vector<dataPartitionP> ptn_list;

    ReversedFileDistribution(int entities_per_proc = ENTITIES_PER_PROC,
                              MPI_Comm comm = MPI_COMM_WORLD) {
      int p, r;
      MPI_Comm_size(comm, &p);
      MPI_Comm_rank(comm, &r);

      int n = entities_per_proc;
      int global_start = r * n;
      local_entities = entitySet(interval(0, n - 1));

      l2g.allocate(local_entities);
      l2f.allocate(local_entities);
      key_domain.allocate(local_entities);

      FORALL(local_entities, i) {
        l2g[i] = global_start + i;
        l2f[i] = global_start + (n - 1 - i);
        key_domain[i] = 0;
      } ENDFORALL;

      vector<entitySet> global_ptn(p);
      for (int i = 0; i < p; ++i)
        global_ptn[i] = entitySet(interval(i * n, (i + 1) * n - 1));

      dataPartitionP dp = std::make_shared<dataPartitionGeneral>(
          global_ptn, comm);
      ptn_list.push_back(dp);

      ptn_info = std::make_shared<entityPartitionInfo>(
          local_entities, l2g, key_domain, l2f, ptn_list, comm);
    }
  };

  // Build a distribution where local, global, and file numberings are
  // all different.  This is useful for testing even on a single processor.
  //
  // Local numbering:   [0, N)
  // Global numbering:  [100 + r*N, 100 + (r+1)*N)  (offset by 100)
  // File numbering:    stride-2 starting from r*N*2: {r*N*2, r*N*2+2, ...}
  struct ShuffledDistribution {
    entitySet local_entities;
    Map l2g;
    Map l2f;
    store<unsigned char> key_domain;
    entityPartitionInfoP ptn_info;
    vector<dataPartitionP> ptn_list;

    ShuffledDistribution(int entities_per_proc = ENTITIES_PER_PROC,
                         MPI_Comm comm = MPI_COMM_WORLD) {
      int p, r;
      MPI_Comm_size(comm, &p);
      MPI_Comm_rank(comm, &r);

      int n = entities_per_proc;
      int global_start = 100 + r * n;
      int file_start = r * n * 2;
      local_entities = entitySet(interval(0, n - 1));

      l2g.allocate(local_entities);
      l2f.allocate(local_entities);
      key_domain.allocate(local_entities);

      FORALL(local_entities, i) {
        l2g[i] = global_start + i;
        l2f[i] = file_start + i * 2;
        key_domain[i] = 0;
      } ENDFORALL;

      vector<entitySet> global_ptn(p);
      for (int i = 0; i < p; ++i)
        global_ptn[i] = entitySet(interval(100 + i * n,
                                           100 + (i + 1) * n - 1));

      dataPartitionP dp = std::make_shared<dataPartitionGeneral>(
          global_ptn, comm);
      ptn_list.push_back(dp);

      ptn_info = std::make_shared<entityPartitionInfo>(
          local_entities, l2g, key_domain, l2f, ptn_list, comm);
    }
  };

  // Build a distribution with multiple key domains.
  // Key domain 0: first half of entities per proc
  // Key domain 1: second half of entities per proc
  struct MultiKeyDomainDistribution {
    entitySet local_entities;
    Map l2g;
    Map l2f;
    store<unsigned char> key_domain;
    entityPartitionInfoP ptn_info;
    vector<dataPartitionP> ptn_list;

    MultiKeyDomainDistribution(int entities_per_proc = ENTITIES_PER_PROC,
                                MPI_Comm comm = MPI_COMM_WORLD) {
      int p, r;
      MPI_Comm_size(comm, &p);
      MPI_Comm_rank(comm, &r);

      int n = entities_per_proc;
      int half = n / 2;
      int global_start = r * n;
      local_entities = entitySet(interval(0, n - 1));

      l2g.allocate(local_entities);
      l2f.allocate(local_entities);
      key_domain.allocate(local_entities);

      FORALL(local_entities, i) {
        l2g[i] = global_start + i;
        l2f[i] = global_start + i;
        key_domain[i] = (i < half) ? 0 : 1;
      } ENDFORALL;

      vector<entitySet> global_ptn0(p), global_ptn1(p);
      for (int i = 0; i < p; ++i) {
        int gs = i * n;
        global_ptn0[i] = entitySet(interval(gs, gs + half - 1));
        global_ptn1[i] = entitySet(interval(gs + half, gs + n - 1));
      }

      dataPartitionP dp0 = std::make_shared<dataPartitionGeneral>(
          global_ptn0, comm);
      dataPartitionP dp1 = std::make_shared<dataPartitionGeneral>(
          global_ptn1, comm);
      ptn_list.push_back(dp0);
      ptn_list.push_back(dp1);

      ptn_info = std::make_shared<entityPartitionInfo>(
          local_entities, l2g, key_domain, l2f, ptn_list, comm);
    }
  };

} // namespace


// ---------------------------------------------------------------------------
// Basic accessor tests
// ---------------------------------------------------------------------------

TEST_CASE("myEntities returns the local entity set") {
  TestDistribution td;
  CHECK(td.ptn_info->myEntities() == td.local_entities);
}

TEST_CASE("getCommunicator returns the communicator used at construction") {
  TestDistribution td;
  int result;
  MPI_Comm_compare(td.ptn_info->getCommunicator(), MPI_COMM_WORLD, &result);
  CHECK((result == MPI_IDENT || result == MPI_CONGRUENT));
}

TEST_CASE("getLocal2GlobalMap covers the local domain") {
  TestDistribution td;
  storeRepP l2g_rep = td.ptn_info->getLocal2GlobalMap();
  CHECK(l2g_rep != (storeRep*)0);
  entitySet l2g_dom = l2g_rep->domain();
  CHECK((td.local_entities - l2g_dom) == EMPTY);
}

TEST_CASE("getLocal2FileMap covers the local domain") {
  TestDistribution td;
  storeRepP l2f_rep = td.ptn_info->getLocal2FileMap();
  CHECK(l2f_rep != (storeRep*)0);
  entitySet l2f_dom = l2f_rep->domain();
  CHECK((td.local_entities - l2f_dom) == EMPTY);
}


// ---------------------------------------------------------------------------
// getKeyDomain tests
// ---------------------------------------------------------------------------

TEST_CASE("getKeyDomain returns 0 for a single key domain distribution") {
  TestDistribution td;
  int kd = td.ptn_info->getKeyDomain(td.local_entities);
  CHECK(kd == 0);
}

TEST_CASE("getKeyDomain distinguishes key domains in a multi-domain setup") {
  MultiKeyDomainDistribution md;
  int n = ENTITIES_PER_PROC;
  int half = n / 2;

  entitySet kd0_entities = entitySet(interval(0, half - 1));
  entitySet kd1_entities = entitySet(interval(half, n - 1));

  int kd0 = md.ptn_info->getKeyDomain(kd0_entities);
  int kd1 = md.ptn_info->getKeyDomain(kd1_entities);

  CHECK(kd0 == 0);
  CHECK(kd1 == 1);
}

TEST_CASE("getKeyDomain returns -1 for a mixed key domain set") {
  MultiKeyDomainDistribution md;
  int kd = md.ptn_info->getKeyDomain(md.local_entities);
  CHECK(kd == -1);
}


// ---------------------------------------------------------------------------
// Local2FileOrder -> File2LocalOrder round-trip (identity)
// ---------------------------------------------------------------------------

TEST_CASE("Local2FileOrder followed by File2LocalOrder is identity for store<int>") {
  TestDistribution td;

  store<int> values;
  values.allocate(td.local_entities);
  int global_start = MPI_rank * ENTITIES_PER_PROC;
  FORALL(td.local_entities, i) {
    values[i] = (global_start + i) * 100 + 7;
  } ENDFORALL;

  int offset = 0;
  storeRepP file_order = td.ptn_info->Local2FileOrder(
      values.Rep(), td.local_entities, offset);

  CHECK(file_order != (storeRep*)0);
  CHECK(file_order->domain().size() > 0);

  store<int> restored;
  restored.allocate(td.local_entities);
  storeRepP restored_rep = restored.Rep();
  td.ptn_info->File2LocalOrder(restored_rep, td.local_entities,
                                file_order, offset);

  FORALL(td.local_entities, i) {
    CAPTURE(i);
    CHECK(const_store<int>(restored_rep)[i] == values[i]);
  } ENDFORALL;
}

TEST_CASE("Local2FileOrder -> File2LocalOrder identity for store<double>") {
  TestDistribution td;

  store<double> values;
  values.allocate(td.local_entities);
  FORALL(td.local_entities, i) {
    values[i] = 3.14159 * (MPI_rank * ENTITIES_PER_PROC + i);
  } ENDFORALL;

  int offset = 0;
  storeRepP file_order = td.ptn_info->Local2FileOrder(
      values.Rep(), td.local_entities, offset);

  store<double> restored;
  restored.allocate(td.local_entities);
  storeRepP restored_rep = restored.Rep();
  td.ptn_info->File2LocalOrder(restored_rep, td.local_entities,
                                file_order, offset);

  FORALL(td.local_entities, i) {
    CAPTURE(i);
    CHECK(const_store<double>(restored_rep)[i] ==
          doctest::Approx(values[i]));
  } ENDFORALL;
}

TEST_CASE("Local2FileOrder -> File2LocalOrder identity for store<vector<int>>") {
  TestDistribution td;

  store<vector<int> > values;
  values.allocate(td.local_entities);
  FORALL(td.local_entities, i) {
    int g = MPI_rank * ENTITIES_PER_PROC + i;
    values[i].push_back(g);
    values[i].push_back(g * 2);
    values[i].push_back(g * 3);
  } ENDFORALL;

  int offset = 0;
  storeRepP file_order = td.ptn_info->Local2FileOrder(
      values.Rep(), td.local_entities, offset);

  store<vector<int> > restored;
  restored.allocate(td.local_entities);
  storeRepP restored_rep = restored.Rep();
  td.ptn_info->File2LocalOrder(restored_rep, td.local_entities,
                                file_order, offset);

  FORALL(td.local_entities, i) {
    CAPTURE(i);
    CHECK(const_store<vector<int> >(restored_rep)[i] == values[i]);
  } ENDFORALL;
}


// ---------------------------------------------------------------------------
// File2LocalOrder -> Local2FileOrder round-trip (identity)
// ---------------------------------------------------------------------------

TEST_CASE("File2LocalOrder followed by Local2FileOrder is identity for store<int>") {
  TestDistribution td;

  // Build a container already in "file order": each processor holds a
  // contiguous slice of the global file numbering.
  int n = ENTITIES_PER_PROC;
  int p = MPI_processes;
  int total = n * p;

  // Create a file-ordered partition: distribute [0, total) across procs
  vector<entitySet> file_ptn = simplePartition(0, total - 1, MPI_COMM_WORLD);
  entitySet my_file_set = file_ptn[MPI_rank];

  int file_offset = my_file_set.Min();
  entitySet shifted_file_set = my_file_set;
  shifted_file_set <<= file_offset;

  store<int> file_values;
  file_values.allocate(shifted_file_set);
  FORALL(shifted_file_set, i) {
    file_values[i] = (i + file_offset) * 11 + 3;
  } ENDFORALL;

  // File2LocalOrder: convert file-ordered data to local numbering
  store<int> local_vals;
  local_vals.allocate(td.local_entities);
  storeRepP local_rep = local_vals.Rep();
  td.ptn_info->File2LocalOrder(local_rep, td.local_entities,
                                file_values.Rep(), file_offset);

  // Verify the local values match what we expect
  FORALL(td.local_entities, i) {
    int global_id = MPI_rank * n + i;
    CAPTURE(i);
    CHECK(const_store<int>(local_rep)[i] == global_id * 11 + 3);
  } ENDFORALL;

  // Now convert back: Local2FileOrder
  int offset2 = 0;
  storeRepP roundtrip = td.ptn_info->Local2FileOrder(
      local_rep, td.local_entities, offset2);

  // Verify the file-ordered data matches the original
  entitySet rt_dom = roundtrip->domain();
  const_store<int> rt_vals(roundtrip);
  FORALL(rt_dom, i) {
    int file_num = i + offset2;
    CAPTURE(file_num);
    CHECK(rt_vals[i] == file_num * 11 + 3);
  } ENDFORALL;
}


// ---------------------------------------------------------------------------
// Round-trip with reversed file numbering
// ---------------------------------------------------------------------------

TEST_CASE("round-trip identity with non-trivial file numbering") {
  ReversedFileDistribution rd;

  store<int> values;
  values.allocate(rd.local_entities);
  int n = ENTITIES_PER_PROC;
  FORALL(rd.local_entities, i) {
    values[i] = (MPI_rank * n + i) * 7 + 13;
  } ENDFORALL;

  int offset = 0;
  storeRepP file_order = rd.ptn_info->Local2FileOrder(
      values.Rep(), rd.local_entities, offset);

  store<int> restored;
  restored.allocate(rd.local_entities);
  storeRepP restored_rep = restored.Rep();
  rd.ptn_info->File2LocalOrder(restored_rep, rd.local_entities,
                                file_order, offset);

  FORALL(rd.local_entities, i) {
    CAPTURE(i);
    CHECK(const_store<int>(restored_rep)[i] == values[i]);
  } ENDFORALL;
}


// ---------------------------------------------------------------------------
// Subset round-trip: only convert a subset of entities
// ---------------------------------------------------------------------------

TEST_CASE("Local2FileOrder -> File2LocalOrder works on a strict subset") {
  TestDistribution td;
  int n = ENTITIES_PER_PROC;

  // Use only the first half of local entities
  int half = n / 2;
  entitySet subset = entitySet(interval(0, half - 1));

  store<int> values;
  values.allocate(td.local_entities);
  FORALL(td.local_entities, i) {
    values[i] = i * 1000 + MPI_rank;
  } ENDFORALL;

  int offset = 0;
  storeRepP file_order = td.ptn_info->Local2FileOrder(
      values.Rep(), subset, offset);

  store<int> restored;
  restored.allocate(subset);
  storeRepP restored_rep = restored.Rep();
  td.ptn_info->File2LocalOrder(restored_rep, subset,
                                file_order, offset);

  FORALL(subset, i) {
    CAPTURE(i);
    CHECK(const_store<int>(restored_rep)[i] == values[i]);
  } ENDFORALL;
}


// ---------------------------------------------------------------------------
// File ordering produces globally consistent numbering
// ---------------------------------------------------------------------------

TEST_CASE("Local2FileOrder produces a contiguous global file domain") {
  TestDistribution td;

  store<int> values;
  values.allocate(td.local_entities);
  FORALL(td.local_entities, i) {
    values[i] = i;
  } ENDFORALL;

  int offset = 0;
  storeRepP file_order = td.ptn_info->Local2FileOrder(
      values.Rep(), td.local_entities, offset);

  entitySet local_file_dom = file_order->domain();
  int local_min = local_file_dom.Min() + offset;
  int local_max = local_file_dom.Max() + offset;
  int local_size = local_file_dom.size();

  // Gather sizes to verify total
  int total_size = 0;
  MPI_Allreduce(&local_size, &total_size, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  CHECK(total_size == ENTITIES_PER_PROC * MPI_processes);

  // Verify non-overlapping ranges by gathering all min/max pairs
  vector<int> all_mins(MPI_processes), all_maxs(MPI_processes);
  MPI_Allgather(&local_min, 1, MPI_INT, all_mins.data(), 1, MPI_INT,
                MPI_COMM_WORLD);
  MPI_Allgather(&local_max, 1, MPI_INT, all_maxs.data(), 1, MPI_INT,
                MPI_COMM_WORLD);

  // Sort by min to check for overlaps
  vector<pair<int, int> > ranges(MPI_processes);
  for (int i = 0; i < MPI_processes; ++i)
    ranges[i] = make_pair(all_mins[i], all_maxs[i]);
  std::sort(ranges.begin(), ranges.end());

  bool no_overlap = true;
  for (int i = 1; i < MPI_processes; ++i) {
    if (ranges[i].first <= ranges[i - 1].second) {
      no_overlap = false;
      break;
    }
  }
  CHECK(no_overlap);
}


// ---------------------------------------------------------------------------
// Empty domain handling
// ---------------------------------------------------------------------------

TEST_CASE("Local2FileOrder handles empty domain gracefully") {
  TestDistribution td;

  store<int> values;
  values.allocate(td.local_entities);
  FORALL(td.local_entities, i) {
    values[i] = i;
  } ENDFORALL;

  int offset = 0;
  storeRepP file_order = td.ptn_info->Local2FileOrder(
      values.Rep(), EMPTY, offset);

  // The result domain should be empty since we passed an empty domain
  // (intersected with myEntities gives empty).  Gather globally to check
  // that this is collectively empty.
  int local_sz = file_order->domain().size();
  int total_sz = 0;
  MPI_Allreduce(&local_sz, &total_sz, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  CHECK(total_sz == 0);
}

TEST_CASE("File2LocalOrder into empty result set produces no change") {
  TestDistribution td;

  store<int> file_vals;
  file_vals.allocate(EMPTY);

  store<int> local_vals;
  local_vals.allocate(EMPTY);
  storeRepP local_rep = local_vals.Rep();

  td.ptn_info->File2LocalOrder(local_rep, EMPTY, file_vals.Rep(), 0);

  CHECK(local_rep->domain() == EMPTY);
}


// ---------------------------------------------------------------------------
// dataPartition integration: partitionEntitySet and getAllocation
// ---------------------------------------------------------------------------

TEST_CASE("dataPartitionGeneral partitionEntitySet covers all entities") {
  int n = ENTITIES_PER_PROC;
  int p = MPI_processes;
  int total = n * p;

  vector<entitySet> ptn(p);
  for (int i = 0; i < p; ++i)
    ptn[i] = entitySet(interval(i * n, (i + 1) * n - 1));

  dataPartitionP dp = std::make_shared<dataPartitionGeneral>(
      ptn, MPI_COMM_WORLD);

  // Partition the entire range
  entitySet all = entitySet(interval(0, total - 1));
  vector<pair<int, entitySet> > result = dp->partitionEntitySet(all);

  entitySet covered;
  for (size_t i = 0; i < result.size(); ++i)
    covered += result[i].second;
  CHECK(covered == all);
}

TEST_CASE("dataPartitionGeneral getAllocation matches construction") {
  int n = ENTITIES_PER_PROC;
  int p = MPI_processes;

  vector<entitySet> ptn(p);
  for (int i = 0; i < p; ++i)
    ptn[i] = entitySet(interval(i * n, (i + 1) * n - 1));

  dataPartitionP dp = std::make_shared<dataPartitionGeneral>(
      ptn, MPI_COMM_WORLD);

  for (int i = 0; i < p; ++i)
    CHECK(dp->getAllocation(i) == ptn[i]);
}

TEST_CASE("dataPartitionGeneral processorLookup returns correct owners") {
  int n = ENTITIES_PER_PROC;
  int p = MPI_processes;

  vector<entitySet> ptn(p);
  for (int i = 0; i < p; ++i)
    ptn[i] = entitySet(interval(i * n, (i + 1) * n - 1));

  dataPartitionP dp = std::make_shared<dataPartitionGeneral>(
      ptn, MPI_COMM_WORLD);

  vector<int> entities;
  for (int i = 0; i < p; ++i) {
    entities.push_back(i * n);
    entities.push_back(i * n + n / 2);
    entities.push_back((i + 1) * n - 1);
  }

  vector<int> procs;
  dp->processorLookup(procs, entities);

  REQUIRE(procs.size() == entities.size());
  for (size_t j = 0; j < entities.size(); ++j) {
    int expected_proc = entities[j] / n;
    CAPTURE(entities[j]);
    CHECK(procs[j] == expected_proc);
  }
}


// ---------------------------------------------------------------------------
// dataPartitionSplits: verify split-based partition
// ---------------------------------------------------------------------------

TEST_CASE("dataPartitionSplits partitions correctly") {
  int n = ENTITIES_PER_PROC;
  int p = MPI_processes;

  vector<int> splits(p);
  for (int i = 0; i < p; ++i)
    splits[i] = i * n;

  dataPartitionP dp = createPartition(splits, MPI_COMM_WORLD);

  // Note: dataPartitionSplits overrides splits[0] with INT_MIN+1,
  // so processor 0 owns [INT_MIN+1, splits[1]-1].
  // Processors 1..p-1 own [splits[i], splits[i+1]-1].
  for (int i = 1; i < p; ++i) {
    entitySet alloc = dp->getAllocation(i);
    CHECK(alloc.size() > 0);
    CHECK(alloc.Min() == i * n);
  }

  // Verify processorLookup for entities clearly in range
  vector<int> entities;
  for (int i = 0; i < p; ++i)
    entities.push_back(i * n + n / 2);

  vector<int> procs;
  dp->processorLookup(procs, entities);

  for (int i = 0; i < p; ++i) {
    CAPTURE(i);
    CHECK(procs[i] == i);
  }
}


// ---------------------------------------------------------------------------
// dataPartitionComputed: verify algorithmic partition
// ---------------------------------------------------------------------------

TEST_CASE("dataPartitionComputed partitions uniformly") {
  int n = ENTITIES_PER_PROC;
  int p = MPI_processes;

  dataPartitionP dp = createPartition(0, n, MPI_COMM_WORLD);

  for (int i = 0; i < p; ++i) {
    entitySet alloc = dp->getAllocation(i);
    CHECK(alloc == entitySet(interval(i * n, (i + 1) * n - 1)));
  }

  vector<int> test_entities;
  for (int i = 0; i < p; ++i) {
    test_entities.push_back(i * n);
    test_entities.push_back((i + 1) * n - 1);
  }

  vector<int> procs;
  dp->processorLookup(procs, test_entities);

  for (int i = 0; i < p; ++i) {
    CAPTURE(i);
    CHECK(procs[2 * i] == i);
    CHECK(procs[2 * i + 1] == i);
  }
}

// Note: dataPartitionComputed::partitionEntitySet has a known bug
// (infinite loop due to commented-out set reduction), so we only test
// getAllocation and processorLookup here.


// ---------------------------------------------------------------------------
// Round-trip with fully shuffled numberings (exercises single-proc path)
// ---------------------------------------------------------------------------

TEST_CASE("round-trip identity with shuffled local/global/file numberings") {
  ShuffledDistribution sd;

  store<int> values;
  values.allocate(sd.local_entities);
  FORALL(sd.local_entities, i) {
    values[i] = i * 37 + MPI_rank * 1000 + 5;
  } ENDFORALL;

  int offset = 0;
  storeRepP file_order = sd.ptn_info->Local2FileOrder(
      values.Rep(), sd.local_entities, offset);

  store<int> restored;
  restored.allocate(sd.local_entities);
  storeRepP restored_rep = restored.Rep();
  sd.ptn_info->File2LocalOrder(restored_rep, sd.local_entities,
                                file_order, offset);

  FORALL(sd.local_entities, i) {
    CAPTURE(i);
    CHECK(const_store<int>(restored_rep)[i] == values[i]);
  } ENDFORALL;
}


// ---------------------------------------------------------------------------
// Round-trip with a larger and more scattered data set
// ---------------------------------------------------------------------------

TEST_CASE("round-trip identity with storeVec<double>") {
  TestDistribution td;
  int vec_size = 3;

  storeVec<double> values;
  values.setVecSize(vec_size);
  values.allocate(td.local_entities);
  FORALL(td.local_entities, i) {
    int g = MPI_rank * ENTITIES_PER_PROC + i;
    for (int j = 0; j < vec_size; ++j)
      values[i][j] = g * 10.0 + j * 0.1;
  } ENDFORALL;

  int offset = 0;
  storeRepP file_order = td.ptn_info->Local2FileOrder(
      values.Rep(), td.local_entities, offset);

  storeVec<double> restored;
  restored.setVecSize(vec_size);
  restored.allocate(td.local_entities);
  storeRepP restored_rep = restored.Rep();
  td.ptn_info->File2LocalOrder(restored_rep, td.local_entities,
                                file_order, offset);

  bool all_match = true;
  FORALL(td.local_entities, i) {
    for (int j = 0; j < vec_size; ++j) {
      if (std::fabs(restored[i][j] - values[i][j]) > 1e-12) {
        all_match = false;
      }
    }
  } ENDFORALL;
  CHECK(all_match);
}


// ---------------------------------------------------------------------------
// Single-entity-per-proc edge case
// ---------------------------------------------------------------------------

TEST_CASE("round-trip with one entity per processor") {
  int p, r;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &r);

  entitySet local = entitySet(interval(0, 0));
  Map l2g, l2f;
  store<unsigned char> kd;
  l2g.allocate(local);
  l2f.allocate(local);
  kd.allocate(local);
  l2g[0] = r;
  l2f[0] = r;
  kd[0] = 0;

  vector<entitySet> gptn(p);
  for (int i = 0; i < p; ++i)
    gptn[i] = entitySet(interval(i, i));

  dataPartitionP dp = std::make_shared<dataPartitionGeneral>(
      gptn, MPI_COMM_WORLD);
  vector<dataPartitionP> ptn_list;
  ptn_list.push_back(dp);

  entityPartitionInfoP ptn_info = std::make_shared<entityPartitionInfo>(
      local, l2g, kd, l2f, ptn_list, MPI_COMM_WORLD);

  store<int> values;
  values.allocate(local);
  values[0] = r * 42 + 7;

  int offset = 0;
  storeRepP file_order = ptn_info->Local2FileOrder(
      values.Rep(), local, offset);

  store<int> restored;
  restored.allocate(local);
  storeRepP restored_rep = restored.Rep();
  ptn_info->File2LocalOrder(restored_rep, local, file_order, offset);

  CHECK(const_store<int>(restored_rep)[0] == r * 42 + 7);
}


// ---------------------------------------------------------------------------
// Consistency: two conversions to file order produce the same result
// ---------------------------------------------------------------------------

TEST_CASE("Local2FileOrder is deterministic across repeated calls") {
  TestDistribution td;

  store<int> values;
  values.allocate(td.local_entities);
  FORALL(td.local_entities, i) {
    values[i] = i + MPI_rank * 1000;
  } ENDFORALL;

  int offset1 = 0, offset2 = 0;
  storeRepP file_order1 = td.ptn_info->Local2FileOrder(
      values.Rep(), td.local_entities, offset1);
  storeRepP file_order2 = td.ptn_info->Local2FileOrder(
      values.Rep(), td.local_entities, offset2);

  CHECK(offset1 == offset2);
  CHECK(file_order1->domain() == file_order2->domain());

  bool all_match = true;
  const_store<int> fo1(file_order1), fo2(file_order2);
  FORALL(file_order1->domain(), i) {
    if (fo1[i] != fo2[i]) {
      all_match = false;
    }
  } ENDFORALL;
  CHECK(all_match);
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
