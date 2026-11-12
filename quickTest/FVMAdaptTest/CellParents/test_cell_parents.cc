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

#include <algorithm>
#include <limits>
#include <set>
#include <vector>
#include <Loci.h>
#include <FVMAdapt/gridInterface.h>
#include <FVMAdapt/dataxferDB.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

namespace {
  using CellPair=std::pair<int,int>; // One-based new and preceding file cell numbers.

  void check_mapping(int old_count,int new_count,const std::set<CellPair> &expected,
                     const CPTR<refinedGridData> &grid) {
    int actual_count=0;
    for(const auto &part:grid->local_cells) actual_count+=part.size();
    CHECK(actual_count==new_count);
    const_store<CellPair> mapping=DataXFER_DB.getItem("c2p");
    std::vector<int> local;
    for(Entity entry:mapping.domain()) {
      local.push_back(mapping[entry].first); local.push_back(mapping[entry].second);
    }
    const int size=local.size(); std::vector<int> sizes(MPI_processes),offsets(MPI_processes,0);
    MPI_Allgather(&size,1,MPI_INT,sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
    for(int rank=1;rank<MPI_processes;++rank) offsets[rank]=offsets[rank-1]+sizes[rank-1];
    std::vector<int> all(offsets.back()+sizes.back());
    MPI_Allgatherv(local.data(),size,MPI_INT,all.data(),sizes.data(),offsets.data(),MPI_INT,MPI_COMM_WORLD);
    std::set<CellPair> actual; std::set<int> old_cells,new_cells;
    for(std::size_t i=0;i<all.size();i+=2) {
      actual.emplace(all[i],all[i+1]); new_cells.insert(all[i]); old_cells.insert(all[i+1]);
    }
    CHECK(actual.size()*2==all.size());
    CHECK(actual==expected);
    CHECK(old_cells.size()==old_count); CHECK(new_cells.size()==new_count);
  }
}

/// A tetrahedron has four vertex-centered children. Refining one makes eleven
/// leaves; coarsening restores each complete family and accounts for every donor.
TEST_CASE("cell donors cover retention refinement and coarsening") {
  rule_db rules; rules.add_rules(global_rule_list); load_module("fvmadapt",rules);
  RefinementOptions options; options.splitFamily="corner"; options.centerMethod="reference";
  RefinementState state; CPTR<refinedGridData> grid; std::string error;
  REQUIRE_MESSAGE(onlineRefineMesh(grid,rules,state,0,"tet",options,error),error);
  check_mapping(1,1,{{1,1}},grid);
  for(int stage=0;stage<5;++stage) {
    CAPTURE(stage);
    store<char> tags; tags.allocate(grid->local_cells[MPI_rank]);
    int first=std::numeric_limits<int>::max();
    for(const auto &part:grid->local_cells) if(part!=EMPTY) first=std::min(first,part.Min());
    for(Entity cell:tags.domain()) tags[cell]=stage==0 ? 1:stage==1 ? 0:stage==2 ? (cell==first ? 1:0):2;
    REQUIRE_MESSAGE(onlineRefineMesh(grid,rules,state,tags.Rep(),"tet",options,error),error);
    std::set<CellPair> expected;
    if(stage==0) {
      for(int cell=1;cell<=4;++cell) expected.emplace(cell,1);
      check_mapping(1,4,expected,grid);
    } else if(stage==1) {
      for(int cell=1;cell<=4;++cell) expected.emplace(cell,cell);
      check_mapping(4,4,expected,grid);
    } else if(stage==2) {
      for(int cell=1;cell<=8;++cell) expected.emplace(cell,1);
      for(int cell=2;cell<=4;++cell) expected.emplace(cell+7,cell);
      check_mapping(4,11,expected,grid);
    } else if(stage==3) {
      for(int cell=1;cell<=8;++cell) expected.emplace(1,cell);
      for(int cell=2;cell<=4;++cell) expected.emplace(cell,cell+7);
      check_mapping(11,4,expected,grid);
    } else {
      for(int cell=1;cell<=4;++cell) expected.emplace(1,cell);
      check_mapping(4,1,expected,grid);
    }
  }
}

int main(int argc,char** argv) {
  Loci::Init(&argc,&argv);
  doctest::Context context; context.applyCommandLine(argc,argv);
  context.setOption("abort-after",0);
  const int local_failure=context.run()==0 ? 0:1; int failure=0;
  MPI_Allreduce(&local_failure,&failure,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  Loci::Finalize();
  return failure;
}
