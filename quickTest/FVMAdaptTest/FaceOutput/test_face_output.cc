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

#include <cstdio>
#include <iostream>
#include <string>

#include <Loci.h>
#include <FVMAdapt/defines.h>
#include <hdf5.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci;

/// Writing faces with empty MPI ranks must not add phantom faces or cells.
/// Supply the writer's inputs directly so grid reading cannot mask this bug.
TEST_CASE("face output preserves the mesh counts when two ranks own no faces") {
  const char* filename = "faces.vog";
  // Seed just the file metadata consumed by face_output. Node generation and
  // the mesh reader are outside this test of the face-writing stage.
  int header_ok = 1;
  if(MPI_rank == 0) {
    hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t group = H5Gcreate2(file, "file_info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(group, "numNodes", H5T_NATIVE_LLONG, space,
                          H5P_DEFAULT, H5P_DEFAULT);
    const long long num_nodes = 12;
    header_ok = H5Awrite(attr, H5T_NATIVE_LLONG, &num_nodes) >= 0;
    H5Aclose(attr);
    H5Sclose(space);
    H5Gclose(group);
    Loci::writeVOGClose(file);
  }
  MPI_Bcast(&header_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  REQUIRE(header_ok);

  fact_db facts;
  store<FineFaces> fine_faces, fine_faces_cell;
  fine_faces.allocate(MPI_rank == 0 ? entitySet(interval(0, 11)) : EMPTY);
  fine_faces_cell.allocate(EMPTY);

  // Three disconnected tetrahedra, all their boundary faces on rank zero.
  // Three cells keep cell coloring nonempty on each rank, isolating empty
  // face ownership from the separate fewer-cells-than-ranks problem.
  if(MPI_rank == 0) {
    for(int cell = 0; cell < 3; ++cell) {
      const int n = 4 * cell;
      // Each record contains: left cell, boundary ID, then node IDs.
      fine_faces[n]     = {{cell, -1, n, n+2, n+1}};
      fine_faces[n + 1] = {{cell, -1, n, n+1, n+3}};
      fine_faces[n + 2] = {{cell, -1, n+1, n+2, n+3}};
      fine_faces[n + 3] = {{cell, -1, n+2, n, n+3}};
    }
  }
  facts.create_fact("fine_faces", fine_faces);
  facts.create_fact("fine_faces_cell", fine_faces_cell);
  param<std::string> outfile;
  *outfile = filename;
  facts.create_fact("outfile_par", outfile);
  // Query the installed module; let Loci select and execute the writer.
  rule_db rules;
  Loci::load_module("fvmadapt", rules);
  REQUIRE(Loci::makeQuery(rules, facts, "face_output"));

  MPI_Barrier(MPI_COMM_WORLD);
  if(MPI_rank == 0) {
    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t group = H5Gopen2(file, "file_info", H5P_DEFAULT);
    for(const char* name : {"numFaces", "numCells"}) {
      CAPTURE(name);
      hid_t attr = H5Aopen(group, name, H5P_DEFAULT);
      long long count = 0;
      CHECK(H5Aread(attr, H5T_NATIVE_LLONG, &count) >= 0);
      CHECK(count == (std::string(name) == "numFaces" ? 12 : 3));
      H5Aclose(attr);
    }
    H5Gclose(group);
    H5Fclose(file);
  }
}

int main(int argc, char** argv) {
  Loci::Init(&argc, &argv);
  if(MPI_processes != 3) {
    if(MPI_rank == 0)
      std::cerr << "test_face_output requires three MPI ranks\n";
    Loci::Finalize();
    return 1;
  }
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  context.setOption("abort-after", 0);
  const int local_failure = context.run() == 0 ? 0 : 1;
  int failure = 0;
  MPI_Allreduce(&local_failure, &failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if(failure == 0 && MPI_rank == 0)
    std::remove("faces.vog");
  Loci::Finalize();
  return failure;
}
