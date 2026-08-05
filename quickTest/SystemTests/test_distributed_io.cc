#include <Loci.h>
#include <distribute_io.h>
#include <mpi_containerIO.h>

#include <cstdio>
#include <iostream>
#include <vector>

using namespace Loci;

namespace {

  /// Build the minimal metadata needed to map an owned entity to its file
  /// number during distributed container I/O.
  fact_db::distribute_infoP make_test_distribution(const entitySet &owned,
                                                    int file_number) {
    fact_db::distribute_infoP dist = new fact_db::distribute_info;
    dist->myid = MPI_rank;
    dist->isDistributed = 1;
    dist->my_entities = owned;
    dist->comp_entities = owned;
    dist->copy_total_size = 0;
    dist->xmit_total_size = 0;

    Map l2g;
    Map l2f;
    store<unsigned char> key_domain;
    l2g.allocate(owned);
    l2f.allocate(owned);
    key_domain.allocate(owned);

    FORALL(owned, entity) {
      l2g[entity] = entity;
      l2f[entity] = file_number;
      key_domain[entity] = 0;
    } ENDFORALL;

    dist->l2g = l2g.Rep();
    dist->l2f = l2f.Rep();
    dist->key_domain = key_domain.Rep();
    return dist;
  }

} // namespace

int main(int argc, char **argv) {
  Loci::Init(&argc, &argv);

  if(MPI_processes != 3) {
    if(MPI_rank == 0)
      std::cerr << "test_distributed_io requires exactly three MPI ranks"
                << std::endl;
    Loci::Finalize();
    return 1;
  }

  int local_failures = 0;

  /// Reproduce the empty-rank case with one stored value across three MPI
  /// ranks. The checks cover partitioning, serial-HDF5 output, a globally
  /// empty store, and redistribution back to the original owner.

  /// Dividing one file number among three ranks should produce one populated
  /// partition and two empty partitions without inventing file numbers.
  const std::vector<entitySet> partitions =
    Loci::simplePartition(42, 42, MPI_COMM_WORLD);
  entitySet partitioned = EMPTY;
  int partition_entries = 0;
  int empty_partitions = 0;
  for(std::vector<entitySet>::const_iterator ptn = partitions.begin();
      ptn != partitions.end(); ++ptn) {
    partitioned += *ptn;
    partition_entries += ptn->size();
    if(*ptn == EMPTY)
      ++empty_partitions;
  }
  if(partitions.size() != 3 ||
     partitioned != entitySet(interval(42, 42)) ||
     partition_entries != 1 || empty_partitions != 2) {
    std::cerr << "simplePartition did not preserve empty MPI partitions"
              << std::endl;
    local_failures = 1;
  }

  /// Start the only value on rank two and assign it a different file number so
  /// serial-HDF5 output must redistribute it while two ranks begin empty.
  Loci::use_parallel_io = false;
  const char *filename = "test_distributed_io.h5";
  const int source_rank = 2;
  const int local_entity = 17;
  const int file_number = 42;
  const entitySet owned = MPI_rank == source_rank
                            ? entitySet(interval(local_entity, local_entity))
                            : EMPTY;
  const std::vector<char> expected = {'p', 'l', 'a', 'n'};

  store<std::vector<char> > values;
  values.allocate(owned);
  if(MPI_rank == source_rank)
    values[local_entity] = expected;

  fact_db facts;
  facts.put_distribute_info(make_test_distribution(owned, file_number));

  hid_t file_id = Loci::hdf5CreateFile(filename, H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT);
  Loci::writeContainer(file_id, "values", values.Rep(), facts);

  /// A globally empty store should remain empty on disk rather than becoming
  /// a reversed interval.
  store<std::vector<char> > empty_values;
  empty_values.allocate(EMPTY);
  Loci::writeContainer(file_id, "empty_values", empty_values.Rep(), facts);

  Loci::hdf5CloseFile(file_id);
  MPI_Barrier(MPI_COMM_WORLD);

  /// Inspect the file-numbered value and empty domain directly on rank zero.
  if(MPI_rank == 0) {
    store<std::vector<char> > written;
    file_id = Loci::hdf5OpenFile(filename, H5F_ACC_RDONLY, H5P_DEFAULT,
                                MPI_COMM_SELF);
    Loci::readContainerRAW(file_id, "values", written.Rep(), MPI_COMM_SELF);

    if(written.domain() != entitySet(interval(file_number, file_number)) ||
       written[file_number] != expected) {
      std::cerr << "distributed container output did not preserve its value"
                << std::endl;
      local_failures = 1;
    }

    hid_t group_id = H5Gopen(file_id, "empty_values", H5P_DEFAULT);
    if(group_id < 0) {
      std::cerr << "empty distributed container group was not written"
                << std::endl;
      local_failures = 1;
    } else {
      entitySet empty_domain;
      Loci::HDF5_ReadDomain(group_id, empty_domain);
      H5Gclose(group_id);
      if(empty_domain != EMPTY) {
        std::cerr << "empty distributed container has a nonempty file domain"
                  << std::endl;
        local_failures = 1;
      }
    }

    Loci::hdf5CloseFile(file_id, MPI_COMM_SELF);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  /// Map the file-numbered value back to its original entity on rank two.
  store<std::vector<char> > restored;
  restored.allocate(owned);
  file_id = Loci::hdf5OpenFile(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  Loci::readContainer(file_id, "values", restored.Rep(), owned, facts);
  Loci::hdf5CloseFile(file_id);

  if(restored.domain() != owned ||
     (MPI_rank == source_rank && restored[local_entity] != expected)) {
    std::cerr << "distributed container input did not restore its value on rank "
              << MPI_rank << std::endl;
    local_failures = 1;
  }

  int failures = 0;
  MPI_Allreduce(&local_failures, &failures, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if(MPI_rank == 0)
    std::remove(filename);

  if(failures == 0 && MPI_rank == 0)
    std::cout << "SUCCESS!" << std::endl;

  Loci::Finalize();
  return failures == 0 ? 0 : 1;
}
