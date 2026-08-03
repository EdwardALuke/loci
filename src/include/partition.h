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
#ifndef PARTITION_H
#define PARTITION_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>
#include <data_traits.h>
#include <mpi.h>
#include <memory>
#include <vector>
#include <entitySet.h>
#include <Tools/debug.h>

namespace Loci {
  /// Abstract base class to define data partioner interface
  class dataPartition {
  public:
    /// MPI Communicator that the data is partitioned over
    MPI_Comm comm ;
    /// Constructor to initialize communicator
    dataPartition(const MPI_Comm &icomm) : comm(icomm) {}
    /// virtual distructor
    virtual ~dataPartition() {} ;

    /// Method that will partition an entity set to the owning processor.
    /// Returns a list of pairs of processor number and set that is owned
    /// by that processor
    virtual std::vector<std::pair<int,entitySet> > partitionEntitySet(entitySet set) = 0 ;
    /// Return set that is owned by processor i
    virtual entitySet getAllocation(int i) = 0 ;
  } ;

  /// The most general partition where the assignment of sets to processors is
  /// just an array of sets.  This will require p intersections to peerform
  /// the partitionEntitySet operation
  class dataPartitionGeneral: public dataPartition {
    /// Partition set array of size p
    std::vector<entitySet> ptn ;
  public:
    dataPartitionGeneral(const std::vector<entitySet> &iptn, const MPI_Comm &icomm):
    dataPartition(icomm), ptn(iptn) {
#ifdef DEBUG
      int p = 1 ;
      MPI_Comm_size(comm,&p) ;
      fatal(int(ptn.size()) != p) ;
      Loci::debugout << "ptn = " << endl ;
      for(int i=0;i<p;++i)
        Loci::debugout << i << " - " << ptn[i] << endl ;
#endif
    }
    std::vector<std::pair<int,entitySet> > partitionEntitySet(entitySet set) ;
    entitySet getAllocation(int i) ;
  } ;


  /// A partition that is provided by a set of split values.  It is assumed
  /// that the entities are assigned to processors by splitting a continuous
  /// interval at (p-1) split locations.  The constructor takes an array of
  /// p values which are the initial value that is owned by that processor
  class dataPartitionSplits: public dataPartition {
    std::vector<int> splits ;
  public:
    dataPartitionSplits(const std::vector<int> &isplits,const MPI_Comm &icomm):
    dataPartition(icomm),splits(isplits) {
      int p = 1 ;
      MPI_Comm_size(comm,&p) ;
      fatal(splits.size() != size_t(p)) ;
      splits[0] = std::numeric_limits<int>::lowest()+1 ;
      splits.push_back(std::numeric_limits<int>::max()-1) ;
#ifdef DEBUG
      Loci::debugout << "splits=" ;
      for(int i=0;i<=p;++i)
        Loci::debugout << " " << splits[i] ;
      Loci::debugout << endl ;
#endif
    }
    std::vector<std::pair<int,entitySet> > partitionEntitySet(entitySet set) ;
    entitySet getAllocation(int i) ;
  } ;


  /// The most efficient partition which assumes that ownership is assigned to
  /// processors with splits that are all equally spaced.  This allows for
  /// the intersection to be computed algorithmically. 
  class dataPartitionComputed: public dataPartition {
    int start, delta ;
  public:
    dataPartitionComputed(int istart,int idelta,const MPI_Comm &icomm) :
    dataPartition(icomm),start(istart),delta(idelta) {
#ifdef DEBUG
      Loci::debugout << "start=" << start << ", delta=" << delta << endl ;
#endif
    }
    std::vector<std::pair<int,entitySet> > partitionEntitySet(entitySet set) ;
    entitySet getAllocation(int i) ;
  } ;    



  /// A convenience pointer to the data partitioner interface
  typedef std::shared_ptr<dataPartition> dataPartitionP ;

  /// Factory functions for creating the general partition object, the
  /// input is a vector of non-overlapping entitySets for each processor
  /// where it is assumed that ptn[i] contains all entities owned by
  /// processor i.
  inline dataPartitionP createPartition(const std::vector<entitySet>  &ptn,
                                        const MPI_Comm &comm) {
    return std::make_shared<dataPartitionGeneral>(ptn,comm) ;
  }

  /// Factory function for creating the splitter based partition object, where
  /// the input is a vector of p values.  The splits[i] contains the first
  /// value that is owned by processor i and ownership of entities is
  /// described by the interval [splits[i],splits[i+1])
  inline dataPartitionP createPartition(const std::vector<int>  &splits,
                                        const MPI_Comm &comm) {
    return std::make_shared<dataPartitionSplits>(splits,comm) ;
  }

  /// Factory function for creating an algorithmic partition which is described
  /// by the starting number of the entity set and a delta which is the number
  /// of entities assigned to each processor.  This is the most constrained
  /// partition, but also the most efficent for partitioning a general set.
  inline dataPartitionP createPartition(int start, int delta, 
                                        const MPI_Comm &comm) {
    return std::make_shared<dataPartitionComputed>(start,delta,comm) ;
  }


}
#endif
