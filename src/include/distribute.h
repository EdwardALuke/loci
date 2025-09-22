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
#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <typeinfo>
#include <vector>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <string>
#include <memory>

#include <mpi.h>

#include <Map.h>
#include <DMap.h>
#include <store_rep.h>
#include <store_def.h>
#include <storeVec_def.h>
namespace Loci {
  class joiner ;
  extern std::ofstream debugout ;
  extern int MPI_processes;
  extern int MPI_rank ;
  
  void Init(int* argc, char*** argv) ;
  void Finalize() ; 
  void Abort() ;
  size_t MPI_process_mem_avail() ;
  


  /// Abstract base class to define data partioner interface
  class dataPartition {
  public:
    /// MPI Communicator that the data is partitioned over
    MPI_Comm comm ;
    /// Constructor to initialize communicator
    dataPartition(const MPI_Comm &icomm) : comm(icomm) {}
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


  /// Communication Schedule Generator for gathering data for a set of
  /// entities to each processor
  class gatherCommSchedule {
  public:
    MPI_Comm comm ;
    int p,r ;
    std::vector<int> send_entities ;
    std::vector<int> send_sizes, send_offsets ;
    std::vector<int> recv_sizes, recv_offsets ;
    // -----------------------------------------------------------------------
    /// @brief generateSchedule generates a communication schedule for
    /// gathering the entities in [request] distributed according to
    /// data partition [ptn]
    ///
    /// @param [request] Set of entities that is requested by each processor
    /// @param [ptn] Partition of entities across processors
    void generateSchedule(entitySet request, dataPartitionP ptn) ;
    // -----------------------------------------------------------------------
    /// @brief gatherData gathers the distributed data that will be accessed
    /// by each processor from a store that is partitioned according to ptn,
    /// 
    /// @param[result] store that will contain the gathered data for
    /// each processor, this will be ordered as the entitySet provided
    /// @param[data] The input const_store container that is distributed across
    /// processors
    template<class T> void gatherData(store<T> &result,
                                      const_store<T> &data) const {
      int result_size = 0 ;
      for(int i=0;i<p;++i)
        result_size += recv_sizes[i] ;

      entitySet dom  ;
      if(result_size > 0)
        dom = interval(0,result_size-1) ;
      result.allocate(dom) ;

      // If single processor run, just copy data directly to result
      if(p == 1) {
        for(size_t i=0;i<send_entities.size();++i) 
          result[i] = data[send_entities[i]] ;
        return ;
      }

      // Allocate buffer for data to be sent
      int send_sz = send_entities.size() ;
      std::vector<T> send_data(send_sz) ;

      // copy data to send buffer
      for(int i=0;i<send_sz;++i) 
        send_data[i] = data[send_entities[i]] ;

      // copy self data to result buffer first (don't send using MPI)
      for(int i=0;i<send_sizes[r];++i)
        result[recv_offsets[r]+i] =
          send_data[send_offsets[r]+i] ;

      // Count how many sends and recvs
      int reqs = 0 ;
      for(int i=0;i<p;++i) {
        if(i!=r && send_sizes[i] > 0)
          reqs++ ;
        if(i!=r && recv_sizes[i] > 0)
          reqs++ ;
      }
      if(reqs == 0) // if no comms, we can return
        return ;

      // requests for outstanding communication
      std::vector<MPI_Request> requests(reqs) ;
      int cnt = 0 ; // Count of communication operations (send/recv)

      MPI_Datatype bytearray ;
      MPI_Type_vector(sizeof(T),1,1,MPI_BYTE,&bytearray) ;
      MPI_Type_commit(&bytearray) ;
      // Perform non-blocking recvs
      for(int i=0;i<p;++i)
        // If I am not recving from myself and I have something to recv
        if(i!=r && recv_sizes[i] > 0) {
          MPI_Irecv(&result[recv_offsets[i]],
                    recv_sizes[i],
                    bytearray,i,99,comm,&requests[cnt]) ;
          cnt++ ;
        }
      // Perform non-blocking sends
      for(int i=0;i<p;++i)
        if(i!=r && send_sizes[i] > 0) {
          // If I am not sending to myself and I have something to send
          MPI_Isend(&send_data[send_offsets[i]],
                    send_sizes[i],
                    bytearray,i,99,comm,&requests[cnt]) ;
          cnt++ ;
        }

      // wait for communication to finish
      std::vector<MPI_Status> status_list(cnt) ;
      MPI_Waitall(cnt,&requests[0],&status_list[0]) ;

      MPI_Type_free(&bytearray) ;
    }

    // -----------------------------------------------------------------------
    /// @brief gatherData compatability conversion that allows const references
    /// to store to work the same as const_store

    template<class T> inline void gatherData(store<T> &result,
                                             const store<T> &data) const {
      const_store<T> data_const ;
      data_const.setRep(data.Rep()) ;
      gatherData(result,data_const) ;
    }
        
    // -----------------------------------------------------------------------
    /// @brief scatterData scatters the local buffer to the processor that
    /// owns the associated entity
    ///
    /// @param[buffer] const_store that will contain the data to be 
    /// scattered, the inverse of the gather operation
    /// @param[data] The output store container that is distributed across
    /// processors
    template<class T> void scatterData(const_store<T> &buffer,
                                       store<T> &data) const {

      // If single processor run, just compy data directly to result
      if(p == 1) {
        warn(buffer.domain().size() != send_entities.size()) ;
        for(size_t i=0;i<send_entities.size();++i) 
          data[send_entities[i]] = buffer[i] ;
        return ;
      }

      // Allocate buffer for data to be sent
      int recv_sz = send_entities.size() ;
      std::vector<T> recv_data(recv_sz) ;

      // copy self data to result buffer first (don't send using MPI)
      for(size_t i=0;i<send_sizes[r];++i)
        recv_data[send_offsets[r]+i] = buffer[recv_offsets[r]+i] ;

      // Count how many sends and recvs
      int reqs = 0 ;
      for(int i=0;i<p;++i) {
        if(send_sizes[i] > 0)
          reqs++ ;
        if(recv_sizes[i] > 0)
          reqs++ ;
      }
      // If this processor does not commuicate we are done
      if(reqs == 0) {
        for(size_t i=0;i<send_entities.size();++i) 
          data[send_entities[i]] = recv_data[i] ;
        return ;
      }

      
      MPI_Datatype bytearray ;
      MPI_Type_vector(sizeof(T),1,1,MPI_BYTE,&bytearray) ;
      MPI_Type_commit(&bytearray) ;
      
      // requests for outstanding communication
      std::vector<MPI_Request> requests(reqs) ;

      int cnt = 0 ; // Count of communication operations (send/recv)

      // Perform non-blocking recvs
      for(int i=0;i<p;++i)
        if(i!=r && send_sizes[i] > 0) {
          // If I am not sending to myself and I have something to send
          MPI_Irecv(&recv_data[send_offsets[i]],
                    send_sizes[i],
                    bytearray,i,99,comm,&requests[cnt]) ;
          cnt++ ;
        }
      // Perform non-blocking sends
      for(int i=0;i<p;++i)
        // If I am not recving from myself and I have something to recv
        if(i!=r && recv_sizes[i] > 0) {
          MPI_Isend(&buffer[recv_offsets[i]],
                    recv_sizes[i],
                    bytearray,i,99,comm,&requests[cnt]) ;
          cnt++ ;
        }

      // otherwise wait for communication to finish
      std::vector<MPI_Status> status_list(cnt) ;
      MPI_Waitall(cnt,&requests[0],&status_list[0]) ;
      MPI_Type_free(&bytearray) ;

      // copy data from from recv buffer
      for(size_t i=0;i<send_entities.size();++i) 
        data[send_entities[i]] = recv_data[i] ;
    }
    // -----------------------------------------------------------------------
    /// @brief scatterData scatters compatibility routine that allows const
    /// references to buffer to work the same as const_store
    template<class T> inline void scatterData(const store<T> &buffer,
                                              store<T> &data) const {
      const_store<T> buffer_const ;
      buffer_const.setRep(buffer.Rep()) ;
      scatterData(buffer_const,data) ;
    }

    // -----------------------------------------------------------------------
    /// @brief gatherData gathers the distributed data that will be accessed
    /// by each processor from a store that is partitioned according to ptn,
    /// 
    /// @param[result] store that will contain the gathered data for
    /// each processor, this will be ordered as the entitySet provided
    /// @param[data] The input const_store container that is distributed across
    /// processors
    template<class T> void gatherData(storeVec<T> &result,
                                      const_storeVec<T> &data) const {
      int result_size = 0 ;
      for(int i=0;i<p;++i)
        result_size += recv_sizes[i] ;

      const int vs = data.vecSize() ;
      result.setVecSize(vs) ;
      entitySet result_set ;
      if(result_size > 0)
        result_set = interval(0,result_size-1) ;
      result.allocate(result_set) ;
      
      // If single processor run, just copy data directly to result
      if(p == 1) {
        fatal(result_size != (int)send_entities.size()) ;
        //        fatal(result_size != data.domain().size()) ;
        for(size_t i=0;i<send_entities.size();++i)
          result[i] = data[send_entities[i]] ;
        return ;
      }

      // Allocate buffer for data to be sent
      int send_sz = send_entities.size() ;
      std::vector<T> send_data(send_sz*vs) ;

      // copy data to send buffer
      for(int i=0;i<send_sz;++i)
        for(int j=0;j<vs;++j)
          send_data[i*vs+j] = data[send_entities[i]][j] ;

      // copy self data to result buffer first (don't send using MPI)
      for(size_t i=0;i<send_sizes[r];++i)
        for(int j=0;j<vs;++j) 
          result[recv_offsets[r]+i][j] =
            send_data[(send_offsets[r]+i)*vs+j] ;

      // Count how many sends and recvs
      int reqs = 0 ;
      for(int i=0;i<p;++i) {
        if(i!=r && send_sizes[i] > 0)
          reqs++ ;
        if(i!=r && recv_sizes[i] > 0)
          reqs++ ;
      }
      if(reqs == 0) // if no comms, we can return
        return ;

      // requests for outstanding communication
      std::vector<MPI_Request> requests(reqs) ;
      int cnt = 0 ; // Count of communication operations (send/recv)

      MPI_Datatype bytearray ;
      MPI_Type_vector(sizeof(T)*vs,1,1,MPI_BYTE,&bytearray) ;
      MPI_Type_commit(&bytearray) ;
      // Perform non-blocking recvs
      for(int i=0;i<p;++i)
        if(i!=r && recv_sizes[i] > 0) { // if recving from another processor
          MPI_Irecv(&result[recv_offsets[i]][0],
                    recv_sizes[i],
                    bytearray,i,99,comm,&requests[cnt]) ;
          cnt++ ;
        }
      // Perform non-blocking sends
      for(int i=0;i<p;++i)
        if(i!=r && send_sizes[i] > 0) { // if sending to another processor
          MPI_Isend(&send_data[send_offsets[i]*vs],
                    send_sizes[i],
                    bytearray,i,99,comm,&requests[cnt]) ;
          cnt++ ;
        }

      // wait for communication to finish
      std::vector<MPI_Status> status_list(cnt) ;
      MPI_Waitall(cnt,&requests[0],&status_list[0]) ;

      MPI_Type_free(&bytearray) ;
    }

    // -----------------------------------------------------------------------
    /// @brief gatherData compatability conversion that allows const references
    /// to storeVec to work the same as const_storeVec

    template<class T> inline void gatherData(storeVec<T> &result,
                                             const storeVec<T> &data) const {
      const_storeVec<T> data_const ;
      data_const.setRep(data.Rep()) ;
      gatherData(result,data_const) ;
    }
        
    // -----------------------------------------------------------------------
    /// @brief scatterData scatters the local buffer to the processor that
    /// owns the associated entity
    ///
    /// @param[buffer] const_storeVec that will contain the data to be 
    /// scattered, the inverse of the gather operation
    /// @param[data] The output storeVec container that is distributed across
    /// processors
    template<class T> void scatterData(const_storeVec<T> &buffer,
                                       storeVec<T> &data) const {
      warn(buffer.vecSize() != data.vecSize()) ;
      const int vs = buffer.vecSize() ;
      // If single processor run, just compy data directly to result
      if(p == 1) {
        warn(buffer.domain().size() != send_entities.size()) ;
        for(size_t i=0;i<send_entities.size();++i)
          data[send_entities[i]] = buffer[i] ;
        return ;
      }

      // Allocate buffer for data to be sent
      int recv_sz = send_entities.size() ;
      std::vector<T> recv_data(recv_sz*vs) ;

      // copy self data to result buffer first (don't send using MPI)
      for(size_t i=0;i<send_sizes[r];++i)
        for(int j=0;j<vs;++j)
          recv_data[(send_offsets[r]+i)*vs+j] = buffer[recv_offsets[r]+i][j] ;

      // Count how many sends and recvs
      int reqs = 0 ;
      for(int i=0;i<p;++i) {
        if(send_sizes[i] > 0)
          reqs++ ;
        if(recv_sizes[i] > 0)
          reqs++ ;
      }

      if(reqs == 0) { // No communication case
        // copy data from buffer to output store
        for(size_t i=0;i<send_entities.size();++i)
          for(int j=0;j<vs;++j)
            data[send_entities[i]][j] = recv_data[i*vs+j] ;
        return ;
      }

      
      MPI_Datatype bytearray ;
      MPI_Type_vector(sizeof(T)*vs,1,1,MPI_BYTE,&bytearray) ;
      MPI_Type_commit(&bytearray) ;
      
      // requests for outstanding communication
      std::vector<MPI_Request> requests(reqs) ;

      int cnt = 0 ; // Count of communication operations (send/recv)

      // Perform non-blocking recvs
      for(int i=0;i<p;++i)
        if(i!=r && send_sizes[i] > 0) {
          // If I am not sending to myself and I have something to send
          MPI_Irecv(&recv_data[send_offsets[i]*vs],
                    send_sizes[i],
                    bytearray,i,99,comm,&requests[cnt]) ;
          cnt++ ;
        }
      // Perform non-blocking sends
      for(int i=0;i<p;++i)
        // If I am not recving from myself and I have something to recv
        if(i!=r && recv_sizes[i] > 0) {
          MPI_Isend(&buffer[recv_offsets[i]][0],
                    recv_sizes[i],
                    bytearray,i,99,comm,&requests[cnt]) ;
          cnt++ ;
        }

      // otherwise wait for communication to finish
      std::vector<MPI_Status> status_list(cnt) ;
      MPI_Waitall(cnt,&requests[0],&status_list[0]) ;
      MPI_Type_free(&bytearray) ;

      // copy data from from recv buffer
      for(size_t i=0;i<send_entities.size();++i) 
        for(int j=0;j<vs;++j)
          data[send_entities[i]][j] = recv_data[i*vs+j] ;
    }
    // -----------------------------------------------------------------------
    /// @brief scatterData scatters compatibility routine that allows const
    /// references to buffer to work the same as const_storeVec
    template<class T> inline void scatterData(const storeVec<T> &buffer,
                                              storeVec<T> &data) const {
      const_storeVec<T> buffer_const ;
      buffer_const.setRep(buffer.Rep()) ;
      scatterData(buffer_const,data) ;
    }
  } ;


  dMap send_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;

  std::vector<dMap> send_global_map(Map &attrib_data, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  void fill_clone(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  
  storeRepP send_clone_non(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  std::vector<storeRepP> send_global_clone_non(storeRepP &sp , entitySet &out_of_dom,  std::vector<entitySet> &init_ptn) ;
  
  std::vector<std::pair<int, entitySet> >
  transpose_entitySet(const std::vector<std::pair<int,entitySet> > &in,
                      MPI_Comm comm) ;
  
  std::vector<entitySet>
  transpose_entitySet(const std::vector<entitySet>& in, MPI_Comm comm) ;
  std::vector<sequence>
  transpose_sequence(const std::vector<sequence>& in, MPI_Comm comm) ;

  std::vector<entitySet> all_collect_vectors(entitySet &e,MPI_Comm comm) ;
  std::vector<entitySet> all_collect_vectors(entitySet &e) ;

  entitySet distribute_entitySet(entitySet e,const std::vector<entitySet> &ptn) ;

  entitySet all_collect_entitySet(const entitySet &e) ;

  // This is equivalent to but more efficient than
  // collectSet(entitySet iset) {
  // dset = all_collect_entitySet(kset) ;
  // return dset & domain
  entitySet collectSet(const entitySet iset, const entitySet domain,
		       MPI_Comm comm) ;
  
  int GLOBAL_OR(int b, MPI_Comm comm=MPI_COMM_WORLD) ;
  int GLOBAL_AND(int b, MPI_Comm comm=MPI_COMM_WORLD) ;
  int GLOBAL_MAX(int b, MPI_Comm comm=MPI_COMM_WORLD) ;
  int GLOBAL_MIN(int b, MPI_Comm comm=MPI_COMM_WORLD) ;
  
  // We've added these back as they seem to be used
  // in the fuel cell program//from distribute.h //////
  Map distribute_global_map(Map &m, const std::vector<entitySet> &vset) ;
  Map distribute_gmap(Map &m, const std::vector<entitySet> &vset) ;

  // a data-structure used in communcating data among processes
  struct P2pCommInfo {
    int proc ;                  // process id to talk/listen to
    entitySet global_dom ;      // the set of entity to send/recv
                                // (in a pre-agreed global numbering)
    entitySet local_dom ;       // remapped dom (if indeed remapped,
                                // if not remapped, equals "dom"
    P2pCommInfo(int p, const entitySet& ge, const entitySet& le)
      :proc(p),global_dom(ge),local_dom(le) {}
  } ;
  // a utility function to remap an entityset
  entitySet
  remap_entitySet(const entitySet& es, const Map& remap) ;
  entitySet
  remap_entitySet(const entitySet& es, const dMap& remap) ;
  // remap of a sequence
  sequence
  remap_sequence(const sequence& seq, const dMap& remap) ;
  sequence
  remap_sequence(const sequence& seq, const Map& remap) ;

  // given a dom partition over an mpi communicator and a requesting
  // entity set on each process, this function returns the point 2 point
  // communication structure that fulfills the request on each process.
  // since the request is in a pre-agreed global numbering system,
  // the send_remap is used to remap the send entities from the
  // pre-agreed numbering scheme to another numbering system that the
  // senders can understand. if it is null, then no remapping is done.
  // similarly the recv_remap is used remap the recv entities set from
  // the pre-agreed numbering scheme to another number scheme so that
  // the receiving side can understand. if it is null, then no remapping
  // is done.
  void
  get_p2p_comm(const std::vector<entitySet>& ptn,
               const entitySet& request,
               const dMap* send_remap,
               const dMap* recv_remap,
               MPI_Comm comm,
               std::vector<P2pCommInfo>& send,
               std::vector<P2pCommInfo>& recv) ;

#ifdef DYNAMICSCHEDULING
  // given a communication structure on every process,
  // this function fulfills the data communication.
  // each process fills the "dst" store with the data that
  // are from the "src" store. the "dst" and "src" CAN be the same.
  // however it is required that the "dst" store will have the
  // space to receive the requested data.
  // all the remaps are provides to convert respective local
  // numbering to/from the pre-agreed global numbering. if they
  // are NULL (or any of them is NULL), then no corresponding
  // remapping is performed.
  //
  // Technically this function *should* work with any of the
  // Loci container type that we have. However some rarely used
  // ones (such as multiStore, etc.) may have some problems due
  // to illy-designed/programmed routines. we might want to
  // investigate and fix those in the future if they do cause
  // problems
  void
  fill_store(storeRepP src, const Map* src_pack,
             storeRepP dst, const dMap* dst_unpack,
             const std::vector<P2pCommInfo>& send,
             const std::vector<P2pCommInfo>& recv,
             MPI_Comm comm) ;
  // this version only remaps the domain in the communication
  // process, this is mainly for maps, their image is not remapped
  // it also returns the actually filled domain
  entitySet
  fill_store_omd(storeRepP src, const Map* src_pack,
                 storeRepP dst, const dMap* dst_unpack,
                 const std::vector<P2pCommInfo>& send,
                 const std::vector<P2pCommInfo>& recv,
                 MPI_Comm comm) ;
  // this version uses the fill_store to expand the dst store
  // given a request (in a pre-agreed global numbering) and a
  // domain partition, this
  // function will expand the dst store to include the requested
  // domains (the data are retrieved from the src store)
  void
  expand_store(storeRepP src,
               // src also needs unpack because the "request"
               // is made in the global numbering scheme and
               // the src needs to understand that in its
               // own local number
               const Map* src_pack, const dMap* src_unpack,
               storeRepP dst,
               // dst does not need pack (obviously because it
               // just receive data)
               const dMap* dst_unpack,
               const entitySet& request,
               const std::vector<entitySet>& src_ptn, MPI_Comm comm) ;

  // this version uses the pack_size(e,packed) version to ensure
  // the requested communicate domain is valid on each process
  // this function returns the actually filled domain
  entitySet
  fill_store2(storeRepP src, const Map* src_pack,
              storeRepP dst, const dMap* dst_unpack,
              const std::vector<P2pCommInfo>& send,
              const std::vector<P2pCommInfo>& recv,
              MPI_Comm comm) ;
  // this version fills a vector of dstS from corresponding srcS
  std::vector<entitySet>
  fill_store2(std::vector<storeRepP>& src, const Map* src_pack,
              std::vector<storeRepP>& dst, const dMap* dst_unpack,
              const std::vector<P2pCommInfo>& send,
              const std::vector<P2pCommInfo>& recv, MPI_Comm comm) ;
  entitySet
  expand_store2(storeRepP src,
                // src also needs unpack because the "request"
                // is made in the global numbering scheme and
                // the src needs to understand that in its
                // own local number
                const Map* src_pack, const dMap* src_unpack,
                storeRepP dst,
                // dst does not need pack (obviously because it
                // just receive data)
                const dMap* dst_unpack,
                const entitySet& request,
                const std::vector<entitySet>& src_ptn, MPI_Comm comm) ;
  
  std::vector<entitySet>
  expand_store2(std::vector<storeRepP>& src,
                // src also needs unpack because the "request"
                // is made in the global numbering scheme and
                // the src needs to understand that in its
                // own local number
                const Map* src_pack, const dMap* src_unpack,
                std::vector<storeRepP>& dst,
                // dst does not need pack (obviously because it
                // just receive data)
                const dMap* dst_unpack,
                const entitySet& request,
                const std::vector<entitySet>& src_ptn, MPI_Comm comm) ;
  
  // we also need to push_store function to push the
  // contents to their originating process, to be done later
  //void
  //push_store(...

  // reduction of store, everything is the same as in the above
  // "expand_store" function. the join_op is the associative
  // operator used for the reduction part. and the preconditions
  // are that the "dst" has the space allocated and has been
  // initialized properly.
  void
  reduce_store(storeRepP src, const Map* src_pack,
               storeRepP dst, const dMap* dst_unpack,
               const std::vector<P2pCommInfo>& send,
               const std::vector<P2pCommInfo>& recv,
               CPTR<joiner> join_op, MPI_Comm comm) ;
  void
  reduce_store(storeRepP src,
               // src also needs unpack because the "request"
               // is made in the global numbering scheme and
               // the src needs to understand that in its
               // own local number
               const Map* src_pack, const dMap* src_unpack,
               storeRepP dst,
               // dst does not need pack (obviously because it
               // just receive data)
               const dMap* dst_unpack,
               const entitySet& request,
               CPTR<joiner> join_op,
               const std::vector<entitySet>& dst_ptn, MPI_Comm comm) ;

#endif
  // this function provides a way to determine if an execution 
  // thread is the leading execution unit in the system.  for threads,
  // this is similar to determine if the calling process is ranked 0
  // in the MPI communication world.
  bool is_leading_execution();

  // these two functions are intended to provide a global atomic
  // region of execution (the ideal use of them is to use the 
  // Loci preprocessor to hide these calls from the users)
  void global_atomic_region_begin();
  void global_atomic_region_end();

  class atomic_region_helper {
  public:
    atomic_region_helper() { global_atomic_region_begin() ; }
    ~atomic_region_helper() { global_atomic_region_end() ; }
  } ;

  // this function generates a unique name on each process/thread
  std::string gen_local_name(const std::string& prefix,
                             const std::string& suffix="");


}





#endif
 
