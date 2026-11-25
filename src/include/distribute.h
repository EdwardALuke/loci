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
#include <partition.h>
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
      if(cnt > 0)
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
      for(int i=0;i<send_sizes[r];++i)
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
      if(cnt > 0)
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
      for(int i=0;i<send_sizes[r];++i)
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
      if(cnt > 0)
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
      for(int i=0;i<send_sizes[r];++i)
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
      if(cnt > 0)
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


  // send_map and send_global_map 
  inline  dMap send_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn)  {
    cerr << "send_map is deprecated" << endl ;
    Loci::Abort() ;
    dMap m ;
    return m ;
  }

  inline  std::vector<dMap> send_global_map(Map &attrib_data, entitySet &out_of_dom, std::vector<entitySet> &init_ptn)  {
    cerr << "send_global_map is deprecated" << endl ;
    Loci::Abort() ;
    std::vector<dMap> m ;
    return m ;
  }
  inline storeRepP send_clone_non(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    cerr << "send_clone_non is deprecated" << endl ;
    Loci::Abort() ;
    return 0 ;
  }
  
  void fill_clone(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  
  
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


  /// This class describes the distribution of a container that exists in
  /// file numbering, before a partitioning has been established and before
  /// global numbering exists.  Containers in this state belong to a single
  /// keyspace.  The file partition describes how file numbers are distributed
  /// across processors.
  class filePartitionInfo {
    /// Entities owned by this processor in file numbering
    entitySet my_entities ;
    /// Partition describing how file numbering is distributed
    dataPartitionP filePartition ;
    /// Key space this partition describes
    int keyspace ;
    /// MPI communicator
    MPI_Comm comm ;
  public:
    filePartitionInfo(entitySet my_file_entities,
                      dataPartitionP filePtn,
                      int keyspace_in,
                      MPI_Comm comm_in) ;

    /// Return entities this processor owns in file numbering
    entitySet myEntities() const { return my_entities ; }
    /// Return the MPI communicator
    MPI_Comm getCommunicator() const { return comm ; }
    /// Return the key space this partition describes
    int getKeySpace() const { return keyspace ; }
    /// Return the file numbering partition
    dataPartitionP getPartition() const { return filePartition ; }

    /// Redistribute container from its current file-numbering distribution
    /// to a new distribution described by the target partition.
    /// This is a collective operation on the stored MPI communicator.
    ///
    /// @param[in] sp store representation for the input container
    /// @param[in] dom file numbering of domain to redistribute
    /// @param[in] target_ptn partition describing the desired distribution
    storeRepP redistribute(storeRepP sp, entitySet dom,
                           dataPartitionP target_ptn) ;
  } ;


  /// This class describes the distribution of a container that exists in
  /// global numbering, before local numbering has been created (i.e. before
  /// category analysis).  Containers in this state belong to a single
  /// keyspace.  The global-to-file map allows conversion between global
  /// and file numbering.  Contains a filePartitionInfo describing the
  /// underlying file numbering distribution for this keyspace.
  class globalPartitionInfo {
    /// Entities owned by this processor in global numbering
    entitySet my_entities ;
    /// Global-to-file numbering map for this keyspace
    dMap g2f ;
    /// Partition describing how global numbering is distributed
    dataPartitionP globalPartition ;
    /// Contained file partition info for this keyspace
    filePartitionInfoP file_info ;
    /// Key space this partition describes
    int keyspace ;
    /// MPI communicator
    MPI_Comm comm ;
  public:
    globalPartitionInfo(entitySet my_global_entities,
                        dMap &g2f_in,
                        dataPartitionP globalPtn,
                        filePartitionInfoP file_info_in,
                        int keyspace_in,
                        MPI_Comm comm_in) ;

    /// Return entities this processor owns in global numbering
    entitySet myEntities() const { return my_entities ; }
    /// Return the MPI communicator
    MPI_Comm getCommunicator() const { return comm ; }
    /// Return the key space this partition describes
    int getKeySpace() const { return keyspace ; }
    /// Return the global numbering partition
    dataPartitionP getPartition() const { return globalPartition ; }
    /// Return the global-to-file map
    storeRepP getGlobal2FileMap() const { return g2f.Rep() ; }
    /// Return the contained file partition info
    filePartitionInfoP getFilePartitionInfo() const { return file_info ; }

    /// Convert container from global numbering to file numbering.
    /// This is a collective operation on the stored MPI communicator.
    ///
    /// @param[in] sp store representation for the input container
    /// @param[in] dom global numbering of domain to convert
    /// @param[out] offset offset in file numbering (each processor will
    ///   allocate from zero, add offset to domain to get actual file numbering)
    storeRepP Global2FileOrder(storeRepP sp, entitySet dom, int &offset) ;

    /// Convert container from file numbering to global numbering.
    /// This is a collective operation on the stored MPI communicator.
    ///
    /// @param[out] result resulting container in global numbering
    /// @param[in] resultSet set in global numbering that will be populated
    /// @param[in] input input container in file numbering
    /// @param[in] offset add this to input domain to obtain true file
    ///   numbering
    void File2GlobalOrder(storeRepP &result, entitySet resultSet,
                          storeRepP input, int offset) ;
  } ;


  /// This class is used to describe how entities are partitioned and
  /// distributed to processors.  It is assumed that the data will be
  /// stored in a local numbering, and a mapping from the local numbering
  /// to the global numbering is provided.  The file numbering is accessed
  /// through the contained globalPartitionInfo for each keyspace.
  /// Note that containers in local numbering may span multiple keyspaces.
  class  entityPartitionInfo {
    entityPartitionInfo() {}
    /// Entities owned by this processor in local numbering
    entitySet my_entities ;
    /// Map from local to global numbering
    Map l2g ;
    /// Key space for each local entity
    store<unsigned char> key_domain ;
    /// Map from local to file numbering
    Map l2f ;
    /// Distribution info for each key space (contains global and file info)
    std::vector<globalPartitionInfoP> keyspaceInfo ;
    /// MPI communicator for distribution
    MPI_Comm comm ;
  public:
    /// Construct from contained globalPartitionInfo objects per keyspace
    entityPartitionInfo(entitySet &scope,
                        Map &l2g_in,
                        store<unsigned char> &key_domain_in,
                        Map &l2f_in,
                        std::vector<globalPartitionInfoP> &keyspace_info_in,
                        MPI_Comm comm_in) ;

    /// Backward-compatible constructor: builds contained objects from raw
    /// maps and partitions.
    entityPartitionInfo(entitySet &scope,
                        Map &l2g_in,
                        store<unsigned char> &key_domain_in,
                        Map &l2f_in,
                        std::vector<dataPartitionP> &ptnlist,
                        MPI_Comm comm_in) ;

    /// return entities this processor owns in local numbering
    entitySet myEntities() const { return my_entities ; }
    /// return the MPI communicator that this distribution is over
    MPI_Comm getCommunicator() const { return comm ; }
    /// return map from local numbering to global numbering
    storeRepP getLocal2GlobalMap() const { return l2g.Rep() ; }
    /// return map from local numbering to file numbering
    storeRepP getLocal2FileMap() const { return l2f.Rep() ; }
    /// return the key domain for the input set in local numbering
    ///
    /// @param[in] dom The input set in local numbering
    int getKeyDomain(entitySet dom) const ;

    /// return the globalPartitionInfo for a given key space
    ///
    /// @param[in] kd the key space index
    globalPartitionInfoP getKeySpaceInfo(int kd) const {
      return keyspaceInfo[kd] ;
    }

    /// return the partition for a given key space in global numbering
    ///
    /// @param[in] kd the key space index
    dataPartitionP getKeySpacePartition(int kd) const {
      return keyspaceInfo[kd]->getPartition() ;
    }

    /// Convert container from local numbering to file numbering, return
    /// converted container.  Note, this is a collective operation on the
    /// stored MPI communicator
    ///
    /// @param[in] sp, store representation for the input store container
    /// @param[in] dom, local numbering of domain to convert to file numbering
    /// @param[out] offset, offset in file numbering (each processor will
    /// allocate from zero, add offset to domain to get actual file numbering)
    storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset) ;

    /// Convert container from file numbering to local numbering
    ///
    /// @param[out] result, resulting container in local numbering
    /// @param[in] resultSet, set in local numbering that will be converted
    /// @param[in] input, input container in file numbering
    /// @param[in] offset, add this to input domain to obtain true file
    /// numbering
    void File2LocalOrder(storeRepP &result,entitySet resultSet,
                         storeRepP input, int offset) ;

    /// Convert container from local numbering to global numbering, return
    /// converted container.  Note, this is a collective operation on the
    /// stored MPI communicator.  Domain must be in a single keyspace.
    ///
    /// @param[in] sp, store representation for the input store container
    /// @param[in] dom, local numbering of domain to convert to global numbering
    /// @param[out] offset, offset in global numbering (each processor will
    /// allocate from zero, add offset to domain to get actual global numbering)
    storeRepP Local2GlobalOrder(storeRepP sp, entitySet dom, int &offset) ;

    /// Convert container from global numbering to local numbering
    ///
    /// @param[out] result, resulting container in local numbering
    /// @param[in] resultSet, set in local numbering that will be populated
    /// @param[in] input, input container in global numbering
    /// @param[in] offset, add this to input domain to obtain true global
    /// numbering
    void Global2LocalOrder(storeRepP &result, entitySet resultSet,
                           storeRepP input, int offset) ;
  } ;
  

  storeRepP redistribute_container(storeRepP sp, std::vector<entitySet> &dom_ptn,
				  MPI_Comm comm) ;
}

#endif
 
