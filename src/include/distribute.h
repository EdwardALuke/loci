//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#include <mpi.h>

#include <Map.h>
#include <DMap.h>
#include <store_rep.h>

namespace Loci {
  class joiner ;
  extern std::ofstream debugout ;
  extern int MPI_processes;
  extern int MPI_rank ;
  
  void Init(int* argc, char*** argv) ;
  void Finalize() ; 
  void Abort() ;
  size_t MPI_process_mem_avail() ;
  
  dMap send_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;

  std::vector<dMap> send_global_map(Map &attrib_data, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  void fill_clone(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  
  storeRepP send_clone_non(storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  std::vector<storeRepP> send_global_clone_non(storeRepP &sp , entitySet &out_of_dom,  std::vector<entitySet> &init_ptn) ;
  
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
 
