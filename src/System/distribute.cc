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
#include <constraint.h>
#include "dist_tools.h"
#include <Tools/debug.h>
#include <entitySet.h>

#include "Tools/debugger.h"
#include "Tools/hash_map.h"

#include <rule.h>

#include <vector>
using std::vector ;

#include <set>
using std::set ;

#include <mpi.h>

#include <iostream>
using std::endl ;

#include <algorithm>
using std::sort ;

#include "execute.h"

namespace Loci {

  std::vector<std::pair<int,entitySet> >
  dataPartitionGeneral::partitionEntitySet(entitySet set) {
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    if(p==1) {
      vector<std::pair<int,entitySet> > res ;
      res.push_back(make_pair(0,set)) ;
      return res ;
    }
    fatal(ptn.size() != p) ;
    std::vector<std::pair<int,entitySet> > splitlist ;
    for(int i=0;i<p;++i) {
      entitySet part = set&ptn[i] ;
      if(part != EMPTY) {
        splitlist.push_back(std::make_pair(i,part)) ;
      }
    }
#ifdef DEBUG
    entitySet sum ;
    for(auto &elem : splitlist)
      sum += elem.second ;
    fatal(sum!=set) ;
#endif
    return splitlist ;
  }

  entitySet dataPartitionGeneral::getAllocation(int i) {
#ifdef DEBUG
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    fatal(i<0 || i >= p) ;
#endif
    return ptn[i] ;
  }
  
  std::vector<std::pair<int,entitySet> >
  dataPartitionSplits::partitionEntitySet(entitySet set) {
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    if(p==1) {
      vector<std::pair<int,entitySet> > res ;
      res.push_back(make_pair(0,set)) ;
      return res ;
    }

    fatal(splits.size() != p+1) ;
    // remove any entity not defined by the partition function
    //    set -= interval(splits[0],splits[p]-1) ;
    entitySet filter = interval(splits[0],splits[p]-1) ;
    fatal((set - filter) != EMPTY) ;
    std::vector<std::pair<int,entitySet> > splitlist ;
    for(int i=0;i<p;++i) {
      if(set == EMPTY)
        break ;
      if(set.Min() < splits[i+1]) {
        interval ival(splits[i],splits[i+1]-1) ;
        entitySet part = set&ival ;
        //        set -= ival ;
        splitlist.push_back(std::make_pair(i,part)) ;
      }
    }
#ifdef DEBUG
    entitySet sum ;
    for(auto &elem : splitlist)
      sum += elem.second ;
    fatal(sum!=set) ;
#endif
    return splitlist ;
  }

  entitySet dataPartitionSplits::getAllocation(int i) {
#ifdef DEBUG
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    fatal(i<0 || i >= p) ;
#endif
    return entitySet(interval(splits[i],splits[i+1]-1)) ;
  }

  std::vector<std::pair<int,entitySet> >
  dataPartitionComputed::partitionEntitySet(entitySet set) {
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    if(p==1) {
      vector<std::pair<int,entitySet> > res ;
      res.push_back(make_pair(0,set)) ;
      return res ;
    }
    
    // remove any entity not defined by the partition function
    //    set -= interval(start,start+delta*p-1) ;
    fatal((set - interval(start,start+delta*p-1)) != EMPTY) ;
    std::vector<std::pair<int,entitySet> > splitlist ;
    while(set != EMPTY) {
      int i = (set.Min()-start)/delta ;
      int s1 = start+i*delta ;
      interval ival(s1,s1+delta-1) ;
      entitySet part = set&ival ;
      //      set -= ival ;
      splitlist.push_back(std::make_pair(i,part)) ;
    }
#ifdef DEBUG
    entitySet sum ;
    for(auto &elem : splitlist)
      sum += elem.second ;
    fatal(sum!=set) ;
#endif
    return splitlist ;
  }

  entitySet dataPartitionComputed::getAllocation(int i) {
#ifdef DEBUG
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    fatal(i<0 || i >= p) ;
#endif
    int s1 = start+i*delta ;
    return entitySet(interval(s1,s1+delta-1)) ;
  }

  void gatherCommSchedule::generateSchedule(entitySet request,
                                            dataPartitionP ptn) {
    comm = ptn->comm ;
    MPI_Comm_size(comm,&p) ;
    MPI_Comm_rank(comm,&r) ;
    // Divide set into parts based on processor that owns it
    Loci::debugout << "request = " << request << endl ;

    vector<pair<int,entitySet> > recv_sets = ptn->partitionEntitySet(request) ;
    Loci::debugout << "recv_sets = " << endl ;
    for(size_t i = 0 ;i< recv_sets.size() ;++i)
      Loci::debugout << recv_sets[i].first << " - " << recv_sets[i].second
                     << endl ;
    // Transpose this to get sets that we need to send to gather request
    vector<pair<int,entitySet> > send_sets =
      transpose_entitySet(recv_sets, comm) ;


    // Setup Information for Sending Side
    int send_sz = 0 ;
    { vector<int> tmp(p,0) ; send_sizes.swap(tmp) ; }
    for(auto i = send_sets.begin();i!=send_sets.end();++i) {
      send_sz += i->second.size() ;
      send_sizes[i->first] = i->second.size() ;
    }

    // Entities to gather from sending side
    { vector<int> tmp(send_sz) ; send_entities.swap(tmp) ; }
    int j = 0 ;
    for(auto i = send_sets.begin();i!=send_sets.end();++i) {
      FORALL(i->second,ii) {
        send_entities[j] = ii ;
        j++ ;
      } ENDFORALL ;
    }

    // Recv data sizes
    { vector<int> tmp(p,0) ; recv_sizes.swap(tmp) ; }
    for(auto i = recv_sets.begin();i!=recv_sets.end();++i)
      recv_sizes[i->first] = i->second.size() ;

    // Compute offsets into send and recv buffers
    { vector<int> tmp(p,0) ; send_offsets.swap(tmp) ; }
    { vector<int> tmp(p,0) ; recv_offsets.swap(tmp) ; }
    for(int i=1;i<p;++i) {
      send_offsets[i] = send_offsets[i-1]+send_sizes[i-1] ;
      recv_offsets[i] = recv_offsets[i-1]+recv_sizes[i-1] ;
    }

#ifdef DEBUG
    // Sanity Check
    size_t buff_size = 0 ;
    for(int i=0;i<p;++i)
        buff_size += send_sizes[i] ;
    warn(buff_size != send_entities.size()) ;
#endif

  }

  std::vector<std::pair<int, entitySet> >
  transpose_entitySet(const std::vector<std::pair<int,entitySet> > &in,
                      MPI_Comm comm) {
    int np ;
    MPI_Comm_size(comm, &np) ;
    
    
    vector<bool> pack_interval(np,false) ;
    // first compute the send count and displacement
    vector<int> send_counts(np,0) ;
    vector<int> send_loc(np,-1) ;
    for(size_t i=0;i<in.size();++i) {
      // since if sending intervals, the send size will
      // be 2 * the number of intervals, then if that is
      // larger than the total element size, then we choose
      // to communicate elements directly, otherwise, we
      // choose to send the intervals
      int num_intervals = in[i].second.num_intervals() ;
      num_intervals*=2 ;
      int size = in[i].second.size() ;
      int j = in[i].first ;
      send_loc[j] = i ;
      if(num_intervals >= size) {
        pack_interval[j] = false ;
        // +1 since we include a flag in the head to indicate
        // whether this is interval, or elements packed
        send_counts[j] = size + 1 ;
      } else {
        pack_interval[j] = true ;
        send_counts[j] = num_intervals + 1 ;
      }
    }

    vector<int> send_displs(np,0) ;
    send_displs[0] = 0 ;
    for(int i=1;i<np;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
    // then communicate this get the recv info.
    vector<int> recv_counts(np,0) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, comm) ;
    vector<int> recv_displs(np,0) ;
    recv_displs[0] = 0 ;
    for(int i=1;i<np;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    // all info. gathered, ready to do MPI_Alltoallv
    // first pack data into a raw buffer.
    int buf_size = send_counts[np-1] +
      send_displs[np-1] ;
    
    vector<int> send_buf(buf_size) ;
    int buf_idx = 0 ;
    for(int i=0;i<np;++i) {
      // only pack non-empty entitySet
      if(send_counts[i] == 0)
        continue ;

      const entitySet& eset = in[send_loc[i]].second ;
      if(pack_interval[i]) {
        // packing intervals
        send_buf[buf_idx++] = 1 ; // set flag to indicate
                                  // that this is intervals packed
        for(size_t k=0;k<eset.num_intervals();++k) {
          send_buf[buf_idx++] = eset[k].first ;
          send_buf[buf_idx++] = eset[k].second ;
        }
      } else {
        send_buf[buf_idx++] = 0 ; // elements directly packed
        for(entitySet::const_iterator ei=eset.begin();
            ei!=eset.end();++ei,++buf_idx)
          send_buf[buf_idx] = *ei ;
      }
    }
    // allocate receive buffer
    int recv_size = recv_displs[np-1] +
      recv_counts[np-1] ;

    vector<int> recv_buf(recv_size) ;
    // communicate
    MPI_Alltoallv(&send_buf[0], &send_counts[0],
                  &send_displs[0], MPI_INT,
                  &recv_buf[0], &recv_counts[0],
                  &recv_displs[0], MPI_INT, comm) ;
    // release buffers that are not needed
    vector<int>().swap(send_counts) ;
    vector<int>().swap(send_displs) ;
    vector<int>().swap(send_buf) ;

    // unpack recv buffer into a vector of entitySet
    vector<pair<int,entitySet> > out ;

    for(int i=0;i<np;++i) {
      int b = recv_displs[i] ;
      int e = b + recv_counts[i] ;
      if(b == e)
        continue ;              // empty buffer

      entitySet eset ;
      int flag = recv_buf[b] ;
      ++b ;
      if(flag == 1) {
        // if packed with interval
        for(;b<e;b+=2) {
          int l = recv_buf[b] ;
          int u = recv_buf[b+1] ;
          eset += interval(l,u) ;
        }
      } else {
        // packed with elements
        for(;b<e;++b) {
          eset += recv_buf[b] ;
        }
      }
      out.push_back(make_pair(i,eset)) ;
    }

    return out ;
    // the end...
  }

  
  std::vector<entitySet>
  transpose_entitySet(const std::vector<entitySet>& in, MPI_Comm comm) {
    int np ;
    MPI_Comm_size(comm, &np) ;
    
    
    vector<bool> pack_interval(np,false) ;
    // first compute the send count and displacement
    vector<int> send_counts(np,0) ;
    for(int i=0;i<np;++i) {
      // we don't pack empty entitySet
      if(in[i] == EMPTY) {
        send_counts[i] = 0 ;
        continue ;
      }
      // since if sending intervals, the send size will
      // be 2 * the number of intervals, then if that is
      // larger than the total element size, then we choose
      // to communicate elements directly, otherwise, we
      // choose to send the intervals
      int num_intervals = in[i].num_intervals() ;
      num_intervals*=2 ;
      int size = in[i].size() ;
      if(num_intervals >= size) {
        pack_interval[i] = false ;
        // +1 since we include a flag in the head to indicate
        // whether this is interval, or elements packed
        send_counts[i] = size + 1 ;
      } else {
        pack_interval[i] = true ;
        send_counts[i] = num_intervals + 1 ;
      }
    }
    vector<int> send_displs(np,0) ;
    send_displs[0] = 0 ;
    for(int i=1;i<np;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
    // then communicate this get the recv info.
    vector<int> recv_counts(np,0) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, comm) ;
    vector<int> recv_displs(np,0) ;
    recv_displs[0] = 0 ;
    for(int i=1;i<np;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    // all info. gathered, ready to do MPI_Alltoallv
    // first pack data into a raw buffer.
    int buf_size = send_counts[np-1] +
      send_displs[np-1] ;
    
    vector<int> send_buf(buf_size) ;
    int buf_idx = 0 ;
    for(int i=0;i<np;++i) {
      // only pack non-empty entitySet
      if(send_counts[i] == 0)
        continue ;
      
      const entitySet& eset = in[i] ;
      if(pack_interval[i]) {
        // packing intervals
        send_buf[buf_idx++] = 1 ; // set flag to indicate
                                  // that this is intervals packed
        for(size_t k=0;k<eset.num_intervals();++k) {
          send_buf[buf_idx++] = eset[k].first ;
          send_buf[buf_idx++] = eset[k].second ;
        }
      } else {
        send_buf[buf_idx++] = 0 ; // elements directly packed
        for(entitySet::const_iterator ei=eset.begin();
            ei!=eset.end();++ei,++buf_idx)
          send_buf[buf_idx] = *ei ;
      }
    }
    // allocate receive buffer
    int recv_size = recv_displs[np-1] +
      recv_counts[np-1] ;

    vector<int> recv_buf(recv_size) ;
    // communicate
    MPI_Alltoallv(&send_buf[0], &send_counts[0],
                  &send_displs[0], MPI_INT,
                  &recv_buf[0], &recv_counts[0],
                  &recv_displs[0], MPI_INT, comm) ;
    // release buffers that are not needed
    vector<int>().swap(send_counts) ;
    vector<int>().swap(send_displs) ;
    vector<int>().swap(send_buf) ;

    // unpack recv buffer into a vector of entitySet
    vector<entitySet> out(np, EMPTY) ;

    for(int i=0;i<np;++i) {
      int b = recv_displs[i] ;
      int e = b + recv_counts[i] ;
      if(b == e)
        continue ;              // empty buffer

      entitySet& eset = out[i] ;
      int flag = recv_buf[b] ;
      ++b ;
      if(flag == 1) {
        // if packed with interval
        for(;b<e;b+=2) {
          int l = recv_buf[b] ;
          int u = recv_buf[b+1] ;
          eset += interval(l,u) ;
        }
      } else {
        // packed with elements
        for(;b<e;++b) {
          eset += recv_buf[b] ;
        }
      }
    }

    return out ;
    // the end...
  }

  std::vector<sequence>
  transpose_sequence(const std::vector<sequence>& in, MPI_Comm comm) {
    int np ;
    MPI_Comm_size(comm, &np) ;
    
    vector<bool> pack_interval(np,false) ;
    // first compute the send count and displacement
    vector<int> send_counts(np,0) ;
    for(int i=0;i<np;++i) {
      if(in[i] == EMPTY) {
        send_counts[i] = 0 ;
        continue ;
      }
      // since if sending intervals, the send size will
      // be 2 * the number of intervals, then if that is
      // larger than the total element size, then we choose
      // to communicate elements directly, otherwise, we
      // choose to send the intervals
      int num_intervals = in[i].num_intervals() ;
      num_intervals*=2 ;
      int size = in[i].size() ;
      if(num_intervals >= size) {
        pack_interval[i] = false ;
        // +1 since we include a flag in the head to indicate
        // whether this is interval, or elements packed
        send_counts[i] = size + 1 ;
      } else {
        pack_interval[i] = true ;
        send_counts[i] = num_intervals + 1 ;
      }
    }
    vector<int> send_displs(np,0) ;
    send_displs[0] = 0 ;
    for(int i=1;i<np;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
    // then communicate this get the recv info.
    vector<int> recv_counts(np,0) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, comm) ;
    vector<int> recv_displs(np,0) ;
    recv_displs[0] = 0 ;
    for(int i=1;i<np;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    // all info. gathered, ready to do MPI_Alltoallv
    // first pack data into a raw buffer.
    int buf_size = send_counts[np-1] +
      send_displs[np-1] ;
    
    vector<int> send_buf(buf_size) ;
    int buf_idx = 0 ;
    for(int i=0;i<np;++i) {
      if(send_counts[i] == 0)
        continue ;              // do not pack for empty message
      
      const sequence& seq = in[i] ;
      if(pack_interval[i]) {
        // packing intervals
        send_buf[buf_idx++] = 1 ; // set flag to indicate
                                  // that this is intervals packed
        for(size_t k=0;k<seq.num_intervals();++k) {
          send_buf[buf_idx++] = seq[k].first ;
          send_buf[buf_idx++] = seq[k].second ;
        }
      } else {
        send_buf[buf_idx++] = 0 ; // elements directly packed
        for(sequence::const_iterator si=seq.begin();
            si!=seq.end();++si,++buf_idx)
          send_buf[buf_idx] = *si ;
      }
    }
    // allocate receive buffer
    int recv_size = recv_displs[np-1] +
      recv_counts[np-1] ;

    vector<int> recv_buf(recv_size) ;
    // communicate
    MPI_Alltoallv(&send_buf[0], &send_counts[0],
                  &send_displs[0], MPI_INT,
                  &recv_buf[0], &recv_counts[0],
                  &recv_displs[0], MPI_INT, comm) ;
    // release buffers that are not needed
    vector<int>().swap(send_counts) ;
    vector<int>().swap(send_displs) ;
    vector<int>().swap(send_buf) ;

    // unpack recv buffer into a vector of entitySet
    vector<sequence> out(np, EMPTY) ;
    
    for(int i=0;i<np;++i) {
      int b = recv_displs[i] ;
      int e = b + recv_counts[i] ;
      if(b==e)
        continue ;
      
      sequence& seq = out[i] ;

      int flag = recv_buf[b] ;
      ++b ;
      if(flag == 1) {
        // if packed with interval
        for(;b<e;b+=2) {
          int f = recv_buf[b] ;
          int s = recv_buf[b+1] ;
          seq += interval(f,s) ;
        }
      } else {
        // packed with elements
        for(;b<e;++b) {
          seq += recv_buf[b] ;
        }
      }
    }

    return out ;
    // the end...
  }

  using std::vector ;
  void fill_clone( storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {

    vector<entitySet> recv_req(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      if(i!=MPI_rank) 
        recv_req[i] = out_of_dom & init_ptn[i] ;

    // send the recieve requests
    vector<int> recv_count(MPI_processes,0) ;
    vector<int> send_count(MPI_processes,0) ;
    vector<int> send_displacement(MPI_processes,0) ;
    vector<int> recv_displacement(MPI_processes,0) ;
    for(int i=0;i<MPI_processes;++i)
      send_count[i] = recv_req[i].num_intervals() * 2 ;
    MPI_Alltoall(&send_count[0],1,MPI_INT, &recv_count[0], 1, MPI_INT,MPI_COMM_WORLD) ;

    
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i=1;i<MPI_processes;++i) {
      send_displacement[i] = send_displacement[i-1]+send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1]+recv_count[i-1] ;
    }
    int mp = MPI_processes-1 ;
    int send_sizes = send_displacement[mp]+send_count[mp] ;
    int recv_sizes = recv_displacement[mp]+recv_count[mp] ;


    vector<int> send_set_alloc(send_sizes+recv_sizes,0) ;
    int * send_set_buf = &send_set_alloc[0] ; 
    int * recv_set_buf = &send_set_alloc[send_sizes] ; 
    for(int i=0;i<MPI_processes;++i) {
      for(size_t j=0;j<recv_req[i].num_intervals();++j) {
        send_set_buf[send_displacement[i]+j*2  ] = recv_req[i][j].first ;
        send_set_buf[send_displacement[i]+j*2+1] = recv_req[i][j].second ;
      }
    }

    MPI_Alltoallv(send_set_buf, &send_count[0], &send_displacement[0] , MPI_INT,
		  recv_set_buf, &recv_count[0], &recv_displacement[0], MPI_INT,
		  MPI_COMM_WORLD) ;

    vector<entitySet> send_set(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i) {
      for(int j=0;j<recv_count[i]/2;++j) {
        int i1 = recv_set_buf[recv_displacement[i]+j*2  ] ;
        int i2 = recv_set_buf[recv_displacement[i]+j*2+1] ;
        send_set[i] += interval(i1,i2) ;
      }
    }
    
#ifdef DEBUG
    // Sanity check, no send set should be outside of entities we own
    entitySet mydom = init_ptn[MPI_rank] ;
    for(int i=0;i<MPI_processes;++i) {
      if((send_set[i] & mydom) != send_set[i]) {
        cerr << "problem with partitioning in fill_clone!" ;
        debugout << "send_set["<< i << "] = " << send_set[i]
                 << "not owned = " << (send_set[i]-mydom) << endl ;
      }
    }
#endif
    
    // Now that we know what we are sending and receiving (from recv_req and
    // send_set) we can communicate the actual information...

    // Compute sizes of sending buffers
    for(int i=0;i<MPI_processes;++i) 
      send_count[i] =  sp->pack_size(send_set[i]) ;

    // Get sizes needed for receiving buffers
    MPI_Alltoall(&send_count[0],1,MPI_INT, &recv_count[0], 1, MPI_INT,MPI_COMM_WORLD) ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i=1;i<MPI_processes;++i) {
      send_displacement[i] = send_displacement[i-1]+send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1]+recv_count[i-1] ;
    }

    send_sizes = send_displacement[mp]+send_count[mp] ;
    recv_sizes = recv_displacement[mp]+recv_count[mp] ;

    vector<unsigned char> store_comm_buf(send_sizes+1024+recv_sizes,0) ;
    unsigned char *send_store = &store_comm_buf[0] ; 
    unsigned char *recv_store = &store_comm_buf[send_sizes+1024] ; 

    for(int i=0;i<send_sizes;++i)
      send_store[i] = 0 ;
    

    for(int i = 0; i <  MPI_processes; ++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_displacement[i]], loc_pack, send_count[i],
               send_set[i]) ;
    }
    
    MPI_Alltoallv(send_store, &send_count[0], &send_displacement[0],
                  MPI_PACKED,
		  recv_store, &recv_count[0], &recv_displacement[0],
                  MPI_PACKED,
		  MPI_COMM_WORLD) ;

    for(int i = 0; i <  MPI_processes; ++i) {
      int loc_pack = 0 ;
      sp->unpack(&recv_store[recv_displacement[i]], loc_pack,recv_count[i],
                 sequence(recv_req[i])) ; 
    }
  }
  
  storeRepP send_clone_non( storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet tmp_dom ;
    for(int i = 0; i <  MPI_processes; ++i)
      tmp_dom += entitySet(recv_dom[i]) ;
    storeRepP tmp_sp = sp->new_store(tmp_dom) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      recv_count[i] =  tmp_sp->pack_size(entitySet(recv_dom[i])) ;
      size_recv += recv_count[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    }
    
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      tmp_sp->unpack(recv_store, loc_pack, size_recv, recv_dom[i]) ; 
    }
    
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_sp ;
  }

  std::vector<storeRepP> send_global_clone_non(storeRepP &sp , entitySet &out_of_dom,  std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    std::vector< sequence> recv_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet total_dom ;
    std::vector<entitySet> e_vec( MPI_processes) ;
    std::vector< storeRepP> tmp_sp( MPI_processes) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      e_vec[i] = entitySet(recv_dom[i]) ;
      tmp_sp[i] = sp->new_store(e_vec[i]) ;
      recv_count[i] =  tmp_sp[i]->pack_size(e_vec[i]) ;
      size_recv += recv_count[i] ;
      total_dom += e_vec[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    }
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  MPI_COMM_WORLD) ;  
    
    loc_pack = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      tmp_sp[i]->unpack(recv_store, loc_pack, size_recv, recv_dom[i]) ; 
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_store ;
    delete [] recv_store ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_sp ;
  }

  // This only used by the code for the min2noslip computation!  Is it really
  // needed?
#define MINNOSLIP
#ifdef MINNOSLIP
  
  dMap send_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    HASH_MAP(int, int) attrib_data ;
    entitySet dm_dom = dm.domain() ;
    for(ei = dm_dom.begin(); ei != dm_dom.end(); ++ei)
      attrib_data[*ei] = dm[*ei] ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    
    std::vector<HASH_MAP(int, int) > map_entities( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(attrib_data.find(*vi) != attrib_data.end())
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    HASH_MAP(int, int) hm ;
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; j+=2) {
	hm[recv_map[j]] = recv_map[j+1];
      }
    }
    dMap tmp_dm ;
    for(HASH_MAP(int, int)::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi) {
      tmp_dm[hmi->first] = hmi->second ;
    }
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_map ;
    delete [] recv_map ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return tmp_dm ;
  }
  
  std::vector<dMap> send_global_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn) {
    int *recv_count = new int[ MPI_processes] ;
    int *send_count = new int[ MPI_processes] ;
    int *send_displacement = new int[ MPI_processes] ;
    int *recv_displacement = new int[ MPI_processes] ;
    entitySet::const_iterator ei ;
    HASH_MAP(int, int) attrib_data ;
    entitySet dm_dom = dm.domain() ;
    for(ei = dm_dom.begin(); ei != dm_dom.end(); ++ei)
      attrib_data[*ei] = dm[*ei] ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy( MPI_processes), send_clone( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i <  MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    
    std::vector<HASH_MAP(int, int) > map_entities( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(attrib_data.find(*vi) != attrib_data.end())
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  MPI_processes; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    std::vector<HASH_MAP(int, int) > hm( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      HASH_MAP(int, int) tmp_hm ;
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; j+=2) {
	tmp_hm[recv_map[j]] = recv_map[j+1];
      }
      hm[i] = tmp_hm ;
    }
    std::vector<dMap> v_dm( MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i) {
      HASH_MAP(int, int) tmp_hm = hm[i] ; 
      dMap tmp_dm ;
      for(HASH_MAP(int, int)::const_iterator hmi = tmp_hm.begin(); hmi != tmp_hm.end(); ++hmi) 
	tmp_dm[hmi->first] = hmi->second ;
      v_dm[i] = tmp_dm.Rep() ;
    }
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_map ;
    delete [] recv_map ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return v_dm ;
  }
#endif

  // Collect entitities to a unified entitySet that is distributed across
  // processors according to the partition ptn.
  entitySet dist_collect_entitySet(entitySet inSet, const vector<entitySet> &ptn) {
    const int p = MPI_processes ;
    const int r = MPI_rank ;
#ifdef DEBUG
    if(r == 0) 
      cerr << "dist_collect_entitySet is depreciated" << endl ;
#endif
    entitySet retval = inSet & ptn[r] ;
    // Check for empty and universal set
    int sbits = ((inSet != EMPTY)?1:0)| ((retval != ptn[r])?2:0) ;
    int rbits = sbits ;
    MPI_Allreduce(&sbits, &rbits, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD) ;
    if((rbits & 1) == 0) // EMPTY set
      return EMPTY ;
    if((rbits & 2) == 0) // UNIVERSE set
      return ptn[r] ;
    inSet -= ptn[r] ;

    stopWatch s ;
    s.start() ;
    
    // Communicate the largest interval to eliminate redundant communication
    interval iset(1,-1) ;
    for(size_t i=0;i<retval.num_intervals();++i)
      if((iset.second-iset.first) < (retval[i].second-retval[i].first))
        iset = retval[i] ;
    vector<interval> rset(p) ;
    MPI_Allgather(&iset.first,1,MPI_2INT,&rset[0].first,1,MPI_2INT,MPI_COMM_WORLD) ;
    
    
    // Remove redundant info
    for(int i=0;i<p;++i)
      if(rset[i].second >= rset[i].first)
        inSet -= rset[i] ;

    
    s.start() ;
    // If the remainder is empty, then we are finished
    if(GLOBAL_AND((inSet==EMPTY),MPI_COMM_WORLD)) {
      return retval ;
    }

    // Remove second largest interval
    interval iset2(1,-1) ;
    for(size_t i=0;i<retval.num_intervals();++i)
      if(iset != retval[i] &&
         (iset2.second-iset2.first) < (retval[i].second-retval[i].first))
        iset2 = retval[i] ;
    
    MPI_Allgather(&iset2.first,1,MPI_2INT,&rset[0].first,1,MPI_2INT,MPI_COMM_WORLD) ;
    
    
    // Remove redundant info
    for(int i=0;i<p;++i)
      if(rset[i].second >= rset[i].first)
        inSet -= rset[i] ;

    
    // If the remainder is empty, then we are finished
    if(GLOBAL_AND((inSet==EMPTY),MPI_COMM_WORLD)) {
      return retval ;
    }
    

    // Otherwise, we need to communicate residual entities to owning processor
    vector<entitySet> distSet(p) ;
    for(int i=0;i<p;++i) 
      distSet[i] = inSet & ptn[i] ;

    vector<int> send_sz(p) ;
    for(int i=0;i<p;++i)
      send_sz[i] = distSet[i].num_intervals() ;
    
    vector<int> recv_sz(p) ;
    MPI_Alltoall(&send_sz[0],1,MPI_INT,
                 &recv_sz[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;

    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<p;++i) {
      size_send += send_sz[i] ;
      size_recv += recv_sz[i] ;
    }

    vector<interval> send_store(size_send) ;
    vector<interval> recv_store(size_recv) ;
    vector<int> send_displacement(p) ;
    vector<int> recv_displacement(p) ;
    
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sz[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sz[i-1] ;
    }
    for(int i = 0; i <  p; ++i)
      for(size_t j=0;j<distSet[i].num_intervals();++j) {
        send_store[send_displacement[i]+j] = distSet[i][j] ;
      }
    
    MPI_Alltoallv(&send_store[0],&send_sz[0], &send_displacement[0], MPI_2INT,
    		  &recv_store[0],&recv_sz[0], &recv_displacement[0], MPI_2INT,
    		  MPI_COMM_WORLD) ;

    // Sort received intervals to increase performance of union code
    sort(recv_store.begin(),recv_store.end()) ;

    for(size_t i=0;i<recv_store.size();++i)
      retval += recv_store[i] ;
    
    return retval ;
  }

  // -----------------------------------------------------------------------
  /// @brief dist_expand_entitySet() is given an input set and a set of
  /// entities that are in the clone region.  The ownership of entities
  /// are given by the partition, it will recieve the data from the clone
  /// regions so that the expanded set is valid representation of the
  /// distributed set in the provided clone regions
  ///
  /// @param [inSet] input set that represents the set distributed across
  /// processors
  /// @param [clone] input set of cloned entities that we wish to fill in
  /// with set values from the remote procesors
  /// @param [ptn]   partition function that describes which processor owns
  /// which entities.
  entitySet dist_expand_entitySet(entitySet inSet, entitySet clone,
                                  const vector<entitySet> &ptn) {
    // We don't need to send data to ourselves
    clone -= inSet ;

    // Send an interval that contains the cloned region for each processor
    // If there is no entities on the processor, mark by an invalid interval
    vector<int> send_req(ptn.size()*2) ;
    for(size_t i=0;i < ptn.size();++i) {
      entitySet clone_i = clone&ptn[i] ;
      if(clone_i != EMPTY) {
        send_req[i*2+0] = clone_i.Min() ;
        send_req[i*2+1] = clone_i.Max() ;
      } else {
        send_req[i*2+0] = 1 ;
        send_req[i*2+1] = 0 ;
      }
    }

    // Send clone requests to partners
    vector<int> recv_req(ptn.size()*2) ;
    MPI_Alltoall(&send_req[0],2,MPI_INT, &recv_req[0],2,MPI_INT,
                 MPI_COMM_WORLD) ;

    // From the clone requests, find a set that we need to send to the
    // corresponding processor
    vector<int> send_sizes(ptn.size()) ;
    vector<entitySet> send_data(ptn.size()) ;
    for(size_t i=0;i < ptn.size();++i) {
      if(recv_req[i*2+0]<=recv_req[i*2+1]) 
        send_data[i] = inSet & interval(recv_req[i*2+0],recv_req[i*2+1]) ;
      send_sizes[i] = send_data[i].num_intervals() ;
    }

    // Share the size information with partners
    vector<int> recv_sizes(ptn.size()) ;
    MPI_Alltoall(&send_sizes[0],1,MPI_INT, &recv_sizes[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;

    // Compress all of the recieves to a single buffer, compute offsets
    // within each buffer for the recieves
    vector<int> recv_offsets(ptn.size()+1) ;
    recv_offsets[0] = 0 ;
    int numRecvs = 0 ;
    for(size_t i=0;i<ptn.size();++i) {
      recv_offsets[i+1] = recv_offsets[i]+recv_sizes[i] ;
      if(recv_sizes[i] > 0)
        numRecvs++ ;
    }
    vector<int> recv_buffer(recv_offsets[ptn.size()]*2,0) ;

    // Allocate request buffers for Irecvs
    vector<MPI_Request> recv_Requests(numRecvs) ;
    int cnt = 0 ;
    for(size_t i=0;i < ptn.size();++i) {
      if(recv_sizes[i] > 0) {
        int recv_size = recv_sizes[i]*2 ;
        int recv_offset = recv_offsets[i]*2 ;
        
        MPI_Irecv(&(recv_buffer[recv_offset]),recv_size,MPI_INT,i,3,
                  MPI_COMM_WORLD,&recv_Requests[cnt++]) ;
      }
    }

    // Send the data to partner
    for(size_t i=0;i < ptn.size();++i) {
      if(send_sizes[i] != 0) {
        MPI_Send(&(send_data[i][0]),send_sizes[i]*2,MPI_INT,i,3,MPI_COMM_WORLD) ;
      }
    }

    // Wait for communication to complete
    vector<MPI_Status> recv_Status(recv_Requests.size()) ;
    MPI_Waitall(numRecvs,&recv_Requests[0],&recv_Status[0]) ;

    // update set to add results from clone regions of other processors
    for(size_t i=0;i<recv_buffer.size();i=i+2)
      inSet += interval(recv_buffer[i+0],recv_buffer[i+1]) ;

    // return modified set
    return inSet ;
  }

  
  inline bool spec_ival_compare(const interval &i1,
                            const interval &i2) {
    if(i1.first < i2.first)
      return true ;
    if(i1.first == i2.first && i1.second > i2.second)
      return true ;
    return false ;
  }
      
  // Return union of all entitySets from all processors
  entitySet all_gather_entitySet(const entitySet &e) {
    const int p = MPI_processes ;
    if(p == 1)
      return e ;
    
    int send_count = 2*e.num_intervals() ;
    vector<int> recv_count(p) ;
    MPI_Allgather(&send_count,1,MPI_INT,&recv_count[0],1,MPI_INT,
                  MPI_COMM_WORLD) ;
    vector<int> recv_disp(p) ;
    recv_disp[0] = 0 ;
    for(int i=1;i<p;++i)
      recv_disp[i] = recv_disp[i-1]+recv_count[i-1] ;
    int tot = recv_disp[p-1]+recv_count[p-1] ;
    if(tot == 0)
      return EMPTY ;
    vector<interval> ivl_list(tot/2) ;
    vector<interval> snd_list(send_count/2) ;
    for(int i=0;i<send_count/2;++i)
      snd_list[i] = e[i] ;
    MPI_Allgatherv(&(snd_list[0]),send_count,MPI_INT,
                   &(ivl_list[0]),&(recv_count[0]), &(recv_disp[0]), MPI_INT,
                   MPI_COMM_WORLD) ;
    sort(ivl_list.begin(),ivl_list.end(),spec_ival_compare) ;
    entitySet tmp = ivl_list[0] ;
    for(size_t i=1;i<ivl_list.size();++i)
      tmp += ivl_list[i] ;
    return tmp ;
  }

  // -----------------------------------------------------------------------
  /// @brief distribute_entitySet() will send data to the owning processor
  /// as described by the ptn argument.  The entities from each processor are
  /// unioned to form the final entitySet that is partititoned to processors
  ///
  /// @param [e]  input set that needs to be distributed to owning processor
  /// @param [ptn]   partition function that describes which processor owns
  /// which entities.
  entitySet distribute_entitySet(entitySet e,const vector<entitySet> &ptn) {

    const int p = MPI_processes ;
    // Single processor, do nothing
    if(p == 1)
      return e ;

    // if any processor contains universal set, all processors hold the
    // universal set
    if(GLOBAL_OR(e == ~EMPTY)) 
      return ~EMPTY ;

    // if no communication need, then just return the input set
    if(GLOBAL_AND((e & ptn[MPI_rank]) == e)) 
      return e ;

    // Use recursive splitting to get set to owning processor.
    const int r = MPI_rank ;


    int nump = p ; // number of processors in the current group
    int start = 0 ; // starting processor number for the current group
    while(nump > 1) {
      int split=(nump+1)/2 ; // how we will split the current group

      // The start of the interval after the split
      int nstart = start ;
      // The number of processors after the split
      int nnump = split ;
      if(r >= start+split) {
        // If this processor is in the upper half, then adjust the
        // interval to be the upper half.
        nnump = nump-split ;
        nstart = start+split ;
      }
      // The rank of the last processor if nump is odd
      int oddproc = -1 ;
      if((nump&1) == 1) {
        oddproc = start+split-1 ;
      }

      // Find what part of the partition represents our subset
      // of processors
      entitySet myset ;
      for(int i=0;i<nnump;++i)
        myset += ptn[i+nstart] ;
      // divide the set into parts we keep and parts that belong
      // to the other half.
      entitySet keep = e&myset ;
      entitySet give = e-keep ;
      // If the split doesn't perfectly divide, then send remainder to
      // next lower processor so it gets shared with the other partition
      if(r == oddproc) {
        int set_ivls = give.num_intervals() ;
        MPI_Send(&set_ivls,1,MPI_INT,r-1,11,MPI_COMM_WORLD);
        
        if(set_ivls > 0) {
          vector<interval> ivl_list(set_ivls) ;
          for(int i=0;i<set_ivls;++i)
            ivl_list[i] = give[i] ;
          MPI_Send(&ivl_list[0],2*set_ivls,MPI_INT,r-1,12,MPI_COMM_WORLD) ;
        }
      }
      if(r+1==oddproc) {
        int set_ivls = 0 ;
        MPI_Status status ;
        MPI_Recv(&set_ivls,1,MPI_INT,r+1,11,MPI_COMM_WORLD,&status) ;

        if(set_ivls > 0) {
          vector<interval> ivl_list(set_ivls) ;
          MPI_Recv(&ivl_list[0],2*set_ivls,MPI_INT,r+1,12,MPI_COMM_WORLD,
                   &status) ;

          for(int i=0;i<set_ivls;++i)
            give += ivl_list[i] ;
        }
      }

      if(r != oddproc) { // exclude odd processor from this comm step
        
        
        if(r < start+split) { // Lower Partition
          int partner = r+split ;
          // Lower partition sends first,  then recieves
          int set_ivls = give.num_intervals() ;
          MPI_Send(&set_ivls,1,MPI_INT,partner,11,MPI_COMM_WORLD);

          if(set_ivls > 0) {
            vector<interval> ivl_list(set_ivls) ;
            for(int i=0;i<set_ivls;++i)
              ivl_list[i] = give[i] ;
            MPI_Send(&ivl_list[0],2*set_ivls,MPI_INT,partner,12,
                     MPI_COMM_WORLD) ;
          }
          MPI_Status status ;
          MPI_Recv(&set_ivls,1,MPI_INT,partner,11,MPI_COMM_WORLD,&status) ;

          if(set_ivls > 0) {
            vector<interval> ivl_list(set_ivls) ;
            MPI_Recv(&ivl_list[0],2*set_ivls,MPI_INT,partner,12,
                     MPI_COMM_WORLD, &status) ;
            for(int i=0;i<set_ivls;++i)
              keep += ivl_list[i] ;
          }
        } else { // Upper Partition
          int partner = r-split ;
          int set_ivls = 0 ;
          MPI_Status status ;
          MPI_Recv(&set_ivls,1,MPI_INT,partner,11,MPI_COMM_WORLD,&status) ;

          if(set_ivls > 0) {
            vector<interval> ivl_list(set_ivls) ;
            MPI_Recv(&ivl_list[0],2*set_ivls,MPI_INT,partner,12,MPI_COMM_WORLD,
                     &status) ;
            for(int i=0;i<set_ivls;++i)
              keep += ivl_list[i] ;
          }
          set_ivls = give.num_intervals() ;
          MPI_Send(&set_ivls,1,MPI_INT,partner,11,MPI_COMM_WORLD);

          if(set_ivls > 0) {
            vector<interval> ivl_list(set_ivls) ;
            for(int i=0;i<set_ivls;++i)
              ivl_list[i] = give[i] ;

            MPI_Send(&ivl_list[0],2*set_ivls,MPI_INT,partner,12,
                     MPI_COMM_WORLD) ;
          }
        }
      }
      e = keep ;
      nump = nnump ;
      start = nstart ;
    }

    return e & ptn[MPI_rank] ;
  }


  // Collect largest interval of entitySet from all processors
  entitySet collectLargest(const entitySet &e) {
    const int p = MPI_processes ;
    // Else we compute set
    //First get largest interval
    interval ivl_large(1,-1) ;

    const int esz = e.num_intervals() ;
    for(int i=0;i<esz;++i)
      if((ivl_large.second-ivl_large.first) < (e[i].second-e[i].first))
        ivl_large = e[i] ;

    vector<interval> ivl_large_p(p) ;

    MPI_Allgather(&ivl_large,2,MPI_INT,&(ivl_large_p[0]),2,MPI_INT,
                  MPI_COMM_WORLD) ;

    entitySet lset ;
    for(int i=0;i<p;++i)
      if(ivl_large_p[i].first <= ivl_large_p[i].second)
        lset += ivl_large_p[i] ;

    return lset ;
  }


  // Conceptually this will be equivalent to
  // entitySet all = all_collect_entitySet(set) ;
  // return domain & all ;
  
  entitySet collectSet(const entitySet iset, const entitySet domain,
		       MPI_Comm comm) {
    int np ;
    MPI_Comm_size(comm, &np) ;
    if(np == 1)
      return iset&domain ;
    int lmn = iset.Min() ;
    int lmx = iset.Max() ;
    int gmn = lmn ;
    int gmx = lmx ;
    MPI_Allreduce(&lmn,&gmn,1,MPI_INT,MPI_MIN,comm) ;
    MPI_Allreduce(&lmx,&gmx,1,MPI_INT,MPI_MAX,comm) ;
    if(gmn > gmx) // This means that all of the sets are empty
      return EMPTY ;
    
    // now just get a simple partition (note this is assuming the set is dense)
    int delta = (gmx-gmn+np)/np ;
    vector<entitySet> setsplit(np),domsplit(np) ;

    for(int i=0;i<np;++i) {
      setsplit[i] = iset & interval(gmn+i*delta,gmn+(i+1)*delta) ;
      domsplit[i] = domain & interval(gmn+i*delta,gmn+(i+1)*delta) ;
    }
    vector<entitySet> setgather =  transpose_entitySet(setsplit,comm)  ;
    entitySet dset=EMPTY ; // we are now redistributing set to processors ;
    for(int i=0;i<np;++i)
      dset += setgather[i] ;
    vector<entitySet> domgather = transpose_entitySet(domsplit,comm) ;
    // Now filter domain to set
    for(int i=0;i<np;++i)
      domgather[i] &= dset ;

    vector<entitySet> domreturn = transpose_entitySet(domgather,comm) ;
    entitySet rset=EMPTY ; // we are now redistributing set to processors ;
    for(int i=0;i<np;++i)
      rset += domreturn[i] ;
    return rset ;
  }
  
  
  entitySet all_collect_entitySet(const entitySet &e) {
    const int p = MPI_processes ;
    // no operation for single processor
    if(p == 1)
      return e ;
    

    // Check to see if the result should be EMPTY or UNIVERSE
    int code = 0 ;
    if(e != EMPTY)
      code = 1 ;
    if(e == ~EMPTY)
      code = 2 ;
    code = GLOBAL_MAX(code,MPI_COMM_WORLD) ;
    
    if(code == 0) // All empty, so return empty
      return EMPTY ;
    if(code == 2) // at least one UNIVERSE, so return UNIVERSE
      return ~EMPTY ;

    // First collect largest intervals, most of the time this will get
    // us the final set.  Compute remainder to determine if we have
    // more work to do.  Repeat collecting the largest interval for a
    // fixed number of tries, then finally just gather the remainder
#ifdef VERBOSE
    stopWatch s ;
    s.start() ;
#endif
    entitySet lset = collectLargest(e) ;
    entitySet rem = e-lset ;
    for(int i=0;i<4;++i) {
      int remsz = rem.num_intervals() ;
      int szmx = GLOBAL_MAX(remsz,MPI_COMM_WORLD) ;
      if(szmx == 0) {
#ifdef VERBOSE
        debugout << "time to get lset = " << s.stop() << endl ;
#endif
        return lset ;
      }
      lset += collectLargest(rem) ;
      rem -= lset ;
    }
#ifdef VERBOSE
    debugout << "time to get lset = " << s.stop() << endl ;
    s.start() ;
    
    debugout << "e="<< e.num_intervals() << ",rem=" << rem.num_intervals()
             << ",lset=" << lset.num_intervals() << endl ;
#endif
    entitySet remtot = all_gather_entitySet(rem) ;
#ifdef VERBOSE
    debugout << "time to gather rem = " << s.stop() << endl ;
#endif
    return lset + remtot ;
  }
  
  
  std::vector<entitySet> all_collect_vectors(entitySet &e,MPI_Comm comm) {
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    vector<entitySet> vset(p) ;
    if(p == 1) {
      vset[0] = e ;
      return vset ;
    }
    
    int send_count = 2*e.num_intervals() ;
    vector<int> recv_count(p) ;
    MPI_Allgather(&send_count,1,MPI_INT,&recv_count[0],1,MPI_INT,
                  comm) ;
    
    vector<int> recv_disp(p) ;
    recv_disp[0] = 0 ;
    for(int i=1;i<p;++i)
      recv_disp[i] = recv_disp[i-1]+recv_count[i-1] ;

    int tot = recv_disp[p-1]+recv_count[p-1] ;
    if(tot == 0)
      return vset ;
    vector<int> ivl_list(tot) ;
    vector<int> snd_list(send_count) ;
    for(int i=0;i<send_count/2;++i) {
      snd_list[i*2] = e[i].first ;
      snd_list[i*2+1] = e[i].second ;
    }
    MPI_Allgatherv(&(snd_list[0]),send_count,MPI_INT,
                   &(ivl_list[0]),&(recv_count[0]), &(recv_disp[0]), MPI_INT,
                   comm) ;

    for(int i = 0; i < p ; ++i) {
      int ind = recv_disp[i] ;
      for(int j=0;j<recv_count[i]/2;++j) {
        vset[i] += interval(ivl_list[ind+j*2],ivl_list[ind+j*2+1]) ;
      }
    }
    return vset ;
  }
  std::vector<entitySet> all_collect_vectors(entitySet &e) {
    return all_collect_vectors(e,MPI_COMM_WORLD) ;
  }

  int GLOBAL_OR(int b, MPI_Comm comm) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_LOR, comm) ;
    return result ;
  }
  
  int GLOBAL_AND(int b, MPI_Comm comm) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_LAND, comm) ;
    return result ;
  }
  
  int GLOBAL_MAX(int b, MPI_Comm comm) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_MAX, comm) ;
    return result ;
  }

  int GLOBAL_MIN(int b, MPI_Comm comm) {
    int result ;
    MPI_Allreduce(&b, &result, 1, MPI_INT, MPI_MIN, comm) ;
    return result ;
  }

  // remap of entitySet
  entitySet
  remap_entitySet(const entitySet& es, const Map& remap) {
    entitySet re ;
    for(entitySet::const_iterator ei=es.begin();ei!=es.end();++ei)
      re += remap[*ei] ;

    return re ;
  }
  
  entitySet
  remap_entitySet(const entitySet& es, const dMap& remap) {
    entitySet re ;
    for(entitySet::const_iterator ei=es.begin();ei!=es.end();++ei)
      re += remap[*ei] ;

    return re ;
  }

  // remap of a sequence
  sequence
  remap_sequence(const sequence& seq, const dMap& remap) {
    sequence re ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      re += remap[*si] ;

    return re ;
  }
  sequence
  remap_sequence(const sequence& seq, const Map& remap) {
    sequence re ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      re += remap[*si] ;

    return re ;
  }

  void
  get_p2p_comm(const std::vector<entitySet>& ptn,
               const entitySet& request,
               const dMap* send_remap,
               const dMap* recv_remap,
               MPI_Comm comm,
               std::vector<P2pCommInfo>& send,
               std::vector<P2pCommInfo>& recv) {
    send.clear() ;
    recv.clear() ;
    int np ;
    //MPI_Comm_rank(comm, &rank) ;
    MPI_Comm_size(comm, &np) ;
    FATAL(ptn.size() != size_t(np)) ;
    // obtain entity set partition for the request
    vector<entitySet> dist(ptn.size()) ;
    for(size_t i=0;i<ptn.size();++i)
      dist[i] = request & ptn[i] ;
    // than obtain the transpose of the dist
    vector<entitySet> dist_t = transpose_entitySet(dist,comm) ;

    for(size_t i=0;i<dist.size();++i) {
      const entitySet& r = dist[i] ;
      entitySet r_remap = r ;
      if(recv_remap != 0)
        r_remap = remap_entitySet(r,*recv_remap) ;
      if(r != EMPTY)
        recv.push_back(P2pCommInfo(i,r,r_remap)) ;

      const entitySet& s = dist_t[i] ;
      entitySet s_remap = s ;
      if(send_remap != 0)
        s_remap = remap_entitySet(s,*send_remap) ;
      if(s != EMPTY)
        send.push_back(P2pCommInfo(i,s,s_remap)) ;
    }
  }

#ifdef DYNAMICSCHEDULING
  void
  fill_store(storeRepP src, const Map* src_pack,
             storeRepP dst, const dMap* dst_unpack,
             const std::vector<P2pCommInfo>& send,
             const std::vector<P2pCommInfo>& recv,
             MPI_Comm comm) {
    int rank ;
    MPI_Comm_rank(comm, &rank) ;
    // we will need to communicate the size of the send buffer
    // to the receiving process so that they can allocate buffer.
    // we also need to communicate the pack entitySet in sequence
    // to the receiving processes so that they can properly unpack
    // the buffer.
    // normall, we need to send 3 messages: 1) the size of the total
    // packed buffer to send to a particular process, 2) the size
    // of the sequence to send, 3) the sequence itself. In order to
    // save message start-up time, we combine 1) and 2) messages
    // together into one message since they are both integer type.
    vector<int> pack_size(send.size(),0) ;
    vector<int> seq_size(send.size(),0) ;
    vector<bool> pack_interval(send.size(),false) ;
    vector<sequence> pack_seq(send.size()) ;
    // compute the pack size first
    for(size_t i=0;i<send.size();++i)
      pack_size[i] = src->pack_size(send[i].local_dom) ;
    // compute the packing sequence in global numbering
    // and also the way to send the sequence (in intervals
    // or direct elements), and the sequence send size.
    for(size_t i=0;i<send.size();++i) {
      // NOTE, we cannot use send[i].global_dom
      // as the pack sequence because if in deed
      // the src uses local numbering, then the
      // global numbering converted to a sequence
      // is not guaranteed to match the exact order
      // of the pack function, we therefore need
      // to map the sequence from the local numbering
      // directly
      // (error) pack_seq[i] = sequence(send[i].global_dom) ;
      sequence& ps = pack_seq[i] ;
      if(src_pack) {
        for(entitySet::const_iterator ei=send[i].local_dom.begin();
            ei!=send[i].local_dom.end();++ei)
          ps += (*src_pack)[*ei] ;
      } else {
        ps = sequence(send[i].local_dom) ;
      }
      int interval_size = pack_seq[i].num_intervals() ;
      int elem_size = pack_seq[i].size() ;
      if(2*interval_size < elem_size) {
        pack_interval[i] = true ;
        seq_size[i] = 2*interval_size + 1 ;
      } else {
        pack_interval[i] = false ;
        seq_size[i] = elem_size + 1 ;
      }
    }
    // now send the total pack size and seq size first
    vector<int> recv_msg_size(recv.size()*2, 0) ;

    vector<MPI_Request> requests(recv.size()) ;
    // first post recv requests to avoid deadlock
    int req = 0 ;
    // this index is used to optimize in the case of sending/receiving
    // messages to itself, we instead would just do a local copy
    int self_msg_buffer_idx = -1 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // remember the buffer index
        self_msg_buffer_idx = i ;
      } else {
        MPI_Irecv(&recv_msg_size[i*2], 2, MPI_INT,
                  recv[i].proc, recv[i].proc, comm, &requests[req++]) ;
      }
    }
    // then post send requests
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc == rank) {
        // just do a copy
        recv_msg_size[2*self_msg_buffer_idx] = pack_size[i] ;
        // this is another optimization, we do not need to
        // communicate the packing sequence for myself.
        // by setting the sequence size to be 0
        recv_msg_size[2*self_msg_buffer_idx+1] = 0 ;
      } else {
        int tmp[2] ;
        // first one is the total pack size (for the data)
        tmp[0] = pack_size[i] ;
        // second one is the pack sequence size
        tmp[1] = seq_size[i] ;
        MPI_Send(tmp, 2, MPI_INT, send[i].proc, rank, comm) ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then we actually need to communicate the packing
    // sequence to all the receiving processes

    // allocate recv buffer first
    int total_recv_size = 0 ;
    for(size_t i=0;i<recv.size();++i)
      total_recv_size += recv_msg_size[2*i+1] ;
    int* unpack_seq_recv_buffer = new int[total_recv_size] ;
    int** unpack_seq_recv_buffer_ptr = new int*[recv.size()] ;

    if(!recv.empty()) {
      unpack_seq_recv_buffer_ptr[0] = unpack_seq_recv_buffer ;
      for(size_t i=1;i<recv.size();++i)
        unpack_seq_recv_buffer_ptr[i] = recv_msg_size[2*(i-1)+1] +
          unpack_seq_recv_buffer_ptr[i-1] ;
    }

    // post recv requests (to receive the sequence)
    req = 0 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc != rank)
        MPI_Irecv(unpack_seq_recv_buffer_ptr[i],
                  recv_msg_size[2*i+1],
                  MPI_INT, recv[i].proc, recv[i].proc,
                  comm, &requests[req++]) ;
    }
    // send the sequence
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc != rank) {
        // allocate a send buffer
        vector<int> buf(seq_size[i]) ;
        // pack it
        if(pack_interval[i]) {
          // pack intervals
          buf[0] = 1 ;    // indicate following contents are intervals
          int count = 1 ;
          const sequence& seq = pack_seq[i] ;
          for(size_t k=0;k<seq.num_intervals();++k) {
            buf[count++] = seq[k].first ;
            buf[count++] = seq[k].second ;
          }
        } else {
          buf[0] = 0 ;          // we are packing elements
          int count = 1 ;
          const sequence& seq = pack_seq[i] ;
          for(sequence::const_iterator si=seq.begin();
              si!=seq.end();++si,++count)
            buf[count] = *si ;
        }
        // send the buffer
        MPI_Send(&buf[0], seq_size[i], MPI_INT, send[i].proc, rank, comm) ;
      } else {
        // remember the buffer index for later retrieval
        self_msg_buffer_idx = i ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then unpack the sequence buffer
    vector<sequence> unpack_seq(recv.size()) ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // just copy
        unpack_seq[i] = pack_seq[self_msg_buffer_idx] ;
      } else {
        // extract the first integer to see if the
        // packed stuff is interval or element
        sequence& seq = unpack_seq[i] ;
        int* p = unpack_seq_recv_buffer_ptr[i] ;
        int size = recv_msg_size[2*i+1] - 1 ;
        if(*p == 1) {
          // extract intervals
          ++p ;
          for(int k=0;k<size;k+=2) {
            int b = *p ; ++p ;
            int e = *p ; ++p ;
            seq += interval(b,e) ;
          }
        } else {
          // extract elements
          ++p ;
          for(int k=0;k<size;++k,++p)
            seq += *p ;
        }
      }
    }
    // release all unnecessary buffers
    vector<int>().swap(seq_size) ;
    vector<bool>().swap(pack_interval) ;
    vector<sequence>().swap(pack_seq) ;
    delete[] unpack_seq_recv_buffer_ptr ;
    delete[] unpack_seq_recv_buffer ;

    // remap the unpack sequence to the dst numbering
    if(dst_unpack != 0) {
      for(size_t i=0;i<recv.size();++i)
        unpack_seq[i] = remap_sequence(unpack_seq[i], *dst_unpack) ;
    }

    // now it is time for us to send/recv the data

    // allocate a recv buffer first
    total_recv_size = 0 ;
    for(size_t i=0;i<recv.size();++i)
      total_recv_size += recv_msg_size[2*i] ;
    unsigned char* unpack_buffer = new unsigned char[total_recv_size] ;
    unsigned char** unpack_buffer_ptr = new unsigned char*[recv.size()] ;

    if(!recv.empty()) {
      unpack_buffer_ptr[0] = unpack_buffer ;
      for(size_t i=1;i<recv.size();++i)
        unpack_buffer_ptr[i] = recv_msg_size[2*(i-1)] +
          unpack_buffer_ptr[i-1] ;
    }
    // post recv requests
    req = 0 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // remember the idx
        self_msg_buffer_idx = i ;
      } else {
        MPI_Irecv(unpack_buffer_ptr[i], recv_msg_size[2*i],
                  MPI_PACKED, recv[i].proc, recv[i].proc,
                  comm, &requests[req++]) ;
      }
    }
    // actually pack and send
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc == rank) {
        // no need to communicate, directly pack
        // into the receiving buffer
        unsigned char* pack_buffer =
          unpack_buffer_ptr[self_msg_buffer_idx] ;
        int position = 0 ;
        if(src_pack != 0)
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom, *src_pack) ;
        else
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom) ;
      } else {
        // first allocate the pack buffer
        unsigned char* pack_buffer = new unsigned char[pack_size[i]] ;
        // the do the pack
        int position = 0 ;
        if(src_pack != 0)
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom, *src_pack) ;
        else
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom) ;
        // send
        MPI_Send(pack_buffer, pack_size[i],
                 MPI_PACKED, send[i].proc, rank, comm) ;
        delete[] pack_buffer ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then unpack the receiving buffer
    for(size_t i=0;i<recv.size();++i) {
      int position = 0 ;
      int unpack_buffer_size = recv_msg_size[2*i] ;
      if(dst_unpack != 0)
        dst->unpack(unpack_buffer_ptr[i], position,
                    unpack_buffer_size, unpack_seq[i], *dst_unpack) ;
      else
        dst->unpack(unpack_buffer_ptr[i], position,
                    unpack_buffer_size, unpack_seq[i]) ;
    }
    // release recv buffer and finishing up
    delete[] unpack_buffer_ptr ;
    delete[] unpack_buffer ;
    // end of function fill_store
  }

  void
  expand_store(storeRepP src,
               const Map* src_pack, const dMap* src_unpack,
               storeRepP dst,
               const dMap* dst_unpack,
               const entitySet& request,
               const std::vector<entitySet>& src_ptn, MPI_Comm comm) {
    int rank ;
    MPI_Comm_rank(comm, &rank) ;
    // first we'll need to generate a send/recv structure
    vector<P2pCommInfo> send, recv ;
    get_p2p_comm(src_ptn, request,
                 src_unpack, dst_unpack, comm, send, recv) ;

    fill_store(src, src_pack, dst, dst_unpack, send, recv, comm) ;
  }
  
  void
  reduce_store(storeRepP src, const Map* src_pack,
               storeRepP dst, const dMap* dst_unpack,
               const std::vector<P2pCommInfo>& send,
               const std::vector<P2pCommInfo>& recv,
               CPTR<joiner> join_op, MPI_Comm comm) {
    int rank ;
    MPI_Comm_rank(comm, &rank) ;
    // we will need to communicate the size of the send buffer
    // to the receiving process so that they can allocate buffer.
    // we also need to communicate the pack entitySet in sequence
    // to the receiving processes so that they can properly unpack
    // the buffer.
    // normall, we need to send 3 messages: 1) the size of the total
    // packed buffer to send to a particular process, 2) the size
    // of the sequence to send, 3) the sequence itself. In order to
    // save message start-up time, we combine 1) and 2) messages
    // together into one message since they are both integer type.
    vector<int> pack_size(send.size(),0) ;
    vector<int> seq_size(send.size(),0) ;
    vector<bool> pack_interval(send.size(),false) ;
    vector<sequence> pack_seq(send.size()) ;
    // compute the pack size first
    for(size_t i=0;i<send.size();++i)
      pack_size[i] = src->pack_size(send[i].local_dom) ;
    // compute the packing sequence in global numbering
    // and also the way to send the sequence (in intervals
    // or direct elements), and the sequence send size.
    for(size_t i=0;i<send.size();++i) {
      // NOTE, we cannot use send[i].global_dom
      // as the pack sequence because if in deed
      // the src uses local numbering, then the
      // global numbering converted to a sequence
      // is not guaranteed to match the exact order
      // of the pack function, we therefore need
      // to map the sequence from the local numbering
      // directly
      // (error) pack_seq[i] = sequence(send[i].global_dom) ;
      sequence& ps = pack_seq[i] ;
      if(src_pack) {
        for(entitySet::const_iterator ei=send[i].local_dom.begin();
            ei!=send[i].local_dom.end();++ei)
          ps += (*src_pack)[*ei] ;
      } else {
        ps = sequence(send[i].local_dom) ;
      }
      int interval_size = pack_seq[i].num_intervals() ;
      int elem_size = pack_seq[i].size() ;
      if(2*interval_size < elem_size) {
        pack_interval[i] = true ;
        seq_size[i] = 2*interval_size + 1 ;
      } else {
        pack_interval[i] = false ;
        seq_size[i] = elem_size + 1 ;
      }
    }
    // now send the total pack size and seq size first
    vector<int> recv_msg_size(recv.size()*2, 0) ;

    vector<MPI_Request> requests(recv.size()) ;
    // first post recv requests to avoid deadlock
    int req = 0 ;
    // this index is used to optimize in the case of sending/receiving
    // messages to itself, we instead would just do a local copy
    int self_msg_buffer_idx = -1 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // remember the buffer index
        self_msg_buffer_idx = i ;
      } else {
        MPI_Irecv(&recv_msg_size[i*2], 2, MPI_INT,
                  recv[i].proc, recv[i].proc, comm, &requests[req++]) ;
      }
    }
    // then post send requests
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc == rank) {
        // just do a copy
        recv_msg_size[2*self_msg_buffer_idx] = pack_size[i] ;
        // this is another optimization, we do not need to
        // communicate the packing sequence for myself.
        // by setting the sequence size to be 0
        recv_msg_size[2*self_msg_buffer_idx+1] = 0 ;
      } else {
        int tmp[2] ;
        // first one is the total pack size (for the data)
        tmp[0] = pack_size[i] ;
        // second one is the pack sequence size
        tmp[1] = seq_size[i] ;
        MPI_Send(tmp, 2, MPI_INT, send[i].proc, rank, comm) ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then we actually need to communicate the packing
    // sequence to all the receiving processes

    // allocate recv buffer first
    int total_recv_size = 0 ;
    for(size_t i=0;i<recv.size();++i)
      total_recv_size += recv_msg_size[2*i+1] ;
    int* unpack_seq_recv_buffer = new int[total_recv_size] ;
    int** unpack_seq_recv_buffer_ptr = new int*[recv.size()] ;

    if(!recv.empty()) {
      unpack_seq_recv_buffer_ptr[0] = unpack_seq_recv_buffer ;
      for(size_t i=1;i<recv.size();++i)
        unpack_seq_recv_buffer_ptr[i] = recv_msg_size[2*(i-1)+1] +
          unpack_seq_recv_buffer_ptr[i-1] ;
    }

    // post recv requests (to receive the sequence)
    req = 0 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc != rank)
        MPI_Irecv(unpack_seq_recv_buffer_ptr[i],
                  recv_msg_size[2*i+1],
                  MPI_INT, recv[i].proc, recv[i].proc,
                  comm, &requests[req++]) ;
    }
    // send the sequence
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc != rank) {
        // allocate a send buffer
        vector<int> buf(seq_size[i]) ;
        // pack it
        if(pack_interval[i]) {
          // pack intervals
          buf[0] = 1 ;    // indicate following contents are intervals
          int count = 1 ;
          const sequence& seq = pack_seq[i] ;
          for(size_t k=0;k<seq.num_intervals();++k) {
            buf[count++] = seq[k].first ;
            buf[count++] = seq[k].second ;
          }
        } else {
          buf[0] = 0 ;          // we are packing elements
          int count = 1 ;
          const sequence& seq = pack_seq[i] ;
          for(sequence::const_iterator si=seq.begin();
              si!=seq.end();++si,++count)
            buf[count] = *si ;
        }
        // send the buffer
        MPI_Send(&buf[0], seq_size[i], MPI_INT, send[i].proc, rank, comm) ;
      } else {
        // remember the buffer index for later retrieval
        self_msg_buffer_idx = i ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then unpack the sequence buffer
    vector<sequence> unpack_seq(recv.size()) ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // just copy
        unpack_seq[i] = pack_seq[self_msg_buffer_idx] ;
      } else {
        // extract the first integer to see if the
        // packed stuff is interval or element
        sequence& seq = unpack_seq[i] ;
        int* p = unpack_seq_recv_buffer_ptr[i] ;
        int size = recv_msg_size[2*i+1] - 1 ;
        if(*p == 1) {
          // extract intervals
          ++p ;
          for(int k=0;k<size;k+=2) {
            int b = *p ; ++p ;
            int e = *p ; ++p ;
            seq += interval(b,e) ;
          }
        } else {
          // extract elements
          ++p ;
          for(int k=0;k<size;++k,++p)
            seq += *p ;
        }
      }
    }
    // release all unnecessary buffers
    vector<int>().swap(seq_size) ;
    vector<bool>().swap(pack_interval) ;
    vector<sequence>().swap(pack_seq) ;
    delete[] unpack_seq_recv_buffer_ptr ;
    delete[] unpack_seq_recv_buffer ;

    // remap the unpack sequence to the dst numbering
    if(dst_unpack != 0) {
      for(size_t i=0;i<recv.size();++i)
        unpack_seq[i] = remap_sequence(unpack_seq[i], *dst_unpack) ;
    }

    // now it is time for us to send/recv the data

    // allocate a recv buffer first
    total_recv_size = 0 ;
    for(size_t i=0;i<recv.size();++i)
      total_recv_size += recv_msg_size[2*i] ;
    unsigned char* unpack_buffer = new unsigned char[total_recv_size] ;
    unsigned char** unpack_buffer_ptr = new unsigned char*[recv.size()] ;

    if(!recv.empty()) {
      unpack_buffer_ptr[0] = unpack_buffer ;
      for(size_t i=1;i<recv.size();++i)
        unpack_buffer_ptr[i] = recv_msg_size[2*(i-1)] +
          unpack_buffer_ptr[i-1] ;
    }
    // post recv requests
    req = 0 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // remember the idx
        self_msg_buffer_idx = i ;
      } else {
        MPI_Irecv(unpack_buffer_ptr[i], recv_msg_size[2*i],
                  MPI_PACKED, recv[i].proc, recv[i].proc,
                  comm, &requests[req++]) ;
      }
    }
    // actually pack and send
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc == rank) {
        // no need to communicate, directly pack
        // into the receiving buffer
        unsigned char* pack_buffer =
          unpack_buffer_ptr[self_msg_buffer_idx] ;
        int position = 0 ;
        if(src_pack != 0)
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom, *src_pack) ;
        else
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom) ;
      } else {
        // first allocate the pack buffer
        unsigned char* pack_buffer = new unsigned char[pack_size[i]] ;
        // the do the pack
        int position = 0 ;
        if(src_pack != 0)
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom, *src_pack) ;
        else
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom) ;
        // send
        MPI_Send(pack_buffer, pack_size[i],
                 MPI_PACKED, send[i].proc, rank, comm) ;
        delete[] pack_buffer ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then unpack the receiving buffer
    entitySet recv_domain ;
    for(size_t i=0;i<recv.size();++i)
      recv_domain += recv[i].global_dom ;
    if(dst_unpack != 0)
      recv_domain = remap_entitySet(recv_domain, *dst_unpack) ;

    // create a tmp one for reduction purpose.
    storeRepP tmp = dst->new_store(recv_domain) ;
    
    for(size_t i=0;i<recv.size();++i) {
      int position = 0 ;
      int unpack_buffer_size = recv_msg_size[2*i] ;
      if(dst_unpack != 0) {
        // unpack to tmp first
        tmp->unpack(unpack_buffer_ptr[i], position,
                    unpack_buffer_size, unpack_seq[i], *dst_unpack) ;
      } else {
        tmp->unpack(unpack_buffer_ptr[i], position,
                    unpack_buffer_size, unpack_seq[i]) ;
      }
      // then join the contents to dst
      join_op->SetArgs(dst,   /*target rep*/
                       tmp    /*source rep*/) ;
      join_op->Join(unpack_seq[i]) ;
    }
    // release recv buffer and finishing up
    delete[] unpack_buffer_ptr ;
    delete[] unpack_buffer ;
    // end of function reduce_store
  }
  
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
               const std::vector<entitySet>& dst_ptn, MPI_Comm comm) {
    int rank, np ;
    MPI_Comm_rank(comm, &rank) ;
    MPI_Comm_size(comm, &np) ;
    // first we'll need to generate a send/recv structure
    vector<P2pCommInfo> send, recv ;
    // since this is a reduce, so we actually need to
    // reverse the send/recv pair as in the expand operation.
    get_p2p_comm(dst_ptn, request,
                 dst_unpack, src_unpack, comm, recv, send) ;

    reduce_store(src, src_pack, dst, dst_unpack,
                 send, recv, join_op, comm) ;
    // end of function reduce_store
  }

  // this version uses the pack_size(e,packed) inside
  entitySet
  fill_store2(storeRepP src, const Map* src_pack,
              storeRepP dst, const dMap* dst_unpack,
              const std::vector<P2pCommInfo>& send,
              const std::vector<P2pCommInfo>& recv,
              MPI_Comm comm) {
    int rank ;
    MPI_Comm_rank(comm, &rank) ;
    // we will need to communicate the size of the send buffer
    // to the receiving process so that they can allocate buffer.
    // we also need to communicate the pack entitySet in sequence
    // to the receiving processes so that they can properly unpack
    // the buffer.
    // normall, we need to send 3 messages: 1) the size of the total
    // packed buffer to send to a particular process, 2) the size
    // of the sequence to send, 3) the sequence itself. In order to
    // save message start-up time, we combine 1) and 2) messages
    // together into one message since they are both integer type.
    vector<int> pack_size(send.size(),0) ;
    vector<int> seq_size(send.size(),0) ;
    vector<bool> pack_interval(send.size(),false) ;
    vector<sequence> pack_seq(send.size()) ;
    // compute the pack size first
    for(size_t i=0;i<send.size();++i) {
      entitySet packed ;
      pack_size[i] = src->pack_size(send[i].local_dom,packed) ;

      // compute the packing sequence in global numbering
      // and also the way to send the sequence (in intervals
      // or direct elements), and the sequence send size.
      
      // NOTE, we cannot use send[i].global_dom
      // as the pack sequence because if in deed
      // the src uses local numbering, then the
      // global numbering converted to a sequence
      // is not guaranteed to match the exact order
      // of the pack function, we therefore need
      // to map the sequence from the local numbering
      // directly
      // (error) pack_seq[i] = sequence(send[i].global_dom) ;
      sequence& ps = pack_seq[i] ;
      if(src_pack) {
        for(entitySet::const_iterator ei=packed.begin();
            ei!=packed.end();++ei)
          ps += (*src_pack)[*ei] ;
      } else {
        ps = sequence(packed) ;
      }
      int interval_size = pack_seq[i].num_intervals() ;
      int elem_size = pack_seq[i].size() ;
      if(2*interval_size < elem_size) {
        pack_interval[i] = true ;
        seq_size[i] = 2*interval_size + 1 ;
      } else {
        pack_interval[i] = false ;
        seq_size[i] = elem_size + 1 ;
      }
    }
    // now send the total pack size and seq size first
    vector<int> recv_msg_size(recv.size()*2, 0) ;

    vector<MPI_Request> requests(recv.size()) ;
    // first post recv requests to avoid deadlock
    int req = 0 ;
    // this index is used to optimize in the case of sending/receiving
    // messages to itself, we instead would just do a local copy
    int self_msg_buffer_idx = -1 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // remember the buffer index
        self_msg_buffer_idx = i ;
      } else {
        MPI_Irecv(&recv_msg_size[i*2], 2, MPI_INT,
                  recv[i].proc, recv[i].proc, comm, &requests[req++]) ;
      }
    }
    // then post send requests
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc == rank) {
        // just do a copy
        recv_msg_size[2*self_msg_buffer_idx] = pack_size[i] ;
        // this is another optimization, we do not need to
        // communicate the packing sequence for myself.
        // by setting the sequence size to be 0
        recv_msg_size[2*self_msg_buffer_idx+1] = 0 ;
      } else {
        int tmp[2] ;
        // first one is the total pack size (for the data)
        tmp[0] = pack_size[i] ;
        // second one is the pack sequence size
        tmp[1] = seq_size[i] ;
        MPI_Send(tmp, 2, MPI_INT, send[i].proc, rank, comm) ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then we actually need to communicate the packing
    // sequence to all the receiving processes

    // allocate recv buffer first
    int total_recv_size = 0 ;
    for(size_t i=0;i<recv.size();++i)
      total_recv_size += recv_msg_size[2*i+1] ;
    int* unpack_seq_recv_buffer = new int[total_recv_size] ;
    int** unpack_seq_recv_buffer_ptr = new int*[recv.size()] ;

    if(!recv.empty()) {
      unpack_seq_recv_buffer_ptr[0] = unpack_seq_recv_buffer ;
      for(size_t i=1;i<recv.size();++i)
        unpack_seq_recv_buffer_ptr[i] = recv_msg_size[2*(i-1)+1] +
          unpack_seq_recv_buffer_ptr[i-1] ;
    }

    // post recv requests (to receive the sequence)
    req = 0 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc != rank)
        MPI_Irecv(unpack_seq_recv_buffer_ptr[i],
                  recv_msg_size[2*i+1],
                  MPI_INT, recv[i].proc, recv[i].proc,
                  comm, &requests[req++]) ;
    }
    // send the sequence
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc != rank) {
        // allocate a send buffer
        vector<int> buf(seq_size[i]) ;
        // pack it
        if(pack_interval[i]) {
          // pack intervals
          buf[0] = 1 ;    // indicate following contents are intervals
          int count = 1 ;
          const sequence& seq = pack_seq[i] ;
          for(size_t k=0;k<seq.num_intervals();++k) {
            buf[count++] = seq[k].first ;
            buf[count++] = seq[k].second ;
          }
        } else {
          buf[0] = 0 ;          // we are packing elements
          int count = 1 ;
          const sequence& seq = pack_seq[i] ;
          for(sequence::const_iterator si=seq.begin();
              si!=seq.end();++si,++count)
            buf[count] = *si ;
        }
        // send the buffer
        MPI_Send(&buf[0], seq_size[i], MPI_INT, send[i].proc, rank, comm) ;
      } else {
        // remember the buffer index for later retrieval
        self_msg_buffer_idx = i ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then unpack the sequence buffer
    vector<sequence> unpack_seq(recv.size()) ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // just copy
        unpack_seq[i] = pack_seq[self_msg_buffer_idx] ;
      } else {
        // extract the first integer to see if the
        // packed stuff is interval or element
        sequence& seq = unpack_seq[i] ;
        int* p = unpack_seq_recv_buffer_ptr[i] ;
        int size = recv_msg_size[2*i+1] - 1 ;
        if(*p == 1) {
          // extract intervals
          ++p ;
          for(int k=0;k<size;k+=2) {
            int b = *p ; ++p ;
            int e = *p ; ++p ;
            seq += interval(b,e) ;
          }
        } else {
          // extract elements
          ++p ;
          for(int k=0;k<size;++k,++p)
            seq += *p ;
        }
      }
    }
    // release all unnecessary buffers
    vector<int>().swap(seq_size) ;
    vector<bool>().swap(pack_interval) ;
    vector<sequence>().swap(pack_seq) ;
    delete[] unpack_seq_recv_buffer_ptr ;
    delete[] unpack_seq_recv_buffer ;

    // remap the unpack sequence to the dst numbering
    if(dst_unpack != 0) {
      for(size_t i=0;i<recv.size();++i)
        unpack_seq[i] = remap_sequence(unpack_seq[i], *dst_unpack) ;
    }

    // now it is time for us to send/recv the data

    // allocate a recv buffer first
    total_recv_size = 0 ;
    for(size_t i=0;i<recv.size();++i)
      total_recv_size += recv_msg_size[2*i] ;
    unsigned char* unpack_buffer = new unsigned char[total_recv_size] ;
    unsigned char** unpack_buffer_ptr = new unsigned char*[recv.size()] ;

    if(!recv.empty()) {
      unpack_buffer_ptr[0] = unpack_buffer ;
      for(size_t i=1;i<recv.size();++i)
        unpack_buffer_ptr[i] = recv_msg_size[2*(i-1)] +
          unpack_buffer_ptr[i-1] ;
    }
    // post recv requests
    req = 0 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // remember the idx
        self_msg_buffer_idx = i ;
      } else {
        MPI_Irecv(unpack_buffer_ptr[i], recv_msg_size[2*i],
                  MPI_PACKED, recv[i].proc, recv[i].proc,
                  comm, &requests[req++]) ;
      }
    }
    // actually pack and send
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc == rank) {
        // no need to communicate, directly pack
        // into the receiving buffer
        unsigned char* pack_buffer =
          unpack_buffer_ptr[self_msg_buffer_idx] ;
        int position = 0 ;
        if(src_pack != 0)
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom, *src_pack) ;
        else
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom) ;
      } else {
        // first allocate the pack buffer
        unsigned char* pack_buffer = new unsigned char[pack_size[i]] ;
        // the do the pack
        int position = 0 ;
        if(src_pack != 0)
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom, *src_pack) ;
        else
          src->pack(pack_buffer, position,
                    pack_size[i], send[i].local_dom) ;
        // send
        MPI_Send(pack_buffer, pack_size[i],
                 MPI_PACKED, send[i].proc, rank, comm) ;
        delete[] pack_buffer ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // first compute the union of all unpack sequence
    entitySet unpack_domain ;
    for(size_t i=0;i<unpack_seq.size();++i) {
      unpack_domain += entitySet(unpack_seq[i]) ;
    }
    dst->guarantee_domain(unpack_domain) ;

    // then unpack the receiving buffer
    for(size_t i=0;i<recv.size();++i) {
      int position = 0 ;
      int unpack_buffer_size = recv_msg_size[2*i] ;
      if(dst_unpack != 0)
        dst->unpack(unpack_buffer_ptr[i], position,
                    unpack_buffer_size, unpack_seq[i], *dst_unpack) ;
      else
        dst->unpack(unpack_buffer_ptr[i], position,
                    unpack_buffer_size, unpack_seq[i]) ;
    }
    // release recv buffer and finishing up
    delete[] unpack_buffer_ptr ;
    delete[] unpack_buffer ;

    return unpack_domain ;
    // end of function fill_store2
  }


  entitySet
  expand_store2(storeRepP src,
               const Map* src_pack, const dMap* src_unpack,
               storeRepP dst,
               const dMap* dst_unpack,
               const entitySet& request,
               const std::vector<entitySet>& src_ptn, MPI_Comm comm) {
    int rank ;
    MPI_Comm_rank(comm, &rank) ;
    // first we'll need to generate a send/recv structure
    vector<P2pCommInfo> send, recv ;
    get_p2p_comm(src_ptn, request,
                 src_unpack, dst_unpack, comm, send, recv) ;

    return
      fill_store2(src, src_pack, dst, dst_unpack, send, recv, comm) ;
  }

  // this version fills a vector of dstS from corresponding srcS
  vector<entitySet>
  fill_store2(std::vector<storeRepP>& src, const Map* src_pack,
              std::vector<storeRepP>& dst, const dMap* dst_unpack,
              const std::vector<P2pCommInfo>& send,
              const std::vector<P2pCommInfo>& recv, MPI_Comm comm) {
    FATAL(src.size() != dst.size()) ;
    int rank ;
    MPI_Comm_rank(comm, &rank) ;
    // we will need to communicate the size of the send buffer
    // to the receiving process so that they can allocate buffer.
    // we also need to communicate the pack entitySet in sequence
    // to the receiving processes so that they can properly unpack
    // the buffer.
    // normall, we need to send 3 messages: 1) the size of the total
    // packed buffer to send to a particular process, 2) the size
    // of the sequence to send, 3) the sequence itself. In order to
    // save message start-up time, we combine 1) and 2) messages
    // together into one message since they are both integer type.
    vector<int> pack_size(send.size(),0) ;
    vector<int> seq_size(send.size()*src.size(), 0) ;
    vector<bool> pack_interval(send.size()*src.size(), false) ;
    vector<entitySet> packed_dom(send.size()*src.size()) ;
    vector<sequence> pack_seq(send.size()*src.size()) ;
    vector<int> pack_seq_size(send.size(), 0) ;
    // compute the pack size first
    size_t ik = 0 ;
    for(size_t i=0;i<send.size();++i) {
      for(size_t k=0;k<src.size();++k,++ik) {
        entitySet& packed = packed_dom[ik] ;
        pack_size[i] += src[k]->pack_size(send[i].local_dom,packed) ;

        // compute the packing sequence in global numbering
        // and also the way to send the sequence (in intervals
        // or direct elements), and the sequence send size.
        
        // NOTE, we cannot use send[i].global_dom
        // as the pack sequence because if in deed
        // the src uses local numbering, then the
        // global numbering converted to a sequence
        // is not guaranteed to match the exact order
        // of the pack function, we therefore need
        // to map the sequence from the local numbering
        // directly
        // (error) pack_seq[i] = sequence(send[i].global_dom) ;
        sequence& ps = pack_seq[ik] ;
        if(src_pack) {
          for(entitySet::const_iterator ei=packed.begin();
              ei!=packed.end();++ei)
            ps += (*src_pack)[*ei] ;
        } else {
          ps = sequence(packed) ;
        }
        int interval_size = ps.num_intervals() ;
        int elem_size = ps.size() ;
        if(2*interval_size < elem_size) {
          pack_interval[ik] = true ;
          seq_size[ik] = 2*interval_size + 1 ;
          pack_seq_size[i] += seq_size[ik] ;
        } else {
          pack_interval[ik] = false ;
          seq_size[ik] = elem_size + 1 ;
          pack_seq_size[i] += seq_size[ik] ;
        }
      }
    }
    // now send the total pack size and seq size first
    int recv_msg_len = 1+dst.size() ;
    vector<int> recv_msg_size(recv.size()*recv_msg_len, 0) ;

    vector<MPI_Request> requests(recv.size()) ;
    // first post recv requests to avoid deadlock
    int req = 0 ;
    // this index is used to optimize in the case of sending/receiving
    // messages to itself, we instead would just do a local copy
    int self_msg_buffer_idx = -1 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // remember the buffer index
        self_msg_buffer_idx = i ;
      } else {
        MPI_Irecv(&recv_msg_size[i*recv_msg_len], recv_msg_len, MPI_INT,
                  recv[i].proc, recv[i].proc, comm, &requests[req++]) ;
      }
    }
    // then post send requests
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc == rank) {
        // just do a copy
        recv_msg_size[self_msg_buffer_idx*recv_msg_len] = pack_size[i] ;
        // this is another optimization, we do not need to
        // communicate the packing sequence for myself.
        // by setting the sequence size to be 0
        for(size_t k=0;k<dst.size();++k)
          recv_msg_size[self_msg_buffer_idx*recv_msg_len+k+1] = 0 ;
      } else {
        vector<int> tmp(recv_msg_len) ;
        // first one is the total pack size (for the data)
        tmp[0] = pack_size[i] ;
        // followings are the pack sequence size
        for(size_t k=0;k<dst.size();++k)
          tmp[k+1] = seq_size[i*src.size()+k] ;
        MPI_Send(&tmp[0], recv_msg_len, MPI_INT, send[i].proc, rank, comm) ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then we actually need to communicate the packing
    // sequence to all the receiving processes

    // allocate recv buffer first
    int total_recv_size = 0 ;
    vector<int> recv_seq_size(recv.size(), 0) ;
    for(size_t i=0;i<recv.size();++i) {
      for(size_t j=0;j<dst.size();++j) {
        recv_seq_size[i] += recv_msg_size[i*recv_msg_len+1+j] ;
      }
      total_recv_size += recv_seq_size[i] ;
    }
    int* unpack_seq_recv_buffer = new int[total_recv_size] ;
    int** unpack_seq_recv_buffer_ptr = new int*[recv.size()] ;

    if(!recv.empty()) {
      unpack_seq_recv_buffer_ptr[0] = unpack_seq_recv_buffer ;
      for(size_t i=1;i<recv.size();++i)
        unpack_seq_recv_buffer_ptr[i] = recv_seq_size[i-1] +
          unpack_seq_recv_buffer_ptr[i-1] ;
    }

    // post recv requests (to receive the sequence)
    req = 0 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc != rank)
        MPI_Irecv(unpack_seq_recv_buffer_ptr[i],
                  recv_seq_size[i],
                  MPI_INT, recv[i].proc, recv[i].proc,
                  comm, &requests[req++]) ;
    }
    // send the sequence
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc != rank) {
        // allocate a send buffer
        vector<int> buf(pack_seq_size[i]) ;
        // pack it
        size_t buf_idx = 0 ;
        for(size_t k=0;k<src.size();++k) {
          ik = i*src.size() + k ;
          if(pack_interval[ik]) {
            // pack intervals
            buf[buf_idx++] = 1 ;// indicate following contents are intervals
            const sequence& seq = pack_seq[ik] ;
            for(size_t j=0;j<seq.num_intervals();++j) {
              buf[buf_idx++] = seq[j].first ;
              buf[buf_idx++] = seq[j].second ;
            }
          } else {
            buf[buf_idx++] = 0 ;          // we are packing elements
            const sequence& seq = pack_seq[ik] ;
            for(sequence::const_iterator si=seq.begin();
                si!=seq.end();++si)
              buf[buf_idx++] = *si ;
          }
        }
        // send the buffer
        MPI_Send(&buf[0], pack_seq_size[i],
                 MPI_INT, send[i].proc, rank, comm) ;
      } else {
        // remember the buffer index for later retrieval
        self_msg_buffer_idx = i ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then unpack the sequence buffer
    vector<sequence> unpack_seq(recv.size() * dst.size()) ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // just copy
        for(size_t k=0;k<dst.size();++k) {
          unpack_seq[i*dst.size()+k]
            = pack_seq[self_msg_buffer_idx*src.size()+k] ;
        }
      } else {
        // extract the first integer to see if the
        // packed stuff is interval or element
        int* p = unpack_seq_recv_buffer_ptr[i] ;
        for(size_t k=0;k<dst.size();++k) {
          sequence& seq = unpack_seq[i*dst.size()+k] ;
          int size = recv_msg_size[i*recv_msg_len+1+k] - 1 ;
          if(*p == 1) {
            // extract intervals
            ++p ;
            for(int j=0;j<size;j+=2) {
              int b = *p ; ++p ;
              int e = *p ; ++p ;
              seq += interval(b,e) ;
            }
          } else {
            // extract elements
            ++p ;
            for(int j=0;j<size;++j,++p)
              seq += *p ;
          }
        }
      }
    }
    // release all unnecessary buffers
    vector<int>().swap(seq_size) ;
    vector<bool>().swap(pack_interval) ;
    vector<sequence>().swap(pack_seq) ;
    vector<int>().swap(pack_seq_size) ;
    delete[] unpack_seq_recv_buffer_ptr ;
    delete[] unpack_seq_recv_buffer ;

    // remap the unpack sequence to the dst numbering
    if(dst_unpack != 0) {
      for(size_t i=0;i<recv.size();++i)
        for(size_t k=0;k<dst.size();++k) {
          sequence& seq = unpack_seq[i*dst.size()+k] ;
          seq = remap_sequence(seq, *dst_unpack) ;
        }
    }

    // now it is time for us to send/recv the data

    // allocate a recv buffer first
    total_recv_size = 0 ;
    for(size_t i=0;i<recv.size();++i)
      total_recv_size += recv_msg_size[i*recv_msg_len] ;
    unsigned char* unpack_buffer = new unsigned char[total_recv_size] ;
    unsigned char** unpack_buffer_ptr = new unsigned char*[recv.size()] ;

    if(!recv.empty()) {
      unpack_buffer_ptr[0] = unpack_buffer ;
      for(size_t i=1;i<recv.size();++i)
        unpack_buffer_ptr[i] = recv_msg_size[(i-1)*recv_msg_len] +
          unpack_buffer_ptr[i-1] ;
    }
    // post recv requests
    req = 0 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // remember the idx
        self_msg_buffer_idx = i ;
      } else {
        MPI_Irecv(unpack_buffer_ptr[i], recv_msg_size[i*recv_msg_len],
                  MPI_PACKED, recv[i].proc, recv[i].proc,
                  comm, &requests[req++]) ;
      }
    }
    // actually pack and send
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc == rank) {
        // no need to communicate, directly pack
        // into the receiving buffer
        unsigned char* pack_buffer =
          unpack_buffer_ptr[self_msg_buffer_idx] ;
        int position = 0 ;
        if(src_pack != 0) {
          for(size_t k=0;k<src.size();++k)
            src[k]->pack(pack_buffer, position, pack_size[i],
                         packed_dom[i*src.size()+k], *src_pack) ;
        } else {
          for(size_t k=0;k<src.size();++k)
            src[k]->pack(pack_buffer, position,
                         pack_size[i], packed_dom[i*src.size()+k]) ;
        }
      } else {
        // first allocate the pack buffer
        unsigned char* pack_buffer = new unsigned char[pack_size[i]] ;
        // the do the pack
        int position = 0 ;
        if(src_pack != 0) {
          for(size_t k=0;k<src.size();++k)
            src[k]->pack(pack_buffer, position, pack_size[i],
                         packed_dom[i*src.size()+k], *src_pack) ;
        } else {
          for(size_t k=0;k<src.size();++k)
            src[k]->pack(pack_buffer, position,
                         pack_size[i], packed_dom[i*src.size()+k]) ;
        }
        // send
        MPI_Send(pack_buffer, pack_size[i],
                 MPI_PACKED, send[i].proc, rank, comm) ;
        delete[] pack_buffer ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // first compute the union of all unpack sequence
    vector<entitySet> unpack_domain(dst.size()) ;
    for(size_t i=0;i<dst.size();++i) {
      for(size_t k=0;k<recv.size();++k)
        unpack_domain[i] += entitySet(unpack_seq[k*dst.size()+i]) ;
      dst[i]->guarantee_domain(unpack_domain[i]) ;
    }

    // then unpack the receiving buffer
    for(size_t i=0;i<recv.size();++i) {
      int position = 0 ;
      int unpack_buffer_size = recv_msg_size[i*recv_msg_len] ;
      if(dst_unpack != 0) {
        for(size_t k=0;k<dst.size();++k)
          dst[k]->unpack(unpack_buffer_ptr[i], position,
                         unpack_buffer_size, unpack_seq[i*dst.size()+k],
                         *dst_unpack) ;
      } else {
        for(size_t k=0;k<dst.size();++k)
          dst[k]->unpack(unpack_buffer_ptr[i], position,
                         unpack_buffer_size, unpack_seq[i*dst.size()+k]) ;
      }
    }
    // release recv buffer and finishing up
    delete[] unpack_buffer_ptr ;
    delete[] unpack_buffer ;

    return unpack_domain ;
    // end of function fill_store2 (vector version)
  }
  
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
                const std::vector<entitySet>& src_ptn, MPI_Comm comm) {
    int rank ;
    MPI_Comm_rank(comm, &rank) ;
    // first we'll need to generate a send/recv structure
    vector<P2pCommInfo> send, recv ;
    get_p2p_comm(src_ptn, request,
                 src_unpack, dst_unpack, comm, send, recv) ;

    return
      fill_store2(src, src_pack, dst, dst_unpack, send, recv, comm) ;
  }
  
  entitySet
  fill_store_omd(storeRepP src, const Map* src_pack,
                 storeRepP dst, const dMap* dst_unpack,
                 const std::vector<P2pCommInfo>& send,
                 const std::vector<P2pCommInfo>& recv,
                 MPI_Comm comm) {
    int rank ;
    MPI_Comm_rank(comm, &rank) ;
    // we will need to communicate the size of the send buffer
    // to the receiving process so that they can allocate buffer.
    // we also need to communicate the pack entitySet in sequence
    // to the receiving processes so that they can properly unpack
    // the buffer.
    // normall, we need to send 3 messages: 1) the size of the total
    // packed buffer to send to a particular process, 2) the size
    // of the sequence to send, 3) the sequence itself. In order to
    // save message start-up time, we combine 1) and 2) messages
    // together into one message since they are both integer type.
    vector<int> pack_size(send.size(),0) ;
    vector<int> seq_size(send.size(),0) ;
    vector<bool> pack_interval(send.size(),false) ;
    vector<sequence> pack_seq(send.size()) ;
    // compute the pack size first
    for(size_t i=0;i<send.size();++i) {
      entitySet packed ;
      pack_size[i] = src->pack_size(send[i].local_dom,packed) ;

      // compute the packing sequence in global numbering
      // and also the way to send the sequence (in intervals
      // or direct elements), and the sequence send size.
      
      // NOTE, we cannot use send[i].global_dom
      // as the pack sequence because if in deed
      // the src uses local numbering, then the
      // global numbering converted to a sequence
      // is not guaranteed to match the exact order
      // of the pack function, we therefore need
      // to map the sequence from the local numbering
      // directly
      // (error) pack_seq[i] = sequence(send[i].global_dom) ;
      sequence& ps = pack_seq[i] ;
      if(src_pack) {
        for(entitySet::const_iterator ei=packed.begin();
            ei!=packed.end();++ei)
          ps += (*src_pack)[*ei] ;
      } else {
        ps = sequence(packed) ;
      }
      int interval_size = pack_seq[i].num_intervals() ;
      int elem_size = pack_seq[i].size() ;
      if(2*interval_size < elem_size) {
        pack_interval[i] = true ;
        seq_size[i] = 2*interval_size + 1 ;
      } else {
        pack_interval[i] = false ;
        seq_size[i] = elem_size + 1 ;
      }
    }
    // now send the total pack size and seq size first
    vector<int> recv_msg_size(recv.size()*2, 0) ;

    vector<MPI_Request> requests(recv.size()) ;
    // first post recv requests to avoid deadlock
    int req = 0 ;
    // this index is used to optimize in the case of sending/receiving
    // messages to itself, we instead would just do a local copy
    int self_msg_buffer_idx = -1 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // remember the buffer index
        self_msg_buffer_idx = i ;
      } else {
        MPI_Irecv(&recv_msg_size[i*2], 2, MPI_INT,
                  recv[i].proc, recv[i].proc, comm, &requests[req++]) ;
      }
    }
    // then post send requests
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc == rank) {
        // just do a copy
        recv_msg_size[2*self_msg_buffer_idx] = pack_size[i] ;
        // this is another optimization, we do not need to
        // communicate the packing sequence for myself.
        // by setting the sequence size to be 0
        recv_msg_size[2*self_msg_buffer_idx+1] = 0 ;
      } else {
        int tmp[2] ;
        // first one is the total pack size (for the data)
        tmp[0] = pack_size[i] ;
        // second one is the pack sequence size
        tmp[1] = seq_size[i] ;
        MPI_Send(tmp, 2, MPI_INT, send[i].proc, rank, comm) ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then we actually need to communicate the packing
    // sequence to all the receiving processes

    // allocate recv buffer first
    int total_recv_size = 0 ;
    for(size_t i=0;i<recv.size();++i)
      total_recv_size += recv_msg_size[2*i+1] ;
    int* unpack_seq_recv_buffer = new int[total_recv_size] ;
    int** unpack_seq_recv_buffer_ptr = new int*[recv.size()] ;

    if(!recv.empty()) {
      unpack_seq_recv_buffer_ptr[0] = unpack_seq_recv_buffer ;
      for(size_t i=1;i<recv.size();++i)
        unpack_seq_recv_buffer_ptr[i] = recv_msg_size[2*(i-1)+1] +
          unpack_seq_recv_buffer_ptr[i-1] ;
    }

    // post recv requests (to receive the sequence)
    req = 0 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc != rank)
        MPI_Irecv(unpack_seq_recv_buffer_ptr[i],
                  recv_msg_size[2*i+1],
                  MPI_INT, recv[i].proc, recv[i].proc,
                  comm, &requests[req++]) ;
    }
    // send the sequence
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc != rank) {
        // allocate a send buffer
        vector<int> buf(seq_size[i]) ;
        // pack it
        if(pack_interval[i]) {
          // pack intervals
          buf[0] = 1 ;    // indicate following contents are intervals
          int count = 1 ;
          const sequence& seq = pack_seq[i] ;
          for(size_t k=0;k<seq.num_intervals();++k) {
            buf[count++] = seq[k].first ;
            buf[count++] = seq[k].second ;
          }
        } else {
          buf[0] = 0 ;          // we are packing elements
          int count = 1 ;
          const sequence& seq = pack_seq[i] ;
          for(sequence::const_iterator si=seq.begin();
              si!=seq.end();++si,++count)
            buf[count] = *si ;
        }
        // send the buffer
        MPI_Send(&buf[0], seq_size[i], MPI_INT, send[i].proc, rank, comm) ;
      } else {
        // remember the buffer index for later retrieval
        self_msg_buffer_idx = i ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // then unpack the sequence buffer
    vector<sequence> unpack_seq(recv.size()) ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // just copy
        unpack_seq[i] = pack_seq[self_msg_buffer_idx] ;
      } else {
        // extract the first integer to see if the
        // packed stuff is interval or element
        sequence& seq = unpack_seq[i] ;
        int* p = unpack_seq_recv_buffer_ptr[i] ;
        int size = recv_msg_size[2*i+1] - 1 ;
        if(*p == 1) {
          // extract intervals
          ++p ;
          for(int k=0;k<size;k+=2) {
            int b = *p ; ++p ;
            int e = *p ; ++p ;
            seq += interval(b,e) ;
          }
        } else {
          // extract elements
          ++p ;
          for(int k=0;k<size;++k,++p)
            seq += *p ;
        }
      }
    }
    // release all unnecessary buffers
    vector<int>().swap(seq_size) ;
    vector<bool>().swap(pack_interval) ;
    vector<sequence>().swap(pack_seq) ;
    delete[] unpack_seq_recv_buffer_ptr ;
    delete[] unpack_seq_recv_buffer ;

    // remap the unpack sequence to the dst numbering
    if(dst_unpack != 0) {
      for(size_t i=0;i<recv.size();++i)
        unpack_seq[i] = remap_sequence(unpack_seq[i], *dst_unpack) ;
    }

    // now it is time for us to send/recv the data

    // allocate a recv buffer first
    total_recv_size = 0 ;
    for(size_t i=0;i<recv.size();++i)
      total_recv_size += recv_msg_size[2*i] ;
    unsigned char* unpack_buffer = new unsigned char[total_recv_size] ;
    unsigned char** unpack_buffer_ptr = new unsigned char*[recv.size()] ;

    if(!recv.empty()) {
      unpack_buffer_ptr[0] = unpack_buffer ;
      for(size_t i=1;i<recv.size();++i)
        unpack_buffer_ptr[i] = recv_msg_size[2*(i-1)] +
          unpack_buffer_ptr[i-1] ;
    }
    // post recv requests
    req = 0 ;
    for(size_t i=0;i<recv.size();++i) {
      if(recv[i].proc == rank) {
        // remember the idx
        self_msg_buffer_idx = i ;
      } else {
        MPI_Irecv(unpack_buffer_ptr[i], recv_msg_size[2*i],
                  MPI_PACKED, recv[i].proc, recv[i].proc,
                  comm, &requests[req++]) ;
      }
    }
    // actually pack and send
    for(size_t i=0;i<send.size();++i) {
      if(send[i].proc == rank) {
        // no need to communicate, directly pack
        // into the receiving buffer
        unsigned char* pack_buffer =
          unpack_buffer_ptr[self_msg_buffer_idx] ;
        int position = 0 ;
        src->pack(pack_buffer, position,
                  pack_size[i], send[i].local_dom) ;
      } else {
        // first allocate the pack buffer
        unsigned char* pack_buffer = new unsigned char[pack_size[i]] ;
        // the do the pack
        int position = 0 ;
        src->pack(pack_buffer, position,
                  pack_size[i], send[i].local_dom) ;
        // send
        MPI_Send(pack_buffer, pack_size[i],
                 MPI_PACKED, send[i].proc, rank, comm) ;
        delete[] pack_buffer ;
      }
    }
    // wait all Irecv to finish
    if(req > 0)
      MPI_Waitall(req, &requests[0], MPI_STATUSES_IGNORE) ;

    // first compute the union of all unpack sequence
    entitySet unpack_domain ;
    for(size_t i=0;i<unpack_seq.size();++i) {
      unpack_domain += entitySet(unpack_seq[i]) ;
    }
    dst->guarantee_domain(unpack_domain) ;

    // then unpack the receiving buffer
    for(size_t i=0;i<recv.size();++i) {
      int position = 0 ;
      int unpack_buffer_size = recv_msg_size[2*i] ;
      dst->unpack(unpack_buffer_ptr[i], position,
                  unpack_buffer_size, unpack_seq[i]) ;
    }
    // release recv buffer and finishing up
    delete[] unpack_buffer_ptr ;
    delete[] unpack_buffer ;

    return unpack_domain ;
    // end of function fill_store_omd
  }
#endif  
  // ... the end of file ...
} 
