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
    fatal(int(ptn.size()) != p) ;
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

    fatal(int(splits.size()) != p+1) ;
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
  void fill_clone( storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn, MPI_Comm comm) {
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    int r = 0 ; MPI_Comm_rank(comm, &r) ;

    vector<entitySet> recv_req(p) ;
    for(int i=0;i<p;++i)
      if(i!=r) 
        recv_req[i] = out_of_dom & init_ptn[i] ;

    // send the recieve requests
    vector<int> recv_count(p,0) ;
    vector<int> send_count(p,0) ;
    vector<int> send_displacement(p,0) ;
    vector<int> recv_displacement(p,0) ;
    for(int i=0;i<p;++i)
      send_count[i] = recv_req[i].num_intervals() * 2 ;
    MPI_Alltoall(&send_count[0],1,MPI_INT, &recv_count[0], 1, MPI_INT,comm) ;

    
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_displacement[i] = send_displacement[i-1]+send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1]+recv_count[i-1] ;
    }
    int mp = p-1 ;
    int send_sizes = send_displacement[mp]+send_count[mp] ;
    int recv_sizes = recv_displacement[mp]+recv_count[mp] ;


    vector<int> send_set_alloc(send_sizes+recv_sizes,0) ;
    int * send_set_buf = &send_set_alloc[0] ; 
    int * recv_set_buf = &send_set_alloc[send_sizes] ; 
    for(int i=0;i<p;++i) {
      for(size_t j=0;j<recv_req[i].num_intervals();++j) {
        send_set_buf[send_displacement[i]+j*2  ] = recv_req[i][j].first ;
        send_set_buf[send_displacement[i]+j*2+1] = recv_req[i][j].second ;
      }
    }

    MPI_Alltoallv(send_set_buf, &send_count[0], &send_displacement[0] , MPI_INT,
		  recv_set_buf, &recv_count[0], &recv_displacement[0], MPI_INT,
		  comm) ;

    vector<entitySet> send_set(p) ;
    for(int i=0;i<p;++i) {
      for(int j=0;j<recv_count[i]/2;++j) {
        int i1 = recv_set_buf[recv_displacement[i]+j*2  ] ;
        int i2 = recv_set_buf[recv_displacement[i]+j*2+1] ;
        send_set[i] += interval(i1,i2) ;
      }
    }
    
#ifdef DEBUG
    // Sanity check, no send set should be outside of entities we own
    entitySet mydom = init_ptn[r] ;
    for(int i=0;i<p;++i) {
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
    for(int i=0;i<p;++i) 
      send_count[i] =  sp->pack_size(send_set[i]) ;

    // Get sizes needed for receiving buffers
    MPI_Alltoall(&send_count[0],1,MPI_INT, &recv_count[0], 1, MPI_INT,comm) ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i=1;i<p;++i) {
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
    

    for(int i = 0; i <  p; ++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_displacement[i]], loc_pack, send_count[i],
               send_set[i]) ;
    }
    
    MPI_Alltoallv(send_store, &send_count[0], &send_displacement[0],
                  MPI_PACKED,
		  recv_store, &recv_count[0], &recv_displacement[0],
                  MPI_PACKED,
		  comm) ;

    for(int i = 0; i <  p; ++i) {
      int loc_pack = 0 ;
      sp->unpack(&recv_store[recv_displacement[i]], loc_pack,recv_count[i],
                 sequence(recv_req[i])) ; 
    }
  }
  
  storeRepP send_clone_non( storeRepP& sp, entitySet &out_of_dom, std::vector<entitySet> &init_ptn, MPI_Comm comm) {
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    int *recv_count = new int[p] ;
    int *send_count = new int[p] ;
    int *send_displacement = new int[p] ;
    int *recv_displacement = new int[p] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy(p), send_clone(p) ;
    for(int i = 0; i <  p; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i <  p; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  p; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    for(int i = 0; i <  p; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    std::vector< sequence> recv_dom(p) ;
    for(int i = 0; i <  p; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom(p) ;
    for(int i = 0; i <  p; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  p; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet tmp_dom ;
    for(int i = 0; i <  p; ++i)
      tmp_dom += entitySet(recv_dom[i]) ;
    storeRepP tmp_sp = sp->new_store(tmp_dom) ;
    for(int i = 0; i <  p; ++i) {
      recv_count[i] =  tmp_sp->pack_size(entitySet(recv_dom[i])) ;
      size_recv += recv_count[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  p; ++i) {
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    }
    
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  comm) ;  
    loc_pack = 0 ;
    for(int i = 0; i <  p; ++i) {
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

  std::vector<storeRepP> send_global_clone_non(storeRepP &sp , entitySet &out_of_dom,  std::vector<entitySet> &init_ptn, MPI_Comm comm) {
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    int *recv_count = new int[p] ;
    int *send_count = new int[p] ;
    int *send_displacement = new int[p] ;
    int *recv_displacement = new int[p] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy(p), send_clone(p) ;
    for(int i = 0; i <  p; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i <  p; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  p; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    for(int i = 0; i <  p; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    std::vector< sequence> recv_dom(p) ;
    for(int i = 0; i <  p; ++i) 
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) 
	recv_dom[i] += *vi ;
    std::vector<entitySet> send_dom(p) ;
    for(int i = 0; i <  p; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	send_dom[i] += *vi ;
    
    size_send = 0 ;
    for(int i = 0; i <  p; ++i) {
      send_count[i] =  sp->pack_size(send_dom[i]) ;
      size_send += send_count[i] ;
    } 
    unsigned char *send_store = new unsigned char[size_send] ;
    int size_recv = 0 ;
    entitySet total_dom ;
    std::vector<entitySet> e_vec(p) ;
    std::vector< storeRepP> tmp_sp(p) ;  
    for(int i = 0; i <  p; ++i) {
      e_vec[i] = entitySet(recv_dom[i]) ;
      tmp_sp[i] = sp->new_store(e_vec[i]) ;
      recv_count[i] =  tmp_sp[i]->pack_size(e_vec[i]) ;
      size_recv += recv_count[i] ;
      total_dom += e_vec[i] ;
    } 
    unsigned char *recv_store = new unsigned char[size_recv] ;
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    int loc_pack = 0 ;
    for(int i = 0; i <  p; ++i) {
      sp->pack(send_store, loc_pack, size_send, send_dom[i]) ;
    }
    MPI_Alltoallv(send_store,send_count, send_displacement , MPI_PACKED,
		  recv_store, recv_count, recv_displacement, MPI_PACKED,
		  comm) ;  
    
    loc_pack = 0 ;
    for(int i = 0; i <  p; ++i) 
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
  
  dMap send_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn, MPI_Comm comm) {
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    int *recv_count = new int[p] ;
    int *send_count = new int[p] ;
    int *send_displacement = new int[p] ;
    int *recv_displacement = new int[p] ;
    entitySet::const_iterator ei ;
    HASH_MAP(int, int) attrib_data ;
    entitySet dm_dom = dm.domain() ;
    for(ei = dm_dom.begin(); ei != dm_dom.end(); ++ei)
      attrib_data[*ei] = dm[*ei] ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy(p), send_clone(p) ;
    for(int i = 0; i <  p; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i <  p; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  p; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    for(int i = 0; i <  p; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    
    std::vector<HASH_MAP(int, int) > map_entities(p) ;
    for(int i = 0; i <  p; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(attrib_data.find(*vi) != attrib_data.end())
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i <  p; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i <  p; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  p; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    HASH_MAP(int, int) hm ;
    for(int i = 0; i <  p; ++i) {
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
  
  std::vector<dMap> send_global_map(Map &dm, entitySet &out_of_dom, std::vector<entitySet> &init_ptn, MPI_Comm comm) {
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    int *recv_count = new int[p] ;
    int *send_count = new int[p] ;
    int *send_displacement = new int[p] ;
    int *recv_displacement = new int[p] ;
    entitySet::const_iterator ei ;
    HASH_MAP(int, int) attrib_data ;
    entitySet dm_dom = dm.domain() ;
    for(ei = dm_dom.begin(); ei != dm_dom.end(); ++ei)
      attrib_data[*ei] = dm[*ei] ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy(p), send_clone(p) ;
    for(int i = 0; i <  p; ++i) {
      entitySet tmp = out_of_dom & init_ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	send_clone[i].push_back(*ei) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
      send_count[i] = send_clone[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i <  p; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  p; ++i)
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    for(int i = 0; i <  p; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	copy[i].push_back(recv_buf[j]) ;
      sort(copy[i].begin(), copy[i].end()) ;
    }
    
    std::vector<HASH_MAP(int, int) > map_entities(p) ;
    for(int i = 0; i <  p; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(attrib_data.find(*vi) != attrib_data.end())
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i <  p; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i <  p; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i <  p; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    std::vector<HASH_MAP(int, int) > hm(p) ;
    for(int i = 0; i <  p; ++i) {
      HASH_MAP(int, int) tmp_hm ;
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; j+=2) {
	tmp_hm[recv_map[j]] = recv_map[j+1];
      }
      hm[i] = tmp_hm ;
    }
    std::vector<dMap> v_dm(p) ;
    for(int i = 0; i <  p; ++i) {
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
  entitySet dist_collect_entitySet(entitySet inSet, const vector<entitySet> &ptn, MPI_Comm comm) {
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    int r = 0 ; MPI_Comm_rank(comm, &r) ;
#ifdef DEBUG
    if(r == 0) 
      cerr << "dist_collect_entitySet is depreciated" << endl ;
#endif
    entitySet retval = inSet & ptn[r] ;
    // Check for empty and universal set
    int sbits = ((inSet != EMPTY)?1:0)| ((retval != ptn[r])?2:0) ;
    int rbits = sbits ;
    MPI_Allreduce(&sbits, &rbits, 1, MPI_INT, MPI_BOR, comm) ;
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
    MPI_Allgather(&iset.first,1,MPI_2INT,&rset[0].first,1,MPI_2INT,comm) ;
    
    
    // Remove redundant info
    for(int i=0;i<p;++i)
      if(rset[i].second >= rset[i].first)
        inSet -= rset[i] ;

    
    s.start() ;
    // If the remainder is empty, then we are finished
    if(GLOBAL_AND((inSet==EMPTY),comm)) {
      return retval ;
    }

    // Remove second largest interval
    interval iset2(1,-1) ;
    for(size_t i=0;i<retval.num_intervals();++i)
      if(iset != retval[i] &&
         (iset2.second-iset2.first) < (retval[i].second-retval[i].first))
        iset2 = retval[i] ;
    
    MPI_Allgather(&iset2.first,1,MPI_2INT,&rset[0].first,1,MPI_2INT,comm) ;
    
    
    // Remove redundant info
    for(int i=0;i<p;++i)
      if(rset[i].second >= rset[i].first)
        inSet -= rset[i] ;

    
    // If the remainder is empty, then we are finished
    if(GLOBAL_AND((inSet==EMPTY),comm)) {
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
                 comm) ;

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
    		  comm) ;

    // Sort received intervals to increase performance of union code
    sort(recv_store.begin(),recv_store.end()) ;

    for(size_t i=0;i<recv_store.size();++i)
      retval += recv_store[i] ;
    
    return retval ;
  }
#ifndef LOCI_STRICT_COMM
  entitySet dist_collect_entitySet(entitySet inSet, const vector<entitySet> &ptn) {
    dist_collect_entitySet(inSet,ptn,MPI_COMM_WORLD) ;
  }
#endif
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
                                  const vector<entitySet> &ptn, MPI_Comm comm) {
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
                 comm) ;

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
                 comm) ;

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
                  comm,&recv_Requests[cnt++]) ;
      }
    }

    // Send the data to partner
    for(size_t i=0;i < ptn.size();++i) {
      if(send_sizes[i] != 0) {
        MPI_Send(&(send_data[i][0]),send_sizes[i]*2,MPI_INT,i,3,comm) ;
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

#ifndef LOCI_STRICT_COMM
  entitySet dist_expand_entitySet(entitySet inSet, entitySet clone,
                                  const vector<entitySet> &ptn) {
    dist_expand_entitySet(inSet,clone,ptn,MPI_COMM_WORLD) ;
  }
#endif
  inline bool spec_ival_compare(const interval &i1,
                            const interval &i2) {
    if(i1.first < i2.first)
      return true ;
    if(i1.first == i2.first && i1.second > i2.second)
      return true ;
    return false ;
  }
      
  // Return union of all entitySets from all processors
  entitySet all_gather_entitySet(const entitySet &e, MPI_Comm comm) {
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    if(p == 1)
      return e ;
    
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
      return EMPTY ;
    vector<interval> ivl_list(tot/2) ;
    vector<interval> snd_list(send_count/2) ;
    for(int i=0;i<send_count/2;++i)
      snd_list[i] = e[i] ;
    MPI_Allgatherv(&(snd_list[0]),send_count,MPI_INT,
                   &(ivl_list[0]),&(recv_count[0]), &(recv_disp[0]), MPI_INT,
                   comm) ;
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
  entitySet distribute_entitySet(entitySet e,const vector<entitySet> &ptn, MPI_Comm comm) {
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    // Single processor, do nothing
    if(p == 1)
      return e ;

    // if any processor contains universal set, all processors hold the
    // universal set
    if(GLOBAL_OR(e == ~EMPTY, comm)) 
      return ~EMPTY ;

    int r = 0 ; MPI_Comm_rank(comm, &r) ;
    // if no communication need, then just return the input set
    if(GLOBAL_AND((e & ptn[r]) == e, comm)) 
      return e ;


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
        MPI_Send(&set_ivls,1,MPI_INT,r-1,11,comm);
        
        if(set_ivls > 0) {
          vector<interval> ivl_list(set_ivls) ;
          for(int i=0;i<set_ivls;++i)
            ivl_list[i] = give[i] ;
          MPI_Send(&ivl_list[0],2*set_ivls,MPI_INT,r-1,12,comm) ;
        }
      }
      if(r+1==oddproc) {
        int set_ivls = 0 ;
        MPI_Status status ;
        MPI_Recv(&set_ivls,1,MPI_INT,r+1,11,comm,&status) ;

        if(set_ivls > 0) {
          vector<interval> ivl_list(set_ivls) ;
          MPI_Recv(&ivl_list[0],2*set_ivls,MPI_INT,r+1,12,comm,
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
          MPI_Send(&set_ivls,1,MPI_INT,partner,11,comm);

          if(set_ivls > 0) {
            vector<interval> ivl_list(set_ivls) ;
            for(int i=0;i<set_ivls;++i)
              ivl_list[i] = give[i] ;
            MPI_Send(&ivl_list[0],2*set_ivls,MPI_INT,partner,12,
                     comm) ;
          }
          MPI_Status status ;
          MPI_Recv(&set_ivls,1,MPI_INT,partner,11,comm,&status) ;

          if(set_ivls > 0) {
            vector<interval> ivl_list(set_ivls) ;
            MPI_Recv(&ivl_list[0],2*set_ivls,MPI_INT,partner,12,
                     comm, &status) ;
            for(int i=0;i<set_ivls;++i)
              keep += ivl_list[i] ;
          }
        } else { // Upper Partition
          int partner = r-split ;
          int set_ivls = 0 ;
          MPI_Status status ;
          MPI_Recv(&set_ivls,1,MPI_INT,partner,11,comm,&status) ;

          if(set_ivls > 0) {
            vector<interval> ivl_list(set_ivls) ;
            MPI_Recv(&ivl_list[0],2*set_ivls,MPI_INT,partner,12,comm,
                     &status) ;
            for(int i=0;i<set_ivls;++i)
              keep += ivl_list[i] ;
          }
          set_ivls = give.num_intervals() ;
          MPI_Send(&set_ivls,1,MPI_INT,partner,11,comm);

          if(set_ivls > 0) {
            vector<interval> ivl_list(set_ivls) ;
            for(int i=0;i<set_ivls;++i)
              ivl_list[i] = give[i] ;

            MPI_Send(&ivl_list[0],2*set_ivls,MPI_INT,partner,12,
                     comm) ;
          }
        }
      }
      e = keep ;
      nump = nnump ;
      start = nstart ;
    }

    return e & ptn[r] ;
  }


  // Collect largest interval of entitySet from all processors
  entitySet collectLargest(const entitySet &e, MPI_Comm comm) {
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    // Else we compute set
    //First get largest interval
    interval ivl_large(1,-1) ;

    const int esz = e.num_intervals() ;
    for(int i=0;i<esz;++i)
      if((ivl_large.second-ivl_large.first) < (e[i].second-e[i].first))
        ivl_large = e[i] ;

    vector<interval> ivl_large_p(p) ;

    MPI_Allgather(&ivl_large,2,MPI_INT,&(ivl_large_p[0]),2,MPI_INT,
                  comm) ;

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
  
  
  entitySet all_collect_entitySet(const entitySet &e, MPI_Comm comm) {
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    // no operation for single processor
    if(p == 1)
      return e ;
    

    // Check to see if the result should be EMPTY or UNIVERSE
    int code = 0 ;
    if(e != EMPTY)
      code = 1 ;
    if(e == ~EMPTY)
      code = 2 ;
    code = GLOBAL_MAX(code,comm) ;
    
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
    entitySet lset = collectLargest(e, comm) ;
    entitySet rem = e-lset ;
    for(int i=0;i<4;++i) {
      int remsz = rem.num_intervals() ;
      int szmx = GLOBAL_MAX(remsz,comm) ;
      if(szmx == 0) {
#ifdef VERBOSE
        debugout << "time to get lset = " << s.stop() << endl ;
#endif
        return lset ;
      }
      lset += collectLargest(rem, comm) ;
      rem -= lset ;
    }
#ifdef VERBOSE
    debugout << "time to get lset = " << s.stop() << endl ;
    s.start() ;
    
    debugout << "e="<< e.num_intervals() << ",rem=" << rem.num_intervals()
             << ",lset=" << lset.num_intervals() << endl ;
#endif
    entitySet remtot = all_gather_entitySet(rem, comm) ;
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
  //  std::vector<entitySet> all_collect_vectors(entitySet &e) {
  //    return all_collect_vectors(e,get_exec_comm()) ;
  //  }

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

  // ... the end of file ...
} 
