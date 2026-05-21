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
#include <vector>
using std::vector;
#include <set>
using std::set;

#include <mpi.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include "dist_tools.h"
#include <rule.h>
#include <fact_db.h>
#include <constraint.h>
#include <multiMap.h>
#include <algorithm>
using std::sort ;

namespace Loci {

  // Pass in a set of entitySets one for each processor to broadcast
  // to other processors.  Return with a vector of entitySets as sent to
  // you by each other processor.
  vector<entitySet> Alltoall_entitySet(vector<entitySet> v) {
    int p = get_exec_size() ;
    WARN(int(v.size()) != p) ;

    if(p == 1)
      return v ;
    vector<int> ivals(p) ;
    for(int i=0;i<p;++i)
      ivals[i] = v[i].num_intervals() ;
    vector<int> ivalr(p) ;
    MPI_Comm comm = get_exec_comm() ;
    MPI_Alltoall(&(ivals[0]), 1, MPI_INT,
                 &(ivalr[0]), 1, MPI_INT, comm) ;
    vector<int> sdispls(p),rdispls(p) ;
    sdispls[0] = 0 ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      sdispls[i] = sdispls[i-1]+ivals[i-1]*2 ;
      rdispls[i] = rdispls[i-1]+ivalr[i-1]*2 ;
    }
    vector<int> sbuf(sdispls[p-1]+ivals[p-1]*2) ;
    vector<int> rbuf(rdispls[p-1]+ivalr[p-1]*2) ;
    vector<int> scounts(p),rcounts(p) ;
    for(int i=0;i<p;++i) {
      scounts[i] = 2*ivals[i] ;
      rcounts[i] = 2*ivalr[i] ;
      for(int j=0;j<ivals[i];++j) {
        sbuf[sdispls[i]+j*2] = v[i][j].first ;
        sbuf[sdispls[i]+j*2+1] = v[i][j].second ;
      }
    }
    MPI_Alltoallv(&(sbuf[0]), &(scounts[0]), &(sdispls[0]), MPI_INT,
                  &(rbuf[0]), &(rcounts[0]), &(rdispls[0]), MPI_INT,
                  comm) ;

    vector<entitySet> retv(p) ;
    for(int i=0;i<p;++i) {
      for(int j=0;j<ivalr[i];++j) {
        retv[i] += interval(rbuf[rdispls[i]+j*2],rbuf[rdispls[i]+j*2+1]) ;
      }
    }
    return retv ;
  }

  dMap distribute_dMap(dMap m, const std::vector<entitySet> &init_ptn) {
    int p = get_exec_size() ;
    if(p == 1)
      return m ;
    entitySet dom = m.domain() ;
    std::vector<entitySet> send_slices(p) ;
    for(int i=0;i<p;++i) {
      send_slices[i] = dom & init_ptn[i] ;
      dom -= init_ptn[i] ;
    }
    WARN(dom != EMPTY) ;

    vector<entitySet> recv_slices = Alltoall_entitySet(send_slices) ;

    vector<int> scounts(p),rcounts(p), sdispls(p),rdispls(p) ;
    for(int i=0;i<p;++i) {
      scounts[i] = send_slices[i].size() ;
      rcounts[i] = recv_slices[i].size() ;
    }
    sdispls[0] = 0 ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }

    vector<int> sbuf(sdispls[p-1]+scounts[p-1]) ;
    vector<int> rbuf(rdispls[p-1]+rcounts[p-1]) ;
    for(int i=0;i<p;++i) {
      int k =0 ;
      for(entitySet::const_iterator ei = send_slices[i].begin();
          ei != send_slices[i].end();++ei,++k) {
        sbuf[sdispls[i]+k] = m[*ei] ;
      }
    }

    MPI_Comm comm = get_exec_comm() ;
    MPI_Alltoallv(&(sbuf[0]), &(scounts[0]), &(sdispls[0]), MPI_INT,
                  &(rbuf[0]), &(rcounts[0]), &(rdispls[0]), MPI_INT,
                  comm) ;

    dMap ret_map ;

    for(int i=0;i<p;++i) {
      int k =0 ;
      for(entitySet::const_iterator ei = recv_slices[i].begin();
          ei != recv_slices[i].end();++ei,++k) {
        ret_map[*ei] = rbuf[rdispls[i]+k] ;
      }
    }

    return ret_map ;

  }

  namespace {
    inline bool fieldSort2(const std::pair<Entity,Entity> &p1,
                          const std::pair<Entity,Entity> &p2) {
      return p1.second < p2.second ;
    }
  }

  void distributed_inverseMap(multiMap &result,
                              vector<pair<Entity,Entity> > &input,
                              entitySet input_image,
                              entitySet input_preimage,
                              const std::vector<entitySet> &init_ptn) {
    MPI_Comm comm = get_exec_comm() ;
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    int r = 0 ; MPI_Comm_rank(comm, &r) ;
    // Sort input according to second field
    sort(input.begin(),input.end(),fieldSort2) ;

    // Now count what we will be sending
    vector<int> send_sz(p) ;
    for(int i=0;i<p;++i)
      send_sz[i] = 0 ;
    int current_p = 0 ;
    for(size_t i=0;i<input.size();++i) {
      int to = input[i].second ;

      if(!init_ptn[current_p].inSet(to))
        for(int j=0;j<p;++j)
          if(init_ptn[j].inSet(to)) {
            current_p = j ;
            break ;
          }
      send_sz[current_p] += 2 ;
    }

    // transfer recv sizes
    vector<int> recv_sz(p) ;
    MPI_Alltoall(&send_sz[0], 1, MPI_INT,
                 &recv_sz[0], 1, MPI_INT, comm) ;
    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<p;++i) {
      size_send += send_sz[i] ;
      size_recv += recv_sz[i] ;
    }

    int *send_store = new int[size_send] ;
    int *recv_store = new int[size_recv] ;
    int *send_displacement = new int[p] ;
    int *recv_displacement = new int[p] ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sz[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sz[i-1] ;
    }

    current_p = 0 ;
    vector<int> offsets(p) ;
    for(int i=0;i<p;++i)
      offsets[i] = 0 ;
    for(size_t i=0;i<input.size();++i) {
      int to = input[i].second ;
      int from = input[i].first ;
      if(!init_ptn[current_p].inSet(to))
        for(int j=0;j<p;++j)
          if(init_ptn[j].inSet(to)) {
            current_p = j ;
            break ;
          }
      send_store[send_displacement[current_p]+offsets[current_p]++] = to ;
      send_store[send_displacement[current_p]+offsets[current_p]++] = from ;
    }
    MPI_Alltoallv(send_store, &send_sz[0], send_displacement,
		  MPI_INT, recv_store, &recv_sz[0],
		  recv_displacement, MPI_INT, comm) ;

    entitySet local_input_image = input_image ;
    local_input_image &= init_ptn[r] ;
    store<int> sizes ;

    sizes.allocate(local_input_image) ;
    FORALL(local_input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    for(int i=0;i<size_recv;++i) {
      int indx = recv_store[i++] ;
      sizes[indx]++ ;
    }
    result.allocate(sizes) ;
    for(int i=0;i<size_recv;++i) {
      int indx = recv_store[i++] ;
      int refer = recv_store[i] ;
#ifdef DEBUG
      if(input_image.inSet(indx)) {
#endif
        sizes[indx]-- ;
        result[indx][sizes[indx]] = refer ;
#ifdef DEBUG
      } else {
	debugout << "indx=" << indx << " not in input_image=" << input_image
		 << endl ;
	Loci::Abort() ;
      }
#endif
    }
#ifdef DEBUG
    FORALL(local_input_image,i) {
      WARN(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
    
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;
  }

  void distributed_inverseMap(dmultiMap &result,
                              vector<pair<Entity,Entity> > &input,
                              entitySet input_image,
                              entitySet input_preimage,
                              const std::vector<entitySet> &init_ptn) {
    MPI_Comm comm = get_exec_comm() ;
    int p = 0 ; MPI_Comm_size(comm, &p) ;
    int r = 0 ; MPI_Comm_rank(comm, &r) ;
    // Sort input according to second field
    sort(input.begin(),input.end(),fieldSort2) ;

    // Now count what we will be sending
    vector<int> send_sz(p) ;
    for(int i=0;i<p;++i)
      send_sz[i] = 0 ;
    int current_p = 0 ;
    for(size_t i=0;i<input.size();++i) {
      int to = input[i].second ;

      if(!init_ptn[current_p].inSet(to))
        for(int j=0;j<p;++j)
          if(init_ptn[j].inSet(to)) {
            current_p = j ;
            break ;
          }
      send_sz[current_p] += 2 ;
    }

    // transfer recv sizes
    vector<int> recv_sz(p) ;
    MPI_Alltoall(&send_sz[0], 1, MPI_INT,
                 &recv_sz[0], 1, MPI_INT, comm) ;
    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<p;++i) {
      size_send += send_sz[i] ;
      size_recv += recv_sz[i] ;
    }

    int *send_store = new int[size_send] ;
    int *recv_store = new int[size_recv] ;
    int *send_displacement = new int[p] ;
    int *recv_displacement = new int[p] ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sz[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sz[i-1] ;
    }

    current_p = 0 ;
    vector<int> offsets(p) ;
    for(int i=0;i<p;++i)
      offsets[i] = 0 ;
    for(size_t i=0;i<input.size();++i) {
      int to = input[i].second ;
      int from = input[i].first ;
      if(!init_ptn[current_p].inSet(to))
        for(int j=0;j<p;++j)
          if(init_ptn[j].inSet(to)) {
            current_p = j ;
            break ;
          }
      send_store[send_displacement[current_p]+offsets[current_p]++] = to ;
      send_store[send_displacement[current_p]+offsets[current_p]++] = from ;
    }
    MPI_Alltoallv(send_store, &send_sz[0], send_displacement,
		  MPI_INT, recv_store, &recv_sz[0],
		  recv_displacement, MPI_INT, comm) ;

    entitySet local_input_image = input_image ;
    local_input_image &= init_ptn[r] ;
    store<int> sizes ;
    sizes.allocate(local_input_image) ;
    FORALL(local_input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    for(int i=0;i<size_recv;++i) {
      int indx = recv_store[i++] ;
      //      if(input_image.inSet(indx))
        sizes[indx]++ ;
    }
    
    FORALL(local_input_image,i) {
      vector<int,malloc_alloc<int> > tmp(sizes[i]) ;
      tmp.swap(result[i]) ;
    }ENDFORALL ;

    for(int i=0;i<size_recv;++i) {
      int indx = recv_store[i++] ;
      int refer = recv_store[i] ;
      sizes[indx]-- ;
      result[indx][sizes[indx]] = refer ;
    }
#ifdef DEBUG
    FORALL(local_input_image,i) {
      WARN(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
    
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;
  }
 
  void distributed_inverseMap(dmultiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {

    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
        inlist.push_back(pair<int,int>(i,elem)) ;
    } ENDFORALL ;
    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  void distributed_inverseMap(dmultiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {

    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
        inlist.push_back(pair<int,int>(i,elem)) ;
    } ENDFORALL ;
    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  void distributed_inverseMap(dmultiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      for(size_t k=0;k<input_map[i].size();++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) 
          inlist.push_back(pair<int,int>(i,elem)) ;
      }
    } ENDFORALL ;

    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  
  
  void distributed_inverseMap(dmultiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int sz = input_map.end(i)-input_map.begin(i) ;
      for(int k=0;k<sz;++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) 
          inlist.push_back(pair<int,int>(i,elem)) ;
      }
    } ENDFORALL ;

    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }
  // non dynamic multiMap
  void distributed_inverseMap(multiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {

    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
        inlist.push_back(pair<int,int>(i,elem)) ;
    } ENDFORALL ;
    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  void distributed_inverseMap(multiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {

    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) 
        inlist.push_back(pair<int,int>(i,elem)) ;
    } ENDFORALL ;
    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  void distributed_inverseMap(multiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      for(size_t k=0;k<input_map[i].size();++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) 
          inlist.push_back(pair<int,int>(i,elem)) ;
      }
    } ENDFORALL ;

    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }

  
  
  void distributed_inverseMap(multiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn) {
    vector<pair<int,int> > inlist ;
    entitySet preloop = input_preimage & input_map.domain() ;
    inlist.reserve(preloop.size()) ;
    FORALL(preloop,i) {
      int sz = input_map.end(i)-input_map.begin(i) ;
      for(int k=0;k<sz;++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) 
          inlist.push_back(pair<int,int>(i,elem)) ;
      }
    } ENDFORALL ;

    distributed_inverseMap(result,inlist,input_image,input_preimage,init_ptn) ;
  }
  


}
