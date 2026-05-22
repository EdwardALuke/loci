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
#ifndef DISTRIBUTE_CONTAINER_H
#define DISTRIBUTE_CONTAINER_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <vector>
#include <algorithm>

#include <DMap.h>
#include <store_rep.h>
#include <DMultiMap.h>
#include <fact_db.h>
#include <Map.h>
#include <multiMap.h>
#include <store.h>
#include <distribute_io.h>

namespace Loci {
  
  void distributed_inverseMap(multiMap &result,
                              std::vector<std::pair<Entity,Entity> > &input,
                              entitySet input_image,
                              entitySet input_preimage,
                              const std::vector<entitySet> &init_ptn,
                              MPI_Comm comm LOCI_DEFAULT_COMM) ;

  dMap distribute_dMap(dMap m, const std::vector<entitySet> &init_ptn, MPI_Comm comm LOCI_DEFAULT_COMM) ;

  void distributed_inverseMap(dmultiMap &result, const dMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn,
                              MPI_Comm comm LOCI_DEFAULT_COMM) ;
  
  void distributed_inverseMap(dmultiMap &result, const Map &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn,
                              MPI_Comm comm LOCI_DEFAULT_COMM) ;
  
  void distributed_inverseMap(dmultiMap &result,
                              const dmultiMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn,
                              MPI_Comm comm LOCI_DEFAULT_COMM) ;
  
  void distributed_inverseMap(dmultiMap &result,
                              const multiMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn,
                              MPI_Comm comm LOCI_DEFAULT_COMM); 

  void distributed_inverseMap(dmultiMap &result,
                              const const_dMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result,
                              const const_Map &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result,
                              const const_dmultiMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn) ;
  
  void distributed_inverseMap(dmultiMap &result,
                              const const_multiMap &input_map,
                              const entitySet &input_image,
                              const entitySet &input_preimage,
                              std::vector<entitySet> &init_ptn); 

  inline void distributed_inverseMap(dmultiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts, size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }
  
  inline void distributed_inverseMap(dmultiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts, size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }      

  inline void distributed_inverseMap(dmultiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts,size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }

  inline void distributed_inverseMap(dmultiMap &result,
                                     const multiMap &input_map,
                                     const entitySet &input_image,
                                     const entitySet &input_preimage,
                                     fact_db &facts,
				     size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn) ;
    } else {
      inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }

  void distributed_inverseMap(multiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn, MPI_Comm comm LOCI_DEFAULT_COMM) ;
  
  void distributed_inverseMap(multiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn, MPI_Comm comm LOCI_DEFAULT_COMM) ;
  
  void distributed_inverseMap(multiMap &result, const dmultiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn, MPI_Comm comm LOCI_DEFAULT_COMM) ;
  
  void distributed_inverseMap(multiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, std::vector<entitySet> &init_ptn, MPI_Comm comm LOCI_DEFAULT_COMM);
  

  inline void distributed_inverseMap(multiMap &result, const Map &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts, size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn,facts.get_comm()) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }      

  inline void distributed_inverseMap(multiMap &result, const multiMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts,size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn,facts.get_comm()) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }      

  inline void distributed_inverseMap(multiMap &result, const dMap &input_map, const entitySet &input_image, const entitySet &input_preimage, fact_db &facts, size_t kd) {
    if(facts.is_distributed_start()) {
      std::vector<entitySet> init_ptn = facts.get_init_ptn(kd) ;
      Loci::distributed_inverseMap(result, input_map, input_image, input_preimage, init_ptn,facts.get_comm()) ;
    } else {
      Loci::inverseMap(result,input_map,input_image,input_preimage) ;
    }
  }


  typedef std::vector<std::pair<int,int> > protoMap;
  
  inline bool equiFF(const std::pair<int,int> &v1,
                     const std::pair<int,int> &v2) {
    return v1.first < v2.first ;
  }

  inline void equiJoinFF(protoMap &in1, protoMap &in2, protoMap &out,
                         MPI_Comm comm LOCI_DEFAULT_COMM) {
    std::sort(in1.begin(),in1.end(),equiFF) ;
    std::sort(in2.begin(),in2.end(),equiFF) ;
    
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    // Sort inputs using same splitters (this will make sure that
    // data that needs to be on the same processor ends up on the
    // same processor
    if(p != 1) {
      std::vector<std::pair<int,int> > splitters ;
      parGetSplitters(splitters,in1,equiFF,comm) ;

      parSplitSort(in1,splitters,equiFF,comm) ;
      parSplitSort(in2,splitters,equiFF,comm) ;
    }

    // Find pairs where first entry are the same and create joined protomap
    out.clear() ;
    size_t j = 0 ;
    for(size_t i=0;i<in1.size();++i) {
      while(j<in2.size() && equiFF(in2[j],in1[i]))
        ++j ;
      size_t k=j ;
      while(k<in2.size() && in2[k].first == in1[i].first) {
        out.push_back(std::pair<int,int>(in1[i].second,in2[k].second)) ;
        k++ ;
      }
    }

    // Remove duplicates from protomap
    parSampleSort(out,equiFF,comm) ;
    std::sort(out.begin(),out.end()) ;
    out.erase(std::unique(out.begin(),out.end()),out.end()) ;
  }

  inline void removeIdentity(protoMap &mp) {
      // Remove self references
    protoMap::iterator ii,ij ;
    for(ii=ij=mp.begin();ij!=mp.end();++ij) {
      if(ij->first != ij->second) {
        *ii = *ij ;
        ii++ ;
      }
    }
    mp.erase(ii,mp.end()) ;
  }
  inline void addToProtoMap(multiMap &m, protoMap &mp) {
    entitySet dom = m.domain() ;
    FORALL(dom,e) {
      int sz = m[e].size() ;
      for(int i=0;i<sz;++i)
        mp.push_back(std::pair<int,int>(e,m[e][i])) ;
    } ENDFORALL ;
  }

  inline void addToProtoMap(Map &m, protoMap &mp) {
    entitySet dom = m.domain() ;
    FORALL(dom,e) {
      mp.push_back(std::pair<int,int>(e,m[e])) ;
    } ENDFORALL ;
  }

    // -----------------------------------------------------------------------
  /// @brief getLocalContextMap returns a global to local mapping that
  /// can be used to translate global numbers to the local numbering
  /// provided by a data access gathering set.  The local numbering is just
  /// the ordering of the data in the get set.
  ///
  /// @param [g2l] The resulting global to local map
  /// @param [get] The sets that each processor will be accessing
  template<class T> void getLocalContextMap(T &g2l, entitySet get) {
    int cnt = 0 ;
    for(size_t i=0;i<get.num_intervals();++i)
      for(Entity ii=get[i].first;ii<=get[i].second;++ii)
        g2l[ii] = cnt++ ;
  }

  // -----------------------------------------------------------------------
  /// @brief gatherData gathers the distributed data that will be accessed
  /// by each processor from a store that is partitioned according to ptn,
  ///
  /// @param[result] std::vector that will contain the gathered data for
  /// each processor, this will be ordered as the entitySet provided
  /// @param[data] The input store container that is distributed across
  /// processors
  /// @param[get] The requested set of entities that will be accessed by each
  /// processor
  /// @param[ptn] The partition of the data stored in container [data]
  template <class T>  void gatherData(std::vector<T> &result,
                                      store<T> &data,
                                      entitySet get,
                                      std::vector<entitySet> &ptn,
                                      MPI_Comm comm LOCI_DEFAULT_COMM) {
    using std::vector ;
    { vector<T> tmp(get.size()) ; result.swap(tmp) ;}

    const int r = Loci::MPI_rank ;
    const int p = Loci::MPI_processes ;
    FATAL(p!=int(ptn.size())) ;

    vector<entitySet> recv_sets(p) ;
    for(int i=0;i<p;++i)
      recv_sets[i] = get&ptn[i] ;

    vector<entitySet> send_sets =
      transpose_entitySet(recv_sets, comm) ;
    int send_sz = 0 ;
    int reqs = 0 ;
    for(int i=0;i<p;++i) {
      send_sz += send_sets[i].size() ;
      if(i!=r) {
        if(send_sets[i] != EMPTY)
          reqs++ ;
        if(recv_sets[i] != EMPTY)
          reqs++ ;
      }
    }
    std::vector<T> send_data(send_sz) ;
    int j = 0 ;
    for(int i=0;i<p;++i) {
      FORALL(send_sets[i],ii) {
        send_data[j] = data[ii] ;
        j++ ;
      } ENDFORALL ;
    }
    vector<int> send_offsets(p,0) ;
    vector<int> recv_offsets(p,0) ;
    for(int i=1;i<p;++i) {
      send_offsets[i] = send_offsets[i-1]+send_sets[i-1].size() ;
      recv_offsets[i] = recv_offsets[i-1]+recv_sets[i-1].size() ;
    }

    // copy self data first
    for(size_t i=0;i<send_sets[r].size();++i)
      result[recv_offsets[r]+i] = send_data[send_offsets[r]+i] ;

    // Single processor we are done
    if(p == 1)
      return ;
    
    std::vector<MPI_Request> requests(reqs) ;
    int cnt = 0 ;
    for(int i=0;i<p;++i)
      // If I am not recving from myself and I have something to recv
      if(i!=r && recv_sets[i] != EMPTY) {
        MPI_Irecv(&result[recv_offsets[i]],
                  recv_sets[i].size()*sizeof(T),
                  MPI_BYTE,i,99,comm,&requests[cnt]) ;
        cnt++ ;
      }
    for(int i=0;i<p;++i)
      if(i!=r && send_sets[i] != EMPTY) {
        // If I am not sending to myself and I have something to send
        MPI_Isend(&send_data[send_offsets[i]],
                  send_sets[i].size()*sizeof(T),
                  MPI_BYTE,i,99,comm,&requests[cnt]) ;
        cnt++ ;
      }

    // If this processor does not commuicate we are done
    if(cnt == 0)
      return ;

    // otherwise wait for communication to finish
    vector<MPI_Status> status_list(cnt) ;
    MPI_Waitall(cnt,&requests[0],&status_list[0]) ;
  }

}

#endif
