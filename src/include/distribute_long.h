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
#ifndef DISTRIBUTE_LONG_H
#define DISTRIBUTE_LONG_H

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



namespace Loci {

  
  extern std::ofstream debugout ;
  extern int MPI_processes;
  extern int MPI_rank ;
  
  template<class T> inline bool g_spec_ival_compare(const std::pair<T, T> &i1,
                                                  const std::pair<T, T> &i2) {
    if(i1.first < i2.first)
      return true ;
    if(i1.first == i2.first && i1.second > i2.second)
      return true ;
    return false ;
  }

  template<class T> T g_GLOBAL_MAX(T b, MPI_Comm comm=MPI_COMM_WORLD) {
    T result ;
    MPI_Datatype MPI_T_type = MPI_traits<T>::get_MPI_type() ;
    MPI_Allreduce(&b, &result, 1, MPI_T_type, MPI_MAX, comm) ;
    return result ;
  }
  
  template<class T> T g_GLOBAL_MIN(T b, MPI_Comm comm=MPI_COMM_WORLD) {
    T result ;
    MPI_Datatype MPI_T_type = MPI_traits<T>::get_MPI_type() ;
    MPI_Allreduce(&b, &result, 1, MPI_T_type, MPI_MIN, comm) ;
    return result ;
  }
  
  // Collect largest interval of entitySet from all processors
  template<class T> genIntervalSet<T>  g_collectLargest(const  genIntervalSet<T>&e,MPI_Comm comm=MPI_COMM_WORLD) {
    int p = 1;
    MPI_Comm_size(comm,&p) ;
    
    // Else we compute set
    //First get largest interval
    std::pair<T, T> ivl_large(1,-1) ;
    
    const long esz = e.num_intervals() ;
    for(long i=0;i<esz;++i)
      if((ivl_large.second-ivl_large.first) < (e[i].second-e[i].first))
        ivl_large = e[i] ;
    
    std::vector<std::pair<T, T> > ivl_large_p(p) ;
    MPI_Datatype MPI_T_type = MPI_traits<T>::get_MPI_type() ;
    MPI_Allgather(&ivl_large,2,MPI_T_type,&(ivl_large_p[0]),2,MPI_T_type, comm) ;

    genIntervalSet<T> lset ;
    for(int i=0;i<p;++i)
      if(ivl_large_p[i].first <= ivl_large_p[i].second)
        lset += ivl_large_p[i] ;
    
    return lset ;
  }

  // Return union of all entitySets from all processors, the actual user interface is g_all_collect_entitySet()
  template<class T> genIntervalSet<T>  g_all_gather_entitySet(const  genIntervalSet<T> &e,MPI_Comm comm = MPI_COMM_WORLD ) {
    int p = 1;
    MPI_Comm_size(comm,&p) ;
    if(p == 1)
      return e ;
    
    int send_count = 2*e.num_intervals() ; //must be int*
    std::vector<int> recv_count(p) ;//must be int*
    MPI_Allgather(&send_count,1,MPI_INT,&recv_count[0],1,MPI_INT, comm) ;

    std::vector<int> recv_disp(p) ;//must be int*
    recv_disp[0] = 0 ;
    for(int i=1;i<p;++i)
      recv_disp[i] = recv_disp[i-1]+recv_count[i-1] ;
    int tot = recv_disp[p-1]+recv_count[p-1] ;
    if(tot == 0)
      return  genIntervalSet<T>::EMPTY ;
    std::vector<std::pair<T, T> > ivl_list(tot/2) ;
    std::vector<std::pair<T, T> > snd_list(send_count/2) ;
    for(long i=0;i<send_count/2;++i)
      snd_list[i] = e[i] ;

    MPI_Datatype MPI_T_type = MPI_traits<T>::get_MPI_type() ;
    MPI_Allgatherv(&(snd_list[0]),send_count,MPI_T_type,
		   &(ivl_list[0]),&(recv_count[0]), &(recv_disp[0]), MPI_T_type,
		   comm) ;

    std::sort(ivl_list.begin(),ivl_list.end(),g_spec_ival_compare<T>) ;
    genIntervalSet<T> tmp = ivl_list[0] ;
    for(size_t i=1;i<ivl_list.size();++i)
      tmp += ivl_list[i] ;
    return tmp ;
  }
  
  //Return union of all entitySets from all processors,
  template<class T> genIntervalSet<T>  g_all_collect_entitySet(const genIntervalSet<T> &e,MPI_Comm comm = MPI_COMM_WORLD ) {
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    // no operation for single processor
    if(p == 1)
      return e ;
    // Check to see if the result should be EMPTY or UNIVERSE
    int code = 0 ;
    if(e != genIntervalSet<T>::EMPTY)
      code = 1 ;
    if(e == ~(genIntervalSet<T>::EMPTY))
      code = 2 ;
    
    int result = 0;
    MPI_Allreduce(&code, &result, 1, MPI_INT, MPI_MAX, comm) ;
    code = result;
    
    //code = g_GLOBAL_MAX<int>(code,comm) ;
    
    if(code == 0) // All empty, so return empty
      return genIntervalSet<T>::EMPTY ;
    if(code == 2) // at least one UNIVERSE, so return UNIVERSE
      return ~genIntervalSet<T>::EMPTY ;
    
    // First collect largest intervals, most of the time this will get
    // us the final set.  Compute remainder to determine if we have
    // more work to do.  Repeat collecting the largest interval for a
    // fixed number of tries, then finally just gather the remainder
#ifdef VERBOSE
    stopWatch s ;
    s.start() ;
#endif
    
    genIntervalSet<T> lset = g_collectLargest<T>(e) ;
    genIntervalSet<T> rem = e-lset ;
    for(int i=0;i<4;++i) {
      T remsz = rem.num_intervals() ;
      T szmx = g_GLOBAL_MAX<T>(remsz,comm) ;
      if(szmx == 0) {
#ifdef VERBOSE
        debugout << "time to get lset = " << s.stop() << endl ;
#endif
        return lset ;
      }
      lset += g_collectLargest<T>(rem) ;
      rem -= lset ;
    }
#ifdef VERBOSE
    debugout << "time to get lset = " << s.stop() << endl ;
    s.start() ;
    
    debugout << "e="<< e.num_intervals() << ",rem=" << rem.num_intervals()
             << ",lset=" << lset.num_intervals() << endl ;
#endif
    genIntervalSet<T> remtot = g_all_gather_entitySet<T>(rem, comm) ;
#ifdef VERBOSE
    debugout << "time to gather rem = " << s.stop() << endl ;
#endif
    return lset + remtot ;
  }
  //Return union of all entitySets from all processors that belongs to this processor
  template<class T> genIntervalSet<T>  g_dist_collect_entitySet(const genIntervalSet<T> &inSet,
                                                                const std::vector<genIntervalSet<T> > &ptn,
                                                                MPI_Comm comm = MPI_COMM_WORLD ) {
    const int r = MPI_rank ;
    genIntervalSet<T> retval = inSet & ptn[r] ; //set from me
    // Check for empty and universal set
    int sbits = ((inSet != genIntervalSet<T>::EMPTY)?1:0)| ((retval != ptn[r])?2:0) ;
    int rbits = sbits ;
    MPI_Allreduce(&sbits, &rbits, 1, MPI_INT, MPI_BOR, comm) ;
    if((rbits & 1) == 0) // EMPTY set
      return genIntervalSet<T>::EMPTY ;
    if((rbits & 2) == 0) // UNIVERSE set
      return ptn[r] ;
    
    genIntervalSet<T> collect_set = inSet- ptn[r] ; //set for others


    genIntervalSet<T> all_set = g_all_collect_entitySet<T>(collect_set, comm); //set from all others
   
    return ((all_set&ptn[r]) + retval); //only return the set that belongs to me
  }

  
  //apply to both entitySet and gEntitySet
  template<class T>  std::vector<genIntervalSet<T> > g_all_collect_vectors(genIntervalSet<T> &e,MPI_Comm comm){
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    std::vector<genIntervalSet<T> > vset(p) ;
    if(p == 1) {
      vset[0] = e ;
      return vset ;
    }
    
    int send_count = 2*e.num_intervals() ;
    std::vector<int> recv_count(p) ;
    MPI_Allgather(&send_count,1,MPI_INT,&recv_count[0],1,MPI_INT,
                  comm) ;
    
    std::vector<int> recv_disp(p) ;
    recv_disp[0] = 0 ;
    for(int i=1;i<p;++i)
      recv_disp[i] = recv_disp[i-1]+recv_count[i-1] ;

    int tot = recv_disp[p-1]+recv_count[p-1] ;
    if(tot == 0)
      return vset ;
    std::vector<T> ivl_list(tot) ;
    std::vector<T> snd_list(send_count) ;
    for(int i=0;i<send_count/2;++i) {
      snd_list[i*2] = e[i].first ;
      snd_list[i*2+1] = e[i].second ;
    }
    MPI_Datatype MPI_T_type = MPI_traits<T>::get_MPI_type() ;
    MPI_Allgatherv(&(snd_list[0]),send_count,MPI_T_type,
		   &(ivl_list[0]),&(recv_count[0]), &(recv_disp[0]), MPI_T_type,
		   comm) ;
    
    for(int i = 0; i < p ; ++i) {
      int ind = recv_disp[i] ;
      for(int j=0;j<recv_count[i]/2;++j) {
        vset[i] += std::pair<T, T>(ivl_list[ind+j*2],ivl_list[ind+j*2+1]) ;
      }
    }
    return vset ;
  }

  
  template<class T>  std::vector<genIntervalSet<T> > g_all_collect_vectors(genIntervalSet<T> &e) {
    return g_all_collect_vectors<T>(e,MPI_COMM_WORLD) ;
  }
 
}

#endif
 
