//#############################################################################
//#
//# Copyright 2008, 2025, Mississippi State University
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
#ifndef SIMPLE_PARTITION_LONG_H
#define SIMPLE_PARTITION_LONG_H
#include <Tools/intervalSet.h>
#include <mpi.h>

using std::vector;

namespace Loci{

  //return num of entities of each processor
  template<class T> vector<T> g_simple_partition_dist(T mn, T mx, int p) {
    vector<T> nums(p) ;
    T n = mx-mn+1 ;
    
    if(p==1){
      nums[0] = n;
    }else{
      T dn = n/p ; // divisor
      T rn = n%p ; // remainder
  
      for(int i=0;i<p;++i) {
        int j = p - i - 1 ;
        nums[i]= dn+((j<rn)?1:0) ;
      }
    }
    return nums ;
  }
  
  
  //return the start number of each processor  
  template<class T> vector<T> g_simple_partition_vec(T mn, T mx, int p) {
    vector<T> nums = g_simple_partition_dist<T>(mn, mx, p);
    vector<T> starts(p+1) ;
  
    T start = mn ;
    starts[0] = start ;
    for(int i=0;i<p;++i) {
      start += nums[i];
      starts[i+1] = start ;
    }
    FATAL(start != mx+1) ;
    return starts ;
  }
  


  //return the entitiySet of each processor
  template<class T> vector<genIntervalSet<T> > g_simple_partition(T mn, T mx, MPI_Comm comm) {
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    vector<T> pl = g_simple_partition_vec(mn,mx,p) ;
    vector<genIntervalSet<T> > ptn(p) ;
    for(int i=0;i<p;++i)
      ptn[i] = std::pair<T, T>(pl[i],pl[i+1]-1) ;
    return ptn ;
  }

  //return the entitiySet of each processor
  template<class T> vector<genIntervalSet<T> > g_simple_partition(T num_entities, MPI_Comm comm) {
    T mn = 0;
    T mx = num_entities -1;
    return g_simple_partition(mn,mx,comm) ;
  }
  //for non continuous keys
  template<class T> vector<genIntervalSet<T> > g_simple_partition( genIntervalSet<T> all_keys ,MPI_Comm comm ) {
    int p = 1 ;
    MPI_Comm_size(comm,&p);
    vector<genIntervalSet<T> > ptn(p) ;
    size_t num_keys = all_keys.size();
    vector<genIntervalSet<T> > result(p);
    //start and end index of each process
    vector<T> nums = g_simple_partition_vec<gEntity>(0, num_keys-1, p);
    int index = 0;
    int pid = 0;
    for(typename genIntervalSet<T>::const_iterator ei = all_keys.begin(); ei != all_keys.end(); ei++){
      if(index < nums[pid+1]) result[pid] += *ei;
      else{
        pid++;
        result[pid] += *ei;
      }
      index++;
    }
    return result;
  }
  
 
}
#endif
