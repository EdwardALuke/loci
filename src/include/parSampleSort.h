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
#ifndef PAR_SAMPLE_SORT_H
#define PAR_SAMPLE_SORT_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <distribute.h>
#include <mpi.h>

#include <vector>
#include <algorithm>

namespace Loci {
   template <class T> void balanceDistribution(std::vector<T> &list,
                                              MPI_Comm comm) {
    
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    if(p == 1) // if serial run, we are finished
      return ;

    MPI_Datatype bytearray ;
    MPI_Type_vector(sizeof(T),1,1,MPI_BYTE,&bytearray) ;
    MPI_Type_commit(&bytearray) ;

    int sz = list.size() ;
    std::vector<int> sizes(p) ;
    MPI_Allgather(&sz,1,MPI_INT,&sizes[0],1,MPI_INT,comm) ;
    long long tot_size = 0  ;
    for(int i=0;i<p;++i) {
      tot_size += sizes[i] ;
    }
    std::vector<int> bsizes(p) ;
    int psz = tot_size/p ;
    int rem = tot_size%p ;
    for(int i=0;i<p;++i) 
      bsizes[i] = psz + (i<rem?1:0) ;
#ifdef DEBUG
    long long sum = 0 ;
    for(int i=0;i<p;++i) 
      sum += bsizes[i] ;
    
    FATAL(sum != tot_size) ;
#endif
    int r=0 ;
    MPI_Comm_rank(comm,&r) ;
    std::vector<T> recv(bsizes[r]) ;
    std::vector<long long> sdist(p+1),rdist(p+1) ;
    sdist[0] = 0 ;
    rdist[0] = 0 ;
    for(int i=0;i<p;++i) {
      sdist[i+1] = sizes[i]+sdist[i] ;
      rdist[i+1] = bsizes[i]+rdist[i] ;
    }
    
    FATAL(sdist[p] != rdist[p]) ;
#ifdef DEBUG
    if(r == 0) {
      debugout << "balanceDistribution: ssplits=" ;
      for(int i=0;i<p+1;++i)  
        debugout << " " << sdist[i] ;
      debugout << endl ;
      debugout << "balanceDistribution: rsplits=" ;
      for(int i=0;i<p+1;++i)  
        debugout << " " << rdist[i] ;
      debugout << endl ;
    }
#endif
    
    // Perform Irecvs
    std::vector<MPI_Request> requests(p) ;
    int req = 0 ;
    long long i1 = rdist[r] ;
    long long i2 = rdist[r+1] ;
    for(int i=0;i<p;++i) {
      if(sdist[i]<i2 && sdist[i+1]>i1) { // intersection
        long long li = std::max(i1,sdist[i]) ;
        long long ri = std::min(i2,sdist[i+1]) ;
        int len = ri-li ;
        FATAL(len < 0) ;
        int s2 = li-rdist[r] ;
        FATAL(s2 < 0) ;
        if(i == r) { // local copy
          int s1 = li-sdist[r] ;
          FATAL(s1 < 0 || s2 < 0) ;
          FATAL(s2+len>bsizes[r]) ;
          FATAL(size_t(s1+len)>list.size()) ;
          for(int j=0;j<len;++j)
            recv[s2+j] = list[s1+j] ;
        } else {
          WARN(s2+len>bsizes[r]) ;
          if(len > 0)
            MPI_Irecv(&recv[s2],len,bytearray,i,55,comm,
                      &requests[req++]) ;
#ifdef DEBUG
          debugout << "balanceDistribution: recving " << len << " items from " << i << endl ;
#endif
        }
      }
    }

    // Perform sends
    i1 = sdist[r] ;
    i2 = sdist[r+1] ;
    for(int i=0;i<p;++i) {
      if(i != r && rdist[i]<i2 && rdist[i+1]>i1) { // intersection
        long long li = std::max(i1,rdist[i]) ;
        long long ri = std::min(i2,rdist[i+1]) ;
        int len = ri-li ;
        int s1 = li-sdist[r] ;
        FATAL(s1 < 0) ;
        FATAL(len < 0) ;
        FATAL(size_t(s1+len) > list.size()) ;
        if(len > 0) 
          MPI_Send(&list[s1],len,bytearray,i,55,comm) ;
#ifdef DEBUG
        debugout << "balanceDistribution: sending " << len << " items to " << i << endl ;
#endif
      }
    }

    if(req > 0) {
      std::vector<MPI_Status> status(p) ;
      MPI_Waitall(req,&requests[0],&status[0]) ;
    }
    list.swap(recv) ;
    MPI_Type_free(&bytearray) ;
  }
  
  // With user supplied comparison operator

  // Utility routine for sample sort
  template <class T,class Cmp>
  void parGetSplitters(std::vector<T> &splitters,
                       const std::vector<T> &input,
                       Cmp cmp,
                       MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    MPI_Datatype bytearray ;
    MPI_Type_vector(sizeof(T),1,1,MPI_BYTE,&bytearray) ;
    MPI_Type_commit(&bytearray) ;

    splitters = std::vector<T>(p-1) ;
    std::vector<T> allsplits(p*(p-1)) ;

    int nlocal = input.size() ;
    if(nlocal == 0) {
      std::cerr << "unable to sort with 0 elements on a processor"
                << std::endl ;
      Loci::Abort() ;
    }

    if(nlocal < p) {
      // If there aren't enough splitters, replicate the first entry
      // to fill in the rest.  May not be efficient, but the sort
      // should still work.
      for(int i=0;i<nlocal;++i)
        splitters[i] = input[i] ;
      for(int i=nlocal;i<p-1;++i)
        splitters[i] = input[0] ;
    } else {
      for(int i=1;i<p;++i) 
        splitters[i-1] = input[i*(nlocal/p)] ;
    }

    MPI_Allgather(&splitters[0],(p-1),bytearray,
                  &allsplits[0],(p-1),bytearray,comm) ;
    
    std::sort(allsplits.begin(),allsplits.end(),cmp) ;
    for(int i=1;i<p;++i)
      splitters[i-1] = allsplits[i*(p-1)] ;

    MPI_Type_free(&bytearray) ;
    //    splitters[p-1] = std::numeric_limits<T>::max() ;
    return ;
  }


  template <class T, class Cmp>
  void parSplitSort(std::vector<T> &list, std::vector<T> &splitters,
                    Cmp cmp, MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1) // if serial run, we are finished
      return ;
    if(int(splitters.size()) != p-1) {
      cerr << "parSplitSort passed invalid splitter" << endl ;
      Loci::Abort() ;
    }

    MPI_Datatype bytearray ;
    MPI_Type_vector(sizeof(T),1,1,MPI_BYTE,&bytearray) ;
    MPI_Type_commit(&bytearray) ;

    int s=0 ;
    std::vector<int> scounts(p,0) ;
    for(size_t i=0;i<list.size();++i)
      if(s == p-1 || cmp(list[i] , splitters[s]) ) 
        scounts[s]++ ;
      else {
        while((s!=p-1) && !cmp(list[i],splitters[s]))
          ++s ;
        scounts[s]++ ;
      }

    std::vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;

    std::vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;

    std::vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }
  
    int result_size = (rdispls[p-1]+rcounts[p-1]) ;
    std::vector<T> sorted_pnts(result_size) ;

    MPI_Alltoallv(&list[0],&scounts[0],&sdispls[0],bytearray,
                  &sorted_pnts[0],&rcounts[0],&rdispls[0],bytearray,
                  comm) ;

    list.swap(sorted_pnts) ;
    std::sort(list.begin(),list.end(),cmp) ;
    MPI_Type_free(&bytearray) ;
    return ;
  }

  template <class T, class Cmp>
  void parSampleSort(std::vector<T> &list, Cmp cmp, MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1) {
      std::sort(list.begin(),list.end(),cmp) ;
      return ;
    }
      

    long long lsz = int(list.size()) ;
    long long tsz = 0 ;
    MPI_Allreduce(&lsz,&tsz,1,MPI_LONG_LONG,MPI_SUM,comm) ;
    // Compute target number of processors needed to sort this set of numbers
    // ideally this will have p^2 numbers per processor
    if(tsz == 0)
      return ;

    int target_p = max(1,int(floor(pow(double(tsz),1./3.)))) ;

    if(target_p < p) { // reduce to subset of processors
#ifdef DEBUG
      debugout << "target_p=" << target_p << endl ;
#endif
      int r = 0 ;
      MPI_Comm_rank(comm,&r) ;
      int color = 1 ;
      MPI_Datatype bytearray ;
      MPI_Type_vector(sizeof(T),1,1,MPI_BYTE,&bytearray) ;
      MPI_Type_commit(&bytearray) ;
      if(r < target_p) {
        color = 0 ;
        size_t loc_size = int(lsz) ;
        std::vector<int> recv_sizes ;
        for(int i=r+target_p;i<p;i+=target_p) {
          int tmp = 0 ;
          MPI_Status stat ;
          MPI_Recv(&tmp,1,MPI_INT,i,1,comm,&stat) ;
          recv_sizes.push_back(tmp) ;
          loc_size += tmp ;
        }

        std::vector<T> nlist(loc_size) ;
        for(long long i=0;i<lsz;++i)
          nlist[i] = list[i] ;

        size_t loc = lsz ;
        int cnt = 0 ;
        for(int i=r+target_p;i<p;i+=target_p) {

          MPI_Status stat ;
          MPI_Recv(&nlist[loc],recv_sizes[cnt],bytearray,i,2,comm,&stat) ;
          loc += recv_sizes[cnt] ;
          cnt++ ;
        }
        list.swap(nlist) ;
      } else {
        int dest = r % target_p ;
        int sz = list.size() ;
        MPI_Send(&sz,1,MPI_INT,dest,1,comm) ;
        MPI_Send(&list[0],sz,bytearray,dest,2,comm) ;
        std::vector<T> nlist ;
        list.swap(nlist) ;
      }
      MPI_Type_free(&bytearray) ;
      // Split off smaller group of processors
      MPI_Comm subset ;
      MPI_Comm_split(comm,color,r,&subset) ;

      if(color == 0) {
        // sort over subset of processors
        parSampleSort(list, cmp, subset) ;
      }
      // Free communicatiors
      MPI_Comm_free(&subset) ;
    } else {
      int sz = list.size() ;
      int msz = 0 ;
      MPI_Allreduce(&sz,&msz,1,MPI_INT,MPI_MIN,comm) ;
      if(msz < p)
        balanceDistribution(list,comm) ;

      // First sort list locally
      std::sort(list.begin(),list.end(),cmp) ;
      
      if(p == 1) // if serial run, we are finished
        return ;
      
      std::vector<T> splitters ;
      parGetSplitters(splitters,list,cmp,comm) ;
      parSplitSort(list,splitters,cmp,comm) ;
    }

  }
  // Now same code with default operator

  // Utility routine for sample sort
  template <class T> void parGetSplitters(std::vector<T> &splitters,
                                          const std::vector<T> &input,
                                          MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    splitters = std::vector<T>(p-1) ;
    std::vector<T> allsplits(p*(p-1)) ;
    
    int nlocal = input.size() ;

    if(nlocal == 0) {
      cerr << "unable to sort with 0 elements on a processor" << endl ;
      Loci::Abort() ;
    }
    if(nlocal < p) {
      // If there aren't enough splitters, replicate the first entry
      // to fill in the rest.  May not be efficient, but the sort
      // should still work.
      for(int i=0;i<nlocal;++i)
        splitters[i] = input[i] ;
      for(int i=nlocal;i<p;++i)
        splitters[i] = input[0] ;
    } else 
      for(int i=1;i<p;++i) 
        splitters[i-1] = input[i*(nlocal/p)] ;

    MPI_Datatype bytearray ;
    MPI_Type_vector(sizeof(T),1,1,MPI_BYTE,&bytearray) ;
    MPI_Type_commit(&bytearray) ;

    MPI_Allgather(&splitters[0],(p-1),bytearray,
                  &allsplits[0],(p-1),bytearray,comm) ;

    MPI_Type_free(&bytearray) ;

    std::sort(allsplits.begin(),allsplits.end()) ;

    for(int i=1;i<p;++i)
      splitters[i-1] = allsplits[i*(p-1)] ;
    //    splitters[p-1] = std::numeric_limits<T>::max() ;
    return ;
  }


  template <class T> void parSplitSort(std::vector<T> &list,
                                       std::vector<T> &splitters,
                                       MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1) // if serial run, we are finished
      return ;

    int s=0 ;
    std::vector<int> scounts(p,0) ;
    for(size_t i=0;i<list.size();++i)
      if(s == p-1 || (list[i] < splitters[s]) ) 
        scounts[s]++ ;
      else {
        while((s!=p-1) && !(list[i] < splitters[s]))
          ++s ;
        scounts[s]++ ;
      }

    std::vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;

    std::vector<int> rcounts(p) ;


    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;

    std::vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }
  
    int result_size = (rdispls[p-1]+rcounts[p-1]) ;

    std::vector<T> sorted_pnts(result_size) ;
    
    MPI_Datatype bytearray ;
    MPI_Type_vector(sizeof(T),1,1,MPI_BYTE,&bytearray) ;
    MPI_Type_commit(&bytearray) ;

    MPI_Alltoallv(&list[0],&scounts[0],&sdispls[0],bytearray,
                  &sorted_pnts[0],&rcounts[0],&rdispls[0],bytearray,
                  comm) ;
    
    MPI_Type_free(&bytearray) ;
    list.swap(sorted_pnts) ;
    std::sort(list.begin(),list.end()) ;
    return ;
  }

  template <class T> void parSampleSort(std::vector<T> &list, MPI_Comm comm) {
    return parSampleSort(list, std::less<T>(), comm) ;
  }

}

#endif
