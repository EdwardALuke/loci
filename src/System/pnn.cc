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
#include "pnn.h"
using std::vector ;
#include <iostream>
using std::cout ;
using std::cerr ;
using std::endl ;
#include <sstream>
#include <algorithm>
using std::max ;
using std::min ;
using std::sort ;
using std::pair ;

namespace Loci {
  using namespace Loci::kdTree ;

  void ParGetBounds(kd_tree::bounds &bnd,
                    const vector<kd_tree::coord_info> &pnts,
                    MPI_Comm comm) {
    kd_tree::bounds bndi ;
    for(int d=0;d<3;++d) {
      bndi.minc[d] = .25*std::numeric_limits<double>::max() ;
      bndi.maxc[d] = -.25*std::numeric_limits<double>::max() ;
    }
    for(size_t i=0;i<pnts.size();++i) 
      for(int d=0;d<3;++d) {
        bndi.maxc[d] = max(bndi.maxc[d],pnts[i].coords[d]) ;
        bndi.minc[d] = min(bndi.minc[d],pnts[i].coords[d]) ;
      }
    MPI_Allreduce(&bndi.maxc[0],&bnd.maxc[0],3,MPI_DOUBLE,MPI_MAX,comm) ;
    MPI_Allreduce(&bndi.minc[0],&bnd.minc[0],3,MPI_DOUBLE,MPI_MIN,comm) ;
  }

  namespace {
    bool comparex(const kd_tree::coord_info &v1,
                  const kd_tree::coord_info &v2) {
      return v1.coords[0] < v2.coords[0] ;
    }
    bool comparey(const kd_tree::coord_info &v1,
                  const kd_tree::coord_info &v2) {
      return v1.coords[1] < v2.coords[1] ;
    }
    bool comparez(const kd_tree::coord_info &v1,
                  const kd_tree::coord_info &v2) {
      return v1.coords[2] < v2.coords[2] ;
    }
  }
              

  void ParGetSplitters(vector<double> &splitters,
                       const vector<kd_tree::coord_info> &tpnts,
                       const vector<kd_tree::coord_info> &spnts,
                       int dim,
                       MPI_Comm comm) {
    const size_t tsz = tpnts.size() ;
    const size_t ssz = spnts.size() ;
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;

    vector<double> vals(tsz+ssz) ;
    size_t k=0 ;
    for(size_t i=0;i<tsz;++i)
      vals[k++] = tpnts[i].coords[dim] ;
    for(size_t i=0;i<ssz;++i)
      vals[k++] = spnts[i].coords[dim] ;
    sort(vals.begin(),vals.end()) ;
    splitters = vector<double>(p,-1.) ;
    vector<double> allsplits(p*(p-1),-1.) ;

    int nlocal = tsz+ssz ;
    for(int i=1;i<p;++i) 
      splitters[i-1] = vals[(i*nlocal)/p] ;

    MPI_Allgather(&splitters[0],p-1,MPI_DOUBLE,
                  &allsplits[0],p-1,MPI_DOUBLE,comm) ;

    sort(allsplits.begin(),allsplits.end()) ;
    for(int i=1;i<p;++i)
      splitters[i-1] = allsplits[i*(p-1)] ;
    splitters[p-1] = std::numeric_limits<double>::max() ;
    return ;
  }

  void ParSplitCoord(vector<kd_tree::coord_info> &pnts,
                     const vector<double> &splitters,
                     int dim,
                     MPI_Comm comm) {
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1)
      return ;
    int me = 0 ;
    MPI_Comm_rank(comm,&me) ;

    switch(dim) {
    case 0:
      sort(pnts.begin(),pnts.end(),comparex) ;
      break ;
    case 1:
      sort(pnts.begin(),pnts.end(),comparey) ;
      break ;
    case 2:
      sort(pnts.begin(),pnts.end(),comparez) ;
      break ;
    }
    int s=0 ;
    vector<int> scounts(p,0) ;
    for(size_t i=0;i<pnts.size();++i)
      if(pnts[i].coords[dim] < splitters[s]) 
        scounts[s]++ ;
      else {
        while(pnts[i].coords[dim] >= splitters[s])
          ++s ;
        scounts[s]++ ;
      }

    for(size_t i=0;i<scounts.size();++i) 
      scounts[i]*=sizeof(kd_tree::coord_info) ;

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;

    vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }
  
    int result_size = (rdispls[p-1]+rcounts[p-1])/sizeof(kd_tree::coord_info) ;

    vector<kd_tree::coord_info> sorted_pnts(result_size) ;

    MPI_Alltoallv(&pnts[0],&scounts[0],&sdispls[0],MPI_BYTE,
                  &sorted_pnts[0],&rcounts[0],&rdispls[0],MPI_BYTE,
                  comm) ;

    pnts.swap(sorted_pnts) ;
    return ;
  }

  void ParSortCoord(vector<kd_tree::coord_info> &tpnts,
                    vector<kd_tree::coord_info> &spnts,
                    int dim,
                    MPI_Comm comm) {
  
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1)
      return ;
    int me = 0 ;
    MPI_Comm_rank(comm,&me) ;

    vector<double> splitters(p,-1.) ;

    ParGetSplitters(splitters,tpnts,spnts,dim,comm) ;
    ParSplitCoord(tpnts,splitters,dim,comm) ;
    ParSplitCoord(spnts,splitters,dim,comm) ;
  }

  void ORB_Partition(vector<kd_tree::coord_info> &tpnts,
                     vector<kd_tree::coord_info> &spnts,
                     MPI_Comm comm) {

    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    if(p == 1)
      return ;
    int myid = 0;
    MPI_Comm_rank(comm,&myid) ;
  
    kd_tree::bounds bnd1 ;
    ParGetBounds(bnd1,tpnts,comm) ;
    kd_tree::bounds bnd2 ;
    ParGetBounds(bnd2,spnts,comm) ;
    kd_tree::bounds bnd ;
    for(int d=0;d<3;++d) {
      bnd.minc[d] = min(bnd1.minc[d],bnd2.minc[d]) ;
      bnd.maxc[d] = max(bnd1.maxc[d],bnd2.maxc[d]) ;
    }
  
    int dim = 0 ;
    double dx = bnd.maxc[0]-bnd.minc[0] ;
    double dy = bnd.maxc[1]-bnd.minc[1] ;
    double dz = bnd.maxc[2]-bnd.minc[2] ;
    if(dy >= dx && dy >= dz)
      dim = 1 ;
    if(dz >= dx && dz >= dy)
      dim = 2 ;

    ParSortCoord(tpnts,spnts,dim,comm) ;

    if(p <= 3)
      return ;
    MPI_Comm sub_com ;
    int color = myid >= (p/2) ;
    MPI_Comm_split(comm,color,myid,&sub_com) ;
    ORB_Partition(tpnts,spnts,sub_com) ;
    MPI_Comm_free(&sub_com) ;
  }

  inline bool pt_in_box(coord3d v,const kd_tree::bounds &b1) {
    return
      (v[0] >= b1.minc[0] && v[0] <= b1.maxc[0]) &&
      (v[1] >= b1.minc[1] && v[1] <= b1.maxc[1]) &&
      (v[2] >= b1.minc[2] && v[2] <= b1.maxc[2]) ;
  }
  inline coord3d pselect(const kd_tree::bounds &b, int sel) {
    coord3d v ;
    v[0] = ((sel&1)==0)?b.maxc[0]:b.minc[0] ;
    v[1] = ((sel&2)==0)?b.maxc[1]:b.minc[1] ;
    v[2] = ((sel&4)==0)?b.maxc[2]:b.minc[2] ;
    return v ;
  }
    
  bool box_intersect(const kd_tree::bounds &b1, const kd_tree::bounds &b2) {
    bool test = false ;
    for(int i=0;i<8;++i)
      test = test || pt_in_box(pselect(b1,i),b2) || pt_in_box(pselect(b2,i),b1) ;
    return test ;
  }

  void parallelNN_DistributedSearch(const vector<coord3d> &target_pnts,
                                    const vector<int> &target_ids,
                                    const vector<coord3d> &search_pnts,
                                    vector<int> &close_pt,
                                    MPI_Comm comm) {

    int p = 0 ;
    int myid = 0 ;
    /* Get the number of processors and my processor identification */
    MPI_Comm_size(comm,&p) ;
    MPI_Comm_rank(comm,&myid) ;

  
    // Make local copies of target points for partitioning
    vector<kd_tree::coord_info> tcopy(target_pnts.size()) ;
    for(size_t i=0;i<target_pnts.size();++i) {
      tcopy[i].coords[0] = target_pnts[i][0] ;
      tcopy[i].coords[1] = target_pnts[i][1] ;
      tcopy[i].coords[2] = target_pnts[i][2] ;
      tcopy[i].id = target_ids[i] ;
    }

    // Make local copies of search points for partitioning
    // Give the search points id's so we can return them to the correct
    // processor when the search is complete
    vector<int> search_sizes(p) ;
    int lM = search_pnts.size() ;

    MPI_Allgather(&lM,1,MPI_INT,&search_sizes[0],1,MPI_INT,comm) ;
    vector<int> search_ids(p+1) ;
    search_ids[0] = 0 ;
    for(int i=1;i<p;++i)
      search_ids[i] = search_ids[i-1]+search_sizes[i-1] ;
    search_ids[p] = std::numeric_limits<int>::max() ;
  
    vector<kd_tree::coord_info> scopy(search_pnts.size()) ;
    for(size_t i=0;i<search_pnts.size();++i) {
      scopy[i].coords[0] = search_pnts[i][0] ;
      scopy[i].coords[1] = search_pnts[i][1] ;
      scopy[i].coords[2] = search_pnts[i][2] ;
      scopy[i].id = i+search_ids[myid] ;
    }

    // Use Orthogonal Recursive Bisection to partition target and search points
    ORB_Partition(tcopy,scopy,comm) ;

    int lN = tcopy.size() ;
    lM = scopy.size() ;

    kd_tree::bounds bndi ;
    for(int d=0;d<3;++d) {
      bndi.minc[d] = .25*std::numeric_limits<double>::max() ;
      bndi.maxc[d] = -.25*std::numeric_limits<double>::max() ;
    }
    for(int i=0;i<lN;++i) 
      for(int d=0;d<3;++d) {
        bndi.maxc[d] = max(bndi.maxc[d],tcopy[i].coords[d]) ;
        bndi.minc[d] = min(bndi.minc[d],tcopy[i].coords[d]) ;
      }
    for(int i=0;i<lM;++i) 
      for(int d=0;d<3;++d) {
        bndi.maxc[d] = max(bndi.maxc[d],scopy[i].coords[d]) ;
        bndi.minc[d] = min(bndi.minc[d],scopy[i].coords[d]) ;
      }

    vector<kd_tree::bounds> pbnds(p) ;
    MPI_Allgather(&bndi,6,MPI_DOUBLE,&pbnds[0],6,MPI_DOUBLE,comm) ;
  
    int check = lN ;

    vector<int> lNp(p) ;
    MPI_Allgather(&lN,1,MPI_INT,&lNp[0],1,MPI_INT,comm) ;
    for(int i=0;i<p;++i)
      check = min(check,lNp[i]) ;
    bool no_tpnts = (lN == 0) ;
    if(check == 0) {
      // Some of the processors have no search points, we need to give them some
      int fromp = -1 ;
      if(lN == 0) {
        double close_dist = 1e100 ;
        coord3d me ;
        for(int d=0;d<3;++d)
          me[d] = .5*( pbnds[myid].minc[d]+pbnds[myid].maxc[d] ) ;
      

        for(int i=0;i<p;++i) 
          if(lNp[i] != 0) {
            coord3d v ;
            for(int d=0;d<3;++d)
              v[d] = .5*( pbnds[i].minc[d]+pbnds[i].maxc[d] ) ;
            if(dist2(me,v) < close_dist) {
              close_dist = dist2(me,v) ;
              fromp = i ;
            }
          }
      }
      vector<int> preq(p) ;
      MPI_Allgather(&fromp,1,MPI_INT,&preq[0],1,MPI_INT,comm) ;
      if(lN == 0) {
        int sz = lNp[fromp] ;
        vector<kd_tree::coord_info> tmp(sz) ;
        sz *= sizeof(kd_tree::coord_info) ;
        MPI_Status status ;
        MPI_Recv(&tmp[0],sz,MPI_BYTE,fromp,0,comm,&status) ;
        tcopy.swap(tmp) ;
        lN = tcopy.size() ;
      } else {
        for(int i=0;i<p;++i)
          if(preq[i] == myid) {
            int sz = lN*sizeof(kd_tree::coord_info) ;
            MPI_Send(&tcopy[0],sz,MPI_BYTE,i,0,comm) ;
          }
      }
    }  

    vector<double> rmin(lM,1e65) ;
    vector<pair<int,int> > close_pair(lM) ;


    // First search with the points that we have
    vector<kd_tree::coord_info> targ_send = tcopy ;

    kd_tree kd(targ_send) ;
    
    for(int i=0;i<lM;++i) {
      close_pair[i].first=scopy[i].id ;
      close_pair[i].second = kd.find_closest(scopy[i].coords,rmin[i]) ;
    }

    if(p > 1) { // If there are other processors, communicate points to
      // complete the search.

      // Compute bounding box of space to be searched
      for(int d=0;d<3;++d) {
        bndi.minc[d] = .25*std::numeric_limits<double>::max() ;
        bndi.maxc[d] = -.25*std::numeric_limits<double>::max() ;
      }
      vector<int> workingSet ; // Points that still need to search are in the
      // working set.
      for(int i=0;i<lM;++i) {
        double r = sqrt(rmin[i]) ;
        coord3d v = scopy[i].coords ;
        bool inbounds =
          (v[0] + r < pbnds[myid].maxc[0] &&
           v[0] - r > pbnds[myid].minc[0]) &&
          (v[1] + r < pbnds[myid].maxc[1] &&
           v[1] - r > pbnds[myid].minc[1]) &&
          (v[2] + r < pbnds[myid].maxc[2] &&
           v[2] - r > pbnds[myid].minc[2])  ;
        if(!inbounds) {
          workingSet.push_back(i) ;
          for(int d=0;d<3;++d) {
            bndi.minc[d] = min(bndi.minc[d],v[d]-r) ;
            bndi.maxc[d] = max(bndi.maxc[d],v[d]+r) ;
          }
        }
      }

      // Communicate bounds request to other processors
      vector<kd_tree::bounds> bnd_req(p) ;
      MPI_Allgather(&bndi,6,MPI_DOUBLE,&bnd_req[0],6,MPI_DOUBLE,comm) ;


      vector<int> scounts(p,0) ;
      vector<kd_tree::coord_info> pntlist ;
      if(!no_tpnts) 
        for(int i=0;i<p;++i) {
          int bsz = pntlist.size() ;
          kd.find_box(pntlist,bnd_req[i]) ;
          scounts[i] = pntlist.size()-bsz ;
        }
      
      for(size_t i=0;i<scounts.size();++i) 
        scounts[i]*=sizeof(kd_tree::coord_info) ;
    
      vector<int> sdispls(p) ;
      sdispls[0] = 0 ;
      for(int i=1;i<p;++i)
        sdispls[i] = sdispls[i-1]+scounts[i-1] ;
    
      vector<int> rcounts(p) ;
      MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;
    
      vector<int> rdispls(p) ;
      rdispls[0] = 0 ;
      for(int i=1;i<p;++i) {
        rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
      }
    
      int result_size = (rdispls[p-1]+rcounts[p-1])/sizeof(kd_tree::coord_info) ;
    
      vector<kd_tree::coord_info> pcollect(result_size) ;
    
      MPI_Alltoallv(&pntlist[0],&scounts[0],&sdispls[0],MPI_BYTE,
                    &pcollect[0],&rcounts[0],&rdispls[0],MPI_BYTE,
                    comm) ;
    
      kd_tree kd_xfer(pcollect) ;
    
      for(size_t ii=0;ii<workingSet.size();++ii) {
        int i=workingSet[ii] ;
        int id = kd_xfer.find_closest(scopy[i].coords,rmin[i]) ;
        if(id != std::numeric_limits<int>::min()) 
          close_pair[i].second = id ;
      }
    }
    

    sort(close_pair.begin(),close_pair.end()) ;

    vector<int> scounts(p,0) ;
    int s = 0 ;
    for(size_t i=0;i<close_pair.size();++i)
      if(close_pair[i].first < search_ids[s+1])
        scounts[s]++ ;
      else {
        while(close_pair[i].first >= search_ids[s+1])
          ++s ;
        scounts[s]++ ;
      }

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i) 
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;

    vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }

    int rsize = rdispls[p-1]+rcounts[p-1] ;

    if(size_t(rsize)!= close_pt.size()) {
      cerr << "rsize = " << rsize << ", close_pt.size() = " << close_pt.size()
           << endl ;
      cerr << "confused close_pt redistribution" << endl ;
    }

    vector<pair<int,int> > rclose(rsize) ;
    MPI_Alltoallv(&close_pair[0],&scounts[0],&sdispls[0],MPI_2INT,
                  &rclose[0],&rcounts[0],&rdispls[0],MPI_2INT,
                  comm) ;
  
    // Redistribute close_pt
    for(int i=0;i<rsize;++i) {
      close_pt[rclose[i].first-search_ids[myid]] = rclose[i].second ;
    }
  }

  void parallelNN_GatherSearch(const vector<coord3d> &target_pnts,
                               const vector<int> &target_ids,
                               const vector<coord3d> &search_pnts,
                               vector<int> &close_pt,
                               MPI_Comm comm) {

    int p = 0 ;
    int myid = 0 ;
    /* Get the number of processors and my processor identification */
    MPI_Comm_size(comm,&p) ;
    MPI_Comm_rank(comm,&myid) ;

    int tsz = target_pnts.size() ;
    vector<int> rcounts(p) ;
    MPI_Allgather(&tsz,1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;
    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }
    int rsize = rdispls[p-1]+rcounts[p-1] ;

    vector<coord3d> tpnts(rsize);
    vector<int> tids(rsize) ;
    MPI_Allgatherv((void *)&target_ids[0],tsz,MPI_INT,
                   &tids[0],&rcounts[0],&rdispls[0],
                   MPI_INT,comm) ;
    for(int i=0;i<p;++i)
      rcounts[i] *= 3 ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }
    rsize = rdispls[p-1]+rcounts[p-1] ;


    MPI_Allgatherv((void *)&target_pnts[0],tsz*3,MPI_DOUBLE,
                   &tpnts[0],&rcounts[0],&rdispls[0],
                   MPI_DOUBLE,comm) ;
    kd_tree kd(tpnts,tids) ;

    int lM = search_pnts.size() ;
    for(int i=0;i<lM;++i) 
      close_pt[i] = kd.find_closest(search_pnts[i]) ;
  }

  void parallelNN_SampleSearch(const vector<coord3d> &target_pnts,
                               const vector<int> &target_ids,
                               const vector<coord3d> &search_pnts,
                               vector<int> &close_pt,
                               MPI_Comm comm) {

    int p = 0 ;
    int myid = 0 ;
    /* Get the number of processors and my processor identification */
    MPI_Comm_size(comm,&p) ;
    MPI_Comm_rank(comm,&myid) ;

    int lM = search_pnts.size() ;

    if(p==1) {
      kd_tree kd(target_pnts,target_ids) ;
      for(int i=0;i<lM;++i) 
        close_pt[i] = kd.find_closest(search_pnts[i]) ;
      return ;
    }

    // Sample 100000 points total max, but at least 5 per processor.
    // If target_pnts is smaller than this sample all the points
    const int sample_size = max(min(max(100000/p,5),(int)target_pnts.size()),1);
    int samp_freq = max(int(target_pnts.size()/sample_size),1) ;
    int nsamples = target_pnts.size()/samp_freq ;

    vector<double> rmin(lM,1e65) ;

    {// First sample the targets
      vector<coord3d> sample_pts(nsamples) ;
      vector<int> sample_ids(nsamples) ;
      for(int i=0;i<nsamples;++i) {
        sample_pts[i] = target_pnts[i*samp_freq] ;
        sample_ids[i] = target_ids[i*samp_freq] ;
      }

      int tsz = nsamples ;
      vector<int> rcounts(p) ;
      MPI_Allgather(&tsz,1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;
      vector<int> rdispls(p) ;
      rdispls[0] = 0 ;
      for(int i=1;i<p;++i) {
        rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
      }
      int rsize = rdispls[p-1]+rcounts[p-1] ;
      
      vector<coord3d> tpnts(rsize);
      vector<int> tids(rsize) ;
      MPI_Allgatherv((void *)&sample_ids[0],tsz,MPI_INT,
                     &tids[0],&rcounts[0],&rdispls[0],
                     MPI_INT,comm) ;
      for(int i=0;i<p;++i)
        rcounts[i] *= 3 ;
      rdispls[0] = 0 ;
      for(int i=1;i<p;++i) {
        rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
      }
      rsize = rdispls[p-1]+rcounts[p-1] ;
      
      
      MPI_Allgatherv((void *)&sample_pts[0],tsz*3,MPI_DOUBLE,
                     &tpnts[0],&rcounts[0],&rdispls[0],
                     MPI_DOUBLE,comm) ;
      
      kd_tree kd(tpnts,tids) ;

      for(int i=0;i<lM;++i) 
        close_pt[i] = kd.find_closest(search_pnts[i],rmin[i]) ;
    }

    // Compute remaining target points
    vector<kd_tree::coord_info> remains(target_pnts.size()-nsamples) ;
    int pt = 0 ;
    int cnt = 0 ;
    const int max_pt = nsamples*samp_freq ;
    for(int i=0;i<int(target_pnts.size());++i) 
      if(i==pt && pt < max_pt) {
        pt += samp_freq ;
      } else {
        remains[cnt].coords[0] = target_pnts[i][0] ;
        remains[cnt].coords[1] = target_pnts[i][1] ;
        remains[cnt].coords[2] = target_pnts[i][2] ;
        remains[cnt].id = target_ids[i] ;
        cnt++ ;
      }


    // Compute bounding box of space to be searched
    kd_tree::bounds bndi ;
    for(int d=0;d<3;++d) {
      bndi.minc[d] = .25*std::numeric_limits<double>::max() ;
      bndi.maxc[d] = -.25*std::numeric_limits<double>::max() ;
    }
    for(int i=0;i<lM;++i) {
      double r = sqrt(rmin[i]) ;
      coord3d v = search_pnts[i] ;
      for(int d=0;d<3;++d) {
        bndi.minc[d] = min(bndi.minc[d],v[d]-r) ;
        bndi.maxc[d] = max(bndi.maxc[d],v[d]+r) ;
      }
    }

    // Communicate bounds request to other processors
    vector<kd_tree::bounds> bnd_req(p) ;
    MPI_Allgather(&bndi,6,MPI_DOUBLE,&bnd_req[0],6,MPI_DOUBLE,comm) ;
    

    // Now communicate points that processors need to complete NN search
    vector<int> scounts(p,0) ;
    vector<kd_tree::coord_info> pntlist ;
    if(remains.size() != 0) {
      kd_tree kdr(remains) ;
      for(int i=0;i<p;++i) { // Find points in i'th processor bounding box
        int bsz = pntlist.size() ;
        kdr.find_box(pntlist,bnd_req[i]) ;
        scounts[i] = pntlist.size()-bsz ;
      }
    }
    for(size_t i=0;i<scounts.size();++i) 
      scounts[i]*=sizeof(kd_tree::coord_info) ;
    
    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;
    
    vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;
    
    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }
    
    int result_size = (rdispls[p-1]+rcounts[p-1])/sizeof(kd_tree::coord_info) ;

    vector<kd_tree::coord_info> pcollect(result_size) ;
    
    MPI_Alltoallv(&pntlist[0],&scounts[0],&sdispls[0],MPI_BYTE,
                  &pcollect[0],&rcounts[0],&rdispls[0],MPI_BYTE,
                  comm) ;

    // Now we have collected all remaining points we need to complete the
    // search
    kd_tree kd_xfer(pcollect) ;

    // Check to see if we have a updated minimum
    for(int i=0;i<lM;++i) {
      int id = kd_xfer.find_closest(search_pnts[i],rmin[i]) ;
      if(id != std::numeric_limits<int>::min()) 
        close_pt[i] = id ;
    }

  }

  void parallelNearestNeighbors(const vector<coord3d> &target_pnts,
                                const vector<int> &target_ids,
                                const vector<coord3d> &search_pnts,
                                vector<int> &close_pt,
                                MPI_Comm comm, bool rebalance) {

    int num_targets = 0 ;
    int lsz = target_pnts.size() ;
    MPI_Allreduce(&lsz,&num_targets,1,MPI_INT,MPI_SUM,comm) ;
    
    if(num_targets <= 100000) {
      parallelNN_GatherSearch(target_pnts,target_ids,search_pnts,close_pt,
                              comm) ;
    } else {

      if(rebalance) {
        // Maybe a re-balanced sample search would be better
        parallelNN_DistributedSearch(target_pnts,target_ids,search_pnts,
                                     close_pt, comm) ;
      } else {
        parallelNN_SampleSearch(target_pnts,target_ids,search_pnts,
                                close_pt, comm) ;
      }
    }

  }
}
