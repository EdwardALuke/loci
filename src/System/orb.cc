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
#include <Loci.h>
#include <parSampleSort.h>
#include <Loci_types.h>
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

using std::max ;
using std::min ;

namespace Loci {
  typedef pair<vector3d<float>, int> vpair ;

  inline bool vpairxcmp(const vpair &x1, const vpair &x2) {
    return x1.first.x < x2.first.x ;
  }
  inline bool vpairycmp(const vpair &x1, const vpair &x2) {
    return x1.first.y < x2.first.y ;
  }
  inline bool vpairzcmp(const vpair &x1, const vpair &x2) {
    return x1.first.z < x2.first.z ;
  }
  
  void ORBSort(vector<vpair> &pnts,MPI_Comm comm) {
    vector3d<float> mx(-1e33,-1e33,-1e33),mn(1e33,1e33,1e33) ;
    int sz = pnts.size() ;
    for(int i=0;i<sz;++i) {
      mx.x = max(mx.x,pnts[i].first.x) ;
      mx.y = max(mx.y,pnts[i].first.y) ;
      mx.z = max(mx.z,pnts[i].first.z) ;
      mn.x = min(mx.x,pnts[i].first.x) ;
      mn.y = min(mx.y,pnts[i].first.y) ;
      mn.z = min(mx.z,pnts[i].first.z) ;
    }
    vector3d<float> gmx,gmn ;
    MPI_Allreduce(&mx.x,&gmx.x,3,MPI_FLOAT,MPI_MAX,comm) ;
    MPI_Allreduce(&mn.x,&gmn.x,3,MPI_FLOAT,MPI_MIN,comm) ;
    vector3d<float> dx = gmx-gmn ;
    float d1,d2 ;
    if(dx.x >= dx.y && dx.x >= dx.z) { // X split
      d1 = dx.x ;
      d2 = max(dx.y,dx.z) ;
      parSampleSort(pnts,vpairxcmp,comm) ;
    } else if(dx.y >= dx.x && dx.y >= dx.z) { // Y split
      d1 = dx.y ;
      d2 = max(dx.x,dx.z) ;
      parSampleSort(pnts,vpairycmp,comm) ;
    } else { // Z split
      d1 = dx.z ;
      d2 = max(dx.x,dx.y) ;
      parSampleSort(pnts,vpairzcmp,comm) ;
    }
    balanceDistribution(pnts,comm) ;

    int p=0 ;
    MPI_Comm_size(comm,&p) ;
    int r=0 ;
    MPI_Comm_rank(comm,&r) ;

    // How much will we divide in this direction
    int ratio = max(4,int(min(ceil(d1/(d2+1e-33)),20.0))) ;

    // Number of processors per subdivision
    int nump = p/ratio ;
    if(nump<2) // If the subdivisions too small, stop partitioning
      return ;

    // Compute color for this processor to partition into ratio parts
    int color = min(ratio-1,(r*ratio)/p) ;

    MPI_Comm sub_com ;
    MPI_Comm_split(comm,color,r,&sub_com) ;
    ORBSort(pnts,sub_com) ;
    MPI_Comm_free(&sub_com) ;
  }

  inline bool firstCompare(const pair<int,int> &i1, const pair<int,int> &i2) {
    return i1.first < i2.first ;
  }

  void ORBPartition(const vector<vector3d<float> > &pnts,
                    vector<int> &procid,
                    MPI_Comm comm) {
    int sz = pnts.size() ;
    int p = 0;
    MPI_Comm_size(comm,&p) ;

    int r = 0 ;
    MPI_Comm_rank(comm,&r) ;
    
    if(p == 1) {
      vector<int> pid(sz) ;
      for(int i=0;i<sz;++i)
        pid[i] = 0 ;
      procid.swap(pid) ;
      return ;
    }

    vector<int> sizes(p) ;

    MPI_Allgather(&sz,1,MPI_INT,&sizes[0],1,MPI_INT,comm) ;

    vector<int> offsets(p+1) ;
    offsets[0] = 0 ;
    for(int i=0;i<p;++i)
      offsets[i+1] = offsets[i] + sizes[i] ;

    
    vector<vpair> info(sz) ;
    for(int i=0;i<sz;++i)
      info[i] = vpair(pnts[i],offsets[r]+i) ;

    ORBSort(info,comm) ;

    vector<pair<int,int> > proc_pairs(info.size()) ;
    for(size_t i=0;i<info.size();++i)
      proc_pairs[i] = pair<int,int>(info[i].second,r) ;
    vector<pair<int,int> > splitters(p-1) ;
    for(int i=0;i<p-1;++i) {
      splitters[i].first  = offsets[i+1] ;
      splitters[i].second = -1 ;
    }
    sort(proc_pairs.begin(),proc_pairs.end(),firstCompare) ;
    parSplitSort(proc_pairs,splitters,firstCompare,comm) ;
    
    vector<int> pid(sz) ;
    if(proc_pairs.size() != pid.size())
      cerr << "something wrong in matching proc_pairs " << endl ;
    for(int i=0;i<sz;++i)
      pid[i] = proc_pairs[i].second ;

    procid.swap(pid) ;
  }
}
