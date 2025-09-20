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
#include "multiStoreIO.h"

namespace Loci {
 
  
  void findFileDist(int &fmin, int &fstart, int &fsz, int &delta,
		    entitySet read_set, fact_db &facts) {
    // Now create data structure distributed over read set to
    // contain mapping information from file to local
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    read_set = read_set & df->my_entities ;
    const_Map l2f ;
    l2f = df->l2f.Rep() ;
    // Find bounds of file number
    int fmin_local = l2f[read_set.Min()] ;
    int fmax_local = fmin_local ;
    FORALL(read_set,ii) {
      int fid = l2f[ii] ;
      fmin_local = min(fmin_local,fid) ;
      fmax_local = max(fmax_local,fid) ;
    } ENDFORALL ;
    int fmin_global = fmin_local ;
    MPI_Allreduce(&fmin_local,&fmin_global,1,MPI_INT,
		  MPI_MIN,MPI_COMM_WORLD) ;
    fmin = fmin_global ;
    int fmax_global = fmax_local ;
    MPI_Allreduce(&fmax_local,&fmax_global,1,MPI_INT,
		  MPI_MAX,MPI_COMM_WORLD) ;

    // Find distribution of file numbers to processors
    int nfilenums = fmax_global-fmin_global+1 ;
    const int p = MPI_processes ;
    delta = (nfilenums+p-1)/MPI_processes ;
    const int r = MPI_rank ;
    fstart = delta*r ; // start of file number on this processor
    fsz = delta ; // size of file numbers on this processor
    if(delta*(r+1) > nfilenums)
      fsz = max(nfilenums-delta*r,0) ;
  }    

  void gatherFileMapping(std::vector<int> &map_local,
			 std::vector<int> &map_proc,
			 int fmin, int fstart, int fsz, int delta,
			 entitySet read_set,
			 fact_db &facts) {
    // Now create data structure distributed over read set to
    // contain mapping information from file to local
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    read_set = read_set & df->my_entities ;
    const_Map l2f ;
    l2f = df->l2f.Rep() ;

    const int p = MPI_processes ;
    // Collect information about mapping between file number and
    // global-processor numbering
    // Find file number to processor mapping
    std::vector<int> sendto(p,0) ;
    FORALL(read_set,ii) {
      int fid = l2f[ii]-fmin ;
      int ps = fid/delta ;
      sendto[ps]++ ;
    } ENDFORALL ;
    std::vector<int> recvfrom(p,0) ;
    MPI_Alltoall(&sendto[0],1,MPI_INT,&recvfrom[0],1,MPI_INT,
		 MPI_COMM_WORLD) ;
    const int r = MPI_rank ;

    std::vector<int> sendbuf(read_set.size()*3) ;
    std::vector<int> offsets(p,0) ;
    std::vector<int> send_offsets(p+1,0) ;
    std::vector<int> recv_offsets(p+1,0) ;
    for(int i=0;i<p;++i) {
      send_offsets[i+1] = send_offsets[i]+sendto[i] ;
      recv_offsets[i+1] = recv_offsets[i]+recvfrom[i] ;
    }
    FORALL(read_set,ii) {
      int fid = l2f[ii]-fmin ;
      int p = fid/delta ;
      int offset = offsets[p]+send_offsets[p] ; ;
      offsets[p]++ ;
      sendbuf[offset*3+0] = fid ; // File Id
      sendbuf[offset*3+1] = ii ; // Local Id
      sendbuf[offset*3+2] = r ; // local processor
    } ENDFORALL ;
    int rsz = 0 ;
    for(int i=0;i<p;++i)
      rsz += recvfrom[i] ;
    std::vector<int> recvbuf(rsz*3) ;
    int nreq = 0 ;
    for(int i=0;i<p;++i) {
      if(sendto[i] > 0)
	nreq++ ;
      if(recvfrom[i] > 0)
	nreq++ ;
    }
    std::vector<MPI_Request> recv_Requests(nreq) ;
    int req = 0 ;
    for(int i=0;i<p;++i)
      if(recvfrom[i] > 0) {
	MPI_Irecv(&recvbuf[recv_offsets[i]*3],recvfrom[i]*3,MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    for(int i=0;i<p;++i)
      if(sendto[i] > 0) {
	MPI_Isend(&sendbuf[send_offsets[i]*3],sendto[i]*3,MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    std::vector<MPI_Status> statuslist(nreq) ;
    MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;

    for(int i=0;i<rsz;++i) {
      int lf =recvbuf[i*3+0] - fstart ; // Find file number to local filenumber
      map_local[lf] = recvbuf[i*3+1] ;
      map_proc[lf] = recvbuf[i*3+2] ;
    }
    
  }
  
  // Find mapping from fileID space to local number/procid
  void findMapping(std::vector<int> &local_num,
		   std::vector<int> &procID,
		   entitySet read_set,
		   const std::vector<int> &fileID,
		   fact_db &facts) {

    // Find lowest number in fileID list since we are assuming that
    // the read_sets match
    int sz = fileID.size() ;
    int local_min = fileID[0] ;
    for(int i=0;i<sz;++i) {
      local_min = min(fileID[i],local_min) ;
    }

    // Handle single processor case
    const int p = MPI_processes ;
    if(p == 1) {
      int imin = read_set.Min() ;
      for(int i=0;i<sz;++i) {
	local_num[i] = fileID[i]-local_min+imin ;
	procID[i] = 0 ;
      }
      return ;
    }

    // Find minimum id from fileID input
    int global_min = local_min ;
    MPI_Allreduce(&local_min,&global_min,1,MPI_INT,
		  MPI_MIN,MPI_COMM_WORLD) ;

    int fmin=0, fstart=0, fsz=0,delta=0 ;
    // fmin is minimum file number in read_set
    // fstart is start of file number partition on this processor
    // fsz is the size of the partition on this processor
    // delta is the delta for mapping between file number and processor
    findFileDist(fmin, fstart, fsz, delta, read_set, facts) ;
    
    std::vector<int> map_local(fsz) ;
    std::vector<int> map_proc(fsz,-1) ;

    // Gather local number and processor to file distribution
    gatherFileMapping(map_local,map_proc,fmin, fstart, fsz, delta,
		      read_set,facts) ;

    // Now go through fileID and send requests to each processor to fill in
    // mapping information from fileID domain.
    std::vector<int> send_req(p,0) ;
    for(int i=0;i<sz;++i) {
      int fid = fileID[i]-global_min ;
      int rp = fid/delta ;
      send_req[rp]++ ;
    }
    std::vector<int> recv_req(p,0) ;
    MPI_Alltoall(&send_req[0],1,MPI_INT,&recv_req[0],1,MPI_INT,
		 MPI_COMM_WORLD) ;
    std::vector<int> sendbuffer(fileID.size()) ;
    std::vector<int> soffsets(p+1,0) ;
    std::vector<int> roffsets(p+1,0) ;
    int rsz = 0 ;
    for(int i=0;i<p;++i) {
      rsz += recv_req[i] ;
      soffsets[i+1] = soffsets[i]+send_req[i] ;
      roffsets[i+1] = roffsets[i]+recv_req[i] ;
    }
    std::vector<int> scounts(p,0) ;
    for(int i=0;i<sz;++i) {
      int fid = fileID[i]-global_min ;
      int rp = fid/delta ;
      sendbuffer[soffsets[rp]+scounts[rp]] = fid ;
      scounts[rp]++ ;
    }
    std::vector<int> recvbuffer(rsz) ;
    
    int nreq = 0 ;
    for(int i=0;i<p;++i) {
      if(send_req[i] > 0)
	nreq++ ;
      if(recv_req[i] > 0)
	nreq++ ;
    }
    std::vector<MPI_Request> recv_Requests(nreq) ;
    int req = 0 ;
    for(int i=0;i<p;++i)
      if(recv_req[i] > 0) {
	MPI_Irecv(&recvbuffer[roffsets[i]],recv_req[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    for(int i=0;i<p;++i)
      if(send_req[i] > 0) {
	MPI_Isend(&sendbuffer[soffsets[i]],send_req[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    std::vector<MPI_Status> statuslist(nreq) ;
    MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;
    // Now we go the requests, get information to return
    std::vector<int> datasend(rsz*2) ;
    int map_error = false ;
    for(int i=0;i<rsz;++i) {
      int fid = recvbuffer[i] ;
      int flocal = fid-fstart ;
      if(flocal >= fsz) {
	map_error = true ;
	datasend[i*2+0] = -1 ;
	datasend[i*2+1] = -1 ;
      } else {
	datasend[i*2+0] = map_local[flocal] ;
	datasend[i*2+1] = map_proc[flocal] ;
      }
    }
    if(map_error)
      cerr << "error in mapping, file number out of bounds" << endl ;

    std::vector<int> datarecv(fileID.size()*2) ;
    
    req = 0 ;
    for(int i=0;i<p;++i)
      if(send_req[i] > 0) {
	MPI_Irecv(&datarecv[soffsets[i]*2],send_req[i]*2,MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    for(int i=0;i<p;++i)
      if(recv_req[i] > 0) {
	MPI_Isend(&datasend[roffsets[i]*2],recv_req[i]*2,MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;
    bool error = false ;

    std::vector<int> sc(p,0) ;
    for(int i=0;i<sz;++i) {
      int fid = fileID[i]-global_min ;
      int rp = fid/delta ;
      local_num[i] = datarecv[(soffsets[rp]+sc[rp])*2+0] ;
      procID[i] = datarecv[(soffsets[rp]+sc[rp])*2+1] ;
      sc[rp]++ ;
      if(procID[i] < 0) {
	error = true ;
      }
    }

    if(error)
      cerr << "unable to find processor mapping" << endl ;      
  }

  void distributeMapMultiStore(std::vector<int> &send_sz,
			       std::vector<int> &recv_sz,
			       std::vector<int> &recv_local_num,
			       const std::vector<int> &local_num,
			       const std::vector<int> &procID) {
    const int p = MPI_processes ;
    for(int i=0;i<p;++i)
      send_sz[i] = 0 ;
    for(size_t i=0;i<procID.size();++i)
      send_sz[procID[i]]++ ;
    MPI_Alltoall(&send_sz[0],1,MPI_INT, &recv_sz[0],1,MPI_INT,
		 MPI_COMM_WORLD) ;
    std::vector<int> send_data(procID.size()) ;
    std::vector<int> soffsets(p+1,0) ;
    for(int i=0;i<p;++i)
      soffsets[i+1] = soffsets[i]+send_sz[i] ;
    std::vector<int> scounts(p,0) ;
    for(size_t i=0;i<procID.size();++i) {
      int pc = procID[i] ;
      send_data[soffsets[pc]+scounts[pc]] = local_num[i] ;
      scounts[pc]++ ;
    }
    std::vector<int> roffsets(p+1,0) ;
    for(int i=0;i<p;++i)
      roffsets[i+1] = roffsets[i]+recv_sz[i] ;
    std::vector<int> recv_data(roffsets[p]) ;
    int nreq = 0 ;
    for(int i=0;i<p;++i) {
      if(send_sz[i] > 0)
	nreq++ ;
      if(recv_sz[i] > 0)
	nreq++ ;
    }
    std::vector<MPI_Request> recv_Requests(nreq) ;
    int req = 0 ;
    for(int i=0;i<p;++i)
      if(recv_sz[i] > 0) {
	MPI_Irecv(&recv_data[roffsets[i]],recv_sz[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    for(int i=0;i<p;++i)
      if(send_sz[i] > 0) {
	MPI_Isend(&send_data[soffsets[i]],send_sz[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    std::vector<MPI_Status> statuslist(nreq) ;
    MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;
    recv_local_num.swap(recv_data) ;
  }    

  void sendCounts(std::vector<int>&recv_count,
		  const std::vector<int> &send_sz,
		  const std::vector<int> &recv_sz,
		  const std::vector<int> &recv_local_num,
		  const std::vector<int> &counts,
		  const std::vector<int> &procID) {
    const int p = MPI_processes ;
    std::vector<int> soffsets(p+1,0) ;
    std::vector<int> roffsets(p+1,0) ;
    for(int i=0;i<p;++i) {
      soffsets[i+1] = soffsets[i] + send_sz[i] ;
      roffsets[i+1] = roffsets[i] + recv_sz[i] ;
    }
    std::vector<int> so(p,0) ;
    std::vector<int> sbuf(counts.size()) ;
    for(size_t i=0;i<counts.size();++i) {
      int ps = procID[i] ;
      sbuf[soffsets[ps]+so[ps]] = counts[i] ;
      so[ps]++ ;
    }

    std::vector<int> rbuf(roffsets[p]) ;
    int nreq = 0 ;
    for(int i=0;i<p;++i) {
      if(send_sz[i] > 0)
	nreq++ ;
      if(recv_sz[i] > 0)
	nreq++ ;
    }
    std::vector<MPI_Request> recv_Requests(nreq) ;
    int req = 0 ;
    for(int i=0;i<p;++i)
      if(recv_sz[i] > 0) {
	MPI_Irecv(&rbuf[roffsets[i]],recv_sz[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    for(int i=0;i<p;++i)
      if(send_sz[i] > 0) {
	MPI_Isend(&sbuf[soffsets[i]],send_sz[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    std::vector<MPI_Status> statuslist(nreq) ;
    MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;
    recv_count.swap(rbuf) ;
    
  }
}
