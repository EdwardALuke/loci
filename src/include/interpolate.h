//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
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
#ifndef INTERPOLATE_H
#define INTERPOLATE_H
#include <kd_tree.h>
#include <store.h>
#include <storeVec.h>
#include <Loci_types.h>
#include <LociGridReaders.h>
#include <vector>

#ifdef MEMORY_DEBUG
#define MEMORY_PROFILE(X) memSpace(string(# X) +";" + string(__FILE__) + ";"+i2s(__LINE__))
#else
#define MEMORY_PROFILE(X)
#endif

namespace Loci {

  std::vector<int> get_stencil(const Loci::kdTree::KDTree<float> &kd,
				 kdTree::coord3df ccenter,
				 double delta) ;
  std::vector<int> get_stencil(const Loci::kdTree::KDTree<double> &kd,
				 kdTree::coord3d ccenter,
				 double delta) ;

  /*  template<class T> inline std::vector<int>
  get_stencil(const Loci::kdTree::KDTree<float> &kd,
	      vector3d<T> pnt,
	      double deltai) {
    kdTree::coord3df ccenter ;
    ccenter[0] = realToFloat(pnt.x) ;
    ccenter[1] = realToFloat(pnt.y) ;
    ccenter[2] = realToFloat(pnt.z) ;
    double delta = deltai ;
    return get_stencil_i(kd,ccenter,delta) ;
  }
    
 
  template<class T> inline std::vector<int>
  get_stencil(const Loci::kdTree::KDTree<double> &kd,
	      vector3d<T> pnt,
	      double deltai) {
    kdTree::coord3d ccenter ;
    ccenter[0] = realToDouble(pnt.x) ;
    ccenter[1] = realToDouble(pnt.y) ;
    ccenter[2] = realToDouble(pnt.z) ;
    double delta = deltai ;
    return get_stencil_i(kd,ccenter,delta) ;
  }
  */
  
  void stencil_weights(std::vector<double> &w,
                       std::vector<int> &neighbors,
                       const store<vector3d<double> > &loc,
                       vector3d<double> ipnt) ;

  int collectPointsSizes(const Loci::kdTree::KDTree<float> &kd,
			 Loci::kdTree::KDTree<float>::bounds bnd) ;
  
  void getStencilBoundingBox2(Loci::kdTree::KDTree<float>::bounds &bnd,
			      double &delta,
			      const Loci::kdTree::KDTree<float> &kd,
			      const vector3d<double> pnts[],
                              int start, int end) ;
  void getStencilBoundingBox(Loci::kdTree::KDTree<float>::bounds &bnd,
                             double &delta,
                             const const_store<vector3d<double> > &pnts,
                             entitySet dom) ;
  void collectPoints(std::vector<Loci::kdTree::KDTree<float>::coord_info> &pout,
                     const Loci::kdTree::KDTree<float> &kd,
                     Loci::kdTree::KDTree<float>::bounds bnd) ;
  
  void getCommSchedFromStencil(std::vector<int> &send_info_out,
                               std::vector<int> &req_sizes_out,
                               std::vector<int> &snd_sizes_out,
                               std::vector<int> &access_out,
                               const std::vector<Array<int,4> > &stencil,
                               const store<int> &ids,
                               const std::vector<int> &distribution) ;
  void remapStencil(std::vector<Array<int,4> > &stencils,
                    const std::vector<int> &access,
                    const store<int> &ids) ;


  template <class T>
  void sendStencilData(store<T> &stencilData,
                       const_store<T> &sourceData,
                       const std::vector<int> &send_info,
                       const std::vector<int> &req_sizes_in,
                       const std::vector<int> &snd_sizes_in) {
    std::vector<T> databuf(send_info.size()) ;
    for(size_t i = 0;i<send_info.size();++i) {
      int id = send_info[i] ;
      databuf[i] = sourceData[id] ;
    }

    int p = Loci::MPI_processes ;
    std::vector<int> req_sizes(p),snd_sizes(p) ;
    const int data_size_b = sizeof(T) ;
    for(int i=0;i<p;++i) {
      req_sizes[i] = req_sizes_in[i]*data_size_b ;
      snd_sizes[i] = snd_sizes_in[i]*data_size_b ;
    }

    std::vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+req_sizes[i-1] ;

    std::vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+snd_sizes[i-1] ;
    }

    int loc_size = 0 ;
    for(int i=0;i<p;++i)
      loc_size += req_sizes_in[i] ;

    stencilData.allocate(entitySet(interval(0,loc_size-1))) ;

    MPI_Alltoallv(&databuf[0],&snd_sizes[0],&rdispls[0],MPI_BYTE,
                  &stencilData[0],&req_sizes[0],&sdispls[0],MPI_BYTE,
                  MPI_COMM_WORLD) ;
  }
  
  template <class T>
  void sendStencilData(storeVec<T> &stencilData,
                       const_storeVec<T> &sourceData,
                       const vector<int> &send_info,
                       const vector<int> &req_sizes_in,
                       const vector<int> &snd_sizes_in) {
    MEMORY_PROFILE(sendStencilDataStartv) ;

#ifdef VERBOSE
    entitySet dom = sourceData.domain() ;
#endif
    int vec_size = sourceData.vecSize() ;
    vector<T> databuf(send_info.size()*vec_size) ;
    for(size_t i = 0;i<send_info.size();++i) {
      int id = send_info[i] ;
#ifdef VERBOSE
      if(!dom.inSet(id)) {
        debugout << "id=" <<id << " out of domain " << dom << endl ;
        id = dom.Min() ;
      }

#endif
      for(int j=0;j<vec_size;++j) {
        databuf[i*vec_size+j] = sourceData[id][j] ;
      }
    }

    int p = MPI_processes ;
    vector<int> req_sizes(p),snd_sizes(p) ;
    for(int i=0;i<p;++i) {
      req_sizes[i] = req_sizes_in[i]*vec_size*sizeof(T) ;
      snd_sizes[i] = snd_sizes_in[i]*vec_size*sizeof(T) ;
    }

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+req_sizes[i-1] ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+snd_sizes[i-1] ;
    }

    int loc_size = 0 ;
    for(int i=0;i<p;++i)
      loc_size += req_sizes_in[i] ;

    stencilData.allocate(entitySet(interval(0,loc_size-1))) ;
    stencilData.setVecSize(vec_size) ;

    MEMORY_PROFILE(sendStencilDataStartall2all) ;
    
    MPI_Alltoallv(&databuf[0],&snd_sizes[0],&rdispls[0],MPI_BYTE,
                  &stencilData[0][0],&req_sizes[0],&sdispls[0],MPI_BYTE,
                  MPI_COMM_WORLD) ;
    MEMORY_PROFILE(sendStencilDataStartEnd3dv) ;
  }

  template<class T> void scatterVector(std::vector<T> &v,
                                  const std::vector<int> &procid,
                                  MPI_Comm comm) {
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;

    if(p==1)
      return ;
    // count sends and recvs
    std::vector<int> recv_count(p,0),send_count(p,0) ;
    
    for(size_t i=0;i<procid.size();++i)
      send_count[procid[i]]++ ;

    MPI_Alltoall(&send_count[0],1,MPI_INT,&recv_count[0],1,MPI_INT,comm) ;

    int final_size = 0 ;
    for(int i=0;i<p;++i)
      final_size += recv_count[i] ;

    std::vector<T> sendbuf(procid.size()) ;
    std::vector<T> recvbuf(final_size) ;
    std::vector<int> soffsets(p,0) ;
    for(int i=1;i<p;++i)
      soffsets[i] = send_count[i-1] + soffsets[i-1] ;
    for(size_t i=0;i<procid.size();++i) {
      sendbuf[soffsets[procid[i]]] = v[i] ;
      soffsets[procid[i]]++ ;
    }

    for(int i=0;i<p;++i) {
      send_count[i] *= sizeof(T) ;
      recv_count[i] *= sizeof(T) ;
    }

    std::vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+send_count[i-1] ;

    std::vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+recv_count[i-1] ;
    }

    MPI_Alltoallv(&sendbuf[0],&send_count[0],&sdispls[0],MPI_BYTE,
                  &recvbuf[0],&recv_count[0],&rdispls[0],MPI_BYTE,
                  comm) ;

    v.swap(recvbuf) ;
  }
    
  template<class T> void gatherVector(std::vector<T> &v,
                                 const std::vector<int> &procid,
                                 MPI_Comm comm) {
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    if(p==1)
      return ;
    
    // count sends and recvs
    std::vector<int> recv_count(p,0),send_count(p,0) ;
    
    for(size_t i=0;i<procid.size();++i)
      recv_count[procid[i]]++ ;

    MPI_Alltoall(&recv_count[0],1,MPI_INT,&send_count[0],1,MPI_INT,comm) ;

    std::vector<T> recvbuf(procid.size()) ;

    std::vector<int> roffsets(p,0) ;
    for(int i=1;i<p;++i)
      roffsets[i] = recv_count[i-1] + roffsets[i-1] ;

    for(int i=0;i<p;++i) {
      send_count[i] *= sizeof(T) ;
      recv_count[i] *= sizeof(T) ;
    }

    std::vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+send_count[i-1] ;

    std::vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+recv_count[i-1] ;
    }

    MPI_Alltoallv(&v[0],&send_count[0],&sdispls[0],MPI_BYTE,
                  &recvbuf[0],&recv_count[0],&rdispls[0],MPI_BYTE,
                  comm) ;

    // unscramble recvbuf
    // deallocate v
    { std::vector<T> tmp; v.swap(tmp) ;}
    // reallocate v
    { std::vector<T> tmp(recvbuf.size()) ; v.swap(tmp) ;} 

    for(size_t i=0;i<procid.size();++i) {
      v[i] = recvbuf[roffsets[procid[i]]] ;
      roffsets[procid[i]]++ ;
    }
  }
  template <class T> void allGatherVec(std::vector<T> &r, std::vector<T> &in,
                                       const MPI_Comm &comm) {
    MEMORY_PROFILE(allGatherVecBegin) ;
    int sz = in.size() ;
    int p ;
    MPI_Comm_size(comm,&p) ;
    std::vector<int> psz(p) ;
    MPI_Allgather(&sz,1,MPI_INT,&psz[0],1,MPI_INT,comm) ;
    int tot = psz[0] ;
    std::vector<int> displ(p) ;
    displ[0] = 0 ;
    psz[0] *= sizeof(T) ;
    for(int i=1;i<p;++i) {
      displ[i] = displ[i-1]+psz[i-1] ;
      tot += psz[i] ;
      psz[i] *= sizeof(T) ;
    }
    std::vector<T> buf(tot) ;
    MPI_Allgatherv(&in[0],sz*sizeof(T),MPI_BYTE,
                   &buf[0],&psz[0],&displ[0],MPI_BYTE,comm) ;
    MEMORY_PROFILE(allGatherVecEnd) ;
    r.swap(buf) ;
  }
  
  struct bound_info {
    Loci::kdTree::KDTree<float>::bounds bnd ;
    int start,stop ;
  } ;
  
  template<class T> 
  void communicateBoxPoints(std::vector<std::vector<T>  >
                            &loc_pbox, 
                            const std::vector<int> &box_sp,
                            const std::vector<int> &box_tp,
                            const std::vector<int> &b_sizes,
                            const std::vector<int> &send_block_id,
                            const std::vector<int> &recv_block_id,
                            std::vector<T> &pbox,
                            const std::vector<bound_info> &boxes,
                            const MPI_Comm &comm) {

    MEMORY_PROFILE(communicateBoxPointsBegin) ;
#ifdef REPORT_TIMES
    Loci::stopWatch s ;
    s.start() ;
#endif    
    int p ;
    MPI_Comm_size(comm,&p) ;
    int r ;
    MPI_Comm_rank(comm,&r) ;
    
    std::vector<std::vector<T>  > tmpx(recv_block_id.size()) ;
    loc_pbox.swap(tmpx) ;
                            
    std::vector<MPI_Request> req_queue(recv_block_id.size()) ;

    // post receives
    int coord_size = sizeof(T) ;
    for(size_t i=0;i<recv_block_id.size();++i) {
      int bk = recv_block_id[i] ;
#ifdef VERBOSE
      Loci::debugout <<"bsize_list["<<bk << "]= " << b_sizes[bk]  << endl ;
#endif
      std::vector<T> tmp(b_sizes[bk]) ;
      loc_pbox[i].swap(tmp) ;
      MPI_Irecv(&loc_pbox[i][0],b_sizes[bk]*coord_size,MPI_BYTE,box_sp[bk],bk,comm,&req_queue[i]) ;
    }

    // send data
    for(size_t i=0;i<send_block_id.size();++i) {
      int bk = send_block_id[i] ;
      int to = box_tp[bk] ;
      MPI_Send(&pbox[boxes[i].start],b_sizes[bk]*coord_size,MPI_BYTE,
               to,bk,comm) ;
    }

#ifdef VERBOSE
    Loci::debugout << "waiting " << recv_block_id.size() << endl ;
#endif
    // wait for comm to complete
    if(recv_block_id.size() != 0) {
      std::vector<MPI_Status> stat_queue(recv_block_id.size()) ;
      MPI_Waitall(recv_block_id.size(),&req_queue[0],&stat_queue[0]) ;
    }

    MEMORY_PROFILE(communicateBoxPointsEnd) ;
#ifdef REPORT_TIMES
    Loci::debugout << "time to collect boxes: " << s.stop() << endl ;
    s.start() ;
#endif    
  }                            
  template <class T>
  void returnBlockData(std::vector<T> &recv_data,
                       std::vector<std::vector<T> > &data_to_send,
                       int group_size,
                       const std::vector<bound_info> &boxes,
                       const std::vector<int> &box_sp,
                       const std::vector<int> &box_tp,
                       const std::vector<int> &b_sizes,
                       const std::vector<int> &send_block_id,
                       const std::vector<int> &recv_block_id,
                       const MPI_Comm &comm) {
    //    MEMORY_PROFILE(returnBlockDataBegin) ;
    std::vector<MPI_Request> req_queue(send_block_id.size()) ;
    // allocate recv data
    int rcv_sz = 0 ;
    for(size_t i=0;i<send_block_id.size();++i) {
      int bk = send_block_id[i] ;
      rcv_sz += b_sizes[bk]*group_size ;
    }
      
    {
      std::vector<T> rcv_data_x(rcv_sz) ;
      recv_data.swap(rcv_data_x) ;
    }

    // Send stencil data
    for(size_t i=0;i<send_block_id.size();++i) {
      int bk = send_block_id[i] ;
      int to = box_tp[bk] ;
      MPI_Irecv(&recv_data[boxes[i].start*group_size],b_sizes[bk]*sizeof(T)*group_size,
                MPI_BYTE,to,bk,comm,&req_queue[i]) ;
    }
    for(size_t i=0;i<recv_block_id.size();++i) {
      int bk = recv_block_id[i] ;
      int to = box_sp[bk] ;
      MPI_Send(&data_to_send[i][0],data_to_send[i].size()*sizeof(T),
               MPI_BYTE,to,bk,comm) ;
    }

    if(send_block_id.size() != 0) {
      std::vector<MPI_Status> stat_queue(send_block_id.size()) ;
      MPI_Waitall(send_block_id.size(),&req_queue[0],&stat_queue[0]) ;
    }
    //    MEMORY_PROFILE(returnBlockDataEnd) ;
  }

  // Group points into a collection of bounding boxes (no more than 2^levels)
  // For each bounding box, find the smallest reciprocal point spacing
  // and assign this spacing to the centroid of the bounding box.
  void collectGroups(std::vector<Loci::vpair> &results,
                     std::vector<Loci::vpair> &inputs, int start, int end,
                     int levels, int min_pnts) ;
  
  void getBoundingBoxDecomp(std::vector<kdTree::KDTree<float>::coord_info> &pnts, 
                            std::vector<bound_info> &boxes, 
                            double split, int dim, int levels,
                            int start, int stop, double delta_lim,
                            int split_lim) ;
  void getBoundingBoxes(std::vector<kdTree::KDTree<float>::coord_info> &pnts, 
                        std::vector<bound_info> &boxes,
                        int levels, double delta_lim,
                        int split_lim) ;
  
  void boundingBoxDistribution(std::vector<int> &box_sp,
                               std::vector<int> &box_tp,
                               std::vector<int> &b_sizes,
                               std::vector<int> &send_block_id,
                               std::vector<int> &recv_block_id,
                               const std::vector<bound_info> &boxes,
                               const MPI_Comm &comm) ;
  
  void recieveTargetPoints(std::vector<std::vector<kdTree::KDTree<float>::coord_info>  > &recv_targets,
                           const Loci::kdTree::KDTree<float> *kd,
                           const std::vector<bound_info> &boxes_g,
                           const std::vector<int> &box_sp,
                           const std::vector<int> &box_tp,
                           const std::vector<int> &b_sizes,
                           const std::vector<int> &send_block_id,
                           const std::vector<int> &recv_block_id,
                           const MPI_Comm &comm) ;
}

#endif
