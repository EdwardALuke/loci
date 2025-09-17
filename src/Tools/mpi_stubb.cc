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

#ifdef MPI_STUBB
#include "mpi.h"

#include <string.h>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
using std::cerr ;
using std::endl ;
const int MPI_TYPE_SIZE[] =
  {
    0,
    sizeof(char), sizeof(short), sizeof(int), sizeof(long),
    sizeof(unsigned char), sizeof(unsigned short), sizeof(unsigned int),
    sizeof(unsigned long), sizeof(float), sizeof(double),
    sizeof(long double), sizeof(long long),	sizeof(int),
    sizeof(float), sizeof(double), 2*sizeof(float), 2*sizeof(double),
    sizeof(bool), sizeof(char), 1, 2, 4, 8, 4, 8, 16, 1, 1, 1, 1,
    sizeof(float)+sizeof(int), sizeof(double)+sizeof(int),
    sizeof(long)+sizeof(int), 2*sizeof(int)
  } ;

extern "C" {
  
  int MPI_GET_TYPE_SIZE(MPI_Datatype type)
  {
    
    if(type >= sizeof(MPI_TYPE_SIZE)/sizeof(int))
      return 1 ;
    else
      return MPI_TYPE_SIZE[type] ;
  }
}

/** Added by Kenny Moser krm104 **/
struct Comm_Buffer
{
  int size;
  int ab_size;

  void *temp_buff;

  void *attached_buff;

  bool buffered;

  Comm_Buffer()
  {
    size = 0;
    ab_size = 0;
    temp_buff = NULL;
    attached_buff = NULL;
    buffered = false;
  }

  Comm_Buffer(int n_Size, MPI_Datatype datatype)
  {
    size = n_Size;
    ab_size = 0;
    temp_buff = malloc(MPI_GET_TYPE_SIZE(datatype) * n_Size);
    buffered = false;
    attached_buff = NULL;
  }

  int copy(const void *buf, int n_Size, MPI_Datatype datatype)
  {
    size = n_Size;
    if(temp_buff != NULL) {
      free(temp_buff) ;
    }
    temp_buff = malloc(MPI_GET_TYPE_SIZE(datatype) * n_Size);
    memcpy(temp_buff, buf, MPI_GET_TYPE_SIZE(datatype) * n_Size);
    return 1;
  }

  int copy_buff(const void *buf, int n_Size, MPI_Datatype datatype)
  {
    if (buffered)
      {
	if (ab_size >= (MPI_GET_TYPE_SIZE(datatype) * n_Size))
	  {
	    memcpy(attached_buff, buf, MPI_GET_TYPE_SIZE(datatype) * n_Size);
	    return 1;
	  }
	else
	  {
	    return 0;
	  }
      }
    else
      {
	size = n_Size;
	if(temp_buff != NULL) {
	  free(temp_buff) ;
	}

	temp_buff = malloc(MPI_GET_TYPE_SIZE(datatype) * n_Size);

	memcpy(temp_buff, buf, MPI_GET_TYPE_SIZE(datatype) * n_Size);
	return 1;
      }
  }

  int attach_buff(void *buff, int n_size)
  {
    buffered = true;

    ab_size = n_size;
		 
    attached_buff = buff;

    return 1;
  }

  int dettach_buff(void *buff)
  {
    if (buffered && buff == attached_buff)
      {
	size = 0;
	attached_buff = NULL;
	buffered = false;
	return 1;
      }
    else
      return 0;
  }

};

extern "C" {
  

  /*************************************/
  /* MPI-1 bindings, sorted by chapter */
  /*************************************/


  /* 3.2 */

  void err_report()
  {
    std::cerr << "MPI Function stub!" << std::endl ;

#ifndef NO_CSTDLIB
    std::abort() ;
#else
    abort() ;
#endif
  }

  //Added by Kenny Moser krm104
  Comm_Buffer Temp_Comm_Buffer;

  int  MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
  {
    Temp_Comm_Buffer.copy(buf, count, datatype);
    return MPI_SUCCESS;
  }

  int  MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, 
		MPI_Comm comm, MPI_Status *status) 
  {
    memcpy(buf, Temp_Comm_Buffer.temp_buff, MPI_GET_TYPE_SIZE(datatype) * count);
    status->MPI_SOURCE = 0;
    status->MPI_TAG = 0;
    status->size = count;
    return MPI_SUCCESS;
  }

  int  MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count)
  {
    if (MPI_GET_TYPE_SIZE(datatype) == 0)
      {
	*count = 0 ;
	return MPI_ERR_TYPE;
      }
    else
      {
	*count = status->size;
	return MPI_SUCCESS;
      }	 
  }

  /* 3.4 */

  int  MPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
  {
    Temp_Comm_Buffer.copy_buff(buf, count, datatype);
    return MPI_SUCCESS;
    //err_report(); return -1;
  }

  int  MPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
  {
    Temp_Comm_Buffer.copy(buf, count, datatype);
    return MPI_SUCCESS;
    //err_report(); return -1;
  }

  int  MPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
  {
    Temp_Comm_Buffer.copy(buf, count, datatype);
    return MPI_SUCCESS;
    //err_report(); return -1;
  }

  /* 3.6 */

  int  MPI_Buffer_attach(void *buf, int count)
  {
    Temp_Comm_Buffer.attach_buff(buf, count);
    return MPI_SUCCESS;
    //err_report(); return -1;
  }

  int  MPI_Buffer_detach(void *buf, int *count)
  {
    //err_report(); return -1;
    if (Temp_Comm_Buffer.dettach_buff(buf))
      return MPI_SUCCESS;
    else
      return MPI_ERR_BUFFER;
  }

  /* 3.7 */

  int  MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
  {
    Temp_Comm_Buffer.copy(buf, count, datatype);
    return MPI_SUCCESS;
    //err_report(); return -1;
  }

  int  MPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
  {
    Temp_Comm_Buffer.copy_buff(buf, count, datatype);
    return MPI_SUCCESS;
    //err_report(); return -1;
  }

  int  MPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
  {
    Temp_Comm_Buffer.copy(buf, count, datatype);
    return MPI_SUCCESS;
    //err_report(); return -1;
  }

  int  MPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
  {
    Temp_Comm_Buffer.copy(buf, count, datatype);
    return MPI_SUCCESS;
    //err_report(); return -1;
  }

  int  MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
  {
    memcpy(buf, Temp_Comm_Buffer.temp_buff, MPI_GET_TYPE_SIZE(datatype) * count);
    //status->MPI_SOURCE = 0;
    //status->MPI_TAG = 0;
    //tatus->size = count;
    return MPI_SUCCESS;
    //err_report(); return -1;
  }

  int  MPI_Wait(MPI_Request *, MPI_Status *)
  {
    cerr << "MPI_Wait" << endl ;
    err_report(); return -1;
    //return MPI_SUCCESS;
  }

  int  MPI_Test(MPI_Request *, int *, MPI_Status *)
  {
    cerr << "MPI_Test" << endl ;
    err_report(); return -1;
    //return MPI_SUCCESS;
  }

  int  MPI_Request_free(MPI_Request *)
  {
    cerr << "MPI_Request_free" << endl ;
    err_report(); return -1;
  }

  int  MPI_Waitany(int, MPI_Request *, int *, MPI_Status *)
  {
    cerr << "MPI_Waitany" << endl ;
    err_report(); return -1;
  }

  int  MPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *)
  {
    cerr << "MPI_Testany" << endl ;
    err_report(); return -1;
  }

  int  MPI_Waitall(int, MPI_Request *, MPI_Status *)
  {
    //    cerr << "MPI_Waitall" << endl ;
    //    err_report(); return -1;
    return 0 ;
  }

  int  MPI_Testall(int, MPI_Request *, int *, MPI_Status *)
  {
    cerr << "MPI_Testall" << endl ;
    err_report(); return -1;
  }

  int  MPI_Waitsome(int, MPI_Request *, int *, int *, MPI_Status *)
  {
    cerr << "MPI_Waitsome" << endl ;
    err_report(); return -1;
  }

  int  MPI_Testsome(int, MPI_Request *, int *, int *, MPI_Status *)
  {
    cerr << "MPI_Testsome" << endl ;
    err_report(); return -1;
  }

  /* 3.8 */

  int  MPI_Iprobe(int, int, MPI_Comm, int *, MPI_Status *)
  {
    cerr << "MPI_Iprobe" << endl ;
    err_report(); return -1;
  }

  int  MPI_Probe(int, int, MPI_Comm, MPI_Status *)
  {
    cerr << "MPI_Probe" << endl ;
    err_report(); return -1;
  }

  int  MPI_Cancel(MPI_Request *)
  {
    cerr << "MPI_Cancel" << endl ;
    err_report(); return -1;
  }

  int  MPI_Test_cancelled(const MPI_Status *, int *)
  {
    cerr << "MPI_Test_cancelled" << endl ;
    err_report(); return -1;
  }

  /* 3.9 */

  int  MPI_Send_init(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *)
  {
    cerr << "MPI_Send_init" << endl ;
    err_report(); return -1;
  }

  int  MPI_Bsend_init(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *)
  {
    cerr << "MPI_Bsend_init" << endl ;
  
    err_report(); return -1;
  }

  int  MPI_Ssend_init(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *)
  {
    cerr << "MPI_Ssend_init" << endl ;
    err_report(); return -1;
  }

  int  MPI_Rsend_init(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *)
  {
    cerr << "MPI_Rsend_init" << endl ;
    err_report(); return -1;
  }

  int  MPI_Recv_init(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *)
  {
    cerr << "MPI_Recv_init" << endl ;
    err_report(); return -1;
  }

  int  MPI_Start(MPI_Request *)
  {
    cerr << "MPI_Start" << endl ;
    err_report(); return -1;
  }

  int  MPI_Startall(int, MPI_Request *)
  {
    cerr << "MPI_Startall" << endl ;
    err_report(); return -1;
  }

  /* 3.10 */

  int  MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag,
		    void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag,
		    MPI_Comm comm, MPI_Status *status)
  {
    if (recvcount <= sendcount)
      {
	memcpy(recvbuf, sendbuf, MPI_GET_TYPE_SIZE(sendtype) * sendcount);
	return MPI_SUCCESS;
      }
    return MPI_ERR_COUNT;
  }

  int  MPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest, int sendtag,
			    int source, int recvtag, MPI_Comm comm, MPI_Status *status)
  {
    return MPI_SUCCESS;
  }

  /* 3.12 */

  int  MPI_Type_contiguous(int, MPI_Datatype, MPI_Datatype *)
  {
    cerr << "MPI_Type_contiguous" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_vector(int, int, int, MPI_Datatype, MPI_Datatype *)
  {
    cerr << "MPI_Type_vector" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *)
  {
    cerr << "MPI_Type_hvector" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_indexed(int, const int *,const int *, MPI_Datatype, MPI_Datatype *)
  {
    cerr << "MPI_Type_indexed" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_hindexed(int, int *, MPI_Aint *, MPI_Datatype, MPI_Datatype *)
  {
    cerr << "MPI_Type_hindexed" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *)
  {
    cerr << "MPI_Type_struct" << endl ;
    err_report(); return -1;
  }

  int  MPI_Address(void *location, MPI_Aint *address)
  {
    cerr << "MPI_Address" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_extent(MPI_Datatype, MPI_Aint *)
  {
    cerr << "MPI_Type_extent" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_size(MPI_Datatype, int *)
  {
    cerr << "MPI_Type_size" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_lb(MPI_Datatype, MPI_Aint *)
  {
    cerr << "MPI_Type_lb" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_ub(MPI_Datatype, MPI_Aint *)
  {
    cerr << "MPI_Type_ub" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_commit(MPI_Datatype *)
  {
    cerr << "MPI_Type_commit" << endl ;
    err_report(); return -1;
  }

  int  MPI_Type_free(MPI_Datatype *)
  {
    cerr << "MPI_Type_free" << endl ;
    err_report(); return -1;
  }

  int MPI_Type_create_resized(MPI_Datatype oldtype,
			      MPI_Aint lb,
			      MPI_Aint extent,
			      MPI_Datatype *newtype) {
    cerr << "MPI_Type_create_resized" << endl ;
    err_report(); return -1;
  }    
  
  int  MPI_Get_elements(const MPI_Status *, MPI_Datatype, int *)
  {
    cerr << "MPI_Get_elements" << endl ;
    err_report(); return -1;
  }

  /* 3.13 */

  int  MPI_Pack(const void *buf, int count, MPI_Datatype datatype,
		void *packbuf, int packsize, int *packpos,
		MPI_Comm comm)
  {
    int size = MPI_GET_TYPE_SIZE(datatype) ;
    size *= count ;
    memcpy(((char *)packbuf)+*packpos,buf,size) ;
    *packpos += size;
    return 0 ;
  }

  int  MPI_Unpack(void *packbuf, int packsize, int *packpos,
		  void *buf, int count, MPI_Datatype datatype,
		  MPI_Comm comm)
  {
    int size = MPI_GET_TYPE_SIZE(datatype) ;
    size *= count ;

    memcpy(buf,((char *)packbuf)+*packpos,size) ;
    *packpos += size;
    return 0 ;
  }

  int  MPI_Pack_size(int incount, MPI_Datatype datatype , MPI_Comm comm, int *size)
  {
    *size = MPI_GET_TYPE_SIZE(datatype)*incount ;
    return 0;
  }

  /* 4.3 */

  int  MPI_Barrier(MPI_Comm)
  {
    return MPI_SUCCESS;
  }

  /* 4.4 */

  int  MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
  {
    if (root == 0)
      return MPI_SUCCESS;
    else
      return MPI_ERR_ROOT;
  }

  /* 4.5 */

  int  MPI_Gather(const void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		  void *recvbuf, int recvcnt, MPI_Datatype recvtype,
		  int root , MPI_Comm comm)
  {
    if (root == 0)
      {
	memcpy(recvbuf, sendbuf, MPI_GET_TYPE_SIZE(recvtype)*recvcnt) ;
	return MPI_SUCCESS ;
      }
    else
      return MPI_ERR_ROOT;
  }


  int  MPI_Gatherv(const void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf,
		   const int *recvcnt, const int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm)
  {
    if (root == 0)
      {
	const int sz = MPI_GET_TYPE_SIZE(recvtype) ;
	memcpy((char *)recvbuf + sz*displs[0], sendbuf, sz * recvcnt[0]);
	return MPI_SUCCESS;
      }
    else 
      return MPI_ERR_ROOT;
    //err_report(); return -1;
  }

  /* 4.6 */

  int  MPI_Scatter(const void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf,
		   int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm)
  {
    if (root == 0)
      {
	memcpy(recvbuf, sendbuf, MPI_GET_TYPE_SIZE(recvtype)*recvcnt);
	return MPI_SUCCESS;
      }
    else
      return MPI_ERR_ROOT;
    //err_report(); return -1;
  }

  int  MPI_Scatterv(const void *sendbuf, const int *sendcnts, const  int *displs, MPI_Datatype sendtype,
		    void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm)
  {
    if (root == 0)
      {
	const int sz = MPI_GET_TYPE_SIZE(recvtype) ;
	memcpy(recvbuf, (char *)sendbuf + sz*displs[0], sz*recvcnt);
	return MPI_SUCCESS;
      }
    else
      return MPI_ERR_ROOT;
    //err_report(); return -1;
  }

  /* 4.7 */

  int  MPI_Allgather(const void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		     void *recvbuf, int recvcnt, MPI_Datatype recvtype,
		     MPI_Comm comm )
  {
    memcpy(recvbuf, sendbuf, MPI_GET_TYPE_SIZE(recvtype)*recvcnt) ;
    return MPI_SUCCESS ; 
  }

  int  MPI_Allgatherv(const void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		      void *recvbuf, const int *recvcounts, const int *displs,
		      MPI_Datatype recvtype, MPI_Comm comm)
  { 
    const int sz = MPI_GET_TYPE_SIZE(recvtype) ;
    memcpy(((char *)recvbuf)+displs[0]*sz, sendbuf, sz*recvcounts[0]);
    return MPI_SUCCESS ; 
  } 

  /* 4.8 */

  int  MPI_Alltoall(const void *send, int ssz, MPI_Datatype dts, void *recv, int rsz,
		    MPI_Datatype dtr, MPI_Comm comm)
  {
    memcpy(recv, send, MPI_GET_TYPE_SIZE(dts)*rsz) ;
    return MPI_SUCCESS;
  }

  int  MPI_Alltoallv(const void *sendbuf, const int *sendcnts, const int *sdispls,
		     MPI_Datatype sendtype,
		     void *recvbuf, const int *recvcnts, const int *rdispls,
		     MPI_Datatype recvtype, MPI_Comm comm)
  { 
    const int sz = MPI_GET_TYPE_SIZE(recvtype) ;
    memcpy(((char *)recvbuf)+rdispls[0]*sz, ((char *)sendbuf)+sdispls[0]*sz,
	   recvcnts[0]*sz) ;
    return MPI_SUCCESS ;
  }

  /* 4.9 */

  int  MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
		  MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm)
  {
    if (root == 0)
      {
	memcpy(recvbuf, sendbuf, MPI_GET_TYPE_SIZE(type)*count) ;
	return MPI_SUCCESS;
      }
    else
      return MPI_ERR_ROOT;
  } 

  int  MPI_Op_create(MPI_User_function *, int, MPI_Op *)
  {
    return MPI_SUCCESS;
  }

  int  MPI_Op_free(MPI_Op *)
  {
    return MPI_SUCCESS;
  }

  int  MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
		     MPI_Datatype type, MPI_Op op , MPI_Comm comm)
  { 
    memcpy(recvbuf, sendbuf, MPI_GET_TYPE_SIZE(type)*count) ;
    return MPI_SUCCESS;
  } 

  /* 4.10 */

  int  MPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int *recvcnts,
			  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
  {
    memcpy(recvbuf, sendbuf, MPI_GET_TYPE_SIZE(datatype)*recvcnts[0]) ;
    return MPI_SUCCESS;
  }

  /* 4.11 */

  int  MPI_Scan(const void *sendbuf, void *recvbuf, int recvcnts, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
  {
    memcpy(recvbuf, sendbuf, MPI_GET_TYPE_SIZE(datatype)*recvcnts) ;
    return MPI_SUCCESS;
  }

  /* 5.3 */

  int  MPI_Group_size(MPI_Group, int *)
  {

    cerr << "MPI_Group_size" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_rank(MPI_Group, int *)
  {
    cerr << "MPI_Group_rank" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_translate_ranks(MPI_Group, int, int *, MPI_Group, int *)
  {
    cerr << "MPI_Group_translate_ranks" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_compare(MPI_Group, MPI_Group, int *)
  {
    cerr << "MPI_Group_compare" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Comm_group(MPI_Comm, MPI_Group *)
  {
    cerr << "MPI_Comm_group" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_union(MPI_Group, MPI_Group, MPI_Group *)
  {
    cerr << "MPI_Group_union" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_intersection(MPI_Group, MPI_Group, MPI_Group *)
  {
    cerr << "MPI_Group_intersection" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_difference(MPI_Group, MPI_Group, MPI_Group *)
  {
    cerr << "MPI_Group_difference" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_incl(MPI_Group, int, int *, MPI_Group *)
  {
    cerr << "MPI_Group_incl" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_excl(MPI_Group, int, int *, MPI_Group *)
  {
    cerr << "MPI_Group_excl" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_range_incl(MPI_Group, int, int [][3], MPI_Group *)
  {
    cerr << "MPI_Group_range_incl" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_range_excl(MPI_Group, int, int [][3], MPI_Group *)
  {
    cerr << "MPI_Group_range_excl" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Group_free(MPI_Group *)
  {
    cerr << "MPI_Group_free" << endl ;
    err_report();
    return -1;
  }

  /* 5.4 */

  int  MPI_Comm_size(MPI_Comm, int *val)
  {
    *val = 1 ;
    return 0;
  }

  int  MPI_Comm_rank(MPI_Comm, int *val)
  {
    *val = 0 ;
    return 0;
  }

  int  MPI_Comm_compare(MPI_Comm, MPI_Comm, int *)
  {
    cerr << "MPI_Comm_compare" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Comm_dup(MPI_Comm, MPI_Comm *)
  {
    cerr << "MPI_Comm_dup" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Comm_create(MPI_Comm, MPI_Group, MPI_Comm *)
  {
    cerr << "MPI_Comm_create" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Comm_split(MPI_Comm, int, int, MPI_Comm *)
  {
    cerr << "MPI_Comm_split" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Comm_free(MPI_Comm *)
  {
    cerr << "MPI_Comm_free" << endl ;
    err_report();
    return -1;
  }

  /* 5.6 */

  int  MPI_Comm_test_inter(MPI_Comm, int *)
  {
    cerr << "MPI_Comm_test_inter" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Comm_remote_size(MPI_Comm, int *)
  {
    cerr << "MPI_Comm_remote_size" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Comm_remote_group(MPI_Comm, MPI_Group *)
  {
    cerr << "MPI_Comm_remote_group" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Intercomm_create(MPI_Comm, int, MPI_Comm, int, int, MPI_Comm *)
  {
    cerr << "MPI_Intercomm_cerate" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Intercomm_merge(MPI_Comm, int, MPI_Comm *)
  {
    cerr << "MPI_Intercomm_merge" << endl ;
    err_report();
    return -1;
  }

  /* 5.7 */

  int  MPI_Keyval_create(MPI_Copy_function *, MPI_Delete_function *, int *, void *)
  {
    cerr << "MPI_Keyval_create" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Keyval_free(int *)
  {
    cerr << "MPI_Keyval_free" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Attr_put(MPI_Comm, int, void *)
  {
    cerr << "MPI_Keyval_put" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Attr_get(MPI_Comm, int, void *, int *)
  {
    cerr << "MPI_Attr_get" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Attr_delete(MPI_Comm, int)
  {
    cerr << "MPI_Attr_delete" << endl ;
    err_report();
    return -1;
  }

  /* 6.5 */

  int  MPI_Cart_create(MPI_Comm, int, const int *, const int *, int, MPI_Comm *)
  {
    cerr << "MPI_Cart_create" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Dims_create(int, int, int *)
  {
    cerr << "MPI_Dims_create" << endl ;
    err_report();
    return -1;
  }

  int  MPI_Graph_create(MPI_Comm, int, int *, int *, int, MPI_Comm *)
  {
    err_report();
    return -1;
  }

  int  MPI_Topo_test(MPI_Comm, int *)
  {
    err_report();
    return -1;
  }

  int  MPI_Graphdims_get(MPI_Comm, int *, int *)
  {
    err_report();
    return -1;
  }

  int  MPI_Graph_get(MPI_Comm, int, int, int *, int *)
  {
    err_report();
    return -1;
  }

  int  MPI_Cartdim_get(MPI_Comm, int *)
  {
    err_report();
    return -1;
  }

  int  MPI_Cart_get(MPI_Comm, int, int *, int *, int *){
    err_report();
    return -1;
  }

  int  MPI_Cart_rank(MPI_Comm, int *, int *)
  {
    err_report();
    return -1;
  }

  int  MPI_Cart_coords(MPI_Comm, int, int, int *)
  {
    err_report();
    return -1;
  }

  int  MPI_Graph_neighbors_count(MPI_Comm, int, int *)
  {
    err_report();
    return -1;
  }

  int  MPI_Graph_neighbors(MPI_Comm, int, int, int *)
  {
    err_report();
    return -1;
  }

  int  MPI_Cart_shift(MPI_Comm, int, int, int *, int *)
  {
    err_report();
    return -1;
  }

  int  MPI_Cart_sub(MPI_Comm, const int *, MPI_Comm *)
  {
    err_report();
    return -1;
  }

  int  MPI_Cart_map(MPI_Comm, int, const int *, const int *, int *)
  {
    err_report();
    return -1;
  }

  int  MPI_Graph_map(MPI_Comm, int, int *, int *, int *)
  {
    err_report();
    return -1;
  }

  /* 7.1 */

  int  MPI_Get_processor_name(char *, int *)
  {
    err_report();
    return -1;
  }

  /* 7.2 */

  int  MPI_Errhandler_create(MPI_Handler_function *, MPI_Errhandler *)
  {
    return 0 ;
  }

  int  MPI_Errhandler_set(MPI_Comm, MPI_Errhandler)
  {
    return 0;
  }

  int  MPI_Errhandler_get(MPI_Comm, MPI_Errhandler *)
  {
    err_report();
    return -1;
  }

  int  MPI_Errhandler_free(MPI_Errhandler *)
  {
    err_report();
    return -1;
  }

  int  MPI_Error_string(int, char *, int *)
  {
    err_report();
    return -1;
  }

  /* 7.3 */

  int  MPI_Error_class(int, int *)
  {
    err_report();
    return -1;
  }

  /* 7.4 */

  double  MPI_Wtime(void)
  {
    timeval tv ;
    struct timezone tz ;
    gettimeofday(&tv,&tz) ;
    return double(tv.tv_sec)+1e-6*double(tv.tv_usec) ;
  }

  double  MPI_Wtick(void)
  {
    return 1e-6 ;
  }

  /* 7.5 */

  int  MPI_Init(int *, char ***)
  { 
    return 0;
  }

  int  MPI_Finalize(void)
  {
    return 0;
  }

  int  MPI_Initialized(int *)
  {
    return 0;
  }

  int  MPI_Abort(MPI_Comm, int)
  {
#ifndef NO_CSTDLIB
    std::exit(-1);
#else
    exit(-1) ;
#endif
    return -1;
  }

  /* 8.3 */

  int  MPI_Pcontrol(int, ...)
  {
    err_report();
    return -1;
  }

  /********************/
  /* MPI-1.2 bindings */
  /********************/

  int  MPI_Get_version(int *, int *)
  {
    return 1;
  }


  /*************************************/
  /* MPI-2 bindings, sorted by chapter */
  /*************************************/

  /* 4.10 */

  int  MPI_Info_create(MPI_Info *)
  {
    //    err_report();
    return -1;
  }

  int MPI_Info_delete(MPI_Info, char *)
  {
    //    err_report();
    return -1;
  }

  int MPI_Info_dup(MPI_Info, MPI_Info *)
  {
    //    err_report();
    return -1;
  }

  int MPI_Info_free(MPI_Info *)
  {
    //    err_report();
    return -1;
  }

  int MPI_Info_get(MPI_Info, char *, int, char *, int *)
  {
    err_report();
    return -1;
  }

  int MPI_Info_get_nkeys(MPI_Info, int *)
  {
    err_report();
    return -1;
  }

  int MPI_Info_get_nthkey(MPI_Info, int, char *)
  {
    err_report();
    return -1;
  }

  int MPI_Info_get_valuelen(MPI_Info, char *, int *, int *)
  {
    err_report();
    return -1;
  }

  int MPI_Info_set(MPI_Info, char *, char *)
  {
    err_report();
    return -1;
  }

  /* 4.11 */

  int MPI_Alloc_mem(MPI_Aint,MPI_Info,void *)
  {
    err_report();
    return -1;
  }

  int MPI_Free_mem(void *)
  {
    err_report();
    return -1;
  }

  /* 4.12 */

  MPI_Fint MPI_Info_c2f(MPI_Info)
  {
    err_report();
    return -1;
  }

  MPI_Info MPI_Info_f2c(MPI_Fint)
  {
    err_report();
    return 0;
  }

  MPI_Fint MPI_Comm_c2f(MPI_Comm)
  {
    err_report();
    return 0;
  }

  MPI_Comm MPI_Comm_f2c(MPI_Fint)
  {
    err_report();
    return 0;
  }

  MPI_Fint MPI_Type_c2f(MPI_Datatype)
  {
    err_report();
    return 0;
  }

  MPI_Datatype MPI_Type_f2c(MPI_Fint)
  {
    err_report();
    return 0;
  }

  MPI_Fint MPI_Group_c2f(MPI_Group)
  {
    err_report();
    return 0;
  }

  MPI_Group MPI_Group_f2c(MPI_Fint)
  {
    err_report();
    return 0;
  }

  MPI_Fint MPI_Request_c2f(MPI_Request)
  {
    err_report();
    return 0;
  }

  MPI_Request MPI_Request_f2c(MPI_Fint)
  {
    err_report();
    return 0;
  }

  MPI_Fint MPI_Op_c2f(MPI_Op)
  {
    err_report();
    return 0;
  }

  MPI_Op MPI_Op_f2c(MPI_Fint)
  {
    err_report();
    return 0;
  }


  /* 4.14 */

  int  MPI_Type_create_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *)
  {
    err_report();
    return 0;
  }

  int  MPI_Type_create_hindexed(int, const int *, const MPI_Aint *, MPI_Datatype, MPI_Datatype *)
  {
    err_report();
    return 0;
  }

  int  MPI_Type_create_struct(int, const int *, const MPI_Aint *,const MPI_Datatype *, MPI_Datatype *)
  {
    err_report();
    return 0;
  }

  int MPI_Get_address(void *, MPI_Aint *)
  {
    err_report();
    return 0;
  }

  /* 5.3 */

  int MPI_Comm_spawn(const char *, char **, int, MPI_Info, int, MPI_Comm, MPI_Comm *, int *)
  {
    err_report();
    return 0;
  }

  int MPI_Comm_spawn_multiple(int , char **, char ***, const int *, const MPI_Info *,
			      int , MPI_Comm, MPI_Comm *, int *)
  {
    err_report();
    return 0;
  }

  int MPI_Comm_get_parent(MPI_Comm *)
  {
    err_report();
    return 0;
  }

  /* referenced in 4.12 of MPI-2 standard with other transfer of handles functions */

  MPI_Fint MPI_Win_c2f(MPI_Win)
  {
    err_report(); 
    return 0;
  }

  MPI_Win MPI_Win_f2c(MPI_Fint)
  {
    err_report();
    return 0;
  }

  /* 6.2 */

  int  MPI_Win_create(void *, MPI_Aint, int, MPI_Info, MPI_Comm, MPI_Win *)
  {
    err_report();
    return 0;
  }

  int  MPI_Win_fence(int, MPI_Win)
  {
    err_report();
    return 0;
  }

  int  MPI_Win_free(MPI_Win *)
  {
    err_report();
    return 0;
  }

  /* 6.3 */


  int  MPI_Put(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
	       MPI_Win)
  {
    err_report();
    return 0;
  }

  int  MPI_Get(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
	       MPI_Win)
  {
    err_report();
    return 0;
  }

  int  MPI_Accumulate(const void *, int, MPI_Datatype, int, MPI_Aint, int,
		      MPI_Datatype, MPI_Op, MPI_Win)
  {
    err_report();
    return 0;
  }

  /* 7.5 */

  int  MPI_Type_get_envelope(MPI_Datatype, int *, int *, int *, int *)
  {
    err_report();
    return 0;
  }

  int  MPI_Type_get_contents(MPI_Datatype, int, int, int, int *, MPI_Aint *, MPI_Datatype *)
  {
    err_report();
    return 0;
  }

  /* 7.8 */

  int  MPI_Type_dup(MPI_Datatype, MPI_Datatype *)
  {
    err_report();
    return 0;
  }

  /* 8.7 */

  int  MPI_Init_thread(int *, char ***, int, int *)
  {
    err_report();
    return 0;
  }

  int  MPI_Query_thread(int *)
  {
    err_report();
    return 0;
  }

  int  MPI_Is_thread_main(int *)
  {
    err_report();
    return 0;
  }

  /* 9.6 */

  int  MPI_Finalized(int *)
  {
    err_report();
    return 0;
  }
  int  MPI_Comm_create_errhandler(MPI_Handler_function *, MPI_Errhandler *) {
    return 0 ;
  }
 
  int  MPI_Comm_set_errhandler(MPI_Comm, MPI_Errhandler) {
    return 0 ;
  }

  // MPI_IO
  int MPI_File_open(MPI_Comm comm, const char *filename, int amode,
		    MPI_Info info, MPI_File *fh)
  {
    err_report();
    return 0;
  }

  int MPI_File_close(MPI_File *fh)
  {
    err_report();
    return 0;
  }
    
  int MPI_File_read(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
		    MPI_Status *status)
  {
    err_report();
    return 0;
  }
    
  int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
			   int count, MPI_Datatype datatype, 
			   MPI_Status *status)
  {
    err_report();
    return 0;
  }
    
  int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf,
		       int count, MPI_Datatype datatype, MPI_Status *status)
  {
    err_report();
    return 0;
  }
    
  int MPI_File_write(MPI_File fh, const void *buf, int count,
		     MPI_Datatype datatype, MPI_Status *status)
  {
    err_report();
    return 0;
  }
    
  int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf,
			int count, MPI_Datatype datatype, MPI_Status *status)
  {
    err_report();
    return 0;
  }
    
  int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf,
			    int count, MPI_Datatype datatype, 
			    MPI_Status *status)
  {
    err_report();
    return 0;
  }
    
  
}
#endif
