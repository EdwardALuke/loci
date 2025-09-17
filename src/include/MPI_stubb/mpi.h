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
#ifndef MPI_H_STUBB
#define MPI_H_STUBB



#ifndef MPI_STUBB
#define MPI_STUBB
#endif

#ifdef __cplusplus
extern "C" {
#endif

/******************/
/* MPI-1 bindings */
/******************/


#define MPI_BOTTOM		((MPI_Aint)0)
#define MPI_STATUSES_IGNORE     (0)
#define MPI_STATUS_IGNORE       (0)
  
typedef long				MPI_Aint;
typedef unsigned int		MPI_Request;
typedef unsigned int		MPI_Group;
typedef unsigned int		MPI_Comm;
typedef unsigned int		MPI_Errhandler;
typedef unsigned int		MPI_Op;
typedef unsigned int		MPI_Datatype;

typedef struct { 
	int MPI_SOURCE;
	int MPI_TAG;
	int MPI_ERROR;
	int size;
	int reserved[2];
} MPI_Status;

enum {
	MPI_COMM_NULL		= 0,
	MPI_COMM_WORLD		= 1,
	MPI_COMM_SELF		= 2
};

enum {
	MPI_ERRHANDLER_NULL	= 0,
	MPI_ERRORS_ARE_FATAL	= 1,
	MPI_ERRORS_RETURN	= 2
};

enum {
	MPI_GROUP_NULL		= 0,
	MPI_GROUP_EMPTY		= 1
};

enum {
	MPI_REQUEST_NULL	= 0
};

enum {
	MPI_OP_NULL		= 0,
	MPI_MAX			= 1,
	MPI_MIN			= 2,
	MPI_SUM			= 3,
	MPI_PROD		= 4,
	MPI_LAND		= 5,
	MPI_BAND 		= 6,
	MPI_LOR			= 7,
	MPI_BOR			= 8,
	MPI_LXOR		= 9,
	MPI_BXOR		= 10,
	MPI_MAXLOC		= 11,
	MPI_MINLOC		= 12
};

enum {
	MPI_DATATYPE_NULL	= 0,

	MPI_CHAR		= 1,
	MPI_SHORT		= 2,
	MPI_INT			= 3,
	MPI_LONG		= 4,
	MPI_UNSIGNED_CHAR	= 5,
	MPI_UNSIGNED_SHORT	= 6,
	MPI_UNSIGNED		= 7,
	MPI_UNSIGNED_LONG	= 8,
	MPI_FLOAT		= 9,
	MPI_DOUBLE		= 10,
	MPI_LONG_DOUBLE		= 11,
	MPI_LONG_LONG		= 12,

	MPI_INTEGER		= 13,
	MPI_REAL		= 14,
	MPI_DOUBLE_PRECISION	= 15,
	MPI_COMPLEX		= 16,
	MPI_DOUBLE_COMPLEX	= 17,
	MPI_LOGICAL		= 18,
	MPI_CHARACTER		= 19,
	MPI_INTEGER1		= 20,
	MPI_INTEGER2		= 21,
	MPI_INTEGER4		= 22,
	MPI_INTEGER8		= 23,
	MPI_REAL4		= 24,
	MPI_REAL8		= 25,
	MPI_REAL16		= 26,

	MPI_BYTE		= 27,
	MPI_PACKED		= 28,
	MPI_UB			= 29,
	MPI_LB			= 30,

	MPI_FLOAT_INT		= 31,
	MPI_DOUBLE_INT		= 32,
	MPI_LONG_INT		= 33,
	MPI_2INT		= 34,
	MPI_SHORT_INT		= 35,
	MPI_LONG_DOUBLE_INT	= 36,

	MPI_2REAL		= 37,
	MPI_2DOUBLE_PRECISION	= 38,
	MPI_2INTEGER		= 39,
	_MPI_SGI_TYPE_LAST
};

#define MPI_INT32_T MPI_INTEGER
#define MPI_INT64_T MPI_LONG_LONG
  
int MPI_GET_TYPE_SIZE(MPI_Datatype type);


#define MPI_LONG_LONG_INT	MPI_LONG_LONG


enum {
	MPI_SUCCESS			= 0,

	/* These 19 error codes are specified by the MPI-1 standard */

	MPI_ERR_BUFFER			= 1,
	MPI_ERR_COUNT			= 2,
	MPI_ERR_TYPE			= 3,
	MPI_ERR_TAG			= 4,
	MPI_ERR_COMM			= 5,
	MPI_ERR_RANK			= 6,
	MPI_ERR_REQUEST			= 7,
	MPI_ERR_ROOT			= 8,
	MPI_ERR_GROUP			= 9,
	MPI_ERR_OP			= 10,
	MPI_ERR_TOPOLOGY		= 11,
	MPI_ERR_DIMS			= 12,
	MPI_ERR_ARG			= 13,
	MPI_ERR_UNKNOWN			= 14,
	MPI_ERR_TRUNCATE		= 15,
	MPI_ERR_OTHER			= 16,
	MPI_ERR_INTERN			= 17,
	MPI_ERR_IN_STATUS		= 18,
	MPI_ERR_PENDING			= 19,


	/* Error codes 20-27 used by MPI on T3E systems. */

	/* These 34 error codes are specified by the MPI-2 standard */

	MPI_ERR_ACCESS			= 28,
	MPI_ERR_AMODE			= 29,
	MPI_ERR_ASSERT			= 30,
	MPI_ERR_BAD_FILE		= 31,
	MPI_ERR_BASE			= 32,
	MPI_ERR_CONVERSION		= 33,
	MPI_ERR_DISP			= 34,
	MPI_ERR_DUP_DATAREP		= 35,
	MPI_ERR_FILE_EXISTS		= 36,
	MPI_ERR_FILE_IN_USE		= 37,
	MPI_ERR_FILE			= 38,
	MPI_ERR_INFO_KEY		= 39,
	MPI_ERR_INFO_NOKEY		= 40,
	MPI_ERR_INFO_VALUE		= 41,
	MPI_ERR_INFO			= 42,
	MPI_ERR_IO			= 43,
	MPI_ERR_KEYVAL			= 44,
	MPI_ERR_LOCKTYPE		= 45,
	MPI_ERR_NAME			= 46,
	MPI_ERR_NO_MEM			= 47,
	MPI_ERR_NOT_SAME		= 48,
	MPI_ERR_NO_SPACE		= 49,
	MPI_ERR_NO_SUCH_FILE		= 50,
	MPI_ERR_PORT			= 51,
	MPI_ERR_QUOTA			= 52,
	MPI_ERR_READ_ONLY		= 53,
	MPI_ERR_RMA_CONFLICT		= 54,
	MPI_ERR_RMA_SYNC		= 55,
	MPI_ERR_SERVICE			= 56,
	MPI_ERR_SIZE			= 57,
	MPI_ERR_SPAWN			= 58,
	MPI_ERR_UNSUPPORTED_DATAREP	= 59,
	MPI_ERR_UNSUPPORTED_OPERATION	= 60,
	MPI_ERR_WIN			= 61,
	MPI_ERR_LASTCODE		= 100	/* last built-in error code */
};

enum {
	MPI_KEYVAL_INVALID	= 0,
	MPI_TAG_UB		= 1,
	MPI_HOST		= 2,
	MPI_IO			= 3,
	MPI_WTIME_IS_GLOBAL	= 4,
	MPI_UNIVERSE_SIZE	= 9,
	MPI_APPNUM		= 11
};

enum {
	MPI_IDENT		= 0,
	MPI_CONGRUENT		= 1,
	MPI_SIMILAR		= 2,
	MPI_UNEQUAL		= 3
};

enum {
	MPI_GRAPH		= 1,
	MPI_CART		= 2
};

enum {
	MPI_UNDEFINED		= -3,
	MPI_ANY_SOURCE		= -2,
	MPI_PROC_NULL		= -1
};

enum {
	MPI_ANY_TAG		= -1
};

enum {
	MPI_BSEND_OVERHEAD	= 32
};

enum {
	MPI_MAX_PROCESSOR_NAME	= 256
};

enum {
	MPI_MAX_ERROR_STRING	= 256
};


typedef int MPI_Copy_function(MPI_Comm, int, void *, void *, void *, int *);
typedef int MPI_Delete_function(MPI_Comm, int, void *, void *);
typedef void MPI_Handler_function(MPI_Comm *, int *, ...);
typedef void MPI_User_function(void *, void *, int *, MPI_Datatype *); 

MPI_Copy_function		MPI_NULL_COPY_FN, MPI_DUP_FN;
MPI_Delete_function		MPI_NULL_DELETE_FN;


/*************************************/
/* MPI-1 bindings, sorted by chapter */
/*************************************/


/* 3.2 */

void err_report();

 int  MPI_Send(const void *, int, MPI_Datatype, int, int, MPI_Comm);

 int  MPI_Recv(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);

 int  MPI_Get_count(const MPI_Status *, MPI_Datatype, int *);

/* 3.4 */

 int  MPI_Bsend(const void *, int, MPI_Datatype, int, int, MPI_Comm);

 int  MPI_Ssend(const void *, int, MPI_Datatype, int, int, MPI_Comm);

 int  MPI_Rsend(const void *, int, MPI_Datatype, int, int, MPI_Comm);

/* 3.6 */

 int  MPI_Buffer_attach(void *, int);

 int  MPI_Buffer_detach(void *, int *);

/* 3.7 */

 int  MPI_Isend(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

 int  MPI_Ibsend(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

 int  MPI_Issend(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

 int  MPI_Irsend(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

 int  MPI_Irecv(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

 int  MPI_Wait(MPI_Request *, MPI_Status *);

 int  MPI_Test(MPI_Request *, int *, MPI_Status *);

 int  MPI_Request_free(MPI_Request *);

 int  MPI_Waitany(int, MPI_Request *, int *, MPI_Status *);

 int  MPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *);

 int  MPI_Waitall(int, MPI_Request *, MPI_Status *);

 int  MPI_Testall(int, MPI_Request *, int *, MPI_Status *);

 int  MPI_Waitsome(int, MPI_Request *, int *, int *, MPI_Status *);

 int  MPI_Testsome(int, MPI_Request *, int *, int *, MPI_Status *);

/* 3.8 */

 int  MPI_Iprobe(int, int, MPI_Comm, int *, MPI_Status *);

 int  MPI_Probe(int, int, MPI_Comm, MPI_Status *);

 int  MPI_Cancel(MPI_Request *);

 int  MPI_Test_cancelled(const MPI_Status *, int *);

/* 3.9 */

 int  MPI_Send_init(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

 int  MPI_Bsend_init(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

 int  MPI_Ssend_init(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

 int  MPI_Rsend_init(const void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

 int  MPI_Recv_init(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);

 int  MPI_Start(MPI_Request *);

 int  MPI_Startall(int, MPI_Request *);

/* 3.10 */

 int  MPI_Sendrecv(const void *, int, MPI_Datatype, int, int, void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);

 int  MPI_Sendrecv_replace(void *, int, MPI_Datatype, int, int, int, int, MPI_Comm, MPI_Status *);

/* 3.12 */

 int  MPI_Type_contiguous(int, MPI_Datatype, MPI_Datatype *);

 int  MPI_Type_vector(int, int, int, MPI_Datatype, MPI_Datatype *);

 int  MPI_Type_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);

 int  MPI_Type_indexed(int, const int *, const int *, MPI_Datatype, MPI_Datatype *);

 int  MPI_Type_hindexed(int, int *, MPI_Aint *, MPI_Datatype, MPI_Datatype *);

 int  MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);

 int  MPI_Address(void *, MPI_Aint *);

 int  MPI_Type_extent(MPI_Datatype, MPI_Aint *);

 int  MPI_Type_size(MPI_Datatype, int *);

 int  MPI_Type_lb(MPI_Datatype, MPI_Aint *);

 int  MPI_Type_ub(MPI_Datatype, MPI_Aint *);

 int  MPI_Type_commit(MPI_Datatype *);

 int  MPI_Type_create_resized(MPI_Datatype, MPI_Aint,MPI_Aint,MPI_Datatype *) ;

 int  MPI_Type_free(MPI_Datatype *);

 int  MPI_Get_elements(const MPI_Status *, MPI_Datatype, int *);

/* 3.13 */

 int  MPI_Pack(const void *buf, int count, MPI_Datatype datatype,
                     void *packbuf, int packsize, int *packpos,
                     MPI_Comm comm);

 int  MPI_Unpack(void *packbuf, int packsize, int *packpos,
                       void *buf, int count, MPI_Datatype datatype,
                       MPI_Comm comm);

 int  MPI_Pack_size(int incount, MPI_Datatype datatype , MPI_Comm comm, int *size);

/* 4.3 */

 int  MPI_Barrier(MPI_Comm);

/* 4.4 */

 int  MPI_Bcast(void *, int, MPI_Datatype, int, MPI_Comm);

/* 4.5 */

 int  MPI_Gather(const void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                       void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                       int root , MPI_Comm comm);

 int  MPI_Gatherv(const void *, int, MPI_Datatype, void *, const int *, const int *, MPI_Datatype, int, MPI_Comm);

/* 4.6 */

 int  MPI_Scatter(const void *, int, MPI_Datatype, void *, int, MPI_Datatype, int, MPI_Comm);

 int  MPI_Scatterv(const void *, const int *, const int *, MPI_Datatype, void *, int, MPI_Datatype, int, MPI_Comm);

/* 4.7 */

 int  MPI_Allgather(const void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		    void *recvbuf, int recvcnt, MPI_Datatype recvtype,
		    MPI_Comm comm );

 int  MPI_Allgatherv(const void *sendbuf, int sendcnt, MPI_Datatype sendtype,
		     void *recvbuf, const int *recvcounts, const int *displs,
		     MPI_Datatype recvtype, MPI_Comm comm);
/* 4.8 */

 int  MPI_Alltoall(const void *send, int ssz, MPI_Datatype dts, void *recv, int rsz, MPI_Datatype dtr, MPI_Comm comm);

  int  MPI_Alltoallv(const void *sendbuf,
		     const int *sendcnts, const int *sdispls,
		     MPI_Datatype sendtype,
		     void *recvbuf, const int *recvcnts, const int *rdispls,
		     MPI_Datatype recvtype, MPI_Comm comm);

/* 4.9 */

 int  MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
		 MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm);

 int  MPI_Op_create(MPI_User_function *, int, MPI_Op *);

 int  MPI_Op_free(MPI_Op *);

 int  MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
		    MPI_Datatype type, MPI_Op op , MPI_Comm comm);

/* 4.10 */

 int  MPI_Reduce_scatter(const void *, void *, const int *, MPI_Datatype, MPI_Op, MPI_Comm);

/* 4.11 */

 int  MPI_Scan(const void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm);
/* 5.3 */

 int  MPI_Group_size(MPI_Group, int *);

 int  MPI_Group_rank(MPI_Group, int *);

 int  MPI_Group_translate_ranks(MPI_Group, int, int *, MPI_Group, int *);

 int  MPI_Group_compare(MPI_Group, MPI_Group, int *);

 int  MPI_Comm_group(MPI_Comm, MPI_Group *);

 int  MPI_Group_union(MPI_Group, MPI_Group, MPI_Group *);

 int  MPI_Group_intersection(MPI_Group, MPI_Group, MPI_Group *);

 int  MPI_Group_difference(MPI_Group, MPI_Group, MPI_Group *);

 int  MPI_Group_incl(MPI_Group, int, int *, MPI_Group *);

 int  MPI_Group_excl(MPI_Group, int, int *, MPI_Group *);

 int  MPI_Group_range_incl(MPI_Group, int, int [][3], MPI_Group *);

 int  MPI_Group_range_excl(MPI_Group, int, int [][3], MPI_Group *);

 int  MPI_Group_free(MPI_Group *);

/* 5.4 */

 int  MPI_Comm_size(MPI_Comm, int *val);

 int  MPI_Comm_rank(MPI_Comm, int *val);

 int  MPI_Comm_compare(MPI_Comm, MPI_Comm, int *);

 int  MPI_Comm_dup(MPI_Comm, MPI_Comm *);

 int  MPI_Comm_create(MPI_Comm, MPI_Group, MPI_Comm *);

 int  MPI_Comm_split(MPI_Comm, int, int, MPI_Comm *);

 int  MPI_Comm_free(MPI_Comm *);

/* 5.6 */

 int  MPI_Comm_test_inter(MPI_Comm, int *);

 int  MPI_Comm_remote_size(MPI_Comm, int *);

 int  MPI_Comm_remote_group(MPI_Comm, MPI_Group *);

 int  MPI_Intercomm_create(MPI_Comm, int, MPI_Comm, int, int, MPI_Comm *);

 int  MPI_Intercomm_merge(MPI_Comm, int, MPI_Comm *);

/* 5.7 */

 int  MPI_Keyval_create(MPI_Copy_function *, MPI_Delete_function *, int *, void *);

 int  MPI_Keyval_free(int *);

 int  MPI_Attr_put(MPI_Comm, int, void *);

 int  MPI_Attr_get(MPI_Comm, int, void *, int *);

 int  MPI_Attr_delete(MPI_Comm, int);

/* 6.5 */

 int  MPI_Cart_create(MPI_Comm, int, const int *, const int *, int, MPI_Comm *);

 int  MPI_Dims_create(int, int, int *);

 int  MPI_Graph_create(MPI_Comm, int, int *, int *, int, MPI_Comm *);

 int  MPI_Topo_test(MPI_Comm, int *);

 int  MPI_Graphdims_get(MPI_Comm, int *, int *);

 int  MPI_Graph_get(MPI_Comm, int, int, int *, int *);

 int  MPI_Cartdim_get(MPI_Comm, int *);

 int  MPI_Cart_get(MPI_Comm, int, int *, int *, int *);

 int  MPI_Cart_rank(MPI_Comm, int *, int *);

 int  MPI_Cart_coords(MPI_Comm, int, int, int *);

 int  MPI_Graph_neighbors_count(MPI_Comm, int, int *);

 int  MPI_Graph_neighbors(MPI_Comm, int, int, int *);

 int  MPI_Cart_shift(MPI_Comm, int, int, int *, int *);

 int  MPI_Cart_sub(MPI_Comm, const int *, MPI_Comm *);

 int  MPI_Cart_map(MPI_Comm, int,const int *,const int *, int *);

 int  MPI_Graph_map(MPI_Comm, int, int *, int *, int *);

/* 7.1 */

 int  MPI_Get_processor_name(char *, int *);

/* 7.2 */

 int  MPI_Comm_create_errhandler(MPI_Handler_function *, MPI_Errhandler *);
 
 int  MPI_Comm_set_errhandler(MPI_Comm, MPI_Errhandler);

 int  MPI_Errhandler_create(MPI_Handler_function *, MPI_Errhandler *);

 int  MPI_Errhandler_set(MPI_Comm, MPI_Errhandler);

 int  MPI_Errhandler_get(MPI_Comm, MPI_Errhandler *);

 int  MPI_Errhandler_free(MPI_Errhandler *);

 int  MPI_Error_string(int, char *, int *);

/* 7.3 */

 int  MPI_Error_class(int, int *);

/* 7.4 */

 double  MPI_Wtime(void);

 double  MPI_Wtick(void);

/* 7.5 */

 int  MPI_Init(int *, char ***);

 int  MPI_Finalize(void);

 int  MPI_Initialized(int *);

 int  MPI_Abort(MPI_Comm, int);
/* 8.3 */

 int  MPI_Pcontrol(int, ...);

/********************/
/* MPI-1.2 bindings */
/********************/

#define MPI_VERSION		1
#define MPI_SUBVERSION		2

 int  MPI_Get_version(int *, int *);


/*************************************/
/* MPI-2 bindings, sorted by chapter */
/*************************************/

/* 4.10 */

typedef unsigned int          MPI_Info;

enum {
	MPI_INFO_NULL		= 0,
	MPI_MAX_INFO_KEY        = 255,
	MPI_MAX_INFO_VAL        = 1024
};

enum {
	MPI_FUNDAMENTAL		= -1
};

 int  MPI_Info_create(MPI_Info *);

 int MPI_Info_delete(MPI_Info, char *);

 int MPI_Info_dup(MPI_Info, MPI_Info *);

 int MPI_Info_free(MPI_Info *);

 int MPI_Info_get(MPI_Info, char *, int, char *, int *);

 int MPI_Info_get_nkeys(MPI_Info, int *);

 int MPI_Info_get_nthkey(MPI_Info, int, char *);

 int MPI_Info_get_valuelen(MPI_Info, char *, int *, int *);

 int MPI_Info_set(MPI_Info, char *, char *);

/* 4.11 */

 int MPI_Alloc_mem(MPI_Aint,MPI_Info,void *);

 int MPI_Free_mem(void *);

/* 4.12 */

typedef int MPI_Fint;

 MPI_Fint MPI_Info_c2f(MPI_Info);

 MPI_Info MPI_Info_f2c(MPI_Fint);

 MPI_Fint MPI_Comm_c2f(MPI_Comm);

 MPI_Comm MPI_Comm_f2c(MPI_Fint);

 MPI_Fint MPI_Type_c2f(MPI_Datatype);

 MPI_Datatype MPI_Type_f2c(MPI_Fint);

 MPI_Fint MPI_Group_c2f(MPI_Group);

 MPI_Group MPI_Group_f2c(MPI_Fint);

 MPI_Fint MPI_Request_c2f(MPI_Request);

 MPI_Request MPI_Request_f2c(MPI_Fint);

 MPI_Fint MPI_Op_c2f(MPI_Op);

 MPI_Op MPI_Op_f2c(MPI_Fint);


/* 4.14 */

 int  MPI_Type_create_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);

 int  MPI_Type_create_hindexed(int, const int *, const MPI_Aint *, MPI_Datatype, MPI_Datatype *);

 int  MPI_Type_create_struct(int, const int *,const MPI_Aint *, const MPI_Datatype *, MPI_Datatype *);

 int MPI_Get_address(void *, MPI_Aint *);

/* 5.3 */

#define MPI_ARGV_NULL              ((char **)NULL)
#define MPI_ARGVS_NULL              ((char ***)NULL)
#define MPI_ERRCODES_IGNORE        ((int *)NULL)

 int MPI_Comm_spawn(const char *, char **, int, MPI_Info, int, MPI_Comm, MPI_Comm *, int *);

 int MPI_Comm_spawn_multiple(int , char **, char ***, const int *, const MPI_Info *,
                            int , MPI_Comm, MPI_Comm *, int *);

 int MPI_Comm_get_parent(MPI_Comm *);

/* 6 */

/* MPI one-sided is supported only under ABI 64.  */
#if 	!_ABIN32

typedef unsigned int MPI_Win;

enum {
	MPI_WIN_NULL 	= 0
};

/* referenced in 4.12 of MPI-2 standard with other transfer of handles functions */

 MPI_Fint MPI_Win_c2f(MPI_Win);

 MPI_Win MPI_Win_f2c(MPI_Fint);

/* 6.2 */

 int  MPI_Win_create(void *, MPI_Aint, int, MPI_Info, MPI_Comm, MPI_Win *);

 int  MPI_Win_fence(int, MPI_Win);

 int  MPI_Win_free(MPI_Win *);

/* 6.3 */

 int  MPI_Put(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
	MPI_Win);

 int  MPI_Get(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
	MPI_Win);

 int  MPI_Accumulate(const void *, int, MPI_Datatype, int, MPI_Aint, int,
	MPI_Datatype, MPI_Op, MPI_Win);

/* 6.4 */

enum {
	MPI_MODE_NOCHECK	= 1,
	MPI_MODE_NOSTORE	= 2,
	MPI_MODE_NOPUT		= 4,
	MPI_MODE_NOPRECEDE	= 8,
	MPI_MODE_NOSUCCEED 	= 16
};


#endif	/* !_ABIN32 */


/* 7.5 */

 int  MPI_Type_get_envelope(MPI_Datatype, int *, int *, int *, int *);

 int  MPI_Type_get_contents(MPI_Datatype, int, int, int, int *, MPI_Aint *, MPI_Datatype *);

/* 7.8 */

 int  MPI_Type_dup(MPI_Datatype, MPI_Datatype *);

/* 8.6 */

enum {
	MPI_COMBINER_NAMED      	= (0),
	MPI_COMBINER_CONTIGUOUS 	= 0,
	MPI_COMBINER_VECTOR     	= 1,
	MPI_COMBINER_HVECTOR    	= 2,
	MPI_COMBINER_INDEXED    	= 3,
	MPI_COMBINER_HINDEXED   	= 4,
	MPI_COMBINER_STRUCT     	= 5,
	MPI_COMBINER_DARRAY		= 6,
	MPI_COMBINER_DUP		= 7,
	MPI_COMBINER_F90_COMPLES	= 8,
	MPI_COMBINER_F90_INTEGER	= 9,
	MPI_COMBINER_F90_REAL		= 10,
	MPI_COMBINER_HINDEXED_INTEGER	= 11,
	MPI_COMBINER_HVECTOR_INTEGER	= 12,
	MPI_COMBINER_INDEXED_BLOCK	= 13,
	MPI_COMBINER_RESIZED		= 14,
	MPI_COMBINER_STRUCT_INTEGER	= 15,
	MPI_COMBINER_SUBARRAY		= 16
};

/* 8.7 */

 int  MPI_Init_thread(int *, char ***, int, int *);
 int  MPI_Query_thread(int *);
 int  MPI_Is_thread_main(int *);

enum {
	MPI_THREAD_SINGLE		= 0,
	MPI_THREAD_FUNNELED		= 1,
	MPI_THREAD_SERIALIZED		= 2,
	MPI_THREAD_MULTIPLE		= 3
};

/* 9.6 */

 int  MPI_Finalized(int *);


  // IO
  typedef int MPI_File ;
  typedef int MPI_Offset ;

  enum {
	MPI_MODE_WRONLY                 = 1,
	MPI_MODE_RDONLY                 = 2,
	MPI_MODE_CREATE                 = 4,
	
	
  };

  int MPI_File_open(MPI_Comm comm, const char *filename, int amode, MPI_Info info,
		    MPI_File *fh) ;
  int MPI_File_close(MPI_File *fh) ;
  int MPI_File_read(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
		    MPI_Status *status) ;
  int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
			   int count, MPI_Datatype datatype, 
			   MPI_Status *status)  ;
  int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf,
		       int count, MPI_Datatype datatype, MPI_Status *status) ;
  int MPI_File_write(MPI_File fh, const void *buf, int count,
		     MPI_Datatype datatype, MPI_Status *status) ;
  int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf,
			int count, MPI_Datatype datatype, MPI_Status *status) ;
  int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf,
			    int count, MPI_Datatype datatype, 
			    MPI_Status *status) ;

#ifdef __cplusplus
}
#endif
#endif	/* MPI_H_INCLUDED */

