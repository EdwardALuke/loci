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
#include <hdf5_readwrite.h>
#include <distribute_io.h>
#include <fact_db.h>
#include <string>
using std::string;
namespace Loci {


  /*
   * Create the appropriate File access property list
   */
  namespace hdf5_const {
    const int PPN(1);
    const int facc_type(FACC_MPIO);		/*Test file access type */
    const int dxfer_coll_type(DXFER_COLLECTIVE_IO); /*use collective IO*/
    // extern const int dxfer_coll_type(DXFER_INDEPENDENT_IO); /*use independent IO*/
  }

  /*turn on parallel io, must be called by all processes collectively */
  void set_parallel_io(bool io_type){
    if(io_type){
#ifdef H5_HAVE_PARALLEL
      if(MPI_processes > 1){ 
	use_parallel_io = true ;
      }else {
	use_parallel_io = false ;
	if(MPI_rank == 0) cerr << "Parallel I/O not used when running on a single CPU!" << endl ;
      }
#else
      if(MPI_rank == 0)cerr << "Cannot enable parallel I/O, Loci linked with serial HDF5 library!" << endl ;
      use_parallel_io = false ;     
#endif
    }else{
      use_parallel_io = false ;
    }
  }
  
  int  create_mpi_info(MPI_Info *info){//create a info and set the
    // *info = MPI_INFO_NULL;
    //return 0;
    MPI_Info_create(info);
    
#ifdef H5_HAVE_PARALLEL //optional setup

    //MPI predefined hints
    MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN) ;

    if(MPI_Info_set(*info, "striping_unit", "8388608")!=MPI_SUCCESS)
      cerr << "failed to set MPI info for 'striping_unit'" << endl ;

    if(MPI_Info_set(*info,"stripping_factor","4")!=MPI_SUCCESS)
      cerr << "failed to set MPI info for 'stripping_factor'" << endl ;

    // This was causing problems
#ifdef HACK
    if(MPI_Info_set(*info,"cb_nodes","4")!=MPI_SUCCESS)
      cerr << "failed to set MPI info for 'cb_nodes'" << endl ;

    if(MPI_Info_set(*info,"cb_buffer_size","4194304")!=MPI_SUCCESS)
      cerr << "failed to set MPI info for 'cb_buffer_size'" << endl ;
#endif
      
    //new algorithm parameters
    if(MPI_Info_set(*info,"ind_rd_buffer_size","41943040")!=MPI_SUCCESS)
      cerr << "failed to set MPI info for 'ind_rd_buffer_size'" << endl ;

    if(MPI_Info_set(*info,"ind_wr_buffer_size","5242880")!=MPI_SUCCESS)
      cerr << "failed to set MPI info for 'ind_wr_buffer_size'" << endl ;

    //platform-specific hints
    if(MPI_Info_set(*info,"IBM_largeblock_io","true")!=MPI_SUCCESS)
      cerr << "failed to set MPI info for 'IBM_largeblock_io'" << endl ;

    //      if(MPI_Info_set(*info,"H5F_ACS_CORE_WRITE_TRACKING_PAGE_SIZE_DEF","524288")!=MPI_SUCCESS)
    //	cerr << "failed to set MPI info for 'H5F_ACS_CORE_WRITE_TRACKING_PAGE_SIZE_DEF'" << endl ;

    //dissble datat sieving in read
    if(MPI_Info_set(*info,"romio_ds_read","disable")!=MPI_SUCCESS)
      cerr << "failed to set MPI info for 'romio_ds_read'" << endl ;

    // disable data sieving in write
    if(MPI_Info_set(*info,"romio_ds_write","disable")!=MPI_SUCCESS)
      cerr << "failed to set MPI info for 'romio_ds_write'" << endl ;

    //enable aggregation
    if(MPI_Info_set(*info,"romio_cb_write","enable")!=MPI_SUCCESS)
      cerr << "failed to set MPI info for 'romio_cb_write'" << endl ;
	
    //Setting the environment variable
    //MPICH_MPIIO_HINTS_DISPLAY=1 to print out available I/O hints and their values
    //aggreation: processors with fast connection perform io for others
    MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_ARE_FATAL) ;

#endif
    return 0;
  }
  
  hid_t
  create_faccess_plist(MPI_Comm comm, MPI_Info info, int l_facc_type)
  {
#ifdef H5_HAVE_PARALLEL
    hid_t ret_pl = -1;
    ret_pl = H5Pcreate (H5P_FILE_ACCESS);
    WARN(ret_pl<0) ;
    if (l_facc_type == FACC_DEFAULT){
      return (ret_pl);
    }
    if (l_facc_type == FACC_MPIO){
      /* set Parallel access with communicator */
      H5Pset_fapl_mpio(ret_pl, comm, info);

      //unleased collective metadata read/write options to avoid read storm
      //one process read small data and broadcast it to avoid I/O access
      // ret = H5Pset_all_coll_metadata_ops(ret_pl, TRUE);
      //WARN(ret < 0);

      //ret = H5Pset_coll_metadata_write(ret_pl, TRUE);
      //WARN(ret<0);
      return(ret_pl);
    }
    return(ret_pl);
#else
    return H5P_DEFAULT ;
#endif
  }

  hid_t
  create_xfer_plist(int l_xfer_type){
#ifdef H5_HAVE_PARALLEL
    hid_t xfer_plist = H5Pcreate (H5P_DATASET_XFER);
    if(l_xfer_type == DXFER_COLLECTIVE_IO)
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    else
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_INDEPENDENT);
    return xfer_plist;
#else
    return H5P_DEFAULT ;
#endif
  }

  using std::vector ;
  using std::ostringstream ;
 

 
  
  /*serial/parallel io
    next several functions, io data is the same to all processes
    so in serial version, only process 0 call these functions, and it is caller's responsibility to broadcast the data after reading
    if use_parallel_io:
    for parallel writing, since all structure modification procedures such as creating groups, dataset, attributes etc. need to be collective,
    all processes call these functions, only process 0 need perform writing, using independent data transfer 
    for parallel reading, to avoid read storm, only process 0 call these functions, and it is caller's responsibility to broadcast the data
  */
  
  //**************************************************************************/
  //in case of parallel io, all processes call it, but only process 0 perform writing
  void HDF5_WriteDomain(hid_t group_id, const entitySet &en , MPI_Comm comm)
  { 
    hsize_t dimension=0;

    int rank = 1;     
    int num_intervals = en.num_intervals(); 

    hid_t dataspace, dataset;
    hid_t datatype  = H5T_NATIVE_INT;
    int prank = 0 ;
    //int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    //MPI_Comm_size(comm,&np) ;

    if( num_intervals < 1) return;

    dimension = num_intervals*2; //size of 1D Array
    if(dimension == 0) return;
    dataspace = H5Screate_simple(rank, &dimension, NULL);
    dataset   = H5Dcreate(group_id, "Interval Set", datatype, dataspace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    if(prank == 0){ //only process 0 perform writing
      interval *it = new interval[num_intervals];
  
      int *data = new int[num_intervals*2];      

      for(int i=0;i<num_intervals;i++){
	it[i]       = en[i];
	data[i*2]   = it[i].first;
	data[i*2+1] = it[i].second;
      }
  
      //no xfer_plist here, must be independently, which is default
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
   

      delete [] data;
      delete [] it;
    }

    H5Sclose(dataspace);
    H5Dclose(dataset);

  }
  
    

  //**************************************************************************/
  //   never been called, assume only process 0 will call it
  void HDF5_Local2Global( hid_t group_id, entitySet &eset, Map &l_g)//never been used
  {

    int        indx = 0;
    hsize_t    dimension;
    hid_t      dataset, dataspace;

    dataset    = H5Dopen(group_id, "Map",H5P_DEFAULT);

    if( dataset > 0) {
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);

      int *data = new int[dimension];
        
      hid_t xfer_plist = H5P_DEFAULT;
     
      H5Dread( dataset, H5T_NATIVE_INT, H5S_ALL, dataspace,
	       xfer_plist, data);

   
      l_g.allocate(eset);
      entitySet :: const_iterator ci;
      for( ci = eset.begin(); ci != eset.end(); ++ci) 
	l_g[*ci] = data[indx++];

      delete [] data;
      if(xfer_plist != H5P_DEFAULT) H5Pclose(xfer_plist);
      H5Sclose(dataspace);
      H5Dclose(dataset);
        
    }

  }

  //**************************************************************************/
  //assume only process 0 will call it, and then broadcast the data
  void HDF5_ReadDomain( hid_t group_id, entitySet &eset)
  {
    hsize_t    dimension;
    hid_t      dataset, dataspace;

    H5Eset_auto (H5E_DEFAULT,NULL, NULL);

    eset = EMPTY;
    dataset  = H5Dopen(group_id, "Interval Set",H5P_DEFAULT);

    if( dataset > 0 ) {
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);

      int *data = new int[dimension];
      hid_t xfer_plist = H5P_DEFAULT ;
      
      H5Dread( dataset, H5T_NATIVE_INT, H5S_ALL, dataspace,
	       xfer_plist, data);

      eset = EMPTY;
      for(size_t i=0;i< dimension;i+=2){
	eset |= interval(data[i],data[i+1]);
      }
      delete [] data;
      if(xfer_plist != H5P_DEFAULT)
	H5Pclose(xfer_plist);
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }

   
  //**************************************************************************/

  //in case of  parallel io, all processes call it , only process 0 perform writing 
  void HDF5_WriteVecSize(hid_t group_id, const int &size,  MPI_Comm comm  )//FOR dynamic map/store, not used currently
  {
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    hsize_t dimension = 1;
    int     rank = 1;
    if(dimension == 0) return;
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_INT;

    hid_t vDataset   = H5Dcreate(group_id, "VecSize", vDatatype, vDataspace,
				 H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    //no xfer_plist here, must be independently, which is default
    if(prank == 0) 
      H5Dwrite(vDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &size);

    H5Sclose( vDataspace );
    H5Dclose( vDataset   );

  }

  
  //**************************************************************************/

  //to avoid read storm, assume only process 0 will call it, and broadcast the data
  void HDF5_ReadVecSize(hid_t group_id, int *size )//FOR dynamic map/store, not used currently
  {

    int     rank=1;
    hsize_t dimension = 1;

    H5Eset_auto (H5E_DEFAULT,NULL, NULL);

    *size = 0;
  
    hid_t vDataset   = H5Dopen( group_id, "VecSize", H5P_DEFAULT);

    if( vDataset > 0) {
      if(dimension == 0) return;
      hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
      hid_t vDatatype  = H5T_NATIVE_INT;
      hid_t xfer_plist = H5P_DEFAULT;
      H5Dread( vDataset, vDatatype, H5S_ALL, vDataspace, xfer_plist, size);
      if(xfer_plist != H5P_DEFAULT) H5Pclose(xfer_plist);
      H5Sclose( vDataspace );
      H5Dclose( vDataset   );
    }

  }
}

