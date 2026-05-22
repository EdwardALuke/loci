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
#ifndef DISTRIBUTE_IO_H
#define DISTRIBUTE_IO_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <parSampleSort.h>

#include <store_rep.h>
#include <DMap.h>
#include <fact_db.h>
#include <mpi.h>
#include <store.h>
#include <MapVec.h>
#include <multiStore.h>

namespace Loci {

  namespace pio {
    /*
      namespace pio is for functions that have separate serial io version and parallel io version.
      these low-level functions are generally used by Loci-developers instead of application-developers
    */
    void write_storeP(hid_t group_id, storeRepP qrep, entitySet dom, int offset, MPI_Comm comm);
    void write_storeS(hid_t group_id, storeRepP qrep, entitySet dom, int offset, MPI_Comm comm);
    void read_storeS(hid_t group_id, storeRepP qrep, int &offset, MPI_Comm comm) ;
    void read_storeP(hid_t group_id, storeRepP qrep, int &offset, MPI_Comm comm);
  }
  
  extern bool use_parallel_io ;

  //-----------------------------------------------------------------------
 
  void redistribute_write_container(hid_t file_id, std::string vname,
                                    Loci::storeRepP var, fact_db &facts) ;
 
  void read_container_redistribute(hid_t file_id, std::string vname,
                                   Loci::storeRepP var, entitySet read_set,
                                   fact_db &facts) ;
  
  //-----------------------------------------------------------------------
 
  hid_t hdf5CreateFile(const char *name, unsigned flags, hid_t create_id, hid_t access_id, MPI_Comm comm, size_t file_size_estimate=0);
 

#ifndef LOCI_STRICT_COMM
  inline hid_t hdf5CreateFile(const char *name, unsigned flags, hid_t create_id, hid_t access_id, size_t file_size_estimate = 0) {
    return hdf5CreateFile(name,flags,create_id,access_id, MPI_COMM_WORLD,file_size_estimate) ;
  }    
#endif

  //-----------------------------------------------------------------------
  hid_t hdf5OpenFile(const char *name, unsigned flags, hid_t access_id,
		     MPI_Comm comm) ;


  //-----------------------------------------------------------------------
#ifndef LOCI_STRICT_COMM
  inline hid_t hdf5OpenFile(const char *name, unsigned flags, hid_t access_id) {
    return hdf5OpenFile(name,flags,access_id,MPI_COMM_WORLD) ;
  }
#endif

  //-----------------------------------------------------------------------
#ifndef LOCI_STRICT_COMM
  inline herr_t hdf5CloseFile(hid_t file_id) {
    if(use_parallel_io || Loci::MPI_rank==0)
      return H5Fclose(file_id) ;
    else
      return 0 ;
  }
#endif

  //-----------------------------------------------------------------------
  inline herr_t hdf5CloseFile(hid_t file_id, MPI_Comm comm) {
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    if(use_parallel_io || rank==0)
      return H5Fclose(file_id) ;
    else
      return 0 ;
  }

  //-----------------------------------------------------------------------
  //general way to open a file for writing
  hid_t writeVOGOpen(std::string filename, MPI_Comm comm LOCI_DEFAULT_COMM) ;
  hid_t readVOGOpen(std::string filename, MPI_Comm comm LOCI_DEFAULT_COMM) ;
  void writeVOGClose(hid_t file_id) ;
  //-----------------------------------------------------------------------  
  hid_t createUnorderedFile(const char * filename, entitySet set,
                            fact_db &facts) ;
  
  //-----------------------------------------------------------------------  
  inline hid_t createUnorderedFile(const char * filename, entitySet set) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::createUnorderedFile()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
      return -1 ;
    } else
      return createUnorderedFile(filename, set, *Loci::exec_current_fact_db) ;
  }
  
  //-----------------------------------------------------------------------  
  void closeUnorderedFile(hid_t file_id, MPI_Comm comm LOCI_DEFAULT_COMM) ;
  




  
   
  //-----------------------------------------------------------------------
  inline  void read_store(hid_t group_id, storeRepP qrep, int &offset, MPI_Comm comm){
    if(use_parallel_io)  pio::read_storeP(group_id, qrep, offset,  comm);
    else pio::read_storeS(group_id, qrep, offset,  comm);
  }

  //-----------------------------------------------------------------------
  inline void write_store(hid_t group_id, storeRepP qrep, entitySet dom, int offset, MPI_Comm comm){
    if(use_parallel_io)
      pio::write_storeP(group_id,qrep,dom,offset,comm) ;
    else
      pio::write_storeS(group_id,qrep,dom,offset,comm) ;
  }

  
  //-----------------------------------------------------------------------
  inline void writeContainer(hid_t file_id,std::string vname, Loci::storeRepP var, fact_db &facts) {
    redistribute_write_container(file_id,vname,var,facts) ;
  }
 
  //-----------------------------------------------------------------------
  inline void readContainer(hid_t file_id, std::string vname, Loci::storeRepP var, entitySet readSet, fact_db &facts) {
    read_container_redistribute(file_id,vname,var,readSet, facts) ;
  }

  //-----------------------------------------------------------------------
  inline void writeContainer(hid_t file_id,std::string vname, Loci::storeRepP var) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::writeContainer()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
    } else
      redistribute_write_container(file_id,vname,var,
                                   *Loci::exec_current_fact_db) ;
  }
 
  //-----------------------------------------------------------------------
  inline void readContainer(hid_t file_id, std::string vname, Loci::storeRepP var, entitySet readSet) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::readContainer()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
    } else
     
      read_container_redistribute(file_id,vname,var,readSet,
                                  *Loci::exec_current_fact_db) ;
  }
 
  void read_parameter(hid_t group_id, storeRepP qrep, MPI_Comm comm) ;
  //-----------------------------------------------------------------------
  void writeContainerRAW(hid_t file_id, std::string vname,
                         storeRepP var, MPI_Comm comm) ;

  //----------------------------------------------------------------------- 
  void readContainerRAW(hid_t file_id, std::string vname,
                        storeRepP var, MPI_Comm comm ) ;


  /*
    
    this function is when MPI_processes == 1 or when MPI_rank == 0
    if MPI_rank == 0, the caller will broadcast the data
    since no structure modification to the file, parallel version is not needed
  */
  template<class T> void readVectorSerial(hid_t group_id,
                                          const char *element_name,
                                          std::vector<T> &v ) {
    hsize_t    dimension;
    hid_t      dataset, dataspace;

    H5Eset_auto (H5E_DEFAULT,NULL, NULL);
    dataset  = H5Dopen(group_id, element_name,H5P_DEFAULT);
    if( dataset > 0 ) {
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);
      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;
      
      hid_t datatype = dp->get_hdf5_type() ;
      
      {
        std::vector<T> tmp(dimension) ;
        v.swap(tmp) ;
      }
      
      H5Dread(dataset, datatype, H5S_ALL, dataspace, H5P_DEFAULT, &v[0]);
      
      H5Sclose(dataspace);
      H5Dclose(dataset);
      H5Tclose(datatype) ;
    } else {
      std::vector<T> tmp ;
      v.swap(tmp) ;
    }
  }

  
  namespace pio {
    /*
      namespace pio is for functions that have separate serial io version and parallel io version.
      these low-level functions are generally used by Loci-developers instead of application-developers
    */


    //serial io version, called only if MPI_processes == 1 or MPI_rank==0
    template<class T> void writeVectorSerialS(hid_t group_id,
                                              const char *element_name,
                                              std::vector<T> &v ) {
      hsize_t array_size_combined = v.size() ;
      if(array_size_combined == 0)
        return ;

      int rank = 1 ;
      hsize_t dimension = array_size_combined ;
      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;

      hsize_t start = 0 ;
      hsize_t stride = 1 ;
      hsize_t count = v.size() ;
      hid_t datatype = dp->get_hdf5_type() ;

      hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                                dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;

      if(count != 0) {

        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            &start, &stride, &count, NULL) ;

        hid_t memspace = H5Screate_simple(rank, &count, NULL) ;

        H5Dwrite(dataset,datatype,memspace,dataspace,
                 H5P_DEFAULT, &v[0]) ;

        H5Sclose(memspace) ;

      }
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      H5Tclose(datatype) ;
      return;
    }

    //Since file_id and group_id is open for prime_comm, even just one process
    //write out data in serial, the dataset should be created by all processes in
    //prim_comm, otherwise the file won't be able to close properly
    template<class T> void writeVectorSerialP(hid_t group_id,
                                              const char *element_name,
                                              std::vector<T> &v ,
                                              MPI_Comm prime_comm) {
#ifndef H5_HAVE_PARALLEL
      writeVectorSerialS(group_id,element_name,v,prime_comm) ;
#else

      hsize_t array_size_combined = v.size() ;

      if(array_size_combined == 0)
        return ;

      if(MPI_COMM_NULL != prime_comm){
        int rank = 1 ;
        hsize_t dimension = array_size_combined ;
        hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

        typedef data_schema_traits<T> traits_type ;
        DatatypeP dp = traits_type::get_type() ;

        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = v.size() ;
        hid_t datatype = dp->get_hdf5_type() ;

        hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                                  dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;

        int prime_rank = -1;

        MPI_Comm_rank(prime_comm, &prime_rank);
        if(prime_rank==0){
          if(count != 0) {
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                &start, &stride, &count, NULL) ;
            hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
            H5Dwrite(dataset,datatype,memspace,dataspace,
                     H5P_DEFAULT, &v[0]) ;
            H5Sclose(memspace) ;
          }
        }
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
        H5Tclose(datatype) ;
      }
      return;
#endif
    }















    
    //-----------------------------------------------------------------------  
    template<class T> void writeUnorderedVectorS(hid_t group_id,
                                                 const char *element_name,
                                                 std::vector<T> &v,
                                                 MPI_Comm comm) {
      int my_rank = 0 ;
      MPI_Comm_rank(comm,&my_rank) ;
      int procs = 1 ;
      MPI_Comm_size(comm,&procs) ;

      //serial version
      if(procs==1){
      
        hsize_t array_size_combined = v.size() ;
        if(array_size_combined == 0)
          return ;
     
        int rank = 1 ;
        hsize_t dimension = array_size_combined ;
        hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
     
        typedef data_schema_traits<T> traits_type ;
        DatatypeP dp = traits_type::get_type() ;
      
        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = v.size() ;
        hid_t datatype = dp->get_hdf5_type() ;
     
        hid_t dataset = H5Dcreate(group_id,element_name,datatype,dataspace,
				  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
     
        if(count != 0) {
     
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start, &stride, &count, NULL) ;
     
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
     
          H5Dwrite(dataset,datatype,memspace,dataspace,
                   H5P_DEFAULT, &v[0]) ;
     
          H5Sclose(memspace) ;
     
        }
     
      
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
        H5Tclose(datatype) ;
        return;
      }

      size_t local_size = v.size() ;
      std::vector<size_t> recv_sizes(procs) ;
      MPI_Gather(&local_size,sizeof(size_t),MPI_BYTE,
                 &recv_sizes[0],sizeof(size_t),MPI_BYTE,0,comm) ;

      if(my_rank == 0) {
        hsize_t array_size = 0 ;
        for(int i=0;i<procs;++i)
          array_size += recv_sizes[i] ;
        if(array_size == 0)
          return ;
        int rank = 1 ;
        hsize_t dimension = array_size ;

        hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

        typedef data_schema_traits<T> traits_type ;
        DatatypeP dp = traits_type::get_type() ;

        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = recv_sizes[0] ;
        hid_t datatype = dp->get_hdf5_type() ;

        hid_t dataset = H5Dcreate(group_id,element_name,datatype,dataspace,
				  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;

        if(count != 0) {
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,datatype,memspace,dataspace,
                   H5P_DEFAULT, &v[0]) ;
          H5Sclose(memspace) ;
        }
        for(int i=1;i<procs;++i) {
          start += recv_sizes[i-1] ;
          if(recv_sizes[i] == 0)
            continue ;
          int flag = 0 ;
          MPI_Send(&flag,1,MPI_INT,i,0,comm) ;
          std::vector<T> rv(recv_sizes[i]) ;
          MPI_Status mstat ;
          MPI_Recv(&rv[0],sizeof(T)*recv_sizes[i],MPI_BYTE,i,0,comm,
                   &mstat) ;
          count = recv_sizes[i] ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;
          H5Dwrite(dataset,datatype,memspace,dataspace,
                   H5P_DEFAULT, &rv[0]) ;
          H5Sclose(memspace) ;
        }

        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
        H5Tclose(datatype) ;
      } else {
        if(local_size == 0)
          return ;

        int flag = 0;
        MPI_Status mstat ;
        MPI_Recv(&flag,1,MPI_INT,0,0,comm,&mstat) ;
        MPI_Send(&v[0],sizeof(T)*local_size,MPI_BYTE,0,0,comm) ;
      }
    }
 
    //-----------------------------------------------------------------------  
#ifndef LOCI_STRICT_COMM
    template<class T> void writeUnorderedVectorS(hid_t group_id,
                                                 const char *element_name,
                                                 std::vector<T> &v) {
      writeUnorderedVectorS(group_id,element_name,v,MPI_COMM_WORLD) ;
    }
#endif
  
    //-----------------------------------------------------------------------  
    template<class T> void writeUnorderedVectorP(hid_t group_id,
                                                 const char *element_name,
                                                 std::vector<T> &v,
                                                 MPI_Comm prime_comm) {
      /*
        in serial version,  process 0 gathers all the data and write it out
        here in parallel version, each process performs parallel writing directly
        this function is called when use_parallel_io is true
      */
      int mpi_size;
      int mpi_rank;
      MPI_Comm_rank(prime_comm, &mpi_rank);
      MPI_Comm_size(prime_comm, &mpi_size);

      //serial version
      if(mpi_size==1){//this code probably never used
        hsize_t array_size_combined = v.size() ;
        if(array_size_combined == 0)
          return ;

        int rank = 1 ;
        hsize_t dimension = array_size_combined ;
        hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;

        typedef data_schema_traits<T> traits_type ;
        DatatypeP dp = traits_type::get_type() ;

        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = v.size() ;
        hid_t datatype = dp->get_hdf5_type() ;

        hid_t dataset = H5Dcreate(group_id,element_name,datatype,
                                  dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;

        if(count != 0) {
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start, &stride, &count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;

          H5Dwrite(dataset,datatype,memspace,dataspace,H5P_DEFAULT, &v[0]) ;
          H5Sclose(memspace) ;
        }
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
        H5Tclose(datatype) ;
        return;
      }

      //    Loci::stopWatch s;
      //    s.start();
    
      //each process figure out the size and the start position to write
      int local_size = v.size() ; //my size to write
      hsize_t rsize = local_size; //individual size for each process 
      std::vector<hsize_t> prime_count(mpi_size);
      MPI_Allgather(&rsize,sizeof(hsize_t),MPI_BYTE,
                    &prime_count[0],sizeof(hsize_t),MPI_BYTE,prime_comm) ;

      std::vector<hsize_t> pdispls(mpi_size) ; //the start point of each process 
      pdispls[0] = 0 ;
      for(int i = 1; i < mpi_size; i++) {
        pdispls[i] = pdispls[i-1]+prime_count[i-1] ;
      }
      hsize_t array_size = pdispls[mpi_size-1]+prime_count[mpi_size-1] ; //size of the whole dataset

      if(array_size == 0)
        return ;


      //create dataset collectively
      int rank = 1 ;
      hsize_t dimension = array_size ;
      hid_t dataspace = H5Screate_simple(rank,&dimension,NULL) ;
      typedef data_schema_traits<T> traits_type ;
      DatatypeP dp = traits_type::get_type() ;
      hsize_t start = pdispls[mpi_rank] ;
      hsize_t stride = 1 ;
     
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t dataset = H5Dcreate(group_id,element_name,datatype,
				dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
      //choose a hyperslab
      herr_t ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                       &start,&stride,&rsize, NULL) ;
      if(ret<0) {
        cerr << "H5Sselect_hyperslab failed in writeUnorderedVector" << endl ;
      }
      WARN(ret < 0) ;

      // create a memory dataspace  
      hid_t memspace = H5Screate_simple (rank, &rsize, NULL);
      WARN(memspace < 0) ;

      //write data
      hid_t xfer_plist = create_xfer_plist(hdf5_const::dxfer_coll_type);
      H5Dwrite(dataset,datatype,memspace,dataspace, xfer_plist,  &v[0]) ;

      //close everything
      H5Pclose(xfer_plist) ;
      H5Sclose(memspace) ;
      H5Sclose(dataspace) ;
      H5Tclose(datatype) ;
      H5Dclose(dataset) ;

    }

    //-----------------------------------------------------------------------  
#ifndef LOCI_STRICT_COMM
    template<class T> void writeUnorderedVectorP(hid_t group_id,
                                                 const char *element_name,
                                                 std::vector<T> &v
                                                 ) {
      writeUnorderedVectorP(group_id,
                            element_name,
                            v,
                            MPI_COMM_WORLD );
    }
#endif



    template<class T> void readUnorderedVectorS(hid_t group_id,
                                                const char *element_name,
                                                std::vector<T> &v,
                                                MPI_Comm comm) {
      int my_rank = 0 ;
      MPI_Comm_rank(comm,&my_rank) ;
      int procs = 1 ;
      MPI_Comm_size(comm,&procs) ;

      if(procs == 1)
        return readVectorSerial(group_id,element_name,v) ;


      hsize_t    dimension=0;
      hid_t      dataset=0, dataspace=0 ;

      if(my_rank == 0) {
        H5Eset_auto (H5E_DEFAULT,NULL, NULL);
        dataset  = H5Dopen(group_id, element_name,H5P_DEFAULT);
      }
      int check = 0 ;
      if(my_rank == 0 && dataset > 0)
        check = 1 ;
      MPI_Bcast(&check,1,MPI_INT,0,comm) ;

      if( check ) {
        if(my_rank == 0) {
          dataspace  = H5Dget_space(dataset);
          H5Sget_simple_extent_dims (dataspace, &dimension, NULL);
        }
      }
      int array_size = dimension ;
      MPI_Bcast(&array_size,1,MPI_INT,0,comm) ;
      int delta = array_size/procs ;
      int rem = array_size%procs ;
      int loc_array_size = delta + ((my_rank<rem)?1:0);
      {
        std::vector<T> tmp(loc_array_size) ;
        v.swap(tmp) ;
      }

      if( check ) {
        if(my_rank == 0) {
          typedef data_schema_traits<T> traits_type ;
          DatatypeP dp = traits_type::get_type() ;

          hid_t datatype = dp->get_hdf5_type() ;
          int rank = 1 ;
          hsize_t start = 0 ;
          hsize_t stride = 1 ;
          hsize_t count = loc_array_size ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              &start,&stride,&count, NULL) ;
          hid_t memspace = H5Screate_simple(rank, &count, NULL) ;

          H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &v[0]);
          H5Sclose(memspace) ;
          start += count ;
          for(int i=1;i<procs;++i) {
            count = delta + ((i<rem)?1:0) ;
            std::vector<T> tmp(count) ;
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                &start,&stride,&count, NULL) ;
            memspace = H5Screate_simple(rank, &count, NULL) ;
            H5Dread(dataset,datatype,memspace,dataspace,H5P_DEFAULT,&tmp[0]) ;
            MPI_Send(&tmp[0],count*sizeof(T),MPI_BYTE,i,0,comm) ;
            H5Sclose(memspace) ;
            start += count ;
          }
          H5Tclose(datatype) ;
        } else {
          MPI_Status mstat ;
          MPI_Recv(&v[0],loc_array_size*sizeof(T),MPI_BYTE,0,0,comm,&mstat) ;
        }
      } else {
        std::vector<T> tmp ;
        v.swap(tmp) ;
      }
      if(my_rank == 0) {
        H5Sclose(dataspace);
        H5Dclose(dataset);
      }
    }

    template<class T> void readUnorderedVectorP(hid_t group_id,
                                                const char *element_name,
                                                std::vector<T> &v,
                                                MPI_Comm prime_comm
                                                ) {


#ifndef H5_HAVE_PARALLEL
      readUnorderedVectorS(group_id,element_name,v,prime_comm) ;
#else
      int procs = 1 ;
      MPI_Comm_size(MPI_COMM_WORLD,&procs) ;

      if(procs == 1)
        return readVectorSerial(group_id,element_name,v) ;


      int mpi_rank = -1;
      int mpi_size = -1;
      hsize_t    dimension=0;
      hid_t      dataset=0, dataspace=0 ;
      if (MPI_COMM_NULL != prime_comm){
        MPI_Comm_rank(prime_comm, &mpi_rank);
        MPI_Comm_size(prime_comm, &mpi_size);
        H5Eset_auto (H5E_DEFAULT,NULL, NULL);
        dataset  = H5Dopen(group_id, element_name,H5P_DEFAULT);
        WARN(dataset<0) ;

        dataspace  = H5Dget_space(dataset);
        H5Sget_simple_extent_dims (dataspace, &dimension, NULL);

        int array_size = dimension ;
        int delta = array_size/mpi_size ;
        int rem = array_size%mpi_size ;
        int loc_array_size = delta + ((mpi_rank<rem)?1:0);

        {
          std::vector<T> tmp(loc_array_size) ;
          v.swap(tmp) ;
        }

        //parallel read here
        typedef data_schema_traits<T> traits_type ;
        DatatypeP dp = traits_type::get_type() ;
        hid_t datatype = dp->get_hdf5_type() ;
        int rank = 1 ;
        hsize_t start = (mpi_rank< rem)?mpi_rank*(delta+1):(delta+1)*rem+delta*(mpi_rank-rem) ;
        hsize_t stride = 1 ;
        hsize_t count = loc_array_size ;

        herr_t ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                         &start,&stride,&count, NULL) ;
        if(ret<0) {
          std::cerr << "warning H5Sselect_hyperslab failed." << std::endl;
        }
        WARN(ret<0) ;
        hid_t memspace = H5Screate_simple(rank, &count, NULL);
        WARN(memspace<0) ;

        H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &v[0]);

        H5Sclose(memspace) ;
        H5Sclose(dataspace) ;
        H5Tclose(datatype) ;
        H5Dclose(dataset) ;
      }
#endif
    }

    
  }


  //-----------------------------------------------------------------------  
  /*this is the general version for application developers*/ 
  template<class T> void writeUnorderedVector(hid_t group_id,
                                              const char *element_name,
                                              std::vector<T> &v,
                                              MPI_Comm prime_comm) {
   
    if(use_parallel_io) pio::writeUnorderedVectorP( group_id,
                                                    element_name,
                                                    v,
                                                    prime_comm);
    else pio::writeUnorderedVectorS( group_id,
                                     element_name,
                                     v,
                                     prime_comm);
  }


  //-----------------------------------------------------------------------  
#ifndef LOCI_STRICT_COMM
  template<class T> void writeUnorderedVector(hid_t group_id,
                                              const char *element_name,
                                              std::vector<T> &v
                                              ) {
    writeUnorderedVector(group_id,
                         element_name,
                         v,
                         MPI_COMM_WORLD );
  }
#endif
  
  //-----------------------------------------------------------------------  
  template<class T> void writeUnorderedStore(hid_t file_id,
                                             const_store<T> &s, entitySet set,
                                             const char *name,
                                             MPI_Comm comm LOCI_DEFAULT_COMM) {
    std::vector<T> v(set.size()) ;
    size_t c = 0 ;
    FORALL(set,ii) {
      v[c++] = s[ii] ;
    } ENDFORALL ;
    writeUnorderedVector(file_id,name,v,comm) ;
  }
  
  //-----------------------------------------------------------------------  
  void writeSetIds(hid_t file_id, entitySet local_set, fact_db &facts) ;
  
 
  //-----------------------------------------------------------------------  
  void parallelWriteGridTopology(const char *filename,
                                 storeRepP upperRep,
                                 storeRepP lowerRep,
                                 storeRepP boundary_mapRep,
                                 storeRepP face2nodeRep,
                                 storeRepP refRep,
                                 storeRepP bnamesRep,
                                 storeRepP posRep,
                                 entitySet localCells,
                                 fact_db &facts) ;

  //-----------------------------------------------------------------------  
  inline
  void parallelWriteGridTopology(const char *filename,
                                 storeRepP upperRep,
                                 storeRepP lowerRep,
                                 storeRepP boundary_mapRep,
                                 storeRepP face2nodeRep,
                                 storeRepP refRep,
                                 storeRepP bnamesRep,
                                 storeRepP posRep,
                                 entitySet localCells) {
    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::parallelWriteGridTopology()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
    } else
      parallelWriteGridTopology(filename, upperRep, lowerRep, boundary_mapRep,
                                face2nodeRep, refRep, bnamesRep, posRep,
                                localCells, *Loci::exec_current_fact_db) ;
  }
  
  
  //open /output/$bc_name/$file_name
  hid_t open_boundary_file(std::string bc_name,
                           std::string file_name
                           );
 
  //get boundary faces that belong to a boundary surface current_bc
  entitySet get_boundary_faces(std::string current_bc,//boundary name
                               storeRepP refRep, // ref map
                               storeRepP bnamesRep,//bounadry name store
                               entitySet fset //all boundary faces 
                               );

  //get boundary nodes that belong to a boundary surface current_bc
  entitySet get_boundary_nodes(std::string current_bc,//boundary name
                               storeRepP face2nodeRep,
                               storeRepP refRep, // ref map
                               storeRepP bnamesRep,//bounadry name store
                               entitySet fset, //all boundary faces 
                               fact_db &facts );

  void writeBoundaryTopo(hid_t file_id, //file_id of this boudnary surface
                         storeRepP face2nodeRep,
                         entitySet bfaces, //boundary faces belong to this surface 
                         fact_db &facts ); 
  


 
  //return value at cut position 
  template<class T>  inline T interpolate_val( double t, //weight of first node
                                               T a, //value at first node
                                               T b //value at second node
                                               ){
    
    return t*a + (1.0 - t)*b ;
  }
 
  
 
  
  struct CutPlane {
    storeRepP edgesWeight; //it's a store<double> allocated on edgesCut, containing the weight for interpoplation for each edge in edgesCut
    std::vector<std::vector<int > > faceLoops;  //loops formed, the values stored are edge ids, which is either local edge entity or negated index to inner_edges
    
    CutPlane( storeRepP ew, 
              std::vector<std::vector<int > >& fl){
      edgesWeight = ew;
      faceLoops = fl;
    }
    
    CutPlane(const CutPlane& cp){
      edgesWeight = cp.edgesWeight ;
      faceLoops =  cp.faceLoops;
    }
    
    CutPlane(){}
    
  };
  
  CutPlane getCutPlane(storeRepP upperRep,
                       storeRepP lowerRep,
                       storeRepP boundary_mapRep,
                       storeRepP face2nodeRep,
                       storeRepP face2edgeRep,
                       storeRepP edge2nodeRep,
                       storeRepP posRep,
                       entitySet localCells,//all geom_cells
                       fact_db &facts);

  inline
  CutPlane getCutPlane(storeRepP upperRep,
                       storeRepP lowerRep,
                       storeRepP boundary_mapRep,
                       storeRepP face2nodeRep,
                       storeRepP face2edgeRep,
                       storeRepP edge2nodeRep,
                       storeRepP posRep,
                       entitySet localCells//all geom_cells
                       ) {

    if(Loci::exec_current_fact_db == 0) {
      std::cerr << "Loci::getCutPlane()" ;
      std::cerr << "this routine needs a fact database argument when called outside of a rule!" << endl ;
      Loci::Abort() ;
      return CutPlane() ;
    } else return getCutPlane( upperRep,
                               lowerRep,
                               boundary_mapRep,
                               face2nodeRep,
                               face2edgeRep,
                               edge2nodeRep,
                               posRep,
                               localCells,//all geom_cells
                               *Loci::exec_current_fact_db);
  }

  // Convert container from local numbering to output file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // fact_db pointer  (facts)
  // MPI Communicator
  storeRepP Local2FileOrder_output(storeRepP sp, entitySet dom,
                                   fact_db& facts, MPI_Comm comm);

  //serial/parallel io
  template<class T>   void writeCutPlaneNodalVal(hid_t file_id,
                                                 std::string element_name,
                                                 storeRepP face2nodeRep,
                                                 storeRepP edge2nodeRep,
                                                 const_store<T> & pos,
                                                 const Loci::CutPlane &cp,
                                                 fact_db &facts){
    
   
    
    const_multiMap face2node(face2nodeRep) ;
    const_MapVec<2> edge2node(edge2nodeRep);
    const_store<double> edgesWeight(cp.edgesWeight); //the weight for interpoplation for each edgesCut, allocated on edgesCut
    entitySet edgesCut = edgesWeight.domain();
        
    
    //check the domain
    if((edgesCut-edge2node.domain())!=EMPTY){
      debugout<< "ERROR: the domain of edge2node is smaller than cp.edgesCut"<<endl;
    }
    
    //compute the cutting positions of edges 
    store<T> edge_pos;
    edge_pos.allocate(edgesCut);
    FORALL(edgesCut, e){
      double w =edgesWeight[e];
      T a = pos[edge2node[e][0]];
      T b = pos[edge2node[e][1]];
      T p = interpolate_val(w, a, b);
      edge_pos[e] = p;
    }ENDFORALL;
    
   
    //transform the store into output order
    store<T> gedge_pos;
    storeRepP geposRep =  Local2FileOrder_output(edge_pos.Rep(),  edgesCut, 
                                                 facts, MPI_COMM_WORLD);
       
    if(geposRep == NULL){
      gedge_pos .allocate(EMPTY);
    }else{
      gedge_pos = geposRep;
    }
   
   
    //get positions std::vector
    entitySet local_edges_cut = gedge_pos.domain();
    int num_edge_nodes = local_edges_cut.size();
    std::vector<T>  vpos(num_edge_nodes);
    int cnt = 0 ;
    entitySet::const_iterator ei ;
    for(ei=local_edges_cut.begin();ei!=local_edges_cut.end();++ei)
      vpos[cnt++] = gedge_pos[*ei];
      
    //write out the vector
    writeUnorderedVector(file_id, element_name.c_str(), vpos) ;
  }
  
 

  template<class T>   void writeCutPlaneNodalVal(hid_t file_id,
                                                 std::string element_name,
                                                 storeRepP face2nodeRep,
                                                 storeRepP edge2nodeRep,
                                                 const_store<T> & pos,
                                                 const Loci::CutPlane &cp
                                                 ){
    
    writeCutPlaneNodalVal( file_id,
                           element_name,
                           face2nodeRep,
                           edge2nodeRep,
                           pos,
                           cp,
                           *Loci::exec_current_fact_db );
  }
  
  
                          
  void writeCutPlaneTopo(hid_t bc_id,
                         const CutPlane& cp,
                         fact_db &facts) ;

 
   
  // Updated container communication code
  class partitionFunctionType: public CPTR_type  {
  public:
    virtual int numRanks() const = 0 ;
    virtual void mapKeyToRank(int keyRank[],
			      const Entity inputKeys[], size_t sz) const = 0 ;
  } ;

  class algorithmicPartition : public partitionFunctionType {
    Entity mn, delta ;
    int nRanks ;
  public:
    algorithmicPartition(Entity mnl, Entity mxl, int nRanksl) {
      mn = mnl ;
      nRanks = nRanksl ;
      int sz = mxl-mnl+1 ;
      delta = (sz+nRanks-1)/nRanks ;
    }
    int numRanks() const { return nRanks; }
    
    void mapKeyToRank(int keyRank[], const Entity inputKeys[], size_t sz) const {
      for(size_t i=0;i<sz;++i) {
	keyRank[i] = max(min(int((inputKeys[i]-mn)/delta),nRanks-1),0) ;
      }
    }
  } ;

  class generalPartition: public partitionFunctionType {
    std::vector<std::pair<interval, int> > splits ;
    int nRanks ;
  public:
    generalPartition(const std::vector<Entity> &splits_in) {
      nRanks = splits_in.size()+1 ;
      std::vector<std::pair<interval, int> > tsplit(splits_in.size()+1) ;
      splits.swap(tsplit) ;
      int cx = std::numeric_limits<Entity>::min() ;
      for(size_t i=0;i<splits_in.size();++i) {
	splits[i].first.first = cx ;
	splits[i].first.second = splits_in[i] ;
	splits[i].second = i ;
	cx = splits_in[i]+1 ;
      }
      splits[splits_in.size()].first.first = cx ;
      splits[splits_in.size()].first.second = std::numeric_limits<Entity>::max() ;
      splits[splits_in.size()].second = splits_in.size() ;
    }
    generalPartition(const std::vector<entitySet> &init_ptn) {
      entitySet totset ;
      nRanks = init_ptn.size() ;
      for(size_t i=0;i<init_ptn.size();++i) {
	totset += init_ptn[i] ;
	for(size_t j=0;j<init_ptn[i].num_intervals();++j) {
	  splits.push_back(std::pair<interval,int>(init_ptn[i][j],i)) ;
	}
      }
      entitySet negspace = ~totset ;
      for(size_t j=0;j<negspace.num_intervals();++j) {
	splits.push_back(std::pair<interval,int>(negspace[j],-1)) ;
      }
      std::sort(splits.begin(),splits.end()) ;
      for(size_t i=0;i<splits.size()-1;++i)
	if(splits[i].first.second+1 != splits[i+1].first.first) {
	  cerr << "set array input does not form partition" << endl ;
	}
    }

    int numRanks() const { return nRanks; }

    void mapKeyToRank(int keyRank[], const Entity inputKeys[],
		      size_t sz) const {
      int lastr = 0 ;
      for(size_t i=0;i<sz;++i) {
	Entity key = inputKeys[i] ;
	int low = 0, high = splits.size()-1 ;
	while(key < splits[lastr].first.first ||
	      key > splits[lastr].first.second) {
	  if(key<splits[lastr].first.first)  
	    high = lastr-1 ; // lastr is to high
	  else  
	    low = lastr+1 ;  // lastr is to low
	  lastr = (low+high)/2 ;
	}
	keyRank[i] = splits[lastr].second ;
      }

    }
  } ;
  void generalMPIComm(Loci::storeRepP op,
		      Loci::storeRepP sp,
		      const std::vector<Loci::entitySet> &sendSets,
		      const std::vector<Loci::sequence> &recvSeqs,
		      MPI_Comm comm) ;
  storeRepP
  generalCommStore(// input store
		   Loci::storeRepP sp,
		   // first: from entity (in container ordering),
		   // second: to global partitioned entity map
		   const std::vector<std::pair<Entity,Entity> > &commMap,
		   // To entity partition
		   Loci::CPTR<Loci::partitionFunctionType> partition,
		   // mapping from global number to local numbering
		   const std::vector<std::pair<Entity,Entity> > &global2local,
		   // If this is null, create new container, otherwise
		   // assume it is allocated already
		   Loci::storeRepP op,
		   MPI_Comm comm) ;
  /// Helper Function to Redistribute store to a new numbering with a partition
  /// specified by the provided splits
  /// @param[out] op is the output store pointer which will hold the
  ///                redistributed store
  /// @param[in] sp is the input store that will be redistributed
  /// @param[in] dom is the domain of the input store that will be communicated
  /// @param[in] toNumbering is the mapping to the new numbering that will
  ///            be used to number the redistributed store
  /// @param[in] splits is an array of values that define how the renumbered
  ///            store will be distributed across processors
  /// @param[in] comm is the MPI communicator of the processes participating in
  ///            the redistribution
  void redistributeStore(Loci::storeRepP &op,
                         Loci::storeRepP sp,
                         entitySet dom,
                         const Map &toNumbering,
                         const std::vector<int> &splits,
                         MPI_Comm comm) ;

  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm) ;

  int getKeyDomain(entitySet dom, fact_db::distribute_infoP dist, MPI_Comm comm) ;
  entitySet
  getF2G(Map &f2g, Loci::entitySet fdom, dMap &g2f, MPI_Comm comm) ;
  entitySet
  getF2G2(Map &f2g, entitySet fdom, Map &l2g, Map &l2f, MPI_Comm comm) ;
  void File2LocalOrderGeneral(storeRepP &result, entitySet resultSet,
                              storeRepP input, int offset,
                              fact_db::distribute_infoP dist,
                              MPI_Comm comm) ;
  void getL2FMap(Map &l2f, entitySet dom, fact_db::distribute_infoP dist, MPI_Comm comm LOCI_DEFAULT_COMM) ;
  void FindSimpleDistribution(entitySet dom, const Map &l2f,
                              std::vector<int> &splits, MPI_Comm comm) ;
  void memoryBalancedDistribution(std::vector<int> &splits_out,
                                  const store<int> &countl,
                                  entitySet dom,
                                  const Map &toNumbering,
                                  MPI_Comm comm) ;
  storeRepP gatherStore(// Input Store
                        storeRepP sp,
                        // EntitySet of input to reorder
                        const std::vector<int> &commPattern,
                        // Splits for partition
                        const std::vector<int> &splits,
                        MPI_Comm comm) ;
  storeRepP gatherMultiStore(// Input Store
                             storeRepP sp,
                             // EntitySet of input to reorder
                             const std::vector<int> &commPattern,
                             // Splits for partition
                             const std::vector<int> &splits,
                             MPI_Comm comm) ;




}
#endif
