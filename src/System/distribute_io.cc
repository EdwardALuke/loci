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
//#define io_performance
#include <vector>
using std::vector;
#include <string>
using std::string;

#include <iostream>
using std::cerr;
using std::endl;
using std::ostringstream;

#include <algorithm>
using std::sort;
using std::pair ;

#include <mpi.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include "dist_tools.h"
#include <fact_db.h>
#include <constraint.h>
#include <execute.h>

using std::cout;

namespace Loci {
  
  extern string PFS_Script ;
  extern bool use_parallel_io ;
  void read_parameter(hid_t group_id, storeRepP qrep, MPI_Comm comm);
  entitySet BcastEntitySet(entitySet set, int root, MPI_Comm comm) {

    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    // Now lets share the domain with all other processors ;
    int sz = set.num_intervals() ;
    MPI_Bcast(&sz,1,MPI_INT,root,comm) ;
    if(sz == 0)
      return EMPTY ;
    vector<interval> vlist(sz) ;
    if(rank == 0) {
      for(int i=0;i<sz;++i)
        vlist[i] = set[i] ;
    }
    MPI_Bcast(&vlist[0],sz*2,MPI_INT,root,comm) ;
    set = EMPTY ;
    for(int i = 0;i<sz;++i)
      set += vlist[i] ;
    return set ;
  }

  vector<int> simplePartitionVec(int mn, int mx, int p) {
    vector<int> nums(p+1) ;
    int n = mx-mn+1 ;
    int dn = n/p ; // divisor
    int rn = n%p ; // remainder
    int start = mn ;
    nums[0] = start ;
    for(int i=0;i<p;++i) {
      start += dn+((i<rn)?1:0) ;
      nums[i+1] = start ;
    }
    FATAL(start != mx+1) ;
    return nums ;
  }

  vector<entitySet> simplePartition(int mn, int mx, MPI_Comm comm) {
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    vector<int> pl = simplePartitionVec(mn,mx,p) ;
    vector<entitySet> ptn(p) ;
    for(int i=0;i<p;++i)
      ptn[i] = interval(pl[i],pl[i+1]-1) ;
    return ptn ;
  }

  //serial/parallel io
  hid_t hdf5CreateFile(const char *name, unsigned flags, hid_t create_id, hid_t access_id, MPI_Comm comm, size_t file_size_estimate)
  {
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    if(use_parallel_io) {
      if(rank == 0 && file_size_estimate > 0 && PFS_Script != "") {
        ostringstream oss ;
        oss << PFS_Script << " " << name << " " << file_size_estimate ;
        string script = oss.str() ;
        int ret =system(script.c_str()) ;
        if(ret !=0)
          cerr << "Error, script '" << script << "' failed!"
               << endl ;
      }
      string filename = name ;
      hid_t file_id = -1;
      hid_t  acc_plist;
      // open collectively by all processor in MPI_COMM_WORLD,
      acc_plist = Loci::create_faccess_plist(MPI_COMM_WORLD,
                                             Loci::PHDF5_MPI_Info,
                                             Loci::hdf5_const::facc_type); 
      file_id = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,acc_plist) ;
      H5Pclose(acc_plist);
      if(file_id < 0) {
        if(rank==0) cerr << "unable to open file " << filename << endl ;
        Loci::Abort() ;
      }
      return file_id ;
    } else if(rank==0)
      return H5Fcreate(name,flags,create_id,access_id) ;
    else
      return 0 ;
  }

  
 

  hid_t writeVOGOpen(string filename) {
    if(use_parallel_io){    
      hid_t file_id = 0;
      hid_t  acc_plist;
      
      acc_plist = Loci::create_faccess_plist(MPI_COMM_WORLD,
                                             Loci::PHDF5_MPI_Info,
                                             Loci::hdf5_const::facc_type); 
      file_id = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,acc_plist) ;
      H5Pclose(acc_plist);
      if(file_id == 0) {
        if(MPI_rank==0) cerr << "unable to open file " << filename << endl ;
        Loci::Abort() ;
      }
      return file_id ;
    }else{
      hid_t file_id = 0 ;
      if(MPI_rank==0) 
        file_id = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;
      return file_id ;
    }
  }
 

  hid_t hdf5OpenFile(const char *name, unsigned flags, hid_t access_id,
                     MPI_Comm comm) {
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    if(!use_parallel_io) {
      if(rank == 0)
        return H5Fopen(name,flags,access_id) ;
      else
        return 0 ;
    }

    hid_t file_id = 0;
    hid_t  acc_plist;
    // open collectively by all processor in MPI_COMM_WORLD,
    acc_plist = Loci::create_faccess_plist(MPI_COMM_WORLD,
                                           Loci::PHDF5_MPI_Info,
                                           Loci::hdf5_const::facc_type); 
    file_id = H5Fopen(name,H5F_ACC_RDONLY, acc_plist) ;
    H5Pclose(acc_plist);

    return file_id ;
    
  }

  hid_t readVOGOpen(string filename) {
    
    if(use_parallel_io){
    
      hid_t file_id = -1;
    
      hid_t  acc_plist;
      // open collectively by all processor in MPI_COMM_WORLD,
      acc_plist = Loci::create_faccess_plist(MPI_COMM_WORLD,
                                             Loci::PHDF5_MPI_Info,
                                             Loci::hdf5_const::facc_type); 
      file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY, acc_plist) ;
      H5Pclose(acc_plist);
      if(file_id < 0) {
        if(MPI_rank==0) cerr << "unable to open file " << filename << endl ;
        Loci::Abort() ;
      }
      return file_id ;
    }else{
      hid_t file_id = 0 ;
      if(MPI_rank==0) 
        file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY, H5P_DEFAULT) ;
      return file_id ;
    }
    
  }


  namespace pio{
    void write_frame_info_paramS(hid_t group_id, frame_info &fi, MPI_Comm comm) {
      int prank = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      // Write out is_stat and vector size
      if(prank==0) {
        hsize_t dimension = 1 ;
        int rank = 1 ;
        hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
        hid_t datatype = H5T_NATIVE_INT ;

        hid_t dataset = H5Dcreate(group_id, "is_stat", datatype, dataspace,
				  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.is_stat) ;
        H5Dclose(dataset) ;

        dataset = H5Dcreate(group_id, "vec_size", datatype, dataspace,
			    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.size) ;
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
        if(fi.is_stat != 0) {
          rank = 1 ;
          dimension = 1 ;
          dataspace = H5Screate_simple(rank,&dimension,NULL) ;

          hid_t dataset = H5Dcreate(group_id,"second_level",H5T_NATIVE_INT,
                                    dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

          hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
          H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &fi.second_level[0]) ;
          H5Sclose(memspace) ;
          H5Dclose(dataset) ;
          H5Sclose(dataspace) ;
        }
      }
    }
  
  
    //H5DWrite is independent, only p0 performs it
    void write_frame_info_paramP(hid_t group_id, frame_info &fi, MPI_Comm comm) {
#ifndef H5_HAVE_PARALLEL
      write_frame_info_paramS(group_id, fi, comm);
#else
      int prank = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      // Write out is_stat and vector size
    
      //each process create dataset
    
      hsize_t dimension = 1 ;
      int rank = 1 ;
      hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
      hid_t datatype = H5T_NATIVE_INT ;

      hid_t dataset = H5Dcreate(group_id, "is_stat", datatype, dataspace,
				H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

      
      if(prank==0) {
        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.is_stat) ;
      }
      H5Dclose(dataset) ;

      dataset = H5Dcreate(group_id, "vec_size", datatype, dataspace,
			  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

      if(prank==0) {
        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.size) ;
      }
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
      if(fi.is_stat != 0) {
        rank = 1 ;
        dimension = 1 ;
        dataspace = H5Screate_simple(rank,&dimension,NULL) ;

        hid_t dataset = H5Dcreate(group_id,"second_level",H5T_NATIVE_INT,
                                  dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

        if(prank==0) {
          hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
          H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &fi.second_level[0]) ;
    
          H5Sclose(memspace) ;
        }
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
      }
  
#endif
    }

  
    void write_frame_infoS(hid_t group_id, frame_info &fi, MPI_Comm comm) {
      int prank = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      // Write out is_stat and vector size
      if(prank==0) {
        hsize_t dimension = 1 ;
        int rank = 1 ;
        hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
        hid_t datatype = H5T_NATIVE_INT ;

        hid_t dataset = H5Dcreate(group_id, "is_stat", datatype, dataspace,
				  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.is_stat) ;
        H5Dclose(dataset) ;

        dataset = H5Dcreate(group_id, "vec_size", datatype, dataspace,
			    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.size) ;
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
      }
      int is_stat = fi.is_stat ;
      int size = fi.size ;
      MPI_Bcast(&is_stat,1,MPI_INT,0,comm) ;
      MPI_Bcast(&size,1,MPI_INT,0,comm) ;
      if(size == 0) { // Two level framing
        write_vector_intS(group_id, "first_level", fi.first_level,comm) ;
        if(is_stat != 0) {
          write_vector_intS(group_id, "second_level", fi.second_level,comm) ;
        }
      } else {
        if(is_stat != 0) {
          write_vector_intS(group_id,"second_level",fi.second_level,comm) ;
        }
      }
    }

    //single int write independntly, only p0 perform writing
    //vector<int> write collectively/independently
    void write_frame_infoP(hid_t group_id, frame_info &fi, MPI_Comm comm) {
#ifndef H5_HAVE_PARALLEL
      write_frame_infoS(group_id, fi, comm);
#else
    
      int prank = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      // Write out is_stat and vector size

      //all process create/close, only process 0 write is_stat and vec_size
      hsize_t dimension = 1 ;
      int rank = 1 ;
      hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
      hid_t datatype = H5T_NATIVE_INT ;

      hid_t dataset = H5Dcreate(group_id, "is_stat", datatype, dataspace,
				H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

      if(prank==0) {
        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.is_stat) ;
      }
      H5Dclose(dataset) ;

      dataset = H5Dcreate(group_id, "vec_size", datatype, dataspace,
			  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

      if(prank==0) {
        H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fi.size) ;
      }
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
  

    
      int is_stat = fi.is_stat ;
      int size = fi.size ;
      MPI_Bcast(&is_stat,1,MPI_INT,0,comm) ;
      MPI_Bcast(&size,1,MPI_INT,0,comm) ;
      if(size == 0) { // Two level framing
        write_vector_intP(group_id, "first_level", fi.first_level,comm) ;
        if(is_stat != 0) {
          write_vector_intP(group_id, "second_level", fi.second_level,comm) ;
        }
      } else {
        if(is_stat != 0) {
          write_vector_intP(group_id,"second_level",fi.second_level,comm) ;
        }
      }
#endif
    }

  
    //only process 0 read in fi 
    frame_info read_frame_info_paramS(hid_t group_id,  int dom_size, MPI_Comm comm) {
      int prank = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      int is_stat = 0 ;
      int sz = 0 ;
      // Read in is_stat and vector size
      frame_info fi ;
      if(prank == 0) {
        hid_t datatype = H5T_NATIVE_INT ;

        hid_t dataset = H5Dopen(group_id, "is_stat",H5P_DEFAULT) ;


        H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &is_stat) ;
        H5Dclose(dataset) ;

        dataset = H5Dopen(group_id, "vec_size",H5P_DEFAULT) ;

        H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &sz) ;
        H5Dclose(dataset) ;
        fi.is_stat = is_stat ;
        fi.size = sz ;
        if(is_stat != 0) {

          hid_t dataset = H5Dopen(group_id, "second_level",H5P_DEFAULT) ;
          hid_t dataspace = H5Dget_space(dataset) ;
          hsize_t dimension = 1 ;
          H5Sget_simple_extent_dims(dataspace, &dimension, NULL) ;
          std::vector<int> vint ;
          int tmp ;
          H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT, &tmp) ;
          vint.push_back(tmp) ;
          fi.second_level = vint ;
          H5Dclose(dataset) ;
          H5Sclose(dataspace) ;
        }
      }
      return fi ;
    }
    
   
  

    //process 0 read in is_stat and vec_size, then broad cast it, then read data in fi
    frame_info read_frame_infoS(hid_t group_id,  int dom_size, MPI_Comm comm) {
      int prank = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      int is_stat = 0 ;
      int sz = 0 ;
      // read in is_stat and vector size
      if(prank == 0) {
        hid_t datatype = H5T_NATIVE_INT ;

        hid_t dataset = H5Dopen(group_id, "is_stat",H5P_DEFAULT) ;
        H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &is_stat) ;
        H5Dclose(dataset) ;

        dataset = H5Dopen(group_id, "vec_size",H5P_DEFAULT) ;
        H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &sz) ;
        H5Dclose(dataset) ;
      }
      int dim[2] ;
      dim[0] = is_stat ;
      dim[1] = sz ;
      MPI_Bcast(&dim, 2, MPI_INT, 0, comm) ;
      frame_info fi = frame_info(dim[0], dim[1]) ;

      if(dim[1] == 0) { // level 1 framing
        read_vector_intS(group_id, "first_level", fi.first_level, dom_size,comm) ;
        int total_size = fi.first_level.size() ;
        if(dim[0] != 0) {
          int dims = 0 ;
          for(int i = 0; i < total_size; ++i)
            dims += (fi.first_level)[i] ;
          // read_multi_vector_int(group_id, "second_level", dims, fi.second_level,comm) ;
          read_vector_intS(group_id, "second_level", fi.second_level, dims, comm) ;
        }
      } else { // level 2 framing only
        if(dim[0] !=0) {
          read_vector_intS(group_id, "second_level", fi.second_level, dom_size,comm) ;
        }
      }
      return fi ;
    }

    //process 0 read in is_stat and vec_size, broadcast it, then all process read in first_level and second_level in parallel 
    frame_info read_frame_infoP(hid_t group_id,  int dom_size, MPI_Comm comm) {
#ifndef H5_HAVE_PARALLEL
      return read_frame_infoS( group_id,   dom_size, comm);
#else
      int prank = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      int is_stat = 0 ;
      int sz = 0 ;
    
      //to avoid read storm
      //read in is_stat and vec_size using MPI_IO
      if(prank == 0) {
        hid_t datatype = H5T_NATIVE_INT ;

        hid_t dataset = H5Dopen(group_id, "is_stat",H5P_DEFAULT) ;
        H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &is_stat) ;
        H5Dclose(dataset) ;

        dataset = H5Dopen(group_id, "vec_size",H5P_DEFAULT) ;
        H5Dread(dataset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT, &sz) ;
        H5Dclose(dataset) ;
      }
     
      int dim[2] ;
      dim[0] = is_stat ;
      dim[1] = sz ;
      MPI_Bcast(&dim, 2, MPI_INT, 0, comm) ;
      frame_info fi = frame_info(dim[0], dim[1]) ;

      //start parallel here
      if(dim[1] == 0) { // level 1 framing
        read_vector_intP(group_id, "first_level", fi.first_level, dom_size,comm) ;
        int total_size = fi.first_level.size() ;
        if(dim[0] != 0) {
          int dims = 0 ;
          for(int i = 0; i < total_size; ++i)
            dims += (fi.first_level)[i] ;
          //read_multi_vector_intP(group_id, "second_level", dims, fi.second_level,comm) ;
          read_vector_intP(group_id, "second_level",  fi.second_level, dims, comm) ;
        }
      } else { // level 2 framing only
        if(dim[0] !=0) {
          read_vector_intP(group_id, "second_level", fi.second_level, dom_size,comm) ;
        }
      }
      return fi ;
#endif
    }
  
    void write_parameterS(hid_t group_id, storeRepP qrep, MPI_Comm comm) {
      int prank = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      frame_info fi = qrep->get_frame_info() ;
      write_frame_info_paramS(group_id,fi,comm) ;

      if(prank==0) {
        int array_size = 0 ;
        if(fi.size)
          if(fi.is_stat)
            for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
              array_size += *vi ;
          else
            array_size = fi.size  ;
        else
          if(fi.is_stat)
            for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
              array_size += *vi ;
          else
            for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
              array_size += *fvi ;

        if(array_size == 0)
          array_size = 1 ;
        hsize_t dimension = array_size ;
        int rank = 1 ;

        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = array_size ;

        hid_t dataspace =  H5Screate_simple(rank, &dimension, NULL) ;
        DatatypeP dp = qrep->getType() ;
        hid_t datatype = dp->get_hdf5_type() ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        dimension = count ;
        start += dimension ;

        hid_t dataset = H5Dcreate(group_id, "data", datatype, dataspace,
				  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

        entitySet dom = ~EMPTY ;
        qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", dom) ;
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
      }
    }

    //even though only process 0 perform write, all process must do create* collectively
    void write_parameterP(hid_t group_id, storeRepP qrep, MPI_Comm comm) {
#ifndef H5_HAVE_PARALLEL
      write_parameterS( group_id, qrep, comm);
#else
      int prank = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      frame_info fi = qrep->get_frame_info() ;
      write_frame_info_paramP(group_id,fi,comm) ;

      int array_size = 0 ;
      if(fi.size)
        if(fi.is_stat)
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        else
          array_size = fi.size  ;
      else
        if(fi.is_stat)
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        else
          for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
            array_size += *fvi ;

      if(array_size == 0)
        array_size = 1 ;
      hsize_t dimension = array_size ;
      int rank = 1 ;
      hsize_t start = 0 ;
      hsize_t stride = 1 ;
      hsize_t count = array_size ;

      hid_t dataspace =  H5Screate_simple(rank, &dimension, NULL) ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
      dimension = count ;
      start += dimension ;

      hid_t dataset = H5Dcreate(group_id, "data", datatype, dataspace,
				H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

      entitySet dom = ~EMPTY ;
      //only process 0 do the writing, and use non-parallel
      if(prank==0)  qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", dom) ;
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
    
#endif
    }

    void write_storeS(hid_t group_id, storeRepP qrep, entitySet dom, int offset, MPI_Comm comm) {
#ifdef io_performance
      MPI_Barrier(MPI_COMM_WORLD);
      Loci::stopWatch s;
      s.start();
#endif

      int prank = 0 ;
      int np = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      MPI_Comm_size(comm,&np) ;

      // Shift domain by offset
      entitySet dom_file = dom >> offset ;

      // Compute overall domain across processors
      std::vector<entitySet> dom_vector = all_collect_vectors(dom_file,comm);
      entitySet q_dom;
      for(int i = 0; i < np; i++)
        q_dom += dom_vector[i];

      if(prank == 0)
        HDF5_WriteDomain(group_id, q_dom, comm);

      // If nothing to write, don't proceed
      if(q_dom == EMPTY)
        return ;

      // Allocate buffer for largest processor buffer size
      std::vector<int> sort_max ;
      int local_size = 1 ;
      if(prank > 0)
        local_size = qrep->pack_size(dom) ;
      sort_max = all_collect_sizes(local_size,comm) ;
      int total_size = *std::max_element(sort_max.begin(), sort_max.end() );
      vector<unsigned char> tmp_send_buf(total_size) ;


      frame_info fi = qrep->get_frame_info() ;
      write_frame_infoS(group_id,fi,comm) ;


      int array_size = 0 ;
      if(fi.size)
        if(fi.is_stat)
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        else
          array_size = fi.size * dom.size() ;
      else
        if(fi.is_stat)
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        else
          for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
            array_size += *fvi ;
      std::vector<int> arr_sizes = all_collect_sizes(array_size,comm) ;
      size_t tot_arr_size = 0 ;
      for(int i = 0; i < np; ++i)
        tot_arr_size += size_t(max(0,arr_sizes[i])) ;


      if(prank != 0) {
        // Client processor code, pack data into send buffer
        MPI_Status status ;
        int send_size_buf ;
        send_size_buf = qrep->pack_size(dom) ;
        int tot_size = send_size_buf ;
        int loc_pack = 0 ;
        qrep->pack(&tmp_send_buf[0], loc_pack, total_size, dom) ;
        // Wait for signal to send message to root processor
        int flag = 0 ;
        MPI_Recv(&flag,1, MPI_INT, 0, 10, comm, &status) ;
        // received token to send, so send message
        if(flag) {
          MPI_Send(&tot_size, 1, MPI_INT, 0, 11, comm) ;
          MPI_Send(&tmp_send_buf[0], tot_size, MPI_PACKED, 0, 12, comm) ;
        }
      } else {
        // Begin writing array
        int rank = 1 ;
        hsize_t dimension = 1 ;
        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = arr_sizes[0] ;
        dimension =  tot_arr_size ;
        if(dimension != 0) {
          // First write local data
          hid_t dataspace =  H5Screate_simple(rank, &dimension, NULL) ;
          DatatypeP dp = qrep->getType() ;
          hid_t datatype = dp->get_hdf5_type() ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          dimension = count ;
          start += dimension ;

          hid_t dataset = H5Dcreate(group_id, "data", datatype, dataspace,
				    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

          qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", dom) ;
          H5Dclose(dataset) ;

          // Now write remaining vectors

          for(int i = 1; i < np; ++i) {
            MPI_Status status ;
            int recv_total_size ;
            // Allocate over 0-size-1, this allows for greater scalability when
            // sets data exceeds 2gig
            entitySet tmpset = interval(0,dom_vector[i].size()-1);

            storeRepP t_qrep = qrep->new_store(tmpset) ;

            int loc_unpack = 0 ;
            int flag = 1 ;
            MPI_Send(&flag, 1, MPI_INT, i, 10, comm) ;
            MPI_Recv(&recv_total_size, 1, MPI_INT, i, 11, comm, &status) ;
            MPI_Recv(&tmp_send_buf[0], recv_total_size, MPI_PACKED, i, 12, comm, &status) ;

            sequence tmp_seq = sequence(tmpset) ;
            t_qrep->unpack(&tmp_send_buf[0], loc_unpack, total_size, tmp_seq) ;
            dimension = arr_sizes[i] ;
            count = dimension ;

            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
            start += count ;

            dataset = H5Dopen(group_id, "data",H5P_DEFAULT) ;

            t_qrep->writehdf5(group_id, dataspace, dataset, dimension, "data", tmpset) ;
            t_qrep->allocate(EMPTY) ;

            H5Dclose(dataset) ;
          }
          H5Sclose(dataspace) ;
          H5Tclose(datatype) ;
        }
        //add else part by Qiuhan to avoid MPI communication get stuck
        else{
          for(int i = 1; i < np; ++i) {
            MPI_Status status ;
            int recv_total_size ;
            int flag = 1 ;
            MPI_Send(&flag, 1, MPI_INT, i, 10, comm) ;
            MPI_Recv(&recv_total_size, 1, MPI_INT, i, 11, comm, &status) ;
            MPI_Recv(&tmp_send_buf[0], recv_total_size, MPI_PACKED, i, 12, comm, &status) ;
          }
        }
      }
#ifdef io_performance
      MPI_Barrier(MPI_COMM_WORLD);
      double wall_time = s.stop();
      if(prank == 0) std::cerr << "                                                    serial time to write_store: "  << "  " << wall_time << endl; 
#endif
    }
  
    void write_storeP(hid_t group_id, storeRepP qrep, entitySet dom, int offset, MPI_Comm comm) {
    
#ifndef H5_HAVE_PARALLEL
      write_storeS(group_id, qrep, dom, offset, comm) ;
#else
#ifdef io_performance
      MPI_Barrier(MPI_COMM_WORLD);
      Loci::stopWatch s;
      s.start();
#endif
      int prank = 0 ;
      int np = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      MPI_Comm_size(comm,&np) ;
      if(np==1){
        write_storeS(group_id, qrep, dom, offset, comm) ;
        return;
      }
      // Shift domain by offset
      entitySet dom_file = dom >> offset ;

      // Compute overall domain across processors
      std::vector<entitySet> dom_vector = all_collect_vectors(dom_file,comm);
      entitySet q_dom;
      for(int i = 0; i < np; i++)
        q_dom += dom_vector[i];

    
    
      HDF5_WriteDomain(group_id, q_dom, comm); //no comm info here ?
      
      // If nothing to write, don't proceed
      if(q_dom == EMPTY)
        return ;

   

      //write frame info
      frame_info fi = qrep->get_frame_info() ;
      write_frame_infoP(group_id,fi,comm) ;

      //get arr_sizes and total_arr_size
      int array_size = 0 ;
      if(fi.size)
        if(fi.is_stat)
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        else
          array_size = fi.size * dom.size() ;
      else
        if(fi.is_stat)
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        else
          for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
            array_size += *fvi ;
      std::vector<int> arr_sizes = all_collect_sizes(array_size,comm) ;
      size_t tot_arr_size = 0 ;
      for(int i = 0; i < np; ++i)
        tot_arr_size += size_t(max(0,arr_sizes[i])) ;


   
      int rank = 1 ;
      hsize_t dimension = 1 ;
      hsize_t start = 0 ;
      hsize_t stride = 1 ;
      hsize_t count = arr_sizes[prank] ;

      for(int i = 0; i < prank; i++)start += arr_sizes[i];
      
      dimension =  tot_arr_size ;
      if(dimension != 0) {
        // First write local data
        hid_t dataspace =  H5Screate_simple(rank, &dimension, NULL) ;
        DatatypeP dp = qrep->getType() ;
        hid_t datatype = dp->get_hdf5_type() ;

        hid_t dataset = H5Dcreate(group_id, "data", datatype, dataspace,
				  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
        
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        dimension = count ;
        hid_t xfer_plist = create_xfer_plist(Loci::hdf5_const::dxfer_coll_type);
        qrep->writehdf5P(group_id, dataspace, dataset, dimension, "data", dom, xfer_plist) ;
        H5Pclose(xfer_plist);
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
        H5Tclose(datatype) ;
      }
#ifdef io_performance
      MPI_Barrier(MPI_COMM_WORLD);
      double wall_time = s.stop();
      if(prank == 0) std::cerr << "                                                    parallel time to write_store "  << "  " << wall_time << endl;
#endif
#endif
    }

    void write_containerS(hid_t group_id, storeRepP qrep) {
      int offset = 0 ;
      if(qrep->RepType() == PARAMETER)
        write_parameterS(group_id,qrep,MPI_COMM_WORLD) ;
      else
        write_storeS(group_id,qrep,qrep->domain(),offset,MPI_COMM_WORLD) ;
    }
  
    void write_containerP(hid_t group_id, storeRepP qrep) {
#ifndef H5_HAVE_PARALLEL
      write_containerS(group_id, qrep);
#else
      int offset = 0 ;
      if(qrep->RepType() == PARAMETER)
        write_parameterP(group_id,qrep,MPI_COMM_WORLD) ;
      else
        write_storeP(group_id,qrep,qrep->domain(),offset,MPI_COMM_WORLD) ;
#endif
    }

  
     

    void read_storeS(hid_t group_id, storeRepP qrep, int &offset, MPI_Comm comm) {
#ifdef io_performance   
      MPI_Barrier(MPI_COMM_WORLD);    
      Loci::stopWatch s;
      s.start();
#endif
      offset = 0 ;
      int prank = 0 ;
      int np = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      MPI_Comm_size(comm,&np) ;

      // Here we read in a store container.  First lets read in the domain
      entitySet q_dom ;
      if(prank == 0)
        HDF5_ReadDomain(group_id, q_dom) ;

      // Now lets share the domain with all other processors ;
      q_dom = BcastEntitySet(q_dom,0,comm) ;

      unsigned char* tmp_buf = 0;
      std::vector<int> interval_sizes ;
      entitySet dom ;
      if(q_dom != EMPTY) {
        vector<entitySet> ptn = simplePartition(q_dom.Min(),q_dom.Max(),comm) ;
        for(int i=0;i<np;++i) {
          entitySet qset = ptn[i] &q_dom ;
          interval_sizes.push_back(qset.size()) ;
        }
        dom = ptn[prank] &q_dom ;
      } else
        for(int i=0;i<np;++i)
          interval_sizes.push_back(0) ;

      if(q_dom==EMPTY) {
        qrep->allocate(q_dom) ;
        return ;
      }
      offset = dom.Min() ;
      dom <<= offset ;
      qrep->allocate(dom) ;

      frame_info fi = read_frame_infoS(group_id,dom.size(),comm) ;
      int array_size = 0 ;
      int vec_size = 0 ;

      if(fi.size) {
        if(fi.is_stat) {
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
          vec_size = fi.second_level.size() ;
        } else {
          if(fi.size > 1)
            qrep->set_elem_size(fi.size) ;
          array_size = fi.size * dom.size() ;
        }
      } else {
        if(fi.is_stat) {
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
          vec_size = fi.second_level.size() + dom.size() ;
        } else {
          for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
            array_size += *fvi ;
          vec_size = dom.size() ;
        }
      }

      std::vector<int> tmp_sizes = all_collect_sizes(vec_size,comm) ;
      int max_tmp_size = *std::max_element(tmp_sizes.begin(), tmp_sizes.end()) ;
      int max_eset_size = *std::max_element(interval_sizes.begin(), interval_sizes.end()) ;
      int* tmp_int  ;
      tmp_int = new int[max_tmp_size] ;
      std::vector<int> arr_sizes = all_collect_sizes(array_size,comm) ;
      //      size_t tot_arr_size = 0 ;
      //      for(int i = 0; i < np; ++i)
      //        tot_arr_size += size_t(max(0,arr_sizes[i])) ;
      MPI_Status status ;
      if(prank != 0) {
        int t = 0 ;
        if(fi.size) {
          if(fi.size > 1)
            qrep->set_elem_size(fi.size) ;
          if(fi.is_stat)
            for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
              tmp_int[t++] = *vi ;
        } else {
          if(fi.is_stat) {
            for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
              tmp_int[t++] = *fvi ;

            for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
              tmp_int[t++] = *vi ;
          } else
            for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
              tmp_int[t++] = *fvi ;
        }
        if(tmp_sizes[prank])
          MPI_Send(tmp_int, tmp_sizes[prank], MPI_INT, 0, 10, comm) ;
        int total_size = 0 ;
        MPI_Recv(&total_size, 1, MPI_INT, 0, 11,comm, &status) ;
        tmp_buf = new unsigned char[total_size] ;
        MPI_Recv(tmp_buf, total_size, MPI_PACKED, 0, 12, comm, &status) ;
        sequence tmp_seq = sequence(dom) ;
        int loc_unpack = 0 ;
        qrep->unpack(tmp_buf, loc_unpack, total_size, tmp_seq) ;
      } else {
        // processor zero

        hid_t dataset =  H5Dopen(group_id, "data",H5P_DEFAULT) ;
        hid_t dataspace = H5Dget_space(dataset) ;

        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = 0 ;
        int curr_indx = 0 ;
        int total_size = 0 ;
        int tmp_total_size = 0 ;
        entitySet max_set;
        if(max_eset_size > 0)
          max_set = interval(0, max_eset_size-1) ;
        if(np == 1) {
          qrep->allocate(dom) ;
          if(fi.size)
            if(fi.size > 1)
              qrep->set_elem_size(fi.size) ;
          hsize_t dimension = arr_sizes[0] ;
          count = dimension ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
          qrep->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, dom) ;
        } else {
          storeRepP tmp_sp ;
          if(fi.size)
            tmp_sp = qrep->new_store(max_set) ;
          for(int p = 0; p < np; ++p) {
            entitySet local_set;
            if(interval_sizes[p] > 0)
              local_set = entitySet(interval(curr_indx, interval_sizes[p]+curr_indx-1)) ;
            curr_indx += interval_sizes[p] ;
            hsize_t dimension = arr_sizes[p] ;
            count = dimension ;
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
            entitySet tmp_set;
            if(local_set.size())
              tmp_set = interval(0, local_set.size()-1) ;
            if(p && tmp_sizes[p]) {
              MPI_Recv(tmp_int, tmp_sizes[p], MPI_INT, p, 10, comm, &status) ;
              std::vector<int> vint, fvint ;
              int t = 0 ;
              if(fi.size) {
                if(fi.is_stat) {
                  for(int i = 0; i < tmp_sizes[p]; ++i)
                    vint.push_back(tmp_int[t++]) ;
                  fi.second_level = vint ;
                }
              } else {
                for(size_t i = 0; i < local_set.size(); ++i)
                  fvint.push_back(tmp_int[t++]) ;
                for(size_t i = 0; i < tmp_sizes[p]-local_set.size(); ++i)
                  vint.push_back(tmp_int[t++]) ;
                fi.first_level = fvint ;
                fi.second_level = vint ;
              }
            }
            storeRepP t_sp ;
            int t = 0 ;
            if(p == 0)
              if(!fi.size)
                for(std::vector<int>::const_iterator vi = fi.first_level.begin(); vi != fi.first_level.end(); ++vi)
                  tmp_int[t++] = *vi ;
            if(fi.size) {
              tmp_sp->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, tmp_set) ;
              tmp_total_size = tmp_sp->pack_size(tmp_set) ;
            } else {
              t_sp = qrep->new_store(tmp_set, tmp_int) ;
              t_sp->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, tmp_set) ;
              tmp_total_size = t_sp->pack_size(tmp_set) ;
            }

            if(tmp_total_size > total_size) {
              total_size = tmp_total_size ;
              if(p)
                delete [] tmp_buf ;
              tmp_buf = new unsigned char[total_size] ;
            }
            start += count ;
            int loc = 0 ;
            if(fi.size)
              tmp_sp->pack(tmp_buf, loc, total_size, tmp_set) ;
            else
              t_sp->pack(tmp_buf, loc, total_size, tmp_set) ;
            if(p == 0) {
              int loc_unpack = 0 ;

              if(fi.size)
                if(fi.size > 1)
                  qrep->set_elem_size(fi.size) ;

              sequence tmp_seq = sequence(dom) ;
              qrep->allocate(dom) ;
              qrep->unpack(tmp_buf, loc_unpack, total_size, tmp_seq) ;
            } else {
              MPI_Send(&total_size, 1, MPI_INT, p, 11, comm) ;
              MPI_Send(tmp_buf, total_size, MPI_PACKED, p, 12, comm) ;
            }
          }
        }
        H5Dclose(dataset) ;
        H5Sclose(dataspace) ;
      }
      delete [] tmp_buf ;
      delete [] tmp_int ;
#ifdef io_performance
      MPI_Barrier(MPI_COMM_WORLD);
      double wall_time = s.stop();
      if(prank == 0) std::cerr << "                                                     serial time to read_store "  << "  " << wall_time << endl; 
#endif
    }

    //process 0 read in domain and broadcast it, data is read in in parallel
    void read_storeP(hid_t group_id, storeRepP qrep, int &offset, MPI_Comm comm) {
#ifndef H5_HAVE_PARALLEL
      read_storeS(group_id, qrep, offset, comm);
#else
#ifdef io_performance
      MPI_Barrier(MPI_COMM_WORLD);
      Loci::stopWatch s;
      s.start();
#endif
      offset = 0 ;
      int prank = 0 ;
      int np = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      MPI_Comm_size(comm,&np) ;

      

      //to avoid read storm, process 0 read in domain and broadcast it
      entitySet q_dom ;
      if(prank == 0) HDF5_ReadDomain(group_id, q_dom) ;
      // Now lets share the domain with all other processors ;
      q_dom = BcastEntitySet(q_dom,0,comm) ;

      
      //each process do partition and compute dom
      std::vector<int> interval_sizes ;
      entitySet dom ;
      if(q_dom != EMPTY) {
        vector<entitySet> ptn = simplePartition(q_dom.Min(),q_dom.Max(),comm) ;
        for(int i=0;i<np;++i) {
          entitySet qset = ptn[i] &q_dom ;
          interval_sizes.push_back(qset.size()) ;
        }
        dom = ptn[prank] &q_dom ;
      } else
        for(int i=0;i<np;++i)
          interval_sizes.push_back(0) ;

      if(q_dom==EMPTY) {
        qrep->allocate(q_dom) ;
        return ;
      }
      offset = dom.Min() ;
      dom <<= offset ;
      qrep->allocate(dom) ;
    

      //each process read in frame_info and compute vez_size and array_size
      frame_info fi = read_frame_infoP(group_id,dom.size(),comm) ;
      int array_size = 0 ;
    
      if(fi.size) {
        if(fi.is_stat) {
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        } else {
          if(fi.size > 1)
            qrep->set_elem_size(fi.size) ;
          array_size = fi.size * dom.size() ;
        }
      } else {
        if(fi.is_stat) {
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        } else {
          for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
            array_size += *fvi ;
        }
      }
   
      //all_collect, now each process also know array_size and vec_size of others
    
      std::vector<int> arr_sizes = all_collect_sizes(array_size,comm) ;
         
      //  each process open dataset and read in data
      hid_t dataset =  H5Dopen(group_id, "data",H5P_DEFAULT) ;

      hid_t dataspace = H5Dget_space(dataset) ;

      hsize_t start = 0 ;
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      
      if(np == 1) {
        qrep->allocate(dom) ;
        if(fi.size)
          if(fi.size > 1)
            qrep->set_elem_size(fi.size) ;
        hsize_t dimension = arr_sizes[0] ;
        count = dimension ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        qrep->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, dom) ;
      } else {
        qrep->allocate(dom) ;
        if(fi.size)
          if(fi.size > 1)
            qrep->set_elem_size(fi.size) ;
        
        //first compute  start
        for(int p = 0; p < prank; ++p){
          start +=  arr_sizes[p] ;
        }
        hsize_t dimension = arr_sizes[prank] ;
        count = dimension ;
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        hid_t xfer_plist = create_xfer_plist(Loci::hdf5_const::dxfer_coll_type);
        qrep->readhdf5P(group_id, dataspace, dataset, dimension, "data", fi, dom, xfer_plist) ;
        H5Pclose(xfer_plist);
      }          
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;
#ifdef io_performance
      MPI_Barrier(MPI_COMM_WORLD);
      double wall_time = s.stop();
      if(prank == 0) std::cerr << "                                                    parallel time to read_store " << "  " << wall_time << endl; 
#endif
#endif
    }
  }

  //to avoid read storm, no read_parameterP
  void read_parameter(hid_t group_id, storeRepP qrep, MPI_Comm comm) {
    int pack_size = 0 ;
    int prank = 0 ;
    int np = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&np) ;

    frame_info fi = pio::read_frame_info_paramS(group_id,1,comm) ;
    if(prank==0) {
      int array_size = 0 ;
      if(fi.size)
        if(fi.is_stat) {
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        } else {
          if(fi.size > 1)
            qrep->set_elem_size(fi.size) ;
          array_size = fi.size ;
        } else {
        if(fi.is_stat) {
          for(std::vector<int>::const_iterator vi = fi.second_level.begin(); vi != fi.second_level.end(); ++vi)
            array_size += *vi ;
        } else {
          for(std::vector<int>::const_iterator fvi = fi.first_level.begin(); fvi != fi.first_level.end(); ++fvi)
            array_size += *fvi ;
        }
      }

      hid_t dimension = array_size ;
      hid_t dataset =  H5Dopen(group_id, "data",H5P_DEFAULT) ;
      hid_t dataspace = H5Dget_space(dataset) ;
      entitySet dom = ~EMPTY ;
      qrep->readhdf5(group_id, dataspace, dataset, dimension, "data", fi, dom) ;
      H5Dclose(dataset) ;
      H5Sclose(dataspace) ;

      pack_size = qrep->pack_size(dom) ;

    }

    // Now broadcast the result to other processors
    if(np > 1) {
      MPI_Bcast(&pack_size,1,MPI_INT,0,comm) ;
      unsigned char *pack_buf = new unsigned char[pack_size] ;
      if(prank == 0) {
        int loc_pack = 0 ;
        int sz = pack_size ;
        entitySet dom = ~EMPTY ;
        qrep->pack(pack_buf,loc_pack,sz,dom) ;
      }
      MPI_Bcast(pack_buf,pack_size,MPI_PACKED,0,comm) ;
      if(prank != 0) {
        int loc_pack = 0 ;
        int sz = pack_size ;
        entitySet dom = ~EMPTY ;
        qrep->unpack(pack_buf,loc_pack,sz,dom) ;
      }
      delete[] pack_buf ;
    }
  }


 
  entitySet findBoundingSet(entitySet in) {
    Entity max_val = in.Max() ;
    Entity min_val = in.Min() ;
    Entity gmin_val = min_val;
    Entity gmax_val = max_val ;
    MPI_Allreduce(&min_val,&gmin_val,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
    MPI_Allreduce(&max_val,&gmax_val,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;

    return entitySet(interval(gmin_val,gmax_val)) ;
  }

  vector<sequence> transposeSeq(const vector<sequence> sv) {
    vector<int> send_sz(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      send_sz[i] = sv[i].num_intervals()*2 ;
    vector<int> recv_sz(MPI_processes) ;
    MPI_Alltoall(&send_sz[0],1,MPI_INT,
                 &recv_sz[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;
    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      size_send += send_sz[i] ;
      size_recv += recv_sz[i] ;
    }
    //    outRep->allocate(new_alloc) ;
    int *send_store = new int[size_send] ;
    int *recv_store = new int[size_recv] ;
    int *send_displacement = new int[MPI_processes] ;
    int *recv_displacement = new int[MPI_processes] ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sz[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sz[i-1] ;
    }
    for(int i = 0; i <  MPI_processes; ++i)
      for(size_t j=0;j<sv[i].num_intervals();++j) {
        send_store[send_displacement[i]+j*2] = sv[i][j].first ;
        send_store[send_displacement[i]+j*2+1] = sv[i][j].second ;
      }


    MPI_Alltoallv(send_store,&send_sz[0], send_displacement , MPI_INT,
                  recv_store, &recv_sz[0], recv_displacement, MPI_INT,
                  MPI_COMM_WORLD) ;

    vector<sequence> sv_t(MPI_processes) ;
    for(int i = 0; i <  MPI_processes; ++i)
      for(int j=0;j<recv_sz[i]/2;++j) {
        int i1 = recv_store[recv_displacement[i]+j*2]  ;
        int i2 = recv_store[recv_displacement[i]+j*2+1] ;
        sv_t[i] += interval(i1,i2) ;
      }
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;

    return sv_t ;
  }

  // convert domain in local numbering into key space
  int getKeyDomain(entitySet dom, fact_db::distribute_infoP dist, MPI_Comm comm) {
    int kdl = 0 ;
    entitySet searchdom = dom & dist->key_domain.domain() ;
    FORALL(searchdom,i) {
      int key = dist->key_domain[i] ;
      kdl = std::max<int>(key,kdl) ;
    } ENDFORALL ;
    int kd=-1 ;
    MPI_Allreduce(&kdl,&kd,1,MPI_INT,MPI_MAX,comm) ;

    bool failure = false ;
    FORALL(searchdom,i) {
      int key = dist->key_domain[i] ;
      if(kd != key){
        failure = true ;
      }
    } ENDFORALL ;


    kdl = failure?-1:kd ;
    MPI_Allreduce(&kdl,&kd,1,MPI_INT,MPI_MIN,comm) ;
    return kd ;
  }


  // Convert container from local numbering to file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // return offset in file numbering (each processor will allocate from zero,
  // add offset to domain to get actual file numbering)
  // distribution info pointer (dist)
  // MPI Communicator
  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm) {


    // Get local numbering of entities owned by this processor, only write
    // out these entities.
    constraint my_entities ;
    my_entities = dist->my_entities ;
    dom = *my_entities & dom ;

    // Get mapping from local to global numbering
    Map l2g ;
    l2g = dist->l2g.Rep() ;
    // Compute domain in global numbering
    entitySet dom_global = l2g.image(dom) ;

    // This shouldn't happen
    FATAL(dom.size() != dom_global.size()) ;

    int kd =  getKeyDomain(dom, dist, comm) ;
    if(kd< 0) {
      cerr << "Local2FileOrder not in single keyspace!" << endl ;
      kd = 0 ;
    }
    // Now get global to file numbering
    Map l2f ;
    l2f = dist->l2f.Rep() ;

    // Compute map from local numbering to file numbering
    Map newnum ;
    newnum.allocate(dom) ;
    FORALL(dom,i) {
      newnum[i] = l2f[i] ; 
    } ENDFORALL ;

    int imx = std::numeric_limits<int>::min() ;
    int imn = std::numeric_limits<int>::max() ;

    // Find bounds in file numbering from this processor
    FORALL(dom,i) {
      imx = max(newnum[i],imx) ;
      imn = min(newnum[i],imn) ;
    } ENDFORALL ;

    // Find overall bounds
    imx = GLOBAL_MAX(imx) ;
    imn = GLOBAL_MIN(imn) ;

    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    // Get partitioning of file numbers across processors
    vector<entitySet> out_ptn = simplePartition(imn,imx,comm) ;

    // Now compute where to send data to put in file ordering
    vector<entitySet> send_sets(p) ;
    vector<sequence> send_seqs(p) ;

    // Loop over processors and compute sets of entities to send
    // To efficiently compute this mapping, first sort the transpose
    // of the newnum map to quickly find the set of entities to send
    // without searching entire newnum map for each processor
    vector<pair<int,int> > file2num(dom.size()) ;
    size_t cnt = 0 ;
    FORALL(dom,ii) {
      file2num[cnt].first = newnum[ii] ;
      file2num[cnt].second = ii ;
      cnt++ ;
    } ENDFORALL ;
    sort(file2num.begin(),file2num.end()) ;

    // Check each processor, find out which sets to send
    cnt = 0 ;
    for(int i=0;i<p;++i) {
      int mxi = out_ptn[i].Max() ;
      while(cnt < file2num.size() && file2num[cnt].first <= mxi) {
        send_sets[i] += file2num[cnt].second ;
        cnt++ ;
      }
      sequence s ;
      FORALL(send_sets[i],j) {
        s+= newnum[j] ;
      } ENDFORALL ;
      send_seqs[i] = s ;
    }

    //Get the sequences of where we place the data when we receive it
    vector<sequence> recv_seqs = transposeSeq(send_seqs) ;


    // shift by the offset
    offset = out_ptn[prank].Min() ;
    for(int i=0;i<p;++i)
      recv_seqs[i] <<= offset ;

    // Compute allocation domain
    entitySet file_dom ;
    for(int i=0;i<p;++i)
      file_dom += entitySet(recv_seqs[i]) ;

    // allocate store over shifted domain
    storeRepP qcol_rep ;
    qcol_rep = sp->new_store(file_dom) ;

    // Now communicate the container
    vector<int> send_sizes(p),recv_sizes(p) ;

    for(int i=0;i<p;++i)
      send_sizes[i] = sp->pack_size(send_sets[i]) ;

    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 comm) ;

    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz) ;
    vector<unsigned char> recv_store(recv_sz) ;


    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
               send_sets[i]) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
                  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
                  comm) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      qcol_rep->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                       recv_seqs[i]) ;
    }
    return qcol_rep ;
  }

  // Convert container from local numbering to output file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // fact_db pointer  (facts)
  // MPI Communicator
  storeRepP Local2FileOrder_output(storeRepP sp, entitySet dom,
                                   fact_db& facts, MPI_Comm comm) {

    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;

    if(p==1) return sp;


    fact_db::distribute_infoP dist = facts.get_distribute_info() ;
    int kd =  getKeyDomain(dom, dist, comm) ;
    if(kd < 0) {
      cerr << "unable to finde key domain in File2LocalOrderOutput"
           << endl ;
      kd = 0 ;
    }
    vector<entitySet> out_ptn = facts.get_init_ptn(kd) ;
    // Get mapping from local to global numbering
    Map l2g ;
    l2g = dist->l2g.Rep() ;
    // Compute domain in global numbering
    entitySet dom_global = l2g.image(dom) ;
    // This shouldn't happen
    FATAL(dom.size() != dom_global.size()) ;

    // Now compute where to send data to put in output ordering
    vector<entitySet> send_sets(p) ;//local numbering
    vector<sequence> send_seqs(p) ;//global numbering

    // Loop over processors and compute sets of entities to send
    // To efficiently compute this mapping, first sort the transpose
    // of the newnum map to quickly find the set of entities to send
    // without searching entire newnum map for each processor
    vector<pair<int,int> > file2num(dom.size()) ;//global2local
    size_t cnt = 0 ;
    FORALL(dom,ii) {
      file2num[cnt].first = l2g[ii] ;
      file2num[cnt].second = ii ;
      cnt++ ;
    } ENDFORALL ;

    //sort according to global numbering
    sort(file2num.begin(),file2num.end()) ;

    // Check each processor, find out which sets to send
    cnt = 0 ;
    for(int i=0;i<p;++i) {
      int mxi = out_ptn[i].Max() ;
      while(cnt < file2num.size() && file2num[cnt].first <= mxi) {
        send_sets[i] += file2num[cnt].second ;
        cnt++ ;
      }
      sequence s ;
      FORALL(send_sets[i],j) {
        s+= l2g[j] ;
      } ENDFORALL ;
      send_seqs[i] = s ;
    }

    //Get the sequences of where we place the data when we receive it
    vector<sequence> recv_seqs = transposeSeq(send_seqs) ;

    // don't need shift by the offset because global number is unique
    int offset = out_ptn[prank].Min() ;
    for(int i=0;i<p;++i){
      recv_seqs[i] <<= offset ;
    }


    // Compute allocation domain
    entitySet file_dom =EMPTY;
    for(int i=0;i<p;++i){
      file_dom += entitySet(recv_seqs[i]) ;
    }


    // Now communicate the container
    vector<int> send_sizes(p),recv_sizes(p) ;

    for(int i=0;i<p;++i){
      if(send_sets[i] != EMPTY) send_sizes[i] = sp->pack_size(send_sets[i]) ;
      else send_sizes[i] = 0;
    }
    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 comm) ;


    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz) ;
    vector<unsigned char> recv_store(recv_sz) ;


    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
               send_sets[i]) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
                  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
                  comm) ;


    if(file_dom == EMPTY) return NULL;

    // allocate store
    storeRepP result_rep;
    result_rep = sp->new_store(file_dom) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      result_rep->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                         recv_seqs[i]) ;
    }

    return result_rep ;
  }

  void File2LocalOrder(storeRepP &result, entitySet resultSet,
                       storeRepP input, int offset,
                       fact_db::distribute_infoP dist,
                       MPI_Comm comm) {
    Map newnum ;
    newnum.allocate(resultSet) ;

    if(dist !=0 ) {
      int kd =  getKeyDomain(resultSet, dist, comm) ;

      if(kd < 0) {
        cerr << "File2LocalOrder not in single keyspace!" << endl ;
        kd = 0 ;
      }

      Map l2f ;
      l2f = dist->l2f.Rep() ;
      FORALL(resultSet,i) {
        newnum[i] = l2f[i] ; 
      } ENDFORALL ;
    } else {
      result->copy(input,resultSet) ;
      return ;
    }

    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int mn = input->domain().Min() ;
    int mx = input->domain().Max() ;
    if(input->domain() != EMPTY) {
      mn += offset ;
      mx += offset ;
    }
    vector<int> allmx(p) ;
    vector<int> allmn(p) ;
    MPI_Allgather(&mx,1,MPI_INT,&allmx[0],1,MPI_INT,comm) ;
    MPI_Allgather(&mn,1,MPI_INT,&allmn[0],1,MPI_INT,comm) ;

    vector<pair<int,int> > file_requests ;
    FORALL(resultSet,i) {
      file_requests.push_back(pair<int,int>(newnum[i],i)) ;
    } ENDFORALL ;
    sort(file_requests.begin(),file_requests.end()) ;
    // Get distribution plan
    vector<vector<pair<int,int> > > dist_plan(p) ;

    int proc = 0 ;
    for(size_t i=0;i<file_requests.size();++i) {
      int fn = file_requests[i].first ;
      while(proc < p && (fn < allmn[proc] || fn > allmx[proc]))
        proc++ ;
      if(fn < allmn[proc] || fn > allmx[proc]) {
        cerr << "Unable to find processor that contains entity!" << endl ;
        Abort() ;
      }
      dist_plan[proc].push_back(pair<int,int>(fn,file_requests[i].second)) ;
    }

    // Compute recv requests from distribution plan
    vector<sequence> recv_seq(p),send_req(p) ;
    for(int i=0;i<p;++i) {
      sort(dist_plan[i].begin(),dist_plan[i].end()) ;
      sequence s1,s2 ;
      int psz = dist_plan[i].size() ;
      for(int j=0;j<psz;++j) {
        s1 +=dist_plan[i][j].first ;
        s2 +=dist_plan[i][j].second ;
      }
      send_req[i] = s1 ;
      recv_seq[i] = s2 ;
    }

    // Transpose the send requests to get the sending sequences
    // from this processor
    vector<sequence> send_seq = transposeSeq(send_req) ;
    vector<entitySet> send_sets(p) ;
    for(int i=0;i<p;++i) {
      send_seq[i] <<= offset ;
      send_sets[i] = entitySet(send_seq[i]) ;
    }

    vector<int> send_sizes(p), recv_sizes(p) ;


    for(int i=0;i<p;++i)
      send_sizes[i] = input->pack_size(send_sets[i]) ;

    MPI_Alltoall(&send_sizes[0],1,MPI_INT,&recv_sizes[0],1,MPI_INT, comm) ;

    vector<int> send_dspl(p), recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz), recv_store(recv_sz) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      input->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
                  send_sets[i]) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
                  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
                  comm) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      result->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                     recv_seq[i]) ;
    }
  }


  void writeContainerRAW(hid_t file_id, std::string vname,
                         storeRepP var, MPI_Comm comm) {
    hid_t group_id = 0 ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    if(use_parallel_io || prank == 0)
      group_id = H5Gcreate(file_id, vname.c_str(),
			   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;

    if(var->RepType() != PARAMETER) {
      int offset = 0 ;
      if(use_parallel_io)
        pio::write_storeP(group_id,var,var->domain(),offset,comm) ;
      else
        pio::write_storeS(group_id,var,var->domain(),offset,comm) ;
    } else {
      if(use_parallel_io)
        pio::write_parameterP(group_id, var,comm) ;
      else 
        pio::write_parameterS(group_id, var,comm) ;
    }
    if(use_parallel_io || prank == 0)
      H5Gclose(group_id) ;
  }
  
 
  void readContainerRAW(hid_t file_id, std::string vname,
                        storeRepP var,
                        MPI_Comm comm ) {
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    hid_t group_id = 0;
    if(prank == 0||use_parallel_io)
      group_id = H5Gopen(file_id, vname.c_str(),H5P_DEFAULT) ;

    if(var->RepType() == PARAMETER) {
      read_parameter(group_id, var, comm) ;
      if(prank == 0 || use_parallel_io)
        H5Gclose(group_id) ;
      return ;
    }

    // Read in store in file numbering
    int offset = 0 ;
    if(use_parallel_io) pio::read_storeP( group_id, var,offset,comm) ;
    else pio::read_storeS( group_id, var,offset,comm) ;
    var->shift(offset) ;

    if(prank == 0 || use_parallel_io)
      H5Gclose(group_id) ;
  }
 
 
  void redistribute_write_container(hid_t file_id, std::string vname,
                                    storeRepP var, fact_db &facts) {
    fact_db::distribute_infoP dist = facts.get_distribute_info() ;
    if(dist == 0) {
      writeContainerRAW(file_id,vname,var,MPI_COMM_WORLD) ;
      return ;
    }

    hid_t group_id = 0 ;

    if(MPI_rank == 0||use_parallel_io)
      group_id = H5Gcreate(file_id, vname.c_str(),
			   H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT) ;

    // Redistribute container to map from local to global numbering
    if(var->RepType() != PARAMETER && MPI_processes != 1) {
      // parallel store write.. reorder to file numbering then write out
      // reorder from local to file numbering
      int offset = 0 ;
      entitySet dom = var->domain() ;
      // Create container vardist that is ordered across processors in the
      // file numbering, the domain of this container shifted by offset
      // is the actual file numbering.
      storeRepP vardist = Local2FileOrder(var,dom,offset,dist,MPI_COMM_WORLD) ;
      // Write out container that has been distributed in the file numbering
      if(use_parallel_io)
        pio::write_storeP(group_id,vardist,vardist->domain(),offset,MPI_COMM_WORLD) ;
      else
        pio::write_storeS(group_id,vardist,vardist->domain(),offset,MPI_COMM_WORLD) ;
    } else {
      // No need to reorder container if parameter or only one
      // processor, so just write out the container.
      if(use_parallel_io)
        pio::write_containerP(group_id, var) ;
      else
        pio::write_containerS(group_id, var) ;
    }

    if(MPI_rank == 0||use_parallel_io)
      H5Gclose(group_id) ;
  }


  int getMinFileNumberFromLocal(entitySet read_set,
                                fact_db::distribute_infoP dist ) {

    int minIDfl = std::numeric_limits<int>::max() ;
    int kd =  getKeyDomain(read_set, dist, MPI_COMM_WORLD) ;
    if(kd< 0) {
      cerr << "read_set not in single keyspace!" << endl ;
      kd = 0 ;
    }
    // local 2 file
    Map l2f ;
    l2f = dist->l2f.Rep() ;

    // Compute map from local numbering to file numbering
    FORALL(read_set,ii) {
      minIDfl = min(minIDfl,l2f[ii]) ;
    } ENDFORALL ;
    int minIDf = minIDfl ;
    MPI_Allreduce(&minIDfl,&minIDf,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
    return minIDf ;
  }

  //serial/parallel io
  void read_container_redistribute(hid_t file_id, std::string vname,
                                   storeRepP var, entitySet read_set,
                                   fact_db &facts) {
    hid_t group_id = 0;
    if(use_parallel_io || MPI_rank == 0)
      group_id = H5Gopen(file_id, vname.c_str(),H5P_DEFAULT) ;

    if(var->RepType() == PARAMETER) {
      read_parameter(group_id, var, MPI_COMM_WORLD) ;
      if(use_parallel_io || MPI_rank == 0)
        H5Gclose(group_id) ;
      return ;
    }

    // Read in store in file numbering
    int offset = 0 ;
    storeRepP new_store = var->new_store(EMPTY) ;
    if(use_parallel_io )
      pio::read_storeP( group_id, new_store,offset,MPI_COMM_WORLD) ;
    else pio::read_storeS( group_id, new_store,offset,MPI_COMM_WORLD) ;
    
    if(use_parallel_io || MPI_rank == 0)
      H5Gclose(group_id) ;

    // map from file number to local numbering
    fact_db::distribute_infoP dist = facts.get_distribute_info() ;
    if(dist != 0) {
      // Correct offset if file numbering changes.  Assume read_set is being
      // read in over the same set
      int minID = offset ;
      MPI_Bcast(&minID,1,MPI_INT,0,MPI_COMM_WORLD) ;
      const int minIDf = getMinFileNumberFromLocal(read_set,dist) ;
      const int correct = minIDf - minID ;
      offset += correct  ;

      // Allocate space for reordered container
      storeRepP result = var->new_store(read_set) ;
      File2LocalOrder(result,read_set,new_store,offset,dist,MPI_COMM_WORLD) ;
      // Copy results into container
      if(read_set == EMPTY) {
        read_set = result->domain() ;
        var->allocate(read_set) ;
      }
      var->copy(result,read_set) ;
    } else {
      if(read_set == EMPTY) {
        read_set = new_store->domain() ;
        var->allocate(read_set) ;
      } else {
        offset = read_set.Min() - new_store->domain().Min() ;
        if(offset != 0) {
          // shift new store by offset to correct alignment
          new_store->shift(offset) ;
        }
      }
      var->copy(new_store,read_set) ;
    }

   
  }

 
  void writeSetIds(hid_t file_id, entitySet local_set, fact_db &facts) {
    vector<int> ids(local_set.size()) ;

    int c = 0 ;
    if(MPI_processes > 1) {
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      int kd =  getKeyDomain(local_set, df, MPI_COMM_WORLD) ;
      if(kd < 0) {
        cerr << "unable to find distribute info in writeSetIds" << endl;
        kd = 0 ;
      }

      Map l2f ;
      l2f = df->l2f.Rep() ;

      FORALL(local_set,ii) {
        ids[c++] = l2f[ii] ; 
      } ENDFORALL ;
    } else {
      // Note, this means that a single processor run will not be compatible
      // with a parallel processor run
      FORALL(local_set,ii) {
        ids[c++] = ii ;
      } ENDFORALL ;
    }
        
    writeUnorderedVector(file_id,"entityIds",ids) ;

  }
  
  
  

  hid_t createUnorderedFile(const char * filename, entitySet set, fact_db &facts) {
    hid_t file_id = 0;
    hid_t group_id = 0 ;
   
    // file_id = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;//error

    file_id = writeVOGOpen(filename);
    if(MPI_rank == 0 || use_parallel_io) {
      group_id = H5Gcreate(file_id,"dataInfo",
			   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    }
    writeSetIds(group_id,set,facts) ;
    if(MPI_rank == 0 || use_parallel_io)
      H5Gclose(group_id) ;
    return file_id ;
  }

  
  void closeUnorderedFile(hid_t file_id) {
    if(MPI_rank == 0 || use_parallel_io)
      H5Fclose(file_id) ;
  }
 

  void generalMPIComm(storeRepP op,
                      storeRepP sp,
                      const vector<entitySet> &sendSets,
                      const vector<sequence> &recvSeqs,
                      MPI_Comm comm) {
    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;

    // Now get the send sizes
    vector<int> send_sizes(p,0),recv_sizes(p,0) ;

    for(int i=0;i<p;++i) {
      if(sendSets[i].num_intervals()>0)
        send_sizes[i] = sp->pack_size(sendSets[i]) ;
      else
        send_sizes[i] = 0 ;
    }

    MPI_Alltoall(&send_sizes[0],1,MPI_INT, &recv_sizes[0],1,MPI_INT, comm) ;
    int maxbuf_send = send_sizes[0] ;
    int maxbuf_recv = recv_sizes[0] ;
    for(int i=1;i<p;++i) {
      maxbuf_send = max(maxbuf_send,send_sizes[i]) ;
      maxbuf_recv = max(maxbuf_recv,recv_sizes[i]) ;
    }
    vector<unsigned char> send_buf(maxbuf_send) ;
    vector<unsigned char> recv_buf(maxbuf_recv) ;

    // First, selfcopy
    if(send_sizes[prank] > 0) {
      int loc_pack = 0 ;
      sp->pack(&send_buf[0],loc_pack, send_sizes[prank], sendSets[prank]) ;
      loc_pack = 0 ;
      op->unpack(&send_buf[0],loc_pack,recv_sizes[prank],  recvSeqs[prank]) ;
    }
    for(int i=1;i<p;++i) { // no loop to communicate
      int ps = (prank+i)%p ; // sending partner
      int pr = (prank-i+p)%p ; // receiving partner
      MPI_Request request ;
      if(recv_sizes[pr] > 0) {
        MPI_Irecv(&recv_buf[0],recv_sizes[pr],MPI_PACKED,pr,901,comm,&request) ;
      }
      if(send_sizes[ps] > 0) {
        int loc_pack = 0 ;
        sp->pack(&send_buf[0],loc_pack, send_sizes[ps], sendSets[ps]) ;
        MPI_Send(&send_buf[0],send_sizes[ps],MPI_PACKED,ps,901,comm) ;
      }
      if(recv_sizes[pr]> 0) {MPI_Status status ;
        MPI_Wait(&request,&status) ;
        int loc_pack = 0 ;
        op->unpack(&recv_buf[0],loc_pack,recv_sizes[pr],  recvSeqs[pr]) ;
      }
    }

  }

  storeRepP
  generalCommStore(// input store
                   storeRepP sp,
                   // first: from entity (in container ordering),
                   // second: to global partitioned entity map
                   const vector<std::pair<Entity,Entity> > &commMap,
                   // To entity partition
                   CPTR<partitionFunctionType> partition,
                   // mapping from global number to local numbering
                   const vector<std::pair<Entity,Entity> > &global2local,
                   // If this is null, create new container, otherwise
                   // assume it is allocated already
                   storeRepP op,
                   MPI_Comm comm) {
    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;

    // Check each processor, find out which sets to recv
    vector<Entity> sourceEntity(commMap.size()) ;
    vector<int> sourceRank(commMap.size()) ;
    for(size_t i=0;i<commMap.size();++i)
      sourceEntity[i] = commMap[i].second ;
    partition->mapKeyToRank(&sourceRank[0],&sourceEntity[0],commMap.size()) ;
    vector<pair<int,pair<Entity,Entity> > > procsort(commMap.size()) ;
    for(size_t i=0;i<commMap.size();++i) {
      procsort[i].first = sourceRank[i] ;
      procsort[i].second=commMap[i] ;
    }
    std::sort(procsort.begin(),procsort.end()) ;
    vector<int> skipsz(p,0) ;

    // Now compute where to send data to put in file ordering
    vector<entitySet> send_sets(p) ; // sets of entities to send
    vector<sequence> send_seqs(p) ;

    for(size_t i=0;i<commMap.size();++i) {
      const int p = procsort[i].first ;

      skipsz[p]++ ;
      send_sets[p] += procsort[i].second.first ;
      send_seqs[p] += procsort[i].second.second ;
    }
    vector<int> offsets(p+1,0) ;
    for(int i = 0;i<p;++i)
      offsets[i+1] = offsets[i]+skipsz[i] ;

    //Get the sequences of where we place the data when we receive it
    vector<sequence> recv_seqs = transposeSeq(send_seqs) ;

    if(global2local.size() > 0) { // now map the recv seqs to local ordering
      int Imn = 0 ;
      int Imx = global2local.size()-1 ;
      int Vmn = global2local[Imn].first ;
      int Vmx = global2local[Imx].first ;
      int delta = max((Imx-Imn+1)/(Vmx-Vmn+1),1) ;
      for(int i=0;i<p;++i) {
        sequence s = recv_seqs[i] ;
        sequence cs ;
        for(sequence::const_iterator ii=s.begin();ii!=s.end();++ii) {
          int v = *ii ;
          int imn = Imn ;
          int imx = Imx ;
          int vmn = Vmn ;
          int vmx = Vmx ;
          int is = imn + (imx-imn)*delta ;
          while(true) {
            int vs = global2local[is].first ;
            if(v == vs || imn == imx) {
              break ;
            } else if(v>vs) {
              imn = is ;
              vmn = vs ;
            } else {
              imx = is ;
              vmx = vs ;
            }
            delta = max((imx-imn+1)/(vmx-vmn+1),1) ;
            is = imn + (v-vmn)*delta ;
            if(is < imn || is > imx)
              is = (imn+imx)/2 ;
          }
          if(global2local[is].first == v)
            cs += global2local[is].second ;
          else
            cs += v ;
        }
        recv_seqs[i] = cs ;
      }
    }

    // recv_allocation
    entitySet dom ;
    for(int i=0;i<p;++i)
      dom += entitySet(recv_seqs[i]) ;

    // allocate store over shifted domain
    //    storeRepP op = 0 ; //output container
    if(op == 0) { // allocate output container
      // now allocate the container
      frame_info spfi = sp->get_frame_info() ;
      if(spfi.size == 0) { // This is a multistore
        // we need to get the counts
        // First get counts for multiStore
        store<int> counts ;
        entitySet sdom = sp->domain() ;
        counts.allocate(sdom) ;
        int loc = 0 ;
        FORALL(sdom,ii) {
          counts[ii] = spfi.first_level[loc] ;
          loc++ ;
        } ENDFORALL ;
        store<int> counts_target ;
        counts_target.allocate(dom) ;
        generalMPIComm(counts_target.Rep(),counts.Rep(),
                       send_sets,recv_seqs,comm) ;
        // Allocate multistore using framing information
        // If not one interval in dom, then we need to serialize counts
        if(dom.num_intervals() != 1) {
          vector<int> countsA(dom.size()) ;
          int cnt = 0 ;
          FORALL(dom,ii) {
            countsA[cnt++] = counts_target[ii] ;
          } ENDFORALL ;
          op = sp->new_store(dom,&countsA[0]) ;
        } else {
          op = sp->new_store(dom,&counts_target[dom.Min()]) ;
        }
      } else {
        // this is a store or storeVec
        op = sp->new_store(dom) ;
        op->set_elem_size(spfi.size) ;
      }
    }
    // Now do communication
    generalMPIComm(op, sp, send_sets, recv_seqs, comm) ;

    return op ;
  }

  void redistributeStore(Loci::storeRepP &op,
                         Loci::storeRepP sp,
                         entitySet dom,
                         const Map &toNumbering,
                         const std::vector<int> &splits,
                         MPI_Comm comm) {
    vector<std::pair<Entity,Entity> > commMap(dom.size());
    int i = 0 ;
    FORALL(dom,ii) {
      commMap[i].first = ii ;
      commMap[i].second = toNumbering[ii] ;
      i++ ;
    } ENDFORALL ;
    
    Loci::CPTR<Loci::partitionFunctionType> partition = new Loci::generalPartition(splits) ;
    vector<std::pair<Entity,Entity> > global2local ;

    op = generalCommStore(// input store
			  sp,
			  // first: from entity (in container ordering),
			  // second: to global partitioned entity map
			  commMap,
			  // To entity partition
			  partition,
			  // mapping from global number to local numbering
			  global2local,
			  op,
			  comm) ;
  }

  // collect file to global map
  entitySet
  getF2G(Map &f2g, entitySet fdom, dMap &g2f, MPI_Comm comm) {
    // First find out the distribution of the input
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int mn = fdom.Min() ;
    int mx = fdom.Max() ;
    vector<int> allmx(p) ;
    vector<int> allmn(p) ;
    MPI_Allgather(&mx,1,MPI_INT,&allmx[0],1,MPI_INT,comm) ;
    MPI_Allgather(&mn,1,MPI_INT,&allmn[0],1,MPI_INT,comm) ;
    int mnv = allmn[0] ;
    int mxv = allmx[0] ;
    vector<int> splits(p) ;
    int last = allmx[0] ;
    for(int i=0;i<p-1;++i) {
      mnv = min(mnv,allmn[i+1]) ;
      mxv = max(mxv,allmx[i+1]) ;
      splits[i]=allmx[i] ;
      if(allmx[i] < allmn[i])
        splits[i] = last ;
      else
        last = splits[i] ;
    }
    splits[p-1] = mxv+1 ;

    MapRepP g2fP = g2f ;
    entitySet FileScope = interval(mnv,mxv) ;
    entitySet dom = g2fP->preimage(FileScope).first ;

    // return maps to requestors
    vector<pair<int,int> > datalist(dom.size()) ;
    int cnt = 0 ;
    FORALL(dom,ii) {
      datalist[cnt].second = ii ; // global number
      datalist[cnt].first = g2f[ii] ; // global number
      cnt++ ;
    } ENDFORALL ;
    sort(datalist.begin(),datalist.end()) ;
    vector<int> sendszs(p,0) ;
    cnt = 0 ;
    for(size_t i=0;i<datalist.size();++i) {
      int f = datalist[i].first ;
      while(f > splits[cnt] && cnt < p)
        cnt++ ;
      sendszs[cnt]++ ;
    }
    vector<int> recvszs(p,0) ;
    MPI_Alltoall(&sendszs[0],1,MPI_INT, &recvszs[0],1,MPI_INT, comm) ;

    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + sendszs[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recvszs[i-1] ;
    }
    int recv_sz = recv_dspl[p-1] + recvszs[p-1] ;

    size_t send_sz = send_dspl[p-1] + sendszs[p-1] ;
    if(send_sz != datalist.size()) {
      cerr << "internal error in getF2G()!" << endl ;
    }
    vector<pair<int,int> > datarecv(recv_sz) ;
    for(int i=0;i<p;++i) {
      int scale = sizeof(pair<int,int>) ;
      sendszs[i] *= scale ;
      recvszs[i] *= scale ;
      send_dspl[i] *= scale ;
      recv_dspl[i] *= scale ;
    }
    MPI_Alltoallv(&datalist[0], &sendszs[0], &send_dspl[0], MPI_BYTE,
                  &datarecv[0], &recvszs[0], &recv_dspl[0], MPI_BYTE,
                  comm) ;
    entitySet rdom ;
    for(int i=0;i<recv_sz;++i)
      rdom += datarecv[i].first ;
    f2g.allocate(rdom) ;
    for(int i=0;i<recv_sz;++i)
      f2g[datarecv[i].first] = datarecv[i].second ;
    return dom ;
  }

  // collect file to global map
  entitySet
  getF2G2(Map &f2g, entitySet fdom, Map &l2g, Map &l2f, MPI_Comm comm) {
    // First find out the distribution of the input
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int mn = fdom.Min() ;
    int mx = fdom.Max() ;
    vector<int> allmx(p) ;
    vector<int> allmn(p) ;
    MPI_Allgather(&mx,1,MPI_INT,&allmx[0],1,MPI_INT,comm) ;
    MPI_Allgather(&mn,1,MPI_INT,&allmn[0],1,MPI_INT,comm) ;
    int mnv = allmn[0] ;
    int mxv = allmx[0] ;
    vector<int> splits(p) ;
    int last = allmx[0] ;
    for(int i=0;i<p-1;++i) {
      mnv = min(mnv,allmn[i+1]) ;
      mxv = max(mxv,allmx[i+1]) ;
      splits[i]=allmx[i] ;
      if(allmx[i] < allmn[i])
        splits[i] = last ;
      else
        last = splits[i] ;
    }
    splits[p-1] = mxv+1 ;

    MapRepP l2fP = l2f ;
    entitySet FileScope = interval(mnv,mxv) ;
    entitySet dom = l2fP->preimage(FileScope).first ;

    // return maps to requestors
    vector<pair<int,int> > datalist(dom.size()) ;
    int cnt = 0 ;
    FORALL(dom,ii) {
      datalist[cnt].second = l2g[ii] ; // global number
      datalist[cnt].first = l2f[ii] ; // global number
      cnt++ ;
    } ENDFORALL ;
    sort(datalist.begin(),datalist.end()) ;
    vector<int> sendszs(p,0) ;
    cnt = 0 ;
    for(size_t i=0;i<datalist.size();++i) {
      int f = datalist[i].first ;
      while(f > splits[cnt] && cnt < p)
        cnt++ ;
      sendszs[cnt]++ ;
    }
    vector<int> recvszs(p,0) ;
    MPI_Alltoall(&sendszs[0],1,MPI_INT, &recvszs[0],1,MPI_INT, comm) ;

    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + sendszs[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recvszs[i-1] ;
    }
    int recv_sz = recv_dspl[p-1] + recvszs[p-1] ;

    size_t send_sz = send_dspl[p-1] + sendszs[p-1] ;
    if(send_sz != datalist.size()) {
      cerr << "internal error in getF2G()!" << endl ;
    }
    vector<pair<int,int> > datarecv(recv_sz) ;
    for(int i=0;i<p;++i) {
      int scale = sizeof(pair<int,int>) ;
      sendszs[i] *= scale ;
      recvszs[i] *= scale ;
      send_dspl[i] *= scale ;
      recv_dspl[i] *= scale ;
    }
    MPI_Alltoallv(&datalist[0], &sendszs[0], &send_dspl[0], MPI_BYTE,
                  &datarecv[0], &recvszs[0], &recv_dspl[0], MPI_BYTE,
                  comm) ;
    entitySet rdom ;
    for(int i=0;i<recv_sz;++i)
      rdom += datarecv[i].first ;
    f2g.allocate(rdom) ;
    for(int i=0;i<recv_sz;++i)
      f2g[datarecv[i].first] = datarecv[i].second ;
    return dom ;
  }

  void File2LocalOrderGeneral(storeRepP &result, entitySet resultSet,
                              storeRepP input, int offset,
                              fact_db::distribute_infoP dist,
                              MPI_Comm comm) {
    using namespace Loci ;

    if(dist ==0 ) {
      result->copy(input,resultSet) ;
      return ;
    }
    int kd =  getKeyDomain(resultSet, dist, comm) ;

    if(kd < 0) {
      cerr << "File2LocalOrder not in single keyspace!" << endl ;
      kd = 0 ;
    }

    Map l2f ;
    l2f = dist->l2f.Rep() ;
    Map l2g ;
    l2g = dist->l2g.Rep() ;

    Map F2G ;
    entitySet fdom = input->domain() ;
    fdom = fdom >> offset ;
    // gdom is the global entitys on source processor.
    entitySet gdom = getF2G2(F2G,fdom, l2g, l2f, comm) ;

    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    vector<int> splits(p) ;
    int gmx = gdom.Max() ;
    MPI_Allgather(&gmx, 1,MPI_INT,&splits[0],1,MPI_INT,comm) ;
    for(int i=1;i<p;++i)
      if(splits[i-1]>splits[i])
        splits[i] = splits[i-1] ;

    CPTR<partitionFunctionType> partition =
      new generalPartition(splits) ;


    entitySet fsdom = F2G.domain() ;

    vector<std::pair<Entity,Entity> > commMap(fsdom.size()) ;
    int cnt = 0 ;
    FORALL(fsdom,ii) {
      commMap[cnt].first = ii-offset ;
      commMap[cnt].second = F2G[ii] ;
      cnt++ ;
    } ENDFORALL ;

    entitySet ldom = MapRepP(l2g)->preimage(gdom).first ;
    vector<std::pair<Entity,Entity> > global2local ;
    FORALL(ldom,ii) {
      global2local.push_back(pair<Entity,Entity>(l2g[ii],ii)) ;
    } ENDFORALL ;
    sort(global2local.begin(),global2local.end()) ;
    result = generalCommStore(// input store
                              input,
                              // first: from entity (in container ordering),
                              // second: to global partitioned entity map
                              commMap,
                              // To entity partition
                              partition,
                              // mapping from global number to local numbering
                              global2local,
                              // input container, if zero then allocate
                              0,
                              comm) ;

  }


  void getL2FMap(Map &l2f, entitySet dom, fact_db::distribute_infoP dist) {
    if(dist == 0) {
      l2f.allocate(dom) ;
      FORALL(dom,ii) {
        l2f[ii] = ii-dom.Min() ;
      } ENDFORALL ;
    } else {
      l2f.allocate(dom) ;

      int mnl = std::numeric_limits<int>::max() ;
      FORALL(dom,ii) {
        l2f[ii] = dist->l2f[ii] ; 
        mnl = min(mnl,l2f[ii]) ;
      } ENDFORALL ;
      int mn=mnl ;
      MPI_Allreduce(&mnl,&mn,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD) ;
      FORALL(dom,ii) {
        l2f[ii] -= mn ;
      } ENDFORALL ;
    }
  }

  void FindSimpleDistribution(entitySet dom, const Map &l2f,
                              vector<int> &splits, MPI_Comm comm) {
    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    vector<int> splits_l(p-1,0) ;
    if(p > 1) {
      int mxl = std::numeric_limits<int>::min() ;
      int mnl = std::numeric_limits<int>::max() ;
      FORALL(dom,ii) {
        mxl = max(mxl,l2f[ii]) ;
        mnl = min(mnl,l2f[ii]) ;
      } ENDFORALL ;
      int mx=mxl,mn=mnl ;
      MPI_Allreduce(&mxl,&mx,1,MPI_INT,MPI_MAX,comm) ;
      MPI_Allreduce(&mnl,&mn,1,MPI_INT,MPI_MIN,comm) ;
      //      FORALL(dom,ii) {
      // l2f[ii] -= mn ;
      //      } ENDFORALL ;
      //int dx = (mx-mn)/p+1 ;
      int dx = (mx-mn)/p ; // reh -- remove +1 rec by ed
      for(int i=0;i<p-1;++i) {
        splits_l[i] = (i+1)*dx ;
      }
      splits.swap(splits_l) ;
    }
  }

  
  void memoryBalancedDistribution(vector<int> &splits_out,
                                  const store<int> &countl,
                                  entitySet dom,
                                  const Map &toNumbering,
                                  MPI_Comm comm) {


    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int r = 0 ;
    MPI_Comm_rank(comm,&r) ;
    vector<int> splits ;
    if(p==1)
      return ;

    FindSimpleDistribution(dom, toNumbering, splits,comm) ;
    store<int> countf ;
    // Redistribute to existing partition
    vector<std::pair<Entity,Entity> > commMap(dom.size());
    int i = 0 ;
    FORALL(dom,ii) {
      commMap[i].first = ii ;
      commMap[i].second = toNumbering[ii] ;
      i++ ;
    } ENDFORALL ;

    CPTR<partitionFunctionType> partition = new generalPartition(splits) ;
    vector<std::pair<Entity,Entity> > global2local ;


    storeRepP vardist = 0 ;
    vardist = generalCommStore(// input store
                               countl.Rep(),
                               // first: from entity (in container ordering),
                               // second: to global partitioned entity map
                               commMap,
                               // To entity partition
                               partition,
                               // mapping from global number to local numbering
                               global2local,
                               vardist,
                               comm) ;

    //    RdistributeStore(vardist,countl.Rep(),dom,toNumbering,splits,comm) ;

    countf.setRep(vardist) ;
    int sumc_local = 0 ;
    entitySet locdom = countf.domain() ;
    FORALL(locdom,ii) {
      sumc_local += countf[ii] ;
    } ENDFORALL ;
    vector<int> sumc_all(p) ;
    MPI_Allgather(&sumc_local,1,MPI_INT,&sumc_all[0],1,MPI_INT,comm) ;

    vector<int> count_offsets(p,0);
    for(int i=1;i<p;++i)
      count_offsets[i] = count_offsets[i-1]+sumc_all[i-1] ;
    int sumc = count_offsets[p-1]+sumc_all[p-1] ;
    double cpp = double(sumc)/double(p) ;
    int csum_r = count_offsets[r] ;
    int p_prev = max(0,min(p-1,int(floor(double(csum_r)/cpp))));

    vector<int> lsplits ;
    FORALL(locdom,ii) {
      csum_r += countf[ii] ;
      int p_next = max(0,min(p-1,int(floor(double(csum_r)/cpp))));
      if(p_next != p_prev)
        lsplits.push_back(ii) ;
      p_prev = p_next ;
    } ENDFORALL ;
    int lsz = lsplits.size() ;
    vector<int> splitsz(p) ;
    MPI_Allgather(&lsz,1,MPI_INT,&splitsz[0],1,MPI_INT,comm) ;
    vector<int> displs(p,0) ;
    for(int i=1;i<p;++i)
      displs[i] = displs[i-1]+splitsz[i-1] ;
    if(splitsz[p-1]+displs[p-1] != p-1) {
      splits_out = splits ;
    } else {
      vector<int> splits_balanced(p-1) ;
      MPI_Allgatherv(&lsplits[0],lsplits.size(),MPI_INT,
                     &splits_balanced[0],&splitsz[0],&displs[0],
                     MPI_INT,comm) ;
      splits_out = splits_balanced ;
    }
  }

  // gather data in container
  storeRepP gatherStore(// Input Store
                        storeRepP sp,
                        // EntitySet of input to reorder
                        const std::vector<int> &commPattern,
                        // Splits for partition
                        const std::vector<int> &splits,
                        MPI_Comm comm) {

    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;

    // Now compute where to send data to put in file ordering
    vector<entitySet> recv_sets(p) ;
    vector<sequence> recv_seqs(p) ;

    // Check each processor, find out which sets to recv
    int cnt = 0 ;
    for(int i=0;i<p;++i) {
      int mxi = i<p-1?splits[i]:std::numeric_limits<int>::max() ;
      while(cnt < int(commPattern.size()) && commPattern[cnt] <= mxi) {
        recv_seqs[i] += commPattern[cnt] ;
        recv_sets[i] += cnt ;
        cnt++ ;
      }
    }

    //Get the sequences of where we place the data when we receive it
    vector<sequence> send_seqs = transposeSeq(recv_seqs) ;


    // Compute allocation domain
    entitySet file_dom = interval(0,commPattern.size()-1) ; ;
    if(commPattern.size() == 0)
      file_dom = EMPTY ;
    // allocate store over shifted domain
    storeRepP op = sp->new_store(file_dom) ;
    // Now communicate the container
    vector<int> send_sizes(p),recv_sizes(p) ;

    for(int i=0;i<p;++i)
      if(send_seqs[i].num_intervals() > 0)
        send_sizes[i] = sp->pack_size(entitySet(send_seqs[i])) ;
      else
        send_sizes[i] = 0 ;

    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 comm) ;

    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz) ;
    vector<unsigned char> recv_store(recv_sz) ;


    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      if(send_seqs[i].num_intervals() > 0)
        sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
                 entitySet(send_seqs[i])) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
                  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
                  comm) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      if(recv_sets[i] != EMPTY)
        op->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                   sequence(recv_sets[i])) ;
    }

    return op ;
  }


  // With allocation for multistore
  storeRepP gatherMultiStore(// Input Store
                             storeRepP sp,
                             // EntitySet of input to reorder
                             const std::vector<int> &commPattern,
                             // Splits for partition
                             const std::vector<int> &splits,
                             MPI_Comm comm) {

    // Get number of processors
    int p = 0 ;
    MPI_Comm_size(comm,&p) ;
    int prank = 0 ;
    MPI_Comm_rank(comm,&prank) ;


    // First get counts for multiStore
    store<int> counts ;
    entitySet dom = sp->domain() ;
    counts.allocate(dom) ;
    frame_info spfi = sp->get_frame_info() ;
    int loc = 0 ;
    FORALL(dom,ii) {
      counts[ii] = spfi.first_level[loc] ;
      loc++ ;
    } ENDFORALL ;
    storeRepP csp = gatherStore(counts.Rep(),commPattern,splits,comm) ;
    counts.setRep(csp) ;

    // Now compute where to send data to put in file ordering
    vector<entitySet> recv_sets(p) ;
    vector<sequence> recv_seqs(p) ;

    // Check each processor, find out which sets to recv
    int cnt = 0 ;
    for(int i=0;i<p;++i) {
      int mxi = i<p-1?splits[i]:std::numeric_limits<int>::max() ;
      while(cnt < int(commPattern.size()) && commPattern[cnt] <= mxi) {
        recv_seqs[i] += commPattern[cnt] ;
        recv_sets[i] += cnt ;
        cnt++ ;
      }
    }

    //Get the sequences of where we place the data when we receive it
    vector<sequence> send_seqs = transposeSeq(recv_seqs) ;


    // Compute allocation domain
    entitySet file_dom = interval(0,commPattern.size()-1) ; ;
    if(commPattern.size() == 0)
      file_dom = EMPTY ;

    // allocate store over shifted domain
    storeRepP op = sp->new_store(file_dom,&counts[0]) ;

    // Now communicate the container
    vector<int> send_sizes(p),recv_sizes(p) ;

    for(int i=0;i<p;++i)
      if(send_seqs[i].num_intervals()>0)
        send_sizes[i] = sp->pack_size(entitySet(send_seqs[i])) ;
      else
        send_sizes[i] = 0 ;

    MPI_Alltoall(&send_sizes[0],1,MPI_INT,
                 &recv_sizes[0],1,MPI_INT,
                 comm) ;

    vector<int> send_dspl(p),recv_dspl(p) ;
    send_dspl[0] = 0 ;
    recv_dspl[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_dspl[i] = send_dspl[i-1] + send_sizes[i-1] ;
      recv_dspl[i] = recv_dspl[i-1] + recv_sizes[i-1] ;
    }
    int send_sz = send_dspl[p-1] + send_sizes[p-1] ;
    int recv_sz = recv_dspl[p-1] + recv_sizes[p-1] ;

    vector<unsigned char> send_store(send_sz) ;
    vector<unsigned char> recv_store(recv_sz) ;


    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      if(send_seqs[i].num_intervals()>0)
        sp->pack(&send_store[send_dspl[i]],loc_pack, send_sizes[i],
                 entitySet(send_seqs[i])) ;
    }

    MPI_Alltoallv(&send_store[0], &send_sizes[0], &send_dspl[0], MPI_PACKED,
                  &recv_store[0], &recv_sizes[0], &recv_dspl[0], MPI_PACKED,
                  comm) ;

    for(int i=0;i<p;++i) {
      int loc_pack = 0 ;
      if(recv_sets[i] != EMPTY)
        op->unpack(&recv_store[recv_dspl[i]],loc_pack,recv_sizes[i],
                   sequence(recv_sets[i])) ;
    }

    return op ;
  }

}
