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
#include <vector>
using std::vector ;

#include <mpi.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include "dist_tools.h"
#include "dist_internal.h"
#include <rule.h>
#include <fact_db.h>
#include <constraint.h>
#include <multiMap.h>


namespace Loci {
  namespace pio{
    void read_vector_intS(hid_t group_id, const char* name, std::vector<int>& vint, int dom_size, MPI_Comm comm) {
      
      int prank=0,pnum = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      MPI_Comm_size(comm,&pnum) ;
      
      hsize_t dimension = 0 ;
      hid_t dataset = 0;
      hid_t dataspace = 0;
      if(prank == 0) {
        dataset = H5Dopen(group_id, name,H5P_DEFAULT) ;
        dataspace = H5Dget_space(dataset) ;
        H5Sget_simple_extent_dims(dataspace, &dimension, NULL) ;
      }
      int dim = dimension ;
      MPI_Bcast(&dim, 1, MPI_INT, 0, comm) ;
      int rank = 1 ;
      hid_t datatype = H5T_NATIVE_INT ;
      int total_size = dim / pnum ;
      std::vector<int> sizes = all_collect_sizes(dom_size,comm) ;
      total_size = *std::max_element(sizes.begin(), sizes.end() );
      int *tmp_int = new int[total_size] ;
      MPI_Status status ;
      if(prank != 0) {
        MPI_Recv(tmp_int, sizes[prank], MPI_INT, 0, 12, comm, &status) ;  
        for(int i = 0; i < sizes[prank]; ++i)
          vint.push_back(tmp_int[i]) ;
      } else {
        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = 0 ;
        for(int p = 0; p < pnum; ++p) {
          dimension = sizes[p] ;
          count = dimension ;
          if(dimension != 0) {
            hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
            hid_t err = H5Dread(dataset, datatype, memspace, dataspace,
                                H5P_DEFAULT, tmp_int) ;
            if(err < 0) {
              cerr << "H5Dread() failed" << endl ;
            }
            H5Sclose(memspace) ;
          }
          start += count ;
          if(p == 0) {
            for(int i = 0; i < sizes[p]; ++i) 
              vint.push_back(tmp_int[i]) ;
          } else 
            MPI_Send(tmp_int, sizes[p], MPI_INT, p, 12, comm) ;
        }
        H5Sclose(dataspace) ;
        H5Dclose(dataset) ;
      }
      delete [] tmp_int ; 
    }

    void read_vector_intP(hid_t group_id, const char* name, std::vector<int>& vint, int dom_size, MPI_Comm comm) {
#ifndef H5_HAVE_PARALLEL
      return read_vector_intS( group_id,  name, vint,  dom_size, comm);
#else
      int prank=0,pnum = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      MPI_Comm_size(comm,&pnum) ;
    
    
      hsize_t dimension = 0 ;
      hid_t dataset = 0;
      hid_t dataspace = 0;
  
      dataset = H5Dopen(group_id, name,H5P_DEFAULT) ;
      dataspace = H5Dget_space(dataset) ;
    
      int rank = 1 ;
      hid_t datatype = H5T_NATIVE_INT ;
    
      std::vector<int> sizes = all_collect_sizes(dom_size,comm) ;

      int *tmp_int = new int[dom_size] ;

   
      hsize_t start = 0 ;
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      for(int p = 0; p < prank; ++p) {
        start += sizes[p] ; 
      }
      dimension = sizes[prank] ;
      count = dimension ;
    
      //read data
    
      {
        //create memespace independently, it is OK for dimension to be 0 
        hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
        if(dimension != 0)
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        else H5Sselect_none(dataspace);
      
        hid_t xfer_plist = create_xfer_plist(Loci::hdf5_const::dxfer_coll_type);
     
        hid_t err = H5Dread(dataset, datatype, memspace, dataspace,
                            xfer_plist, tmp_int) ;
        if(err < 0) {
          cerr << "H5Dread() failed" << endl ;
        }
        H5Pclose(xfer_plist);
        H5Sclose(memspace) ;
      }
	
	
      for(int i = 0; i < sizes[prank]; ++i) 
        vint.push_back(tmp_int[i]) ;

      H5Sclose(dataspace) ;
      H5Dclose(dataset) ;
  
      delete [] tmp_int ;
#endif
    }

    void read_multi_vector_intS(hid_t group_id, const char* name, int dim,  std::vector<int>& vint, MPI_Comm comm) {
      int prank=0,pnum = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      MPI_Comm_size(comm,&pnum) ;

      hsize_t dimension = 0 ;
      hid_t dataset = 0;
      hid_t dataspace = 0;
      if(prank == 0) {
        dataset = H5Dopen(group_id, name,H5P_DEFAULT) ;
        dataspace = H5Dget_space(dataset) ;
        H5Sget_simple_extent_dims(dataspace, &dimension, NULL) ;
      }
      int rank = 1 ;
      hid_t datatype = H5T_NATIVE_INT ;
      int total_size = 0 ;
      std::vector<int> sizes = all_collect_sizes(dim,comm) ;
      total_size = *std::max_element(sizes.begin(), sizes.end());
      int *tmp_int = new int[total_size] ;
      MPI_Status status ;
      if(prank != 0) {
        MPI_Recv(tmp_int, sizes[prank], MPI_INT, 0, 12, comm, &status) ;  
        for(int i = 0; i < sizes[prank]; ++i)
          vint.push_back(tmp_int[i]) ;
      } else { 
        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = 0 ;
        for(int p = 0; p < pnum; ++p) {
          dimension = sizes[p] ;
          count = dimension ;
          if(dimension != 0) {
            hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
            hid_t err = H5Dread(dataset, datatype, memspace, dataspace,
                                H5P_DEFAULT, tmp_int) ;
            if(err < 0) {
              cerr << "H5Dread() failed" << endl ;
            }
            H5Sclose(memspace) ;
          }
          start += count ;
          if(p == 0) {
            for(int i = 0; i < sizes[p]; ++i) 
              vint.push_back(tmp_int[i]) ;
          } else 
            MPI_Send(tmp_int, sizes[p], MPI_INT, p, 12, comm) ;
        }
        H5Sclose(dataspace) ;
        H5Dclose(dataset) ;
      }
      delete [] tmp_int ; 
    }
  
    void write_vector_intS(hid_t group_id, const char* name, std::vector<int>& vint,MPI_Comm comm) {
      int prank=0,pnum = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      MPI_Comm_size(comm,&pnum) ;
      std::vector<int> sort_max(pnum) ;
      int tot_entities = vint.size() ;
      sort_max = all_collect_sizes(tot_entities,comm) ;
      std::vector<int> sizes = sort_max ;
      std::sort(sort_max.begin(), sort_max.end()) ;
      tot_entities = 0 ;
      for(int i = 0; i < pnum; ++i)
        tot_entities += sort_max[i] ;
      int *tmp_int = new int[sort_max[pnum-1]] ;
      int tmp = 0 ;
      for(std::vector<int>::iterator vi = vint.begin(); vi != vint.end(); ++vi)
        tmp_int[tmp++] = *vi ;
      if(prank != 0) {
        MPI_Status status ;
        int flag = 0 ;
        MPI_Recv(&flag,1, MPI_INT, 0, 11, comm, &status) ;
        if(flag)
          MPI_Send(tmp_int, sizes[prank], MPI_INT, 0, 12, comm) ;
      } else {
        hid_t datatype = H5T_NATIVE_INT ;
        int rank = 1 ;
        hsize_t start = 0 ;
        hsize_t stride = 1 ;
        hsize_t count = 0 ;
        hsize_t dimension = tot_entities ;
        if(dimension != 0) {
          hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
          count = sizes[0] ;
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
      
          dimension = sizes[0] ;
          start += dimension ;
          hid_t dataset = H5Dcreate(group_id, name , datatype, dataspace,
				    H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
          if(dimension != 0) {
            hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
            H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_int) ;
            H5Sclose(memspace) ;
          }
          H5Dclose(dataset) ;
		  
          for(int i = 1; i < pnum; ++i) {
            MPI_Status status ;
            int flag = 1 ;
            MPI_Send(&flag, 1, MPI_INT, i, 11, comm) ;
            MPI_Recv(tmp_int, sizes[i], MPI_INT, i, 12, comm, &status) ;
            dimension = sizes[i] ;
            count = dimension ;
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ; 
            start += count ;
            if(dimension != 0) {
              hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
              dataset = H5Dopen(group_id, name,H5P_DEFAULT) ;
              H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_int) ;
              H5Dclose(dataset) ;
              H5Sclose(memspace) ;
            }
          }
          H5Sclose(dataspace) ;
        }
      }
      delete [] tmp_int ;
    }

    void write_vector_intP(hid_t group_id, const char* name, std::vector<int>& vint,MPI_Comm comm) {
#ifndef H5_HAVE_PARALLEL
      write_vector_intS(group_id,  name,  vint, comm);
#else
    
      int prank=0,pnum = 0 ;
      MPI_Comm_rank(comm,&prank) ;
      MPI_Comm_size(comm,&pnum) ;
      std::vector<int> sort_max(pnum) ;
      int tot_entities = vint.size() ;
      sort_max = all_collect_sizes(tot_entities,comm) ;
      std::vector<int> sizes = all_collect_sizes(vint.size());
      std::sort(sort_max.begin(), sort_max.end()) ;
      tot_entities = 0 ;
      for(int i = 0; i < pnum; ++i)
        tot_entities += sort_max[i] ;
      int *tmp_int = new int[sort_max[pnum-1]] ;
      int tmp = 0 ;
      for(std::vector<int>::iterator vi = vint.begin(); vi != vint.end(); ++vi)
        tmp_int[tmp++] = *vi ;
    
    
    
      //every process write
      hid_t datatype = H5T_NATIVE_INT ;
      int rank = 1 ;
      hsize_t start = 0 ;
      hsize_t stride = 1 ;
      hsize_t count = 0 ;
      hsize_t dimension = tot_entities ;
      if(dimension != 0) {
        //create dataset collectively
        hid_t dataspace = H5Screate_simple(rank, &dimension, NULL) ;
        count = sizes[prank] ;
        for(int i = 0; i < prank; i++)start += sizes[i] ;
        if(count != 0)H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL) ;
        else H5Sselect_none(dataspace);
      
        hid_t dataset = H5Dcreate(group_id, name , datatype, dataspace,
				  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
        dimension = sizes[prank] ;
           
        //it is ok for dimension to be 0
        hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
        hid_t xfer_plist = create_xfer_plist(Loci::hdf5_const::dxfer_coll_type);
        H5Dwrite(dataset, datatype, memspace, dataspace, xfer_plist, tmp_int) ;
        H5Sclose(memspace) ;
        H5Pclose(xfer_plist);
        
        H5Dclose(dataset) ;
        H5Sclose(dataspace);
      }	  
      delete [] tmp_int ;
      
#endif
    }
  }   
  std::vector<int> all_collect_sizes(int size,MPI_Comm comm) {
    int prank=0,pnum = 0 ;
    MPI_Comm_rank(comm,&prank) ;
    MPI_Comm_size(comm,&pnum) ;
    std::vector<int> vset( pnum) ;
    if(pnum > 1) {
      int *recv_count = new int[ pnum ] ;
      int *send_count = new int[ pnum ] ;
      for(int i = 0; i <  pnum; ++i) 
	send_count[i] = size ;
      
      MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,comm) ;
      for(int i = 0; i <  pnum ; ++i)
	vset[i] = recv_count[i] ;
      
      delete [] send_count ;
      delete [] recv_count ;
    }
    else
      vset[0] = size ;
    return vset ;
  }

}
