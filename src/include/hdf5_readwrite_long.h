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
#ifndef HDF5_READWRITE_LONG_H
#define HDF5_READWRITE_LONG_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>
#include <vector>
#include <entitySet.h>
#include <Map.h>
#include <data_traits.h>
#include <ostream>
#include <fstream>
#include <iostream>
namespace Loci {
  extern std::ofstream debugout ;
  //writing will check sizeof(T) and decide datatype 
  template<class T> void g_HDF5_WriteDomain(hid_t group_id, const genIntervalSet<T> &en ){
    hsize_t dimension=0;
    
    int rank = 1;     
    size_t num_intervals = en.num_intervals(); 

    hid_t dataspace, dataset;
    hid_t datatype  = H5T_NATIVE_INT;
    if(sizeof(T) != sizeof(int)) datatype = HDF5_GENTITY_TYPE;
    if( num_intervals < 1) return;

    dimension = num_intervals*2; //size of 1D Array
    if(dimension == 0) return;
    dataspace = H5Screate_simple(rank, &dimension, NULL);
    dataset   = H5Dcreate(group_id, "Interval Set", datatype, dataspace,
			  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    std::vector<std::pair<T, T> > it(num_intervals);

    std::vector<T> data(num_intervals*2);      

    for(size_t i=0;i<num_intervals;i++){
      it[i]       = en[i];
      data[i*2]   = it[i].first;
      data[i*2+1] = it[i].second;
    }

    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
  // template<class T> void HDF5_Local2Global( hid_t group_id, genIntervalSet<T> &eset, Map &l_g) ;
  
  //reading will check datatype first instead of assuming it is H5T_NATIVE_INT,
  //if types do not match, read into int_type vector and type cast into T
  template<class T> void g_HDF5_ReadDomain( hid_t group_id, genIntervalSet<T> &eset){
    hsize_t    dimension;
    hid_t      dataset, dataspace;

    H5Eset_auto (H5E_DEFAULT,NULL, NULL);

    eset = genIntervalSet<T>::EMPTY;
    dataset  = H5Dopen(group_id, "Interval Set",H5P_DEFAULT);
    if( dataset > 0 ) {
      dataspace  = H5Dget_space(dataset);
      H5Sget_simple_extent_dims (dataspace, &dimension, NULL);
     hid_t datatype = H5Dget_type (dataset);  /* datatype identifier */
     
    
     
     if(datatype == H5T_NATIVE_INT){
       std::vector<int_type> tmp_data(dimension);
       H5Dread( dataset, datatype, H5S_ALL, dataspace,
                H5P_DEFAULT, &tmp_data[0]);
        
       eset =genIntervalSet<T>::EMPTY;
       for(size_t i=0;i< dimension;i+=2){
         eset |= std::pair<T, T>(T(tmp_data[i]),T(tmp_data[i+1]));
       }
     }else{//assume type match
       std::vector<T> data(dimension);
       H5Dread( dataset, datatype, H5S_ALL, dataspace,
                H5P_DEFAULT, &data[0]);
       eset =genIntervalSet<T>::EMPTY;
       for(size_t i=0;i< dimension;i+=2){
         eset |= std::pair<T, T>(data[i],data[i+1]);
       }
     }
     H5Sclose(dataspace);
     H5Dclose(dataset);
    }
  }
  //writing will check sizeof(T) and decide datatype 
  template<class T> void g_HDF5_WriteVecSize(hid_t group_id, const T &size ){
   hsize_t dimension = 1;
   int     rank = 1;
   if(dimension == 0) return;
   hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
   hid_t vDatatype  = H5T_NATIVE_INT;
   if(sizeof(T) != sizeof(int)) vDatatype = HDF5_GENTITY_TYPE;
   hid_t vDataset   = H5Dcreate(group_id, "VecSize", vDatatype, vDataspace,
                                H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
   
   H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &size);
   H5Sclose( vDataspace );
   H5Dclose( vDataset   );
  }
  
  template<class T> void g_HDF5_ReadVecSize(hid_t group_id, T *size ){
    int     rank=1;
    hsize_t dimension = 1;

    H5Eset_auto (H5E_DEFAULT,NULL, NULL);

    *size = 0;
  
    hid_t vDataset   = H5Dopen( group_id, "VecSize",H5P_DEFAULT);

    if( vDataset > 0) {
      if(dimension == 0) return;
      hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
      hid_t vDatatype  = H5Dget_type (vDataset);  /* datatype identifier */
      
      if( vDatatype == H5T_NATIVE_INT){
        int_type tmp_size = 0;
        H5Dread( vDataset, vDatatype, H5S_ALL, vDataspace, H5P_DEFAULT, &tmp_size);
        *size = T(tmp_size);
      }else{
        H5Dread( vDataset, vDatatype, H5S_ALL, vDataspace, H5P_DEFAULT, size);
      }
      H5Sclose( vDataspace );
      H5Dclose( vDataset   );
    }
  }

}
#endif
