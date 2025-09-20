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
#ifndef HDF5_READWRITE_H
#define HDF5_READWRITE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <entitySet.h>
#include <Map.h>
#include <data_traits.h>

#define DXFER_COLLECTIVE_IO 0x1  /* Collective data transfer */
#define DXFER_INDEPENDENT_IO 0x2 /* Independent data transfer */
//#define dxfer_coll_type 0x1 /* define the data transfer type as Collective IO*/

//from partest/testphdf5.h,version 1.10.2
#define FACC_DEFAULT    0x0     /* default file access type*/
#define FACC_MPIO       0x1     /* MPIO file accesss type*/

namespace Loci {
  extern MPI_Info PHDF5_MPI_Info ;

  namespace hdf5_const {
    extern const int PPN; //processes per node
    extern const int facc_type; //file access type
    extern const int dxfer_coll_type;/* define the data transfer type as Collective IO or Independent IO collectively*/
  }



  
  /*serial/parallel io
    next several functions, io data is the same to all processes
    so in serial version, only process 0 call these functions, and it is caller's responsibility to broadcast the data after reading
    if use_parallel_io:
    for parallel writing, since all structure modification procedures such as creating groups, dataset, attributes etc. need to be collective,
    all processes call these functions, only process 0 need perform writing, using independent data transfer 
    for parallel reading, to avoid read storm, only process 0 call these functions, and it is caller's responsibility to broadcast the data
  */
  
  void HDF5_WriteDomain(hid_t group_id, const entitySet &en,  MPI_Comm comm ) ;
  void HDF5_Local2Global( hid_t group_id, entitySet &eset, Map &l_g) ;
  void HDF5_ReadDomain( hid_t group_id, entitySet &eset) ;
  void HDF5_WriteVecSize(hid_t group_id, const int &size,  MPI_Comm comm  ) ;//FOR dynamic map/store, not used currently
  void HDF5_ReadVecSize(hid_t group_id, int *size ) ;//FOR dynamic map/store, not used currently
   
  

   
  /*
   * Create the appropriate File access property list
   */
  hid_t create_faccess_plist(MPI_Comm comm, MPI_Info info, int l_facc_type);

  /* create data transfer property list*/
  hid_t  create_xfer_plist(int l_xfer_type);

  /* cretae MPI_Info*/
  int create_mpi_info(MPI_Info *info);

  /*turn on parallel io, must be called by all processes collectively */
  void set_parallel_io(bool io_type);

}
#endif
