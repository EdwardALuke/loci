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
#define io_performance
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi_containerIO.h>
#include <execute.h>
#include <dist_internal.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include <distribute_io.h>
#include "dist_tools.h"
#include <fact_db.h>
#include <constraint.h>
#include <string>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;


#ifndef MPI_STUBB

namespace Loci {

  static void handle_error(int errcode, string str)
  {
    char msg[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string(errcode, msg, &resultlen);
    cerr <<  str << " : " << msg << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
#define MPI_CHECK(fn) { int errcode; errcode = (fn); if (errcode != MPI_SUCCESS) handle_error(errcode, #fn ); }
  
  //only process 0 will call this function to read in header and domain
  //Notice:: this function need to be called right after file is opened, fh not updated
  //MPI_file_read() is used, it uses individual file pointer, blocking, non-collective
  void pmpi_ReadHeader( MPI_File fh,  store_header& header,  entitySet &eset)
  {
    // read in  header 
    MPI_CHECK(MPI_File_read(fh, &header, sizeof(header), MPI_BYTE, MPI_STATUS_IGNORE));
    eset = EMPTY;
  
    unsigned long    dimension = header.dom_size;
    if(dimension == 0 || !header.ordered)return; //if unordered or domain is empty, done

    //read in domain
    
    int *data = new int[dimension];
    MPI_CHECK(MPI_File_read(fh, data, dimension, MPI_INT, MPI_STATUS_IGNORE));
    eset = EMPTY;
    for(size_t i=0;i< dimension;i+=2){
      eset |= interval(data[i],data[i+1]);
    }
    delete [] data;
  }


  //only process 0 will call this function to write header and domain,
  //Notice:: this function need to be called right after file is opened, fh not updated
  //MPI_file_write() is used, it uses individual file pointer, blocking, non-collective
  void pmpi_WriteHeader(MPI_File fh, const entitySet &eset, const store_header& header )
  {
    //Domain is written as MPI_INT, may need change later

    //write header
    MPI_CHECK(MPI_File_write(fh,
                             &header, sizeof(header), MPI_BYTE,
                             MPI_STATUS_IGNORE) );
    
    if(!header.ordered) return; //if write unordered, done 

    //if write in file numbering order, write out domain
    int num_intervals = eset.num_intervals(); 
    if( num_intervals < 1) return;
    int  dimension = num_intervals*2; //size of 1D Array
      
    interval *it = new interval[num_intervals];
      
    int *data = new int[num_intervals*2];      
      
    for(int i=0;i<num_intervals;i++){
      it[i]       = eset[i];
      data[i*2]   = it[i].first;
      data[i*2+1] = it[i].second;
    }
    MPI_CHECK(MPI_File_write(fh,
                             data, dimension, MPI_INT,
                             MPI_STATUS_IGNORE) );
    delete [] data;
    delete [] it;
  }

}


#endif





