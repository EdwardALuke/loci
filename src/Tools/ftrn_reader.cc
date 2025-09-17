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
#include <Tools/debug.h>
#include <Tools/ftrn_reader.h>

#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>

#include <typeinfo>

namespace Loci {

fortran_binary_file::fortran_binary_file(const char *filename)
{
    error = 0 ;
    record_size = 0 ;
    if((fd = open(filename,O_RDONLY,0666)) < 0) {
        char *nerror = new char[512] ;
        snprintf(nerror,512,"can't open '%s'",filename) ;
        error = nerror ;
    }
}

fortran_binary_file::~fortran_binary_file()
{
    close(fd) ;
}

int fortran_binary_file::get_array(int *num_array,unsigned long sz)
{
    if(record_size == 0) {
        if(read(fd,(char *)&record_size,sizeof(record_size))
           != sizeof(record_size)) {
            warn(true) ;
            error = "premature EOF" ;
            return -1 ;
        }
    }
    
    sz = sz*sizeof(*num_array) ;
    warn(record_size<sz) ;
    if(record_size < sz) {
        error = "record sizes don't match request" ;
        return -1 ;
    }


    if(read(fd,(char *)num_array,sz) != (int)sz) {
        warn(true) ;
        error = "premature EOF reading array" ;
        return -1 ;
    }
    record_size -= sz ;
    if(record_size == 0)
      if(read(fd,(char *)&sz,sizeof(sz))!=sizeof(sz)) {
          warn(true) ;
          error = "premature EOF" ;
          return -1 ;
      }
    return 0 ;
}

int fortran_binary_file::get_array(float *num_array,unsigned long sz)
{
    if(record_size == 0) {
        if(read(fd,(char *)&record_size,sizeof(record_size))
           != sizeof(record_size)) {
            warn(true) ;
            error = "premature EOF" ;
            return -1 ;
        }
    }
    
    sz = sz*sizeof(*num_array) ;
    warn(record_size < sz) ;
    if(record_size < sz) {
        error = "record sizes don't match request" ;
        return -1 ;
    }

    if(read(fd,(char *)num_array,sz) != (int)sz) {
        warn(true) ;
        error = "premature EOF reading array" ;
        return -1 ;
    }
    record_size -= sz ;
    if(record_size == 0)
      if(read(fd,(char *)&sz,sizeof(sz))!=sizeof(sz)) {
          warn(true) ;
          error = "premature EOF" ;
          return -1 ;
      }
    return 0 ;
}

int fortran_binary_file::get_array(double *num_array,unsigned long sz)
{
    if(record_size == 0) {
        if(read(fd,(char *)&record_size,sizeof(record_size))
           != sizeof(record_size)) {
            warn(true) ;
            error = "premature EOF" ;
            return -1 ;
        }
    }
    
    sz = sz*sizeof(*num_array) ;
    warn(record_size < sz) ;
    if(record_size < sz) {
        error = "record sizes don't match request" ;
        return -1 ;
    }

    if(read(fd,(char *)num_array,sz) != (int)sz) {
        warn(true) ;
        error = "premature EOF reading array" ;
        return -1 ;
    }
    record_size -= sz ;
    if(record_size == 0)
      if(read(fd,(char *)&sz,sizeof(sz))!=sizeof(sz)) {
          warn(true) ;
          error = "premature EOF" ;
          return -1 ;
      }
    return 0 ;
}

}
