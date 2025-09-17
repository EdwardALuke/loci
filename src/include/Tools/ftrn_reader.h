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
#ifndef FTRN_READER_H
#define FTRN_READER_H 1

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

namespace Loci {
class fortran_binary_file {
    int fd ;
    unsigned long record_size ;
    const char *error ;
  public:
    fortran_binary_file(const char *filename);
    ~fortran_binary_file() ;
    const char *error_msg() { return error; }
    int get_array(int *num_array, unsigned long  sz) ;
    int get_array(float *num_array,unsigned long  sz) ;
    int get_array(double *num_array,unsigned long  sz) ;
} ;
}

#endif
