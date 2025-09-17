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
#ifndef FRAME_INFO_H
#define FRAME_INFO_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

namespace Loci {
  //frame info is used in hdf5 file_io of containers 
  struct frame_info {
    //is the data user defined type with variable size,
    //0: No, IDENTITY_CONVERTER, size = 1
    //1: Yes, USER_DEFINED_CONVERTER, size = 1
    int is_stat ; 
    int size ;
    std::vector<int> first_level ;
    std::vector<int> second_level ;
    frame_info() {
      is_stat = 0 ;
      size = 0 ;
    }
    frame_info(int a , int b) {
      is_stat = a ;
      size = b ;
      
    }
    frame_info(int a , int b, std::vector<int> c, std::vector<int> d) {
      is_stat = a ;
      size = b ;
      
      first_level = c ;
      second_level = d ;
    } 
    frame_info(const frame_info &fi) { 
      is_stat = fi.is_stat ;
      size = fi.size ;
      if(!size)
	first_level = fi.first_level ;
      if(is_stat) 
	second_level = fi.second_level ;
    }
    
    frame_info &operator = (const frame_info &fi) { 
      is_stat = fi.is_stat ;
      size = fi.size ;
      if(!size) 
	first_level = fi.first_level ;
      if(is_stat)
	second_level = fi.second_level ;
      
      return *this ;
    }
  } ;

  
 
}

#endif
