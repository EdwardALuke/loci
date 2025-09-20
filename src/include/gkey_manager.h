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
#ifndef GKEY_MANAGER_H
#define GKEY_MANAGER_H

#include <vector>
#include <Tools/intervalSet.h>
#include <Tools/cptr.h>
#include <Tools/simple_partition_long.h>
#include <gmap.h>

namespace Loci {
  class gKeyManager;
  typedef CPTR<gKeyManager> gKeyManagerP;
  
  class gKeyManager: public CPTR_type {
  public:
    virtual ~gKeyManager() {}
    // generate a set of gEntitys, whose size is passed in as an argument
    virtual gEntitySet
    generate_key(gEntity size) = 0 ;
    virtual gKeyManagerP
    clone() const = 0 ;
  };

 
  
  class distKeyManager:public gKeyManager{
    //this is a global maximum value, all process has the same value
    gEntity max_allocated; 
    MPI_Comm comm;
    int p;
    int rank;
  public:
    distKeyManager(): max_allocated(-1), comm(MPI_COMM_WORLD){
      MPI_Comm_size(comm,&p) ;
      MPI_Comm_rank(comm, &rank);
    }
    
    // generate a set of keys distributedly,
    // size: local number  of keys to generate
    //return value: local keys generated
    gEntitySet
    generate_key(gEntity size) ;
      
    gKeyManagerP
    clone() const { return new distKeyManager(*this) ;
    }
   
  };
    

} // end of namespace Loci

#endif
