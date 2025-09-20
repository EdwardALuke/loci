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


#include "Config/conf.h"
#include <gkey_manager.h>

namespace Loci {
  gEntitySet
  distKeyManager::generate_key(gEntity size) {
    if(p>1){
      vector<gEntity> recv_buf(p) ;
      MPI_Datatype MPI_T_type = MPI_traits<gEntity>::get_MPI_type() ;
      MPI_Allgather(&size, 1, MPI_T_type, &recv_buf[0], 1, MPI_T_type, comm) ;
      gEntity local_min = max_allocated+1 ;
      for(int i = 0; i < rank; ++i)
        local_min += recv_buf[i] ;
      gEntity total_size = 0 ;
      for(int i = 0; i < p; ++i)
        total_size += recv_buf[i] ;
      gEntity local_max = local_min + size -1;
      max_allocated += total_size;
      return gEntitySet(std::pair<gEntity, gEntity>(local_min, local_max));
    }else{
      gEntity local_min = max_allocated+1 ;
      gEntity local_max = local_min + size -1;
      max_allocated += size;
      return gEntitySet(std::pair<gEntity, gEntity>(local_min, local_max));
    }
  }




}
