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
#ifndef DISTRIBUTED_INVERSE_LONG_H
#define DISTRIBUTED_INVERSE_LONG_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <typeinfo>
#include <vector>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <string>
#include <mpi.h>
#include <field_sort.h>
#include <fact_db.h>
namespace Loci {
  
  template<class T> void distributed_inverse_map(std::vector<T> &recv_store,
                                                 std::vector<std::pair<T,T> > &input, 
                                                 const std::vector<genIntervalSet<T> > &init_ptn){
    MPI_Comm comm = get_exec_comm() ;
   
    // Sort input according to second field
    sort(input.begin(),input.end(),field_sort2<T>) ;

    // Now count what we will be sending
    vector<int> send_sz(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      send_sz[i] = 0 ;
    int current_p = 0 ;
    for(size_t i=0;i<input.size();++i) {
      T to = input[i].second ;

      if(!init_ptn[current_p].inSet(to))
        for(int j=0;j<MPI_processes;++j)
          if(init_ptn[j].inSet(to)) {
            current_p = j ;
            break ;
          }
      send_sz[current_p] += 2 ;
    }

 
    // transfer recv sizes
    vector<int> recv_sz(MPI_processes) ;
    MPI_Alltoall(&send_sz[0],1,MPI_INT,
                 &recv_sz[0],1,MPI_INT,
                 comm) ;
    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<MPI_processes;++i) {
      size_send += send_sz[i] ;
      size_recv += recv_sz[i] ;
    }
 
    T *send_store = new T[size_send] ;
    recv_store.resize(size_recv) ;
    int *send_displacement = new int[MPI_processes] ;
    int *recv_displacement = new int[MPI_processes] ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sz[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sz[i-1] ;
    }
 
    current_p = 0 ;
    vector<int> offsets(MPI_processes) ;
    for(int i=0;i<MPI_processes;++i)
      offsets[i] = 0 ;
    for(size_t i=0;i<input.size();++i) {
      T to = input[i].second ;
      T from = input[i].first ;
      if(!init_ptn[current_p].inSet(to))
        for(int j=0;j<MPI_processes;++j)
          if(init_ptn[j].inSet(to)) {
            current_p = j ;
            break ;
          }
      send_store[send_displacement[current_p]+offsets[current_p]++] = to ;
      send_store[send_displacement[current_p]+offsets[current_p]++] = from ;
    }
 
    MPI_Datatype MPI_T_type = MPI_traits<T>::get_MPI_type() ;
    MPI_Alltoallv(&send_store[0], &send_sz[0], send_displacement, MPI_T_type,
		  &recv_store[0], &recv_sz[0], recv_displacement, MPI_T_type,
		  comm) ;
 
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] send_store ;
  }

}

#endif
