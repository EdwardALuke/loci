//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef FLOWTYPES_H
#define FLOWTYPES_H

#include <Loci>

#include <iostream>
#include <string>
#include <sstream>
#include <ctype.h>

#include <Tools/tools.h>
#include <Tools/parse.h>
#include <Tools/unit_type.h>


namespace flowPsi {
  
  typedef Loci::real_t real ;
  typedef float realF ;

  typedef unsigned char byte_t ;
  
  
  using Loci::vector3d ;
  using Loci::tensor3d ;
  using Loci::norm ;
  using Loci::dot ;
  using Loci::cross ;
  using Loci::options_list ;

  using Loci::rigid_transform ;
  using Loci::periodic_info ;
  typedef Loci::vector3d<real> vect3d ;
  typedef Loci::tensor3d<real> tens3d ;

  // Used to do a priority join using pairs
  template <class T> struct priority_joiner {
    void operator()(T &r, const T &s) {
      if(r.first < s.first)
        r = s ;
    }
  } ;
    

  

  template <class T> class tmp_array {
    int sz ;
    T data[25] ;
    T * p ;
    void alloc(int size) {
      sz = size ;
      p = data ;
      if(sz > 25)
        p = new T[sz] ;
    }
    void free() {
      if(sz > 25)
        delete[] p ;
    }
    tmp_array() { alloc(0) ; }
  public:
    tmp_array(int size) {
      alloc(size) ;
    }
    tmp_array(const tmp_array &ta) {
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
    }
    tmp_array &operator=(const tmp_array &ta) {
      free() ;
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
      return *this ;
    }
    ~tmp_array() { free(); }
    T & operator[](int i) { return p[i] ; }
    T & operator[](int i) const { return p[i] ; }
    operator T *() { return p ; }
    operator const T *() const { return p ; }
  } ;
      
}

#endif
