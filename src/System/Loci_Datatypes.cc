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
#include <Loci_Datatypes.h>
#include <Tools/parse.h>

#include <iostream>
#include <string>

using std::string ;
using std::cerr ;
using std::endl ;

namespace Loci {
  hid_t   AtomicType::get_hdf5_type() const {
    switch( atom ) {
    case BOOL:
      return (H5Tcopy(H5T_NATIVE_HBOOL));
    case FLOAT:
      return (H5Tcopy(H5T_NATIVE_FLOAT));
    case DOUBLE:
      return (H5Tcopy(H5T_NATIVE_DOUBLE));
    case LONG_DOUBLE:
      return (H5Tcopy(H5T_NATIVE_LDOUBLE));
    case INT:
      return (H5Tcopy(H5T_NATIVE_INT));
    case UNSIGNED:
      return (H5Tcopy(H5T_NATIVE_UINT));
    case CHAR:
      return (H5Tcopy(H5T_NATIVE_CHAR));
    case UNSIGNED_CHAR:
      return (H5Tcopy(H5T_NATIVE_UCHAR));
    case SHORT:
      return (H5Tcopy(H5T_NATIVE_SHORT));
    case UNSIGNED_SHORT:
      return (H5Tcopy(H5T_NATIVE_USHORT));
    case LONG:
      return (H5Tcopy(H5T_NATIVE_LONG));
    case UNSIGNED_LONG:
      return (H5Tcopy(H5T_NATIVE_ULONG));
    case LONGLONG:
      return (H5Tcopy(H5T_NATIVE_LLONG)) ;
    case UNSIGNED_LONGLONG:
      return (H5Tcopy(H5T_NATIVE_ULLONG)) ;
    default:
      cerr << "Unknown Basic datatype  " << atom << endl;
      abort();
    }
    return 0;
  }
  
  std::ostream  &AtomicType::output(std::ostream &s, const void *p) const {
    switch( atom ) {
    case BOOL:
      s << *(reinterpret_cast<const bool *>(p)) ;
      return s ;
    case FLOAT:
      s << *(reinterpret_cast<const float *>(p)) ;
      return s ;
    case DOUBLE:
      s << *(reinterpret_cast<const double *>(p)) ;
      return s ;
    case LONG_DOUBLE:
      {double d = *(reinterpret_cast<const long double *>(p)) ;
      s << d ;
      }
      return s ;
    case INT:
      s << *(reinterpret_cast<const int *>(p)) ;
      return s ;
    case UNSIGNED:
      s << *(reinterpret_cast<const unsigned int *>(p)) ;
      return s ;
    case CHAR:
      s << *(reinterpret_cast<const char *>(p)) ;
      return s ;
    case UNSIGNED_CHAR:
      s << *(reinterpret_cast<const unsigned char *>(p)) ;
      return s ;
    case SHORT:
      s << *(reinterpret_cast<const short *>(p)) ;
      return s ;
    case UNSIGNED_SHORT:
      s << *(reinterpret_cast<const unsigned short *>(p)) ;
      return s ;
    case LONG:
      s << *(reinterpret_cast<const long *>(p)) ;
      return s ;
    case UNSIGNED_LONG:
      s << *(reinterpret_cast<const unsigned long *>(p)) ;
      return s ;
    default:
      cerr << "Unknown Basic datatype  " << atom << endl;
      abort();
    }
    return s ;
  }

  std::istream  &AtomicType::input(std::istream &s, void *p) const {
    switch( atom ) {
    case BOOL:
      s >> *(reinterpret_cast<bool *>(p)) ;
      return s ;
    case FLOAT:
      s >> *(reinterpret_cast<float *>(p)) ;
      return s ;
    case DOUBLE:
      s >> *(reinterpret_cast<double *>(p)) ;
      return s ;
    case LONG_DOUBLE:
      {double d ;
      s >> d ;
      *(reinterpret_cast<long double *>(p)) = d ;
      return s ;
      }
    case INT:
      s >> *(reinterpret_cast<int *>(p)) ;
      return s ;
    case UNSIGNED:
      s >> *(reinterpret_cast<unsigned int *>(p)) ;
      return s ;
    case CHAR:
      s >> *(reinterpret_cast<char *>(p)) ;
      return s ;
    case UNSIGNED_CHAR:
      s >> *(reinterpret_cast<unsigned char *>(p)) ;
      return s ;
    case SHORT:
      s >> *(reinterpret_cast<short *>(p)) ;
      return s ;
    case UNSIGNED_SHORT:
      s >> *(reinterpret_cast<unsigned short *>(p)) ;
      return s ;
    case LONG:
      s >> *(reinterpret_cast<long *>(p)) ;
      return s ;
    case UNSIGNED_LONG:
      s >> *(reinterpret_cast<unsigned long *>(p)) ;
      return s ;
    default:
      cerr << "Unknown Basic datatype  " << atom << endl;
      abort();
    }
    return s ;
  }

  int   AtomicType::bytesize() const {
    switch( atom ) {
    case BOOL:
      return sizeof(bool) ;
    case FLOAT:
      return sizeof(float) ;
    case DOUBLE:
      return sizeof(double) ;
    case LONG_DOUBLE:
      return sizeof(long double) ;
    case INT:
      return sizeof(int) ;
    case UNSIGNED:
      return sizeof(unsigned int) ;
    case CHAR:
      return sizeof(char) ;
    case UNSIGNED_CHAR:
      return sizeof(unsigned char) ;
    case SHORT:
      return sizeof(short) ;
    case UNSIGNED_SHORT:
      return sizeof(unsigned short) ;
    case LONG:
      return sizeof(long) ;
    case UNSIGNED_LONG:
      return sizeof(unsigned long) ;
    default:
      cerr << "Unknown Basic datatype  " << atom << endl;
      abort();
    }
    return 0 ;
  }

  hid_t  ArrayType::get_hdf5_type() const {
    hid_t hdf5T = type_data->get_hdf5_type() ;
    hsize_t array_dims[10];

    for(int k = 0; k < rank; k++)
      array_dims[k]  = dimension[k];


    hid_t rtype = H5Tarray_create(hdf5T,rank,array_dims) ;

    H5Tclose(hdf5T) ;
    return  rtype ;
  }

  std::ostream &ArrayType::output(std::ostream &s, const void *p) const {
    int total = 1 ;
    for(int k=0;k<rank;++k)
      total *= dimension[k] ;
    int sz = type_data->bytesize() ;

    for(int i=0;i<total;++i) {
      const void *np = reinterpret_cast<const char *>(p)+sz*i ;
      type_data->output(s,np) ;
      if(i<total-1)
        s << ' ' ;
    }
    s << endl ;
    return s ;
  }
  std::istream &ArrayType::input(std::istream &s, void *p) const {
    int total = 1 ;
    for(int k=0;k<rank;++k)
      total *= dimension[k] ;
    int sz = type_data->bytesize() ;

    for(int i=0;i<total;++i) {
      void *np = reinterpret_cast<char *>(p)+sz*i ;
      type_data->input(s,np) ;
    }
    return s ;
  }

  int ArrayType::bytesize() const  {
    int total = 1 ;
    for(int k=0;k<rank;++k)
      total *= dimension[k] ;

    return total * (type_data->bytesize()) ;
  }

  hid_t CompoundType::get_hdf5_type() const
  {
    hid_t vDatatype  = H5Tcreate( H5T_COMPOUND, numBytes);
    hid_t hdf5T ;

    for(size_t i = 0; i < type_list.size(); i++) {
      hdf5T = type_list[i].type_data->get_hdf5_type() ;
      H5Tinsert(vDatatype,type_list[i].name.c_str(), type_list[i].offset,hdf5T) ;
      H5Tclose(hdf5T) ;
    }
    return vDatatype ;
  }

  std::ostream &CompoundType::output(std::ostream &s, const void *p) const {
    s << '<' ;
    for(size_t i = 0; i < type_list.size(); i++) {
      s << type_list[i].name << '=' ;
      const void *np = reinterpret_cast<const char *>(p)+type_list[i].offset ;
      type_list[i].type_data->output(s,np) ;
      if(i+1<type_list.size())
        s << ',' ;
    }
    s << '>' ;
    return s ;
  }

  std::istream &CompoundType::input(std::istream &s, void *p) const {
    parse::kill_white_space(s) ;
    size_t i = 0 ;
    if(s.peek() == '<') {
      s.get() ;
      while(true) {
        string name = parse::get_name(s) ;
        size_t j =0;
        while(type_list[i].name != name && j!=type_list.size()) {
          ++i;
          i=i%type_list.size() ;
          ++j;
        }
        if(type_list[i].name != name) {
          cerr << "unable to find " << name << " in typelist for CompoundType::input" << endl ;
        }
        parse::kill_white_space(s) ;
        if(s.peek() != '=') {
          cerr << "syntax error in CompoundType::input()" << endl ;
          return s ;
        }
        s.get() ;
        parse::kill_white_space(s) ;
        void *np = reinterpret_cast<char *>(p)+type_list[i].offset ;
        type_list[i].type_data->input(s,np) ;
        parse::kill_white_space(s) ;
        
        if(s.peek() == ',') {
          s.get() ;
          continue ;
        }
        if(s.peek() == '>') {
          s.get();
          return s;
        }
        cerr << "Syntax error in CompoundType::input()" << endl ;
        return s ;
      }
    } else {
      cerr << "Syntax error in CompoundType::input()" << endl ;
      return s ;
    }
    //    return s ;
  }

  int CompoundType::bytesize() const {
    return numBytes ;
  }
}


