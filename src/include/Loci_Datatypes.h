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
#ifndef LOCI_DATATYPE_H
#define LOCI_DATATYPE_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <Tools/cptr.h>

#include <set>
#include <string>
#include <vector>
#include <sstream>
#include <stddef.h>

#ifdef RCSID
#undef RCSID
#endif
#include <hdf5.h>
#ifdef inline
#undef inline
#endif

#if(H5_VERS_MAJOR>1)
#define H5_INTERFACE_1_6_4
#define H5_INTERFACE_1_8
#endif
#if(H5_VERS_MAJOR==1 && H5_VERS_MINOR > 6)
#define H5_INTERFACE_1_6_4
#endif
#if(H5_VERS_MAJOR==1 && H5_VERS_MINOR == 6 && H5_VERS_RELEASE >3)
#define H5_INTERFACE_1_6_4
#endif
#if(H5_VERS_MAJOR==1 && H5_VERS_MINOR > 7)
#define H5_INTERFACE_1_8
#endif

namespace Loci {
  
  enum AtomType {
    BOOL=1, CHAR, UNSIGNED_CHAR, BYTE,
    SHORT, UNSIGNED_SHORT, INT, UNSIGNED, LONG, 
    UNSIGNED_LONG, LONGLONG, UNSIGNED_LONGLONG,
    FLOAT, DOUBLE, LONG_DOUBLE
  };

  class  AbstractDatatype: public CPTR_type {
  public:
    virtual hid_t  get_hdf5_type() const = 0;
    virtual std::ostream &output(std::ostream &s, const void *p) const = 0 ;
    virtual std::istream &input(std::istream &s, void *p) const = 0 ;
    virtual int bytesize() const = 0 ;
  };

  typedef CPTR<AbstractDatatype> DatatypeP ;

  class AtomicType : public AbstractDatatype {
    AtomType     atom;
  public:
    AtomicType(const AtomType &a) {atom = a;}
    hid_t        get_hdf5_type() const ;
    std::ostream &output(std::ostream &s, const void *p) const ;
    std::istream &input(std::istream &s, void *p) const ;
    int bytesize() const ;
  };

  typedef CPTR<AtomicType> AtomicDatatypeP ;

  class ArrayType : public AbstractDatatype {
  public:
    int numBytes ;
    DatatypeP type_data ;
    int rank ;
    int dimension[10] ;

    ArrayType(DatatypeP p, int sz, int in_rank, int *in_dim) {
      numBytes = sz ;
      rank = in_rank ;
      type_data = p ;
      for(int i=0;i<rank;++i)
        dimension[i] = in_dim[i] ;
    }
    hid_t  get_hdf5_type() const ;
    std::ostream &output(std::ostream &s, const void *p) const ;
    std::istream &input(std::istream &s, void *p) const ;
    int bytesize() const ;
  } ;
  
  typedef CPTR<ArrayType> ArrayDatatypeP ;

  class CompoundType : public AbstractDatatype  {
    int numBytes ;
  public:
    struct CompoundField  {
      DatatypeP      type_data ;
      std::string    name;
      size_t         offset;
    } ;
    CompoundType(int sz) { numBytes = sz; }

  private:
    std::vector<CompoundField>     type_list;

  public:
    // Member function ...
    void    insert( const std::string &name, size_t  offset,
                    const DatatypeP &p) {
      CompoundField f ;
      f.name = name ;
      f.offset = offset ;
      f.type_data = p ;
      type_list.push_back(f) ;
    }
    int     getNumFields() const { return type_list.size() ; }
    const CompoundField &getField(int i) const { return type_list[i]; }
    virtual hid_t  get_hdf5_type() const ;
    std::ostream &output(std::ostream &s, const void *p) const ;
    std::istream &input(std::istream &s, void *p) const ;
    int bytesize() const ;
  };
  typedef CPTR<CompoundType> CompoundDatatypeP ;

  template<class T> inline CompoundDatatypeP CompoundFactory(T in) {
    return new CompoundType(sizeof(T)) ;
  }

  template<class T> inline AtomicDatatypeP AtomicFactory(T in) {
    return new AtomicType(in) ;
  }

  inline ArrayDatatypeP ArrayFactory(DatatypeP dtype, int sz, 
                                          int in_rank, int *in_dim) {
   return new ArrayType(dtype, sz, in_rank, in_dim);
 }

}



#endif
