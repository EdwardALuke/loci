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
#ifndef DATA_TRAITS_H_
#define DATA_TRAITS_H_

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#ifndef NO_OFFSETOF
#include <stddef.h>
#endif

#include <vector>
#include <utility>
#include <string>

#include <entitySet.h>

#include <Loci_Datatypes.h>
#include <Tools/expr.h>
#include <mpi.h>

namespace Loci {

  //-----------STD pair-------------------------------//
  template<class T1,class T2> inline std::ostream &
    operator<<(std::ostream &s, const std::pair<T1,T2> &v) {
    s<<"["<<v.first<<","<<v.second<<"]";
    return s;
  }

  template<class T1,class T2> inline std::istream &
    operator>>(std::istream &s, std::pair<T1,T2> &i) {
    char ch ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!='[') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a '[' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> i.first ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=',') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a ',' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> i.second ;
    
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=']') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a ']' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    return s;
  }

  //-----------STD pair unsigned int specializations---------------------//
  template<> inline std::ostream &
    operator<<(std::ostream &s, const std::pair<unsigned int,unsigned int> &v) {
    s<<"["<<v.first<<","<<v.second<<"]";
    return s;
  }

  template<> inline std::istream &
    operator>>(std::istream &s, std::pair<unsigned int,unsigned int> &i) {
    char ch ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!='[') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a '[' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> i.first ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=',') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a ',' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> i.second ;
    
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=']') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a ']' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    return s;
  }

  template<class T> inline std::ostream &operator << (std::ostream &s,
                                                      const std::vector<T> &v){
    s << v.size() ;
    for(typename std::vector<T>::const_iterator i=v.begin();i!=v.end(); ++i) 
      s << ' ' << *i ;
    s << std::endl ;
    return s ;
  }

  template<class T> inline std::istream &operator >> (std::istream &s,
                                                      std::vector<T> &v){
    v.clear() ;
    int sz ;
    s >> sz ;
    for(int i=0;i<sz;++i) {
      T val ;
      s >> val ;
      v.push_back(val) ;
    }
    return s ;
  }


  template <class T> 
  class data_schema_traits { };

  template<class T> class StringStreamConverter {
    T &ref ;
  public:
    StringStreamConverter(T &iref) : ref(iref) {}
    int getSize() const {
      std::ostringstream oss ;
      oss << ref ;
      return oss.str().size() ;
    }
    void getState(char *buf, int &size) const {
      std::ostringstream oss ;
      oss << ref ;
      const std::string &s = oss.str() ;
      size = s.size() ;
      for(int i=0;i<size;++i)
        buf[i] = s[i] ;
    }
    void setState(const char *buf, int size) {
      std::string s ;
      for(int i=0;i<size;++i)
        s += buf[i] ;
      std::istringstream iss(s) ;
      iss >> ref ;
    }
  } ;
    
  template<class T> class UndefinedConverter {
    T &ref ;
  public:
    UndefinedConverter(T &iref) : ref(iref) {}
    int getSize() const {
      std::cerr << "undefined converter" << std::endl ;
      return 0 ;
    }
    void getState(char *buf, int &size) const {
      std::cerr << "undefined converter" << std::endl ;
    }
    void setState(const char *buf, int size) {
      std::cerr << "undefined converter" << std::endl ;
    }
  } ;
    
    
  struct IDENTITY_CONVERTER{};
  struct USER_DEFINED_CONVERTER{};

  template<class T> inline std::istream &
  streaminput_SEL(IDENTITY_CONVERTER,T *v, int sz, std::istream &s) {
    typedef data_schema_traits<T> traits_type;

    DatatypeP dtype = traits_type::get_type();
    for(int i=0;i<sz;++i) {
      dtype->input(s,(void *)(v+i)) ;
    }
    return s ;
  }

  template<class T> inline std::istream &
  streaminput_SEL(USER_DEFINED_CONVERTER,T *v, int sz, std::istream &s) {
    for(int i=0;i<sz;++i) {
      s >> v[i] ; 
    }
    return s ;
  }
  
  template<class T> inline std::istream &
    streaminput(T *v,int sz, std::istream &s) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return streaminput_SEL(schema_converter(), v, sz, s) ;
  }

  template<class T> inline std::ostream &
  streamoutput_SEL(IDENTITY_CONVERTER,const T *v, int sz, std::ostream &s) {
    typedef data_schema_traits<T> traits_type;

    DatatypeP dtype = traits_type::get_type();
    for(int i=0;i<sz;++i) {
      dtype->output(s,(void *)(v+i)) ;
      s << std::endl ;
    }
    return s ;
  }

  template<class T> inline std::ostream &
  streamoutput_SEL(USER_DEFINED_CONVERTER,const T *v, int sz, std::ostream &s) {
    for(int i=0;i<sz;++i) {
      s << v[i] ;
      s << std::endl ;
    }
    return s ;
  }
  
  template<class T> inline std::ostream &
    streamoutput(const T *v,int sz, std::ostream &s) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return streamoutput_SEL(schema_converter(), v, sz, s) ;
  }

  //**************************************************************************

  template <class T> 
  class data_schema_converter_traits { };

  template <>
  struct data_schema_traits<char> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(CHAR)) ;
    }
  };

  template <>
  struct data_schema_traits<unsigned char> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED_CHAR)) ;
    }
  };

  template <>
  struct data_schema_traits<short> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(SHORT)) ;
    }
  };

  template <>
  struct data_schema_traits<int> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(INT)) ;
    }
  } ;
    
  template <>
  struct data_schema_traits<unsigned> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED)) ;
    }
  };

  template <>
  struct data_schema_traits<unsigned short> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED_SHORT)) ;
    }
  };

  template <>
  struct data_schema_traits<long> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED)) ;
    }      
  };

  template <>
  struct data_schema_traits<unsigned long> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED_LONG)) ;
    }      
  } ;

  template <>
  struct data_schema_traits<long long> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(LONGLONG)) ;
    }      
  } ;

  template <>
  struct data_schema_traits<unsigned long long> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED_LONGLONG)) ;
    }      
  } ;

  template <>
  struct data_schema_traits<float> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(FLOAT)) ;
    }      
  };

  template <>
  struct data_schema_traits<double> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(DOUBLE)) ;
    }      
  };

  template <>
  struct data_schema_traits<long double> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(LONG_DOUBLE)) ;
    }      
  };

  template <>
  struct data_schema_traits<bool> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(BOOL)) ;
    }      
  };


  class entitySetSchemaConverter {
    // For the schema converter, we always store a reference to the object
    // we are converting schmata for.
    entitySet &eref ;
  public:
    explicit entitySetSchemaConverter(entitySet &iset): eref(iset) {}
    int getSize() const {
      return eref.size() ;
    }
    void getState(int *buf, int &size) const {
      size = getSize() ;
      int ii=0; 
      for(entitySet::const_iterator i=eref.begin();i!=eref.end();++i)
        buf[ii++] = *i ;
    }
    void setState(const int *buf, int size) {
      eref = EMPTY ;
      for(int i=0;i<size;++i)
        eref += buf[i] ;
    }
  } ;
  
  template <>
  struct data_schema_traits<entitySet> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef int Converter_Base_Type ;
    typedef entitySetSchemaConverter Converter_Type ;
  };

  template<class T> class vectorSchemaConverter {
    // For the schema converter, we always store a reference to the object
    // we are converting schmata for.
    std::vector<T> &eref ;
  public:
    explicit vectorSchemaConverter(std::vector<T> &iset): eref(iset) {}
    int getSize() const {
      return eref.size() ;
    }
    void getState(T *buf, int &size) const {
      size = getSize() ;
      for(int i=0;i<size;++i)
        buf[i] = eref[i] ;
    }
    void setState(const T *buf, int size) {
      eref.clear() ;
      eref.reserve(size) ;
      for(int i=0;i<size;++i)
        eref.push_back(buf[i]) ;
    }
  } ;

  template<class T> struct data_schema_traits<std::vector<T> > {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef T Converter_Base_Type ;
    typedef vectorSchemaConverter<T> Converter_Type ;
  } ;

  class std_string_schema_converter {
    std::string &ref ;
  public:
    explicit std_string_schema_converter(std::string &iref): ref(iref) {}
    int getSize() const {
      return ref.size() ;
    }
    void getState(char *buf, int &size) const {
      size = getSize() ;
      for(int i=0;i<size;++i)
        buf[i] = ref[i] ;
    }
    void setState(const char *buf, int size) {
      ref = "" ;
      for(int i=0;i<size;++i)
        ref += buf[i] ;
    }
  } ;

  template<> struct data_schema_traits<std::string> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef std_string_schema_converter Converter_Type ;
  } ;

  template<> struct data_schema_traits<exprP> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<exprP> Converter_Type ;
  } ;

  template<class T> inline DatatypeP getLociType(const T &in) {
    return data_schema_traits<T>::get_type() ;
  }

#ifdef NO_OFFSETOF
#define LOCI_INSERT_TYPE(ct,type,variable) \
{ type X ; size_t offset = reinterpret_cast<char *>(&(X.variable)) - reinterpret_cast<char *>(&X) ;\
   ct->insert(# variable, offset,getLociType(X.variable));}
#else
#define LOCI_INSERT_TYPE(ct,type,variable) \
{ type X ; size_t offset = offsetof(type,variable) ;\
   ct->insert(# variable, offset,getLociType(X.variable));}
#endif

  template<class T, class S> struct data_schema_traits<std::pair<T,S> > {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      typedef std::pair<T,S> pair ;
      pair p ;
      CompoundDatatypeP ct = CompoundFactory(p) ;
      LOCI_INSERT_TYPE(ct,pair,first) ;
      LOCI_INSERT_TYPE(ct,pair,second) ;
      return DatatypeP(ct) ;
    }
  } ;

  template <class T> class MPI_traits {} ;


  template<> class MPI_traits <unsigned char> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_UNSIGNED_CHAR ; } ;
  } ;

  template<> class MPI_traits <char> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_BYTE ; } ;
  } ;

  template<> class MPI_traits <short> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_SHORT ; } ;
  } ;

  template<> class MPI_traits <unsigned int> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_UNSIGNED ; } ;
  } ;



  template<> class MPI_traits <int> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_INT ; } ;
  } ;

  template<> class MPI_traits <long> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_LONG ; } ;
  } ;

  template<> class MPI_traits <unsigned long> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_UNSIGNED_LONG ; } ;
  } ;

  template<> class MPI_traits <long long> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_LONG_LONG_INT; } ;
  } ;

  template<> class MPI_traits <float> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_FLOAT ; } ;
  } ;

  template<> class MPI_traits <double> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_DOUBLE ; } ;
  } ;

  template<> class MPI_traits <long double> {
  public:
    static MPI_Datatype get_MPI_type() { return MPI_LONG_DOUBLE ; } ;
  } ;



  template <class T> class HDF5_traits {} ;

  template<> class HDF5_traits <char> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_CHAR ; } ;
  } ;

  template<> class HDF5_traits <short> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_SHORT ; } ;
  } ;

  template<> class HDF5_traits <int> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_INT ; } ;
  } ;

  template<> class HDF5_traits <long> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_LONG ; } ;
  } ;

  template<> class HDF5_traits <long long> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_LLONG; } ;
  } ;


  template<> class HDF5_traits <unsigned char> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_UCHAR ; } ;
  } ;

  template<> class HDF5_traits <unsigned short> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_USHORT ; } ;
  } ;

  template<> class HDF5_traits <unsigned int> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_UINT ; } ;
  } ;

  template<> class HDF5_traits <unsigned long> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_ULONG ; } ;
  } ;

  template<> class HDF5_traits <unsigned long long> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_ULLONG; } ;
  } ;

  template<> class HDF5_traits <float> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_FLOAT ; } ;
  } ;

  template<> class HDF5_traits <double> {
  public:
    static hid_t get_HDF5_type() { return H5T_NATIVE_DOUBLE ; } ;
  } ;
    
}

#endif
