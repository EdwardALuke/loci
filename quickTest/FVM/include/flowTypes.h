//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
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

#ifdef LOCI_V5
using Loci::gKeySpaceP ;
using Loci::GEMPTY ;
using Loci::gParam ;
#endif

namespace flowPsi {
  
  typedef Loci::real_t real ;
#ifdef USE_AUTODIFF
  typedef Loci::real_t realF ;
#else
  typedef float realF ;
#endif
  typedef realF real_fj ;

  using Loci::realToDouble ;
  using Loci::realToFloat ;

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
    
  class list_input {
  public:
    std::string namelist ;
  } ;

  inline std::ostream & operator <<(std::ostream &s, const list_input &n)
  {
    s << n.namelist << std::endl ;
    return s ;
  }

  inline std::istream &operator>>(std::istream &s, list_input &n)
  {
    n.namelist = string() ;
    while(s.peek() != '\n' && s.peek() != '\r' && s.peek() != EOF) {
      char c = s.get() ;
      if(c != ' ' && c != '\t') {
        n.namelist += c ;
      }
    }
    return s ;
  }


}

namespace Loci {

   class list_input_schema_converter {
    flowPsi::list_input &ref ;
  public:
    explicit list_input_schema_converter(flowPsi::list_input &iref): ref(iref) {}
    int getSize() const {
      return ref.namelist.size() ;
    }
    void getState(char *buf, int &size) {
      size = getSize() ;
      for(int i=0;i<size;++i)
        buf[i] = ref.namelist[i] ;
    }
    void setState(char *buf, int size) {
      ref.namelist = "" ;
      for(int i=0;i<size;++i)
        ref.namelist += buf[i] ;
    }
  } ;

  template<> struct data_schema_traits<flowPsi::list_input> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef list_input_schema_converter Converter_Type ;
  } ;
}

namespace flowPsi {

  
  template <unsigned int n> class vec : public Loci::Array<real,n> {
  } ;

  struct symmetricTensor {
    real xx, xy, xz, yy, yz, zz ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, const symmetricTensor &ss)
  {
    s << ss.xx << ' ' << ss.xy << ' ' << ss.xz << ' ' <<
      ss.yy << ' ' << ss.yz << ' ' << ss.zz << std::endl ; 
    return s ;
  }

  inline std::istream &operator>>(std::istream &s, symmetricTensor &ss)
    {
      s >> ss.xx >> ss.xy >> ss.xz >> ss.yy >> ss.yz >> ss.zz ;
      return s ;
    }

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
      
  struct Quaternion {
    real x,y,z,w ;
    Quaternion() {}
    Quaternion(real xi, real yi, real zi, real wi):x(xi),y(yi),z(zi),w(wi) {}
    Quaternion(vect3d axis, real angle) {
      real sinAngle; 
      angle *= 0.5; 
      axis *= 1.0/(norm(axis)+1e-30) ; 
      sinAngle = sin(angle);
      x = (axis.x * sinAngle); 
      y = (axis.y * sinAngle); 
      z = (axis.z * sinAngle); 
      w = cos(angle);
    }
    Quaternion operator*(const Quaternion &q) const {
      vect3d vector1(x,y,z), vector2(q.x,q.y,q.z); 

      const real angle = ((w * q.w) - (dot(vector1, vector2))); 
      const vect3d across = cross(vector1, vector2);
      vector1 *= q.w ;
      vector2 *= w ;
      Quaternion result; 
      result.x = (vector1.x + vector2.x + across.x); 
      result.y = (vector1.y + vector2.y + across.y); 
      result.z = (vector1.z + vector2.z + across.z); 
      result.w = angle;
      return result ;
    }
    Quaternion &Normalize() {
      // reciprocal of the l2 norm 
      const real rl2 = 1.0 / sqrt((x*x) + (y*y) + (z*z) + (w*w)) ;
      x*=rl2 ;
      y*=rl2 ;
      z*=rl2 ;
      w*=rl2 ;
      return *this ;
    }
    
    Quaternion Inverse() const {
      Quaternion result = *this ;
      result.x *= -1 ; 
      result.y *= -1 ; 
      result.z *= -1 ; 
      return result ;
    }
    vect3d operator*(const vect3d &v) const {
      const Quaternion Q = *this ;
      const Quaternion Qinv = Inverse() ;
      const Quaternion vQ(v.x,v.y,v.z,0) ;
      Quaternion result = vQ*Qinv ;
      result = Q*result ;
      return vect3d(result.x,result.y,result.z) ;
    }
  } ;
}

namespace Loci {

  template<unsigned int n> struct data_schema_traits<flowPsi::vec<n> > {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      return getLociType(Loci::Array<flowPsi::real,n>()) ;
    }
  } ;

  template<> struct data_schema_traits<flowPsi::symmetricTensor> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::symmetricTensor()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::symmetricTensor,xx) ;
      LOCI_INSERT_TYPE(ct,flowPsi::symmetricTensor,xy) ;
      LOCI_INSERT_TYPE(ct,flowPsi::symmetricTensor,xz) ;
      LOCI_INSERT_TYPE(ct,flowPsi::symmetricTensor,yy) ;
      LOCI_INSERT_TYPE(ct,flowPsi::symmetricTensor,yz) ;
      LOCI_INSERT_TYPE(ct,flowPsi::symmetricTensor,zz) ;
      return DatatypeP(ct) ;
    }
  } ;

  template<> struct data_schema_traits<flowPsi::Quaternion> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::Quaternion()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Quaternion,x) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Quaternion,y) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Quaternion,z) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Quaternion,w) ;
      return DatatypeP(ct) ;
    }
  } ;

}


namespace flowPsi {
  inline double get_val_in_units(std::istream &s, const char *units) {
    double val ;
    s >> val ;
    while(s.peek() == ' ' || s.peek() == '\t')
      s.get() ;
    if(isalpha(s.peek())) {
      std::string units_in ;
      while(s.peek() != EOF && s.peek() != '\n' && s.peek() !='\r')
        units_in += s.get() ;
      Loci::UNIT_type units_value(Loci::UNIT_type::MKS,"general",val,units_in) ;
      if(!units_value.is_compatible(units)) {
        std::ostringstream oss ;
        oss << "Expecting units compatible with " << units
            << "!" << std::endl
            << "incompatible with input of '"<<units<<"'" << std::endl ;
        throw Loci::StringError(oss.str()) ;
      }
      val = units_value.get_value_in(units) ;
    }
    return val ;
  }

  struct TimeValue {
    double val ;
    operator double() const { return val ; }
    TimeValue &operator=(const TimeValue &v) { val = v.val ; return *this ;}
    TimeValue &operator=(double v) { val = v ; return *this; }
    TimeValue &operator=(float v) { val = v ; return *this; }
    TimeValue &operator=(int v) { val = v ; return *this; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const TimeValue &ss)
  {
    s << ss.val << " second" << std::endl ;
    return s ;
  }
    
  inline std::istream &operator>>(std::istream &s, TimeValue &ss) {
    ss.val = get_val_in_units(s,"second") ;
    return s ;
  }        

  struct TemperatureValue {
    double val ;
    operator double() const { return val ; }
    TemperatureValue &operator=(const TemperatureValue  &v) { val = v.val ; return *this; }
    TemperatureValue &operator=(double v) { val = v ; return *this; }
    TemperatureValue &operator=(float v) { val = v ; return *this; }
    TemperatureValue &operator=(int v) { val = v ; return *this; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const TemperatureValue &ss)
  {
    s << ss.val << " kelvin" << std::endl ;
    return s ;
  }
    
  inline std::istream &operator>>(std::istream &s, TemperatureValue &ss) {
    ss.val = get_val_in_units(s,"kelvin") ;
    return s ;
  }        

  struct PressureValue {
    double val ;
    operator double() const { return val ; }
    PressureValue &operator=(const PressureValue &v) { val = v.val ; return *this; }
    PressureValue &operator=(double v) { val = v ; return *this; }
    PressureValue &operator=(float v) { val = v ; return *this; }
    PressureValue &operator=(int v) { val = v ; return *this; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const PressureValue &ss)
  {
    s << ss.val << " pascal" << std::endl ;
    return s ;
  }
    
  inline std::istream &operator>>(std::istream &s, PressureValue &ss) {
    ss.val = get_val_in_units(s,"pascal") ;
    return s ;
  }        
  
}

namespace Loci {

  template<class T> class ValConverter {
    T & ref ;
  public:
    ValConverter(T &iref) : ref(iref) {}
    int getSize() const {
      return 1 ;
    }
    void getState(double *buf, int &size) {
      buf[0] = ref.val ;
      size = 1 ;
    }
    void setState(double *buf, int size) {
      ref.val = buf[0] ;
    }
  } ;
  
  template<> struct data_schema_traits<flowPsi::TimeValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef ValConverter<flowPsi::TimeValue> Converter_Type ;
  } ;
  template<> struct data_schema_traits<flowPsi::TemperatureValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef ValConverter<flowPsi::TemperatureValue> Converter_Type ;
  } ;
  template<> struct data_schema_traits<flowPsi::PressureValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef ValConverter<flowPsi::PressureValue> Converter_Type ;
  } ;
}

namespace flowPsi {
  struct residual {
    real rtrms,etrms,mtrms,mtmax,etmax,rtmax ;
    vect3d ccmaxe_center ;
    vect3d ccmaxr_center ;
    unsigned long bcErrorCode ;
  } ;

}


namespace Loci {
  template<> struct data_schema_traits<flowPsi::residual> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::residual()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::residual,rtrms) ;
      LOCI_INSERT_TYPE(ct,flowPsi::residual,etrms) ;
      LOCI_INSERT_TYPE(ct,flowPsi::residual,mtrms) ;

      LOCI_INSERT_TYPE(ct,flowPsi::residual,etmax) ;
      LOCI_INSERT_TYPE(ct,flowPsi::residual,rtmax) ;

      LOCI_INSERT_TYPE(ct,flowPsi::residual,ccmaxe_center) ;
      LOCI_INSERT_TYPE(ct,flowPsi::residual,ccmaxr_center) ;

      LOCI_INSERT_TYPE(ct,flowPsi::residual,bcErrorCode) ;
      return DatatypeP(ct) ;
    }
  } ;
}

#endif
