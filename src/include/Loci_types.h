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
#ifndef LOCI_TYPES_H
#define LOCI_TYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <data_traits.h>

#include <Tools/basic_types.h>
#include <Tools/vtypes.h>
#include <Tools/options_list.h>

namespace Loci {

  // Make special type manager for FADdType
  class FADdType : public AbstractDatatype {
  public:
    FADdType() {}
    hid_t get_hdf5_type() const {
      FADd tmp ;
      hid_t vDatatype = H5Tcreate( H5T_COMPOUND, sizeof(tmp) ) ;
#ifdef NO_OFFSETOF
      size_t offset1 = reinterpret_cast<char *>(&(tmp.value)) - reinterpret_cast<char *>(&tmp) ;
#else
      size_t offset1 = offsetof(FADd,value) ;
#endif
      H5Tinsert(vDatatype,"value",offset1,H5T_NATIVE_DOUBLE) ;

#ifdef NO_OFFSETOF
      size_t offset2 = reinterpret_cast<char *>(&(tmp.grad)) - reinterpret_cast<char *>(&tmp) ;
#else
      size_t offset2 = offsetof(FADd,grad) ;
#endif
      H5Tinsert(vDatatype,"grad",offset2,H5T_NATIVE_DOUBLE) ;
      return vDatatype ;
    }
    std::ostream &output(std::ostream &s, const void *p) const
      { s << *(reinterpret_cast<const FADd *>(p)) ; return s ;}
    std::istream &input(std::istream &s, void *p) const
      { s >> *(reinterpret_cast<FADd *>(p)) ; return s ; }
    int bytesize() const
    { return sizeof(FADd) ; }
  } ;

  template<> struct data_schema_traits<FADd> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      return new FADdType() ;
    }
  };

  // Make special type manager for FAD2dType
  class FAD2dType : public AbstractDatatype {
  public:
    FAD2dType() {}
    hid_t get_hdf5_type() const {
      FAD2d tmp ;
      hid_t vDatatype = H5Tcreate( H5T_COMPOUND, sizeof(tmp) ) ;
#ifdef NO_OFFSETOF
      size_t offset1 = reinterpret_cast<char *>(&(tmp.value)) - reinterpret_cast<char *>(&tmp) ;
#else
      size_t offset1 = offsetof(FAD2d,value) ;
#endif
      H5Tinsert(vDatatype,"value",offset1,H5T_NATIVE_DOUBLE) ;

#ifdef NO_OFFSETOF
      size_t offset2 = reinterpret_cast<char *>(&(tmp.grad)) - reinterpret_cast<char *>(&tmp) ;
#else
      size_t offset2 = offsetof(FAD2d,grad) ;
#endif
      H5Tinsert(vDatatype,"grad",offset2,H5T_NATIVE_DOUBLE) ;

#ifdef NO_OFFSETOF
      size_t offset3 = reinterpret_cast<char *>(&(tmp.grad2)) - reinterpret_cast<char *>(&tmp) ;
#else
      size_t offset3 = offsetof(FAD2d,grad2) ;
#endif
      H5Tinsert(vDatatype,"grad2",offset3,H5T_NATIVE_DOUBLE) ;
      return vDatatype ;
    }
    std::ostream &output(std::ostream &s, const void *p) const
      { s << *(reinterpret_cast<const FAD2d *>(p)) ; return s ;}
    std::istream &input(std::istream &s, void *p) const
      { s >> *(reinterpret_cast<FAD2d *>(p)) ; return s ; }
    int bytesize() const
    { return sizeof(FAD2d) ; }
  } ;

  template<> struct data_schema_traits<FAD2d> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      return new FAD2dType() ;
    }
  };

  // Make special type manager for MFADdType
  class MFADdType : public AbstractDatatype {
  public:
    MFADdType() {}
    hid_t get_hdf5_type() const {
      MFADd tmp ;
      hid_t vDatatype = H5Tcreate( H5T_COMPOUND, sizeof(tmp) ) ;
#ifdef NO_OFFSETOF
      size_t offset1 = reinterpret_cast<char *>(&(tmp.value)) - reinterpret_cast<char *>(&tmp) ;
#else
      size_t offset1 = offsetof(MFADd,value) ;
#endif
      H5Tinsert(vDatatype,"value",offset1,H5T_NATIVE_DOUBLE) ;

      // note, not an array, this needs to be fixed
#ifdef NO_OFFSETOF
      size_t offset2 = reinterpret_cast<char *>(&(tmp.grad[0])) - reinterpret_cast<char *>(&tmp) ;
#else
      size_t offset2 = offsetof(MFADd,grad) ;
#endif
      H5Tinsert(vDatatype,"grad",offset2,H5T_NATIVE_DOUBLE) ;

      //#ifdef NO_OFFSETOF
      //      size_t offset3 = reinterpret_cast<char *>(&(tmp.maxN)) - reinterpret_cast<char *>(&tmp) ;
      //#else
      //      size_t offset3 = offsetof(MFADd,maxN) ;
      //#endif
      //      H5Tinsert(vDatatype,"maxN",offset3,H5T_NATIVE_INT) ;
      return vDatatype ;
    }
    std::ostream &output(std::ostream &s, const void *p) const
      { s << *(reinterpret_cast<const MFADd *>(p)) ; return s ;}
    std::istream &input(std::istream &s, void *p) const
      { s >> *(reinterpret_cast<MFADd *>(p)) ; return s ; }
    int bytesize() const
    { return sizeof(MFADd) ; }
  } ;

  template<> struct data_schema_traits<MFADd> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      return new MFADdType() ;
    }
  };

  // Make special type manager for VFADType
  class VFADType : public AbstractDatatype {
  public:
    VFADType() {}
    hid_t get_hdf5_type() const {
      VFAD tmp ;
      hid_t vDatatype = H5Tcreate( H5T_COMPOUND, sizeof(tmp) ) ;
#ifdef NO_OFFSETOF
      size_t offset1 = reinterpret_cast<char *>(&(tmp.data.value)) - reinterpret_cast<char *>(&tmp) ;
#else
      size_t offset1 = offsetof(VFAD,data.value) ;
#endif
      H5Tinsert(vDatatype,"value",offset1,H5T_NATIVE_DOUBLE) ;

      // note, not an array, this needs to be fixed
#ifdef NO_OFFSETOF
      size_t offset2 = reinterpret_cast<char *>(&(tmp.data.grad[0])) - reinterpret_cast<char *>(&tmp) ;
#else
      size_t offset2 = offsetof(VFAD,data.grad) ;
#endif
      int dim = VFAD::maxN ;
      DatatypeP atype = new ArrayType(getLociType(float()),sizeof(tmp.data.grad),1,&dim) ;
      H5Tinsert(vDatatype,"grad",offset2,atype->get_hdf5_type()) ;

      return vDatatype ;
    }
    std::ostream &output(std::ostream &s, const void *p) const
      { s << *(reinterpret_cast<const VFAD *>(p)) ; return s ;}
    std::istream &input(std::istream &s, void *p) const
      { s >> *(reinterpret_cast<VFAD *>(p)) ; return s ; }
    int bytesize() const
    { return sizeof(VFAD) ; }
  } ;

  template<> struct data_schema_traits<VFAD> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      return new VFADType() ;
    }
  };


  template <class T>
    struct data_schema_traits< vector3d<T> > {
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        CompoundDatatypeP ct = CompoundFactory(vector3d<T>() ) ;

        LOCI_INSERT_TYPE(ct,vector3d<T>,x) ;
        LOCI_INSERT_TYPE(ct,vector3d<T>,y) ;
        LOCI_INSERT_TYPE(ct,vector3d<T>,z) ;
        return DatatypeP(ct) ;
      }
    };


  template <class T>
    struct data_schema_traits< tensor3d<T> > {
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        CompoundDatatypeP ct = CompoundFactory(tensor3d<T>()) ;

        LOCI_INSERT_TYPE(ct,tensor3d<T>,x) ;
        LOCI_INSERT_TYPE(ct,tensor3d<T>,y) ;
        LOCI_INSERT_TYPE(ct,tensor3d<T>,z) ;
        return DatatypeP(ct) ;
      }
    };



  template <class T>
    struct  data_schema_traits< vector2d<T> > {
    public:
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        CompoundDatatypeP ct = CompoundFactory(vector2d<T>()) ;
        LOCI_INSERT_TYPE(ct,vector2d<T>,x) ;
        LOCI_INSERT_TYPE(ct,vector2d<T>,y) ;
        return DatatypeP(ct) ;
      }
    };




  template <class T,size_t n> inline std::ostream &
    operator<<(std::ostream &s, const Array<T,n> &v) {
    for(size_t i=0;i<n;++i)
      s << v[i] << ' ' ;
    return s ;
  }

  template <class T,size_t n> inline std::istream &
    operator>>(std::istream &s, Array<T,n> &v) {
    for(size_t i=0;i<n;++i)
      s >> v[i] ;
    return s ;
  }


  template <class T,size_t n>
  class data_schema_traits< Array<T,n> > {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      int dim = n ;
      return new ArrayType(getLociType(T()),sizeof(Array<T,n>),1,&dim) ;
    }
  };

  template <class T,size_t n> inline std::ostream &
    operator<<(std::ostream &s, const vtype<T,n> &v) {
    for(size_t i=0;i<n;++i)
      s << v[i] << ' ' ;
    return s ;
  }

  template <class T,size_t n> inline std::istream &
    operator>>(std::istream &s, vtype<T,n> &v) {
    for(size_t i=0;i<n;++i)
      s >> v[i] ;
    return s ;
  }


  template <class T,size_t n>
  class data_schema_traits< vtype<T,n> > {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      int dim = n ;
      return new ArrayType(getLociType(T()),sizeof(vtype<T,n>),1,&dim) ;
    }
  };

  template <>
  class data_schema_traits< vtype<float,8> > {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      int dim = 8 ;
      return new ArrayType(getLociType(float()),
                           sizeof(vtype<float,8>),1,&dim) ;
    }
  };

  template <>
  class data_schema_traits< vtype<double,4> > {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      int dim = 4 ;
      return new ArrayType(getLociType(double()),
                           sizeof(vtype<double,4>),1,&dim) ;
    }
  };

  template <>
  class data_schema_traits< vtype<int32_t,8> > {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      int dim = 8 ;
      return new ArrayType(getLociType(int32_t()),
                           sizeof(vtype<int32_t,8>),1,&dim) ;
    }
  };

  template <> struct data_schema_traits<options_list> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<options_list> Converter_Type ;
  } ;

  template<class T> class UnitValueConverter {
    T & ref ;
  public:
    UnitValueConverter(T &iref) : ref(iref) {}
    int getSize() const {
      return ref.getSize() ;
    }
    void getState(double *buf, int &size) const {
      ref.getState(buf,size) ;
    }
    void setState(double *buf, int size) {
      ref.setState(buf,size) ;
    }
  } ;

  template<> struct data_schema_traits<NonDimensionalValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<NonDimensionalValue> Converter_Type ;
  } ;
  template<> struct data_schema_traits<TimeValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<TimeValue> Converter_Type ;
  } ;
  template<> struct data_schema_traits<LengthValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<LengthValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<MassValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<MassValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<AreaValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<AreaValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<DensityValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<DensityValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<EnergyValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<EnergyValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<PowerValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<PowerValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<FlowValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<FlowValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<SpeedValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<SpeedValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<DynamicViscosityValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<DynamicViscosityValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<KinematicViscosityValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<KinematicViscosityValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<ThermalConductivityValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<ThermalConductivityValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<VolumeValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<VolumeValue> Converter_Type ;
  } ;

  template<> struct data_schema_traits<TemperatureValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<TemperatureValue> Converter_Type ;
  } ;
  template<> struct data_schema_traits<PressureValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef UnitValueConverter<PressureValue> Converter_Type ;
  } ;

  /// Class to allow for simple comma-separated lists in vars file values.
  class list_input {
  public:
    // Raw list expression from the vars option (whitespace removed).
    std::string namelist;

    // Output operator for list_input, which prints the raw list expression.
    friend std::ostream& operator<<(std::ostream& s, const list_input& n) {
      s << n.namelist << std::endl;
      return s;
    }

    // Input operator for list_input, which reads in characters until a newline, removing whitespace.
    friend std::istream& operator>>(std::istream& s, list_input& n) {
      n.namelist = std::string();
      while (s.peek() != '\n' && s.peek() != '\r' && s.peek() != EOF) {
        char c = s.get();
        if (c != ' ' && c != '\t') n.namelist += c;
      }
      return s;
    }
};

  /// Class to serialize list_input.
  class list_input_schema_converter {
    list_input &ref ;
  public:
    explicit list_input_schema_converter(list_input &iref): ref(iref) {}
    int getSize() const {
      return ref.namelist.size() ;
    }
    void getState(char *buf, int &size) {
      size = getSize() ;
      for(int i=0;i<size;++i) {
        buf[i] = ref.namelist[i] ;
      }
    }
    void setState(char *buf, int size) {
      ref.namelist = "" ;
      for(int i=0;i<size;++i) {
        ref.namelist += buf[i] ;
      }
    }
  } ;

  /// Schema converter for list_input.
  template<> struct data_schema_traits<list_input> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef list_input_schema_converter Converter_Type ;
  } ;


}

#endif
