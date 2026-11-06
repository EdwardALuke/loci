// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_QUALITY_H
#define LOCI_FVMADAPT_QUALITY_H
#include <Loci.h>
#include "parametric/mesh_quality.h"

namespace Loci {
namespace FVMAdapt {
  /// One emitted polygon on an original face. Indices belong to the incident roots.
  struct QualityFace { FaceMeasure geometry ; int left=-1, right=-1 ; } ;
  inline std::ostream &operator<<(std::ostream &s,const FaceMeasure &v) {
    for(auto p:{v.center,v.area,v.wireCenter}) for(auto x:p) s<<x<<' ' ;
    for(auto x:v.moment) s<<x<<' ' ;
    return s<<v.warpage ;
  }
  inline std::istream &operator>>(std::istream &s,FaceMeasure &v) {
    for(auto p:{&v.center,&v.area,&v.wireCenter}) for(auto &x:*p) s>>x ;
    for(auto &x:v.moment) s>>x ;
    return s>>v.warpage ;
  }
  inline std::ostream &operator<<(std::ostream &s,const CellMeasure &v) { return s<<v.center[0]<<' '<<v.center[1]<<' '<<v.center[2]<<' '<<v.volume ; }
  inline std::istream &operator>>(std::istream &s,CellMeasure &v) { return s>>v.center[0]>>v.center[1]>>v.center[2]>>v.volume ; }
  inline std::ostream &operator<<(std::ostream &s,const QualityFace &v) { return s<<v.geometry<<' '<<v.left<<' '<<v.right ; }
  inline std::istream &operator>>(std::istream &s,QualityFace &v) { return s>>v.geometry>>v.left>>v.right ; }
  inline std::ostream &operator<<(std::ostream &s,const QualityValue &v) { return s<<v.value<<' '<<v.root<<' '<<v.leaf<<' '<<v.face ; }
  inline std::istream &operator>>(std::istream &s,QualityValue &v) { return s>>v.value>>v.root>>v.leaf>>v.face ; }
  inline std::ostream &operator<<(std::ostream &s,const MeshQualityReport &v) {
    return s<<v.nonOrthogonality<<' '<<v.faceWarpage<<' '<<v.childVolumeRatio<<' '<<v.invalidRoot<<' '<<v.invalidLeaf<<' '<<v.invalidFace ;
  }
  inline std::istream &operator>>(std::istream &s,MeshQualityReport &v) {
    return s>>v.nonOrthogonality>>v.faceWarpage>>v.childVolumeRatio>>v.invalidRoot>>v.invalidLeaf>>v.invalidFace ;
  }
}
  template<> struct data_schema_traits<FVMAdapt::FaceMeasure> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      using T=FVMAdapt::FaceMeasure ; T value ; auto ct=CompoundFactory(value) ;
      const auto add=[&](const std::string &name,const double &field) {
        ct->insert(name,std::size_t(reinterpret_cast<const char *>(&field)-reinterpret_cast<const char *>(&value)),getLociType(field)) ;
      } ;
      for(unsigned d=0;d<3;++d) { add("center"+std::to_string(d),value.center[d]) ; add("area"+std::to_string(d),value.area[d]) ; add("wireCenter"+std::to_string(d),value.wireCenter[d]) ; }
      for(unsigned i=0;i<9;++i) add("moment"+std::to_string(i),value.moment[i]) ;
      add("warpage",value.warpage) ; return DatatypeP(ct) ;
    }
  } ;
  template<> struct data_schema_traits<FVMAdapt::CellMeasure> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      using T=FVMAdapt::CellMeasure ; T value ; auto ct=CompoundFactory(value) ;
      for(unsigned d=0;d<3;++d) ct->insert("center"+std::to_string(d),offsetof(T,center)+d*sizeof(double),getLociType(double())) ;
      LOCI_INSERT_TYPE(ct,T,volume) ; return DatatypeP(ct) ;
    }
  } ;
  template<> struct data_schema_traits<FVMAdapt::QualityFace> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      using T=FVMAdapt::QualityFace ; auto ct=CompoundFactory(T()) ;
      LOCI_INSERT_TYPE(ct,T,geometry) ; LOCI_INSERT_TYPE(ct,T,left) ; LOCI_INSERT_TYPE(ct,T,right) ; return DatatypeP(ct) ;
    }
  } ;
  template<> struct data_schema_traits<FVMAdapt::QualityValue> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      using T=FVMAdapt::QualityValue ; auto ct=CompoundFactory(T()) ;
      LOCI_INSERT_TYPE(ct,T,value) ; LOCI_INSERT_TYPE(ct,T,root) ; LOCI_INSERT_TYPE(ct,T,leaf) ; LOCI_INSERT_TYPE(ct,T,face) ; return DatatypeP(ct) ;
    }
  } ;
  template<> struct data_schema_traits<FVMAdapt::MeshQualityReport> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      using T=FVMAdapt::MeshQualityReport ; auto ct=CompoundFactory(T()) ;
      LOCI_INSERT_TYPE(ct,T,nonOrthogonality) ; LOCI_INSERT_TYPE(ct,T,faceWarpage) ; LOCI_INSERT_TYPE(ct,T,childVolumeRatio) ;
      LOCI_INSERT_TYPE(ct,T,invalidRoot) ; LOCI_INSERT_TYPE(ct,T,invalidLeaf) ; LOCI_INSERT_TYPE(ct,T,invalidFace) ; return DatatypeP(ct) ;
    }
  } ;
}
#endif
