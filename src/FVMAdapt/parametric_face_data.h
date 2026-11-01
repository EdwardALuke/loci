// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_PARAMETRIC_FACE_DATA_H
#define LOCI_FVMADAPT_PARAMETRIC_FACE_DATA_H

#include <Loci.h>
#include <FVMAdapt/parametric/face_plan.h>

namespace Loci {
namespace FVMAdapt {
  /// A face reduction keeps accepted divisions or a failure for collective rollback.
  struct FacePlanData {
    FacePlan plan ;
    bool failed = false ;
    void merge(const FacePlanData &other) {
      if(this == &other) return ;
      if(failed || other.failed) { failed = true ; plan.clear() ; return ; }
      try { plan = mergeFacePlans(plan,other.plan) ; }
      catch(const std::exception &) { failed = true ; plan.clear() ; }
    }
  } ;

  class FacePlanDataConverter {
    FacePlanData &value_ ;
    static constexpr int maxBytes = 64*1024*1024 ;
  public:
    explicit FacePlanDataConverter(FacePlanData &value) : value_(value) {}
    int getSize() const {
      if(value_.plan.size() > maxBytes) throw std::length_error("face plan record exceeds module limit") ;
      return 1+int(value_.plan.size()) ;
    }
    void getState(unsigned char *buffer, int &size) const {
      size = getSize() ; buffer[0] = value_.failed ;
      std::copy(value_.plan.begin(),value_.plan.end(),buffer+1) ;
    }
    void setState(const unsigned char *buffer, int size) {
      value_ = FacePlanData{} ;
      if(size < 1 || size > 1+maxBytes || buffer[0] > 1) { value_.failed = true ; return ; }
      value_.failed = buffer[0] ;
      if(!value_.failed) value_.plan.assign(buffer+1,buffer+size) ;
    }
  } ;
  inline std::ostream &operator<<(std::ostream &out, const FacePlanData &value) {
    return value.failed ? out << "invalid or conflicting face plan" : out << "face plan (" << value.plan.size() << " bytes)" ;
  }
  inline std::istream &operator>>(std::istream &in, FacePlanData &) { in.setstate(std::ios::failbit) ; return in ; }
}
  template<> struct data_schema_traits<FVMAdapt::FacePlanData> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef unsigned char Converter_Base_Type ;
    typedef FVMAdapt::FacePlanDataConverter Converter_Type ;
  } ;
}
#endif
