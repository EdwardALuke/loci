// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_PARAMETRIC_FACE_DATA_H
#define LOCI_FVMADAPT_PARAMETRIC_FACE_DATA_H

#include <Loci.h>
#include <FVMAdapt/parametric/face_plan.h>
#include <cmath>

namespace Loci {
namespace FVMAdapt {
  /// Canonical source-edge controls agreed by incident cells; empty controls mean straight.
  struct EdgeGeometryData {
    std::vector<Point> controls ;
    /// Unassigned is the reduction identity; a conflict remains failed through every merge.
    bool assigned=false, failed=false ;
    void merge(const EdgeGeometryData &other) {
      if(this==&other) return ;
      if(failed || other.failed || (assigned && other.assigned && controls!=other.controls)) {
        failed=true ; controls.clear() ; return ;
      }
      if(!assigned && other.assigned) { controls=other.controls ; assigned=true ; }
    }
    Point position(const Point &first,const Point &last,double t) const {
      if(failed || !assigned || (!controls.empty() && controls.size()!=2))
        throw std::invalid_argument("invalid or conflicting source edge geometry") ;
      Point result ; const double s=1-t ;
      for(unsigned d=0;d<3;++d) result[d]=first[d]+(controls.empty() ? t*(last[d]-first[d]):
        3*s*s*t*(controls[0][d]-first[d])+3*s*t*t*(controls[1][d]-first[d])+t*t*t*(last[d]-first[d])) ;
      return result ;
    }
  } ;
  class EdgeGeometryDataConverter {
    EdgeGeometryData &value_ ;
  public:
    explicit EdgeGeometryDataConverter(EdgeGeometryData &value):value_(value) {}
    int getSize() const { return !value_.failed && value_.assigned && !value_.controls.empty() ? 7:1 ; }
    void getState(double *buffer,int &size) const {
      // Tags: unassigned, straight, cubic, failed; cubic adds two XYZ controls.
      size=getSize() ; buffer[0]=value_.failed ? 3 : !value_.assigned ? 0 : value_.controls.empty() ? 1:2 ;
      if(size==7) for(unsigned i=0;i<2;++i) for(unsigned d=0;d<3;++d) buffer[1+3*i+d]=value_.controls.at(i)[d] ;
    }
    void setState(const double *buffer,int size) {
      value_=EdgeGeometryData{} ;
      if(size<1 || !std::isfinite(buffer[0]) || buffer[0]<0 || buffer[0]>3 || buffer[0]!=std::floor(buffer[0]) ||
         size!=(buffer[0]==2 ? 7:1)) { value_.failed=true ; return ; }
      value_.assigned=buffer[0]==1 || buffer[0]==2 ; value_.failed=buffer[0]==3 ;
      if(buffer[0]==2) {
        value_.controls.resize(2) ;
        for(unsigned i=0;i<2;++i) for(unsigned d=0;d<3;++d) {
          const double coordinate=buffer[1+3*i+d] ;
          if(!std::isfinite(coordinate)) value_.failed=true ;
          value_.controls[i][d]=coordinate ;
        }
        if(value_.failed) value_.controls.clear() ;
      }
    }
  } ;
  inline std::ostream &operator<<(std::ostream &out,const EdgeGeometryData &value) {
    return out<<(value.failed ? "conflicting source edge geometry":value.controls.empty() ? "straight edge":"cubic edge") ;
  }
  inline std::istream &operator>>(std::istream &in,EdgeGeometryData &) { in.setstate(std::ios::failbit) ; return in ; }

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
  template<> struct data_schema_traits<FVMAdapt::EdgeGeometryData> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef FVMAdapt::EdgeGeometryDataConverter Converter_Type ;
  } ;
  template<> struct data_schema_traits<FVMAdapt::FacePlanData> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef unsigned char Converter_Base_Type ;
    typedef FVMAdapt::FacePlanDataConverter Converter_Type ;
  } ;
}
#endif
