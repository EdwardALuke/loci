// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#include <FVMAdapt/spacing_source.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace Loci {
namespace FVMAdapt {
  namespace {
    double distance(const Point &a,const Point &b) { return std::hypot(a[0]-b[0],a[1]-b[1],a[2]-b[2]) ; }
    Point average(const std::vector<Point> &points) {
      if(points.empty()) throw std::invalid_argument("spacing requests require cell vertices") ;
      Point center{} ;
      for(const auto &p:points) for(unsigned d=0;d<3;++d) center[d]+=p[d]/points.size() ;
      return center ;
    }
  }

  std::vector<SpacingSource> readSpacingSources(const std::string &filename) {
    std::ifstream input(filename) ;
    if(!input) throw std::invalid_argument("unable to read spacing source file: "+filename) ;
    int count=0 ;
    if(!(input>>count) || count<=0) throw std::invalid_argument("spacing source count must be positive") ;
    std::vector<SpacingSource> sources ;
    for(int i=0;i<count;++i) {
      SpacingSource s{} ;
      for(auto &x:s.start) input>>x ;
      for(auto &x:s.end) input>>x ;
      input>>s.innerRadius>>s.innerSpacing>>s.outerRadius>>s.outerSpacing>>s.exponent ;
      if(!input) throw std::invalid_argument("incomplete spacing source "+std::to_string(i+1)) ;
      bool finite=true ;
      for(double x:s.start) finite=finite && std::isfinite(x) ;
      for(double x:s.end) finite=finite && std::isfinite(x) ;
      for(double x:{s.innerRadius,s.innerSpacing,s.outerRadius,s.outerSpacing,s.exponent}) finite=finite && std::isfinite(x) ;
      if(!finite || s.innerRadius<0 || s.outerRadius<=0 || s.innerRadius>s.outerRadius ||
         s.innerSpacing<=0 || s.outerSpacing<=0 || s.exponent<0)
        throw std::invalid_argument("invalid radii, spacings or exponent in spacing source "+std::to_string(i+1)) ;
      sources.push_back(s) ;
    }
    std::string extra ;
    if(input>>extra) throw std::invalid_argument("unexpected data after spacing sources") ;
    return sources ;
  }

  double sourceSpacing(const Point &p,const SpacingSource &s) {
    const double length=distance(s.start,s.end) ;
    double along=0 ; Point closest=s.start ;
    if(length>0) {
      Point direction{} ;
      for(unsigned d=0;d<3;++d) { direction[d]=(s.end[d]-s.start[d])/length ; along+=(p[d]-s.start[d])*direction[d] ; }
      for(unsigned d=0;d<3;++d) closest[d]+=std::clamp(along,0.,length)*direction[d] ;
    }
    const double radius=distance(p,closest) ;
    if(along>0 && along<length && radius<=s.innerRadius) return s.innerSpacing ;
    if(radius<=s.outerRadius || s.exponent==0) return s.outerSpacing ;
    // A distant source may request effectively unlimited spacing without overflowing.
    const double growth=std::log(radius)-std::log(s.outerRadius),base=std::log(s.outerSpacing) ;
    if(growth<=0) return s.outerSpacing ;
    const double maximum=std::numeric_limits<double>::max() ;
    if(s.exponent>=(std::log(maximum)-base)/growth) return maximum ;
    return std::exp(base+s.exponent*growth) ;
  }

  double minimumSourceSpacing(const std::vector<Point> &points,const std::vector<SpacingSource> &sources) {
    const auto center=average(points) ; double spacing=std::numeric_limits<double>::max() ;
    for(const auto &source:sources) {
      for(const auto &p:points) spacing=std::min(spacing,sourceSpacing(p,source)) ;
      spacing=std::min(spacing,sourceSpacing(center,source)) ;
    }
    return spacing ;
  }

  bool refineForSpacing(const std::vector<Point> &points,const std::vector<SpacingSource> &sources,double shortest) {
    const auto center=average(points) ; Point lo=points.front(),hi=lo ;
    for(const auto &p:points) for(unsigned d=0;d<3;++d) { lo[d]=std::min(lo[d],p[d]) ; hi[d]=std::max(hi[d],p[d]) ; }
    const auto inside=[&](const Point &p) { for(unsigned d=0;d<3;++d) if(p[d]<lo[d] || p[d]>hi[d]) return false ; return true ; } ;
    for(const auto &source:sources) if(shortest>source.outerRadius/2) {
      Point midpoint{} ; for(unsigned d=0;d<3;++d) midpoint[d]=source.start[d]/2+source.end[d]/2 ;
      if(inside(source.start) || inside(source.end) || distance(center,midpoint)<shortest+source.outerRadius) return true ;
    }
    return shortest/2>=minimumSourceSpacing(points,sources) ;
  }
}
}
