//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
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
#ifndef FITFUNC_H
#define FITFUNC_H
#include <vector>
#include <map> 
#include <string>
#include <algorithm>
#include <iostream>
#include <Tools/basic_types.h>

namespace fitfunction {
  using std::vector ;
  using std::map ;
  using std::string ;
  using std::ostream ;
  using std::istream ;
  using std::endl ;
  using std::cerr ;
  using std::min ;
  using std::max ;
  using Loci::realToDouble ;

  struct spline_segments {
    double y0,y1,a,b ;
    template<class T> T ft(T t) const { 
      T omt = 1.-t ; 
      return y0*omt+t*y1+t*omt*(a*omt+b*t) ;
    }
    template<class T> T ft(T t, T &fp) const {
      const T omt = 1.-t ; 
      const T ab=(a*omt+b*t) ;
      fp = y1-y0+omt*ab-t*(ab+omt*(a-b)) ;
      return y0*omt+t*y1+t*omt*(ab) ;
    }
    template<class T> T ftp(T t) const {
      T omt = 1.-t ; 
      return y1-y0+omt*(a*omt+b*t)-t*(a*omt+b*t+omt*(a-b)) ;
    }
  } ;

  // FitFunction is a class that stores and evaluates the fitted function.
  class FitFunction {
    double start, delta,rdelta  ;
    vector<spline_segments> fit ;
  public:
    FitFunction() { start = 0; delta = 0; rdelta = 1.0; }
    FitFunction(double constval) {  // Constant function
      start=-1e-30 ;
      delta = 2e30 ;
      rdelta = 1./delta ;
      fit = vector<spline_segments>(1) ;
      fit[0].y0 = constval ;
      fit[0].y1 = constval ;
      fit[0].a = 0 ;
      fit[0].b = 0 ;
    }
    // Line Segment
    FitFunction(double x0, double y0, double x1, double y1) {
      start =  x0 ;
      delta =  x1 - x0 ;
      rdelta = 1./delta ;
      fit = vector<spline_segments>(1) ;
      fit[0].y0 = y0 ;
      fit[0].y1 = y1 ;
      fit[0].a = 0 ;
      fit[0].b = 0 ;
    }

    FitFunction(double strt, double end, const vector<spline_segments> &spline) {
      fit = spline ;
      start = strt ;
      delta = (end-start)/double(fit.size()) ;
      rdelta = 1./delta ;
    }

    // Get minimum input for eval
    double minInputValue() const { return start ; }
    // Get maximum output for eval
    double maxInputValue() const { return start + delta*double(fit.size()); }
    // Evaluate function
    template<class T> T eval(T x) const {
      T c = (x-start)*rdelta ;
      int i = max(0,min(int(realToDouble(c)),int(fit.size())-1)) ;
      T t = max<T>(0.0,min<T>(c-T(i),1.0)) ;
      return fit[i].ft(t) ;
    }
    // Evaluate function and its derivative
    template<class T> T evalp(T x, T &fp) const {
      const T c = (x-start)*rdelta ;
      const int i = max(0,min(int(realToDouble(c)),int(fit.size())-1)) ;
      const T t = max<T>(0.0,min<T>(c-T(i),1.0)) ;
      const T f = fit[i].ft(t,fp) ;
      //convert derivative to global coordinates (or zero if at endpoints)
      fp = (c<0.0 || c > T(fit.size()))?0.0:fp*rdelta ;

      return f ;
    }
    // Evaluate inverse function
    template<class T> T evalInverse(T f, T tol = 1e-5, T xguess = 0) const {
      T xmin = minInputValue() ;
      T xmax = maxInputValue() ;
      T x = xguess ;
      if(xguess < xmin || xguess > xmax)
	x = 0.5*(xmin+xmax) ;

      int MAXITER = 50 ;
      int iter = 0;
      // Use bracketed newton method to find solution
      for(iter=0;iter<MAXITER;++iter) {
	T dfdx ;
	T fdiff = evalp(x,dfdx)-f ;
	if(fabs(fdiff)< tol*fabs(f))
	  return x ;
	if(xmax-xmin < tol*xmax)
	  return x ;
	if(fdiff > 0.0)
	  xmax = x ;
	if(fdiff < 0.0)
	  xmin = x ;
	x -= fdiff/dfdx ;
	if(x > xmax || x < xmin)
	  x = 0.5*(xmax+xmin) ;
      }
      cerr << "evalInverse failed to converge" << endl ;
      return x ;
    }

    int SerializeSize() {
      return fit.size()*4+3 ;
    }
    void injectBuffer(double *db, int size) {
      if(size != SerializeSize()) {
	std::cerr << "size error in injectBuffer" << std::endl ;
      }
      db[0] = start ;
      db[1] = delta ;
      db[2] = double(fit.size()) ;
      for(size_t i=0;i<fit.size();++i) {
	int off = 3+i*4 ;
	db[off+0] = fit[i].y0 ;
	db[off+1] = fit[i].y1 ;
	db[off+2] = fit[i].a ;
	db[off+3] = fit[i].b ;
      }
    }
    void extractBuffer(double *db,int size) {
      start = db[0] ;
      delta = db[1] ;
      rdelta = 1./delta ;
      int nfit = int(db[2]);
      fit = vector<spline_segments>(nfit) ;
      for(int i=0;i<nfit;++i) {
	int off = 3+i*4 ;
	fit[i].y0 = db[off+0] ;
	fit[i].y1 = db[off+1] ;
	fit[i].a = db[off+2] ;
	fit[i].b = db[off+3] ;
      }
    }    
    size_t size() const { return fit.size() ; }

    std::ostream &Print(std::ostream &s) const {
      s.precision(14) ;
      s << start << ' ' << delta << ' ' << fit.size() << std::endl ;
      int n = fit.size() ;
      s << fit[0].y0 << std::endl ;

      for(int i=0;i<n;++i) 
	s << fit[i].y1 << std::endl ;

      for(int i=0;i<n;++i) 
	s << fit[i].a << ' ' << fit[i].b << std::endl ;

      return s ;
    }
    std::istream &Input(std::istream &s) {
      double val = 0.0;
      //      s >> start ;
      s >> val ;
      start = val ;
      //      s >> delta ;
      s >> val ;
      delta = val ;
      rdelta = 1./delta ;
      int n = 0 ;
      s >> n ;
      fit = vector<spline_segments>(n) ;
      //      s >> fit[0].y0 ;
      s >> val ;
      fit[0].y0 = val ;
      for(int i=0;i<n;++i) {
	//	s >> fit[i].y1 ;
	s >> val ;
	fit[i].y1 = val ;
      }

      for(int i=0;i<n-1;++i) 
	fit[i+1].y0 = fit[i].y1 ;

      for(int i=0;i<n;++i) {
	//	s >> fit[i].a >> fit[i].b  ;
	s >> val ;
	fit[i].a = val ;
	s >> val ;
	fit[i].b = val ;
      }
      return s ;
    }
  } ;
  
  inline std::ostream &operator<<(std::ostream &s, const FitFunction &fit) {
    return fit.Print(s) ;
  }
  inline std::istream &operator>>(std::istream &s, FitFunction &fit) {
    return fit.Input(s) ;
  }


  // This class creates a database of fit functions that can be read or
  // written to an ascii file
  class FitFunctionDB {
    map <string,FitFunction> fitsDB ;
  public:
    void clearDB() { fitsDB.clear() ; }

    void insertDB(const string &name, const FitFunction &fit) {
      fitsDB[name] = fit ;
    }
    bool nameExists(const string &name) const {
      map<string,FitFunction>::const_iterator mi = fitsDB.find(name) ;
      return (mi != fitsDB.end()) ;
    }
    FitFunction getEntry(const string &name) const {
      map<string,FitFunction>::const_iterator mi = fitsDB.find(name) ;
      if(mi != fitsDB.end()) {
	return mi->second ;
      }
      return FitFunction() ;
    }
    std::ostream &Print(std::ostream &s) const {
      s << fitsDB.size() << endl ;
      map<string,FitFunction>::const_iterator mi ;
      for(mi=fitsDB.begin();mi!=fitsDB.end();++mi) {
	s << mi->first << endl ;
	s << mi->second ;
      }
      return s ;
    }
    std::istream &Input(std::istream &s) {
      int nents = 0 ;
      s >> nents ;
      for(int i=0;i<nents;++i) {
	string name ;
	FitFunction Fit ;
	s >> name ;
	s >> Fit ;
	insertDB(name,Fit) ;
      }
      return s ;
    }
  } ;

  inline std::ostream &operator<<(std::ostream &s, const FitFunctionDB &fit) {
    return fit.Print(s) ;
  }
  inline std::istream &operator>>(std::istream &s, FitFunctionDB &fit) {
    return fit.Input(s) ;
  }

  // Find a cubic fit for a function using n segments.  Establish C1 continuity
  // with limiting to minimize oscillations in the fit.
  template <class FUNC> 
  vector<spline_segments> fit_func(FUNC &f,double start, double end, int n) {
    vector<spline_segments> set1(n) ;
    double delta = (end - start)/double(n) ;
    // First fit the function within each interval
    for(int i=0;i<n;++i) {
      double istart = start+delta*double(i) ;
      double y0 = f(istart) ;
      double y1 = f(istart+delta) ;
      double ya = f(istart+delta/3.0) ;
      double yb = f(istart+2.*delta/3.0) ;
      double ca = 9.*(3.*ya-2.*y0-y1)/2. ;
      double cb = 9.*(3.*yb-y0-2.*y1)/2. ;

      set1[i].y0 = y0 ;
      set1[i].y1 = y1 ;
      set1[i].a = (2.*ca-cb)/3.0 ;
      set1[i].b = (2.*cb-ca)/3.0 ;
    }

    // Now we fix continuity with limiting
    for(int i=0;i<n-1;++i) {
      double fpi = set1[i].ftp(1.0) ;
      double fpn = set1[i+1].ftp(0.0) ;
      if(fabs(fpi) < fabs(fpn))
	set1[i+1].a = fpi-(set1[i+1].y1-set1[i+1].y0) ;
      else
	set1[i].b = (set1[i].y1-set1[i].y0)-fpn ;
    }  
    return set1 ;
  }

  // Fit a function to a specified relative and absolute tolerance.  Relative
  // tolerance is a RMS value.  This routine iteratively adjusts the number
  // of intervals to meet target error.  By default, limits number of segments
  // to 1024.
  template <class FUNC> 
  FitFunction fit_function(FUNC &f,double start, double end, double reltol, 
			   int minn = 3, int maxn=1024, double abstol=1e-20) {
    int n = minn ;
    double err = 0 ;
    double err8 = 0 ;
    vector<spline_segments> spline ;
    double delta = (end-start)/double(n) ;
    do {
      spline=fit_func(f,start,end,n) ;
      delta = (end-start)/double(n) ;
    
      err = 0 ;
      err8 = 0 ;
      for(int j=0;j<n;++j) {
	for(int k=1;k<4;++k) {
	  double t = double(k)*0.25 ;
	  double fapprox = spline[j].ft(t) ;
	  double fexact = f(start+(double(j)+t)*delta) ;
	  double delta = (fexact-fapprox)/max(fabs(fexact),abstol) ;
	  err8 = max(err8,fabs(delta)) ;
	  err += delta*delta/3.0 ;
	}
      }
      err = sqrt(err/n) ;
      if((err8 > 0.05 || err > reltol) && n < maxn) {
	double factor=pow(err/reltol,1./3.)*.5 ;
	if(factor < 1.0)
	  factor *= 2.0 ;
	n = min(max(int(double(n)*factor),n+1),maxn) ;
      } else {
	return FitFunction(start,end,spline) ;
      }
    } while(true) ;
    return FitFunction(start,end,spline) ;
  }

  
}

#endif
