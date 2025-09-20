//#############################################################################
//#
//# Copyright 2017-2025, Mississippi State University
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
#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <iostream> 
#include <Tools/tools.h>
#include <algorithm>
#include <math.h>

namespace Loci {
  
  //---------------------- second derivative type

  class FAD2d {
  public:

    double value;
    double grad;
    double grad2;
    /* MPGCOMMENT [05-12-2017 13:49] ---> constructor */
#define gradcheck() //if(grad!=0) debugger_() ; //fatal(grad!=0.0) ;
    template< class T1>
      FAD2d(T1 a0): value(a0), grad(0.0),grad2(0.0) { gradcheck();}
  FAD2d(double a0, double b0, double c0): value(a0), grad(b0), grad2(c0) { gradcheck();}
  FAD2d(const FAD2d &u) : value(u.value), grad(u.grad), grad2(u.grad2) { gradcheck();}

  FAD2d() : value(0.0), grad(0.0),grad2(0.0) { gradcheck();}
// --------------------------------------------------------------------
// REH ADDED BELOW - 20220202 - Missing code from initial 2017 dev mods
// --------------------------------------------------------------------
 FAD2d(int a0         ,int b0)        : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(int a0         ,float b0)        : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(int a0         ,double b0)        : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(int a0         ,long double b0)        : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(float a0      ,int b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(float a0      ,float b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(float a0      ,double b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(float a0      ,long double b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(double a0      ,int b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(double a0      ,float b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(double a0      ,double b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(double a0      ,long double b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(long double a0      ,int b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(long double a0      ,float b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(long double a0      ,double b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}
 FAD2d(long double a0      ,long double b0)       : value(a0), grad(b0), grad2(0.0) { gradcheck();}

 FAD2d(double a0      ,int b0        , int c0)       : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,int b0        , float c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,int b0        , double c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,int b0        , long double c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,float b0      , int c0)       : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,float b0      , float c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,float b0      , double c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,float b0      , long double c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,double b0     , int c0)       : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,double b0     , float c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 //FAD2d(double a0      ,double b0     , double c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,double b0     , long double c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,long double b0, int c0)       : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,long double b0, float c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,long double b0, double c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(double a0      ,long double b0, long double c0)     : value(a0), grad(b0), grad2(c0) { gradcheck();}
 FAD2d(bool &u)        : value((u)?1.0:0.0), grad(0.0), grad2(0.0) { gradcheck();}
 //FAD2d(const FAD2d &u) : value(u.value), grad(u.grad), grad2(u.grad2) { gradcheck();}
 FAD2d(int a0)        : value(a0), grad(0.0), grad2(0.0) { gradcheck();}
 FAD2d(float a0)        : value(a0), grad(0.0), grad2(0.0) { gradcheck();}
 FAD2d(double a0)        : value(a0), grad(0.0), grad2(0.0) { gradcheck();}
 FAD2d(long double a0)        : value(a0), grad(0.0), grad2(0.0) { gradcheck();}
// --------------------------------------------------------------------
    FAD2d& operator =(const FAD2d &u) {
      value = u.value;
      grad = u.grad;
      grad2 = u.grad2 ;
      return *this;
    }
    FAD2d& operator =(const int &u) {
      value = u;
      grad = 0.0;
      grad2 = 0.0 ;
      return *this;
    }
    FAD2d& operator =(const float &u) {
      value = u;
      grad = 0.0;
      grad2 = 0.0;
      return *this;
    }
    FAD2d& operator =(const double &u) {
      value = u;
      grad = 0.0;
      grad2 = 0.0;
      return *this;
    }
    FAD2d& operator =(const long double &u) {
      value = u;
      grad = 0.0;
      grad2 = 0.0;
      return *this;
    }


    FAD2d operator +(const FAD2d &v) const{
      return FAD2d(value+v.value, grad+v.grad, grad2+v.grad2);
    }
    FAD2d operator +() const{
      return FAD2d(value, grad, grad2);
    }
    FAD2d operator +(const int &v) const{
      return FAD2d(value+(double)v, grad, grad2);
    }
    FAD2d operator +(const float &v) const{
      return FAD2d(value+(double)v, grad, grad2);
    }
    FAD2d operator +(const double &v) const{
      return FAD2d(value+(double)v, grad, grad2);
    }
    FAD2d operator +(const long double &v) const{
      return FAD2d(double(value+v), grad, grad2);
    }

    FAD2d& operator +=(const FAD2d &u) {
      value+=u.value;
      grad+=u.grad;
      grad2+=u.grad2 ;
      return *this;
    }
    FAD2d& operator +=(const int &u) {
      value+=u;
      return *this;
    }
    FAD2d& operator +=(const float &u) {
      value+=u;
      return *this;
    }
    FAD2d& operator +=(const double &u) {
      value+=u;
      return *this;
    }
    FAD2d& operator +=(const long double &u) {
      value+=u;
      return *this;
    }

    FAD2d operator -(const FAD2d &v) const{
      return FAD2d(value-v.value, grad-v.grad,grad2-v.grad2);
    }
    FAD2d operator -() const {
      return FAD2d(-value, -grad, -grad2);
    }
    FAD2d operator -(const int &v) const{
      return FAD2d(value-(double)v, grad, grad2);
    }
    FAD2d operator -(const float &v) const{
      return FAD2d(value-(double)v, grad, grad2);
    }
    FAD2d operator -(const double &v) const{
      return FAD2d(value-(double)v, grad, grad2);
    }
    FAD2d operator -(const long double &v) const{
      return FAD2d(double(value-v), grad, grad2);
    }

    FAD2d& operator -=(const FAD2d &u) {
      value-=u.value;
      grad-=u.grad;
      grad2 -= u.grad2 ;
      return *this;
    }
    FAD2d& operator -=(const int &u) {
      value-=u;
      return *this;
    }
    FAD2d& operator -=(const float &u) {
      value-=u;
      return *this;
    }
    FAD2d& operator -=(const double &u) {
      value-=u;
      return *this;
    }
    FAD2d& operator -=(const long double &u) {
      value-=u;
      return *this;
    }

    FAD2d operator *(const FAD2d &v) const {
      return FAD2d(value*v.value,  grad*v.value + v.grad*value, 
		   grad2*v.value + value*v.grad2 + 2.*grad*v.grad);
    }
    FAD2d operator *(const int &v) const{
      return FAD2d(value*(double)v,  grad*(double)v, grad2*(double) v);
    }
    FAD2d operator *(const float &v) const{
      return FAD2d(value*(double)v,  grad*(double)v, grad2*(double) v);
    }
    FAD2d operator *(const double &v) const{
      return FAD2d(value*(double)v,  grad*(double)v, grad2*(double) v);
    }
    FAD2d operator *(const long double &v) const{
      return FAD2d(double(value*v),  double(grad*v), double(grad2*v));
    }
    //    template<class T> FAD2d& operator *=(const T &u) {
    //      value*=(double)u;
    //      return *this;
    //    }
    FAD2d& operator *=(const FAD2d &u) {
      grad2 = grad2*u.value + value*u.grad2 + 2.*grad*u.grad ;
      grad = grad*u.value + u.grad * value;
      value*=u.value;
      return *this;
    }
    FAD2d& operator *=(const int &u) {
      value*=u;
      grad *=u;
      grad2 *= u ;
      return *this;
    }
    FAD2d& operator *=(const float &u) {
      value*=u;
      grad *=u;
      grad2 *= u;
      return *this;
    }
    FAD2d& operator *=(const double &u) {
      value*=u;
      grad *=u;
      grad2 *= u ;
      return *this;
    }
    FAD2d& operator *=(const long double &u) {
      value*=u;
      grad *=u;
      grad2 *= u;
      return *this;
    }

    FAD2d operator /(const FAD2d &v) const{
      double d = v.value ;
      double d2 = d*d ;
      double d3 = d*d2 ;
      return FAD2d(value/d, 
		   (grad*v.value - value*v.grad)/d2,
		   grad2/d - grad*v.grad/d2 +
		   2.*v.grad*v.grad*value/d3 -
		   (v.grad*grad+v.grad2*value)/d2) ;
    }
    FAD2d operator /(const int &v) const{
      return FAD2d(value/v, (grad/(double)v), (grad2/double(v)));
    }
    FAD2d operator /(const float &v) const{
      return FAD2d(value/v, (grad/(double)v), (grad2/double(v)));
    }
    FAD2d operator /(const double &v) const{
      return FAD2d(value/v, (grad/(double)v), (grad2/double(v)));
    }
    FAD2d operator /(const long double &v) const{
      return FAD2d(double(value/v), double(grad/v), double(grad2/v));
    }

    FAD2d& operator /=(const FAD2d &u) {
      double d = u.value ;
      double d2 = d*d ;
      double d3 = d*d2 ;
      grad2 = (grad2/d - grad*u.grad/d2 +2.*u.grad*u.grad*value/d3 -
	       (u.grad*grad+u.grad2*value)/d2) ;
      grad = (grad*u.value - value * u.grad) / d2 ;
      value/=u.value;
      return *this;
    }
    FAD2d& operator /=(const int &u) {
      value/=u;
      grad /=u;
      grad2 /= u ;
      return *this;
    }
    FAD2d& operator /=(const float &u) {
      value/=u;
      grad /=u;
      grad2 /=u ;
      return *this;
    }
    FAD2d& operator /=(const double &u) {
      value/=u;
      grad /=u;
      grad2 /= u ;
      return *this;
    }
    FAD2d& operator /=(const long double &u) {
      value/=u;
      grad /=u;
      grad2 /= u ;
      return *this;
    }


    bool operator ==(const FAD2d &u)const  {
      return ((value==u.value)?true:false);
    }
    bool operator !=(const FAD2d &u) const {
      return ((value!=u.value)?true:false);
    }
    bool operator >(const FAD2d &u) const {
      return ((value>u.value)?true:false);
    }
    bool operator <(const FAD2d &u) const {
      return ((value<u.value)?true:false);
    }
    bool operator >=(const FAD2d &u) const {
      return ((value>=u.value)?true:false);
    }
    bool operator <=(const FAD2d &u) const {
      return ((value<=u.value)?true:false);
    }
    bool operator ==(const int &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const int &u) const {
      return ((value!=u)?true:false);
    }
    bool operator >(const int &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const int &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const int &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const int &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const float &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const float &u) const {
      return ((value!=u)?true:false);
    }
    bool operator >(const float &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const float &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const float &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const float &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const double &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const double &u) const {
      return ((value!=u)?true:false);
    }
    bool operator  >(const double &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const double &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const double &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const double &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const long double &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const long double &u) const {
      return ((value!=u)?true:false);
    }
    bool operator  >(const long double &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const long double &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const long double &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const long double &u) const {
      return ((value<=u)?true:false);
    }

  };

  inline std::ostream &operator <<(std::ostream& stream, const FAD2d &u) {
    if(u.grad != 0 || u.grad2 != 0)
      stream << u.value << '^' << u.grad << "^^" << u.grad2 ;
    else
      stream << u.value ;
    return stream;
  }
  inline std::istream &operator >> (std::istream& stream, FAD2d &u) {
    u.grad2 = 0 ;
    u.grad = 0 ;
    stream >> u.value;
    if(stream.peek() == '^') {
      stream.get() ;
      stream >> u.grad ;
      if(stream.peek() == '^') {
	stream.get() ;
	if(stream.peek() == '^') 
	  stream.get() ;
	stream >> u.grad2 ;
      }
    }

    return stream;
  }

#ifdef NO_CMATH
  inline double ceil(const FAD2d u) {
    return ::ceil(u.value) ;
  }
  inline double floor(const FAD2d u) {
    return ::floor(u.value) ;
  }

  inline FAD2d erf(const FAD2d u) {
    double ef = ::erf(u.value) ;
    double efp = (2./sqrt(M_PI))*::exp(-u.value*u.value) ;
    double efpp = -2.*u.value*efp ;
    return FAD2d(ef,u.grad*efp,u.grad*u.grad*efpp+efp*u.grad2) ;
  }

  inline FAD2d sin(const FAD2d u) {
    double su = ::sin(u.value) ;
    double cu = ::cos(u.value) ;
    return FAD2d(su, u.grad*cu, cu*u.grad2 - u.grad*u.grad*su);
  }
  inline FAD2d cos(const FAD2d u) {
    double su = ::sin(u.value) ;
    double cu = ::cos(u.value) ;
    return FAD2d(cu, -u.grad*su, -su*u.grad2 - u.grad*u.grad*cu);
  }
  inline FAD2d tan(const FAD2d u) {
    double tanu = ::tan(u.value) ;
    double secu = 1./::cos(u.value) ;
    double secu2 = secu*secu ;
    return FAD2d(tanu, u.grad*secu2,
		 u.grad2*secu2+2.*u.grad*u.grad*secu2*tanu) ;
  }


  inline FAD2d exp(const FAD2d u) {
    double eu = ::exp(u.value) ;
    return FAD2d(eu, u.grad*eu, eu*u.grad2+u.grad*u.grad*eu);
  }
  inline FAD2d log(const FAD2d u) {
    return FAD2d(::log(u.value), u.grad/u.value, 
		 u.grad2/u.value - u.grad*u.grad/(u.value*u.value));
  }
  inline FAD2d log10(const FAD2d u) {
    return FAD2d(::log(u.value), u.grad/u.value,
		 u.grad2/u.value - u.grad*u.grad/(u.value*u.value))/::log(10.0);
  }

  inline FAD2d fabs(FAD2d u) {
    return FAD2d(::fabs(u.value), 
		 u.grad*( (u.value<0.0)?-1.0:1.0 ),
		 u.grad2*( (u.value<0.0)?-1.0:1.0 )) ;
  }
  inline FAD2d abs(FAD2d u) {
    return FAD2d(::fabs(u.value), 
		 u.grad*( (u.value<0.0)?-1.0:1.0 ),
		 u.grad2*( (u.value<0.0)?-1.0:1.0 )) ;
  }
  /* MPGCOMMENT [05-12-2017 15:04] ---> POW */
  inline FAD2d pow(const FAD2d u, const int k) {
    return FAD2d(::pow(u.value, k), 
		 (double)k * ::pow(u.value, k-1)*u.grad,
		 k*((k - 1)*(::pow(u.value, k - 2)*::pow(u.grad, 2)) + ::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const float k) {
    return FAD2d(::pow(u.value, k), 
		 (double)k * ::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(::pow(u.value, k - 2)*::pow(u.grad, 2)) + ::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const double k) {
    return FAD2d(::pow(u.value, k), 
		 (double)k * ::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(::pow(u.value, k - 2)*::pow(u.grad, 2)) + ::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const long double k) {
    return FAD2d(::pow(u.value, k), 
		 (double)k * ::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(::pow(u.value, k - 2)*::pow(u.grad, 2)) + ::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d sqrt(const FAD2d u) {
    double su = ::sqrt(u.value) ;
    return FAD2d(su,0.5*u.grad/max(su,1e-30),
		 0.5*u.grad2/(max(su,1e-30)) - 0.25*u.grad*u.grad/(max(su*su*su,1e-30))) ;
  }
  inline FAD2d pow(const FAD2d k, const FAD2d u) {
    double kpu = ::pow(k.value,u.value) ;
    double kpu1 = ::pow(k.value,u.value-1.) ;
    double kpu2 = ::pow(k.value,u.value-2.) ;
    double lxu = ::log(k.value) ;
    return FAD2d(kpu,kpu1*k.grad*u.value + kpu*lxu*u.grad,
		 u.value*(k.grad*(kpu2*(k.grad*(u.value - 1)) + 
				  kpu1*(lxu*u.grad)) + kpu1*k.grad2) +
		 2*(kpu1*(k.grad*u.grad)) + 
		 lxu*(u.grad*(kpu1*(k.grad*u.value) + 
			      kpu*(lxu*u.grad)) + kpu*u.grad2)) ;

  }
  inline FAD2d pow(const int k, const FAD2d u) {
    double kpu = ::pow(k,u.value) ;
    double lnk = ::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const float k, const FAD2d u) {
    double kpu = ::pow(k,u.value) ;
    double lnk = ::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const double k, const FAD2d u) {
    double kpu = ::pow(k,u.value) ;
    double lnk = ::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const long double k, const FAD2d u) {
    double kpu = ::pow(k,u.value) ;
    double lnk = ::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }

  inline FAD2d sinh(const FAD2d u) {
    double expu = ::exp(u.value) ;
    double expmu = ::exp(-u.value) ;
    return FAD2d(::sinh(u.value),
		 0.5*u.grad*(expu + expmu),
		 0.5*(expmu*(u.grad2-u.grad*u.grad)+expu*(u.grad2+u.grad*u.grad))) ;
  }
  inline FAD2d cosh(const FAD2d u) {
    double expu = ::exp(u.value) ;
    double expmu = ::exp(-u.value) ;
    return FAD2d(::cosh(u.value),
		 0.5*u.grad*(expu - expmu),
		 0.5*(expmu*(u.grad*u.grad-u.grad2)+expu*(u.grad2+u.grad*u.grad))) ;
  }

  inline FAD2d tanh(const FAD2d u) {
    double ex = ::exp(min(u.value,350.0)) ;
    double exm = ::exp(min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    double rex = dex/sex ;
    double rex2 = rex*rex ;
    return FAD2d(::tanh(u.value),u.grad*(1.-rex2),
		 u.grad*u.grad*2.*(rex2-1.)*rex+
		 u.grad2*(1.-rex2)) ;
  }
  inline FAD2d asin(const FAD2d u) {
    double rsrt = 1./::sqrt(1.-u.value*u.value) ;
    return FAD2d(::asin(u.value), u.grad*rsrt,
		 u.grad2*rsrt + u.grad*u.grad*u.value*rsrt*rsrt*rsrt);
  }
  inline FAD2d acos(const FAD2d u) {
    double rsrt = 1./::sqrt(1.-u.value*u.value) ;
    return FAD2d(::acos(u.value), -u.grad*rsrt,
		 -u.grad2*rsrt - u.grad*u.grad*u.value*rsrt*rsrt*rsrt);
  }
  inline FAD2d atan(const FAD2d u) {
    double rval = 1./(1.+u.value*u.value) ;
    return FAD2d(::atan(u.value), u.grad*rval,
		 u.grad2*rval-2.*u.grad*u.grad*u.value*rval*rval);
  }
  //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // This will not work in general
  //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  inline FAD2d atan2(const FAD2d u, const FAD2d v) {  
    return atan(u/v);
  }

  inline FAD2d asinh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double strm = ::sqrt(uv2+1.) ;
    const double uvstrm = uv+strm ;
    const double uvstrm2 = uvstrm*uvstrm ;
    return FAD2d(::asinh(uv), u.grad*(uv/strm+1.)/uvstrm,
		 ((1. / strm - uv2 / ::pow(strm, 3)) / uvstrm 
		  - ::pow(uv / strm + 1., 2) / uvstrm2)*u.grad*u.grad
		 + u.grad2*(uv / strm + 1.) / uvstrm) ;
  }
  inline FAD2d acosh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double strm = ::sqrt(uv*uv-1.) ;
    const double uvstrm = uv+strm ;
    const double uvstrm2 = uvstrm*uvstrm ;
    return FAD2d(::acosh(uv), u.grad*(uv/strm +1.)/uvstrm,
		 ((1 / strm - uv2 / ::pow(strm, 3)) / uvstrm 
		  - ::pow(uv / strm + 1., 2) / uvstrm2)*u.grad*u.grad 
		 + u.grad2*(uv / strm + 1.) / uvstrm) ;	 
  }
  inline FAD2d atanh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double factor = .5*(1./(1.+uv)+1./(1.-uv)) ;
    return FAD2d(::atanh(uv), .5*u.grad*factor,
		 0.5*(1./ (1. - 2.*uv + uv2) - 1./ (1. + 2.*uv + uv2))*u.grad*u.grad + u.grad2*factor) ;
  }


#else
  inline double ceil(const FAD2d u) {
    return std::ceil(u.value) ;
  }
  inline double floor(const FAD2d u) {
    return std::floor(u.value) ;
  }

  inline FAD2d erf(const FAD2d u) {
    double ef = ::erf(u.value) ;
    double efp = (2./sqrt(M_PI))*std::exp(-u.value*u.value) ;
    double efpp = -2.*u.value*efp ;
    return FAD2d(ef,u.grad*efp,u.grad*u.grad*efpp+efp*u.grad2) ;
  }
  inline FAD2d sin(const FAD2d u) {
    double su = std::sin(u.value) ;
    double cu = std::cos(u.value) ;
    return FAD2d(su, u.grad*cu, cu*u.grad2 - u.grad*u.grad*su);
  }
  inline FAD2d cos(const FAD2d u) {
    double su = std::sin(u.value) ;
    double cu = std::cos(u.value) ;
    return FAD2d(cu, -u.grad*su, -su*u.grad2 - u.grad*u.grad*cu);
  }
  inline FAD2d tan(const FAD2d u) {
    double tanu = std::tan(u.value) ;
    double secu = 1./std::cos(u.value) ;
    double secu2 = secu*secu ;
    return FAD2d(tanu, u.grad*secu2,
		 u.grad2*secu2+2.*u.grad*u.grad*secu2*tanu) ;
  }

  inline FAD2d exp(const FAD2d u) {
    double eu = std::exp(u.value) ;
    return FAD2d(eu, u.grad*eu, eu*u.grad2+u.grad*u.grad*eu);
  }
  inline FAD2d log(const FAD2d u) {
    return FAD2d(std::log(u.value), u.grad/u.value, 
		 u.grad2/u.value - u.grad*u.grad/(u.value*u.value));
  }
  inline FAD2d log10(const FAD2d u) {
    return FAD2d(std::log(u.value), u.grad/u.value,
		 u.grad2/u.value - u.grad*u.grad/(u.value*u.value))/std::log(10.0);
  }

  inline FAD2d fabs(FAD2d u) {
    return FAD2d(std::fabs(u.value), 
		 u.grad*( (u.value<0.0)?-1.0:1.0 ),
		 u.grad2*( (u.value<0.0)?-1.0:1.0 )) ;
  }
  inline FAD2d abs(FAD2d u) {
    return FAD2d(std::fabs(u.value), 
		 u.grad*( (u.value<0.0)?-1.0:1.0 ),
		 u.grad2*( (u.value<0.0)?-1.0:1.0 )) ;
  }
  /* MPGCOMMENT [05-12-2017 15:04] ---> POW */
  inline FAD2d pow(const FAD2d u, const int k) {
    return FAD2d(std::pow(u.value, k), 
		 (double)k * std::pow(u.value, k-1)*u.grad,
		 k*((k - 1)*(std::pow(u.value, k - 2)*std::pow(u.grad, 2)) + std::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const float k) {
    return FAD2d(std::pow(u.value, double(k)), 
		 (double)k * std::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(std::pow(u.value, k-2.0)*std::pow(u.grad, 2)) + std::pow(u.value, k-1.0)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const double k) {
    return FAD2d(std::pow(u.value, k), 
		 (double)k * std::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(std::pow(u.value, k-2.0)*std::pow(u.grad, 2)) + std::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d pow(const FAD2d u, const long double ki) {
    double k = double(ki) ;
    return FAD2d(std::pow(u.value, k), 
		 (double)k * std::pow(u.value, k-1.0)*u.grad,
		 k*((k - 1)*(std::pow(u.value, k - 2)*std::pow(u.grad, 2)) + std::pow(u.value, k - 1)*u.grad2)) ;
  }
  inline FAD2d sqrt(const FAD2d u) {

    double su = std::sqrt(u.value) ;
    return FAD2d(su,0.5*u.grad/max(su,1e-30),
		 0.5*u.grad2/(max(su,1e-30)) - 0.25*u.grad*u.grad/(max(su*su*su,1e-30))) ;
  }
  inline FAD2d pow(const FAD2d k, const FAD2d u) {
    double kpu = std::pow(k.value,u.value) ;
    double kpu1 = std::pow(k.value,u.value-1.) ;
    double kpu2 = std::pow(k.value,u.value-2.) ;
    double lxu = std::log(k.value) ;
    return FAD2d(kpu,kpu1*k.grad*u.value + kpu*lxu*u.grad,
		 u.value*(k.grad*(kpu2*(k.grad*(u.value - 1)) + 
				  kpu1*(lxu*u.grad)) + kpu1*k.grad2) +
		 2*(kpu1*(k.grad*u.grad)) + 
		 lxu*(u.grad*(kpu1*(k.grad*u.value) + 
			      kpu*(lxu*u.grad)) + kpu*u.grad2)) ;

  }
  inline FAD2d pow(const int k, const FAD2d u) {
    double kpu = std::pow(k,u.value) ;
    double lnk = std::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const float k, const FAD2d u) {
    double kpu = std::pow(double(k),u.value) ;
    double lnk = std::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const double k, const FAD2d u) {
    double kpu = std::pow(k,u.value) ;
    double lnk = std::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }
  inline FAD2d pow(const long double ki, const FAD2d u) {
    double k = double(ki) ;
    double kpu = std::pow(k,u.value) ;
    double lnk = std::log(double(k)) ;
    return FAD2d(kpu,kpu*lnk*u.grad,  
		 kpu*((lnk*lnk)*(u.grad*u.grad)) + kpu*(lnk*u.grad2)) ;
  }

  inline FAD2d sinh(const FAD2d u) {
    double expu = std::exp(u.value) ;
    double expmu = std::exp(-u.value) ;
    return FAD2d(std::sinh(u.value),
		 0.5*u.grad*(expu + expmu),
		 0.5*(expmu*(u.grad2-u.grad*u.grad)+expu*(u.grad2+u.grad*u.grad))) ;
  }
  inline FAD2d cosh(const FAD2d u) {
    double expu = std::exp(u.value) ;
    double expmu = std::exp(-u.value) ;
    return FAD2d(std::cosh(u.value),
		 0.5*u.grad*(expu - expmu),
		 0.5*(expmu*(u.grad*u.grad-u.grad2)+
		      expu*(u.grad2+u.grad*u.grad))) ;
  }

  inline FAD2d tanh(const FAD2d u) {
    double ex = std::exp(min(u.value,350.0)) ;
    double exm = std::exp(min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    double rex = dex/sex ;
    double rex2 = rex*rex ;
    return FAD2d(std::tanh(u.value),u.grad*(1.-rex2),
		 u.grad*u.grad*2.*(rex2-1.)*rex+
		 u.grad2*(1.-rex2)) ;
  }
  inline FAD2d asin(const FAD2d u) {
    double rsrt = 1./std::sqrt(1.-u.value*u.value) ;
    return FAD2d(std::asin(u.value), u.grad*rsrt,
		 u.grad2*rsrt + u.grad*u.grad*u.value*rsrt*rsrt*rsrt);
  }
  inline FAD2d acos(const FAD2d u) {
    double rsrt = 1./std::sqrt(1.-u.value*u.value) ;
    return FAD2d(std::acos(u.value), -u.grad*rsrt,
		 -u.grad2*rsrt - u.grad*u.grad*u.value*rsrt*rsrt*rsrt);
  }
  inline FAD2d atan(const FAD2d u) {
    double rval = 1./(1.+u.value*u.value) ;
    return FAD2d(std::atan(u.value), u.grad*rval,
		 u.grad2*rval-2.*u.grad*u.grad*u.value*rval*rval);
  }
  //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // This will not work in general
  //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  inline FAD2d atan2(const FAD2d u, const FAD2d v) {  
    return atan(u/v);
  }

  inline FAD2d asinh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double strm = std::sqrt(uv2+1.) ;
    const double uvstrm = uv+strm ;
    const double uvstrm2 = uvstrm*uvstrm ;
    return FAD2d(::asinh(uv), u.grad*(uv/strm+1.)/uvstrm,
		 ((1. / strm - uv2 / std::pow(strm, 3)) / uvstrm 
		  - std::pow(uv / strm + 1., 2) / uvstrm2)*u.grad*u.grad
		 + u.grad2*(uv / strm + 1.) / uvstrm) ;
  }
  inline FAD2d acosh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double strm = std::sqrt(uv*uv-1.) ;
    const double uvstrm = uv+strm ;
    const double uvstrm2 = uvstrm*uvstrm ;
    return FAD2d(::acosh(uv), u.grad*(uv/strm +1.)/uvstrm,
		 ((1 / strm - uv2 / std::pow(strm, 3)) / uvstrm 
		  - std::pow(uv / strm + 1., 2) / uvstrm2)*u.grad*u.grad 
		 + u.grad2*(uv / strm + 1.) / uvstrm) ;	 
  }
  inline FAD2d atanh(const FAD2d u) {
    const double uv = u.value ;
    const double uv2 = uv*uv ;
    const double factor = .5*(1./(1.+uv)+1./(1.-uv)) ;
    return FAD2d(::atanh(uv), .5*u.grad*factor,
		 0.5*(1./ (1. - 2.*uv + uv2) -
		      1./ (1. + 2.*uv + uv2))*u.grad*u.grad + u.grad2*factor) ;
  }


#endif
  

  inline FAD2d operator +(const double u,const FAD2d v)
  { return FAD2d(u+v.value, v.grad, v.grad2); }
  inline FAD2d operator -(const double u,const FAD2d v) 
  { return FAD2d(u-v.value, -v.grad, -v.grad2); }
  inline FAD2d operator *(const double u,const FAD2d v) 
    { return FAD2d(u*v.value,  v.grad*u, v.grad2*u ); }
  inline FAD2d operator /(const double &u,const FAD2d &v) 
  { return FAD2d(u/v.value, -u*v.grad/v.value/v.value,
		 2.*u*v.grad*v.grad/(v.value*v.value*v.value) -
		 u*v.grad2/(v.value*v.value)); }

  class FADd {
  public:

    double value;
    double grad;
#ifdef LIVING_DANGEROUS
    /* MPGCOMMENT [05-12-2017 13:50] ---> type-cast */
    explicit operator char() {return (char)value;}
    explicit operator bool() {return (bool)value;}
    explicit operator int() {return (int)value;}
    explicit operator float() {return (float)value;}
    operator double() {return (double)value;}
    explicit operator FADd() {return FADd(value,grad);}
    explicit operator long double() const {return (long double)value;}

    //    operator double*()  {return &value;}/* MPGCOMMENT [05-15-2017 10:51] ---> only one used, others set for completeness */
    //    operator const double*()  {return (const double*)&value;}/* MPGCOMMENT [05-15-2017 10:51] ---> only one used, others set for completeness */
#endif
    
    /* MPGCOMMENT [05-12-2017 13:49] ---> constructor */
#define gradcheck() //if(grad!=0) debugger_() ; //fatal(grad!=0.0) ;
    template< class T1>
      FADd(T1 a0): value(a0), grad(0.0){ gradcheck();}/* MPGCOMMENT [05-15-2017 14:46] ---> used in chem, chem::Temperature */
  FADd(int a0         ,int b0)        : value(a0), grad(b0) { gradcheck();}
  FADd(int a0         ,float b0)        : value(a0), grad(b0) { gradcheck();}
  FADd(int a0         ,double b0)        : value(a0), grad(b0) { gradcheck();}
  FADd(int a0         ,long double b0)        : value(a0), grad(b0) { gradcheck();}
  FADd(float a0      ,int b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(float a0      ,float b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(float a0      ,double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(float a0      ,long double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(double a0      ,int b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(double a0      ,float b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(double a0      ,double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(double a0      ,long double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(long double a0      ,int b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(long double a0      ,float b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(long double a0      ,double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(long double a0      ,long double b0)       : value(a0), grad(b0) { gradcheck();}
  FADd(const FAD2d &u): value(u.value), grad(u.grad) { gradcheck();}
    //    explicit FADd(bool &u) : value((u)?1.0:0.0), grad(0.0) { gradcheck();}
    FADd(const FADd &u) : value(u.value), grad(u.grad) { gradcheck();}
    //  explicit FADd(int a0)        : value(a0), grad(0.0) { gradcheck();}
    //  explicit FADd(float a0)        : value(a0), grad(0.0) { gradcheck();}
    //  explicit FADd(double a0)        : value(a0), grad(0.0) { gradcheck();}
    //  explicit FADd(long double a0)        : value(a0), grad(0.0) { gradcheck();}
  FADd() : value(0.0), grad(0.0) { gradcheck();}


    //    friend void swap(double& a, FADd& b) {
    //      using std::swap;
    //      swap(a, b.value);
    //    }

    /**************************************************************/
    /* MPGCOMMENT [05-12-2017 10:38] ---> FADd --- FADd OPERATORS */

    /* MPGCOMMENT [05-12-2017 10:35] ---> ASSIGNMENT */ 
    /* template<class T> Loci::const_store<Loci::vector3d<FADd> >& operator =(Loci::const_store<Loci::vector3d<T> > &u) { */
    /*   value = u; */
    /*   grad = 0.0; */
    /*   return *this; */
    /* } */

    //    template<class T> FADd& operator =(const T &u) {
    //      value = u;
    //      grad = 0.0;
    //      return *this;
    //    }
    FADd& operator =(const FADd &u) {
      value = u.value;
      grad = u.grad;
      return *this;
    }
    FADd& operator =(const int &u) {
      value = u;
      grad = 0.0;
      return *this;
    }
    FADd& operator =(const float &u) {
      value = u;
      grad = 0.0;
      return *this;
    }
    FADd& operator =(const double &u) {
      value = u;
      grad = 0.0;
      return *this;
    }
    FADd& operator =(const long double &u) {
      value = u;
      grad = 0.0;
      return *this;
    }


    /* MPGCOMMENT [05-12-2017 10:35] ---> ADDITION */
    //    template<class T> FADd operator +(const T &v) const{
    //      return FADd(value+(double)v, grad);
    //    }
    FADd operator +(const FADd &v) const{
      return FADd(value+v.value, grad+v.grad);
    }
    FADd operator +() const{
      return FADd(value, grad);
    }
    FADd operator +(const int &v) const{
      return FADd(value+(double)v, grad);
    }
    FADd operator +(const float &v) const{
      return FADd(value+(double)v, grad);
    }
    FADd operator +(const double &v) const{
      return FADd(value+(double)v, grad);
    }
    FADd operator +(const long double &v) const{
      return FADd(value+(long double)v, grad);
    }
    //    template<class T> FADd& operator +=(const T &u) {
    //      value+=(double)u;
    //      return *this;
    //    }
    FADd& operator +=(const FADd &u) {
      value+=u.value;
      grad+=u.grad;
      return *this;
    }
    FADd& operator +=(const int &u) {
      value+=u;
      return *this;
    }
    FADd& operator +=(const float &u) {
      value+=u;
      return *this;
    }
    FADd& operator +=(const double &u) {
      value+=u;
      return *this;
    }
    FADd& operator +=(const long double &u) {
      value+=u;
      return *this;
    }


    /* MPGCOMMENT [05-12-2017 10:35] ---> SUBTRACTION */
    //    template<class T> FADd operator -(const T &v) const{
    //      return FADd(value-(double)v, grad);
    //    }
    FADd operator -(const FADd &v) const{
      return FADd(value-v.value, grad-v.grad);
    }
    FADd operator -() const {
      return FADd(-value, -grad);
    }
    FADd operator -(const int &v) const{
      return FADd(value-(double)v, grad);
    }
    FADd operator -(const float &v) const{
      return FADd(value-(double)v, grad);
    }
    FADd operator -(const double &v) const{
      return FADd(value-(double)v, grad);
    }
    FADd operator -(const long double &v) const{
      return FADd(value-(long double)v, grad);
    }
    //    template<class T> FADd& operator -=(const T &u) {
    //      value-=(double)u;
    //      return *this;
    //    }
    FADd& operator -=(const FADd &u) {
      value-=u.value;
      grad-=u.grad;
      return *this;
    }
    FADd& operator -=(const int &u) {
      value-=u;
      return *this;
    }
    FADd& operator -=(const float &u) {
      value-=u;
      return *this;
    }
    FADd& operator -=(const double &u) {
      value-=u;
      return *this;
    }
    FADd& operator -=(const long double &u) {
      value-=u;
      return *this;
    }

    /* MPGCOMMENT [05-12-2017 10:36] ---> MULTIPLICATION */
    //    template<class T> FADd operator *(const T &v) const{
    //      return FADd(value*(double)v, grad);
    //    }
    FADd operator *(const FADd &v) const {
      return FADd(value*v.value,  grad*v.value + v.grad*value);
    }
    FADd operator *(const int &v) const{
      return FADd(value*(double)v,  grad*(double)v);
    }
    FADd operator *(const float &v) const{
      return FADd(value*(double)v,  grad*(double)v);
    }
    FADd operator *(const double &v) const{
      return FADd(value*(double)v,  grad*(double)v);
    }
    FADd operator *(const long double &v) const{
      return FADd(value*(long double)v,  grad*(long double)v);
    }
    //    template<class T> FADd& operator *=(const T &u) {
    //      value*=(double)u;
    //      return *this;
    //    }
    FADd& operator *=(const FADd &u) {
      grad = grad*u.value + u.grad * value;
      value*=u.value;
      return *this;
    }
    FADd& operator *=(const int &u) {
      value*=u;
      grad *=u;// + (u.grad is zero) * value;
      return *this;
    }
    FADd& operator *=(const float &u) {
      value*=u;
      grad *=u;// + (u.grad is zero) * value;
      return *this;
    }
    FADd& operator *=(const double &u) {
      value*=u;
      grad *=u;// + (u.grad is zero) * value;
      return *this;
    }
    FADd& operator *=(const long double &u) {
      value*=u;
      grad *=u;// + (u.grad is zero) * value;
      return *this;
    }

    /* MPGCOMMENT [05-12-2017 10:36] ---> DIVISION */
    //    template<class T> FADd operator /(const T &v) const{
    //      return FADd(value/(double)v, grad);
    //    }
    FADd operator /(const FADd &v) const{
      return FADd(value/v.value, (grad*v.value - value*v.grad)/v.value/v.value);
    }
    FADd operator /(const int &v) const{
      return FADd(value/v, (grad/(double)v));
    }
    FADd operator /(const float &v) const{
      return FADd(value/v, (grad/(double)v));
    }
    FADd operator /(const double &v) const{
      return FADd(value/v, (grad/(double)v));
    }
    FADd operator /(const long double &v) const{
      return FADd(value/v, (grad/(long double)v));
    }
    //    template<class T> FADd& operator /=(const T &u) {
    //      value/=(double)u;
    //      return *this;
    //    }
    FADd& operator /=(const FADd &u) {
      grad = (grad*u.value - value * u.grad) / u.value / u.value;
      value/=u.value;
      return *this;
    }
    FADd& operator /=(const int &u) {
      value/=u;
      grad /=u;// (grad*u.value /*- value * u.grad=0*/) / u / u;
      return *this;
    }
    FADd& operator /=(const float &u) {
      value/=u;
      grad /=u;// (grad*u.value /*- value * u.grad=0*/) / u / u;
      return *this;
    }
    FADd& operator /=(const double &u) {
      value/=u;
      grad /=u;// (grad*u.value /*- value * u.grad=0*/) / u / u;
      return *this;
    }
    FADd& operator /=(const long double &u) {
      value/=u;
      grad /=u;// (grad*u.value /*- value * u.grad=0*/) / u / u;
      return *this;
    }


    /* MPGCOMMENT [05-12-2017 10:37] ---> COMPARISON OPERATORS */
    //    template<class T> bool operator ==(const T &u) const {
    //      return ((value==(double)u)?true:false);
    //    }
    //    template<class T> bool operator !=(const T &u) const {
    //      return ((value!=(double)u)?true:false);
    //    }
    //    template<class T> bool operator >(const T &u) const {
      //      return ((value>(double)u)?true:false);
    //    }
    //    template<class T> bool operator <(const T &u) const {
    //      return ((value<(double)u)?true:false);
    //    }
    //    template<class T> bool operator >=(const T &u) const {
    //      return ((value>=(double)u)?true:false);
    //    }
    //    template<class T> bool operator <=(const T &u) const {
    //      return ((value<=(double)u)?true:false);
    //    }
    bool operator ==(const FADd &u)const  {
      return ((value==u.value)?true:false);
    }
    bool operator !=(const FADd &u) const {
      return ((value!=u.value)?true:false);
    }
    bool operator >(const FADd &u) const {
      return ((value>u.value)?true:false);
    }
    bool operator <(const FADd &u) const {
      return ((value<u.value)?true:false);
    }
    bool operator >=(const FADd &u) const {
      return ((value>=u.value)?true:false);
    }
    bool operator <=(const FADd &u) const {
      return ((value<=u.value)?true:false);
    }
    bool operator ==(const int &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const int &u) const {
      return ((value!=u)?true:false);
    }
    bool operator >(const int &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const int &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const int &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const int &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const float &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const float &u) const {
      return ((value!=u)?true:false);
    }
    bool operator >(const float &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const float &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const float &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const float &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const double &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const double &u) const {
      return ((value!=u)?true:false);
    }
    bool operator  >(const double &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const double &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const double &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const double &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const long double &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const long double &u) const {
      return ((value!=u)?true:false);
    }
    bool operator  >(const long double &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const long double &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const long double &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const long double &u) const {
      return ((value<=u)?true:false);
    }



    /* Array<T,n> &operator +=(const Array<T,n> &v) */
    /*   { for(size_t i=0;i<n;++i) x[i] += v.x[i] ; return *this ; } */
    /* Array<T,n> &operator -=(const Array<T,n> &v) */
    /*   { for(size_t i=0;i<n;++i) x[i] -= v.x[i] ; return *this ; } */
    /* Array<T,n> &operator *=(const Array<T,n> &v) */
    /*   { for(size_t i=0;i<n;++i) x[i] *= v.x[i] ; return *this ; } */
    /* Array<T,n> &operator /=(const Array<T,n> &v) */
    /*   { for(size_t i=0;i<n;++i) x[i] /= v.x[i] ; return *this ;} */


  };

  inline std::ostream &operator <<(std::ostream& stream, const FADd &u) {
    if(u.grad != 0) 
      stream << u.value << '^' << u.grad ;
    else
      stream << u.value ;
      
    return stream;
  }
  inline std::istream &operator >> (std::istream& stream, FADd &u) {
    u.grad = 0 ;
    stream >> u.value;
    if(stream.peek() == '^') {
      stream.get() ;
      stream >> u.grad ;
    }
    return stream;
  }
  
#ifdef NO_CMATH
  inline double ceil(const FADd u) {
    return ::ceil(u.value) ;
  }
  inline double floor(const FADd u) {
    return ::floor(u.value) ;
  }

  inline FADd erf(const FADd u) {
    double ef = ::erf(u.value) ;
    double efp = (2./sqrt(M_PI))*::exp(-u.value*u.value) ;
    return FADd(ef,u.grad*efp) ;
  }
  inline FADd sin(const FADd u) {
    return FADd(::sin(u.value), u.grad*::cos(u.value));
  }
  inline FADd cos(const FADd u) {
    return FADd(::cos(u.value), -u.grad*::sin(u.value));
  }
  inline FADd tan(const FADd u) {
    return FADd(::tan(u.value), u.grad/pow(::cos(u.value),2)) ;
  }
  inline FADd exp(const FADd u) {
    return FADd(::exp(u.value), u.grad*::exp(u.value));
  }
  inline FADd log(const FADd u) {
    return FADd(::log(u.value), u.grad/u.value);
  }
  inline FADd log10(const FADd u) {
    return FADd(::log(u.value), u.grad/u.value)/::log(10.0);
  }

  inline FADd fabs(FADd u) {
    return FADd(::fabs(u.value), 
		u.grad*( (u.value<0.0)?-1.0:1.0 ) );
  }
  inline FADd abs(FADd u) {
    return FADd(::fabs(u.value), 
		u.grad*( (u.value<0.0)?-1.0:1.0 ) );
  }
  /* MPGCOMMENT [05-12-2017 15:04] ---> POW */
  inline FADd pow(const FADd u, const int k) {
    return FADd(::pow(u.value, k), 
		(double)k * ::pow(u.value, (double)k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const float k) {
    return FADd(::pow(u.value, k), 
		(double)k * ::pow(u.value, (double)k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const double k) {
    return FADd(::pow(u.value, k), 
		k * ::pow(u.value, k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const long double k) {
    return FADd(::pow(u.value, k), 
		(double)k * ::pow(u.value,  (double)k-1.0)*u.grad );
  }
  inline FADd sqrt(const FADd u) {
    double su = ::sqrt(u.value) ;
    return FADd(su,0.5*u.grad/max(su,1e-30)) ;
  }
  inline FADd pow(const FADd k, const FADd u) {
    double kpu = ::pow(k.value,u.value) ;
    return FADd(kpu,::pow(k.value,u.value-1.)*k.grad*u.value +
		kpu*::log(k.value)*u.grad) ;
  }
  inline FADd pow(const int k, const FADd u) {
    double kpu = ::pow(k,u.value) ;
    return FADd(kpu,kpu*::log(double(k))*u.grad) ;
  }
  inline FADd pow(const float k, const FADd u) {
    double kpu = ::pow(k,u.value) ;
    return FADd(kpu,kpu*::log(double(k))*u.grad) ;
  }
  inline FADd pow(const double k, const FADd u) {
    double kpu = ::pow(k,u.value) ;
    return FADd(kpu,kpu*::log(k)*u.grad) ;
  }
  inline FADd pow(const long double k, const FADd u) {
    double kpu = ::pow(k,u.value) ;
    return FADd(kpu,kpu*::log(k)*u.grad) ;
  }


  inline FADd sinh(const FADd u) {
    return FADd(::sinh(u.value),
		0.5*u.grad*(::exp(u.value) + ::exp(-(u.value)))) ;
  }
  inline FADd cosh(const FADd u) {
    return FADd(::cosh(u.value),
		0.5*u.grad*(::exp(u.value) - ::exp(-(u.value)))) ;
  }
  inline FADd tanh(const FADd u) {
    double ex = ::exp(min(u.value,350.0)) ;
    double exm = ::exp(min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    return FADd(::tanh(u.value),u.grad*(1.-dex*dex/(sex*sex))) ;
  }
  inline FADd asin(const FADd u) {
    return FADd(::asin(u.value), u.grad/::sqrt(1.0-u.value*u.value) );
  }
  inline FADd acos(const FADd u) {
    return FADd(::acos(u.value), -u.grad/::sqrt(1.0-u.value*u.value) );
  }
  inline FADd atan(const FADd u) {
    return FADd(::atan(u.value), u.grad/(1.0+u.value*u.value) );
  }
  // This will not work in general
  inline FADd atan2(const FADd u, const FADd v) {  
    return atan(u/v);
  }

  inline FADd asinh(const FADd u) {
    double strm = ::sqrt(u.value*u.value+1.) ;
    double uvstrm = u.value+strm ;
    return FADd(::asinh(u.value), u.grad*(u.value/strm+1.)/uvstrm) ;
  }
  inline FADd acosh(const FADd u) {
    double strm = ::sqrt(u.value*u.value-1.) ;
    double uvstrm = u.value+strm ;
    return FADd(::acosh(u.value), u.grad*(u.value/strm +1.)/uvstrm) ;
  }
  inline FADd atanh(const FADd u) {
    const double uv = u.value ;
    return FADd(::atanh(uv), .5*u.grad*(1./(1.+uv)+1./(1.-uv))) ;
  }


#else
  inline double ceil(const FADd u) {
    return std::ceil(u.value) ;
  }
  inline double floor(const FADd u) {
    return std::floor(u.value) ;
  }

  inline FADd erf(const FADd u) {
    double ef = ::erf(u.value) ;
    double efp = (2./sqrt(M_PI))*std::exp(-u.value*u.value) ;
    return FADd(ef,u.grad*efp) ;
  }
  inline FADd sin(const FADd u) {
    return FADd(std::sin(u.value), u.grad*std::cos(u.value));
  }
  inline FADd cos(const FADd u) {
    return FADd(std::cos(u.value), -u.grad*std::sin(u.value));
  }
  inline FADd tan(const FADd u) {
    return FADd(std::tan(u.value), u.grad/pow(std::cos(u.value),2)) ;
  }

  inline FADd exp(const FADd u) {
    return FADd(std::exp(u.value), u.grad*std::exp(u.value));
  }
  inline FADd log(const FADd u) {
    return FADd(std::log(u.value), u.grad/u.value);
  }
  inline FADd log10(const FADd u) {
    return FADd(std::log(u.value), u.grad/u.value)/std::log(10.0);
  }
  inline FADd abs(const FADd u) {
    return FADd(std::fabs(u.value), 
		u.grad*( (u.value<0.0)?-1.0:1.0 ) );
  }
  inline FADd fabs(const FADd u) {
    return FADd(std::fabs(u.value), 
		u.grad*( (u.value<0.0)?-1.0:1.0 ) );
  }
  /* MPGCOMMENT [05-12-2017 15:04] ---> POW */
  inline FADd pow(const FADd u, const int k) {
    return FADd(std::pow(u.value, k), 
		(double)k * std::pow(u.value, (double)k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const float k) {
    return FADd(std::pow(u.value, double(k)), 
		(double)k * std::pow(u.value, (double)k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const double k) {
    return FADd(std::pow(u.value, k), 
		k * std::pow(u.value, k-1.0)*u.grad );
  }
  inline FADd pow(const FADd u, const long double k) {
    return FADd(std::pow(u.value,  double(k)), 
		(double)k * std::pow(u.value,  (double)k-1.0)*u.grad );
  }
  inline FADd sqrt(const FADd u) {
    double su = std::sqrt(u.value) ;
    return FADd(su,0.5*u.grad/max(su,1e-30)) ;
  }
  inline FADd pow(const FADd k, const FADd u) {
    double kpu = std::pow(k.value,u.value) ;
    return FADd(kpu,std::pow(k.value,u.value-1.)*k.grad*u.value +
		kpu*std::log(k.value)*u.grad) ;
  }
  inline FADd pow(const int k, const FADd u) {
    double kpu = std::pow(k,u.value) ;
    return FADd(kpu,kpu*std::log(double(k))*u.grad) ;
  }
  inline FADd pow(const float k, const FADd u) {
    double kpu = std::pow(double(k),u.value) ;
    return FADd(kpu,kpu*std::log(double(k))*u.grad) ;
  }
  inline FADd pow(const double k, const FADd u) {
    double kpu = std::pow(k,u.value) ;
    return FADd(kpu,kpu*std::log(k)*u.grad) ;
  }
  inline FADd pow(const long double k, const FADd u) {
    double kpu = std::pow(double(k),u.value) ;
    return FADd(kpu,kpu*std::log(k)*u.grad) ;
  }


  inline FADd sinh(const FADd u) {
    return FADd(std::sinh(u.value),
		0.5*u.grad*(std::exp(u.value) + std::exp(-u.value))) ;
  }
  inline FADd cosh(const FADd u) {
    return FADd(std::cosh(u.value),
		0.5*u.grad*(std::exp(u.value) - std::exp(-u.value))) ;
  }
  inline FADd tanh(const FADd u) {
    double ex = std::exp(min(u.value,350.0)) ;
    double exm = std::exp(min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    return FADd(std::tanh(u.value),
		u.grad*(1.-dex*dex/(sex*sex))) ;
  }

  inline FADd asin(const FADd u) {
    return FADd(std::asin(u.value), u.grad/std::sqrt(1.0-u.value*u.value) );
  }
  inline FADd acos(const FADd u) {
    return FADd(std::acos(u.value), -u.grad/std::sqrt(1.0-u.value*u.value) );
  }
  inline FADd atan(const FADd u) {
    return FADd(std::atan(u.value), u.grad/(1.0+u.value*u.value) );
  }
  // This will not work in general
  inline FADd atan2(const FADd u, const FADd v) {  
    return atan(u/v);
  }

  inline FADd asinh(const FADd u) {
    double strm = std::sqrt(u.value*u.value+1.) ;
    double uvstrm = u.value+strm ;
    return FADd(::asinh(u.value), u.grad*(u.value/strm+1.)/uvstrm) ;
  }
  inline FADd acosh(const FADd u) {
    double strm = std::sqrt(u.value*u.value-1.) ;
    double uvstrm = u.value+strm ;
    return FADd(::acosh(u.value), u.grad*(u.value/strm +1.)/uvstrm) ;
  }
  inline FADd atanh(const FADd u) {
    const double uv = u.value ;
    return FADd(::atanh(uv), .5*u.grad*(1./(1.+uv)+1./(1.-uv))) ;
  }


#endif
  

  inline FADd operator +(const double u,const FADd v)
  { return FADd(u+v.value, v.grad); }
  inline FADd operator -(const double u,const FADd v) 
  { return FADd(u-v.value, -v.grad); }
  inline FADd operator *(const double u,const FADd v) 
  { return FADd(u*v.value,  v.grad*u); }
  inline FADd operator /(const double &u,const FADd &v) 
  { return FADd(u/v.value, -u*v.grad/v.value/v.value); }

#ifndef MFAD_SIZE
#define MFAD_SIZE 1
#endif

  class MFADd {
  public:

    double value;
    double grad[MFAD_SIZE];
    static const size_t maxN=MFAD_SIZE;

    //    void setVecSize(size_t s) {
    //      this->maxN = s;
    //    }
    template< class T1>      
      MFADd(T1 a0) : value(a0), grad()
    { for(int i=0;i<MFAD_SIZE;++i) grad[i] =0 ; }
    MFADd() : value(0.0), grad()
    { for(int i=0;i<MFAD_SIZE;++i) grad[i] =0 ; }
    //explicit MFADd(int u) : value(u), grad(){ }
    //explicit MFADd(float u) : value(u), grad() { }
    //explicit MFADd(double u) : value(u), grad() { }
    //explicit MFADd(long double u) : value(u), grad() { }
    MFADd(double u, size_t n) : value(u), grad()
    { for(int i=0;i<MFAD_SIZE;++i) grad[i] =0 ; }
    //explicit MFADd(bool &u)  : value((u)?1.0:0.0), grad() { }
    MFADd(const double u, const double * g, const int n): value(u), grad() {
      for (size_t i=0; i<maxN; i++) this->grad[i]=g[i];      
    }
    MFADd(const MFADd &u, size_t n): value(u.value), grad() {
      for (size_t i=0; i<maxN; i++) this->grad[i]=u.grad[i];      
    }
    MFADd(const MFADd &u): value(u.value), grad() {
      for (size_t i=0; i<maxN; i++) this->grad[i]=u.grad[i];      
    }

    MFADd& operator =(const MFADd &u) {    
      this->value = u.value;
      for (size_t i=0; i<maxN; i++) this->grad[i]=u.grad[i];
      return *this;
    }
    MFADd& operator =(const int &u) {    
      this->value = u ;
      for (size_t i=0; i<maxN; i++) this->grad[i] = 0.0 ;
      return *this;
    }
    MFADd& operator =(const float &u) {    
      this->value = u ;
      for (size_t i=0; i<maxN; i++) this->grad[i] = 0.0 ;
      return *this;
    }
    MFADd& operator =(const double &u) {    
      this->value = u ;
      for (size_t i=0; i<maxN; i++) this->grad[i] = 0.0 ;
      return *this;
    }
    MFADd& operator =(const long double &u) {    
      this->value = u ;
      for (size_t i=0; i<maxN; i++) this->grad[i] = 0.0 ;
      return *this;
    }

    MFADd operator +(const MFADd &v) const{
      size_t s = maxN;
      MFADd out(this->value+v.value,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i]+v.grad[i];
      return out;
    }
    MFADd operator +() const { 
      size_t s = maxN;
      MFADd out(this->value,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i];
      return out; 
    }
    MFADd operator +(const int &v) const{
      size_t s = maxN;
      MFADd out(this->value+(double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i];
      return out;
    }
    MFADd operator +(const float &v) const{
      size_t s = maxN;
      MFADd out(this->value+(double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i];
      return out;
    }
    MFADd operator +(const double &v) const{
      size_t s = maxN;
      MFADd out(this->value+(double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i];
      return out;
    }
    MFADd operator +(const long double &v) const{
      size_t s = maxN;
      MFADd out(this->value+(long double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i];
      return out;
    }
    
    MFADd& operator +=(const MFADd &u) {
      this->value+=u.value;
      for (size_t i=0; i<maxN; i++) this->grad[i]+=u.grad[i];
      return *this;
    }
    MFADd& operator +=(const int &u) {
      this->value+=u;
      return *this;
    }
    MFADd& operator +=(const float &u) {
      this->value+=u;
      return *this;
    }
    MFADd& operator +=(const double &u) {
      this->value+=u;
      return *this;
    }
    MFADd& operator +=(const long double &u) {
      this->value+=u;
      return *this;
    }
    
    MFADd operator -(const MFADd &v) const{ 
      size_t s = maxN ;
      MFADd out(this->value-v.value,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i]-v.grad[i];
      return out;
    }
    MFADd operator -() const {
      size_t s = maxN;
      MFADd out(-this->value,s);
      for (size_t i=0; i<s; i++) out.grad[i] = -this->grad[i];
      return out;
    }
    MFADd operator -(const int &v) const{ 
      size_t s = maxN;
      MFADd out(this->value-(double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i];
      return out;
    }
    MFADd operator -(const float &v) const{ 
      size_t s = maxN;
      MFADd out(this->value-(double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i];
      return out;
    }
    MFADd operator -(const double &v) const{ 
      size_t s = maxN;
      MFADd out(this->value-(double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i];
      return out;
    }
    MFADd operator -(const long double &v) const{ 
      size_t s = maxN;
      MFADd out(this->value-(long double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i];
      return out;
    }
    
    MFADd& operator -=(const MFADd &u) {
      this->value -= u.value;
      for (size_t i=0; i<maxN; i++) this->grad[i]-=u.grad[i];
      return *this;
    }
    MFADd& operator -=(const int &u) {
      this->value -= u;
      return *this;
    }
    MFADd& operator -=(const float &u) {
      this->value -= u;
      return *this;
    }
    MFADd& operator -=(const double &u) {
      this->value -= u;
      return *this;
    }
    MFADd& operator -=(const long double &u) {
      this->value -= u;
      return *this;
    }
    
    MFADd operator *(const MFADd &v) const {
      size_t s = maxN ; 
      MFADd out(this->value*v.value,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i]*v.value + v.grad[i]*this->value;
      return out;
    }
    MFADd operator *(const int &v) const {
      size_t s = maxN;
      MFADd out(this->value*(double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i] * (double)v;
      return out;
    }
    MFADd operator *(const float &v) const {
      size_t s = maxN;
      MFADd out(this->value*(double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i] * (double)v;
      return out;
    }
    MFADd operator *(const double &v) const {
      size_t s = maxN;
      MFADd out(this->value*(double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i] * (double)v;
      return out;
    }
    MFADd operator *(const long double &v) const {
      size_t s = maxN;
      MFADd out(this->value*(long double)v,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i] * (long double)v;
      return out;
    }
     
    MFADd& operator *=(const MFADd &u) {
      size_t s = maxN;
      for (size_t i=0; i<s; i++) this->grad[i] = this->grad[i]*u.value + u.grad[i] * this->value;
      this->value*=u.value;
      return *this;
    }
    MFADd& operator *=(const int &u) {
      size_t s = maxN;
      for (size_t i=0; i<s; i++) this->grad[i] *= u;
      this->value *=u ;
      return *this;
    }
    MFADd& operator *=(const float &u) {
      size_t s = maxN;
      for (size_t i=0; i<s; i++) this->grad[i] *= u;
      this->value *=u ;
      return *this;
    }
    MFADd& operator *=(const double &u) {
      size_t s = maxN;
      for (size_t i=0; i<s; i++) this->grad[i] *= u;
      this->value *=u ;
      return *this;
    }
    MFADd& operator *=(const long double &u) {
      size_t s = maxN;
      for (size_t i=0; i<s; i++) this->grad[i] *= u;
      this->value *=u ;
      return *this;
    }
    
    MFADd operator /(const MFADd &v) const{ 
      size_t s = maxN ;
      double div = this->value/v.value;
      MFADd out(div,s);
      for (size_t i=0; i<s; i++) out.grad[i] = (this->grad[i] - div*v.grad[i])/v.value;
      return out;
    }
    MFADd operator /(const int &v) const{ 
      size_t s = maxN;
      double div = this->value/(double)v;
      MFADd out(div,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i]/(double)v;
      return out;
    }
    MFADd operator /(const float &v) const{ 
      size_t s = maxN;
      double div = this->value/(double)v;
      MFADd out(div,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i]/(double)v;
      return out;
    }
    MFADd operator /(const double &v) const{ 
      size_t s = maxN;
      double div = this->value/(double)v;
      MFADd out(div,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i]/(double)v;
      return out;
    }
    MFADd operator /(const long double &v) const{ 
      size_t s = maxN;
      double div = this->value/(long double)v;
      MFADd out(div,s);
      for (size_t i=0; i<s; i++) out.grad[i] = this->grad[i]/(long double)v;
      return out;
    }

    MFADd& operator /=(const MFADd &u) {
      size_t s = maxN;
      double div = this->value/u.value;
      for (size_t i=0; i<s; i++) this->grad[i] = (this->grad[i] - div*u.grad[i])/u.value;
      this->value = div ;
      return *this;
    }
    MFADd& operator /=(const int &u) {
      size_t s = maxN;
      for (size_t i=0; i<s; i++) this->grad[i] /= u;
      this->value /= u ;
      return *this;
    }
    MFADd& operator /=(const float &u) {
      size_t s = maxN;
      for (size_t i=0; i<s; i++) this->grad[i] /= u;
      this->value /= u ;
      return *this;
    }
    MFADd& operator /=(const double &u) {
      size_t s = maxN;
      for (size_t i=0; i<s; i++) this->grad[i] /= u;
      this->value /= u ;
      return *this;
    }
    MFADd& operator /=(const long double &u) {
      size_t s = maxN;
      for (size_t i=0; i<s; i++) this->grad[i] /= u;
      this->value /= u ;
      return *this;
    }
    
    bool operator ==(const MFADd &u)const  {
      return ((this->value==u.value)?true:false);
    }
    bool operator !=(const MFADd &u) const {
      return ((this->value!=u.value)?true:false);
    }
    bool operator >(const MFADd &u) const {
      return ((this->value>u.value)?true:false);
    }
    bool operator <(const MFADd &u) const {
      return ((this->value<u.value)?true:false);
    }
    bool operator >=(const MFADd &u) const {
      return ((this->value>=u.value)?true:false);
    }
    bool operator <=(const MFADd &u) const {
      return ((this->value<=u.value)?true:false);
    }
    bool operator ==(const int &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const int &u) const {
      return ((value!=u)?true:false);
    }
    bool operator >(const int &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const int &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const int &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const int &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const float &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const float &u) const {
      return ((value!=u)?true:false);
    }
    bool operator >(const float &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const float &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const float &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const float &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const double &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const double &u) const {
      return ((value!=u)?true:false);
    }
    bool operator  >(const double &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const double &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const double &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const double &u) const {
      return ((value<=u)?true:false);
    }
    bool operator ==(const long double &u) const {
      return ((value==u)?true:false);
    }
    bool operator !=(const long double &u) const {
      return ((value!=u)?true:false);
    }
    bool operator  >(const long double &u) const {
      return ((value>u)?true:false);
    }
    bool operator <(const long double &u) const {
      return ((value<u)?true:false);
    }
    bool operator >=(const long double &u) const {
      return ((value>=u)?true:false);
    }
    bool operator <=(const long double &u) const {
      return ((value<=u)?true:false);
    }
  } ;
  
  inline std::ostream &operator <<(std::ostream& stream, const MFADd &u) {
    stream << u.value ;
    for (size_t i=0; i<MFADd::maxN; i++ ) stream << " " << u.grad[i];
    return stream;
  }
  inline std::istream &operator >> (std::istream& stream, MFADd &u) {
    u.value = 0. ;
    for (size_t i=0;i<MFADd::maxN;i++) u.grad[i] = 0. ;
    stream >> u.value;
    if(stream.peek() == '^') {
      stream.get() ;
      stream >> u.grad[0] ;
      for (size_t i=1;i<MFADd::maxN;i++) {
        if(stream.peek()=='^') {
	  while(stream.peek()=='^')
	    stream.get() ;
 	  stream >> u.grad[i] ;
        }
      }
    }
    return stream;
  }
  
  using std::sqrt;
  using std::sin;
  using std::cos;
  using std::tan;
  using std::exp;
  using std::log;
  using std::fabs;
  //using std::cbrt;
  using std::sinh;
  using std::cosh;
  using std::tanh;
  using std::asin;
  using std::acos;
  using std::atan;
  //using std::asinh;
  //using std::acosh;
  //using std::atanh;
  using std::floor;
  using std::ceil;
  using std::pow;
  using std::max;
  using std::min;

#ifdef NO_CMATH
  inline double ceil(const MFADd u) {
    return ::ceil(u.value) ;
  }
  inline double floor(const MFADd u) {
    return ::floor(u.value) ;
  }

  inline MFADd erf(const MFADd u) {
    double ef = ::erf(u.value) ;
    size_t s = MFADd::maxN ;
    MFADd out(ef,s) ;
    double efp = (2./sqrt(M_PI))*std::exp(-u.value*u.value) ;
    for(size_t i=0;i<s;++i)
      out.grad[i] = u.grad[i]*efp ; 
    return out ;
  }

  inline MFADd sin(const MFADd &u) {
    size_t s = MFADd::maxN;
    MFADd out(::sin(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*::cos(u.value);
    return out;
  }
  inline MFADd cos(const MFADd &u) {
    size_t s = MFADd::maxN;
    MFADd out(::cos(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = -u.grad[i]*::sin(u.value);
    return out;
  }
  inline MFADd tan(const MFADd &u) {
    size_t s = MFADd::maxN;
    MFADd out(::tan(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]/pow(::cos(u.value),2) ;
    return out;
  }
  inline MFADd exp(const MFADd &u) {
    size_t s = MFADd::maxN;
    MFADd out(::exp(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*::exp(u.value);
    return out;
  }
  inline MFADd log(const MFADd &u) {
    size_t s = MFADd::maxN;
    MFADd out(::log(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]/u.value ;
    return out;
  }
  inline MFADd log10(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::log(u.value)/::log(10.0),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]/u.value/::log(10.0);
    return out;
  }

  inline MFADd abs(MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::fabs(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*( (u.value<0.0)?-1.0:1.0 ); 
    return out;
  }
  inline MFADd fabs(MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::fabs(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*( (u.value<0.0)?-1.0:1.0 );
    return out;
  }
  
  inline MFADd pow(const MFADd u, const int k) {
    size_t s = MFADd::maxN;
    MFADd out(::pow(u.value, k),s);
    for (size_t i=0; i<s; i++) out.grad[i] = (double)k * ::pow(u.value, (double)k-1.0)*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const MFADd u, const float k) {
    size_t s = MFADd::maxN;
    MFADd out(::pow(u.value, k),s);
    for (size_t i=0; i<s; i++) out.grad[i] = (double)k * ::pow(u.value, (double)k-1.0)*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const MFADd u, const double k) {
    size_t s = MFADd::maxN;
    MFADd out(::pow(u.value, k),s);
    for (size_t i=0; i<s; i++) out.grad[i] = k * ::pow(u.value, k-1.0)*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const MFADd u, const long double k) {
    size_t s = MFADd::maxN;
    MFADd out(::pow(u.value, k),s);
    for (size_t i=0; i<s; i++) out.grad[i] = (double)k * ::pow(u.value, (double)k-1.0)*u.grad[i] ;
    return out ;
  }
  inline MFADd sqrt(const MFADd &u) { 
    size_t s = MFADd::maxN;
    double sq = ::sqrt(u.value);
    MFADd out(sq,s) ;
    for (size_t i=0; i<s; i++) out.grad[i] = 0.5*u.grad[i]/::max(sq,1e-30) ;
    return out;
  }
  inline MFADd pow(const MFADd k, const MFADd u) {
    size_t s = MFADd::maxN;
    double kpu = ::pow(k.value,u.value) ;
    MFADd out(kpu,s);
    for (size_t i=0; i<s; i++) out.grad[i] = ::pow(k.value,u.value-1.)*k.grad[i]*u.value + kpu*::log(k.value)*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const int k, const MFADd u) {
    size_t s = MFADd::maxN;
    double kpu = ::pow(k,u.value) ;
    MFADd out(kpu,s);
    for (size_t i=0; i<s; i++) out.grad[i] = kpu*::log(double(k))*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const float k, const MFADd u) {
    size_t s = MFADd::maxN;
    double kpu = ::pow(k,u.value) ;
    MFADd out(kpu,s);
    for (size_t i=0; i<s; i++) out.grad[i] = kpu*::log(double(k))*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const double k, const MFADd u) {
    size_t s = MFADd::maxN;
    double kpu = ::pow(k,u.value) ;
    MFADd out(kpu,s);
    for (size_t i=0; i<s; i++) out.grad[i] = kpu*::log(k)*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const long double k, const MFADd u) {
    size_t s = MFADd::maxN;
    double kpu = ::pow(k,u.value) ;
    MFADd out(kpu,s);
    for (size_t i=0; i<s; i++) out.grad[i] = kpu*::log(k)*u.grad[i] ;
    return out ;
  }

  inline MFADd sinh(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::sinh(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = 0.5*u.grad[i]*(::exp(u.value) + ::exp(-(u.value))) ;
    return out;
  }
  inline MFADd cosh(const MFADd u) {
    size_t s = FMADd::maxN;
    MFADd out(::cosh(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = 0.5*u.grad[i]*(::exp(u.value) - ::exp(-(u.value))) ;
    return out;
  }
  inline MFADd tanh(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::tanh(u.value),s);
    double ex = ::exp(::min(u.value,350.0)) ;
    double exm = ::exp(::min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*(1.-dex*dex/(sex*sex)) ;
    return out;
  }
  inline MFADd asin(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::asin(u.value),s);
    double val = ::sqrt(1.0-u.value*u.value);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]/val;
    return out;
  }
  inline MFADd acos(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::acos(u.value),s);
    double val = ::sqrt(1.0-u.value*u.value);
    for (size_t i=0; i<s; i++) out.grad[i] = -u.grad[i]/val;
    return out;
  }
  inline MFADd atan(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::atan(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]/(u.value*u.value + 1.0);
    return out;
  }
  inline MFADd atan2(const MFADd u, const MFADd v) {  
    return atan(u/v);
  }
  inline MFADd asinh(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::asinh(u.value),s);
    double strm = ::sqrt(u.value*u.value+1.) ;
    double uvstrm = u.value+strm ;
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*(u.value/strm+1.)/uvstrm ;
    return out;
  }
  inline MFADd acosh(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::acosh(u.value),s);
    double strm = ::sqrt(u.value*u.value-1.) ;
    double uvstrm = u.value+strm ;
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*(u.value/strm +1.)/uvstrm ;
    return out;
  }
  inline MFADd atanh(const MFADd &u) {
    size_t s = MFADd::maxN;
    const double uv = u.value ;
    MFADd out(::atanh(uv),s);
    for (size_t i=0; i<s; i++) out.grad[i] = .5*u.grad[i]*(1./(1.+uv)+1./(1.-uv)) ;
    return out;
  }
#else
  inline double ceil(const MFADd u) {
    return std::ceil(u.value) ;
  }
  inline double floor(const MFADd u) {
    return std::floor(u.value) ;
  }
  inline MFADd erf(const MFADd u) {
    double ef = ::erf(u.value) ;
    size_t s = MFADd::maxN ;
    MFADd out(ef,s) ;
    double efp = (2./sqrt(M_PI))*std::exp(-u.value*u.value) ;
    for(size_t i=0;i<s;++i)
      out.grad[i] = u.grad[i]*efp ; 
    return out ;
  }
  inline MFADd sin(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::sin(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*std::cos(u.value);
    return out;
  }
  inline MFADd cos(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::cos(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = -u.grad[i]*std::sin(u.value);
    return out;
  }
  inline MFADd tan(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::tan(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]/pow(std::cos(u.value),2) ;
    return out;
  }

  inline MFADd exp(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::exp(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*std::exp(u.value);
    return out;
  }
  inline MFADd log(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::log(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]/u.value ;
    return out;
  }
  inline MFADd log10(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::log(u.value)/std::log(10.0),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]/u.value/std::log(10.0);
    return out;
  }
  inline MFADd abs(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::fabs(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*( (u.value<0.0)?-1.0:1.0 ); 
    return out;
  }
  inline MFADd fabs(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::fabs(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*( (u.value<0.0)?-1.0:1.0 );
    return out;
  }
  
  inline MFADd pow(const MFADd u, const int k) {
    size_t s = MFADd::maxN;
    MFADd out(std::pow(u.value, k),s);
    for (size_t i=0; i<s; i++) out.grad[i] = (double)k * std::pow(u.value, (double)k-1.0)*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const MFADd u, const float k) {
    size_t s = MFADd::maxN;
    MFADd out(std::pow(u.value, double(k)),s);
    for (size_t i=0; i<s; i++) out.grad[i] = (double)k * std::pow(u.value, (double)k-1.0)*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const MFADd u, const double k) {
    size_t s = MFADd::maxN;
    MFADd out(std::pow(u.value, k),s);
    for (size_t i=0; i<s; i++) out.grad[i] = k * std::pow(u.value, k-1.0)*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const MFADd u, const long double k) {
    size_t s = MFADd::maxN;
    MFADd out(std::pow(u.value, double(k)),s);
    for (size_t i=0; i<s; i++) out.grad[i] = (double)k * std::pow(u.value, (double)k-1.0)*u.grad[i] ;
    return out ;
  }
  inline MFADd sqrt(const MFADd u) { 
    size_t s = MFADd::maxN;
    double sq = std::sqrt(u.value);
    MFADd out(sq,s) ;
    for (size_t i=0; i<s; i++) out.grad[i] = 0.5*u.grad[i]/::max(sq,1e-30) ;
    return out;
  }
  inline MFADd pow(const MFADd k, const MFADd u) {
    size_t s = MFADd::maxN;
    double kpu = std::pow(k.value,u.value) ;
    MFADd out(kpu,s);
    for (size_t i=0; i<s; i++) out.grad[i] = std::pow(k.value,u.value-1.)*k.grad[i]*u.value + kpu*std::log(k.value)*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const int k, const MFADd u) {
    size_t s = MFADd::maxN;
    double kpu = std::pow(k,u.value) ;
    MFADd out(kpu,s);
    for (size_t i=0; i<s; i++) out.grad[i] = kpu*std::log(double(k))*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const float k, const MFADd u) {
    size_t s = MFADd::maxN;
    double kpu = std::pow(k,u.value) ;
    MFADd out(kpu,s);
    for (size_t i=0; i<s; i++) out.grad[i] = kpu*std::log(double(k))*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const double k, const MFADd u) {
    size_t s = MFADd::maxN;
    double kpu = std::pow(k,u.value) ;
    MFADd out(kpu,s);
    for (size_t i=0; i<s; i++) out.grad[i] = kpu*std::log(k)*u.grad[i] ;
    return out ;
  }
  inline MFADd pow(const long double k, const MFADd u) {
    size_t s = MFADd::maxN;
    double kpu = std::pow(double(k),u.value) ;
    MFADd out(kpu,s);
    for (size_t i=0; i<s; i++) out.grad[i] = kpu*std::log(k)*u.grad[i] ;
    return out ;
  }
  
  inline MFADd max (const MFADd &a, const MFADd& b) { 
    size_t s = MFADd::maxN;
    MFADd out((a.value<b.value)?b.value:a.value);
    for (size_t i=0; i<s; i++) out.grad[i] = (a.value<b.value)?b.grad[i]:a.grad[i];
    return out;
  }
  inline MFADd max (const MFADd &a, const int& b) {
    return max(a,MFADd((double)b));
  }
  inline MFADd max (const MFADd &a, const float& b) {
    return max(a,MFADd((double)b));
  }
  inline MFADd max (const MFADd &a, const double & b) {
    return max(a,MFADd((double)b));
  }
  inline MFADd max (const MFADd &a, const long double & b) {
    return max(a,MFADd((double)b));
  }
  inline MFADd min (const MFADd &a, const MFADd& b) { 
    size_t s = MFADd::maxN;
    MFADd out((a.value<b.value)?a.value:b.value);
    for (size_t i=0; i<s; i++) out.grad[i] = (a.value<b.value)?a.grad[i]:b.grad[i];
    return out;
  }
  inline MFADd min (const MFADd &a, const int& b) {
    return min(a, MFADd((double)b));
  }
  inline MFADd min (const MFADd &a, const float& b) {
    return min(a, MFADd((double)b));
  }
  inline MFADd min (const MFADd &a, const double& b) {
    return min(a, MFADd((double)b));
  }
  inline MFADd min (const MFADd &a, const long double& b) {
    return min(a, MFADd((double)b));
  }
  
  inline MFADd sinh(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::sinh(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = 0.5*u.grad[i]*(std::exp(u.value) + std::exp(-(u.value))) ;
    return out;
  }
  inline MFADd cosh(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::cosh(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = 0.5*u.grad[i]*(std::exp(u.value) - std::exp(-(u.value))) ;
    return out;
  }
  inline MFADd tanh(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::tanh(u.value),s);
    double ex = std::exp(::min(u.value,350.0)) ;
    double exm = std::exp(::min(-u.value,350.0)) ;
    double dex = ex-exm ;
    double sex = ex+exm ;
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*(1.-dex*dex/(sex*sex)) ;
    return out;
  }

  inline MFADd asin(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::asin(u.value),s);
    double val = std::sqrt(1.0-u.value*u.value);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]/val;
    return out;
  }
  inline MFADd acos(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::acos(u.value),s);
    double val = std::sqrt(1.0-u.value*u.value);
    for (size_t i=0; i<s; i++) out.grad[i] = -u.grad[i]/val;
    return out;
  }
  inline MFADd atan(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(std::atan(u.value),s);
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]/(u.value*u.value + 1.0);
    return out;
  }
  inline MFADd atan2(const MFADd u, const MFADd v) { return atan(u/v); }

  inline MFADd asinh(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::asinh(u.value),s);
    double strm = std::sqrt(u.value*u.value+1.) ;
    double uvstrm = u.value+strm ;
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*(u.value/strm+1.)/uvstrm ;
    return out;
  }
  inline MFADd acosh(const MFADd u) {
    size_t s = MFADd::maxN;
    MFADd out(::acosh(u.value),s);
    double strm = std::sqrt(u.value*u.value-1.) ;
    double uvstrm = u.value+strm ;
    for (size_t i=0; i<s; i++) out.grad[i] = u.grad[i]*(u.value/strm +1.)/uvstrm ;
    return out;
  }
  inline MFADd atanh(const MFADd u) {
    size_t s = MFADd::maxN;
    const double uv = u.value ;
    MFADd out(::atanh(uv),s);
    for (size_t i=0; i<s; i++) out.grad[i] = .5*u.grad[i]*(1./(1.+uv)+1./(1.-uv)) ;
    return out;
  }
#endif
  inline MFADd operator +(const double u,const MFADd v) { 
    size_t s = MFADd::maxN;
    MFADd out(u+v.value,s);
    for (size_t i=0; i<s; i++) out.grad[i] = v.grad[i];
    return out;
  }
  inline MFADd operator -(const double u,const MFADd v) { 
    size_t s = MFADd::maxN;
    MFADd out(u-v.value,s);
    for (size_t i=0; i<s; i++) out.grad[i] = v.grad[i];
    return out;
  }
  inline MFADd operator *(const double u,const MFADd v) { 
    size_t s = MFADd::maxN;
    MFADd out(u*v.value,s);
    for (size_t i=0; i<s; i++) out.grad[i] = u * v.grad[i];
    return out;
  }
  inline MFADd operator /(const double &u,const MFADd &v) { 
    size_t s = MFADd::maxN;
    MFADd out(u/v.value,s);
    for (size_t i=0; i<s; i++) out.grad[i] = -u*v.grad[i]/v.value/v.value;
    return out;
  }
  //----------------helpers
  inline float realToFloat(double v) { return float(v) ; }
  inline double realToDouble(double v) { return v ; } 

  inline float realToFloat(const FADd &v) { return float(v.value) ; }
  inline double realToDouble(const FADd &v) { return v.value ; }

  inline float realToFloat(const FAD2d &v) { return float(v.value) ; }
  inline double realToDouble(const FAD2d &v) { return v.value ; }

  inline float realToFloat(const MFADd &v) { return float(v.value) ; }
  inline double realToDouble(const MFADd &v) { return v.value ; }

}

#endif
