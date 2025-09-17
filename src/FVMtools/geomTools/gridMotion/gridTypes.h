//#############################################################################
//#
//# Copyright 2014-2025, Mississippi State University
//#
//# The GridMover is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The GridMover software is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the GridMover software. If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef GRIDMOTIONTYPES_H
#define GRIDMOTIONTYPES_H
#ifdef USE_LOCI_TYPES
#include <Loci>
#include <Tools/tools.h>
namespace gridMotion {
  using Loci::vector3d ;
  using Loci::tensor3d ;
  using Loci::dot ;
  using Loci::cross ;
  using Loci::norm ;
  using Loci::Array ;
  using Loci::tmp_array ;
  using Loci::Scalar ;
  using Loci::Vect ;
  using Loci::const_Vect ;
  using Loci::dotprod_accum ;
  using Loci::solve_lu ;
  using Loci::pivot_type ;
  using Loci::solve_lu_pivot ;
  using Loci::const_Mat ;
  using Loci::Mat ;
  using Loci::stopWatch ;
}
#else
#include <iostream>

#ifdef restrict
#undef restrict
#endif

#if defined(__GNUC__)

#if defined(__INTEL_COMPILER)
// Intel Compiler

#define HAVE_IVDEP

#define restrict __restrict
// else INTEL compiler
#else

#define restrict __restrict__
// else on intel compiler
#endif
#else
// not gcc compiler
#if defined(__IBMCPP__)
/* IBM C++ Compiler */
#define restrict __restrict
#else
// not ibm compiler or gcc compiler
#if defined(__sgi)
// SGI Compiler
#define NO_CSTDLIB
#define NO_CMATH
#define NO_SIGNBIT
#define restrict __restrict
#else
// not ibm, gcc, or sgi compiler
#if defined(__PGI)
/* Portland group compiler*/
#define restrict __restrict
#define HAVE_IVDEP
#define NO_SIGNBIT
#else
// not gcc, ibm, sgi, or pgi compiler
#if defined(__SUNPRO_CC)
/* Sun CC compiler */
#define restrict 
#else
// not gcc or ibm or sgi or pgi or sunpro compiler
#define restrict

#endif
// not gcc or ibm or sgi or pgi compiler
#endif

// not gcc or ibm or sgi compiler
#endif
// not gcc or ibm  compiler
#endif
// not gcc compiler
#endif



#ifdef NO_CMATH
#include <math.h>
#else
#include <cmath>
// trigonemetric
using std::acos ;
using std::asin ;
using std::atan ;
using std::atan2 ;
using std::cos ;
using std::sin ;
using std::tan ;
// hyperbolic
using std::sinh ;
using std::cosh ;
using std::tanh ;
// exponetial and logrithmic
using std::exp ;
using std::log ;
using std::log10 ;
using std::sqrt ;
using std::pow ;
using std::fabs ;
// misc
using std::ceil ;
using std::floor ;
using std::fmod ;
#endif

#include <mpi.h>

namespace gridMotion {

  //---------------------Array----------------------//
  template <class T,size_t n> class Array {
    T x[n] ;
  public:
    typedef T * iterator ;
    typedef const T * const_iterator ;
    
    //    Array() {} ;
    //    Array(const Array<T,n> &v)
    //    { for(size_t i=0;i<n;++i) x[i] = v.x[i] ; } 
    //    Array<T,n> &operator=(const Array<T,n> &v)
    //    { for(size_t i=0;i<n;++i) x[i] = v.x[i] ; return *this ; } 

    Array<T,n> &operator +=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] += v.x[i] ; return *this ; }
    Array<T,n> &operator -=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] -= v.x[i] ; return *this ; }
    Array<T,n> &operator *=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] *= v.x[i] ; return *this ; }
    Array<T,n> &operator /=(const Array<T,n> &v)
    { for(size_t i=0;i<n;++i) x[i] /= v.x[i] ; return *this ; }

    T &operator[](size_t indx) { return x[indx]; }
    const T &operator[](size_t indx) const { return x[indx] ; }

    iterator begin() { return &x[0] ; }
    iterator end() { return begin()+n ; }
    const_iterator begin() const { return &x[0] ; }
    const_iterator end() const { return begin()+n ; }

    size_t size() const  { return n ; }
  } ;

  //---------------------vector3d------------------//
  template <class T> 
    struct vector3d {
      T x,y,z ;
      vector3d() {} 
      vector3d(T xx,T yy, T zz) : x(xx),y(yy),z(zz) {}
      vector3d(const vector3d &v) {x=v.x;y=v.y;z=v.z;}
      template <class S> vector3d(const vector3d<S> &v) {x=v.x;y=v.y;z=v.z;}
      template <class S> vector3d(const Array<S,3> &a) {x=a[0];y=a[1];z=a[2];}
      template <class S> operator Array<S,3>() {
	Array<S,3> a ;
	a[0] = x ;
	a[1] = y ;
	a[2] = z ;
	return a;
      }
    } ;
  
  template <class T> inline std::ostream & operator<<(std::ostream &s, const vector3d<T> &v)
    {
      s << v.x << ' ' << v.y << ' ' << v.z << ' ' ;
      return s ;
    }

  template <class T> inline std::istream &operator>>(std::istream &s, vector3d<T> &v)
    {
      s >> v.x >> v.y >> v.z ;
      return s ;
    }

  template <class T> inline T dot(const vector3d<T> &v1, const vector3d<T> &v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z ;
  }

  template <class T> inline T norm(const vector3d<T> &v) {
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z) ;
  }

  template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.y*v2.z-v1.z*v2.y,
                       v1.z*v2.x-v1.x*v2.z,
                       v1.x*v2.y-v1.y*v2.x) ;
  }

  template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const T ra2[]) {
    return vector3d<T>(v1.y*ra2[2]-v1.z*ra2[1],
                       v1.z*ra2[0]-v1.x*ra2[2],
                       v1.x*ra2[1]-v1.y*ra2[0]) ;
  }
  template<class T> inline vector3d<T> cross(const T ra1[], const vector3d<T> &v2) {
    return vector3d<T>(ra1[1]*v2.z-ra1[2]*v2.y,
                       ra1[2]*v2.x-ra1[0]*v2.z,
                       ra1[0]*v2.y-ra1[1]*v2.x) ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, float val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, float val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, double val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, double val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, long double val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
  }

  template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, long double val) {
    target.x /= val ;
    target.y /= val ;
    target.z /= val ;
    return target ;
  }

  template<class T> inline vector3d<T> operator+=(vector3d<T> &target, const vector3d<T> &val) {
    target.x += val.x ;
    target.y += val.y ;
    target.z += val.z ;
    return target ;
  }

  template<class T> inline vector3d<T> operator-=(vector3d<T> &target, const vector3d<T> &val) {
    target.x -= val.x ;
    target.y -= val.y ;
    target.z -= val.z ;
    return target ;
  }

  template<class T> inline vector3d<T> operator+(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z) ;
  }

  template<class T> inline vector3d<T> operator-(const vector3d<T> &v1, const vector3d<T> &v2) {
    return vector3d<T>(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, float r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(float r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, float r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, double r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(double r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, double r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, long double r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

  template<class T> inline vector3d<T> operator*(long double r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
  }

  template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, long double r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
  }

  template<class T>  struct tensor3d : public vector3d<vector3d< T > > {
    tensor3d() {}
    tensor3d(vector3d<T> xx,vector3d<T> yy, vector3d<T> zz)
      : vector3d<vector3d< T> > (xx,yy,zz) {}
    tensor3d(const tensor3d &v) : vector3d<vector3d< T> >(v) {}
  } ;

  template<class T> inline vector3d<T> dot(const tensor3d<T> &t,
                                           const vector3d<T> &v) {
    return vector3d<T>(dot(t.x,v),dot(t.y,v),dot(t.z,v)) ;
  }

  template<class T> inline tensor3d<T> product(const tensor3d<T> &t1,
                                               const tensor3d<T> &t2) {
    tensor3d<T> temp ;
    temp.x.x = t1.x.x*t2.x.x+t1.x.y*t2.y.x+t1.x.z*t2.z.x ;
    temp.y.x = t1.y.x*t2.x.x+t1.y.y*t2.y.x+t1.y.z*t2.z.x ;
    temp.z.x = t1.z.x*t2.x.x+t1.z.y*t2.y.x+t1.z.z*t2.z.x ;

    temp.x.y = t1.x.x*t2.x.y+t1.x.y*t2.y.y+t1.x.z*t2.z.y ;
    temp.y.y = t1.y.x*t2.x.y+t1.y.y*t2.y.y+t1.y.z*t2.z.y ;
    temp.z.y = t1.z.x*t2.x.y+t1.z.y*t2.y.y+t1.z.z*t2.z.y ;

    temp.x.z = t1.x.x*t2.x.z+t1.x.y*t2.y.z+t1.x.z*t2.z.z ;
    temp.y.z = t1.y.x*t2.x.z+t1.y.y*t2.y.z+t1.y.z*t2.z.z ;
    temp.z.z = t1.z.x*t2.x.z+t1.z.y*t2.y.z+t1.z.z*t2.z.z ;

    return temp ;
  }

 template <class T,size_t n> inline std::ostream &
    operator<<(std::ostream &s, const Array<T,n> &v) {
    for(size_t i=0;i<n;++i)
      s << v[i] << ' ' ;
    return s ;
  }

  template <class T,size_t n> inline std::istream &
    operator>>(std::istream &s, Array<T,n> &v) {
    for(int i=0;i<n;++i)
      s >> v[i] ;
    return s ;
  }

  
  // For allocating temporary arrays of small size
  const int tmp_array_internal_SIZE=25 ;
  template <class T> class tmp_array {
    int sz ;
    T data[tmp_array_internal_SIZE] ;
    T * p ;
    void alloc(int size) {
      sz = size ;
      p = data ;
      if(sz > tmp_array_internal_SIZE)
        p = new T[sz] ;
    }
    void free() {
      if(sz > tmp_array_internal_SIZE)
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


  class stopWatch {
#ifdef USE_PAPI
    long_long start_time ;
#else
    double start_time ;
#endif
  public:
    void start() { // This method resets the clock
#ifdef USE_PAPI
      start_time = PAPI_get_real_usec();
#else
      start_time = MPI_Wtime() ;
#endif
    }
    double stop() { // This method returns time since last start call
#ifdef USE_PAPI
      return 1e-6*double(PAPI_get_real_usec()-start_time) ;
#else
      return MPI_Wtime()-start_time ;
#endif
    }
  } ;
  


  //*******************************************************************/
  template <class T> struct Scalar {
    T val ;
    Scalar(T v) : val(v) { }
  } ;
  template <class T> Scalar<T> mk_Scalar(T v) { return Scalar<T>(v) ;} 
  //*******************************************************************/
  template <class T> class Vect ;
  template <class T> class const_Vect {
  public:
    friend class Vect<T> ;
    const T * restrict ptr ;
#ifdef BOUNDS_CHECK
    int size ;
#endif
  public:
    const_Vect(const T *p ,int sz) {
      ptr = p ;
#ifdef BOUNDS_CHECK
      size = sz ;
#endif
    }
    const T &restrict operator[](int idx) const restrict {
#ifdef BOUNDS_CHECK
      fatal(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    operator const T *restrict () const restrict {
      return ptr ;
    }
  } ;
  //******************************************************************/
  template <class T> class Vect {
  public:
    T *ptr ;
    int size ;
  public:
    Vect() {};
    void setSize( int s ) {
      size = s;
    }
    int getSize() { return size; }
    
    Vect(T *p ,int sz) {
      ptr = p ;
      size = sz ;
    }

    void operator=(const Vect<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ = *p2++ ;
    }


    void operator=(const const_Vect<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ = *p2++ ;
    }


    template <class S> void operator=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ = s.val ;
    }

    template <class S> void operator+=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ += s.val ;
    }
      
    template <class S> void operator*=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ *= s.val ;
    }
      
    template <class S> void operator-=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ -= s.val ;
    }
      
    template <class S> void operator/=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ /= s.val ;
    }

    template <class S> void operator+=(const Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ += *p2++ ;
    }

    template <class S> void operator+=(const const_Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ += *p2++ ;
    }

    template <class S> void operator-=(const Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ -= *p2++ ;
    }

    template <class S> void operator-=(const const_Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ -= *p2++ ;
    }

    T &operator[](int idx) {
#ifdef BOUNDS_CHECK
      fatal(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    const T &operator[](int idx) const {
#ifdef BOUNDS_CHECK
      fatal(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    operator T*() {
      return ptr ;
    }

    operator const T *() const {
      return ptr ;
    }
  } ;

  typedef int pivot_type ;

  template< class T1, class T2, class T3 >
  inline void dotprod_accum(const T1 * const restrict A,
                const T2 * const restrict vin,
                T3 * restrict vo,
                int size) {
    for(int j=0;j<size;++j) {
      const T2 in = vin[j] ;
      const T1 *restrict Aj = A+size*j ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int i=0;i<size;++i)
          vo[i] += (Aj[i])*in ;
      }
    }

  template<class T1, class T2, class T3> 
  inline void solve_lu(const T1 * restrict A,
                       const T2 *restrict b,
                       T3 *restrict x,
                       int size) {
    // Perform forward solve Ly = b, note b becomes y after this step
    for(int i=0;i<size;++i) {
      T3 tmp = b[i] ;
      const T1 * restrict Aij = A+i ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int j=0;j<i;++j,Aij+=size)
        tmp -= *Aij*x[j] ;
      x[i] = tmp ;


    }
    // Do back solver Ux = y
    const T1 *restrict Ai = A + size*(size-1) ;
    for(int i=size-1;i>=0;--i,Ai-=size) {
      const T1 *restrict Aj = Ai + size ;
      T3 tmp = x[i] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int j=i+1;j<size;++j,Aj+=size)
        tmp -= Aj[i]*x[j] ;
      x[i] = tmp/Ai[i] ;
    }
  }


  template<class T1, class T2, class T3> 
  inline void solve_lu_pivot(const T1 * restrict A,
                             const T2 * restrict b,
                             T3 * restrict x,
                             const pivot_type * restrict pivot,
                             int size) {
      // Perform forward solve Ly = b, note b becomes x after this step
      for(int i=0;i<size;++i) {
        T3 xi = b[pivot[i]] ;
        const T1 * restrict Aj = A ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=0;j<i;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi ;
      }
      // Do back solver Ux = y
      const T1 * restrict Ai = A + size*(size-1) ;
      for(int i=size-1;i>=0;--i,Ai-=size) {
        const T1 *restrict Aj = Ai + size ;
        T3 xi = x[i] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int j=i+1;j<size;++j,Aj+=size)
          xi -= Aj[i]*x[j] ;
        x[i] = xi/Ai[i] ;
      }
    }
  
  template <class T> class Mat ;

  //**************************************************************************
  
  template <class T> class const_Mat_partial {
    const T * restrict ptr ;
    int size ;
    int i ;
  public:
    const_Mat_partial(const T *p,int sz, int ii) :ptr(p),size(sz),i(ii) {}
    const T & restrict operator[](int j) {
#ifdef BOUNDS_CHECK
      warn(j>=size || j < 0) ;
#endif
      return ptr[j*size + i] ;
    }
  } ;

  //**************************************************************************

  
  template <class T> class const_Mat {
  public:
    const T *restrict ptr ;
    int size ;
    friend class Mat<T> ;
  public:
    const_Mat(const T *p ,int sz) : ptr(p),size(sz){}

    const_Mat_partial<T> operator[](int idx) const {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return const_Mat_partial<T>(ptr,size,idx) ;
    }

    //************************************************************************

    template<class S1, class S2> 
    void solve_lu(const S1 *b, S2 *x) const restrict {
      gridMotion::solve_lu(ptr,b,x,size) ;
    }

    //*************************************************************************

    template<class T1,class T2> 
    void solve_lu(const_Vect<T1> b, T2 * x) const  restrict {
      gridMotion::solve_lu(ptr,&b[0],x,size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu(const_Vect<T1> b, Vect<T2> x) const  restrict {
      gridMotion::solve_lu(ptr,&b[0],&x[0],size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu(const T1 *b, Vect<T2> x) const restrict {
      gridMotion::solve_lu(ptr,b,&x[0],size) ;
    }

    //************************************************************************

    template<class S> 
    void solve_lu_pivot(const S *b, S * x,const pivot_type * pivot) const restrict {
      gridMotion::solve_lu_pivot(ptr,b,x,pivot,size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu_pivot( const_Vect<T1> b, Vect<T2> x, 
                         const_Vect<pivot_type> pivot) const restrict {
      gridMotion::solve_lu_pivot(ptr,&b[0],&x[0],&pivot[0],size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu_pivot(const T1 *b, Vect<T2> x,
                        const pivot_type *pivot) const restrict {
      gridMotion::solve_lu_pivot(ptr,b,&x[0],pivot,size) ;
    }

    //************************************************************************

    template<class Tin,class Tout>
    void dotprod_accum(const Tin *vin, Tout *vout) const restrict {
      gridMotion::dotprod_accum(ptr,vin,vout,size) ;
    } 

    //************************************************************************

    template<class Tin,class Tout>
    void dotprod_accum(const Vect<Tin> vin, Tout *vout) const restrict {
      gridMotion::dotprod_accum(ptr,&vin[0],vout,size) ;
    }

    //************************************************************************

    template<class Tin,class Tout>
    void dotprod_accum(const const_Vect<Tin> vin, Tout *vout) const restrict {
      gridMotion::dotprod_accum(ptr,&vin[0],vout,size) ;
    }

    template<class Tin,class Tout>
    void dotprod_accum(const const_Vect<Tin> vin, Vect<Tout> vout) const restrict {
      gridMotion::dotprod_accum(ptr,&vin[0],&vout[0],size) ;
    }

  } ;

  //************************************************************************

  template <class T> class Mat_partial {
    T *restrict ptr ;
    int size ;
    int i ;
  public:
    Mat_partial(T *p,int sz, int ii) : ptr(p),size(sz),i(ii) {}
    T &operator[](int j) {
#ifdef BOUNDS_CHECK
      warn(j>=size || j < 0) ;
#endif
      return ptr[j*size + i] ;
    }
  } ;

  //**************************************************************************
  
  template <class T> class Mat {
  public:
    T *ptr ;
    int size ;
  public:
    Mat(T *p ,int sz) : ptr(p),size(sz) {}

    //------------------------------------------------------------------------
    void operator=(const Mat<T> &t) restrict {
      T *restrict p1 = ptr ;
      const T *restrict p2 = t.ptr ;
      const int sz2 = size*size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] = p2[i] ;
    }
    //------------------------------------------------------------------------
    void operator=(const const_Mat<T> &t) restrict {
      T *restrict p1 = ptr ;
      const T *restrict p2 = t.ptr ;
      const int sz2 = size*size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] = p2[i] ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator=(const Scalar<S> &restrict s) restrict {
      T *restrict p1 = ptr ;
      const int sz2 = size*size ;
      S val = s.val ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] = val ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator+=(const Scalar<S> &restrict s) restrict {
      T *restrict p1 = ptr ;
      const int sz2 = size*size ; ;
      S val = s.val ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] += val ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator*=(const Scalar<S> &restrict s) restrict {
      T *restrict p1 = ptr ;
      const int sz2 = size*size ;
      S val = s.val ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] *= val ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator-=(const Scalar<S> &restrict s) restrict {
      T *restrict p1 = ptr ;
      S val = s.val ;
      const int sz2 = size*size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] -= val ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator/=(const Scalar<S> &restrict s) restrict {
      T *restrict p1 = ptr ;
      S val = s.val ;
      const int sz2 = size*size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] /= val ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator+=(const Mat<S> &restrict t) restrict {
      T *restrict p1 = ptr ;
      const S *restrict p2 = t.ptr ;
      const int sz2 = size*size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] += p2[i] ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator+=(const const_Mat<S> &restrict t) restrict {
      T *restrict p1 = ptr ;
      const S *restrict p2 = t.ptr ;
      const int sz2 = size*size ;
      
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] += p2[i] ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator-=(const Mat<S> &restrict t) restrict {
      T *restrict p1 = ptr ;
      const S *restrict p2 = t.ptr ;
      const int sz2 = size*size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] -= p2[i] ;
    }
    //------------------------------------------------------------------------
    template <class S> 
    void operator-=(const const_Mat<S> &restrict t) restrict {
      T *restrict p1 = ptr ;
      const S *restrict p2 = t.ptr ;
      const int sz2 = size*size ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int i=0;i<sz2;++i)
        p1[i] -= p2[i] ;
    }
    //------------------------------------------------------------------------
    
    Mat_partial<T> operator[](int idx) restrict {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    //------------------------------------------------------------------------
    Mat_partial<T> operator[](int idx) const restrict {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    //------------------------------------------------------------------------

    void decompose_lu() restrict {
      // GAXPY from of LU decomposition algorithm 
      T *restrict Aj = ptr ;
      for(int j=0;j<size;++j,Aj += size) {
        const T *restrict Ak = ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
          for(int i=k+1;i<j;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }
        Ak = (T * restrict) ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
          for(int i=j;i<size;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }

        const T Ajjr = 1./Aj[j] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
        for(int i=j+1;i<size;++i)
          Aj[i] *= Ajjr ;
      }
    }
    //------------------------------------------------------------------------
    
    void decompose_lu_pivot(pivot_type *restrict pivot) restrict {
      pivot_type piv[256] ;  // Maximum matrix size for pivoting
      for(int i=0;i<size;++i)
        pivot[i] = i ;
      T *restrict Aj = ptr ;
      for(int j=0;j<size;++j,Aj+= size) {
        for(int k=0;k<j;++k)
          if(k!=piv[k])
            std::swap(Aj[k],Aj[piv[k]]) ;
        T *restrict Ak = ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
          for(int i=k+1;i<j;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }
        Ak = (T * restrict) ptr ;
        for(int k=0;k<j;++k,Ak += size) {
          const T Ajk = Aj[k] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
          for(int i=j;i<size;++i)
            Aj[i] -= Ak[i]*Ajk ;
        }
        int mu = j ;
        for(int k=j+1;k<size;++k)
          if(fabs(Aj[mu]) < fabs(Aj[k]))
            mu = k ;
        piv[j] = mu ;
        if(j!= mu)
          std::swap(pivot[j],pivot[mu]) ;
        Ak = (T * restrict) ptr ;
        for(int k=0;k<j+1;++k,Ak += size)
          if(j != piv[j])
            std::swap(Ak[j],Ak[piv[j]]) ;
        if(Aj[j] != 0) {
          T ajjr = 1./Aj[j] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
          for(int i=j+1;i<size;++i)
            Aj[i] *= ajjr ;
        }
      }
    }
    //------------------------------------------------------------------------

    template<class S1, class S2> 
    void solve_lu(const S1 * b, S2 *x) const restrict {
      gridMotion::solve_lu(ptr,b,x,size) ;
    }

    //------------------------------------------------------------------------

    template<class S1, class S2> 
    void solve_lu(const_Vect<S1> b, S2 * x) const restrict {
      gridMotion::solve_lu(ptr,&b[0],x,size) ;
    }

    //------------------------------------------------------------------------

    template<class S1, class S2> 
    void solve_lu(const S1 *b, Vect<S2> x) const restrict {
      gridMotion::solve_lu(ptr,b,&x[0],size) ;
    }

    //------------------------------------------------------------------------
    template<class S1, class S2> 
    void solve_lu(const_Vect<S1> b, Vect<S2> x) const restrict {
      gridMotion::solve_lu(ptr,&b[0],&x[0],size) ;
    }

    //------------------------------------------------------------------------

    template<class T1, class T2> 
    void solve_lu_pivot(const T1 *b, T2 *x,
                        const pivot_type *pivot) const restrict {
      gridMotion::solve_lu_pivot(ptr,b,x,pivot,size) ;
    }

    //------------------------------------------------------------------------

    template<class T1, class T2> 
    void solve_lu_pivot(const_Vect<T1> &b, T2 *x,
                        const pivot_type *restrict pivot) const restrict {
      gridMotion::solve_lu_pivot(ptr,&b[0],x,pivot,size) ;
    }
    //------------------------------------------------------------------------

    template<class T1, class T2> 
    void solve_lu_pivot(const T1 *b, Vect<T2> &x,
                        const pivot_type *pivot) const restrict {
      gridMotion::solve_lu_pivot(ptr,b,&x[0],pivot,size) ;
    }

    //------------------------------------------------------------------------

  } ;
}


#endif
namespace gridMotion {
  struct Quaternion {
    typedef vector3d<double> vect3d;
    typedef tensor3d<double> tens3d;
    typedef double real;
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
#endif
