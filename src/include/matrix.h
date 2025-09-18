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
#ifndef MATRIX_H
#define MATRIX_H 

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>

using std::pair ;
using std::make_pair ;

namespace Loci {
  typedef unsigned char pivot_type ;

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
  

  template<class T1, class T2, class T3> 
  inline void solve_qr(const T1 &A, 
                      const T2 *restrict b, 
                      T3 *restrict x, 
                      const pivot_type *restrict pivot, 
                      int size) {
    const int s = size;

    T2 b_temp[s];
    for(int i = 0; i < s; ++i) {
        b_temp[i] = b[i];
    }
    
    // Apply Q^T to b
    for (int k = 0; k < s - 1; ++k) {
        const T1 *restrict v_ptr = A.ptr + k + k * s;
        T1 dot_prod_v_b = 0.0;
        for (int i = k; i < s; ++i) {
            dot_prod_v_b += v_ptr[i] * b_temp[i];
        }
        T1 temp_val = 2.0 * dot_prod_v_b;
        for (int i = k; i < s; ++i) {
            b_temp[i] -= temp_val * v_ptr[i];
        }
    }
    
    // Back substitution to solve Ry = Q^T b
    for (int i = s - 1; i >= 0; --i) {
        T3 sum = 0.0;
        for (int j = i + 1; j < s; ++j) {
            sum += A.ptr[i + j * s] * x[j];
        }
        x[i] = (b_temp[i] - sum) / A.ptr[i + i * s];
    }
  }


  template<class T1, class T2, class T3> 
  inline void solve_qr_pivot(const T1 &A, 
                            const T2 *restrict b, 
                            T3 *restrict x, 
                            const pivot_type *restrict pivot, 
                            int size) 
  {
    const int s = size;

    T2 b_temp[s];
    for (int i = 0; i < s; ++i) {
        b_temp[i] = b[i];
    }
    
    // Apply Q^T to b
    for (int k = 0; k < s - 1; ++k) {
        const T1 *restrict v_ptr = A.ptr + k + k * s;
        T1 dot_prod_v_b = 0.0;
        for (int i = k; i < s; ++i) {
            dot_prod_v_b += v_ptr[i] * b_temp[i];
        }
        T1 temp_val = 2.0 * dot_prod_v_b;
        for (int i = k; i < s; ++i) {
            b_temp[i] -= temp_val * v_ptr[i];
        }
    }
    
    // Back substitution to solve Ry = Q^T b
    T3 y[s];
    for (int i = s - 1; i >= 0; --i) {
        T3 sum = 0.0;
        for (int j = i + 1; j < s; ++j) {
            sum += A_factored.ptr[i + j * s] * y[j];
        }
        y[i] = (b_temp[i] - sum) / A_factored.ptr[i + i * s];
    }
    

    // Permute the solution vector y to get the final solution x
    for (int i = 0; i < s; ++i) {
        x[pivot[i]] = y[i];
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
    const T & restrict operator[](size_t j) {
#ifdef BOUNDS_CHECK
      warn(j>=size || j < 0) ;
#endif
      return ptr[j*size + i] ;
    }
    const T & restrict operator[](unsigned int j) {
#ifdef BOUNDS_CHECK
      warn(j>=size || j < 0) ;
#endif
      return ptr[j*size + i] ;
    }
    const T & restrict operator[](unsigned char j) {
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
    const_Mat_partial<T> operator[](size_t idx) const {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return const_Mat_partial<T>(ptr,size,idx) ;
    }
    const_Mat_partial<T> operator[](unsigned int idx) const {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return const_Mat_partial<T>(ptr,size,idx) ;
    }
    const_Mat_partial<T> operator[](unsigned char idx) const {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return const_Mat_partial<T>(ptr,size,idx) ;
    }

    //************************************************************************

    template<class S1, class S2> 
    void solve_lu(const S1 *b, S2 *x) const restrict {
      Loci::solve_lu(ptr,b,x,size) ;
    }

    //*************************************************************************

    template<class T1,class T2> 
    void solve_lu(const_Vect<T1> b, T2 * x) const  restrict {
      Loci::solve_lu(ptr,&b[0],x,size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu(const_Vect<T1> b, Vect<T2> x) const  restrict {
      Loci::solve_lu(ptr,&b[0],&x[0],size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu(const T1 *b, Vect<T2> x) const restrict {
      Loci::solve_lu(ptr,b,&x[0],size) ;
    }

    //************************************************************************

    template<class S> 
    void solve_lu_pivot(const S *b, S * x,const pivot_type * pivot) const restrict {
      Loci::solve_lu_pivot(ptr,b,x,pivot,size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu_pivot( const_Vect<T1> b, Vect<T2> x, 
                         const_Vect<pivot_type> pivot) const restrict {
      Loci::solve_lu_pivot(ptr,&b[0],&x[0],&pivot[0],size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_lu_pivot(const T1 *b, Vect<T2> x,
                        const pivot_type *pivot) const restrict {
      Loci::solve_lu_pivot(ptr,b,&x[0],pivot,size) ;
    }
 // QR
     //************************************************************************

    template<class S1, class S2> 
    void solve_qr(const S1 *b, S2 *x) const restrict {
      Loci::solve_qr(ptr,b,x,size) ;
    }

    //*************************************************************************

    template<class T1,class T2> 
    void solve_qr(const_Vect<T1> b, T2 * x) const  restrict {
      Loci::solve_qr(ptr,&b[0],x,size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_qr(const_Vect<T1> b, Vect<T2> x) const  restrict {
      Loci::solve_qr(ptr,&b[0],&x[0],size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_qr(const T1 *b, Vect<T2> x) const restrict {
      Loci::solve_qr(ptr,b,&x[0],size) ;
    }

    //************************************************************************

    template<class S> 
    void solve_qr_pivot(const S *b, S * x,const pivot_type * pivot) const restrict {
      Loci::solve_qr_pivot(ptr,b,x,pivot,size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_qr_pivot( const_Vect<T1> b, Vect<T2> x, 
                         const_Vect<pivot_type> pivot) const restrict {
      Loci::solve_qr_pivot(ptr,&b[0],&x[0],&pivot[0],size) ;
    }

    //************************************************************************

    template<class T1,class T2> 
    void solve_qr_pivot(const T1 *b, Vect<T2> x,
                        const pivot_type *pivot) const restrict {
      Loci::solve_qr_pivot(ptr,b,&x[0],pivot,size) ;
    }
    
    //************************************************************************

    template<class Tin,class Tout>
    void dotprod_accum(const Tin *vin, Tout *vout) const restrict {
      Loci::dotprod_accum(ptr,vin,vout,size) ;
    } 

    //************************************************************************

    template<class Tin,class Tout>
    void dotprod_accum(const Vect<Tin> vin, Tout *vout) const restrict {
      Loci::dotprod_accum(ptr,&vin[0],vout,size) ;
    }

    //************************************************************************

    template<class Tin,class Tout>
    void dotprod_accum(const const_Vect<Tin> vin, Tout *vout) const restrict {
      Loci::dotprod_accum(ptr,&vin[0],vout,size) ;
    }

    template<class Tin,class Tout>
    void dotprod_accum(const const_Vect<Tin> vin, Vect<Tout> vout) const restrict {
      Loci::dotprod_accum(ptr,&vin[0],&vout[0],size) ;
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
    T &operator[](size_t j) {
#ifdef BOUNDS_CHECK
      warn(j>=size || j < 0) ;
#endif
      return ptr[j*size + i] ;
    }
    T &operator[](unsigned int j) {
#ifdef BOUNDS_CHECK
      warn(j>=size || j < 0) ;
#endif
      return ptr[j*size + i] ;
    }
    T &operator[](unsigned char j) {
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

    Mat_partial<T> operator[](int idx) const restrict {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    Mat_partial<T> operator[](size_t idx) restrict {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    Mat_partial<T> operator[](size_t idx) const restrict {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    Mat_partial<T> operator[](unsigned int idx) restrict {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    Mat_partial<T> operator[](unsigned int idx) const restrict {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    Mat_partial<T> operator[](unsigned char idx) restrict {
#ifdef BOUNDS_CHECK
      warn(idx >= size || idx < 0) ;
#endif
      return Mat_partial<T>(ptr,size,idx) ;
    }

    Mat_partial<T> operator[](unsigned char idx) const restrict {
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
        Ak = ptr ;
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
        Ak = ptr ;
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
          if(abs(Aj[mu]) < abs(Aj[k]))
            mu = k ;
        piv[j] = mu ;
        if(j!= mu)
          std::swap(pivot[j],pivot[mu]) ;
        Ak = ptr ;
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
    void decompose_qr() restrict {
        const int s = this->size;

        for (int k = 0; k < s - 1; ++k) {
            T *restrict a_k = this->ptr + k + k * s;
            T norm_a_k_sq = 0.0;
            for (int i = k; i < s; ++i) {
                norm_a_k_sq += a_k[i] * a_k[i];
            }
            T norm_a_k = std::sqrt(norm_a_k_sq);

            if (norm_a_k < 1e-12) continue;

            T alpha = (a_k[0] >= 0) ? norm_a_k : -norm_a_k;
            a_k[0] += alpha;
            T v_norm_sq = 0.0;
            for (int i = k; i < s; ++i) {
                v_norm_sq += a_k[i] * a_k[i];
            }
            T inv_v_norm = 1.0 / std::sqrt(v_norm_sq);

            for (int i = k; i < s; ++i) {
                a_k[i] *= inv_v_norm;
            }

            for (int j = k + 1; j < s; ++j) {
                T dot_prod_a_v = 0.0;
                for (int i = k; i < s; ++i) {
                    dot_prod_a_v += this->ptr[i + j * s] * this->ptr[i + k * s];
                }
                T temp_val = 2.0 * dot_prod_a_v;
                for (int i = k; i < s; ++i) {
                    this->ptr[i + j * s] -= temp_val * this->ptr[i + k * s];
                }
            }
        }
    }
    //------------------------------------------------------------------------
    
    void decompose_qr_pivot(pivot_type *restrict pivot) restrict {
      pivot_type piv[256] ;  // Maximum matrix size for pivoting
      const int s = this->size;
        
        // Initialize pivot array
        std::iota(pivot, pivot + s, 0);

        T col_norms[s];
        for (int j = 0; j < s; ++j) {
            T sum_sq = 0.0;
            for (int i = 0; i < s; ++i) {
                sum_sq += this->ptr[i + j * s] * this->ptr[i + j * s];
            }
            col_norms[j] = std::sqrt(sum_sq);
        }

        for (int k = 0; k < s - 1; ++k) {
            int pivot_col = k;
            T max_norm = col_norms[k];
            for (int j = k + 1; j < s; ++j) {
                if (col_norms[j] > max_norm) {
                    max_norm = col_norms[j];
                    pivot_col = j;
                }
            }
            
            if (pivot_col != k) {
                // Swap columns in matrix
                for (int i = 0; i < s; ++i) {
                    std::swap(this->ptr[i + k * s], this->ptr[i + pivot_col * s]);
                }
                // Swap norms and indices
                std::swap(col_norms[k], col_norms[pivot_col]);
                std::swap(pivot[k], pivot[pivot_col]);
            }

            T *restrict a_k = this->ptr + k + k * s;
            T norm_a_k_sq = 0.0;
            for (int i = k; i < s; ++i) {
                norm_a_k_sq += a_k[i] * a_k[i];
            }
            T norm_a_k = std::sqrt(norm_a_k_sq);

            if (norm_a_k < 1e-12) continue;
            
            T alpha = (a_k[0] >= 0) ? norm_a_k : -norm_a_k;
            a_k[0] += alpha;
            T v_norm_sq = 0.0;
            for (int i = k; i < s; ++i) {
                v_norm_sq += a_k[i] * a_k[i];
            }
            T inv_v_norm = 1.0 / std::sqrt(v_norm_sq);

            for (int i = k; i < s; ++i) {
                a_k[i] *= inv_v_norm;
            }
            
            for (int j = k + 1; j < s; ++j) {
                T dot_prod_a_v = 0.0;
                for (int i = k; i < s; ++i) {
                    dot_prod_a_v += this->ptr[i + j * s] * this->ptr[i + k * s];
                }
                T temp_val = 2.0 * dot_prod_a_v;
                for (int i = k; i < s; ++i) {
                    this->ptr[i + j * s] -= temp_val * this->ptr[i + k * s];
                }
            }
        }
    }

    //------------------------------------------------------------------------

    template<class S1, class S2> 
    void solve_lu(const S1 * b, S2 *x) const restrict {
      Loci::solve_lu(ptr,b,x,size) ;
    }

    //------------------------------------------------------------------------

    template<class S1, class S2> 
    void solve_lu(const_Vect<S1> b, S2 * x) const restrict {
      Loci::solve_lu(ptr,&b[0],x,size) ;
    }

    //------------------------------------------------------------------------

    template<class S1, class S2> 
    void solve_lu(const S1 *b, Vect<S2> x) const restrict {
      Loci::solve_lu(ptr,b,&x[0],size) ;
    }

    //------------------------------------------------------------------------
    template<class S1, class S2> 
    void solve_lu(const_Vect<S1> b, Vect<S2> x) const restrict {
      Loci::solve_lu(ptr,&b[0],&x[0],size) ;
    }

    //------------------------------------------------------------------------

    template<class T1, class T2> 
    void solve_lu_pivot(const T1 *b, T2 *x,
                        const pivot_type *pivot) const restrict {
      Loci::solve_lu_pivot(ptr,b,x,pivot,size) ;
    }

    //------------------------------------------------------------------------

    template<class T1, class T2> 
    void solve_lu_pivot(const_Vect<T1> &b, T2 *x,
                        const pivot_type *restrict pivot) const restrict {
      Loci::solve_lu_pivot(ptr,&b[0],x,pivot,size) ;
    }
    //------------------------------------------------------------------------

    template<class T1, class T2> 
    void solve_lu_pivot(const T1 *b, Vect<T2> &x,
                        const pivot_type *pivot) const restrict {
      Loci::solve_lu_pivot(ptr,b,&x[0],pivot,size) ;
    }

    //========================================================================
    
    template<class S1, class S2> 
    void solve_qr(const S1 * b, S2 *x) const restrict {
      Loci::solve_qr(ptr,b,x,size) ;
    }

    //------------------------------------------------------------------------

    template<class S1, class S2> 
    void solve_qr(const_Vect<S1> b, S2 * x) const restrict {
      Loci::solve_qr(ptr,&b[0],x,size) ;
    }

    //------------------------------------------------------------------------

    template<class S1, class S2> 
    void solve_qr(const S1 *b, Vect<S2> x) const restrict {
      Loci::solve_qr(ptr,b,&x[0],size) ;
    }

    //------------------------------------------------------------------------
    template<class S1, class S2> 
    void solve_qr(const_Vect<S1> b, Vect<S2> x) const restrict {
      Loci::solve_qr(ptr,&b[0],&x[0],size) ;
    }

    //------------------------------------------------------------------------

    template<class T1, class T2> 
    void solve_qr_pivot(const T1 *b, T2 *x,
                        const pivot_type *pivot) const restrict {
      Loci::solve_qr_pivot(ptr,b,x,pivot,size) ;
    }

    //------------------------------------------------------------------------

    template<class T1, class T2> 
    void solve_qr_pivot(const_Vect<T1> &b, T2 *x,
                        const pivot_type *restrict pivot) const restrict {
      Loci::solve_qr_pivot(ptr,&b[0],x,pivot,size) ;
    }
    //------------------------------------------------------------------------

    template<class T1, class T2> 
    void solve_qr_pivot(const T1 *b, Vect<T2> &x,
                        const pivot_type *pivot) const restrict {
      Loci::solve_qr_pivot(ptr,b,&x[0],pivot,size) ;
    }
  } ;

}

#endif

  
