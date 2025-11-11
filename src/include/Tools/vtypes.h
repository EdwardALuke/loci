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
#ifndef LOCI_VTYPES_H
#define LOCI_VTYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#if defined(__x86_64__)

#if defined(__AVX2__)
#include <immintrin.h>
#else 
#warning "X86_64 architecture without AVX2 support, vtype not vectorized!, include -mavx2 in compiler options"
#define VTYPE_SCALAR_ONLY
#endif

#endif

#if !defined(__x86_64__)
// compatibility for non x86_64 architecures
#define VTYPE_SCALAR_ONLY
#endif

#include <cstdint>
#include <cmath>
#include <initializer_list>


// TODO
// 1. What should be the order of registers in a multi-register member of the union in a vtype?
// 2. Is it useful to make vtype types as literals? (requires c++11)
// 3. Check the assembly code generated from vtype expressions.
// 4. Write unit test suites, comparing SIMD results with scalar results.
// 5. Write performance tests (latency and throughput) for vtype operations.

namespace Loci {
  // This is a wapper around 1D statically sized arrays in C (similar to
  // std::array, equivalent to Loci::Array). This has been reinvented
  // here since it may need to be ported to GPGPU, but this needs to be
  // replaced with Loci::Array once Loci::Array is ported to GPGPU.
  template<typename T, int N>
  class alignas(alignof(T)) StaticSizedArray {
  private:
    T data[N];

  public:
    StaticSizedArray() = default ;

    StaticSizedArray(T e) {
      for(int i = 0; i < N; ++i) {
        data[i] = e;
      }
    }

    StaticSizedArray(std::initializer_list<T> l) {
      auto iter = l.begin();
      for(int i = 0; i < N && iter != l.end(); ++i) {
        data[i] = *iter++;
      }
    }

    StaticSizedArray(StaticSizedArray<T, N> const & src) {
      for(int i = 0; i < N; ++i) {
        data[i] = src.data[i];
      }
    }

    StaticSizedArray & operator=(StaticSizedArray<T, N> const & src) {
      for(int i = 0; i < N; ++i) {
        data[i] = src.data[i];
      }
      return *this;
    }

    StaticSizedArray & operator=(std::initializer_list<T> l) {
      auto iter = l.begin();
      for(int i = 0; i < N && iter != l.end(); ++i) {
        data[i] = *iter++;
      }
      return *this;
    }

    T & operator[](unsigned int index) {
      return data[index];
    }

    T const & operator[](unsigned int index) const {
      return data[index];
    }
  };

  /**
   * Generic vectorized type.
   */
  template<typename T, int W>
  struct vtype {
  };

  /**
   * Vectorized 8x int32
   */
  template<>
  struct vtype<int32_t, 8> {
    static constexpr unsigned int W = 8;

    union {
      StaticSizedArray<int32_t, 8> data;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      __m256i avxReg;
#endif

    };

    vtype() = default ;

    vtype(
          int32_t e0, int32_t e1, int32_t e2, int32_t e3,
          int32_t e4, int32_t e5, int32_t e6, int32_t e7
          )
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : avxReg(_mm256_set_epi32(e7, e6, e5, e4, e3, e2, e1, e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data({e0, e1, e2, e3, e4, e5, e6, e7}) { }
#endif

    vtype(int32_t e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : avxReg(_mm256_set1_epi32(e)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(e) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    explicit vtype(__m256i reg) : avxReg(reg) { }
#endif

    vtype(vtype<int32_t, 8> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : avxReg(rhs.avxReg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(rhs.data) { }
#endif

    vtype<int32_t, 8> & operator=(vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      avxReg = rhs.avxReg;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
      data = rhs.data;
#endif
      return *this;
    }

    int32_t & operator[](unsigned int index) {
      return data[index];
    }

    int32_t const & operator[](unsigned int index) const {
      return data[index];
    }

    vtype<int32_t, 8> & operator+=(vtype<int32_t, 8> const & rhs);
    vtype<int32_t, 8> & operator-=(vtype<int32_t, 8> const & rhs);
    vtype<int32_t, 8> & operator*=(vtype<int32_t, 8> const & rhs);
    vtype<int32_t, 8> & operator/=(vtype<int32_t, 8> const & rhs);
    vtype<int32_t, 8> & operator%=(vtype<int32_t, 8> const & rhs);

    vtype<int32_t, 8> & operator+=(int32_t rhs);
    vtype<int32_t, 8> & operator-=(int32_t rhs);
    vtype<int32_t, 8> & operator*=(int32_t rhs);
    vtype<int32_t, 8> & operator/=(int32_t rhs);
    vtype<int32_t, 8> & operator%=(int32_t rhs);

    // variable bit manipulation operations
    vtype<int32_t, 8> & operator<<=(vtype<int32_t, 8> const & rhs);
    vtype<int32_t, 8> & operator>>=(vtype<int32_t, 8> const & rhs);

    // uniform bit manipulation operations
    vtype<int32_t, 8> & operator<<=(int32_t shift);
    vtype<int32_t, 8> & operator>>=(int32_t shift);

    // constant bit manipulation operations
    // constexpr vtype<int32_t, 8> & operator<<=(int32_t const shift);
    // constexpr vtype<int32_t, 8> & operator>>=(int32_t const shift);
  };


  inline vtype<int32_t, 8> & vtype<int32_t, 8>::operator+=(vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_add_epi32(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
    data[3] += rhs.data[3];
    data[4] += rhs.data[4];
    data[5] += rhs.data[5];
    data[6] += rhs.data[6];
    data[7] += rhs.data[7];
#endif
    return *this;
  }

  inline vtype<int32_t, 8> & vtype<int32_t, 8>::operator-=(vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_sub_epi32(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] -= rhs.data[0];
    data[1] -= rhs.data[1];
    data[2] -= rhs.data[2];
    data[3] -= rhs.data[3];
    data[4] -= rhs.data[4];
    data[5] -= rhs.data[5];
    data[6] -= rhs.data[6];
    data[7] -= rhs.data[7];
#endif
    return *this;
  }

  inline vtype<int32_t, 8> & vtype<int32_t, 8>::operator*=(vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_mullo_epi32(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] *= rhs.data[0];
    data[1] *= rhs.data[1];
    data[2] *= rhs.data[2];
    data[3] *= rhs.data[3];
    data[4] *= rhs.data[4];
    data[5] *= rhs.data[5];
    data[6] *= rhs.data[6];
    data[7] *= rhs.data[7];
#endif
    return *this;
  }
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
  inline void avx2_si32_si32_quorem(
                                    __m256i lhs, __m256i rhs,
                                    __m256i & quotient, __m256i & remainder
                                    ) {
#  if defined(__INTEL_COMPILER)
    quotient = _mm256_idivrem_epi32(&remainder, lhs, rhs);
#  else
    quotient = _mm256_setzero_si256();
    remainder = _mm256_setzero_si256();
    __m256i dividend = _mm256_abs_epi32(lhs);
    __m256i divisor = _mm256_abs_epi32(rhs);
    __m256i signDividend = _mm256_cmpeq_epi32(dividend, lhs);
    __m256i signDivisor = _mm256_cmpeq_epi32(divisor, rhs);
    __m256i sign = _mm256_xor_si256(signDividend, signDivisor);
    for(int i = 30; i >= 0; --i) {
      quotient = _mm256_slli_epi32(quotient, 1);
      remainder = _mm256_slli_epi32(remainder, 1);
      __m256i shift = _mm256_set1_epi32(i);
      __m256i mask = _mm256_sllv_epi32(_mm256_set1_epi32(1), shift);
      remainder = _mm256_or_si256(remainder, _mm256_srlv_epi32(_mm256_and_si256(dividend, mask), shift));
      __m256i mask1 = _mm256_or_si256(_mm256_cmpgt_epi32(remainder, divisor), _mm256_cmpeq_epi32(remainder, divisor));
      remainder = _mm256_sub_epi32(remainder, _mm256_and_si256(divisor, mask1));
      quotient = _mm256_or_si256(quotient, _mm256_and_si256(_mm256_set1_epi32(1), mask1));
    }
    remainder = _mm256_add_epi32(_mm256_xor_si256(remainder, signDividend), _mm256_and_si256(signDividend, _mm256_set1_epi32(1)));
    quotient = _mm256_add_epi32(_mm256_xor_si256(quotient, sign), _mm256_and_si256(sign, _mm256_set1_epi32(1)));
#  endif
  }
#endif

  inline vtype<int32_t, 8> & vtype<int32_t, 8>::operator/=(vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256i r;
    avx2_si32_si32_quorem(avxReg, rhs.avxReg, avxReg, r);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] /= rhs.data[0];
    data[1] /= rhs.data[1];
    data[2] /= rhs.data[2];
    data[3] /= rhs.data[3];
    data[4] /= rhs.data[4];
    data[5] /= rhs.data[5];
    data[6] /= rhs.data[6];
    data[7] /= rhs.data[7];
#endif
    return *this;
  }

  inline vtype<int32_t, 8> & vtype<int32_t, 8>::operator%=(vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256i q;
    avx2_si32_si32_quorem(avxReg, rhs.avxReg, q, avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
#endif
    data[0] %= rhs.data[0];
    data[1] %= rhs.data[1];
    data[2] %= rhs.data[2];
    data[3] %= rhs.data[3];
    data[4] %= rhs.data[4];
    data[5] %= rhs.data[5];
    data[6] %= rhs.data[6];
    data[7] %= rhs.data[7];
#endif
    return *this;
  }

  inline vtype<int32_t, 8> & vtype<int32_t, 8>::operator<<=(vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_sllv_epi32(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] <<= rhs.data[0];
    data[1] <<= rhs.data[1];
    data[2] <<= rhs.data[2];
    data[3] <<= rhs.data[3];
    data[4] <<= rhs.data[4];
    data[5] <<= rhs.data[5];
    data[6] <<= rhs.data[6];
    data[7] <<= rhs.data[7];
#endif
    return *this;
  }

  inline vtype<int32_t, 8> & vtype<int32_t, 8>::operator>>=(vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_srlv_epi32(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] >>= rhs.data[0];
    data[1] >>= rhs.data[1];
    data[2] >>= rhs.data[2];
    data[3] >>= rhs.data[3];
    data[4] >>= rhs.data[4];
    data[5] >>= rhs.data[5];
    data[6] >>= rhs.data[6];
    data[7] >>= rhs.data[7];
#endif
    return *this;
  }

  inline vtype<int32_t, 8> & vtype<int32_t, 8>::operator<<=(int32_t shift) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m128i sh = _mm_set1_epi64x(shift);
    avxReg = _mm256_sll_epi32(avxReg, sh);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] <<= shift;
    data[1] <<= shift;
    data[2] <<= shift;
    data[3] <<= shift;
    data[4] <<= shift;
    data[5] <<= shift;
    data[6] <<= shift;
    data[7] <<= shift;
#endif
    return *this;
  }

  inline vtype<int32_t, 8> & vtype<int32_t, 8>::operator>>=(int32_t shift) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m128i sh = _mm_set1_epi64x(shift);
    avxReg = _mm256_srl_epi32(avxReg, sh);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    for(int i = 0; i < 8; ++i) {
      data[i] >>= shift;
    }
#endif
    return *this;
  }



  inline vtype<int32_t, 8> operator+(vtype<int32_t, 8> const & lhs, vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<int32_t, 8>(_mm256_add_epi32(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<int32_t, 8>(
                             lhs.data[0] + rhs.data[0],
                             lhs.data[1] + rhs.data[1],
                             lhs.data[2] + rhs.data[2],
                             lhs.data[3] + rhs.data[3],
                             lhs.data[4] + rhs.data[4],
                             lhs.data[5] + rhs.data[5],
                             lhs.data[6] + rhs.data[6],
                             lhs.data[7] + rhs.data[7]
                             );
#endif
  }

  inline vtype<int32_t, 8> operator-(vtype<int32_t, 8> const & lhs, vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<int32_t, 8>(_mm256_sub_epi32(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<int32_t, 8>(
                             lhs.data[0] - rhs.data[0],
                             lhs.data[1] - rhs.data[1],
                             lhs.data[2] - rhs.data[2],
                             lhs.data[3] - rhs.data[3],
                             lhs.data[4] - rhs.data[4],
                             lhs.data[5] - rhs.data[5],
                             lhs.data[6] - rhs.data[6],
                             lhs.data[7] - rhs.data[7]
                             );
#endif
  }

  inline vtype<int32_t, 8> operator*(vtype<int32_t, 8> const & lhs, vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<int32_t, 8>(_mm256_mullo_epi32(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<int32_t, 8>(
                             lhs.data[0] * rhs.data[0],
                             lhs.data[1] * rhs.data[1],
                             lhs.data[2] * rhs.data[2],
                             lhs.data[3] * rhs.data[3],
                             lhs.data[4] * rhs.data[4],
                             lhs.data[5] * rhs.data[5],
                             lhs.data[6] * rhs.data[6],
                             lhs.data[7] * rhs.data[7]
                             );
#endif
  }

  inline vtype<int32_t, 8> operator/(vtype<int32_t, 8> const & lhs, vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256i q, r;
    avx2_si32_si32_quorem(lhs.avxReg, rhs.avxReg, q, r);
    return vtype<int32_t, 8>(q);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<int32_t, 8>(
                             lhs.data[0] / rhs.data[0],
                             lhs.data[1] / rhs.data[1],
                             lhs.data[2] / rhs.data[2],
                             lhs.data[3] / rhs.data[3],
                             lhs.data[4] / rhs.data[4],
                             lhs.data[5] / rhs.data[5],
                             lhs.data[6] / rhs.data[6],
                             lhs.data[7] / rhs.data[7]
                             );
#endif
  }


  inline vtype<int32_t, 8> operator%(vtype<int32_t, 8> const & lhs, vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256i q, r;
    avx2_si32_si32_quorem(lhs.avxReg, rhs.avxReg, q, r);
    return vtype<int32_t, 8>(r);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<int32_t, 8>(
                             lhs.data[0] % rhs.data[0],
                             lhs.data[1] % rhs.data[1],
                             lhs.data[2] % rhs.data[2],
                             lhs.data[3] % rhs.data[3],
                             lhs.data[4] % rhs.data[4],
                             lhs.data[5] % rhs.data[5],
                             lhs.data[6] % rhs.data[6],
                             lhs.data[7] % rhs.data[7]
                             );
#endif
  }


  /**
   * Vectorized 4x float.
   */
  template<>
  struct vtype<float, 4> {
    static constexpr unsigned int W = 4;

    union {
      StaticSizedArray<float, 4> data;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      __m128 sseReg;
#endif
    };

    vtype() = default ;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    explicit vtype(__m128 reg) : sseReg(reg) { }
#endif

    vtype(vtype<float, 4> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : sseReg(rhs.sseReg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(rhs.data) { }
#endif

    vtype(float e0, float e1, float e2, float e3)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : sseReg(_mm_set_ps(e3, e2, e1, e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data{e0, e1, e2, e3} { }
#endif

    explicit vtype(float e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : sseReg(_mm_set1_ps(e)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(e) { }
#endif

    vtype<float, 4> & operator=(vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      sseReg = rhs.sseReg;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
      data = rhs.data;
#endif
      return *this;
    }

    float & operator[](unsigned int index) {
      return data[index];
    }

    float const & operator[](unsigned int index) const {
      return data[index];
    }

    vtype<float, 4> & operator+=(vtype<float, 4> const & rhs);
    vtype<float, 4> & operator-=(vtype<float, 4> const & rhs);
    vtype<float, 4> & operator*=(vtype<float, 4> const & rhs);
    vtype<float, 4> & operator/=(vtype<float, 4> const & rhs);
  };

  inline vtype<float, 4> & vtype<float, 4>::operator+=(vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    sseReg = _mm_add_ps(sseReg, rhs.sseReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
    data[3] += rhs.data[3];
#endif
    return *this;
  }

  inline vtype<float, 4> & vtype<float, 4>::operator-=(vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    sseReg = _mm_sub_ps(sseReg, rhs.sseReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] -= rhs.data[0];
    data[1] -= rhs.data[1];
    data[2] -= rhs.data[2];
    data[3] -= rhs.data[3];
#endif
    return *this;
  }

  inline vtype<float, 4> & vtype<float, 4>::operator*=(vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    sseReg = _mm_mul_ps(sseReg, rhs.sseReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] *= rhs.data[0];
    data[1] *= rhs.data[1];
    data[2] *= rhs.data[2];
    data[3] *= rhs.data[3];
#endif
    return *this;
  }

  inline vtype<float, 4> & vtype<float, 4>::operator/=(vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    sseReg = _mm_div_ps(sseReg, rhs.sseReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] /= rhs.data[0];
    data[1] /= rhs.data[1];
    data[2] /= rhs.data[2];
    data[3] /= rhs.data[3];
#endif
    return *this;
  }

  inline vtype<float, 4> operator+(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_add_ps(lhs.sseReg, rhs.sseReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(
                           lhs.data[0] + rhs.data[0],
                           lhs.data[1] + rhs.data[1],
                           lhs.data[2] + rhs.data[2],
                           lhs.data[3] + rhs.data[3]
                           );
#endif
  }

  inline vtype<float, 4> operator-(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_sub_ps(lhs.sseReg, rhs.sseReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(
                           lhs.data[0] - rhs.data[0],
                           lhs.data[1] - rhs.data[1],
                           lhs.data[2] - rhs.data[2],
                           lhs.data[3] - rhs.data[3]
                           );
#endif
  }

  inline vtype<float, 4> operator*(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_mul_ps(lhs.sseReg, rhs.sseReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(
                           lhs.data[0] * rhs.data[0],
                           lhs.data[1] * rhs.data[1],
                           lhs.data[2] * rhs.data[2],
                           lhs.data[3] * rhs.data[3]
                           );
#endif
  }

  inline vtype<float, 4> operator/(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_div_ps(lhs.sseReg, rhs.sseReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(
                           lhs.data[0] / rhs.data[0],
                           lhs.data[1] / rhs.data[1],
                           lhs.data[2] / rhs.data[2],
                           lhs.data[3] / rhs.data[3]
                           );
#endif
  }

  inline vtype<float, 4> ceil(vtype<float, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_ceil_ps(arg.sseReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(
                           std::ceil(arg.data[0]),
                           std::ceil(arg.data[1]),
                           std::ceil(arg.data[2]),
                           std::ceil(arg.data[3])
                           );
#endif
  }

  inline vtype<float, 4> floor(vtype<float, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_floor_ps(arg.sseReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(
                           std::floor(arg.data[0]),
                           std::floor(arg.data[1]),
                           std::floor(arg.data[2]),
                           std::floor(arg.data[3])
                           );
#endif
  }

  inline vtype<float, 4> round(vtype<float, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_round_ps(arg.sseReg, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(
                           std::round(arg.data[0]),
                           std::round(arg.data[1]),
                           std::round(arg.data[2]),
                           std::round(arg.data[3])
                           );
#endif
  }

  inline vtype<float, 4> sqrt(vtype<float, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_sqrt_ps(arg.sseReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(
                           std::sqrt(arg.data[0]),
                           std::sqrt(arg.data[1]),
                           std::sqrt(arg.data[2]),
                           std::sqrt(arg.data[3])
                           ) ;
#endif
      }

  inline vtype<float, 4> rsqrt(vtype<float, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_rsqrt_ps(arg.sseReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(
                           1.0f/std::sqrt(arg.data[0]),
                           1.0f/std::sqrt(arg.data[1]),
                           1.0f/std::sqrt(arg.data[2]),
                           1.0f/std::sqrt(arg.data[3])
                           );
#endif
  }

  /**
   * Vectorized 8x float.
   */
  template<>
  struct vtype<float, 8> {
    static constexpr unsigned int W = 8;

    union {
      StaticSizedArray<float, 8> data;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      __m256 avxReg;
#endif

    };

    vtype() = default ;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    explicit vtype(__m256 reg) : avxReg(reg) { }
#endif


    vtype(vtype<float, 8> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : avxReg(rhs.avxReg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(rhs.data) { }
#endif

    vtype(float e0, float e1, float e2, float e3, float e4, float e5, float e6, float e7)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : avxReg(_mm256_set_ps(e7, e6, e5, e4, e3, e2, e1, e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    : data{e0, e1, e2, e3, e4, e5, e6, e7} { }
#endif

    explicit vtype(float e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : avxReg(_mm256_set1_ps(e)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(e) { }
#endif

    vtype<float, 8> & operator=(vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      avxReg = rhs.avxReg;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data = rhs.data;
#endif
      return *this;
    }

    float & operator[](unsigned int index) {
      return data[index];
    }

    float const & operator[](unsigned int index) const {
      return data[index];
    }

    vtype<float, 8> & operator+=(vtype<float, 8> const & rhs);
    vtype<float, 8> & operator-=(vtype<float, 8> const & rhs);
    vtype<float, 8> & operator*=(vtype<float, 8> const & rhs);
    vtype<float, 8> & operator/=(vtype<float, 8> const & rhs);
  };

  inline vtype<float, 8> & vtype<float, 8>::operator+=(vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_add_ps(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
    data[3] += rhs.data[3];
    data[4] += rhs.data[4];
    data[5] += rhs.data[5];
    data[6] += rhs.data[6];
    data[7] += rhs.data[7];
#endif
    return *this;
  }

  inline vtype<float, 8> & vtype<float, 8>::operator-=(vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_sub_ps(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] -= rhs.data[0];
    data[1] -= rhs.data[1];
    data[2] -= rhs.data[2];
    data[3] -= rhs.data[3];
    data[4] -= rhs.data[4];
    data[5] -= rhs.data[5];
    data[6] -= rhs.data[6];
    data[7] -= rhs.data[7];
#endif
    return *this;
  }

  inline vtype<float, 8> & vtype<float, 8>::operator*=(vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_mul_ps(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] *= rhs.data[0];
    data[1] *= rhs.data[1];
    data[2] *= rhs.data[2];
    data[3] *= rhs.data[3];
    data[4] *= rhs.data[4];
    data[5] *= rhs.data[5];
    data[6] *= rhs.data[6];
    data[7] *= rhs.data[7];
#endif
    return *this;
  }

  inline vtype<float, 8> & vtype<float, 8>::operator/=(vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_div_ps(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] /= rhs.data[0];
    data[1] /= rhs.data[1];
    data[2] /= rhs.data[2];
    data[3] /= rhs.data[3];
    data[4] /= rhs.data[4];
    data[5] /= rhs.data[5];
    data[6] /= rhs.data[6];
    data[7] /= rhs.data[7];
#endif
    return *this;
  }

  inline vtype<float, 8> operator+(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_add_ps(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           lhs.data[0] + rhs.data[0],
                           lhs.data[1] + rhs.data[1],
                           lhs.data[2] + rhs.data[2],
                           lhs.data[3] + rhs.data[3],
                           lhs.data[4] + rhs.data[4],
                           lhs.data[5] + rhs.data[5],
                           lhs.data[6] + rhs.data[6],
                           lhs.data[7] + rhs.data[7]
                           );
#endif
  }

  inline vtype<float, 8> operator-(vtype<float, 8> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256 sign_bit = _mm256_set1_ps(-0.0) ;
    return vtype<float, 8>(_mm256_xor_ps(sign_bit, v.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           - v.data[0],
                           - v.data[1],
                           - v.data[2],
                           - v.data[3],
                           - v.data[4],
                           - v.data[5],
                           - v.data[6],
                           - v.data[7]
                           );
#endif
  }

  inline vtype<float, 8> operator-(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_sub_ps(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           lhs.data[0] - rhs.data[0],
                           lhs.data[1] - rhs.data[1],
                           lhs.data[2] - rhs.data[2],
                           lhs.data[3] - rhs.data[3],
                           lhs.data[4] - rhs.data[4],
                           lhs.data[5] - rhs.data[5],
                           lhs.data[6] - rhs.data[6],
                           lhs.data[7] - rhs.data[7]
                           );
#endif
  }

  inline vtype<float, 8> operator*(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_mul_ps(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           lhs.data[0] * rhs.data[0],
                           lhs.data[1] * rhs.data[1],
                           lhs.data[2] * rhs.data[2],
                           lhs.data[3] * rhs.data[3],
                           lhs.data[4] * rhs.data[4],
                           lhs.data[5] * rhs.data[5],
                           lhs.data[6] * rhs.data[6],
                           lhs.data[7] * rhs.data[7]
                           );
#endif
  }

  inline vtype<float, 8> operator*(vtype<float, 8> const & lhs, float const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256 C = _mm256_set1_ps(rhs) ;
    return vtype<float, 8>(_mm256_mul_ps(lhs.avxReg, C));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           lhs.data[0] * rhs,
                           lhs.data[1] * rhs,
                           lhs.data[2] * rhs,
                           lhs.data[3] * rhs,
                           lhs.data[4] * rhs,
                           lhs.data[5] * rhs,
                           lhs.data[6] * rhs,
                           lhs.data[7] * rhs
                           );
#endif
  }

  inline vtype<float, 8> operator*(float const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256 C = _mm256_set1_ps(lhs) ;
    return vtype<float, 8>(_mm256_mul_ps(C, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           lhs * rhs.data[0],
                           lhs * rhs.data[1],
                           lhs * rhs.data[2],
                           lhs * rhs.data[3],
                           lhs * rhs.data[4],
                           lhs * rhs.data[5],
                           lhs * rhs.data[6],
                           lhs * rhs.data[7]
                           );
#endif
  }

  inline vtype<float, 8> operator*(vtype<float, 8> const & lhs, double const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256 C = _mm256_set1_ps(float(rhs)) ;
    return vtype<float, 8>(_mm256_mul_ps(lhs.avxReg, C));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           lhs.data[0] * rhs,
                           lhs.data[1] * rhs,
                           lhs.data[2] * rhs,
                           lhs.data[3] * rhs,
                           lhs.data[4] * rhs,
                           lhs.data[5] * rhs,
                           lhs.data[6] * rhs,
                           lhs.data[7] * rhs
                           );
#endif
  }

  inline vtype<float, 8> operator*(double const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256 C = _mm256_set1_ps(float(lhs)) ;
    return vtype<float, 8>(_mm256_mul_ps(C, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           lhs * rhs.data[0],
                           lhs * rhs.data[1],
                           lhs * rhs.data[2],
                           lhs * rhs.data[3],
                           lhs * rhs.data[4],
                           lhs * rhs.data[5],
                           lhs * rhs.data[6],
                           lhs * rhs.data[7]
                           );
#endif
  }
  
  inline vtype<float, 8> operator/(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_div_ps(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           lhs.data[0] / rhs.data[0],
                           lhs.data[1] / rhs.data[1],
                           lhs.data[2] / rhs.data[2],
                           lhs.data[3] / rhs.data[3],
                           lhs.data[4] / rhs.data[4],
                           lhs.data[5] / rhs.data[5],
                           lhs.data[6] / rhs.data[6],
                           lhs.data[7] / rhs.data[7]
                           );
#endif
  }

  inline vtype<float, 8> operator/(vtype<float, 8> const & lhs, float const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256 C = _mm256_set1_ps(rhs) ;
    return vtype<float, 8>(_mm256_div_ps(lhs.avxReg, C));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           lhs.data[0] / rhs,
                           lhs.data[1] / rhs,
                           lhs.data[2] / rhs,
                           lhs.data[3] / rhs,
                           lhs.data[4] / rhs,
                           lhs.data[5] / rhs,
                           lhs.data[6] / rhs,
                           lhs.data[7] / rhs
                           );
#endif
  }

  inline vtype<float, 8> ceil(vtype<float, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_ceil_ps(arg.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           std::ceil(arg.data[0]),
                           std::ceil(arg.data[1]),
                           std::ceil(arg.data[2]),
                           std::ceil(arg.data[3]),
                           std::ceil(arg.data[4]),
                           std::ceil(arg.data[5]),
                           std::ceil(arg.data[6]),
                           std::ceil(arg.data[7])
                           );
#endif
  }

  inline vtype<float, 8> floor(vtype<float, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_floor_ps(arg.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           std::floor(arg.data[0]),
                           std::floor(arg.data[1]),
                           std::floor(arg.data[2]),
                           std::floor(arg.data[3]),
                           std::floor(arg.data[4]),
                           std::floor(arg.data[5]),
                           std::floor(arg.data[6]),
                           std::floor(arg.data[7])
                           );
#endif
  }

  inline vtype<float, 8> round(vtype<float, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_round_ps(arg.avxReg, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           std::round(arg.data[0]),
                           std::round(arg.data[1]),
                           std::round(arg.data[2]),
                           std::round(arg.data[3]),
                           std::round(arg.data[4]),
                           std::round(arg.data[5]),
                           std::round(arg.data[6]),
                           std::round(arg.data[7])
                           );
#endif
  }

  inline vtype<float, 8> sqrt(vtype<float, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_sqrt_ps(arg.avxReg));
#else
#if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           std::sqrt(arg.data[0]),
                           std::sqrt(arg.data[1]),
                           std::sqrt(arg.data[2]),
                           std::sqrt(arg.data[3]),
                           std::sqrt(arg.data[4]),
                           std::sqrt(arg.data[5]),
                           std::sqrt(arg.data[6]),
                           std::sqrt(arg.data[7])
                           );
#endif
  }

  inline vtype<float, 8> rsqrt(vtype<float, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_rsqrt_ps(arg.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(
                           1.0f/std::sqrt(arg.data[0]),
                           1.0f/std::sqrt(arg.data[1]),
                           1.0f/std::sqrt(arg.data[2]),
                           1.0f/std::sqrt(arg.data[3]),
                           1.0f/std::sqrt(arg.data[4]),
                           1.0f/std::sqrt(arg.data[5]),
                           1.0f/std::sqrt(arg.data[6]),
                           1.0f/std::sqrt(arg.data[7])
                           );
#endif
  }

  /**
   * Vectorized 4x double.
   */
  template<>
  struct vtype<double, 4> {
    static constexpr unsigned int W = 4;

    union {
      StaticSizedArray<double, 4> data;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      __m256d avxReg;
#endif

    };

    vtype() = default ;

    vtype(
          double e0, double e1, double e2, double e3
          )
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : avxReg(_mm256_set_pd(e3, e2, e1, e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    : data({e0, e1, e2, e3}) { }
#endif

    explicit vtype(double e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : avxReg(_mm256_set1_pd(e)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(e) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    explicit vtype(__m256d avxReg) : avxReg(avxReg) { }
#endif

    vtype(vtype<double, 4> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : avxReg(rhs.avxReg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    : data(rhs.data) { }
#endif

    vtype<double, 4> & operator=(vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      avxReg = rhs.avxReg;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data = rhs.data;
#endif
      return *this;
    }

    double & operator[](unsigned int index) {
      return data[index];
    }

    double const & operator[](unsigned int index) const {
      return data[index];
    }

    vtype<double, 4> & operator+=(vtype<double, 4> const & rhs);
    vtype<double, 4> & operator-=(vtype<double, 4> const & rhs);
    vtype<double, 4> & operator*=(vtype<double, 4> const & rhs);
    vtype<double, 4> & operator/=(vtype<double, 4> const & rhs);
  };

  inline vtype<double, 4> & vtype<double, 4>::operator+=(vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_add_pd(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
    data[3] += rhs.data[3];
#endif
    return *this;
  }

  inline vtype<double, 4> & vtype<double, 4>::operator-=(vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_sub_pd(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] -= rhs.data[0];
    data[1] -= rhs.data[1];
    data[2] -= rhs.data[2];
    data[3] -= rhs.data[3];
#endif
    return *this;
  }

  inline vtype<double, 4> & vtype<double, 4>::operator*=(vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_mul_pd(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] *= rhs.data[0];
    data[1] *= rhs.data[1];
    data[2] *= rhs.data[2];
    data[3] *= rhs.data[3];
#endif
    return *this;
  }

  inline vtype<double, 4> & vtype<double, 4>::operator/=(vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    avxReg = _mm256_div_pd(avxReg, rhs.avxReg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    data[0] /= rhs.data[0];
    data[1] /= rhs.data[1];
    data[2] /= rhs.data[2];
    data[3] /= rhs.data[3];
#endif
    return *this;
  }

  inline vtype<double, 4> operator+(vtype<double, 4> const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_add_pd(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(
                            lhs.data[0] + rhs.data[0],
                            lhs.data[1] + rhs.data[1],
                            lhs.data[2] + rhs.data[2],
                            lhs.data[3] + rhs.data[3]
                            );
#endif
  }

  inline vtype<double,4> operator-(vtype<double,4> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256d sign_bit = _mm256_set1_pd(-0.0) ;
    return vtype<double,4>(_mm256_xor_pd(sign_bit,v.avxReg)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(-v.data[0],-v.data[1],-v.data[2],-v.data[3]) ;
#endif
  }

  inline vtype<double, 4> operator-(vtype<double, 4> const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_sub_pd(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(
                            lhs.data[0] - rhs.data[0],
                            lhs.data[1] - rhs.data[1],
                            lhs.data[2] - rhs.data[2],
                            lhs.data[3] - rhs.data[3]
                            );
#endif
  }

  inline vtype<double, 4> operator*(vtype<double, 4> const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_mul_pd(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(
                            lhs.data[0] * rhs.data[0],
                            lhs.data[1] * rhs.data[1],
                            lhs.data[2] * rhs.data[2],
                            lhs.data[3] * rhs.data[3]
                            );
#endif
  }

  inline vtype<double, 4> operator/(vtype<double, 4> const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_div_pd(lhs.avxReg, rhs.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(
                            lhs.data[0] / rhs.data[0],
                            lhs.data[1] / rhs.data[1],
                            lhs.data[2] / rhs.data[2],
                            lhs.data[3] / rhs.data[3]
                            );
#endif
  }

  inline vtype<double,4> fabs(const vtype<double,4> &v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256d sign_bit = _mm256_set1_pd(-0.0) ;
    return vtype<double,4>(_mm256_andnot_pd(sign_bit,v.avxReg)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(std::fabs(v.data[0]),
                            std::fabs(v.data[1]),
                            std::fabs(v.data[2]),
                            std::fabs(v.data[3])) ;

#endif
  }
  inline vtype<double,4> abs(const vtype<double,4> &v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256d sign_bit = _mm256_set1_pd(-0.0) ;
    return vtype<double,4>(_mm256_andnot_pd(sign_bit,v.avxReg)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(std::fabs(v.data[0]),
                            std::fabs(v.data[1]),
                            std::fabs(v.data[2]),
                            std::fabs(v.data[3])) ;

#endif
  }

  inline vtype<double,4> max(const vtype<double,4> &m1, const vtype<double,4> &m2) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double,4>(_mm256_max_pd(m1.avxReg,m2.avxReg)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(std::max(m1.data[0],m2.data[0]),
                            std::max(m1.data[1],m2.data[1]),
                            std::max(m1.data[2],m2.data[2]),
                            std::max(m1.data[3],m2.data[3])) ;
#endif
  }
  inline vtype<double,4> min(const vtype<double,4> &m1, const vtype<double,4> &m2) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double,4>(_mm256_min_pd(m1.avxReg,m2.avxReg)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(std::min(m1.data[0],m2.data[0]),
                            std::min(m1.data[1],m2.data[1]),
                            std::min(m1.data[2],m2.data[2]),
                            std::min(m1.data[3],m2.data[3]));
#endif

  }

  // boolean comparison
  inline vtype<double, 4> operator<(vtype<double,4>  const &d1, vtype<double,4> const &d2) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double,4>(_mm256_cmp_pd(d1.avxReg,d2.avxReg,_CMP_LT_OQ)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(d1.data[0]<d2.data[0],
                            d1.data[1]<d2.data[1],
                            d1.data[2]<d2.data[2],
                            d1.data[3]<d2.data[3]) ;

#endif
  }

  inline vtype<double, 4> operator>(vtype<double,4>  const &d1, vtype<double,4> const &d2) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double,4>(_mm256_cmp_pd(d1.avxReg,d2.avxReg,_CMP_GT_OQ)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(d1.data[0]<d2.data[0],
                            d1.data[1]<d2.data[1],
                            d1.data[2]<d2.data[2],
                            d1.data[3]<d2.data[3]) ;

#endif
  }

  inline vtype<double, 4> operator==(vtype<double,4>  const &d1, vtype<double,4> const &d2) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double,4>(_mm256_cmp_pd(d1.avxReg,d2.avxReg,_CMP_EQ_OQ)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(d1.data[0]==d2.data[0],
                            d1.data[1]==d2.data[1],
                            d1.data[2]==d2.data[2],
                            d1.data[3]==d2.data[3]) ;

#endif
  }

  inline vtype<double, 4> select(vtype<double,4> const c, vtype<double,4>  const &d1, vtype<double,4> const &d2) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double,4>(_mm256_blendv_pd(d2.avxReg,d1.avxReg,c.avxReg)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>((c.data[0]!=0.)?d1.data[0]:d2.data[0],
                            (c.data[1]!=0.)?d1.data[1]:d2.data[1],
                            (c.data[2]!=0.)?d1.data[2]:d2.data[2],
                            (c.data[3]!=0.)?d1.data[3]:d2.data[3]) ;
#endif
  }

  inline vtype<double, 4> ceil(vtype<double, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_ceil_pd(arg.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(
                            std::ceil(arg.data[0]),
                            std::ceil(arg.data[1]),
                            std::ceil(arg.data[2]),
                            std::ceil(arg.data[3])
                            );
#endif
  }

  inline vtype<double, 4> floor(vtype<double, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_floor_pd(arg.avxReg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(
                            std::floor(arg.data[0]),
                            std::floor(arg.data[1]),
                            std::floor(arg.data[2]),
                            std::floor(arg.data[3])
                            );
#endif
  }

  inline vtype<double, 4> round(vtype<double, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_round_pd(arg.avxReg, _MM_FROUND_TO_NEAREST_INT| _MM_FROUND_NO_EXC));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(
                            std::round(arg.data[0]),
                            std::round(arg.data[1]),
                            std::round(arg.data[2]),
                            std::round(arg.data[3])
                            );
#endif
  }
  inline vtype<double, 4> sqrt(vtype<double, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX__)
    return vtype<double, 4>(
                            _mm256_sqrt_pd(arg.avxReg)
                            );
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(
                            std::sqrt(arg.data[0]),
                            std::sqrt(arg.data[1]),
                            std::sqrt(arg.data[2]),
                            std::sqrt(arg.data[3])
                            );
#endif
  }
}


#endif
