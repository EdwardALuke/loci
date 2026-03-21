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

#if defined(__AVX2__) || defined(__AVX512F__)
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
#include <type_traits>
#include <initializer_list>


// TODO
// 1. What should be the order of registers in a multi-register member of the union in a vtype?
// 2. Is it useful to make vtype types as literals? (requires c++11)
// 3. Check the assembly code generated from vtype expressions.
// 4. Write unit test suites, comparing SIMD results with scalar results.
// 5. Write performance tests (latency and throughput) for vtype operations.
// 6. Implementent transcendental functions. Potential solutions:
//    a. Sleef (potential licensing conflicts, easy of integration in Loci etc.), portable
//    b. Arm Performance Libraries (ArmPL) for ARM processors
//    c. SVML for Intel compiler on x86 processors
//    d. Roll your own implementation
// 7. Extend to use NEON and SVE instructions for Arm CPUs
//    a. Arm C Language Extension (ACLE) provides <arm_neon.h> and <arm_sve.h>

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

    template<typename... Args>
    StaticSizedArray(Args... args) : data{args...} {
      static_assert(
        (std::is_same_v<Args, T> && ...),
        "all arguments must be of the same type as StaticSizedArray"
      );
      static_assert(sizeof...(Args) == N, "number of arguments must be the same size as StaticSizedArray");
    }

    StaticSizedArray(std::initializer_list<T> l) {
      auto iter = l.begin();
      for(int i = 0; i < N && iter != l.end(); ++i) {
        data[i] = *iter++;
      }
    }

    StaticSizedArray(T const * mem) {
      for(int i = 0; i < N; ++i) {
        data[i] = mem[i];
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

  template<int B, int W>
  struct vmask {
  };

  /**
   * Vectorized mask for 32-bit 4x values
   */
  template<>
  struct vmask<4, 4> {
    static constexpr unsigned int W = 4;

    typedef uint8_t bitmask_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    typedef uint8_t mask_t;
    typedef __mmask8 arch_mask_t;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef StaticSizedArray<uint32_t, 4> mask_t;
    typedef __m128i arch_mask_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    typedef uint8_t mask_t;
    typedef uint8_t arch_mask_t;
#endif

    union {
      mask_t data;
      arch_mask_t reg;
    };

    vmask() = default;

    explicit vmask(bitmask_t m)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      : reg(_cvtu32_mask8(m)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(m) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      {
        __m128i shifts = _mm_set_epi32(28, 29, 30, 31);
        __m128i mask = _mm_set1_epi32(m);
        __m128i shifted_mask = _mm_sllv_epi32(mask, shifts);
        reg = _mm_cmpgt_epi32(_mm_setzero_si128(), shifted_mask);
      }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(m) { }
#endif

    explicit vmask(bool e0, bool e1, bool e2, bool e3)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      : reg(_cvtu32_mask8((uint32_t)e0 << 0 |
                          (uint32_t)e1 << 1 |
                          (uint32_t)e2 << 2 |
                          (uint32_t)e3 << 3)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(uint8_t((uint32_t)e0 << 0 |
                    (uint32_t)e1 << 1 |
                    (uint32_t)e2 << 2 |
                    (uint32_t)e3 << 3)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm_set_epi32(-(uint32_t)e3,
                          -(uint32_t)e2,
                          -(uint32_t)e1,
                          -(uint32_t)e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data((uint32_t)e0 << 0 |
             (uint32_t)e1 << 1 |
             (uint32_t)e2 << 2 |
             (uint32_t)e3 << 3) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__)
    explicit vmask(__m128i reg) : reg(reg) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__)
    explicit vmask(__m128 reg) : reg(_mm_castps_si128(reg)) { }
#endif

    vmask(vmask<4, 4> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(rhs.reg) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(rhs.data) { }
#endif

    vmask<4, 4> & operator=(vmask<4, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg = rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data = rhs.data;
#endif
      return *this;
    }

    explicit operator bitmask_t() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return data;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return _mm_movemask_ps(_mm_castsi128_ps(reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return data;
#endif
    }

    bool operator[](unsigned int index) const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return (data >> index) & 1u;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return data[index] != 0;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return (data >> index) & 1u;
#endif
    }

    vmask<4, 4> & operator|=(vmask<4, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kor_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg |= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = _mm_or_si128(reg, rhs.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data |= rhs.data;
#endif
      return *this;
    }

    vmask<4, 4> & operator&=(vmask<4, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kand_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg &= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = _mm_and_si128(reg, rhs.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data &= rhs.data;
#endif
      return *this;
    }

    vmask<4, 4> & operator^=(vmask<4, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kxor_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg ^= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = _mm_xor_si128(reg, rhs.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data ^= rhs.data;
#endif
      return *this;
    }

    vmask<4, 4> operator~() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<4, 4>(_knot_mask8(reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<4, 4>(~reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return vmask<4, 4>(_mm_xor_si128(reg, _mm_set1_epi32(0xffffffff)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<4, 4>(~data);
#endif
    }

    vmask<4, 4> operator!() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<4, 4>(_knot_mask8(reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<4, 4>(~reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return vmask<4, 4>(_mm_xor_si128(reg, _mm_set1_epi32(0xffffffff)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<4, 4>(~data);
#endif
    }
  };

  inline vmask<4, 4> operator|(vmask<4, 4> const & lhs, vmask<4, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 4>(_kor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 4>(lhs.reg | rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 4>(_mm_or_si128(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 4>(lhs.data | rhs.data);
#endif
  }

  inline vmask<4, 4> operator||(vmask<4, 4> const & lhs, vmask<4, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 4>(_kor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 4>(lhs.reg | rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 4>(_mm_or_si128(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 4>(lhs.data | rhs.data);
#endif
  }

  inline vmask<4, 4> operator&(vmask<4, 4> const & lhs, vmask<4, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 4>(_kand_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 4>(lhs.reg & rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 4>(_mm_and_si128(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 4>(lhs.data & rhs.data);
#endif
  }

  inline vmask<4, 4> operator&&(vmask<4, 4> const & lhs, vmask<4, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 4>(_kand_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 4>(lhs.reg & rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 4>(_mm_and_si128(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 4>(lhs.data & rhs.data);
#endif
  }

  inline vmask<4, 4> operator^(vmask<4, 4> const & lhs, vmask<4, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 4>(_kxor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 4>(lhs.reg ^ rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 4>(_mm_xor_si128(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 4>(lhs.data ^ rhs.data);
#endif
  }

  /**
   * Vectorized mask for 32-bit 8x values
   */
  template<>
  struct vmask<4, 8> {
    static constexpr unsigned int W = 8;

    typedef uint8_t bitmask_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    typedef uint8_t mask_t;
    typedef __mmask8 arch_mask_t;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef StaticSizedArray<uint32_t, 8> mask_t;
    typedef __m256i arch_mask_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    typedef uint8_t mask_t;
    typedef uint8_t arch_mask_t;
#endif

    union {
      mask_t data;
      arch_mask_t reg;
    };

    vmask() = default;

    explicit vmask(bitmask_t m)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      : reg(_cvtu32_mask8(m)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(m) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      {
        __m256i shifts = _mm256_set_epi32(24, 25, 26, 27, 28, 29, 30, 31);
        __m256i mask = _mm256_set1_epi32(m);
        __m256i shifted_mask = _mm256_sllv_epi32(mask, shifts);
        reg = _mm256_cmpgt_epi32(_mm256_setzero_si256(), shifted_mask);
      }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(m) { }
#endif

    explicit vmask(bool e0, bool e1, bool e2, bool e3,
                   bool e4, bool e5, bool e6, bool e7)
#if !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512DQ__)
      : reg(_cvtu32_mask8((uint32_t)e0 << 0 |
                          (uint32_t)e1 << 1 |
                          (uint32_t)e2 << 2 |
                          (uint32_t)e3 << 3 |
                          (uint32_t)e4 << 4 |
                          (uint32_t)e5 << 5 |
                          (uint32_t)e6 << 6 |
                          (uint32_t)e7 << 7)) { }
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512F__)
      : reg(uint8_t((uint32_t)e0 << 0 |
                    (uint32_t)e1 << 1 |
                    (uint32_t)e2 << 2 |
                    (uint32_t)e3 << 3 |
                    (uint32_t)e4 << 4 |
                    (uint32_t)e5 << 5 |
                    (uint32_t)e6 << 6 |
                    (uint32_t)e7 << 7)) { }
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX2__)
      : reg(_mm256_set_epi32(-(uint32_t)e7,
                             -(uint32_t)e6,
                             -(uint32_t)e5,
                             -(uint32_t)e4,
                             -(uint32_t)e3,
                             -(uint32_t)e2,
                             -(uint32_t)e1,
                             -(uint32_t)e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data((uint32_t)e0 << 0 |
             (uint32_t)e1 << 1 |
             (uint32_t)e2 << 2 |
             (uint32_t)e3 << 3 |
             (uint32_t)e4 << 4 |
             (uint32_t)e5 << 5 |
             (uint32_t)e6 << 6 |
             (uint32_t)e7 << 7) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__)
    explicit vmask(__m256i reg) : reg(reg) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__) && !defined(__AVX512DQ__)
    explicit vmask(__m256 reg) : reg(_mm256_castps_si256(reg)) { }
#endif

    vmask(vmask<4, 8> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(rhs.reg) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(rhs.data) { }
#endif

    vmask<4, 8> & operator=(vmask<4, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg = rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data = rhs.data;
#endif
      return *this;
    }

    explicit operator bitmask_t() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return data;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return _mm256_movemask_ps(_mm256_castsi256_ps(reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return data;
#endif
    }

    bool operator[](unsigned int index) const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return (data >> index) & 1u;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return data[index] != 0;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return (data >> index) & 1u;
#endif
    }

    vmask<4, 8> & operator|=(vmask<4, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kor_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg |= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = _mm256_or_si256(reg, rhs.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data |= rhs.data;
#endif
      return *this;
    }

    vmask<4, 8> & operator&=(vmask<4, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX5122DQ__)
      reg = _kand_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg &= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = _mm256_and_si256(reg, rhs.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data &= rhs.data;
#endif
      return *this;
    }

    vmask<4, 8> & operator^=(vmask<4, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kxor_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg ^= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = _mm256_xor_si256(reg, rhs.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data ^= rhs.data;
#endif
      return *this;
    }

    vmask<4, 8> operator~() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<4, 8>(_knot_mask8(reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<4, 8>(~reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return vmask<4, 8>(_mm256_xor_si256(reg, _mm256_set1_epi32(0xffffffff)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<4, 8>(~data);
#endif
    }

    vmask<4, 8> operator!() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<4, 8>(_knot_mask8(reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<4, 8>(~reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return vmask<4, 8>(_mm256_xor_si256(reg, _mm256_set1_epi32(0xffffffff)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<4, 8>(~data);
#endif
    }
  };

  inline vmask<4, 8> operator|(vmask<4, 8> const & lhs, vmask<4, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 8>(_kor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 8>(lhs.reg | rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 8>(_mm256_or_si256(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 8>(lhs.data | rhs.data);
#endif
  }

  inline vmask<4, 8> operator||(vmask<4, 8> const & lhs, vmask<4, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 8>(_kor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 8>(lhs.reg | rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 8>(_mm256_or_si256(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 8>(lhs.data | rhs.data);
#endif
  }

  inline vmask<4, 8> operator&(vmask<4, 8> const & lhs, vmask<4, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 8>(_kand_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 8>(lhs.reg & rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 8>(_mm256_and_si256(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 8>(lhs.data & rhs.data);
#endif
  }

  inline vmask<4, 8> operator&&(vmask<4, 8> const & lhs, vmask<4, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 8>(_kand_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 8>(lhs.reg & rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 8>(_mm256_and_si256(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 8>(lhs.data & rhs.data);
#endif
  }

  inline vmask<4, 8> operator^(vmask<4, 8> const & lhs, vmask<4, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 8>(_kxor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 8>(lhs.reg ^ rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 8>(_mm256_xor_si256(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 8>(lhs.data ^ rhs.data);
#endif
  }

  template<>
  struct vmask<4, 16> {
    static constexpr unsigned int W = 16;

    typedef uint16_t bitmask_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    typedef uint16_t mask_t;
    typedef __mmask16 arch_mask_t;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef StaticSizedArray<uint32_t, 16> mask_t;
    typedef struct {
      __m256i x, y;
    } arch_mask_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    typedef uint16_t mask_t;
    typedef uint16_t arch_mask_t;
#endif

    union {
      mask_t data;
      arch_mask_t reg;
    };

    vmask() = default;

    explicit vmask(bitmask_t m)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(_cvtu32_mask16(m)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      {
        __m256i shifts1 = _mm256_set_epi32(24, 25, 26, 27, 28, 29, 30, 31);
        __m256i shifts2 = _mm256_set_epi32(16, 17, 18, 19, 20, 21, 22, 23);
        __m256i mask = _mm256_set1_epi32(m);
        __m256i shifted_mask1 = _mm256_sllv_epi32(mask, shifts1);
        __m256i shifted_mask2 = _mm256_sllv_epi32(mask, shifts2);
        __m256i zero = _mm256_setzero_si256();
        reg.x = _mm256_cmpgt_epi32(zero, shifted_mask1);
        reg.y = _mm256_cmpgt_epi32(zero, shifted_mask2);
      }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(m) { }
#endif

    explicit vmask(bool e0, bool e1, bool e2, bool e3,
                   bool e4, bool e5, bool e6, bool e7,
                   bool e8, bool e9, bool e10, bool e11,
                   bool e12, bool e13, bool e14, bool e15)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(_cvtu32_mask16((uint32_t)e0 << 0 |
                           (uint32_t)e1 << 1 |
                           (uint32_t)e2 << 2 |
                           (uint32_t)e3 << 3 |
                           (uint32_t)e4 << 4 |
                           (uint32_t)e5 << 5 |
                           (uint32_t)e6 << 6 |
                           (uint32_t)e7 << 7 |
                           (uint32_t)e8 << 8 |
                           (uint32_t)e9 << 9 |
                           (uint32_t)e10 << 10 |
                           (uint32_t)e11 << 11 |
                           (uint32_t)e12 << 12 |
                           (uint32_t)e13 << 13 |
                           (uint32_t)e14 << 14 |
                           (uint32_t)e15 << 15)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_set_epi32(-(uint32_t)e7,
                             -(uint32_t)e6,
                             -(uint32_t)e5,
                             -(uint32_t)e4,
                             -(uint32_t)e3,
                             -(uint32_t)e2,
                             -(uint32_t)e1,
                             -(uint32_t)e0),
            _mm256_set_epi32(-(uint32_t)e15,
                             -(uint32_t)e14,
                             -(uint32_t)e13,
                             -(uint32_t)e12,
                             -(uint32_t)e11,
                             -(uint32_t)e10,
                             -(uint32_t)e9,
                             -(uint32_t)e8)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data((uint32_t)e0 << 0 |
             (uint32_t)e1 << 1 |
             (uint32_t)e2 << 2 |
             (uint32_t)e3 << 3 |
             (uint32_t)e4 << 4 |
             (uint32_t)e5 << 5 |
             (uint32_t)e6 << 6 |
             (uint32_t)e7 << 7 |
             (uint32_t)e8 << 8 |
             (uint32_t)e9 << 9 |
             (uint32_t)e10 << 10 |
             (uint32_t)e11 << 11 |
             (uint32_t)e12 << 12 |
             (uint32_t)e13 << 13 |
             (uint32_t)e14 << 14 |
             (uint32_t)e15 << 15) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__)
    explicit vmask(__m256i reg0, __m256i reg1) : reg{reg0, reg1} { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__)
    explicit vmask(__m256 reg0, __m256 reg1)
      : reg{_mm256_castps_si256(reg0), _mm256_castps_si256(reg1)} { }
#endif

    vmask(vmask<4, 16> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(rhs.reg) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(rhs.data) { }
#endif

    vmask<4, 16> & operator=(vmask<4, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg = rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data = rhs.data;
#endif
      return *this;
    }

    explicit operator bitmask_t() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return data;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      uint16_t m0 = _mm256_movemask_ps(_mm256_castsi256_ps(reg.x));
      uint16_t m1 = _mm256_movemask_ps(_mm256_castsi256_ps(reg.y));
      return m0 | (m1 << 8);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return data;
#endif
    }

    bool operator[](unsigned int index) const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return (data >> index) & 1u;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return data[index] != 0;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return (data >> index) & 1u;
#endif
    }

    vmask<4, 16> & operator|=(vmask<4, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kor_mask16(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg |= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg.x = _mm256_or_si256(reg.x, rhs.reg.x);
      reg.y = _mm256_or_si256(reg.y, rhs.reg.y);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data |= rhs.data;
#endif
      return *this;
    }

    vmask<4, 16> & operator&=(vmask<4, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kand_mask16(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg &= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg.x = _mm256_and_si256(reg.x, rhs.reg.x);
      reg.y = _mm256_and_si256(reg.y, rhs.reg.y);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data &= rhs.data;
#endif
      return *this;
    }

    vmask<4, 16> & operator^=(vmask<4, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kxor_mask16(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg ^= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg.x = _mm256_xor_si256(reg.x, rhs.reg.x);
      reg.y = _mm256_xor_si256(reg.y, rhs.reg.y);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data ^= rhs.data;
#endif
      return *this;
    }

    vmask<4, 16> operator~() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<4, 16>(_knot_mask16(reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<4, 16>(~reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      __m256i ones = _mm256_set1_epi32(0xffffffff);
      return vmask<4, 16>(_mm256_xor_si256(reg.x, ones),
                          _mm256_xor_si256(reg.y, ones));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<4, 16>(~data);
#endif
    }

    vmask<4, 16> operator!() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<4, 16>(_knot_mask16(reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<4, 16>(~reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      __m256i ones = _mm256_set1_epi32(0xffffffff);
      return vmask<4, 16>(_mm256_xor_si256(reg.x, ones),
                          _mm256_xor_si256(reg.y, ones));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<4, 16>(~data);
#endif
    }
  };

  inline vmask<4, 16> operator|(vmask<4, 16> const & lhs, vmask<4, 16> const & rhs) {
#if !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512DQ__)
    return vmask<4, 16>(_kor_mask16(lhs.reg, rhs.reg));
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512F__)
    return vmask<4, 16>(lhs.reg | rhs.reg);
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX2__)
    return vmask<4, 16>(_mm256_or_si256(lhs.reg.x, rhs.reg.x),
                        _mm256_or_si256(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 16>(lhs.data | rhs.data);
#endif
  }

  inline vmask<4, 16> operator||(vmask<4, 16> const & lhs, vmask<4, 16> const & rhs) {
#if !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512DQ__)
    return vmask<4, 16>(_kor_mask16(lhs.reg, rhs.reg));
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512F__)
    return vmask<4, 16>(lhs.reg | rhs.reg);
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX2__)
    return vmask<4, 16>(_mm256_or_si256(lhs.reg.x, rhs.reg.x),
                        _mm256_or_si256(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 16>(lhs.data | rhs.data);
#endif
  }

  inline vmask<4, 16> operator&(vmask<4, 16> const & lhs, vmask<4, 16> const & rhs) {
#if !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512DQ__)
    return vmask<4, 16>(_kand_mask16(lhs.reg, rhs.reg));
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512F__)
    return vmask<4, 16>(lhs.reg & rhs.reg);
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX2__)
    return vmask<4, 16>(_mm256_and_si256(lhs.reg.x, rhs.reg.x),
                        _mm256_and_si256(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 16>(lhs.data & rhs.data);
#endif
  }

  inline vmask<4, 16> operator&&(vmask<4, 16> const & lhs, vmask<4, 16> const & rhs) {
#if !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512DQ__)
    return vmask<4, 16>(_kand_mask16(lhs.reg, rhs.reg));
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512F__)
    return vmask<4, 16>(lhs.reg & rhs.reg);
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX2__)
    return vmask<4, 16>(_mm256_and_si256(lhs.reg.x, rhs.reg.x),
                        _mm256_and_si256(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 16>(lhs.data & rhs.data);
#endif
  }

  inline vmask<4, 16> operator^(vmask<4, 16> const & lhs, vmask<4, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<4, 16>(_kxor_mask16(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<4, 16>(lhs.reg ^ rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<4, 16>(_mm256_xor_si256(lhs.reg.x, rhs.reg.x),
                        _mm256_xor_si256(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<4, 16>(lhs.data ^ rhs.data);
#endif
  }

  /**
   * Vectorized mask for 64-bit 4x values
   */
  template<>
  struct vmask<8, 4> {
    static constexpr unsigned int W = 4;

    typedef uint8_t bitmask_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    typedef uint8_t mask_t;
    typedef __mmask8 arch_mask_t;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef StaticSizedArray<uint64_t, 4> mask_t;
    typedef __m256i arch_mask_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    typedef uint8_t mask_t;
    typedef uint8_t arch_mask_t;
#endif

    union {
      mask_t data;
      arch_mask_t reg;
    };

    vmask() = default;

    explicit vmask(bitmask_t m)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      : reg(_cvtu32_mask8(m)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(m) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      {
        __m256i shifts = _mm256_set_epi64x(60, 61, 62, 63);
        __m256i mask = _mm256_set1_epi64x(m);
        __m256i shifted_mask = _mm256_sllv_epi64(mask, shifts);
        reg = _mm256_cmpgt_epi64(_mm256_setzero_si256(), shifted_mask);
      }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(m) { }
#endif

    explicit vmask(bool e0, bool e1, bool e2, bool e3)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      : reg(_cvtu32_mask8((uint8_t)e0 << 0 |
                          (uint8_t)e1 << 1 |
                          (uint8_t)e2 << 2 |
                          (uint8_t)e3 << 3)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(uint8_t((uint8_t)e0 << 0 |
                    (uint8_t)e1 << 1 |
                    (uint8_t)e2 << 2 |
                    (uint8_t)e3 << 3)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm256_set_epi64x(-(uint64_t)e3,
                              -(uint64_t)e2,
                              -(uint64_t)e1,
                              -(uint64_t)e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data((uint8_t)e0 << 0 |
             (uint8_t)e1 << 1 |
             (uint8_t)e2 << 2 |
             (uint8_t)e3 << 3) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__) && !defined(__AVX512DQ__)
    explicit vmask(__m256i reg) : reg(reg) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__) && !defined(__AVX512DQ__)
    explicit vmask(__m256d reg) : reg(_mm256_castpd_si256(reg)) { }
#endif

    vmask(vmask<8, 4> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(rhs.reg) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(rhs.data) { }
#endif

    vmask<8, 4> & operator=(vmask<8, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg = rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data = rhs.data;
#endif
      return *this;
    }

    explicit operator bitmask_t() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return data;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return _mm256_movemask_pd(_mm256_castsi256_pd(reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return data;
#endif
    }

    bool operator[](unsigned int index) const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return (data >> index) & 1u;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return data[index] != 0;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return (data >> index) & 1u;
#endif
    }

    vmask<8, 4> & operator|=(vmask<8, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kor_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg |= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = _mm256_or_si256(reg, rhs.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data |= rhs.data;
#endif
      return *this;
    }

    vmask<8, 4> & operator&=(vmask<8, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kand_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg &= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = _mm256_and_si256(reg, rhs.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data &= rhs.data;
#endif
      return *this;
    }

    vmask<8, 4> & operator^=(vmask<8, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kxor_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg ^= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = _mm256_xor_si256(reg, rhs.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data ^= rhs.data;
#endif
      return *this;
    }

    vmask<8, 4> operator~() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<8, 4>(_knot_mask8(reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<8, 4>(~reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return vmask<8, 4>(_mm256_xor_si256(reg, _mm256_set1_epi64x(0xffffffffffffffff)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<8, 4>(~data);
#endif
    }

    vmask<8, 4> operator!() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<8, 4>(_knot_mask8(reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<8, 4>(~reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return vmask<8, 4>(_mm256_xor_si256(reg, _mm256_set1_epi64x(0xffffffffffffffff)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<8, 4>(~data);
#endif
    }
  };

  inline vmask<8, 4> operator|(vmask<8, 4> const & lhs, vmask<8, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<8, 4>(_kor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<8, 4>(lhs.reg | rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<8, 4>(_mm256_or_si256(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 4>(lhs.data | rhs.data);
#endif
  }

  inline vmask<8, 4> operator||(vmask<8, 4> const & lhs, vmask<8, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<8, 4>(_kor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<8, 4>(lhs.reg | rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<8, 4>(_mm256_or_si256(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 4>(lhs.data | rhs.data);
#endif
  }

  inline vmask<8, 4> operator&(vmask<8, 4> const & lhs, vmask<8, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<8, 4>(_kand_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<8, 4>(lhs.reg & rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<8, 4>(_mm256_and_si256(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 4>(lhs.data & rhs.data);
#endif
  }

  inline vmask<8, 4> operator&&(vmask<8, 4> const & lhs, vmask<8, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<8, 4>(_kand_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<8, 4>(lhs.reg & rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<8, 4>(_mm256_and_si256(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 4>(lhs.data & rhs.data);
#endif
  }

  inline vmask<8, 4> operator^(vmask<8, 4> const & lhs, vmask<8, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<8, 4>(_kxor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<8, 4>(lhs.reg ^ rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<8, 4>(_mm256_xor_si256(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 4>(lhs.data ^ rhs.data);
#endif
  }

  /**
   * Vectorized mask for 64-bit 8x values
   */
  template<>
  struct vmask<8, 8> {
    static constexpr unsigned int W = 8;

    typedef uint8_t bitmask_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    typedef uint8_t mask_t;
    typedef __mmask8 arch_mask_t;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef StaticSizedArray<uint64_t, 8> mask_t;
    typedef struct {
      __m256i x, y;
    } arch_mask_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    typedef uint8_t mask_t;
    typedef uint8_t arch_mask_t;
#endif

    union {
      mask_t data;
      arch_mask_t reg;
    };

    vmask() = default;

    vmask(bool e0, bool e1, bool e2, bool e3,
          bool e4, bool e5, bool e6, bool e7)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      : reg(_cvtu32_mask8((uint32_t)e0 << 0 |
                          (uint32_t)e1 << 1 |
                          (uint32_t)e2 << 2 |
                          (uint32_t)e3 << 3 |
                          (uint32_t)e4 << 4 |
                          (uint32_t)e5 << 5 |
                          (uint32_t)e6 << 6 |
                          (uint32_t)e7 << 7)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(uint8_t((uint32_t)e0 << 0 |
                    (uint32_t)e1 << 1 |
                    (uint32_t)e2 << 2 |
                    (uint32_t)e3 << 3 |
                    (uint32_t)e4 << 4 |
                    (uint32_t)e5 << 5 |
                    (uint32_t)e6 << 6 |
                    (uint32_t)e7 << 7)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_set_epi64x(-(uint64_t)e3,
                              -(uint64_t)e2,
                              -(uint64_t)e1,
                              -(uint64_t)e0),
            _mm256_set_epi64x(-(uint64_t)e7,
                              -(uint64_t)e6,
                              -(uint64_t)e5,
                              -(uint64_t)e4)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data((uint32_t)e0 << 0 |
             (uint32_t)e1 << 1 |
             (uint32_t)e2 << 2 |
             (uint32_t)e3 << 3 |
             (uint32_t)e4 << 4 |
             (uint32_t)e5 << 5 |
             (uint32_t)e6 << 6 |
             (uint32_t)e7 << 7) { }
#endif

    explicit vmask(bitmask_t m)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      : reg(_cvtu32_mask8(m)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(m) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      {
        __m256i shifts1 = _mm256_set_epi64x(60, 61, 62, 63);
        __m256i shifts2 = _mm256_set_epi64x(56, 57, 58, 59);
        __m256i mask = _mm256_set1_epi64x(m);
        __m256i shifted_mask1 = _mm256_sllv_epi64(mask, shifts1);
        __m256i shifted_mask2 = _mm256_sllv_epi64(mask, shifts2);
        __m256i zero = _mm256_setzero_si256();
        reg.x = _mm256_cmpgt_epi64(zero, shifted_mask1);
        reg.y = _mm256_cmpgt_epi64(zero, shifted_mask2);
      }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(m) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__) && !defined(__AVX512DQ__)
    explicit vmask(__m256i reg0, __m256i reg1) : reg{reg0, reg1} { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__) && !defined(__AVX512DQ__)
    explicit vmask(__m256d reg0, __m256d reg1) : reg{_mm256_castpd_si256(reg0), _mm256_castpd_si256(reg1)} { }
#endif

    vmask(vmask<8, 8> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(rhs.reg) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(rhs.data) { }
#endif

    vmask<8, 8> & operator=(vmask<8, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg = rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data = rhs.data;
#endif
      return *this;
    }

    explicit operator bitmask_t() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return data;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      bitmask_t m0 = _mm256_movemask_pd(_mm256_castsi256_pd(reg.x));
      bitmask_t m1 = _mm256_movemask_pd(_mm256_castsi256_pd(reg.y));
      return m0 | (m1 << 4);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return data;
#endif
    }

    bool operator[](unsigned int index) const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return (data >> index) & 1u;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return data[index] != 0;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return (data >> index) & 1u;
#endif
    }

    vmask<8, 8> & operator|=(vmask<8, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kor_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg |= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg.x = _mm256_or_si256(reg.x, rhs.reg.x);
      reg.y = _mm256_or_si256(reg.y, rhs.reg.y);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data |= rhs.data;
#endif
      return *this;
    }

    vmask<8, 8> & operator&=(vmask<8, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kand_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg &= rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg.x = _mm256_and_si256(reg.x, rhs.reg.x);
      reg.y = _mm256_and_si256(reg.y, rhs.reg.y);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data &= rhs.data;
#endif
      return *this;
    }

    vmask<8, 8> & operator^=(vmask<8, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg = _kxor_mask8(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg ^= rhs.reg;
#elif!defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg.x = _mm256_xor_si256(reg.x, rhs.reg.x);
      reg.y = _mm256_xor_si256(reg.y, rhs.reg.y);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data ^= rhs.data;
#endif
      return *this;
    }

    vmask<8, 8> operator~() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<8, 8>(_knot_mask8(reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<8, 8>(~reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      __m256i allone = _mm256_set1_epi64x(0xffffffffffffffff);
      return vmask<8, 8>(_mm256_xor_si256(reg.x, allone),
                         _mm256_xor_si256(reg.y, allone));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<8, 8>(~data);
#endif
    }

    vmask<8, 8> operator!() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<8, 8>(_knot_mask8(reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<8, 8>(~reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      __m256i allone = _mm256_set1_epi64x(0xffffffffffffffff);
      return vmask<8, 8>(_mm256_xor_si256(reg.x, allone),
                         _mm256_xor_si256(reg.y, allone));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<8, 8>(~data);
#endif
    }
  };

  inline vmask<8, 8> operator|(vmask<8, 8> const & lhs, vmask<8, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<8, 8>(_kor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<8, 8>(lhs.reg | rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<8, 8>(_mm256_or_si256(lhs.reg.x, rhs.reg.x),
                       _mm256_or_si256(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 8>(lhs.data | rhs.data);
#endif
  }

  inline vmask<8, 8> operator||(vmask<8, 8> const & lhs, vmask<8, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<8, 8>(_kor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<8, 8>(lhs.reg | rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<8, 8>(_mm256_or_si256(lhs.reg.x, rhs.reg.x),
                       _mm256_or_si256(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 8>(lhs.data | rhs.data);
#endif
  }

  inline vmask<8, 8> operator&(vmask<8, 8> const & lhs, vmask<8, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<8, 8>(_kand_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<8, 8>(lhs.reg & rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<8, 8>(_mm256_and_si256(lhs.reg.x, rhs.reg.x),
                       _mm256_and_si256(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 8>(lhs.data & rhs.data);
#endif
  }

  inline vmask<8, 8> operator&&(vmask<8, 8> const & lhs, vmask<8, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<8, 8>(_kand_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<8, 8>(lhs.reg & rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<8, 8>(_mm256_and_si256(lhs.reg.x, rhs.reg.x),
                       _mm256_and_si256(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 8>(lhs.data & rhs.data);
#endif
  }

  inline vmask<8, 8> operator^(vmask<8, 8> const & lhs, vmask<8, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    return vmask<8, 8>(_kxor_mask8(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<8, 8>(lhs.reg ^ rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<8, 8>(_mm256_xor_si256(lhs.reg.x, rhs.reg.x),
                       _mm256_xor_si256(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 8>(lhs.data ^ rhs.data);
#endif
  }

  template<>
  struct vmask<8, 16> {
    static constexpr unsigned int W = 16;

    typedef uint16_t bitmask_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    typedef uint16_t mask_t;
    typedef struct {
      __mmask8 x, y;
    } arch_mask_t;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef StaticSizedArray<uint64_t, 16> mask_t;
    typedef struct {
      __m256i x, y, z, w;
    } arch_mask_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    typedef uint16_t mask_t;
    typedef uint16_t arch_mask_t;
#endif

    union {
      mask_t data;
      arch_mask_t reg;
    };

    vmask() = default;

    explicit vmask(bitmask_t m)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512FDQ__)
      : reg{_cvtu32_mask8(m), _cvtu32_mask8(m >> 8)} { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg{(__mmask8)(m & 0xff), (__mmask8)((m >> 8) & 0xff)} { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      {
        __m256i shifts1 = _mm256_set_epi64x(60, 61, 62, 63);
        __m256i shifts2 = _mm256_set_epi64x(56, 57, 58, 59);
        __m256i shifts3 = _mm256_set_epi64x(52, 53, 54, 55);
        __m256i shifts4 = _mm256_set_epi64x(48, 49, 50, 51);
        __m256i mask = _mm256_set1_epi64x(m);
        __m256i shifted_mask1 = _mm256_sllv_epi64(mask, shifts1);
        __m256i shifted_mask2 = _mm256_sllv_epi64(mask, shifts2);
        __m256i shifted_mask3 = _mm256_sllv_epi64(mask, shifts3);
        __m256i shifted_mask4 = _mm256_sllv_epi64(mask, shifts4);
        __m256i zero = _mm256_setzero_si256();
        reg.x = _mm256_cmpgt_epi64(zero, shifted_mask1);
        reg.y = _mm256_cmpgt_epi64(zero, shifted_mask2);
        reg.z = _mm256_cmpgt_epi64(zero, shifted_mask3);
        reg.w = _mm256_cmpgt_epi64(zero, shifted_mask4);
      }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(m) { }
#endif

    explicit vmask(bool e0, bool e1, bool e2, bool e3,
                   bool e4, bool e5, bool e6, bool e7,
                   bool e8, bool e9, bool e10, bool e11,
                   bool e12, bool e13, bool e14, bool e15)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      : reg{_cvtu32_mask8((uint32_t)e0 << 0 |
                          (uint32_t)e1 << 1 |
                          (uint32_t)e2 << 2 |
                          (uint32_t)e3 << 3 |
                          (uint32_t)e4 << 4 |
                          (uint32_t)e5 << 5 |
                          (uint32_t)e6 << 6 |
                          (uint32_t)e7 << 7),
            _cvtu32_mask8((uint32_t)e8 << 0 |
                          (uint32_t)e9 << 1 |
                          (uint32_t)e10 << 2 |
                          (uint32_t)e11 << 3 |
                          (uint32_t)e12 << 4 |
                          (uint32_t)e13 << 5 |
                          (uint32_t)e14 << 6 |
                          (uint32_t)e15 << 7)} { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg{(__mmask8)((uint32_t)e0 << 0 |
                       (uint32_t)e1 << 1 |
                       (uint32_t)e2 << 2 |
                       (uint32_t)e3 << 3 |
                       (uint32_t)e4 << 4 |
                       (uint32_t)e5 << 5 |
                       (uint32_t)e6 << 6 |
                       (uint32_t)e7 << 7),
            (__mmask8)((uint32_t)e8 << 0 |
                       (uint32_t)e9 << 1 |
                       (uint32_t)e10 << 2 |
                       (uint32_t)e11 << 3 |
                       (uint32_t)e12 << 4 |
                       (uint32_t)e13 << 5 |
                       (uint32_t)e14 << 6 |
                       (uint32_t)e15 << 7)} { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_set_epi64x(-(uint64_t)e3,
                              -(uint64_t)e2,
                              -(uint64_t)e1,
                              -(uint64_t)e0),
            _mm256_set_epi64x(-(uint64_t)e7,
                              -(uint64_t)e6,
                              -(uint64_t)e5,
                              -(uint64_t)e4),
            _mm256_set_epi64x(-(uint64_t)e11,
                              -(uint64_t)e10,
                              -(uint64_t)e9,
                              -(uint64_t)e8),
            _mm256_set_epi64x(-(uint64_t)e15,
                              -(uint64_t)e14,
                              -(uint64_t)e13,
                              -(uint64_t)e12)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data((uint32_t)e0 << 0 |
             (uint32_t)e1 << 1 |
             (uint32_t)e2 << 2 |
             (uint32_t)e3 << 3 |
             (uint32_t)e4 << 4 |
             (uint32_t)e5 << 5 |
             (uint32_t)e6 << 6 |
             (uint32_t)e7 << 7 |
             (uint32_t)e8 << 8 |
             (uint32_t)e9 << 9 |
             (uint32_t)e10 << 10 |
             (uint32_t)e11 << 11 |
             (uint32_t)e12 << 12 |
             (uint32_t)e13 << 13 |
             (uint32_t)e14 << 14 |
             (uint32_t)e15 << 15) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    explicit vmask(__mmask8 reg0, __mmask8 reg1) : reg{reg0, reg1} { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__)
    explicit vmask(__m256i reg0, __m256i reg1, __m256i reg2, __m256i reg3)
      : reg{reg0, reg1, reg2, reg3} { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__)
    explicit vmask(__m256d reg0, __m256d reg1, __m256d reg2, __m256d reg3)
      : reg{_mm256_castpd_si256(reg0), _mm256_castpd_si256(reg1),
            _mm256_castpd_si256(reg2), _mm256_castpd_si256(reg3)} { }
#endif

    vmask(vmask<8, 16> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(rhs.reg) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(rhs.data) { }
#endif

    vmask<8, 16> & operator=(vmask<8, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg = rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data = rhs.data;
#endif
      return *this;
    }

    explicit operator bitmask_t() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return data;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      unsigned int m0 = _mm256_movemask_pd(_mm256_castsi256_pd(reg.x));
      unsigned int m1 = _mm256_movemask_pd(_mm256_castsi256_pd(reg.y));
      unsigned int m2 = _mm256_movemask_pd(_mm256_castsi256_pd(reg.z));
      unsigned int m3 = _mm256_movemask_pd(_mm256_castsi256_pd(reg.w));
      return m0 | (m1 << 4) | (m2 << 8) | (m3 << 12);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return data;
#endif
    }

    bool operator[](unsigned int index) const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return (data >> index) & 1u;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      return data[index] != 0;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return (data >> index) & 1u;
#endif
    }

    vmask<8, 16> & operator|=(vmask<8, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ)
      reg.x = _kor_mask8(reg.x, rhs.reg.x);
      reg.y = _kor_mask8(reg.y, rhs.reg.y);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg.x |= rhs.reg.x;
      reg.y |= rhs.reg.y;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg.x = _mm256_or_si256(reg.x, rhs.reg.x);
      reg.y = _mm256_or_si256(reg.y, rhs.reg.y);
      reg.z = _mm256_or_si256(reg.z, rhs.reg.z);
      reg.w = _mm256_or_si256(reg.w, rhs.reg.w);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data |= rhs.data;
#endif
      return *this;
    }

    vmask<8, 16> & operator&=(vmask<8, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg.x = _kand_mask8(reg.x, rhs.reg.x);
      reg.y = _kand_mask8(reg.y, rhs.reg.y);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg.x &= rhs.reg.x;
      reg.y &= rhs.reg.y;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg.x = _mm256_and_si256(reg.x, rhs.reg.x);
      reg.y = _mm256_and_si256(reg.y, rhs.reg.y);
      reg.z = _mm256_and_si256(reg.z, rhs.reg.z);
      reg.w = _mm256_and_si256(reg.w, rhs.reg.w);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data &= rhs.data;
#endif
      return *this;
    }

    vmask<8, 16> & operator^=(vmask<8, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      reg.x = _kxor_mask8(reg.x, rhs.reg.x);
      reg.y = _kxor_mask8(reg.y, rhs.reg.y);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg.x ^= rhs.reg.x;
      reg.y ^= rhs.reg.y;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg.x = _mm256_xor_si256(reg.x, rhs.reg.x);
      reg.y = _mm256_xor_si256(reg.y, rhs.reg.y);
      reg.z = _mm256_xor_si256(reg.z, rhs.reg.z);
      reg.w = _mm256_xor_si256(reg.w, rhs.reg.w);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      data ^= rhs.data;
#endif
      return *this;
    }

    vmask<8, 16> operator~() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<8, 16>(_knot_mask8(reg.x),
                          _knot_mask8(reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<8, 16>(~reg.x, ~reg.y);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      __m256i allone = _mm256_set1_epi64x(0xffffffffffffffff);
      return vmask<8, 16>(_mm256_xor_si256(reg.x, allone),
                          _mm256_xor_si256(reg.y, allone),
                          _mm256_xor_si256(reg.z, allone),
                          _mm256_xor_si256(reg.w, allone));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<8, 16>(~data);
#endif
    }

    vmask<8, 16> operator!() const {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
      return vmask<8, 16>(_knot_mask8(reg.x),
                          _knot_mask8(reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      return vmask<8, 16>(~reg.x, ~reg.y);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      __m256i allone = _mm256_set1_epi64x(0xffffffffffffffff);
      return vmask<8, 16>(_mm256_xor_si256(reg.x, allone),
                          _mm256_xor_si256(reg.y, allone),
                          _mm256_xor_si256(reg.z, allone),
                          _mm256_xor_si256(reg.w, allone));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      return vmask<8, 16>(~data);
#endif
    }
  };

  inline vmask<8, 16> operator|(vmask<8, 16> const & lhs, vmask<8, 16> const & rhs) {
#if !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512DQ__)
    return vmask<8, 16>(_kor_mask8(lhs.reg.x, rhs.reg.x),
                        _kor_mask8(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512F__)
    return vmask<8, 16>(lhs.reg.x | rhs.reg.x,
                        lhs.reg.y | rhs.reg.y);
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX2__)
    return vmask<8, 16>(_mm256_or_si256(lhs.reg.x, rhs.reg.x),
                        _mm256_or_si256(lhs.reg.y, rhs.reg.y),
                        _mm256_or_si256(lhs.reg.z, rhs.reg.z),
                        _mm256_or_si256(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 16>(lhs.data | rhs.data);
#endif
  }

  inline vmask<8, 16> operator||(vmask<8, 16> const & lhs, vmask<8, 16> const & rhs) {
#if !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512DQ__)
    return vmask<8, 16>(_kor_mask8(lhs.reg.x, rhs.reg.x),
                        _kor_mask8(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512F__)
    return vmask<8, 16>(lhs.reg.x | rhs.reg.x,
                        lhs.reg.y | rhs.reg.y);
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX2__)
    return vmask<8, 16>(_mm256_or_si256(lhs.reg.x, rhs.reg.x),
                        _mm256_or_si256(lhs.reg.y, rhs.reg.y),
                        _mm256_or_si256(lhs.reg.z, rhs.reg.z),
                        _mm256_or_si256(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 16>(lhs.data | rhs.data);
#endif
  }

  inline vmask<8, 16> operator&(vmask<8, 16> const & lhs, vmask<8, 16> const & rhs) {
#if !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX5122DQ__)
    return vmask<8, 16>(_kand_mask8(lhs.reg.x, rhs.reg.x),
                        _kand_mask8(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512F__)
    return vmask<8, 16>(lhs.reg.x & rhs.reg.x,
                        lhs.reg.y & rhs.reg.y);
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX2__)
    return vmask<8, 16>(_mm256_and_si256(lhs.reg.x, rhs.reg.x),
                        _mm256_and_si256(lhs.reg.y, rhs.reg.y),
                        _mm256_and_si256(lhs.reg.z, rhs.reg.z),
                        _mm256_and_si256(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 16>(lhs.data & rhs.data);
#endif
  }

  inline vmask<8, 16> operator&&(vmask<8, 16> const & lhs, vmask<8, 16> const & rhs) {
#if !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX5122DQ__)
    return vmask<8, 16>(_kand_mask8(lhs.reg.x, rhs.reg.x),
                        _kand_mask8(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512F__)
    return vmask<8, 16>(lhs.reg.x & rhs.reg.x,
                        lhs.reg.y & rhs.reg.y);
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX2__)
    return vmask<8, 16>(_mm256_and_si256(lhs.reg.x, rhs.reg.x),
                        _mm256_and_si256(lhs.reg.y, rhs.reg.y),
                        _mm256_and_si256(lhs.reg.z, rhs.reg.z),
                        _mm256_and_si256(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 16>(lhs.data & rhs.data);
#endif
  }

  inline vmask<8, 16> operator^(vmask<8, 16> const & lhs, vmask<8, 16> const & rhs) {
#if !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512DQ__)
    return vmask<8, 16>(_kxor_mask8(lhs.reg.x, rhs.reg.x),
                        _kxor_mask8(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX512F__)
    return vmask<8, 16>(lhs.reg.x ^ rhs.reg.x,
                        lhs.reg.y ^ rhs.reg.y);
#elif !defined(VTYPE_ENABLE_SCALAR) && defined(__AVX2__)
    return vmask<8, 16>(_mm256_xor_si256(lhs.reg.x, rhs.reg.x),
                        _mm256_xor_si256(lhs.reg.y, rhs.reg.y),
                        _mm256_xor_si256(lhs.reg.z, rhs.reg.z),
                        _mm256_xor_si256(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<8, 16>(lhs.data ^ rhs.data);
#endif
  }

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

    typedef int32_t scalar_t;
    typedef vmask<sizeof(int32_t), 8> mask_t;
    typedef StaticSizedArray<int32_t, 8> value_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef __m256i register_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    typedef StaticSizedArray<int32_t, 8> register_t;
#endif

    union {
      value_t data;
      register_t reg;
    };

    vtype() = default ;

    vtype(int32_t e0, int32_t e1, int32_t e2, int32_t e3,
          int32_t e4, int32_t e5, int32_t e6, int32_t e7)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm256_set_epi32(e7, e6, e5, e4, e3, e2, e1, e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(e0, e1, e2, e3, e4, e5, e6, e7) { }
#endif

    vtype(int32_t e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm256_set1_epi32(e)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(e) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    explicit vtype(__m256i reg) : reg(reg) { }
#endif

    explicit vtype(int32_t const * mem)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
      : reg(_mm256_loadu_epi32(mem)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm256_loadu_si256((__m256i const *)mem)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
      : data(mem) { }
#endif

    vtype(vtype<int32_t, 8> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
      : data(rhs.data) { }
#endif

    vtype<int32_t, 8> & operator=(vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
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
    reg = _mm256_add_epi32(reg, rhs.reg);
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
    reg = _mm256_sub_epi32(reg, rhs.reg);
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
    reg = _mm256_mullo_epi32(reg, rhs.reg);
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
    avx2_si32_si32_quorem(reg, rhs.reg, reg, r);
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
    avx2_si32_si32_quorem(reg, rhs.reg, q, reg);
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
    reg = _mm256_sllv_epi32(reg, rhs.reg);
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
    // TODO: srav is arithmetic shift, that is, fills MSBs with sign bit
    reg = _mm256_srav_epi32(reg, rhs.reg);
    // TODO: srlv is logical shift, that is fills MSBs with zero
    //reg = _mm256_srav_epi32(reg, rhs.reg);
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
    reg = _mm256_sll_epi32(reg, sh);
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
    // TODO: sra is arithmetic shift, that is, fills MSBs with sign bit
    reg = _mm256_sra_epi32(reg, sh);
    // TODO: srl is logical shift, that is fills MSBs with zero
    //reg = _mm256_srl_epi32(reg, sh);
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
    return vtype<int32_t, 8>(_mm256_add_epi32(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<int32_t, 8>(lhs.data[0] + rhs.data[0],
                             lhs.data[1] + rhs.data[1],
                             lhs.data[2] + rhs.data[2],
                             lhs.data[3] + rhs.data[3],
                             lhs.data[4] + rhs.data[4],
                             lhs.data[5] + rhs.data[5],
                             lhs.data[6] + rhs.data[6],
                             lhs.data[7] + rhs.data[7]);
#endif
  }

  inline vtype<int32_t, 8> operator-(vtype<int32_t, 8> const & lhs, vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<int32_t, 8>(_mm256_sub_epi32(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<int32_t, 8>(lhs.data[0] - rhs.data[0],
                             lhs.data[1] - rhs.data[1],
                             lhs.data[2] - rhs.data[2],
                             lhs.data[3] - rhs.data[3],
                             lhs.data[4] - rhs.data[4],
                             lhs.data[5] - rhs.data[5],
                             lhs.data[6] - rhs.data[6],
                             lhs.data[7] - rhs.data[7]);
#endif
  }

  inline vtype<int32_t, 8> operator*(vtype<int32_t, 8> const & lhs, vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<int32_t, 8>(_mm256_mullo_epi32(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<int32_t, 8>(lhs.data[0] * rhs.data[0],
                             lhs.data[1] * rhs.data[1],
                             lhs.data[2] * rhs.data[2],
                             lhs.data[3] * rhs.data[3],
                             lhs.data[4] * rhs.data[4],
                             lhs.data[5] * rhs.data[5],
                             lhs.data[6] * rhs.data[6],
                             lhs.data[7] * rhs.data[7]);
#endif
  }

  inline vtype<int32_t, 8> operator/(vtype<int32_t, 8> const & lhs, vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256i q, r;
    avx2_si32_si32_quorem(lhs.reg, rhs.reg, q, r);
    return vtype<int32_t, 8>(q);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<int32_t, 8>(lhs.data[0] / rhs.data[0],
                             lhs.data[1] / rhs.data[1],
                             lhs.data[2] / rhs.data[2],
                             lhs.data[3] / rhs.data[3],
                             lhs.data[4] / rhs.data[4],
                             lhs.data[5] / rhs.data[5],
                             lhs.data[6] / rhs.data[6],
                             lhs.data[7] / rhs.data[7]);
#endif
  }


  inline vtype<int32_t, 8> operator%(vtype<int32_t, 8> const & lhs, vtype<int32_t, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256i q, r;
    avx2_si32_si32_quorem(lhs.reg, rhs.reg, q, r);
    return vtype<int32_t, 8>(r);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<int32_t, 8>(lhs.data[0] % rhs.data[0],
                             lhs.data[1] % rhs.data[1],
                             lhs.data[2] % rhs.data[2],
                             lhs.data[3] % rhs.data[3],
                             lhs.data[4] % rhs.data[4],
                             lhs.data[5] % rhs.data[5],
                             lhs.data[6] % rhs.data[6],
                             lhs.data[7] % rhs.data[7]);
#endif
  }

  /**
   * Vectorized 4x float.
   */
  template<>
  struct vtype<float, 4> {
    static constexpr unsigned int W = 4;

    typedef float scalar_t;
    typedef vmask<sizeof(float), 4> mask_t;
    typedef StaticSizedArray<float, 4> value_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef __m128 register_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    typedef StaticSizedArray<float, 4> register_t;
#endif

    union {
      value_t data;
      register_t reg;
    };

    vtype() = default ;

    vtype(float e0, float e1, float e2, float e3)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm_set_ps(e3, e2, e1, e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data{e0, e1, e2, e3} { }
#endif

    explicit vtype(float e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm_set1_ps(e)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(e) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    explicit vtype(__m128 reg) : reg(reg) { }
#endif

    explicit vtype(float const * mem)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm_loadu_ps(mem)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(mem) { }
#endif

    vtype(vtype<float, 4> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(rhs.data) { }
#endif

    vtype<float, 4> & operator=(vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
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
    reg = _mm_add_ps(reg, rhs.reg);
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
    reg = _mm_sub_ps(reg, rhs.reg);
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
    reg = _mm_mul_ps(reg, rhs.reg);
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
    reg = _mm_div_ps(reg, rhs.reg);
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
    return vtype<float, 4>(_mm_add_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(lhs.data[0] + rhs.data[0],
                           lhs.data[1] + rhs.data[1],
                           lhs.data[2] + rhs.data[2],
                           lhs.data[3] + rhs.data[3]);
#endif
  }

  inline vtype<float, 4> operator-(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_sub_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(lhs.data[0] - rhs.data[0],
                           lhs.data[1] - rhs.data[1],
                           lhs.data[2] - rhs.data[2],
                           lhs.data[3] - rhs.data[3]);
#endif
  }

  inline vtype<float, 4> operator*(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_mul_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(lhs.data[0] * rhs.data[0],
                           lhs.data[1] * rhs.data[1],
                           lhs.data[2] * rhs.data[2],
                           lhs.data[3] * rhs.data[3]);
#endif
  }

  inline vtype<float, 4> operator/(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_div_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(lhs.data[0] / rhs.data[0],
                           lhs.data[1] / rhs.data[1],
                           lhs.data[2] / rhs.data[2],
                           lhs.data[3] / rhs.data[3]);
#endif
  }

  inline vtype<float, 4> operator-(vtype<float, 4> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m128 sign_mask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    return vtype<float, 4>(_mm_xor_ps(v.reg, sign_mask));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(-v.data[0],
                           -v.data[1],
                           -v.data[2],
                           -v.data[3]);
#endif
  }

  inline vmask<sizeof(float), 4> operator<(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_LT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 4>((__mmask8)_mm512_cmp_ps_mask(
      _mm512_castps128_ps512(lhs.reg), _mm512_castps128_ps512(rhs.reg), _CMP_LT_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps(lhs.reg, rhs.reg, _CMP_LT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 4>(lhs.data[0] < rhs.data[0],
                                   lhs.data[1] < rhs.data[1],
                                   lhs.data[2] < rhs.data[2],
                                   lhs.data[3] < rhs.data[3]);
#endif
  }

  inline vmask<sizeof(float), 4> operator<=(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_LE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 4>((__mmask8)_mm512_cmp_ps_mask(
      _mm512_castps128_ps512(lhs.reg), _mm512_castps128_ps512(rhs.reg), _CMP_LE_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps(lhs.reg, rhs.reg, _CMP_LE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 4>(lhs.data[0] <= rhs.data[0],
                                   lhs.data[1] <= rhs.data[1],
                                   lhs.data[2] <= rhs.data[2],
                                   lhs.data[3] <= rhs.data[3]);
#endif
  }

  inline vmask<sizeof(float), 4> operator>(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_GT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 4>((__mmask8)_mm512_cmp_ps_mask(
      _mm512_castps128_ps512(lhs.reg), _mm512_castps128_ps512(rhs.reg), _CMP_GT_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps(lhs.reg, rhs.reg, _CMP_GT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 4>(lhs.data[0] > rhs.data[0],
                                   lhs.data[1] > rhs.data[1],
                                   lhs.data[2] > rhs.data[2],
                                   lhs.data[3] > rhs.data[3]);
#endif
  }

  inline vmask<sizeof(float), 4> operator>=(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_GE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 4>((__mmask8)_mm512_cmp_ps_mask(
      _mm512_castps128_ps512(lhs.reg), _mm512_castps128_ps512(rhs.reg), _CMP_GE_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps(lhs.reg, rhs.reg, _CMP_GE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 4>(lhs.data[0] >= rhs.data[0],
                                   lhs.data[1] >= rhs.data[1],
                                   lhs.data[2] >= rhs.data[2],
                                   lhs.data[3] >= rhs.data[3]);
#endif
  }

  inline vmask<sizeof(float), 4> operator==(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_EQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 4>((__mmask8)_mm512_cmp_ps_mask(
      _mm512_castps128_ps512(lhs.reg), _mm512_castps128_ps512(rhs.reg), _CMP_EQ_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps(lhs.reg, rhs.reg, _CMP_EQ_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 4>(lhs.data[0] == rhs.data[0],
                                   lhs.data[1] == rhs.data[1],
                                   lhs.data[2] == rhs.data[2],
                                   lhs.data[3] == rhs.data[3]);
#endif
  }

  inline vmask<sizeof(float), 4> operator!=(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_NEQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 4>((__mmask8)_mm512_cmp_ps_mask(
      _mm512_castps128_ps512(lhs.reg), _mm512_castps128_ps512(rhs.reg), _CMP_NEQ_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 4>(_mm_cmp_ps(lhs.reg, rhs.reg, _CMP_NEQ_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 4>(lhs.data[0] != rhs.data[0],
                                   lhs.data[1] != rhs.data[1],
                                   lhs.data[2] != rhs.data[2],
                                   lhs.data[3] != rhs.data[3]);
#endif
  }

  inline void mem_store(float * mem, vtype<float, 4> const & op) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    _mm_storeu_ps(mem, op.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    mem[0] = op.data[0];
    mem[1] = op.data[1];
    mem[2] = op.data[2];
    mem[3] = op.data[3];
#endif
  }

  inline vtype<float, 4> select(vmask<sizeof(float), 4> const & c,
                                vtype<float, 4> const & t,
                                vtype<float, 4> const & f) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512DQ__) && defined(__AVX512VL__)
    return vtype<float, 4>(_mm_mask_blend_ps(c.reg, f.reg, t.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    __m512 ta = _mm512_castps128_ps512(t.reg);
    __m512 fa = _mm512_castps128_ps512(f.reg);
    __m512 res = _mm512_mask_blend_ps((__mmask16)c.reg, fa, ta);
    return vtype<float, 4>(_mm512_castps512_ps128(res));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_blendv_ps(f.reg, t.reg, _mm_castsi128_ps(c.reg)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>((c.data & 1u) ? t.data[0] : f.data[0],
                           (c.data & 2u) ? t.data[1] : f.data[1],
                           (c.data & 4u) ? t.data[2] : f.data[2],
                           (c.data & 8u) ? t.data[3] : f.data[3]);
#endif
  }

  inline vtype<float, 4> fabs(vtype<float, 4> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m128 sign_mask = _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF));
    return vtype<float, 4>(_mm_and_ps(sign_mask, v.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(std::fabs(v.data[0]),
                           std::fabs(v.data[1]),
                           std::fabs(v.data[2]),
                           std::fabs(v.data[3]));
#endif
  }

  inline vtype<float, 4> abs(vtype<float, 4> const & v) {
    return fabs(v);
  }

  inline vtype<float, 4> max(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_max_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(std::max(lhs.data[0], rhs.data[0]),
                           std::max(lhs.data[1], rhs.data[1]),
                           std::max(lhs.data[2], rhs.data[2]),
                           std::max(lhs.data[3], rhs.data[3]));
#endif
  }

  inline vtype<float, 4> min(vtype<float, 4> const & lhs, vtype<float, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_min_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(std::min(lhs.data[0], rhs.data[0]),
                           std::min(lhs.data[1], rhs.data[1]),
                           std::min(lhs.data[2], rhs.data[2]),
                           std::min(lhs.data[3], rhs.data[3]));
#endif
  }

  inline vtype<float, 4> ceil(vtype<float, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_ceil_ps(arg.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(std::ceil(arg.data[0]),
                           std::ceil(arg.data[1]),
                           std::ceil(arg.data[2]),
                           std::ceil(arg.data[3]));
#endif
  }

  inline vtype<float, 4> floor(vtype<float, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_floor_ps(arg.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(std::floor(arg.data[0]),
                           std::floor(arg.data[1]),
                           std::floor(arg.data[2]),
                           std::floor(arg.data[3]));
#endif
  }

  inline vtype<float, 4> round(vtype<float, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_round_ps(arg.reg, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(std::round(arg.data[0]),
                           std::round(arg.data[1]),
                           std::round(arg.data[2]),
                           std::round(arg.data[3]));
#endif
  }

  inline vtype<float, 4> sqrt(vtype<float, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_sqrt_ps(arg.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(std::sqrt(arg.data[0]),
                           std::sqrt(arg.data[1]),
                           std::sqrt(arg.data[2]),
                           std::sqrt(arg.data[3])) ;
#endif
      }

  inline vtype<float, 4> rsqrt(vtype<float, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vtype<float, 4>(_mm_rsqrt14_ps(arg.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 4>(_mm_rsqrt_ps(arg.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 4>(1.0f/std::sqrt(arg.data[0]),
                           1.0f/std::sqrt(arg.data[1]),
                           1.0f/std::sqrt(arg.data[2]),
                           1.0f/std::sqrt(arg.data[3]));
#endif
  }

  /**
   * Vectorized 8x float.
   */
  template<>
  struct vtype<float, 8> {
    static constexpr unsigned int W = 8;

    typedef float scalar_t;
    typedef vmask<sizeof(float), 8> mask_t;
    typedef StaticSizedArray<float, 8> value_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef __m256 register_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    typedef StaticSizedArray<float, 8> register_t;
#endif

    union {
      value_t data;
      register_t reg;
    };

    vtype() = default ;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    explicit vtype(__m256 reg) : reg(reg) { }
#endif

    vtype(float e0, float e1, float e2, float e3, float e4, float e5, float e6, float e7)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm256_set_ps(e7, e6, e5, e4, e3, e2, e1, e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    : data{e0, e1, e2, e3, e4, e5, e6, e7} { }
#endif

    explicit vtype(float e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm256_set1_ps(e)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(e) { }
#endif

    explicit vtype(float const * mem)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm256_loadu_ps(mem)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
      : data(mem) { }
#endif

    vtype(vtype<float, 8> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(rhs.data) { }
#endif

    vtype<float, 8> & operator=(vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
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
    reg = _mm256_add_ps(reg, rhs.reg);
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
    reg = _mm256_sub_ps(reg, rhs.reg);
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
    reg = _mm256_mul_ps(reg, rhs.reg);
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
    reg = _mm256_div_ps(reg, rhs.reg);
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
    return vtype<float, 8>(_mm256_add_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(lhs.data[0] + rhs.data[0],
                           lhs.data[1] + rhs.data[1],
                           lhs.data[2] + rhs.data[2],
                           lhs.data[3] + rhs.data[3],
                           lhs.data[4] + rhs.data[4],
                           lhs.data[5] + rhs.data[5],
                           lhs.data[6] + rhs.data[6],
                           lhs.data[7] + rhs.data[7]);
#endif
  }

  inline vtype<float, 8> operator-(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_sub_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(lhs.data[0] - rhs.data[0],
                           lhs.data[1] - rhs.data[1],
                           lhs.data[2] - rhs.data[2],
                           lhs.data[3] - rhs.data[3],
                           lhs.data[4] - rhs.data[4],
                           lhs.data[5] - rhs.data[5],
                           lhs.data[6] - rhs.data[6],
                           lhs.data[7] - rhs.data[7]);
#endif
  }

  inline vtype<float, 8> operator*(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_mul_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(lhs.data[0] * rhs.data[0],
                           lhs.data[1] * rhs.data[1],
                           lhs.data[2] * rhs.data[2],
                           lhs.data[3] * rhs.data[3],
                           lhs.data[4] * rhs.data[4],
                           lhs.data[5] * rhs.data[5],
                           lhs.data[6] * rhs.data[6],
                           lhs.data[7] * rhs.data[7]);
#endif
  }

  inline vtype<float, 8> operator*(vtype<float, 8> const & lhs, float const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_mul_ps(lhs.reg, _mm256_set1_ps(rhs)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(lhs.data[0] * rhs,
                           lhs.data[1] * rhs,
                           lhs.data[2] * rhs,
                           lhs.data[3] * rhs,
                           lhs.data[4] * rhs,
                           lhs.data[5] * rhs,
                           lhs.data[6] * rhs,
                           lhs.data[7] * rhs);
#endif
  }

  inline vtype<float, 8> operator*(float const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_mul_ps(_mm256_set1_ps(lhs), rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(lhs * rhs.data[0],
                           lhs * rhs.data[1],
                           lhs * rhs.data[2],
                           lhs * rhs.data[3],
                           lhs * rhs.data[4],
                           lhs * rhs.data[5],
                           lhs * rhs.data[6],
                           lhs * rhs.data[7]);
#endif
  }

  inline vtype<float, 8> operator*(vtype<float, 8> const & lhs, double const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_mul_ps(lhs.reg, _mm256_set1_ps(float(rhs))));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(lhs.data[0] * rhs,
                           lhs.data[1] * rhs,
                           lhs.data[2] * rhs,
                           lhs.data[3] * rhs,
                           lhs.data[4] * rhs,
                           lhs.data[5] * rhs,
                           lhs.data[6] * rhs,
                           lhs.data[7] * rhs);
#endif
  }

  inline vtype<float, 8> operator*(double const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_mul_ps(_mm256_set1_ps(float(lhs)), rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(lhs * rhs.data[0],
                           lhs * rhs.data[1],
                           lhs * rhs.data[2],
                           lhs * rhs.data[3],
                           lhs * rhs.data[4],
                           lhs * rhs.data[5],
                           lhs * rhs.data[6],
                           lhs * rhs.data[7]);
#endif
  }
  
  inline vtype<float, 8> operator/(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_div_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(lhs.data[0] / rhs.data[0],
                           lhs.data[1] / rhs.data[1],
                           lhs.data[2] / rhs.data[2],
                           lhs.data[3] / rhs.data[3],
                           lhs.data[4] / rhs.data[4],
                           lhs.data[5] / rhs.data[5],
                           lhs.data[6] / rhs.data[6],
                           lhs.data[7] / rhs.data[7]);
#endif
  }

  inline vtype<float, 8> operator/(vtype<float, 8> const & lhs, float const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_div_ps(lhs.reg, _mm256_set1_ps(rhs)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(lhs.data[0] / rhs,
                           lhs.data[1] / rhs,
                           lhs.data[2] / rhs,
                           lhs.data[3] / rhs,
                           lhs.data[4] / rhs,
                           lhs.data[5] / rhs,
                           lhs.data[6] / rhs,
                           lhs.data[7] / rhs);
#endif
  }

  inline vtype<float, 8> operator-(vtype<float, 8> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256 sign_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
    return vtype<float, 8>(_mm256_xor_ps(v.reg, sign_mask));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(-v.data[0],
                           -v.data[1],
                           -v.data[2],
                           -v.data[3],
                           -v.data[4],
                           -v.data[5],
                           -v.data[6],
                           -v.data[7]);
#endif
  }

  inline vmask<sizeof(float), 8> operator<(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_LT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 8>(
      (__mmask8)_mm512_cmp_ps_mask(
        _mm512_castps256_ps512(lhs.reg),
        _mm512_castps256_ps512(rhs.reg),
        _CMP_LT_OQ
      )
    );
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps(lhs.reg, rhs.reg, _CMP_LT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 8>(lhs.data[0] < rhs.data[0],
                                   lhs.data[1] < rhs.data[1],
                                   lhs.data[2] < rhs.data[2],
                                   lhs.data[3] < rhs.data[3],
                                   lhs.data[4] < rhs.data[4],
                                   lhs.data[5] < rhs.data[5],
                                   lhs.data[6] < rhs.data[6],
                                   lhs.data[7] < rhs.data[7]);
#endif
  }

  inline vmask<sizeof(float), 8> operator<=(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_LE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 8>(
      (__mmask8)_mm512_cmp_ps_mask(
        _mm512_castps256_ps512(lhs.reg),
        _mm512_castps256_ps512(rhs.reg),
        _CMP_LE_OQ
      )
    );
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps(lhs.reg, rhs.reg, _CMP_LE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 8>(lhs.data[0] <= rhs.data[0],
                                   lhs.data[1] <= rhs.data[1],
                                   lhs.data[2] <= rhs.data[2],
                                   lhs.data[3] <= rhs.data[3],
                                   lhs.data[4] <= rhs.data[4],
                                   lhs.data[5] <= rhs.data[5],
                                   lhs.data[6] <= rhs.data[6],
                                   lhs.data[7] <= rhs.data[7]);
#endif
  }

  inline vmask<sizeof(float), 8> operator>(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_GT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 8>(
      (__mmask8)_mm512_cmp_ps_mask(
        _mm512_castps256_ps512(lhs.reg),
        _mm512_castps256_ps512(rhs.reg),
        _CMP_GT_OQ
      )
    );
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps(lhs.reg, rhs.reg, _CMP_GT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 8>(lhs.data[0] > rhs.data[0],
                                   lhs.data[1] > rhs.data[1],
                                   lhs.data[2] > rhs.data[2],
                                   lhs.data[3] > rhs.data[3],
                                   lhs.data[4] > rhs.data[4],
                                   lhs.data[5] > rhs.data[5],
                                   lhs.data[6] > rhs.data[6],
                                   lhs.data[7] > rhs.data[7]);
#endif
  }

  inline vmask<sizeof(float), 8> operator>=(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_GE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 8>(
      (__mmask8)_mm512_cmp_ps_mask(
        _mm512_castps256_ps512(lhs.reg),
        _mm512_castps256_ps512(rhs.reg),
        _CMP_GE_OQ
      )
    );
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps(lhs.reg, rhs.reg, _CMP_GE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 8>(lhs.data[0] >= rhs.data[0],
                                   lhs.data[1] >= rhs.data[1],
                                   lhs.data[2] >= rhs.data[2],
                                   lhs.data[3] >= rhs.data[3],
                                   lhs.data[4] >= rhs.data[4],
                                   lhs.data[5] >= rhs.data[5],
                                   lhs.data[6] >= rhs.data[6],
                                   lhs.data[7] >= rhs.data[7]);
#endif
  }

  inline vmask<sizeof(float), 8> operator==(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_EQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 8>(
      (__mmask8)_mm512_cmp_ps_mask(
        _mm512_castps256_ps512(lhs.reg),
        _mm512_castps256_ps512(rhs.reg),
        _CMP_EQ_OQ
      )
    );
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps(lhs.reg, rhs.reg, _CMP_EQ_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 8>(lhs.data[0] == rhs.data[0],
                                   lhs.data[1] == rhs.data[1],
                                   lhs.data[2] == rhs.data[2],
                                   lhs.data[3] == rhs.data[3],
                                   lhs.data[4] == rhs.data[4],
                                   lhs.data[5] == rhs.data[5],
                                   lhs.data[6] == rhs.data[6],
                                   lhs.data[7] == rhs.data[7]);
#endif
  }

  inline vmask<sizeof(float), 8> operator!=(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_NEQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 8>(
      (__mmask8)_mm512_cmp_ps_mask(
        _mm512_castps256_ps512(lhs.reg),
        _mm512_castps256_ps512(rhs.reg),
        _CMP_NEQ_OQ
      )
    );
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 8>(_mm256_cmp_ps(lhs.reg, rhs.reg, _CMP_NEQ_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 8>(lhs.data[0] != rhs.data[0],
                                   lhs.data[1] != rhs.data[1],
                                   lhs.data[2] != rhs.data[2],
                                   lhs.data[3] != rhs.data[3],
                                   lhs.data[4] != rhs.data[4],
                                   lhs.data[5] != rhs.data[5],
                                   lhs.data[6] != rhs.data[6],
                                   lhs.data[7] != rhs.data[7]);
#endif
  }

  inline void mem_store(float * mem, vtype<float, 8> const & op) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    _mm256_storeu_ps(mem, op.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    mem[0] = op.data[0];
    mem[1] = op.data[1];
    mem[2] = op.data[2];
    mem[3] = op.data[3];
    mem[4] = op.data[4];
    mem[5] = op.data[5];
    mem[6] = op.data[6];
    mem[7] = op.data[7];
#endif
  }

  inline vtype<float, 8> select(vmask<sizeof(float), 8> const & c,
                                vtype<float, 8> const & t,
                                vtype<float, 8> const & f) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    __m512 ta = _mm512_castps256_ps512(t.reg);
    __m512 fa = _mm512_castps256_ps512(f.reg);
    __m512 res = _mm512_mask_blend_ps((__mmask16)c.reg, fa, ta);
    return vtype<float, 8>(_mm512_castps512_ps256(res));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_blendv_ps(f.reg, t.reg, _mm256_castsi256_ps(c.reg)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>((c.data & 1u) ? t.data[0] : f.data[0],
                           (c.data & 2u) ? t.data[1] : f.data[1],
                           (c.data & 4u) ? t.data[2] : f.data[2],
                           (c.data & 8u) ? t.data[3] : f.data[3],
                           (c.data & 16u) ? t.data[4] : f.data[4],
                           (c.data & 32u) ? t.data[5] : f.data[5],
                           (c.data & 64u) ? t.data[6] : f.data[6],
                           (c.data & 128u) ? t.data[7] : f.data[7]);
#endif
  }

  inline vtype<float, 8> fabs(vtype<float, 8> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256 sign_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
    return vtype<float, 8>(_mm256_and_ps(sign_mask, v.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(std::fabs(v.data[0]),
                           std::fabs(v.data[1]),
                           std::fabs(v.data[2]),
                           std::fabs(v.data[3]),
                           std::fabs(v.data[4]),
                           std::fabs(v.data[5]),
                           std::fabs(v.data[6]),
                           std::fabs(v.data[7]));
#endif
  }

  inline vtype<float, 8> abs(vtype<float, 8> const & v) {
    return fabs(v);
  }

  inline vtype<float, 8> max(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_max_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(std::max(lhs.data[0], rhs.data[0]),
                           std::max(lhs.data[1], rhs.data[1]),
                           std::max(lhs.data[2], rhs.data[2]),
                           std::max(lhs.data[3], rhs.data[3]),
                           std::max(lhs.data[4], rhs.data[4]),
                           std::max(lhs.data[5], rhs.data[5]),
                           std::max(lhs.data[6], rhs.data[6]),
                           std::max(lhs.data[7], rhs.data[7]));
#endif
  }

  inline vtype<float, 8> min(vtype<float, 8> const & lhs, vtype<float, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_min_ps(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(std::min(lhs.data[0], rhs.data[0]),
                           std::min(lhs.data[1], rhs.data[1]),
                           std::min(lhs.data[2], rhs.data[2]),
                           std::min(lhs.data[3], rhs.data[3]),
                           std::min(lhs.data[4], rhs.data[4]),
                           std::min(lhs.data[5], rhs.data[5]),
                           std::min(lhs.data[6], rhs.data[6]),
                           std::min(lhs.data[7], rhs.data[7]));
#endif
  }

  inline vtype<float, 8> ceil(vtype<float, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_ceil_ps(arg.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(std::ceil(arg.data[0]),
                           std::ceil(arg.data[1]),
                           std::ceil(arg.data[2]),
                           std::ceil(arg.data[3]),
                           std::ceil(arg.data[4]),
                           std::ceil(arg.data[5]),
                           std::ceil(arg.data[6]),
                           std::ceil(arg.data[7]));
#endif
  }

  inline vtype<float, 8> floor(vtype<float, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_floor_ps(arg.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(std::floor(arg.data[0]),
                           std::floor(arg.data[1]),
                           std::floor(arg.data[2]),
                           std::floor(arg.data[3]),
                           std::floor(arg.data[4]),
                           std::floor(arg.data[5]),
                           std::floor(arg.data[6]),
                           std::floor(arg.data[7]));
#endif
  }

  inline vtype<float, 8> round(vtype<float, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_round_ps(arg.reg, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(std::round(arg.data[0]),
                           std::round(arg.data[1]),
                           std::round(arg.data[2]),
                           std::round(arg.data[3]),
                           std::round(arg.data[4]),
                           std::round(arg.data[5]),
                           std::round(arg.data[6]),
                           std::round(arg.data[7]));
#endif
  }

  inline vtype<float, 8> sqrt(vtype<float, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_sqrt_ps(arg.reg));
#else
#if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(std::sqrt(arg.data[0]),
                           std::sqrt(arg.data[1]),
                           std::sqrt(arg.data[2]),
                           std::sqrt(arg.data[3]),
                           std::sqrt(arg.data[4]),
                           std::sqrt(arg.data[5]),
                           std::sqrt(arg.data[6]),
                           std::sqrt(arg.data[7]));
#endif
  }

  inline vtype<float, 8> rsqrt(vtype<float, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vtype<float, 8>(_mm256_rsqrt14_ps(arg.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    __m512 zero = _mm512_setzero_ps();
    __m512 v = _mm512_mask_blend_ps(0x00FF, zero, _mm512_castps256_ps512(arg.reg));
    return vtype<float, 8>(_mm512_castps512_ps256(_mm512_rsqrt14_ps(v)));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 8>(_mm256_rsqrt_ps(arg.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 8>(1.0f/std::sqrt(arg.data[0]),
                           1.0f/std::sqrt(arg.data[1]),
                           1.0f/std::sqrt(arg.data[2]),
                           1.0f/std::sqrt(arg.data[3]),
                           1.0f/std::sqrt(arg.data[4]),
                           1.0f/std::sqrt(arg.data[5]),
                           1.0f/std::sqrt(arg.data[6]),
                           1.0f/std::sqrt(arg.data[7]));
#endif
  }

  /**
   * Vectorized 16x float
   */
  template<>
  struct vtype<float, 16> {
    static constexpr unsigned int W = 16;

    typedef float scalar_t;
    typedef vmask<sizeof(float), 16> mask_t;
    typedef StaticSizedArray<float, 16> value_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    typedef __m512 register_t;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef struct {
      __m256 x, y;
    } register_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    typedef StaticSizedArray<float, 16> register_t;
#endif

    union {
      value_t data;
      register_t reg;
    };

    vtype() = default;

    vtype(float e0, float e1, float e2, float e3,
          float e4, float e5, float e6, float e7,
          float e8, float e9, float e10, float e11,
          float e12, float e13, float e14, float e15)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(_mm512_set_ps(e15, e14, e13, e12,
                          e11, e10, e9, e8,
                          e7, e6, e5, e4,
                          e3, e2, e1, e0)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_set_ps(e7, e6, e5, e4, e3, e2, e1, e0),
            _mm256_set_ps(e15, e14, e13, e12, e11, e10, e9, e8)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(e0, e1, e2, e3, e4, e5, e6, e7,
             e8, e9, e10, e11, e12, e13, e14, e15) { }
#endif

    explicit vtype(float e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(_mm512_set1_ps(e)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_set1_ps(e), _mm256_set1_ps(e)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(e) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__)
    explicit vtype(__m256 reg0, __m256 reg1) : reg{reg0, reg1} { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    explicit vtype(__m512 reg) : reg(reg) { }
#endif

    explicit vtype(float const * mem)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(_mm512_loadu_ps(mem)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_loadu_ps(mem), _mm256_loadu_ps(mem+8)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(mem) { }
#endif

    vtype(vtype<float, 16> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(rhs.reg) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(rhs.data) { }
#endif

    vtype<float, 16> & operator=(vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg = rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AXV2__)
      reg = rhs.reg
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

    vtype<float, 16> & operator+=(vtype<float, 16> const & rhs);
    vtype<float, 16> & operator-=(vtype<float, 16> const & rhs);
    vtype<float, 16> & operator*=(vtype<float, 16> const & rhs);
    vtype<float, 16> & operator/=(vtype<float, 16> const & rhs);
  };

  inline vtype<float, 16> & vtype<float, 16>::operator+=(vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg = _mm512_add_ps(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_add_ps(reg.x, rhs.reg.x);
    reg.y = _mm256_add_ps(reg.y, rhs.reg.y);
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
    data[8] += rhs.data[8];
    data[9] += rhs.data[9];
    data[10] += rhs.data[10];
    data[11] += rhs.data[11];
    data[12] += rhs.data[12];
    data[13] += rhs.data[13];
    data[14] += rhs.data[14];
    data[15] += rhs.data[15];
#endif
    return *this;
  }

  inline vtype<float, 16> & vtype<float, 16>::operator-=(vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg = _mm512_sub_ps(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_sub_ps(reg.x, rhs.reg.x);
    reg.y = _mm256_sub_ps(reg.y, rhs.reg.y);
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
    data[8] -= rhs.data[8];
    data[9] -= rhs.data[9];
    data[10] -= rhs.data[10];
    data[11] -= rhs.data[11];
    data[12] -= rhs.data[12];
    data[13] -= rhs.data[13];
    data[14] -= rhs.data[14];
    data[15] -= rhs.data[15];
#endif
    return *this;
  }

  inline vtype<float, 16> & vtype<float, 16>::operator*=(vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg = _mm512_mul_ps(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_mul_ps(reg.x, rhs.reg.x);
    reg.y = _mm256_mul_ps(reg.y, rhs.reg.y);
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
    data[8] *= rhs.data[8];
    data[9] *= rhs.data[9];
    data[10] *= rhs.data[10];
    data[11] *= rhs.data[11];
    data[12] *= rhs.data[12];
    data[13] *= rhs.data[13];
    data[14] *= rhs.data[14];
    data[15] *= rhs.data[15];
#endif
    return *this;
  }

  inline vtype<float, 16> & vtype<float, 16>::operator/=(vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg = _mm512_div_ps(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_div_ps(reg.x, rhs.reg.x);
    reg.y = _mm256_div_ps(reg.y, rhs.reg.y);
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
    data[8] /= rhs.data[8];
    data[9] /= rhs.data[9];
    data[10] /= rhs.data[10];
    data[11] /= rhs.data[11];
    data[12] /= rhs.data[12];
    data[13] /= rhs.data[13];
    data[14] /= rhs.data[14];
    data[15] /= rhs.data[15];
#endif
    return *this;
  }

  inline vtype<float, 16> operator+(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_add_ps(lhs.reg, rhs.reg));
#elif defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_add_ps(lhs.reg.x, rhs.reg.x),
                            _mm256_add_ps(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(lhs.data[0] + rhs.data[0],
                            lhs.data[1] + rhs.data[1],
                            lhs.data[2] + rhs.data[2],
                            lhs.data[3] + rhs.data[3],
                            lhs.data[4] + rhs.data[4],
                            lhs.data[5] + rhs.data[5],
                            lhs.data[6] + rhs.data[6],
                            lhs.data[7] + rhs.data[7],
                            lhs.data[8] + rhs.data[8],
                            lhs.data[9] + rhs.data[9],
                            lhs.data[10] + rhs.data[10],
                            lhs.data[11] + rhs.data[11],
                            lhs.data[12] + rhs.data[12],
                            lhs.data[13] + rhs.data[13],
                            lhs.data[14] + rhs.data[14],
                            lhs.data[15] + rhs.data[15]);
#endif
  }

  inline vtype<float, 16> operator-(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_sub_ps(lhs.reg, rhs.reg));
#elif defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_sub_ps(lhs.reg.x, rhs.reg.x),
                            _mm256_sub_ps(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(lhs.data[0] - rhs.data[0],
                            lhs.data[1] - rhs.data[1],
                            lhs.data[2] - rhs.data[2],
                            lhs.data[3] - rhs.data[3],
                            lhs.data[4] - rhs.data[4],
                            lhs.data[5] - rhs.data[5],
                            lhs.data[6] - rhs.data[6],
                            lhs.data[7] - rhs.data[7],
                            lhs.data[8] - rhs.data[8],
                            lhs.data[9] - rhs.data[9],
                            lhs.data[10] - rhs.data[10],
                            lhs.data[11] - rhs.data[11],
                            lhs.data[12] - rhs.data[12],
                            lhs.data[13] - rhs.data[13],
                            lhs.data[14] - rhs.data[14],
                            lhs.data[15] - rhs.data[15]);
#endif
  }

  inline vtype<float, 16> operator*(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_mul_ps(lhs.reg, rhs.reg));
#elif defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_mul_ps(lhs.reg.x, rhs.reg.x),
                            _mm256_mul_ps(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(lhs.data[0] * rhs.data[0],
                            lhs.data[1] * rhs.data[1],
                            lhs.data[2] * rhs.data[2],
                            lhs.data[3] * rhs.data[3],
                            lhs.data[4] * rhs.data[4],
                            lhs.data[5] * rhs.data[5],
                            lhs.data[6] * rhs.data[6],
                            lhs.data[7] * rhs.data[7],
                            lhs.data[8] * rhs.data[8],
                            lhs.data[9] * rhs.data[9],
                            lhs.data[10] * rhs.data[10],
                            lhs.data[11] * rhs.data[11],
                            lhs.data[12] * rhs.data[12],
                            lhs.data[13] * rhs.data[13],
                            lhs.data[14] * rhs.data[14],
                            lhs.data[15] * rhs.data[15]);
#endif
  }

  inline vtype<float, 16> operator/(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_div_ps(lhs.reg, rhs.reg));
#elif defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_div_ps(lhs.reg.x, rhs.reg.x),
                            _mm256_div_ps(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(lhs.data[0] / rhs.data[0],
                            lhs.data[1] / rhs.data[1],
                            lhs.data[2] / rhs.data[2],
                            lhs.data[3] / rhs.data[3],
                            lhs.data[4] / rhs.data[4],
                            lhs.data[5] / rhs.data[5],
                            lhs.data[6] / rhs.data[6],
                            lhs.data[7] / rhs.data[7],
                            lhs.data[8] / rhs.data[8],
                            lhs.data[9] / rhs.data[9],
                            lhs.data[10] / rhs.data[10],
                            lhs.data[11] / rhs.data[11],
                            lhs.data[12] / rhs.data[12],
                            lhs.data[13] / rhs.data[13],
                            lhs.data[14] / rhs.data[14],
                            lhs.data[15] / rhs.data[15]);
#endif
  }

  inline vtype<float, 16> operator-(vtype<float, 16> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    __m512 sign_mask = _mm512_castsi512_ps(_mm512_set1_epi32(0x80000000));
    return vtype<float, 16>(_mm512_xor_ps(sign_mask, v.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    __m512i sign_mask = _mm512_set1_epi32(0x80000000);
    return vtype<float, 16>(_mm512_castsi512_ps(
      _mm512_xor_si512(sign_mask, _mm512_castps_si512(v.reg))
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256 sign_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
    return vtype<float, 16>(_mm256_xor_ps(sign_mask, v.reg.x),
                            _mm256_xor_ps(sign_mask, v.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(-v.data[0],
                            -v.data[1],
                            -v.data[2],
                            -v.data[3],
                            -v.data[4],
                            -v.data[5],
                            -v.data[6],
                            -v.data[7],
                            -v.data[8],
                            -v.data[9],
                            -v.data[10],
                            -v.data[11],
                            -v.data[12],
                            -v.data[13],
                            -v.data[14],
                            -v.data[15]);
#endif
  }

  inline vmask<sizeof(float), 16> operator<(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 16>(_mm512_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_LT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 16>(_mm256_cmp_ps(lhs.reg.x, rhs.reg.x, _CMP_LT_OQ),
                                    _mm256_cmp_ps(lhs.reg.y, rhs.reg.y, _CMP_LT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 16>(lhs.data[0] < rhs.data[0],
                                    lhs.data[1] < rhs.data[1],
                                    lhs.data[2] < rhs.data[2],
                                    lhs.data[3] < rhs.data[3],
                                    lhs.data[4] < rhs.data[4],
                                    lhs.data[5] < rhs.data[5],
                                    lhs.data[6] < rhs.data[6],
                                    lhs.data[7] < rhs.data[7],
                                    lhs.data[8] < rhs.data[8],
                                    lhs.data[9] < rhs.data[9],
                                    lhs.data[10] < rhs.data[10],
                                    lhs.data[11] < rhs.data[11],
                                    lhs.data[12] < rhs.data[12],
                                    lhs.data[13] < rhs.data[13],
                                    lhs.data[14] < rhs.data[14],
                                    lhs.data[15] < rhs.data[15]);
#endif
  }

  inline vmask<sizeof(float), 16> operator<=(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 16>(_mm512_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_LE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 16>(_mm256_cmp_ps(lhs.reg.x, rhs.reg.x, _CMP_LE_OQ),
                                    _mm256_cmp_ps(lhs.reg.y, rhs.reg.y, _CMP_LE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 16>(lhs.data[0] <= rhs.data[0],
                                    lhs.data[1] <= rhs.data[1],
                                    lhs.data[2] <= rhs.data[2],
                                    lhs.data[3] <= rhs.data[3],
                                    lhs.data[4] <= rhs.data[4],
                                    lhs.data[5] <= rhs.data[5],
                                    lhs.data[6] <= rhs.data[6],
                                    lhs.data[7] <= rhs.data[7],
                                    lhs.data[8] <= rhs.data[8],
                                    lhs.data[9] <= rhs.data[9],
                                    lhs.data[10] <= rhs.data[10],
                                    lhs.data[11] <= rhs.data[11],
                                    lhs.data[12] <= rhs.data[12],
                                    lhs.data[13] <= rhs.data[13],
                                    lhs.data[14] <= rhs.data[14],
                                    lhs.data[15] <= rhs.data[15]);
#endif
  }

  inline vmask<sizeof(float), 16> operator>(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 16>(_mm512_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_GT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 16>(_mm256_cmp_ps(lhs.reg.x, rhs.reg.x, _CMP_GT_OQ),
                                    _mm256_cmp_ps(lhs.reg.y, rhs.reg.y, _CMP_GT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 16>(lhs.data[0] > rhs.data[0],
                                    lhs.data[1] > rhs.data[1],
                                    lhs.data[2] > rhs.data[2],
                                    lhs.data[3] > rhs.data[3],
                                    lhs.data[4] > rhs.data[4],
                                    lhs.data[5] > rhs.data[5],
                                    lhs.data[6] > rhs.data[6],
                                    lhs.data[7] > rhs.data[7],
                                    lhs.data[8] > rhs.data[8],
                                    lhs.data[9] > rhs.data[9],
                                    lhs.data[10] > rhs.data[10],
                                    lhs.data[11] > rhs.data[11],
                                    lhs.data[12] > rhs.data[12],
                                    lhs.data[13] > rhs.data[13],
                                    lhs.data[14] > rhs.data[14],
                                    lhs.data[15] > rhs.data[15]);
#endif
  }

  inline vmask<sizeof(float), 16> operator>=(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 16>(_mm512_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_GE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 16>(_mm256_cmp_ps(lhs.reg.x, rhs.reg.x, _CMP_GE_OQ),
                                    _mm256_cmp_ps(lhs.reg.y, rhs.reg.y, _CMP_GE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 16>(lhs.data[0] >= rhs.data[0],
                                    lhs.data[1] >= rhs.data[1],
                                    lhs.data[2] >= rhs.data[2],
                                    lhs.data[3] >= rhs.data[3],
                                    lhs.data[4] >= rhs.data[4],
                                    lhs.data[5] >= rhs.data[5],
                                    lhs.data[6] >= rhs.data[6],
                                    lhs.data[7] >= rhs.data[7],
                                    lhs.data[8] >= rhs.data[8],
                                    lhs.data[9] >= rhs.data[9],
                                    lhs.data[10] >= rhs.data[10],
                                    lhs.data[11] >= rhs.data[11],
                                    lhs.data[12] >= rhs.data[12],
                                    lhs.data[13] >= rhs.data[13],
                                    lhs.data[14] >= rhs.data[14],
                                    lhs.data[15] >= rhs.data[15]);
#endif
  }

  inline vmask<sizeof(float), 16> operator==(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 16>(_mm512_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_EQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 16>(_mm256_cmp_ps(lhs.reg.x, rhs.reg.x, _CMP_EQ_OQ),
                                    _mm256_cmp_ps(lhs.reg.y, rhs.reg.y, _CMP_EQ_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 16>(lhs.data[0] == rhs.data[0],
                                    lhs.data[1] == rhs.data[1],
                                    lhs.data[2] == rhs.data[2],
                                    lhs.data[3] == rhs.data[3],
                                    lhs.data[4] == rhs.data[4],
                                    lhs.data[5] == rhs.data[5],
                                    lhs.data[6] == rhs.data[6],
                                    lhs.data[7] == rhs.data[7],
                                    lhs.data[8] == rhs.data[8],
                                    lhs.data[9] == rhs.data[9],
                                    lhs.data[10] == rhs.data[10],
                                    lhs.data[11] == rhs.data[11],
                                    lhs.data[12] == rhs.data[12],
                                    lhs.data[13] == rhs.data[13],
                                    lhs.data[14] == rhs.data[14],
                                    lhs.data[15] == rhs.data[15]);
#endif
  }

  inline vmask<sizeof(float), 16> operator!=(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(float), 16>(_mm512_cmp_ps_mask(lhs.reg, rhs.reg, _CMP_NEQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(float), 16>(_mm256_cmp_ps(lhs.reg.x, rhs.reg.x, _CMP_NEQ_OQ),
                                    _mm256_cmp_ps(lhs.reg.y, rhs.reg.y, _CMP_NEQ_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(float), 16>(lhs.data[0] != rhs.data[0],
                                    lhs.data[1] != rhs.data[1],
                                    lhs.data[2] != rhs.data[2],
                                    lhs.data[3] != rhs.data[3],
                                    lhs.data[4] != rhs.data[4],
                                    lhs.data[5] != rhs.data[5],
                                    lhs.data[6] != rhs.data[6],
                                    lhs.data[7] != rhs.data[7],
                                    lhs.data[8] != rhs.data[8],
                                    lhs.data[9] != rhs.data[9],
                                    lhs.data[10] != rhs.data[10],
                                    lhs.data[11] != rhs.data[11],
                                    lhs.data[12] != rhs.data[12],
                                    lhs.data[13] != rhs.data[13],
                                    lhs.data[14] != rhs.data[14],
                                    lhs.data[15] != rhs.data[15]);
#endif
  }

  inline void mem_store(float * mem, vtype<float, 16> const & op) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    _mm512_storeu_ps(mem, op.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    _mm256_storeu_ps(mem, op.reg.x);
    _mm256_storeu_ps(mem+8, op.reg.y);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    mem[0] = op.data[0];
    mem[1] = op.data[1];
    mem[2] = op.data[2];
    mem[3] = op.data[3];
    mem[4] = op.data[4];
    mem[5] = op.data[5];
    mem[6] = op.data[6];
    mem[7] = op.data[7];
    mem[8] = op.data[8];
    mem[9] = op.data[9];
    mem[10] = op.data[10];
    mem[11] = op.data[11];
    mem[12] = op.data[12];
    mem[13] = op.data[13];
    mem[14] = op.data[14];
    mem[15] = op.data[15];
#endif
  }

  inline vtype<float, 16> select(vmask<sizeof(float), 16> const & c,
                                 vtype<float, 16> const & t,
                                 vtype<float, 16> const & f) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_mask_blend_ps(c.reg, f.reg, t.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_blendv_ps(f.reg.x, t.reg.x, _mm256_castsi256_ps(c.reg.x)),
                            _mm256_blendv_ps(f.reg.y, t.reg.y, _mm256_castsi256_ps(c.reg.y)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>((c.data & 1u) ? t.data[0] : f.data[0],
                            (c.data & 2u) ? t.data[1] : f.data[1],
                            (c.data & 4u) ? t.data[2] : f.data[2],
                            (c.data & 8u) ? t.data[3] : f.data[3],
                            (c.data & 16u) ? t.data[4] : f.data[4],
                            (c.data & 32u) ? t.data[5] : f.data[5],
                            (c.data & 64u) ? t.data[6] : f.data[6],
                            (c.data & 128u) ? t.data[7] : f.data[7],
                            (c.data & 256u) ? t.data[8] : f.data[8],
                            (c.data & 512u) ? t.data[9] : f.data[9],
                            (c.data & 1024u) ? t.data[10] : f.data[10],
                            (c.data & 2048u) ? t.data[11] : f.data[11],
                            (c.data & 4096u) ? t.data[12] : f.data[12],
                            (c.data & 8192u) ? t.data[13] : f.data[13],
                            (c.data & 16384u) ? t.data[14] : f.data[14],
                            (c.data & 32768u) ? t.data[15] : f.data[15]);
#endif
  }

  inline vtype<float, 16> fabs(vtype<float, 16> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    __m512 sign_mask = _mm512_castsi512_ps(_mm512_set1_epi32(0x7FFFFFFF));
    return vtype<float, 16>(_mm512_and_ps(sign_mask, v.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    __m512i sign_mask = _mm512_set1_epi32(0x7FFFFFFF);
    return vtype<float, 16>(_mm512_castsi512_ps(
      _mm512_and_si512(sign_mask, _mm512_castps_si512(v.reg))
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256 sign_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
    return vtype<float, 16>(_mm256_and_ps(sign_mask, v.reg.x),
                            _mm256_and_ps(sign_mask, v.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(std::fabs(v.data[0]),
                            std::fabs(v.data[1]),
                            std::fabs(v.data[2]),
                            std::fabs(v.data[3]),
                            std::fabs(v.data[4]),
                            std::fabs(v.data[5]),
                            std::fabs(v.data[6]),
                            std::fabs(v.data[7]),
                            std::fabs(v.data[8]),
                            std::fabs(v.data[9]),
                            std::fabs(v.data[10]),
                            std::fabs(v.data[11]),
                            std::fabs(v.data[12]),
                            std::fabs(v.data[13]),
                            std::fabs(v.data[14]),
                            std::fabs(v.data[15]));
#endif
  }

  inline vtype<float, 16> abs(vtype<float, 16> const & v) {
    return fabs(v);
  }

  inline vtype<float, 16> max(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_max_ps(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_max_ps(lhs.reg.x, rhs.reg.x),
                            _mm256_max_ps(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(std::max(lhs.data[0], rhs.data[0]),
                            std::max(lhs.data[1], rhs.data[1]),
                            std::max(lhs.data[2], rhs.data[2]),
                            std::max(lhs.data[3], rhs.data[3]),
                            std::max(lhs.data[4], rhs.data[4]),
                            std::max(lhs.data[5], rhs.data[5]),
                            std::max(lhs.data[6], rhs.data[6]),
                            std::max(lhs.data[7], rhs.data[7]),
                            std::max(lhs.data[8], rhs.data[8]),
                            std::max(lhs.data[9], rhs.data[9]),
                            std::max(lhs.data[10], rhs.data[10]),
                            std::max(lhs.data[11], rhs.data[11]),
                            std::max(lhs.data[12], rhs.data[12]),
                            std::max(lhs.data[13], rhs.data[13]),
                            std::max(lhs.data[14], rhs.data[14]),
                            std::max(lhs.data[15], rhs.data[15]));
#endif
  }

  inline vtype<float, 16> min(vtype<float, 16> const & lhs, vtype<float, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_min_ps(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_min_ps(lhs.reg.x, rhs.reg.x),
                            _mm256_min_ps(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(std::min(lhs.data[0], rhs.data[0]),
                            std::min(lhs.data[1], rhs.data[1]),
                            std::min(lhs.data[2], rhs.data[2]),
                            std::min(lhs.data[3], rhs.data[3]),
                            std::min(lhs.data[4], rhs.data[4]),
                            std::min(lhs.data[5], rhs.data[5]),
                            std::min(lhs.data[6], rhs.data[6]),
                            std::min(lhs.data[7], rhs.data[7]),
                            std::min(lhs.data[8], rhs.data[8]),
                            std::min(lhs.data[9], rhs.data[9]),
                            std::min(lhs.data[10], rhs.data[10]),
                            std::min(lhs.data[11], rhs.data[11]),
                            std::min(lhs.data[12], rhs.data[12]),
                            std::min(lhs.data[13], rhs.data[13]),
                            std::min(lhs.data[14], rhs.data[14]),
                            std::min(lhs.data[15], rhs.data[15]));
#endif
  }

  inline vtype<float, 16> ceil(vtype<float, 16> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_roundscale_ps(arg.reg, _MM_FROUND_TO_POS_INF|_MM_FROUND_NO_EXC));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_ceil_ps(arg.reg.x),
                            _mm256_ceil_ps(arg.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(std::ceil(arg.data[0]),
                            std::ceil(arg.data[1]),
                            std::ceil(arg.data[2]),
                            std::ceil(arg.data[3]),
                            std::ceil(arg.data[4]),
                            std::ceil(arg.data[5]),
                            std::ceil(arg.data[6]),
                            std::ceil(arg.data[7]),
                            std::ceil(arg.data[8]),
                            std::ceil(arg.data[9]),
                            std::ceil(arg.data[10]),
                            std::ceil(arg.data[11]),
                            std::ceil(arg.data[12]),
                            std::ceil(arg.data[13]),
                            std::ceil(arg.data[14]),
                            std::ceil(arg.data[15]));
#endif
  }

  inline vtype<float, 16> floor(vtype<float, 16> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_roundscale_ps(arg.reg, _MM_FROUND_TO_NEG_INF|_MM_FROUND_NO_EXC));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_floor_ps(arg.reg.x),
                            _mm256_floor_ps(arg.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(std::floor(arg.data[0]),
                            std::floor(arg.data[1]),
                            std::floor(arg.data[2]),
                            std::floor(arg.data[3]),
                            std::floor(arg.data[4]),
                            std::floor(arg.data[5]),
                            std::floor(arg.data[6]),
                            std::floor(arg.data[7]),
                            std::floor(arg.data[8]),
                            std::floor(arg.data[9]),
                            std::floor(arg.data[10]),
                            std::floor(arg.data[11]),
                            std::floor(arg.data[12]),
                            std::floor(arg.data[13]),
                            std::floor(arg.data[14]),
                            std::floor(arg.data[15]));
#endif
  }

  inline vtype<float, 16> round(vtype<float, 16> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_roundscale_ps(arg.reg, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_round_ps(arg.reg.x, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC),
                            _mm256_round_ps(arg.reg.y, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(std::round(arg.data[0]),
                            std::round(arg.data[1]),
                            std::round(arg.data[2]),
                            std::round(arg.data[3]),
                            std::round(arg.data[4]),
                            std::round(arg.data[5]),
                            std::round(arg.data[6]),
                            std::round(arg.data[7]),
                            std::round(arg.data[8]),
                            std::round(arg.data[9]),
                            std::round(arg.data[10]),
                            std::round(arg.data[11]),
                            std::round(arg.data[12]),
                            std::round(arg.data[13]),
                            std::round(arg.data[14]),
                            std::round(arg.data[15]));
#endif
  }

  inline vtype<float, 16> sqrt(vtype<float, 16> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_sqrt_ps(arg.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_sqrt_ps(arg.reg.x),
                            _mm256_sqrt_ps(arg.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(std::sqrt(arg.data[0]),
                            std::sqrt(arg.data[1]),
                            std::sqrt(arg.data[2]),
                            std::sqrt(arg.data[3]),
                            std::sqrt(arg.data[4]),
                            std::sqrt(arg.data[5]),
                            std::sqrt(arg.data[6]),
                            std::sqrt(arg.data[7]),
                            std::sqrt(arg.data[8]),
                            std::sqrt(arg.data[9]),
                            std::sqrt(arg.data[10]),
                            std::sqrt(arg.data[11]),
                            std::sqrt(arg.data[12]),
                            std::sqrt(arg.data[13]),
                            std::sqrt(arg.data[14]),
                            std::sqrt(arg.data[15]));
#endif
  }

  inline vtype<float, 16> rsqrt(vtype<float, 16> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<float, 16>(_mm512_rsqrt14_ps(arg.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<float, 16>(_mm256_rsqrt_ps(arg.reg.x),
                            _mm256_rsqrt_ps(arg.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<float, 16>(1.0f/std::sqrt(arg.data[0]),
                            1.0f/std::sqrt(arg.data[1]),
                            1.0f/std::sqrt(arg.data[2]),
                            1.0f/std::sqrt(arg.data[3]),
                            1.0f/std::sqrt(arg.data[4]),
                            1.0f/std::sqrt(arg.data[5]),
                            1.0f/std::sqrt(arg.data[6]),
                            1.0f/std::sqrt(arg.data[7]),
                            1.0f/std::sqrt(arg.data[8]),
                            1.0f/std::sqrt(arg.data[9]),
                            1.0f/std::sqrt(arg.data[10]),
                            1.0f/std::sqrt(arg.data[11]),
                            1.0f/std::sqrt(arg.data[12]),
                            1.0f/std::sqrt(arg.data[13]),
                            1.0f/std::sqrt(arg.data[14]),
                            1.0f/std::sqrt(arg.data[15]));
#endif
  }

  /**
   * Vectorized 4x double.
   */
  template<>
  struct vtype<double, 4> {
    static constexpr unsigned int W = 4;

    typedef double scalar_t;
    typedef vmask<sizeof(double), 4> mask_t;
    typedef StaticSizedArray<double, 4> value_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef __m256d register_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    typedef StaticSizedArray<double, 4> register_t;
#endif

    union {
      value_t data;
      register_t reg;
    };

    vtype() = default ;

    vtype(double e0, double e1, double e2, double e3)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm256_set_pd(e3, e2, e1, e0)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    : data(e0, e1, e2, e3) { }
#endif

    explicit vtype(double e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm256_set1_pd(e)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "serial implementation"
# endif
    : data(e) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    explicit vtype(__m256d reg) : reg(reg) { }
#endif

    explicit vtype(double const * mem)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(_mm256_loadu_pd(mem)) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(mem) { }
#endif

    vtype(vtype<double, 4> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    : data(rhs.data) { }
#endif

    vtype<double, 4> & operator=(vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
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
    reg = _mm256_add_pd(reg, rhs.reg);
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
    reg = _mm256_sub_pd(reg, rhs.reg);
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
    reg = _mm256_mul_pd(reg, rhs.reg);
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
    reg = _mm256_div_pd(reg, rhs.reg);
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
    return vtype<double, 4>(_mm256_add_pd(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(lhs.data[0] + rhs.data[0],
                            lhs.data[1] + rhs.data[1],
                            lhs.data[2] + rhs.data[2],
                            lhs.data[3] + rhs.data[3]);
#endif
  }

  inline vtype<double, 4> operator-(vtype<double, 4> const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_sub_pd(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(lhs.data[0] - rhs.data[0],
                            lhs.data[1] - rhs.data[1],
                            lhs.data[2] - rhs.data[2],
                            lhs.data[3] - rhs.data[3]);
#endif
  }

  inline vtype<double, 4> operator*(vtype<double, 4> const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_mul_pd(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(lhs.data[0] * rhs.data[0],
                            lhs.data[1] * rhs.data[1],
                            lhs.data[2] * rhs.data[2],
                            lhs.data[3] * rhs.data[3]);
#endif
  }

  inline vtype<double, 4> operator/(vtype<double, 4> const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_div_pd(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(lhs.data[0] / rhs.data[0],
                            lhs.data[1] / rhs.data[1],
                            lhs.data[2] / rhs.data[2],
                            lhs.data[3] / rhs.data[3]);
#endif
  }

  inline vtype<double, 4> operator-(vtype<double, 4> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256d sign_mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000));
    return vtype<double, 4>(_mm256_xor_pd(v.reg, sign_mask)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(-v.data[0],
                            -v.data[1],
                            -v.data[2],
                            -v.data[3]);
#endif
  }

  // boolean comparison
  inline vmask<sizeof(double), 4> operator<(vtype<double, 4> const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(double), 4>(_mm256_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_LT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 4>((__mmask8)_mm512_cmp_pd_mask(
      _mm512_castpd256_pd512(lhs.reg), _mm512_castpd256_pd512(rhs.reg), _CMP_LT_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 4>(_mm256_cmp_pd(lhs.reg, rhs.reg, _CMP_LT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 4>(lhs.data[0] < rhs.data[0],
                                    lhs.data[1] < rhs.data[1],
                                    lhs.data[2] < rhs.data[2],
                                    lhs.data[3] < rhs.data[3]);
#endif
  }

  inline vmask<sizeof(double), 4> operator<=(vtype<double, 4>  const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(double), 4>(_mm256_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_LE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 4>((__mmask8)_mm512_cmp_pd_mask(
      _mm512_castpd256_pd512(lhs.reg), _mm512_castpd256_pd512(rhs.reg), _CMP_LE_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 4>(_mm256_cmp_pd(lhs.reg, rhs.reg, _CMP_LE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 4>(lhs.data[0] <= rhs.data[0],
                                    lhs.data[1] <= rhs.data[1],
                                    lhs.data[2] <= rhs.data[2],
                                    lhs.data[3] <= rhs.data[3]);
#endif
  }

  inline vmask<sizeof(double), 4> operator>(vtype<double,4>  const & lhs, vtype<double,4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(double), 4>(_mm256_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_GT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 4>((__mmask8)_mm512_cmp_pd_mask(
      _mm512_castpd256_pd512(lhs.reg), _mm512_castpd256_pd512(rhs.reg), _CMP_GT_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double),4>(_mm256_cmp_pd(lhs.reg, rhs.reg, _CMP_GT_OQ)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 4>(lhs.data[0] > rhs.data[0],
                                    lhs.data[1] > rhs.data[1],
                                    lhs.data[2] > rhs.data[2],
                                    lhs.data[3] > rhs.data[3]) ;
#endif
  }

  inline vmask<sizeof(double), 4> operator>=(vtype<double,4>  const & lhs, vtype<double,4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(double), 4>(_mm256_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_GE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 4>((__mmask8)_mm512_cmp_pd_mask(
      _mm512_castpd256_pd512(lhs.reg), _mm512_castpd256_pd512(rhs.reg), _CMP_GE_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double),4>(_mm256_cmp_pd(lhs.reg, rhs.reg, _CMP_GE_OQ)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 4>(lhs.data[0] >= rhs.data[0],
                                    lhs.data[1] >= rhs.data[1],
                                    lhs.data[2] >= rhs.data[2],
                                    lhs.data[3] >= rhs.data[3]) ;
#endif
  }

  inline vmask<sizeof(double), 4> operator==(vtype<double, 4>  const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(double), 4>(_mm256_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_EQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 4>((__mmask8)_mm512_cmp_pd_mask(
      _mm512_castpd256_pd512(lhs.reg), _mm512_castpd256_pd512(rhs.reg), _CMP_EQ_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double),4>(_mm256_cmp_pd(lhs.reg, rhs.reg, _CMP_EQ_OQ)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 4>(lhs.data[0] == rhs.data[0],
                                    lhs.data[1] == rhs.data[1],
                                    lhs.data[2] == rhs.data[2],
                                    lhs.data[3] == rhs.data[3]) ;
#endif
  }

  inline vmask<sizeof(double), 4> operator!=(vtype<double, 4>  const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vmask<sizeof(double), 4>(_mm256_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_NEQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 4>((__mmask8)_mm512_cmp_pd_mask(
      _mm512_castpd256_pd512(lhs.reg), _mm512_castpd256_pd512(rhs.reg), _CMP_NEQ_OQ
    ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 4>(_mm256_cmp_pd(lhs.reg, rhs.reg, _CMP_NEQ_OQ)) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 4>(lhs.data[0] != rhs.data[0],
                                    lhs.data[1] != rhs.data[1],
                                    lhs.data[2] != rhs.data[2],
                                    lhs.data[3] != rhs.data[3]) ;
#endif
  }

  inline void mem_store(double * mem, vtype<double, 4> const & op) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    _mm256_storeu_pd(mem, op.reg);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    mem[0] = op.data[0];
    mem[1] = op.data[1];
    mem[2] = op.data[2];
    mem[3] = op.data[3];
#endif
  }

  inline vtype<double, 4> select(vmask<sizeof(double), 4> const & c,
                                 vtype<double, 4> const & t,
                                 vtype<double, 4> const & f) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vtype<double, 4>(_mm256_mask_blend_pd(c.reg, f.reg, t.reg)) ;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 4>(
      _mm512_castpd512_pd256(_mm512_mask_blend_pd(
        c.reg,
        _mm512_castpd256_pd512(f.reg),
        _mm512_castpd256_pd512(t.reg)
      ))
    );
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_blendv_pd(f.reg, t.reg, _mm256_castsi256_pd(c.reg))) ;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>((c.data & 1u) ? t.data[0] : f.data[0],
                            (c.data & 2u) ? t.data[1] : f.data[1],
                            (c.data & 4u) ? t.data[2] : f.data[2],
                            (c.data & 8u) ? t.data[3] : f.data[3]) ;
#endif
  }

  inline vtype<double,4> fabs(vtype<double,4> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256d sign_mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
    return vtype<double, 4>(_mm256_and_pd(v.reg, sign_mask)) ;
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

  inline vtype<double,4> abs(vtype<double,4> const & v) {
    return fabs(v);
  }

  inline vtype<double,4> max(vtype<double, 4> const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_max_pd(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(std::max(lhs.data[0], rhs.data[0]),
                            std::max(lhs.data[1], rhs.data[1]),
                            std::max(lhs.data[2], rhs.data[2]),
                            std::max(lhs.data[3], rhs.data[3]));
#endif
  }

  inline vtype<double, 4> min(vtype<double, 4> const & lhs, vtype<double, 4> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_min_pd(lhs.reg, rhs.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(std::min(lhs.data[0], rhs.data[0]),
                            std::min(lhs.data[1], rhs.data[1]),
                            std::min(lhs.data[2], rhs.data[2]),
                            std::min(lhs.data[3], rhs.data[3]));
#endif
  }

  inline vtype<double, 4> ceil(vtype<double, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_ceil_pd(arg.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(std::ceil(arg.data[0]),
                            std::ceil(arg.data[1]),
                            std::ceil(arg.data[2]),
                            std::ceil(arg.data[3]));
#endif
  }

  inline vtype<double, 4> floor(vtype<double, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_floor_pd(arg.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(std::floor(arg.data[0]),
                            std::floor(arg.data[1]),
                            std::floor(arg.data[2]),
                            std::floor(arg.data[3]));
#endif
  }

  inline vtype<double, 4> round(vtype<double, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_round_pd(arg.reg, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(std::round(arg.data[0]),
                            std::round(arg.data[1]),
                            std::round(arg.data[2]),
                            std::round(arg.data[3]));
#endif
  }

  inline vtype<double, 4> sqrt(vtype<double, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 4>(_mm256_sqrt_pd(arg.reg));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(std::sqrt(arg.data[0]),
                            std::sqrt(arg.data[1]),
                            std::sqrt(arg.data[2]),
                            std::sqrt(arg.data[3]));
#endif
  }

  inline vtype<double, 4> rsqrt(vtype<double, 4> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__) && defined(__AVX512VL__)
    return vtype<double, 4>(_mm256_rsqrt14_pd(arg.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    // TODO: This might have high latency. Need to find better approach.
    __m256d sqrtx = _mm256_sqrt_pd(arg.reg);
    __m256d one = _mm256_set1_pd(1.0);
    return vtype<double, 4>(_mm256_div_pd(one, sqrtx));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 4>(1.0/std::sqrt(arg.data[0]),
                            1.0/std::sqrt(arg.data[1]),
                            1.0/std::sqrt(arg.data[2]),
                            1.0/std::sqrt(arg.data[3]));
#endif
  }

  /**
   * Vectorized 8x double
   */
  template<>
  struct vtype<double, 8> {
    static constexpr unsigned int W = 8;

    typedef double scalar_t;
    typedef vmask<sizeof(double), 8> mask_t;
    typedef StaticSizedArray<double, 8> value_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    typedef __m512d register_t;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef struct {
      __m256d x, y;
    } register_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    typedef StaticSizedArray<double, 8> register_t;
#endif

    union {
      value_t data;
      register_t reg;
    };

    vtype() = default;

    vtype(double e0, double e1, double e2, double e3,
          double e4, double e5, double e6, double e7)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(_mm512_set_pd(e7, e6, e5, e4, e3, e2, e1, e0)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_set_pd(e3, e2, e1, e0), _mm256_set_pd(e7, e6, e5, e4)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(e0, e1, e2, e3, e4, e5, e6, e7) { }
#endif

    explicit vtype(double e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(_mm512_set1_pd(e)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_set1_pd(e), _mm256_set1_pd(e)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(e) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__)
    explicit vtype(__m256d reg0, __m256d reg1) : reg{reg0, reg1} { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    explicit vtype(__m512d reg) : reg(reg) { }
#endif

    explicit vtype(double const * mem)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(_mm512_loadu_pd(mem)) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_loadu_pd(mem), _mm256_loadu_pd(mem+4)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(mem) { }
#endif

    vtype(vtype<double, 8> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(rhs.reg) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(rhs.data) { }
#endif

    vtype<double, 8> & operator=(vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg = rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
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

    vtype<double, 8> & operator+=(vtype<double, 8> const & rhs);
    vtype<double, 8> & operator-=(vtype<double, 8> const & rhs);
    vtype<double, 8> & operator*=(vtype<double, 8> const & rhs);
    vtype<double, 8> & operator/=(vtype<double, 8> const & rhs);
  };

  inline vtype<double, 8> & vtype<double, 8>::operator+=(vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg = _mm512_add_pd(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_add_pd(reg.x, rhs.reg.x);
    reg.y = _mm256_add_pd(reg.y, rhs.reg.y);
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

  inline vtype<double, 8> & vtype<double, 8>::operator-=(vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg = _mm512_sub_pd(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_sub_pd(reg.x, rhs.reg.x);
    reg.y = _mm256_sub_pd(reg.y, rhs.reg.y);
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
  
  inline vtype<double, 8> & vtype<double, 8>::operator*=(vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg = _mm512_mul_pd(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_mul_pd(reg.x, rhs.reg.x);
    reg.y = _mm256_mul_pd(reg.y, rhs.reg.y);
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
  
  inline vtype<double, 8> & vtype<double, 8>::operator/=(vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg = _mm512_div_pd(reg, rhs.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_div_pd(reg.x, rhs.reg.x);
    reg.y = _mm256_div_pd(reg.y, rhs.reg.y);
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

  inline vtype<double, 8> operator+(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_add_pd(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_add_pd(lhs.reg.x, rhs.reg.x),
                            _mm256_add_pd(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(lhs.data[0] + rhs.data[0],
                            lhs.data[1] + rhs.data[1],
                            lhs.data[2] + rhs.data[2],
                            lhs.data[3] + rhs.data[3],
                            lhs.data[4] + rhs.data[4],
                            lhs.data[5] + rhs.data[5],
                            lhs.data[6] + rhs.data[6],
                            lhs.data[7] + rhs.data[7]);
#endif
  }

  inline vtype<double, 8> operator-(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_sub_pd(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_sub_pd(lhs.reg.x, rhs.reg.x),
                            _mm256_sub_pd(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(lhs.data[0] - rhs.data[0],
                            lhs.data[1] - rhs.data[1],
                            lhs.data[2] - rhs.data[2],
                            lhs.data[3] - rhs.data[3],
                            lhs.data[4] - rhs.data[4],
                            lhs.data[5] - rhs.data[5],
                            lhs.data[6] - rhs.data[6],
                            lhs.data[7] - rhs.data[7]);
#endif
  }

  inline vtype<double, 8> operator*(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_mul_pd(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_mul_pd(lhs.reg.x, rhs.reg.x),
                            _mm256_mul_pd(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(lhs.data[0] * rhs.data[0],
                            lhs.data[1] * rhs.data[1],
                            lhs.data[2] * rhs.data[2],
                            lhs.data[3] * rhs.data[3],
                            lhs.data[4] * rhs.data[4],
                            lhs.data[5] * rhs.data[5],
                            lhs.data[6] * rhs.data[6],
                            lhs.data[7] * rhs.data[7]);
#endif
  }

  inline vtype<double, 8> operator/(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_div_pd(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_div_pd(lhs.reg.x, rhs.reg.x),
                            _mm256_div_pd(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(lhs.data[0] / rhs.data[0],
                            lhs.data[1] / rhs.data[1],
                            lhs.data[2] / rhs.data[2],
                            lhs.data[3] / rhs.data[3],
                            lhs.data[4] / rhs.data[4],
                            lhs.data[5] / rhs.data[5],
                            lhs.data[6] / rhs.data[6],
                            lhs.data[7] / rhs.data[7]);
#endif
  }

  inline vtype<double, 8> operator-(vtype<double, 8> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    __m512d sign_mask = _mm512_castsi512_pd(_mm512_set1_epi64(0x8000000000000000));
    return vtype<double, 8>(_mm512_xor_pd(v.reg, sign_mask));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    __m512i sign_mask = _mm512_set1_epi64(0x8000000000000000);
    return vtype<double, 8>(_mm512_castsi512_pd(_mm512_xor_epi64(_mm512_castpd_si512(v.reg), sign_mask)));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256d sign_mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000));
    return vtype<double, 8>(_mm256_xor_pd(v.reg.x, sign_mask),
                            _mm256_xor_pd(v.reg.y, sign_mask));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(-v.data[0],
                            -v.data[1],
                            -v.data[2],
                            -v.data[3],
                            -v.data[4],
                            -v.data[5],
                            -v.data[6],
                            -v.data[7]);
#endif
  }

  inline vmask<sizeof(double), 8> operator<(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 8>(_mm512_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_LT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 8>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_LT_OQ),
                                    _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_LT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 8>(lhs.data[0] < rhs.data[0],
                                    lhs.data[1] < rhs.data[1],
                                    lhs.data[2] < rhs.data[2],
                                    lhs.data[3] < rhs.data[3],
                                    lhs.data[4] < rhs.data[4],
                                    lhs.data[5] < rhs.data[5],
                                    lhs.data[6] < rhs.data[6],
                                    lhs.data[7] < rhs.data[7]);
#endif
  }

  inline vmask<sizeof(double), 8> operator<=(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 8>(_mm512_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_LE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 8>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_LE_OQ),
                                    _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_LE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 8>(lhs.data[0] <= rhs.data[0],
                                    lhs.data[1] <= rhs.data[1],
                                    lhs.data[2] <= rhs.data[2],
                                    lhs.data[3] <= rhs.data[3],
                                    lhs.data[4] <= rhs.data[4],
                                    lhs.data[5] <= rhs.data[5],
                                    lhs.data[6] <= rhs.data[6],
                                    lhs.data[7] <= rhs.data[7]);
#endif
  }

  inline vmask<sizeof(double), 8> operator>(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 8>(_mm512_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_GT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 8>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_GT_OQ),
                                    _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_GT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 8>(lhs.data[0] > rhs.data[0],
                                    lhs.data[1] > rhs.data[1],
                                    lhs.data[2] > rhs.data[2],
                                    lhs.data[3] > rhs.data[3],
                                    lhs.data[4] > rhs.data[4],
                                    lhs.data[5] > rhs.data[5],
                                    lhs.data[6] > rhs.data[6],
                                    lhs.data[7] > rhs.data[7]);
#endif
  }

  inline vmask<sizeof(double), 8> operator>=(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 8>(_mm512_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_GE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 8>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_GE_OQ),
                                    _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_GE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 8>(lhs.data[0] >= rhs.data[0],
                                    lhs.data[1] >= rhs.data[1],
                                    lhs.data[2] >= rhs.data[2],
                                    lhs.data[3] >= rhs.data[3],
                                    lhs.data[4] >= rhs.data[4],
                                    lhs.data[5] >= rhs.data[5],
                                    lhs.data[6] >= rhs.data[6],
                                    lhs.data[7] >= rhs.data[7]);
#endif
  }

  inline vmask<sizeof(double), 8> operator==(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 8>(_mm512_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_EQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 8>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_EQ_OQ),
                                    _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_EQ_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 8>(lhs.data[0] == rhs.data[0],
                                    lhs.data[1] == rhs.data[1],
                                    lhs.data[2] == rhs.data[2],
                                    lhs.data[3] == rhs.data[3],
                                    lhs.data[4] == rhs.data[4],
                                    lhs.data[5] == rhs.data[5],
                                    lhs.data[6] == rhs.data[6],
                                    lhs.data[7] == rhs.data[7]);
#endif
  }

  inline vmask<sizeof(double), 8> operator!=(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 8>(_mm512_cmp_pd_mask(lhs.reg, rhs.reg, _CMP_NEQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 8>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_NEQ_OQ),
                                    _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_NEQ_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 8>(lhs.data[0] != rhs.data[0],
                                    lhs.data[1] != rhs.data[1],
                                    lhs.data[2] != rhs.data[2],
                                    lhs.data[3] != rhs.data[3],
                                    lhs.data[4] != rhs.data[4],
                                    lhs.data[5] != rhs.data[5],
                                    lhs.data[6] != rhs.data[6],
                                    lhs.data[7] != rhs.data[7]);
#endif
  }

  inline void mem_store(double * mem, vtype<double, 8> const & op) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    _mm512_storeu_pd(mem, op.reg);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    _mm256_storeu_pd(mem, op.reg.x);
    _mm256_storeu_pd(mem+4, op.reg.y);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    mem[0] = op.data[0];
    mem[1] = op.data[1];
    mem[2] = op.data[2];
    mem[3] = op.data[3];
    mem[4] = op.data[4];
    mem[5] = op.data[5];
    mem[6] = op.data[6];
    mem[7] = op.data[7];
#endif
  }

  inline vtype<double, 8> select(vmask<sizeof(double), 8> const & c,
                                 vtype<double, 8> const & t,
                                 vtype<double, 8> const & f) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_mask_blend_pd(c.reg, f.reg, t.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_blendv_pd(f.reg.x, t.reg.x, _mm256_castsi256_pd(c.reg.x)),
                            _mm256_blendv_pd(f.reg.y, t.reg.y, _mm256_castsi256_pd(c.reg.y)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>((c.data & 1u) ? t.data[0] : f.data[0],
                            (c.data & 2u) ? t.data[1] : f.data[1],
                            (c.data & 4u) ? t.data[2] : f.data[2],
                            (c.data & 8u) ? t.data[3] : f.data[3],
                            (c.data & 16u) ? t.data[4] : f.data[4],
                            (c.data & 32u) ? t.data[5] : f.data[5],
                            (c.data & 64u) ? t.data[6] : f.data[6],
                            (c.data & 128u) ? t.data[7] : f.data[7]);
#endif
  }

  inline vtype<double, 8> fabs(vtype<double, 8> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    __m512d sign_mask = _mm512_castsi512_pd(_mm512_set1_epi64(0x7FFFFFFFFFFFFFFF));
    return vtype<double, 8>(_mm512_and_pd(v.reg, sign_mask));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    __m512i sign_mask = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFF);
    return vtype<double, 8>(_mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(v.reg), sign_mask)));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256d sign_mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
    return vtype<double, 8>(_mm256_and_pd(v.reg.x, sign_mask),
                            _mm256_and_pd(v.reg.y, sign_mask));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(std::fabs(v.data[0]),
                            std::fabs(v.data[1]),
                            std::fabs(v.data[2]),
                            std::fabs(v.data[3]),
                            std::fabs(v.data[4]),
                            std::fabs(v.data[5]),
                            std::fabs(v.data[6]),
                            std::fabs(v.data[7]));
#endif
  }

  inline vtype<double, 8> abs(vtype<double, 8> const & v) {
    return fabs(v);
  }

  inline vtype<double, 8> max(vtype<double, 8>  const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_max_pd(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_max_pd(lhs.reg.x, rhs.reg.x),
                            _mm256_max_pd(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(std::max(lhs.data[0], rhs.data[0]),
                            std::max(lhs.data[1], rhs.data[1]),
                            std::max(lhs.data[2], rhs.data[2]),
                            std::max(lhs.data[3], rhs.data[3]),
                            std::max(lhs.data[4], rhs.data[4]),
                            std::max(lhs.data[5], rhs.data[5]),
                            std::max(lhs.data[6], rhs.data[6]),
                            std::max(lhs.data[7], rhs.data[7]));
#endif
  }

  inline vtype<double, 8> min(vtype<double, 8> const & lhs, vtype<double, 8> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_min_pd(lhs.reg, rhs.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_min_pd(lhs.reg.x, rhs.reg.x),
                            _mm256_min_pd(lhs.reg.y, rhs.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(std::min(lhs.data[0], rhs.data[0]),
                            std::min(lhs.data[1], rhs.data[1]),
                            std::min(lhs.data[2], rhs.data[2]),
                            std::min(lhs.data[3], rhs.data[3]),
                            std::min(lhs.data[4], rhs.data[4]),
                            std::min(lhs.data[5], rhs.data[5]),
                            std::min(lhs.data[6], rhs.data[6]),
                            std::min(lhs.data[7], rhs.data[7]));
#endif
  }

  inline vtype<double, 8> ceil(vtype<double, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_roundscale_pd(arg.reg, _MM_FROUND_TO_POS_INF|_MM_FROUND_NO_EXC));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_ceil_pd(arg.reg.x),
                            _mm256_ceil_pd(arg.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(std::ceil(arg.data[0]),
                            std::ceil(arg.data[1]),
                            std::ceil(arg.data[2]),
                            std::ceil(arg.data[3]),
                            std::ceil(arg.data[4]),
                            std::ceil(arg.data[5]),
                            std::ceil(arg.data[6]),
                            std::ceil(arg.data[7]));
#endif
  }

  inline vtype<double, 8> floor(vtype<double, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_roundscale_pd(arg.reg, _MM_FROUND_TO_NEG_INF|_MM_FROUND_NO_EXC));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_floor_pd(arg.reg.x),
                            _mm256_floor_pd(arg.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(std::floor(arg.data[0]),
                            std::floor(arg.data[1]),
                            std::floor(arg.data[2]),
                            std::floor(arg.data[3]),
                            std::floor(arg.data[4]),
                            std::floor(arg.data[5]),
                            std::floor(arg.data[6]),
                            std::floor(arg.data[7]));
#endif
  }

  inline vtype<double, 8> round(vtype<double, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_roundscale_pd(arg.reg, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_round_pd(arg.reg.x, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC),
                            _mm256_round_pd(arg.reg.y, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(std::round(arg.data[0]),
                            std::round(arg.data[1]),
                            std::round(arg.data[2]),
                            std::round(arg.data[3]),
                            std::round(arg.data[4]),
                            std::round(arg.data[5]),
                            std::round(arg.data[6]),
                            std::round(arg.data[7]));
#endif
  }

  inline vtype<double, 8> sqrt(vtype<double, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_sqrt_pd(arg.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 8>(_mm256_sqrt_pd(arg.reg.x),
                            _mm256_sqrt_pd(arg.reg.y));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(std::sqrt(arg.data[0]),
                            std::sqrt(arg.data[1]),
                            std::sqrt(arg.data[2]),
                            std::sqrt(arg.data[3]),
                            std::sqrt(arg.data[4]),
                            std::sqrt(arg.data[5]),
                            std::sqrt(arg.data[6]),
                            std::sqrt(arg.data[7]));
#endif
  }

  inline vtype<double, 8> rsqrt(vtype<double, 8> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 8>(_mm512_rsqrt14_pd(arg.reg));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    // TODO: This might have high latency. Need to find better approach.
    __m256d sqrtx = _mm256_sqrt_pd(arg.reg.x);
    __m256d sqrty = _mm256_sqrt_pd(arg.reg.y);
    __m256d one = _mm256_set1_pd(1.0);
    return vtype<double, 8>(_mm256_div_pd(one, sqrtx),
                            _mm256_div_pd(one, sqrty));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 8>(1.0/std::sqrt(arg.data[0]),
                            1.0/std::sqrt(arg.data[1]),
                            1.0/std::sqrt(arg.data[2]),
                            1.0/std::sqrt(arg.data[3]),
                            1.0/std::sqrt(arg.data[4]),
                            1.0/std::sqrt(arg.data[5]),
                            1.0/std::sqrt(arg.data[6]),
                            1.0/std::sqrt(arg.data[7]));
#endif
  }

  /**
   * Vectorized 16x double
   */
  template<>
  struct vtype<double, 16> {
    static constexpr unsigned int W = 16;

    typedef double scalar_t;
    typedef vmask<sizeof(double), 16> mask_t;
    typedef StaticSizedArray<double, 16> value_t;

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    typedef struct {
      __m512d x, y;
    } register_t;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    typedef struct {
      __m256d x, y, z, w;
    } register_t;
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    typedef StaticSizedArray<double, 16> register_t;
#endif

    union {
      value_t data;
      register_t reg;
    };

    vtype() = default;

    vtype(double e0, double e1, double e2, double e3,
          double e4, double e5, double e6, double e7,
          double e8, double e9, double e10, double e11,
          double e12, double e13, double e14, double e15)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg{_mm512_set_pd(e7, e6, e5, e4, e3, e2, e1, e0),
            _mm512_set_pd(e15, e14, e13, e12, e11, e10, e9, e8)} { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_set_pd(e3, e2, e1, e0),
            _mm256_set_pd(e7, e6, e5, e4),
            _mm256_set_pd(e11, e10, e9, e8),
            _mm256_set_pd(e15, e14, e13, e12)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(e0, e1, e2, e3, e4, e5, e6, e7,
             e8, e9, e10, e11, e12, e13, e14, e15) { }
#endif

    explicit vtype(double e)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg{_mm512_set1_pd(e), _mm512_set1_pd(e)} { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_set1_pd(e), _mm256_set1_pd(e),
               _mm256_set1_pd(e), _mm256_set1_pd(e)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(e) { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__) && !defined(__AVX512F__)
    explicit vtype(__m256d reg0, __m256d reg1, __m256d reg2, __m256d reg3)
      : reg{reg0, reg1, reg2, reg3} { }
#endif

#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    explicit vtype(__m512d reg0, __m512d reg1) : reg{reg0, reg1} { }
#endif

    explicit vtype(double const * mem)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg{_mm512_loadu_pd(mem), _mm512_loadu_pd(mem+8)} { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg{_mm256_loadu_pd(mem), _mm256_loadu_pd(mem+4),
            _mm256_loadu_pd(mem+8), _mm256_loadu_pd(mem+12)} { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(mem) { }
#endif

    vtype(vtype<double, 16> const & rhs)
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      : reg(rhs.reg) { }
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      : reg(rhs.reg) { }
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
      : data(rhs.data) { }
#endif

    vtype<double, 16> & operator=(vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
      reg = rhs.reg;
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
      reg = rhs.reg;
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

    vtype<double, 16> & operator+=(vtype<double, 16> const & rhs);
    vtype<double, 16> & operator-=(vtype<double, 16> const & rhs);
    vtype<double, 16> & operator*=(vtype<double, 16> const & rhs);
    vtype<double, 16> & operator/=(vtype<double, 16> const & rhs);
  };

  inline vtype<double, 16> & vtype<double, 16>::operator+=(vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg.x = _mm512_add_pd(reg.x, rhs.reg.x);
    reg.y = _mm512_add_pd(reg.y, rhs.reg.y);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_add_pd(reg.x, rhs.reg.x);
    reg.y = _mm256_add_pd(reg.y, rhs.reg.y);
    reg.z = _mm256_add_pd(reg.z, rhs.reg.z);
    reg.w = _mm256_add_pd(reg.w, rhs.reg.w);
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
    data[8] += rhs.data[8];
    data[9] += rhs.data[9];
    data[10] += rhs.data[10];
    data[11] += rhs.data[11];
    data[12] += rhs.data[12];
    data[13] += rhs.data[13];
    data[14] += rhs.data[14];
    data[15] += rhs.data[15];
#endif
    return *this;
  }

  inline vtype<double, 16> & vtype<double, 16>::operator-=(vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg.x = _mm512_sub_pd(reg.x, rhs.reg.x);
    reg.y = _mm512_sub_pd(reg.y, rhs.reg.y);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_sub_pd(reg.x, rhs.reg.x);
    reg.y = _mm256_sub_pd(reg.y, rhs.reg.y);
    reg.z = _mm256_sub_pd(reg.z, rhs.reg.z);
    reg.w = _mm256_sub_pd(reg.w, rhs.reg.w);
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
    data[8] -= rhs.data[8];
    data[9] -= rhs.data[9];
    data[10] -= rhs.data[10];
    data[11] -= rhs.data[11];
    data[12] -= rhs.data[12];
    data[13] -= rhs.data[13];
    data[14] -= rhs.data[14];
    data[15] -= rhs.data[15];
#endif
    return *this;
  }

  inline vtype<double, 16> & vtype<double, 16>::operator*=(vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg.x = _mm512_mul_pd(reg.x, rhs.reg.x);
    reg.y = _mm512_mul_pd(reg.y, rhs.reg.y);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_mul_pd(reg.x, rhs.reg.x);
    reg.y = _mm256_mul_pd(reg.y, rhs.reg.y);
    reg.z = _mm256_mul_pd(reg.z, rhs.reg.z);
    reg.w = _mm256_mul_pd(reg.w, rhs.reg.w);
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
    data[8] *= rhs.data[8];
    data[9] *= rhs.data[9];
    data[10] *= rhs.data[10];
    data[11] *= rhs.data[11];
    data[12] *= rhs.data[12];
    data[13] *= rhs.data[13];
    data[14] *= rhs.data[14];
    data[15] *= rhs.data[15];
#endif
    return *this;
  }

  inline vtype<double, 16> & vtype<double, 16>::operator/=(vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    reg.x = _mm512_div_pd(reg.x, rhs.reg.x);
    reg.y = _mm512_div_pd(reg.y, rhs.reg.y);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    reg.x = _mm256_div_pd(reg.x, rhs.reg.x);
    reg.y = _mm256_div_pd(reg.y, rhs.reg.y);
    reg.z = _mm256_div_pd(reg.z, rhs.reg.z);
    reg.w = _mm256_div_pd(reg.w, rhs.reg.w);
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
    data[8] /= rhs.data[8];
    data[9] /= rhs.data[9];
    data[10] /= rhs.data[10];
    data[11] /= rhs.data[11];
    data[12] /= rhs.data[12];
    data[13] /= rhs.data[13];
    data[14] /= rhs.data[14];
    data[15] /= rhs.data[15];
#endif
    return *this;
  }

  inline vtype<double, 16> operator+(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_add_pd(lhs.reg.x, rhs.reg.x),
                             _mm512_add_pd(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_add_pd(lhs.reg.x, rhs.reg.x),
                             _mm256_add_pd(lhs.reg.y, rhs.reg.y),
                             _mm256_add_pd(lhs.reg.z, rhs.reg.z),
                             _mm256_add_pd(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(lhs.data[0] + rhs.data[0],
                             lhs.data[1] + rhs.data[1],
                             lhs.data[2] + rhs.data[2],
                             lhs.data[3] + rhs.data[3],
                             lhs.data[4] + rhs.data[4],
                             lhs.data[5] + rhs.data[5],
                             lhs.data[6] + rhs.data[6],
                             lhs.data[7] + rhs.data[7],
                             lhs.data[8] + rhs.data[8],
                             lhs.data[9] + rhs.data[9],
                             lhs.data[10] + rhs.data[10],
                             lhs.data[11] + rhs.data[11],
                             lhs.data[12] + rhs.data[12],
                             lhs.data[13] + rhs.data[13],
                             lhs.data[14] + rhs.data[14],
                             lhs.data[15] + rhs.data[15]);
#endif
  }

  inline vtype<double, 16> operator-(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_sub_pd(lhs.reg.x, rhs.reg.x),
                             _mm512_sub_pd(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_sub_pd(lhs.reg.x, rhs.reg.x),
                             _mm256_sub_pd(lhs.reg.y, rhs.reg.y),
                             _mm256_sub_pd(lhs.reg.z, rhs.reg.z),
                             _mm256_sub_pd(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(lhs.data[0] - rhs.data[0],
                             lhs.data[1] - rhs.data[1],
                             lhs.data[2] - rhs.data[2],
                             lhs.data[3] - rhs.data[3],
                             lhs.data[4] - rhs.data[4],
                             lhs.data[5] - rhs.data[5],
                             lhs.data[6] - rhs.data[6],
                             lhs.data[7] - rhs.data[7],
                             lhs.data[8] - rhs.data[8],
                             lhs.data[9] - rhs.data[9],
                             lhs.data[10] - rhs.data[10],
                             lhs.data[11] - rhs.data[11],
                             lhs.data[12] - rhs.data[12],
                             lhs.data[13] - rhs.data[13],
                             lhs.data[14] - rhs.data[14],
                             lhs.data[15] - rhs.data[15]);
#endif
  }

  inline vtype<double, 16> operator*(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_mul_pd(lhs.reg.x, rhs.reg.x),
                             _mm512_mul_pd(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_mul_pd(lhs.reg.x, rhs.reg.x),
                             _mm256_mul_pd(lhs.reg.y, rhs.reg.y),
                             _mm256_mul_pd(lhs.reg.z, rhs.reg.z),
                             _mm256_mul_pd(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(lhs.data[0] * rhs.data[0],
                             lhs.data[1] * rhs.data[1],
                             lhs.data[2] * rhs.data[2],
                             lhs.data[3] * rhs.data[3],
                             lhs.data[4] * rhs.data[4],
                             lhs.data[5] * rhs.data[5],
                             lhs.data[6] * rhs.data[6],
                             lhs.data[7] * rhs.data[7],
                             lhs.data[8] * rhs.data[8],
                             lhs.data[9] * rhs.data[9],
                             lhs.data[10] * rhs.data[10],
                             lhs.data[11] * rhs.data[11],
                             lhs.data[12] * rhs.data[12],
                             lhs.data[13] * rhs.data[13],
                             lhs.data[14] * rhs.data[14],
                             lhs.data[15] * rhs.data[15]);
#endif
  }

  inline vtype<double, 16> operator/(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_div_pd(lhs.reg.x, rhs.reg.x),
                             _mm512_div_pd(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_div_pd(lhs.reg.x, rhs.reg.x),
                             _mm256_div_pd(lhs.reg.y, rhs.reg.y),
                             _mm256_div_pd(lhs.reg.z, rhs.reg.z),
                             _mm256_div_pd(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(lhs.data[0] / rhs.data[0],
                             lhs.data[1] / rhs.data[1],
                             lhs.data[2] / rhs.data[2],
                             lhs.data[3] / rhs.data[3],
                             lhs.data[4] / rhs.data[4],
                             lhs.data[5] / rhs.data[5],
                             lhs.data[6] / rhs.data[6],
                             lhs.data[7] / rhs.data[7],
                             lhs.data[8] / rhs.data[8],
                             lhs.data[9] / rhs.data[9],
                             lhs.data[10] / rhs.data[10],
                             lhs.data[11] / rhs.data[11],
                             lhs.data[12] / rhs.data[12],
                             lhs.data[13] / rhs.data[13],
                             lhs.data[14] / rhs.data[14],
                             lhs.data[15] / rhs.data[15]);
#endif
  }

  inline vtype<double, 16> operator-(vtype<double, 16> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    __m512d sign_mask = _mm512_castsi512_pd(_mm512_set1_epi64(0x8000000000000000));
    return vtype<double, 16>(_mm512_xor_pd(v.reg.x, sign_mask),
                             _mm512_xor_pd(v.reg.y, sign_mask));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    __m512i sign_mask = _mm512_set1_epi64(0x8000000000000000);
    return vtype<double, 16>(_mm512_castsi512_pd(_mm512_xor_epi64(_mm512_castpd_si512(v.reg.x), sign_mask)),
                             _mm512_castsi512_pd(_mm512_xor_epi64(_mm512_castpd_si512(v.reg.y), sign_mask)));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256d sign_mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000));
    return vtype<double, 16>(_mm256_xor_pd(v.reg.x, sign_mask),
                             _mm256_xor_pd(v.reg.y, sign_mask),
                             _mm256_xor_pd(v.reg.z, sign_mask),
                             _mm256_xor_pd(v.reg.w, sign_mask));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(-v.data[0],
                             -v.data[1],
                             -v.data[2],
                             -v.data[3],
                             -v.data[4],
                             -v.data[5],
                             -v.data[6],
                             -v.data[7],
                             -v.data[8],
                             -v.data[9],
                             -v.data[10],
                             -v.data[11],
                             -v.data[12],
                             -v.data[13],
                             -v.data[14],
                             -v.data[15]);
#endif
  }

  inline void mem_store(double * mem, vtype<double, 16> const & op) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    _mm512_storeu_pd(mem, op.reg.x);
    _mm512_storeu_pd(mem+8, op.reg.y);
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    _mm256_storeu_pd(mem, op.reg.x);
    _mm256_storeu_pd(mem+4, op.reg.y);
    _mm256_storeu_pd(mem+8, op.reg.z);
    _mm256_storeu_pd(mem+12, op.reg.w);
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNIGNS)
#warning "scalar implementation"
# endif
    mem[0] = op.data[0];
    mem[1] = op.data[1];
    mem[2] = op.data[2];
    mem[3] = op.data[3];
    mem[4] = op.data[4];
    mem[5] = op.data[5];
    mem[6] = op.data[6];
    mem[7] = op.data[7];
    mem[8] = op.data[8];
    mem[9] = op.data[9];
    mem[10] = op.data[10];
    mem[11] = op.data[11];
    mem[12] = op.data[12];
    mem[13] = op.data[13];
    mem[14] = op.data[14];
    mem[15] = op.data[15];
#endif
  }

  inline vtype<double, 16> fabs(vtype<double, 16> const & v) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512DQ__)
    __m512d sign_mask = _mm512_castsi512_pd(_mm512_set1_epi64(0x7FFFFFFFFFFFFFFF));
    return vtype<double, 16>(_mm512_and_pd(v.reg.x, sign_mask),
                             _mm512_and_pd(v.reg.y, sign_mask));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    __m512i sign_mask = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFF);
    return vtype<double, 16>(_mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(v.reg.x), sign_mask)),
                             _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(v.reg.y), sign_mask)));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256d sign_mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
    return vtype<double, 16>(_mm256_and_pd(v.reg.x, sign_mask),
                             _mm256_and_pd(v.reg.y, sign_mask),
                             _mm256_and_pd(v.reg.z, sign_mask),
                             _mm256_and_pd(v.reg.w, sign_mask));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(std::fabs(v.data[0]),
                             std::fabs(v.data[1]),
                             std::fabs(v.data[2]),
                             std::fabs(v.data[3]),
                             std::fabs(v.data[4]),
                             std::fabs(v.data[5]),
                             std::fabs(v.data[6]),
                             std::fabs(v.data[7]),
                             std::fabs(v.data[8]),
                             std::fabs(v.data[9]),
                             std::fabs(v.data[10]),
                             std::fabs(v.data[11]),
                             std::fabs(v.data[12]),
                             std::fabs(v.data[13]),
                             std::fabs(v.data[14]),
                             std::fabs(v.data[15]));
#endif
  }

  inline vtype<double, 16> abs(vtype<double, 16> const & v) {
    return fabs(v);
  }

  inline vtype<double, 16> max(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_max_pd(lhs.reg.x, rhs.reg.x),
                             _mm512_max_pd(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_max_pd(lhs.reg.x, rhs.reg.x),
                             _mm256_max_pd(lhs.reg.y, rhs.reg.y),
                             _mm256_max_pd(lhs.reg.z, rhs.reg.z),
                             _mm256_max_pd(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(std::max(lhs.data[0], rhs.data[0]),
                             std::max(lhs.data[1], rhs.data[1]),
                             std::max(lhs.data[2], rhs.data[2]),
                             std::max(lhs.data[3], rhs.data[3]),
                             std::max(lhs.data[4], rhs.data[4]),
                             std::max(lhs.data[5], rhs.data[5]),
                             std::max(lhs.data[6], rhs.data[6]),
                             std::max(lhs.data[7], rhs.data[7]),
                             std::max(lhs.data[8], rhs.data[8]),
                             std::max(lhs.data[9], rhs.data[9]),
                             std::max(lhs.data[10], rhs.data[10]),
                             std::max(lhs.data[11], rhs.data[11]),
                             std::max(lhs.data[12], rhs.data[12]),
                             std::max(lhs.data[13], rhs.data[13]),
                             std::max(lhs.data[14], rhs.data[14]),
                             std::max(lhs.data[15], rhs.data[15]));
#endif
  }

  inline vtype<double, 16> min(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_min_pd(lhs.reg.x, rhs.reg.x),
                             _mm512_min_pd(lhs.reg.y, rhs.reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_min_pd(lhs.reg.x, rhs.reg.x),
                             _mm256_min_pd(lhs.reg.y, rhs.reg.y),
                             _mm256_min_pd(lhs.reg.z, rhs.reg.z),
                             _mm256_min_pd(lhs.reg.w, rhs.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(std::min(lhs.data[0], rhs.data[0]),
                             std::min(lhs.data[1], rhs.data[1]),
                             std::min(lhs.data[2], rhs.data[2]),
                             std::min(lhs.data[3], rhs.data[3]),
                             std::min(lhs.data[4], rhs.data[4]),
                             std::min(lhs.data[5], rhs.data[5]),
                             std::min(lhs.data[6], rhs.data[6]),
                             std::min(lhs.data[7], rhs.data[7]),
                             std::min(lhs.data[8], rhs.data[8]),
                             std::min(lhs.data[9], rhs.data[9]),
                             std::min(lhs.data[10], rhs.data[10]),
                             std::min(lhs.data[11], rhs.data[11]),
                             std::min(lhs.data[12], rhs.data[12]),
                             std::min(lhs.data[13], rhs.data[13]),
                             std::min(lhs.data[14], rhs.data[14]),
                             std::min(lhs.data[15], rhs.data[15]));
#endif
  }

  inline vmask<sizeof(double), 16> operator<(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 16>(_mm512_cmp_pd_mask(lhs.reg.x, rhs.reg.x, _CMP_LT_OQ),
                                     _mm512_cmp_pd_mask(lhs.reg.y, rhs.reg.y, _CMP_LT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 16>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_LT_OQ),
                                     _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_LT_OQ),
                                     _mm256_cmp_pd(lhs.reg.z, rhs.reg.z, _CMP_LT_OQ),
                                     _mm256_cmp_pd(lhs.reg.w, rhs.reg.w, _CMP_LT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 16>(lhs.data[0] < rhs.data[0],
                                     lhs.data[1] < rhs.data[1],
                                     lhs.data[2] < rhs.data[2],
                                     lhs.data[3] < rhs.data[3],
                                     lhs.data[4] < rhs.data[4],
                                     lhs.data[5] < rhs.data[5],
                                     lhs.data[6] < rhs.data[6],
                                     lhs.data[7] < rhs.data[7],
                                     lhs.data[8] < rhs.data[8],
                                     lhs.data[9] < rhs.data[9],
                                     lhs.data[10] < rhs.data[10],
                                     lhs.data[11] < rhs.data[11],
                                     lhs.data[12] < rhs.data[12],
                                     lhs.data[13] < rhs.data[13],
                                     lhs.data[14] < rhs.data[14],
                                     lhs.data[15] < rhs.data[15]);
#endif
  }

  inline vmask<sizeof(double), 16> operator<=(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 16>(_mm512_cmp_pd_mask(lhs.reg.x, rhs.reg.x, _CMP_LE_OQ),
                                     _mm512_cmp_pd_mask(lhs.reg.y, rhs.reg.y, _CMP_LE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 16>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_LE_OQ),
                                     _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_LE_OQ),
                                     _mm256_cmp_pd(lhs.reg.z, rhs.reg.z, _CMP_LE_OQ),
                                     _mm256_cmp_pd(lhs.reg.w, rhs.reg.w, _CMP_LE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 16>(lhs.data[0] <= rhs.data[0],
                                     lhs.data[1] <= rhs.data[1],
                                     lhs.data[2] <= rhs.data[2],
                                     lhs.data[3] <= rhs.data[3],
                                     lhs.data[4] <= rhs.data[4],
                                     lhs.data[5] <= rhs.data[5],
                                     lhs.data[6] <= rhs.data[6],
                                     lhs.data[7] <= rhs.data[7],
                                     lhs.data[8] <= rhs.data[8],
                                     lhs.data[9] <= rhs.data[9],
                                     lhs.data[10] <= rhs.data[10],
                                     lhs.data[11] <= rhs.data[11],
                                     lhs.data[12] <= rhs.data[12],
                                     lhs.data[13] <= rhs.data[13],
                                     lhs.data[14] <= rhs.data[14],
                                     lhs.data[15] <= rhs.data[15]);
#endif
  }

  inline vmask<sizeof(double), 16> operator>(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 16>(_mm512_cmp_pd_mask(lhs.reg.x, rhs.reg.x, _CMP_GT_OQ),
                                     _mm512_cmp_pd_mask(lhs.reg.y, rhs.reg.y, _CMP_GT_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 16>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_GT_OQ),
                                     _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_GT_OQ),
                                     _mm256_cmp_pd(lhs.reg.z, rhs.reg.z, _CMP_GT_OQ),
                                     _mm256_cmp_pd(lhs.reg.w, rhs.reg.w, _CMP_GT_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 16>(lhs.data[0] > rhs.data[0],
                                     lhs.data[1] > rhs.data[1],
                                     lhs.data[2] > rhs.data[2],
                                     lhs.data[3] > rhs.data[3],
                                     lhs.data[4] > rhs.data[4],
                                     lhs.data[5] > rhs.data[5],
                                     lhs.data[6] > rhs.data[6],
                                     lhs.data[7] > rhs.data[7],
                                     lhs.data[8] > rhs.data[8],
                                     lhs.data[9] > rhs.data[9],
                                     lhs.data[10] > rhs.data[10],
                                     lhs.data[11] > rhs.data[11],
                                     lhs.data[12] > rhs.data[12],
                                     lhs.data[13] > rhs.data[13],
                                     lhs.data[14] > rhs.data[14],
                                     lhs.data[15] > rhs.data[15]);
#endif
  }

  inline vmask<sizeof(double), 16> operator>=(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 16>(_mm512_cmp_pd_mask(lhs.reg.x, rhs.reg.x, _CMP_GE_OQ),
                                     _mm512_cmp_pd_mask(lhs.reg.y, rhs.reg.y, _CMP_GE_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 16>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_GE_OQ),
                                     _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_GE_OQ),
                                     _mm256_cmp_pd(lhs.reg.z, rhs.reg.z, _CMP_GE_OQ),
                                     _mm256_cmp_pd(lhs.reg.w, rhs.reg.w, _CMP_GE_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 16>(lhs.data[0] >= rhs.data[0],
                                     lhs.data[1] >= rhs.data[1],
                                     lhs.data[2] >= rhs.data[2],
                                     lhs.data[3] >= rhs.data[3],
                                     lhs.data[4] >= rhs.data[4],
                                     lhs.data[5] >= rhs.data[5],
                                     lhs.data[6] >= rhs.data[6],
                                     lhs.data[7] >= rhs.data[7],
                                     lhs.data[8] >= rhs.data[8],
                                     lhs.data[9] >= rhs.data[9],
                                     lhs.data[10] >= rhs.data[10],
                                     lhs.data[11] >= rhs.data[11],
                                     lhs.data[12] >= rhs.data[12],
                                     lhs.data[13] >= rhs.data[13],
                                     lhs.data[14] >= rhs.data[14],
                                     lhs.data[15] >= rhs.data[15]);
#endif
  }

  inline vmask<sizeof(double), 16> operator==(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 16>(_mm512_cmp_pd_mask(lhs.reg.x, rhs.reg.x, _CMP_EQ_OQ),
                                     _mm512_cmp_pd_mask(lhs.reg.y, rhs.reg.y, _CMP_EQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 16>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_EQ_OQ),
                                     _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_EQ_OQ),
                                     _mm256_cmp_pd(lhs.reg.z, rhs.reg.z, _CMP_EQ_OQ),
                                     _mm256_cmp_pd(lhs.reg.w, rhs.reg.w, _CMP_EQ_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 16>(lhs.data[0] == rhs.data[0],
                                     lhs.data[1] == rhs.data[1],
                                     lhs.data[2] == rhs.data[2],
                                     lhs.data[3] == rhs.data[3],
                                     lhs.data[4] == rhs.data[4],
                                     lhs.data[5] == rhs.data[5],
                                     lhs.data[6] == rhs.data[6],
                                     lhs.data[7] == rhs.data[7],
                                     lhs.data[8] == rhs.data[8],
                                     lhs.data[9] == rhs.data[9],
                                     lhs.data[10] == rhs.data[10],
                                     lhs.data[11] == rhs.data[11],
                                     lhs.data[12] == rhs.data[12],
                                     lhs.data[13] == rhs.data[13],
                                     lhs.data[14] == rhs.data[14],
                                     lhs.data[15] == rhs.data[15]);
#endif
  }

  inline vmask<sizeof(double), 16> operator!=(vtype<double, 16> const & lhs, vtype<double, 16> const & rhs) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vmask<sizeof(double), 16>(_mm512_cmp_pd_mask(lhs.reg.x, rhs.reg.x, _CMP_NEQ_OQ),
                                     _mm512_cmp_pd_mask(lhs.reg.y, rhs.reg.y, _CMP_NEQ_OQ));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vmask<sizeof(double), 16>(_mm256_cmp_pd(lhs.reg.x, rhs.reg.x, _CMP_NEQ_OQ),
                                     _mm256_cmp_pd(lhs.reg.y, rhs.reg.y, _CMP_NEQ_OQ),
                                     _mm256_cmp_pd(lhs.reg.z, rhs.reg.z, _CMP_NEQ_OQ),
                                     _mm256_cmp_pd(lhs.reg.w, rhs.reg.w, _CMP_NEQ_OQ));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vmask<sizeof(double), 16>(lhs.data[0] != rhs.data[0],
                                     lhs.data[1] != rhs.data[1],
                                     lhs.data[2] != rhs.data[2],
                                     lhs.data[3] != rhs.data[3],
                                     lhs.data[4] != rhs.data[4],
                                     lhs.data[5] != rhs.data[5],
                                     lhs.data[6] != rhs.data[6],
                                     lhs.data[7] != rhs.data[7],
                                     lhs.data[8] != rhs.data[8],
                                     lhs.data[9] != rhs.data[9],
                                     lhs.data[10] != rhs.data[10],
                                     lhs.data[11] != rhs.data[11],
                                     lhs.data[12] != rhs.data[12],
                                     lhs.data[13] != rhs.data[13],
                                     lhs.data[14] != rhs.data[14],
                                     lhs.data[15] != rhs.data[15]);
#endif
  }

  inline vtype<double, 16> select(vmask<sizeof(double), 16> const & c,
                                  vtype<double, 16> const & t,
                                  vtype<double, 16> const & f) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_mask_blend_pd(c.reg.x, f.reg.x, t.reg.x),
                             _mm512_mask_blend_pd(c.reg.y, f.reg.y, t.reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_blendv_pd(f.reg.x, t.reg.x, _mm256_castsi256_pd(c.reg.x)),
                             _mm256_blendv_pd(f.reg.y, t.reg.y, _mm256_castsi256_pd(c.reg.y)),
                             _mm256_blendv_pd(f.reg.z, t.reg.z, _mm256_castsi256_pd(c.reg.z)),
                             _mm256_blendv_pd(f.reg.w, t.reg.w, _mm256_castsi256_pd(c.reg.w)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>((c.data & 1u) ? t.data[0] : f.data[0],
                             (c.data & 2u) ? t.data[1] : f.data[1],
                             (c.data & 4u) ? t.data[2] : f.data[2],
                             (c.data & 8u) ? t.data[3] : f.data[3],
                             (c.data & 16u) ? t.data[4] : f.data[4],
                             (c.data & 32u) ? t.data[5] : f.data[5],
                             (c.data & 64u) ? t.data[6] : f.data[6],
                             (c.data & 128u) ? t.data[7] : f.data[7],
                             (c.data & 256u) ? t.data[8] : f.data[8],
                             (c.data & 512u) ? t.data[9] : f.data[9],
                             (c.data & 1024u) ? t.data[10] : f.data[10],
                             (c.data & 2048u) ? t.data[11] : f.data[11],
                             (c.data & 4096u) ? t.data[12] : f.data[12],
                             (c.data & 8192u) ? t.data[13] : f.data[13],
                             (c.data & 16384u) ? t.data[14] : f.data[14],
                             (c.data & 32768u) ? t.data[15] : f.data[15]);
#endif
  }

  inline vtype<double, 16> ceil(vtype<double, 16> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_roundscale_pd(arg.reg.x, _MM_FROUND_TO_POS_INF|_MM_FROUND_NO_EXC),
                             _mm512_roundscale_pd(arg.reg.y, _MM_FROUND_TO_POS_INF|_MM_FROUND_NO_EXC));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_ceil_pd(arg.reg.x),
                             _mm256_ceil_pd(arg.reg.y),
                             _mm256_ceil_pd(arg.reg.z),
                             _mm256_ceil_pd(arg.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(std::ceil(arg.data[0]),
                             std::ceil(arg.data[1]),
                             std::ceil(arg.data[2]),
                             std::ceil(arg.data[3]),
                             std::ceil(arg.data[4]),
                             std::ceil(arg.data[5]),
                             std::ceil(arg.data[6]),
                             std::ceil(arg.data[7]),
                             std::ceil(arg.data[8]),
                             std::ceil(arg.data[9]),
                             std::ceil(arg.data[10]),
                             std::ceil(arg.data[11]),
                             std::ceil(arg.data[12]),
                             std::ceil(arg.data[13]),
                             std::ceil(arg.data[14]),
                             std::ceil(arg.data[15]));
#endif
  }

  inline vtype<double, 16> floor(vtype<double, 16> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_roundscale_pd(arg.reg.x, _MM_FROUND_TO_NEG_INF|_MM_FROUND_NO_EXC),
                             _mm512_roundscale_pd(arg.reg.y, _MM_FROUND_TO_NEG_INF|_MM_FROUND_NO_EXC));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_floor_pd(arg.reg.x),
                             _mm256_floor_pd(arg.reg.y),
                             _mm256_floor_pd(arg.reg.z),
                             _mm256_floor_pd(arg.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(std::floor(arg.data[0]),
                             std::floor(arg.data[1]),
                             std::floor(arg.data[2]),
                             std::floor(arg.data[3]),
                             std::floor(arg.data[4]),
                             std::floor(arg.data[5]),
                             std::floor(arg.data[6]),
                             std::floor(arg.data[7]),
                             std::floor(arg.data[8]),
                             std::floor(arg.data[9]),
                             std::floor(arg.data[10]),
                             std::floor(arg.data[11]),
                             std::floor(arg.data[12]),
                             std::floor(arg.data[13]),
                             std::floor(arg.data[14]),
                             std::floor(arg.data[15]));
#endif
  }

  inline vtype<double, 16> round(vtype<double, 16> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_roundscale_pd(arg.reg.x, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC),
                             _mm512_roundscale_pd(arg.reg.y, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_round_pd(arg.reg.x, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC),
                             _mm256_round_pd(arg.reg.y, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC),
                             _mm256_round_pd(arg.reg.z, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC),
                             _mm256_round_pd(arg.reg.w, _MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(std::round(arg.data[0]),
                             std::round(arg.data[1]),
                             std::round(arg.data[2]),
                             std::round(arg.data[3]),
                             std::round(arg.data[4]),
                             std::round(arg.data[5]),
                             std::round(arg.data[6]),
                             std::round(arg.data[7]),
                             std::round(arg.data[8]),
                             std::round(arg.data[9]),
                             std::round(arg.data[10]),
                             std::round(arg.data[11]),
                             std::round(arg.data[12]),
                             std::round(arg.data[13]),
                             std::round(arg.data[14]),
                             std::round(arg.data[15]));
#endif
  }

  inline vtype<double, 16> sqrt(vtype<double, 16> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_sqrt_pd(arg.reg.x),
                             _mm512_sqrt_pd(arg.reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    return vtype<double, 16>(_mm256_sqrt_pd(arg.reg.x),
                             _mm256_sqrt_pd(arg.reg.y),
                             _mm256_sqrt_pd(arg.reg.z),
                             _mm256_sqrt_pd(arg.reg.w));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(std::sqrt(arg.data[0]),
                             std::sqrt(arg.data[1]),
                             std::sqrt(arg.data[2]),
                             std::sqrt(arg.data[3]),
                             std::sqrt(arg.data[4]),
                             std::sqrt(arg.data[5]),
                             std::sqrt(arg.data[6]),
                             std::sqrt(arg.data[7]),
                             std::sqrt(arg.data[8]),
                             std::sqrt(arg.data[9]),
                             std::sqrt(arg.data[10]),
                             std::sqrt(arg.data[11]),
                             std::sqrt(arg.data[12]),
                             std::sqrt(arg.data[13]),
                             std::sqrt(arg.data[14]),
                             std::sqrt(arg.data[15]));
#endif
  }

  inline vtype<double, 16> rsqrt(vtype<double, 16> const & arg) {
#if !defined(VTYPE_SCALAR_ONLY) && defined(__AVX512F__)
    return vtype<double, 16>(_mm512_rsqrt14_pd(arg.reg.x),
                             _mm512_rsqrt14_pd(arg.reg.y));
#elif !defined(VTYPE_SCALAR_ONLY) && defined(__AVX2__)
    __m256d one = _mm256_set1_pd(1.0);
    return vtype<double, 16>(_mm256_div_pd(one, _mm256_sqrt_pd(arg.reg.x)),
                             _mm256_div_pd(one, _mm256_sqrt_pd(arg.reg.y)),
                             _mm256_div_pd(one, _mm256_sqrt_pd(arg.reg.z)),
                             _mm256_div_pd(one, _mm256_sqrt_pd(arg.reg.w)));
#else
# if defined(VTYPE_ENABLE_SCALAR_WARNINGS)
#warning "scalar implementation"
# endif
    return vtype<double, 16>(1.0/std::sqrt(arg.data[0]),
                             1.0/std::sqrt(arg.data[1]),
                             1.0/std::sqrt(arg.data[2]),
                             1.0/std::sqrt(arg.data[3]),
                             1.0/std::sqrt(arg.data[4]),
                             1.0/std::sqrt(arg.data[5]),
                             1.0/std::sqrt(arg.data[6]),
                             1.0/std::sqrt(arg.data[7]),
                             1.0/std::sqrt(arg.data[8]),
                             1.0/std::sqrt(arg.data[9]),
                             1.0/std::sqrt(arg.data[10]),
                             1.0/std::sqrt(arg.data[11]),
                             1.0/std::sqrt(arg.data[12]),
                             1.0/std::sqrt(arg.data[13]),
                             1.0/std::sqrt(arg.data[14]),
                             1.0/std::sqrt(arg.data[15]));
#endif
  }

} // end: namespace Loci

#endif
