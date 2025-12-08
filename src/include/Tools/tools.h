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
#ifndef TOOLS_H
#define TOOLS_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

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
//using std::isnan ;

namespace Loci {
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
  using std::isnan ;
}

#ifdef NO_CSTDLIB
#include <stdlib.h>
#ifdef NO_ABS
inline float abs(float f1) { return f1>0?f1:-f1; }
inline double abs(double f1) { return f1>0?f1:-f1; }
inline long double abs(long double f1) { return f1>0?f1:-f1; }
#endif
#else
#include <cstdlib>
using std::abs ;
namespace Loci {
  using std::abs ;
}
#endif

#include <algorithm>
using std::max ;
using std::min ;

#include <utility>
#include <iosfwd>

// Hack to define signbit if one is not provided by the compiler
#ifdef NO_SIGNBIT
#undef signbit
#define signbit(X) signbit_replace(X)

union double_long_helper {
  double d ;
  long long ll;
};

union float_long_helper {
  float f ;
  long l ;
} ;

inline int signbit_replace(double x) {
  double_long_helper ul ;
  ul.d = x ;
  return (ul.ll < 0)?1:0 ;
}

inline int signbit_replace(float x) {
  float_long_helper ul ;
  ul.f = x ;
  return ((ul.l & 0x80000000) != 0)?1:0 ;
}
#else
using std::signbit ;
#endif

template <class T> inline T square(const T &a)   { return a*a ; }

template <class T> inline T sign(T a, T b) { return (b>=0)?abs(a):-abs(a) ; }

namespace Loci {
  void register_closing_function(void (*fptr)(int code)) ;
}

#endif
