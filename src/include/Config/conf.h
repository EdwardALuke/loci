//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#ifndef LOCI_CONFIGURATION_INCLUDE
#define LOCI_CONFIGURATION_INCLUDE

#include "Config/version_info.h"

#define DYNAMICSCHEDULING

#ifdef restrict
#undef restrict
#endif

#if defined(__GNUC__)

#if defined(__clang__) 
// clang compiler
#define USE_MAP_FOR_HASH_MAP
#define restrict __restrict__
#define NO_OFFSETOF

#else
#if defined(__INTEL_COMPILER)
// Intel Compiler

//#define NO_FENV

#define USE_MAP_FOR_HASH_MAP
//#define EXT_HASH_MAP
//#define EXT_NAMESPACE std
#define HAVE_IVDEP
#define NO_SIGNBIT

#define restrict __restrict
#else
#define restrict __restrict__
// GCC Compilers
#if (__GNUC__==3) && (__GNUC_MINOR__==0)
#define EXT_HASH_MAP
#endif

#if (__GNUC__==3) && (__GNUC_MINOR__!=0)
#define EXT_HASH_MAP
#define EXT_NAMESPACE __gnu_cxx
#define NO_OFFSETOF
#endif

#if (__GNUC__==4)
#if (__GNUC_MINOR__ < 3)
#define EXT_HASH_MAP
#define EXT_NAMESPACE __gnu_cxx
#define NO_OFFSETOF
#else
#define USE_MAP_FOR_HASH_MAP
#define NO_OFFSETOF
#endif
#endif
#if (__GNUC__ > 4)
#define USE_MAP_FOR_HASH_MAP
#define NO_OFFSETOF
#endif

#endif

#endif
#else
#if defined(__IBMCPP__)
/* IBM C++ Compiler */
#define USE_MAP_FOR_HASH_MAP
#define NO_OFFSETOF
#define restrict __restrict
#else
#if defined(__sgi)
// SGI Compiler
#define NO_CSTDLIB
#define NO_CMATH
#define NO_SIGNBIT
//#define MPI_NO_CPPBIND
#define restrict __restrict
#else
#if defined(__PGI)
/* Portland group compiler*/
#define restrict __restrict
#define HAVE_IVDEP
#define NO_FENV
#define NO_SIGNBIT
#define USE_MAP_FOR_HASH_MAP
#define __thread
#else
#if defined(__SUNPRO_CC)
/* Sun CC compiler */
#define USE_MAP_FOR_HASH_MAP
#define SUN_CC
#define restrict 
#else

#define restrict

#endif

#endif

#endif

#endif

#endif

#if defined(LOCI_SYS_IRIX64)
#undef SGI

#define SGI
#endif

#if defined(LOCI_SYS_Linux)
#undef LINUX
#define LINUX

#endif

#if defined(LOCI_ARCH_ia64)
#define SYSTEM_ITANIUM64
#endif

#if defined(LOCI_ARCH_arm64)
#define SYSTEM_ARM64
#endif

#if defined(LOCI_SYS_SunOS)
#undef SPARC
#define SUN
#define NO_SIGNBIT
#define SPARC
#endif

#if defined(LOCI_SYS_Darwin)
#ifndef BSD
#define BSD
#endif
#define DARWIN
#define NO_THREAD_MEMORY
#endif

#ifdef __CYGWIN__
#define NO_XDR_CPP_PROTOTYPES
#endif

#ifndef MPI_NO_CPPBIND
#define MPI_NO_CPPBIND
#endif

#endif
