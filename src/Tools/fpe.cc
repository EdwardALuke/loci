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
#include <Tools/fpe.h>
#include <Tools/debugger.h>
#include <stdlib.h>

#ifdef DARWIN
#ifdef SYSTEM_ARM64
#include <fenv.h>
#include <iostream>
#include <signal.h>
#include <math.h>
#else
#include <xmmintrin.h>
#endif
#endif

#ifdef SPARC

extern "C" {
  extern void SPARC_set_fpe_abort() ;
  void fpe_debugger_() { Loci::debugger_() ; }
}

namespace Loci {
  void set_fpe_abort() {
    SPARC_set_fpe_abort() ;
  }
}
#else

#ifdef SGI

#include <iostream>
#include <stdio.h>
#include <sigfpe.h>
extern "C" {
  void sgi_fpe_user_(unsigned int v[5], int y[2]) {
    std::cerr << "Floating Point Exception" << std::endl ;
    Loci::debugger_() ;
  }
  void sgi_fpe_abort_(unsigned int **p) {
    std::cerr << "Floating Point Exception" << std::endl ;
    Loci::debugger_() ;
  }
}

namespace Loci {
  void set_fpe_abort() {
    //    sigfpe_[_DIVZERO].abort = 1 ;
    //    sigfpe_[_OVERFL].abort = 1 ;
    //    sigfpe_[_INVALID].abort = 1 ;
    handle_sigfpes(_ON,_EN_OVERFL|_EN_DIVZERO|_EN_INVALID,
                   sgi_fpe_user_,_ABORT_ON_ERROR,sgi_fpe_abort_) ;
                   //                   0,_ABORT_ON_ERROR,0) ;
  }
}
                      
#else
#ifdef LINUX

#include <stdio.h>
#include <signal.h>
#include <iostream>

#ifndef NO_FENV
// Floating point exception environment
#ifndef _GNU_SOURCE
#define _GNU_SOURCE     /* needed to get non POSIX extensions. */
#endif
#include <fenv.h>
#endif

using std::cerr ;
using std::endl ;

extern "C" {
  void fpe_debugger_(int i)
  {
    std::cerr << "floating point exception" << std::endl ;
    Loci::debugger_() ;
  }
}

namespace Loci {
  
  void set_fpe_abort()
  {
#ifndef NO_FENV

#ifdef SYSTEM_ITANIUM64
    
    if(feenableexcept((FE_DIVBYZERO)) == -1) {
      std::cerr << "feenableexcept had error, floating point exceptions not caught" << std::endl ;
    } else {
      signal(SIGFPE,fpe_debugger_) ;
    }
#else
    if(feenableexcept((FE_DIVBYZERO|FE_OVERFLOW|FE_INVALID)) == -1) {
      std::cerr << "feenableexcept had error, floating point exceptions not caught" << std::endl ;
    } else {
      signal(SIGFPE,fpe_debugger_) ;
    }
#endif

#endif
  }
}

#else

namespace Loci {
#ifdef DARWIN
#ifdef SYSTEM_ARM64
  static void
  fpe_signal_handler( int sig, siginfo_t *sip, void *scp )
  {
    int fe_code = sip->si_code;

    if (fe_code == ILL_ILLTRP)
      std::cerr << "Illegal trap detected" << std::endl ;
    else
      std::cerr << "Code detected : " << fe_code << std::endl ;

    debugger_();
  }
#endif
  void set_fpe_abort()
  {
#ifdef SYSTEM_ARM64
        fenv_t env;
    fegetenv(&env);

    env.__fpcr = env.__fpcr | __fpcr_trap_invalid;
    fesetenv(&env);

    struct sigaction act;
    act.sa_sigaction = fpe_signal_handler ;
    sigemptyset (&act.sa_mask);
    act.sa_flags = SA_SIGINFO;
    sigaction(SIGILL, &act, NULL);
#else
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID) ;
#endif
  }
#else
  void set_fpe_abort()
  {
    
  }
#endif
}
#endif
#endif
#endif

