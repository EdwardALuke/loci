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
#include <Config/conf.h>
#ifdef SPARC
#ifdef NOSUNMATH
void SPARC_set_fpe_abort() {
}
#else

#include <stdio.h>
#include <stdlib.h>
#include <floatingpoint.h>
#include <siginfo.h>
#include <ucontext.h>
#include <sunmath.h>

extern void fpe_debugger_() ;

void SPARC_ieee_abort(int sig, siginfo_t *sip, ucontext_t *uap) {
  const char *label;
  
  switch (sip->si_code) {
  case FPE_FLTINV: label = "invalid operand"; break;
  case FPE_FLTRES: label = "inexact"; break;
  case FPE_FLTDIV: label = "division-by-zero"; break;
  case FPE_FLTUND: label = "underflow"; break;
  case FPE_FLTOVF: label = "overflow"; break;
  default: label = "???"; break;
  }
    
  fprintf(stderr, "FP exception %s (0x%x) occurred at address %p.\n",
          label, sip->si_code, (void *) sip->si_addr);
  fpe_debugger_() ;
  abort();
}

void SPARC_set_fpe_abort() {
  sigfpe_handler_type hdl ;
  hdl = (sigfpe_handler_type) SPARC_ieee_abort ;
  ieee_handler("set","common",hdl) ;
}
#endif
#else
void SPARC_set_fpe_abort() {
}
#endif

