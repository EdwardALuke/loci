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

#include <stdio.h>
#include <signal.h>
#ifdef SUN
#include <floatingpoint.h>
#include <siginfo.h>
#endif

#ifndef __CYGWIN__
#ifndef DARWIN
#include <ucontext.h>
#endif
#endif

#ifdef SUN
#include <sys/machsig.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <stdlib.h>
#include <unistd.h>

#include "Tools/debugger.h"

#include <iostream>

using std::cerr ;
using std::endl ;

#define HOST_ID "localhost"

namespace Loci {

  long getmaxrss() {
    rusage usage ;
    getrusage(RUSAGE_SELF,&usage) ;
    return usage.ru_maxrss ;
  }
  
  const char *debug_hostname = HOST_ID ;
  const char *debug_execname = "a.out" ;
#ifdef SGI
  const char *debug_program = "dbx" ;
#else
  const char *debug_program = "gdb" ;
#endif
  
  bool debugger_setup = false ;

  void (*debug_callback)() = 0 ;

  void set_debug_callback(void (*dcb)()) {
    debug_callback = dcb ;
  }

  
  void setup_debugger(const char *execname, const char *gdb, const char *hostname) {
    debug_execname = execname ;
    debug_program = gdb ;
    debug_hostname = hostname ;
    debugger_setup = true ;
  }

  struct closing_function_list {
    struct closing_function_list *next ;
    void (*fptr)(int code) ;
  } ;

  closing_function_list *closing_funcs_list = 0 ;

  void register_closing_function(void (*fptr)(int code)) {
    closing_function_list *fl = new closing_function_list ;
    fl->next = closing_funcs_list ;
    fl->fptr = fptr ;
    closing_funcs_list = fl ;
  }

  void call_closing_functions(int code) {
    // loop over registed closing functions
    closing_function_list *fl = closing_funcs_list ;
    while(fl != 0) {
      (*(fl->fptr))(code) ;
      fl = fl->next ;
    }
  }

  void debugger_()
  {
    call_closing_functions(-1) ;

    void (*dbc)() = debug_callback ;
    debug_callback = 0 ;

    if(dbc == 0) {
      call_closing_functions(-2) ;
      exit(-1) ;
    }

    if(dbc != 0)
      (*dbc)() ;

    if(!debugger_setup) {
      call_closing_functions(-2) ;
      exit(-1) ;
    }
    
#ifdef BSD
    exit(-1) ;
#else
    char buf[512] ;
    int pid = getpid() ;
    const char *xtermpath = "xterm" ;
#ifndef SGI
    snprintf(buf,512,"%s  -display %s -e %s %s %d &",
            xtermpath,
            debug_hostname,debug_program,debug_execname,pid) ;
#else
    snprintf(buf,512,"%s  -display %s -e %s -p %d %s &",
            xtermpath,
            debug_hostname,debug_program,pid,debug_execname) ;
#endif
    cerr << buf << endl ;
    int err = system(buf) ;
    if(err != 0)
      cerr << "system call failed with on '"<< buf << "'" << endl ;


    sleep(100) ; /* Wait for debugger to attach */
#endif
  }

}

extern "C" {
  void program_trap(int sig) /*,code,scp,addr)*/
  {
    const char *sigtype = "(undefined)" ;
    
    switch(sig) {
    case SIGBUS:
      sigtype = "a Bus Error" ;
      break ;
    case SIGSEGV:
      sigtype = "a Segmentation Violation" ;
      break ;
    case SIGILL:
      sigtype = "an Illegal Instruction Call" ;
      break ;
    case SIGSYS:
      sigtype = "an Illegal System Call" ;
      break ;
    case SIGFPE:
      sigtype = "a Floating Point Exception" ;
      break ;
    }
    fprintf(stderr,"ERROR: Program terminated due to %s\n",sigtype) ;
    Loci::debugger_() ;
    Loci::call_closing_functions(-2) ;
    exit(-1) ;
  }
}


namespace Loci {
  

  void chopsigs_()
  {
    signal(SIGBUS,program_trap) ;
    signal(SIGSEGV,program_trap) ;
    signal(SIGILL,program_trap) ;
    signal(SIGSYS,program_trap) ;
    signal(SIGFPE,program_trap) ;
    signal(SIGABRT,program_trap) ;
  }
}  
