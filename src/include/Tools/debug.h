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
#ifndef DEBUG_H
#define DEBUG_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/debugger.h>

#include <ostream>
#include <iostream>
#ifdef NO_CSTDLIB
#include <stdlib.h>
#else
#include <cstdlib>
using std::abort ;
#endif

#ifdef DEBUG

/* 
   When DEBUG is defined, warn(bool) and fatal(bool) print warning messages 
   if the boolean argument they are passed evaluates to a true value.  
   fatal(bool) terminates the program after printing the message.  
   If DEBUG is not defined, these macros do not insert any code.
   
   These routines are used to provide some level of assurance about program
   correctness and should not be interpreted as a means of providing warnings
   or error messages to the user.
   */

#define DEBUG_TYPE_NAME(v) (typeid(v).name())

#define warn(cond) {if(cond){ std::cerr << "Warning: (" << # cond << \
                                ") true -- file \"" << __FILE__ << \
                                    "\", line " << __LINE__ << \
                    ", Class '" << DEBUG_TYPE_NAME(*this) <<"'" << std::endl;  }}

#define fatal(cond) {if(cond){ std::cerr << "FATAL Error: (" # cond << \
                                ") true -- file\"" << __FILE__ << \
                                  "\", line " << __LINE__ <<  \
                 ", Class '" << DEBUG_TYPE_NAME(*this) << "'" << std::endl; \
                                 Loci::debugger_() ; \
                                 abort() ; }}
 
#define WARN(cond) {if(cond){ std::cerr << "Warning: (" << # cond << \
                                ") true -- file \"" << __FILE__ << \
                                    "\", line " << __LINE__ << std::endl; }}
 
#define FATAL(cond) {if(cond){ std::cerr << "FATAL Error: (" # cond << \
                                ") true -- file\"" << __FILE__ << \
                                  "\", line " << __LINE__ << std::endl; \
                                 Loci::debugger_() ; \
                                 abort() ; }}


#else  // DEBUG

#define warn(cond)
#define fatal(cond)

#define WARN(cond)
#define FATAL(cond)


#endif // DEBUG

namespace Loci {
  long getmaxrss() ;
  // memory profile function
  inline double currentMem(void) { return ::Loci::getmaxrss()*1024.0 ;  }    
}

#ifdef MEMDEBUG
#define REPORTMEM() { Loci::debugout << "MEM: file " << __FILE__ << \
                                    "\", line " << __LINE__ << \
      "- Max RSS = " << ::Loci::getmaxrss()  << std::endl;  }
#else
#define REPORTMEM() 
#endif

#endif // DEBUG_H
