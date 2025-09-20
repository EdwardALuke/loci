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
#include "comp_tools.h"

// flags to enable(disable) writing out the
// these internals to the schedule file
#define WRITE_OUT
//#define DONT_WRITE_OUT

namespace Loci {

  void promote_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
  }

  void promote_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
  }

  executeP promote_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
#ifdef DONT_WRITE_OUT
    return executeP(0) ;
#endif
#ifdef WRITE_OUT
    std::ostringstream oss ;
    oss << r ;
    executeP execute = executeP(new execute_msg(oss.str())) ;
    return execute;	
#endif
  }

  void generalize_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
  }

  void generalize_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
  }

  executeP generalize_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
#ifdef DONT_WRITE_OUT
    return executeP(0) ;
#endif
#ifdef WRITE_OUT
    std::ostringstream oss ;
    oss << r ;
    executeP execute = executeP(new execute_msg(oss.str())) ;
    return execute;
#endif
  }

  void priority_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
  }

  void priority_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
  }

  executeP priority_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
#ifdef DONT_WRITE_OUT
    return executeP(0) ;
#endif
#ifdef WRITE_OUT
    std::ostringstream oss ;
    oss << r ;
    executeP execute = executeP(new execute_msg(oss.str())) ;
    return execute;
#endif
  }
}
