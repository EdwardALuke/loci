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
#include <ostream>
#include <execute.h>
#include "loci_globs.h"
#include <distribute.h>

using std::cerr ;
using std::endl ;
using std::ostream ;


namespace Loci {
  struct exec_info  {
    Loci::executeP exec_routine;
    Loci::fact_db *current_fact_db ;
    exec_info() {} ;
    exec_info(Loci::executeP &ep, Loci::fact_db &facts) {
      exec_routine = ep ; current_fact_db = &facts ;
    }
  } ;

  Loci::fact_db *current_fact_db ;

  void execute_list::execute(fact_db &facts, sched_db& scheds) {
    std::vector<executeP>::iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->execute(facts, scheds) ;
  }

  void execute_list::Print(std::ostream &s) const {
    std::vector<executeP>::const_iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->Print(s) ;
  }

  void execute_list::dataCollate(collectData &data_collector) const {
    std::vector<executeP>::const_iterator eli ;

    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->dataCollate(data_collector) ;
  }


  void execute_sequence::execute(fact_db &facts, sched_db& scheds) {
    std::vector<executeP>::iterator eli ;

    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->execute(facts, scheds) ;
  }

  void execute_sequence::Print(std::ostream &s) const {
    std::vector<executeP>::const_iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->Print(s) ;
  }

  void execute_sequence::dataCollate(collectData &data_collector) const {
    std::vector<executeP>::const_iterator eli ;

    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->dataCollate(data_collector) ;
  }

}
