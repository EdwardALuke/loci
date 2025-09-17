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
#include <vector>

using std::vector ;

#include <sys/stat.h>

#include "dist_tools.h"
#include <parameter.h>

#include <map>
using std::map ;

#include "visitorabs.h"

namespace Loci {
  
  dag_compiler::dag_compiler(rulecomp_map &rule_process, digraph dag, int id):cid(id) {

    ////////////////////
    // store the dag structure and the relevant rulecompiler map
    dag_gr = dag ;
    ruleSet allrules = extract_rules(dag_gr.get_all_vertices()) ;
    
    ruleSet::const_iterator ri ;
    for(ri=allrules.begin();ri!=allrules.end();++ri) {
      rulecomp_map::const_iterator rmi ;
      rmi = rule_process.find(*ri) ;
      FATAL(rmi == rule_process.end()) ;
      rule_compiler_map[*ri] = rmi->second ;
    }
    ////////////////////
  }

  void dag_compiler::accept(visitor& v) {
    v.visit(*this) ;
  }

  void dag_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i)
      (*i)->set_var_existence(facts, scheds) ;
  }
  
  void dag_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=dag_comp.rbegin();ri!=dag_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;
  }
  
  executeP dag_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    CPTR<execute_list> elp = new execute_list ;
    
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i) {
      elp->append_list((*i)->create_execution_schedule(facts, scheds)) ;
    }

    return executeP(elp) ;
  }
}

