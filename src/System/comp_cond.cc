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
#include <sstream>
using std::vector ;

using std::ostream ;
using std::endl ;
using std::ostringstream ;

#include "visitorabs.h"
#include "loci_globs.h"

namespace Loci {

  class execute_conditional : public execute_modules {
    executeP conditional ;
    variable cvar ;
  public:
    execute_conditional(const executeP &cond, const variable &cv) :
      conditional(cond),cvar(cv) { warn(cond==0) ; }
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
	virtual string getName() { return "execute_conditional";};
    virtual void dataCollate(collectData &data_collector) const  ;
  } ;
  
  void execute_conditional::execute(fact_db &facts, sched_db& scheds) {
    param<bool> test ;
    test = facts.get_variable(cvar) ;

    if(*test) {
      conditional->execute(facts, scheds) ;
    }
  }

  void execute_conditional::Print(ostream &s) const {
    printIndent(s) ;
    s << "if(" << cvar << ") {" << endl ;
    printLevel++ ;
    conditional->Print(s) ;
    printLevel-- ;
    printIndent(s) ;
    s << "} // if("<< cvar <<")" << endl ;
  }

  void execute_conditional::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "conditional("<<cvar<<")" ;
    int group = data_collector.openGroup(oss.str()) ;
    conditional->dataCollate(data_collector) ;
    data_collector.closeGroup(group) ;
  }

  conditional_compiler::conditional_compiler(rulecomp_map &rule_process,
					     digraph dag,
                                             variable conditional,
                                             int id):cid(id)
  {
    cond_var = conditional ;
    
    cond_gr = dag ;
    ruleSet allrules = extract_rules(cond_gr.get_all_vertices()) ;
    
    ruleSet::const_iterator ri ;
    for(ri=allrules.begin();ri!=allrules.end();++ri) {
      rulecomp_map::const_iterator rmi ;
      rmi = rule_process.find(*ri) ;
      FATAL(rmi == rule_process.end()) ;
      rule_compiler_map[*ri] = rmi->second ;
    }    
  }
  
  void conditional_compiler::accept(visitor& v) {
    v.visit(*this) ;
  }
  
  void conditional_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i)
      (*i)->set_var_existence(facts, scheds) ;
  }
  
  void conditional_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=dag_comp.rbegin();ri!=dag_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;
  }
  
  executeP conditional_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    CPTR<execute_list> elp = new execute_list ;
    
    std::vector<rule_compilerP>::iterator i ;
    for(i=dag_comp.begin();i!=dag_comp.end();++i) {
      elp->append_list((*i)->create_execution_schedule(facts, scheds)) ;
    }
	
    executeP execute = new execute_conditional(executeP(elp),cond_var);
    return  execute;
  }
}
