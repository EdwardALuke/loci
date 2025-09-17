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
#include "dist_tools.h"
using std::ostream ;
using std::endl ;
#include <sstream>
using std::ostringstream ;

namespace Loci {

  extern int current_rule_id ;

  class execute_constraint_rule: public execute_modules {
    rule_implP rp ;
    rule rule_tag ;
    sequence exec_seq ;
    timeAccumulator timer ;
  public:
    execute_constraint_rule(rule fi, sequence seq,
                            fact_db &facts, sched_db &scheds) ;
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() { return "execute_constraint_rule";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  execute_constraint_rule::
  execute_constraint_rule(rule fi, sequence seq,
                          fact_db& facts, sched_db& scheds) {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
  }
  
  void execute_constraint_rule::execute(fact_db &facts, sched_db &scheds) {
    stopWatch s ;
    s.start() ;
    current_rule_id = rule_tag.ident() ;
    rp->compute(exec_seq) ;
    current_rule_id = 0 ;
    timer.addTime(s.stop(),1) ;
  }
  
  void execute_constraint_rule::Print(ostream &s) const {
    s << rule_tag << "  over sequence " << exec_seq << endl ;
  }
  
  void execute_constraint_rule::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "constraint rule: " << rule_tag;

    data_collector.accumulateTime(timer,EXEC_CONTROL,oss.str()) ;
  }

  void constraint_compiler::set_var_existence(fact_db& facts,
                                              sched_db& scheds)
  {
    existential_rule_analysis(constraint_rule,facts, scheds) ;
  }
  
  void constraint_compiler::process_var_requests(fact_db& facts,
                                                 sched_db& scheds) {
    variableSet var_requests = constraint_rule.targets() ;
    variableSet::const_iterator vi ;
    // all constraints would need to request for everything
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      scheds.variable_request(*vi,~EMPTY) ;//scheds.variable_existence(*vi)) ;
    }
    
    entitySet exec_seq = process_rule_requests(constraint_rule,facts, scheds) ;
    scheds.update_exec_seq(constraint_rule, exec_seq);
  }
  
  executeP constraint_compiler::create_execution_schedule(fact_db& facts,
                                                          sched_db& scheds) {
    entitySet exec_seq = scheds.get_exec_seq(constraint_rule);
    executeP execute = new execute_constraint_rule(constraint_rule,
                                                   exec_seq, facts,scheds) ;
    return execute;
  }

}
