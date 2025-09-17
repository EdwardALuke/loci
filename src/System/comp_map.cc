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
  
  class execute_map_rule: public execute_modules {
    rule_implP rp ;
    rule rule_tag ;
    sequence exec_seq ;
    variableSet sources ; // source vars of the map rule
    variableSet targets ; // target vars of the map rule
    timeAccumulator timer ;
  public:
    // tsv is the target and source variables of this map rule
    execute_map_rule(rule fi, sequence seq,
                     const variableSet& svars,
                     const variableSet& tvars,
                     fact_db &facts, sched_db &scheds) ;
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() { return "execute_map_rule";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;
  
  execute_map_rule::
  execute_map_rule(rule fi, sequence seq,
                   const variableSet& svars,
                   const variableSet& tvars,
                   fact_db& facts, sched_db& scheds) {
    rp = fi.get_rule_implP() ;
    sources = svars ;
    targets = tvars ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
  }
  
  void execute_map_rule::execute(fact_db &facts, sched_db &scheds) {
    stopWatch s ;
    s.start() ;
    
    current_rule_id = rule_tag.ident() ;
    // before the execute begins, we need to restore
    // all the facts associated with this map rule
    // to their global numbering scheme
    if(facts.is_distributed_start()) {
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      dMap dl2g ;
      dl2g = MapRepP(df->l2g.Rep())->thaw() ;
      for(variableSet::const_iterator vi=sources.begin();
          vi!=sources.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        // Note, we can use facts.replace_fact here
        // But doing so will require rebind the rule
        // to the fact_db
        facts.update_fact(*vi,srp->remap(dl2g)) ;
        //facts.replace_fact(*vi,srp->remap(l2g)) ;
      }
    }
    rp->compute(exec_seq) ;
    // and then after computing, we'll convert all
    // all the facts associated with this map rule
    // to their local numbering scheme
    // But since we have generated new maps
    // we could possibly end up having new entities
    // included in the clone region on local process.
    // Therefore we should only transform all the source
    // facts since transforming the generated maps
    // at last will likely to lose data.
    if(facts.is_distributed_start()) {
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      dMap g2l ;
      entitySet ldom = df->l2g.domain() ;
      FORALL(ldom,i) {
        g2l[df->l2g[i]] = i ;
      } ENDFORALL ;
      for(variableSet::const_iterator vi=sources.begin();
          vi!=sources.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        facts.update_fact(*vi,srp->remap(g2l)) ;
      }
    }
    current_rule_id = 0 ;
    timer.addTime(s.stop(),1) ;
  }
  
  void execute_map_rule::Print(ostream &s) const {
    s << rule_tag << "  over sequence " << exec_seq << endl ;
  }
  
  void execute_map_rule::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "map rule: " << rule_tag ;
    data_collector.accumulateTime(timer,EXEC_CONTROL,oss.str()) ;
  }

  void map_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    existential_rule_analysis(map_impl, facts, scheds) ;
  }
  
  void map_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    entitySet exec_seq = process_rule_requests(map_impl, facts, scheds) ;
    scheds.update_exec_seq(map_impl, exec_seq);
  }
  
  executeP map_compiler::create_execution_schedule(fact_db& facts, sched_db& scheds) {

    entitySet exec_seq = scheds.get_exec_seq(map_impl);
    executeP execute = new execute_map_rule(map_impl, exec_seq,
					    map_impl.sources(),
                                            map_impl.targets(), facts,scheds);
    return execute;
  }
  
}
