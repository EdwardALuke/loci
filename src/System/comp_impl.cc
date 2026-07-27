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
#include "loci_globs.h"
#include "thread.h"
#include <sstream>
#include <set>
#include <map>

using std::ostream ;
using std::endl ;
using std::set ;
using std::map ;
using std::ostringstream ;

namespace Loci {
  //#define DYNAMIC_TIMING ;
#ifdef DYNAMIC_TIMING
#include <execute.h>
  // performance measure
  stopWatch sw_expand ;
  stopWatch sw_expand2 ;
  stopWatch sw_expand_start ;
  stopWatch sw_expand_cache ;
  stopWatch sw_expand_collect_img ;
  stopWatch sw_expand_missing ;
  stopWatch sw_context ;
  stopWatch sw_context_nonepd ;
  stopWatch sw_context_nonepd_domt ;
  stopWatch sw_context_nonepd_ints ;
  stopWatch sw_context_epdend ;
  stopWatch sw_context_epdmid ;
  stopWatch sw_context_epdsta ;
  stopWatch sw_output_oh ;
  stopWatch sw_compute ;
  stopWatch sw_erase ;
  stopWatch sw_invalidate ;
  stopWatch sw_keyremoval ;
  stopWatch sw_insertion ;
  stopWatch sw_keyinsert ;
  stopWatch sw_key_dist ;
  stopWatch sw_dist_renumber ;
  stopWatch sw_dist ;
  stopWatch sw_push ;
  stopWatch sw_expand_comm ;
  stopWatch sw_pw_push ;
  stopWatch sw_param_push ;
  stopWatch sw_param_pack ;
  stopWatch sw_param_unpack ;
  stopWatch sw_param_reduce ;
  stopWatch sw_record_erase ;
  stopWatch sw_dctrl ;

  timeAccumulator ta_expand ;
  timeAccumulator ta_expand2 ;
  timeAccumulator ta_expand_start ;
  timeAccumulator ta_expand_cache ;
  timeAccumulator ta_expand_collect_img ;
  timeAccumulator ta_expand_missing ;
  timeAccumulator ta_context ;
  timeAccumulator ta_context_nonepd ;
  timeAccumulator ta_context_nonepd_domt ;
  timeAccumulator ta_context_nonepd_ints ;
  timeAccumulator ta_context_epdend ;
  timeAccumulator ta_context_epdmid ;
  timeAccumulator ta_context_epdsta ;
  timeAccumulator ta_output_oh ;
  timeAccumulator ta_compute ;
  timeAccumulator ta_erase ;
  timeAccumulator ta_invalidate ;
  timeAccumulator ta_keyremoval ;
  timeAccumulator ta_insertion ;
  timeAccumulator ta_keyinsert ;
  timeAccumulator ta_key_dist ;
  timeAccumulator ta_dist_renumber ;
  timeAccumulator ta_dist ;
  int             ta_dist_number = 0 ;
  timeAccumulator ta_push ;
  timeAccumulator ta_expand_comm ;
  timeAccumulator ta_pw_push ;
  timeAccumulator ta_param_push ;
  timeAccumulator ta_param_pack ;
  timeAccumulator ta_param_unpack ;
  timeAccumulator ta_param_reduce ;
  int             ta_param_reduce_num = 0 ;
  timeAccumulator ta_record_erase ;
  timeAccumulator ta_dctrl ;
  int             ta_drule_executes = 0 ;
  int             ta_dctrl_executes = 0 ;
#endif

  extern bool in_internal_query;
  extern bool threading_pointwise;
  extern int num_threaded_pointwise;
  extern int num_total_pointwise;

  int current_rule_id = 0 ;
  int rule_count = 0;

  
  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts, const sched_db &scheds)  {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    exec_seq = seq ;
    exec_size = seq.size() ;
    //
#ifdef PAPI_DEBUG
    papi_events[0] = PAPI_L1_DCM;
    papi_events[1] = PAPI_L2_DCM;
    papi_values[0] = papi_values[1] = 0;
    l1_dcm = l2_dcm = 0;
#endif
  }

  execute_rule::execute_rule(rule fi, sequence seq, fact_db &facts,
                             variable v, const storeRepP &p, const sched_db &scheds) {
    rp = fi.get_rule_implP() ;
    rule_tag = fi ;
    rp->initialize(facts) ;
    rp->set_store(v,p) ;
    exec_seq = seq ;
    exec_size = seq.size() ;
#ifdef PAPI_DEBUG
    papi_events[0] = PAPI_L1_DCM;
    papi_events[1] = PAPI_L2_DCM;
    papi_values[0] = papi_values[1] = 0;
    l1_dcm = l2_dcm = 0;
#endif
  }

  void execute_rule::execute(fact_db &facts, sched_db &scheds) {
#ifdef PAPI_DEBUG
    if( (PAPI_start_counters(papi_events,2)) != PAPI_OK) {
      cerr << "PAPI failed to start counters" << endl;
      Loci::Abort();
    }
#endif

    stopWatch s ;
    s.start() ;
    current_rule_id = rule_tag.ident() ;
#ifdef VERBOSE
    Loci::debugout << "executing " << rule_tag << endl ;
#endif
    //rp->compute(exec_seq);
    execute_prelude(exec_seq);
    execute_kernel(exec_seq);
    execute_postlude(exec_seq);
    current_rule_id = 0 ;
    timer.addTime(s.stop(),exec_size) ;
#ifdef PAPI_DEBUG
    if( (PAPI_stop_counters(papi_values,2)) != PAPI_OK) {
      cerr << "PAPI failed to read counters" << endl;
      Loci::Abort();
    }
    l1_dcm += papi_values[0];
    l2_dcm += papi_values[1];
#endif
  }

  inline void 
  execute_rule::execute_kernel(const sequence& s)
  { rp->compute(s); }
  
  inline void
  execute_rule::execute_prelude(const sequence& s)
  { rp->prelude(s); }

  inline void
  execute_rule::execute_postlude(const sequence& s)
  { rp->postlude(s); }

  void execute_rule::Print(ostream &s) const {
    printIndent(s) ;
    s << rule_tag << " over sequence " ;
    if(verbose || exec_seq.num_intervals() < 4) {
      s << exec_seq << endl ;
    } else {
      s << "[ ... ], l=" << exec_seq.size() << endl ;
    }
  }

  void execute_rule::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "rule: "<<rule_tag ;

    data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;

#ifdef PAPI_DEBUG
    data_collector.accumulateSchedMemory(oss.str(),
                                         8*exec_seq.num_intervals());

    data_collector.accumulateDCM(oss.str(),
                                 papi_values[0], papi_values[1]);
#endif
  }

  void impl_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    existential_rule_analysis(impl,facts, scheds) ;
  }

  void impl_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    entitySet exec_seq = process_rule_requests(impl,facts, scheds) ;
    scheds.update_exec_seq(impl, exec_seq);
  }

  executeP impl_compiler::
  create_execution_schedule(fact_db &facts,sched_db &scheds) {
    //    if(GLOBAL_AND(exec_seq.size()==0)) {
    //      return executeP(0) ;
    //    }
    variableSet targets = impl.targets() ;
    WARN(targets.size() == 0) ;

    extern int method ;
    entitySet exec_seq = scheds.get_exec_seq(impl);
    if (impl.get_info().rule_impl->dynamic_schedule_rule() &&
        use_dynamic_scheduling) {
      executeP execute_dynamic =
        new dynamic_schedule_rule(impl,exec_seq,facts, scheds,method) ;

      return execute_dynamic;
    }

#ifdef PTHREADS
    ++num_total_pointwise;
    bool threadable = 
      impl.get_info().rule_impl->thread_rule() &&
      (impl.get_info().rule_impl->get_rule_class() 
                       == rule_impl::POINTWISE ||
       impl.get_info().rule_impl->get_rule_class()
                       == rule_impl::UNIT);
    rule_implP ti = impl.get_rule_implP() ;
    for (variableSet::const_iterator vi=targets.begin();
        vi!=targets.end();++vi) {
      storeRepP tr = ti->get_store(*vi);
      if (isPARAMETER(tr)) {
        threadable = false;
        break;
      }
    }
    if(!in_internal_query && threading_pointwise && threadable) {
      int tnum = thread_control->num_threads();
      int minw = thread_control->min_work_per_thread();
      // if a rule is not for threading, then generate a normal module,
      // also no multithreading if the execution sequence is too small
      if(exec_seq.size() < (size_t)tnum*minw)
        // normal case
        return new execute_rule(impl,sequence(exec_seq),facts, scheds);
      else {
        // generate multithreaded execution module
        ++num_threaded_pointwise;
        return new Threaded_execute_rule(impl, exec_seq, facts, scheds);
      }
    }
#endif
    // normal case
    executeP exec_rule =
      new execute_rule(impl,sequence(exec_seq),facts, scheds);

    return exec_rule;
  }


  
  // blackbox_compiler code
  void
  blackbox_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    //    existential_blackboxrule_analysis(impl, facts, scheds) ;
    // set UNIVERSE existence for all targets
     variableSet targets = impl.targets() ;
     for(variableSet::const_iterator vi=targets.begin();vi!=targets.end();++vi)
       scheds.set_existential_info(*vi, impl, ~EMPTY) ;
  }

  void blackbox_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    //    exec_seq = process_blackboxrule_requests(impl, facts, scheds) ;

    //everyone will need to request for their existence
    variableSet targets = impl.targets() ;
    for(variableSet::const_iterator vi = targets.begin(); vi != targets.end(); ++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;

    variableSet sources = impl.sources() ;
    for(variableSet::const_iterator vi=sources.begin(); vi != sources.end(); ++vi)
      scheds.variable_request(*vi, scheds.variable_existence(*vi)) ;
    entitySet exec_seq = ~EMPTY ;
    scheds.update_exec_seq(impl, exec_seq);
  }

  executeP blackbox_compiler::
  create_execution_schedule(fact_db& facts, sched_db& scheds) {
    entitySet exec_seq = scheds.get_exec_seq(impl);
    executeP execute = new execute_rule(impl,
                                        sequence(exec_seq), facts, scheds);
    return execute ;
  }

  // superRule_compiler code
  void
  superRule_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    CPTR<super_rule> rp(impl.get_rule_implP()) ;
    rp->process_existential(impl,facts,scheds) ;
  }

  void superRule_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    CPTR<super_rule> rp(impl.get_rule_implP()) ;
    rp->process_requests(impl,facts,scheds) ;
  }

  executeP superRule_compiler::create_execution_schedule(fact_db& facts, sched_db& scheds) {
    executeP execute = new execute_rule(impl, ~EMPTY, facts, scheds);
    return execute;
  }


  ///////////////////////////////////////////////
  // Lets set up some common super rule functions
  namespace {
    const char *systemVarDocs[] = {
      "SYSTEM\000NOT(X)\000param<bool>\000Exists over complement of X",
      "SYSTEM\000OR(X,Y)\000param<bool>\000Exists over union of X and Y",
      "SYSTEM\000OR(X,Y,Z)\000param<bool>\000Exists over union of X, Y and Z",
      "SYSTEM\000OR(W,X,Y,Z)\000param<bool>\000Exists over union of W, X, Y and Z",
      "SYSTEM\000AND(X,Y)\000param<bool>\000Exists over intersection of X and Y",
      "SYSTEM\000AND(X,Y,Z)\000param<bool>\000Exists over intersection of X, Y and Z",
      "SYSTEM\000AND(W,X,Y,Z)\000param<bool>\000Exists over intersection of W, X, Y and Z",
      "\000\000\000"
    } ;
  }
  class NOT_rule : public super_rule {
    param<bool> NOT ;
  public:
    NOT_rule() {
      comments("Internal System Rule.  Use with care to prevent infinite regression!") ;
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("NOT(X)",0) ;
      name_store("NOT(X)",NOT) ;
      constraint("X") ;
      output("NOT(X)") ;
    }
    void compute(const sequence &seq) {
      *NOT = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = ~EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints &= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      //      debugout << "constraints = " << constraints << endl ;
      // Now complement (entities we own)
      constraints = (~constraints) & my_entities ;
      //      debugout << "constraints & my_entities =" << constraints << endl ;

      variableSet output = r.targets() ;
      //      debugout << "setting " << output << endl ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<NOT_rule> register_NOT_rule ;

  class OR_rule : public super_rule {
    param<bool> OR ;
  public:
    OR_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("OR(X,Y)",1) ;
      name_store("OR(X,Y)",OR) ;
      constraint("X,Y") ;
      output("OR(X,Y)") ;
    }
    void compute(const sequence &seq) {
      *OR = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints |= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      // Now complement (entities we own)
      constraints = constraints & my_entities ;

      variableSet output = r.targets() ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<OR_rule> register_OR_rule ;

  class OR3_rule : public super_rule {
    param<bool> OR ;
  public:
    OR3_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("OR(X,Y,Z)",2) ;
      name_store("OR(X,Y,Z)",OR) ;
      constraint("X,Y,Z") ;
      output("OR(X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *OR = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints |= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      // Now complement (entities we own)
      constraints = constraints & my_entities ;

      variableSet output = r.targets() ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<OR3_rule> register_OR3_rule ;

  class OR4_rule : public super_rule {
    param<bool> OR ;
  public:
    OR4_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("OR(W,X,Y,Z)",3) ;
      name_store("OR(W,X,Y,Z)",OR) ;
      constraint("W,X,Y,Z") ;
      output("OR(W,X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *OR = true ;
    }
    void process_existential(rule r, fact_db &facts, sched_db &scheds) {
      const rule_impl::info &rinfo = r.get_info().desc ;
      set<vmap_info>::const_iterator si ;

      entitySet constraints = EMPTY ;
      for(si=rinfo.constraints.begin();si!=rinfo.constraints.end();++si)
        constraints |= vmap_source_exist(*si,facts, scheds) ;

      entitySet my_entities = ~EMPTY ;

      if(facts.isDistributed()) {
      // For the distributed memory case we restrict the sources and
      // constraints to be within my_entities.
        fact_db::distribute_infoP d = facts.get_distribute_info() ;
        my_entities = d->my_entities ;
      }
      
      // Now complement (entities we own)
      constraints = constraints & my_entities ;

      variableSet output = r.targets() ;
      for(variableSet::const_iterator vi=output.begin();vi!=output.end();++vi)
        scheds.set_existential_info(*vi, r, constraints) ;
    }
    void process_requests(rule r, fact_db &facts, sched_db &scheds) {
    }
  } ;

  register_rule<OR4_rule> register_OR4_rule ;

  class AND2_rule : public singleton_rule {
    param<bool> AND ;
  public:
    AND2_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("AND(X,Y)",4) ;
      name_store("AND(X,Y)",AND) ;
      constraint("X,Y") ;
      output("AND(X,Y)") ;
    }
    void compute(const sequence &seq) {
      *AND = true ;
    }
  } ;

  register_rule<AND2_rule> register_AND2_rule ;

  class AND3_rule : public singleton_rule {
    param<bool> AND ;
  public:
    AND3_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("AND(X,Y,Z)",5) ;
      name_store("AND(X,Y,Z)",AND) ;
      constraint("X,Y,Z") ;
      output("AND(X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *AND = true ;
    }
  } ;

  register_rule<AND3_rule> register_AND3_rule ;

  class AND4_rule : public singleton_rule {
    param<bool> AND ;
  public:
    AND4_rule() {
      setvardoc(systemVarDocs) ;
      set_file( __FILE__ ) ;
      store_info_id("AND(W,X,Y,Z)",6) ;
      name_store("AND(W,X,Y,Z)",AND) ;
      constraint("W,X,Y,Z") ;
      output("AND(W,X,Y,Z)") ;
    }
    void compute(const sequence &seq) {
      *AND = true ;
    }
  } ;

  register_rule<AND4_rule> register_AND4_rule ;
    
  
}

// ... the end ...
