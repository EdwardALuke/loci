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
#include "Config/conf.h"
#include "dist_tools.h" // for use of Loci::debugout
#include "comp_tools.h"
#include "visitorabs.h"
#include "loci_globs.h"
#include "sched_tools.h"
#include "thread.h"
#include <vector>
using std::vector ;
#include <deque>
using std::deque ;
#include <map>
using std::map ;
using std::set;
#include <sstream>
using std::ostringstream ;


namespace Loci {
  extern int current_rule_id ;
  extern int chomping_size ;
  extern bool profile_memory_usage ;
  extern bool chomp_verbose;

  extern double LociAppPeakMemory ;
  extern double LociAppAllocRequestBeanCounting ;
  extern double LociAppFreeRequestBeanCounting ;
  extern double LociAppPeakMemoryBeanCounting ;
  extern double LociAppPMTemp ;

  extern bool in_internal_query;
  extern bool threading_chomping;
  extern int num_total_chomping;
  extern int num_threaded_chomping;

  execute_chomp::execute_chomp
  (const entitySet& td,
   const vector<pair<rule,rule_compilerP> >& comp,
   const deque<entitySet>& seq,
   const variableSet& cv,
   fact_db& facts):
    total_domain(td),chomp_comp(comp),rule_seq(seq),
    chomp_vars(cv),chomp_size(0),chomp_iter(0),
    D_cache_size(0), execute_times(0)
  {  
    comp_timers = vector<timeAccumulator>(comp.size()) ;
    for(vector<pair<rule,rule_compilerP> >::const_iterator vi=comp.begin();
        vi!=comp.end();++vi)
      chomp_compP.push_back(pair<int,rule_implP>
                            ((vi->first).ident(),
                             (vi->first).get_rule_implP())) ;
    
    for(vector<pair<int,rule_implP> >::iterator vi=chomp_compP.begin();
        vi!=chomp_compP.end();++vi)
      vi->second->initialize(facts) ;
    
    for(variableSet::const_iterator vi=chomp_vars.begin();
        vi!=chomp_vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      chomp_vars_rep.push_back(srp) ;
    }
    D_cache_size = chomping_size * 1024 ; // assume 128KB now
    set_seq_table();
  }
  
  void execute_chomp::set_seq_table() {
    // we'll need to set up the chomp_size
    // and the chomping sequence table
    entitySet test_alloc = interval(1,1) ;
    int_type total_obj_size = 0 ;
    for(vector<storeRepP>::iterator vi=chomp_vars_rep.begin();
        vi!=chomp_vars_rep.end();++vi) {
      int_type my_size = (*vi)->pack_size(test_alloc) ;
      total_obj_size += my_size ;
    }

    double scale = (double)D_cache_size / (double)total_obj_size ;
    
    if(scale <= 1.0)
      chomp_size = 1 ;
    else
      chomp_size = (size_t)(scale+.5) ;
    
    // this is to check that the allocated chomping domain
    // cannot exceed the original total domain of these rules
    // this is to prevent bugs in the extreme cases where
    // the container sizes are only accurately known at the
    // runtime, thus causing the analysis here to predict
    // a very large chomping allocation domain and leads to
    // extremely large allocations at the runtime.
    // for example, a storeMat<double> may have a large matrix
    // associated with each index (say 100KB), but here at
    // compile time, it only says that each index has an
    // allocation size of 4 bytes. Then if we don't check
    // the original total_domain size, we may end up computing
    // a very large chomping domain and fail the program
    if(total_domain.size() < chomp_size)
      chomp_size = total_domain.size() ;
    
    entitySet copy_total_domain = total_domain ;
    chomp_offset.clear() ;
    seq_table.clear();
    int_type low_pos ;
    int_type high_pos ;
    while(copy_total_domain != EMPTY) {
      low_pos = copy_total_domain.Min() ;
      high_pos = low_pos + chomp_size - 1 ;
      entitySet seg = interval(low_pos, high_pos) ;
      vector<entitySet> seq_vec ;
      for(deque<entitySet>::const_iterator vi=rule_seq.begin();
          vi!=rule_seq.end();++vi)
        seq_vec.push_back(*vi & seg) ;
      
      seq_table.push_back(seq_vec) ;
        
      copy_total_domain &= interval(high_pos+1,Loci::UNIVERSE_MAX) ;

      if(copy_total_domain != EMPTY)
        chomp_offset.push_back(copy_total_domain.Min() - low_pos) ;
      else
        chomp_offset.push_back(0) ;
    }

    chomp_iter = seq_table.size() ;
  }
  
    
  void execute_chomp::execute(fact_db& facts, sched_db &scheds) {
    stopWatch stot ;
    stot.start() ;

    {
      entitySet first_alloc =
        entitySet(interval(total_domain.Min(),
                           total_domain.Min()+chomp_size-1)
                  ) ;
      // first we need to allocate the chunks of variables
      for(vector<storeRepP>::iterator vi=chomp_vars_rep.begin();
          vi!=chomp_vars_rep.end();++vi) {
        (*vi)->allocate(first_alloc) ;
      }
    }
    // do memory profiling
    if(profile_memory_usage) {
      double currmen = currentMem() ;
      if(currmen > LociAppPeakMemory)
        LociAppPeakMemory = currmen ;
      LociAppAllocRequestBeanCounting += D_cache_size ;
      LociAppPMTemp += D_cache_size ;
      if(LociAppPMTemp > LociAppPeakMemoryBeanCounting)
        LociAppPeakMemoryBeanCounting = LociAppPMTemp ;
    }
   
    // begin execution, the loop number would be seq_table.size()

    // first call prelude
    vector<pair<int,rule_implP> >::iterator vri ;
    for(vri=chomp_compP.begin();vri!=chomp_compP.end();++vri) {
      vri->second->prelude(sequence(EMPTY)) ;
    }
    
    vector<vector<entitySet> >::const_iterator vvi ;
    int count = 0 ;
    for(vvi=seq_table.begin();vvi!=seq_table.end();++vvi,++count) {
      for(size_t i=0;i<chomp_compP.size();++i) {
        stopWatch s ;
        s.start() ;
        current_rule_id = chomp_compP[i].first ;
        chomp_compP[i].second->compute(sequence((*vvi)[i])) ;;
        current_rule_id = 0 ;
        comp_timers[i].addTime(s.stop(),((*vvi)[i]).size()) ;
      }
      // we shift the alloc domain for each chomp_vars_repS
      // first get the offset
      int_type offset = chomp_offset[count] ;
      if(offset != 0) {
        for(vector<storeRepP>::iterator vsi=chomp_vars_rep.begin();
            vsi!=chomp_vars_rep.end();++vsi)
          (*vsi)->shift(offset) ;
      }
    }

    //first time execute, reset_seq_table
    if(execute_times==0){
      set_seq_table();
     
    }
      
    // at last, we deallocate all the chomp_vars
    for(vector<storeRepP>::iterator vi=chomp_vars_rep.begin();
        vi!=chomp_vars_rep.end();++vi) {
      (*vi)->allocate(EMPTY) ;
    }
    
    // do memory profiling
    if(profile_memory_usage) {
      LociAppFreeRequestBeanCounting += D_cache_size ;
      LociAppPMTemp -= D_cache_size ;
      if(LociAppPMTemp < 0)
        std::cout << "MEMORY PROFILING WARNING: negative memory size"
                  << endl ;
    }
    timer.addTime(stot.stop(),1) ;
    execute_times++;
    
  }

  void execute_chomp::Print(std::ostream& s) const {
    printIndent(s) ;
    s << "--Start chomping: (chomping interval size: "
      << chomp_size << ", iter number: " << chomp_iter
      << ", total domain: " << total_domain
      << ")" << endl ;
    printIndent(s) ;
    s << "--Perform chomping for the following rule sequence: " << endl ;
    deque<entitySet>::const_iterator di=rule_seq.begin();
    for(vector<pair<rule,rule_compilerP> >::const_iterator
          vi=chomp_comp.begin();vi!=chomp_comp.end();++vi,++di) {
      printIndent(s) ;
      s << "-- " << vi->first;
      if(chomp_verbose)
        s << ", over: " << *di;
      s << endl ;
    }
    printIndent(s) ;
    s << "--End chomping" << endl ;
  }

  void execute_chomp::dataCollate(collectData &data_collector) const {
    int group = data_collector.openGroup("chomp") ;

    double tot = 0 ;
    for(size_t i=0;i<chomp_comp.size();++i) {
      ostringstream oss ;
      oss << "rule: " << chomp_comp[i].first ;

      data_collector.accumulateTime(comp_timers[i],EXEC_COMPUTATION,oss.str()) ;
      tot += comp_timers[i].getTime() ;
    }
    data_collector.closeGroup(group) ;
    group = data_collector.openGroup("chomp-overhead") ;
    timeAccumulator ov ;
    ov.addTime(timer.getTime()-tot,timer.getEvents()) ;

    ostringstream oss ;
    oss <<"chomp:" ;
    for(size_t i=0;i<chomp_comp.size();++i) 
      oss << "[" << chomp_comp[i].first  << "]";
    data_collector.accumulateTime(ov,EXEC_CONTROL,oss.str()) ;
    data_collector.closeGroup(group) ;
  }

  chomp_compiler::chomp_compiler(const digraph& cgraph,
                                 const variableSet& cvars,
                                 const map<rule,rule>& a2u)
    :chomp_graph(cgraph),chomp_vars(cvars),apply2unit(a2u) {
  }

  void chomp_compiler::set_var_existence(fact_db& facts, sched_db& scheds) {
    
    barrier_sets.clear();
    for(unsigned int i = 0; i < old_barrier_sets.size(); i++)
      barrier_sets.push_back(old_barrier_sets[i]); 
    
    vector<pair<rule,rule_compilerP> >::iterator i ;
    variableSet tvars;
    for(i=old_chomp_comp.begin();i!=old_chomp_comp.end();++i) {
      rule r = i->first ;
      rule_compilerP bc = i->second ;
      // first we check if r is a fake rule
      if(r.type() == rule::INTERNAL)
        if(r.get_info().qualifier() == "FAKE") {
          bc->set_var_existence(facts,scheds) ;
          continue ;
        }

      if(r.get_info().rule_impl->get_rule_class() != rule_impl::APPLY)
        existential_rule_analysis(r,facts,scheds) ;
      else
        existential_applyrule_analysis(r,facts,scheds) ;

      if(duplicate_work) {
	const rule_impl::info &rinfo = r.get_info().desc;
	set<vmap_info>::const_iterator si ;
	for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
	  tvars += si->var;
	}
      }
    }

    if(duplicate_work) {
      variableSet all_barrier_vars;
      for(unsigned int i = 0; i < barrier_sets.size(); i++) 
	all_barrier_vars += barrier_sets[i];

      //If work duplication is selected then making sure that model based approach
      //is not used inside chomp compiler.
      //Mainly because all variables in chomp chain should be consistent because
      //there is no communication allowed inside the chain.
      for(variableSet::const_iterator vi = all_barrier_vars.begin();
	  vi != all_barrier_vars.end(); vi++)
	scheds.add_policy(*vi, sched_db::ALWAYS);
   
      tvars -= all_barrier_vars;
      barrier_sets.push_back(tvars);
      all_barrier_vars += tvars;

      //To find which are variables to connected through rules
      //Each element of vector represent a set of variables which are connected
      //through rules of chomp compiler.
      //These variables are either barrier_variables or target variables
      vector<variableSet> variable_dependancy; 
      map<variable, bool> duplication_map; 
      for(int i = barrier_sets.size() - 1; i >= 0; i--) {
 	std::map<variable, ruleSet> var2rules;
 	for(variableSet::const_iterator vi = barrier_sets[i].begin();
 	    vi != barrier_sets[i].end(); vi++) {
 	  var2rules[*vi] = scheds.get_existential_rules(*vi);
 	}

	for(variableSet::const_iterator vi = barrier_sets[i].begin(); vi != barrier_sets[i].end(); vi++) {
	  if(!tvars.inSet(*vi))
	    duplication_map[*vi] = process_policy_duplication(*vi, scheds, facts);
	  bool found = false;  //To check if variable already exist in variable_dependancy
	  for(unsigned int j = 0; j < variable_dependancy.size(); j++) {
	    if(variable_dependancy[j].inSet(*vi)) {
	      found = true;
	      for(ruleSet::const_iterator rsi = var2rules[*vi].begin();
		  rsi != var2rules[*vi].end(); ++rsi) {
		variableSet temp = input_variables(*rsi);
		for(variableSet::const_iterator vsi = temp.begin(); vsi != temp.end(); vsi++)
		  variable_dependancy[j] += scheds.get_synonyms(*vsi) & all_barrier_vars;
		  
	      }
	    }
	  }

	  if(!found) {
	    variableSet tmpVars;
	    tmpVars += *vi;
	    for(ruleSet::const_iterator rsi = var2rules[*vi].begin();
		rsi != var2rules[*vi].end(); ++rsi) {
	      variableSet temp = input_variables(*rsi);
	      for(variableSet::const_iterator vsi = temp.begin(); vsi != temp.end(); vsi++)
		tmpVars += scheds.get_synonyms(*vsi) & all_barrier_vars;
	      
	    }

	    variable_dependancy.push_back(tmpVars);
	  }
	}
      }
      for(unsigned int i = 0; i < variable_dependancy.size(); i++) {
	bool alltrue = true;  //To check if each variable in chomp chain is capable of duplication 
	for(variableSet::const_iterator vi = variable_dependancy[i].begin();
	    vi != variable_dependancy[i].end(); vi++) 
	  if(!tvars.inSet(*vi))
	    alltrue &= duplication_map[*vi];

	//If any of the variables in chomp chain cannot be duplicated then
	//all of the variables should not be duplicated including the target variables
	//which are using chomp variables as input
	if(!alltrue)
	  for(variableSet::const_iterator vi = variable_dependancy[i].begin();
	      vi != variable_dependancy[i].end(); vi++)
	    scheds.add_policy(*vi, sched_db::NEVER);
      }
    }
  }
  
  void chomp_compiler::process_var_requests(fact_db& facts, sched_db& scheds) {
    rule_seq.clear();
    deque<pair<rule,rule_compilerP> > new_chomp_comp ;
    vector<pair<rule,rule_compilerP> >::reverse_iterator ri ;
    for(ri=old_chomp_comp.rbegin();ri!=old_chomp_comp.rend();++ri) {
      rule r = ri->first ;
      rule_compilerP bc = ri->second ;

      // first we check if r is a fake rule
      if(r.type() == rule::INTERNAL)
        if(r.get_info().qualifier() == "FAKE") {
          bc->process_var_requests(facts,scheds) ;
          continue ;
        }
      
      entitySet exec_seq ;
      if(r.get_info().rule_impl->get_rule_class() != rule_impl::APPLY)
        exec_seq = process_rule_requests(r,facts,scheds) ;
      else {
        map<rule,rule>::const_iterator mi ;
        mi = apply2unit.find(r) ;
        FATAL(mi == apply2unit.end()) ;
	//We just need this variable to pass to the process_applyrule_requests. 
	bool output_mapping; 
        exec_seq = process_applyrule_requests(r,mi->second,output_mapping, facts,scheds) ;
      }

      new_chomp_comp.push_front(*ri) ;
      rule_seq.push_front(exec_seq) ;
    }
    // finally update the real chomp_comp list
    chomp_comp.clear() ;
    for(deque<pair<rule,rule_compilerP> >::const_iterator
          di=new_chomp_comp.begin();di!=new_chomp_comp.end();++di)
      chomp_comp.push_back(*di) ;
  }

  executeP chomp_compiler::create_execution_schedule(fact_db& facts,
                                                     sched_db& scheds) {
    // we first union all the rule sequence
    entitySet total ;
    
    deque<entitySet>::size_type i ;
    
    for(i=0;i<rule_seq.size();++i)
      total += rule_seq[i] ;

#ifdef PTHREADS
    ++num_total_chomping;
    if(!in_internal_query && threading_chomping) {
      int tnum = thread_control->num_threads();
      int minw = thread_control->min_work_per_thread();
      // no multithreading if the execution sequence is too small
      if(total.size() < (size_t)tnum*minw)
        // normal case
        return
          new execute_chomp(total,chomp_comp,rule_seq,chomp_vars,facts);
      else {
        // generate multithreaded execution module
        ++num_threaded_chomping;
        return new
          Threaded_execute_chomp(sequence(total),chomp_comp,rule_seq,
                                 chomp_vars,facts,scheds);
      }      
    } else {
#endif
      executeP execute =
        new execute_chomp(total,chomp_comp,rule_seq,chomp_vars,facts);
      return  execute;
#ifdef PTHREADS
    }
#endif
  }
  
} // end of namespace Loci
