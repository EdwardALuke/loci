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
//#define VERBOSE
#include "comp_tools.h"
#include <vector>
#include <sstream>
using std::vector ;
using std::list ;
using std::map ;

using std::ostream ;
using std::endl ;
using std::ostringstream ;

#include "visitorabs.h"
#include "loci_globs.h"
namespace Loci {
  int printLevel = 0 ;
  
  class execute_loop : public execute_modules {
    executeP collapse, advance ;
    variable cvar ;
    variable tvar ;
    time_ident tlevel ;
    list<list<variable> > rotate_lists ;
  public:
    execute_loop(const variable &cv,
                 const executeP &col, const executeP &adv,
                 const time_ident &tl, 
                 list<list<variable> > &rl) :
      collapse(col),advance(adv),
      cvar(cv),
      tlevel(tl),rotate_lists(rl) {
      warn(col==0 || advance==0) ; tvar = variable(tlevel) ;
    }
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() { return "execute_loop";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;
  
  void execute_loop::execute(fact_db &facts, sched_db &scheds) {
    param<bool> test ;
    test = facts.get_variable(cvar) ;
    // initialize conditional variables to true
    //    *test = true ;
    
    param<int> time_var ;
    time_var = facts.get_variable(tvar) ;
    // Start iteration by setting iteration variable to zero 


    *time_var = 0 ;

    for(;;) { // Begin looping
      //First evaluate the collapse condition
      collapse->execute(facts, scheds) ;

#ifdef VERBOSE
      debugout << cvar << " = " << *test << endl ;
      debugout << tvar << " = " << *time_var << endl ;
#endif

      if(*test) {// If collapse condition satisfied, were finished
        return ;
      }
      // Advance to the next iterate
      advance->execute(facts, scheds) ;

      // We are finished with this iteration, so rotate variables and
      // add one to the iterate
      list<list<variable> >::const_iterator rli ;
      for(rli = rotate_lists.begin();rli!=rotate_lists.end();++rli) {
#ifdef VERBOSE
        debugout << "rotating [ " ;
        for(list<variable>::const_iterator rlii=rli->begin();rlii!=rli->end();
            ++rlii) {
          debugout << *rlii << ' ' ;
        }
        debugout << "]" << endl ;
#endif
#ifdef DYNAMICSCHEDULING
        // before rotation, we need to adjust the
        // history variables (if necessary)
        facts.adjust_rotation_vars(*rli) ;
#endif
        // then rotate them
        facts.rotate_vars(*rli) ;
      }
      *time_var += 1 ;
    }
  }
  
  void execute_loop::Print(ostream &s) const {
    printIndent(s) ;
    s << "Iteration Loop{"<< tlevel << "} {" << endl ;
    printLevel+=2 ;
    //    printIndent(s) ;
    //    s << "--compute collapse of {" <<tlevel << "} , conditional on " << cvar << endl ;
    //    printIndent(s) ;
    //    s << "--collapse iteration {" << tlevel << "}"<< endl ;
    collapse->Print(s) ;
    s<< endl ;
    printIndent(s) ;
    s << "-------------- Exit of Loop{" << tlevel << "}"<< endl ;
    printIndent(s) ;
    s << "if(" << cvar << ") break ;" << endl ;
    s << endl ;
    //    printIndent(s) ;
    //    s << "--advance iteration {" << tlevel << "}" << endl ;
    advance->Print(s) ;
    printLevel-=2;
    printIndent(s) ;
    s << "} // {" << tlevel  << "}" << endl ;
  }

  void execute_loop::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "iteration("<<tlevel<<")"  ;
    int group = data_collector.openGroup(oss.str()) ;
    advance->dataCollate(data_collector) ;
    collapse->dataCollate(data_collector) ;
    data_collector.closeGroup(group) ;
  }
  
  inline bool offset_sort(const variable &v1, const variable &v2)
  { return v1.get_info().offset > v2.get_info().offset ; }
  
  loop_compiler::loop_compiler(rulecomp_map &rule_process, digraph dag, int id):cid(id) {
    ////////////////////
    // store the graph structure and the relevant rulecompiler map
    loop_gr = dag ;
    ruleSet allrules = extract_rules(loop_gr.get_all_vertices()) ;
    
    for(ruleSet::const_iterator ri=allrules.begin();ri!=allrules.end();++ri) {
      rulecomp_map::const_iterator rmi ;
      rmi = rule_process.find(*ri) ;
      FATAL(rmi == rule_process.end()) ;
      rule_compiler_map[*ri] = rmi->second ;
    }
    ////////////////////
    ruleSet collapse_rules, all_rules, loopset ;

    all_rules = extract_rules(loop_gr.get_all_vertices()) ;
    for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) {
      if(ri->get_info().qualifier() == "looping")
        loopset += ri->ident() ;
      if(ri->target_time().before(ri->source_time()))
        collapse_rules += *ri ;
    }

    if(loopset.size() != 1) {
      cerr << "internal consistency error, loopset size != 1" << endl ;
      cerr << "looping rules = " << loopset << endl ;
      Loci::Abort() ;
    }

    tlevel = loopset.begin()->source_time() ;
    loop_gr.remove_vertex((*loopset.begin()).ident()) ;
    loop_gr.remove_dangling_vertices() ;

    
    variable OUTPUT(variable("OUTPUT"),tlevel) ;

    digraph loop_grt = loop_gr.transpose() ;
  
    if(collapse_rules.size() != 1 ||
       collapse_rules.begin()->get_info().desc.conditionals.size() != 1 ) {
      cerr << "collapse for loop at iteration level " << tlevel << " ill-formed"
           << endl << "error is not recoverable" << endl ;
      cerr << "rules = " << collapse_rules << endl ;
      Loci::Abort() ;
    }

    variableSet input_vars ;
    variableSet all_vars = extract_vars(loop_gr.get_all_vertices()) ;
    for(variableSet::const_iterator vi=all_vars.begin();
        vi!=all_vars.end();++vi ) {

        if(vi->get_info().offset == 1)
          advance_vars += *vi ;
        if(loop_grt[vi->ident()] == EMPTY)
          input_vars += *vi ;
    }
        
  
    collapse_vars = collapse_rules.begin()->targets() ;
    cond_var = *(collapse_rules.begin()->get_info().desc.conditionals.begin());

    variableSet outcond ;
    ruleSet outrules = extract_rules(loop_grt[OUTPUT.ident()]) ;
    for(ruleSet::const_iterator ri=outrules.begin();ri!=outrules.end();++ri) 
      outcond += ri->get_info().desc.conditionals ;

#ifdef VERBOSE
    debugout << "loop: collapse vars = " << collapse_vars << endl ;
    debugout << "loop: input_vars = " << input_vars << endl ;
    debugout << "loop: outcond = " << outcond << endl ;
#endif
    digraph::vertexSet collapse_search = collapse_vars ;

    collapse_search += input_vars ;
    //    collapse_search += outcond ;
    
    digraph::vertexSet collapse_part = visit_vertices(loop_grt,collapse_search) ;
    collapse_part += collapse_search ;
    // for the collapse part, we also need to add all the rule targets
    // into the graph. Because the collapse part is formed by searching
    // back from the collapse_search, it is possible that we won't include
    // a rule's targets if the rule's targets do not connect to other parts
    // of the collapse graph
    ruleSet collapse_rule_set = extract_rules(collapse_part) ;
    for(ruleSet::const_iterator ri=collapse_rule_set.begin();
        ri!=collapse_rule_set.end();++ri) {
      collapse_part += extract_vars(loop_gr[ri->ident()]) ;
    }

#ifdef MOVE_OUTPUT_TO_COLLAPSE
    digraph::vertexSet output_set ;

    for(ruleSet::const_iterator ri=outrules.begin();ri!=outrules.end();++ri)
      if((ri->sources() - collapse_part) == EMPTY) {
          collapse_search += ri->ident() ;
      } else {
          // Later make this verbose
#ifdef VERBOSE
	debugout << *ri << "put in advance because of " << extract_vars(ri->sources()-collapse_part) << endl ;
#endif
          digraph::vertexSet tmp_search  ;
          tmp_search += ri->ident() ;
          output_set += visit_vertices(loop_grt,tmp_search)+tmp_search ;
      }
    //    output_set -= visit_vertices(loop_grt,collapse_search)+collapse_search ;

    //    cerr << "collapse group = " << extract_rules(output_set) << endl ;
#endif
    
    collapse_gr = loop_gr.subgraph(collapse_part) ;
    collapse_gr.remove_dangling_vertices() ;
    
    
    digraph::vertexSet collapse_rulesV = extract_rules(collapse_gr.get_all_vertices()) ;
    digraph::vertexSet advance_subset = loop_gr.get_all_vertices() - collapse_rulesV ;
    advance_gr = loop_gr.subgraph(advance_subset) ;
    advance_gr.remove_dangling_vertices() ;

    all_loop_vars = all_vars ;

  }

  void loop_compiler::accept(visitor& v) {
    v.visit(*this) ;
  }

  void loop_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    if(duplicate_work) {
      for(variableSet::const_iterator vi = advance_vars.begin();
	  vi != advance_vars.end(); vi++) {
	scheds.add_policy(*vi, sched_db::NEVER);
	variable tmp_var = vi->new_offset(vi->get_info().offset - 1);
	scheds.add_policy(tmp_var, sched_db::NEVER);
      }
    
      list<list<variable> >::const_iterator rli ;
      for(rli = rotate_lists.begin();rli!=rotate_lists.end();++rli) {
        variableSet rvar ;
	list<variable>::const_iterator li ;
	for(li=rli->begin();li!=rli->end();++li) {
          rvar += *li ;
	  scheds.add_policy(*li, sched_db::NEVER);
	}
        scheds.set_variable_rotations(rvar) ;
      }
    }

    list<list<variable> >::const_iterator rli ;
    for(rli = rotate_lists.begin();rli!=rotate_lists.end();++rli) {
      variableSet rvar ;
      list<variable>::const_iterator li ;
      for(li=rli->begin();li!=rli->end();++li) {
        rvar += *li ;
      }
      scheds.set_variable_rotations(rvar) ;
    }

    std::vector<rule_compilerP>::iterator i ;
    for(i=collapse_comp.begin();i!=collapse_comp.end();++i)
      (*i)->set_var_existence(facts, scheds) ;
    for(i=advance_comp.begin();i!=advance_comp.end();++i)
      (*i)->set_var_existence(facts, scheds) ;
  }
  
  void loop_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    variableSet var_requests = advance_vars ;
    variableSet::const_iterator vi ;
    
    Loci::fact_db::distribute_infoP d ;
    if(facts.isDistributed()) {
      d = facts.get_distribute_info() ;
    }
    
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      entitySet vexist = scheds.variable_existence(*vi) ;
      if(facts.isDistributed()) 
	vexist &= d->my_entities;
      scheds.variable_request(*vi,vexist) ;
    }
    
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=advance_comp.rbegin();ri!=advance_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;

    for(ri=collapse_comp.rbegin();ri!=collapse_comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;

    // spread requests over all time iterations
    list<list<variable> >::const_iterator rli ;
    for(rli = rotate_lists.begin();rli!=rotate_lists.end();++rli) {
      list<variable>::const_iterator li ;
      entitySet tot_request ;
      for(li=rli->begin();li!=rli->end();++li)
        tot_request += scheds.get_variable_requests(*li) ;

      for(li=rli->begin();li!=rli->end();++li)
        scheds.variable_request(*li,tot_request) ;
    }


    // then we need to remove all the non-store variables from
    // var_requests to prevent any parallel communications
    variableSet non_stores ;
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      storeRepP sp = facts.get_variable(*vi) ;
      if(!isSTORE(sp))
        non_stores += *vi ;
    }
    var_requests -= non_stores ;
    

    if(facts.isDistributed()) {
      // Communication of n+1 variables should match n=0 variables,
      // so change variable set accordingly
      variableSet adjust_var_requests ;
      for(variableSet::const_iterator vi=var_requests.begin();
          vi!= var_requests.end();++vi) {
        if(vi->time() == tlevel) {
          variable::info vinfo = vi->get_info() ;
          vinfo.assign=true ;
          vinfo.offset=0 ;
	  vinfo.priority = std::vector<std::string>() ;
          adjust_var_requests += variable(vinfo) ;
        }
      }

      list<comm_info> advance_variables_barrier = barrier_process_rule_requests(adjust_var_requests, facts, scheds);
      advance_variables_barrier = sort_comm(advance_variables_barrier, facts);
      scheds.update_comm_info_list(advance_variables_barrier, sched_db::LOOP_ADVANCE_LIST);
      
    }
  }
  
  executeP loop_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    
    CPTR<execute_list> col = new execute_list ;
    
    std::vector<rule_compilerP>::iterator i ;
    for(i=collapse_comp.begin();i!=collapse_comp.end();++i) {
      col->append_list((*i)->create_execution_schedule(facts, scheds)) ;
    }

    CPTR<execute_list> adv = new execute_list ;
    
    for(i=advance_comp.begin();i!=advance_comp.end();++i) {
      adv->append_list((*i)->create_execution_schedule(facts, scheds)) ;
    }

    if(facts.isDistributed()) {
      // Communication of n+1 variables should match n=0 variables,
      // so change variable set accordingly
      variableSet adjust_advance_vars ;
      for(variableSet::const_iterator vi=advance_vars.begin();
          vi!= advance_vars.end();++vi) {
        if(vi->time() == tlevel) {
          variable::info vinfo = vi->get_info() ;
          vinfo.assign=true ;
          vinfo.offset=0 ;
          adjust_advance_vars += variable(vinfo) ;
        }
      }
      std::list<comm_info> advance_variables_barrier = scheds.get_comm_info_list(adjust_advance_vars, facts, sched_db::LOOP_ADVANCE_LIST);
      
      // Now change variable names back before creating communication routine
      // in schedule
      std::list<comm_info>::iterator li ;
      for(li=advance_variables_barrier.begin();li!=advance_variables_barrier.end();++li) {
        
        variable::info vinfo = li->v.get_info() ;
        vinfo.assign=false ;
        vinfo.offset=1 ;
        li->v = variable(vinfo) ;
      }
      execute_comm2::inc_comm_step() ;
      //executeP exec_comm =
      //new execute_comm(advance_variables_barrier, facts);
      executeP exec_comm2 =
        new execute_comm2(advance_variables_barrier, facts);
      adv->append_list(exec_comm2);
      //adv->append_list(exec_comm);
    }
    
    executeP execute = new execute_loop(cond_var,executeP(col),executeP(adv),tlevel,rotate_lists) ;
    return execute;
  }

}
