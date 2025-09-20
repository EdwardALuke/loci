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

#include <Tools/stream.h>

using std::vector ;

#include <sys/stat.h>

#include "dist_tools.h"
#include <parameter.h>

#include <map>
using std::map ;

#include "visit_tools.h"

namespace Loci {
  
  dynamic_compiler::dynamic_compiler(rulecomp_map &rule_process,
                                     const digraph& g, int id):cid(id) {
    ////////////////////
    // store the dag structure and the relevant rulecompiler map
    gr = g ;
    ruleSet allrules = extract_rules(gr.get_all_vertices()) ;
    
    ruleSet::const_iterator ri ;
    for(ri=allrules.begin();ri!=allrules.end();++ri) {
      rulecomp_map::const_iterator rmi ;
      rmi = rule_process.find(*ri) ;
      FATAL(rmi == rule_process.end()) ;
      rule_compiler_map[*ri] = rmi->second ;
    }
    ////////////////////
  }

  void dynamic_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    std::vector<rule_compilerP>::iterator i ;
    for(i=comp.begin();i!=comp.end();++i)
      (*i)->set_var_existence(facts, scheds) ;
  }
  
  void dynamic_compiler::process_var_requests(fact_db &facts,
                                              sched_db &scheds) {
    std::vector<rule_compilerP>::reverse_iterator ri ;
    for(ri=comp.rbegin();ri!=comp.rend();++ri)
      (*ri)->process_var_requests(facts, scheds) ;
  }
  
  executeP dynamic_compiler::create_execution_schedule(fact_db &facts,
                                                       sched_db &scheds) {
    CPTR<execute_list> elp = new execute_list ;
    std::vector<rule_compilerP>::iterator i ;
    for(i=comp.begin();i!=comp.end();++i) {
      elp->append_list((*i)->create_execution_schedule(facts, scheds)) ;
    }

    return executeP(elp) ;
  }

  void dynamic_compiler::collect_reduce_info() {
    digraph grt = gr.transpose() ;
    variableSet allvars = extract_vars(gr.get_all_vertices()) ;
    map<variable,ruleSet> info ;
    
    for(variableSet::const_iterator vi=allvars.begin();
        vi!=allvars.end();++vi) {
      ruleSet pre_rules = extract_rules(grt[vi->ident()]) ;
      ruleSet use_rules ;
      bool reduction = false ;
      bool unit_rule_exists = false ;
      for(ruleSet::const_iterator ri=pre_rules.begin();
          ri!=pre_rules.end();++ri)
        if(ri->get_info().rule_class != rule::INTERNAL) {
          use_rules += *ri ;
          rule_implP rimp = ri->get_rule_implP() ;
          if(rimp->get_rule_class() == rule_impl::UNIT ||
             rimp->get_rule_class() == rule_impl::APPLY)
            reduction = true ;
          if(rimp->get_rule_class() == rule_impl::UNIT)
            unit_rule_exists = true ;          
        }
      if(reduction) {
        info[*vi] = use_rules ;
        if(reduction && unit_rule_exists)
          all_reduce_vars += *vi ;
      }
    }
    if(all_reduce_vars != EMPTY) {
      map<variable,ruleSet>::const_iterator xi ;
      for(xi=info.begin();xi!=info.end();++xi) {
        rule unit_rule ;
        ruleSet apply_rules ;
        CPTR<joiner> join_op = CPTR<joiner>(0) ;
        ruleSet::const_iterator ri ;
        for(ri=xi->second.begin();ri!=xi->second.end();++ri) {
          if(ri->get_rule_implP()->get_rule_class() == rule_impl::UNIT)
            unit_rule = *ri ;
          else if(ri->get_rule_implP()->get_rule_class() == rule_impl::APPLY){
            apply_rules += *ri ;
            join_op = ri->get_rule_implP()->get_joiner() ;
#ifdef TYPEINFO_CHECK
            else if(typeid(*join_op) !=
                    typeid(*(ri->get_rule_implP()->get_joiner()))) {
              cerr << "Warning:  Not all apply rules for variable "
                   << xi->first << " have identical join operations!"
                   << endl ;
            }
#endif
          } else {
            cerr << "Warning: reduction variable " << xi->first
                 << " has a non-reduction rule contributing\
 to its computation,"
                 << endl << "offending rule is " << *ri << endl ;
          }
        }
        if(join_op == 0) {
          cerr << "unable to find any apply rules to complete\
 the reduction defined by rule"
               << endl
               << unit_rule << endl ;
        }
        FATAL(join_op == 0) ;
        reduce_info[xi->first] = make_pair(unit_rule,join_op) ;
      }
    }
  }

  void dynamic_compiler::schedule() {
    digraph grt = gr.transpose() ;
    
    // First schedule any vertices that have no edges leading into them
    //and have not been scheduled previously (in start vertices)
    digraph::vertexSet working = gr.get_source_vertices() -
      gr.get_target_vertices() ;
    if(working != EMPTY) 
      sched.push_back(working) ;
    
    // visited vertices are all vertices that have already been scheduled
    digraph::vertexSet visited_vertices = working ;
    // In the beginning our working set are all scheduled vertices
    working = visited_vertices ;
    while(working != EMPTY) {
      // While we have vertices to work on, compute additional vertices that
      // can be scheduled
      digraph::vertexSet new_vertices ;
      digraph::vertexSet::const_iterator ni ;
      // loop over working set and create a list of candidate vertices
      for(ni=working.begin();ni != working.end(); ++ni) 
        new_vertices += gr[*ni] ;

      // If a vertex has already been scheduled it can't be scheduled again,
      // so remove visited vertices
      new_vertices = new_vertices - visited_vertices    ;
      working = new_vertices ;
      new_vertices = EMPTY ;
      // Find any vertex from this working set that has had all
      // vertices leading to it scheduled
      for(ni=working.begin();ni != working.end(); ++ni) 
        if((grt[*ni] & visited_vertices) == grt[*ni])
          new_vertices += *ni ;
      working = new_vertices ;
      // add these new vertices to the schedule
      if(new_vertices != EMPTY)
        sched.push_back(new_vertices) ;
      // update visited vertices set to include scheduled vertices
      visited_vertices += new_vertices ;
    }
  }

  void dynamic_compiler::compile() {
    digraph grt = gr.transpose() ;
    if(sched.size() == 0)
      return ;
    for(size_t i=0;i<sched.size();++i) {
      variableSet vars = extract_vars(sched[i]) ;
      ruleSet rules = extract_rules(sched[i]) ;
      if(rules == EMPTY && i+1<sched.size()) {
        ++i ;
        vars += extract_vars(sched[i]) ;
        rules = extract_rules(sched[i]) ;
      }

      variableSet barrier_vars, reduce_vars,singleton_vars,all_vars ;
      variableSet::const_iterator vi ;
      
      for(vi=vars.begin();vi!=vars.end();++vi) {
        ruleSet var_rules = extract_rules(grt[(*vi).ident()]) ;
        ruleSet::const_iterator ri ;
        ruleSet use_rules ;
        bool reduction = false ;
        bool pointwise = false ;
        bool singleton = false ;
        bool recursive = false ;
        bool unit_rule_exists = false ;
        bool priority_rule = false ;

        for(ri=var_rules.begin();ri!=var_rules.end();++ri) {
          if(!is_virtual_rule(*ri) ||
             (ri->get_info().rule_class == rule::INTERNAL &&
              ri->get_info().qualifier() == "priority")) {
            use_rules += *ri ;

            // Check for a priority rule
            if(ri->get_info().rule_class == rule::INTERNAL &&
               ri->get_info().qualifier() == "priority")
              priority_rule = true ;

            if(ri->get_info().rule_class == rule::INTERNAL)
              if(is_chomp_node(ri)) {
                // we need to actually look into the chomping
                // graph to find out the rules that generate
                // this variable and dertermin the types of
                // rules
                rulecomp_map::const_iterator rmi ;
                rmi = rule_compiler_map.find(*ri) ;
                FATAL(rmi == rule_compiler_map.end()) ;
                rule_compilerP rcp = rmi->second ;
                chomp_compiler* chc =
                  dynamic_cast<chomp_compiler*>(&(*(rcp))) ;
                FATAL(chc == 0) ;
                digraph cgrt = (chc->chomp_graph).transpose() ;
                ruleSet chomp_var_rules = extract_rules(cgrt[vi->ident()]) ;
                for(ruleSet::const_iterator cri=chomp_var_rules.begin() ;
                    cri!=chomp_var_rules.end();++cri) {
                  
                  rule_implP crimp = cri->get_rule_implP() ;
                  if(crimp->get_rule_class() == rule_impl::POINTWISE)
                    pointwise = true ;
                  
                  if(crimp->get_rule_class() == rule_impl::UNIT ||
                     crimp->get_rule_class() == rule_impl::APPLY)
                    reduction = true ;
                  
                  if(crimp->get_rule_class() == rule_impl::UNIT)
                    unit_rule_exists = true ;
                  
                  if(crimp->get_rule_class() == rule_impl::SINGLETON)
                    singleton = true ;
                }
                continue ;
              }
          
            rule_implP rimp = ri->get_rule_implP() ;
            if(rimp->get_rule_class() == rule_impl::POINTWISE)
              pointwise = true ;
            
            if(rimp->get_rule_class() == rule_impl::UNIT ||
               rimp->get_rule_class() == rule_impl::APPLY)
              reduction = true ;

            if(rimp->get_rule_class() == rule_impl::UNIT)
              unit_rule_exists = true ;
            
            if(rimp->get_rule_class() == rule_impl::SINGLETON)
              singleton = true ;
          } else {
            if((ri->sources() & ri->targets()) != EMPTY)
              recursive = true ;
          }
        }
        
        WARN((reduction && pointwise) || (pointwise && singleton) ||
             (reduction && singleton)) ;

        if((use_rules != EMPTY)) {
          if((priority_rule || pointwise) && !recursive && (vi->get_info().name != "OUTPUT")){
            // Don't use the priority variables for variable barriers
            if(vi->get_info().priority.size() == 0)
              barrier_vars += *vi ;
          }
          if((priority_rule || pointwise) && recursive && (vi->get_info().name != "OUTPUT")) {
            if(vi->get_info().priority.size() == 0)
              all_vars += *vi ;
          }

          if(reduction && unit_rule_exists)
            reduce_vars += *vi ;
          
          if(singleton) {
            singleton_vars += *vi ;
          }
        }

      }

      all_vars += barrier_vars ;
      if(barrier_vars != EMPTY)
        comp.push_back(new barrier_compiler(barrier_vars)) ;
      
      all_vars += singleton_vars ;

      if(singleton_vars != EMPTY)
        comp.push_back(new singleton_var_compiler(singleton_vars)) ;

      all_vars += reduce_vars;

      vector<CPTR<joiner> > join_op_vector ;
      vector<rule> unit_rule_vector ;
      vector<variable> reduce_var_vector ;
      
      for(variableSet::const_iterator rdvi=reduce_vars.begin();
          rdvi!=reduce_vars.end();++rdvi) {
        map<variable,pair<rule,CPTR<joiner> > >::const_iterator xi ;
        xi = reduce_info.find(*rdvi) ;
        FATAL(xi == reduce_info.end()) ;
        rule unit_rule = (xi->second).first ;
        CPTR<joiner> join_op = (xi->second).second ;

        if(join_op != 0) {
          storeRepP sp = join_op->getTargetRep() ;
          if(isPARAMETER(sp)) {
            reduce_var_vector.push_back(xi->first) ;
            unit_rule_vector.push_back(unit_rule) ;
            join_op_vector.push_back(join_op) ;
          }
          else if (isBLACKBOX(sp)) {
            cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
          }else {
            comp.push_back(new reduce_store_compiler(xi->first,unit_rule,
                                                         join_op)) ;
          }
        }
      }
      if(reduce_var_vector.size() != 0)  
        comp.push_back(new reduce_param_compiler(reduce_var_vector, unit_rule_vector, join_op_vector));
      
      
      if(rules != EMPTY) {
        ruleSet::const_iterator ri ;
        for(ri=rules.begin();ri!=rules.end();++ri) {
          rulecomp_map::const_iterator rmi ;
          rmi = rule_compiler_map.find(*ri) ;
          FATAL(rmi == rule_compiler_map.end()) ;
          comp.push_back(rmi->second) ;
        }
      }
    }    
  }
  
}

