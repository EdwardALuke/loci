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
#include <ostream>
#include <iostream>
#include "sched_tools.h"
#include "comp_tools.h"
#include "dist_tools.h"
#include "visitor.h"
#include "loci_globs.h"
#include "thread.h"

#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <list>
using std::list ;
#include <set>
using std::set ;
#include <string>
using std::string ;
#include <sstream>
using std::stringstream ;
#include <utility>
using std::pair ;


using std::cout ;
using std::ostream ;
using std::cerr ;
using std::endl ;

//#define HACK ; 

namespace Loci {
  extern void create_digraph_dot_file(const digraph&, const char*) ;
  // Loci option flags
  extern double total_memory_usage ;
  extern bool show_decoration ;
  extern bool use_dynamic_memory ;
  extern bool dmm_no_deallocation ;
  extern bool show_dmm_verbose ;
  extern bool use_chomp ;
  extern bool show_chomp ;
  extern bool chomp_verbose ;
  extern int chomping_size ;
  extern bool profile_memory_usage ;
  extern bool memory_greedy_schedule ;
  extern bool randomized_memory_greedy_schedule;
  // memory profiling counters
  extern double LociAppPeakMemory ;
  extern double LociAppAllocRequestBeanCounting ;
  extern double LociAppFreeRequestBeanCounting ;
  extern double LociAppPeakMemoryBeanCounting ;
  extern double LociAppLargestAlloc ;
  extern variable LociAppLargestAllocVar ;
  extern double LociAppLargestFree ;
  extern variable LociAppLargestFreeVar ;
  extern double LociAppPMTemp ;
  extern double LociInputVarsSize ;
  // threads
  extern bool threading_pointwise;
  extern bool threading_global_reduction;
  extern bool threading_local_reduction;
  extern bool threading_chomping;
  extern int num_threads;
  bool in_internal_query = false ;
  namespace {
    // used to pre-process preallocation memory profiling
    variableSet LociRecurrenceVarsRealloc ;
    variableSet LociPreallocationReallocVars ;
    variableSet LociAppGivenVars ;
    double LociPreallocationReallocSize = 0 ;
    // all chomped variables
    //variableSet all_chomped_vars = variableSet(EMPTY) ;
    variableSet all_chomped_vars ;
  }
  
  class error_compiler : public rule_compiler {
  public:
    error_compiler() {}
    virtual void accept(visitor& v)
    { cerr << "Internal consistency error" << endl ; Loci::Abort();}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds)
    { cerr << "Internal consistency error" << endl ; Loci::Abort();}
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) 
    { cerr << "Internal consistency error" << endl ; Loci::Abort();}
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds)
    { cerr << "Internal consistency error" << endl ; Loci::Abort();
    return executeP(0);}
  } ;

  int get_supernode_num(int rid) {
    rule r(rid) ;
    string rqualifier = r.get_info().qualifier() ;
    string head = rqualifier.substr(0,2) ;
    if(head != "SN") return -1 ;
    
    string number = rqualifier.substr(2,rqualifier.size()-2) ;
    stringstream ss ;
    ss << number ;
    int ret ;
    ss >> ret ;

    return ret ;
  }
 
  graph_compiler::graph_compiler(decomposed_graph &deco,variableSet
				 initial_vars) {
    fact_db_comm = new barrier_compiler(initial_vars) ;
    multiLevelGraph &mlg = deco.mlg ;
    vector<int> levels ;
    baserule = rule(mlg.toplevel) ;
    digraph::vertexSet working ;
    working += mlg.toplevel ;
    digraph::vertexSet::const_iterator vi ;
    while(working != EMPTY) {
      digraph::vertexSet next_level ;
      for(vi=working.begin();vi!=working.end();++vi) {
        levels.push_back(*vi) ;
        next_level += mlg.find(*vi)->graph_v & mlg.subgraphs ;
      }
      working = next_level ;
    }
    for(vector<int>::reverse_iterator ii=levels.rbegin();
        ii!=levels.rend();
        ++ii) {
      multiLevelGraph::subGraph *p = mlg.find(*ii) ;
      ruleSet level_rules = extract_rules(p->graph_v-mlg.subgraphs) ;
      digraph gr = p->gr ;

      // Remove any rules in that are not in graph_v but are in gr
      ruleSet errs = extract_rules(gr.get_all_vertices()-p->graph_v) ;
      gr.remove_vertices(errs) ;
      gr.remove_dangling_vertices() ;

      digraph grt = gr.transpose() ;

      // Collect information on unit rules on this level
      map<rule,rule> apply_to_unit ;
      for(ruleSet::const_iterator ri=level_rules.begin();
          ri!=level_rules.end();
          ++ri) {
        if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
          variable unit_var = *(ri->targets().begin()) ;
          ruleSet apply_rules = extract_rules(grt[unit_var.ident()]) ;
          for(ruleSet::const_iterator rii=apply_rules.begin();
              rii!=apply_rules.end();
              ++rii) {
            apply_to_unit[*rii] = *ri ;
          }
        }
      }

      // Assign compilers to rules on this level
      for(ruleSet::const_iterator ri=level_rules.begin();
          ri!=level_rules.end();
          ++ri) {
        if(ri->type() == rule::INTERNAL) {
          if(ri->get_info().qualifier() == "promote")
            rule_process[*ri] = new promote_compiler(*ri) ;
          else if(ri->get_info().qualifier() == "generalize")
            rule_process[*ri] = new generalize_compiler(*ri) ;
          else if(ri->get_info().qualifier() == "priority")
            rule_process[*ri] = new priority_compiler(*ri) ;
	  else if(ri->get_info().qualifier() == "gpu2cpu") 
	    rule_process[*ri] = new gpu2cpu_compiler(*ri) ;
	  else if(ri->get_info().qualifier() == "cpu2gpu")
	    rule_process[*ri] = new cpu2gpu_compiler(*ri) ;
	  else if(ri->get_info().qualifier() == "map2gpu")
	    rule_process[*ri] = new map2gpu_compiler(*ri) ;
          else if(ri->get_info().qualifier() == "ALLOCATE")
            rule_process[*ri] = new allocate_var_compiler(ri->targets()) ;
          else
            rule_process[*ri] = new error_compiler ;
        } else {
#ifdef DYNAMICSCHEDULING
          if(ri->get_info().rule_impl->get_rule_class()
             == rule_impl::INSERTION) // an insertion rule
            rule_process[*ri] = new insertion_rule_compiler(*ri) ;
          else if(ri->get_info().rule_impl->get_rule_class()
                  == rule_impl::DELETION) // a deletion rule
            rule_process[*ri] = new deletion_rule_compiler(*ri) ;
          else if(ri->get_info().rule_impl->get_rule_class()
                  == rule_impl::ERASE) // a erase rule
            rule_process[*ri] = new erase_rule_compiler(*ri) ;
          else
#endif
	    if(ri->get_info().rule_impl->get_rule_class()
                  == rule_impl::CONSTRAINT_RULE) // a constraint rule
            rule_process[*ri] = new constraint_compiler(*ri) ;
          else if(ri->get_info().rule_impl->get_rule_class()
                  == rule_impl::MAP_RULE) // a map rule (to be revised)
            rule_process[*ri] = new map_compiler(*ri) ;
          else if(ri->get_info().rule_impl->get_rule_class()
                  == rule_impl::BLACKBOX_RULE) // a blackbox rule
            rule_process[*ri] = new blackbox_compiler(*ri) ;
          else if(ri->get_info().rule_impl->get_rule_class()
                  == rule_impl::SUPER_RULE) // a blackbox rule
            rule_process[*ri] = new superRule_compiler(*ri) ;
          else if(ri->get_info().rule_impl->get_rule_class()
                  == rule_impl::APPLY) // an apply rule
            rule_process[*ri] = new apply_compiler(*ri,apply_to_unit[*ri]) ;
          else                  // a normal concrete rule
            rule_process[*ri] = new impl_compiler(*ri) ;
        }
      }

      // We have created compilers for every rule on this level, now
      // we create a compiler for this supernode

      rule snrule = rule(*ii) ;
      // push each snrule into the vector super_rules for later use
      super_rules.push_back(snrule) ;
      int id = get_supernode_num(*ii) ;

      if(deco.loops.inSet(*ii)) {
        // Looping supernode
        rule_process[snrule] = new loop_compiler(rule_process,gr,id) ;
      } else if(deco.recursive.inSet(*ii)) {
        // recursive supernode
        ruleSet recurse_rules = extract_rules(p->graph_v) ;
        if(recurse_rules.size() == 1 &&
           recurse_rules.begin()->get_info().desc.constraints.size() == 0
           && ((*recurse_rules.begin()).get_rule_implP()->is_relaxed()))
          // Single rule recursion
          rule_process[snrule] =
            new impl_recurse_compiler(*(recurse_rules.begin()),id) ;
        else
          // Multi-rule recursion
          rule_process[snrule] =
            new recurse_compiler(rule_process, recurse_rules, id) ;
      } else if(deco.conditional.inSet(*ii)) {
        // conditional supernode
        variable cond_var = *(snrule.get_info().desc.conditionals.begin()) ;
        rule_process[snrule] =
          new conditional_compiler(rule_process,gr,cond_var,id) ;
      } else {
        // DAG supernode
        rule_process[snrule] = new dag_compiler(rule_process,gr,id) ;
      }
    } 
  }

  void graph_compiler::top_down_visit(visitor& v) {
    for(vector<rule>::reverse_iterator ii=super_rules.rbegin();
        ii!=super_rules.rend();++ii) {
      rulecomp_map::iterator ri ;
      ri = rule_process.find(*ii) ;
      FATAL(ri == rule_process.end()) ;
      (ri->second)->accept(v) ;
    }
  }
  
  void graph_compiler::bottom_up_visit(visitor& v) {
    for(vector<rule>::const_iterator ii=super_rules.begin();
        ii!=super_rules.end();++ii) {
      rulecomp_map::iterator ri ;
      ri = rule_process.find(*ii) ;
      FATAL(ri == rule_process.end()) ;
      (ri->second)->accept(v) ;
    }
  }
  
  fact_db *exec_current_fact_db = 0 ;

  //#define COMPILE_PROGRESS
  void graph_compiler::compile(fact_db& facts,sched_db& scheds,
                               const variableSet& given,
                               const variableSet& target) {
#ifdef COMPILE_PROGRESS
    if(Loci::MPI_rank==0)
      cerr << "[Graph Compile Phase] Start!" << endl ;
#endif
    /***********************************
    1. unordered visitor phase (1..n)
         top -> down
         bottom -> up
       graph editing phase (1..n)

    2. ordering phase (1..n)

    3. assembly phase.
    ************************************/
    // timing variables
    double dst1=0,det1=0,dst2=0,det2=0,cst=0,cet=0,schedst=0,schedet=0 ;
    dst1 = MPI_Wtime() ;

    // get all the recurrence information
    // i.e. the generalize, promote, priority and rename info
    recurInfoVisitor recv ;
    top_down_visit(recv) ;

    // must do this visitation
    // the loop rotate_lists is computed here
    rotateListVisitor rotlv(scheds,
                            recv.get_rename_s2t(),
                            recv.get_rename_t2s()) ;
    top_down_visit(rotlv) ;

    // set the reallocated recurrence vars for
    // memory profiling preallocation
    if(profile_memory_usage) {
      variableSet sources = variableSet(recv.get_recur_source_vars()
                                        - recv.get_recur_target_vars()) ;
      map<variable,variableSet> s2t_table = recv.get_recur_vars_s2t() ;
      LociRecurrenceVarsRealloc +=
        get_recur_target_for_vars(sources,s2t_table) ;
      LociAppGivenVars = given ;

      for(variableSet::const_iterator vi=given.begin();
          vi!=given.end();++vi) {
        storeRepP srp = facts.get_variable(*vi) ;
        LociInputVarsSize += srp->pack_size(srp->domain()) ;
      }
    }
    
    // get the input variables
    // get all the variables in the multilevel graph
    getAllVarVisitor allvarV ;
    top_down_visit(allvarV) ;
    variableSet input = variableSet(given & allvarV.get_all_vars_ingraph()) ;
    variableSet all_vars = allvarV.get_all_vars_ingraph();
    // check preconditions of recurrence
    check_recur_precondition(recv,input) ;

    // get the unit apply info (reduce info)
    unitApplyMapVisitor reduceV ;
    top_down_visit(reduceV) ;

#ifdef DYNAMICSCHEDULING
    // replace all the dynamic rules with corresponding compilers
    DynamicKeyspaceVisitor dkv(facts,scheds,
                               reduceV.get_apply2unit(),
                               rotlv.get_overlap_rotvars()) ;
    top_down_visit(dkv) ;
#endif
    
    // insert rule that redistribute dynamic keyspaces

#ifdef COMPILE_PROGRESS
    if(Loci::MPI_rank==0)
      cerr << "[Graph Compile Phase] Passed Information Collection!" << endl ;
#endif

    det1 = MPI_Wtime() ;

    if(use_dynamic_memory) {
      if(Loci::MPI_rank == 0)
        if(!in_internal_query)
          cout << "USING DYNAMIC MEMORY MANAGEMENT" << endl ;
    }

    // chomping searching and graph editing
    if(use_chomp) {
      if(Loci::MPI_rank == 0)
        if(!in_internal_query)
          cout << "USING CHOMPING"
               << " (chomping size: "
               << chomping_size << "KB)"
               << endl ;
      
      //	  timer_token chomping_subgraph_searching_timer = new timer_token;
      cst = MPI_Wtime() ;

      chompPPVisitor cppv(facts,
                          rotlv.get_rotate_vars_table(),
                          rotlv.get_loop_shared_table(),
                          variableSet(recv.get_rename_source_vars() +
                                      recv.get_rename_target_vars())
                          ) ;
      top_down_visit(cppv) ;

      chompRuleVisitor crv(cppv.get_good_vars(),
                           cppv.get_bad_vars(),
                           reduceV.get_apply2unit()) ;
      top_down_visit(crv) ;
      cet = MPI_Wtime() ;
		
      if(show_chomp)
        crv.visualize(cout) ;
      if(chomp_verbose)
        crv.summary(cout) ;
      
      dagCheckVisitor dagcV1(true) ;
      top_down_visit(dagcV1) ;

      // we set all chomped variables here
      all_chomped_vars = crv.get_all_chomped_vars() ;
    }

    // dynamic memory management graph decoration
    if(use_dynamic_memory) {
      dst2 = MPI_Wtime() ;

      // get inter/intra supernode information
      snInfoVisitor snv ;
      top_down_visit(snv) ;
      
      // get all untyped variables
      unTypedVarVisitor untypevarV(recv.get_recur_vars_s2t(),
                                   recv.get_recur_vars_t2s(),
                                   input) ;
      //top_down_visit(untypevarV) ;

      // get the representitive variable for each rename target cluster
      variableSet only_targets = recv.get_rename_target_vars() ;
      only_targets -= recv.get_rename_source_vars() ;

      variableSet cluster_remaining ;
      cluster_remaining = pick_rename_target(recv.get_rename_s2t(),
                                             recv.get_rename_t2s(),
                                             only_targets) ;

      // get representitive for each promoted variable cluster
      promotePPVisitor pppV(recv.get_promote_t2s(),
                            recv.get_promote_s2t(),
                            recv.get_promote_source_vars(),
                            recv.get_promote_target_vars(),
                            snv.get_graph_sn(),
                            snv.get_cond_sn(),input) ;
      top_down_visit(pppV) ;

      // reserved variableSet for the deleteInfoVisitor
      variableSet delInfoV_reserved ;

      variableSet promoted_rep =
        variableSet(pppV.get_rep() - recv.get_rename_source_vars()) ;

      delInfoV_reserved += recv.get_recur_source_vars() ;
      delInfoV_reserved += pppV.get_remaining() ;
      delInfoV_reserved -= promoted_rep ;
      delInfoV_reserved += input ;
      delInfoV_reserved += get_recur_target_for_vars(input,
                                                     recv.get_recur_vars_s2t()
                                                     ) ;
      delInfoV_reserved += target ;
      // any recurrence variables connected to target should not
      // be deallocated
      delInfoV_reserved += recv.get_reachable(target) ;

      delInfoV_reserved += untypevarV.get_untyped_vars() ;
      delInfoV_reserved += cluster_remaining ;

      // compute how to do allocation
      allocInfoVisitor aiv(snv.get_graph_sn(),
                           recv.get_recur_vars_s2t(),
                           recv.get_recur_vars_t2s(),
                           variableSet(recv.get_recur_source_vars()+
                                       recv.get_recur_target_vars()),
                           snv.get_loop_sn(),
                           rotlv.get_rotate_vars_table(),
                           rotlv.get_loop_shared_table(),
                           // we also don't want to allocate
                           // all the dynamic rule targets because
                           // they are handled dynamically
#ifdef DYNAMICSCHEDULING
                           variableSet(untypevarV.get_untyped_vars() +
                                       dkv.get_dynamic_targets())) ;
#else
                           variableSet(untypevarV.get_untyped_vars())) ;
#endif
      top_down_visit(aiv) ;

      // compute how to do deletion
      deleteInfoVisitor div(get_loop_alloc_table(aiv.get_alloc_table(),
                                                 snv.get_subnode_table(),
                                                 snv.get_loop_sn(),
                                                 snv.get_graph_sn(),
                                                 recv.get_recur_vars_s2t()),
                            recv.get_recur_vars_t2s(),
                            recv.get_recur_vars_s2t(),
                            snv.get_graph_sn(),
                            get_parentnode_table(snv.get_subnode_table()),
                            snv.get_loop_col_table(),
                            snv.get_loop_sn(),
                            snv.get_cond_sn(),
                            rotlv.get_rotate_vars_table(),
                            rotlv.get_loop_shared_table(),
                            promoted_rep,delInfoV_reserved
                            ) ;
      if(!dmm_no_deallocation) {
        top_down_visit(div) ;
      }
      
      // decorate the graph to insert allocation rules
      allocGraphVisitor agv(aiv.get_alloc_table(),
                            snv.get_loop_sn(),
                            rotlv.get_rotate_vars_table(),
                            recv.get_priority_s2t(),
                            recv.get_priority_source_vars()
                            ) ;
      top_down_visit(agv) ;

      if(profile_memory_usage) {
        // decorate the graph to include
        // the memory profiling alloc compiler
        memProfileAllocDecoVisitor mpadv(aiv.get_alloc_table(),
                                         rotlv.get_rotate_vars_table()
                                         ) ;
        top_down_visit(mpadv) ;
      }
      
      // decorate the graph to insert deletion rules
      deleteGraphVisitor dgv(div.get_delete_table(),
                             div.get_recur_source_other_rules()) ;
      if(!dmm_no_deallocation)
        top_down_visit(dgv) ;

      det2 = MPI_Wtime() ;
	  
#ifdef COMPILE_PROGRESS
    if(Loci::MPI_rank==0)
      cerr << "[Graph Compile Phase] Passed Memory "
           << "Management Decoration!" << endl ;
#endif
      
      // check if the decorated graphs are acyclic
      dagCheckVisitor dagcV(true) ;
      top_down_visit(dagcV) ;

#ifdef COMPILE_PROGRESS
    if(Loci::MPI_rank==0)
      cerr << "[Graph Compile Phase] Passed Graph Cycle Check!" << endl ;
#endif

      if(show_dmm_verbose && !in_internal_query) {
        // the output stream
        ostream& os = cout ;
        //print out the loop rotate list variables
        os << endl ;
        os << "---------loop rotate list variables-----------" << endl ;
        os << rotlv.get_rotate_vars_table() ;
        //print out the loop shared variables
        os << endl ;
        os << "---------loop shared variables-----------" << endl ;
        os << rotlv.get_loop_shared_table() ;
        os << endl ;
        //print out the untyped variables
        os << "----------untyped variables--------------" << endl ;
        os << untypevarV.get_untyped_vars() ;
        os << endl ;

        os << endl ;
        // get all typed variables in the fact database
        variableSet vars = facts.get_typed_variables() ;
        variableSet allgraphv = allvarV.get_all_vars_ingraph() ;
        // the difference
        variableSet typed_not_in_graph = variableSet(vars - allgraphv) ;
        if(typed_not_in_graph != EMPTY) {
          os << "These are typed variables in fact database, "
             << "but are not in the multilevel graph: "
             << typed_not_in_graph << endl << endl ;
        }
        
        variableSet allocated_vars ;
        map<int,variableSet> alloc_table = aiv.get_alloc_table() ;
        map<int,variableSet>::const_iterator miter ;
        for(miter=alloc_table.begin();miter!=alloc_table.end();++miter)
          allocated_vars += miter->second ;
        vars -= allocated_vars ;
        vars -= recv.get_recur_target_vars() ;
        // input variables
        rulecomp_map::iterator ri ;
        ri = rule_process.find(*super_rules.rbegin()) ;
        FATAL(ri == rule_process.end()) ;
        dag_compiler* dc = dynamic_cast<dag_compiler*>(&(*(ri->second))) ;
        digraph gr = dc->dag_gr ;
        digraph grt = gr.transpose() ;
        
        digraph::vertexSet input_vars ;
        digraph::vertexSet all_vertices = grt.get_all_vertices() ;  
        for(digraph::vertexSet::const_iterator vi=all_vertices.begin();
            vi!=all_vertices.end();++vi) {
          if(grt[*vi] == EMPTY)
            input_vars += *vi ;
        }
        vars -= extract_vars(input_vars) ;
        vars &= allgraphv ;
        
        if(vars != EMPTY) {
          os << "These are typed variables in the multilevel "
             << "graph, but are not scheduled to be allocated: "
             << vars << endl << endl ;
        }
        
        // report some allocation and deletion info
        allocDeleteStat adstat(aiv.get_alloc_table(),
                               div.get_delete_table(),
                               recv.get_recur_vars_t2s(),
                               recv.get_recur_vars_s2t(),
                               recv.get_recur_source_vars(),
                               recv.get_recur_target_vars()
                               ) ;
        adstat.report(os) ;
        os << endl ;

        allocDelNumReportVisitor adnrv(cout) ;
        top_down_visit(adnrv) ;
        //os<<recv.get_recur_vars_s2t()<<endl ;
        os<<recv.get_rename_s2t()<<endl ;
        os << endl ;
      } // end of if(show_dmm_verbose)
      
    } // end of if(use_dynamic_memory)
    
    // visualize each decorated graph
    if(show_decoration) {
      graphVisualizeVisitor visV ;
      top_down_visit(visV) ;
    }
    
    if(use_chomp) {
      // compile all the chomp_compilers.
      // Note that all the chomping compilers
      // are added into the multilevel graph
      // later after we built the multilevel
      // graph structure. Therefore all the
      // chomping compilers are NOT in the
      // graph structure that we built. So
      // we need to schedule and compile
      // all the chomping compilers separately.
      // We do this here. The scheduling
      // and compilation of all other compilers
      // are done at below by graphSchedulerVisitor
      // and the assembleVisitor.
      compChompVisitor compchompv(reduceV.get_reduceInfo()) ;
      top_down_visit(compchompv) ;
#ifdef COMPILE_PROGRESS
    if(Loci::MPI_rank==0)
      cerr << "[Graph Compile Phase] Passed Chomping Compilation!" << endl ;
#endif    
    }

    //orderVisitor ov ;    
    //bottom_up_visit(ov) ;
	
    schedst = MPI_Wtime() ;
   
    if(randomized_memory_greedy_schedule && (!in_internal_query)){ 
      if(Loci::MPI_rank == 0)
        if(!in_internal_query)
          cout << "graph scheduling... (randomized memory greedy)" << endl ;
      
      sched_db local_scheds =  scheds;
      fact_db local_facts = facts ;
      
      compGreedyPrio cgp2 ;
      graphSchedulerVisitor gsv2(cgp2) ;
      top_down_visit(gsv2) ;
#ifdef DYNAMICSCHEDULING
      assembleVisitor av2(local_facts, local_scheds,
                          reduceV.get_all_reduce_vars(),
                          reduceV.get_reduceInfo(),
                          dkv.get_dynamic_targets(),
                          dkv.get_self_clone(),
                          dkv.get_shadow_clone(),
                          dkv.get_drule_ctrl()) ;
      bottom_up_visit(av2) ;
#else
      assembleVisitor av2(reduceV.get_all_reduce_vars(),
                          reduceV.get_reduceInfo()) ;
      bottom_up_visit(av2) ;
#endif
      existential_analysis(local_facts, local_scheds);
    
      map<variable, double> var_info;
      int num_vars = all_vars.size();
      double* local_sizes = new double[num_vars];
      double* global_sizes = new double[num_vars];
      
      //collect the info
      int id = 0;
      for(variableSet::const_iterator
            mi=all_vars.begin();mi!=all_vars.end();++mi) {
        entitySet eset = local_scheds.get_variable_requests(*mi);
        storeRepP st = facts.get_variable(*mi);
        local_sizes[id++] = st->estimated_pack_size(eset);
      }
      
      for(id = 0; id < num_vars; id++)global_sizes[id] = 0;
       MPI_Allreduce(local_sizes,
                     global_sizes,
                     num_vars, MPI_DOUBLE,
                     MPI_MAX,MPI_COMM_WORLD) ;
       id = 0;
       for(variableSet::const_iterator
             mi=all_vars.begin();mi!=all_vars.end();++mi) {
         var_info[*mi]= global_sizes[id++];
       }
       delete [] local_sizes;
       delete [] global_sizes;
       
       SchedClearVisitor scv;
       top_down_visit(scv) ;
       MemGreedyScheduler mgs(facts,
                              scheds,
                              recv.get_recur_vars_s2t(),
                              recv.get_recur_vars_t2s(),
                              var_info);
      
       top_down_visit(mgs) ;
      
    }else if(!memory_greedy_schedule) {
      if(Loci::MPI_rank == 0)
        if(!in_internal_query)
          cout << "graph scheduling... (computation greedy)" << endl ;
      
      compGreedyPrio cgp ;
      graphSchedulerVisitor gsv(cgp) ;
      top_down_visit(gsv) ;
      
      
    } else {
      if(Loci::MPI_rank == 0)
        if(!in_internal_query)
          cout << "graph scheduling... (memory greedy)" << endl ;
      
      memGreedyPrio mgp(facts) ;
      graphSchedulerVisitor gsv(mgp) ;
      top_down_visit(gsv) ;
    }
#ifdef COMPILE_PROGRESS
    if(Loci::MPI_rank==0)
      cerr << "[Graph Compile Phase] Passed Graph Scheduling!" << endl ;
#endif

    schedet = MPI_Wtime() ;
    
#ifdef DYNAMICSCHEDULING
    assembleVisitor av(facts, scheds,
                       reduceV.get_all_reduce_vars(),
                       reduceV.get_reduceInfo(),
                       dkv.get_dynamic_targets(),
                       dkv.get_self_clone(),
                       dkv.get_shadow_clone(),
                       dkv.get_drule_ctrl()) ;
    bottom_up_visit(av) ;
#else
  assembleVisitor av(reduceV.get_all_reduce_vars(),
                     reduceV.get_reduceInfo()) ;
  bottom_up_visit(av) ;
#endif



#ifdef COMPILE_PROGRESS
    if(Loci::MPI_rank==0)
      cerr << "[Graph Compile Phase] Passed Schedule Assembly!" << endl ;
#endif

    // setting this external pointer
    exec_current_fact_db = &facts ;

    if(!in_internal_query)
      if(use_dynamic_memory)
        Loci::debugout << "Time taken for dmm graph decoration = "
                       << (det1-dst1) + (det2-dst2) << " seconds " << endl ;
    if(!in_internal_query)
      if(use_chomp)
        Loci::debugout << "Time taken for chomping subgraph searching = "
                       << cet-cst << " seconds " << endl ;
    if(!in_internal_query)
      Loci::debugout << "Time taken for graph scheduling = "
                     << schedet-schedst << " sceonds " << endl ;
#ifdef COMPILE_PROGRESS
    if(Loci::MPI_rank==0)
      cerr << "[Graph Compile Phase] Graph Compile Phase End!" << endl ;
#endif
  }
  
  void graph_compiler::existential_analysis(fact_db &facts, sched_db &scheds) {
    fact_db_comm->set_var_existence(facts, scheds) ;
    (rule_process[baserule])->set_var_existence(facts, scheds) ;
    variableSet var_requests = baserule.targets() ;
    variableSet::const_iterator vi ;
    entitySet my_entities = ~EMPTY ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      my_entities = d->my_entities ;
    }
    for(vi=var_requests.begin();vi!=var_requests.end();++vi) {
      entitySet vexist = scheds.variable_existence(*vi)&my_entities ;
      scheds.variable_request(*vi,vexist) ;
    }
    scheds.add_possible_duplicate_vars(var_requests);
    (rule_process[baserule])->process_var_requests(facts, scheds) ;
    fact_db_comm->process_var_requests(facts, scheds) ;
  }
  
  class allocate_all_vars : public execute_modules {
    variableSet vars ;
    bool is_alloc_all ;
    std::map<variable,entitySet> v_requests, v_existence ;
  public:
    allocate_all_vars() { }
    allocate_all_vars(fact_db &facts, sched_db &scheds,
                      const variableSet& alloc,
                      bool is_alloc_all) ;
    void fill_in_requests(fact_db &facts, sched_db &scheds) ;
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "allocate_all_vars";};
    virtual void dataCollate(collectData &data_collector) const {} ;
  } ;
  
  
  allocate_all_vars::allocate_all_vars(fact_db &facts, sched_db &scheds,
                                       const variableSet& alloc,
                                       bool is_alloc_all)
    :vars(alloc),is_alloc_all(is_alloc_all) {
    fill_in_requests(facts, scheds) ;
  }

  void allocate_all_vars::fill_in_requests(fact_db &facts, sched_db &scheds) {
    //    variableSet vars = facts.get_typed_variables() ;
    variableSet::const_iterator vi,vii ;
    
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      variableSet aliases = scheds.get_aliases(*vi) ;
      entitySet requests, existence ;
      for(vii=aliases.begin();vii!=aliases.end();++vii) {
	existence += scheds.variable_existence(*vii) ;
	requests += scheds.get_variable_requests(*vii) ;
      }
      v_requests[*vi] = requests ;
      v_existence[*vi] = existence ;
    }
  }
  
  void allocate_all_vars::execute(fact_db &facts, sched_db& scheds) {
    //exec_current_fact_db = &facts ;
    //variableSet vars = facts.get_typed_variables() ;
    variableSet::const_iterator vi ;  
    double total_size = 0 ;
    entitySet dom, total, unused ;
    double total_wasted = 0 ;
    
    //#define DIAGNOSTICS
#ifdef DIAGNOSTICS
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      if(isSTORE(srp)) {
	dom = v_requests[*vi] ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	total_size += srp->pack_size(dom) ; 
	total_wasted += srp->pack_size(unused) ;
      }
    }
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
      if(isSTORE(srp)) {
	dom = v_requests[*vi] ;
	double size, wasted_space ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	size = srp->pack_size(dom) ;
	wasted_space = srp->pack_size(unused) ;
	Loci::debugout << " ****************************************************" << endl ;
	Loci::debugout << " Total_size = " << total_size << endl ;
	Loci::debugout << "Variable = "  << *vi << endl ;
	Loci::debugout << "Domain = " << dom << endl ;
	Loci::debugout << "Size allocated = " << size << endl ; 
	if(facts.isDistributed() )  {
	  Loci::fact_db::distribute_infoP d ;
	  d   = facts.get_distribute_info() ;
	  entitySet my_entities = d->my_entities ; 
	  entitySet clone = dom - my_entities ;
	  double clone_size = srp->pack_size(clone) ;
	  Loci::debugout << "----------------------------------------------------" << endl;
	  Loci::debugout << " My_entities = " << my_entities << endl ;
	  Loci::debugout << " Clone entities = " << clone << endl ;
	  Loci::debugout << "Memory required for the  clone region  = " << clone_size << endl ;
	  Loci::debugout << "Percentage of clone memory required (of size allocated)  = " << double(double(100*clone_size) / size)<< endl ;
	  Loci::debugout << "Percentage of clone memory required (of total size allocated)  = " << double(double(100*clone_size) / total_size)<< endl ;
	  Loci::debugout << "----------------------------------------------------" << endl;
	}
	Loci::debugout << "Percentage of total memory allocated  = " << double(double(100*size) / (total_size+total_wasted)) << endl ;
	Loci::debugout << "----------------------------------------------------" << endl;
	Loci::debugout << "Total wasted size = " << total_wasted << endl ;
	Loci::debugout << "Unused entities = " << unused << endl ;
	Loci::debugout << "Wasted space = " << wasted_space << endl ;
	Loci::debugout << "Percentage of total memory wasted  = " << double(double(100*wasted_space) / (total_size + total_wasted)) << endl ;
	Loci::debugout << " ***************************************************" << endl << endl << endl ;
      }
    }
#endif

    //#define HACK
    for(vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP srp = facts.get_variable(*vi) ;
#ifdef HACK
      entitySet alloc_dom = v_existence[*vi]  + srp->domain() ;
#else
      entitySet alloc_dom = v_requests[*vi] + srp->domain() ;
#endif
      if(srp->domain() == EMPTY) {
	srp->allocate(alloc_dom) ;
      }
      else {
        if(profile_memory_usage) {
          // take off any space that's been reallocated
          if(use_dynamic_memory) {
            int packsize = srp->pack_size(srp->domain()) ;
            LociAppPMTemp -= packsize ;
          }else {//collect reallocation information
            entitySet old_domain = srp->domain() ;
            if(alloc_dom != old_domain) {
              int packsize = srp->pack_size(old_domain) ;
              LociPreallocationReallocSize += packsize ;
              LociPreallocationReallocVars += *vi ;
            }
          }
        }
	if(isSTORE(srp)) {
	  entitySet tmp = interval(alloc_dom.Min(), alloc_dom.Max()) ;
	  if(verbose && tmp.size() >= 2*srp->domain().size())
	    Loci::debugout << "Variable = " << *vi << "  more than twice the space allocated :  allocated over " << alloc_dom << " size = " << tmp.size()  << "  while domain is only  " << srp->domain() << " size = " << srp->domain().size() << endl ;
	  if(alloc_dom != srp->domain()) {
	    if(verbose)
	      Loci::debugout << "reallocating " << *vi << "  over  " << alloc_dom << " initially it was over  " << srp->domain() << endl ;
	    srp->allocate(alloc_dom) ;
	  }
	} 
      }
      // do memory profiling
      if(profile_memory_usage) {
        // we only do profiling for parallel start for dmm
        // since only in parallel and dmm, input vars are resized
        if(use_dynamic_memory) {
          int packsize = srp->pack_size(alloc_dom) ;
          LociAppAllocRequestBeanCounting += packsize ;
          LociAppPMTemp += packsize ;
          if(LociAppPMTemp > LociAppPeakMemoryBeanCounting)
            LociAppPeakMemoryBeanCounting = LociAppPMTemp ;
          if(packsize > LociAppLargestAlloc) {
            LociAppLargestAlloc = packsize ;
            LociAppLargestAllocVar = *vi ;
          }
        }
      }
      
    }
    total_memory_usage = total_size + total_wasted ;
  }
  
  void allocate_all_vars::Print(std::ostream &s) const {
    if(is_alloc_all)
      s << "allocate all variables" << endl ;
    else {
      s << "reallocate all given variables" << endl ;
      if(profile_memory_usage)
        if(use_dynamic_memory)
          s << "memory profiling check point" << endl ;
    }
  }


  executeP graph_compiler::
  execution_schedule(fact_db &facts, sched_db &scheds,
                     const variableSet& alloc) {

    CPTR<execute_list> schedule = new execute_list ;

    if(!use_dynamic_memory) {
      variableSet alloc_vars = facts.get_typed_variables() ;
      if(use_chomp) {
        // we used chomp with no dmm, we don't need to allocate
        // all those chomped variables
        alloc_vars -= all_chomped_vars ;
      }
      schedule->append_list(new allocate_all_vars(facts,scheds,alloc_vars,true)) ;
    }
    else
      if(facts.is_distributed_start())
        if((MPI_processes > 1))
          schedule->append_list(new allocate_all_vars(facts,scheds,alloc,false)) ;

    // we first create a execution module to initialize all keyspaces
#ifdef DYNAMICSCHEDULING
    if(!in_internal_query)
      schedule->append_list(new execute_init_keyspace(facts,scheds)) ;
#endif
    schedule->append_list(fact_db_comm->create_execution_schedule(facts, scheds));
    executeP top_level_schedule = 0;

#ifdef PTHREADS
    // no threading schedule generated for internal query
    if(!in_internal_query) {
      if(threading_pointwise || threading_global_reduction
          || threading_local_reduction || threading_chomping) {
        thread_control =
          new ThreadControl_pthread(num_threads,facts,scheds);
        schedule->append_list(new StartThreads());
      }
    }
#endif
    
    top_level_schedule = (rule_process[baserule])->
      create_execution_schedule(facts, scheds) ;

    if(top_level_schedule == 0) 
      return executeP(0) ;

    schedule->append_list(top_level_schedule) ;

#ifdef PTHREADS
    if(!in_internal_query) {
      if(threading_pointwise || threading_global_reduction
          || threading_local_reduction || threading_chomping) {
        schedule->append_list(new ShutDownThreads());
      }
    }
#endif
    
    if(!use_dynamic_memory)
      if(profile_memory_usage) {
        variableSet profile_vars = facts.get_typed_variables() ;
        // given variables don't need to be profiled
        profile_vars -= LociAppGivenVars ;
        // we also need to exclude those recurrence variables
        profile_vars -= LociRecurrenceVarsRealloc ;
        // but we need to profile those reallocated given vars
        profile_vars += LociPreallocationReallocVars ;
        LociAppPMTemp -= LociPreallocationReallocSize ;
        if(use_chomp) {
          // if use chomping, then chomped vars don't need to be profiled
          profile_vars -= all_chomped_vars ;
        }
        
        schedule->append_list(new execute_memProfileAlloc(profile_vars)) ;
      }
    
    return executeP(schedule) ;
  }

  void clean_empties(digraph &gr,variableSet given,
                     variableSet target, variableSet empties) {
    variableSet work = empties ;

    ruleSet del_rules ;
    for(variableSet::const_iterator vi=work.begin();vi!=work.end();++vi)
      del_rules += extract_rules(gr[vi->ident()]) ;
    variableSet candidates ;
    for(ruleSet::const_iterator ri=del_rules.begin();ri!=del_rules.end();++ri)
      candidates += ri->targets() ;
    
    digraph::vertexSet killvertices = digraph::vertexSet(del_rules) ;
    killvertices += digraph::vertexSet(work) ;
    gr.remove_vertices(killvertices) ;
    //digraph grt = gr.transpose() ;
  }
  
  
  // this is the function to compute all maps
  // in the stationary time level

  // NOTE: the user_query argument is what the user issued
  // for query in the makeQuery function. This is passed in
  // bacause we don't need to issue internal queries for
  // these user_query facts if they happen to be relations also.
#define RENUMBER
  void stationary_relation_gen(rule_db& par_rdb,
                               fact_db& facts,
                               const variableSet& user_query) {
    // we'll first generate the dependency graph
    // according to the initial facts and rules database
    variableSet given = facts.get_typed_variables() ;
    digraph gr ;
    given -= variable("EMPTY") ;
    gr = dependency_graph2(par_rdb,given,user_query).get_graph() ;

    //create_digraph_dot_file(gr,"dependgr.dot") ;
    //std::string cmd = "dotty dependgr.dot" ;
    //system(cmd.c_str()) ;
    
    // If graph is empty, we just return without any further actions
    if(gr.get_target_vertices() == EMPTY)
      return ;

    // we need to collect which relations are to be
    // generated in the stationary time level
    // but we'll need to query the relations according
    // to their topological order since the new
    // relations may change the existential analysis
    vector<digraph::vertexSet> order =
      component_sort(gr).get_components() ;

    // this variable is used to keep track of which set
    // of relations have been reached, we don't need to
    // issue queries for these relations later. This is
    // necessary because it is possible that a relation
    // is generated by multiple rules and the topological
    // sorting of the dependency graph serialized these
    // rules so that we might see the relation in multiple
    // places.
    bool has_map = false ;
    variableSet seen_relations = user_query ;
    vector<pair<ruleSet,variableSet> > relations ;
    for(vector<digraph::vertexSet>::const_iterator vi=order.begin();
        vi!=order.end();++vi) {
      ruleSet rules = extract_rules(*vi) ;

      ruleSet step_rules ;
      variableSet step_relations ;
      for(ruleSet::const_iterator ri=rules.begin();
          ri!=rules.end();++ri) {
        variableSet local_relations ;
        if(ri->type() != rule::INTERNAL) {
          // collect mapping rule and constraint rule target info.
          if( (ri->get_info().rule_impl->get_rule_class() ==
               rule_impl::MAP_RULE) ||
              (ri->get_info().rule_impl->get_rule_class() ==
               rule_impl::CONSTRAINT_RULE)
              ) {
            variableSet targets = ri->targets() ;
            for(variableSet::const_iterator vi=targets.begin();
                vi!=targets.end();++vi) {
              if(vi->time() != time_ident()) {
                cerr << "\tDynamic scheduling of non-stationary time relation"
                     << " is not currently supported. Aborting..." << endl ;
                cerr << "\tAborted because of: " << *vi << endl ;
                Loci::Abort() ;
              }
              if(!has_map) {
                storeRepP srp ;
                srp = ri->get_info().rule_impl->get_store(*vi) ;
                if(isMAP(srp))
                  has_map = true ;
              }
            }
            // this is to ensure that the targets are either
            // not user queried facts or the by-products of
            // user queried facts
            if( (targets & user_query) == EMPTY)
              local_relations += targets ;
          }
        }
        local_relations -= seen_relations ;
        if(local_relations != EMPTY) {
          step_rules += *ri ;
          step_relations += local_relations ;
          seen_relations += local_relations ;
        }
      }
      if(step_relations != EMPTY)
        relations.push_back(make_pair(step_rules,step_relations)) ;
    }
    // if we don't find any relation, we then quit.
    if(relations.empty()) {
      return ;
    }
    // Okay we have gathered all the relations to be generated
    // in the stationary time level. We can now issue the
    // internal queries to compute these relations.

    // we need to collect all the empty constraints if any
    // at the end, we need to remove any rule that involves
    // the empty constraints
    in_internal_query = true ;
    variableSet empty_constraints ;
    // if we don't see any maps being generated, then we
    // can optimize the internal query. i.e., we can group
    // the queries together, we don't need to query according
    // to the topological order as there will not be any
    // changes in the existential and global->local numbering
    // changes during the internal queries.
    if(!has_map) {
      ruleSet all_rules ;
      variableSet all_queries ;
      for(vector<pair<ruleSet,variableSet> >::const_iterator
            vi=relations.begin();vi!=relations.end();++vi) {
        all_rules += vi->first ;
        all_queries += vi->second ;
      }
      fact_db clone(facts) ;
      // before each internal query, we need to
      // perform the global -> local renumbering
      // since the fact_db facts is in global numbering state
#ifdef RENUMBER
      if(clone.is_distributed_start()) {
        if((MPI_processes > 1)) 
          get_clone(clone, par_rdb) ;
        else
          Loci::serial_freeze(clone) ; 
      } else {
        Loci::serial_freeze(clone) ;
      }
#else
      Loci::serial_freeze(clone) ;
#endif
      if(!internalQuery(par_rdb, clone, all_queries)) {
        cerr << "Internal Query Failed, Aborting..." << endl ;
        Loci::Abort() ;
      }
      // then we remove in the rule database the rules
      // that generate the relations
      par_rdb.remove_rules(all_rules) ;
      
      // Okay, now we need to put back the computed relations
      // to the original fact_db and restore the global numbering
#ifdef RENUMBER
      if(clone.is_distributed_start()) {
        fact_db::distribute_infoP df = clone.get_distribute_info() ;
        dMap dl2g ;
        dl2g = MapRepP(df->l2g.Rep())->thaw() ;
        for(variableSet::const_iterator vi2=all_queries.begin();
            vi2!=all_queries.end();++vi2) {
          storeRepP srp = clone.get_variable(*vi2) ;
          if(srp->RepType() == Loci::CONSTRAINT)
            if(GLOBAL_AND(srp->domain()==EMPTY)) {
              empty_constraints += *vi2 ;
              // we don't create empty constraints in
              // the global fact database
              continue ;
            }
          facts.create_intensional_fact(*vi2,(srp->remap(dl2g))->freeze()) ;
        }
      }else{
        for(variableSet::const_iterator vi2=all_queries.begin();
            vi2!=all_queries.end();++vi2) {
          storeRepP srp = clone.get_variable(*vi2) ;
          if(srp->RepType() == Loci::CONSTRAINT)
            if(GLOBAL_AND(srp->domain()==EMPTY)) {
              empty_constraints += *vi2 ;
              continue ;
            }
          facts.create_intensional_fact(*vi2,srp) ;
        }
      }
#else
      for(variableSet::const_iterator vi2=all_queries.begin();
          vi2!=all_queries.end();++vi2) {
        storeRepP srp = clone.get_variable(*vi2) ;
        if(srp->RepType() == Loci::CONSTRAINT)
          if(GLOBAL_AND(srp->domain()==EMPTY)) {
            empty_constraints += *vi2 ;
            continue ;
          }
        facts.create_intensional_fact(*vi2,srp) ;
      }
#endif      
    }else { // if(has_map), then we need to make successive queries
      for(vector<pair<ruleSet,variableSet> >::const_iterator
            vi=relations.begin();vi!=relations.end();++vi) {
        fact_db clone(facts) ;
        // before each internal query, we need to
        // perform the global -> local renumbering
        // since the fact_db facts is in global numbering state
#ifdef RENUMBER
        if(clone.is_distributed_start()) {
          if((MPI_processes > 1)) 
            get_clone(clone, par_rdb) ;
          else
            Loci::serial_freeze(clone) ; 
        } else {
          Loci::serial_freeze(clone) ;
        }
#else
        Loci::serial_freeze(clone) ;
#endif
        // get the relations that we need to query
        ruleSet relationRules = vi->first ;
        variableSet queries = vi->second ;
        if(!internalQuery(par_rdb, clone, queries)) {
          cerr << "Internal Query Failed, Aborting..." << endl ;
          Loci::Abort() ;
        }
        // then we remove in the rule database the rules
        // that generate the relations
        par_rdb.remove_rules(relationRules) ;
      
        // Okay, now we need to put back the computed relations
        // to the original fact_db and restore the global numbering
#ifdef RENUMBER
        if(clone.is_distributed_start()) {
          fact_db::distribute_infoP df = clone.get_distribute_info() ;
          dMap dl2g ;
          dl2g = MapRepP(df->l2g.Rep())->thaw() ;
          for(variableSet::const_iterator vi2=queries.begin();
              vi2!=queries.end();++vi2) {
            storeRepP srp = clone.get_variable(*vi2) ;
            if(srp->RepType() == Loci::CONSTRAINT)
              if(GLOBAL_AND(srp->domain()==EMPTY)) {
                empty_constraints += *vi2 ;
                // we don't create empty constraints in
                // the global fact database
                continue ;
              }
            // we do not need to remap any maps from local -> global
            // because they are generated in the global numbering
            if(isMAP(srp))
              facts.create_intensional_fact(*vi2,srp) ;
            else
              facts.create_intensional_fact(*vi2,srp->remap(dl2g)) ;
          }
        }else{
          for(variableSet::const_iterator vi2=queries.begin();
              vi2!=queries.end();++vi2) {
            storeRepP srp = clone.get_variable(*vi2) ;
            if(srp->RepType() == Loci::CONSTRAINT)
              if(GLOBAL_AND(srp->domain()==EMPTY)) {
                empty_constraints += *vi2 ;
                continue ;
              }
            facts.create_intensional_fact(*vi2,srp) ;
          }
        }
#else
        for(variableSet::const_iterator vi2=queries.begin();
            vi2!=queries.end();++vi2) {
          storeRepP srp = clone.get_variable(*vi2) ;
          if(srp->RepType() == Loci::CONSTRAINT)
            if(GLOBAL_AND(srp->domain()==EMPTY)) {
              empty_constraints += *vi2 ;
              continue ;
            }
          facts.create_intensional_fact(*vi2,srp) ;
        }
#endif        
      }
    } // end of if(!has_map)
    // finally we remove all the rules that derived from
    // any empty constraints
    //    debugout << "empty constraints: " << empty_constraints << endl ;

    // Delete rules dependent on empty constraints if they are not super rules
    // that don't have the AND semantic
    ruleSet del_rules ;
    for(variableSet::const_iterator vi=empty_constraints.begin();
        vi!=empty_constraints.end();++vi) {
      ruleSet empty_rules = extract_rules(gr[vi->ident()]) ;
      for(ruleSet::const_iterator ri=empty_rules.begin();
	  ri != empty_rules.end();++ri) {
	if(!(ri->type() != rule::INTERNAL 
	     &&ri->get_rule_implP()->get_rule_class()==rule_impl::SUPER_RULE)) {
	  del_rules += *ri ;
	}
      }
      //      del_rules += extract_rules(gr[vi->ident()]) ;
    }
    
    par_rdb.remove_rules(del_rules) ;

    in_internal_query = false ;
  }//end of stationary_map_gen


}
