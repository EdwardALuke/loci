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
#include "visitor.h"
#include "visit_tools.h"
#include <sys/stat.h>
#include "comp_tools.h"

#include "dist_tools.h"

#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <set>
using std::set ;
#include <string>
using std::string ;

#include <iostream>
using std::cout ;
using std::cerr;
using std::endl;
using std::ios ;

#include <fstream>
using std::ostringstream ;
using std::ofstream ;

#include <list>
using std::list ;
#include <stack>
using std::stack ;
#include <deque>
using std::deque ;
#include <queue>
using std::queue ;
using std::priority_queue ;

#include <algorithm>
using std::remove_if ;
using std::sort ;
using std::stable_sort ;
//using std::random_shuffle ;

#include <sys/time.h>
//#define VERBOSE

namespace Loci {

  ///////////////////////////////////////////////////////////////////
  // orderVisitor
  ///////////////////////////////////////////////////////////////////
  // Create a schedule for traversing a directed acyclic graph.  This schedule
  // may be concurrent, or many vertices of the graph may be visited at each
  // step of the schedule  If the graph contains cycles, the schedule may
  // not include all of the vertices in the graph.
  vector<digraph::vertexSet>
  orderVisitor::order_dag(const digraph &g,
                          digraph::vertexSet start_vertices,
                          digraph::vertexSet only_vertices) {
    digraph gt = g.transpose() ;
    
    vector<digraph::vertexSet> schedule ; 
    // First schedule any vertices that have no edges leading into them
    //and have not been scheduled previously (in start vertices)

    digraph::vertexSet working = g.get_source_vertices() -
      (g.get_target_vertices()+start_vertices) ;
    working &= only_vertices ;
    if(working != EMPTY) 
      schedule.push_back(working) ;
    
    // visited vertices are all vertices that have already been scheduled
    digraph::vertexSet visited_vertices = start_vertices + working ;
    // In the beginning our working set are all scheduled vertices
    working = visited_vertices ;
    while(working != EMPTY) {
      // While we have vertices to work on, compute additional vertices that
      // can be scheduled
      digraph::vertexSet new_vertices ;
      digraph::vertexSet::const_iterator ni ;
      // loop over working set and create a list of candidate vertices
      for(ni=working.begin();ni != working.end(); ++ni) 
        new_vertices += g[*ni] ;

      // If a vertex has already been scheduled it can't be scheduled again,
      // so remove visited vertices
      new_vertices = new_vertices - visited_vertices    ;
      // We only schedule vertices that are also in the only_vertices set
      working = new_vertices & only_vertices ;
      new_vertices = EMPTY ;
      // Find any vertex from this working set that has had all
      // vertices leading to it scheduled
      for(ni=working.begin();ni != working.end(); ++ni) 
        if((gt[*ni] & visited_vertices) == gt[*ni])
          new_vertices += *ni ;
      working = new_vertices ;
      // and these new vertices to the schedule
      if(new_vertices != EMPTY)
        schedule.push_back(new_vertices) ;
      // update visited vertices set to include scheduled vertices
      visited_vertices += new_vertices ;
    }
    return schedule ;
  }  

  void orderVisitor::visit(loop_compiler& lc) {
    
    lc.collapse_sched = order_dag(lc.collapse_gr) ;
    lc.advance_sched = order_dag(lc.advance_gr) ;
  }

  void orderVisitor::visit(dag_compiler& dc) {
    dc.dag_sched = order_dag(dc.dag_gr) ;
  }

  void orderVisitor::visit(conditional_compiler& cc) {
    cc.dag_sched = order_dag(cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // assembleVisitor
  /////////////////////////////////////////////////////////////////
  namespace {
    class check_dump_on_startup {
      bool checked_dump ;
      bool do_dump ;
    public :
      check_dump_on_startup() {
        do_dump = false ;
        checked_dump = false ;
      }
      bool ok() {
        if(!checked_dump) {
          do_dump=true ;
          struct stat statbuf ;
          if(stat("dump_vars",&statbuf))
            do_dump = false ;
          else if(!S_ISDIR(statbuf.st_mode)) 
            do_dump = false ;
          checked_dump = true ;
        } 
        return do_dump ;
      }
      
    } ;
    
    check_dump_on_startup check_dump_vars ;
    
    class execute_dump_var : public execute_modules {
      variableSet dump_vars ;
      timeAccumulator timer ;
    public:
      execute_dump_var(variableSet vars) : dump_vars(vars) {}
      virtual void execute(fact_db &facts, sched_db &scheds) ;
      virtual void Print(std::ostream &s) const ;
      virtual string getName() { return "execute_dump_var";};
      virtual void dataCollate(collectData &data_collector) const ;
    } ;
    
    map<variable, int> dump_var_lookup ;
    
    void execute_dump_var::execute(fact_db &facts, sched_db &scheds) {
      stopWatch s ;
      s.start() ;
      for(variableSet::const_iterator vi=dump_vars.begin();
          vi!=dump_vars.end();++vi) {
        variable v = *vi ;
        debugout << "dumping variable " << v << endl ;
        storeRepP st = facts.get_variable(v) ;
        storeRepP sc = st ;
        if(isSTORE(st)) {
          ostringstream oss ;
          oss << "dump_vars/"<<v ;
          if(dump_var_lookup.find(v) ==dump_var_lookup.end())
            dump_var_lookup[v] = 0 ;
          if(dump_var_lookup[v] != 0)
            oss << "_"<<dump_var_lookup[v] ;
          dump_var_lookup[v]++ ;
          oss <<".hdf5" ;
          string filename = oss.str() ;
          hid_t file_id=0 ;
          file_id = hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                   H5P_DEFAULT, H5P_DEFAULT) ;
          writeContainer(file_id,v.get_info().name,st,facts) ;
          hdf5CloseFile(file_id) ;
        } else {
          if(MPI_rank == 0) {
            ostringstream oss ;
            oss << "dump_vars/"<<v ;
            if(dump_var_lookup.find(v) ==dump_var_lookup.end())
              dump_var_lookup[v] = 0 ;
            if(dump_var_lookup[v] != 0)
              oss << "."<<dump_var_lookup[v] ;
            dump_var_lookup[v]++ ;
          
            string filename = oss.str() ;
            ofstream ofile(filename.c_str(),ios::out) ;
            ofile.precision(6) ;
            sc->Print(ofile) ;
          }
        }
      }
      timer.addTime(s.stop(),1) ;
    }
    
    void execute_dump_var::Print(std::ostream &s) const {
      s << "dumping variables " << dump_vars << endl ;
    }
    
    void execute_dump_var::dataCollate(collectData &data_collector) const {
      string name ;
      name = "dump vars" ;
      data_collector.accumulateTime(timer,EXEC_CONTROL,name) ;
    }
    
    class dump_vars_compiler : public rule_compiler {
      variableSet dump_vars ;
    public:
      dump_vars_compiler(variableSet &vars) : dump_vars(vars) {}
      virtual void compile() {}
      virtual void accept(visitor& v) {}
      virtual void set_var_existence(fact_db &facts, sched_db &scheds)  {}
      virtual void process_var_requests(fact_db &facts, sched_db &scheds) {}
      virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) {
        return executeP(new execute_dump_var(dump_vars)) ;
      }
    } ;
    
  } // end of namespace

  assembleVisitor::
  assembleVisitor(
#ifdef DYNAMICSCHEDULING
                  fact_db& fd, sched_db& sd,
#endif
                  const variableSet& arv,
                  const std::map<variable,
                  std::pair<rule,CPTR<joiner> > >& ri
#ifdef DYNAMICSCHEDULING
                  ,const variableSet& dt,
                  const map<variable,string>& sc,
                  const map<variable,set<string> >& sac,
                  const map<variable,set<string> >& drc
#endif
                  ):
#ifdef DYNAMICSCHEDULING
    facts(fd),scheds(sd),
#endif
    all_reduce_vars(arv),reduceInfo(ri)
#ifdef DYNAMICSCHEDULING
    ,dynamic_targets(dt),
    self_clone(sc),shadow_clone(sac),drule_ctrl(drc)
#endif
  {
#ifdef DYNAMICSCHEDULING
    // get all the keyspace critical vars
    map<string,KeySpaceP>::iterator ki ;
    for(ki=facts.keyspace.begin();ki!=facts.keyspace.end();++ki) {
      string space_name = ki->second->get_name() ;
      
      variableSet critical_vars = ki->second->get_critical_vars() ;
      variableSet all_synonyms ;

      for(variableSet::const_iterator vi=critical_vars.begin();
          vi!=critical_vars.end();++vi) {
        all_synonyms += scheds.try_get_synonyms(*vi) ;
        all_synonyms += scheds.try_get_aliases(*vi) ;
        all_synonyms += scheds.try_get_rotations(*vi) ;
      }

      critical_vars += all_synonyms ;
      for(variableSet::const_iterator vi=critical_vars.begin();
            vi!=critical_vars.end();++vi)
        var2keyspace[*vi] = space_name ;

      keyspace_critical_vars += critical_vars ;
    } // finish all keyspace init
    // build the dynamic clone data
    map<variable,string>::const_iterator sci ;
    for(sci=self_clone.begin();sci!=self_clone.end();++sci)
      self_clone_vars += sci->first ;
    map<variable,set<string> >::const_iterator sai ;
    for(sai=shadow_clone.begin();sai!=shadow_clone.end();++sai)
      shadow_clone_vars += sai->first ;

    clone_vars = self_clone_vars + shadow_clone_vars ;

    // build the dynamic rule control data
    for(sai=drule_ctrl.begin();sai!=drule_ctrl.end();++sai)
      drule_inputs += sai->first ;
#endif
  }
  
  void assembleVisitor::compile_dag_sched
  (std::vector<rule_compilerP> &dag_comp,
   const std::vector<digraph::vertexSet> &dag_sched,
   const rulecomp_map &rcm,
   const digraph &dag) {

    digraph dagt = dag.transpose() ;
    if(dag_sched.size() == 0)
      return ;
#ifdef VERBOSE
    debugout << "in compile_dag_sched for dag with variables "
             << extract_vars(dag.get_all_vertices()) << endl
             << "dag rules = " << endl
             << extract_rules(dag.get_all_vertices()) << endl ;
#endif
    for(size_t i=0;i<dag_sched.size();++i) {
#ifdef VERBOSE
      debugout << " in comp_dag.cc dag_sched[i] = " << dag_sched[i] << endl ;
#endif
      variableSet vars = extract_vars(dag_sched[i]) ;
#ifdef VERBOSE
      debugout << " in comp_dag.cc vars = " << vars  << endl ;
#endif
      ruleSet rules = extract_rules(dag_sched[i]) ;

      if(rules == EMPTY && i+1<dag_sched.size()) {
        ++i ;
        vars += extract_vars(dag_sched[i]) ;
        rules = extract_rules(dag_sched[i]) ;
      }

      variableSet barrier_vars,reduce_vars,
        singleton_vars,dide_vars,all_vars,dvars ;
      variableSet::const_iterator vi ;
      
      for(vi=vars.begin();vi!=vars.end();++vi) {
        ruleSet var_rules = extract_rules(dagt[(*vi).ident()]) ;
        ruleSet::const_iterator ri ;
        ruleSet use_rules ;
        bool reduction = false ;
        bool pointwise = false ;
        bool singleton = false ;
        bool recursive = false ;
        bool dynamic_ide = false ;
        bool unit_rule_exists = false ;
        bool priority_rule = false ;

        for(ri=var_rules.begin();ri!=var_rules.end();++ri) {
          if(!is_super_node(ri))
            dvars += *vi ;
          
          if(!is_virtual_rule(*ri) ||
             (ri->get_info().rule_class == rule::INTERNAL &&
              ri->get_info().qualifier() == "priority")) {
            use_rules += *ri ;

            // Check for a priority rule
            if(ri->get_info().rule_class == rule::INTERNAL &&
               ri->get_info().qualifier() == "priority")
              priority_rule = true ;

            if(ri->get_info().rule_class == rule::INTERNAL) {
              if(is_chomp_node(ri)) {
                // we need to actually look into the chomping
                // graph to find out the rules that generate
                // this variable and dertermin the types of
                // rules
                rulecomp_map::const_iterator rmi ;
                rmi = rcm.find(*ri) ;
                FATAL(rmi == rcm.end()) ;
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

                  if(crimp->get_rule_class() == rule_impl::INSERTION ||
                     crimp->get_rule_class() == rule_impl::DELETION ||
                     crimp->get_rule_class() == rule_impl::ERASE)
                    dynamic_ide = true ;
                }
                //continue ;
              }
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

            if(rimp->get_rule_class() == rule_impl::INSERTION ||
               rimp->get_rule_class() == rule_impl::DELETION ||
               rimp->get_rule_class() == rule_impl::ERASE)
              dynamic_ide = true ;
            
          } else {
            if((ri->sources() & ri->targets()) != EMPTY)
              recursive = true ;
          }
        }

        WARN(((reduction && pointwise) || (pointwise && singleton) ||
              (reduction && singleton)) && vi->get_info().name != "OUTPUT") ;
        if(((reduction && pointwise) || (pointwise && singleton) ||
            (reduction && singleton)) && vi->get_info().name != "OUTPUT") {
          cerr << "Warning: invalid mix of rule types for variable " << *vi << endl ;
        }

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

          /*
          if(all_reduce_vars.inSet(*vi))
          reduce_vars += *vi ;
          */
          if(reduction && unit_rule_exists)
            reduce_vars += *vi ;
          
          if(singleton) {
            singleton_vars += *vi ;
          }

          if(dynamic_ide) {
            dide_vars += *vi ;
          }
        }

      }

      all_vars += barrier_vars ;

#ifdef DYNAMICSCHEDULING
      // remove dynamic targets
      barrier_vars -= dynamic_targets ;
#endif

      if(barrier_vars != EMPTY) {
        dag_comp.push_back(new barrier_compiler(barrier_vars)) ;
      }
      
      all_vars += singleton_vars ;

#ifdef DYNAMICSCHEDULING
      singleton_vars -= dynamic_targets ;
#endif
      if(singleton_vars != EMPTY)
        dag_comp.push_back(new singleton_var_compiler(singleton_vars)) ;

      all_vars += reduce_vars;

      all_vars += dide_vars ;

#ifdef DYNAMICSCHEDULING
      reduce_vars -= dynamic_targets ;
#endif
#ifdef VERBOSE
      debugout << "all_vars = " << all_vars << endl ;
      debugout << "singleton_vars =" <<singleton_vars << endl ;
      debugout << "barrier_vars = " << barrier_vars << endl ;
      debugout << "reduce_vars = " << reduce_vars << endl ;
#endif

      vector<CPTR<joiner> > join_op_vector ;
      vector<rule> unit_rule_vector ;
      vector<variable> reduce_var_vector ;
      
      for(variableSet::const_iterator rdvi=reduce_vars.begin();
          rdvi!=reduce_vars.end();++rdvi) {
        map<variable,pair<rule,CPTR<joiner> > >::const_iterator xi ;
        xi = reduceInfo.find(*rdvi) ;
        FATAL(xi == reduceInfo.end()) ;
        rule unit_rule = (xi->second).first ;
        CPTR<joiner> join_op = (xi->second).second ;

        if(join_op != 0) {
          storeRepP sp = join_op->getTargetRep() ;
          /* old code
             if(isPARAMETER(sp)) {
             dag_comp.push_back(new reduce_param_compiler(xi->first,unit_rule,
             join_op)) ;
             }
          */
          if(sp!=0) {
            if(isPARAMETER(sp)) {
              reduce_var_vector.push_back(xi->first) ;
              unit_rule_vector.push_back(unit_rule) ;
              join_op_vector.push_back(join_op) ;
            } else {
              warn(!isSTORE(sp)) ;
#ifdef VERBOSE
              debugout << "reduce_store_compiler("
                       << xi->first << ","
                       << unit_rule << ")" << endl ;
#endif
              dag_comp.push_back(new reduce_store_compiler(xi->first,unit_rule,
                                                           join_op)) ;
            }
          }
        }
      }
      if(reduce_var_vector.size() != 0)  
        dag_comp.push_back(new reduce_param_compiler(reduce_var_vector, unit_rule_vector, join_op_vector));

#ifdef DYNAMICSCHEDULING
      // check for dynamic clone invalidator first
      for(variableSet::const_iterator dci=dvars.begin();
          dci!=dvars.end();++dci) {
        variable t = facts.remove_synonym(*dci) ;
        if(clone_vars.inSet(t)) {
          // create a compiler for the task
          KeySpaceP self_clone_kp = 0 ;
          vector<KeySpaceP> shadow_clone_kp ;

          map<variable,string>::const_iterator sci ;
          sci = self_clone.find(t) ;
          if(sci != self_clone.end()) {
            self_clone_kp = facts.get_keyspace(sci->second) ;
            if(self_clone_kp == 0) {
              cerr << "Error: keyspace: " << sci->second
                   << " does not exist, error from assembleVisitor"
                   << endl ;
              Loci::Abort() ;
            }
          }

          map<variable,set<string> >::const_iterator sai ;
          sai = shadow_clone.find(t) ;
          if(sai != shadow_clone.end()) {
            const set<string>& space_names = sai->second ;
            for(set<string>::const_iterator si=space_names.begin();
                si!=space_names.end();++si) {
              KeySpaceP kp = facts.get_keyspace(*si) ;
              if(kp == 0) {
                cerr << "Error: keyspace: " << *si
                     << " does not exist, error from assembleVisitor"
                     << endl ;
                Loci::Abort() ;
              }
              shadow_clone_kp.push_back(kp) ;
            }
          }

          dag_comp.push_back
            (new dclone_invalidate_compiler(*dci, t,
                                            self_clone_kp, shadow_clone_kp)) ;
        }
        // then check for dynamic rule control relations
        if(drule_inputs.inSet(*dci)) {
          // create a compiler for the task
          vector<KeySpaceP> register_kp ;
          map<variable,set<string> >::const_iterator msi ;
          msi = drule_ctrl.find(*dci) ;
          if(msi != drule_ctrl.end()) {
            const set<string>& space_names = msi->second ;
            dag_comp.push_back
              (new dcontrol_reset_compiler(*dci, t, space_names)) ;
          }
        }
      }

      // check for automatic keyspace distribution
      variableSet level_critical =
        variableSet(all_vars & keyspace_critical_vars) ;
      set<string> keyspaces ;
      for(variableSet::const_iterator lcvi=level_critical.begin();
          lcvi!=level_critical.end();++lcvi) {
        map<variable,string>::const_iterator mi =
          var2keyspace.find(*lcvi) ;
        FATAL(mi == var2keyspace.end()) ;
        keyspaces.insert(mi->second) ;
      }
      if(!keyspaces.empty()) {
        vector<string> keyspace_dist ;
        for(set<string>::const_iterator si=keyspaces.begin();
            si!=keyspaces.end();++si)
          keyspace_dist.push_back(*si) ;

        dag_comp.push_back(new keyspace_dist_compiler(keyspace_dist)) ;
      }
#endif
      if(check_dump_vars.ok()) {
        if(all_vars != EMPTY) 
          dag_comp.push_back(new dump_vars_compiler(all_vars)) ;
      }
      
      if(rules != EMPTY) {
        ruleSet::const_iterator ri ;
        for(ri=rules.begin();ri!=rules.end();++ri) {
          rulecomp_map::const_iterator rmi ;
          rmi = rcm.find(*ri) ;
          if(rmi == rcm.end()) {
            cerr << "could not find rule compiler for " <<*ri << endl ;
          }
          FATAL(rmi == rcm.end()) ;
          dag_comp.push_back(rmi->second) ;
        }

#ifdef DYNAMICSCHEDULING
        // check for manual keyspace distribution
        for(ri=rules.begin();ri!=rules.end();++ri) {
          if(is_internal_rule(ri))
            continue ;
          rule_implP rp = ri->get_rule_implP() ;
          if(rp->affect_keyspace_dist())
            dag_comp.push_back
              (new keyspace_dist_compiler(rp->gather_keyspace_dist())) ;
        }
#endif
      }
    }
  }

  void assembleVisitor::visit(loop_compiler& lc) {
#ifdef VERBOSE
    debugout << "collapse comp dag scheduler:" << endl ;
#endif
    compile_dag_sched(lc.collapse_comp,lc.collapse_sched,
                      lc.rule_compiler_map,lc.collapse_gr) ;
#ifdef VERBOSE
    debugout << "advance comp dag scheduler:" << endl ;
#endif
    compile_dag_sched(lc.advance_comp,lc.advance_sched,
                      lc.rule_compiler_map,lc.advance_gr) ;
  }
  
  void assembleVisitor::visit(dag_compiler& dc) {
    compile_dag_sched(dc.dag_comp,dc.dag_sched,
                      dc.rule_compiler_map,dc.dag_gr) ;
  }
  
  void assembleVisitor::visit(conditional_compiler& cc) {
    compile_dag_sched(cc.dag_comp,cc.dag_sched,
                      cc.rule_compiler_map,cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // utility functions
  /////////////////////////////////////////////////////////////////
  namespace {

    vector<digraph::vertexSet> insertAlloc2Sched
    (const digraph &gr,
     const vector<digraph::vertexSet>& firstSched) {
      // then we schedule those alloc rules
      digraph grt = gr.transpose() ;
      vector<digraph::vertexSet> finalSched ;
      
      digraph::vertexSet last ;
      digraph::vertexSet sched_alloc ;
      for(vector<digraph::vertexSet>::const_iterator vi=firstSched.begin();
          vi!=firstSched.end();++vi) {
        digraph::vertexSet step = *vi ;
        ruleSet step_rules = extract_rules(step) ;
        
        digraph::vertexSet pre_vertices ;
        for(ruleSet::const_iterator ri=step_rules.begin();
            ri!=step_rules.end();++ri)
          pre_vertices += grt[ri->ident()] ;
        
        ruleSet pre_rules = extract_rules(pre_vertices) ;
        digraph::vertexSet needed_alloc ;
        for(ruleSet::const_iterator ri=pre_rules.begin();
            ri!=pre_rules.end();++ri)
          if(ri->get_info().qualifier() == "ALLOCATE")
            needed_alloc += ri->ident() ;
        
        needed_alloc -= sched_alloc ;
        sched_alloc += needed_alloc ;
        
        digraph::vertexSet sum = needed_alloc + last ;
        if(sum != EMPTY)
          finalSched.push_back(sum) ;
        
        last = step ;
      }
      if(last != EMPTY)
        finalSched.push_back(last) ;
      
      return finalSched ;
    }

    
  } // end of namespace

  /////////////////////////////////////////////////////////////////
  // simLazyAllocSchedVisitor
  /////////////////////////////////////////////////////////////////
  std::vector<digraph::vertexSet>
  simLazyAllocSchedVisitor::get_firstSched(const digraph& gr) const {
    // first we get all vertices except all the allocate rules
    ruleSet rules = extract_rules(gr.get_all_vertices()) ;
    
    digraph::vertexSet alloc_rules_vertices ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      if(ri->get_info().qualifier() == "ALLOCATE")
        alloc_rules_vertices += ri->ident() ;
    
    vector<digraph::vertexSet> first_sched ;
    
    first_sched =
      orderVisitor::order_dag(gr,alloc_rules_vertices,
                              gr.get_all_vertices()-alloc_rules_vertices) ;

    return first_sched ;
  }

  void simLazyAllocSchedVisitor::visit(loop_compiler& lc) {
    vector<digraph::vertexSet> first_sched ;

    first_sched = get_firstSched(lc.collapse_gr) ;
    lc.collapse_sched = insertAlloc2Sched(lc.collapse_gr,first_sched) ;

    first_sched = get_firstSched(lc.advance_gr) ;
    lc.advance_sched = insertAlloc2Sched(lc.advance_gr,first_sched) ;
  }

  void simLazyAllocSchedVisitor::visit(dag_compiler& dc) {
    vector<digraph::vertexSet> first_sched = get_firstSched(dc.dag_gr) ;
    dc.dag_sched = insertAlloc2Sched(dc.dag_gr,first_sched) ;
  }

  void simLazyAllocSchedVisitor::visit(conditional_compiler& cc) {
    vector<digraph::vertexSet> first_sched = get_firstSched(cc.cond_gr) ;
    cc.dag_sched = insertAlloc2Sched(cc.cond_gr,first_sched) ;
  }

  /////////////////////////////////////////////////////////////////
  // memGreedySchedVisitor
  /////////////////////////////////////////////////////////////////
  std::vector<digraph::vertexSet>
  memGreedySchedVisitor::get_firstSched(const digraph& gr) {
    ruleSet rules = extract_rules(gr.get_all_vertices()) ;    
    digraph::vertexSet alloc_rules_vertices ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      if(ri->get_info().qualifier() == "ALLOCATE")
        alloc_rules_vertices += ri->ident() ;

    digraph grt = gr.transpose() ;
    
    vector<digraph::vertexSet> schedule ; 
    digraph::vertexSet waiting = gr.get_source_vertices() -
      gr.get_target_vertices() - alloc_rules_vertices ;
    
    // visited vertices are all vertices that have already been scheduled
    digraph::vertexSet visited_vertices = alloc_rules_vertices ;
    while(waiting != EMPTY) {
      // we need to schedule the waiting set of vertices
      // but we assign priority to vertices insice the waiting set
      // we first schedule all variables, delete rules, rules that cause
      // deletion, or rules that don't change memory status.
      digraph::vertexSet step_sched ;
      for(digraph::vertexSet::const_iterator vi=waiting.begin();
          vi!=waiting.end();++vi) {
        if(*vi >= 0)
          step_sched += *vi ;
        else {
          rule r(*vi) ;
          if(r.get_info().qualifier() == "DELETE")
            step_sched += *vi ;
          else {
            // look back
            ruleSet pre_rules = extract_rules(grt[*vi]) ;
            bool has_alloc = false ;
            for(ruleSet::const_iterator ri=pre_rules.begin();
                ri!=pre_rules.end();++ri)
              if(ri->get_info().qualifier() == "ALLOCATE") {
                variableSet targets = ri->targets() ;
                for(variableSet::const_iterator vi2=targets.begin();
                    vi2!=targets.end();++vi2) {
                  storeRepP srp = facts.get_variable(*vi2) ;
                  if(isSTORE(srp)) {
                    has_alloc = true ;
                    break ;
                  }
                }
              }
            if(!has_alloc)
              step_sched += *vi ;
          }
        }
      }
      step_sched -= visited_vertices ;
      digraph::vertexSet valid_sched ;
      for(digraph::vertexSet::const_iterator vi=step_sched.begin();
          vi!=step_sched.end();++vi)
        if( (grt[*vi] - visited_vertices) == EMPTY)
          valid_sched += *vi ;
      
      step_sched = valid_sched ;
      // if any thing inside, we schedule this set
      if(step_sched != EMPTY) {
        digraph::vertexSet new_vertices ;
        for(digraph::vertexSet::const_iterator vi=step_sched.begin();
            vi!=step_sched.end();++vi)
          new_vertices += gr[*vi] ;
        waiting -= step_sched ;
        waiting += new_vertices ;
        visited_vertices += step_sched ;
        schedule.push_back(step_sched) ;
        continue ;
      }

      // now, step_sched must be empty,
      // we pick a rule from the waiting set
      // all rules left now must lead to allocation, we first
      // check to see if there's any rule that lead to deletion also.
      // we schedule the one that lead to most deletions
      int num_del = 0 ;
      int num_alloc = UNIVERSE_MAX ;
      int_type sched_v = 0 ;
      for(digraph::vertexSet::const_iterator vi=waiting.begin();
          vi!=waiting.end();++vi) {
        if(*vi < 0) { // a rule
          if( (grt[*vi] - visited_vertices) != EMPTY)
            continue ;
          int local_num_del = 0 ;
          ruleSet next = extract_rules(gr[*vi]) ;
          for(ruleSet::const_iterator ri=next.begin();
              ri!=next.end();++ri) {
            if(ri->get_info().qualifier() == "DELETE") {
              variableSet targets = ri->targets() ;
              for(variableSet::const_iterator vi2=targets.begin();
                  vi2!=targets.end();++vi2) {
                storeRepP srp = facts.get_variable(*vi2) ;
                if(isSTORE(srp)) {
                  ++local_num_del ;
                  break ; // !!!!!!!!!!!!!!modified!!!!!!!!!!!!!!!
                }
              }
            }
          }
          int local_num_alloc = 0 ;
          ruleSet pre = extract_rules(grt[*vi]) ;
          for(ruleSet::const_iterator ri=pre.begin();
              ri!=pre.end();++ri) {
            if(ri->get_info().qualifier() == "ALLOCATE") {
              variableSet targets = ri->targets() ;
              for(variableSet::const_iterator vi2=targets.begin();
                  vi2!=targets.end();++vi2) {
                storeRepP srp = facts.get_variable(*vi2) ;
                if(isSTORE(srp)) {
                  ++local_num_alloc ;
                  break ; // !!!!!!!!!!!!!!modified!!!!!!!!!!!!!!!
                }
              }
            }
          }

          if(local_num_del > num_del) {
            num_del = local_num_del ;
            sched_v = *vi ;
            num_alloc = local_num_alloc ;
          }else if( (local_num_del != 0) && (local_num_del == num_del))
            //}if(local_num_del == num_del)
            if(local_num_alloc < num_alloc) {
              cout << "I am working" << endl ;
              num_alloc = local_num_alloc ;
              sched_v = *vi ;
            }
        }
      }
      // see if we found any to schedule
      if(sched_v < 0) {
        if(visited_vertices.inSet(sched_v)) {
          cerr << "ERROR in scheduling..." << endl ;
          Loci::Abort() ;
        }
        step_sched += sched_v ;
        waiting -= step_sched ;
        waiting += gr[sched_v] ;
        visited_vertices += step_sched ;
        schedule.push_back(step_sched) ;
        continue ;        
      }
      // otherwise we pick one rule from the rest vertex set
      // we pick one that has fewest outgoing edges from all
      // of its target variables that are store
      int out_edges = UNIVERSE_MAX ;
      for(digraph::vertexSet::const_iterator vi=waiting.begin();
          vi!=waiting.end();++vi) {
        if(*vi < 0) { // a rule
          if( (grt[*vi] - visited_vertices) != EMPTY)
            continue ;
          rule r(*vi) ;
          variableSet targets = r.targets() ;
          int local_out_edges = 0 ;
          for(variableSet::const_iterator vi2=targets.begin();
              vi2!=targets.end();++vi2) {
            storeRepP srp = facts.get_variable(*vi2) ;
            if(isSTORE(srp)) {
              digraph::vertexSet nextv = gr[vi2->ident()] ;
              local_out_edges += nextv.size() ;
            }
          }
          if(local_out_edges < out_edges) {
            out_edges = local_out_edges ;
            sched_v = *vi ;
          }
        }
      }
      // if there's no scheduling candidate, the graph
      // must have problem
      if(sched_v >= 0) {
        cout << "Cannot proceed in graph scheduling!"
             << " The graph may have cycles!" << endl ;
        Loci::Abort() ;
      }
      // we schedule this rule
      if(visited_vertices.inSet(sched_v)) {
        cerr << "ERROR in scheduling..." << endl ;
        Loci::Abort() ;
      }
      step_sched += sched_v ;
      waiting -= step_sched ;
      waiting += gr[sched_v] ;
      visited_vertices += step_sched ;
      schedule.push_back(step_sched) ;
    }
    return schedule ;
  }

  void memGreedySchedVisitor::visit(loop_compiler& lc) {
    vector<digraph::vertexSet> first_sched ;
    
    first_sched = get_firstSched(lc.collapse_gr) ;
    lc.collapse_sched = insertAlloc2Sched(lc.collapse_gr,first_sched) ;

    first_sched = get_firstSched(lc.advance_gr) ;
    lc.advance_sched = insertAlloc2Sched(lc.advance_gr,first_sched) ;
  }

  void memGreedySchedVisitor::visit(dag_compiler& dc) {
    vector<digraph::vertexSet> first_sched = get_firstSched(dc.dag_gr) ;
    dc.dag_sched = insertAlloc2Sched(dc.dag_gr,first_sched) ;
  }

  void memGreedySchedVisitor::visit(conditional_compiler& cc) {
    vector<digraph::vertexSet> first_sched = get_firstSched(cc.cond_gr) ;
    cc.dag_sched = insertAlloc2Sched(cc.cond_gr,first_sched) ;
  }
  


  void SchedClearVisitor::visit(loop_compiler& lc) {
   
    (lc.collapse_sched).clear();
    (lc.advance_sched).clear();
    (lc.collapse_comp).clear();
    (lc.advance_comp).clear();
    
  }

  void SchedClearVisitor::visit(dag_compiler& dc) {
  
    (dc.dag_sched).clear();
    (dc.dag_comp).clear();
  }

  void SchedClearVisitor::visit(conditional_compiler& cc) {
    (cc.dag_sched).clear();
    (cc.dag_comp).clear();
  }
  
  
  /////////////////////////////////////////////////////////////////
  // graphSchedulerVisitor
  /////////////////////////////////////////////////////////////////
  // utility data-structure and functions
  namespace {
    struct sched_item {
      int_type vertex_id ;
      int_type weight ;
      sched_item(int_type v=0,int_type w=0):vertex_id(v),weight(w) {}
    } ;
    struct comp_sched_item {
      bool operator()(const sched_item& si1,const sched_item& si2) const
      { return si1.weight > si2.weight ;}
    } ;
    
    inline void create_queue(const digraph::vertexSet& vs,
                             const map<int,int>& pmap,
                             priority_queue<sched_item,
                               vector<sched_item>,comp_sched_item>& q) {
      for(digraph::vertexSet::const_iterator vi=vs.begin();
          vi!=vs.end();++vi) {
        map<int_type,int_type>::const_iterator mi ;
        mi = pmap.find(*vi) ;
        FATAL(mi == pmap.end()) ;
        sched_item item(*vi,mi->second) ;
        q.push(item) ;
      }
    }

    inline digraph::vertexSet
    pop_queue(priority_queue<sched_item,vector<sched_item>,
              comp_sched_item>& q) {
      if(q.empty()) return EMPTY ;
      
      digraph::vertexSet ret ;
      sched_item si ; // current poped item
      int_type pre_w ; // previous weight

      si = q.top() ;
      q.pop() ;
      ret += si.vertex_id ;
      pre_w = si.weight ;
      while(!q.empty()) {
        si = q.top() ;
        if(si.weight == pre_w) {
          q.pop() ;
          ret += si.vertex_id ;
        } else
          break ;
      }

      return ret ;
    }

    inline digraph::vertexSet
    get_valid_sched(const digraph::vertexSet& vs,
                    const digraph& gt,
                    const digraph::vertexSet& visited) {
      digraph::vertexSet valid ;
      for(digraph::vertexSet::const_iterator vi=vs.begin();
          vi!=vs.end();++vi) {
        digraph::vertexSet pre = gt[*vi] ;
        if( (pre-visited) == EMPTY)
          valid += *vi ;
      }
      // get rid of those already scheduled
      valid -= visited ;
      
      return valid ;
    }

    inline digraph::vertexSet
    get_new_reachable(const digraph::vertexSet& vs,
                      const digraph& g) {
      digraph::vertexSet reachable ;
      for(digraph::vertexSet::const_iterator vi=vs.begin();
          vi!=vs.end();++vi)
        reachable += g[*vi] ;

      return reachable ;
    }

    inline digraph allocFreeGr(const digraph& gr) {
      digraph::vertexSet allv = gr.get_all_vertices() ;
      ruleSet rules = extract_rules(gr.get_all_vertices()) ;    
      digraph::vertexSet alloc_rules_vertices ;
      for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
        if(ri->get_info().qualifier() == "ALLOCATE")
          alloc_rules_vertices += ri->ident() ;

      digraph new_gr = gr.subgraph(allv-alloc_rules_vertices) ;
      return new_gr ;
    }

  } // end of namespace (unnamed)
  
  std::vector<digraph::vertexSet>
  graphSchedulerVisitor::schedule(const digraph& gr) {
    // we schedule the graph without allocate rules first
    digraph gwoa = allocFreeGr(gr) ;
    digraph gwoat = gwoa.transpose() ;
    digraph::vertexSet allv = gwoa.get_all_vertices() ;

    // create a priority map
    map<int_type,int_type> pmap ;
    for(digraph::vertexSet::const_iterator vi=allv.begin();
        vi!=allv.end();++vi)
      pmap[*vi] = 0 ;

    // prioritize the graph
    prio(gr,pmap) ;
    // begin schedule
    vector<digraph::vertexSet> sched ;
    digraph::vertexSet visited ;
    digraph::vertexSet waiting = gwoa.get_source_vertices()
      - gwoa.get_target_vertices() ;
    
    while(waiting != EMPTY) {
      // create a priority queue
      priority_queue<sched_item,vector<sched_item>,comp_sched_item> q ;
      create_queue(waiting,pmap,q) ;
      FATAL(q.size() != (unsigned int)waiting.size()) ;
      digraph::vertexSet step_sched ;
      while( (!q.empty()) && (step_sched == EMPTY)) {
        step_sched = pop_queue(q) ;
        step_sched = get_valid_sched(step_sched,gwoat,visited) ;
      }
      if( (q.empty()) && (step_sched == EMPTY)) {
        cerr << "Loci Graph Scheduler Warning: input graph has cycle(s)"
             << endl ;
        cerr << "waiting rules = " << extract_rules(waiting) << endl ;
        cerr << "waiting variables = " << extract_vars(waiting) << endl ;
        Loci::Abort() ;
      }
      waiting -= step_sched ;
      waiting += get_new_reachable(step_sched,gwoa) ;
      visited += step_sched ;
      sched.push_back(step_sched) ;
    }

    // here we insert all the allocate rule back into
    // the schedule
    sched = insertAlloc2Sched(gr,sched) ;

    return sched ;    
  }

  void graphSchedulerVisitor::visit(loop_compiler& lc) {
    lc.collapse_sched = schedule(lc.collapse_gr) ;
    lc.advance_sched = schedule(lc.advance_gr) ;
  }

  void graphSchedulerVisitor::visit(dag_compiler& dc) {
    dc.dag_sched = schedule(dc.dag_gr) ;
  }

  void graphSchedulerVisitor::visit(conditional_compiler& cc) {
    cc.dag_sched = schedule(cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // computation greedy prioritize
  /////////////////////////////////////////////////////////////////
  void compGreedyPrio::operator()(const digraph& gr,
                                  map<int_type,int_type>& pmap) const {
    digraph::vertexSet allv = gr.get_all_vertices() ;
    for(digraph::vertexSet::const_iterator vi=allv.begin();
        vi!=allv.end();++vi)
      pmap[*vi] = 0 ;
  }

  namespace {
    // utility data-structure and functions
    struct stat_item {
      int a ; // allocation rules attached
      int d ; // deletion rules attached
      int o ; // outgoing edges number for all targets
      int_type id ; // vertex id
      stat_item(int a=0,int d=0,int o=0,int_type id=0):
        a(a),d(d),o(o),id(id) {}
    } ;
    // different comparison functions
    bool ascend_a(const stat_item& s1, const stat_item& s2)
    {return s1.a < s2.a ;}
    bool descend_d(const stat_item& s1, const stat_item& s2)
    {return s1.d > s2.d ;}
    bool ascend_o(const stat_item& s1, const stat_item& s2)
    {return s1.o < s2.o ;}
    // predicate used in remove_if algorithm
    struct comp_rm {
      comp_rm(const digraph::vertexSet& vs): vset(vs) {}
      bool operator()(const stat_item& s)
      {return vset.inSet(s.id) ;}
    private:
      digraph::vertexSet vset ;
    } ;
  } // end of namespace (unnamed)
  
  void memGreedyPrio::operator()(const digraph& gr,
                                 map<int_type,int_type>& pmap) const {
    digraph::vertexSet allv = gr.get_all_vertices() ;
    digraph grt = gr.transpose() ;
    
    vector<stat_item> l ;
    for(digraph::vertexSet::const_iterator vi=allv.begin();
        vi!=allv.end();++vi) {
      if(*vi >= 0) // a variable
        l.push_back(stat_item(0,0,0,*vi)) ;
      else { // a rule
        ruleSet pre = extract_rules(grt[*vi]) ;
        ruleSet next = extract_rules(gr[*vi]) ;
        rule r(*vi) ;
        variableSet targets = r.targets() ;
        
        int a=0,d=0,o=0 ;
        // first count for allocation rules attached
        for(ruleSet::const_iterator ri=pre.begin();
            ri!=pre.end();++ri)
          if(ri->get_info().qualifier() == "ALLOCATE") {
            variableSet tvars = ri->targets() ;
            for(variableSet::const_iterator invi=tvars.begin();
                invi!=tvars.end();++invi) {
              storeRepP srp = facts.get_variable(*invi) ;
              if(isSTORE(srp)) {
                ++a ;
                break ;
              }
            }
          }
        // count the deletion rules attached
        for(ruleSet::const_iterator ri=next.begin() ;
            ri!=next.end();++ri)
          if(ri->get_info().qualifier() == "DELETE") {
            variableSet tvars = ri->targets() ;
            for(variableSet::const_iterator invi=tvars.begin();
                invi!=tvars.end();++invi) {
              storeRepP srp = facts.get_variable(*invi) ;
              if(isSTORE(srp)) {
                ++d ;
                break ;
              }
            }
          }
        // finally count the outgoing edges number
        for(variableSet::const_iterator vari=targets.begin() ;
            vari!=targets.end();++vari) {
          storeRepP srp = facts.get_variable(*vari) ;
          if(isSTORE(srp)) {
            digraph::vertexSet nxv = gr[vari->ident()] ;
            o += nxv.size() ;
          }
        }
        // push to the vector
        l.push_back(stat_item(a,d,o,*vi)) ;
      } // end of else
    } // end of for(allv)

    // now we start to assign weight
    // first pick out those with no allocations
    vector<stat_item>::const_iterator vi ;
    digraph::vertexSet remove ;
    vector<stat_item>::size_type size ;
    vector<stat_item>::iterator old_end, new_end ;
    int weight = 0 ;
    
    for(vi=l.begin();vi!=l.end();++vi) {
      stat_item s = *vi ;
      if(s.a == 0) {
        pmap[s.id] = weight ;
        remove += s.id ;
      }
    }
    // remove the assigned vertices
    size = l.size() ; 
    old_end = l.end() ;
    new_end = remove_if(l.begin(),l.end(),comp_rm(remove)) ;
    l.resize(size-(old_end-new_end)) ;
    if(l.empty())
      return ;

    // we sort the remaining vector
    // first ascending order of a
    sort(l.begin(),l.end(),ascend_a) ;
    // then stable sort according to descending order of d
    stable_sort(l.begin(),l.end(),descend_d) ;
    
    weight = 1 ;
    remove = EMPTY ;
    for(vi=l.begin();vi!=l.end();++vi) {
      stat_item s = *vi ;
      if(s.d > 0) {
        pmap[s.id] = weight++ ;
        remove += s.id ;
      }
    }
    // remove processed vertices
    size = l.size() ; 
    old_end = l.end() ;
    new_end = remove_if(l.begin(),l.end(),comp_rm(remove)) ;
    l.resize(size-(old_end-new_end)) ;
    if(l.empty())
      return ;

    // then sort according to ascending order of o
    sort(l.begin(),l.end(),ascend_o) ;

    for(vi=l.begin();vi!=l.end();++vi) {
      stat_item s = *vi ;
      pmap[s.id] = weight++ ;
    }

  } // end of function

  
  // graph prioritize function that tries
  // to maximize the cache benefit of chomping
  void chompingPrio::operator()(const digraph& gr,
                                map<int_type,int_type>& pmap) const {
    // initial weight
    int_type weight = 0 ;
    // get a depth first order by doing a component sort
    vector<digraph::vertexSet> components =
      component_sort(gr).get_components() ;
    // looping over the depth first order to assign weight
    vector<digraph::vertexSet>::const_iterator vi ;
    for(vi=components.begin();vi!=components.end();++vi) {
      if(vi->size() != 1) {
        cerr << "Chomping graph error: Cycle(s) detected" << endl ;
        Loci::Abort() ;

      }
      pmap[*(vi->begin())] = weight++ ;
    }
  }


  namespace {
    // a small function to generate random integers from
    // a specified range [b,e]
    size_t gen_rand(size_t b, size_t e) {
      double o = drand48() ;
      size_t scale = e - b + 1 ;
      return b + static_cast<size_t>(o*scale) ;
    }
  }
  
  MemGreedyScheduler::
  MemGreedyScheduler(fact_db& fd,
                     sched_db& sd,
                     const std::map<variable,variableSet>& s,
                     const std::map<variable,variableSet>& t,
                     const map<variable, double>& info)
    :s2t(s),t2s(t),var_cluster(5) {
    // initialize random seed
    srand48(s.size()*t.size()) ;

    int max_id = 0 ;

    for(map<variable, double>::const_iterator mi = info.begin(); mi != info.end(); mi++){
      if((mi->first).ident() > max_id)max_id = (mi->first).ident();
    }
      

     
    // construct the var_cluster info
    var_cluster.resize(max_id+1) ;
    for(map<variable,double>::const_iterator
          mi=info.begin();mi!=info.end();++mi) {
      var_cluster.makeSet( (mi->first).ident()) ;
    }
    for(map<variable,variableSet>::const_iterator
          mi=s2t.begin();mi!=s2t.end();++mi) {
      var_cluster.tryMakeSetWithResize( (mi->first).ident()) ;
      for(variableSet::const_iterator vi=(mi->second).begin();
          vi!=(mi->second).end();++vi)
        var_cluster.tryMakeSetWithResize(vi->ident()) ;
    }
    for(map<variable,variableSet>::const_iterator
          mi=t2s.begin();mi!=t2s.end();++mi) {
      var_cluster.tryMakeSetWithResize( (mi->first).ident()) ;
      for(variableSet::const_iterator vi=(mi->second).begin();
          vi!=(mi->second).end();++vi)
        var_cluster.tryMakeSetWithResize(vi->ident()) ;
    }
    for(map<variable,double>::const_iterator
          mi=info.begin();mi!=info.end();++mi) {
      int vid = (mi->first).ident() ;
      variableSet rv = get_all_recur_vars(s2t,mi->first) ;
      rv += get_all_recur_vars(t2s,mi->first) ;
      for(variableSet::const_iterator vi=rv.begin();
          vi!=rv.end();++vi)
        var_cluster.unionSets(vid,vi->ident()) ;
    }
    // fill in the alloc_info map
    for(map<variable,double>::const_iterator
          mi=info.begin();mi!=info.end();++mi) {
      int vid = (mi->first).ident() ;
      int rep = var_cluster.findSet(vid) ;
      alloc_info[variable(rep)] = mi->second ;
    }
  }

  variable MemGreedyScheduler::
  get_rep(const variable& v) {
    int rep = var_cluster.findSet(v.ident()) ;
    if(rep == -1) {
      ostringstream oss ;
      oss << "MemGreedyScheduler::get_rep(" << v << ") failed!" ;
      throw Loci::StringError(oss.str()) ;
    }
    return variable(rep) ;
  }

  variableSet MemGreedyScheduler::
  get_rep(const variableSet& vs) {
    variableSet ret ;
    for(variableSet::const_iterator vi=vs.begin();vi!=vs.end();++vi) {
      ret += get_rep(*vi) ;
    }
    return ret ;
  }

  double MemGreedyScheduler::
  get_alloc_size(const variable& v) {
    map<variable,double>::const_iterator mi = alloc_info.find(v) ;
    if(mi == alloc_info.end()) {
      ostringstream oss ;
      oss << "MemGreedyScheduler::get_alloc_size(" << v << ") failed!" ;
      throw Loci::StringError(oss.str()) ;
    }
    return mi->second ;
  }
  
  void MemGreedyScheduler::
  rule_effects(const rule& r, const variableSet& c,
               vector<pair<variable,double> >& e) {
    variableSet ts = get_rep(r.targets()) ;
    if(r.get_info().qualifier() == "DELETE") {
      for(variableSet::const_iterator vi=ts.begin();vi!=ts.end();++vi) {
        e.push_back(make_pair(*vi,-(get_alloc_size(*vi)))) ;
      }
    } else if(is_internal_rule(r)) {
      for(variableSet::const_iterator vi=ts.begin();vi!=ts.end();++vi) {
        e.push_back(make_pair(*vi,0)) ;
      }
    } else {
      for(variableSet::const_iterator vi=ts.begin();vi!=ts.end();++vi) {
        if(c.inSet(*vi)) {
          // already allocated, set size zero
          e.push_back(make_pair(*vi,0)) ;
        } else {
          e.push_back(make_pair(*vi,get_alloc_size(*vi))) ;
        }
      }
    }
  }

  std::vector<digraph::vertexSet>
  MemGreedyScheduler::schedule(const digraph& dag) {
    digraph::vertexSet dagv = dag.get_all_vertices() ;
    // schedule the graph without allocate rules first
    digraph dag_noalloc = allocFreeGr(dag) ;
    digraph dagt_noalloc = dag_noalloc.transpose() ;
    // and also simplify it to remove all the variables
    digraph g = simplify_graph(dag_noalloc) ;
    digraph::vertexSet gv = g.get_all_vertices() ;
    digraph gt = g.transpose() ;
    // vertices that have been scheduled
    digraph::vertexSet visited ;
    // variables in the graph that has been created by rules
    variableSet created_vars ;
    // candidate list to be scheduled next
    digraph::vertexSet candidates ;

    for(digraph::vertexSet::const_iterator
          vi=gv.begin();vi!=gv.end();++vi) {
      if(gt[*vi]==EMPTY)
        candidates += *vi ;
    }

    const double alpha = 0.1 ;
    vector<digraph::vertexSet> sched ;
    while(candidates != EMPTY) {
      // choose a candidate
      map<int,vector<pair<variable,double> > > effects ;
      map<int,double> score ;
      for(digraph::vertexSet::const_iterator
            vi=candidates.begin();vi!=candidates.end();++vi) {
        rule_effects(rule(*vi), created_vars, effects[*vi]) ;
        const vector<pair<variable,double> >& ve = effects[*vi] ;
        double s = 0 ;
        for(size_t i=0;i!=ve.size();++i)
          s += ve[i].second ;
        score[*vi] = s ;
      }
      double smax = -std::numeric_limits<double>::max() ;
      double smin = std::numeric_limits<double>::max() ;
      for(map<int,double>::const_iterator
            mi=score.begin();mi!=score.end();++mi) {
        if(mi->second > smax)
          smax = mi->second ;
        if(mi->second < smin)
          smin = mi->second ;
      }
      // creates RCL
      vector<int> rcl ;
      for(map<int,double>::const_iterator
            mi=score.begin();mi!=score.end();++mi) {
        double s = mi->second ;
        if(s <= smin+alpha*(smax-smin))
          rcl.push_back(mi->first) ;
      }

      // randomly pick one candidate from "rcl"
      // random_shuffle doesn't seem to generate
      // a unique random sequence in parallel ...
      // thus we are in trouble using it since that
      // creates different schedule on each process ...
      //random_shuffle(rcl.begin(), rcl.end()) ;
      //int c = rcl[0] ;
      // we'll just use drand48 with the same seed on each process
      int c = rcl[gen_rand(0,rcl.size()-1)] ;
      rule r(c) ;
      visited += c ;
      // see if rule "c"'s sources and targets are scheduled
      digraph::vertexSet sources = dagt_noalloc[c] ;
      digraph::vertexSet targets = dag_noalloc[c] ;
      // limit the sources and targets to be variables
      sources = get_vertexSet(extract_vars(sources)) ;
      targets = get_vertexSet(extract_vars(targets)) ;
      // remove scheduled ones
      sources -= visited ;
      targets -= visited ;
      // remove targets whose rules are not all scheduled
      digraph::vertexSet targets_not_ready ;
      for(digraph::vertexSet::const_iterator vi=targets.begin();
          vi!=targets.end();++vi)
        if(dagt_noalloc[*vi]-visited != EMPTY)
          targets_not_ready += *vi ;
      targets -= targets_not_ready ;

      visited += sources ;
      visited += targets ;
      // add to schedule
      if(sources != EMPTY)
        sched.push_back(sources) ;
      digraph::vertexSet main ;
      main += c ;
      sched.push_back(main) ;
      if(targets != EMPTY)
        sched.push_back(targets) ;
      // now add more candidates
      digraph::vertexSet nxt = g[c] ;
      candidates += get_valid_sched(nxt,gt,visited) ;
      candidates -= c ;
      // modify created_vars
      const vector<pair<variable,double> >& e = effects[c] ;
      for(size_t i=0;i!=e.size();++i) {
        if(e[i].second > 0)
          created_vars += e[i].first ;
        else if(e[i].second < 0)
          created_vars -= e[i].first ;
      }
    } // end while(candidates)
    
    // check for topology consistency
    vector<digraph::vertexSet> final = insertAlloc2Sched(dag,sched) ;
    //#define CHECK_TOPO
#ifdef CHECK_TOPO
    if(!final.empty()) {
      if(Loci::MPI_rank == 0)
        cout << "Checking graph schedule consistency..." ;
      // first check if all vertices are scheduled and
      // if we included any vertices that are not in the graph
      digraph::vertexSet allv ;
      for(size_t i=0;i!=final.size();++i)
        allv += final[i] ;
      if(dagv - allv != EMPTY) {
        cerr << endl
             << (dagv-allv).size() << " vertices not scheduled" << endl ;
        Loci::Abort() ;
      }
      if(allv-dagv != EMPTY) {
        cerr << endl
             << (allv-dagv).size() << " extra vertices included" << endl ;
        Loci::Abort() ;
      }
      // then check partial order
      digraph dagt = dag.transpose() ;
      digraph::vertexSet initial = final[0] ;
      for(digraph::vertexSet::const_iterator vi=initial.begin();
          vi!=initial.end();++vi)
        if(dagt[*vi] != EMPTY) {
          cerr << endl
               << "initial vertex error: " << *vi << endl ;
          if(*vi<0)
            cerr << "*vi = " << rule(*vi) << endl ;
          else
            cerr << "*vi = " << variable(*vi) << endl ;
          cerr << "depends on: " << extract_vars(dagt[*vi])
               << " and " << endl
               << extract_rules(dagt[*vi]) << endl ;
          Loci::Abort() ;
        }
      for(size_t i=1;i!=final.size();++i) {
        digraph::vertexSet cur = final[i] ;
        digraph::vertexSet pre ;
        for(size_t j=0;j!=i;++j)
          pre += final[j] ; 
        for(digraph::vertexSet::const_iterator ci=cur.begin();
            ci!=cur.end();++ci) {
          digraph::vertexSet depend = get_reachable(dagt,*ci) ;
          if(depend-pre != EMPTY) {
            cerr << endl << "scheduling error: " << *ci
                 << " violates its dependency, "
                 << (depend-pre) << " not yet scheduled!" << endl ;
            cerr << "vertex " << *ci << " = " ;
            if(*ci<0)
              cerr << rule(*ci) << endl ;
            else
              cerr << variable(*ci) << endl ;
            cerr << "depends on: " << extract_vars(depend-pre)
                 << " and " << endl
                 << extract_rules(depend-pre) << endl ;
            Loci::Abort() ;
          }
        }
      }
      if(Loci::MPI_rank == 0)
        cout << " passed!" << endl ;
    }
#endif
    return final ;
  }

  digraph MemGreedyScheduler::
  simplify_graph(const digraph& g) {
    // build a transposed g
    digraph gtr = g.transpose() ;
    // first copy the input graph
    digraph tmp_g = g ;
    // then strip all of the variables inside while keeping the dependency
    digraph::vertexSet vs = tmp_g.get_all_vertices() ;
    variableSet vars = extract_vars(vs) ;
    digraph::vertexSet varsvs ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      varsvs += vi->ident() ;
      // get the rules connected
      digraph::vertexSet pre = gtr[vi->ident()] ;
      digraph::vertexSet nxt = g[vi->ident()] ;
      for(digraph::vertexSet::const_iterator dvi=pre.begin();
          dvi!=pre.end();++dvi)
        tmp_g.add_edges(*dvi, nxt) ;
    }
    // remove all variables from the graph
    tmp_g = tmp_g.subgraph(vs-varsvs) ;

#ifdef COLLAPSE_LEAF
    // now we will collapse all the leaf nodes (inputs and outputs)
    gtr = tmp_g.transpose() ;
    int start = tmp_g.max_vertex() ; ++start ;
    vs = tmp_g.get_all_vertices() ;
    digraph::vertexSet visited ;
    digraph::vertexSet rm ;
    for(digraph::vertexSet::const_iterator
          dvi=vs.begin();dvi!=vs.end();++dvi) {
      digraph::vertexSet collapse ;
      // first check to collapse the targets of each vertex
      digraph::vertexSet nxt = tmp_g[*dvi] ;
      for(digraph::vertexSet::const_iterator
            dvi2=nxt.begin();dvi2!=nxt.end();++dvi2)
        if(tmp_g[*dvi2]==EMPTY) // leaf vertex
          collapse += *dvi2 ;

      collapse -= visited ;
      if(collapse != EMPTY) {
        int v = start++ ;
        digraph::vertexSet reachable ;
        for(digraph::vertexSet::const_iterator
              dvi2=collapse.begin();dvi2!=collapse.end();++dvi2)
          reachable += gtr[*dvi2] ;
        for(digraph::vertexSet::const_iterator
              dvi2=reachable.begin();dvi2!=reachable.end();++dvi2)
          tmp_g.add_edge(*dvi2, v) ;
        
        rm += collapse ;
        visited += collapse ;
      }

      // then check to collapse the sources
      collapse = EMPTY ;
      nxt = gtr[*dvi] ;
      for(digraph::vertexSet::const_iterator
            dvi2=nxt.begin();dvi2!=nxt.end();++dvi2)
        if(gtr[*dvi2]==EMPTY)
          collapse += *dvi2 ;

      collapse -= visited ;
      if(collapse != EMPTY) {
        int v = start++ ;
        digraph::vertexSet reachable ;
        for(digraph::vertexSet::const_iterator
              dvi2=collapse.begin();dvi2!=collapse.end();++dvi2)
          reachable += tmp_g[*dvi2] ;
        tmp_g.add_edges(v, reachable) ;
        rm += collapse ;
        visited += collapse ;
      }
    }
    // remove all the collapsed leaves
    tmp_g = tmp_g.subgraph(vs-rm) ;
#endif

    return tmp_g ;
  }

 
  


 //#define CLUSTER_TEST
  void MemGreedyScheduler::visit(loop_compiler& lc) {
    lc.collapse_sched = schedule(lc.collapse_gr) ;
    lc.advance_sched = schedule(lc.advance_gr) ;
  }
  
  void MemGreedyScheduler::visit(dag_compiler& dc) {
    dc.dag_sched = schedule(dc.dag_gr) ;
  }

 void MemGreedyScheduler::visit(conditional_compiler& cc) {
    cc.dag_sched = schedule(cc.cond_gr) ;
  }
  






 

} // end of namespace Loci
