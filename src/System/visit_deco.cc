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
#include "visitor.h"
#include "visit_tools.h"
#include "comp_tools.h"


using std::vector ;
using std::map ;
using std::set ;

using std::cerr ;
using std::endl ;

namespace Loci {
  extern bool profile_memory_usage ;

  ////////////////////////////////////////////////////////////////
  // allocInfoVisitor
  ////////////////////////////////////////////////////////////////
  namespace {
    // working is a set of recurrence variables. This function
    // filters the variables that do not need to be allocated.
    // i.e., this function returns the recurrence variables in
    // the "working" set that do NOT need to be allocated.
    variableSet filte_recur_nonalloc(const variableSet& working,
                                     const variableSet& testing,
                                     const map<variable,variableSet>& s2t,
                                     const map<variable,variableSet>& t2s
                                     ) {
      variableSet ret ;
      variableSet processed ;
      for(variableSet::const_iterator vi=working.begin();
          vi!=working.end();++vi) {
        if(processed.inSet(*vi))
          continue ;
        // first get all the source variables from this one
        variableSet sources = get_all_recur_vars(t2s,*vi) ;
        if( (sources == EMPTY) ||
            (variableSet(sources&testing) == EMPTY)
            ) {
          // then this variable is chosen to allocate
          variableSet targets = get_all_recur_vars(s2t,*vi) ;
          ret += sources ;
          ret += targets ;
          processed += *vi ;
          processed += sources ;
          processed += targets ;
        }else {
          ret += *vi ;
          processed += *vi ;
        }
      } // end of for(vi)

      return ret ;
    }

  } // end of unnamed namespace
  
  variableSet allocInfoVisitor::get_start_info(const digraph& gr,int id) {
    variableSet working_vars ;

    // obtain all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;
    
    // looping over all the rules and gather variables
    // need to be processed in this graph for allocation
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      // get the targets of the rule
      working_vars += ruleIter->targets() ;
    }
      
    // we remove recurrence target variables from the working_vars.
    //working_vars -= recur_target_vars ;

    // the following code attempts to fix the allocation
    // problem related to the priority variables. for example:
    // heat::u_f{n,it} vs. u_f{n,it}
    // For any recurrence variable, we allocate the one
    // we first see, and once the first one gets allocated,
    // we mark any other related recurrence variables not
    // to be allocated.
    working_vars -= allocated_vars ;
    variableSet takeoff =
      filte_recur_nonalloc(variableSet(working_vars&all_recur_vars),
                           variableSet
                           (extract_vars
                            (gr.get_all_vertices()) & all_recur_vars),
                           recurs2t,recurt2s) ;
    allocated_vars += takeoff ;
    working_vars -= takeoff ;

    return working_vars ;
  }
  
  variableSet allocInfoVisitor::gather_info(const digraph& gr,
                                            const variableSet& working_vars) {
    // NOTE: even if the passed in working_vars is EMPTY, we cannot return
    // from this function yet, because we need to look for possible
    // loop supernodes in the current graph and if found, allocate its
    // rotation list!
    
    // obtain the transpose of the graph
    digraph grt = gr.transpose() ;

    // variables that allocated in this level (graph)
    variableSet alloc_here ;

    // looping over working_vars and gather allocation info
    for(variableSet::const_iterator varIter=working_vars.begin();
        varIter!=working_vars.end();++varIter) {
      // if this variable has been allocated before
      // we don't need to process it here
      if(allocated_vars.inSet(*varIter))
        continue ;      
      // get the rules that produce this variable
      ruleSet source_rules = extract_rules(grt[varIter->ident()]) ;
      // see if these rules are all super node
      bool all_super_nodes = true ;
      for(ruleSet::const_iterator ruleIter=source_rules.begin();
          ruleIter!=source_rules.end();++ruleIter) {
        if(!is_super_node(ruleIter)) {
          all_super_nodes = false ;
          break ;
        }
        else {
          // we eliminate the recursive super node
          // by look up the graph_sn set
          set<int>::const_iterator found ;
          found = graph_sn.find(get_supernode_num(*ruleIter)) ;
          if(found == graph_sn.end()) {
            all_super_nodes = false ;
            break ;
          }
        }
      }
      // if all the source rules are super node
      // there are two possible cases:
      
      // the first is there is only ONE source super node
      // then we don't allocate it in this level (graph)
      // the second is there are more than one super node
      // that produce this variable, then we have to allocate
      // this variable in this level (graph)
      if( (all_super_nodes) &&
          (source_rules.size() == 1)
          ) {
        if(all_recur_vars.inSet(*varIter)) {
          // if this variable is a recurrence variable,
          // then we push back all its disabled source
          // and target variables. This gives the chance
          // to reselect which one to allocate in the
          // following graphs
          variableSet sources = get_all_recur_vars(recurt2s,*varIter) ;
          variableSet targets = get_all_recur_vars(recurs2t,*varIter) ;
          allocated_vars -= sources ;
          allocated_vars -= targets ;
        }
        continue ;
      }
      // else we allocate it in this level
      // add the variable into the alloc_here set
      alloc_here += *varIter ;
    }
    // we check if there is any loop super node in the graph,
    // if yes, we also allocate the rotate list variables
    // for the loop in this graph.

    // obtain all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    variableSet all_rot_vars ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(is_super_node(ruleIter)) {
        int id = get_supernode_num(*ruleIter) ;
        if(inSet(loop_sn,id)) {
          // get the rotate list variables
          map<int,variableSet>::const_iterator found ;
          found = rotate_vtable.find(id) ;
          FATAL(found == rotate_vtable.end()) ;
          all_rot_vars += found->second ;
        }
      }
    }
    all_rot_vars -= allocated_vars ;
    // if any rotate list variable belongs to
    // recurrence variables cluster, we allocate
    // the first one of them
    variableSet takeoff =
      filte_recur_nonalloc(variableSet(all_rot_vars&all_recur_vars),
                           variableSet
                           (extract_vars
                            (gr.get_all_vertices()) & all_recur_vars),
                           recurs2t,recurt2s) ;
    allocated_vars += takeoff ;
    all_rot_vars -= takeoff ;
    
    alloc_here += all_rot_vars ;
    
    return alloc_here ;
  }

  void allocInfoVisitor::visit(loop_compiler& lc) {
    // we first get the loop shared variables
    // we need to instruct the collpase part allocate
    // them, but the advance part don't need to allocate
    // them. the advance part only deletes them. and upon
    // exit of the loop, we make an additional deletion
    // of these shared variables in the conditional node
    // of this loop compiler
    variableSet shared ;
    map<int,variableSet>::const_iterator found ;
    found = loop_shared_table.find(lc.cid) ;
    FATAL(found == loop_shared_table.end()) ;
    shared += found->second ;

    variableSet tmp,loop_alloc_vars, working_vars ;

    // suppose the id of loop_compiler is x,
    // then -x refers to the collapse part
    // x refers to the advance part
    working_vars = get_start_info(lc.collapse_gr,-lc.cid) ;
    tmp=gather_info(lc.collapse_gr,working_vars) ;
    // edit corresponding variables
    allocated_vars += tmp ;
    alloc_table[-lc.cid] += tmp ;
    loop_alloc_vars += tmp ;

    // the advance part
    working_vars = get_start_info(lc.advance_gr,lc.cid) ;
    working_vars -= shared ;
    tmp=gather_info(lc.advance_gr,working_vars) ;

    // edit corresponding variables
    allocated_vars += tmp ;
    alloc_table[lc.cid] += tmp ;
    loop_alloc_vars += tmp ;

    loop_alloc_table[lc.cid]=loop_alloc_vars ;
  }

  void allocInfoVisitor::visit(dag_compiler& dc) {
    variableSet working_vars = get_start_info(dc.dag_gr,dc.cid) ;
    variableSet ret = gather_info(dc.dag_gr,working_vars) ;

    allocated_vars += ret ;
    alloc_table[dc.cid] += ret ;
  }

  void allocInfoVisitor::visit(conditional_compiler& cc) {
    variableSet working_vars = get_start_info(cc.cond_gr,cc.cid) ;
    variableSet ret = gather_info(cc.cond_gr,working_vars) ;

    allocated_vars += ret ;
    alloc_table[cc.cid] += ret ;
  }

  /////////////////////////////////////////////////////////
  // allocGraphVisitor
  /////////////////////////////////////////////////////////

  // we decorate the graph to include the allocation rules
  void allocGraphVisitor::edit_gr(digraph& gr,rulecomp_map& rcm,int id) {
    // first get the transposed graph
    digraph grt = gr.transpose() ;
    // obtain all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    // variables processed in this graph
    variableSet working_vars ;

    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      working_vars += ruleIter->targets() ;
    }

    // first we get the info in alloc_table, find out which variables
    // need to be allocated in the graph
    map<int,variableSet>::const_iterator found ;
    found = alloc_table.find(id) ;
    FATAL(found == alloc_table.end()) ;
    variableSet alloc_vars = found->second ;

    if(alloc_vars == EMPTY) return ;
    
    for(variableSet::const_iterator varIter=working_vars.begin();
        varIter!=working_vars.end();++varIter) {
      // we check if the variable is in the alloc_vars
      // if it is not in the alloc_vars, we don't allocate it here
      if(!alloc_vars.inSet(*varIter))
        continue ;
      // now we decorate the graph to include the allocation
      // we first obtain the ident of rules that produce the variable
      digraph::vertexSet source_rules = grt[varIter->ident()] ;
      // if source_rules is empty, then the variable is not
      // present in this graph, we need a search to find out
      // the source_rules
      if(source_rules == EMPTY) {
        for(ruleSet::const_iterator ruleIter=rules.begin();
            ruleIter!=rules.end();++ruleIter) {
          variableSet targets = ruleIter->targets() ;
          if(targets.inSet(*varIter))
            source_rules += ruleIter->ident() ;
        }
      }
      // if a variable is a priority source variable,
      // then we need to make the allocation rule
      // point to all the rules that generate other priority
      // target variables to make sure that the variable
      // is allocated before any of the priority variables
      // get computed
      if(prio_sources.inSet(*varIter)) {
        variableSet prio_targets = get_all_recur_vars(prio_s2t,*varIter) ;
        // variable set to mark variables not in the current graph
        // and needed to perform a complete search of rules later.
        // see comments below.
        variableSet empty_psrules_vars = variableSet(EMPTY) ;
        for(variableSet::const_iterator pvi=prio_targets.begin();
            pvi!=prio_targets.end();++pvi) {
          digraph::vertexSet psrules = grt[pvi->ident()] ;
          // however if the psrules (source rules) for this
          // priority variable is empty, then the variable
          // is not in the present graph, but it is still possible
          // that we have rules in the graph that compute it
          // we therefore need to perform a search on all rule
          // targets to see if it actually get generated in
          // the current graph. We mark this variable for
          // such a search later.
          ruleSet rs = extract_rules(psrules) ;
          if(rs == EMPTY) {
            empty_psrules_vars += *pvi ;
          }
          source_rules += get_vertexSet(rs) ;
        } // end for(pvi)
        // perform a complete search if necessary
        if(empty_psrules_vars != EMPTY) {
          for(ruleSet::const_iterator ruleIter=rules.begin();
              ruleIter!=rules.end();++ruleIter) {
            variableSet targets = ruleIter->targets() ;
            if( (targets & empty_psrules_vars) != EMPTY)
              source_rules += ruleIter->ident() ;
          }
        }
      }
      
      // now we create a rule for allocation
      variable sv("CREATE") ;
      rule alloc_rule = create_rule(sv,*varIter,"ALLOCATE") ;
      // edit graph
      gr.add_edges(alloc_rule.ident(),source_rules) ;
      // modify the rulecomp_map
      rcm[alloc_rule] = new allocate_var_compiler(*varIter) ;

      // remove variable in alloc_vars
      alloc_vars -= *varIter ;
    }
    // if alloc_vars is not empty, it must contain loop rotate list
    // variables. we allocate them at here
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      variableSet rotate_vars ;
      variableSet allocated ;
      
      if(is_super_node(ruleIter)) {
        int id = get_supernode_num(*ruleIter) ;
        if(inSet(loop_sn,id)) {
          // get the rotate list variables
          map<int,variableSet>::const_iterator found ;
          found = rotate_vtable.find(id) ;
          FATAL(found == rotate_vtable.end()) ;
          rotate_vars += found->second ;

          for(variableSet::const_iterator vi=alloc_vars.begin();
              vi!=alloc_vars.end();++vi) {
            // rotate list vars
            if(rotate_vars.inSet(*vi)) {
              allocated += *vi ;
              
              // we create a rule for allocation
              variable sv("CREATE") ;
              rule alloc_rule = create_rule(sv,*vi,"ALLOCATE") ;
              // edit graph
              gr.add_edge(alloc_rule.ident(),ruleIter->ident()) ;
              // modify the rulecomp_map
              rcm[alloc_rule] = new allocate_var_compiler(*vi) ;
            }
          }
          alloc_vars -= allocated ;
        }
      }
    }
    
    // if alloc_vars is not empty at this stage, it means
    // there still have some variables need to be allocated in this
    // graph but these variables do not belong to this graph

    // generally it is not possible to have such case at this stage, we
    // leave it for future development. We make a warning here
    if(alloc_vars != EMPTY) {
      cerr << "Warning: allocation not complete in graph: "
           << id << endl
           << "Remaining variables are: " << alloc_vars << endl ;
    }
  }

  void allocGraphVisitor::visit(loop_compiler& lc) {
    edit_gr(lc.collapse_gr,lc.rule_compiler_map,-lc.cid) ;
    edit_gr(lc.advance_gr,lc.rule_compiler_map,lc.cid) ;
  }

  void allocGraphVisitor::visit(dag_compiler& dc) {
    edit_gr(dc.dag_gr,dc.rule_compiler_map,dc.cid) ;
  }

  void allocGraphVisitor::visit(conditional_compiler& cc) {
    edit_gr(cc.cond_gr,cc.rule_compiler_map,cc.cid) ;
  }

  //////////////////////////////////////////////////////////
  // deleteInfoVisitor
  //////////////////////////////////////////////////////////
  deleteInfoVisitor::
  deleteInfoVisitor(const map<int,variableSet>& lat,
                    const map<variable,variableSet>& rvt2s,
                    const map<variable,variableSet>& rvs2t,
                    const set<int>& gsn,
                    const map<int,int>& pnt,
                    const map<int,int>& lct,
                    const set<int>& lsn,
                    const set<int>& csn,
                    const map<int,variableSet>& rot_vt,
                    const map<int,variableSet>& lsharedt,
                    const variableSet& promoted_rep,
                    const variableSet& reserved_vars)
    :recur_vars_t2s(rvt2s), recur_vars_s2t(rvs2t),
     loop_alloc_table(lat), graph_sn(gsn),
     pnode_table(pnt), loop_ctable(lct),
     loop_sn(lsn), cond_sn(csn), rotate_vtable(rot_vt),
     loop_shared_table(lsharedt), promoted_rep(promoted_rep) {
      
    for(std::map<int,variableSet>::const_iterator mi=rotate_vtable.begin();
        mi!=rotate_vtable.end();++mi)
      all_rot_vars += mi->second ;

    deleted_vars += reserved_vars ;
  }

  variableSet deleteInfoVisitor::get_start_info(const digraph& gr,int id) {

    variableSet working_vars ;

    // obtain all the possible variables in the graph
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      working_vars += ri->sources() ;
    }
    
    working_vars += extract_vars(allvertices) ;

    // remove any deleted or reserved variables
    working_vars -= deleted_vars ;

    // get the interface variable
    // flip id if it is negative
    if(id < 0) id = -id ;
    map<int,variableSet>::const_iterator mf ;
    variableSet vars_by_others ;
    for(mf=sn_del_interface.begin();mf!=sn_del_interface.end();++mf) {
      if(mf->first != id)
        vars_by_others += mf->second ;
    }

    working_vars -= vars_by_others ;
    
    return working_vars ;    
  }

  variableSet deleteInfoVisitor::gather_info(const digraph& gr,
                                             const variableSet& working_vars) {
    if(working_vars == EMPTY)
      return variableSet(EMPTY) ;
    // get all vertices in the graph
    digraph::vertexSet allv = gr.get_all_vertices() ;
    // variables that been deleted in this graph
    variableSet delete_here ;
    // looping over all the working variables to gather info
    for(variableSet::const_iterator varIter=working_vars.begin();
        varIter!=working_vars.end();++varIter) {

      // some recurrence variable processing
      // get all the recurrence source variables of *varIter
      variableSet all_recur_sources ;

      all_recur_sources = get_all_recur_vars(recur_vars_t2s,*varIter) ;

      // get the rules that this variable reaches
      ruleSet target_rules = extract_rules(gr[varIter->ident()]) ;
      // if this variable is in the recurrence variable table,
      // we'll also need to get the rules that the recurrence variable
      // reach

      if(all_recur_sources != EMPTY) {
        // get the recurrence source variables that are in this graph
        digraph::vertexSet in_graph_recur_vars ;
        for(variableSet::const_iterator vi=all_recur_sources.begin();
            vi!=all_recur_sources.end();++vi) {
          if(allv.inSet(vi->ident()))
            in_graph_recur_vars += vi->ident() ;
        }

        if(in_graph_recur_vars != EMPTY) {
          ruleSet others ;
          for(digraph::vertexSet::const_iterator
                vi=in_graph_recur_vars.begin();
              vi!=in_graph_recur_vars.end();
              ++vi)
              others = extract_rules(gr[*vi]) ;
          // but we must exclude those "generalize", "promote" and
          // "priority" internal rules from others
          variableSet only_vars = extract_vars(in_graph_recur_vars) ;
          only_vars += *varIter ;
          ruleSet subtract ;
          for(ruleSet::const_iterator ruleIter=others.begin();
              ruleIter!=others.end();++ruleIter) {
            if( (ruleIter->get_info().qualifier() == "promote") ||
                (ruleIter->get_info().qualifier() == "generalize") ||
                (ruleIter->get_info().qualifier() == "priority")
                )
              subtract += *ruleIter ;

            //variable target = *(ruleIter->targets().begin()) ;
            //if(only_vars.inSet(target))
            //subtract += *ruleIter ;
          }

          others -= subtract ;
          if(others != EMPTY) {
            target_rules += others ;
            recur_source_other_rules[*varIter] = others ;
          }
        }
      }

      // check if this variable is a current loop
      // shared variable, if it is, then check if it is a promoted
      // variable, if it is, then we don't need to look at
      // other things, just delete it here
      if(current_lshared_vars.inSet(*varIter)) {
        if(promoted_rep.inSet(*varIter)) {
          delete_here += *varIter ;
          //cout << "hacked variable: " << *varIter << endl ;
	  //cout << "working vars: " << working_vars << endl ;
          continue ;
        }
      }

      // see if we can delete the variable in this graph
      if(let_it_go(gr,target_rules,*varIter)) {
        continue ;
      }

      // otherwise we can delete the variable here
      delete_here += *varIter ;
    }

    // return the delete_here
    return delete_here ;    
  }

  void deleteInfoVisitor::visit(loop_compiler& lc) {
    // first in the loop, we don't delete the time_level variable
    variable tvar = variable(lc.tlevel) ;
    deleted_vars += tvar ;

    // we then get the variables shared between the
    // advance and collapse part of the loop
    map<int,variableSet>::const_iterator found ;
    variableSet shared ;
    found = loop_shared_table.find(lc.cid) ;
    FATAL(found == loop_shared_table.end()) ;
    shared = found->second ;
    current_lshared_vars = shared ;

    // gather variables in the rotate list
    variableSet rotate_vars ;

    found = rotate_vtable.find(lc.cid) ;
    FATAL(found == rotate_vtable.end()) ;
    rotate_vars = found->second ;

    // we then gather deletion information for the collapse part
    variableSet working_vars = get_start_info(lc.collapse_gr,-lc.cid) ;
    // we delete the shared variables until the last iteration
    // finishes and inside the conditional super node
    working_vars -= shared ;
    // we don't delete variabls in rotate lists, we defer that until
    // the last iteration finishes
    working_vars -= rotate_vars ;

    looping_algr(working_vars,lc.collapse_gr,-lc.cid,0) ;
    
    // we gather deletion information for the advance part
    working_vars = get_start_info(lc.advance_gr,lc.cid) ;
    working_vars -= rotate_vars ;
    // we don't need to exclude those shared vars, but we
    // will need to count which shared variable gets
    // deleted in the advance graph, because we want to
    // delete them again later upon exit of the loop
    variableSet before = variableSet(shared - deleted_vars) ;
    looping_algr(working_vars,lc.advance_gr,lc.cid,0) ;
    variableSet after = variableSet(shared - deleted_vars) ;

    // reset to empty
    current_lshared_vars = variableSet(EMPTY) ;

    // we now schedule the deletion of shared and rotate_vars

    // first get the collapse node of this loop
    map<int,int>::const_iterator found2 ;
    found2 = loop_ctable.find(lc.cid) ;
    FATAL(found2 == loop_ctable.end()) ;
    int collapse_id = found2->second ;

    working_vars = rotate_vars ;
    working_vars -= deleted_vars ;
    // we need to add those deleted shared vars
    working_vars += variableSet(before - after) ;
    looping_algr(working_vars,lc.loop_gr,collapse_id,2) ;

  }

  void deleteInfoVisitor::visit(dag_compiler& dc) {
    variableSet working_vars = get_start_info(dc.dag_gr,dc.cid) ;
    looping_algr(working_vars,dc.dag_gr,dc.cid,1) ;    
  }

  
  void deleteInfoVisitor::visit(conditional_compiler& cc) {
    variableSet working_vars = get_start_info(cc.cond_gr,cc.cid) ;
    looping_algr(working_vars,cc.cond_gr,cc.cid,1) ;
  }

  // return the found loop node id, -1 if not found
  // loop_num specifies how many loops to find and include
  // along in the parent node path
  int deleteInfoVisitor::only_loop_alloc(variableSet& working_vars,
                                         int id, int loop_num) {
    // the top level 
    if(id == 0)
      return -1 ;

    if(loop_num < 0) {
      cerr << "ERROR calling only_loop_alloc() with loop_num < 0" << endl ;
      Loci::Abort() ;
    }
    // we check if we need to flip the id, if we are processing
    // the collapse graph of the loop
    if(id < 0) {
      if(inSet(loop_sn,-id))
        id = -id ;
      else {
        cerr << "ERROR calling only_loop_alloc() function,"
             << "with negative id number, only loop compiler "
             << "is allowed to have negative id number for its"
             << "collapse graph, the caller(id): " << id
             << " is not a loop compiler" << endl ;
        Loci::Abort() ;
      }
    }

    // if loop_num == 0, the only possible pass of such argument
    // to this function is inside the loop compiler, we do a check
    if(loop_num == 0) {
      if(!inSet(loop_sn,id)) {
        cerr << "ERROR calling only_loop_alloc() function,"
             << "caller(id): " << id << " is not a loop compiler"
             << " and therefore cannot use loop_num with 0!" << endl ;
        Loci::Abort() ;
      }
    }

    map<int,int>::const_iterator found ;
    map<int,variableSet>::const_iterator found_alloc ;
    int ret_id = -1 ;

    int pnode = id ;
    int found_loop_num = 0 ;

    // we search for loop_num loops
    while( (pnode != 0) && (found_loop_num != loop_num)) {
      found = pnode_table.find(pnode) ;
      FATAL(found == pnode_table.end()) ;
      pnode = found->second ;
      
      if(inSet(loop_sn,pnode))
        ++found_loop_num ;
    }
    
    if(found_loop_num == loop_num) {
      // get the loop's allocation
      found_alloc = loop_alloc_table.find(pnode) ;
      FATAL(found_alloc == loop_alloc_table.end()) ;
      
      working_vars &= found_alloc->second ;
      
      ret_id = pnode ;
    }
    
    return ret_id ;
  }

  // the looping algorithm for schedule deletion of variables
  // in the multilevel graph

  // working_vars is the initial unconstrained variableSet to
  // be evaluated, dg is the graph of a particular compiler
  // id is the compiler's id, start_loop_num is the starting
  // number for loop allocation constrain (0 for loop compiler,
  // 1 for dag & conditional compiler).
  void deleteInfoVisitor::looping_algr(const variableSet& working_vars,
                                       const digraph& dg,int id,
                                       int start_loop_num) {
    // it is important to use "+=" here, because of the loop,
    // we may have something in the delete_table[id] already
    // use "=" only will erase the previous contents.

    // *************************************************** //
    // IN ALL CODE HERE, WHENEVER WE WANT TO USE delete_table
    // AS LEFT VALUE, WE SHOULD ALWAYS USE "+=".
    // *************************************************** //
    
    // we can also skip this checking, it is just slightly
    // faster to do it here
    if(working_vars == EMPTY) {
      delete_table[id] += variableSet(EMPTY) ;
      return ;
    }
    // the looping algorithm
    int loop_num = start_loop_num ;
    int found_loop_id, last_loop_id ;
    variableSet processed_vars = working_vars ;
    variableSet ret ;
    
    found_loop_id = only_loop_alloc(processed_vars,id,loop_num) ;
    ret = gather_info(dg,processed_vars) ;
    last_loop_id = found_loop_id ;

    // edit table
    deleted_vars += ret ;
    delete_table[id] += ret ;

    variableSet remaining_vars = variableSet(working_vars - processed_vars) ;

    while(remaining_vars != EMPTY) {
      last_loop_id = found_loop_id ;
      ++loop_num ;
      found_loop_id = only_loop_alloc(remaining_vars,id,loop_num) ;
      ret = gather_info(dg,remaining_vars) ;

      if(last_loop_id != -1) {
        // put the ret into the collapse node of last_loop_id loop
        if(ret != EMPTY) {
          // we then find out the collapse super node inside this loop node
          map<int,int>::const_iterator found ;
          found = loop_ctable.find(last_loop_id) ;
          FATAL(found == loop_ctable.end()) ;
          int collapse_id = found->second ;
          
          deleted_vars += ret ;
          delete_table[collapse_id] += ret ;
        }
      }
      processed_vars += remaining_vars ;
      remaining_vars = variableSet(working_vars - processed_vars) ;
    }
  }
                                              
  // given a rule set and a graph
  // return true if there is only one super node and
  // all other rules have path to the super node, or if
  // it is only one super node.
  // return false otherwise
  bool deleteInfoVisitor::let_it_go(const digraph& gr, ruleSet rules,
                                    const variable& v) {
    if(rules == EMPTY) {
      // if the target rules are empty, we look back to the rules
      // that generate the variable, if it is only one super node
      // we let it go, meaning we don't delete the variable here
      // this only applies to non loop rotate variables
      digraph grt = gr.transpose() ;
      ruleSet source_rules = extract_rules(grt[v.ident()]) ;

      if(source_rules.size() == 1) {
        if(is_super_node(source_rules.begin())) {
          int sn_id = get_supernode_num(*source_rules.begin()) ;
          if(inSet(graph_sn,sn_id)) {
            if(!all_rot_vars.inSet(v)) {
              // add the escape info to the record
              sn_del_interface[sn_id] += v ;
              return true ;
            }
          }
        }
      }
      return false ;
    }
    
    int supernode_num = 0 ;
    int supernode_id = 1 ; //vertex number (id) for supernode
    rule supernode ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(is_super_node(ruleIter)) {
        // we eliminate the recursive super node
        // by looking up the graph_sn set
        if(inSet(graph_sn,get_supernode_num(*ruleIter))) {
          supernode = *ruleIter ;
          supernode_id = ruleIter->ident() ;
          ++supernode_num ;
        }
      }
    }

    int sn_id = -1 ;
    if(supernode_num == 1)
      sn_id = get_supernode_num(supernode) ;

    if(supernode_num == 1) {
      // check to see if it is a conditional node
      // if it is, we delete it in the graph
      if(inSet(cond_sn,sn_id))
        return false ;
      else if(rules.size() == 1) {
        // add the escape info.
        sn_del_interface[sn_id] += v ;
        return true ;
      }
    }
    if(supernode_num != 1)
      return false ;
    // the remaining case must be there are only one super node
    // and at least one non super node, we check if there is path
    // from each non super node to the super node in the graph
    FATAL(supernode_id > 0) ;
    rules -= rule(supernode_id) ;
    bool ret = true ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(!has_path(gr,ruleIter->ident(),supernode_id)) {
        ret = false ;
        break ;
      }
    }
    if(ret) {
      // add the escape info.
      sn_del_interface[sn_id] += v ;
    }
    return ret ;
  }

  //////////////////////////////////////////////////////////
  // deleteGraphVisitor
  //////////////////////////////////////////////////////////
  void deleteGraphVisitor::edit_gr(digraph& gr,rulecomp_map& rcm,int id) {
    // first we get the info in delete_table, find out which variables
    // need to be deleted in the graph
    map<int,variableSet>::const_iterator found ;
    found = delete_table.find(id) ;
    variableSet delete_vars = variableSet(EMPTY) ;
    if(found != delete_table.end())
      delete_vars = found->second ;

    if(delete_vars == EMPTY) return ;
    
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;
    
    // get all possible variables in the graph
    variableSet working_vars ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      working_vars += ri->sources() ;
    
    working_vars += extract_vars(allvertices) ;

    for(variableSet::const_iterator varIter=working_vars.begin();
        varIter!=working_vars.end();++varIter) {
      // we check if the variable is in the delete_vars
      // if it is not in the delete_vars, we don't delete it here
      if(!delete_vars.inSet(*varIter))
        continue ;
      // now we decorate the graph to include the deletion
      // we first obtain the ident of rules that are reachable
      // from the variable
      digraph::vertexSet target_rules = gr[varIter->ident()] ;

      // if the variable is a recurrence target variable,
      // we should also get possible rules that the recurrence
      // source variable leads to 
      map<variable,ruleSet>::const_iterator found ;
      found = recur_source_other_rules.find(*varIter) ;
      
      if(found != recur_source_other_rules.end()) {
        ruleSet others = found->second ;
        for(ruleSet::const_iterator ruleIter=others.begin();
            ruleIter!=others.end();++ruleIter)
          target_rules += ruleIter->ident() ;
      }

      // if target_rules is empty, then the variable is not
      // present in this graph, we need a search to find out
      // the target_rules
      if(target_rules == EMPTY) {
        for(ruleSet::const_iterator ruleIter=rules.begin();
            ruleIter!=rules.end();++ruleIter) {
          variableSet sources = ruleIter->sources() ;
          if(sources.inSet(*varIter))
            target_rules += ruleIter->ident() ;
        }
      }

      // if target_ruls is still empty, then the variable must
      // be in the graph, but it is a leaf node, we append the
      // deletion rule after the variable itself.
      if(target_rules == EMPTY)
        target_rules += varIter->ident() ;
       
      // we create a rule for deletion
      variable sv("DESTROY") ;
      rule delete_rule = create_rule(sv,*varIter,"DELETE") ;
      // edit graph
      gr.add_edges(target_rules,delete_rule.ident()) ;
      // modify the rulecomp_map
      rcm[delete_rule] = new free_var_compiler(*varIter) ;

      // add memory profiling alloc compiler here
      if(profile_memory_usage) {
        variable sv("PROFILE") ;
        rule profile_delete_rule = create_rule(sv,*varIter,"MEMPROFILED") ;
        gr.add_edges(target_rules,profile_delete_rule.ident()) ;
        // gather profile alloc compilers connect to the target rules
        digraph::vertexSet profileAcompilers ;
        for(digraph::vertexSet::const_iterator paci=target_rules.begin();
            paci!=target_rules.end();++paci) {
          ruleSet next = extract_rules(gr[*paci]) ;
          for(ruleSet::const_iterator pacii=next.begin() ;
              pacii!=next.end();++pacii)
            if(pacii->get_info().qualifier() == "MEMPROFILEA")
              profileAcompilers += pacii->ident() ;
        }
        gr.add_edges(profileAcompilers,profile_delete_rule.ident()) ;
        // this must happen before the actual deletion
        gr.add_edge(profile_delete_rule.ident(),delete_rule.ident()) ;
        rcm[profile_delete_rule] = new memProfileFree_compiler(*varIter) ;
      }
      
      // remove variable in alloc_vars
      delete_vars -= *varIter ;
    }
    // if delete_vars is not empty at this stage, it means
    // there must have some variables need to be deleted in this
    // graph but these variables do not belong to this graph

    if(delete_vars != EMPTY) {
      // first find out all the leaf node in the graph
      digraph::vertexSet leaf_nodes ;
      for(digraph::vertexSet::const_iterator vi=allvertices.begin();
          vi!=allvertices.end();++vi) {
        if(gr[*vi] == EMPTY)
          leaf_nodes += *vi ;
      }
      
      for(variableSet::const_iterator vi=delete_vars.begin();
          vi!=delete_vars.end();++vi) {
        // we create a rule for deletion
        variable sv("DESTROY") ;
        rule delete_rule = create_rule(sv,*vi,"DELETE") ;
        // we make the rule happen at last in the graph
        // we connect all the leaf node in the graph
        // to this deletion rule
        gr.add_edges(leaf_nodes,delete_rule.ident()) ;
        // modify the rulecomp_map
        rcm[delete_rule] = new free_var_compiler(*vi) ;

        // add memory profiling alloc compiler here
        if(profile_memory_usage) {
          variable sv("PROFILE") ;
          rule profile_delete_rule = create_rule(sv,*vi,"MEMPROFILED") ;
          gr.add_edges(leaf_nodes,profile_delete_rule.ident()) ;
          // this must happen before the actual deletion
          gr.add_edge(profile_delete_rule.ident(),delete_rule.ident()) ;
          rcm[profile_delete_rule] = new memProfileFree_compiler(*vi) ;
        }// end if profile_memory_usage
        
      }
    }
  }

  void deleteGraphVisitor::visit(loop_compiler& lc) {
    edit_gr(lc.collapse_gr,lc.rule_compiler_map,-lc.cid) ;
    edit_gr(lc.advance_gr,lc.rule_compiler_map,lc.cid) ;
  }
  
  void deleteGraphVisitor::visit(dag_compiler& dc) {
    edit_gr(dc.dag_gr,dc.rule_compiler_map,dc.cid) ;
  }
  
  void deleteGraphVisitor::visit(conditional_compiler& cc) {
    edit_gr(cc.cond_gr,cc.rule_compiler_map,cc.cid) ;
  }

  //////////////////////////////////////////////////////////
  // memProfileAllocDecoVisitor
  //////////////////////////////////////////////////////////
  memProfileAllocDecoVisitor::
  memProfileAllocDecoVisitor(const std::map<int,variableSet>& t,
                             const std::map<int,variableSet>& rot_vt)
    :alloc_table(t) {
    map<int,variableSet>::const_iterator mi ;
    for(mi=t.begin();mi!=t.end();++mi)
      all_alloc_vars += mi->second ;
    for(mi=rot_vt.begin();mi!=rot_vt.end();++mi)
      loop_rotation_vars += mi->second ;
  }

  void memProfileAllocDecoVisitor::
  edit_gr(digraph& gr,rulecomp_map& rcm,int id) {
    // first get the transposed graph
    digraph grt = gr.transpose() ;
    // obtain all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    // variables processed in this graph
    variableSet working_vars ;

    ruleSet alloc_profile_rules ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if( (ruleIter->get_info().qualifier() == "MEMPROFILEA") ||
          (ruleIter->get_info().qualifier() == "ALLOCATE")
          )
        alloc_profile_rules += *ruleIter ;
    }
    rules -= alloc_profile_rules ;

    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      working_vars += ruleIter->targets() ;
    }


    // first we get the info in alloc_table, find out which variables
    // need to be allocated in the graph
    map<int,variableSet>::const_iterator found ;
    found = alloc_table.find(id) ;
    FATAL(found == alloc_table.end()) ;
    variableSet profile_vars = found->second ;

    // rotation list vars are handled specially
    profile_vars -= loop_rotation_vars ;
    {
      variableSet local_loop_rotation_vars =
        variableSet(loop_rotation_vars & working_vars) ;
      variableSet local_alloc_lrt_vars =
        variableSet(local_loop_rotation_vars & all_alloc_vars) ;
      profile_vars += local_alloc_lrt_vars ;
    }
    for(variableSet::const_iterator varIter=working_vars.begin();
        varIter!=working_vars.end();++varIter) {
      // we look for loop rotation vars and normal vars that
      // need to be profiled
      if(!profile_vars.inSet(*varIter))
        continue ;
      // now we decorate the graph to include the profiling
      // we first obtain the ident of rules that produce the variable
      digraph::vertexSet source_rules = grt[varIter->ident()] ;
      // if source_rules is empty, then the variable is not
      // present in this graph, we need a search to find out
      // the source_rules
      if(source_rules == EMPTY) {
        for(ruleSet::const_iterator ruleIter=rules.begin();
            ruleIter!=rules.end();++ruleIter) {
          variableSet targets = ruleIter->targets() ;
          if(targets.inSet(*varIter))
            source_rules += ruleIter->ident() ;
        }
      }

      digraph::vertexSet target_rules ;
      if(allvertices.inSet(varIter->ident())) {
        target_rules = gr[varIter->ident()] ;
      } else {
        for(digraph::vertexSet::const_iterator veri=source_rules.begin();
            veri!=source_rules.end();++veri)
          target_rules += get_vertexSet(extract_rules(gr[*veri])) ;
      }
      // but we need to exclude any MEMPROFILEA compilers from it
      digraph::vertexSet take_off_from_targets ;
      ruleSet trules = extract_rules(target_rules) ;
      for(ruleSet::const_iterator ri=trules.begin();
          ri!=trules.end();++ri)
        if(ri->get_info().qualifier() == "MEMPROFILEA")
          take_off_from_targets += ri->ident() ;

      target_rules -= take_off_from_targets ;

      variable sv("PROFILE") ;
      rule profile_alloc_rule = create_rule(sv,*varIter,"MEMPROFILEA") ;
      
      gr.add_edges(source_rules,profile_alloc_rule.ident()) ;
      gr.add_edge(varIter->ident(),profile_alloc_rule.ident()) ;
      if(target_rules != EMPTY)
        gr.add_edges(profile_alloc_rule.ident(),target_rules) ;        

      rcm[profile_alloc_rule] = new memProfileAlloc_compiler(*varIter) ;
      
    }
  }

  void memProfileAllocDecoVisitor::visit(loop_compiler& lc) {
    edit_gr(lc.collapse_gr,lc.rule_compiler_map,-lc.cid) ;
    edit_gr(lc.advance_gr,lc.rule_compiler_map,lc.cid) ;
  }

  void memProfileAllocDecoVisitor::visit(dag_compiler& dc) {
    edit_gr(dc.dag_gr,dc.rule_compiler_map,dc.cid) ;
  }

  void memProfileAllocDecoVisitor::visit(conditional_compiler& cc) {
    edit_gr(cc.cond_gr,cc.rule_compiler_map,cc.cid) ;
  }



} // end of namespace Loci
