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
#include "comp_tools.h"


#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include<map>
using std::map ;
#include <list>
using std::list ;
#include <algorithm>
using std::sort ;
using std::find_if ;
using std::copy ;
using std::find ;
#include <functional>
//using std::bind2nd ;
//using std::ptr_fun ;
#include <iterator>
using std::ostream_iterator ;
#include <utility>
using std::pair ;
using std::make_pair ;
using std::min ;

using std::cerr ;
using std::endl ;
using std::ostream ;

namespace Loci {

  //////////////////////////////////////////////////////////////
  // function to get all the recur targets of a variableSet
  // from a given recur mapping table
  //////////////////////////////////////////////////////////////
  variableSet
  get_recur_target_for_vars(const variableSet& vars,
                            const map<variable,variableSet>& t) {
    variableSet ret ;
    for(variableSet::const_iterator vi=vars.begin();
        vi!=vars.end();++vi) {
      ret += get_all_recur_vars(t,*vi) ;
    }

    return ret ;
  }


  //////////////////////////////////////////////////////////////
  // recurInfoVisitor
  //////////////////////////////////////////////////////////////
  void recurInfoVisitor::gather_info2(const ruleSet& rs) {
    for(ruleSet::const_iterator ruleIter=rs.begin();
        ruleIter!=rs.end();++ruleIter) {
      // we check for rename (inplace update rules)
      set<vmap_info>::const_iterator vmsi ;
      for(vmsi=ruleIter->get_info().desc.targets.begin();
          vmsi!=ruleIter->get_info().desc.targets.end(); ++vmsi) {
        if(vmsi->assign.size() != 0) {
          for(size_t i=0;i<vmsi->assign.size();++i) {
            variable new_name = vmsi->assign[i].first ;
            variable old_name = vmsi->assign[i].second ;
            
            recur_vars_t2s[new_name] += old_name ;
            recur_vars_s2t[old_name] += new_name ;
            recur_source_vars += old_name ;
            recur_target_vars += new_name ;
            
            rename_t2s[new_name] += old_name ;
            rename_s2t[old_name] += new_name ;
            rename_source_vars += old_name ;
            rename_target_vars += new_name ;
          }
        }
      }
    }
  }
  
  void recurInfoVisitor::gather_info(const digraph& gr) {
    // obtain all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    // looping over all the rules and gather info
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      // recurrence variables are the target and source of
      // these internal rules
      if(ruleIter->get_info().qualifier() == "generalize") {
        variable target,source ;
        target = *(ruleIter->targets().begin()) ;
        source = *(ruleIter->sources().begin()) ;
        
        recur_vars_t2s[target] += source ;
        recur_vars_s2t[source] += target ;
        recur_source_vars += source ;
        recur_target_vars += target ;
        
        generalize_t2s[target] += source ;
        generalize_s2t[source] += target ;
        generalize_source_vars += source ;
        generalize_target_vars += target ;
      }
      else if(ruleIter->get_info().qualifier() == "promote") {
        variable target,source ;
        target = *(ruleIter->targets().begin()) ;
        source = *(ruleIter->sources().begin()) ;
        
        recur_vars_t2s[target] += source ;
        recur_vars_s2t[source] += target ;
        recur_source_vars += source ;
        recur_target_vars += target ;
        
        promote_t2s[target] += source ;
        promote_s2t[source] += target ;
        promote_source_vars += source ;
        promote_target_vars += target ;
      }
      else if(ruleIter->get_info().qualifier() == "priority") {
        variable target,source ;
        target = *(ruleIter->targets().begin()) ;
        source = *(ruleIter->sources().begin()) ;
        
        recur_vars_t2s[target] += source ;
        recur_vars_s2t[source] += target ;
        recur_source_vars += source ;
        recur_target_vars += target ;
        
        priority_t2s[target] += source ;
        priority_s2t[source] += target ;
        priority_source_vars += source ;
        priority_target_vars += target ;
      }
      else {
        // we check for rename (inplace update rules)
        set<vmap_info>::const_iterator vmsi ;
        for(vmsi=ruleIter->get_info().desc.targets.begin();
            vmsi!=ruleIter->get_info().desc.targets.end(); ++vmsi) {
          if(vmsi->assign.size() != 0) {
            for(size_t i=0;i<vmsi->assign.size();++i) {
              variable new_name = vmsi->assign[i].first ;
              variable old_name = vmsi->assign[i].second ;
              
              recur_vars_t2s[new_name] += old_name ;
              recur_vars_s2t[old_name] += new_name ;
              recur_source_vars += old_name ;
              recur_target_vars += new_name ;
              
              rename_t2s[new_name] += old_name ;
              rename_s2t[old_name] += new_name ;
              rename_source_vars += old_name ;
              rename_target_vars += new_name ;
            }
          }
        }
      }
    }
    
  }
  
  void recurInfoVisitor::visit(loop_compiler& lc) {
    gather_info(lc.collapse_gr) ;
    gather_info(lc.advance_gr) ;
  }

  void recurInfoVisitor::visit(dag_compiler& dc) {
    gather_info(dc.dag_gr) ;
  }

  void recurInfoVisitor::visit(conditional_compiler& cc) {
    gather_info(cc.cond_gr) ;
  }

  void recurInfoVisitor::visit(impl_recurse_compiler& irc) {
    gather_info2(irc.get_rules()) ;
  }

  void recurInfoVisitor::visit(recurse_compiler& rc) {
    gather_info2(rc.get_rules()) ;
  }

  variableSet
  recurInfoVisitor::get_reachable(const variableSet& vs) const {
    variableSet working, visited ;
    map<variable, variableSet>::const_iterator mi ;
    working += vs ;
    while(working != EMPTY) {
      visited += working ;
      variableSet next ;
      for(variableSet::const_iterator vi=working.begin();
          vi!=working.end();++vi) {
        mi = recur_vars_t2s.find(*vi) ;
        if(mi != recur_vars_t2s.end())
          next += mi->second ;
        mi = recur_vars_s2t.find(*vi) ;
        if(mi != recur_vars_s2t.end())
          next += mi->second ;
      }
      next -= visited ;
      working = next ;
    }
    visited -= vs ;
    return visited ;
  }
  
  /////////////////////////////////////////////////////////////
  // snInfoVisitor
  /////////////////////////////////////////////////////////////
  void snInfoVisitor::fill_subnode_table(const digraph& gr, int id) {
    // first get all the rules
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    ruleSet rules = extract_rules(allvertices) ;

    set<int> subnode ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(is_super_node(ruleIter)) {
        subnode.insert(get_supernode_num(*ruleIter)) ;
      }
    }
    subnode_table[id] = subnode ;
  }

  void snInfoVisitor::visit(loop_compiler& lc) {
    graph_sn.insert(lc.cid) ;
    loop_sn.insert(lc.cid) ;

    fill_subnode_table(lc.loop_gr,lc.cid) ;

    rule collapse_node ;
    ruleSet all_rules = extract_rules(lc.collapse_gr.get_all_vertices()) ;
    for(ruleSet::const_iterator ri=all_rules.begin();
        ri!=all_rules.end();++ri) {
      if(ri->target_time().before(ri->source_time())) {
        collapse_node = *ri ;
        break ;
      }
    }
    if(!is_super_node(collapse_node)) {
      if(MPI_rank == 0)
        cerr << "Internal error, the collapse part of loop compiler: "
             << lc.cid << " does not have a collapse rule" << endl ;
      debugout << "Internal error, the collapse part of loop compiler: "
               << lc.cid << " does not have a collapse rule" << endl ;
      
      Loci::Abort() ;
    }
    int collapse_id = get_supernode_num(collapse_node) ;
    if(collapse_id == -1) {
      if(MPI_rank == 0)
        cerr << "Error: conditional node has wrong id number" << endl ;
      debugout << "Error: conditional node has wrong id number" << endl ;
      Loci::Abort() ;
    }

    loop_col_table[lc.cid] = collapse_id ;

  }
  
  void snInfoVisitor::visit(dag_compiler& dc) {
    graph_sn.insert(dc.cid) ;
    fill_subnode_table(dc.dag_gr,dc.cid) ;
  }
  
  void snInfoVisitor::visit(conditional_compiler& cc) {
    graph_sn.insert(cc.cid) ;
    cond_sn.insert(cc.cid) ;
    fill_subnode_table(cc.cond_gr,cc.cid) ;
  }

  ///////////////////////////////////////////////////////////
  // allocDeleteStat
  ///////////////////////////////////////////////////////////
  allocDeleteStat::allocDeleteStat(const map<int,variableSet>& alloc_table,
                                   const map<int,variableSet>& delete_table,
                                   const map<variable,variableSet>& t2s,
                                   const map<variable,variableSet>& s2t,
                                   const variableSet& rsv,
                                   const variableSet& rtv)
    :recur_vars_t2s(t2s),recur_vars_s2t(s2t),
     recur_source_vars(rsv),recur_target_vars(rtv){
    map<int,variableSet>::const_iterator mi ;
    for(mi=alloc_table.begin();mi!=alloc_table.end();++mi)
      allocated_vars += mi->second ;
    for(mi=delete_table.begin();mi!=delete_table.end();++mi)
      deleted_vars += mi->second ;
  }

  ostream& allocDeleteStat::report(ostream& s) const {
    variableSet d_a, a_d ;
    variableSet::const_iterator vi,vii ;

    s << "***** Begin allocation and deletion info report *****" << endl ;
    variableSet miss_alloc ;
    for(vi=allocated_vars.begin();vi!=allocated_vars.end();++vi) {
      if(recur_target_vars.inSet(*vi))
        miss_alloc += *vi ;
    }
    if(miss_alloc != EMPTY) {
      s << endl ;
      s << "These variables should not been allocated, because they are "
        << "recurrence target variables: " << endl << miss_alloc << endl ;
    }


    d_a = deleted_vars - allocated_vars ;

    // some recurrence variable processing
    // NEEDS SOME COMMENTS LATER HERE!!!!!
    variableSet remove ;
    for(vi=d_a.begin();vi!=d_a.end();++vi) {
      variableSet init_sources = get_leaf_recur_vars(recur_vars_t2s,*vi);
      
      if(init_sources != EMPTY) {
        bool all_allocated = true ;
        for(vii=init_sources.begin();vii!=init_sources.end();++vii)
          if(!allocated_vars.inSet(*vii)) {
            all_allocated = false ;
            break ;
          }

        if(all_allocated)
          remove += *vi ;
      }
    }
    d_a -= remove ;

    remove = variableSet(EMPTY) ;
    a_d = allocated_vars - deleted_vars ;
    for(vi=a_d.begin();vi!=a_d.end();++vi) {
      variableSet all_targets = get_all_recur_vars(recur_vars_s2t,*vi) ;
      
      if(all_targets != EMPTY) {
        for(vii=all_targets.begin();vii!=all_targets.end();++vii)
          if(deleted_vars.inSet(*vii)) {
            remove += *vi ;
            break ;
          }
      }
    }
    a_d -= remove ;

    if(a_d == EMPTY) {
      s << endl ;
      s << "All allocated variables are scheduled to be deleted" << endl ;
    }else {
      s << endl ;
      s << "These variables are scheduled to be allocated, but not to be deleted:" << endl
        << a_d << endl ;
    }

    if(d_a != EMPTY) {
      s << endl ;
      s << "These variables are scheduled to be deleted, but they are input variables: " << endl << d_a << endl ;
    }
    s << endl ;
    s << "***** Finish allocation and deletion info report *****" << endl ;

    return s ;
  }

  /////////////////////////////////////////////////////////////////////
  // function to get the parent table
  /////////////////////////////////////////////////////////////////////  
  map<int,int>
  get_parentnode_table(const map<int,set<int> >& subnode_table) {
    map<int,int> ret ;

    // first get the parent node table
    map<int,set<int> >::const_iterator mi ;
    for(mi=subnode_table.begin();mi!=subnode_table.end();++mi) {
      // get all the subnodes of a node
      set<int> subnode = mi->second ;
      for(set<int>::const_iterator si=subnode.begin();
          si!=subnode.end();++si) {
        map<int,int>::const_iterator found ;
        found = ret.find(*si) ;
        if(found != ret.end()) {
          if(MPI_rank==0)
            cerr << "multilevel graph error!" << endl ;
          debugout << "multilevel graph error!" << endl ;
          Loci::Abort() ;
        }
        ret[*si] = mi->first ;
      }
    }
    // super node 0 has no parent node
    ret[0] = -1 ;

    return ret ;
  }

  ///////////////////////////////////////////////////////////////////////
  // function to get the allocation of each loop
  ///////////////////////////////////////////////////////////////////////
  map<int,variableSet>
  get_loop_alloc_table(const map<int,variableSet>& alloc_table,
                       const map<int,set<int> >& subnode_table,
                       const set<int>& loop_sn,
                       const set<int>& graph_sn,
                       const map<variable,variableSet>& rvs2t) {
    map<int,variableSet> ret ;

    for(set<int>::const_iterator si=loop_sn.begin();
        si!=loop_sn.end();++si) {
      // get all the subnode of a loop
      set<int> all_loop_subnode ;
      set<int> working ;

      all_loop_subnode.insert(*si) ;
      working.insert(*si) ;
      
      while(!working.empty()) {
        set<int> tmp ;
        for(set<int>::const_iterator si2=working.begin();
            si2!=working.end();++si2) {
          // skip recursive super node
          if(!inSet(graph_sn,*si2))
            continue ;
          map<int,set<int> >::const_iterator found ;
          found = subnode_table.find(*si2) ;
          FATAL(found == subnode_table.end()) ;
          tmp.insert( (found->second).begin(), (found->second).end()) ;
        }
        working = tmp ;
        all_loop_subnode.insert(tmp.begin(),tmp.end()) ;
      }

      // compute all the allocations of this loop
      variableSet loop_alloc ;
      for(set<int>::const_iterator si3=all_loop_subnode.begin();
          si3!=all_loop_subnode.end();++si3) {
        // skip recursive super node
        if(!inSet(graph_sn,*si3))
          continue ;
        map<int,variableSet>::const_iterator found ;
        if(inSet(loop_sn,*si3)) {
          // collapse part
          found = alloc_table.find(-(*si3)) ;
          FATAL(found == alloc_table.end()) ;
          loop_alloc += found->second ;
          // advance part
          found = alloc_table.find(*si3) ;
          FATAL(found == alloc_table.end()) ;
          loop_alloc += found->second ;
        } else {
          found = alloc_table.find(*si3) ;
          FATAL(found == alloc_table.end()) ;
          loop_alloc += found->second ;
        }
      }
      // we need to include any recurrence variables
      variableSet loop_alloc_recur ;
      for(variableSet::const_iterator vi=loop_alloc.begin();
          vi!=loop_alloc.end();++vi) {
        // get all the recurrence target variables
        variableSet targets = get_all_recur_vars(rvs2t,*vi) ;
        if(targets != EMPTY)
          loop_alloc_recur += targets ;
      }
      loop_alloc += loop_alloc_recur ;
      // fill the ret table
      ret[*si] = loop_alloc ;
    }

    return ret ;
  }

  /////////////////////////////////////////////////////////////////////
  // rotateListVisitor
  /////////////////////////////////////////////////////////////////////
  namespace {
    inline bool offset_sort(const variable& v1, const variable& v2)
    { return v1.get_info().offset > v2.get_info().offset ; }
  }
  void rotateListVisitor::visit(loop_compiler& lc) {
    // compute the rotate list first
    
    map<variable,list<variable> > vlist ;
    for(variableSet::const_iterator vi=lc.all_loop_vars.begin();
        vi!=lc.all_loop_vars.end();++vi) {
      if(vi->time() == lc.tlevel && !vi->assign) {
        variable var_stationary(*vi,time_ident()) ;
        
        list<variable> s ;
        s.push_back(*vi) ;
        vlist[var_stationary].merge(s,offset_sort) ;
      }
    }
    
    map<variable,list<variable> >::const_iterator ii ;
    for(ii=vlist.begin();ii!=vlist.end();++ii) {
      if(ii->second.size() < 2) 
        continue ;
      std::list<variable>::const_iterator jj ;
      bool overlap = false ;
      variableSet vtouch ;
      for(jj=ii->second.begin();jj!=ii->second.end();++jj) {
        variableSet aliases = scheds.get_aliases(*jj) ;
        variableSet as = aliases ;
        variableSet::const_iterator vi ;
        
        for(vi=aliases.begin();vi!=aliases.end();++vi) 
          as += scheds.get_synonyms(*vi) ;
        if((as & vtouch) != EMPTY)
          overlap = true ;
        vtouch += as ;
      }
      if(!overlap) {
        lc.rotate_lists.push_back(ii->second) ;
        // set the scheds rotation information
        variableSet drots ;
        std::list<variable>::const_iterator rotii ;
        for(rotii=ii->second.begin();rotii!=ii->second.end();++rotii)
          drots += *rotii ;
        scheds.set_variable_rotations(drots) ;
      } else {
        if(ii->second.size() !=2) {
          if(MPI_rank == 0) 
            cerr << "unable to have history on variables aliased in time"
                 << endl
                 << "error occured on variable " << ii->first
                 << "{" << lc.tlevel << "}"
                 << endl ;
          
          debugout << "unable to have history on variables aliased in time"
                   << endl
                   << "error occured on variable " << ii->first
                   << "{" << lc.tlevel << "}"
                   << endl ;
          Loci::Abort() ;
        }
        variableSet orv ;
        for(list<variable>::const_iterator iii=ii->second.begin();
            iii!=ii->second.end();++iii)
          orv += *iii ;
        for(variableSet::const_iterator iii=orv.begin();
            iii!=orv.end();++iii)
          overlap_rotvars[*iii] = orv ;
      }
    }

    // then get all the variables in the list and fill the table
    variableSet rotate_vars ;
    /*
    for(list<list<variable> >::const_iterator li=lc.rotate_lists.begin();
        li!=lc.rotate_lists.end();++li) {
      for(list<variable>::const_iterator lii=li->begin();
          lii!=li->end();++lii) {
        rotate_vars += *lii ;
      }
    }
    */
    for(ii=vlist.begin();ii!=vlist.end();++ii) {
      if(ii->second.size() < 2) 
        continue ;
      std::list<variable>::const_iterator jj ;
      for(jj=ii->second.begin();jj!=ii->second.end();++jj) {
        rotate_vars += *jj ;
      }
    }
    // then we include all the recurrence variable in
    variableSet add_rotate_vars ;
    for(variableSet::const_iterator vi=rotate_vars.begin();
        vi!=rotate_vars.end();++vi) {
      add_rotate_vars += get_all_recur_vars(rvs2t,*vi) ;
      add_rotate_vars += get_all_recur_vars(rvt2s,*vi) ;
    }
    rotate_vars += add_rotate_vars ;
    
    rotate_vars_table[lc.cid] = rotate_vars ;

    // and then compute the shared variables between
    // the advance part and the collapse part
    
    // we then get the variables shared between the
    // advance and collapse part of the loop
    digraph::vertexSet collapse_gr_ver = lc.collapse_gr.get_all_vertices() ;
    digraph::vertexSet advance_gr_ver = lc.advance_gr.get_all_vertices() ;
    variableSet shared = extract_vars(collapse_gr_ver & advance_gr_ver) ;
    // exclude time level variable
    variable tvar = variable(lc.tlevel) ;
    shared -= tvar ;

    loop_shared_table[lc.cid] = shared ;
      
  }

  ////////////////////////////////////////////////////////////////
  // dagCheckVisitor
  ////////////////////////////////////////////////////////////////
  // return true if gr has NO cycles
  // return false if gr has cycles
  bool dagCheckVisitor::check_dag(digraph gr) {
    // if graph is empty, it is dag
    if(is_dg_empty(gr))
      return true ;

    /*
    // get all the leaf nodes of this graph
    digraph::vertexSet leaf_nodes ;
    do {
      digraph::vertexSet allv ;

      leaf_nodes = EMPTY ;
      allv = gr.get_all_vertices() ;
      for(digraph::vertexSet::const_iterator vi=allv.begin();
          vi!=allv.end();++vi) {
        if(gr[*vi] == EMPTY)
          leaf_nodes += *vi ;
      }
      gr.remove_vertices(leaf_nodes) ;
      
    }while(leaf_nodes != EMPTY) ;

    cycle = gr ;
    
    return is_dg_empty(gr) ;
    */

    // a different method --- use component sort
    vector<digraph::vertexSet> clusters =
      component_sort(gr).get_components() ;
    
    for(vector<digraph::vertexSet>::size_type i=0;i<clusters.size();++i) {
      digraph::vertexSet potential_cycle_v = clusters[i] ;
      if(potential_cycle_v.size() != 1) {
        cycle = gr.subgraph(potential_cycle_v) ;
        return false ;
      }
    }
    return true ;
  }

  void print_cycles(digraph gr) {
    // if graph is empty, it is dag
    if(is_dg_empty(gr))
      return ;

    vector<digraph::vertexSet> clusters =
      component_sort(gr).get_components() ;
    
    for(vector<digraph::vertexSet>::size_type i=0;i<clusters.size();++i) {
      digraph::vertexSet potential_cycle_v = clusters[i] ;
      if(potential_cycle_v.size() != 1) {
        if(MPI_rank == 0) {
          cerr << "potential cycle contains variables: " << extract_vars(potential_cycle_v) << endl ;
          cerr << "rules:" << endl << extract_rules(potential_cycle_v) << endl ;
        }
        debugout << "potential cycle contains variables: " << extract_vars(potential_cycle_v) << endl ;
        debugout << "rules:" << endl << extract_rules(potential_cycle_v) << endl ;
        
      }

    }
  }

  void dagCheckVisitor::visit(loop_compiler& lc) {
    if(!check_dag(lc.collapse_gr)) {
      if(MPI_rank == 0)
        cerr << "ERROR: the collapse graph of loop super node("
             << lc.cid << ") has cycle(s)" << endl ;
      debugout << "ERROR: the collapse graph of loop super node("
               << lc.cid << ") has cycle(s)" << endl ;
      
      print_cycles(lc.collapse_gr) ;
      if(viz)
        visualize(cerr) ;
      Loci::Abort() ;
    }
    if(!check_dag(lc.advance_gr)) {
      if(MPI_rank == 0)
        cerr << "ERROR: the advance graph of loop super node("
             << lc.cid << ") has cycle(s)" << endl ;
      debugout << "ERROR: the advance graph of loop super node("
               << lc.cid << ") has cycle(s)" << endl ;
      
      print_cycles(lc.advance_gr) ;
      if(viz)
        visualize(cerr) ;
      Loci::Abort() ;
    }
  }
  
  void dagCheckVisitor::visit(dag_compiler& dc) {
    if(!check_dag(dc.dag_gr)) {
      if(MPI_rank == 0)
        cerr << "ERROR: the graph of dag super node("
             << dc.cid << ") has cycle(s)" << endl ;

      debugout << "ERROR: the graph of dag super node("
               << dc.cid << ") has cycle(s)" << endl ;
      
      print_cycles(dc.dag_gr) ;
      if(viz)
        visualize(cerr) ;
      Loci::Abort() ;
    }
  }
  
  void dagCheckVisitor::visit(conditional_compiler& cc) {
    if(!check_dag(cc.cond_gr)) {
      if(MPI_rank == 0)
        cerr << "ERROR: the graph of conditional super node("
             << cc.cid << ") has cycle(s)" << endl ;
      debugout << "ERROR: the graph of conditional super node("
               << cc.cid << ") has cycle(s)" << endl ;

      print_cycles(cc.cond_gr) ;
      if(viz)
        visualize(cerr) ;
      Loci::Abort() ;
    }
  }

  ////////////////////////////////////////////////////////////////
  // unTypedVarVisitor
  ////////////////////////////////////////////////////////////////
  unTypedVarVisitor::
  unTypedVarVisitor(const std::map<variable,variableSet>& s2t,
                    const std::map<variable,variableSet>& t2s,
                    const variableSet& input)
    :recur_vars_s2t(s2t),recur_vars_t2s(t2s) {

    typed_vars += input ;
    variableSet add ;
    for(variableSet::const_iterator vi=typed_vars.begin();
        vi!=typed_vars.end();++vi) {
      variableSet all_targets = get_all_recur_vars(recur_vars_s2t,*vi) ;
      add += all_targets ;
    }
    typed_vars += add ;
  }

  
  void unTypedVarVisitor::discover(const digraph& gr) {
    digraph::vertexSet all_vertices = gr.get_all_vertices() ;
    ruleSet all_rules = extract_rules(all_vertices) ;
    variableSet all_vars ;
    for(ruleSet::const_iterator ri=all_rules.begin();
        ri!=all_rules.end();++ri)
      all_vars += ri->targets() ;

    // first assume all variables in the graph are untyped
    untyped_vars += all_vars ;
    // then we discover those typed, and remove them from untyped_vars

    // first we get all the non virtual rules
    // they are all the concrete rules (and possibly chomps)
    ruleSet concrete_rules ;
    for(ruleSet::const_iterator ri=all_rules.begin();
        ri!=all_rules.end();++ri)
      if(!is_virtual_rule(*ri))
        concrete_rules += *ri ;
    // all their targets are typed
    variableSet typed_here ;
    for(ruleSet::const_iterator ri=concrete_rules.begin();
        ri!=concrete_rules.end();++ri)
      typed_here += ri->targets() ;

    typed_vars += typed_here ;
    // and all the recurrence variables of typed_here are typed
    for(variableSet::const_iterator vi=typed_here.begin();
        vi!=typed_here.end();++vi)
      typed_vars += get_all_recur_vars(recur_vars_s2t,*vi) ;

    for(variableSet::const_iterator vi=typed_here.begin();
        vi!=typed_here.end();++vi)
      typed_vars += get_all_recur_vars(recur_vars_t2s,*vi) ;

    // remove...
    untyped_vars -= typed_vars ;
  }
  
  void unTypedVarVisitor::visit(loop_compiler& lc) {
    discover(lc.collapse_gr) ;
    discover(lc.advance_gr) ;
  }
  
  void unTypedVarVisitor::visit(dag_compiler& dc) {
    discover(dc.dag_gr) ;
  }
  
  void unTypedVarVisitor::visit(conditional_compiler& cc) {
    discover(cc.cond_gr) ;
  }

  ////////////////////////////////////////////////////////////////
  // getAllVarVisitor
  ////////////////////////////////////////////////////////////////

  void getAllVarVisitor::collect_vars(const digraph& gr) {
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    all_vars += extract_vars(allvertices) ;
  }

  void getAllVarVisitor::visit(loop_compiler& lc) {
    collect_vars(lc.loop_gr) ;
  }
  
  void getAllVarVisitor::visit(dag_compiler& dc) {
    collect_vars(dc.dag_gr) ;
  }
  
  void getAllVarVisitor::visit(conditional_compiler& cc) {
    collect_vars(cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // promotePPVisitor
  /////////////////////////////////////////////////////////////////
  promotePPVisitor::
  promotePPVisitor(const map<variable,variableSet>& pt2s,
                   const map<variable,variableSet>& ps2t,
                   const variableSet& psource,
                   const variableSet& ptarget,
                   const set<int>& gsn,
                   const std::set<int>& csn,
                   variableSet& input):
    promote_source_vars(psource), promote_target_vars(ptarget),
    promote_t2s(pt2s), promote_s2t(ps2t), graph_sn(gsn),
    cond_sn(csn) {
    
    reserved_vars += input ;
    reserved_vars += get_recur_target_for_vars(input,promote_s2t) ;
  }

  // given a rule set a variable can reach, if the variable
  // must be deleted in this graph, then it is a rep
  bool promotePPVisitor::is_rep(const digraph& gr, ruleSet rules) {
    if(rules == EMPTY) {
      return true ;
    }
    
    int supernode_num = 0 ;
    int supernode_id = 1 ; // vertex number of the supernode
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
    if(supernode_num == 1) {
      // check to see if it is a conditional node
      // if it is, we delete it in the graph
      if(inSet(cond_sn,get_supernode_num(supernode)))
        return true ;
      else if(rules.size() == 1)
        return false ;
    }
    if( (rules.size() == 1) && is_virtual_rule(*rules.begin()))
      return false ;
    if(supernode_num != 1)
      return true ;
    /*
    if( (rules.size() == 1) && is_virtual_rule(*rules.begin()))
      return false ;
    if(supernode_num != 1)
      return true ;
    */
    // the remaining case must be there is only one super node
    // and at least one non super node, we check if there is path
    // from each non super node to the super node in the graph
    FATAL(supernode_id > 0) ;
    rules -= rule(supernode_id) ;
    bool ret = false ;
    for(ruleSet::const_iterator ruleIter=rules.begin();
        ruleIter!=rules.end();++ruleIter) {
      if(!has_path(gr,ruleIter->ident(),supernode_id)) {
        ret = true ;
        break ;
      }
    }
    
    return ret ;
  }

  void promotePPVisitor::pick_rep(const digraph& gr) {
    // get variables to evaluate
    digraph::vertexSet allv = gr.get_all_vertices() ;
    
    ruleSet rules = extract_rules(allv) ;
    variableSet vars ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      vars += ri->sources() ;

    vars &= variableSet(promote_source_vars + promote_target_vars) ;
    vars -= reserved_vars ;

    for(variableSet::const_iterator vi=vars.begin();
        vi!=vars.end();++vi) {
      if(processed.inSet(*vi))
        continue ;
      // get the final target
      variableSet targets = get_leaf_recur_vars(promote_s2t,*vi) ;

      if(targets.size() == 1) {
        if(allv.inSet(targets.begin()->ident())) {
          processed += *vi ;
          continue ;
        }
      }

      // get the rules all_source can reach
      // we need to get all the rules that all the
      // promoted vars can reach
      digraph::vertexSet promoted =
        get_vertexSet(get_all_recur_vars(promote_s2t,*vi)) ;
      promoted += vi->ident() ;
      promoted &= allv ;
      
      ruleSet others ;
      for(digraph::vertexSet::const_iterator pvi=promoted.begin();
          pvi!=promoted.end();++pvi)
        others += extract_rules(gr[*pvi]) ;

      if(others == EMPTY) {
        cerr << "others should not be empty, promoted = " << extract_vars(promoted) << endl ;
      }
      FATAL(others == EMPTY) ;

      if(is_rep(gr,others)) {
        rep += *vi ;
        variableSet all_source = get_all_recur_vars(promote_t2s,*vi) ;
        variableSet all_target = get_all_recur_vars(promote_s2t,*vi) ;
        processed += *vi ;
        processed += all_source ;
        processed += all_target ;
        remaining += all_source ;
        remaining += all_target ;
      }
    }
  }
  
  void promotePPVisitor::visit(loop_compiler& lc) {
    pick_rep(lc.loop_gr) ;
  }

  void promotePPVisitor::visit(dag_compiler& dc) {
    pick_rep(dc.dag_gr) ;
  }

  void promotePPVisitor::visit(conditional_compiler& cc) {
    pick_rep(cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // function to analyze the renamed target variable cluster
  // and pick out a representitive variable for each cluster,
  // only this variable is allowed to be deleted
  /////////////////////////////////////////////////////////////////
  variableSet pick_rename_target(const map<variable,variableSet>& s2t,
                                 const map<variable,variableSet>& t2s,
                                 const variableSet& allvars) {
    variableSet filtered_vars ;
    
    variableSet processed ;
    for(variableSet::const_iterator vi=allvars.begin();
        vi!=allvars.end();++vi) {
      if(processed.inSet(*vi))
        continue ;

      variableSet all_other_target = get_all_leaf_target(t2s,s2t,*vi) ;

      if(all_other_target == EMPTY)
        continue ;
      
      processed += all_other_target ;
      
      vector<variable> all_targets ;
      all_targets.push_back(*vi) ;
      for(variableSet::const_iterator vi2=all_other_target.begin();
          vi2!=all_other_target.end();++vi2)
        all_targets.push_back(*vi2) ;

      sort(all_targets.begin(),all_targets.end(),time_before) ;

      
      vector<variable>::iterator pos ;
      vector<variable>::iterator old_pos = all_targets.begin() ;
      while(old_pos != all_targets.end()) {
	// old way
	//	pos = find_if(old_pos+1,all_targets.end(),
	//		      bind2nd(ptr_fun(time_after),*old_pos)
	//		      ) ;
        // now replaced with lambda expression
	pos = find_if(old_pos+1,all_targets.end(),
		      [old_pos](variable elem) { return time_after(elem,*old_pos) ; }) ;
        if(pos - old_pos > 1) {
          variableSet time_ident_vars ;
          variableSet real_varset ;
          for(vector<variable>::const_iterator vi2=old_pos;
              vi2!=pos;++vi2) {
            time_ident_vars += *vi2 ;
            variable v = *vi2 ;
            while(v.get_info().priority.size() != 0)
              v = v.drop_priority() ;
            real_varset += v ;
          }
          if(real_varset.size() > 1) {
            if(MPI_rank == 0)
              cerr << "WARNING: These renamed variables coexist in the same time level, and they refer to the same memory location, this is dangerous!: " << time_ident_vars << endl ;
            debugout << "WARNING: These renamed variables coexist in the same time level, and they refer to the same memory location, this is dangerous!: " << time_ident_vars << endl ;
          }
        }
        old_pos = pos ;
      }

      for(vector<variable>::const_iterator vi2=all_targets.begin()+1;
          vi2!=all_targets.end();++vi2)
        filtered_vars += *vi2 ;
    }

    return filtered_vars ;
  }

  ///////////////////////////////////////////////////////////////////
  // function that checks some preconditions that all the recurrence
  // variables should meet.
  // generalize and priority rules should not form branches
  // promoted variables should not been renamed, and input variables
  // also should not been renamed
  ///////////////////////////////////////////////////////////////////
  void check_recur_precondition(const recurInfoVisitor& v,
                                const variableSet& input) {
    variableSet::const_iterator vi ;
    map<variable,variableSet>::const_iterator found ;

    variableSet generalize_source = v.get_generalize_source_vars() ;
    variableSet promote_source = v.get_promote_source_vars() ;
    variableSet priority_source = v.get_priority_source_vars() ;
    variableSet rename_source = v.get_rename_source_vars() ;
    
    map<variable,variableSet> generalize_s2t = v.get_generalize_s2t() ;
    map<variable,variableSet> promote_s2t = v.get_promote_s2t() ;
    map<variable,variableSet> priority_s2t = v.get_priority_s2t() ;
    map<variable,variableSet> rename_s2t = v.get_rename_s2t() ;
    
    // generalize rules should not form branches
    for(vi=generalize_source.begin();vi!=generalize_source.end();++vi) {
      found = generalize_s2t.find(*vi) ;
      FATAL(found == generalize_s2t.end()) ;
      if(found->second.size() > 1) {
        if(MPI_rank == 0)
          cerr << "WARNING: " << *vi << " is in the chain of "
               << "generalize rules, but it forms multiple targets: "
               << found->second << endl ;
        debugout << "WARNING: " << *vi << " is in the chain of "
               << "generalize rules, but it forms multiple targets: "
               << found->second << endl ;
      }
      if(promote_source.inSet(*vi)) {
        if(MPI_rank == 0)
          cerr << "\tit is also in the chain of promote rules " ;
        debugout << "\tit is also in the chain of promote rules " ;
        
        found = promote_s2t.find(*vi) ;
        FATAL(found == promote_s2t.end()) ;
        if(MPI_rank == 0)
          cerr << "forms targets: " << found->second << endl ;
        debugout << "forms targets: " << found->second << endl ;
      }
      if(priority_source.inSet(*vi)) {
        if(MPI_rank == 0)
          cerr << "\tit is also in the chain of priority rules " ;
        debugout << "\tit is also in the chain of priority rules " ;
        
        found = priority_s2t.find(*vi) ;
        FATAL(found == priority_s2t.end()) ;
        if(MPI_rank == 0)
          cerr << "forms targets: " << found->second << endl ;
        debugout << "forms targets: " << found->second << endl ;
      }
      if(rename_source.inSet(*vi)) {
        if(MPI_rank == 0)
          cerr << "\tit is also in the chain of rename rules " ;
        debugout << "\tit is also in the chain of rename rules " ;
        found = rename_s2t.find(*vi) ;
        FATAL(found == rename_s2t.end()) ;
        if(MPI_rank == 0)
          cerr << "forms targets: " << found->second << endl ;
        debugout << "forms targets: " << found->second << endl ;
      }
    }

    // priority rules should not form branches
    for(vi=priority_source.begin();vi!=priority_source.end();++vi) {
      found = priority_s2t.find(*vi) ;
      FATAL(found == priority_s2t.end()) ;
      if(found->second.size() > 1) {
        if(MPI_rank == 0)
          cerr << "WARNING: " << *vi << " is in the chain of "
               << "priority rules, but it forms multiple targets: "
               << found->second << endl ;
        debugout << "WARNING: " << *vi << " is in the chain of "
             << "priority rules, but it forms multiple targets: "
             << found->second << endl ;
      }
      if(generalize_source.inSet(*vi)) {
        if(MPI_rank == 0)
          cerr << "\tit is also in the chain of generalize rules " ;
        debugout << "\tit is also in the chain of generalize rules " ;
        found = generalize_s2t.find(*vi) ;
        FATAL(found == generalize_s2t.end()) ;
        if(MPI_rank == 0)
          cerr << "forms targets: " << found->second << endl ;
        debugout << "forms targets: " << found->second << endl ;
      }
      if(promote_source.inSet(*vi)) {
        if(MPI_rank == 0)
          cerr << "\tit is also in the chain of promote rules " ;
        debugout << "\tit is also in the chain of promote rules " ;
        found = promote_s2t.find(*vi) ;
        FATAL(found == promote_s2t.end()) ;
        if(MPI_rank == 0)
          cerr << "forms targets: " << found->second << endl ;
        debugout << "forms targets: " << found->second << endl ;
      }
      if(rename_source.inSet(*vi)) {
        if(MPI_rank == 0)
          cerr << "\tit is also in the chain of rename rules " ;
        debugout << "\tit is also in the chain of rename rules " ;
        found = rename_s2t.find(*vi) ;
        FATAL(found == rename_s2t.end()) ;
        if(MPI_rank == 0)
          cerr << "forms targets: " << found->second << endl ;
        debugout << "forms targets: " << found->second << endl ;
      }
    }

    // promoted variables should not be renamed
    variableSet promote_vars = v.get_promote_target_vars() ;
    for(vi=promote_vars.begin();vi!=promote_vars.end();++vi) {
      found = rename_s2t.find(*vi) ;
      if(found != rename_s2t.end()) {
        variableSet rt = found->second ;
        if(MPI_rank == 0) 
          cerr << "WARNING: promoted variable: " << *vi
               << " is renamed to: " << rt << endl ;
        debugout << "WARNING: promoted variable: " << *vi
                 << " is renamed to: " << rt << endl ;
        
      }
    }

    // input variables should not be renamed
    for(vi=input.begin();vi!=input.end();++vi) {
      found = rename_s2t.find(*vi) ;
      if(found != rename_s2t.end()) {
        variableSet rt = found->second ;
        if(MPI_rank==0)
          cerr << "WARNING: input variable: " << *vi
               << " is renamed to: " << rt << endl ;
        debugout << "WARNING: input variable: " << *vi
                 << " is renamed to: " << rt << endl ;
      }
      variableSet::const_iterator vi2 ;
      variableSet target ;
      
      target = get_all_recur_vars(promote_s2t,*vi) ;
      for(vi2=target.begin();vi2!=target.end();++vi2) {
        found = rename_s2t.find(*vi) ;
        if(found != rename_s2t.end()) {
          if(MPI_rank == 0)
            cerr << "WARNING: variable " << *vi2
                 << " is promoted from input variable " << *vi
                 << " but is renamed to: " << found->second
                 << endl ;
          debugout << "WARNING: variable " << *vi2
                   << " is promoted from input variable " << *vi
                   << " but is renamed to: " << found->second
                   << endl ;
        }
      }
    }

  }


  /////////////////////////////////////////////////////////////////
  // unitApplyMapVisitor
  /////////////////////////////////////////////////////////////////
  void unitApplyMapVisitor::gather_info(const digraph& gr) {
    digraph grt = gr.transpose() ;
    variableSet allvars = extract_vars(gr.get_all_vertices()) ;
    map<variable,ruleSet> reduce_info ;
    variableSet reduce_vars ;
    
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
        reduce_info[*vi] = use_rules ;
        if(reduction && unit_rule_exists)
          reduce_vars += *vi ;
      }
    }
    if(reduce_vars != EMPTY) {
      all_reduce_vars += reduce_vars ;

      map<variable,ruleSet>::const_iterator xi ;
      for(xi=reduce_info.begin();xi!=reduce_info.end();++xi) {
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
            else
              if(typeid(*join_op) !=
                 typeid(*(ri->get_rule_implP()->get_joiner()))) {
                if(MPI_rank ==0)
                  cerr << "Warning:  Not all apply rules for variable "
                       << xi->first << " have identical join operations!"
                       << endl ;
                debugout << "Warning:  Not all apply rules for variable "
                         << xi->first << " have identical join operations!"
                         << endl ;
              }
#endif
          } else {
            if(MPI_rank == 0)
              cerr << "Warning: reduction variable " << xi->first
                   << " has a non-reduction rule contributing\
 to its computation,"
                   << endl << "offending rule is " << *ri << endl ;
            debugout << "Warning: reduction variable " << xi->first
                     << " has a non-reduction rule contributing\
 to its computation,"
                     << endl << "offending rule is " << *ri << endl ;
          }
        }
        if(join_op == 0) {
          if(MPI_rank == 0)
            cerr << "unable to find any apply rules to complete\
 the reduction defined by rule"
                 << endl
                 << unit_rule << endl ;
          debugout << "unable to find any apply rules to complete\
 the reduction defined by rule"
                   << endl
                   << unit_rule << endl ;
        }
        WARN(join_op == 0) ;
        // fill all the maps
        for(ri=apply_rules.begin();ri!=apply_rules.end();++ri)
          apply2unit[*ri] = unit_rule ;
        reduceInfo[xi->first] = make_pair(unit_rule,join_op) ;
      }
    }
  }
  
  void unitApplyMapVisitor::visit(loop_compiler& lc) {
    gather_info(lc.collapse_gr) ;
    gather_info(lc.advance_gr) ;
  }
  
  void unitApplyMapVisitor::visit(dag_compiler& dc) {
    gather_info(dc.dag_gr) ;
  }
  
  void unitApplyMapVisitor::visit(conditional_compiler& cc) {
    gather_info(cc.cond_gr) ;
  }

  /////////////////////////////////////////////////////////////////
  // allocDelNumReportVisitor
  /////////////////////////////////////////////////////////////////
  void allocDelNumReportVisitor::adNum(const digraph& gr,
                                      int& alloc_num,
                                      int& del_num) {
    ruleSet rules = extract_rules(gr.get_all_vertices()) ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      if(ri->get_info().qualifier() == "ALLOCATE")
        ++alloc_num ;
      if(ri->get_info().qualifier() == "DELETE")
        ++del_num ;
    }
  }
  
  void allocDelNumReportVisitor::visit(loop_compiler& lc) {
    s << "In Loop Node (id = " << lc.cid << ")" << endl ;
    int alloc = 0 ;
    int del = 0 ;
    adNum(lc.collapse_gr,alloc,del) ;
    s << "\tIn the collapse part: " << endl ;
    s << "\t\tThere are: " << alloc << " allocation rule(s)" << endl ;
    s << "\t\tThere are: " << del << " deletion rule(s)" << endl ;

    alloc = 0 ; del = 0 ;
    adNum(lc.advance_gr,alloc,del) ;
    s << "\tIn the advance part: " << endl ;
    s << "\t\tThere are: " << alloc << " allocation rule(s)" << endl ;
    s << "\t\tThere are: " << del << " deletion rule(s)" << endl ;
  }
  
  void allocDelNumReportVisitor::visit(dag_compiler& dc) {
    s << "In Dag Node (id = " << dc.cid << ")" << endl ;
    int alloc = 0 ;
    int del = 0 ;
    adNum(dc.dag_gr,alloc,del) ;
    s << "\tThere are: " << alloc << " allocation rule(s)" << endl ;
    s << "\tThere are: " << del << " deletion rule(s)" << endl ;
  }
  
  void allocDelNumReportVisitor::visit(conditional_compiler& cc) {
    s << "In Conditional Node (id = " << cc.cid << ")" << endl ;
    int alloc = 0 ;
    int del = 0 ;
    adNum(cc.cond_gr,alloc,del) ;
    s << "\tThere are: " << alloc << " allocation rule(s)" << endl ;
    s << "\tThere are: " << del << " deletion rule(s)" << endl ;
  }

  ////////////////////////////////////////////////////////////
  // get_islands
  ////////////////////////////////////////////////////////////
  std::vector<digraph> get_islands(const digraph& gr) {
    digraph::vertexSet vertices = gr.get_all_vertices() ;

    list<digraph::vertexSet> pre_islands ;
    // we get clusters that are reachable by each vertex
    for(digraph::vertexSet::const_iterator vi=vertices.begin();
        vi!=vertices.end();++vi) {
      digraph::vertexSet working ;
      digraph::vertexSet visited ;
      digraph::vertexSet next ;

      working += *vi ;
      while(working != EMPTY) {
        visited += working ;
        for(digraph::vertexSet::const_iterator vi2=working.begin();
            vi2!=working.end();++vi2)
          next += gr[*vi2] ;
        next -= visited ;
        working = next ;
      }

      pre_islands.push_back(visited) ;
    }

    // now we merge them into real islands
    vector<digraph> islands ;
    while(!pre_islands.empty()) {
      digraph::vertexSet island_vertices ;
      island_vertices = pre_islands.front() ;
      pre_islands.pop_front() ;

      bool merge = true ;
      while(merge) {
        merge = false ;
        vector<list<digraph::vertexSet>::iterator> remove ;
        for(list<digraph::vertexSet>::iterator ili=pre_islands.begin();
            ili!=pre_islands.end();++ili) {
          if( (island_vertices & *ili) != EMPTY) {
            island_vertices += *ili ;
            merge = true ;
            remove.push_back(ili) ;
          }
        }
        for(vector<list<digraph::vertexSet>::iterator>::size_type i=0;
            i!=remove.size();++i)
          pre_islands.erase(remove[i]) ;
      }

      digraph one_island = gr.subgraph(island_vertices) ;
      islands.push_back(one_island) ;
    }

    return islands ;
  }

#ifdef DYNAMICSCHEDULING
  /////////////////////////////////////////////////////////////////
  // DynamicKeyspaceVisitor
  /////////////////////////////////////////////////////////////////
  namespace {
    variableSet
    remove_synonym(const variableSet& vs, fact_db& facts) {
      variableSet ret ;
      for(variableSet::const_iterator vi=vs.begin();vi!=vs.end();++vi)
        ret += facts.remove_synonym(*vi) ;

      return ret ;
    }
  }
  void
  DynamicKeyspaceVisitor::replace_compiler(digraph& gr,
                                           rulecomp_map& rcm, int id) {
    // get all the rules inside a graph
    digraph::vertexSet all_vertices = gr.get_all_vertices() ;
    ruleSet all_rules = extract_rules(all_vertices) ;
    // then we'll obtain the input variables to this graph
    // (those variables that are not computed in the graph)
    // and the intra variables to this graph (those variables
    // that are computed in the graph).
    variableSet input_vars = extract_vars(gr.get_source_vertices() -
                                          gr.get_target_vertices()) ;
    variableSet intra_vars = extract_vars(gr.get_source_vertices() &
                                          gr.get_target_vertices()) ;
    // we need to obtain the actual concrete rule impls
    for(ruleSet::const_iterator ri=all_rules.begin();
        ri!=all_rules.end();++ri) {
      if(is_internal_rule(ri))
        continue ;
      // "main" keyspace is currently always setup to
      // be STATIC
      string ks_tag = ri->get_info().rule_impl->get_keyspace_tag() ;
      if(ks_tag == "main")
        continue ;
      // get the corresponding keyspace impl from fact_db
      map<string,KeySpaceP>::const_iterator mi ;
      mi = facts.keyspace.find(ks_tag) ;
      if(mi == facts.keyspace.end()) {
        cerr << "Error: rule: " << *ri << " in Non-exist keyspace: "
             << ks_tag << endl ;
        Loci::Abort() ;
      }
      KeySpaceP kp = mi->second ;
      if(kp->get_dynamism() != DYNAMIC)
        continue ;
      // finally we have a rule in a dynamic keyspace
      // we need to replace the rule compiler with a
      // corresponding dynamic_impl_compiler

      // insertion rule does not need dynamic control reset
      // so we skip it first
      if(ri->get_info().rule_impl->get_rule_class()
         == rule_impl::INSERTION) {
        dynamic_targets += ri->targets() ;
        kp->add_control_vars(ri->targets()) ;
        continue ;
      }

      variableSet rule_sources = ri->sources() ;
      // collect the dynamic rule control map
      for(variableSet::const_iterator vi=rule_sources.begin();
          vi!=rule_sources.end();++vi) {
        drule_ctrl[*vi].insert(ks_tag) ;
        // register the map in the keyspace
        kp->set_dcontrol_map(*vi, *ri) ;
        // get all the rotation variables
        variableSet rot = scheds.try_get_rotations(*vi) ;
        map<variable,variableSet>::const_iterator orvi ;
        orvi = overlap_rotvars.find(*vi) ;
        if(orvi != overlap_rotvars.end())
          rot += orvi->second ;
        // we also need to set rotation variables
        for(variableSet::const_iterator vii=rot.begin();
            vii!=rot.end();++vii) {
          drule_ctrl[*vii].insert(ks_tag) ;
          kp->set_dcontrol_map(*vii, *ri) ;
        }
      }

      
      // deletion and erase rule need to have dynamic
      // control information recorded, but other than
      // this, we don't need to anything to them except
      // for adding their targets to the keyspace
      if(ri->get_info().rule_impl->get_rule_class()
         == rule_impl::DELETION ||
         ri->get_info().rule_impl->get_rule_class()
         == rule_impl::ERASE) {
        dynamic_targets += ri->targets() ;
        // also add targets to the keyspace
        kp->add_control_vars(ri->targets()) ;
        continue ;
      }

      // determine rule's volatile and static sources
      // with regard to this graph
      variableSet static_sources =
        variableSet(rule_sources & input_vars) ;
      variableSet volatile_sources =
        variableSet(rule_sources & intra_vars) ;

      FATAL(rule_sources !=
            variableSet(static_sources + volatile_sources)) ;

      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::APPLY) {
        // get the unit rule associated
        std::map<rule,rule>::const_iterator mi = apply2unit.find(*ri) ;
        if(mi == apply2unit.end()) {
          cerr << "Error: cannot find associated unit rule for "
               << "apply rule: " << *ri << endl ;
          Loci::Abort() ;
        }
        dynamic_apply_compiler* d =
          new dynamic_apply_compiler(*ri, mi->second, kp) ;
        rcm[*ri] = d ;
      } else {
        dynamic_impl_compiler* d =
          new dynamic_impl_compiler(*ri, kp,
                                    volatile_sources, static_sources) ;
        rcm[*ri] = d ;
      }
      // we now collect target variables
      variableSet tunnels = kp->get_tunnels() ;
      tunnels = remove_synonym(tunnels,facts) ;
      
      set<vmap_info>::const_iterator vmsi ;
      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::APPLY) {
        bool cross_space = false ;

        FATAL(ri->targets().size() != 1) ;
        
        vmsi = ri->get_info().desc.targets.begin() ;

        for(size_t i=0;i<vmsi->mapping.size();++i) {
          variableSet m = vmsi->mapping[i] ;
          m = remove_synonym(m,facts) ;

          variableSet k = variableSet(m & tunnels) ;
          if(k != EMPTY) {
            FATAL(k != m) ;
            cross_space = true ;
            break ;
          } 
        }

        if(!cross_space) {
          dynamic_targets += ri->targets() ;
          // also add targets to the keyspace
          kp->add_control_vars(ri->targets()) ;
        }
      } else {
        // we'll need to check to make sure that
        // no target is crossing the keyspace boundary
        for(vmsi=ri->get_info().desc.targets.begin();
            vmsi!=ri->get_info().desc.targets.end();++vmsi) {
          for(size_t i=0;i<vmsi->mapping.size();++i) {
            variableSet m = vmsi->mapping[i] ;
            m = remove_synonym(m,facts) ;
            
            variableSet k = variableSet(m & tunnels) ;
            if(k != EMPTY) {
              FATAL(k != m) ;
              if(Loci::MPI_rank == 0) {
                cerr << "Error: only apply rule is allowed to "
                     << "produce targets in other keyspace!!" << endl ;
                cerr << "Offending rule: " << *ri << endl ;
              }
              Loci::Abort() ;
            } 
          }
          // add the targets
          // however we remove "OUTPUT" from dynamic targets
          variableSet output_vars ;
          for(variableSet::const_iterator vi=vmsi->var.begin();
              vi!=vmsi->var.end();++vi) {
            variable sv(*vi, time_ident()) ;
            if(sv == variable("OUTPUT"))
              output_vars += *vi ;
          }

          variableSet vars = vmsi->var ;
          vars -= output_vars ;
          
          dynamic_targets += vars ;
          kp->add_control_vars(vars) ;
        }
      }
      
      // now we need to collect all the input clone information
      for(vmsi=ri->get_info().desc.sources.begin();
          vmsi!=ri->get_info().desc.sources.end();++vmsi) {
        // no mapping, no clone needed
        if(vmsi->mapping.empty())
          continue ;
        
        bool cross_space = false ;
        // evaluate the first mapping block
        variableSet mappings = vmsi->mapping[0] ;
        variableSet mappings_unique = remove_synonym(mappings,facts) ;
        variableSet cross = variableSet(mappings_unique & tunnels) ;
        if(cross != EMPTY) {
          FATAL(cross != mappings_unique) ;
          cross_space = true ;
        }
        // evaluate the rest of the blocks
        for(size_t i=1;i<vmsi->mapping.size();++i) {
          mappings = vmsi->mapping[i] ;
          mappings_unique = remove_synonym(mappings,facts) ;
          if(cross_space) {
            // shadow clone
            for(variableSet::const_iterator
                  vi=mappings_unique.begin();
                vi!=mappings_unique.end();++vi) {
              // create the shadow first
              storeRepP srp = facts.get_variable(*vi) ;
              srp = srp->new_store(EMPTY)->thaw() ;
              kp->create_shadow(*vi, srp) ;
              shadow_clone[*vi].insert(ks_tag) ;
            }
          } else {
            // self clone
            for(variableSet::const_iterator vi=mappings_unique.begin();
                vi!=mappings_unique.end();++vi) {
              map<variable,string>::const_iterator
                mi = self_clone.find(*vi) ;
              if(mi == self_clone.end()) {
                self_clone[*vi] = ks_tag ;
              } else {
                // see if conflict exist
                if(mi->second != ks_tag) {
                  cerr << "Error: self_clone data error! variable "
                       << *vi << " has self clone in keyspace: "
                       << mi->second << " and keyspace: " << ks_tag
                       << endl ;
                  Loci::Abort() ;
                }
              }
            }
          }
          // evaluate if cross space
          if(!cross_space) {
            cross = variableSet(mappings_unique & tunnels) ;
            if(cross != EMPTY) {
              FATAL(cross != mappings_unique) ;
              cross_space = true ;
            }
          }
        } // end of mapping
        // now evaluate the var
        mappings = vmsi->var ;
        mappings_unique = remove_synonym(mappings,facts) ;
        if(cross_space) {
          // shadow clone
          for(variableSet::const_iterator
                vi=mappings_unique.begin();
              vi!=mappings_unique.end();++vi) {
            // create the shadow first
            storeRepP srp = facts.get_variable(*vi) ;
            srp = srp->new_store(EMPTY)->thaw() ;
            kp->create_shadow(*vi, srp) ;
            shadow_clone[*vi].insert(ks_tag) ;
          }
        } else {
          // self clone
          for(variableSet::const_iterator vi=mappings_unique.begin();
              vi!=mappings_unique.end();++vi) {
            map<variable,string>::const_iterator
              mi = self_clone.find(*vi) ;
            if(mi == self_clone.end()) {
              self_clone[*vi] = ks_tag ;
            } else {
              // see if conflict exist
              if(mi->second != ks_tag) {
                cerr << "Error: self_clone data error! variable "
                     << *vi << " has self clone in keyspace: "
                     << mi->second << " and keyspace: " << ks_tag
                     << endl ;
                Loci::Abort() ;
              }
            }
          }
        }
      } // end of source

      for(vmsi=ri->get_info().desc.constraints.begin();
          vmsi!=ri->get_info().desc.constraints.end();++vmsi) {
        // no mapping, no clone needed
        if(vmsi->mapping.empty())
          continue ;
        
        bool cross_space = false ;
        // evaluate the first mapping block
        variableSet mappings = vmsi->mapping[0] ;
        variableSet mappings_unique = remove_synonym(mappings,facts) ;
        variableSet cross = variableSet(mappings_unique & tunnels) ;
        if(cross != EMPTY) {
          FATAL(cross != mappings_unique) ;
          cross_space = true ;
        }
        // evaluate the rest of the blocks
        for(size_t i=1;i<vmsi->mapping.size();++i) {
          mappings = vmsi->mapping[i] ;
          mappings_unique = remove_synonym(mappings,facts) ;
          if(cross_space) {
            // shadow clone
            for(variableSet::const_iterator
                  vi=mappings_unique.begin();
                vi!=mappings_unique.end();++vi) {
              // create the shadow first
              storeRepP srp = facts.get_variable(*vi) ;
              srp = srp->new_store(EMPTY)->thaw() ;
              kp->create_shadow(*vi, srp) ;
              shadow_clone[*vi].insert(ks_tag) ;
            }
          } else {
            // self clone
            for(variableSet::const_iterator vi=mappings_unique.begin();
                vi!=mappings_unique.end();++vi) {
              map<variable,string>::const_iterator
                mi = self_clone.find(*vi) ;
              if(mi == self_clone.end()) {
                self_clone[*vi] = ks_tag ;
              } else {
                // see if conflict exist
                if(mi->second != ks_tag) {
                  cerr << "Error: self_clone data error! variable "
                       << *vi << " has self clone in keyspace: "
                       << mi->second << " and keyspace: " << ks_tag
                       << endl ;
                  Loci::Abort() ;
                }
              }
            }
          }
          // evaluate if cross space
          if(!cross_space) {
            cross = variableSet(mappings_unique & tunnels) ;
            if(cross != EMPTY) {
              FATAL(cross != mappings_unique) ;
              cross_space = true ;
            }
          }
        } // end of mapping
        // now evaluate the var
        mappings = vmsi->var ;
        mappings_unique = remove_synonym(mappings,facts) ;
        if(cross_space) {
          // shadow clone
          for(variableSet::const_iterator
                vi=mappings_unique.begin();
              vi!=mappings_unique.end();++vi) {
            // create the shadow first
            storeRepP srp = facts.get_variable(*vi) ;
            srp = srp->new_store(EMPTY)->thaw() ;
            kp->create_shadow(*vi, srp) ;
            shadow_clone[*vi].insert(ks_tag) ;
          }
        } else {
          // self clone
          for(variableSet::const_iterator vi=mappings_unique.begin();
              vi!=mappings_unique.end();++vi) {
            map<variable,string>::const_iterator
              mi = self_clone.find(*vi) ;
            if(mi == self_clone.end()) {
              self_clone[*vi] = ks_tag ;
            } else {
              // see if conflict exist
              if(mi->second != ks_tag) {
                cerr << "Error: self_clone data error! variable "
                     << *vi << " has self clone in keyspace: "
                     << mi->second << " and keyspace: " << ks_tag
                     << endl ;
                Loci::Abort() ;
              }
            }
          }
        }
      } // end of constraints

      for(vmsi=ri->get_info().desc.targets.begin();
          vmsi!=ri->get_info().desc.targets.end();++vmsi) {
        // no mapping, no clone needed
        if(vmsi->mapping.empty())
          continue ;
        
        bool cross_space = false ;
        // evaluate the first mapping block
        variableSet mappings = vmsi->mapping[0] ;
        variableSet mappings_unique = remove_synonym(mappings,facts) ;
        variableSet cross = variableSet(mappings_unique & tunnels) ;
        if(cross != EMPTY) {
          FATAL(cross != mappings_unique) ;
          cross_space = true ;
        }
        // evaluate the rest of the blocks
        for(size_t i=1;i<vmsi->mapping.size();++i) {
          mappings = vmsi->mapping[i] ;
          mappings_unique = remove_synonym(mappings,facts) ;
          if(cross_space) {
            // shadow clone
            for(variableSet::const_iterator
                  vi=mappings_unique.begin();
                vi!=mappings_unique.end();++vi) {
              // create the shadow first
              storeRepP srp = facts.get_variable(*vi) ;
              srp = srp->new_store(EMPTY)->thaw() ;
              kp->create_shadow(*vi, srp) ;
              shadow_clone[*vi].insert(ks_tag) ;
            }
          } else {
            // self clone
            for(variableSet::const_iterator vi=mappings_unique.begin();
                vi!=mappings_unique.end();++vi) {
              map<variable,string>::const_iterator
                mi = self_clone.find(*vi) ;
              if(mi == self_clone.end()) {
                self_clone[*vi] = ks_tag ;
              } else {
                // see if conflict exist
                if(mi->second != ks_tag) {
                  cerr << "Error: self_clone data error! variable "
                       << *vi << " has self clone in keyspace: "
                       << mi->second << " and keyspace: " << ks_tag
                       << endl ;
                  Loci::Abort() ;
                }
              }
            }
          }
          // evaluate if cross space
          if(!cross_space) {
            cross = variableSet(mappings_unique & tunnels) ;
            if(cross != EMPTY) {
              FATAL(cross != mappings_unique) ;
              cross_space = true ;
            }
          }
        } // end of mapping
        // now evaluate the var
        mappings = vmsi->var ;
        mappings_unique = remove_synonym(mappings,facts) ;
        if(cross_space) {
          // shadow clone
          for(variableSet::const_iterator
                vi=mappings_unique.begin();
              vi!=mappings_unique.end();++vi) {
            // create the shadow first
            storeRepP srp = facts.get_variable(*vi) ;
            srp = srp->new_store(EMPTY)->thaw() ;
            kp->create_shadow(*vi, srp) ;
            // since we are in the target chain
            // it is therefore NOT necessary to declare
            // this variable as a shadow clone because
            // this variable is the one to be generated
            // so it is not really a clone or so because
            // we are not pulling data remotely to fill
            // this variable, instead we will be pushing
            // data to remote processes for this variable
            // in an apply rule
          }
        } else {
          // no need to do the self clone part
          // because this is not a clone variable in the target chain
        }
        // end of var processing
      } // end of targets
    }
  }

  void
  DynamicKeyspaceVisitor::replace_compiler(const ruleSet& rs) {
    for(ruleSet::const_iterator ri=rs.begin();ri!=rs.end();++ri) {
      if(is_internal_rule(ri)) {
        cerr << "Error: internal rule: " << *ri
             << " in recursive supernode" << endl ;
        Loci::Abort() ;
      }
      string ks_tag = ri->get_info().rule_impl->get_keyspace_tag() ;
      if(ks_tag != "main") {
        cerr << "Error: recursive rule: " << *ri
             << " in keyspace other than \"main\" is currently NOT"
             << " supported!" << endl ;
        Loci::Abort() ;
      }
    }
  }
  
  void
  DynamicKeyspaceVisitor::visit(loop_compiler& lc) {
    replace_compiler(lc.loop_gr, lc.rule_compiler_map, lc.cid) ;
  }

  void
  DynamicKeyspaceVisitor::visit(dag_compiler& dc) {
    replace_compiler(dc.dag_gr, dc.rule_compiler_map, dc.cid) ;
  }

  void
  DynamicKeyspaceVisitor::visit(conditional_compiler& cc) {
    replace_compiler(cc.cond_gr, cc.rule_compiler_map, cc.cid) ;
  }

  void
  DynamicKeyspaceVisitor::visit(impl_recurse_compiler& irc) {
    replace_compiler(irc.get_rules()) ;
  }

  void
  DynamicKeyspaceVisitor::visit(recurse_compiler& rc) {
    replace_compiler(rc.get_rules()) ;
  }
#endif

  ////////////////////////////////////////////////////////////
  //                        The-End                         //
  ////////////////////////////////////////////////////////////
} // end of namespace Loci

