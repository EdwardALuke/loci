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
#include "loci_globs.h"

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include<map>
using std::map ;
#include <list>
using std::list ;
#include <utility>
using std::pair ;
using std::make_pair ;
using std::min ;

using std::cerr ;
using std::endl ;

// all the searching and graph
// editing algorithm for chomping 
namespace Loci {

  // show some debug information for the chomp chain discover algorithm
  //#define CHOMP_DEBUG
  
  //////////////////////////////////////////////////////////////
  // some helper functions 
  //////////////////////////////////////////////////////////////
  namespace {
    // given a contrete rule, return all the sources
    // that have maps
    inline variableSet map_sources(const rule& r) {
      const rule_impl::info& rinfo = r.get_info().desc ;
      std::set<vmap_info>::const_iterator si ;
      variableSet sources ;
      for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
        if(!si->mapping.empty())
          sources += si->var ;
      }
      return sources ;
    }
    
    // given a contrete rule, return all the targetss
    // that have maps
    inline variableSet map_targets(const rule& r) {
      const rule_impl::info& rinfo = r.get_info().desc ;
      std::set<vmap_info>::const_iterator si ;
      variableSet targets ;
      for(si=rinfo.targets.begin();si!=rinfo.targets.end();++si) {
        if(!si->mapping.empty())
          targets += si->var ;
      }
      return targets ;
    }

    // return the vertexSet of the rule's sources, targets
    // and the rule itself
    inline digraph::vertexSet get_rule_vertexSet(const rule& r) {
      digraph::vertexSet ret ;
      variableSet vars ;
      vars += r.sources() ;
      vars += r.targets() ;
      ret = get_vertexSet(vars) ;
      ret += r.ident() ;
      return ret ;
    }

    digraph::vertexSet get_ruleSet_vertexSet(const ruleSet& rules) {
      digraph::vertexSet ret ;
      for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
        ret += get_rule_vertexSet(*ri) ;
      }
      return ret ;
    }
    
    inline bool thread_rule(const rule& r) {
      return (r.get_rule_implP()->thread_rule()&&
              !r.get_rule_implP()->dynamic_schedule_rule()) ;
    }

    // given a rule, return if there are any OUTPUT
    // in the targets
    bool has_output_in_targets(const rule& r) {
      variableSet targets = r.targets() ;
      variable output("OUTPUT") ;
      for(variableSet::const_iterator vi=targets.begin();
          vi!=targets.end();++vi) {
        variable var_stationary(*vi,time_ident()) ;
        if(var_stationary == output)
          return true ;
      }

      return false ;
    }

    // find any chains that contains v and then
    // remove those chains from the list and
    // patch them back to the search list
    // HOWEVER, BECAUSE OF THE RECENT REVISION OF ALGORITHMS,
    // ESPECIALLY WITH THE MORE ROBUST gen_tmp_graph FUNCTION,
    // THIS NOTIFY ACTION MAY NOT BE NECESSARY. We'll revise this
    // in the future.
    void notify_existchains(list<chomp_chain>& result,
                            list<variable>& chomp_vars_order,
                            variableSet& valid_chomp_vars,
                            const variable& v) {
      vector<list<chomp_chain>::iterator> remove ;
      variableSet patch_back ;
      for(list<chomp_chain>::iterator li=result.begin() ;
          li!=result.end();++li) {
        digraph& gr = li->first ;
        variableSet internal_vars =
          extract_vars(gr.get_source_vertices() & gr.get_target_vertices()) ;
        if(internal_vars.inSet(v)) {
          remove.push_back(li) ;
          patch_back += li->second ;
        }
      }
      for(vector<list<chomp_chain>::iterator>::size_type i=0;
          i!=remove.size();++i)
        result.erase(remove[i]) ;

      valid_chomp_vars += patch_back ;
      for(variableSet::const_iterator vi=patch_back.begin();
          vi!=patch_back.end();++vi)
        chomp_vars_order.push_back(*vi) ;
    }

    // get all chomping candidates in sources
    // and targets of a graph
    inline variableSet get_chompCandInST(const digraph& gr,
                                         const variableSet& chompCand) {
      digraph::vertexSet sources = gr.get_source_vertices()
        - gr.get_target_vertices() ;
      digraph::vertexSet targets = gr.get_target_vertices()
        - gr.get_source_vertices() ;
      variableSet stv = extract_vars(sources+targets) ;

      return variableSet(stv & chompCand) ;
    }

    // test for self-cycle and non-chompable
    // internal variables in a graph
    // test_gr is a subgraph of gr
#ifdef CHOMP_DEBUG
    namespace {
      int problem_code = 0 ;
    } 
#endif
    inline bool has_problem(const digraph& gr, const digraph& test_gr,
                            const variableSet& validChompVars,
                            const variableSet& merged_vars) {
      digraph::vertexSet sources = test_gr.get_source_vertices()
        - test_gr.get_target_vertices() ;
      digraph::vertexSet targets = test_gr.get_target_vertices()
        - test_gr.get_source_vertices() ;
      
      digraph::vertexSet internal = test_gr.get_source_vertices()
        & test_gr.get_target_vertices() ;
      
      FATAL(internal !=
            (test_gr.get_all_vertices() - sources - targets)) ;
      
      variableSet internal_vars = extract_vars(internal) ;
      variableSet problem_vars =
        variableSet(internal_vars - validChompVars) ;

      digraph::vertexSet merged = get_vertexSet(merged_vars) ;

#ifdef CHOMP_DEBUG
      if(has_path(gr,targets,sources))
        problem_code = 1 ;
      else if(problem_vars != EMPTY)
        problem_code = 2 ;
#endif
      return (has_path(gr,targets,sources) ||
              (has_path(gr,sources,merged) && has_path(gr,merged,targets)) ||
              (problem_vars != EMPTY)) ;
    }
    
    // this function computes the internal variables
    // (variables that are not the source and target) of a digraph
    inline variableSet get_internal_vars(const digraph& gr) {
      return extract_vars(gr.get_source_vertices() &
                          gr.get_target_vertices()) ;
    }

    // this function generates a temporary graph for chomping
    // it first merges in a chomping candidate and then brings
    // in any necessary rules and variables

    // because when we merge new vertices into the existing graph
    // the source & target vertices of the existing graph could
    // become internal, and thus we need to bring in all the rules
    // that connect to these "new" internal variables. And this
    // is a repeated process, we will stop only if no new internal
    // variables are genrated
    digraph gen_tmp_graph(const digraph& cur_gr,
                          const digraph::vertexSet& new_vertices,
                          const digraph& gr) {
      digraph grt = gr.transpose() ;
      digraph::vertexSet graph_vertices = cur_gr.get_all_vertices() ;
      graph_vertices += new_vertices ;
      variableSet cur_internal_vars = get_internal_vars(cur_gr) ;
      variableSet new_internal_vars ;
      // we repeat until no new internal vertices are generated
      digraph new_gr ;
      while(true) {
        new_gr = gr.subgraph(graph_vertices) ;
        new_internal_vars = get_internal_vars(new_gr) ;
        variableSet diff = variableSet(new_internal_vars - cur_internal_vars) ;
        if(diff == EMPTY)
          break ;
        // then bring in any relevant rules
        ruleSet addon_rules ;
        for(variableSet::const_iterator vi=diff.begin();
            vi!=diff.end();++vi) {
          addon_rules += extract_rules(gr[vi->ident()]) ;
          addon_rules += extract_rules(grt[vi->ident()]) ;
        }
        addon_rules -= extract_rules(new_gr.get_all_vertices()) ;
        graph_vertices += get_ruleSet_vertexSet(addon_rules) ;
        cur_internal_vars = new_internal_vars ;
      }
      return new_gr ;
    }

    // this function merges two chomping graphs and
    // returns the new resulting graph


    //////////////////////////////////////////////////
    // CHOMP_OPT currently has a logical bug inside 
    // More work is needed to fix it in the future. 
    // The bug is that in the first opt pass, when the algorithm
    // merges two chomp chains, it is possible to include
    // chomped variables in other chains as either source
    // or target variables in the merge chain. This would
    // cause problem as it brings the chomped variable in other
    // chains visible.
    //////////////////////////////////////////////////
    //#define CHOMP_OPT
    //#define CHOMP_OPT_MORE

    // given a digraph and all the chomping candidate
    // variables inside it, this function will return
    // all chomping chains by merging variables together
    void get_chomp_chains(list<chomp_chain>& result,
                          const variableSet& chomp_vars,
                          const digraph& gr) {
      if(chomp_vars.size() == 0)
        return ;
      // we first get a topological order
      // of the chomp_vars 
      vector<digraph::vertexSet> gr_order =
        component_sort(gr).get_components() ;
      list<variable> chomp_vars_order ;
      for(vector<digraph::vertexSet>::const_iterator vi=gr_order.begin();
          vi!=gr_order.end();++vi) {
        FATAL(vi->size() != 1) ;
        int_type fv = *(vi->begin()) ;
        if(fv >= 0){
          variable v(fv) ;
          if(chomp_vars.inSet(v))
            chomp_vars_order.push_back(v) ;
        }
      }

      digraph grt = gr.transpose() ;
      variableSet valid_chomp_vars = chomp_vars ;
      variableSet merged_vars ;
      list<chomp_chain> temp_result ;
      
      while(!chomp_vars_order.empty()) {
        // here is the basis, we start by selecting
        // a single chompable variable and form the
        // smallest chomp chain
        variable begin = chomp_vars_order.front() ;
        ruleSet all_rules = extract_rules(gr[begin.ident()] +
                                          grt[begin.ident()]) ;
        variableSet sub_chomp_vars ;
        sub_chomp_vars += begin ;
        // generate an initial graph for this single starting var
        digraph cur_chomp_gr =
          gr.subgraph(get_ruleSet_vertexSet(all_rules)) ;
        // however, we need to make sure that this is a complete graph.
        // that is, we need to make sure any internal variables have
        // all connected vertices included. This is an issue because for
        // complex rules, it may be possible to have multiple targets,
        // and one of these target may be chosen as our "begin" variable.
        // for complex graphs, other targets may become internal
        // variables in the graph. If we don't include all their
        // connected vertices, it would cause problem later.
        while(true) {
          // get current internal variables
          variableSet cur_internal_vars =
            get_internal_vars(cur_chomp_gr) ;
          digraph::vertexSet graph_vertices =
            cur_chomp_gr.get_all_vertices() ;
          // then bring in any relevant rules
          ruleSet add_rules ;
          for(variableSet::const_iterator vi=cur_internal_vars.begin();
              vi!=cur_internal_vars.end();++vi) {
            add_rules += extract_rules(gr[vi->ident()]) ;
            add_rules += extract_rules(grt[vi->ident()]) ;
          }
          add_rules -= extract_rules(graph_vertices) ;
          if(add_rules == EMPTY)
            break ;
          graph_vertices += get_ruleSet_vertexSet(add_rules) ;
          cur_chomp_gr = gr.subgraph(graph_vertices) ;
        }
        
        // remove the first variable from the list
        chomp_vars_order.pop_front() ;

#ifdef CHOMP_DEBUG
        cerr << "Start merging from var: " << begin << endl ;
#endif
        bool merge = true ;
        while(merge) {
          merge = false ;
          vector<list<variable>::iterator> remove ;
          for(list<variable>::iterator lvi=chomp_vars_order.begin();
              lvi!=chomp_vars_order.end();++lvi){
#ifdef CHOMP_DEBUG
            cerr << "--testing merging on var: " << *lvi << endl ;
#endif
            ruleSet rvs_rules = extract_rules(grt[lvi->ident()]) ;
            ruleSet fwd_rules = extract_rules(gr[lvi->ident()]) ;
            ruleSet reachable_rules = ruleSet(rvs_rules + fwd_rules) ;

            if( (all_rules & reachable_rules) != EMPTY) {
              // test for self-cycle
              digraph::vertexSet new_vertices =
                get_ruleSet_vertexSet(reachable_rules) ;
              digraph test_gr =
                gen_tmp_graph(cur_chomp_gr,new_vertices,gr) ;

              if(!has_problem(gr,test_gr,valid_chomp_vars,merged_vars)){
                sub_chomp_vars += *lvi ;
                all_rules += reachable_rules ;
                cur_chomp_gr = test_gr ;
                remove.push_back(lvi) ;
                merge = true ;
#ifdef CHOMP_DEBUG
                cerr << "  --merge accepted" << endl ;
#endif
              }
#ifdef CHOMP_DEBUG
              else {
                cerr << "  --merge rejected (due to constraints" ;
                if(problem_code == 1)
                  cerr << " [path problem]" ;
                else if(problem_code == 2)
                  cerr << " [internval var problem]" ;
                cerr << ")" << endl ;
              }
#endif
            }
#ifdef CHOMP_DEBUG
            else
              cerr << "  --merge rejected (due to no shared vars)" << endl ;
#endif
          } // end-of-for(chomp_vars_order)
          // remove those merged variables from the list
          for(vector<list<variable>::iterator>::size_type i=0;
              i!=remove.size();++i)
            chomp_vars_order.erase(remove[i]) ;
#ifdef CHOMP_DEBUG
          if(merge) {
            cerr << endl ;
            cerr << "--restart merge testing" << endl ;
          }
#endif
        }//end-of-while(merge)

        // if the chain just contains one chompable variable
        // we need to test for self-cycle, if has, it is
        // then discarded
        if(sub_chomp_vars.size() == 1) {
          if(has_problem(gr,cur_chomp_gr,valid_chomp_vars,merged_vars)) {
            variable v = *(sub_chomp_vars.begin()) ;
            valid_chomp_vars -= v ;

#ifdef CHOMP_DEBUG
            cerr << "Valid_chomp_vars removed (single): " << v << endl ;
#endif
            chomp_vars_order.remove(v) ;
            // we need to notify existed chains of this change
            // this variable might be silently included
            // into other chains because it may be target of
            // other chompable variables and thus is already included
            // as an internal variable (as a by-product) in other chain.
            // HOWEVER, BECAUSE OF THE RECENT REVISION OF ALGORITHMS,
            // ESPECIALLY WITH THE MORE ROBUST gen_tmp_graph FUNCTION,
            // THIS NOTIFY ACTION MAY NOT BE NECESSARY. We'll revise this
            // in the future.
            notify_existchains(result,chomp_vars_order,
                               valid_chomp_vars,v) ;
            continue ;
          }
        }
        // we need to remove any chompable variables
        // that are the sources or targets of this
        // sub_gr. They are the cause of the self-cycles.
        // They cannot be merged into other chains because
        // they are now visible outside of this chain
        variableSet intersection =
          get_chompCandInST(cur_chomp_gr,chomp_vars) ;
        if(intersection != EMPTY) {
          for(variableSet::const_iterator rmvi=intersection.begin();
              rmvi!=intersection.end();++rmvi)
            chomp_vars_order.remove(*rmvi) ;
          valid_chomp_vars -= intersection ;

#ifdef CHOMP_DEBUG
          cerr << "Valid_chomp_vars removed (intersection): "
               << intersection << endl ;
#endif
        }

        valid_chomp_vars -= sub_chomp_vars ;
        merged_vars += sub_chomp_vars ;

#ifdef CHOMP_DEBUG
        cerr << "Valid_chomp_vars removed (sub_chomp): "
             << sub_chomp_vars << endl ;
#endif
        
#ifdef CHOMP_OPT
        temp_result.push_back(make_pair(cur_chomp_gr,sub_chomp_vars)) ;
#else
        result.push_back(make_pair(cur_chomp_gr,sub_chomp_vars)) ;
#endif
      } // end-of-while(!chomp_vars_order.empty())
#ifdef CHOMP_OPT      
      // here we do an optimization phase
      // we try to merge the formed chains
      // and hope we can get back some of the chomping
      // variables that are excluded in the above merging
      list<chomp_chain> temp_result2 ;
      while(!temp_result.empty()) {
        chomp_chain begin = temp_result.front() ;
        digraph new_gr = begin.first ;

        temp_result.pop_front() ;

        bool merge = true ;
        while(merge) {
          merge = false ;
          vector<list<chomp_chain>::iterator> remove ;
          for(list<chomp_chain>::iterator li=temp_result.begin();
              li!=temp_result.end();++li) {
            digraph this_gr = li->first ;
            variableSet shared =
              extract_vars(new_gr.get_all_vertices() &
                           this_gr.get_all_vertices()) ;
            if( (shared & chomp_vars) != EMPTY) {
              digraph test_gr = merge_2_graphs(new_gr,this_gr,gr) ;

              if(!has_problem(gr,test_gr,chomp_vars,merged_vars)) {
                new_gr = test_gr ;
                remove.push_back(li) ;
                merge = true ;
              }
            }
          } // end for
          // remove those merged chains from the list
          for(vector<list<variable>::iterator>::size_type i=0;
              i!=remove.size();++i)
            temp_result.erase(remove[i]) ;
        } // end while(merge)

        variableSet chomp_vars_in_chain =
          extract_vars(new_gr.get_source_vertices() &
                       new_gr.get_target_vertices()) ;
        temp_result2.push_back(make_pair(new_gr,chomp_vars_in_chain)) ;
      } // end of while(!temp_result.empty())
      
      for(list<chomp_chain>::iterator li=temp_result2.begin();
          li!=temp_result2.end();++li) {
#ifdef CHOMP_OPT_MORE
        // In this section of chomp optimization, we try to
        // recover any possible chomp vars in a single chomp chain.
        // At the beginning, we know what the chomp candidates are.
        // However at the above steps, they may be discarded due to
        // various problems. Here for any chomp candiates that are
        // in the chain and are discarded and are NOT in any other
        // chains, we try to reclaim them again. If no problems
        // occur, we will accept them.
        
        // it is important to use reference(&) here,
        // because we want to dynamically update temp_result2
        // so that in the following "gcio iter" we can have
        // correct result
        digraph& new_gr = li->first ;
        variableSet& chomp_vars_in_chain = li->second ;
        
        variableSet chompCandInThis =
          get_chompCandInST(li->first,chomp_vars) ;

        variableSet chompCandInOthers ;

        // gcio iter
        for(list<chomp_chain>::iterator li2=temp_result2.begin();
            li2!=temp_result2.end();++li2)
          if(li2 != li)
            chompCandInOthers += get_chompCandInST(li->first, chomp_vars) ;

        for(variableSet::const_iterator vi=chompCandInThis.begin();
            vi!=chompCandInThis.end();++vi) {
          if(!chompCandInOthers.inSet(*vi)) {
            ruleSet rvs_rules = extract_rules(grt[vi->ident()]) ;
            ruleSet fwd_rules = extract_rules(gr[vi->ident()]) ;
            ruleSet reachable_rules = ruleSet(rvs_rules + fwd_rules) ;
            digraph::vertexSet new_vertices =
              get_ruleSet_vertexSet(reachable_rules) ;
            
            digraph test_gr = gr.subgraph(new_vertices +
                                          new_gr.get_all_vertices()) ;
            variableSet new_chomp_vars_in_chain = chomp_vars_in_chain ;
            new_chomp_vars_in_chain += *vi ;
            if(!has_problem(gr,test_gr,chomp_vars,merged_vars)) {
              new_gr = test_gr ;
              chomp_vars_in_chain = new_chomp_vars_in_chain ;
            }
          }
        }
#endif
        // finally we accept the chain
        result.push_back(*li) ;
      }// end of for(temp_result2)
#endif
    }// end of the get_chomp_chains function

  } // end of unamed namespace

  ///////////////////////////////////////////////////////////////
  // chompPPVisitor
  ///////////////////////////////////////////////////////////////
  //#define DISABLE_APPLY
  
  chompPPVisitor::chompPPVisitor(fact_db& fd,
                                 const map<int,variableSet>& rot_vt,
                                 const map<int,variableSet>& lsharedt,
                                 const variableSet& rv)
    :facts(fd),rename_vars(rv) {
    map<int,variableSet>::const_iterator mi ;
    for(mi=rot_vt.begin();mi!=rot_vt.end();++mi)
      rotate_vars += mi->second ;
    for(mi=lsharedt.begin();mi!=lsharedt.end();++mi)
      loop_shared_vars += mi->second ;
  }

  void chompPPVisitor::discover(const digraph& gr) {
    digraph grt = gr.transpose() ;
    // first we need to get all the contrete rules in the graph
    ruleSet rules = extract_rules(gr.get_all_vertices()) ;
    ruleSet remove ;
    // variables that are sources or targets of internal rules
    variableSet intstv ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      if(is_internal_rule(*ri)) {
        remove += *ri ;
        intstv += ri->sources() ;
        intstv += ri->targets() ;
      }

    rules -= remove ;

    // first get all possible variables that can be chomped
    variableSet allvars, remaining_vars, apply_targets ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      allvars += ri->sources() ;
      allvars += ri->targets() ;
      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::APPLY)
        apply_targets += ri->targets() ;
    }
    remaining_vars = allvars ;
    // intstv cannot be chomped for sure, we remove them
    bad_vars += intstv ;
    // targets of apply rule cannot be chomped
    bad_vars += apply_targets ;

    // then we collect relations that have maps involved...
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      bad_vars += map_sources(*ri) ;
      bad_vars += map_targets(*ri) ;
    }

    // remove the bad_vars discovered so far...
    remaining_vars -= bad_vars ;

    // then we get all the remaining variables that are not suitable
    // for chomping, i.e. they are not STORE, they are in the loop
    // rotate list, loop shared variables, and variables that
    // don't have outgoing edges, and variables that connect
    // to any internal rules or rename variables
    // or variables that are generated in more that one levels
    // and variables that connect to unit rules (since unit rule can
    // only be executed once).
    for(variableSet::const_iterator vi=remaining_vars.begin();
        vi!=remaining_vars.end();++vi) {
      
      if(seen_vars.inSet(*vi)) {
        bad_vars += *vi ;
        continue ;
      }
      
      storeRepP srp = facts.get_variable(*vi) ;

      if(!isSTORE(srp)) {
        bad_vars += *vi ;
        continue ;
      }
        
      if(rotate_vars.inSet(*vi)) {
        bad_vars += *vi ;
        continue ;
      }
      if(loop_shared_vars.inSet(*vi)) {
        bad_vars += *vi ;
        continue ;
      }
      
      if(rename_vars.inSet(*vi)) {
        bad_vars += *vi ;
        continue ;
      }
        
      digraph::vertexSet next_vertices = gr[vi->ident()] ;
      if(next_vertices == EMPTY) {
        bad_vars += *vi ;
        continue ;
      }

      digraph::vertexSet next_vertices_t = grt[vi->ident()] ;
      if(next_vertices_t == EMPTY) {
        bad_vars += *vi ;
        continue ;
      }

      ruleSet tmp = extract_rules(next_vertices) ;
      for(ruleSet::const_iterator rii=tmp.begin();
          rii!=tmp.end();++rii)
        if(is_internal_rule(*rii) || !thread_rule(*rii) ||
           has_output_in_targets(*rii) ||
#ifdef DISABLE_APPLY             
           rii->get_info().rule_impl->get_rule_class() == rule_impl::APPLY ||
#endif
           rii->get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
        bad_vars += *vi ;
        break ;
      }
        
      tmp = extract_rules(next_vertices_t) ;
      for(ruleSet::const_iterator rii=tmp.begin();
          rii!=tmp.end();++rii)
        if(is_internal_rule(*rii) || !thread_rule(*rii) ||
           has_output_in_targets(*rii)) {
          bad_vars += *vi ;
          break ;
        }
    }
    remaining_vars -= bad_vars ;

    // fill up...
    good_vars += remaining_vars ;
    seen_vars += allvars ;
  }

  void chompPPVisitor::visit(loop_compiler& lc) {
    discover(lc.collapse_gr) ;
    discover(lc.advance_gr) ;
  }

  void chompPPVisitor::visit(dag_compiler& dc) {
    discover(dc.dag_gr) ;
  }

  void chompPPVisitor::visit(conditional_compiler& cc) {
    discover(cc.cond_gr) ;
  }

  ///////////////////////////////////////////////////////////////
  // chompRuleVisitor
  ///////////////////////////////////////////////////////////////
  list<chomp_chain> chompRuleVisitor::find_chain(const digraph& gr) {
    variableSet allvars = extract_vars(gr.get_all_vertices()) ;
    // now we start the searching, we are first very greedy,
    // we group any potential variables that can be chomped,
    // and then later, we perform checkings and necessary
    // adjustments.
    variableSet cand_vars = variableSet(allvars & good_vars) ;

    list<chomp_chain> chomp_chain_list ;
    get_chomp_chains(chomp_chain_list,cand_vars,gr) ;

    // collect all chomped variables
    for(list<chomp_chain>::const_iterator li=chomp_chain_list.begin();
        li!=chomp_chain_list.end();++li)
      all_chomped_vars += li->second ;

    return chomp_chain_list ;    
  } // end-of-find_chain function

  namespace {
    // here is a small function that is used as a predicate
    // to replace "(){}," chars in a chomp rule signature
    // UNUSED
//     inline bool
//     chomp_sig_replace1(char c) {
//       return (c=='(' || c=='{' || c==',') ;
//     }
//     inline bool
//     chomp_sig_replace2(char c) {
//       return (c==')' || c=='}') ;
//     }
    inline bool
    chomp_sig_replace(char c) {
      return !(isalnum(c) || c=='_') ;
    }
  }
  
  // edit the graph to have the chomp node,
  void chompRuleVisitor::edit_gr(digraph& gr,const list<chomp_chain>& ccin,
                                 rulecomp_map& rcm) {
    // First check to see if the chomp chain will produce cycles if introduced
    // into the graph
    list<chomp_chain> cc ;
    for(list<chomp_chain>::const_iterator li=ccin.begin();li!=ccin.end();++li) {
      // Form set of chomp chain variables and rules
      digraph chomp_graph = li->first ;
      variableSet chomp_vars = li->second ;
      digraph::vertexSet chomp_vars_vertices = get_vertexSet(chomp_vars) ;
      digraph::vertexSet all_vertices = chomp_graph.get_all_vertices() ;
      ruleSet all_rules = extract_rules(all_vertices) ;
      digraph::vertexSet rules_vertices = get_vertexSet(all_rules) ;
      // chomp_set is the set of vertices that form the chomp
      digraph::vertexSet chomp_set = rules_vertices + chomp_vars_vertices ;
      
      // Compute the vertices that are outgoing edges from the chomp chain
      digraph::vertexSet out_vertices ;
      digraph::vertexSet::const_iterator ei ;
      for(ei=chomp_set.begin();ei!=chomp_set.end();++ei)
	out_vertices += gr[*ei] ;
      out_vertices -= chomp_set ;
      
      // Now follow outgoing vertices until no new vertices are found
      digraph::vertexSet visit_set = out_vertices ;
      digraph::vertexSet found_set = out_vertices ;
      do {
        digraph::vertexSet new_set ;
        found_set -= chomp_set ;
	for(ei=found_set.begin();ei!=found_set.end();++ei)
	  new_set += gr[*ei] ;
	found_set = new_set - visit_set ;
	visit_set += found_set ;
      } while(found_set!=EMPTY) ;
      
      // Check to see if any outgoing edges led back to chomp, if so then
      // making this chomp a supernode will cause cycle, otherwise
      // if it doesn't, then it is ok to process
      if((visit_set & chomp_set) == EMPTY)
	cc.push_back(*li) ; // OK chomp, schedule it
      else {
        // This chomp could cause a cycle, print diagnostic in debug file
	debugout << "NOTE:  Removing chomp chain for vars = " << chomp_vars
		 << endl ;
	debugout << "cycle variables = " 
		 <<  extract_vars(visit_set & chomp_set)
		 << endl ;
      }
    }

    if(cc.empty())
      return ;
    
    for(list<chomp_chain>::const_iterator li=cc.begin();li!=cc.end();++li) {
      digraph chomp_graph = li->first ;
      variableSet chomp_vars = li->second ;
      digraph::vertexSet chomp_vars_vertices = get_vertexSet(chomp_vars) ;

      digraph::vertexSet all_vertices = chomp_graph.get_all_vertices() ;
      ruleSet all_rules = extract_rules(all_vertices) ;
      digraph::vertexSet rules_vertices = get_vertexSet(all_rules) ;

      // nodes that lead to the constructed super node
      // and nodes that leave the super node
      digraph::vertexSet
        source_vars_vertices = chomp_graph.get_source_vertices() -
        chomp_graph.get_target_vertices() ;
      digraph::vertexSet
        target_vars_vertices = chomp_graph.get_target_vertices() -
        chomp_graph.get_source_vertices() ;

      // check to see if there's any rules inside them, if it is,
      // then they are removed and their sources or targets are added
      ruleSet ruleinsource = extract_rules(source_vars_vertices) ;
      ruleSet ruleintarget = extract_rules(target_vars_vertices) ;
      if(ruleinsource != EMPTY) {
        digraph grt = gr.transpose() ;
        digraph::vertexSet add ;
        for(ruleSet::const_iterator ri=ruleinsource.begin();
            ri!=ruleinsource.end();++ri)
          add += grt[ri->ident()] ;
        source_vars_vertices -= get_vertexSet(ruleinsource) ;
        source_vars_vertices += add ;
      }
      if(ruleintarget != EMPTY) {
        digraph::vertexSet add ;
        for(ruleSet::const_iterator ri=ruleintarget.begin();
            ri!=ruleintarget.end();++ri)
          add += gr[ri->ident()] ;
        target_vars_vertices -= get_vertexSet(ruleintarget) ;
        target_vars_vertices += add ;
      }

      // make a rule for the chomp_graph

      // NOTE: it is possible to create identical rule signatures
      // for the different chomping node (because their sources
      // and targets are exactly the same). therefore we will
      // need to try to create unique rule signature for each
      // chomping node in order to avoid clash in associating
      // chomping node with rule compilers. we do so by creating
      // a rule qualifier that includes the chomped variables.
      // For example, if a chomp node has sources A,B, and targets
      // C,D, and the actual chomped variable is X,Y, then its
      // rule signature will look like:
      // CHOMP_X_Y_CHOMP:C,D <- A,B
      // this guarantees that each chomp rule signature is
      // different from each other
      std::stringstream chomp_qualifier_ss ;
      chomp_qualifier_ss << "CHOMP_" ;
      for(variableSet::const_iterator vi=chomp_vars.begin();
          vi!=chomp_vars.end();++vi)
        chomp_qualifier_ss << *vi << "_" ;
      chomp_qualifier_ss << "CHOMP" ;
      std::string chomp_qualifier = chomp_qualifier_ss.str() ;
      // NOTE: we need to remove all characters inside the
      // chomp_qualifier that are not alpha-numeric and '_'
      // characters. We do so by replacing them with '_'.
      // The reason for doing so is that the Loci expression
      // parser cannot properly extract a type qualifier string
      // with non alpha-numeric characters.
      std::replace_if(chomp_qualifier.begin(),
                      chomp_qualifier.end(),
                      chomp_sig_replace, '_') ;

      rule chomp_rule = create_rule(extract_vars(source_vars_vertices),
                                    extract_vars(target_vars_vertices),
                                    chomp_qualifier) ;
#ifdef CHOMP_DEBUG
      rulecomp_map::iterator rcm_find = rcm.find(chomp_rule) ;
      if(rcm_find != rcm.end()) {
        cerr << "ERROR: Creating identical chomping super rules. "
             << "This should not happen. Chomping processing failed!"
             << endl ;
        cerr << "--Problem occurred on rule chain: " << all_rules << endl ;
        chomp_compiler* chc =
          dynamic_cast<chomp_compiler*>(&(*(rcm_find->second))) ;
        if(chc != 0) {
          cerr << "--This is what we already have: "
               << extract_rules((chc->chomp_graph).get_all_vertices())
               << endl ;
        }
        Loci::Abort() ;
      } else
        rcm[chomp_rule] = new chomp_compiler(chomp_graph,
                                             chomp_vars,apply2unit) ;
#else
      rcm[chomp_rule] = new chomp_compiler(chomp_graph,
                                           chomp_vars,apply2unit) ;
#endif
      
      // the vertices to be taken out
      digraph::vertexSet takeout_vertices =
        chomp_vars_vertices + rules_vertices ;
      
      if(takeout_vertices != (all_vertices -
                              source_vars_vertices - target_vars_vertices)) {
        cerr << "WARNING: inconsistency in chomping graph editing!" << endl ;
        digraph::vertexSet ast ;
        ast = all_vertices - source_vars_vertices - target_vars_vertices ;
        digraph::vertexSet diff ;
        diff = takeout_vertices - ast ;
        if(diff == EMPTY)
          diff = ast - takeout_vertices ;
        variableSet diff_vars = extract_vars(diff) ;
        ruleSet diff_rules = extract_rules(diff) ;
        if(diff_vars != EMPTY)
          cerr << "These vars should be taken out: " << diff_vars << endl ;
        if(diff_rules != EMPTY)
          cerr << "These rules should be taken out: " << diff_rules << endl ;
        cerr << "The chomped vars in this chain are: " << chomp_vars << endl ;

        Loci::Abort() ;
      }
      
      // get other possible nodes (outside of the chomp_graph)
      // that lead to any internal vertices of the chomp graph

      digraph grt = gr.transpose() ;
      for(digraph::vertexSet::const_iterator vi=takeout_vertices.begin();
          vi!=takeout_vertices.end();++vi) {
        source_vars_vertices += grt[*vi] - all_vertices ;
      }

      // get other possible nodes that any internal vertices of
      // the chomp graph also reach
      for(digraph::vertexSet::const_iterator
            vi=takeout_vertices.begin();vi!=takeout_vertices.end();++vi) {
        target_vars_vertices += gr[*vi] - all_vertices ;
      }

      // graph editing
      // take the chomp graph out from the original graph
      gr = gr.subgraph(gr.get_all_vertices() - takeout_vertices) ;

      // edit the graph again
      gr.add_edges(source_vars_vertices, chomp_rule.ident()) ;
      gr.add_edges(chomp_rule.ident(), target_vars_vertices) ;
    }
  }

  namespace {
    // a utility function for defining += for list<chomp_chain>
    void operator+=(list<chomp_chain>& l1, const list<chomp_chain>& l2) {
      for(list<chomp_chain>::const_iterator li=l2.begin();
          li!=l2.end();++li)
        l1.push_back(*li) ;
    }
  }
  
  void chompRuleVisitor::visit(loop_compiler& lc) {
    list<chomp_chain> c ;

    map<rule,rule_compilerP> tmp ;
    
    c = find_chain(lc.collapse_gr) ;
    if(!c.empty()) {
      // all_chains[-lc.cid] = c ;
      // here we sould not distinguish the collapse graph
      // from the advance graph by setting its graph id to
      // a negative number. Because the purpose of the
      // "all_chains" is to record the whole chomping chains
      // in the multi-level graph for later visualization and
      // summarization. For if we adding the negative number
      // as the index, later we will not be able to get a sorted
      // order to print the information. And thus it is
      // important to use the "+=" in assigning the chains
      // to the record, because otherwise we'll erase the
      // existing records.
      all_chains[lc.cid] += c ;
      edit_gr(lc.collapse_gr,c,lc.rule_compiler_map) ;
      edit_gr(lc.loop_gr,c,tmp) ;
    }
    
    c = find_chain(lc.advance_gr) ;
    if(!c.empty()) {
      all_chains[lc.cid] += c ;
      edit_gr(lc.advance_gr,c,lc.rule_compiler_map) ;
      edit_gr(lc.loop_gr,c,tmp) ;
    }
  }

  void chompRuleVisitor::visit(dag_compiler& dc) {
    list<chomp_chain> c ;
    
    c = find_chain(dc.dag_gr) ;
    if(!c.empty()) {
      all_chains[dc.cid] = c ;
      edit_gr(dc.dag_gr,c,dc.rule_compiler_map) ;
    }
  }

  void chompRuleVisitor::visit(conditional_compiler& cc) {
    list<chomp_chain> c ;
    
    c = find_chain(cc.cond_gr) ;
    if(!c.empty()) {
      all_chains[cc.cid] = c ;
      edit_gr(cc.cond_gr,c,cc.rule_compiler_map) ;
    }
  }

  ///////////////////////////////////////////////////////////////
  // compChompVisitor
  ///////////////////////////////////////////////////////////////
  namespace {
    // scheduling utilities for chomping graph
    chompingPrio cpf ;
    graphSchedulerVisitor cgsv(cpf) ;
  }
  
  void compChompVisitor::schedule_chomp(chomp_compiler& chc) {
    //chc.chomp_sched = orderVisitor::order_dag(chc.chomp_graph) ;
    chc.chomp_sched = cgsv.schedule(chc.chomp_graph) ;
  }
  
  void compChompVisitor::compile_chomp(chomp_compiler& chc,
                                       const rulecomp_map& rcm) {
    for(vector<digraph::vertexSet>::const_iterator vi=chc.chomp_sched.begin();
        vi!=chc.chomp_sched.end();++vi) {
      ruleSet rules = extract_rules(*vi) ;
      variableSet vars = extract_vars(*vi) ;
      if(rules == EMPTY && vi+1<chc.chomp_sched.end()) {
        ++vi ;
        vars += extract_vars(*vi) ;
        rules = extract_rules(*vi) ;
      }
      vars &= chc.chomp_vars ;
      digraph chomp_graph_t = chc.chomp_graph.transpose() ;
      variableSet barrier_vars, reduce_vars,singleton_vars,all_vars ;
      
      for(variableSet::const_iterator vii=vars.begin();vii!=vars.end();++vii) {
        ruleSet var_rules = extract_rules(chomp_graph_t[(*vii).ident()]) ;
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
          if( (priority_rule || pointwise) &&
              !recursive && (vii->get_info().name != "OUTPUT")) {
            // Don't use the priority variables for variable barriers
            if(vii->get_info().priority.size() == 0)
              barrier_vars += *vii ;
          }
          if( (priority_rule || pointwise) &&
              recursive && (vii->get_info().name != "OUTPUT")) {
            if(vii->get_info().priority.size() == 0)
              all_vars += *vii ;
          }

          if(reduction && unit_rule_exists)
            reduce_vars += *vii ;
          
          if(singleton) {
            singleton_vars += *vii ;
          }
        } 
      }

      // create a fake rule for using in the record
      rule fake = create_rule(variable("A"),variable("B"),"FAKE");

      all_vars += barrier_vars ;
      if(barrier_vars != EMPTY) {
        chc.
          old_chomp_comp.
          push_back(make_pair(fake,
                              CPTR<Loci::rule_compiler>(new barrier_compiler(barrier_vars)))) ;
	if(duplicate_work)
	  chc.old_barrier_sets.push_back(barrier_vars);
      }
      
      all_vars += singleton_vars ;

      if(singleton_vars != EMPTY)
        chc.
          old_chomp_comp.
          push_back(make_pair(fake,
                              CPTR<Loci::rule_compiler>(new singleton_var_compiler(singleton_vars)))) ;

      all_vars += reduce_vars;

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
          if(sp!=0) {
            if(isPARAMETER(sp)) {
              reduce_var_vector.push_back(xi->first) ;
              unit_rule_vector.push_back(unit_rule) ;
              join_op_vector.push_back(join_op) ;
            } else {
              WARN(!isSTORE(sp)) ;
              chc.
                old_chomp_comp.
                push_back(make_pair(fake,
                                    CPTR<rule_compiler>(
                       new reduce_store_compiler(xi->first,
                                                 unit_rule,
                                                 join_op)))) ;
	      if(duplicate_work) {
		variableSet temp;
		temp += xi->first;
		chc.old_barrier_sets.push_back(temp);
	      }
            }
          }
        }
      }
      if(reduce_var_vector.size() != 0) {
        chc.
          old_chomp_comp.
          push_back(make_pair(fake,
                              CPTR<rule_compiler>(new reduce_param_compiler(reduce_var_vector,
                                                        unit_rule_vector,
                                                        join_op_vector)))) ;
	if(duplicate_work) {
	  variableSet myVars;
	  for(unsigned int i = 0; i < reduce_var_vector.size(); i++)
	    myVars += reduce_var_vector[i];
	  
	  chc.old_barrier_sets.push_back(myVars);
	}
      }
      
      for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri){
        rulecomp_map::const_iterator rmi ;
        rmi = rcm.find(*ri) ;
        FATAL(rmi == rcm.end()) ;
        //rule_compilerP bc = new barrier_compiler(vars) ;
        // the rule compiler rmi->second is only included
        // as a place holder but servers no purpose in this case
        chc.old_chomp_comp.push_back(make_pair(*ri,rmi->second)) ;
      }
    }
  }

  void compChompVisitor::process_rcm(rulecomp_map& rcm) {
    rulecomp_map::iterator ri ;
    for(ri=rcm.begin();ri!=rcm.end();++ri) {
      chomp_compiler* chc =
        dynamic_cast<chomp_compiler*>(&(*(ri->second))) ;
      if(chc != 0) {
        schedule_chomp(*chc) ;
        compile_chomp(*chc,rcm) ;
      }
    }
  }

  void compChompVisitor::visit(loop_compiler& lc) {
    process_rcm(lc.rule_compiler_map) ;
  }
  
  void compChompVisitor::visit(dag_compiler& dc) {
    process_rcm(dc.rule_compiler_map) ;
  }
  
  void compChompVisitor::visit(conditional_compiler& cc) {
    process_rcm(cc.rule_compiler_map) ;
  }

} // end of namespace

