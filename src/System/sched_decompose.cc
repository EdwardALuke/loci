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
//#define VERBOSE
#include "sched_tools.h"
#include "distribute.h"


using std::map ;
using std::vector ;
using std::set ;
using std::list ;
using std::cout ;
using std::cerr ;
using std::endl ;
using std::pair ;
using std::make_pair ;

namespace Loci {
  void cleanup_component(const digraph &dg, digraph::vertexSet &component) {
    // Remove variables from component if they have references outside
    // the component
    variableSet cvars = extract_vars(component) ;
    digraph::vertexSet remove_vars ;
    digraph gr = dg ;
    digraph grt = gr.transpose() ;
    for(variableSet::const_iterator vi = cvars.begin();
        vi!=cvars.end();
        ++vi) {
      const int id = (*vi).ident() ;
      if((grt[id]-component) != EMPTY
         || (gr[id]-component) != EMPTY) {
        remove_vars += id ; 
        
      }
    }

    component -= remove_vars ;
#ifdef VERBOSE
    variableSet remove_vset = extract_vars(remove_vars) ;
    debugout << "remove_vars = " << remove_vset ;
    debugout <<endl ;
#endif
  }

  inline multiLevelGraph create_mlg(digraph dg,
                                    digraph::vertexSet sources,
                                    digraph::vertexSet targets) {
    digraph::vertexSet grvtx = dg.get_all_vertices() ;
    grvtx -= sources ;
    grvtx -= targets ;
    return multiLevelGraph(dg,grvtx) ;
  }

  /*! prepare_diagraph() adds edges for reduction:
    unit_var->unit_rule, and unit_var->apply_rules
    so that unit rules will stay in components and will be connected to apply rules
  */
  void prepare_digraph(digraph &dg) {
    digraph::vertexSet grvtx = dg.get_all_vertices() ;
    ruleSet all_rules = extract_rules(grvtx) ;
    ruleSet::const_iterator ri ;
    variableSet reduce_vars ;
    for(ri=all_rules.begin();ri!=all_rules.end();++ri) {
      // collect information about unit rules
      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT) {
        if(ri->targets().size() != 1) {
          cerr << "unit rule must have only one variable as target: "
               << endl << *ri << endl ;
          Loci::Abort() ;
        }
        variable unit_var = *(ri->targets().begin()) ;
        if(reduce_vars.inSet(unit_var)) {
          cerr << "only one unit rule may be defined for reduction variable "
               << unit_var << endl ;
          Loci::Abort() ;
        }
        reduce_vars += unit_var ;

        // Make unit rule loop on itself so that it is included in
        // components as we decompose the graph (the unit rule
        // needs to stay with its apply rules.
        dg.add_edge(unit_var.ident(),(*ri).ident()) ;

        // Now we make all apply rules depend explicitly on the unit rule
        // This will make sure that when we begin to schedule that we
        // always call the unit rule before any of the apply rules.
        // We also do a few sanity checks here as well.
        ruleSet apply_rules = extract_rules(dg.transpose()[unit_var.ident()]) ;
        apply_rules -= *ri ;
        for(ruleSet::const_iterator rii=apply_rules.begin();
            rii!=apply_rules.end();
            ++rii) {
          if(rii->get_info().rule_impl->get_rule_class() != rule_impl::APPLY) {
            cerr << "pointwise rule is producing reduction variable: "
                 << unit_var << endl ;
            cerr << "This rule should be an apply rule. (offending rule is):"
                 << endl << *rii << endl ;
            Loci::Abort() ;
          }
          if(rii->targets().size() != 1) {
            cerr << "apply rule should apply to only one reduction variable"
                 << endl << "offending rule is:" << endl << *rii << endl ;
            Loci::Abort() ;
          }
          // Here we add the dependency
          dg.add_edge(unit_var.ident(),(*rii).ident()) ;
        }
      }

    }
    
  }

  
  /*!post_process_mlg() removes edges for reduction:
  unit_var->unit_rule, and unit_var->apply_rules
  adds edges: unit_rule->apply_rules
  and removes dangling vertices
  */
  void post_process_mlg(multiLevelGraph::subGraph &sg) {
    digraph &gr = sg.gr ;
    ruleSet gr_rules = extract_rules(sg.graph_v) ;
    ruleSet::const_iterator ri ;
    for(ri=gr_rules.begin();ri!=gr_rules.end();++ri) {
      if((ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT)) {
        variable unit_var = *(ri->targets().begin()) ;
        gr.remove_edge(unit_var.ident(),ri->ident()) ;
        ruleSet reduce_rules = extract_rules(gr.transpose()[unit_var.ident()]) ;
        reduce_rules -= *ri ;
        ruleSet::const_iterator rri ;
        for(rri = reduce_rules.begin();rri != reduce_rules.end();++rri) {
          gr.remove_edge(unit_var.ident(),rri->ident()) ;
          gr.add_edge(ri->ident(),rri->ident())  ;
        }
      }
    }
    // remove extra dangling nodes
    gr.remove_dangling_vertices() ;
  }
                                             

  digraph::vertexSet sort_time_hierarchy(multiLevelGraph &mlg) {
    digraph::vertexSet new_vertices ;
    
    int toplevel = mlg.toplevel ;
    new_vertices += toplevel ;
    multiLevelGraph::subGraph &sg = *mlg.find(toplevel) ;
    digraph &dg = sg.gr ;
    variableSet vars = extract_vars(sg.graph_v) ;

    map<time_ident,digraph::vertexSet> time_sort_vertices ;

    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi)
      time_sort_vertices[vi->time()] += (*vi).ident() ;
    ruleSet rules = extract_rules(dg.get_all_vertices()) ;

    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri) {
      // sort rules based on time level
      const time_ident rule_time =
        ri->type()==rule::COLLAPSE?ri->source_time():ri->target_time() ;
      time_sort_vertices[rule_time] += (*ri).ident() ;
    }

    // compute a depth first ordering of time levels
    vector<time_ident> time_stack,time_sequence ;
    time_stack.push_back(time_ident()) ;
    for(;!time_stack.empty();) {
      time_ident ctime = time_stack.back() ;
      time_stack.pop_back() ;
      time_sequence.push_back(ctime) ;
      const vector<time_ident> &children = ctime.children()  ;
      for(vector<time_ident>::const_iterator
            i=children.begin();i!=children.end();++i)
        time_stack.push_back(*i) ;
    }

    // proceed from the bottom of the time hierarchy tree, extracting each
    // time level and constructing a super node based on that time level.
    digraph main_graph = dg ;
    for(vector<time_ident>::reverse_iterator ti = time_sequence.rbegin();
        ti!=time_sequence.rend();
        ++ti) {
      cleanup_component(dg,time_sort_vertices[*ti]) ;
      if(time_sort_vertices[*ti] != EMPTY) {
        if(*ti != time_ident()) {
          int new_node = mlg.mksnode(toplevel,time_sort_vertices[*ti]) ;
          multiLevelGraph::subGraph *g = mlg.find(new_node) ;
          if(extract_vars((g->incoming_v & g->outgoing_v)) != EMPTY) {
            cerr << __LINE__ << "recursive supernode found for variables " <<
              extract_vars((g->incoming_v & g->outgoing_v)) << endl ;
          }
          new_vertices += new_node;
          // The new node that we create will be in the level
          // above our time iteration, so insert it in that
          // iteration set.
          time_sort_vertices[ti->parent()] += new_node ;

          ruleSet temp = extract_rules(time_sort_vertices[*ti]) ;
        }
      }
    }
    return new_vertices ;
  }

  digraph::vertexSet partition_loops(multiLevelGraph &mlg, int supernode,
                                     digraph::vertexSet &loops,
                                     digraph::vertexSet &recursive) {
    multiLevelGraph::subGraph &sg = *mlg.find(supernode) ;

    // First search for looping and collapse rules in this subgraph
    digraph::vertexSet looping, collapse, unit, apply, new_nodes ;

    ruleSet all_rules = extract_rules(sg.graph_v) ;
    for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) {
      if(ri->get_info().qualifier() == "looping")
        looping += ri->ident() ;
      if(ri->type() == rule::COLLAPSE) 
        collapse += ri->ident() ;
      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT)
        unit += ri->ident() ;
      if(ri->get_info().rule_impl->get_rule_class() == rule_impl::APPLY)
        apply += ri->ident() ;
    }

    WARN(looping == EMPTY && collapse != EMPTY) ;
    if(looping == EMPTY && collapse != EMPTY) {
      ruleSet collapse_rules = extract_rules(collapse) ;
      cerr << "A collapse rule exists without a well formed loop" << endl ;
      cerr << "collapse = " << collapse_rules << endl ;
    }
    WARN(looping.size() > 1) ;

    if(looping.size() > 1) {
      cerr << extract_rules(looping) << endl ;
    }


    // Now find any recursive rules and form supernodes for these components
    // We first have to temporarily remove the looping rules so that they
    // don't appear to be big collections of recursive rules.  
    digraph tempg = sg.gr ;
    if(looping != EMPTY)
      tempg.remove_vertices(looping) ;


    // We use a component sort to identify recursive rules.
    vector<digraph::vertexSet> recurse_components =
      component_sort(tempg).get_components() ;
    vector<digraph::vertexSet>::const_iterator ci ;
   

    // Here we check all of the components (larger than one vertex).  We
    // eliminate the reduction rules (which can't be part of a recursive block)
    for(ci=recurse_components.begin();ci!=recurse_components.end();++ci)
      if(ci->size() > 1) {
        ruleSet rr = extract_rules(*ci) ;
        digraph::vertexSet sn = rr ;
        if((sn & unit) == EMPTY && (sn & apply) == EMPTY) {
          // Found a recursive rule, so make a supernode for it.
          int newnode = mlg.mksnode(supernode,sn) ;

          new_nodes += newnode ;
          recursive += newnode ;
          // Edit graph so that the resulting supernode is no longer recursive
          digraph::vertexSet recursive_bits = sg.gr[newnode] & sg.gr.transpose()[newnode] ;

          // remove the recursive dependency and replace with a dependency
          // that will make sure that the recursive rule block is executed last.
          digraph::vertexSet::const_iterator cvi ;
          digraph::vertexSet recurse_depend ;
          for(cvi=recursive_bits.begin();cvi!=recursive_bits.end();++cvi) {
            recurse_depend += sg.gr.transpose()[*cvi] ;
            sg.gr.remove_edge(*cvi,newnode) ;
          }
          recurse_depend -= newnode ;
            
          sg.gr.add_edges(recurse_depend,newnode);
        } else {
          sn -= unit ;
          sn -= apply ;
          if(sn != EMPTY) {
            cerr << "reductions cannot be part of a recursive block!" << endl ;
            cerr << "offending rules = " << endl ;
            cerr << rr << endl ;
          }
        }
      }

    if(looping == EMPTY)
      return EMPTY ;

#ifdef VERBOSE
    debugout << "looping = " << extract_rules(looping) << endl ;
#endif
    // Add these edges so that the collapse rules will be
    // included in the loop component.
    sg.gr.add_edges(collapse,*(looping.begin())) ;

    vector<digraph::vertexSet> components =
      component_sort(sg.gr).get_components() ;

    // Identify components, Also identify candidates that
    // may be added to components if needed.
    digraph::vertexSet candidates ;
    vector<digraph::vertexSet> comp ;
    for(ci=components.begin();ci!=components.end();++ci)
      if(ci->size() > 1)
        comp.push_back(*ci) ;
      else {
        int id = *(*ci).begin() ;
        if(id<0) {
          rule r(id) ;
          if(r.type() == rule::INTERNAL) {
            if(r.get_info().qualifier() != "generalize" &&
               r.get_info().qualifier() != "promote")
              candidates += *ci ;
          } else
            candidates += *ci ;
        } else
          candidates += *ci ;
      }

    // Include any candidates that are uniquely referenced by a component
    // into the component.
    digraph grt = sg.gr.transpose() ;
    for(size_t i=0;i<comp.size();++i) {
      digraph::vertexSet cs = comp[i] ;
      digraph::vertexSet imageSet ;
      digraph::vertexSet::const_iterator vi ;
      for(vi=cs.begin();vi!=cs.end();++vi) 
        imageSet += grt[*vi] ;
      // imageSet is the set of all rules/variables that have inputs
      // to this component
      imageSet -= cs ;
      imageSet &= candidates ;

      digraph::vertexSet found ;
      do {
        // Search through the imageSet for components are exclusively used
        // by this component.
        found = EMPTY ;
        for(vi=imageSet.begin();vi!=imageSet.end();++vi) {
          digraph::vertexSet touch = sg.gr[*vi] ;
          if((touch & cs)==touch) // vertex only accessing component, include
            found += *vi ;
        }
        cs += found ; // Add found vertices to component
        for(vi=found.begin();vi!=found.end();++vi) 
          imageSet += grt[*vi] ;
        imageSet -= cs ;
        imageSet -= found ;
        imageSet &= candidates ;
        // Repeat until no new vertices are found
      } while(found != EMPTY) ;


      // remove any variables in the set that access verticies outside
      // of the new candidate component
      variableSet vars = extract_vars(cs) ;
      digraph::vertexSet except ;
      for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
        digraph::vertexSet s = grt[vi->ident()] ;
        if(!((s & cs) == s))
          except += vi->ident() ;
      }
      cs -= except ;

      // Add these vertices of exclusive reference to the component
      // and remove them from the list of candidates
      cs -= comp[i] ;
      candidates -= cs ;
#ifdef VERBOSE
      debugout << "adding variables to loop " << extract_vars(cs) << endl ;
      debugout << "adding rules to loop:" << endl << extract_rules(cs) << endl ;
#endif
      comp[i] += cs ;
    }

    // Make sure any iterating variable that is computed will be included in
    // the component, Otherwise things can get confused as iteration variables
    // rotate
    digraph::vertexSet add_to_looping ;
    for(size_t i=0;i<comp.size();++i) {
      if((comp[i] & looping) != EMPTY) {
        ruleSet loopr = extract_rules(looping) ;
        variableSet loopvars = loopr.begin()->targets() ;
        variableSet::const_iterator vi ;
        for(vi = loopvars.begin();vi!=loopvars.end();++vi) {
          ruleSet rs = extract_rules(grt[vi->ident()]-comp[i]) ;
          ruleSet::const_iterator ri ;
          for(ri = rs.begin();ri != rs.end();++ri) {
            if(ri->type() != rule::INTERNAL) {
              add_to_looping += ri->ident() ;
              
            }
          }
        }
      }
    }
    for(size_t i=0;i<comp.size();++i) {
      if((comp[i] & looping) != EMPTY) {
        comp[i] += add_to_looping ;
      } else {
        comp[i] -= add_to_looping ;
      }
    }
      
    
    // Now make subgraphs for these looping components
    for(ci=comp.begin();ci!=comp.end();++ci)
      if((*ci & looping) != EMPTY) {
        digraph::vertexSet component_parts = *ci ;
        //          cleanup_component(sg.gr,component_parts) ;
        int newnode = mlg.mksnode(supernode,component_parts) ; 
        multiLevelGraph::subGraph *g = mlg.find(newnode) ;
        if(extract_vars((g->incoming_v & g->outgoing_v)) != EMPTY) {
          cerr << __LINE__ << "recursive supernode found for variables " <<
            extract_vars((g->incoming_v & g->outgoing_v)) << endl ;
          variableSet s = extract_vars((g->incoming_v & g->outgoing_v)) ;
          multiLevelGraph::subGraph *t = mlg.find(supernode) ;
          variableSet::const_iterator ii ;
          for(ii=s.begin();ii!=s.end();++ii) {
            ruleSet r = extract_rules(t->gr[ii->ident()]) ;
            digraph grt = t->gr.transpose() ;
            cerr << *ii << " --> rules " << r << endl ;
            ruleSet rt = extract_rules(grt[ii->ident()]) ;
            cerr << *ii << " <-- rules " << rt << endl ;
          }
        
        }
        new_nodes += newnode ;
        loops += newnode ;
        break ;
      }

    return new_nodes ;
  }


  digraph::vertexSet decompose_conditional_rules(multiLevelGraph &mlg,
                                                 int supernode)
  {
    multiLevelGraph::subGraph &sg = *mlg.find(supernode) ;


    map<variable,ruleSet> cond_rules ;
    ruleSet all_rules = extract_rules(sg.graph_v) ;

    digraph::vertexSet looping ;
    for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) {
      if(ri->get_info().qualifier() == "looping")
        looping += ri->ident() ;
      if(ri->get_info().desc.conditionals != EMPTY) {
        variableSet conds = ri->get_info().desc.conditionals ;
        if(conds.size() != 1) {
          cerr << "improper rule: " << *ri << endl
               << "Rule Improperly specifies more than one conditional" << endl ;
          cerr << "Error not recoverable." << endl ;
          Loci::Abort() ;
        }
        variable cond = *(conds.begin()) ;
        cond_rules[cond] += *ri ;
      }
    }
    // If no conditional rules return 
    if(cond_rules.begin() == cond_rules.end()) {
      return EMPTY ;
    }

    digraph::vertexSet new_rules ;

    digraph::vertexSet::const_iterator vi ;

    digraph tmpgr = sg.gr ;

    tmpgr.remove_vertices(looping) ;

    // Remove rule-rule dependencies.  These are not real variable
    // dependencies, but rather created by rename rules.
    ruleSet graph_rules = extract_rules(tmpgr.get_all_vertices()) ;

    for(ruleSet::const_iterator ri = graph_rules.begin();
        ri!=graph_rules.end();
        ++ri) {
      int id = (*ri).ident() ;
      digraph::vertexSet rm_v = tmpgr[id] & interval(UNIVERSE_MIN,-1) ;
      if(rm_v != EMPTY) {
        for(digraph::vertexSet::const_iterator ii = rm_v.begin() ;
            ii != rm_v.end() ;
            ++ii) {
#ifdef VERBOSE
          debugout << "removing edge " << *ri << "," << rule(*ii) << endl  ;
#endif
          tmpgr.remove_edge(id,*ii) ;
        }
      }
    }

    
    // We need to make sure that apply rules work as a group.  We use
    // a component sort to group them together and remove them
    // from the graph.  We will replace them after we've decided
    // what will be in the conditional component.
    vector<digraph::vertexSet> cs = component_sort(tmpgr).get_components() ;
    map<int,digraph::vertexSet> cm ;
    digraph::vertexSet incm ;

    int nid = 0 ;
    for(size_t i=0;i<cs.size();++i)
      if(cs[i].size() > 1) {
        char buf[512] ;
        snprintf(buf,512,"__internal__@Vertex_%d",nid) ;
        nid++ ;
        variable tmp(buf) ;
        int new_vertex = tmp.ident() ; //tmpgr.max_vertex()+1 ;
        digraph::vertexSet in ;
        digraph::vertexSet out ;
        cm[new_vertex] = cs[i] ;
        cs[i] = extract_rules(cs[i]) ;
        for(vi=cs[i].begin();vi!=cs[i].end();++vi) {
          out += tmpgr[*vi] ;
          in += tmpgr.transpose()[*vi] ;
        }
        incm += new_vertex ;

        out -= cs[i] ;
        in -= cs[i] ;

        in -= out ;// Disable recursive loops

        tmpgr.remove_vertices(cs[i]) ;
        tmpgr.add_edges(new_vertex,out) ;
        tmpgr.add_edges(in,new_vertex) ;
      }

#ifdef PRIORITY_FIX_EXPERIMENTAL
    // Find Priority variable clusters
    std::map<variable, variableSet> prio_vars ;

    ruleSet rset = extract_rules(tmpgr.get_all_vertices()) ;

    for(ruleSet::const_iterator ri = rset.begin();ri!=rset.end();++ri) {
      if(ri->get_info().qualifier() == "priority") {
        variable v = *ri->sources().begin() ;
        while(v.get_info().priority.size() != 0)
          v = v.drop_priority() ;
        prio_vars[v] += ri->sources() ;
        prio_vars[v] += ri->targets() ;
      }
    }

    std::map<variable, variableSet>::const_iterator mvi ;
    for(mvi=prio_vars.begin();mvi!=prio_vars.end();++mvi) {
      digraph::vertexSet  prio_cluster ;
      digraph tr = tmpgr.transpose() ;
      variableSet::const_iterator vi ;
      for(vi=mvi->second.begin(); vi!=mvi->second.end();++vi) {
        prio_cluster += vi->ident() ;
        prio_cluster += tr[vi->ident()] ;
      }
      char buf[512] ;
      snprintf(buf,512,"__internal_prio__@Vertex_%d",nid) ;
      nid++ ;
      variable tmp(buf) ;
      int new_vertex = tmp.ident() ; //tmpgr.max_vertex()+1 ;
      digraph::vertexSet in ;
      digraph::vertexSet out ;
      cm[new_vertex] = prio_cluster ;
      prio_cluster = extract_rules(prio_cluster) ;
      digraph::vertexSet::const_iterator Vi ;
      for(Vi=prio_cluster.begin();Vi!=prio_cluster.end();++Vi) {
        out += tmpgr[*Vi] ;
        in += tmpgr.transpose()[*Vi] ;
      }
      incm += new_vertex ;

      out -= prio_cluster ;
      in -= prio_cluster ;

      in -= out ;// Disable recursive loops
      
      tmpgr.remove_vertices(prio_cluster) ;
      tmpgr.add_edges(new_vertex,out) ;
      tmpgr.add_edges(in,new_vertex) ;
    }      

#endif
    // Loop over each conditional variable and find the part of the graph
    // that exclusively connects to the conditional components.  All of
    // these rules can be treated together when evaluating the conditional
    map<variable,ruleSet>::const_iterator mi ;
    for(mi=cond_rules.begin();mi!= cond_rules.end();++mi) {

      digraph tr = tmpgr.transpose() ;
      digraph::vertexSet component ;
      // We start with the conditional rules
      digraph::vertexSet working = mi->second ;

#ifdef VERBOSE
      debugout  << "computing component for conditional " <<  mi->first <<endl ;
#endif
      digraph::vertexSet candidates ;
      while(working != EMPTY) {
        // Repeatedly search for new vertices in the graph that are exclusively
        // associated with the conditional rules.  As we find them add them
        // to the current set of component rules
        component += working ;

        // First find candidates by following all edges that lead to the
        // current working set.
        for(vi=working.begin();vi!=working.end();++vi)
          candidates += tr[*vi] ;

        // Create a new working set by searching through the candidates to
        // find those verticies that only refer to vertices already found
        // in the component.

        // We don't need to try things already in the component
        candidates -= component ;
        // We don't search the conditional variable as it needs to be
        // outside the conditional block.
        candidates -= mi->first.ident() ;
#ifdef VERBOSE
        debugout << "candidate vars = " << extract_vars(candidates) <<endl ;
        debugout << "candidate rules = " << extract_rules(candidates) << endl ;
#endif

        working = EMPTY ;
        for(vi=candidates.begin();vi!=candidates.end();++vi) {
          if((tmpgr[*vi] - component) == EMPTY) {
            working += *vi ;
          } 
        }
      }

      bool clean = false ;
      do {
        variableSet component_vars = extract_vars(component) ;
        component_vars -= incm ;// Don't include virtual node "variables"
        clean = true ;
        digraph::vertexSet cleanSet ;
        for(variableSet::const_iterator v=component_vars.begin();
            v != component_vars.end();
            ++v) {
          if((tmpgr[v->ident()] - component) != EMPTY  ||
             (tr[v->ident()] - component) != EMPTY) {
            // Variable should not be inside component, remove all things that
            // lead to this variable from the component
            clean = false ;
            working = EMPTY ;
            working += v->ident() ;
            cleanSet += v->ident() ;
            digraph::vertexSet vtouch ;
            while(working != EMPTY) {
              vtouch += working ;
              digraph::vertexSet candidates = EMPTY ;
              digraph::vertexSet::const_iterator vv ;
              for(vv=working.begin();vv!=working.end();++vv)
                candidates += tr[*vv] ;
              working = candidates ;
              working -= vtouch ;
            }
            cleanSet += vtouch ;
          }
        }
        ruleSet r = extract_rules(cleanSet & component) ;
#ifdef VERBOSE
        if(r!=EMPTY)
          debugout << "cleaning rules " << r << endl ;
#endif

        component -= cleanSet ;
      } while(!clean) ;


#ifdef VERBOSE
      debugout << "variables for conditional " << mi->first
               << " are " << extract_vars(component) << endl ;
#endif
      
#define NEWWAY
#ifdef NEWWAY

      digraph::vertexSet shared_cond_rules = tmpgr[mi->first.ident()] ;
      shared_cond_rules &= component ;

      vector<digraph::vertexSet> sub_comp ;

      ruleSet cond_rules = extract_rules(shared_cond_rules) ;
      bool collapse = false ;
      for(ruleSet::const_iterator ri=cond_rules.begin();
          ri!=cond_rules.end();++ri) 
        if(ri->type() == rule::COLLAPSE)
          collapse = true ;
        

        
      if(collapse) { // If this is part of the collapse, keep grouped
        digraph::vertexSet compcm = component & incm ;
        component -= incm ;
        for(vi=compcm.begin();vi!=compcm.end();++vi)
          component += cm[*vi] ;

        sub_comp.push_back(component) ;
      } else {
        // Find out if there are multiple independent sub components of the
        // conditional.  If we break them into separate supernodes, then
        // memory allocation should be more efficient.  First identify what
        // part of the component is associated with each conditonal rule:
        digraph sgr = tmpgr.subgraph(component) ;
        digraph sgrt = sgr.transpose() ;
        for(digraph::vertexSet::const_iterator vi=shared_cond_rules.begin();
            vi!=shared_cond_rules.end();++vi) {
          digraph::vertexSet vtmp ;
          vtmp += *vi ;
          digraph::vertexSet vtmp2 = visit_vertices(sgrt,vtmp) + vtmp ;
          digraph::vertexSet compcm = vtmp2 & incm ;
          vtmp2 -= incm ;
          for(digraph::vertexSet::const_iterator vi2=compcm.begin();
              vi2!=compcm.end();++vi2)
            vtmp2 += cm[*vi2] ;
              
          sub_comp.push_back(vtmp2) ;
        }
      } 

      //Find input variables that are computed in this level that are not
      // part of the component.  Add them, because if multiple conditionals
      // share the same input, then it is appropriate for them to be grouped
      // together.  If the computed inputs don't overlap, the components
      // should be independent to improve memory allocation opportunities
      
      digraph grt = sg.gr.transpose() ;
      for(size_t i=0;i<sub_comp.size();++i) {
        variableSet comp_vars ;
        ruleSet rs = extract_rules(sub_comp[i]) ;
        for(ruleSet::const_iterator ri=rs.begin();ri!=rs.end();++ri)
          comp_vars += ri->sources() ;
        comp_vars -= mi->first ;

        variableSet cv ;
        // find computed variables
        for(variableSet::const_iterator ii=comp_vars.begin();
            ii!=comp_vars.end();++ii) {
          ruleSet rs2 = extract_rules(grt[ii->ident()]) ;

          if((rs2&sub_comp[i])==EMPTY) {
            for(ruleSet::const_iterator ri=rs2.begin();
                ri != rs2.end();++ri) {
              if(ri->type() != rule::INTERNAL)
                cv += *ii ;
            }
          }
        }
        sub_comp[i] += cv ;
      }

      // Now determine if there is any overlap between subcomponents:
      digraph::vertexSet overlap,comunion ;
      for(size_t i=0;i<sub_comp.size();++i) {
        overlap += sub_comp[i] & comunion ;
        comunion += sub_comp[i] ;
      }

      // Create conditional components for each non-overlapping component
      // put remaining components here
      digraph::vertexSet remain ;

      for(size_t i=0;i<sub_comp.size();++i)
        if((sub_comp[i]&overlap) == EMPTY) {
          // Make sure that we don't go outside of the current subgraph
          sub_comp[i] &=sg.graph_v ;
          // Remove variables from component if they have references outside
          // the component
          cleanup_component(sg.gr,sub_comp[i]) ;

          int new_node =  mlg.mksnode(supernode,sub_comp[i],mi->first) ;
          multiLevelGraph::subGraph *g = mlg.find(new_node) ;
          if(extract_vars((g->incoming_v & g->outgoing_v)) != EMPTY) {
            cerr << __LINE__ << "recursive supernode found for variables " <<
              extract_vars((g->incoming_v & g->outgoing_v)) << endl ;
          }
          new_rules += new_node ;
        } else
          remain += sub_comp[i] ;


      if(remain != EMPTY) {
        // Make sure that we don't go outside of the current subgraph
        remain &= sg.graph_v ;

        // Remove variables from component if they have references outside
        // the component
        cleanup_component(sg.gr,remain) ;

        int new_node =  mlg.mksnode(supernode,remain,mi->first) ;
        multiLevelGraph::subGraph *g = mlg.find(new_node) ;
        if(extract_vars((g->incoming_v & g->outgoing_v)) != EMPTY) {
          cerr << __LINE__ << "recursive supernode found for variables " <<
            extract_vars((g->incoming_v & g->outgoing_v)) << endl ;
        }
        new_rules += new_node ;
      }

#else
      
      // If any apply rules were included, add them back into the current
      // set of components

      digraph::vertexSet compcm = component & incm ;
      component -= incm ;
      for(vi=compcm.begin();vi!=compcm.end();++vi)
        component += cm[*vi] ;


      // Make sure that we don't go outside of the current subgraph
      component &= sg.graph_v ;

      // Remove variables from component if they have references outside
      // the component
      cleanup_component(sg.gr,component) ;

      int new_node =  mlg.mksnode(supernode,component,mi->first) ;

      multiLevelGraph::subgraph *g = mlg.find(new_node) ;
      if(extract_vars((g->incoming_v & g->outgoing_v)) != EMPTY) {
        cerr << __LINE__ << "recursive supernode found for variables " <<
          extract_vars((g->incoming_v & g->outgoing_v)) << endl ;
      }
      new_rules += new_node ;
#endif
      
    }

    return new_rules ;
  }
  
  decomposed_graph::decomposed_graph(digraph dg,
                                     digraph::vertexSet sources,
                                     digraph::vertexSet targets) {
    prepare_digraph(dg) ;
    mlg = create_mlg(dg,sources,targets) ;

#ifdef DEBUG
    digraph::vertexSet all_vertices = mlg.find(mlg.toplevel)->graph_v ;
#endif

    digraph::vertexSet new_vertices = sort_time_hierarchy(mlg) ;
    digraph::vertexSet::const_iterator vi ;
    digraph::vertexSet looping_components ;
    
    for(vi=new_vertices.begin();vi!=new_vertices.end();++vi) {
      looping_components += partition_loops(mlg,*vi,loops,recursive) ;
    }
    digraph::vertexSet working = mlg.subgraphs ;
    digraph::vertexSet conditional_components;
    for(vi=working.begin();vi!=working.end();++vi) {
      conditional += decompose_conditional_rules(mlg,*vi) ;
    }
    
    for(vi=mlg.subgraphs.begin();vi!=mlg.subgraphs.end();++vi)
      post_process_mlg(*mlg.find(*vi)) ;
  
#ifdef DEBUG
    digraph::vertexSet visit ;
    for(vi=mlg.subgraphs.begin();vi!=mlg.subgraphs.end();++vi) {
      digraph::vertexSet vl = mlg.find(*vi)->graph_v ;
      if(extract_rules(vl & visit) != EMPTY) {
        cerr << "overlapping graphs in mlg" << endl ;
        cerr << "overlapping variables = " << extract_vars(vl&visit)<< endl ;
        cerr << "overlapping rules = " << extract_rules(vl&visit)<<endl ;
      }
      visit += vl ;
    }

    if((all_vertices - visit) != EMPTY) {
      cerr << "some verticies lost in mlg generation" << endl ;
      cerr << "lost variables = " << extract_vars(all_vertices-visit) << endl ;
      cerr << "lost rules = " << extract_rules(all_vertices-visit) << endl ;
    }
#endif
    
  }

}
