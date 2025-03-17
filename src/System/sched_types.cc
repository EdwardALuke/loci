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
#include "sched_tools.h"

using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

using std::cerr ;
using std::endl ;
using std::ofstream ;

namespace Loci {

  extern std::ofstream debugout ;
    // Create variable types in fact database
  // This is a necessary precursor to binding rules to the database
  void set_var_types(fact_db &facts, const digraph &dg, sched_db &scheds) {

    // Get all the variables and rules represented in the graph
    variableSet all_vars = extract_vars(dg.get_all_vertices()) ;
    ruleSet     all_rules = extract_rules(dg.get_all_vertices()) ;
    // We need the transpose of the graph in order to find the rules that
    // generate a particular variable
    digraph dgt = dg.transpose() ;
    
    // We need to account for iteration variables
    // So we can allocate them (No rules explicitly allocate them)
    variableSet iteration_variables ;

    // We will also identify unused variables taht we don't need to setup
    // types for.
    variableSet unused_variables ;
    
    // First we group variables into equivalence classes
    digraph variable_groupings ;

    for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) {
      if(ri->type() == rule::INTERNAL) {
        if(ri->get_info().qualifier() == "promote" ||
           ri->get_info().qualifier() == "priority") {
          FATAL((ri->sources()).size()!=1)  ;
          FATAL((ri->targets()).size()!=1) ;
          variable source = *((ri->sources()).begin()) ;
          variable target = *((ri->targets()).begin()) ;
          variable_groupings.add_edge(source.ident(),target.ident()) ;
          variable_groupings.add_edge(target.ident(),source.ident()) ;
        } else if(ri->get_info().qualifier() == "looping") {
          // create the iteration variable type
          time_ident tl = ri->source_time() ;
          variable v(tl) ;
          iteration_variables += v ;
          // Here we will identify unused variables!
          variableSet candidates = ri->sources() ;
          candidates += ri->targets() ;
          for(variableSet::const_iterator vi=candidates.begin();
              vi!=candidates.end();++vi) {
            ruleSet rs = extract_rules(dg[vi->ident()]+dgt[vi->ident()]) ;
            if(rs.size() <= 1)
              unused_variables += *vi ;
          }
        }
        if(ri->get_info().qualifier() == "generalize") {
          FATAL((ri->sources()).size()!=1)  ;
          FATAL((ri->targets()).size()!=1) ;
          variable source = *((ri->sources()).begin()) ;
          variable target = *((ri->targets()).begin()) ;
          ruleSet rs = extract_rules(dgt[source.ident()]) ;
	  if(rs == EMPTY) {// If a generalize has no rules to generate its
	    // input, then remove it
	    unused_variables += source ;
            rs = extract_rules(dgt[target.ident()]) ;
            // If Nothing else generates the output then we won't need the
            // result of the generalize either
            
	    if(rs.size() == 2) {// looping + generalize 
              unused_variables += target ;
	    }
          } else {
            variable_groupings.add_edge(source.ident(),target.ident()) ;
            variable_groupings.add_edge(target.ident(),source.ident()) ;
          }
        }
      }
    }

    debugout << "unused variables = " << unused_variables << endl ;
    all_vars -= unused_variables ;
    for(variableSet::const_iterator vi = all_vars.begin();
        vi!=all_vars.end();++vi) {
      variable_groupings.add_edge(vi->ident(),vi->ident()) ;
    }
    
    variableSet given_vars = facts.get_typed_variables() ;

    debugout << "iteration_variables = " << iteration_variables << endl ;

    for(variableSet::const_iterator vi = iteration_variables.begin();
        vi!=iteration_variables.end();
        ++vi) {
      param<int> timevar ;
      *timevar = 0 ;
      storeRepP st = timevar.Rep() ;
      scheds.set_variable_type(*vi,st, facts) ;
    }

    variableSet prev_typed_vars = given_vars ;
    prev_typed_vars += iteration_variables ;
    
    map<variable,variable> var_to_cluster ;
    map<variable,variableSet> cluster_map ;
    
    // Find strongly connected components in variable graph,
    // Select the variable out of this group that should represent the
    // storage of the variable.
   
    vector<digraph::vertexSet> var_clusters =
      component_sort(variable_groupings).get_components() ;

    for(size_t i=0;i<var_clusters.size();++i) {
      variableSet cluster = extract_vars(var_clusters[i]) ;
      if(cluster.size() == 1) {
        variable v = *cluster.begin() ;
        cluster_map[v] = cluster ;
	var_to_cluster[v] = v ;
	continue ;
      }
      FATAL(cluster.size() == 0) ;
      
      variable min_time = *cluster.begin() ;
      if((cluster & given_vars) != EMPTY) {
        variableSet test = cluster ;
        test &= given_vars ;
        if(test.size() > 1) {
          cerr << "rules describe synonyms between variables that exist in fact_db" << endl ;
          cerr << "overlap = " << test << endl ;
          abort() ;
        }
        min_time = *test.begin() ;
      }
      
      variableSet::const_iterator vi ;
      for(vi = cluster.begin();vi!=cluster.end();++vi)
        if(vi->time() < min_time.time() ) {
          min_time = *vi ;
        } else if(vi->time() == min_time.time()) {
          if(min_time.get_info().assign)
            min_time = *vi ;
          if(min_time.get_info().priority.size() != 0)
            min_time = *vi ;
        }
      cluster_map[min_time] = cluster ;
      for(vi=cluster.begin();vi!=cluster.end();++vi) {
        var_to_cluster[*vi] = min_time ;
      }
    }

    
    
    // extract the qualified rules, these are rules that are
    // automatically generated by the system.  Since these rules
    // cannot provide type information directly, they are
    // singled out so that they can be handled as a separate case.
    ruleSet qualified_rules ;
    for(ruleSet::const_iterator ri=all_rules.begin();ri!=all_rules.end();++ri) {
      if(ri->type() == rule::INTERNAL)
        qualified_rules += *ri ;
    }

    vector<pair<variable,variable> > rename_vars ;

    // Loop through the variables and type them when there is an appropriate
    // rule specification.  Deal with typing rules that rename their targets
    // separately
    for(variableSet::const_iterator
          vi=all_vars.begin();vi!=all_vars.end();++vi) {
      // only deal with rules given to the system by the user
      ruleSet rs = extract_rules(dgt[(*vi).ident()] - qualified_rules) ;
      if(rs == EMPTY) {
        continue ;
      }
      // storage for rename specification if one is found
      pair<variable,variable> rename_pair ;
      bool rename = false ;
      
      // Check all rules generating this variable, if the rules specify
      // renaming a variable for the target, then by default the target
      // will be an identical type to that of the previous variable.
      for(ruleSet::const_iterator ri=rs.begin();ri!=rs.end();++ri) {
        set<vmap_info>::const_iterator vmsi ;
        for(vmsi=ri->get_info().desc.targets.begin();
            vmsi!=ri->get_info().desc.targets.end(); ++vmsi)
          if(vmsi->assign.size() != 0) 
            for(size_t i=0;i<vmsi->assign.size();++i) {
              variable new_name = vmsi->assign[i].first ;
	      if(new_name == *vi) {
                variable old_name = vmsi->assign[i].second ;
		if(rename && old_name != rename_pair.second) {
                  cerr << "rule " << *ri << endl 
                       << " is not consistent with other rules w/respect to "
                       << "renaming." << endl ;
                  Loci::Abort() ;
                }
                rename_pair = make_pair(new_name,old_name) ;
                rename = true ;
              }
            }
      }
      
      // If a rename variable was found, save this information for later
      // use and continue ;
      if(rename) {
	variable v1 = rename_pair.first ;
        variable v2 = rename_pair.second ;
        variable c1 = v1 ;
        variable c2 = v2 ;
        if(var_to_cluster.find(v1) != var_to_cluster.end())
          c1 = var_to_cluster[v1] ;
        else {
          var_to_cluster[v1] = v1 ;
          cerr << "undefined var_to_cluster for " << v1 << endl ;
        }
        if(var_to_cluster.find(v2) != var_to_cluster.end())
          c2 = var_to_cluster[v2] ;
        else {
          var_to_cluster[v2] = v2 ;
          cerr << "undefined var_to_cluster for " << v2 << endl ;
        }
       	if(c1 == c2) {
          cerr << "Rules that rename variables make variables distinct that are al so equivalent!" << endl ;
          cerr << "variables = " << v1 << "," << v2
               << ", cluster = " << c1 << "," << c2 << endl ;
        }
        rename_vars.push_back(rename_pair) ;
        continue ;
      }
    }


    digraph rename_graph ;
    map<variable,variableSet>::const_iterator ami ;
    for(ami=cluster_map.begin();ami!=cluster_map.end();++ami) {
      variable cv = ami->first ;
      rename_graph.add_edge(cv.ident(),cv.ident()) ;
    }
    for(size_t i=0;i<rename_vars.size();++i) {
      variable new_var = var_to_cluster[rename_vars[i].first] ;
      variable old_var = var_to_cluster[rename_vars[i].second] ;
      rename_graph.add_edge(old_var.ident(),new_var.ident()) ;
    }  
    
    digraph rename_tr = rename_graph.transpose() ;
    // Perform a topological sort on the interal rule graph
    vector<digraph::vertexSet> topo =
      component_sort(rename_graph).get_components() ;
    for(size_t i=0;i<topo.size();++i) {
      variableSet vs = extract_vars(topo[i]) ;
      if(vs.size() != 1) {
        cerr << "recursion in rename graph on variables " << vs
             << endl ;
        abort() ;
      }
      variable v = *vs.begin() ;
      if(rename_tr[v.ident()].size() > 1) {
        variableSet image= extract_vars(rename_tr[v.ident()]) ;
	image -= v ;
        if(image.size() != 1) {
          cerr << "rename shadow has more than one variable!" << endl ;
          cerr << "variables =" <<  image<< ", rename_variable =" <<  v << endl;
          //          abort() ;
        }
        variable orig = *image.begin() ;
        scheds.alias_variable(orig,v,facts) ;
      } else if(!prev_typed_vars.inSet(v)) {
	variable vget = v ;
        
        ruleSet rs = extract_rules(dgt[v.ident()] - qualified_rules) ;
	if(rs == EMPTY) {
          variableSet cluster = cluster_map[v] ;
          for(variableSet::const_iterator vi = cluster.begin();
              vi!=cluster.end();++vi) {
            rs = extract_rules(dgt[vi->ident()] - qualified_rules) ;
            vget = *vi ;
            if(rs != EMPTY)
              break ;
          }
          if(rs == EMPTY) 
            for(variableSet::const_iterator vi = cluster.begin();
                vi!=cluster.end();++vi) {
              rs = extract_rules(dg[vi->ident()] - qualified_rules) ;
              vget = *vi ;
              if(rs != EMPTY)
                break ;
            }
        }

        // Pick the first rule in the set of rules that generates it
        // and query the rule for the type of variable *vi.  Set the variable
        // type appropriately in the fact database.
        if(rs != EMPTY) {
	  rule pick = *rs.begin() ;
	  // scan for rule that can type variable
	  for(auto ri = rs.begin();ri!=rs.end();++ri)
	    if(ri->get_info().rule_impl->get_store(vget)!= 0) {
	      pick = *ri ;
	      break ;
	    }
	  // Now get the type from the picked rule
	  storeRepP st ;
	  st = pick.get_info().rule_impl->get_store(vget) ;
	  if(st == 0) {
	    cerr << "rule " << pick << " unable to provide type for " << vget
		 << endl ;
	    abort() ;
	  }
	  st = st->new_store(EMPTY) ;
	  scheds.set_variable_type(v,st, facts) ;
	} else {

	  cerr << "Warning, variable " << v << " not typed!" << endl;
          
	}
	
      }
      variableSet syns = cluster_map[v] ;
      syns -= v ;
      for(variableSet::const_iterator vi=syns.begin();vi!=syns.end();++vi) {
        scheds.synonym_variable(v,*vi,facts) ;
      }
    }

    // Check to make sure there are no type conflicts, loop over all
    // rules that are not internally generated.
    bool type_error = false ;
    ruleSet::const_iterator ri ;
    ruleSet rs = all_rules ;
    rs -= qualified_rules ;
    variableSet typed_vars = facts.get_typed_variables() ;
    for(ri=rs.begin();ri!=rs.end();++ri) {
      variableSet varcheck ;
      const rule_impl::info &finfo = ri->get_info().desc ;
      // get the keyspace tag the rule has
      bool skip_type_check = false ;
#ifdef DYNAMICSCHEDULING
      string ks_tag = ri->get_info().rule_impl->get_keyspace_tag() ;
      if(ks_tag != "main") {
        map<string,KeySpaceP>::const_iterator mi ;
        mi = facts.keyspace.find(ks_tag) ;
        if(mi == facts.keyspace.end()) {
          cerr << "Error: rule: " << *ri << " in Non-exist keyspace: "
               << ks_tag << endl ;
          type_error = true ;
          continue ;
        }
        if(mi->second->get_dynamism() == DYNAMIC)
          skip_type_check = true ;
      }
#endif
      // Collect all variables for which are actually read or written in the class
      set<vmap_info>::const_iterator i ;
      for(i=finfo.sources.begin();i!=finfo.sources.end();++i) {
        for(size_t j=0;j<i->mapping.size();++j) {
          varcheck += i->mapping[j] ;
	}
        varcheck += i->var ;
      }
      for(i=finfo.targets.begin();i!=finfo.targets.end();++i) {
        for(size_t j=0;j<i->mapping.size();++j) {
          varcheck += i->mapping[j] ;
	}
        varcheck += i->var ;
	for(size_t k=0;k<i->assign.size();++k) {
          varcheck -= i->assign[k].first ;
          varcheck += i->assign[k].second ;
	}
      }
      variableSet::const_iterator vi ;
      for(vi = varcheck.begin();vi!=varcheck.end();++vi) {
        if(ri->get_rule_implP()==0) {
          cerr << "rule " << *ri << " has no rule_implP" <<endl ;
          type_error = true ;
          continue ;
        } else if(ri->get_rule_implP()->get_store(*vi) == 0) {
          cerr << "variable " << *vi << " not found in rule "
               << *ri << endl ;
          type_error = true ;
          continue ;
        }

        storeRepP rule_type = ri->get_rule_implP()->get_store(*vi)->getRep() ;
        if(typed_vars.inSet(*vi)) {
          if(!skip_type_check) {
            storeRepP fact_type = facts.get_variable(*vi)->getRep() ;
	    auto &rule_type_dr = *rule_type ;
	    auto &fact_type_dr = *fact_type ;
            if(typeid(rule_type_dr) != typeid(fact_type_dr)) {
              cerr << "variable type mismatch for variable " << *vi << " in rule "
                   << *ri << endl ;
              cerr << "fact database has type " << typeid(fact_type_dr).name() << endl ;
              cerr << "rule has type " << typeid(rule_type_dr).name() << endl ;
              type_error = true ;
            }
          }
        } else {
          cerr << "Untyped Variable " << *vi << endl ;
          scheds.set_variable_type(*vi,rule_type, facts) ;
          typed_vars += *vi ;
        }
      }
    }
    if(type_error)
      Loci::Abort() ;
  }


}
