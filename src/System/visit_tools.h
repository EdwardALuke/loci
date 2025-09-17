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
#ifndef VISIT_TOOLS_H
#define VISIT_TOOLS_H

#include "visitor.h"
#include "comp_tools.h"

#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include <list>

namespace Loci {

  // some inlined small functions
  
  inline rule create_rule(variable sv, variable tv, std::string qualifier) {
    std::ostringstream oss ;
    oss << "source(" << sv << ')' ;
    oss << ",target(" << tv << ')' ;
    oss << ",qualifier(" << qualifier << ')' ;
    std::string sig = oss.str() ;
    rule r(sig) ;
    return r ;
  }
  
  inline rule create_rule(variableSet source, variableSet target,
                          std::string qualifier) {
    std::ostringstream oss ;
    oss << "source(" << source << ')' ;
    oss << ",target(" << target << ')' ;
    oss << ",qualifier(" << qualifier << ")" ;
    std::string sig = oss.str() ;
    rule r(sig) ;
    return r ;
  }

  inline bool is_super_node(int rid)
  {
    if(rid >= 0) // a variable
      return false ;
    rule r(rid) ;
    std::string rqualifier = r.get_info().qualifier() ;
    
    return (rqualifier.substr(0,2) == "SN") ;
  }
  
  inline bool is_super_node(const ruleSet::const_iterator& ruleIter) {
    std::string rqualifier = ruleIter->get_info().qualifier() ;

    return (rqualifier.substr(0,2) == "SN") ;
  }

  inline bool is_super_node(const rule& r) {
    std::string rqualifier = r.get_info().qualifier() ;

    return ( (!rqualifier.empty()) && rqualifier.substr(0,2) == "SN") ;
  }

  // this one returns if a super rule node is a chomping rule
  // by comparing the first 5 characters to the string "CHOMP"
  inline bool
  is_chomp_node(const rule& r) {
    std::string q = r.get_info().qualifier() ;
    return ( (!q.empty()) && q.substr(0,5) == "CHOMP") ;
  }
  inline bool
  is_chomp_node(const ruleSet::const_iterator& ri) {
    return is_chomp_node(*ri) ;
  }

  inline int get_supernode_num(const rule& r) {
    std::string rqualifier = r.get_info().qualifier() ;
    std::string head = rqualifier.substr(0,2) ;
    if(head != "SN") {
      std::cerr << "get_supernode_num error! pass in rule is not a super node!"
           << std::endl ;
      exit(-1) ;
    }
      
    std::string number = rqualifier.substr(2,rqualifier.size()-2) ;
    std::stringstream ss ;
    ss << number ;
    int ret ;
    ss >> ret ;
      
    return ret ;
  }

  template<typename T>
  inline bool inSet(const std::set<T>& s, const T& elem) {
    typename std::set<T>::const_iterator si ;
    si = s.find(elem) ;
    if(si == s.end())
      return false ;
    return true ;
  }

  // is there a path from source vertex to target vertex in gr?
  inline bool has_path(const digraph& gr, int source, int target) {
    bool path = false ;
    digraph::vertexSet working = gr[source] ;
    while(working != EMPTY) {
      if(working.inSet(target)) {
        path = true ;
        break ;
      }
      digraph::vertexSet children ;
      for(digraph::vertexSet::const_iterator vi=working.begin();
          vi!=working.end();++vi) {
        children += gr[*vi] ;
      }
      working = children ;
    }

    return path ;
  }

  // is there path(s) from sources vertexSet to
  // targets vertexSet in the given gr?
  inline bool has_path(const digraph& gr,
                       const digraph::vertexSet& sources,
                       const digraph::vertexSet& targets) {
    // if sources and targets are overlapped, then path(s) exist
    if( (sources & targets) != EMPTY)
      return true ;
    if(targets == EMPTY)
      return false ;
    digraph::vertexSet working ;
    for(digraph::vertexSet::const_iterator vi=sources.begin();
        vi!=sources.end();++vi)
      working += gr[*vi] ;
    
    digraph::vertexSet visited ;
    while(working != EMPTY) {
      if( (working & targets) != EMPTY)
        return true ;
      visited += working ;
      digraph::vertexSet next ;
      for(digraph::vertexSet::const_iterator vi=working.begin();
          vi!=working.end();++vi)
        next += gr[*vi] ;
      next -= visited ;
      working = next ;
    }

    return false ;
  }
                       

  // given a variableSet, convert them to vertexSet
  inline digraph::vertexSet get_vertexSet(const variableSet& vars) {
    digraph::vertexSet ret ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi)
      ret += vi->ident() ;

    return ret ;
  }

  // given a ruleSet, convert them to vertexSet
  inline digraph::vertexSet get_vertexSet(const ruleSet& rules) {
    digraph::vertexSet ret ;
    for(ruleSet::const_iterator ri=rules.begin();ri!=rules.end();++ri)
      ret += ri->ident() ;

    return ret ;
  }

  // get the init source variableSet from a target variable and
  // the recurrence target to source table, or
  // get the final target variableSet from a source variable and
  // the recurrence source to target table
  inline
  variableSet get_leaf_recur_vars(const std::map<variable,variableSet>& t,
                                  const variable& v) {
    std::map<variable,variableSet>::const_iterator found ;
    variableSet ret ;
    std::list<variable> working ;

    // first check if v is in the table
    found = t.find(v) ;
    // if v is not recurrence variable, return empty set
    if(found == t.end())
      return variableSet(EMPTY) ;

    // do a depth first search
    working.push_front(v) ;
    while(!working.empty()) {
      variable cur = working.front() ;
      working.pop_front() ;
        
      variableSet tmp ;
      found = t.find(cur) ;
        
      if(found != t.end())
        tmp = found->second ;
      else
        ret += cur ;
        
      for(variableSet::const_iterator vi=tmp.begin();
          vi!=tmp.end();++vi)
        working.push_front(*vi) ;
    }

    return ret ;
  }

  // get all the source variables from a target variable and
  // the recurrence target to source table, or
  // get all the target variables from a source variable and
  // the recurrence source to target table
  inline
  variableSet get_all_recur_vars(const std::map<variable,variableSet>& t,
                                 const variable& v) {
    std::map<variable,variableSet>::const_iterator found ;
    variableSet ret, working ;
      
    // first check if v is in the table
    found = t.find(v) ;
    // if v is not recurrence variable, return empty set
    if(found == t.end())
      return variableSet(EMPTY) ;

    // do a breadth first search, and add up all the results
    working += v ;
    while(working != EMPTY) {
      variableSet tmp ;
      for(variableSet::const_iterator vi=working.begin();
          vi!=working.end();++vi) {
        found = t.find(*vi) ;
        if(found != t.end())
          tmp += found->second ;
      }
      ret += tmp ;
      working = tmp ;
    }

    return ret ;
  }

  // given a variable, if it is a recurrence variable, then get
  // all other leaf target variables that all its source can reach
  inline
  variableSet get_all_leaf_target(const std::map<variable,variableSet>& t2s,
                                  const std::map<variable,variableSet>& s2t,
                                  const variable& v) {
    std::map<variable,variableSet>::const_iterator found ;
    variableSet ret, working ;

    working += v ;
    while(working != EMPTY) {
      variableSet tmp ;
      for(variableSet::const_iterator vi=working.begin();
          vi!=working.end();++vi) {
        found = t2s.find(*vi) ;
        if(found != t2s.end()) {
          variableSet cur = found->second ;
          tmp += cur ;
          for(variableSet::const_iterator vi2=cur.begin();
              vi2!=cur.end();++vi2)
            ret += get_leaf_recur_vars(s2t,*vi2) ;
        }
      }
      working = tmp ;
    }

    ret -= v ;
    return ret ;
  }

  inline bool is_dg_empty(const digraph& gr) {
    digraph::vertexSet allv = gr.get_all_vertices() ;
    return (allv == EMPTY) ;
  }

  // return if a rule is a generalize, promote or priority rule
  inline bool is_recur_rule(const ruleSet::const_iterator& ruleIter) {
    if(ruleIter->type() == rule::INTERNAL) {
      if( (ruleIter->get_info().qualifier() == "promote") ||
          (ruleIter->get_info().qualifier() == "generalize") ||
          (ruleIter->get_info().qualifier() == "priority")
          )
        return true ;
    }
    return false ;
  }
  
  inline bool is_internal_rule(const ruleSet::const_iterator& ruleIter) {
    //return (ruleIter->type() == rule::INTERNAL) ;
    return (ruleIter->get_info().rule_class == rule::INTERNAL) ;
  }

  inline bool is_internal_rule(const rule& r) {
    //return (r.type() == rule::INTERNAL) ;
    return (r.get_info().rule_class == rule::INTERNAL) ;
  }
  
  inline bool is_internal_rule(int i) {
    if(i>=0) // is a variable
      return false ;
    rule r(i) ;
    return (r.type() == rule::INTERNAL) ;
  }

  inline bool is_virtual_rule(const rule& r) {
    if(r.type() == rule::INTERNAL) {
      return (
              is_super_node(r) ||
              (r.get_info().qualifier() == "generalize") ||
              (r.get_info().qualifier() == "promote") ||
              (r.get_info().qualifier() == "priority")
              ) ;
    }else
      return false ;
  }

  inline bool time_before(const variable& v1, const variable& v2) {
    return v1.time().before(v2.time()) ;
  }
  
  inline bool time_equal(const variable& v1, const variable& v2) {
    return (!time_before(v1,v2) && !time_before(v2,v1)) ;
  }
  
  inline bool time_after(variable v1, variable v2) {
    return (!time_before(v1,v2) && !time_equal(v1,v2)) ;
  }

  // given a digraph, this function returns a vector
  // of digraph that are separate subgraphs in the
  // given digraph. That means each subgraph is just
  // dangling there in the digraph and does not have
  // connections to each other
  std::vector<digraph> get_islands(const digraph& gr) ;

} // end of namespace Loci

#endif
