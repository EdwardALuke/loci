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
#include "sched_tools.h"
#include "loci_globs.h"

using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;
using std::ostringstream ;

namespace Loci {
// Create a sequential ordering of rules based on the concurrent dag schedule
  void extract_rule_sequence(vector<rule> &rule_seq,
                             const vector<digraph::vertexSet> &v) {
    vector<digraph::vertexSet>::const_iterator i ;
    for(i=v.begin();i!=v.end();++i) {
      ruleSet rules = extract_rules(*i) ;
      ruleSet::const_iterator ri ;
      for(ri=rules.begin();ri!=rules.end();++ri)
        rule_seq.push_back(*ri) ;
    }
  }
  
  rule make_rename_rule(variable new_name, variable old_name) {
    ostringstream oss ;
    oss << "source(" << old_name << "),target(" << new_name
        << "),qualifier(rename)"  ;
    return rule(oss.str()) ;
  }

  entitySet vmap_source_exist(const vmap_info &vmi, fact_db &facts, sched_db &scheds) {
    variableSet::const_iterator vi ;
    entitySet sources = ~EMPTY ;
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi) {
      sources &= scheds.variable_existence(*vi) ;
    }

    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      entitySet working = ~EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!scheds.is_a_Map(*vi)) ;
	working &= scheds.preimage(*vi,sources).first ;
      }
      sources = working ;
    }
    return sources ;
  }


  entitySet vmap_target_exist(const vmap_info &vmi, fact_db &facts,
                              entitySet compute, sched_db &scheds) {
    vector<variableSet>::const_iterator mi ;
    for(mi=vmi.mapping.begin();mi!=vmi.mapping.end();++mi) {
      if(mi->size() == 1) {
        variable v = *(mi->begin()) ;
        FATAL(!scheds.is_a_Map(v)) ;
        compute = scheds.image(v,compute) ;
      } else {
        variableSet::const_iterator vi ;
        entitySet images ;
        for(vi=mi->begin();vi!=mi->end();++vi) {
          variable v = *vi ;
          FATAL(!scheds.is_a_Map(v)) ;
          images |= scheds.image(v,compute) ;
        }
        compute = images ;
      }
    }
    return compute ;
  }

}
