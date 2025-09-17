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
#ifndef DEPEND_GRAPH_H
#define DEPEND_GRAPH_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/digraph.h>
#include <Tools/cptr.h>
#include <Tools/variable.h>
#include <rule.h>

namespace Loci {

  inline ruleSet extract_rules(digraph::vertexSet ns) {
    return ruleSet(ns & interval(UNIVERSE_MIN,-1)) ;
  }

  inline variableSet extract_vars(digraph::vertexSet ns) {
    return variableSet(ns & interval(0,UNIVERSE_MAX)) ;
  }

  // experimental version
  class dependency_graph2 {
    digraph gr ;
  public:
    dependency_graph2(const rule_db& rdb,
                      const variableSet& given,
                      const variableSet& target) ;
    
    digraph get_graph() const { return gr; } 
  };

  digraph partition_iteration(digraph gr) ;

  void clean_graph(digraph &gr,
                   const variableSet& given,
                   const variableSet& target) ;
}

#endif
