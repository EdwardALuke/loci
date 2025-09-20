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
#include "sched_mlg.h"
#include "distribute.h"
using std::ostringstream ;

//#define VERBOSE
namespace Loci {

  extern bool verbose ;

  void multiLevelGraph::subGraph::create(digraph &gin,
                                         digraph::vertexSet &grvtx) {
    digraph gint = gin.transpose() ;

    digraph::vertexSet sourcev, targetv ;
    digraph::vertexSet::const_iterator vi ;

    for(vi=grvtx.begin();vi!=grvtx.end();++vi) {
      sourcev += gint[*vi] ;
      targetv += gin[*vi] ;
    }

    sourcev -= grvtx ;
    targetv -= grvtx ;

    graph_v = grvtx ;
    incoming_v = sourcev ;
    outgoing_v = targetv ;
    digraph::vertexSet all = sourcev+targetv+grvtx ;
    
    gr = gin.subgraph(all) ;
  }
  
  multiLevelGraph::multiLevelGraph(digraph gr, digraph::vertexSet grvtx):
    super_node_number(0) {
    subGraph g(gr,grvtx) ;

    variableSet sv = extract_vars(g.incoming_v) ;
    variableSet tv = extract_vars(g.outgoing_v) ;

    rule r = make_super_rule(sv,tv) ;
    if(verbose) {
      debugout << "----------------------------------------------------------------------------" << endl ;
      debugout << "multiLevelGraph: creating super node: " << r << endl ;
      ruleSet rset = extract_rules(g.graph_v) ;
      debugout << rset << endl ;
      debugout << "----------------------------------------------------------------------------" << endl ;
    }

    toplevel = r.ident() ;

    insert(toplevel,g) ;
  }

  
  multiLevelGraph::subGraph *multiLevelGraph::find(int id) {
    std::map<int,subGraph>::iterator mi ;
    mi = subgrmap.find(id) ;
    if(mi == subgrmap.end())
      return 0 ;
    else
      return &mi->second ;
  }

  rule multiLevelGraph::make_super_rule(variableSet sources,
                                        variableSet targets,
                                        variable cond) {
    if(targets == EMPTY) {
      cerr << "make_super_rule called with targets == EMPTY" << endl ;
      cerr << "sources = " << sources << endl ;
      cerr << "cond = " << cond << endl ;
    }
    
    FATAL(targets == EMPTY) ;
    ostringstream oss ;
    oss << "source("<<sources << "),target(" << targets << ")," ;
    if(cond != variable()) 
      oss<< "conditional(" << cond << ")," ;
    oss << "qualifier(SN" << super_node_number++ << ")" ;
   
    return rule(oss.str()) ;
  }
  
  int multiLevelGraph::mksnode(int gr_id, digraph::vertexSet grvtx,
                               variable cond_var) {
    subGraph *g = find(gr_id) ;
    fatal(g == 0) ;

    subGraph sg(g->gr,grvtx) ;
    variableSet sv = extract_vars(sg.incoming_v) ;
    variableSet tv = extract_vars(sg.outgoing_v) ;
    rule r = make_super_rule(sv,tv,cond_var) ;
    
    if(verbose) {
      debugout << "----------------------------------------------------------------------------" << endl ;
      debugout << "mksnode: creating super node: " << r << endl ;
      ruleSet rset = extract_rules(sg.graph_v) ;
      debugout << rset << endl ;
      debugout << "----------------------------------------------------------------------------" << endl ;
    }
    insert(r.ident(),sg) ;
#ifdef VERBOSE
    if(extract_rules(sg.incoming_v).size() !=0) {
      debugout << "extract_rules(sg.incoming_v) = " << extract_rules(sg.incoming_v)
           << endl ;
    }
    if(extract_rules(sg.outgoing_v).size() != 0) {
      debugout << "extract_rules(sg.outgoing_v) = " << extract_rules(sg.outgoing_v)
           << endl ;
    }
#endif
    g->gr.remove_vertices(grvtx) ;
    g->gr.add_edges(sg.incoming_v,r.ident()) ;
    g->gr.add_edges(r.ident(),sg.outgoing_v) ;
    g->graph_v -= grvtx ;
    g->graph_v += r.ident() ;
    return r.ident() ;
  }
  
}
