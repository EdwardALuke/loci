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
#include <fact_db.h>


#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <list>
using std::list ;

#include <istream>
#include <ostream>
#include <iostream>

using std::ostream ;
using std::istream ;
using std::cerr ;
using std::endl ;

#include <fstream>

using std::ofstream ;

namespace Loci {

  namespace {
    
    void create_digraph_dot_file(const digraph &dg, const char* fname)
    {
      digraph dgt = dg.transpose() ;
      digraph::vertexSet allvertices = dg.get_all_vertices() ;
      digraph::vertexSet::const_iterator ri ;
      ofstream outf(fname) ;
      
      outf<<"digraph G {\n" ;
      //outf<<"size = \"8.5,11\";\n" ;

      // variable's ident() >= 0 rule's ident() < 0
      for(ri=allvertices.begin();ri!=allvertices.end();++ri) {
        digraph::vertexSet outvertices = dg[*ri] ;
        digraph::vertexSet incomevertices = dgt[*ri] ;

        if(*ri < 0) {
          rule r(*ri) ;
          if(r.type() == rule::INTERNAL) {
            outf<<"\""<<pretty_sig(r)<<"\""
                <<"[shape=doubleoctagon,style=filled,color=gold];\n" ;
          }
          else {
            outf<<"\""<<pretty_sig(r)<<"\""<<"[shape=box,style=filled,color=gold];\n" ;
          }
        }
        else {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=green];\n" ;
        }
        if(incomevertices == EMPTY) {
          if(*ri >= 0) {
            variable v(*ri) ;
            outf<<"\""<<v<<"\""<<"[style=filled,color=red];\n" ;
          }else {
            rule r(*ri) ;
            outf<<"\""<<pretty_sig(r)<<"\""<<"[style=filled,color=red];\n" ;
          }
        }
        if(outvertices == EMPTY) {
          if(*ri >= 0) {
            variable v(*ri) ;
            outf<<"\""<<v<<"\""<<"[style=filled,color=blue];\n" ;
          }else {
            rule r(*ri) ;
            outf<<"\""<<pretty_sig(r)<<"\""<<"[style=filled,color=blue];\n" ;
          }
          continue ;
        }
        
        if(*ri < 0) {
          rule r(*ri) ;
          digraph::vertexSet::const_iterator ii ;
          for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
            if(*ii < 0) {
              rule r2(*ii) ;
              outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""<<";\n" ;
            }else {
              variable v(*ii) ;
              outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<v<<"\""<<";\n" ;
            }
          }
        }else {
          variable v(*ri) ;
          digraph::vertexSet::const_iterator ii ;
          for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
            if(*ii < 0) {
              rule r(*ii) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<pretty_sig(r)<<"\""<<";\n" ;
            }else {
              variable v2(*ii) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<v2<<"\""<<";\n" ;
            }
          }        
        }
      }
      outf<<"}\n" ;
      outf.close() ;
    }

  } // end of unnamed namespace

  /////////////////////////////////////////////////////////////
  // graphVisualizeVisitor
  /////////////////////////////////////////////////////////////
  void graphVisualizeVisitor::visit(loop_compiler& lc) {
    if(Loci::get_exec_rank()==0) {
      cerr << "the looping part graph (collapse and advance), id: "
           << lc.cid << endl ;
      create_digraph_dot_file(lc.loop_gr,"visual_graph.dot") ;
      int err = system("dotty visual_graph.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << "dotty ..." << "'" 
	     << endl ;

      err = system("rm -fr visual_graph.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << "rm -fr ..." << "'" 
	     << endl ;

      cerr << "\tthe collapse part graph" << endl ;
      create_digraph_dot_file(lc.collapse_gr,"visual_graph.dot") ;
      err = system("dotty visual_graph.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << "dotty ..." << "'" 
	     << endl ;
      err = system("rm -fr visual_graph.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << "rm ..." << "'" 
	     << endl ;
      
      cerr << "\tthe advance part graph" << endl ;
      create_digraph_dot_file(lc.advance_gr,"visual_graph.dot") ;
      err = system("dotty visual_graph.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << "dotty ..." << "'" 
	     << endl ;
      err = system("rm -fr visual_graph.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << "rm ..." << "'" 
	     << endl ;
    }
  }

  void graphVisualizeVisitor::visit(dag_compiler& dc) {
    if(Loci::get_exec_rank()==0) {
      cerr << "the dag graph, id: " << dc.cid << endl ;
      create_digraph_dot_file(dc.dag_gr,"visual_graph.dot") ;
      int err = system("dotty visual_graph.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << "dotty ..." << "'" 
	     << endl ;
      err = system("rm -fr visual_graph.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << "rm ..." << "'" 
	     << endl ;
      
    }
  }

  void graphVisualizeVisitor::visit(conditional_compiler& cc) {
    if(Loci::get_exec_rank()==0) {
      cerr << "the conditional graph, id: " << cc.cid << endl ;
      create_digraph_dot_file(cc.cond_gr,"visual_graph.dot") ;
      int err = system("dotty visual_graph.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << "dotty ..." << "'" 
	     << endl ;
      err = system("rm -fr visual_graph.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << "rm ..." << "'" 
	     << endl ;

    }
  }

  ///////////////////////////////////////////////////////////////
  // dagCheckVisitor
  ///////////////////////////////////////////////////////////////
  ostream& dagCheckVisitor::visualize(ostream& s) const {
    if(Loci::get_exec_rank() == 0) {
      s << "visualizing detected cycle..." << endl ;
      create_digraph_dot_file(cycle,"cycle_in_dag.dot") ;
      int err = system("dotty cycle_in_dag.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" 
	     << "dotty cycle_in_dag.dot"
	     << "'" << endl ;
      err = system("rm -fr cycle_in_dag.dot") ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" 
	     << "rm -fr cycle_in_dag.dot"
	     << "'" << endl ;
    }
    return s ;
  }

  ///////////////////////////////////////////////////////////////
  // chompRuleVisitor
  ///////////////////////////////////////////////////////////////
  ostream& chompRuleVisitor::visualize(ostream& s) const {
    if(Loci::get_exec_rank() == 0) {
      if(all_chains.empty()) {
        s << "NO chomping chains found!" << endl ;
        return s ;
      }
      map<int,list<chomp_chain> >::const_iterator mi ;
      list<chomp_chain>::const_iterator li ;
      int total_chomp_vars = 0 ;
      for(mi=all_chains.begin();mi!=all_chains.end();++mi) {
        list<chomp_chain> cclist = mi->second ;
        s << "There are " << cclist.size()
          << " chomping rule chain(s) in super node SN"
          << mi->first << ": " << endl ;
        int i ;
        for(li=cclist.begin(),i=1;li!=cclist.end();++li,++i) {
          int total = li->second.size() ;
          total_chomp_vars += total ;
          s << "Visualizing the NO" << i << " chomping chain..." << endl ;
          s << "variables that can be chomped in this chain: "
            << li->second << " total: " << total << endl ;
          
          s << "source variables to the chomp node: " << extract_vars(li->first.get_source_vertices() - li->first.get_target_vertices()) << endl ;
          s << "target variables of the chomp node: " << extract_vars(li->first.get_target_vertices() - li->first.get_source_vertices()) << endl ;
          
          create_digraph_dot_file(li->first,"chomping_rules.dot") ;
          int err = system("dotty chomping_rules.dot") ;
	  if(err != 0)
	    cerr << "system call returned " << err << " in call '" 
		 << "dotty chomping_rules.dot"
		 << "'" << endl ;
          err = system("rm -fr chomping_rules.dot") ;
	  if(err != 0)
	    cerr << "system call returned " << err << " in call '" 
		 << "rm -fr chomping_rules.dot"
		 << "'" << endl ;
        }
      }
      s << "total variables can be chomped in the program: "
        << total_chomp_vars << endl ;
    }
    return s ;
  }

  ostream& chompRuleVisitor::summary(ostream& s) const {
    if(Loci::get_exec_rank() == 0) {
      s << "--------------begin chomping summary--------------" << endl ;
      s << "Theoretically there are: " << good_vars.size()
        << " variables that can be chomped in"
        << " the program (the upper bound)." << endl ;
      if(all_chains.empty()) {
        s << "NO chomping chains found!" << endl ;
        s << "---------------end chomping summary---------------" << endl ;
        return s ;
      }
      map<int,list<chomp_chain> >::const_iterator mi ;
      list<chomp_chain>::const_iterator li ;
      int total_chomp_vars = 0 ;
      for(mi=all_chains.begin();mi!=all_chains.end();++mi) {
        list<chomp_chain> cclist = mi->second ;
        s << "There are: " << cclist.size()
          << " chomping rule chain(s) in super node SN"
          << mi->first << ": " << endl ;
        int i ;
        for(li=cclist.begin(),i=1;li!=cclist.end();++li,++i) {
          int total = li->second.size() ;
          total_chomp_vars += total ;
          s << "NO" << i << ": " ;
          s << "chomped variables: "
            << li->second << " (total: " << total << ")" << endl ;
        }
      }
      s << endl << "Total variables chomped in the program: "
        << total_chomp_vars << endl ;
      s << "---------------end chomping summary---------------" << endl ;
    }

    return s ;
  }

  
} // end of namespace Loci
