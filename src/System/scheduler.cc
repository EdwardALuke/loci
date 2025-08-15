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
#include "loci_globs.h"
#include "sched_tools.h"
#include "dist_tools.h"
#include "param_rule.h"
#include "thread.h"
#include <Tools/except.h>
#include <constraint.h>
#include <new>
using std::bad_alloc ;
using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

using std::ostringstream ;
using std::string ;
using std::endl ;
using std::cout ;
using std::ios ;
using std::ofstream ;
using std::istream ;
using std::ostream ;

///////////////////////////////////
#include <sstream>
#include <algorithm>
#include <sys/time.h> // for gettimeofday function
///////////////////////////////////

#define PROFILE_CODE

#ifdef USE_PAPI
#include "papi.h"
#define N 64*64
#define NCOUNTS 6
#endif

namespace Loci {
  //#define DYNAMIC_TIMING
#ifdef DYNAMIC_TIMING
  extern timeAccumulator ta_expand ;
  extern timeAccumulator ta_expand2 ;
  extern timeAccumulator ta_expand_start ;
  extern timeAccumulator ta_expand_cache ;
  extern timeAccumulator ta_expand_collect_img ;
  extern timeAccumulator ta_expand_missing ;
  extern timeAccumulator ta_context ;
  extern timeAccumulator ta_context_nonepd ;
  extern timeAccumulator ta_context_nonepd_domt ;
  extern timeAccumulator ta_context_nonepd_ints ;
  extern timeAccumulator ta_context_epdend ;
  extern timeAccumulator ta_context_epdmid ;
  extern timeAccumulator ta_context_epdsta ;
  extern timeAccumulator ta_output_oh ;
  extern timeAccumulator ta_compute ;
  extern timeAccumulator ta_erase ;
  extern timeAccumulator ta_invalidate ;
  extern timeAccumulator ta_keyremoval ;
  extern timeAccumulator ta_insertion ;
  extern timeAccumulator ta_keyinsert ;
  extern timeAccumulator ta_key_dist ;
  extern timeAccumulator ta_dist_renumber ;
  extern timeAccumulator ta_dist ;
  extern int             ta_dist_number ;
  extern timeAccumulator ta_push ;
  extern timeAccumulator ta_expand_comm ;
  extern timeAccumulator ta_pw_push ;
  extern timeAccumulator ta_param_push ;
  extern timeAccumulator ta_param_pack ;
  extern timeAccumulator ta_param_unpack ;
  extern timeAccumulator ta_param_reduce ;
  extern int             ta_param_reduce_num ;
  extern timeAccumulator ta_record_erase ;
  extern timeAccumulator ta_dctrl ;
  extern int             ta_drule_executes ;
  extern int             ta_dctrl_executes ;

  stopWatch sw_dist_keys ;
  stopWatch sw_dist_all2all ;
  stopWatch sw_dist_pack ;
  stopWatch sw_dist_unpack ;
  timeAccumulator ta_dist_keys ;
  timeAccumulator ta_dist_all2all ;
  timeAccumulator ta_dist_pack ;
  timeAccumulator ta_dist_unpack ;
#endif

  //#define PARTICLE_PERF
#ifdef PARTICLE_PERF
  stopWatch sw_particle ;
  double pwalltime = 0 ;
#endif

  double LociAppPeakMemory = 0 ;
  double LociAppAllocRequestBeanCounting = 0 ;
  double LociAppFreeRequestBeanCounting = 0 ;
  double LociAppPeakMemoryBeanCounting = 0 ;
  double LociAppLargestAlloc = 0 ;
  variable LociAppLargestAllocVar("EMPTY") ;
  double LociAppLargestFree = 0 ;
  variable LociAppLargestFreeVar("EMPTY") ;
  double LociAppPMTemp = 0 ;

  double LociInputVarsSize = 0 ;

  ////////////////////////////
  // global flags
  extern bool profile_memory_usage ;
  extern bool collect_memory_info ;
  extern bool show_graphs ;
  extern void deco_depend_gr(digraph& gr,const variableSet& given) ;
  extern bool threading_pointwise;
  extern bool threading_global_reduction;
  extern bool threading_local_reduction;
  extern bool threading_chomping;
  extern bool threading_recursion;
  extern int num_threads;  
  extern int num_thread_blocks;
  // 
  int num_threaded_pointwise = 0;
  int num_total_pointwise = 0;

  int num_threaded_global_reduction = 0;
  int num_total_global_reduction = 0;

  int num_threaded_local_reduction = 0;
  int num_total_local_reduction = 0;

  int num_threaded_chomping = 0;
  int num_total_chomping = 0;

  int num_threaded_recursion = 0;
  int num_total_recursion = 0;
  ////////////////////////////

  namespace {
    // pretty printing of a rule's signature
    // i.e. remove namespace info, if any
    inline std::string pretty_sig(const rule& r) {
      //the following line is the simplest, but it does not
      //include the qualify, such like "SN1:" in a super node
      //return r.get_info().desc.rule_identifier() ;
      //ostringstream oss ;
      //oss << "R" << r.ident() << endl ;
      //return oss.str() ;
      std::string name = r.get_info().name() ;
      if(r.type() == rule::INTERNAL) {
        return name ;
      }else {
        std::string::iterator pos ;
        pos = std::find(name.begin(),name.end(),'#') ;
        return std::string( (pos==name.end()?name.begin():pos+1),name.end()) ;
      }
    }
  } // end of unnamed namespace

  ////////////////////////////////////////////////////////////////////
  // the following part are functions to visualize loci internal graphs
  // they create files to be used by "dot", "lefty" and "dotty"
  // programs from AT&T research
  void create_digraph_dot_file(const digraph &dg, const char* fname)
  {
    digraph dgt = dg.transpose() ;
    digraph::vertexSet allvertices = dg.get_all_vertices() ;
    digraph::vertexSet::const_iterator ri ;
    ofstream outf(fname) ;

    outf<<"digraph G {\n" ;
    //outf<<"size = \"8.5,11\";\n" ;

    for(ri=allvertices.begin();ri!=allvertices.end();++ri) {

      digraph::vertexSet outvertices = dg[*ri] ;
      digraph::vertexSet incomevertices = dgt[*ri] ;

      if(*ri < 0) {
        rule r(*ri) ;
        if(r.type() == rule::INTERNAL)
          outf<<"\""<<pretty_sig(r)<<"\""
              <<"[shape=doubleoctagon,style=filled,color=gold];\n" ;
        else
          outf<<"\""<<pretty_sig(r)<<"\""
              <<"[shape=box,style=filled,color=gold];\n" ;
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
            outf<<"\""<<pretty_sig(r)<<"\""<<" -> "
                <<"\""<<pretty_sig(r2)<<"\""<<";\n" ;
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

  bool is_super_rule(int rid)
  {
    if(rid >= 0) // a variable
      return false ;
    rule r(rid) ;
    std::string rqualifier = r.get_info().qualifier() ;

    return (rqualifier.substr(0,2) == "SN") ;
  }

  std::string get_assigned_cluster_name(int rid)
  {
    rule r(rid) ;
    std::string rqualifier = r.get_info().qualifier() ;
    // reads until ':'
    // e.g. SN12: .....
    std::string number = rqualifier.substr(2,rqualifier.size()-2) ;

    // get the cluster name
    std::string clustername = "cluster" ;
    clustername += number ;

    return clustername ;
  }

  // given rid is a supernode, connect it to its targets
  void writeout_super_rule(multiLevelGraph &mlg, int rid,
                           const digraph::vertexSet &targets,
                           std::ofstream &outf)
  {
    int repnode = *(mlg.find(rid)->graph_v-mlg.subgraphs).begin() ;
    // if the picked node is a rule
    if(repnode < 0) {
      rule r(repnode) ;
      digraph::vertexSet::const_iterator ii ;
      for(ii=targets.begin();ii!=targets.end();++ii) {
        if(*ii < 0) {// if the target is a rule
          if(is_super_rule(*ii)) {// if the target rule is a super node
            //then we pick a node in that super node
            int repnode2 = *(mlg.find(*ii)->graph_v-mlg.subgraphs).begin() ;
            if(repnode2 < 0) { // if the picked node is a rule
              rule r2(repnode2) ;
              outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""
                  <<pretty_sig(r2)<<"\""
                  <<"[style=dotted,color=red,ltail="
                  << get_assigned_cluster_name(rid) << ",lhead="
                  << get_assigned_cluster_name(*ii) << "]" << ";\n" ;
            }else { // the picked node is a variable
              variable v2(repnode2) ;
              outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<v2<<"\""
                  <<"[style=dotted,color=red,ltail="
                  << get_assigned_cluster_name(rid) << ",lhead="
                  << get_assigned_cluster_name(*ii) << "]" << ";\n" ;
            }
          }else { // if the target rule is a common single rule
            rule r2(*ii) ;
            outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""
                <<"[style=dotted,color=red,ltail="
                << get_assigned_cluster_name(rid) << "]" << ";\n" ;
          }
        }else { // if the target is variable
          variable v2(*ii) ;
          outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<v2<<"\""
              <<"[style=dotted,color=red,ltail="
              << get_assigned_cluster_name(rid) << "]" << ";\n" ;
        }
      }
    }else { // the picked node is a variable
      variable v(repnode) ;
      digraph::vertexSet::const_iterator ii ;
      for(ii=targets.begin();ii!=targets.end();++ii) {
        if(*ii < 0) {// if the target is a rule
          if(is_super_rule(*ii)) {// if the target rule is a super node
            //then we pick a node in that super node
            int repnode2 = *(mlg.find(*ii)->graph_v-mlg.subgraphs).begin() ;
            if(repnode2 < 0) { // if the picked node is a rule
              rule r2(repnode2) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""
                  <<"[style=dotted,color=red,ltail="
                  << get_assigned_cluster_name(rid) << ",lhead="
                  << get_assigned_cluster_name(*ii) << "]" << ";\n" ;
            }else { // the picked node is a variable
              variable v2(repnode2) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<v2<<"\""
                  <<"[style=dotted,color=red,ltail="
                  << get_assigned_cluster_name(rid) << ",lhead="
                  << get_assigned_cluster_name(*ii) << "]" << ";\n" ;
            }
          }else { // if the target rule is a common single rule
            rule r2(*ii) ;
            outf<<"\""<<v<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""
                <<"[style=dotted,color=red,ltail="
                << get_assigned_cluster_name(rid) << "]" << ";\n" ;
          }
        }else { // if the target is variable
          variable v2(*ii) ;
          outf<<"\""<<v<<"\""<<" -> "<<"\""<<v2<<"\""
              <<"[style=dotted,color=red,ltail="
              << get_assigned_cluster_name(rid) << "]" << ";\n" ;
        }
      }
    }
  }

  // given sid is a NON-super node, tid IS a super node
  // connect sid to tid
  void writeout_super_rule2(multiLevelGraph &mlg, int sid, int tid,
                            std::ofstream &outf)
  {
    if(sid < 0) { // sid is a rule
      rule srule(sid) ;

      // pick a node in tid
      int repnode = *(mlg.find(tid)->graph_v-mlg.subgraphs).begin() ;
      if(repnode < 0) { // if the picked node is a rule
        rule r2(repnode) ;
        outf<<"\""<<srule<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""
            <<"[style=dotted,color=red,lhead="
            << get_assigned_cluster_name(tid) << "]" << ";\n" ;
      }else { // the picked node is a variable
        variable v2(repnode) ;
        outf<<"\""<<srule<<"\""<<" -> "<<"\""<<v2<<"\""
            <<"[style=dotted,color=red,lhead="
            << get_assigned_cluster_name(tid) << "]" << ";\n" ;
      }
    }else { // sid is a variable
      variable sv(sid) ;

      // pick a node in tid
      int repnode = *(mlg.find(tid)->graph_v-mlg.subgraphs).begin() ;
      if(repnode < 0) { // if the picked node is a rule
        rule r2(repnode) ;
        outf<<"\""<<sv<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""
            <<"[style=dotted,color=red,lhead="
            << get_assigned_cluster_name(tid) << "]" << ";\n" ;
      }else { // the picked node is a variable
        variable v2(repnode) ;
        outf<<"\""<<sv<<"\""<<" -> "<<"\""<<v2<<"\""
            <<"[style=dotted,color=red,lhead="
            << get_assigned_cluster_name(tid) << "]" << ";\n" ;
      }
    }
  }

  // recursively layout from the toplevel to the bottom
  // this function just lists which nodes belongs to which cluster(supernode)
  // and does not do connections.
  void layout_super_node(multiLevelGraph &mlg, int sid, std::ofstream &outf)
  {
    static int baseclustercounter = 0 ;
    // colors that fill the clusters
    static vector<string> cluster_color ;
    /*
      cluster_color.push_back("thistle") ;
      cluster_color.push_back("wheat") ;
      cluster_color.push_back("seashell") ;
      cluster_color.push_back("honeydew") ;
    */
    cluster_color.push_back("blueviolet") ;
    cluster_color.push_back("darkorange") ;
    cluster_color.push_back("chartreuse") ;
    cluster_color.push_back("aquamarine") ;
    static int color_counter = 0 ;
    // output cluster header
    std::string clustername = "subgraph cluster" ;
    std::stringstream ss ;
    ss << baseclustercounter++ ;
    clustername += ss.str() ;

    outf << clustername << " {" << '\n' ;
    // label the cluster
    outf << "label = " << "\"" << rule(sid) << "\"" << ";" << '\n' ;
    //    outf << "style=filled;" << '\n' ;
    outf << "color=" << cluster_color[(color_counter++)%4] << ";" << '\n' ;
    outf << "fontcolor=slateblue;" << '\n' ;
    outf << "fontsize=10;" << '\n' ;
    outf << "fontname=Helvetica;" << '\n' ;

    multiLevelGraph::subGraph *p = mlg.find(sid) ;
    ruleSet level_rules = extract_rules(p->graph_v-mlg.subgraphs) ;
    digraph gr = p->gr ;

    // Remove any rules in that are not in graph_v but are in gr
    ruleSet errs = extract_rules(gr.get_all_vertices()-p->graph_v) ;
    gr.remove_vertices(errs) ;

    digraph grt = gr.transpose() ;

    // do the layout of this subgraph
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    digraph::vertexSet::const_iterator ri ;

    // vertex set that holds all the previously placed vertices
    static digraph::vertexSet placed ;

    for(ri=allvertices.begin();ri!=allvertices.end();++ri) {
      digraph::vertexSet outvertices = gr[*ri] ;
      digraph::vertexSet incomevertices = grt[*ri] ;

      // only list the nodes that belongs to this super node
      // exclude any super nodes
      if( (!placed.inSet(*ri)) && (!is_super_rule(*ri))) {
        placed += *ri ;

        if(*ri < 0) {
          rule r(*ri) ;
          if(r.type() == rule::INTERNAL)
            outf<<"\""<<pretty_sig(r)<<"\""
                <<"[shape=doubleoctagon,style=filled,color=gold];\n" ;
          else
            outf<<"\""<<pretty_sig(r)<<"\""
                <<"[shape=box,style=filled,color=gold];\n" ;
        }
        else {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=green];\n" ;
        }
        if(incomevertices == EMPTY) {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=red];\n" ;
        }
        if(outvertices == EMPTY) {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=blue];\n" ;
          continue ;
        }
      }
    }

    for(ri=allvertices.begin();ri!=allvertices.end();++ri) {

      if(*ri < 0) {
        // recursively layout all the super nodes in this node
        if(is_super_rule(*ri)) {
          layout_super_node(mlg,*ri,outf) ;
        }
      }
    }
    // end of the subgraph cluster
    outf << "}" << '\n' ;
  }

  // connect all the nodes together
  void fill_mlg_edges(multiLevelGraph &mlg, int sid, std::ofstream &outf)
  {
    // get all the levels of the multi level graph
    vector<int> levels ;
    digraph::vertexSet working ;
    working += mlg.toplevel ;
    digraph::vertexSet::const_iterator vi ;
    while(working != EMPTY) {
      digraph::vertexSet next_level ;
      for(vi=working.begin();vi!=working.end();++vi) {
        levels.push_back(*vi) ;
        next_level += mlg.find(*vi)->graph_v & mlg.subgraphs ;
      }
      working = next_level ;
    }

    // define the usual edge color
    std::string edge_color = "slategray" ;

    for(vector<int>::const_iterator levelItr=levels.begin();
        levelItr!=levels.end();++levelItr) {
      multiLevelGraph::subGraph *p = mlg.find(*levelItr) ;
      digraph gr = p->gr ;

      // Remove any rules in that are not in graph_v but are in gr
      ruleSet errs = extract_rules(gr.get_all_vertices()-p->graph_v) ;
      gr.remove_vertices(errs) ;
      digraph grt = gr.transpose() ;

      digraph::vertexSet allvertices = gr.get_all_vertices() ;
      digraph::vertexSet::const_iterator ri ;

      for(ri=allvertices.begin();ri!=allvertices.end();++ri) {
        digraph::vertexSet outvertices = gr[*ri] ;
        digraph::vertexSet incomevertices = grt[*ri] ;

        if(*ri < 0) {// a rule
          //if it is a super rule
          if(is_super_rule(*ri)) {
            writeout_super_rule(mlg,*ri,outvertices,outf) ;
          }else {
            rule r(*ri) ;
            digraph::vertexSet::const_iterator ii ;
            for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
              if(*ii < 0) {// the target is a rule
                if(is_super_rule(*ii)) { // if the rule is a super rule
                  writeout_super_rule2(mlg,*ri,*ii,outf) ;
                }else {// the rule is a common single rule
                  rule r2(*ii) ;
                  outf<<"\""<<pretty_sig(r)<<"\""<<" -> "
                      <<"\""<<pretty_sig(r2)<<"\""
                      <<"[style=bold,color="<<edge_color<<"]"<<";\n" ;
                }
              }else {// the target is a variable
                variable v(*ii) ;
                outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<v<<"\""
                    <<"[style=bold,color="<<edge_color<<"]"<<";\n" ;
              }
            }
          }
        }else {// a variable
          variable v(*ri) ;
          digraph::vertexSet::const_iterator ii ;
          for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
            if(*ii < 0) {//the target is a rule
              if(is_super_rule(*ii)) { // the rule is a super rule
                writeout_super_rule2(mlg,*ri,*ii,outf) ;
              }else {// the rule is a common single rule
                rule r(*ii) ;
                outf<<"\""<<v<<"\""<<" -> "<<"\""<<pretty_sig(r)<<"\""
                    <<"[style=bold,color="<<edge_color<<"]"<<";\n" ;
              }
            }else {// the target is a variable
              variable v2(*ii) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<v2<<"\""
                  <<"[style=bold,color="<<edge_color<<"]"<<";\n" ;
            }
          }
        }
      }
    }
  }

  void create_mlg_dot_file(decomposed_graph &deco, const char* fname)
  {
    multiLevelGraph &mlg = deco.mlg ;

    // out put dot file
    ofstream outf(fname) ;
    outf << "digraph G {\n" ;
    // allow edges between clusters
    outf << "compound=true;\n" ;

    // beginning by layout the toplevel of the multi-level-graph
    layout_super_node(mlg,mlg.toplevel,outf) ;
    fill_mlg_edges(mlg,mlg.toplevel,outf) ;
    // end of digraph
    outf << "}" << '\n' ;
    outf.close() ;
  }

  // visualize the content of each super node
  // from the toplevel down to the bottom
  void visualize_mlg(decomposed_graph &deco)
  {
    multiLevelGraph &mlg = deco.mlg ;
    // get all the levels of the multi level graph
    vector<int> levels ;
    digraph::vertexSet working ;
    working += mlg.toplevel ;
    digraph::vertexSet::const_iterator vi ;
    while(working != EMPTY) {
      digraph::vertexSet next_level ;
      for(vi=working.begin();vi!=working.end();++vi) {
        levels.push_back(*vi) ;
        next_level += mlg.find(*vi)->graph_v & mlg.subgraphs ;
      }
      working = next_level ;
    }

    int levelcounter = 1 ;
    for(vector<int>::const_iterator levelItr=levels.begin();
        levelItr!=levels.end();++levelItr) {
      multiLevelGraph::subGraph *p = mlg.find(*levelItr) ;
      digraph gr = p->gr ;

      // Remove any rules in that are not in graph_v but are in gr
      ruleSet errs = extract_rules(gr.get_all_vertices()-p->graph_v) ;
      gr.remove_vertices(errs) ;

      std::string sname = "levels" ;
      std::stringstream ss ;
      ss << levelcounter++ ;
      sname.append(ss.str()) ;
      sname.append(".dot") ;

      create_digraph_dot_file(gr,sname.c_str()) ;
      std::string cmd = "dotty " ;
      cmd += sname ;
      //cmd += "&" ;
      cout << "This is super node: " << "SN" << levelcounter-2 << '\n' ;
      int err = system(cmd.c_str()) ;
      if(err != 0)
	cerr << "system call returned " << err << " in call '" << cmd << "'" 
	     << endl ;

      // rm the generated "dot" file
      cmd = "rm -fr " ;
      cmd += sname ;
      err = system(cmd.c_str()) ;
      if(err != 0)
	cerr << "system call returned " << err << " on system('"
	     << cmd << "')" << endl ;
    }
  }

  void prune_graph(digraph& gr, variableSet& given,
                   const variableSet& target, fact_db& facts) {
    variableSet known_vars = facts.get_typed_variables() ;
    variableSet vars = extract_vars(gr.get_all_vertices()) ;
    variableSet empty_constraints ;
    for(variableSet::const_iterator vi=vars.begin();
        vi!=vars.end();++vi) {
      if(!known_vars.inSet(*vi))
        continue ;
      storeRepP srp = facts.get_variable(*vi) ;
      if(srp->RepType() != Loci::CONSTRAINT)
        continue ;
      if(GLOBAL_AND(srp->domain() == EMPTY))
        empty_constraints += *vi ;
    }
    ruleSet del_rules ;
    for(variableSet::const_iterator vi=empty_constraints.begin();
        vi!=empty_constraints.end();++vi)
      del_rules += extract_rules(gr[vi->ident()]) ;
    gr.remove_vertices(digraph::vertexSet(empty_constraints) +
                       digraph::vertexSet(del_rules)) ;
    given -= empty_constraints ;
    clean_graph(gr,given,target) ;
  }

  extern rule_db rename_gpu_containers(fact_db  &facts,
				       const rule_db &rdb) ;
#define ENABLE_RELATION_GEN
  executeP create_execution_schedule(const rule_db &rdb,
                                     fact_db &facts,
                                     sched_db &scheds,
                                     const variableSet& target,
                                     int nth) {
    variableSet parVars = target ;
    parVars += facts.get_extensional_facts() ;
    rule_db par_rdb ;
    par_rdb = parametric_rdb(rdb,parVars) ;
    // Not sure exactly why we need this, perhaps this is used by the
    // dynamic execution feature?  CHEM passes quicktest without
    // it.
    par_rdb = replace_map_constraints(facts,par_rdb) ;

    // Rename gpu containers and insert
    par_rdb = rename_gpu_containers(facts,par_rdb) ;
    ////////////////decorate the dependency graph/////////////////////
    //if(Loci::MPI_rank==0)
    //cout << "decorating dependency graph to include allocation..." << endl ;
    //deco_depend_gr(gr,given) ;
    //////////////////////////////////////////////////////////////////
#ifdef ENABLE_RELATION_GEN
    if(Loci::MPI_rank==0)
      cout << "Stationary Relation Generation..." << endl ;
    stationary_relation_gen(par_rdb, facts, target) ;
#endif
    // then we need to perform global -> local renumbering
    if(facts.is_distributed_start()) {
      if((MPI_processes > 1))
        get_clone(facts, par_rdb) ;
      else
        Loci::serial_freeze(facts) ;
    } else {
      Loci::serial_freeze(facts) ;
    }

    // then we can generate the dependency graph
    variableSet given = facts.get_typed_variables() ;
    if(Loci::MPI_rank==0)
      cout << "generating dependency graph..." << endl ;

    
    stopWatch sw ;
    sw.start() ;
	
    digraph gr ;

    given -= variable("EMPTY") ;

    gr = dependency_graph2(par_rdb,given,target).get_graph() ;

    // If graph is empty, return a null schedule
    if(gr.get_target_vertices() == EMPTY) {
      if(Loci::MPI_rank == 0)
        cerr << "Warning: empty dependency graph!" << endl ;
      return executeP(0) ;
    }
    ////////////////////////////////////////////////////////////////////////
    std::string dottycmd = "dotty " ;
    if(Loci::MPI_rank==0) {
      if(show_graphs) {
        cout << "creating visualization file for dependency graph..." << endl ;
        create_digraph_dot_file(gr,"dependgr.dot") ;
        std::string cmd = dottycmd + "dependgr.dot" ;
        int err = system(cmd.c_str()) ;
	if(err != 0)
	  cerr << "system call returned " << err << " on system('"
	       << cmd << "')" << endl ;
	
      }
    }
    ////////////////////////////////////////////////////////////////////////

    ////////////////////
    //prune_graph(gr,given,target,facts) ;
    ////////////////////

    scheds.init(facts) ;
    if(Loci::MPI_rank==0)
      cout << "setting up variable types..." << endl ;
    set_var_types(facts,gr,scheds) ;

    //////////////
    //scheds.print_summary(facts,cout) ;
    //////////////

    if(Loci::MPI_rank==0)
      cout << "decomposing graph..." << endl ;
    decomposed_graph decomp(gr,given,target) ;

    //////////////////////////////////////////////////////////////////
    if(Loci::MPI_rank==0) {
      if(show_graphs) {
        cout << "creating visualization file for decomposed graph..." << '\n' ;
        cout << "visualizing decomposed graph..." << '\n' ;
        create_mlg_dot_file(decomp,"decogr.dot") ;
        std::string cmd = dottycmd + "decogr.dot&" ;
        int err = system(cmd.c_str()) ;
	if(err != 0)
	  cerr << "system call returned " << err << " on system('"
	       << cmd << "')" << endl ;
        visualize_mlg(decomp) ;
      }
    }

    //////////////////////////////////////////////////////////////////
    variableSet fact_vars, initial_vars ;
    fact_vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=fact_vars.begin();vi!=fact_vars.end();++vi) {
      storeRepP vp = facts.get_variable(*vi) ;
      //Existence of a map is actually its domain on a processor.
      //It is necessary that existence of a map does not only include the subeset
      //of my_entities, but also includes the clone entities.
      //Mainly because, we have some constraints that are applied over maps.
      //If the map is actually not used in the rule other than its constraints,
      //then that map may not have expanded enough to include necessrary existence.
      //For regular execution it won't affect the schedule but for duplication of work,
      //it is required for saving communication.
      if(isMAP(vp)) {
	if(facts.isDistributed()) {
	  entitySet exist = scheds.variable_existence(*vi);
	  exist = fill_entitySet(exist, facts);
	  scheds.set_variable_existence(*vi, exist);
	}
      }
      if(variable(*vi).time().level_name() == "*" ) {
	if(isSTORE(vp)) {
	  ostringstream oss ;
	  oss << "source(" <<"EMPTY"<<')' ;
	  oss << ",target(" << *vi << ')' ;
	  string sig = oss.str() ;
	  rule r(sig) ;
	  if(par_rdb.rules_by_target(*vi) == EMPTY) {
	    if(facts.isDistributed()) {
	      scheds.set_existential_info(*vi, r, scheds.variable_existence(*vi));
	      initial_vars += *vi ;
	    }
	  }
	}
      }
    }
    Loci::debugout << " initial_vars = " << initial_vars << endl ;

    if(duplicate_work) {
      if(use_duplicate_model) {
	std::ifstream fin(model_file);
	if(!fin) {
	  cerr << "Error: Opening the model file " << model_file << endl;
	  cerr << "Using default duplication policies." << endl;
	  use_duplicate_model = false;
	}
	
	if(use_duplicate_model) {
	      double comm_ts, comm_tw;
	      double comm_ts1, comm_ts2;
	      fin >> comm_ts1 >> comm_ts2 >> comm_tw;
	      comm_ts = comm_ts2;
	      
	      if(comm_tw < 0)
		    comm_tw = 0;

	      unsigned int count;
	      fin >> count;

	      map<rule, pair<double, double> > comp_info;
	      string  rule_name;
	      double ts, tw;
	      double ts1, ts2;
	      for(unsigned int i = 0; i < count; i++) {
		    fin >> rule_name >> ts1 >> ts2 >> tw;
		    ts = ts2;
		    if(tw < 0)
			  tw = 0;
		    
		    pair<double, double> tmpModel(ts, tw);
		    rule myRule = rule::get_rule_by_name(rule_name);
		    if(myRule.get_info().name() == "NO_RULE") {
			  cerr << "Warning (Rule Ignored): " << rule_name << " read from model file is not in rule database" << endl;
		    }
		    else
			  comp_info[myRule] = tmpModel;
	      }
	      scheds.add_model_info(comm_ts, comm_tw, comp_info);	
	}
      }
    }

    graph_compiler compile_graph(decomp, initial_vars) ;
    compile_graph.compile(facts,scheds,given,target) ;
	
    Loci::debugout << "Time taken for graph processing  = "
                   << sw.stop() << "  seconds " << endl ;

    if(Loci::MPI_rank==0)
      cout << "existential analysis..." << endl ;
    sw.start() ;
      
    compile_graph.existential_analysis(facts, scheds) ;
      
    Loci::debugout << "Time taken for existential_analysis  = "
                   << sw.stop() << "  seconds " << endl ;
    ///////////////////////////////////
    /*
      if(Loci::MPI_rank==0) {
      ofstream exinfo(".exis") ;
      scheds.print_summary(facts,exinfo) ;
      exinfo.close() ;
      }
    */
    ///////////////////////////////////
    if(Loci::MPI_rank==0) {
#ifdef PTHREADS
      if(threading_pointwise || threading_global_reduction
         || threading_local_reduction || threading_chomping 
         || threading_recursion) {
        cout << "creating multithreaded execution schedule ("
             << num_threads << " threads per MPI process, "
             << num_thread_blocks << " blocks each)" << endl;
        cout << "--threading suitable ";
        if(threading_pointwise)
          cout << "[pointwise] ";
        if(threading_global_reduction)
          cout << "[global reduction] ";
        if(threading_local_reduction)
          cout << "[local reduction] ";
        if(threading_chomping)
          cout << "[chomping] ";
        if(threading_recursion)
          cout << "[recursive] ";
        cout << "rules" << endl;
      } else 
#endif
        cout << "creating execution schedule..." << endl;
    }
    sw.start() ;
    
    executeP sched =  compile_graph.execution_schedule
      (facts,scheds,initial_vars) ;
    Loci::debugout << "Time taken for create execution schedule = "
                   << sw.stop() << " seconds " << endl ;

    if(GLOBAL_OR(scheds.errors_found())) {
      if(MPI_rank == 0) {
        cerr << "error in generating schedule, dumping schedule files" << endl ;
        if(MPI_processes != 1)
          cerr << "see debug files for more information" << endl ;

      }
      ostringstream oss ;
      oss << "debug/schedule" ;

      if(MPI_processes > 1) {
        oss << "-" << MPI_rank ;
      }

      string sched_filename = oss.str() ;
      ofstream sched_file(sched_filename.c_str(),ios::out) ;
      sched->Print(sched_file) ;
      sched_file.close() ;


      Loci::Abort() ;
    }
    //scheds.print_summary(facts,Loci::debugout) ;

    // setting this external pointer
    exec_current_fact_db = &facts ;
    return sched ;
  }


  // get profiling information from schedule
  class collectTiming : public collectData {
    struct timingData {
      executeEventType eventType ;
      std::string groupName ;
      std::string eventName ;
      timeAccumulator accumTime ;
      // Added for analysis across processors
      double totalTime, maxTime, meanTime ;
      double totalEvents ;
      double maxEvents ;
bool operator <(const timingData &d) const {
        return (max(accumTime.getTime(),maxTime) < 
		max(d.accumTime.getTime(),d.maxTime)) ;
      }
      timingData() : eventType(EXEC_CONTROL),totalTime(0),maxTime(0),meanTime(0),totalEvents(0),maxEvents(0) {}
    } ;
    std::list<timingData>  timing_data ;

    struct schedData {
      std::string eventName;
      double bytes;
      bool operator<(const schedData& s) const
      { return bytes < s.bytes; }
    };
    std::list<schedData> sched_data;
    struct cmpSchedName {
      bool operator()(const schedData& s1, const schedData& s2) const
      { return s1.eventName < s2.eventName; }
    };
    struct cacheData {
      std::string eventName;
      long_long l1_dcm;
      long_long l2_dcm;
      bool operator<(const cacheData& c) const
      { return l1_dcm < c.l1_dcm; }
    };
    std::list<cacheData> cache_data;
  public:

    void accumulateTime(const timeAccumulator &ta, executeEventType t,
                        string eventName) ;
    void accumulateMemory(const std::string &var,
                          allocEventType t,
                          double maxMallocMemory,
                          double maxBeanMemory) ;
    void accumulateSchedMemory(const std::string& eventName,
                               double bytes);
    void accumulateDCM(const std::string& eventName,
                       long_long l1_dcm, long_long l2_dcm);
      
    double getComputeTime() ;
    double getTotalTime() ;
    
    void balanceAnalysis(MPI_Comm comm) ;

    ostream &PrintSummary(ostream &s) ;
    
  } ;
  
  void collectTiming::accumulateTime(const timeAccumulator &ta, executeEventType t, string eventName) {
    timingData td ;
    td.eventType =
      t ;
    td.eventName = eventName ;
    td.accumTime = ta ;
    if(!groups.empty())
      td.groupName = groups[0] ;
    
    for(size_t i=1;i<groups.size();++i) {
      td.groupName += ':' ;
      td.groupName += groups[i] ;
    }
    timing_data.push_back(td) ;
  }

  void collectTiming::accumulateMemory(const std::string &var,
                                       allocEventType t,
                                       double maxMallocMemory,
                                       double beanMemory) {
  }

  void collectTiming::accumulateSchedMemory(const std::string& eventName,
                                            double bytes)
  {
    schedData sd;
    sd.eventName = eventName;
    sd.bytes = bytes;
    sched_data.push_back(sd);
  }

  void collectTiming::accumulateDCM(const std::string& eventName,
                                    long_long l1_dcm, long_long l2_dcm)
  {
    cacheData cd;
    cd.eventName = eventName;
    cd.l1_dcm = l1_dcm;
    cd.l2_dcm = l2_dcm;
    cache_data.push_back(cd);
  }

  void collectTiming::balanceAnalysis(MPI_Comm comm)  {
    int np = 1 ;
    int r =  0; 
    MPI_Comm_size(comm,&np) ;
    MPI_Comm_rank(comm,&r) ;
    map<string,timingData> dataList ;
    map<string,timingData>::iterator lp ;
    std::list<timingData>::iterator  ti = timing_data.begin() ;
    // Collect timing data on this processor
    for(ti=timing_data.begin();ti!=timing_data.end();++ti) {
      
      if((lp = dataList.find(ti->eventName)) != dataList.end()) {
	debugout << "balanceAnalysis: event " << ti->eventName << " is duplicate!" << endl ;
	lp->second.accumTime.addTime(ti->accumTime.getTime(),ti->accumTime.getEvents()) ;
      } else {
	dataList[ti->eventName] = *ti ;
      }
    }
    // Not all processors have the same data so we need to make it consistent
    set<string> unprocessed ;
    set<string>::iterator si ;
    for(lp=dataList.begin();lp!=dataList.end();++lp)  {
      unprocessed.insert(lp->first) ;
    }

    // start with process 0 as the source of entry types
    int psource = 0 ;
    
    // process unprocessed data until all processors have no unprocessed
    // data
    bool processmore = false ;
    do {
      int lsz = unprocessed.size() ;
      MPI_Bcast(&lsz,1,MPI_INT,psource,comm) ;
      vector<int> sizes(lsz) ;
      if(r == psource) {
	int i = 0 ;
	for(si=unprocessed.begin();si!=unprocessed.end();++si)
	  sizes[i++] = si->size() ;
      }
      MPI_Bcast(&sizes[0],lsz,MPI_INT,psource,comm) ;
      int tot = 0 ;
      for(int i=0;i<lsz;++i)
	tot += sizes[i] ;
      vector<char> data(tot) ;
      if(r == psource) {
	int ii = 0 ;
	int i = 0 ;
	for(si=unprocessed.begin();si!=unprocessed.end();++si) {
	  string val = *si ;
	  for(int j=0;j<sizes[i];++j)
	    data[ii++] =val[j] ;
	  i++ ;
	}
      }
      MPI_Bcast(&data[0],tot,MPI_BYTE,psource,comm) ;
      int ii = 0 ;
      for(int i=0;i<lsz;++i) {
	string event ;
	for(int j=0;j<sizes[i];++j)
	  event += data[ii++] ;
	if((si = unprocessed.find(event)) != unprocessed.end()) {
	  unprocessed.erase(si) ;
	} else { // Not in my list so add entry to dataList
	  dataList[event] = timingData() ;
	}
      }
      int pcandidate = np ;
      if(unprocessed.size() != 0)
	pcandidate = r ;
      MPI_Allreduce(&pcandidate,&psource, 1, MPI_INT,MPI_MIN,comm) ;

      processmore = (psource != np) ;
			    
    }  while (processmore) ;

    // Now that dataList is consistent across all processors, lets gather
    // data from all processors to get idea of load imbalances not visible from
    // any one processor
    int datasize = dataList.size() ;
    
    vector<double> localTimes(datasize) ;
    vector<double> localEvents(datasize) ;
    int i = 0 ;
    for(lp=dataList.begin();lp!=dataList.end();++lp)  {
      localTimes[i] = lp->second.accumTime.getTime() ;
      localEvents[i] = lp->second.accumTime.getEvents() ;
      i++ ;
    }
    vector<double> totalTimes(datasize) ;
    vector<double> maxTimes(datasize) ;
    vector<double> totalEvents(datasize) ;
    vector<double> maxEvents(datasize) ;
    MPI_Allreduce(&localTimes[0],&totalTimes[0],datasize,MPI_DOUBLE,MPI_SUM,comm) ;
    MPI_Allreduce(&localTimes[0],&maxTimes[0],datasize,MPI_DOUBLE,MPI_MAX,comm) ;
    MPI_Allreduce(&localEvents[0],&totalEvents[0],datasize,MPI_DOUBLE,MPI_SUM,comm) ;
    MPI_Allreduce(&localEvents[0],&maxEvents[0],datasize,MPI_DOUBLE,MPI_MAX,comm) ;

    // Replace timing_data with new combined data list
    i = 0 ;
    list<timingData> newList ;
    for(lp=dataList.begin();lp!=dataList.end();++lp)  {
      lp->second.totalTime = totalTimes[i] ;
      lp->second.meanTime = totalTimes[i]/np ;
      lp->second.maxTime = maxTimes[i] ;
      lp->second.totalEvents = totalEvents[i] ;
      lp->second.maxEvents = maxEvents[i] ;
      i++ ;
      newList.push_back(lp->second) ;
    }
    timing_data.swap(newList) ;
  }

  ostream &collectTiming::PrintSummary(ostream &s) {
    timing_data.sort() ;
    double totComp = 0 ;
    double totComm = 0 ;
    double totCtrl = 0 ;
    std::list<timingData>::const_iterator  ti = timing_data.begin() ;
    for(ti=timing_data.begin();ti!=timing_data.end();++ti)
      switch(ti->eventType) {
      case EXEC_COMMUNICATION:
        totComm += ti->accumTime.getTime() ;
        break ;
      case EXEC_COMPUTATION:
        totComp += ti->accumTime.getTime() ;
        break ;
      case EXEC_CONTROL:
        totCtrl += ti->accumTime.getTime() ;
        break ;
      }

    s << "------------------------------------------------------------------------------" << endl ;
    s << "Timing Categories:" << endl ;
    s << " -- Computation:   " << totComp << endl ;
    s << " -- Communication: " << totComm << endl ;
    s << " -- Control:       " << totCtrl << endl ;
    double totTime = totComp+totComm+totCtrl+1e-30 ;
    s << " -- totalTime:     " << totTime << endl ;
    s << endl ;
    s << "------------------------------------------------------------------------------" << endl ;

    std::list<timingData>::const_reverse_iterator rti = timing_data.rbegin() ;
    int lcnt = 10 ;
    if(verbose)
      lcnt = timing_data.size() ;
    s << "Top " << lcnt << " Most Expensive Steps:" << endl ;
    lcnt = min(int(timing_data.size()),lcnt) ;
    for(int i =0;i<lcnt;++i,++rti) {
      s << i << "- " << rti->eventName << endl ;
      if(rti->groupName != "")
        s << " --- Group " << rti->groupName << endl ;

      double t = rti->accumTime.getTime() ;
      double e = double(rti->accumTime.getEvents()) ;
      s << " --- Local Time: " << t << " "
        << ceil(1000.0*t/totTime)/10.0 << "% of total," 
        <<  " time per entity: " << t/max(e,1.0)
        << endl ;
      double meanEvents = rti->totalEvents/double(MPI_processes) ;
      s << " === max " << rti->maxTime << ", mean = " << rti->meanTime
	<< ", imbalance = " << 100.0*(rti->maxTime-rti->meanTime)/max(rti->meanTime,1e-10)<<"%" << endl	 ;
      if(rti->eventType == EXEC_COMPUTATION) 
	s << " === partition imbalance =" <<  100.0*(rti->maxEvents-meanEvents)/max(meanEvents,1.0) << "%" 
	  << ", mean time per entity = " << rti->totalTime/max(rti->totalEvents,1.0) << endl;

      // DEBUG
      s << " --- Type: " ;
      switch(rti->eventType) {
      case EXEC_COMMUNICATION:
        s << "Communication" ;
        break ;
      case EXEC_COMPUTATION:
        s << "Computation" ;
        break ;
      case EXEC_CONTROL:
        s << "Control" ;
        break ;
      }
      s << endl ;
      ////////      
      s << "------------------------------------------------------------------------------" << endl ;
    }

    map<string,double> group_times ;
    map<string,double>::const_iterator gi ;
    for(ti=timing_data.begin();ti!=timing_data.end();++ti) {
      string group = ti->groupName ;
      vector<string> groups ;
      string working ;
      for(size_t i=0;i<group.size();++i)
        if(group[i] == ':') {
          if(working != "")
            groups.push_back(working) ;
          working = "" ;
        } else {
          working += group[i] ;
        }
      if(working != "")
        groups.push_back(working) ;
      for(size_t i=0;i<groups.size();++i) 
        group_times[groups[i]] += ti->accumTime.getTime() ;
    }

    s << "----- Time per category:" << endl ;
    for(gi=group_times.begin();gi!=group_times.end();++gi) {
      double t = gi->second ;
      s << "Group " << gi->first << " time = " << gi->second << ", " 
        << ceil(1000.0*t/totTime)/10.0 << "% of total" << endl ;
    }    
    s << "------------------------------------------------------------------------------" << endl ;

    // display the schedule memory consumption
    // first collapse the list by merging names together
    sched_data.sort(cmpSchedName());
    std::list<schedData>::iterator si = sched_data.begin();
    std::list<schedData>::iterator si2 = si; ++si2;
    while(si!=sched_data.end()) {
      if(si2 == sched_data.end())
        break;
      if(si->eventName == si2->eventName) {
        si->bytes += si2->bytes;
        si2 = sched_data.erase(si2);
      } else {
        si = si2;
        ++si2;
      }
    }
    // then sort based on bytes size
    sched_data.sort();
    double total_mem = 0;
    for(si=sched_data.begin();si!=sched_data.end();++si)
      total_mem += si->bytes;
    s << "Total scheduler bean-counting memory: " << total_mem << " bytes"
      << endl;
    
    std::list<schedData>::const_reverse_iterator rsi = sched_data.rbegin();
    lcnt = 10 ;
    if(verbose)
      lcnt = sched_data.size() ;
    s << "Top " << lcnt << " Steps with Most Scheduler Memory:" << endl;
    lcnt = min(int(sched_data.size()),lcnt);
    for(int i =0;i<lcnt;++i,++rsi) {
      s << i << "- " << rsi->eventName
        << " | " << rsi->bytes << " bytes" << endl ;
      s << "------------------------------------------------------------------------------" << endl ;
    }
    
#ifdef PAPI_DEBUG
    // display the cache misses
    cache_data.sort();
    std::list<cacheData>::const_reverse_iterator rci = cache_data.rbegin();
    lcnt = 10 ;
    if(verbose)
      lcnt = cache_data.size() ;
    s << "Top " << lcnt << " Steps with Most Data Cache Misses:" << endl;
    lcnt = min(int(cache_data.size()),lcnt);
    for(int i =0;i<lcnt;++i,++rci) {
      s << i << "- " << rci->eventName
        << " | L1 DCM: " << rci->l1_dcm
        << " | L2 DCM: " << rci->l2_dcm << endl ;
      s << "------------------------------------------------------------------------------" << endl ;
    }
#endif
    
    return s ;
  }
  
  double collectTiming::getComputeTime() {
    double totComp = 0 ;
    std::list<timingData>::const_iterator  ti = timing_data.begin() ;
    for(ti=timing_data.begin();ti!=timing_data.end();++ti)
      if(ti->eventType==EXEC_COMPUTATION || ti->eventType == EXEC_CONTROL) 
        totComp += ti->accumTime.getTime() ;
    return totComp ;
  }  

  double collectTiming::getTotalTime() {
    double totComp = 0 ;
    std::list<timingData>::const_iterator  ti = timing_data.begin() ;
    for(ti=timing_data.begin();ti!=timing_data.end();++ti)
      totComp += ti->accumTime.getTime() ;
    return totComp ;
  }  
  
  // get profiling information from schedule
  class collectMemory : public collectData {
    struct spaceData {
      double max_memory ;
      string event_name ;
      allocEventType type ;
      bool operator <(const spaceData &d) const {
        return max_memory < d.max_memory ;
      }
    } ;
    struct var_size_info {
      double mem_size ;
      double event_tick ;
      var_size_info(){mem_size=0;event_tick=0;}
      var_size_info(double x, double y) : mem_size(x),event_tick(y) {}
    } ;
    struct alloc_event {
      variable alloc_var ;
      variableSet live_set ;
      double live_mem ;
      bool operator <(const alloc_event &d) const {
        return live_mem < d.live_mem ;
      }
      alloc_event() { }
      alloc_event(const variableSet &ls, double mem) : live_set(ls),live_mem(mem)
      {}
    } ;
    std::list<spaceData>  space_data ;
    std::map<variable,var_size_info> variable_size ;
    double event_count ;
  public:
    collectMemory() { event_count = 0 ; }
      
    void accumulateTime(const timeAccumulator &ta, executeEventType t,
                        string eventName) ;
    void accumulateMemory(const std::string &var,
                          allocEventType t,
                          double maxMallocMemory,
                          double maxBeanMemory) ;
      
    void accumulateSchedMemory(const std::string& eventName,
                               double bytes) {}
    void accumulateDCM(const std::string& eventName,
                       long_long l1_dcm, long_long l2_dcm) {}

    double getComputeTime() ;
    double getTotalTime() ;
    ostream &PrintSummary(ostream &s) ;
    
  } ;
  
  void collectMemory::accumulateTime(const timeAccumulator &ta, executeEventType t, string eventName) {
  }
  
  void collectMemory::accumulateMemory(const std::string &var,
                                       allocEventType t,
                                       double maxMallocMemory,
                                       double beanMemory) {
    event_count += 1.0 ;
    spaceData sd ;
    sd.event_name = var ;
    sd.max_memory = maxMallocMemory ;
    sd.type = t ;
    space_data.push_back(sd) ;
    variable v(var) ;
    if(t == ALLOC_CREATE) {
      variable_size[v] = var_size_info(beanMemory,-event_count) ;
    }
    if(t == ALLOC_DELETE) {
      // Update lifetime and memory space
      double create_count = variable_size[v].event_tick ;
      variable_size[v] = var_size_info(beanMemory,event_count+create_count) ;
    }
  }
  
  ostream &collectMemory::PrintSummary(ostream &s) {
    // Now scan through the data to get memory totals
    double tot_memory = 0 ;
    list<spaceData>::const_iterator li ;
    list<alloc_event> l2 ;
    variableSet active_vars ;

    for(li=space_data.begin();li!=space_data.end();++li) {
      variable v(li->event_name) ;
      if(li->type == ALLOC_CREATE) {
        tot_memory += variable_size[v].mem_size ;
        active_vars += v ;
      }
      if(li->type == ALLOC_DELETE) {
        alloc_event ae ;
        ae.live_set = active_vars ;
        ae.live_mem = tot_memory ;
        ae.alloc_var = v ;
        l2.push_back(ae) ;
        tot_memory -= variable_size[v].mem_size ;
        active_vars -= v ;
      }
    }
    l2.sort();

    std::list<alloc_event>::const_reverse_iterator rti = l2.rbegin() ;
    int lcnt = 10 ;
    if(verbose)
      lcnt = l2.size() ;
    s << "Top " << lcnt << " Most Memory Expensive Steps:" << endl ;
    lcnt = min(int(l2.size()),lcnt) ;
    for(int i =0;i<lcnt;++i,++rti) {
      // Output variable that was allocated during max allocation
      s << i << "- allocate var " << rti->alloc_var << " tot mem = " 
        << rti->live_mem/(1024.0*1024.0) << "MB" << endl ;
      // Now sort liveset according to procesor-time product
      vector<pair<double,variable> > live_list ;
      variableSet live_set = rti->live_set ;
      variableSet::const_iterator vi ;
      for(vi=live_set.begin();vi!=live_set.end();++vi) {
        double mt = variable_size[*vi].mem_size * variable_size[*vi].event_tick ;
        live_list.push_back(pair<double,variable>(mt,*vi)) ;
      }
      std::sort(live_list.begin(),live_list.end()) ;
      int sz = live_list.size() ;
      for(int i=sz-1; i>=0;--i) {
        variable v = live_list[i].second ;
        if(verbose || variable_size[v].mem_size > 2048)
          s << "    " << v << ":" << variable_size[v].mem_size/1024.0 << "k,"
            << " lifetime=" << variable_size[v].event_tick << endl ;
      }
    }
    return s ;
  }

  // this function is used to create an execution schedule
  // for internalQuery below. This function and the internalQuery
  // function are mainly intended to be used by the Loci scheduler
  // to compute intermediate values.

  // NOTE: the passed in rule database is the expanded parametric
  // rule database since the expanded rule base is what we needed
  // and this expansion process only needs to be performed once.
  //#define INTERNAL_VERBOSE
  executeP create_internal_execution_schedule(rule_db& par_rdb,
                                              fact_db &facts,
                                              sched_db &scheds,
                                              const variableSet& target,
                                              int nth) {
    // since this function is always executed inside
    // the create_execution_schedule function so the
    // fact database is always in the local number state
    // thus we don't need to perform the global -> local
    // renumbering step

    variableSet given = facts.get_typed_variables() ;
#ifdef INTERNAL_VERBOSE
    if(Loci::MPI_rank==0) {
      cout << "[Internal] generating dependency graph..." << endl ;
    }
#endif
    // as the usual process, we'll need to generate
    // the dependency graph
    digraph gr ;
#ifdef INTERNAL_VERBOSE
    if(Loci::MPI_rank==0)
      cout << "\t[Internal] (recursive backward searching version)" << endl ;
#endif
    given -= variable("EMPTY") ;
    gr = dependency_graph2(par_rdb,given,target).get_graph() ;
    // If graph is empty, return a null schedule
    if(gr.get_target_vertices() == EMPTY)
      return executeP(0) ;

    
//         std::string dottycmd = "dotty " ;
//         if(Loci::MPI_rank==0) {
//           if(show_graphs) {
//             cout << "creating visualization file for dependency graph..." << endl ;
//             create_digraph_dot_file(gr,"dependgr.dot") ;
//             std::string cmd = dottycmd + "dependgr.dot" ;
//             system(cmd.c_str()) ;
//           }
//         }
    

    scheds.init(facts) ;
#ifdef INTERNAL_VERBOSE
    if(Loci::MPI_rank==0)
      cout << "[Internal] setting up variable types..." << endl ;
#endif
    set_var_types(facts,gr,scheds) ;

#ifdef INTERNAL_VERBOSE
    if(Loci::MPI_rank==0)
      cout << "[Internal] decomposing graph..." << endl ;
#endif
    decomposed_graph decomp(gr,given,target) ;

#ifdef INTERNAL_VERBOSE
    if(Loci::MPI_rank==0) {
      cerr << "[Internal] setting initial variables..." << endl ;
    }
#endif
    variableSet fact_vars, initial_vars ;
    fact_vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=fact_vars.begin();
        vi!=fact_vars.end();++vi) {
      storeRepP vp = facts.get_variable(*vi) ;
      //Existence of a map is actually its domain on a processor.
      //It is necessary that existence of a map does not only include
      //the subeset of my_entities, but also includes the clone entities.
      //Mainly because, we have some constraints that are applied over maps.
      //If the map is actually not used in the rule other than its constraints,
      //then that map may not have expanded enough to include necessrary
      //existence. For regular execution it won't affect the schedule but
      //for duplication of work, it is required for saving communication.
      if(isMAP(vp)) {
	if(facts.isDistributed()) {
	  entitySet exist = scheds.variable_existence(*vi);
	  exist = fill_entitySet(exist, facts);
	  scheds.set_variable_existence(*vi, exist);
	}
      }
      if(variable(*vi).time().level_name() == "*" ) {
	if(isSTORE(vp)) {
	  ostringstream oss ;
	  oss << "source(" <<"EMPTY"<<')' ;
	  oss << ",target(" << *vi << ')' ;
	  string sig = oss.str() ;
	  rule r(sig) ;
	  if(par_rdb.rules_by_target(*vi) == EMPTY) {
	    if(facts.isDistributed()) {
	      scheds.set_existential_info(*vi, r,
                                          scheds.variable_existence(*vi));
	      initial_vars += *vi ;
	    }
	  }
	}
      }
    }
    //    Loci::debugout << " initial_vars = " << initial_vars << endl ;

#ifdef INTERNAL_VERBOSE
    if(Loci::MPI_rank==0)
      cerr << "[Internal] compiling graph..." << endl;
#endif
    graph_compiler compile_graph(decomp, initial_vars) ;
    compile_graph.compile(facts,scheds,given,target) ;

#ifdef INTERNAL_VERBOSE
    if(Loci::MPI_rank==0)
      cout << "[Internal] existential analysis..." << endl ;
#endif
    compile_graph.existential_analysis(facts, scheds) ;

#ifdef INTERNAL_VERBOSE
    if(Loci::MPI_rank==0)
      cout << "[Internal] creating execution schedule..." << endl;
#endif
    executeP sched =  compile_graph.execution_schedule
      (facts,scheds,initial_vars) ;

    if(GLOBAL_OR(scheds.errors_found())) {
      if(MPI_rank == 0) {
        cerr << "[Internal] error in generating schedule" << endl ;
        if(MPI_processes != 1)
          cerr << "[Internal] see debug files for more information" << endl ;
        cerr << "[Internal] Aborting..." << endl ;

      }
      Loci::Abort() ;
    }
    exec_current_fact_db = &facts ;
    return sched ;
  }

  // this function is used by the Loci scheduler to issue
  // queries for intermediate relations. It is basically a
  // reduced version of the user function makeQuery
  // NOTE: the passed in rule database is the expanded parametric
  // rule database since the expanded rule base is what we needed
  // and this expansion process only needs to be performed once.
  bool internalQuery(rule_db& par_rdb, fact_db& facts,
                     const variableSet& query) {
    stopWatch sw ;
    sw.start() ;
    
    if(MPI_rank == 0) {
      cout << "[Internal] Quering facts: " << query << endl ;
    }
    // Here, we won't erase the intentional facts since
    // we are using them to provide schedules later

    // But we do need to copy the fact_db before we
    // start to make the query
    // This is because we want to only put the queried facts
    // back into the global fact_db
    fact_db local_facts(facts) ;
    sched_db local_scheds ;

    executeP schedule =
      create_internal_execution_schedule(par_rdb,
                                         local_facts,local_scheds,query) ;
    if(schedule == 0)
      return false ;

    // If a schedule was generated, execute it
#ifdef INTERNAL_VERBOSE
    if(MPI_rank == 0)
      cout << "[Internal] begin query execution" << endl ;
#endif
    exec_current_fact_db = &local_facts ;
    schedule->execute(local_facts, local_scheds) ;

    for(variableSet::const_iterator vi=query.begin();
        vi!=query.end();++vi) {
      storeRepP srp = local_facts.get_variable(*vi) ;
      facts.create_intensional_fact(*vi,srp) ;
    }
    double tlocal = sw.stop() ;
    double tglobal = 0 ;
    MPI_Allreduce(&tlocal,&tglobal, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD) ;
    debugout << "time to execute internal query " << tglobal <<endl;

    return true ;
  }

  bool makeQuery(const rule_db &rdb, fact_db &facts,
                 const std::string& query) {
	/*	
	  #ifdef USE_PAPI
	  int perr,ev_set=PAPI_NULL;
	  int i,ncnt,k;
	  if(PAPI_VER_CURRENT!=(perr=PAPI_library_init(PAPI_VER_CURRENT)))
	  cerr<<"\nerror during initialization\n";
	  unsigned char v[N];
	  long_long counts[NCOUNTS];
	  int evlist[NCOUNTS];
	  char evname[NCOUNTS][PAPI_MAX_STR_LEN];
	  int retval;
	  #endif
	*/
    facts.setupDefaults(rdb) ;
    stopWatch sw ;
    sw.start() ;

    try {
      if(MPI_rank == 0) {
        cout << "Quering facts: " << query << endl ;
      }

      // first check whether queried facts are extensional facts
      // if so, we don't need to actually perform query on
      // these extensional facts
      variableSet target(expression::create(query)) ;
      variableSet efacts = facts.get_extensional_facts() ;
      variableSet remove_query ;
      for(variableSet::const_iterator vi=target.begin();
          vi!=target.end();++vi) {
        if(efacts.inSet(*vi))
          remove_query += *vi ;
      }
      target -= remove_query ;

      if(remove_query != EMPTY)
        if(MPI_rank == 0) {
          cout << "Queried facts: \"" << remove_query << "\" are extensional" ;
          cout << " facts, action not performed on these facts!" << endl ;
        }
      if(target == EMPTY)
        return true ;

      // because we could have multiple queries, we will need to
      // perform a clean up of fact_db at the beginning
      facts.erase_intensional_facts() ;
      // then we need to copy the fact_db
      // This is because we want to only put the queried facts
      // back into the global fact_db
      fact_db local_facts(facts) ;
      sched_db local_scheds ;

      /*
#ifdef USE_PAPI

      if((perr=PAPI_create_eventset(&ev_set)))
        cout<<"\nPAPAI_create_evebtset failed."<<PAPI_strerror(perr)<<"\n";


      if((retval= PAPI_multiplex_init())<PAPI_OK)
        cout<<"\nEvent set multiplexing initialization error\n";

      retval=PAPI_set_multiplex(ev_set);
      if((retval==PAPI_EINVAL) &&(PAPI_get_multiplex(ev_set)>0))
        cout<<"This event set already hs multiplexing enabled";
      else if(retval !=PAPI_OK)  cout<<"\nSet multiplexing error\n";
      else cout<<"\nsuccess\n";

      retval=PAPI_get_multiplex(ev_set);
      if(retval>0) cout<<"This event set is ready for multiplexing";
      if(retval==0)cout<<"This venet set is not enabled for multip0lexig";
      if(retval<0) cout<<"\nerror\n";

      if((perr=PAPI_add_event(ev_set,PAPI_L1_DCH)))
        cout<<__LINE__<<"PAPI_add_event failed."<<PAPI_strerror(perr)<<"\n";


      if((perr=PAPI_add_event(ev_set,PAPI_L1_ICH)))
        cout<<__LINE__<<"PAPI_add_event failed."<<PAPI_strerror(perr)<<"\n";

      if((perr=PAPI_add_event(ev_set,PAPI_L2_DCM)))
        cout<<__LINE__<<"PAPI_add_event failed."<<PAPI_strerror(perr)<<"\n";

      if((perr=PAPI_add_event(ev_set,PAPI_L2_ICM)))
        cout<<__LINE__<<"PAPI_add_event failed."<<PAPI_strerror(perr)<<"\n";


      if((perr=PAPI_list_events(ev_set,evlist,&ncnt)))
        cout<<__LINE__<<"PAPI_list_events failed."<<PAPI_strerror(perr)<<"\n";


      cout<<"\n number of events"<<ncnt<<"\n";
      for(i=0;i<ncnt;i++)
        if((perr=PAPI_event_code_to_name(evlist[i],evname[i])) == PAPI_ENOTPRESET)
          {}
        else if(perr!=PAPI_OK)
          cout<<__LINE__<<" Naming event failed."<<PAPI_strerror(perr)<<"[i="<<i<<" event="<<evlist[i]<<"\n";
      if((perr=PAPI_start(ev_set)))
        cout<<"\nPAPI_start_event failed."<<PAPI_strerror(perr)<<"\n";
#endif
      
      */

      sw.start() ;
      executeP schedule =
        create_execution_schedule(rdb,local_facts,local_scheds,target) ;
      if(schedule == 0)
        throw StringError("makeQuery: query failed!") ;


      double tlocal = sw.stop() ;
      double tglobal = 0 ;
      MPI_Allreduce(&tlocal,&tglobal, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD) ;
      debugout << "time to create schedule " << tglobal  << endl ;
      
      // If a schedule was generated, execute it
      if(MPI_rank == 0)
        cout << "begin execution" << endl ;

      if(MPI_rank == 0) {
	if (threading_pointwise)
	  cout << "--threading " << num_threaded_pointwise
	       << "/" << num_total_pointwise << " pointwise rules" << endl;
	if (threading_global_reduction)
	  cout << "--threading " << num_threaded_global_reduction
	       << "/" << num_total_global_reduction 
	       << " global reduction rules" << endl;
	if (threading_local_reduction)
	  cout << "--threading " << num_threaded_local_reduction
	       << "/" << num_total_local_reduction 
	       << " local reduction rules" << endl;
	if (threading_chomping)
	  cout << "--threading " << num_threaded_chomping
	       << "/" << num_total_chomping << " chomping rules" << endl;
	if (threading_recursion)
	  cout << "--threading " << num_threaded_recursion
	       << "/" << num_total_recursion << " recursive rules" << endl;
      }
      
      if(schedule_output) {
        // Save the schedule in the file schedule for reference
        ostringstream oss ;
        oss << "debug/schedule" ;

        if(MPI_processes > 1) {
          oss << "-" << MPI_rank ;
        }

        string sched_filename = oss.str() ;
        ofstream sched_file(sched_filename.c_str(),ios::out) ;
        schedule->Print(sched_file) ;
        sched_file.close() ;
      }

      // execute schedule
      stopWatch sw ;
      sw.start() ;
      // setting this external pointer
      exec_current_fact_db = &local_facts ;
      schedule->execute(local_facts, local_scheds) ;
#ifdef PARTICLE_PERF
      pwalltime = sw_particle.stop() ;
#endif
      double exec_time = sw.stop() ;

      if(collect_memory_info) {
        debugout << "collect memory info:" << endl ;
        collectMemory memProf ;
        schedule->dataCollate(memProf) ;
        memProf.PrintSummary(debugout) ;
      }
      
      collectTiming timeProf ;
      schedule->dataCollate(timeProf) ;

      timeProf.balanceAnalysis(MPI_COMM_WORLD) ;

      double compute_time_local = timeProf.getComputeTime() ;
      double prof_exec_time = timeProf.getTotalTime() ;
      double compute_time_total = 0 ;
      MPI_Allreduce(&compute_time_local,&compute_time_total, 1, MPI_DOUBLE,
                    MPI_SUM,MPI_COMM_WORLD) ;


      timeProf.PrintSummary(debugout) ;


      /*
#ifdef USE_PAPI
      if((perr=PAPI_read(ev_set,counts)))
        cout<<"PAPI_read failed."<<PAPI_strerror(perr)<<"\n";

      cout<<"Counts registered\n";
      for(i=0;i<ncnt;i++)
        cout<<evname[i]<<"="<<counts[i]<<"\n";
#endif
      */

      Loci::debugout << "Time taken for execution of the schedule = " << exec_time << " seconds " << endl ;

      // put the computed results back to the global facts
      // but we want to restore the facts back to its global
      // numbering scheme if the fact_db was started
      // distributed at the beginning, since if it was
      // started distributed at the beginning, then we've already
      // done the local renumbering step to the facts.
      if(local_facts.is_distributed_start()) {
        fact_db::distribute_infoP df = local_facts.get_distribute_info() ;
        // first get the local to global dMap
        dMap l2g ;
        entitySet dom = df->l2g.domain() ;
        for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei)
          l2g[*ei] = df->l2g[*ei] ;

        for(variableSet::const_iterator vi=target.begin();
            vi!=target.end();++vi) {
          storeRepP srp = local_facts.get_variable(*vi) ;
          // the results are clearly intensional facts
          facts.create_intensional_fact(*vi,srp->remap(l2g)) ;
        }
      }else{
        for(variableSet::const_iterator vi=target.begin();
            vi!=target.end();++vi) {
          storeRepP srp = local_facts.get_variable(*vi) ;
          facts.create_intensional_fact(*vi,srp) ;
        }
      }
	  
		
      if(profile_memory_usage) {
        Loci::debugout << "++++++++Memory Profiling Report++++++++"
                       << endl ;
        Loci::debugout << "Peak Memory used: " << LociAppPeakMemory
                       << " bytes ("
                       << LociAppPeakMemory/(1024*1024)
                       << "MB)" << endl ;

        Loci::debugout << "Total allocation requests: "
                       << LociAppAllocRequestBeanCounting
                       << " bytes ("
                       << LociAppAllocRequestBeanCounting/(1024*1024)
                       << "MB)" << endl ;

        Loci::debugout << "Total recycle requests: "
                       << LociAppFreeRequestBeanCounting
                       << " bytes ("
                       << LociAppFreeRequestBeanCounting/(1024*1024)
                       << "MB)" << endl ;

        Loci::debugout << "Peak Memory in bean counting: "
                       << LociAppPeakMemoryBeanCounting << " bytes ("
                       << LociAppPeakMemoryBeanCounting/(1024*1024)
                       << "MB)" << endl ;
        Loci::debugout << "The largest allocation was: "
                       << LociAppLargestAlloc << " bytes ("
                       << LociAppLargestAlloc/(1024*1024) << "MB)"
                       << " for variable: " << LociAppLargestAllocVar
                       << endl ;
        Loci::debugout << "The largest recycle was: "
                       << LociAppLargestFree << " bytes ("
                       << LociAppLargestFree/(1024*1024) << "MB)"
                       << " for variable: " << LociAppLargestFreeVar
                       << endl ;
        Loci::debugout << "All input variables size: "
                       << LociInputVarsSize << " bytes ("
                       << LociInputVarsSize/(1024*1024) << "MB)"
                       << endl ;
        if(MPI_processes > 1) {
          // code to find out the largest memory bounds on all processes
          double LargestPeakMemory = 0 ;
          MPI_Allreduce(&LociAppPeakMemory,
                        &LargestPeakMemory,
                        1, MPI_DOUBLE,
                        MPI_MAX,MPI_COMM_WORLD) ;

          double totalPeakMemory = 0 ;
          MPI_Allreduce(&LociAppPeakMemory,
                        &totalPeakMemory,
                        1, MPI_DOUBLE,
                        MPI_SUM,MPI_COMM_WORLD) ;
          double avgPeakMemory = totalPeakMemory/MPI_processes ;

          double LargestPeakMemoryBeanCounting = 0 ;
          MPI_Allreduce(&LociAppPeakMemoryBeanCounting,
                        &LargestPeakMemoryBeanCounting,
                        1, MPI_DOUBLE,
                        MPI_MAX,MPI_COMM_WORLD) ;

          double totalPeakMemoryBeanCounting = 0 ;
          MPI_Allreduce(&LociAppPeakMemoryBeanCounting,
                        &totalPeakMemoryBeanCounting,
                        1, MPI_DOUBLE,
                        MPI_SUM,MPI_COMM_WORLD) ;
          double avgPeakMemoryBeanCounting = totalPeakMemoryBeanCounting/MPI_processes ;

          Loci::debugout << endl ;
          Loci::debugout << "The global largest Peak Memory: "
                         << LargestPeakMemory << " bytes ("
                         << LargestPeakMemory/(1024*1024) << "MB)"
                         << endl ;
          Loci::debugout << "The global mean Peak Memory: "
                         << avgPeakMemory << " bytes ("
                         << avgPeakMemory/(1024*1024) << "MB)"
                         << ", imbalance="
                         << 100*(LargestPeakMemory-avgPeakMemory)/avgPeakMemory
                         << "%" << endl ;
          Loci::debugout << "The global largest Peak Memory in bean counting: "
                         << LargestPeakMemoryBeanCounting << " bytes ("
                         << LargestPeakMemoryBeanCounting/(1024*1024) << "MB)"
                         << endl ;
          Loci::debugout << "The global mean Peak Memory in bean counting: "
                         << avgPeakMemoryBeanCounting << " bytes ("
                         << avgPeakMemoryBeanCounting/(1024*1024) << "MB)"
                         << ", imbalance="
                         << 100*(LargestPeakMemoryBeanCounting-avgPeakMemoryBeanCounting)/avgPeakMemoryBeanCounting
                         << "%" << endl ;
          // Reset flags for next run
          LociAppPeakMemory = 0 ;
          LociAppAllocRequestBeanCounting = 0 ;
          LociAppFreeRequestBeanCounting = 0 ;
          LociAppPeakMemoryBeanCounting = 0 ;
          LociAppLargestAlloc = 0 ;
          LociAppLargestAllocVar= variable("EMPTY") ;
          LociAppLargestFree = 0 ;
          LociAppLargestFreeVar= variable("EMPTY") ;
          LociAppPMTemp = 0 ;
          LociInputVarsSize = 0 ;
        }
      }

#ifdef DYNAMIC_TIMING
      {
        double mytime = ta_expand.getTime() ;
        double maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total expand time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_expand2.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total expand time (2): "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_expand_start.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total expand start block time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_expand_cache.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total expand cache management time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_expand_collect_img.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total expand image collection time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_expand_missing.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total expand missing domain comp time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_expand_comm.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total expand comm time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_context.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total context time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_context_nonepd.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total context non-expand block time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_context_nonepd_domt.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total context non-epd block (dom) time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_context_nonepd_ints.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total context non-epd block (int) time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_context_epdend.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total context expand-end block time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_context_epdmid.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total context expand-mid block time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_context_epdsta.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total context expand-sta block time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_output_oh.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total output overhead: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_compute.getTime() ;
        maxtime = 0 ;        
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total compute time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_record_erase.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total record erase time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_erase.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total erase time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_invalidate.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total invalidate time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_keyremoval.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total key removal time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_insertion.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total insertion time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_keyinsert.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total key insert time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_dist.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total distribution time: "
                       << mytime << ", max = " << maxtime << endl ;
        Loci::debugout << "[dynamic] total distribution number: "
                       << ta_dist_number << endl ;

        mytime = ta_key_dist.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total key distribution time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_dist_renumber.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total distribution renumber time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_push.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total push time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_pw_push.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total pointwise push time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_param_push.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total param push time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_param_pack.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total param pack time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_param_unpack.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total param unpack time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_param_reduce.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total param reduce time: "
                       << mytime << ", max = " << maxtime << endl ;
        Loci::debugout << "[dynamic] total param reduction num: "
                       << ta_param_reduce_num << endl ;

        mytime = ta_dist_keys.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total dist keys time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_dist_all2all.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total dist all2all time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_dist_pack.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total dist pack time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_dist_unpack.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total dist unpack time: "
                       << mytime << ", max = " << maxtime << endl ;

        mytime = ta_dctrl.getTime() ;
        maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "[dynamic] total dcontrol reset time: "
                       << mytime << ", max = " << maxtime << endl ;

        Loci::debugout << "[dynamic] total drule executing num: "
                       << ta_drule_executes << endl ;

        Loci::debugout << "[dynamic] total dCTRL executing num: "
                       << ta_dctrl_executes << endl ;

      }
#endif
#ifdef PARTICLE_PERF
      {
        double pmaxtime = 0 ;
        MPI_Allreduce(&pwalltime,&pmaxtime,1,
                      MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "Particle Walltime: " << pwalltime
                       << ", max = " << pmaxtime << endl ;
      }
#endif
      // communicate the execution time
      if(MPI_processes > 1) {
        double mytime = exec_time ;
        double maxtime = 0 ;
        MPI_Allreduce(&mytime,&maxtime,1,
                      MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        Loci::debugout << "Global max time taken for execution"
                       << " of the schedule = "
                       << maxtime << " seconds " << endl ;
        mytime = prof_exec_time ;
        MPI_Allreduce(&mytime,&maxtime,1,
                      MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
        if(MPI_rank == 0 ) {
          double eff = compute_time_total/(double(MPI_processes)*maxtime) ;
          cout << "Schedule execution complete, estimated parallel efficiency = "
               << ceil(1000.0*eff)/10.0 << "%." << endl ;
        }
      }
    } catch(const BasicException &err) {
      cerr << "Loci found an error during MakeQuery" << endl ;
      err.Print(cerr) ;
      Loci::Abort() ;
    } catch (const bad_alloc &x) {
      cerr << "Out of memory: " << x.what() << endl ;
      Loci::Abort() ;
    } catch(...) {
      cerr << "Unknown Exception Caught" << endl ;
      Loci::Abort() ;
    }

    debugout << "Time to execute query for '" << query << "' is " << sw.stop()
             << endl ;

    return true ;

  }

} // end of namespace Loci

