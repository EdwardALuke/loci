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
#ifndef COMP_TOOLS_H
#define COMP_TOOLS_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <fact_db.h>
#include <execute.h>
#include <Tools/digraph.h>

#include "sched_tools.h"
#include <vector>
#include <map>
#include <deque>
#include <gpurep.h>
#include <gpuMap.h>

#include <mpi.h>
using std::vector;

namespace Loci {

  bool rule_has_mapping_in_output(rule r);
  variableSet input_variables_with_mapping(rule r);
  variableSet input_variables(rule r);
  bool is_intensive_rule_output_mapping(rule my_rule, const fact_db &facts);
  bool process_policy_duplication(variable v, sched_db &scheds, fact_db &facts);
  void set_duplication_of_variables(variableSet vlst, sched_db &scheds, fact_db &facts);
  entitySet vmap_source_exist(const vmap_info &vmi, fact_db &facts, sched_db &scheds) ;
  entitySet vmap_target_exist(const vmap_info &vmi, fact_db &facts,
                              entitySet compute, sched_db &scheds) ;

  entitySet vmap_source_exist_apply(const vmap_info &vmi, fact_db &facts,
                                    variable reduce_var, sched_db &scheds) ;
 
  
  void existential_rule_analysis(rule f, fact_db &facts, sched_db &scheds) ;
  entitySet process_rule_requests(rule f, fact_db &facts, sched_db &scheds) ;
  
  void existential_applyrule_analysis(rule apply, fact_db &facts, sched_db &scheds) ;
  entitySet process_applyrule_requests(rule apply, rule unit_tag, bool &output_mapping,fact_db &facts, sched_db &scheds) ;
  
  void existential_blackboxrule_analysis
  (rule f, fact_db &facts, sched_db &scheds) ;

  entitySet process_blackboxrule_requests
  (rule f, fact_db &facts, sched_db &scheds) ;

  std::vector<std::pair<variable,entitySet> >
    barrier_existential_rule_analysis(variableSet vlst, fact_db &facts, sched_db &scheds) ;
  std::vector<std::pair<variable,entitySet> >
    send_ent_for_plist(variableSet vlst, fact_db &facts, sched_db &scheds);
  std::list<comm_info>
  barrier_process_rule_requests(variableSet vars, fact_db &facts, sched_db &scheds) ;

  entitySet send_requests(const entitySet& e, variable v, fact_db &facts,
                          std::list<comm_info> &clist) ;
  
  std::list<comm_info>
  put_precomm_info(std::vector<std::pair<variable,entitySet> > send_entities,
                   fact_db &facts) ;
  
  std::list<comm_info> sort_comm(std::list<comm_info> slist, fact_db &facts) ;
  
  void create_user_function(unsigned char* , unsigned char* , int*,
                            MPI_Datatype* ) ;
  
  typedef std::map<variable,entitySet> vdefmap ;
  entitySet vmap_target_requests(const vmap_info &vmi, const vdefmap &tvarmap,
                                 fact_db &facts, sched_db &scheds, 
				 bool is_request_modification_allowed=true) ;
  entitySet vmap_source_requests(const vmap_info &vmi, fact_db &facts,
                                 entitySet compute, sched_db &scheds) ;

  std::vector<digraph::vertexSet> schedule_dag(const digraph &g,
                           digraph::vertexSet start_vertices = EMPTY,
                           digraph::vertexSet only_vertices =
                                               interval(UNIVERSE_MIN,UNIVERSE_MAX)) ;

  class visitor ;
  
  class loop_compiler : public rule_compiler {
  public:
    int cid ; // id number of this compiler
    std::vector<rule_compilerP> collapse_comp ; 
    variable cond_var ;
    std::vector<rule_compilerP> advance_comp ;
    variableSet advance_vars ;
    time_ident tlevel ;
    variableSet all_loop_vars ;
    ////////////
    variableSet collapse_vars ; // the collapse variables
    std::vector<digraph::vertexSet> collapse_sched ;
    std::vector<digraph::vertexSet> advance_sched ;
    // the rotate list
    std::list<std::list<variable> > rotate_lists ;
    ////////////////
    // the internal loop graph and the rulecomp map
    digraph loop_gr ;
    digraph collapse_gr ;
    digraph advance_gr ;
    rulecomp_map rule_compiler_map ;
    ////////////////

    //To add communication manually for advance variable at the end of advance_comp
    // std::list<comm_info> advance_variables_barrier; 
	
  public:
    loop_compiler(rulecomp_map &rp, digraph gin, int id) ;
    virtual void accept(visitor& v) ;
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  // rule compiler for rule with concrete implementation
  class impl_compiler : public rule_compiler {
    rule impl ;  // rule to implement
    // existential analysis info
    //entitySet exec_seq ;
  public:
    impl_compiler(rule r)  { impl=r;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
#ifdef DYNAMICSCHEDULING
  // rule compiler for rule intended to be evaluated at runtime.
  // example rules are those that may reside in a totally dynamic
  // keyspace
  class dynamic_impl_compiler : public rule_compiler {
    rule impl ;  // rule to implement
    // keyspace that it exists
    KeySpaceP space ;
    // all the rule's sources that are changing
    variableSet volatile_sources ;
    // and the sources that are static
    variableSet static_sources ;
  public:
    dynamic_impl_compiler(rule r, KeySpaceP kp,
                          const variableSet& vs,
                          const variableSet& ss)  {
      impl=r ;
      space = kp ;
      volatile_sources = vs ;
      static_sources = ss ;
    }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
#endif
  
  // rule compiler for single rule recursion
  class impl_recurse_compiler : public rule_compiler {
    rule impl ;  // rule to implement
    struct fcontrol {
      std::list<entitySet> control_list ;
      std::map<variable,entitySet> generated ;
      bool use_constraints ;
      entitySet nr_sources, constraints ;
      struct mapping_info {
        entitySet total ;
        std::vector<MapRepP> mapvec ;
        std::vector<variable> mapvar ;
        variable v ;
      } ;
      std::vector<mapping_info> recursion_maps, target_maps ;
    } ;
    fcontrol  control_set ;
    sequence fastseq ;
    // std::vector<entitySet > par_schedule ;//never used

    // std::vector<std::pair<variable,entitySet> > pre_send_entities ;//assigned, never used
    //std::list<comm_info> clist ;// assigned, but never used
  public:
    int cid ; // id number of this compiler
    impl_recurse_compiler(rule r, int id):cid(id)
      { impl = r ;}
    ruleSet get_rules() const
      {ruleSet rs; rs += impl; return rs ;}
    virtual void accept(visitor& v) ;
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class recurse_compiler : public rule_compiler {
    rulecomp_map &rule_process ;
  
    rule_compilerP calc(const rule &r) 
    {return rule_process[r] ;}

    ruleSet recurse_rules ;
    variableSet recurse_vars ;
    struct fcontrol {
      std::list<entitySet> control_list ;

      std::map<variable,entitySet> generated ;
      bool use_constraints ;
      entitySet nr_sources, constraints ;
      struct mapping_info {
        entitySet total ;
        std::vector<MapRepP> mapvec ;
        variable v ;
      } ;
      std::vector<mapping_info> recursion_maps, target_maps ;
    } ;
    std::map<rule,fcontrol > control_set ;
    std::list<std::vector<std::pair<variable,entitySet> > > recurse_send_entities ;
    std::map<variable,std::vector<std::list<comm_info> > > send_req_var ;
  public:
    int cid ; //id number of this compiler
    recurse_compiler(rulecomp_map &rp, ruleSet rs, int id) : rule_process(rp),cid(id)
    {
      recurse_rules = rs ;
      for(ruleSet::const_iterator ri=rs.begin();ri!=rs.end();++ri) 
        recurse_vars += ri->targets() ;
    }
    ruleSet get_rules() const {return recurse_rules ;}
    virtual void accept(visitor& v) ;
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class dag_compiler : public rule_compiler {
  public:
    int cid ; // id number of this compiler
    std::vector<rule_compilerP> dag_comp ;
    ////////////
    std::vector<digraph::vertexSet> dag_sched ;
    /////////////////
    digraph dag_gr ;
    rulecomp_map rule_compiler_map ;
    /////////////////
  public:
    dag_compiler(rulecomp_map &rp, digraph dag, int id) ;
    /////////////////
    virtual void accept(visitor& v) ;
    /////////////////
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class barrier_compiler : public rule_compiler {
    variableSet barrier_vars ;
   

   //  std::vector<std::pair<variable,entitySet> > send_entities ;
//     std::list<comm_info> clist ;
//     std::list<comm_info> plist ;
  public:
    barrier_compiler(variableSet &vars)
      : barrier_vars(vars) {}
    /////////////
    virtual void accept(visitor& v) {}
    /////////////
    virtual void set_var_existence(fact_db &facts,  sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class singleton_var_compiler : public rule_compiler {
    variableSet barrier_vars ;
  public:	
    singleton_var_compiler(variableSet &vars)
      :  barrier_vars(vars) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
  class reduce_param_compiler : public rule_compiler {
    vector<variable> reduce_vars ;
    vector<rule> unit_rules ;
    vector<CPTR<joiner> > join_ops ;
  public:
    reduce_param_compiler(const vector<variable> &v, const vector<rule> &ur,
                          vector<CPTR<joiner> >&jop) :
      reduce_vars(v), unit_rules(ur), join_ops(jop) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class reduce_store_compiler : public rule_compiler {
    variable reduce_var ;
    rule unit_rule ;
    CPTR<joiner> join_op ;

    // std::list<comm_info> rlist ;  // reduction communication
    //std::list<comm_info> clist ;  // comm from owner to requester
  public:
    reduce_store_compiler(const variable &v, const rule &ur,
                          CPTR<joiner> &jop) :
      reduce_var(v), unit_rule(ur), join_op(jop) {
	std::set<vmap_info> uoutput = unit_rule.get_info().desc.targets ;
	std::set<vmap_info> uinput = unit_rule.get_info().desc.sources ;
	bool outputfail = false ;
	if(uoutput.size() != 1) {
	  outputfail = true ;
	} else if(uoutput.begin()->mapping.size()) {
	  outputfail = true ;
	}
	if(outputfail)
	  cerr << "WARNING: unit rule must have a single varaible as output"
	       << endl ;
	outputfail = false ;
	std::set<vmap_info>::const_iterator si ;
	for(si=uinput.begin();si!=uinput.end();++si) 
	  if(si->mapping.size()) 
	    outputfail = true ;
	if(outputfail) {
	  cerr << "WARNING: unit rule has invalid input, mappings not allowed!"
	       << endl ;
	}
      }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
  class conditional_compiler : public rule_compiler {
  public:
    int cid ; //id number of this compiler
    std::vector<rule_compilerP> dag_comp ;
    variable cond_var ;
    /////////
    std::vector<digraph::vertexSet> dag_sched ;
    ///////////
    digraph cond_gr ;
    rulecomp_map rule_compiler_map ;
  public:
    conditional_compiler(rulecomp_map &rp, digraph gin,
                         variable conditional, int id) ;
    virtual void accept(visitor& v) ;
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  // apply rule compiler 
  class apply_compiler : public rule_compiler {
    rule apply,unit_tag ;  // rule to applyement
    
    // existential analysis info
    //entitySet exec_seq ;
    bool output_mapping ;
  public:
    apply_compiler(rule r, rule ut)
    { apply=r; unit_tag = ut ; }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
#ifdef DYNAMICSCHEDULING
  // dynamic apply rule compiler 
  class dynamic_apply_compiler : public rule_compiler {
    rule apply, unit_tag ;  // rule to apply
    // keyspace that it exist
    KeySpaceP space ;
  public:
    dynamic_apply_compiler(rule a, rule u, KeySpaceP kp) {
      apply = a ;
      unit_tag = u ;
      space = kp ;
    }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  // compiler that invalidates the dynamic clones
  class dclone_invalidate_compiler: public rule_compiler {
    // the variable name, "var_unique" is the name removed synonym
    variable var, var_unique ;
    KeySpaceP self_clone ;
    std::vector<KeySpaceP> shadow_clone ;
  public:
    dclone_invalidate_compiler(const variable& v,
                               const variable& vu,
                               KeySpaceP sc,
                               const std::vector<KeySpaceP>& sac)
      :var(v),var_unique(vu),self_clone(sc),shadow_clone(sac) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts,
                                               sched_db &scheds) ;
  } ;

  // compiler responsible for keyspace redistribution
  class keyspace_dist_compiler: public rule_compiler {
    // keyspaces that needs distribution
    std::vector<std::string> space_names ;
  public:
    keyspace_dist_compiler(const std::vector<std::string>& sn)
      :space_names(sn) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db& facts, sched_db& scheds) ;
    virtual void process_var_requests(fact_db& facts, sched_db& scheds) ;
    virtual executeP create_execution_schedule(fact_db& facts,
                                               sched_db& scheds) ;
  } ;

  class insertion_rule_compiler: public rule_compiler {
    rule rule_tag ;
  public:
    insertion_rule_compiler(const rule& r):rule_tag(r) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db& facts, sched_db& scheds) {}
    virtual void process_var_requests(fact_db& facts, sched_db& scheds) {}
    virtual executeP create_execution_schedule(fact_db& facts,
                                               sched_db& scheds) ;
  } ;
  
  class deletion_rule_compiler: public rule_compiler {
    rule rule_tag ;
  public:
    deletion_rule_compiler(const rule& r):rule_tag(r) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db& facts, sched_db& scheds) {}
    virtual void process_var_requests(fact_db& facts, sched_db& scheds) {}
    virtual executeP create_execution_schedule(fact_db& facts,
                                               sched_db& scheds) ;
  } ;
  
  class erase_rule_compiler: public rule_compiler {
    rule rule_tag ;
  public:
    erase_rule_compiler(const rule& r):rule_tag(r) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db& facts, sched_db& scheds) {}
    virtual void process_var_requests(fact_db& facts, sched_db& scheds) {}
    virtual executeP create_execution_schedule(fact_db& facts,
                                               sched_db& scheds) ;
  } ;
  
  // compiler that resets the dynamic rule control
  class dcontrol_reset_compiler: public rule_compiler {
    // the variable name, "var_unique" is the name removed synonym
    variable var, var_unique ;
    std::set<std::string> register_keyspaces ;
  public:
    dcontrol_reset_compiler(const variable& v,const variable& vu,
                            const std::set<std::string>& rk)
      :var(v),var_unique(vu),register_keyspaces(rk) {}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts,
                                               sched_db &scheds) ;
  } ;
#endif
  
  class promote_compiler : public rule_compiler {
    rule r ;
  public:
    promote_compiler(rule rin)
    { r = rin; }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
  class generalize_compiler : public rule_compiler {
    rule r ;
  public:
    generalize_compiler(rule rin)
    { r = rin; }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
  class priority_compiler : public rule_compiler {
    rule r ;
  public:
    priority_compiler(rule rin)
    { r = rin; }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class gpu2cpu_compiler : public rule_compiler {
    rule r ;
  public:
    gpu2cpu_compiler(rule rin)
    { r = rin; }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class cpu2gpu_compiler : public rule_compiler {
    rule r ;
  public:
    cpu2gpu_compiler(rule rin)
    { r = rin; }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class map2gpu_compiler : public rule_compiler {
    rule r ;
  public:
    map2gpu_compiler(rule rin)
    { r = rin; }
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
  
  class execute_msg : public execute_modules {
    std::string msg ;
  public:
    execute_msg(std::string m) : msg(m) {}
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_msg";};
    virtual void dataCollate(collectData &data_collector) const {}
  } ;

  class execute_gpuSync: public execute_modules {
    variableSet vars ;
  public:
    execute_gpuSync(variableSet v) : vars(v) {}
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() { return "execute_gpuSync" ; }
    virtual void dataCollate(collectData &data_collector) const {}
  } ;
  
  class execute_gpu2cpu_copy: public execute_modules {
    rule r ;
    gpuRepP gpuvar ;
    storeRepP cpuvar ;
    entitySet copyset ;
  public:
    execute_gpu2cpu_copy(rule rin, gpuRepP vgpu, storeRepP vcpu,entitySet s) :
      r(rin),gpuvar(vgpu),cpuvar(vcpu),copyset(s) {}
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_gpu2cpu";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  class execute_cpu2gpu_copy: public execute_modules {
    rule r ;
    gpuRepP gpuvar ;
    storeRepP cpuvar ;
    entitySet copyset ;
  public:
    execute_cpu2gpu_copy(rule rin, gpuRepP vgpu, storeRepP vcpu, entitySet s) :
      r(rin),gpuvar(vgpu),cpuvar(vcpu),copyset(s) {}
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_cpu2gpu";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  class execute_map2gpu_copy: public execute_modules {
    rule r ;
    gpuMapRepP gpuvar ;
    storeRepP cpuvar ;
    entitySet copyset ;
  public:
    execute_map2gpu_copy(rule rin, gpuMapRepP vgpu, storeRepP vcpu, entitySet s) :
      r(rin),gpuvar(vgpu),cpuvar(vcpu),copyset(s) {}
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_map2gpu";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  class execute_comm : public execute_modules {
    std::vector<std::pair<int,std::vector<send_var_info> > > send_info ;
    std::vector<std::pair<int,std::vector<recv_var_info> > > recv_info ;
    std::vector<std::vector<storeRepP> > send_vars ;
    std::vector<std::vector<storeRepP> > recv_vars ;
    int *maxr_size, *maxs_size, *r_size, *s_size, *recv_sizes ;
    unsigned char **recv_ptr , **send_ptr ;
    MPI_Request *request;
    MPI_Status *status ;
    timeAccumulator timer ;
  public:
    execute_comm(std::list<comm_info> &plist, fact_db &facts) ;
    ~execute_comm() ;
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_comm";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  class execute_comm2: public execute_modules {
    struct send_unit {
      variable v ;
      storeRepP rep ;
      entitySet send ;
    } ;
    struct recv_unit {
      variable v ;
      storeRepP rep ;
      sequence recv ;
    } ;
    struct send_proc {
      int proc ;
      int send_size ;           // size of data to be sent to proc
      int max_send_size ;       // maximum recorded buffer size so far
      unsigned char* buf ;      // pointer to send buffer
      std::vector<send_unit> units ;
    } ;
    struct recv_proc {
      int proc ;
      // note, this records the maximum buffer size received
      // so far. the actual message size is extracted from
      // the status objects
      int recv_size ;
      unsigned char* buf ;      // pointer to recv buffer
      std::vector<recv_unit> units ;
    } ;

    std::vector<send_proc> send_info ;
    std::vector<recv_proc> recv_info ;

    // this is the base MPI tag number for all execute_comm2 objects
    static int tag_base ;
    // MPI tags that will be used in a particular comm object
    int tag1, tag2 ;
    
    timeAccumulator timer ;
  public:
    execute_comm2(std::list<comm_info>& plist, fact_db& facts) ;
    ~execute_comm2(){}
    virtual void execute(fact_db& facts, sched_db& scheds) ;
    virtual void Print(std::ostream& s) const ;
    virtual string getName() {return "execute_comm2";};
    virtual void dataCollate(collectData& data_collector) const ;
    static void inc_comm_step() {
      int lt=tag_base, gt=0 ;
      MPI_Allreduce(&lt, &gt, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;
      tag_base = gt ;
      if(tag_base > 32500) {
	tag_base=1500 ; // recycle tags
      }
    }
  } ;
  
  class allocate_var_compiler : public rule_compiler {
    variableSet allocate_vars ;
  public:
    allocate_var_compiler(const variableSet& vars)
      : allocate_vars(vars) {}
    allocate_var_compiler(const variable& var)
    {allocate_vars += var ;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class free_var_compiler : public rule_compiler {
    variableSet free_vars ;
  public:
    free_var_compiler(const variableSet& vars)
      : free_vars(vars) {}
    free_var_compiler(const variable& var)
    {free_vars += var ;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class execute_memProfileAlloc : public execute_modules {
    variableSet vars ;
  public:
    execute_memProfileAlloc(const variableSet& vars): vars(vars) {}
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_memProfileAlloc";};	
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  class execute_memProfileFree : public execute_modules {
    variableSet vars ;
  public:
    execute_memProfileFree(const variableSet& vars) : vars(vars) {}
    virtual void execute(fact_db &facts, sched_db& scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_memProfileFree";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  class memProfileAlloc_compiler : public rule_compiler {
    variableSet vars ;
 public:
    memProfileAlloc_compiler(const variableSet& vars): vars(vars) {}
    memProfileAlloc_compiler(const variable& var)
    {vars += var ;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) {}
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) {}
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class memProfileFree_compiler : public rule_compiler {
    variableSet vars ;
 public:
    memProfileFree_compiler(const variableSet& vars): vars(vars) {}
    memProfileFree_compiler(const variable& var)
    {vars += var ;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) {}
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) {}
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  class chomp_compiler: public rule_compiler {
  public:
    digraph chomp_graph ;
    variableSet chomp_vars ;
    std::vector<digraph::vertexSet> chomp_sched ;
    std::vector<std::pair<rule,rule_compilerP> > chomp_comp ;
    std::deque<entitySet> rule_seq ;
    std::map<rule,rule> apply2unit ;
    std::vector<variableSet> barrier_sets;
    std::vector<variableSet> old_barrier_sets;//generated by visitor
    std::vector<std::pair<rule,rule_compilerP> > old_chomp_comp ;//generated by visitor

    
  public:
    chomp_compiler(const digraph& cgraph,const variableSet& cvars,
                   const std::map<rule,rule>& a2u) ;
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db& facts,sched_db& scheds) ;
    virtual void process_var_requests(fact_db& facts,sched_db& scheds) ;
    virtual executeP create_execution_schedule(fact_db& facts,sched_db& scheds) ;
  } ;

  // compiler to handle the dynamic scheduling stuff
  class dynamic_compiler: public rule_compiler {
  public:
    int cid ;
    std::vector<rule_compilerP> comp ;
    std::vector<digraph::vertexSet> sched ;
    digraph gr ;
    rulecomp_map rule_compiler_map ;
    variableSet all_reduce_vars ;
    std::map<variable,std::pair<rule,CPTR<joiner> > > reduce_info ;
  public:
    dynamic_compiler(rulecomp_map& rp, const digraph& g, int id) ;
    virtual void accept(visitor& v) {}
    void collect_reduce_info() ;
    void schedule() ;
    void compile() ;
    virtual void set_var_existence(fact_db& facts, sched_db& scheds) ;
    virtual void process_var_requests(fact_db& facts, sched_db& scheds) ;
    virtual executeP create_execution_schedule(fact_db& facts,
                                               sched_db& scheds) ;
  } ;
  
  // rule compiler for constraint rules
  class constraint_compiler : public rule_compiler {
    rule constraint_rule ;  // the constraint rule
    // entitySet exec_seq ;
 public:
    constraint_compiler(rule r)  { constraint_rule=r;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db& facts,
                                               sched_db& scheds) ;
  } ;

  // map rule compiler for rule that computes maps
  class map_compiler : public rule_compiler {
    rule map_impl ;  // rule that computes a map
    // existential analysis info
    // entitySet exec_seq ;
  public:
    map_compiler(rule r)  { map_impl=r;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP create_execution_schedule(fact_db &facts,
                                               sched_db &scheds) ;
  } ;

  // rule compiler for blackbox rule
  class blackbox_compiler : public rule_compiler {
    rule impl ;  // rule to implement
    //entitySet exec_seq ;
  public:
    blackbox_compiler(rule r)  { impl=r;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP
    create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;

  // rule compiler for blackbox rule
  class superRule_compiler : public rule_compiler {
    rule impl ;  // rule to implement
  public:
    superRule_compiler(rule r)  { impl=r;}
    virtual void accept(visitor& v) {}
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) ;
    virtual executeP
    create_execution_schedule(fact_db &facts, sched_db &scheds) ;
  } ;
}

#endif
