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
#ifndef SCHED_TOOLS_H
#define SCHED_TOOLS_H

#include <map>
#include <vector>
#include <set>
#include <list>
#include <deque>

#include <Tools/cptr.h>
#include <scheduler.h>
#include <Tools/digraph.h>
#include <fact_db.h>
#include <sched_db.h>
#include <execute.h>
#include <depend_graph.h>
#include <Map.h>

#ifdef PROFILE_CODE
#include <time.h>
#endif
#include "sched_mlg.h"
#include "Loci_types.h"
using std::vector;

#ifdef PAPI_DEBUG
#include <papi.h>
#else
#define long_long long
#endif

namespace Loci {
  void extract_rule_sequence(std::vector<rule> &rule_seq,
                             const std::vector<digraph::vertexSet> &v) ;
  void set_var_types(fact_db &facts, const digraph &dg, sched_db &scheds) ;
  rule make_rename_rule(variable new_name, variable old_name) ;

  class execute_rule : public execute_modules {
  protected:
    rule_implP rp ;
    rule rule_tag ; 
    sequence exec_seq ;
    size_t exec_size ;
    timeAccumulator timer ;
#ifdef PAPI_DEBUG
    int papi_events[2];
    long_long papi_values[2];
    long_long l1_dcm;
    long_long l2_dcm;
#endif
  public:
    execute_rule(rule fi, sequence seq, fact_db &facts, const sched_db &scheds);
    execute_rule(rule fi, sequence seq, fact_db &facts, variable v, const storeRepP &p, const sched_db &scheds);
    // this method executes the prelude (if any), the kernel,
    // and the postlude (if any) of a rule
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    // this method executes the computation kernel of a rule
    virtual void execute_kernel(const sequence&);
    // this method executes the prelude 
    virtual void execute_prelude(const sequence&);
    // this one executes the postlude
    virtual void execute_postlude(const sequence&);
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_rule";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

#ifdef DYNAMICSCHEDULING
  struct ExpandStartUnit {
    variable var ;              // name of the starting var
    storeRepP rep ;             // NOTE we cannot use MapRepP here
                                // because the rep here is actually
                                // a store_ref that contains a
                                // *possible* changing underlying rep
    ExpandStartUnit(const variable& v, storeRepP r):var(v),rep(r) {}
  } ;
  std::ostream& operator<<(std::ostream& s, const ExpandStartUnit& u) ;
  std::ostream& operator<<(std::ostream& s,
                           const std::vector<ExpandStartUnit>& vu) ;
  
  // a type for recording the chains of input
  // first is the source and the second is the destination
  struct ExpandUnit {
    variable var ;              // name of the expanding var
    storeRepP src, dst ;        // same reason, must use storeRepP here
    size_t dst_reserve_size ;
    ExpandUnit(const variable& v, storeRepP s, storeRepP d, size_t r)
      :var(v),src(s),dst(d),dst_reserve_size(r) {}
  } ;

  struct ExpandBlock {
    const std::vector<entitySet>* src_ptn ;
    MPI_Comm src_comm ;
    const Map* src_pack ;
    const dMap* src_unpack ;
    const Map* dst_pack ;
    const dMap* dst_unpack ;
    std::vector<ExpandUnit> units ;
  } ;

  std::ostream& operator<<(std::ostream& s, const ExpandUnit& u) ;
  std::ostream& operator<<(std::ostream& s, const ExpandBlock& b) ;

  struct ExpandChain {
    std::string tag ;
    std::vector<ExpandStartUnit> expand_start ;
    std::vector<ExpandBlock> expand_chain ;
    ExpandBlock expand_end ;
    ExpandChain():tag("") {}
  } ;
  std::ostream& operator<<(std::ostream& s, const ExpandChain& c) ;

  // this function collects the expand chain from
  // a vmap_info pointer to a chain of mappings.
  // if an expand chain is collected, it is then
  // appended to the supplied chains. it also
  // initializes corresponding stores in the rule

  // vmi_tag is to be "source", "constraint", or "target"

  // the last 4 variables are set when the supplied info.
  // is for the target chains, they are used to indicate
  // whether the target var is crossing keyspace, the
  // key partition for the target variable, and the
  // process id for the target variables, and the target
  // keyspace's mpi communicator. if any one is0, then all
  // of them won't be set at all

  // this function returns true if a mapping chain is collected
  // and pushed back to the supplied "chains", returns
  // false otherwise (i.e., no chain is collected)
  bool 
  collect_expand_chain(const vmap_info& vmi,
                       const std::string& vmi_tag,
                       rule_implP rp, KeySpaceP space,
                       fact_db& facts, sched_db& scheds,
                       std::vector<ExpandChain>& chains,
                       bool* output_cross_space,
                       std::string* other_space_name,
                       const std::vector<entitySet>** output_ptn,
                       int* target_process_rank,
                       MPI_Comm* target_comm) ;

  struct NonExpandUnit {
    variable var ;              // name of the variable
    storeRepP var_rep ;         // rep iin the fact_db
    NonExpandUnit(const variable& v, storeRepP r):var(v),var_rep(r) {}
  } ;
  std::ostream& operator<<(std::ostream& s, const NonExpandUnit& u) ;
  std::ostream& operator<<(std::ostream& s,
                           const std::vector<NonExpandUnit>& vu) ;

  // this is a function used to collect non mapping chains,
  // for example, if "A" appears alone in a rule's signature,
  // then it also forms an input chain by itself without any
  // mapping associated, in this case, the chain will be
  // collected by this function. in addition, this function also
  // filters parameters that are not in the local keyspace, i.e.,
  // they don't need to be considered other than setting their
  // proper fact_db storeReps. returns "true" if a chain is
  // collected, "false" is not.
  bool
  collect_nonexpand_chain(const vmap_info& vmi,
                          rule_implP rp, KeySpaceP space,
                          fact_db& facts, sched_db& scheds,
                          std::vector<NonExpandUnit>& chain) ;
                       
  // this function actually expands a chain
  // returns the final image of the chain
  // for input chain, the final block is not
  // considered in the image.
  entitySet
  execute_expand_chain(ExpandChain& chain, entitySet* start) ;

  // this function will return the proper context of a dynamic rule
  // based on its inputs
  entitySet
  dynamic_rule_context(const rule& rule_tag, KeySpaceP space,
                       fact_db& facts, sched_db& scheds,
                       const std::vector<ExpandChain>& input_chains,
                       const std::vector<NonExpandUnit>& input_nes,
                       const std::vector<NonExpandUnit>& input_nec) ;

  // this is a base class for all the dynamic rule execute modules
  class execute_drule_module: public execute_modules {
  public:
    execute_drule_module() {}
    virtual ~execute_drule_module() {}
    virtual void reset_ctrl() = 0 ;
  } ;
  typedef CPTR<execute_drule_module> execute_druleP ;

  class execute_dynamic_rule : public execute_drule_module {
  protected:
    rule_implP rp ;
    rule rule_tag ; 
    timeAccumulator timer ;

    KeySpaceP space ;

    std::vector<ExpandChain> input_chains ;
    std::vector<ExpandChain> output_chains ;

    // this chain records the single non expanding vars
    // e.g.,  A <- B->C->D, E->F, G, then G will appear
    // in the chain

    // single source variables, with non local space variables removed
    std::vector<NonExpandUnit> input_ne_source ;
    // single constraint variables, with non local space variables removed
    std::vector<NonExpandUnit> input_ne_constraint ;

    // small structure used to assist the output target processing
    // for example, we have an output chain A->B->(C,D),
    // then expand is used to indicate whether this chain needs
    // expansion, in this case "yes". if no mapping is involved,
    // then "no". "targets" will be the final block of variables
    // involved, in this case "C,D". "expand chain index will the
    // the index to "output_chains" for expansion purpose.
    // "target_alloc" will be those targets that need to be allocated.
    // "target_exist" is the domain that the targets should be
    // allocated over. "send, recv" is the point to point communication
    // structure if output mapping is involved. Only target_exist and
    // send, recv are changing everytime we execute the rule,
    // all others stay the same.
    struct OutInfo {
      bool expand ;
      size_t expand_index ;
      variableSet targets ;
      variableSet targets_alloc ;
      entitySet target_exist ;
      std::vector<P2pCommInfo> send, recv ;
    } ;
    
    std::vector<OutInfo> output_info ;

    // naturally, these will need to be put into the
    // OutInfo structure for each set of target variables,
    // But since we don't allow pointwise rule to cross
    // keyspace, then every one of them are essentially
    // just the same. we therefore use just one set of
    // these parameters and share them for all output_info[i]
    bool output_cross_space ;
    std::string other_space_name ;
    // process entity partition for the target var
    const std::vector<entitySet>* output_ptn ;
    // the process rank in the target keyspace
    int target_process_rank ;
    MPI_Comm target_comm ;

    size_t large_context_size ;

    // dynamic control flag
    bool dflag ;
    entitySet context ;
  public:
    execute_dynamic_rule(rule r, KeySpaceP kp,
                         fact_db& facts, sched_db& scheds) ;
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_dynamic_rule";};
    virtual void dataCollate(collectData &data_collector) const ;
    virtual void reset_ctrl() { dflag = true ; }
  } ;

  
  class execute_dynamic_applyrule : public execute_drule_module {
  protected:
    rule_implP rp ;
    rule apply_tag, unit_tag ; 
    timeAccumulator timer ;

    KeySpaceP space ;

    std::vector<ExpandChain> input_chains ;
    ExpandChain output_chain ;

    // single source variables, with non local space variables removed
    std::vector<NonExpandUnit> input_ne_source ;
    // single constraint variables, with non local space variables removed
    std::vector<NonExpandUnit> input_ne_constraint ;

    const Map* target_pack ;
    const dMap* target_unpack ;

    // output map
    bool output_map ;

    bool output_cross_space ;
    std::string other_space_name ;
    // process entity partition for the target var
    const std::vector<entitySet>* output_ptn ;
    // the process rank in the target keyspace
    int target_process_rank ;
    MPI_Comm target_comm ;
    variable target_var ;
    storeRepP target_rep_in_facts ;

    size_t large_context_size ;

    bool dflag ;
    entitySet context ;
    entitySet target_var_exist ;
    vector<P2pCommInfo> send, recv ;
  public:
    execute_dynamic_applyrule(rule a, rule u, KeySpaceP kp,
                              fact_db& facts, sched_db& scheds) ;
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_dynamic_applyrule";};
    virtual void dataCollate(collectData &data_collector) const ;
    virtual void reset_ctrl() { dflag = true ; }
  } ;

  // this one is for dynamic param reduction
  class execute_dynamic_applyrule_param : public execute_drule_module {
  protected:
    rule_implP rp ;
    rule apply_tag, unit_tag ; 
    timeAccumulator timer ;

    KeySpaceP space ;

    std::vector<ExpandChain> input_chains ;
     // single source variables, with non local space variables removed
    std::vector<NonExpandUnit> input_ne_source ;
    // single constraint variables, with non local space variables removed
    std::vector<NonExpandUnit> input_ne_constraint ;

    // no output chain needed
    variable target_var ;
    storeRepP target_rep ;

    MPI_Op mpi_reduction_op ;

    size_t large_context_size ;

    static CPTR<joiner> apply_join_op ;

    bool dflag ;
    entitySet context ;
  public:
    execute_dynamic_applyrule_param(rule a, rule u, KeySpaceP kp,
                                    fact_db& facts, sched_db& scheds) ;
    ~execute_dynamic_applyrule_param() ;
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "execute_dynamic_applyrule";};
    virtual void dataCollate(collectData &data_collector) const ;
    static void create_user_function(void* send_ptr, void* result_ptr,
                                     int* size, MPI_Datatype* dptr) ;
    virtual void reset_ctrl() { dflag = true ; }
  } ;

  class execute_dclone_invalidator: public execute_modules {
    variable var, var_unique ;
    KeySpaceP self_clone ;
    std::vector<KeySpaceP> shadow_clone ;

    timeAccumulator timer ;
  public:
    execute_dclone_invalidator(variable v, variable vu,
                               KeySpaceP sc,
                               const std::vector<KeySpaceP>& sac,
                               fact_db& facts, sched_db& scheds) ;
    virtual void execute(fact_db& facts, sched_db& scheds) ;
    virtual void Print(std::ostream& s) const ;
    virtual string getName() {return "execute_dclone_invalidator" ;}
    virtual void dataCollate(collectData& data_collector) const ;
  } ;

  // execution module for keyspace redistribution
  class execute_keyspace_dist: public execute_modules {
    std::vector<KeySpaceP> spaces ;
    std::vector<variableSet> tunnels ;
    std::vector<variableSet> space_vars ;

    std::string space_names ;

    timeAccumulator timer ;
  public:
    execute_keyspace_dist(const std::vector<KeySpaceP>& s,
                          fact_db& facts, sched_db& scheds) ;
    virtual void execute(fact_db& facts, sched_db& scheds) ;
    virtual void Print(std::ostream& s) const ;
    virtual string getName() {return "execute_keyspace_dist" ;}
    virtual void dataCollate(collectData& data_collector) const ;
  } ;

  // this one is responsible for initializing all keyspaces
  class execute_init_keyspace: public execute_modules {
    timeAccumulator timer ;
  public:
    execute_init_keyspace(fact_db& facts, sched_db& scheds) ;
    virtual void execute(fact_db& facts, sched_db& scheds) ;
    virtual void Print(std::ostream& s) const ;
    virtual string getName() {return "execute_init_keyspace" ;}
    virtual void dataCollate(collectData& data_collector) const ;
  } ;
  
  class execute_insertion: public execute_modules {
    rule rule_tag ;
    insertion_rule_interfaceP rp ;
    KeyManagerP key_manager ;
    KeySpaceP space ;

    timeAccumulator timer ;
  public:
    execute_insertion(const rule& r, fact_db& facts, sched_db& scheds) ;
    virtual void execute(fact_db& facts, sched_db& scheds) ;
    virtual void Print(std::ostream& s) const ;
    virtual string getName() {return "execute_insertion" ;}
    virtual void dataCollate(collectData& data_collector) const ;
  } ;

  class execute_key_destruction: public execute_drule_module {
    rule rule_tag ;
    deletion_ruleP rp ;
    KeyManagerP key_manager ;
    KeySpaceP space ;
    std::vector<storeRepP> controlled_stores ;

    std::vector<ExpandChain> input_chains ;
    // single source variables, with non local space variables removed
    std::vector<NonExpandUnit> input_ne_source ;
    // single constraint variables, with non local space variables removed
    std::vector<NonExpandUnit> input_ne_constraint ;

    variable target ;
    storeRepP target_rep ;

    int large_context_size ;

    bool dflag ;
    entitySet context ;
    
    timeAccumulator timer ;
  public:
    execute_key_destruction(const rule& r, fact_db& facts, sched_db& scheds) ;
    virtual void execute(fact_db& facts, sched_db& scheds) ;
    virtual void Print(std::ostream& s) const ;
    virtual string getName() {return "execute_key_destruction" ;}
    virtual void dataCollate(collectData& data_collector) const ;
    virtual void reset_ctrl() { dflag = true ; }
  } ;

  class execute_erase: public execute_drule_module {
    rule rule_tag ;
    erase_ruleP rp ;
    KeySpaceP space ;

    std::vector<ExpandChain> input_chains ;
    // single source variables, with non local space variables removed
    std::vector<NonExpandUnit> input_ne_source ;
    // single constraint variables, with non local space variables removed
    std::vector<NonExpandUnit> input_ne_constraint ;
    // output variables, since erase rule cannot have output mapping,
    // there is no need to output_chains
    std::vector<NonExpandUnit> outputs ;

    variableSet targets ;

    int large_context_size ;

    bool dflag ;
    entitySet context ;
    
    timeAccumulator timer ;
  public:
    execute_erase(const rule& r, fact_db& facts, sched_db& scheds) ;
    virtual void execute(fact_db& facts, sched_db& scheds) ;
    virtual void Print(std::ostream& s) const ;
    virtual string getName() {return "execute_erase" ;}
    virtual void dataCollate(collectData& data_collector) const ;
    virtual void reset_ctrl() { dflag = true ; }
  } ;

  // this execution module resets the control bit of
  // dynamic rule at run-time
  class execute_dcontrol_reset: public execute_modules {
    variable var, var_unique ;
    std::vector<execute_druleP> drules ;

    std::vector<KeySpaceP> register_keyspaces ;

    timeAccumulator timer ;
  public:
    execute_dcontrol_reset(variable v, variable vu,
                           const std::vector<KeySpaceP>& ks,
                           fact_db& facts, sched_db& scheds) ;
    virtual void execute(fact_db& facts, sched_db& scheds) ;
    virtual void Print(std::ostream& s) const ;
    virtual string getName() {return "execute_dcontrol_reset" ;}
    virtual void dataCollate(collectData& data_collector) const ;
    void add_drule(execute_druleP drp) { drules.push_back(drp) ; }
  } ;

#endif
  class execute_rule_null : public execute_modules {
  protected:
    rule rule_tag ; 
  public:
    execute_rule_null(rule fi) : rule_tag(fi) {}
    virtual void execute(fact_db &facts, sched_db &scheds) {}
    virtual void Print(std::ostream &s) const
    {s << rule_tag << " over empty sequence."<< endl ;}
    virtual string getName() {return "execute_rule_null";};
    virtual void dataCollate(collectData &data_collector) const {}
  } ;

  class dynamic_schedule_rule: public execute_modules {
    int LBMethod ;
    
    rule_implP rp,main_comp,local_comp1 ;
    bool compress_set ;
    variableSet inputs, outputs ;
    fact_db backup_facts ;
    fact_db local_facts;
    rule_implP local_compute1;
    rule rule_tag ;
    entitySet given_exec_set ;
    entitySet exec_set ;
    timeAccumulator timer ;
    timeAccumulator comp_timer ;

    fact_db *facts1;
    fact_db *local_facts1;
    
    std::vector<double> workTime ;	// execution time of items belonging
                                  	// to proc i
    std::vector<double> aveWorkTime;	// average work time in proc i
    std::vector<double> ewt;		// expected work time (EWT) of proc i

    int foreMan;			// rank of loop scheduler
    std::vector<int> yMapSave;	// copy of item count, first item
    

    //=======================================================
    // used by IWS
    //=======================================================
    // Chunk information
    struct chunkInfo {
      entitySet chunkDef ;
      float chunkTime[2] ;
    } ;
    // Chunk Commuication data
    struct chunkCommInfo {
      int proc ;
      std::vector<int> chunkList ;
      int send_size,recv_size ;
    } ;

    // This tells us about the chunks formed by this processor
    std::vector<chunkInfo> chunkData ;
    // computation/communication schedule
    // chunks that stay on this processor
    std::vector<int> selfChunks ;
    // note: each processor will either be a sending or recving chunks,
    // not both
    // chunks to send
    std::vector<chunkCommInfo> sendChunks ;
    // chunks to recv 
    std::vector<chunkCommInfo> recvChunks ;
    // number of remote chunks this processor will process
    int numRemoteChunks ;
    // number of balancing steps
    int numBalances ;
    // chunk size
    int iwsSize;			// IWS chunk size
    
      
    int nCalls;			// no. of calls to execute()
    int allItems;			// total  no. of items
    int nProcs;			// no. of procs
    int myRank;			// process rank
    int myItemCount;		// local item count
    int myFirstItem;		// first local item

    double ideal1;		// perfect balance time
    double local1;		// time spent by this proc on own items
    double remote1;		// execution time of migrated items
    double wall1;		// local1 + remote1 + load balancing overheads
    double wall2;		// wall1 + any post-processing (i.e., MPI_Allreduce() )

    
    void ReceiveOutput (int src, int msgSize, int *iters, double *tTime) ;
    void SendOutput (int dest, int tStart, int tSize, double *tTime) ;

    void ReceiveInput (int src, int msgSize, int *tStart, int *tSize) ;

    void SendInput (int dest, int tStart, int tSize) ;
    void
    RecvInfo (int src, int action, int *chunkStart, int *chunkSize,
              int *chunkDest, double *mu) ;
    void
    SendInfo (int dest, int action, int chunkStart, int chunkSize,
              int chunkDest, double mu) ;

    void
    GetChunkSize (int method, int minChunkSize, int source,
                  int *yMap, int *chunkSize, int *batchSize, int *batchRem) ;
    
    void GetBuffer (int size) ;
    void AllocateLBspace (int nItems) ;
    void FreeLBspace () ;
    int inputPackSize (int tStart, int tSize) ;
    int outputPackSize (int tStart, int tSize) ;

    
    void iterative_weighted_static () ;
    void loop_scheduling () ;
    
  public:
    dynamic_schedule_rule(rule fi, entitySet eset, fact_db &facts, sched_db &scheds,int method) ;
    virtual ~dynamic_schedule_rule() ;
  
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() {return "dynamic_schedule_rule";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;
  
  class visitor ;
  
  class rule_compiler : public CPTR_type {
  public:
    ////////////////////
    virtual void accept(visitor& v) = 0 ;//method to accept a visitor
    ////////////////////
    virtual void set_var_existence(fact_db &facts, sched_db &scheds) = 0 ;
    virtual void process_var_requests(fact_db &facts, sched_db &scheds) = 0 ;
    virtual executeP create_execution_schedule(fact_db &facts, sched_db &scheds) = 0;
  } ;

  typedef CPTR<rule_compiler> rule_compilerP ;
  typedef std::map<rule, rule_compilerP > rulecomp_map ;

  struct decomposed_graph {
    multiLevelGraph mlg ;
    digraph::vertexSet loops, recursive, conditional ;
    decomposed_graph(digraph dg, digraph::vertexSet sources,
                     digraph::vertexSet targets) ;
  } ;

  struct graph_compiler {
    rule_compilerP fact_db_comm ;
    rulecomp_map rule_process ;
    rule baserule ;
    /////////////////////
    std::vector<rule> super_rules ;
    /////////////////////
    graph_compiler(decomposed_graph &deco, variableSet initial_vars) ;
    ///////////////////
    // visit order
    void top_down_visit(visitor& v) ;
    void bottom_up_visit(visitor& v) ;
    void compile(fact_db& facts, sched_db& scheds,
                 const variableSet& given, const variableSet& target) ;
    ///////////////////
    void existential_analysis(fact_db &facts, sched_db &scheds) ;
    executeP execution_schedule(fact_db &facts, sched_db &scheds,
                                const variableSet& alloc) ;
  } ;
  
 
   class execute_param_red : public execute_modules {
     vector<variable> reduce_vars ;
     vector<rule> unit_rules ;
     MPI_Op create_join_op ;
     vector<CPTR<joiner> >join_ops ;
     timeAccumulator timer ;
   public:
     execute_param_red(vector<variable> reduce_vars, vector<rule> unit_rules,
                       vector<CPTR<joiner> > join_ops) ; 
     ~execute_param_red() ;
     virtual void execute(fact_db &facts, sched_db &scheds) ;
     virtual void Print(std::ostream &s) const ;
     virtual string getName() {return "execute_param_red";};
     virtual void dataCollate(collectData &data_collector) const ;
  } ;

  class execute_chomp: public execute_modules {
    entitySet total_domain ;
    vector<pair<rule,rule_compilerP> > chomp_comp ;
    vector<pair<int,rule_implP> > chomp_compP ;
    std::deque<entitySet> rule_seq ;
    variableSet chomp_vars ;
    vector<vector<entitySet> > seq_table ;
    size_t chomp_size ;
    int_type chomp_iter ;
    vector<int_type> chomp_offset ;
    vector<storeRepP> chomp_vars_rep ;
    int_type D_cache_size ;
    timeAccumulator timer ;
    vector<timeAccumulator> comp_timers ;
    int execute_times;
  public:
    execute_chomp(const entitySet& td,
                  const vector<pair<rule,rule_compilerP> >& comp,
                  const std::deque<entitySet>& seq,
                  const variableSet& cv,
                  fact_db& facts);
    virtual void set_seq_table();
    virtual void execute(fact_db& facts, sched_db &scheds) ;
    virtual void Print(std::ostream& s) const ;
    virtual string getName() {return "execute_chomp";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  // experimental dynamic mapping generation
  // in the stationary time level
  void stationary_relation_gen(rule_db&, fact_db&, const variableSet&) ;
  
}
#endif

