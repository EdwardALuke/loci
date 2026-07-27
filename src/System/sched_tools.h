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

