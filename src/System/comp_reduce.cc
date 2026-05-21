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
#include "comp_tools.h"
#include "dist_tools.h"

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <map>
using std::map;

using std::ostream ;
using std::endl ;
using std::ostringstream ;

#include <Tools/hash_map.h>
#include "loci_globs.h"
using std::list ;

#include "thread.h"

//#define VERBOSE

namespace Loci {

  extern bool threading_global_reduction;
  extern bool threading_local_reduction;

  extern bool in_internal_query;
  extern int num_threaded_global_reduction;
  extern int num_total_global_reduction;
  extern int num_threaded_local_reduction;
  extern int num_total_local_reduction;

  entitySet vmap_source_exist_apply(const vmap_info &vmi, fact_db &facts,
                                    variable reduce_var, sched_db &scheds) {
    variableSet::const_iterator vi ;
    entitySet sources = ~EMPTY ;
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      if(*vi != reduce_var)
        sources &= scheds.variable_existence(*vi) ;
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

  void apply_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    existential_applyrule_analysis(apply, facts, scheds);
  }

  void apply_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    entitySet exec_seq = process_applyrule_requests(apply, unit_tag, output_mapping, facts, scheds);
    scheds.update_exec_seq(apply, exec_seq);
  }

  executeP apply_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    entitySet exec_seq = scheds.get_exec_seq(apply);
#ifdef PTHREADS
    if(apply.get_info().output_is_parameter) { // global reduction
      ++num_total_global_reduction;
      if(!in_internal_query && threading_global_reduction) {
        int tnum = thread_control->num_threads();
        int minw = thread_control->min_work_per_thread();
        if(!apply.get_info().rule_impl->thread_rule() ||
           exec_seq.size() < (size_t)tnum*minw)
          return new execute_rule(apply,sequence(exec_seq),facts,scheds);
        else {
          variableSet targets = apply.targets() ;
          if (targets.size() != 1) {
            cerr << "Apply rule has more than one target variables!!"
              << " threading schedule fails!!!" << endl;
            Loci::Abort();
          }
          variable t = *(targets.begin()) ;
          ++num_threaded_global_reduction;
          return new Threaded_execute_param_reduction
            (apply, exec_seq, t, facts, scheds);
        }
      } else {
        return new execute_rule(apply,sequence(exec_seq),facts,scheds);
      }
    } else { // local reduction
      ++num_total_local_reduction;
      if(!in_internal_query && threading_local_reduction) {
        int tnum = thread_control->num_threads();
        int minw = thread_control->min_work_per_thread();
        bool threadable = apply.get_info().rule_impl->thread_rule()
          && exec_seq.size() >= (size_t)tnum*minw;
        {
          rule_implP ti = apply.get_rule_implP() ;
          variableSet targets = apply.targets();
          storeRepP tr = ti->get_store(*(targets.begin()));
          if (tr->RepType() == BLACKBOX)
           threadable = false; 
        }
        if(!threadable) {
          return new execute_rule(apply,sequence(exec_seq),facts,scheds);
        } else {
          variableSet targets = apply.targets();
          if (targets.size() != 1) {
            cerr << "Apply rule has more than one target variables!!"
              << " threading schedule fails!!!" << endl;
            Loci::Abort();
          }
          ++num_threaded_local_reduction;
          return new Threaded_execute_local_reduction
            (apply,unit_tag,sequence(exec_seq),facts,scheds);
        }        
      } else
        return new execute_rule(apply,sequence(exec_seq),facts,scheds);
    }
#else
    executeP exec_rule =
      new execute_rule(apply, sequence(exec_seq), facts, scheds);
    return exec_rule ;
#endif
  }

  vector<CPTR<joiner> > global_join_ops ;
  void create_user_function(void *send_ptr, void *result_ptr, int *size, MPI_Datatype* dptr) {
    entitySet e ;
    sequence seq ;
    int unpack_send_position = 0 ;
    int unpack_result_position = 0;
    int pack_result_position = 0;
    for(size_t i = 0; i < global_join_ops.size(); i++) {
      storeRepP sp, tp ;
      sp = global_join_ops[i]->getTargetRep() ;
      tp = global_join_ops[i]->getTargetRep() ;
      sp->unpack(send_ptr, unpack_send_position, *size, seq) ;
      tp->unpack(result_ptr, unpack_result_position, *size, seq) ;
      global_join_ops[i]->SetArgs(tp, sp) ;
      global_join_ops[i]->Join(seq) ;
      tp->pack(result_ptr, pack_result_position, *size, e) ;
    }
  }

  execute_param_red::execute_param_red(vector<variable> red, vector<rule> unit, vector<CPTR<joiner> > j_op) {
    reduce_vars = red ;
    unit_rules = unit ;
    join_ops = j_op ;
    MPI_Op_create(&create_user_function, 0, &create_join_op) ;
  }

  execute_param_red::~execute_param_red() {
    MPI_Op_free(&create_join_op) ;
  }

  void joinHelper(vector<storeRepP> &results,
                  unsigned char *target, int target_size,
                  unsigned char *source, int source_size,
                  vector<CPTR<joiner> > &join_ops) {
    entitySet e = interval(0,0) ;
    sequence seq = sequence(e) ;

    int loc_target = 0 ;
    int loc_source = 0 ;
    
    for(size_t i = 0; i < join_ops.size(); i++) {
      storeRepP tp = join_ops[i]->getTargetRep() ;
      storeRepP sp = join_ops[i]->getTargetRep() ;
      tp->unpack(&target[0], loc_target, target_size, e);
      sp->unpack(&source[0], loc_source, source_size, e);
      join_ops[i]->SetArgs(tp, sp);
      join_ops[i]->Join(seq) ;
      results[i] = tp ;
    }
  }
   /** ************************************************************************
   *
   * @brief groupAllReduce() Performs a reduction on a group of parameters
   * and provides functionality that will work even with user defined 
   * dynamic types such as param<list<int> > 
   *
   * The implementation is based on the hypercube AllReduce algorithm where
   * the largest set of processors that can form a hypercube participate
   * The processors that don't fit into this are non-participating and 
   * contribute their contents to the reduction at the beginning and 
   * recieve the result at the end
   *
   * @param [sp] array of storeRep pointers to the input params
   * @param [join_ops] array of joiners for the communication
   * @param [comm] MPI communicator used to perform the operation
   *
   ** ************************************************************************/

  void groupAllReduce(vector<storeRepP> &sp,vector<CPTR<joiner> > &join_ops,
                      MPI_Comm comm) {
    int r=1, p=1 ;
    MPI_Comm_size(comm,&p) ;
    MPI_Comm_rank(comm,&r) ;

    FATAL(sp.size() != join_ops.size()) ;
    // if one processor, no neeed to continue
    if(p==1)
      return ;

    // find the ceil(log2(p))
    int logP = 0 ;
    while((1<<logP) < p)
      logP++ ;

    int dim = logP ;
    if((1<<logP) != p)
      dim-- ;
    
                          
    entitySet e = interval(0,0) ;
    sequence seq = sequence(e) ;
    int msgsize = 0 ;
    for(size_t i = 0; i < sp.size(); i++) 
      msgsize += sp[i]->pack_size(e);

    // allocate message buffers
    vector<unsigned char> sendbuf(msgsize,0), recvbuf(msgsize,0) ;
    vector<unsigned char> vsendbuf,vrecvbuf ;
    
    // pack initial message buffer
    int position = 0 ;
    for(size_t i = 0; i < sp.size(); i++) 
      sp[i]->pack(&sendbuf[0], position, msgsize, e) ;
    WARN(position != msgsize) ;

    MPI_Status status ;
    int recvsize = 0 ;
    int sendsize = sendbuf.size() ;
    int tag = 1 ;
    // If processors are not perfectly a power of 2, we first
    // reduce to the largest perfect power of 2, then perform
    // the hypercube reduction, then communicate the results
    // back
    bool non_participating = (r & (1<<dim))!=0 ;
    if(non_participating) {
      // These processors are not participating so communicate their contents
      // to their partner
      int partner = r ^(1<<dim) ;
      MPI_Send(&sendsize,1,MPI_INT,partner,tag,comm) ;
      MPI_Send(&sendbuf[0],msgsize,MPI_BYTE,partner,tag,comm) ;
    } else  if((r|(1<<dim)) < p) {
      // partners to the non-participating processors join their contributions
      // before the reduction proceeds
      int partner = r ^(1<<dim) ;
      MPI_Recv(&recvsize,1,MPI_INT,partner,tag,comm,&status)  ;
      if(recvsize >  int(recvbuf.size()))
        recvbuf.resize(recvsize) ;
      MPI_Recv(&recvbuf[0],recvsize,MPI_BYTE,partner,tag,comm,&status) ;
      vector<storeRepP> results(join_ops.size()) ;
      joinHelper(results,&recvbuf[0],recvsize,&sendbuf[0],sendsize,join_ops) ;
      int msgsize = 0 ;
      for(size_t i = 0; i < sp.size(); i++) 
        msgsize += results[i]->pack_size(e);

      if(msgsize > int(sendbuf.size())) 
        sendbuf.resize(msgsize) ;

      sendsize = msgsize ;
      position = 0 ;
      for(size_t i = 0; i < sp.size(); i++) 
        results[i]->pack(&sendbuf[0], position, sendsize, e) ;
      WARN(position != sendsize) ;
    }

    // null loop for non-participating processors
    int ldim = non_participating?0:dim ;
     
    for(int i=0;i<ldim;++i) {
      int partner = r ^ (1<<i) ;
      // If partner is a nonphysical processor, then find the processor
      // emulating the missing processor
      if(partner >= p)
        partner ^= 1<<(logP-1) ;


      MPI_Sendrecv(&sendsize,1,MPI_INT,partner,tag,
                   &recvsize,1,MPI_INT,partner,tag,
                   comm,&status); 
      if(recvsize >  int(recvbuf.size()))
        recvbuf.resize(recvsize) ;
        
      MPI_Sendrecv(&sendbuf[0],sendsize,MPI_BYTE,partner,tag,
                   &recvbuf[0],recvsize,MPI_BYTE,partner,tag,
                   comm,&status) ;


      vector<storeRepP> results(join_ops.size()) ;
      if(partner < r) 
        joinHelper(results,&sendbuf[0],sendsize,&recvbuf[0],recvsize,join_ops) ;
      else 
        joinHelper(results,&recvbuf[0],recvsize,&sendbuf[0],sendsize,join_ops) ;
      
      int msgsize = 0 ;
      for(size_t i = 0; i < sp.size(); i++) 
        msgsize += results[i]->pack_size(e);

      if(msgsize > int(sendbuf.size())) 
        sendbuf.resize(msgsize) ;

      sendsize = msgsize ;
      position = 0 ;
      for(size_t i = 0; i < sp.size(); i++) 
        results[i]->pack(&sendbuf[0], position, sendsize, e) ;
      WARN(position != sendsize) ;

    }
    if(non_participating) {
      // we recieve the combined result into the send buffer for unpacking
      int partner = r ^(1<<dim) ;
      MPI_Recv(&sendsize,1,MPI_INT,partner,tag,comm,&status)  ;
      if(sendsize >  int(sendbuf.size()))
        sendbuf.resize(sendsize) ;
      MPI_Recv(&sendbuf[0],sendsize,MPI_BYTE,partner,tag,comm,&status) ;
    } else  if((r|(1<<dim)) < p) {
      // send the result of the reduction to the non_participating processors
      int partner = r ^(1<<dim) ;
      MPI_Send(&sendsize,1,MPI_INT,partner,tag,comm) ;
      MPI_Send(&sendbuf[0],sendsize,MPI_BYTE,partner,tag,comm) ;
    }      
    position = 0 ;
    for(size_t i = 0; i < sp.size(); i++) 
      sp[i]->unpack(&sendbuf[0], position, sendsize, e) ;
    WARN(position != sendsize) ;
  }
#define GROUP_ALLREDUCE

  void execute_param_red::execute(fact_db &facts, sched_db &scheds) {
    //    debugout << "reduce vars=" ;
    //    for(size_t i = 0; i < reduce_vars.size(); i++) {
    //      debugout << ' ' << reduce_vars[i];
    //    }
    
    stopWatch s ;
    s.start() ;
    vector<storeRepP> sp;
    for(size_t i = 0; i < reduce_vars.size(); i++) {
      sp.push_back(facts.get_variable(reduce_vars[i])) ;
    }

#ifdef GROUP_ALLREDUCE
    groupAllReduce(sp, join_ops, facts.get_comm()) ;
#else
    // Old Code (Works only for fixed sized types)
    int size = 0;
    entitySet e = interval(0,0);
    sequence seq ;
    for(size_t i = 0; i < sp.size(); i++) {
      size += sp[i]->pack_size(e);
    }
    vector<unsigned char> result_ptr(size),send_ptr(size) ;
    int position = 0;
    for(size_t i = 0; i < sp.size(); i++) {
      sp[i]->pack(&send_ptr[0], position, size, e) ;
    }

    global_join_ops = join_ops;
    MPI_Allreduce(&send_ptr[0], &result_ptr[0], size,
                   MPI_PACKED, create_join_op,
                   facts.get_comm()) ;
    position = 0;
    for(size_t i = 0; i < sp.size(); i++) {
      sp[i]->unpack(&result_ptr[0], position, size, seq) ;
    }
#endif
    timer.addTime(s.stop(),1) ;
  }

  void execute_param_red::Print(ostream &s) const {
    if(verbose || get_exec_size() > 1)
      for(size_t i = 0 ; i < reduce_vars.size(); i++) {
        printIndent(s) ;
        s << "param reduction on " << reduce_vars[i] << endl ;
      }
  }

  void execute_param_red::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "param reduce: " ;

    variableSet vars  ;
    for(size_t i = 0 ; i < reduce_vars.size(); i++) 
      vars += reduce_vars[i] ;
    oss << vars ;

    data_collector.accumulateTime(timer,EXEC_COMMUNICATION,oss.str()) ;
  }

  void reduce_param_compiler::set_var_existence(fact_db &facts, sched_db &scheds)  {
    if(facts.isDistributed()) {
      
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      for(size_t i = 0; i < unit_rules.size(); i++) {
    	entitySet targets ;
	targets = scheds.get_existential_info(reduce_vars[i], unit_rules[i]) ;
    	targets += send_entitySet(targets, facts) ;
	targets &= d->my_entities ;
	targets += fill_entitySet(targets, facts) ;
	scheds.set_existential_info(reduce_vars[i],unit_rules[i],targets) ;
      }
    }
  }

  void reduce_param_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {
      for(size_t i = 0; i < unit_rules.size(); i++) {
	entitySet requests = scheds.get_variable_requests(reduce_vars[i]) ;
	requests += send_entitySet(requests, facts) ;
	scheds.variable_request(reduce_vars[i],requests) ;
      }
    }
  }

  executeP reduce_param_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {
      vector<variable> red ;
      vector<rule> ulist ;
      vector<CPTR<joiner> > jop ;

      for(size_t i=0;i<reduce_vars.size();++i) {
        if(GLOBAL_OR(scheds.get_variable_requests(reduce_vars[i])!=EMPTY)) {
          red.push_back(reduce_vars[i]) ;
          ulist.push_back(unit_rules[i]) ;
          jop.push_back(join_ops[i]) ;
        }
      }
      if(red.size() > 0) {
        executeP execute = new execute_param_red(red,ulist,jop);
        return execute;
      } else
        return 0 ;

    }
    if(verbose || facts.get_comm_size() > 1) {
      ostringstream oss ;
      for(size_t i = 0; i < reduce_vars.size(); i++)
        oss << "reduce param " << reduce_vars[i] << std::endl;
      executeP exec_msg = executeP(new execute_msg(oss.str())) ;
      return exec_msg;
    } else
      return 0 ;

  }

  void reduce_store_compiler::set_var_existence(fact_db &facts, sched_db &scheds)  {
  }

  void swap_send_recv(list<comm_info> &cl) {
    list<comm_info>::iterator li ;
    for(li=cl.begin();li!=cl.end();++li) {
      entitySet tmp  = li->send_set ;
      li->send_set = entitySet(li->recv_set) ;
      li->recv_set = sequence(tmp) ;
    }
  }

  //For duplication of work only.
  //It adds information of reduce_proc_able_entities and reduction_output_map to sched_db.
  void set_reduction_info(variableSet vars, sched_db &scheds, fact_db &facts) {
    entitySet reduce_filter = ~EMPTY ;
    fact_db::distribute_infoP d = facts.get_distribute_info() ;
    if(multilevel_duplication)
      //Because we are calling context_for_map_output multiple times,
      //we have made sure that we have full preimage of any output map
      //for comp_entities.
      //That means comp entities can be successfully computed for
      //any reduction rule
      reduce_filter = d->comp_entities ;
    else
      //Otherwise we can gurantee successful computation of only my_entities
      //for reduction rules
      reduce_filter = d->my_entities ;

    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      variable v = *vi ;
      //Find information for reduction variables
      ruleSet r = scheds.get_existential_rules(v);
      bool outputmap = false;

      entitySet reduce_proc_able_entities = ~EMPTY;
      for(ruleSet::const_iterator ri = r.begin();
	  ri != r.end(); ri++) {
	if(rule_has_mapping_in_output(*ri))
	  outputmap = true;
	if(ri->get_info().rule_impl->get_rule_class() == rule_impl::UNIT){
	  reduce_proc_able_entities &= scheds.get_proc_able_entities(v, *ri);
	}
      }
      entitySet tmpSet;
      if(outputmap)
	tmpSet = ~EMPTY;
      else
	tmpSet = EMPTY;

      //if mapping in output, reduce_proc_able_entities are based on minimal approach
      //which finds intersection entities that can be computed by
      //all apply and unit rules
      //if no mapping in output then, then it is union of entities produced by rules
      for(ruleSet::const_iterator ri = r.begin();
	  ri != r.end(); ri++) {
	if(ri->get_info().rule_impl->get_rule_class() == rule_impl::APPLY) {
	  if(!outputmap)
	    tmpSet |= scheds.get_proc_able_entities(v, *ri);
	  else
	    tmpSet &= scheds.get_proc_able_entities(v, *ri);
	}
      }

      reduce_proc_able_entities &= tmpSet;
      //If mapping in output, we may not be able compute entities
      //outside reduce_filter successfully
      if(outputmap)
	reduce_proc_able_entities &= reduce_filter;

      //Information needed for later use
      scheds.set_reduce_proc_able_entities(v, reduce_proc_able_entities);
      scheds.set_reduction_outputmap(v, outputmap);
    }
  }

  void reduce_store_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {
      std::list<comm_info> rlist;
      std::list<comm_info> clist;
      
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      variableSet vars ;
      vars += reduce_var ;

      //Find out duplication of variables that are associated with rules
      //that compute reduce variables
      if(duplicate_work) {
	set_reduction_info(vars, scheds, facts);
	set_duplication_of_variables(vars, scheds, facts);
      }
      list<comm_info> request_comm = barrier_process_rule_requests(vars,facts, scheds) ;
      entitySet requests = scheds.get_variable_requests(reduce_var) ;
      entitySet shadow;
      //Shadow should be empty for the variable which is going to be duplicated
      //Shadow is going to be entities need to be sent to the owner processor
      if(!duplicate_work || !scheds.is_duplicate_variable(reduce_var)) {
	requests += fill_entitySet(requests,facts) ;
	// pass requests across processors so we know to compute something
	// where the result will be sent to another processor
	scheds.variable_request(reduce_var,requests) ;
	shadow = scheds.get_variable_shadow(reduce_var) ;
	shadow &= requests ;
      }

      list<comm_info> slist ;
      //If variable is duplicate variable, we don't need to send request to the
      //other procesors because owner processor can ablways compute requests
      if(!duplicate_work || !scheds.is_duplicate_variable(reduce_var)) {
	entitySet response = send_requests(shadow, reduce_var,facts,slist) ;
	swap_send_recv(slist) ;
	rlist = sort_comm(slist,facts) ;
      }

      clist = sort_comm(request_comm,facts) ;
      
      scheds.update_comm_info_list(rlist, sched_db::REDUCE_RLIST);
      scheds.update_comm_info_list(clist, sched_db::REDUCE_CLIST);
      
#ifdef VERBOSE
      if(shadow != EMPTY) {
	debugout << "shadow = " << shadow << endl ;
	shadow -= d->my_entities ;
	debugout << "shadow/my_entites = " << shadow << endl ;
      }
#endif
    }
  }

  class execute_comm_reduce : public execute_modules {
    vector<pair<int,vector<send_var_info> > > send_info ;
    vector<pair<int,vector<recv_var_info> > > recv_info ;
    std::vector<std::vector<storeRepP> > send_vars ;
    std::vector<std::vector<storeRepP> > recv_vars ;
    CPTR<joiner> join_op ;
    int *maxr_size, *maxs_size, *r_size, *s_size, *recv_sizes ;
    unsigned char **recv_ptr , **send_ptr ;
    MPI_Request *request;
    MPI_Status *status ;
    timeAccumulator timer ;
  public:
    execute_comm_reduce(list<comm_info> &plist, fact_db &facts,
                        CPTR<joiner> jop) ;
    ~execute_comm_reduce() ;
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    virtual void Print(std::ostream &s) const ;
    virtual string getName() { return "execute_comm_reduce";};
    virtual void dataCollate(collectData &data_collector) const ;
  } ;

  execute_comm_reduce::execute_comm_reduce(list<comm_info> &plist,
                                           fact_db &facts,
                                           CPTR<joiner> jop) {

    join_op = jop ;
    HASH_MAP(int,vector<send_var_info> ) send_data ;
    HASH_MAP(int,vector<recv_var_info> ) recv_data ;
    list<comm_info>::const_iterator cli ;
    intervalSet send_procs, recv_procs ;
    for(cli=plist.begin();cli!=plist.end();++cli) {
      variable v = cli->v ;
      if(cli->send_set.size() > 0) {
        int send_proc = cli->processor ;
        send_procs += send_proc ;
        entitySet send_set = cli->send_set ;
        send_data[send_proc].push_back(send_var_info(v,send_set)) ;
      }
      if(cli->recv_set.size() > 0) {
        int recv_proc = cli->processor ;
        sequence recv_seq = cli->recv_set ;
        recv_procs += recv_proc ;
        recv_data[recv_proc].push_back(recv_var_info(v,recv_seq)) ;
      }
    }

    for(intervalSet::const_iterator ii=send_procs.begin();
        ii!=send_procs.end();
        ++ii) {
      send_info.push_back(make_pair(*ii,send_data[*ii])) ;
      send_vars.push_back(std::vector<storeRepP>()) ;
      for(size_t i=0;i<send_data[*ii].size();++i)
        send_vars.back().push_back(facts.get_variable(send_data[*ii][i].v)) ;
    }

    for(intervalSet::const_iterator ii=recv_procs.begin();
        ii!=recv_procs.end();
        ++ii) {
      recv_info.push_back(make_pair(*ii,recv_data[*ii])) ;
      recv_vars.push_back(std::vector<storeRepP>()) ;
      for(size_t i=0;i<recv_data[*ii].size();++i)
        recv_vars.back().push_back(facts.get_variable(recv_data[*ii][i].v)) ;
    }

    int nsend = send_info.size() ;
    int nrecv = recv_info.size() ;
    r_size = new int[nrecv] ;
    recv_sizes = new int[nrecv] ;
    maxr_size = new int[nrecv] ;
    maxs_size = new int[nsend] ;
    s_size = new int[nsend] ;
    for(int i = 0; i < nrecv; ++i) {
      r_size[i] = 0 ;
      maxr_size[i] = sizeof(int) ;
    }
    for(int i = 0; i < nsend; ++i) {
      maxs_size[i] = sizeof(int) ;
      s_size[i] = 0 ;
    }

    recv_ptr = new unsigned char*[max(nrecv,1)] ;
    send_ptr = new unsigned char*[max(nsend,1)] ;
    request =  new MPI_Request[nrecv] ;
    status =  new MPI_Status[nrecv] ;
  }

  execute_comm_reduce::~execute_comm_reduce() {
    delete [] maxr_size ;
    delete [] maxs_size ;
    delete [] recv_sizes ;
    delete [] r_size ;
    delete [] s_size ;
    delete [] recv_ptr ;
    delete [] send_ptr ;
    delete [] request ;
    delete [] status ;
  }

  static unsigned char *recv_ptr_buf = 0;
  static int recv_ptr_buf_size = 0;
  static unsigned char *send_ptr_buf = 0 ;
  static int send_ptr_buf_size = 0 ;

  void execute_comm_reduce::execute(fact_db  &facts, sched_db &scheds) {
    //    debugout << "reduce local vars" << endl ;
    stopWatch s ;
    s.start() ;
    const int nrecv = recv_info.size() ;
    int resend_size = 0, rerecv_size = 0 ;
    std::vector<int> send_index ;
    std::vector<int> recv_index ;
    int total_size = 0 ;
    MPI_Request *re_request = 0;
    MPI_Status *re_status = 0 ;
    for(int i=0;i<nrecv;++i) {
      r_size[i] = maxr_size[i] ;
      total_size += maxr_size[i] ;
    }

    if(recv_ptr_buf_size < total_size) {
      if(recv_ptr_buf)
        delete[] recv_ptr_buf ;
      recv_ptr_buf = new unsigned char[total_size] ;
      recv_ptr_buf_size = total_size ;
    }
    recv_ptr[0] = recv_ptr_buf ;
    for(int i=1;i<nrecv;++i)
      recv_ptr[i] = recv_ptr[i-1] + r_size[i-1] ;

    for(int i=0;i<nrecv;++i) {
      int proc = recv_info[i].first ;
      MPI_Irecv(recv_ptr[i], r_size[i], MPI_PACKED, proc, 1,
                facts.get_comm(), &request[i]) ;
    }
    total_size = 0 ;
    const int nsend = send_info.size() ;
    entitySet resend_procs, rerecv_procs ;
    for(int i=0;i<nsend;++i) {
      s_size[i] = 0 ;
      for(size_t j=0;j<send_info[i].second.size();++j) {
        //facts.get_variable(send_info[i].second[j].v) ;
	storeRepP sp = send_vars[i][j] ;

        s_size[i] += sp->pack_size(send_info[i].second[j].set) ;
      }
      if((s_size[i] > maxs_size[i]) || ( s_size[i] == sizeof(int))) {
	if(s_size[i] > maxs_size[i])
	  maxs_size[i] = s_size[i] ;
	int proc = send_info[i].first ;
	s_size[i] = sizeof(int) ;
	resend_procs += proc ;
	send_index.push_back(i) ;
      }
      total_size += maxs_size[i] ;
    }
    if(send_ptr_buf_size < total_size) {
      if(send_ptr_buf)
        delete[] send_ptr_buf ;
      send_ptr_buf = new unsigned char[total_size] ;
      send_ptr_buf_size = total_size ;
    }
    send_ptr[0] = send_ptr_buf ;
    for(int i = 1; i < nsend; i++)
      send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
    // Pack the buffer for sending
    for(int i=0;i<nsend;++i) {
      int loc_pack = 0 ;
      if(!resend_procs.inSet(send_info[i].first)) {
	for(size_t j=0;j<send_info[i].second.size();++j) {
	  storeRepP sp = send_vars[i][j] ;//facts.get_variable(send_info[i].second[j].v) ;
	  sp->pack(send_ptr[i], loc_pack,s_size[i],send_info[i].second[j].set);
	}
      }
      else
	MPI_Pack(&maxs_size[i], sizeof(int), MPI_BYTE,
                send_ptr[i], s_size[i], &loc_pack,
                facts.get_comm()) ;

    }
    // Send Buffer
    for(int i=0;i<nsend;++i) {
      int proc = send_info[i].first ;
      MPI_Send(send_ptr[i], s_size[i], MPI_PACKED, proc, 1,
               facts.get_comm()) ;
    }
    if(nrecv > 0) {
#ifdef DEBUG
      int err =
#endif
        MPI_Waitall(nrecv, request, status) ;
      FATAL(err != MPI_SUCCESS) ;
      for(int i = 0 ; i < nrecv; i++) {
        int rcv_sizes ;
	MPI_Get_count(&status[i], MPI_BYTE, &rcv_sizes) ;
	if(rcv_sizes == sizeof(int)) {
	  rerecv_procs += recv_info[i].first ;
	  recv_index.push_back(i) ;
	}
      }
    }
    for(int i=0;i<nrecv;++i) {
      int loc_unpack = 0;
      if(rerecv_procs.inSet(recv_info[i].first)) {
	int temp ;
	MPI_Unpack(recv_ptr[i], r_size[i], &loc_unpack,
                   &temp, sizeof(int), MPI_BYTE,
                   facts.get_comm()) ;
	if(temp > maxr_size[i])
	  maxr_size[i] = temp ;
      }
      else
	for(size_t j=0;j<recv_info[i].second.size();++j) {
	  storeRepP sp = recv_vars[i][j] ;//facts.get_variable(recv_info[i].second[j].v) ;
	  storeRepP sr = sp->new_store(EMPTY) ;
	  sr->allocate(entitySet(recv_info[i].second[j].seq)) ;
	  sr->unpack(recv_ptr[i], loc_unpack, r_size[i],
		     recv_info[i].second[j].seq) ;
	  CPTR<joiner> op = join_op->clone() ;
	  op->SetArgs(sp,sr) ;
	  op->Join(recv_info[i].second[j].seq) ;
	}
    }
    rerecv_size = rerecv_procs.size() ;
    resend_size = resend_procs.size() ;

    if(rerecv_size > 0) {
      re_request =  new MPI_Request[rerecv_size] ;
      re_status =  new MPI_Status[rerecv_size] ;
    }
    for(int i = 0; i < rerecv_size; i++) {
      int proc = recv_info[recv_index[i]].first ;
      recv_ptr[recv_index[i]] = new unsigned char[maxr_size[recv_index[i]]] ;
      MPI_Irecv(recv_ptr[recv_index[i]],
                maxr_size[recv_index[i]], MPI_PACKED,
                proc, 2, facts.get_comm(),
                &re_request[i]) ;
    }

    for(int i=0;i<resend_size;++i) {
      int loc_pack = 0 ;
      send_ptr[send_index[i]] = new unsigned char[maxs_size[send_index[i]]] ;
      for(size_t j=0;j<send_info[send_index[i]].second.size();++j) {
	storeRepP sp = send_vars[i][j] ;//facts.get_variable(send_info[send_index[i]].second[j].v) ;
	sp->pack(send_ptr[send_index[i]], loc_pack,maxs_size[send_index[i]],send_info[send_index[i]].second[j].set);
      }
    }

    // Send Buffer
    for(int i=0;i<resend_size;++i) {
      int proc = send_info[send_index[i]].first ;
      MPI_Send(send_ptr[send_index[i]],
               maxs_size[send_index[i]], MPI_PACKED,
               proc, 2, facts.get_comm()) ;
      delete [] send_ptr[send_index[i]] ;
    }
    if(rerecv_size > 0) {
#ifdef DEBUG
      int err =
#endif
        MPI_Waitall(rerecv_size, re_request, re_status) ;
      FATAL(err != MPI_SUCCESS) ;
    }
    for(int i=0;i<rerecv_size;++i) {
      int loc_unpack = 0;
      for(size_t j=0;j<recv_info[recv_index[i]].second.size();++j) {
	storeRepP sp = recv_vars[i][j] ;//facts.get_variable(recv_info[recv_index[i]].second[j].v) ;
	storeRepP sr = sp->new_store(EMPTY) ;
	sr->allocate(entitySet(recv_info[recv_index[i]].second[j].seq)) ;
	sr->unpack(recv_ptr[recv_index[i]], loc_unpack, maxr_size[recv_index[i]],
		   recv_info[recv_index[i]].second[j].seq) ;

	CPTR<joiner> op = join_op->clone() ;
	op->SetArgs(sp,sr) ;
	op->Join(recv_info[recv_index[i]].second[j].seq) ;
      }
      delete [] recv_ptr[recv_index[i]] ;
    }
    if(rerecv_size > 0) {
      delete [] re_status ;
      delete [] re_request ;
    }
    timer.addTime(s.stop(),1) ;
  }

  void execute_comm_reduce::Print(ostream &s) const {
    int sz = 0 ;
    if(send_info.size()+recv_info.size() > 0) {
      printIndent(s) ;
      s << "reduction block {" << endl ;
      printIndent(s) ;
      if(send_info.size() > 0) {
        s << "Send:" << endl ;
        printIndent(s) ;
        for(size_t i=0;i<send_info.size();++i) {
          for(size_t j=0;j<send_info[i].second.size();++j) {
            s << send_info[i].second[j].v << "  " ;
	    sz += (send_info[i].second[j].set).size() ;
	  }
	  s << " to " << send_info[i].first << endl ;
          printIndent(s) ;
        }
	s << " Total entities sent = " << sz << endl ;
      }
      sz = 0 ;
      if(recv_info.size() > 0) {
        printIndent(s) ;
        s << "Recv:" << endl ;
        for(size_t i=0;i<recv_info.size();++i) {
          for(size_t j=0;j<recv_info[i].second.size();++j) {
            s << recv_info[i].second[j].v << "  " ;
	    sz += (recv_info[i].second[j].seq).size() ;
	  }
          s << " from " << recv_info[i].first << endl ;
          printIndent(s) ;
        }
        printIndent(s) ;
	s << " Total entities recieved = " << sz << endl ;
      }
      s << "}" << endl ;
    }
  }

  void execute_comm_reduce::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "comm reduce: " ;

    variableSet vars  ;
    for(size_t i=0;i<send_info.size();++i)
      for(size_t j=0;j<send_info[i].second.size();++j) 
        vars += send_info[i].second[j].v ;

    for(size_t i=0;i<recv_info.size();++i) 
      for(size_t j=0;j<recv_info[i].second.size();++j)
        vars += recv_info[i].second[j].v ;

    oss << vars ;

    data_collector.accumulateTime(timer,EXEC_COMMUNICATION,oss.str()) ;
  }
  
  executeP reduce_store_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {
      variableSet vars;
      vars += reduce_var;
      std::list<comm_info> clist = scheds.get_comm_info_list(vars, facts, sched_db::REDUCE_CLIST);
      std::list<comm_info> rlist = scheds.get_comm_info_list(vars, facts, sched_db::REDUCE_RLIST);

      CPTR<execute_sequence> el = new execute_sequence ;
      if(!rlist.empty()) {
        executeP exec_comm_reduce = new execute_comm_reduce(rlist, facts, join_op);
        el->append_list(exec_comm_reduce);
      }
      execute_comm2::inc_comm_step() ;
      if(!clist.empty()) {
        //executeP exec_comm = new execute_comm(clist, facts);
        executeP exec_comm2 = new execute_comm2(clist, facts);
        el->append_list(exec_comm2) ;
        //el->append_list(exec_comm) ;
      }
      if(verbose || facts.get_comm_size() > 1) {
        ostringstream oss ;
        oss << "reduce store " << reduce_var ;
        executeP exec_msg = new execute_msg(oss.str());
        el->append_list(exec_msg) ;
      }
      return executeP(el) ;
    }

    if(verbose || facts.get_comm_size() > 1) {
      ostringstream oss ;
      oss << "reduce store " << reduce_var ;
      executeP exec_msg2 = new execute_msg(oss.str());
      return executeP(exec_msg2) ;
    }
    return executeP(0) ;
  }
}
