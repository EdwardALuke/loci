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
#ifndef THREAD_H
#define THREAD_H

// disable this file by default for now...
#ifdef PTHREADS

#include <variable.h>
#include <constraint.h>
#include <rule.h>
#include <execute.h>
#include <distribute.h>
#include <Tools/cptr.h>

#include <pthread.h>
#include <semaphore.h>
#include <sched.h>

#include <vector>
#include <map>
#include <set>
#include <deque>
#include <string>
#include <exception>
#include <utility>
#include <cctype>

#include "comp_tools.h"
#include "sched_tools.h"
#include "loci_globs.h"

#ifndef PAPI_DEBUG
#define long_long long
#endif

namespace Loci {

  extern int current_rule_id;
  extern int num_thread_blocks;

  // a simple thread exception
  struct ThreadException : public std::exception {
    ThreadException(const std::string& s="thread exception"):msg(s) {}
    ~ThreadException() throw() {} // make g++ happy
    virtual const char* what() const throw()
    { return msg.c_str(); }
  private:
    std::string msg;
  };

  // this is a convenience type for a thread block ID.  it is just a
  // pair of integers, the first being the thread id, the second being
  // the block number for that particular thread
  typedef std::pair<int,int> TBlockID;
  class ThreadPartition_simple;
  // reduction lock is a special type of thread lock used only inside
  // the partial reduction (or some chomping) rule execution.  It uses
  // an internal conflict table to coordinate the execution of thread
  // blocks.  before a thread executes a block, it will need to acquire
  // the reduction lock by passing in the thread and block IDs.  after
  // the block finishes execution, the thread will release the reduction
  // lock.  the exact working mechanism is not exposed to any working
  // thread.  the reduction lock contains a conflict table internally
  // for all involved threads/blocks, upon a locking request, it will
  // find the corresponding internal lock (or locks) and will try to 
  // lock it (or them) it will return a status to the requesting 
  // thread whether the lock is successfully obtained.
  class ReductionLock {
    friend class ThreadPartition_simple;
    // a static data structure that holds a lock for every block
    // in every thread.  this part is shared globally.  each thread
    // rule will have a different lock pattern apon which the locks
    // are acquired internally
    static pthread_spinlock_t** locks;
    // this is the actual lock pattern, each block in every thread
    // holds a list of the internal locks to acquire, NOTE: the lock
    // list is ORDERED and the lock acquisition must obey this order
    // to avoid thread synchonization problems.
    std::vector<std::vector<std::vector<pthread_spinlock_t*> > >
      lock_pattern;
  public:
    // returns true if lock is successful, false if not
    bool try_lock(int tid, int bid);
    void unlock(int tid, int bid);
    // this is for debugging and statistical info use, returns
    // the number of fine lock entries a thread block has to obtain
    // while calling "try_lock"
    int try_lock_entries(int tid, int bid);
  };

  // a utility function that partitions a sequence into n parts with
  // equal size
  std::vector<sequence> equal_cut_seq(const sequence& s, int n);
  // a utility function that randomizes a passed in sequence
  sequence randomize_seq(const sequence& s);
  // this structure manages the per thread entity partition
  // created during the initial stages. these partitions are created
  // for each thread in order to maximize the locality
  // for later thread execution
  // this is an interface...
  //
  // Here we will have two kinds of partitions: 1) domain --- this is 
  // the total entity set a thread will get, i.e., all the entities
  // a thread is supposed to work on, 2) block --- this is a further
  // partition within each thread domain, the purpose of such blocks
  // are to deal with thread conflicts within the local reduction and
  // chomping scheduling, they may also be useful for future enhancements
  // such as dynamic scheduling with load balancing, etc.
  //
  // Right now, both the domains and blocks are created and maintained 
  // by the ThreadPartition interface.  this has the advantage of keeping
  // just a single partition per MPI process (saving a lot of space).
  // There are several notes: 1) together with the MPI partition, now
  // we have three kinds of partitions: MPI domain, thread domain, and
  // thread blocks. In the future, it would be helpful to generate the
  // entity local numbering considering all these partitions, this would
  // mean a more compact and consistent local number in all levels of
  // partition schemes.  2) in blocked thread execution, each thread 
  // might have different execution sequence than the maintained blocks,
  // e.g., we have blocks a, b, and c for thread T, a=[1,10], b=[11,20],
  // and c=[21,30].  For a particular rule, the execution sequence
  // might be: [5, 25] for example, then the rule will execute over
  // blocks a, b, and c for thread T.  However the exact execution 
  // sequence will be [5,10] for block a, full block b, and [21,25] for
  // block c.  currently to obtain the exact information, we will 
  // compute the exact execution sequence every time before the thread
  // runs, otherwise, we will have to save the information, which will
  // likely be very costly in terms of storage requirement.  we believe
  // the small amount of calculation to figure out the exact sequence
  // each time before execution should not be a performance problem.
  // 3) actually we may have other levels of partition even within the
  // blocks.  for example, chomping rules might further decompose a 
  // single block into multiple chunks so that each chunk is executed 
  // as a unit to better suit the cache performance optimization.  
  // however such further partition is solely the decision of each
  // thread rule and the execution coordination below the block level
  // is maintained by the thread rules themselves.  the block level
  // is mainly used to guarantee the thread synchronization.
  struct ThreadPartition {
    virtual ~ThreadPartition() {}
    // this interface generates the partitioned 
    // domain and block entities for all threads
    virtual void create_thread_domain(fact_db& facts) = 0;
    // this is for generate thread-domain partitioned execution sequence
    // for a threaded rule.
    virtual std::vector<sequence> 
      partition_by_domain(const sequence& s) const = 0;
    // this is to generate blocks for a given sequence on a given thread.
    virtual std::vector<int>
      partition_by_blocks(int t, const sequence& s) const = 0;
    // this is for retrieve the compilete block list for any thread
    virtual const std::vector<entitySet>&
      get_thread_block_list(int t) const = 0;
    // given a thread id and a block id of that thread, and an
    // execution sequence, this method returns the "TRUE" execution
    // sequence contained within the block, i.e., shrink the sequence
    // so that it is entirely contained within the block
    virtual sequence refine_seq_by_block
      (int t, int b, const sequence& s) const = 0;
    virtual entitySet refine_seq_by_block
      (int t, int b, const entitySet& s) const = 0;
    // this method generates a reduction lock for all the
    // passed-in blocks and a chain of maps.
    //
    // the passed in blocks are organized as a 2D array that represents
    // the block IDs for each thread, for example:
    // 0,2,3 | 2,5 | 1,4,6
    // means we have three threads, thread 1 uses blocks 0,2,3
    // thread 2 uses blocks 2 and 5, thread 3 uses blocks 1,4,6.
    // 
    // this method needs to be in the thread partitioner class because
    // the actual blocks contents are only known to the partitioner
    // all the other parts only speak in terms of block IDs
    virtual ReductionLock generate_conflict_lock
      (const std::vector<std::vector<int> >& blocks_id,
       // map_info is the chain of maps contained in a reduction target
       // we take in a vector of these because a thread could execute
       // a series of reduction rules as a group
       const std::vector<vmap_info>& map_info,
       fact_db& facts, sched_db& scheds) const = 0;
    virtual ReductionLock generate_empty_conflict_lock() const = 0;
  };

  // a concrete thread entity set partitioner based on the simple
  // scheme that cuts the cells equally among the threads and then
  // associates the faces and nodes.
  struct ThreadPartition_simple: public ThreadPartition {
    int thread_num; // number of threads
    // number of blocks, defaults to 10 per thread domain
    // this is the number we will try to achieve for each thread,
    int initial_block_num; 
    // this is a distribution of actual block numbers per thread domain
    std::vector<int> real_block_num;
    // this records the total number of blocks
    int total_blocks;
    entitySet whole_entities; // all the entities managed
    std::vector<entitySet> ptn; // the partition table
    std::vector<std::vector<entitySet> > blocks; // the block table
    ThreadPartition_simple(int tn, int bn=num_thread_blocks)
      :thread_num(tn),initial_block_num(bn),
      real_block_num(std::vector<int>(tn,0)),total_blocks(0),
      blocks(std::vector<std::vector<entitySet> >(tn)) {}
    ~ThreadPartition_simple();
    // the implementation
    void create_thread_domain(fact_db& facts);
    std::vector<sequence> partition_by_domain(const sequence& s) const;
    std::vector<int> partition_by_blocks(int t, const sequence& s) const;
    sequence refine_seq_by_block(int t, int b, const sequence& s) const;
    entitySet refine_seq_by_block(int t, int b, const entitySet& s) const;
    const std::vector<entitySet>& get_thread_block_list(int t) const;
    ReductionLock generate_conflict_lock
      (const std::vector<std::vector<int> >& blocks_id,
       const std::vector<vmap_info>& map_info,
       fact_db& facts,sched_db& scheds) const;
    ReductionLock generate_empty_conflict_lock() const;
  private:
    void partition_topology(constraint& geom_cells,
        const Map& cl, const Map& cr, const multiMap& face2node,
        const entitySet& limit, int p, std::vector<entitySet>& ptn);
  };

  // this is the abstract interface for all the thread execution modules
  // all original sequential execution modules that could be made
  // multi-threaded would generate corresponding threaded version as
  // described later in this file
  class ThreadedExecutionModule : public CPTR_type {
  public:
    virtual void execute(fact_db& facts, sched_db& scheds) = 0;
    // print execution sequence for the current thread
    virtual void print_seq(std::ostream& s) const = 0;
    // return the timing information of this execution module
    virtual double get_time() const = 0;
    // get cache misses
    virtual long_long get_l1_dcm() const = 0;
    virtual long_long get_l2_dcm() const = 0;
    virtual ~ThreadedExecutionModule() {}
  };
  typedef CPTR<ThreadedExecutionModule> ThreadedEMP;

  // an abstract interface for thread control.
  // this is used to globally control threads behavior, such as
  // creating/starting/stopping threads, do thread sync, etc.
  //
  // this interface is only used by the main control thread, which
  // does NOT participate in the work performed by all the work
  // threads created through this interface.
  class ThreadControl {
  public:
    // this function creates all the work threads
    // (the number of the threads is determined at the construction time)
    // Each work thread will be in a suspended state
    // upon creation (essentially sits there idly waiting to be fed
    // a thread execution module). 
    virtual void create_threads() = 0;
    // this function terminates all created work threads.
    virtual void shutdown() = 0;
    // this function restarts all suspended threads, feeding each
    // a new thread execution module.
    virtual void restart(std::vector<ThreadedEMP>& tem) = 0;
    // this function causes all the work threads to be sequentialized
    // i.e., thread 0 executes first, and then thread 1 executes,
    // and so on... this is mainly supplied as a debugging tool.
    virtual void sequential_restart(std::vector<ThreadedEMP>& tem) = 0;
    // this function causes the main control thread to suspend
    // waiting for all the work-thread to finish their current
    // execution module, and then wakes up.
    // It also means that after this function returns, all working
    // threads will be in a suspended state.
    virtual void wait_threads() = 0;
    // returns the total number of work threads created so far
    virtual int num_threads() const = 0;
    // sets up the the global fact and sched database for use in execution.
    virtual void setup_facts_scheds(fact_db&, sched_db&) = 0;
    // this just returns the minimum work per thread (as entity size)
    // i.e., if the total work is less than this value times the
    // number of threads, then it is not worth to lauch the work on
    // multiple threads.
    virtual int min_work_per_thread() const = 0;
    // given a thread id, this function will return the reference
    // to a vector of thread id that represents the (global) reduction
    // partner of the given thread
    virtual const std::vector<int>& get_reduction_partners(int id) const = 0;
    // returns the thread id of the reduction tree root
    virtual int get_reduction_root() const = 0;

    // this one determines if a calling thread is the leading thread
    virtual bool is_leading_thread() const = 0;

    // these provide a global lock on all threads managed by 
    // the thread control
    virtual void atomic_begin() = 0;
    virtual void atomic_end() = 0;

    // this provides a barrier for all the threads controlled
    // by the thread controller
    virtual void global_barrier() = 0;

    // this method provides the local ID for each calling thread
    virtual std::string get_local_id() const = 0;

    // this method returns a thread partition for a sequence
    // passed in.  this is intended to be used with the thread
    // entitySet partition manager to generate a locality aware
    // thread entitySet partition.
    virtual std::vector<sequence>
      partition_seq(const sequence& s) const = 0;

    // this method takes a thread ID and a execution sequence, it
    // returns the list of blocks within that thread's domain that
    // intersect the passed in sequence.  this just breaks the 
    // execution sequence into blocks.
    virtual std::vector<int>
      block_thread_seq(int t, const sequence& s) const = 0;

    // shrinks an execute sequence to be within a specific thread block
    virtual sequence refine_seq_by_block
      (int t, int b, const sequence& s) const = 0;
    virtual entitySet refine_seq_by_block
      (int t, int b, const entitySet& s) const = 0;

    // creates a conflict lock for partial reduction
    virtual ReductionLock create_empty_conflict_lock() const = 0;
    virtual ReductionLock create_conflict_lock
      (const std::vector<std::vector<int> >& blocks,
       const vmap_info& map_info,
       fact_db& facts, sched_db& scheds) const = 0;
    virtual ReductionLock create_conflict_lock
      (const std::vector<std::vector<int> >& blocks,
       const std::vector<vmap_info>& map_info,
       fact_db& facts, sched_db& scheds) const = 0;

    virtual ~ThreadControl() {}
  };

  extern ThreadControl* thread_control;

  // this is an implementation of the ThreadControl interface
  // based on the POSIX thread and semaphore implementation.
  class ThreadControl_pthread : public ThreadControl {
  public:
    ThreadControl_pthread(int n,fact_db& facts, sched_db& scheds);
    ~ThreadControl_pthread();
    void create_threads();
    void shutdown();
    void restart(std::vector<ThreadedEMP>&);
    void sequential_restart(std::vector<ThreadedEMP>&);
    void wait_threads();
    int num_threads() const { return tnum; }
    void setup_facts_scheds(fact_db& facts, sched_db& scheds)
    { factsP = &facts; schedsP = &scheds; }
    // 10 is just an arbitrary value for now...
    // we may want to define a better measure in the future...
    int min_work_per_thread() const { return 10; }
    const std::vector<int>& get_reduction_partners(int id) const;
    int get_reduction_root() const { return 0; }
    bool is_leading_thread() const;
    void atomic_begin() { global_spin.lock(); }
    void atomic_end()  { global_spin.unlock(); }
    void global_barrier();
    std::string get_local_id() const;
    std::vector<sequence> partition_seq(const sequence& s) const;
    std::vector<int> block_thread_seq(int t, const sequence& s) const;
    sequence refine_seq_by_block(int t, int b, const sequence& s) const;
    entitySet refine_seq_by_block(int t, int b, const entitySet& s) const;
    ReductionLock create_empty_conflict_lock() const;
    ReductionLock create_conflict_lock
      (const std::vector<std::vector<int> >& blocks,
       const vmap_info& map_info,
       fact_db& facts, sched_db& scheds) const;
    ReductionLock create_conflict_lock
      (const std::vector<std::vector<int> >& blocks,
       const std::vector<vmap_info>& map_info,
       fact_db& facts, sched_db& scheds) const;

  private:
    // disable copy and assignment
    ThreadControl_pthread(const ThreadControl_pthread&);
    ThreadControl_pthread& operator=(const ThreadControl_pthread&);
    // the thread function assigned to eacher worker
    static void* thread_fun(void* arg);
    // whether the work threads are active or not.
    bool active;
    // flag to indicate shutdown
    bool stop;
    // this would point to the actual work for each work thread
    std::vector<ThreadedEMP>* work;
    // pointers to the current fact and sched database
    fact_db* factsP;
    sched_db* schedsP;
    // a small structure to build the argument to the thread function
    struct TArg {
      int tid;                // internal thread id, from 0 to n-1
      pthread_t ptid;           // thread id assigned by the pthread lib
      ThreadControl_pthread* self; // pointer to the class object created
                                   // used to access the object internals
                                   // from each thread
    };
    // records the number of threads created
    int tnum;
    std::vector<TArg> thread_args;
    // the pthreads attributes when creating threads
    pthread_attr_t pattr;
    // semaphores used to block the work threads
    std::vector<sem_t> worker_block;
    // semaphore used to block the main thread
    sem_t main_block;
    // this is used to parallel start all threads
    std::vector<std::vector<int> > start_partners;
    // these types and data structures are used to signal the
    // finishing of all the work threads
    struct Terminate {
      bool done;           // finish flag
      pthread_spinlock_t done_lock; // spin locks that guard "done"
      int p; // the parent node
      std::vector<int> c; // all the child nodes, could be 0, 1, or 2
    };
    std::vector<Terminate> finish_notification;
    std::vector<vector<int> > reduction_partners;
    // the pthread id of the main and the leading work thread
    pthread_t main_pid;
    pthread_t lead_work_pid;
    lmutex global_spin;
    pthread_barrier_t barrier;
    std::map<pthread_t,std::string> id_map;
    void build_term_tree(int pid, int depth, size_t& myid);
    // pointer to the thread entity partition manager
    // this is used to initialize a thread entity partition
    ThreadPartition* tpn;
  };

  // this execution module creates all work threads
  class StartThreads : public execute_modules {
  public:
    StartThreads() {}
    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      thread_control->create_threads();
      timer.addTime(s.stop(),1);
    }
    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << "starting threads" << std::endl;
    }
    std::string getName() { return "StartThreads"; }
    void dataCollate(collectData& data_collector) const
    {data_collector.accumulateTime(timer,EXEC_CONTROL,"start_threads");}
  protected:
    timeAccumulator timer;
  };
  
  // this execution module causes all work threads to be shutdown,
  // and the thread control object to be deleted.
  class ShutDownThreads : public execute_modules {
  public:
    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      thread_control->shutdown();
      delete thread_control;
      thread_control = 0;
      timer.addTime(s.stop(),1);
    }
    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << "shutting down threads" <<std::endl;
    }
    std::string getName() { return "ShutDownThreads"; }
    void dataCollate(collectData& data_collector) const
    {data_collector.accumulateTime(timer,EXEC_CONTROL,"shutdown_threads");}
  protected:
    timeAccumulator timer;
  };

  // we will then create threaded execution modules for the various
  // Loci computations that can be made multithreaded.

  // this is the module for threading pointwise rules
  class ExecuteThreaded_pointwise: public ThreadedExecutionModule {
  public:
    ExecuteThreaded_pointwise(const sequence& s, executeP er)
      :seq(s), exec_size(s.size()), exec_rule(er) {}
    void execute(fact_db& facts, sched_db& scheds)
    { 
      stopWatch s;
      s.start();
      exec_rule->execute_kernel(seq);
      timer.addTime(s.stop(), exec_size);
    }
    void print_seq(std::ostream& s) const
    {
      if(verbose || seq.num_intervals() < 4) {
        s << seq;
      } else {
        s << "{" << seq[0] << "...}";
      }
    }
    double get_time() const { return timer.getTime(); }
    long_long get_l1_dcm() const { return 0; } // to be implemented!!
    long_long get_l2_dcm() const { return 0; } // to be implemented!!
  private:
    // the execution sequence
    sequence seq;
    // the size of the execution sequence
    size_t exec_size;
    // we will just delegate the execution to the underlying concrete rule
    executeP exec_rule;
    // a timer to record the execution time
    timeAccumulator timer;
  };

  // this is a module that hooks up the compiler and the threaded
  // execution module for threaded pointwise rules defined above
  class Threaded_execute_rule: public execute_modules {
  public:
    Threaded_execute_rule(rule r, const sequence& s,
                          fact_db& facts, sched_db& scheds)
      :rule_name(r), seq(s)
    {
      exec_rule = new execute_rule(r, s, facts, scheds);
      int tnum = thread_control->num_threads();
      std::vector<sequence> partition =
        thread_control->partition_seq(seq);
      tot_mem = 0;
      for(int i=0;i<tnum;++i) {
        texec.push_back(new ExecuteThreaded_pointwise(partition[i], 
                                                      exec_rule));
        tot_mem += 16*partition[i].num_intervals();
      }
      tot_mem += 8*seq.num_intervals();
      // DEBUG
      // Loci::debugout << "--------" << r << std::endl;
      // Loci::debugout << "total seq intervals: "
      //                << seq.num_intervals()
      //                << " size: " << seq.size() << std::endl;
      // for(int i=0;i<tnum;++i)
      //   Loci::debugout << "thread " << i << " intervals: "
      //                  << partition[i].num_intervals()
      //                  << " size: " << partition[i].size() << std::endl;
      ////////
    }
    void execute(fact_db& facts, sched_db& scheds)
    {
      current_rule_id = rule_name.ident();
      stopWatch s;
      s.start();

      exec_rule->execute_prelude(seq);
      timer_comp.addTime(s.stop(), seq.size());

      s.start();
      thread_control->setup_facts_scheds(facts,scheds);
      thread_control->restart(texec);
      timer_ctrl.addTime(s.stop(),seq.size());

      thread_control->wait_threads();

      s.start();
      exec_rule->execute_postlude(seq);
      timer_comp.addTime(s.stop(), seq.size());
      current_rule_id = 0;
    }
    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << rule_name << " over sequence ";
      if(verbose || seq.num_intervals() < 4)
        s << seq;
      else
        s << "{" << seq[0] << "...}";
      s << " using " << texec.size() << " threads";
      if(verbose) {
        s << " (";
        s << "thread " << 0 << ": ";
        texec[0]->print_seq(s);
        for(size_t i=1;i<texec.size();++i) {
          s << ", thread " << i << ": " ;
          texec[i]->print_seq(s);
        }
        s << ")";
      }
      s << std::endl;
    }

    std::string getName() { return "Threaded_execute_rule"; }

    void dataCollate(collectData& data_collector) const
    {
      std::ostringstream oss;
      oss << "rule: " << rule_name;
      // control time
      data_collector.accumulateTime(timer_ctrl, EXEC_CONTROL, oss.str());
      // we also need to add the computation time, which is the max
      // time from all threads plus the compute time at this level
      double mt = 0;
      for(size_t i=0;i<texec.size();++i) {
        double lt = texec[i]->get_time();
        if (lt > mt)
          mt = lt;
      }
      timeAccumulator new_comp = timer_comp;
      new_comp.addTime(mt, seq.size());
      data_collector.accumulateTime(new_comp, 
                                    EXEC_COMPUTATION, oss.str());
      // collect scheduling memory footprint size
      data_collector.accumulateSchedMemory(oss.str(),tot_mem);
      // DEBUG (report individual thread timing in debug file
      Loci::debugout << "--------" << rule_name << std::endl;
      for(size_t i=0;i<texec.size();++i) {
        Loci::debugout << "thread " << i << " time: "
                       << texec[i]->get_time() << std::endl;
      }
      ////////
      // collect cache misses
      long_long total_l1_dcm = 0;
      long_long total_l2_dcm = 0;
      for(size_t i=0;i<texec.size();++i) {
        total_l1_dcm += texec[i]->get_l1_dcm();
        total_l2_dcm += texec[i]->get_l2_dcm();
      }
      data_collector.accumulateDCM(oss.str(),total_l1_dcm,total_l2_dcm);
      // report individual thread cache misses in debug file
      Loci::debugout << "********" << rule_name << std::endl;
      for(size_t i=0;i<texec.size();++i) {
        long_long l1 = texec[i]->get_l1_dcm();
        long_long l2 = texec[i]->get_l2_dcm();
        Loci::debugout << "thread " << i << " L1 DCM: " << l1
                       << ", L2 DCM: " << l2 << std::endl;
      }
    }
  protected:
    timeAccumulator timer_comp; // computatation timer
    timeAccumulator timer_ctrl; // control timer
    rule rule_name;
    sequence seq;
    std::vector<ThreadedEMP> texec;
    double tot_mem;             // tabulating the total memory
    executeP exec_rule; // the serial execution module
  };

  // a type used in threaded global reduction
  struct ParamReductionSignal {
    bool done;
    pthread_spinlock_t done_lock;
  };

  // threaded execution module for global reduction computation
  class ExecuteThreaded_param_reduction: public ThreadedExecutionModule {
  public:
    ExecuteThreaded_param_reduction
    (int i, const sequence& s, executeP er,
     CPTR<joiner> jop, const std::vector<int>& rp,
     std::vector<ParamReductionSignal>& ds,
     std::vector<storeRepP>& pr, storeRepP f)
      :tid(i),seq(s),exec_size(s.size()),exec_rule(er),
       join_op(jop),partners(rp),done_signals(ds),partial(pr),
       final_s(f) {}
    
    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      // @PROBLEM:
      // we will need to revist this part later to determine an overall
      // better strategy to control the hidden thread sync cost that is
      // buried in the container designs.
      //
      // it is important that we should copy the unit (initial) value
      // to the local partial result before executing the rule
      //
      // @PROBLEM:
      // it is currently NOT safe to call the "domain()" function on
      // a globally shared container in the thread function. the reason
      // is that the entitySet that represents the domain that gets
      // returned will be copied, which is currently not thread safe.
      // the copying of entitySet in the current implementation relies
      // on an underlying "handle" to share data via a counted pointer.
      // and this reference counted pointer is not protected for thread
      // access.  we do not want to protect it though since that will
      // possibly severely degrade thread performance.
      //
      // another potential source of error lies between the conversion
      // between sequence and entitySet.  currently converting to
      // an entitySet from a sequence is thread safe, however, converting
      // to a sequence from an entitySet is NOT thread safe.  the reason
      // still relates to the way how the internal data sharing is
      // implemented in those data types.
      //
      // so global entitySet and sequence in any current thread 
      // functions should be used with care!!!
      partial[tid]->fast_copy(final_s, EMPTY);//EMPTY);//final_s->domain());//entitySet(seq));
      //partial[tid]->copy(final_s, final_s->domain());//entitySet(seq));
      // we will then execute the rule
      exec_rule->execute_kernel(seq);
      // then we will need to perform the reduction in parallel.
      // we will just need to reverse the termination tree order
      // to perform the reduction
      for(size_t i=0;i<partners.size();++i) {
        int p = partners[i];
        bool done = false;
        while(!done) {
          pthread_spin_lock(&(done_signals[p].done_lock));
          done = done_signals[p].done;
          pthread_spin_unlock(&(done_signals[p].done_lock));
        }
        // we will then reset the done flag in the child
        // we don't need to acquire the lock at this time since
        // this child should have finished by now.
        // the last finished thread's done signal is reset
        // by the main thread when the final result is read
        done_signals[p].done = false;
        // perform the join operation
        join_op->SetArgs(partial[tid], partial[p]);
        join_op->Join(seq);
      }
      // after all join is done, we will then write the local done signal
      pthread_spin_lock(&(done_signals[tid].done_lock));
      done_signals[tid].done = true;
      pthread_spin_unlock(&(done_signals[tid].done_lock));
      timer.addTime(s.stop(),exec_size);
    }
    void print_seq(std::ostream& s) const
    {
      if(verbose || seq.num_intervals() < 4) {
        s << seq;
      } else {
        s << "{" << seq[0] << "...}";
      }
    }
    double get_time() const { return timer.getTime(); }
    long_long get_l1_dcm() const { return 0; } // not implemented yet
    long_long get_l2_dcm() const { return 0; }
  private:
    // thread id
    int tid;
    // the execution sequence (for printing purpose only)
    sequence seq;
    size_t exec_size;
    // we will just delegate the main computation to the underlying
    // concrete rule.
    executeP exec_rule;
    // the join operator
    CPTR<joiner> join_op;
    // the reduction partners of this thread
    const std::vector<int>& partners;
    // this is the done signal pointer
    std::vector<ParamReductionSignal>& done_signals;
    // pointer to the partial results of every thread
    std::vector<storeRepP>& partial;
    // this is the pointer to the store that the final result should go
    // here it is used to copy the initial value
    storeRepP final_s;
    timeAccumulator timer;
  };  

  // this is a module that hooks up the compiler and the
  // threaded execution module for param reduction rule defined above
  class Threaded_execute_param_reduction: public execute_modules {
  public:
    Threaded_execute_param_reduction(rule r, const sequence& s,
                        // target is the name of the final result
                                     const variable& target,
                                     fact_db& facts, sched_db& scheds)
      :rule_name(r), seq(s)
    {
      exec_rule = new execute_rule(r, s, facts, scheds);
      int tnum = thread_control->num_threads();
      done_signals.resize(tnum);
      for(int i=0;i<tnum;++i) {
        done_signals[i].done = false;
        if(pthread_spin_init(&(done_signals[i].done_lock), 0)) {
          // need better error handling later...
          Loci::Abort();
        }
      }

      final_s = facts.get_variable(target);
      std::vector<sequence> partition =
        thread_control->partition_seq(seq);

      rule_implP ar = rule_name.get_info().rule_impl;
      CPTR<joiner> jop = ar->get_joiner();
      
      for(int i=0;i<tnum;++i) {
        storeRepP sp = final_s->new_store(EMPTY);
        //s->allocate(entitySet(partition[i]));
        sp->allocate(final_s->domain());
        partial.push_back(sp);
      }

      tot_mem = 0;
      for(int i=0;i<tnum;++i) {
        texec.push_back
          (new ExecuteThreaded_param_reduction
           (i, partition[i], //exec_rule,
            new execute_rule(rule_name,partition[i],
                             facts,target,partial[i],scheds),
            jop->clone(),thread_control->get_reduction_partners(i),
            done_signals,partial,final_s));
        tot_mem += 16*partition[i].num_intervals();
      }
      tot_mem += 8*seq.num_intervals();

      reduction_root = thread_control->get_reduction_root();
    }

    ~Threaded_execute_param_reduction()
    {
      for(size_t i=0;i<done_signals.size();++i)
        pthread_spin_destroy(&(done_signals[i].done_lock));
      // destroy the local storage
      for(size_t i=0;i<partial.size();++i)
        partial[i]->allocate(EMPTY);
    }
    
    void execute(fact_db& facts, sched_db& scheds)
    {
      current_rule_id = rule_name.ident();
      stopWatch s;
      s.start();

      exec_rule->execute_prelude(seq);
      timer_comp.addTime(s.stop(), seq.size());

      s.start();
      thread_control->setup_facts_scheds(facts,scheds);
      thread_control->restart(texec);
      timer_ctrl.addTime(s.stop(), seq.size());

      thread_control->wait_threads();
      // there are two things to do upon resuming the main thread:
      // 1) reset the done signal of the first work thread so that
      //    it can be reused the next time. (this is not really needed
      //    anyway since no one is going to examine the reduction root's
      //    done flag)
      // 2) store the reduction result into the final place
      s.start();
      done_signals[reduction_root].done = false;
      final_s->copy(partial[reduction_root], final_s->domain());
      // finally execute the postlude
      exec_rule->execute_postlude(seq);
      timer_comp.addTime(s.stop(), seq.size());
      current_rule_id = 0;
    }

    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << rule_name << " over sequence ";
      if(verbose || seq.num_intervals() < 4)
        s << seq;
      else
        s << "{" << seq[0] << "...}";
      s << " using " << texec.size() << " threads";
      if(verbose) {
        s << " (";
        s << "thread " << 0 << ": ";
        texec[0]->print_seq(s);
        for(size_t i=1;i<texec.size();++i) {
          s << ", thread " << i << ": " ;
          texec[i]->print_seq(s);
        }
        s << ")";
      }
      s << std::endl;
    }

    std::string getName() { return "Threaded_execute_param_reduction"; }

    void dataCollate(collectData& data_collector) const
    {
      std::ostringstream oss;
      oss << "rule: " << rule_name;
      // we will need to collect the control time, which is
      // the control time at this level and the maximum control
      // time spent in any thread
      data_collector.accumulateTime(timer_ctrl, EXEC_CONTROL, oss.str());

      double mt = 0;
      for(size_t i=0;i<texec.size();++i) {
        double lt = texec[i]->get_time();
        if (lt > mt)
          mt = lt;
      }
      timeAccumulator new_comp = timer_comp;
      new_comp.addTime(mt, seq.size());
      data_collector.accumulateTime(new_comp,
                                    EXEC_COMPUTATION, oss.str());
      //
      data_collector.accumulateSchedMemory(oss.str(), tot_mem);
    }
  protected:
    timeAccumulator timer_comp;
    timeAccumulator timer_ctrl;
    rule rule_name;
    sequence seq;
    std::vector<ThreadedEMP> texec;
    // here is a vector of partial results for the threads to fill
    // after the thread excution, partial[root] will contain the final result
    std::vector<storeRepP> partial;
    // a vector of reduction signals for the work threads to use
    std::vector<ParamReductionSignal> done_signals;
    // the store in the fact database that the final result will go
    storeRepP final_s;
    // the root thread id of the reduction tree
    int reduction_root;
    double tot_mem;             // tabulate the scheduler memory
    executeP exec_rule;
  };

  extern int chomping_size;
  // this is the module for threading chomped rules with no apply rules
  // inside.  in order to save memory space, all chomping execution
  // sequence for each iteration is calculated on the fly.  this execution
  // module only takes in the number of entities to be computed per
  // chomp iteration.  Another point is that each execution module
  // needs to allocate local storage for all chomped variables and then
  // associate them with the concrete rules that produce them!
  class ExecuteThreaded_chomp: public ThreadedExecutionModule {
  public:
    ExecuteThreaded_chomp(int i, const entitySet* s, const size_t* ces,
       const std::vector<rule>& crs, const std::deque<entitySet>& rss,
       const variableSet& cvr, fact_db& facts) 
      :tid(i), seq(s), chomp_entity_size(ces), chomp_rules(&crs),
       rule_seqs(&rss), chomp_vars(&cvr)
    {
      // first of all we will need to generate local storage
      // for all of the chomping variables.
      for(variableSet::const_iterator
          vi=cvr.begin();vi!=cvr.end();++vi) {
        storeRepP rep = facts.get_variable(*vi);
        // make a local replication
        // allocation will happen at the runtime
        storeRepP new_rep = rep->new_store(EMPTY);
        // add it to the storage map
        chomp_storage[*vi] = new_rep;
      }
      // next we will create local rule impls for execution 
      // and also set up the storage association within the fact_db
      // as well as with the local chomping storage
      for(std::vector<rule>::const_iterator
          ri=crs.begin();ri!=crs.end();++ri) {
        rule_implP rimp = ri->get_rule_implP(); 
        // strangely, we do not need make a copy of the rule impl
        // something like: new_impl = rimp->new_rule_implP()
        // this is because the call to "get_rule_implP()" already
        // makes a copy inside and returns the new copy.  this is
        // not consistent with the variable storage interface, where
        // you have to specifically request a copy to be made
        chomp_impls.push_back(rimp);
        // connect the impl with the fact_db first
        rimp->initialize(facts);
        // then replace those chomped variables with the local storage
        variableSet all_vars;
        all_vars += ri->sources();
        all_vars += ri->targets();
        for(variableSet::const_iterator
            vi=all_vars.begin();vi!=all_vars.end();++vi) {
          // see if the variable is a chomped one
          std::map<variable,storeRepP>::const_iterator f
            = chomp_storage.find(*vi);
          if(f != chomp_storage.end()) {
            // connect to the local storage instead
            rimp->set_store(*vi, f->second); 
            // also sets up the association
            var_in_rules[*vi].push_back(rimp);
          }
        }
      }
      // all set for now...
    }

    // this method is unique for the threaded-chomping execution module.
    // what it does is that it prepares the local chomping storage to
    // be setup properly.  

    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      // first of all, it is IMPORTANT that we run the following 
      // every time before a thread executes the main chomping
      // computation: we will need to regenerate the local storage
      // for the chomping by copying the storage type again from
      // the corresponding fact database.  this is because at runtime,
      // it is possible for the variable storage type to be modified,
      // e.g., through the "prelude" function that could change the
      // vector size of a store_vec, etc.  normally everytime this 
      // happens, it is the same.  and indeed this is the assumption
      // of the older code that only does this once at first execution
      // path.  however, since Loci rules are not pure and side-effects
      // free, it is possible that each rule execution could change 
      // something differently.  therefore to correctly handle this,
      // we will have to reevaulate this every time we run the code.
      // after that, we will also reinitialize the rule impls to
      // reconnect to the new storage.
      for(std::map<variable,storeRepP>::iterator
          mi=chomp_storage.begin();mi!=chomp_storage.end();++mi) {
        storeRepP rep = facts.get_variable(mi->first);
        mi->second = rep->new_store(EMPTY);
        std::vector<rule_implP> rv = var_in_rules[mi->first];
        for(std::vector<rule_implP>::iterator
            ri=rv.begin();ri!=rv.end();++ri)
          (*ri)->set_store(mi->first, mi->second);
      }
      // we will make a local copy of the the rule_seqs restricted
      // only to the allowed total domain for this chomping block,
      // doing this will increase the local storage requirement,
      // but will save time in later entitySet manipulation.
      std::vector<entitySet> restricted_rule_seqs(rule_seqs->size());
      for(size_t i=0;i<rule_seqs->size();++i) {
        restricted_rule_seqs[i] = ((*rule_seqs)[i]) & (*seq);
      }
      // to execute the chomp, we will need to dynamically compute
      // the chomping iteration and the exact execution sequence
      // for each of the chomping rule.
      // 
      // This step requires a bit of explanation:
      //
      // 1) first we will show that all the chomping rules involved
      // must have the same execution sequence and all the chomping
      // variables must have the exact same domain.  the reason for this
      // is because all the chomping variables are strictly internal
      // to the chomping graph, their sole purpose in the program is
      // to compute the targets of the chomping graph.  they are not
      // visible outside of the chomping graph.  also because of that
      // no maps are involved in accessing (either writing or reading)
      // to these chomping variables, together with the existential
      // analysis phases (e.g., the pruning phase), these chomping
      // variables must exist over the same entity domain.
      //
      // 2) to compute by chomping means that we subdivide the domains
      // of these chomping variables and compute each one of them for
      // the whole chomping graph.  here the key step is how to subdivide
      // the domain.  in the older code (in "comp_chomp.cc"), we just
      // work from the smallest entity and each step working on one
      // interval whose size equals the precomputed "chomp size".
      // e.g., if we have the initial domain: [0,101] and our chomp size
      // is 5, then we will break up the domain into: [0,4], [5,9],
      // [10,14],......,[95,99],[100,101].  a total of 21 chomps.
      // however the initial domain may not be contiguous, it is possible
      // we have something like: [0,10],[33,45],[88,123].  then the 
      // older way of partitioning would "waste" the compute power by
      // generating partitions like: [0,4],[5,9],[10,14],[15,19],[20,24]...
      // as some of these chomps corresponds to empty computations.
      //
      // the main reason why the older code is setup like so is because
      // after each chomp, we want to "shift" the domain allocation
      // of each chomping variables so they are ready for the next
      // chomping iteration.  there cannot be gaps inbetween otherwise
      // this will destroy the cache benefit.
      //
      // the ideal partition for chomping for the previous domain will be
      // like this:
      //   1st chomp: [0,4]
      //   2nd chomp: [5,9]
      //   3rd chomp: [10,10],[33,36],
      //   4th chomp: [37,41]
      //   5th chomp: [42,45],[88,88]
      //   6th chomp: [89,93]
      //   ......
      // in this chomp scheme, each iteration will compute exactly 5
      // entities and will not waste time on empty domains.  however
      // the allocation shift for chomping variables will not work
      // when there are gaps.  for example, from allocation [0,4] to
      // [5,9], we can shift the base allocation by 5, but from [5,9]
      // to "[10,10],[33,36]" there is no simple shifting and right 
      // now, static containers must be allocated contiguously.
      // if we change the allocation to [10,36] to accommodate the
      // "[10,10],[33,36]" chomp, then we would have wasted space 
      // as well as allocating perhaps a too large chunk of memory
      // for the cache benefit to take place.
      // note however that this case would probably not occur normally
      // as the static containers have usually been optimized for
      // contiguous memory allocation.  but this could occur at least
      // theoretically and we should deal with it.
      //
      // the method I use here is to remove the gaps by allowing smaller
      // than the planned chomp size.  for example, the previous 
      // initial domain could be divided like this:
      //   1st chomp: [0,4]
      //   2nd chomp: [5,9]
      //   3rd chomp: [10,10]
      //   4th chomp: [33,37]
      //   5th chomp: [38,42]
      //   6th chomp: [43,45]
      //   7th chomp: [88,92]
      //   ......
      // in this scheme, when we reach a gap while not meeting our chomp
      // size, we will just accept the current subdivision,  it may be 
      // smaller than the planned size, but being smaller should not
      // degrade the performance (we may not be fully utilizing the
      // cache capacity though in this case).  the benefit is that
      // we can totally commit to useful computations.  another note is
      // that our shifting offset is no longer fixed in this case.
      // for example, from the 1st chomp to 2nd chomp, we will shift
      // by size 5, but from the 2nd chomp to the 3rd chomp, we need to
      // shift the domain by an offset 23.
      //
      // note: if the gap(s) are not large enough to break the entire
      // chomping size, then we can span the chomping domain over
      // multiple intervals.  for example, we have the entity set:
      // [0,4],[8,20],[50,100] and our chomping size is 10, then the
      // first chomping intervals will be: [0,4],[8,12].  We don't need
      // to stop at the first gap, because the first gap is not wide
      // enough to "consume" all the remaining chomping size.

      // we use several pointers to walk through the initial compute
      // sequence and subdivide it into chomps.
      //
      // this records the total intervals in the sequence
      int seq_size = seq->num_intervals();
      // this points to the current interval index in the sequence
      int cur_seg = 0;
      // this points to the current entity inside the interval
      int cur_ent;
      bool keep_chomping = (seq_size>0);
      if(keep_chomping)
        cur_ent = (*seq)[0].first; // initialize the first entity

      // before we start, we will just allocate the chomping variables
      // first, initially we just wanted to allocate over a fake domain
      // to just get started, later we will shift the domains to
      // the correct ones.
      int pre_dom_start = 0; // the start of previous allocation
      entitySet fake_dom = entitySet(interval(0,*chomp_entity_size-1));
      for(std::map<variable,storeRepP>::iterator
          mi=chomp_storage.begin();mi!=chomp_storage.end();++mi)
        mi->second->allocate(fake_dom);

      while(keep_chomping) {
        int cur_end = cur_ent+*chomp_entity_size-1;
        interval chomp_ivl = interval(cur_ent,cur_end);
        while(cur_seg<seq_size) {
          if(cur_end < (*seq)[cur_seg].first) {
            cur_ent = (*seq)[cur_seg].first;
            break;
          } else if(cur_end < (*seq)[cur_seg].second) {
            cur_ent = cur_end+1;
            break;
          } else {
            ++cur_seg;
          }
        }
        if(cur_seg >= seq_size)
          keep_chomping = false;
        
        // then we will shift the domains of the variables
        int shift = chomp_ivl.first - pre_dom_start;
        pre_dom_start = chomp_ivl.first;
        for(std::map<variable,storeRepP>::iterator
            mi=chomp_storage.begin();mi!=chomp_storage.end();++mi)
          mi->second->shift(shift);
        
        // @PROBLEM:
        // we might want to revist the entitySet/sequence types used
        // in the rules execution.  One problem that stands out is the
        // inconvenience of having to convert between sequence and
        // entitySet in different situations just because certain 
        // operations are not compatible between these two types, or
        // even within the same type (e.g., "-" operators on sequence
        // is not defined, etc).  I think the current rule's execution
        // context is of type "sequence" might not be the most
        // appropriate.  the semantics of most rules say that the 
        // execution order doesn't matter, then why not just use the type
        // "entitySet" as the type for the execution context? there 
        // may be rules whose execution order over entities matter,
        // such as the recursive rules perhaps.  however the advertised
        // semantics of the Loci system is that the rules are 
        // transparent, which includes the execution order over the
        // entities.  so far, many aspects in a Loci rule would break
        // the transparency and as a result, the implementation also
        // has to consider many of these cases, which complicates the
        // matter greatly.  for this particular issue, at least from
        // the threads point of view, a rule's execution context should
        // not care about the orders, and I would prefer this to be
        // reflected in the type of the execution context.  incidentally
        // this change would make many conversion problems go away.
        //
        // we can then compute the rules
        // sequence chomp_seq = sequence(chomp_ivl);
        //
        // NOTE: we could compute the "exec" everytime from scratch by
        // sequence(chomp_ival & (*rule_seqs)[rs] & (*seq)). 
        // in this way we don't have to store anything.  here instead
        // we have chosen to recompute the (*rule_seqs)[rs]&(*seq)
        // part and store that in "restricted_rule_seqs" variable so
        // we can save some time.  when we have many chomp blocks, this
        // could cost us a lot of storage space.  overall we should
        // revisit these issues in the future if this approach (the
        // multi-level partition) is still being used.  together with
        // the the entitySet/sequence issue mentioned above, we should
        // carefully tune the performance and memory cost.
        size_t rs=0;
        for(std::vector<rule_implP>::iterator
            ri=chomp_impls.begin();ri!=chomp_impls.end();++ri,++rs) {
          sequence exec = 
            sequence(chomp_ivl & restricted_rule_seqs[rs]);
          //((*rule_seqs)[rs]) & (*seq));
          (*ri)->compute(exec);
        }
      }

      // finally when we are done, we will deallocate the chomp variables
      for(std::map<variable,storeRepP>::iterator
          mi=chomp_storage.begin();mi!=chomp_storage.end();++mi)
        mi->second->allocate(EMPTY);
    }

    void print_seq(std::ostream& s) const
    {
      if(verbose || seq->num_intervals() < 4) {
        s << *seq;
      } else {
        s << "{" << (*seq)[0] << "...}";
      }
    }
    double get_time() const { return timer.getTime(); }
    long_long get_l1_dcm() const { return 0; }
    long_long get_l2_dcm() const { return 0; }
  private:
    // thread ID
    int tid;
    // the total execution sequence for this execution sequence.
    const entitySet* seq;
    // the number of entities to compute in each chomp iteration
    const size_t* chomp_entity_size;
    // all the rules in the chomping graph (this is taken as a pointer
    // to the control module to save some space)
    const std::vector<rule>* chomp_rules;
    // the execution sequence of each rule (provided by the control
    // module to save space)
    const std::deque<entitySet>* rule_seqs;
    // the concrete rule impl for each chomped rule
    std::vector<rule_implP> chomp_impls;
    // the set of variables in chomping
    const variableSet* chomp_vars;
    // the local storage for chomping
    // this maps each variable to the storage location
    std::map<variable,storeRepP> chomp_storage;
    // this is mapping that tells what chomping variables each rule uses
    std::map<variable,std::vector<rule_implP> > var_in_rules;
    // a timer
    timeAccumulator timer;
  };

  // this class deals with blocked chomping execution.  Each thread
  // creates one of these objects and what this object does is to
  // create for each block a real chomping execution object and then
  // schedule them according to the conflict lock status.
  class ExecuteThreaded_block_chomp: public ThreadedExecutionModule {
  private:
    // a small type used for handling block execution
    struct ChompBlock {
      int bid;  // the block id
      ThreadedEMP exec; // the chomp execution module
      entitySet seq; // the execution sequence of the module
    };
 public:
    ExecuteThreaded_block_chomp(int i, const entitySet& s,
        const std::vector<int>& bs, ReductionLock* l,
        const size_t* ces, const std::vector<rule>& crs,
        const std::deque<entitySet>& rss, const variableSet& cvr,
        fact_db& facts)
      :tid(i), seq(s), exec_size(s.size()), lock(l)
    {
      // the main task in the contructor is to create the real
      // underlying chomping execution modules.
      blocks.resize(bs.size());
      size_t k = 0;
      for(std::list<ChompBlock>::iterator
          li=blocks.begin();li!=blocks.end();++li,++k) {
        li->bid = bs[k];
        li->exec = new ExecuteThreaded_chomp
          (tid, &(li->seq), ces, crs, rss, cvr, facts);
      }
    } 

    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();

      while(!blocks.empty()) {
        std::list<ChompBlock>::iterator candidate = blocks.begin();
        while(candidate != blocks.end()) {
          int bid = candidate->bid;
          if(lock->try_lock(tid, bid)) {
            candidate->seq =
              thread_control->refine_seq_by_block(tid,bid,seq);
            candidate->exec->execute(facts, scheds);
            lock->unlock(tid, bid);
            // NOTE: it is important to reset the sequence back
            // to an empty one after the computation.  the entire
            // purpose of dynamically computing the block sequence
            // is to avoid storing the sequence permanently, thus
            // saving storage space.  if we do not reset the 
            // sequence back to empty at here, then we will end up
            // storing all the sequences generated for all blocks...
            candidate->seq = EMPTY;
            std::list<ChompBlock>::iterator r = candidate;
            ++candidate;
            done_blocks.splice(done_blocks.begin(), blocks, r);
          } else {
            ++candidate;
          }
        }
      }

      blocks.swap(done_blocks);
      timer.addTime(s.stop(), exec_size);
    }

    void print_seq(std::ostream& s) const
    {
      if(verbose || seq.num_intervals() < 4)
        s << seq;
      else 
        s << "{" << seq[0] << "...}";
    }

    double get_time() const {return timer.getTime(); }
    long_long get_l1_dcm() const {return 0;}
    long_long get_l2_dcm() const {return 0;}
  private:
    int tid; // thread id
    const entitySet& seq;
    size_t exec_size;
    std::list<ChompBlock> blocks; // the execution blocks
    std::list<ChompBlock> done_blocks;  // the finished blocks
    ReductionLock* lock; // points to the conflict lock
    timeAccumulator timer;
  };

  // The threaded execution module for chomping rules is the most
  // complicated one as it requires rewrite instead of delegating the
  // work to an underlying normal execution module (as in the pointwise
  // rule).  Several problems require this rewrite: 1) the chomped 
  // variables have to be allocated separately on each thread, also the
  // association of these chomped variables with the respective rule
  // executions needs to be established. We cannot simply reuse the
  // sequential execution module for this purpose.  2) in the normal
  // sequential chomping execution module, we precomputed the execution
  // sequence for each chomping iteration for all the rules and store
  // them in a table.  We could do the same thing here, however since
  // we will be running with potentially large number of threads, such 
  // a design can lead to storage problems as the execution sequences
  // for all chomping iterations on all threads can lead to a large
  // memory footprint.  Instead here we will adopt the same strategy
  // as used in the blocking execution.  We will not precompute and
  // store the sequence for all chomping iterations.  When a thread 
  // executes a chomping iteration, it will then compute the sequence
  // on the fly.  There is a small performance penalty, but it shouldn't
  // be too much.  However this would save us from a potentially huge
  // memory problem.
  class Threaded_execute_chomp: public execute_modules {
  public:
    Threaded_execute_chomp
    (const sequence& s,
     const std::vector<std::pair<rule,rule_compilerP> >& comp,
     const std::deque<entitySet>& rseq,
     const variableSet& cv,fact_db& facts,sched_db& scheds)
      :rule_seqs(rseq),chomp_vars(cv),seq(s),chomp_entity_size(1)
    {
      seq_size = s.size();
      for(std::vector<std::pair<rule,rule_compilerP> >::const_iterator
          vi=comp.begin();vi!=comp.end();++vi) {
        chomp_rules.push_back(vi->first);
        rule_implP rimp = (vi->first).get_rule_implP();
        rimp->initialize(facts);
        chomp_rule_impls.push_back(rimp);
      }
      
      for(std::deque<entitySet>::const_iterator
          di=rseq.begin();di!=rseq.end();++di) {
        rule_seq_seq.push_back(sequence(*di));
      }

      // setting the default data cache size
      D_cache_size = chomping_size * 1024;

      // we will now create the parallel thread execution modules.
      // we do this by splitting the total chomping execution sequence
      // among the threads.  however if the chomping rules contain
      // an apply rule at the final sequence, then we will need to split
      // the sequences into thread blocks instead.

      // we first check to see if a final apply rule is involved:
      // this records the mapping info for all the apply rules (if any)
      std::vector<vmap_info> map_info;
      for(std::vector<rule_implP>::const_iterator
          vi=chomp_rule_impls.begin();vi!=chomp_rule_impls.end();++vi) {
        if( (*vi)->get_rule_class() == rule_impl::APPLY) {
          std::set<vmap_info>::const_iterator vmi = 
            (*vi)->get_info().targets.begin();
          map_info.push_back(*vmi);
        }
      }
      // we will then need to split the rule sequences.
      // Note: we need to partition the entire collection of
      // sequence (the passed in parameter "s") instead of splitting
      // the individual rule sequences since we want to assign
      // each entity to only one thread
      int threads = thread_control->num_threads();

      std::vector<sequence> pt = thread_control->partition_seq(seq);
      partition.resize(pt.size());
      for(int i=0;i<threads;++i)
        partition[i] = entitySet(pt[i]);

      // debugging code: begin()
      // std::vector<std::vector<int> > blocks(threads);
      // for(int i=0;i<threads;++i)
      //   blocks[i] = thread_control->block_thread_seq(i,partition[i]);
      // lock = thread_control->create_empty_conflict_lock();
      // for(int t=0;t<threads;++t)
      //   texec.push_back
      //     (new ExecuteThreaded_block_chomp
      //      (t,partition[t],blocks[t],&lock,&chomp_entity_size,
      //       chomp_rules,rule_seqs,chomp_vars,facts));
      // debugging code: end()

      if(!map_info.empty()) {
        // we will then need to partition the sequence for
        // each thread into blocks
        std::vector<std::vector<int> > blocks(threads);
        for(int i=0;i<threads;++i)
          blocks[i] = thread_control->block_thread_seq(i,pt[i]);
        // then generate a conflict lock
#ifdef SHOW_LOCK_STAT
        std::cerr << "Generating conflict lock for chomping group"
          << " (total rules: " << chomp_rules.size() 
          << ", total reduction: " << map_info.size() 
          << ")--" << std::endl;
        for(size_t i=0;i<chomp_rules.size();++i)
          std::cerr << "--" << chomp_rules[i] << std::endl;
#endif
        lock = thread_control->create_conflict_lock
          (blocks,map_info,facts,scheds);
        // create blocked chomp modules
        for(int t=0;t<threads;++t)
          texec.push_back
            (new ExecuteThreaded_block_chomp
             (t,partition[t],blocks[t],&lock,&chomp_entity_size,
              chomp_rules,rule_seqs,chomp_vars,facts));
      } else {
        // create threaded execute chomp modules
        for(int t=0;t<threads;++t)
          texec.push_back
            (new ExecuteThreaded_chomp
             (t, &(partition[t]), &chomp_entity_size, chomp_rules,
              rule_seqs, chomp_vars, facts));
      }
    }
    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      // first of all, we will need to run the prelude part of all
      // the chomping rules.  we should do it at here since this
      // is the non-threaded part and by definition of the semantics
      // of the prelude, it is sequential.  also by calling the 
      // preludes of all the chomping rules, we can then proceed to 
      // calculate the accurate size of the chomping since some preludes
      // may change the container size.  a slight concern here is that
      // in chomping, most of the targets of the chomping rules are 
      // not yet allocated, will this cause problem with the preludes?
      // for normal rules, the targets are already allocated when we
      // execute the rule and its prelude, and in the normal non-threaded
      // chomping execution, the preludes are called after the chomped
      // variables are allocated, so in those cases, this is not a 
      // concern.  here we have not yet allocated the chomping variables,
      // they will be allocated by the working threads.  but this
      // should not be a problem because the users are not allowed to
      // directly manipulate containers at will, the only thing allowed
      // on the containers for the users is actually to set the 
      // general container size and this would not cause any problem
      // for our case at here.
      for(size_t i=0;i<chomp_rule_impls.size();++i)
        chomp_rule_impls[i]->prelude(rule_seq_seq[i]);
      // next at runtime is that we need to figure out the
      // exact size of an allocation for a single entity in the chomping.
      // this can only be done during runtime since we will generally
      // only "know" the true size at runtime once everything is 
      // calculated and properly allocated.  
      // Note: in the old program, we will do this only once the first 
      // time we ever execute this module, on subsequent execution, we 
      // will not recalculate the chomping size.  we usually use a flag
      // to indicate whether we are in the first run or not.  
      // however now we wanted to calculate the chomping entity size
      // EVERY time we execute the module because there is no guarantee
      // in Loci that the stores won't change inbetween the executions.
      // this is because Loci rules are not pure and side-effects free.
      // thus we cannot guarantee that each execution of a rule will 
      // have the same result
      {
        // first we will test the total size of all the chomping variables
        // for a single entity, this is the unit allocation size
        entitySet test_alloc = interval(1,1);
        int unit_allocation_size = 0;
        for(variableSet::const_iterator
            vi=chomp_vars.begin();vi!=chomp_vars.end();++vi) {
          storeRepP vrep = facts.get_variable(*vi);
          unit_allocation_size += vrep->pack_size(test_alloc);
        }
        // scale is the ratio of the data cache size and the unit
        // allocation, i.e., how many unit allocation can the data
        // cache hold
        double scale = (double)D_cache_size / (double)unit_allocation_size;
        if(scale <= 1.0)
          chomp_entity_size = 1;
        else
          chomp_entity_size = (size_t)(scale+0.5);
        // safty precautions (see "comp_chomp.cc" for detail)
        if((size_t)seq_size < chomp_entity_size)
          chomp_entity_size = seq_size;
      }
      timer_comp.addTime(s.stop(), seq_size);

      s.start();
      // we can now begin to execute the chomping rules, because
      // chomping is super structure, not a single concrete rule,
      // there is no prelude and postlude to execute in this case.
      // we just launch the threaded modules directly:
      thread_control->setup_facts_scheds(facts,scheds);
      thread_control->restart(texec);
      timer_ctrl.addTime(s.stop(), seq_size);

      thread_control->wait_threads();

      // @PROBLEM:
      // we will call the postlude of each chomping rule as a way to
      // execute the postlude of the entire chomping graph.  this is the
      // same strategy as used in the prelude part.  the only potential
      // concern that it might have (the same as in the prelude part) is
      // that normally these postludes are expected to be called one by
      // one after each corresponding rule is run for their sequence.
      // this *might* create an ordering issue since now we have a different
      // order of execution.  because of Loci rules not truely side
      // effects free, different orders can create different results.
      // however since the Loci rules are defined individually without
      // any ordering information specified (i.e., rules are assembled
      // together later by the scheduler according to different application
      // configuration), this is a unique problem introduced by this
      // particular processing strategy.  if we want to ensure a unique
      // and consistent semantics, then we need to address this at a
      // deeper level.  I am adding a tag here so that we may revisit
      // the issue at a later time.
      s.start();
      for(size_t i=0;i<chomp_rule_impls.size();++i)
        chomp_rule_impls[i]->postlude(rule_seq_seq[i]);
      timer_comp.addTime(s.stop(), seq_size);
    }
    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << "--Start chomping over sequence ";
      if(verbose || seq.num_intervals() < 4)
        s<< seq;
      else
        s << "{" << seq[0] << "...}";
      s << " using " << texec.size() << " threads" << std::endl;
      if(verbose) {
        s << " (";
        s << "thread " << 0 << ": ";
        texec[0]->print_seq(s);
        for(size_t i=1;i<texec.size();++i) {
          s << ", thread " << i << ": " ;
          texec[i]->print_seq(s);
        }
        s << ")" << std::endl;
      }
      printIndent(s);
      s << "--Perform chomping on these rules: " << std::endl;
      for(size_t i=0;i<chomp_rules.size();++i) {
        printIndent(s);
        s << "-- " << chomp_rules[i] << std::endl;
      }
      printIndent(s);
      s << "--End chomping" << std::endl;
    }
    std::string getName() { return "Threaded_execute_chomp"; }

    void dataCollate(collectData& data_collector) const
    { // TO BE REVISED...
      std::ostringstream oss ;
      oss <<"chomp: " ;
      for(size_t i=0;i<chomp_rules.size();++i) 
        oss << "[" << chomp_rules[i] << "]";
      data_collector.accumulateTime(timer_comp,EXEC_COMPUTATION,oss.str());
    }
  private:
    timeAccumulator timer_comp;
    timeAccumulator timer_ctrl;
    // rules executing the chomping
    std::vector<rule> chomp_rules;
    // the impls of the rules for chomping, matching the order of
    // the previous variable
    std::vector<rule_implP> chomp_rule_impls;
    // this is the total execution sequence for each chomping rule
    std::deque<entitySet> rule_seqs;
    std::deque<sequence> rule_seq_seq;
    // all of the variables computed in chomping style
    variableSet chomp_vars;
    // this is the overall execution entities for the entire chomping
    sequence seq;               // union of all rules' sequences
    int seq_size;
    // this is the assumed data cache size
    int D_cache_size;
    // this is the number of entities in each chomp
    size_t chomp_entity_size;
    // this is a record of thread domain partitioning of the initial
    // sequence.  we store it here because of the design of the chomping
    // threading modules require a pointer access to the partitioned
    // sequence.
    std::vector<entitySet> partition;
    std::vector<ThreadedEMP> texec;
    // this is the lock used when in the case that the chomping
    // contains an apply rule at the end of the graph
    ReductionLock lock;
  };
  
  // threaded execution module for local reduction computation
  // this is the computation part
  class ExecuteThreaded_local_reduction: public ThreadedExecutionModule {
  public:
    ExecuteThreaded_local_reduction
      (int i, const sequence& s, executeP er,
       const std::vector<int>& bs, ReductionLock* l)
      :tid(i), seq(s), exec_size(s.size()), exec_rule(er),
       blocks(bs.begin(),bs.end()), lock(l) {}

    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();

      // we will iterate over the blocks and execute them one by one,
      // before a block is executed, we will need to first obtain
      // the lock
      while(!blocks.empty()) {
        std::list<int>::iterator candidate = blocks.begin();
        while(candidate != blocks.end()) {
          // try to execute the candidate block
          int bid = *candidate;
          if(lock->try_lock(tid, bid)) {
            // locking successful, execute the block
            sequence block_seq =
              thread_control->refine_seq_by_block(tid,bid,seq);
            exec_rule->execute_kernel(block_seq);
            lock->unlock(tid,bid);
            // now we will just move the block from the execution list
            // to the done blocks list
            std::list<int>::iterator r = candidate;
            ++candidate;
            done_blocks.splice(done_blocks.begin(), blocks, r);
          } else {
            ++candidate;
          }
        }
      }

      // finally swap the block lists
      blocks.swap(done_blocks);
      timer.addTime(s.stop(), exec_size);
    }

    void print_seq(std::ostream& s) const
    {
      if(verbose || seq.num_intervals() < 4)
        s << seq;
      else
        s << "{" << seq[0] << "...}";
    }

    double get_time() const { return timer.getTime(); }
    long_long get_l1_dcm() const { return 0; }
    long_long get_l2_dcm() const { return 0; }
  private:
    int tid; // thread id
    sequence seq; 
    size_t exec_size;
    executeP exec_rule;
    std::list<int> blocks; // the execution blocks
    // this list is used to store blocks that have been finished
    // in this round of execution. at the end, we will need to
    // swap the done_blocks with the blocks so that the next time,
    // we can start over again.
    std::list<int> done_blocks; 
    ReductionLock* lock;
    timeAccumulator timer;
  };

  class Threaded_execute_local_reduction: public execute_modules {
  public:
    Threaded_execute_local_reduction
      (rule r,rule u,const sequence& s,fact_db& facts,sched_db& scheds)
      :rule_name(r), unit_rule(u), seq(s)
      {
        exec_rule = new execute_rule(r, s, facts, scheds);
        int tnum = thread_control->num_threads();
        // partition the execution sequence according to thread domain first
        std::vector<sequence> partition = 
          thread_control->partition_seq(seq);
        // then partition the sequence for each thread into blocks
        std::vector<std::vector<int> > blocks(tnum);
        for(int i=0;i<tnum;++i)
          blocks[i] = thread_control->block_thread_seq(i,partition[i]);
        // finally we will need to generate a conflict lock to be
        // used during the execution of the reduction.
        const rule_impl::info& rinfo = rule_name.get_info().desc;
        std::set<vmap_info>::const_iterator vmi = rinfo.targets.begin();
        lock = 
          thread_control->create_conflict_lock(blocks,*vmi,facts,scheds);
        // finally create all the thread execution modules
        for(int i=0;i<tnum;++i)
          texec.push_back(new ExecuteThreaded_local_reduction
              (i, partition[i], exec_rule, blocks[i], &lock));
      } 

    void execute(fact_db& facts, sched_db& scheds) 
    {
      current_rule_id = rule_name.ident();
      stopWatch s;
      s.start();
      
      exec_rule->execute_prelude(seq);
      timer_comp.addTime(s.stop(), seq.size());

      s.start();
      thread_control->setup_facts_scheds(facts,scheds);
      thread_control->restart(texec);
      timer_ctrl.addTime(s.stop(), seq.size());

      thread_control->wait_threads();

      s.start();
      exec_rule->execute_postlude(seq);
      timer_comp.addTime(s.stop(), seq.size());
      current_rule_id = 0;
    }

    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << rule_name << " over sequence ";
      if(verbose || seq.num_intervals() < 4)
        s << seq;
      else
        s << "{" << seq[0] << "...}";
      s << " using " << texec.size() << " threads";
      if(verbose) {
        s << " (";
        s << "thread " << 0 << ": ";
        texec[0]->print_seq(s);
        for(size_t i=1;i<texec.size();++i) {
          s << ", thread " << i << ": ";
          texec[i]->print_seq(s);
        }
        s << ")";
        s << std::endl;
      }
    }

    std::string getName() { return "Threaded_execute_local_reduction"; }

    void dataCollate(collectData& data_collector) const
    {
    }
  protected:
    timeAccumulator timer_comp;
    timeAccumulator timer_ctrl;
    rule rule_name;
    rule unit_rule;
    sequence seq; // the total execution sequence
    std::vector<ThreadedEMP> texec; // the threaded execution units
    executeP exec_rule; // the underlying rule impl
    ReductionLock lock; // the central lock used in partial reduction
  };

} // end of namespace Loci

#endif  // end #ifdef PTHREADS

#endif
