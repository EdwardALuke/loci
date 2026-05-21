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
//
// Iterative weighted static load balancing and dynamic loop
//   scheduling rule for Loci
// by RL Carino (RLC@HPC.MsState.Edu)
// Iterative weighted static redone by Ed Luke
//
// Ordinarily, Loci executes
//      rp->compute (sequence (exec_set));
// on each processor; however, the exec_sets may have different
// sizes, or the amount of computations required by the items
// may be non-uniform. Therefore, the processor loads may be
// imbalanced, leading to application performance degradation.
//
// This rule adds a load balancing functionality into Loci,
// hopefully to increase performance when the distribution
// of items among processors result in load imbalance. Two
// strategies are implemented: an iterative weighted static
// load balancing strategy, and a dynamic loop scheduling
// strategy.  Descriptions of these are given in the "NOTES"
// below.
//
//

//#define VERBOSE

#include "sched_tools.h"
#include "distribute.h"
#include "Loci_types.h"
#include <sstream>
using std::ostringstream ;


namespace Loci
{

  // set to 1 to activate debugging statements for:
#define DEBUGOUT 0		// load balancing messages
#define SHOWTIMES 0		// summary timings

  // static methods
#define NLB		0	// no load balancing
#define IWS		1	// iterative weighted static load balancing

  // loop scheduling methods
#define FSC		2	// fixed size chunking (same no. of chunks as FAC)
#define FAC		3	// factoring
#define GSS		4	// guided self scheduling

  // weight/contribution of latest time to average time
#define CONTRIB 0.5

  //  skip the JSTART most expensive chunks when deciding transfers
#define JSTART 1

  // message tags
#define TAG_ACK      1		// handshake acknowledgement
#define TAG_INFO     2		// size of succeeding message
#define TAG_INPUT    3		// input data
#define TAG_OUTPUT   4		// output data
#define TAG_TIMES    5          // timing data

  // info[i][] index; info regarding chunk of items to be migrated
#define MSG_LEN      0		// packing size
#define ITEM_COUNT   1		// number of items
#define LOCAL_START  2		// index in facts/local_facts


  //=======================================================
  // NOTES for dynamic loop scheduling (DLS)
  //=======================================================
  //
  // Load balancing by DLS is based on a foreman-worker strategy.
  // The workers execute loads specified by the foreman, and the
  // foreman, in addition to executing loads, is responsible for
  // determining the loads and for detecting termination. A loop
  // scheduling technique is utilized by the foreman to compute
  // loads. A load is basically a chunk of items, which may be
  // local to a worker, or migrated from another worker. The
  // following "#define"s symbolize the actions that can be taken
  // by a processor.

  //=======================================================
  // #define's for DLS
  //=======================================================

  // actions
#define WAIT4_MSG          1
#define TEST4_MSG          2
#define RETRIEVE_MSG       3
#define WORK_LOCAL         4
#define WORK_REMOTE        5
#define SEND_INPUT         6
#define RECV_OUTPUT        7
#define FILL_REQUEST       8
#define QUIT               9
#define EXEC_LOCAL_PART   10
#define EXEC_REMOTE_PART  11
#define EXEC_REMOTE_WHOLE 12

  // The actions FILL_REQUEST and EXEC_REMOTE_PART are unique to
  // the foreman; EXEC_REMOTE_WHOLE is unique to workers; and
  // the rest of the actions are common to the foreman and to
  // workers.

  // =======================================================
  // Outline of algorithm
  // =======================================================

  // =======================================================
  // Initial loads
  // =======================================================
  // if (i_am_foreMan)
  //   Send initial local work info to workers
  //   Set own local work info, nextAction = WORK_LOCAL
  // else // i_am_worker
  //   Set nextAction = WAIT4_MSG
  //
  // =======================================================
  // Main loop
  // =======================================================
  // while (all termination conditions not met)
  //
  //   switch (nextAction)
  //
  // =======================================================
  //     case WAIT4_MSG
  //       MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, ...)
  //       Set nextAction = RETRIEVE_MSG, lastAction = WAIT4_MSG
  // =======================================================
  //
  // =======================================================
  //     case TEST4_MSG
  //       MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, ...)
  //       if (there was a message) Set nextAction = RETRIEVE_MSG
  //       else Set nextAction = lastAction
  // =======================================================
  //
  // =======================================================
  //     case RETRIEVE_MSG
  //       Set nextAction from message status.MPI_TAG (one of the #define)
  //       Retrieve content of message, which may be
  //         - signal to quit
  //         - info on local work to execute
  //         - info on input items to migrate
  //         - a request for load & timing data (for foreman)
  //         - input items
  //         - output items & timing data
  // =======================================================
  //
  // =======================================================
  //     case QUIT
  //       Set one of the termination flags, nextAction = WAIT4_MSG
  // =======================================================
  //
  // =======================================================
  //     case WORK_LOCAL
  //       Update local work info
  //       Compute probing frequency, when to request next work info
  //       Set nextAction = EXEC_LOCAL_PART
  // =======================================================
  //
  // =======================================================
  //     case WORK_REMOTE
  //       Update remote work info
  //       if (i_am_foreman) Set nextAction = EXEC_REMOTE_PART
  //       else Set nextAction = EXEC_REMOTE_WHOLE
  // =======================================================
  //
  // =======================================================
  //     case SEND_INPUT
  //       Send input according to info retrieved from message
  //       Set nextAction = lastAction
  // =======================================================
  //
  // =======================================================
  //     case RECV_OUTPUT
  //       Decrement counter for outputs to receive
  //       Set nextAction = lastAction
  // =======================================================
  //
  // =======================================================
  //     case EXEC_LOCAL_PART
  //       If this is the last subchunk, then send request for next work info
  //       Execute a subchunk from current local work
  //       Update counters and accumulators
  //       if (current local chunk not finished)
  //         Set nextAction = TEST4_MSG, lastAction = EXEC_LOCAL_PART
  //       else if (i_am_foreman) // foreman finished current local chunk
  //         // pretend foreman sent request (to itself)
  //         Set lastAction = WAIT4_MSG, nextAction = FILL_REQUEST
  //       else Set nextAction = WAIT4_MSG
  // =======================================================
  //
  // =======================================================
  //     case EXEC_REMOTE_PART
  //       Foreman executes subchunk of remote items
  //       if (chunk not finished)
  //         Set nextAction = TEST4_MSG, lastAction = EXEC_REMOTE_PART
  //       else // foreman finished current remote chunk
  //         Return outputs, Set lastAction = WAIT4_MSG
  //         // pretend foreman sent a request to self
  //         Set nextAction = FILL_REQUEST
  // =======================================================
  //
  // =======================================================
  //     case EXEC_REMOTE_WHOLE
  //       Worker executes remote chunk, return the outputs,
  //         requests for new work info
  //       Set nextAction = WAIT4_MSG
  // =======================================================
  //
  // =======================================================
  //     case FILL_REQUEST
  //       Foreman takes one of the following actions,
  //         - inform requesting proc that there are no more items to execute
  //         - send new local work info to requesting proc
  //         - inform a heavily loaded worker to send items to requesting proc
  //       Foreman checks for termination
  // =======================================================
  //
  //   end switch
  // end while
  //
  // MPI_Barrier (MPI_COMM_WORLD)
  //


  //=======================================================
  // additional #define's for IWS
  //=======================================================

  // send request for next chunk when the
  // remaining iterations in current chunk reaches PCT_SEND_REQUEST %
  // subject to the minimum of MIN_SEND_REQUEST iterations
#define PCT_SEND_REQUEST 0.10
#define MIN_SEND_REQUEST 50

  // IF THE REQUEST FOR NEXT CHUNK HAS NOT BEEN SENT,
  // probe for messages after executing
  // PCT_PROBE_FREQ % iterations in current chunk
  // subject to the minimum MIN_PROBE_FREQ iterations
#define PCT_PROBE_FREQ 0.05
#define MIN_PROBE_FREQ 10

  //=======================================================
  // common global variables
  //=======================================================

  unsigned char *buf = 0;	// message buffer
  int buf_size = 0;		// size of message buffer
  
  //=======================================================
  // common routines
  //=======================================================

  void dynamic_schedule_rule::GetBuffer (int size) {
    if (size > buf_size) {
      if (buf != 0)
        delete[]buf;
      buf_size = size + size / 5;
      buf = new unsigned char[buf_size];
    }
  }

  void dynamic_schedule_rule::AllocateLBspace (int nItems) {
    for (variableSet::const_iterator vi = inputs.begin ();
	 vi != inputs.end (); ++vi) {
      storeRepP sp = local_facts1->get_variable (*vi);
      sp->allocate (interval (0, nItems - 1));
    }
    for (variableSet::const_iterator vi = outputs.begin ();
	 vi != outputs.end (); ++vi) {
      storeRepP sp = local_facts1->get_variable (*vi);
      sp->allocate (interval (0, nItems - 1));
    }
  }

  void dynamic_schedule_rule::FreeLBspace () {
    for (variableSet::const_iterator vi = inputs.begin ();
	 vi != inputs.end (); ++vi) {
      storeRepP sp = local_facts1->get_variable (*vi);
      sp->allocate (EMPTY);
    }
    for (variableSet::const_iterator vi = outputs.begin ();
	 vi != outputs.end (); ++vi) {
      storeRepP sp = local_facts1->get_variable (*vi);
      sp->allocate (EMPTY);
    }
  }


  // packing size of a chunk of inputs in facts
  int dynamic_schedule_rule::inputPackSize (int tStart, int tSize) {
    int size;
    size = 0;
    for (variableSet::const_iterator vi = inputs.begin ();
	 vi != inputs.end (); ++vi) {
      storeRepP s_ptr = facts1->get_variable (*vi);
      size += s_ptr->pack_size (interval (tStart, tStart + tSize - 1));
    }
    return size;
  }


  // packing size of a chunk of outputs in local_facts
  int dynamic_schedule_rule::outputPackSize (int tStart, int tSize) {
    int size;
    size = 0;
    for (variableSet::const_iterator vi = outputs.begin ();
	 vi != outputs.end (); ++vi) {
      storeRepP s_ptr = local_facts1->get_variable (*vi);
      size += s_ptr->pack_size (interval (tStart, tStart+tSize - 1));
    }
    return size;
  }

  //
  // rule initializer
  //
  dynamic_schedule_rule::dynamic_schedule_rule (rule fi, entitySet eset,
						fact_db & facts,
						sched_db & scheds,
                                                int method) {
    LBMethod = method ;
    int i2[2];

    // Loci initializations
    rp = fi.get_rule_implP ();
    rule_tag = fi;
    given_exec_set = eset;
    exec_set = given_exec_set ;
    compress_set = false ;
    local_compute1 = rp->new_rule_impl ();
    entitySet in = rule_tag.sources ();
    outputs = rule_tag.targets ();
    if(exec_set.num_intervals() > 1) { // If non-contiguous, make contiguous
      // Copy to a compressed set if exec set not already compressed
      int sz = eset.size() ;
      exec_set = interval(0,sz-1) ;
      compress_set = true ;
      main_comp = rp->new_rule_impl() ;
#ifdef VERBOSE
      Loci::debugout << "using compressed set formulation" << endl
                     << "given_set = " << given_exec_set << endl 
                     << "exec_set = " << exec_set << endl ;
#endif
      
      //Setup local facts input variables(types only no allocation)
      for (variableSet::const_iterator vi = in.begin (); vi != in.end (); ++vi) {
        storeRepP store_ptr = rp->get_store (*vi);
        if (isSTORE(store_ptr)) {
          backup_facts.create_fact (*vi, store_ptr->new_store (EMPTY));
        } else {
          backup_facts.create_fact (*vi, facts.get_variable (*vi));
        }
      }

      //Setup local facts output variables
      for (variableSet::const_iterator vi = outputs.begin ();
           vi != outputs.end (); ++vi) {
        storeRepP store_ptr = rp->get_store (*vi);
        backup_facts.create_fact (*vi, store_ptr->new_store (EMPTY));
      }
    } else {
      main_comp = rp ;
    }

    
    if(LBMethod == IWS) {
      // Setup for iterative weighted static
      // first we need to get the size of the iterate space
      int lsz = exec_set.size() ;
      int gsz = 0 ;
      MPI_Allreduce(&lsz, &gsz, 1, MPI_INT, MPI_SUM,
                    facts.get_comm()) ;
      // Now compute chomp size
      int p = facts.get_comm_size() ;
      // Chunk size
      int CHUNKING_FACTOR=256 ;
      int csz = max(gsz/(p*CHUNKING_FACTOR),10) ;
      int nc = max(lsz/csz-1,0) ;
      int first_chunk = lsz - csz*nc ;
#ifdef VERBOSE
      debugout << "chunk_size = " << csz << ", num chunks = " << nc << endl ;
#endif
      int start = exec_set.Min() ;
      int end = start+first_chunk -1 ;
      iwsSize = csz ;
      chunkInfo tmp ;
      tmp.chunkDef = entitySet(interval(start,end))&exec_set ;
      tmp.chunkTime[0] = 0 ;
      tmp.chunkTime[1] = 0 ;
      chunkData.push_back(tmp) ;
      int last = exec_set.Max()+1 ;
      for(start = start+first_chunk;start < last;start+= csz ) {
        tmp.chunkDef = entitySet(interval(start,start+csz-1)) ;
        chunkData.push_back(tmp) ;
      }
      int nchunks = chunkData.size() ;
      for(int i=0;i<nchunks;++i)
        selfChunks.push_back(i) ;
      numBalances = 0 ;
      numRemoteChunks = 0 ;
    }
    
    //Setup local facts input variables(types only no allocation)
    for (variableSet::const_iterator vi = in.begin (); vi != in.end (); ++vi) {
      storeRepP store_ptr = rp->get_store (*vi);
      if (isSTORE(store_ptr)) {
        inputs += *vi;
        local_facts.create_fact (*vi, store_ptr->new_store (EMPTY));
      } else {
        local_facts.create_fact (*vi, facts.get_variable (*vi));
      }
    }

    //Setup local facts output variables
    for (variableSet::const_iterator vi = outputs.begin ();
         vi != outputs.end (); ++vi) {
      storeRepP store_ptr = rp->get_store (*vi);
      local_facts.create_fact (*vi, store_ptr->new_store (EMPTY));
    }

    //Initialize both functions for remote and local execution.
    local_compute1->initialize (local_facts);
    rp->initialize (facts);
    if(compress_set) {
      main_comp->initialize(backup_facts) ;
    } else {
      main_comp->initialize(facts) ;
    }
      

    // ========================================================================
    // Step 0. Initializations
    // ========================================================================

    // proc count, rank
    nProcs = facts.get_comm_size() ;
    
    myRank = facts.get_comm_rank() ;
    foreMan = 0;

    // first items, item counts
    myFirstItem = exec_set.Min ();
    myItemCount = exec_set.size ();

    //Total work
    i2[0] = myFirstItem;
    i2[1] = myItemCount;

    if(yMapSave.size() == 0) {
      vector<int> tmp(nProcs*2) ;
      yMapSave.swap(tmp) ;
    }
    
    MPI_Allgather (i2, 2, MPI_INT, &yMapSave[0], 2, MPI_INT,
                   facts.get_comm()) ;
    allItems = 0;
    for (int proc = 0; proc < nProcs; proc++)
      allItems = allItems + yMapSave[2 * proc + 1];

    // initialize work times
    if(workTime.size() == 0) {
      vector<double> tmp(nProcs) ;
      workTime.swap(tmp) ;
    }
    if(aveWorkTime.size() == 0) {
      vector<double> tmp(nProcs) ;
      aveWorkTime.swap(tmp) ;
    }
    for (int proc = 0; proc < nProcs; proc++) {
      workTime[proc] = 0.0;
      aveWorkTime[proc] = 0.0;
    }

    // number of calls to execute()
    nCalls = 0;


  }

  dynamic_schedule_rule::~dynamic_schedule_rule () {
    if (buf != 0) {
      delete[]buf;
      buf = 0;
      buf_size = 0;
    }
  }
  
  
  
  //---------------------------------------------------------
  // Reworked iterative weighted static code.
  //
  // 
  void dynamic_schedule_rule::iterative_weighted_static () {
    //=======================================================
    // First perform communication schedule
    //=======================================================
    // Note, we will either be sending or receiving, but not both.
    //    MPI_Barrier(MPI_COMM_WORLD) ;
    stopWatch timer ;
    timer.start() ;
    if(sendChunks.size() != 0) {
      for(size_t i=0;i<sendChunks.size();++i) {
        int buf_size = sendChunks[i].send_size ;
        vector<unsigned char> buf(buf_size) ;
        int position = 0 ;
        for(size_t j=0;j<sendChunks[i].chunkList.size();++j) {
          int ch = sendChunks[i].chunkList[j] ;
          entitySet sendSet = chunkData[ch].chunkDef ;
          for (variableSet::const_iterator vi = inputs.begin ();
               vi != inputs.end (); ++vi) {
              storeRepP s_ptr = facts1->get_variable (*vi);
              s_ptr->pack (&buf[0], position, buf_size,
                           sendSet);
          }
        }
        int dest = sendChunks[i].proc ;
#ifdef VERBOSE
        debugout << "Sending data to " << dest << endl ;
#endif
        MPI_Send (&buf[0], buf_size, MPI_PACKED, dest, TAG_INPUT,
                  facts1->get_comm()) ;
      }
    } else if(recvChunks.size() != 0) {
      // First post Irecvs
      int tot_buf_size = 0 ;
      for(size_t i=0;i<recvChunks.size();++i) {
	tot_buf_size += recvChunks[i].recv_size ;
      }
      vector<MPI_Request> req_list(recvChunks.size()) ;

      vector<unsigned char> buf(tot_buf_size) ;
      int offset = 0 ;
      for(size_t i=0;i<recvChunks.size();++i) {
        int buf_size = recvChunks[i].recv_size ;
        int src = recvChunks[i].proc ;
#ifdef VERBOSE
        debugout << "recieving data from " << src << endl ;
#endif
        MPI_Irecv (&buf[offset], buf_size, MPI_PACKED, src, TAG_INPUT,
                    facts1->get_comm(), &req_list[i]) ;
	offset += buf_size ;
      }

      // Allocate space for the data we will receive
      int nchunks = 0 ;
      for(size_t i=0;i<recvChunks.size();++i) {
        nchunks += recvChunks[i].chunkList.size() ;
      }
      entitySet loc_set = interval(0,nchunks*iwsSize-1) ;
      for (variableSet::const_iterator vi = inputs.begin ();
           vi != inputs.end (); ++vi) {
        storeRepP sp = local_facts1->get_variable (*vi);
        sp->allocate (loc_set);
      }
      for (variableSet::const_iterator vi = outputs.begin ();
           vi != outputs.end (); ++vi) {
        storeRepP sp = local_facts1->get_variable (*vi);
        sp->allocate (loc_set);
      }

      // now wait for receives to complete
      int nreqs = req_list.size() ;
      vector<MPI_Status> status_list(nreqs) ;
      MPI_Waitall(nreqs,&req_list[0],&status_list[0]) ;
      // Now unpack the data
      int cnt = 0 ;
      offset = 0 ;
      for(size_t i=0;i<recvChunks.size();++i) {
        int buf_size = recvChunks[i].recv_size ;
        int position = 0 ;
        for(size_t j=0;j<recvChunks[i].chunkList.size();++j) {
          entitySet recvSet = interval(cnt,cnt+iwsSize-1) ;
          for (variableSet::const_iterator vi = inputs.begin ();
               vi != inputs.end (); ++vi) {
            storeRepP s_ptr = local_facts1->get_variable (*vi);
            s_ptr->unpack (&buf[offset], position, buf_size,sequence(recvSet)) ;
          }
          cnt+= iwsSize ;
        }
	offset += buf_size ;
      }
    }

#ifdef VERBOSE
    double comm_time = timer.stop() ;
#endif
    //=======================================================
    // Execute Chunks recording time
    //=======================================================

    vector<float> remote_times(numRemoteChunks,0) ;

    stopWatch exec_clock ;
    exec_clock.start() ;
    // Compute local work
    for(size_t i = 0;i<selfChunks.size();++i) {
      int ch = selfChunks[i] ;
      stopWatch s ;
      s.start() ;
      sequence seq(chunkData[ch].chunkDef) ;
      main_comp->prelude(seq) ;
      main_comp->compute(seq) ;
      double elapsed_time = s.stop() ;
      chunkData[ch].chunkTime[0] += elapsed_time ;
      comp_timer.addTime(elapsed_time,1) ;
    }
    // Now compute remote work
    int cnt = 0 ;
    int skip = iwsSize ;
    for(int i=0;i<numRemoteChunks;++i) {
      stopWatch s ;
      s.start() ;
      sequence seq(interval(cnt,cnt+skip-1)) ;
      local_comp1->prelude(seq) ;
      local_comp1->compute(seq) ;
      double elapsed_time = s.stop() ;
      remote_times[i] = elapsed_time;
      comp_timer.addTime(elapsed_time,1) ;
      cnt += skip ;
    }
    float exec_time = exec_clock.stop() ;
    //=======================================================
    // Return chunks to owning processor
    //=======================================================
    timer.start() ;
    if(sendChunks.size() != 0) {
      // First post Irecvs
      int tot_buf_size = 0 ;
      for(size_t i=0;i<sendChunks.size();++i) {
	tot_buf_size += sendChunks[i].recv_size ;
      }
      vector<unsigned char> buf(tot_buf_size) ;
      int tot_chunks = 0 ;
      for(size_t i=0;i<sendChunks.size();++i) {
        tot_chunks += sendChunks[i].chunkList.size() ;
      }
      vector<float> times(tot_chunks) ;

      vector<MPI_Request> req_list(sendChunks.size()*2) ;

      int offset1 = 0 ;
      int offset2 = 0 ;

      for(size_t i=0;i<sendChunks.size();++i) {
        int buf_size = sendChunks[i].recv_size ;
        int dest = sendChunks[i].proc ;
#ifdef VERBOSE
        debugout << "recieving computed data from " << dest << endl ;
#endif
        MPI_Irecv (&buf[offset1], buf_size, MPI_PACKED, dest, TAG_OUTPUT,
                    facts1->get_comm(), &req_list[i*2]) ;
	offset1 += buf_size ;
        int nchunks = sendChunks[i].chunkList.size() ;
        MPI_Irecv(&times[offset2], nchunks, MPI_FLOAT, dest, TAG_TIMES,
		  facts1->get_comm(), &req_list[i*2+1]) ;
	offset2 += nchunks ;
      }

      // wait on recvs to complete
      int nreqs = req_list.size() ;
      vector<MPI_Status> status_list(nreqs) ;
      MPI_Waitall(nreqs,&req_list[0],&status_list[0]) ;
      
      offset1 = 0 ;
      offset2 = 0 ;
      for(size_t i=0;i<sendChunks.size();++i) {
        int buf_size = sendChunks[i].recv_size ;
        int position = 0 ;
        for(size_t j=0;j<sendChunks[i].chunkList.size();++j) {
          int ch = sendChunks[i].chunkList[j] ;
          entitySet sendSet = chunkData[ch].chunkDef ;
          for (variableSet::const_iterator vi = outputs.begin ();
               vi != outputs.end (); ++vi) {
              storeRepP s_ptr = facts1->get_variable (*vi);
              s_ptr->unpack (&buf[offset1], position, buf_size,
                             sendSet);
          }
        }
	offset1 += buf_size ;
        // Receive times
        int nchunks = sendChunks[i].chunkList.size() ;
        for(int j=0;j<nchunks;++j) {
          int ch = sendChunks[i].chunkList[j] ;
          chunkData[ch].chunkTime[0] += times[offset2+j] ;
        }
	offset2 += nchunks ;
      }
    } else if(recvChunks.size() != 0) {
      // Now send output data back
      int cnt = 0 ;
      int cnt2 = 0 ;
      for(size_t i=0;i<recvChunks.size();++i) {
        int buf_size = recvChunks[i].send_size ;
        vector<unsigned char> buf(buf_size) ;
        int dest = recvChunks[i].proc ;
        
        int position = 0 ;
        for(size_t j=0;j<recvChunks[i].chunkList.size();++j) {
          entitySet sendSet = interval(cnt,cnt+iwsSize-1) ;
          for (variableSet::const_iterator vi = outputs.begin ();
               vi != outputs.end (); ++vi) {
            storeRepP s_ptr = local_facts1->get_variable (*vi);
            s_ptr->pack (&buf[0], position, buf_size,sendSet) ;
          }
          cnt += iwsSize ;
        }
#ifdef VERBOSE
        debugout << "sending computed data to " << dest << endl ;
#endif
        MPI_Send (&buf[0], buf_size, MPI_PACKED, dest, TAG_OUTPUT,
                  facts1->get_comm()) ;
        
        int nchunks = recvChunks[i].chunkList.size() ;
        MPI_Send(&remote_times[cnt2], nchunks, MPI_FLOAT, dest,
                 TAG_TIMES, facts1->get_comm()) ;
        cnt2 += nchunks ;
      }
      // Release space for the remote computation
      for (variableSet::const_iterator vi = inputs.begin ();
           vi != inputs.end (); ++vi) {
        storeRepP sp = local_facts1->get_variable (*vi);
        sp->allocate (EMPTY);
      }
      for (variableSet::const_iterator vi = outputs.begin ();
           vi != outputs.end (); ++vi) {
        storeRepP sp = local_facts1->get_variable (*vi);
        sp->allocate (EMPTY);
      }
    }
#ifdef VERBOSE
    comm_time += timer.stop() ;
    debugout << "done communication, collect times" << endl ;
#endif
    //=======================================================
    // Estimate parallel efficiency
    //=======================================================
    float total_time = 0 ;
    MPI_Allreduce(&exec_time, &total_time, 1, MPI_FLOAT, MPI_SUM,
                  facts1->get_comm()) ;
    float max_time = 0 ;
    MPI_Allreduce(&exec_time, &max_time, 1, MPI_FLOAT, MPI_MAX,
                  facts1->get_comm()) ;
    double ef = total_time/(max_time*float(get_exec_size())) ;

    if(get_exec_rank() == 0) 
      Loci::debugout << "parallel efficiency of IWS scheduled computation is "
                     << ef*100.0 << "%, time =" << max_time << endl ;
    //=======================================================
    // If efficiency too low, regenerate load balance schedule
    //=======================================================
    timer.start() ;
    if((numBalances == 0 && ef < .5) ||
       (numBalances != 0 && ef < 0.89 && (numBalances&0x3) == 0) ) {
      const float eff_tol = 0.01 ;
      // Compute balanced schedule, first compute total time
      float time = 0 ;
      for(size_t i=0;i<chunkData.size();++i)
	time += chunkData[i].chunkTime[0]+chunkData[i].chunkTime[1] ;
      float total_time = 0 ;
#ifdef VERBOSE
      debugout << "reducing " << time << endl ;
#endif
      MPI_Allreduce(&time, &total_time, 1, MPI_FLOAT, MPI_SUM,
                    facts1->get_comm()) ;
      float mean_time = total_time/float(get_exec_size()) ;
      float diff_time = time-mean_time ;
      if(fabs(diff_time) < mean_time*eff_tol) // If close to average time, zero
	diff_time = 0 ;
      
      // Compute chunks that will be sent to other processors
      vector<int> send_list ;
      if(diff_time > 0) { // If this processor is taking too much time, 
	// allocate chunks that get it to the mean_time ;
	vector<pair<float,int> > chunk_times(chunkData.size()-1) ;
	for(size_t i=1;i<chunkData.size();++i) {
	  chunk_times[i-1].first = 
	    chunkData[i].chunkTime[0]+chunkData[i].chunkTime[1] ;
	  chunk_times[i-1].second = i ;
	}
	std::sort(chunk_times.begin(),chunk_times.end()) ;
        
        
#define SELECT_GIVE
#ifdef SELECT_GIVE
        // select chunks to give away
        float time_give = 0 ;
	for(int kk = chunk_times.size()-1;kk>=0;--kk) {
          if(time_give > diff_time)
            break ;
	  if(time_give+chunk_times[kk].first < diff_time+mean_time*eff_tol) {
	    send_list.push_back(chunk_times[kk].second) ;
            time_give += chunk_times[kk].first ;
          }
	}
#endif
#ifdef SELECT_KEEP
        // This code for this processor to keep what it can to meet mean time
	float time_keep = chunkData[0].chunkTime[0]+chunkData[0].chunkTime[1] ;
        // keep the smallest timed chunks up to keep_percent of mean to reduce
        // number of chunks communicated
        float keep_percent = 0.01 ;
        int stop = 0 ;
        while(time_keep+chunk_times[stop].first < keep_percent*mean_time &&
              stop < int(chunk_times.size())) {
          time_keep += chunk_times[stop].first ;
          stop++ ;
        }
          
	int kk ;
	for(kk = chunk_times.size()-1;kk>=stop;--kk) {
          if(time_keep > mean_time)
            break ;
	  if(time_keep+chunk_times[kk].first < mean_time*(1.+eff_tol)) {
	    time_keep += chunk_times[kk].first ;
	  } else
	    send_list.push_back(chunk_times[kk].second) ;
	}
	for(;kk>=0;--kk)
	  send_list.push_back(chunk_times[kk].second) ;
#endif
      }	
      if(send_list.size() != 0) {
	diff_time = 0 ;
	for(size_t i=0;i<send_list.size();++i)
	  diff_time += (chunkData[send_list[i]].chunkTime[0]+
			chunkData[send_list[i]].chunkTime[1]) ;
      }
      vector<float> time_xfers(get_exec_size()) ;
      MPI_Allgather(&diff_time, 1, MPI_FLOAT, &time_xfers[0], 1,
		    MPI_FLOAT, facts1->get_comm()) ;
#ifdef VERBOSE
      if(get_exec_rank() == 0) {
	debugout << "mean_time = " << mean_time << endl ;
	debugout << "time xfers =";
	for(int i=0;i<get_exec_size();++i)
	  debugout << " " << time_xfers[i] ;
	debugout << endl ;
      }
#endif
      vector<pair<float,int> > send_chunks ;
      vector<pair<float,int> > recv_chunks ;
      for(int i=0;i<get_exec_size();++i) {
	if(time_xfers[i] > 0.0) 
	  send_chunks.push_back(pair<float,int>(time_xfers[i],i)) ;
	if(time_xfers[i] < 0.0) 
	  recv_chunks.push_back(pair<float,int>(-time_xfers[i],i)) ;
      }
      sort(send_chunks.begin(),send_chunks.end()) ;
      sort(recv_chunks.begin(),recv_chunks.end()) ;

      // Compute chunk sendto schedule
      vector<pair<int,vector<int> > > sendto ;
      
      for(int i=recv_chunks.size()-1;i>=0;--i) {
	double ptime = recv_chunks[i].first ;
	for(int j=send_chunks.size()-1;j>=0;--j) {
	  if(send_chunks[j].second >=0) {
	    if(ptime - send_chunks[j].first > -mean_time*eff_tol) {
	      ptime -= send_chunks[j].first ;
	      //assign all chunks from this processor
	      if(send_chunks[j].second == get_exec_rank()) {
#ifdef VERBOSE
                debugout << "adding " << send_list.size() << "chunks to sendto"
                         << endl ;
#endif
                if(send_list.size() > 0) 
                  sendto.push_back(pair<int,vector<int> >(recv_chunks[i].second,
                                                          send_list)) ;
		send_list.clear() ;
	      }
	      time_xfers[send_chunks[j].second] = 0 ;
	      send_chunks[j].second = -1 ;
              if(ptime < 0)
                break ;
	    }
	  }
	}
	ptime = -max(ptime,0.0) ;
	if(fabs(ptime) < mean_time*eff_tol)
	  ptime = 0 ;
	time_xfers[recv_chunks[i].second] = ptime ;
      }
      
#ifdef VERBOSE
      if(get_exec_rank() == 0) {
	debugout << "time xfers2 =";
	for(int i=0;i<get_exec_size();++i)
	  debugout << " " << time_xfers[i] ;
	debugout << endl ;
      }
#endif
      
      bool rebalance = false ;
      for(int i=0;i<get_exec_size();++i)
	if(time_xfers[i] > 0)
	  rebalance = true ;
      bool sources = false ;
      for(int i=0;i<get_exec_size();++i)
	if(time_xfers[i] < 0)
	  sources = true ;
      if(rebalance && sources) { // we have residual work to distribute
	int nchunks = send_list.size() ;
	vector<int> chunk_groups(get_exec_size()) ;
	MPI_Allgather(&nchunks, 1, MPI_INT, &chunk_groups[0], 1,
		      MPI_INT, facts1->get_comm()) ;
	vector<float> chunk_times(nchunks) ;
	for(int i=0;i<nchunks;++i) {
	  int ch = send_list[i] ;
	  chunk_times[i] = (chunkData[ch].chunkTime[0]+
			    chunkData[ch].chunkTime[1]) ;
	}
	vector<int> chunk_displ(get_exec_size()) ;
	chunk_displ[0] = 0 ;
	for(int i=1;i<get_exec_size();++i) 
	  chunk_displ[i] = chunk_displ[i-1]+chunk_groups[i-1] ;

	int ntchunks = (chunk_displ[get_exec_size()-1]+
			chunk_groups[get_exec_size()-1]) ;

	vector<float> chunk_time_gather(ntchunks) ;
	
	MPI_Allgatherv(&chunk_times[0], nchunks, MPI_FLOAT,
		       &chunk_time_gather[0],
		       &chunk_groups[0],
		       &chunk_displ[0],
		       MPI_FLOAT,
		       facts1->get_comm()) ;
	vector<pair<float,pair<int,int> > > chunkQueue(ntchunks) ;
        int cnk = 0 ;
	for(int i=0;i<get_exec_size();++i) {
	  for(int j=0;j<chunk_groups[i];++j) {
	    pair<int,int> chunk_info(i,j) ;
	    float chunk_time = chunk_time_gather[chunk_displ[i]+j] ;
	    chunkQueue[cnk] = pair<float,pair<int,int> >(chunk_time,
                                                         chunk_info) ;
            cnk++ ;
	  }
	}
        cnk = 0 ;
	for(int i=0;i<get_exec_size();++i) {
          if(chunk_groups[i] > 1) {
            std::sort(&chunkQueue[cnk],&chunkQueue[cnk+chunk_groups[i]]) ;
          }
          cnk += chunk_groups[i] ;
        }
            
            //	sort(chunkQueue.begin(),chunkQueue.end()) ;
        
	int chunkQueueStart = chunkQueue.size()-1 ;
	for(int i=get_exec_size()-1;i>=0;--i) 
	  if(time_xfers[i] < 0) {
	    vector<int> sendto_p ;
	    float time_x = -time_xfers[i] ;
	    for(int j=chunkQueueStart;j >=0;--j)
	      if(time_x > 0 && (chunkQueue[j].first > 0) &&
		 (time_x - chunkQueue[j].first > -mean_time*eff_tol)) {
		// assign chunk 
		time_x -= chunkQueue[j].first ;
		const int cp = chunkQueue[j].second.first ;
		if(cp == get_exec_rank())
		  sendto_p.push_back(send_list[chunkQueue[j].second.second]) ;
		time_xfers[cp] -= chunkQueue[j].first ;
                if(time_xfers[cp] < mean_time*eff_tol)
                  time_xfers[cp] = 0 ;
		chunkQueue[j].first = -1.0 ; // remove from consideration
		if(time_x < 0)
		  break ;
	      }
	    //skip deleted entries
	    for(;chunkQueueStart>0;--chunkQueueStart)
	      if(chunkQueue[chunkQueueStart].first > 0)
		break ;
	    if(sendto_p.size() > 0) {
	      sendto.push_back(pair<int,vector<int> >(i,sendto_p)) ;
	    }
	    time_x = -time_x ;
	    if(fabs(time_x) < mean_time*eff_tol)
	      time_x = 0 ;
	    time_xfers[i] = time_x ;
	  }

        
        if(get_exec_rank() == 0) {
          bool found = false ;
          for(int j=chunkQueueStart;j >=0;--j) 
            if((chunkQueue[j].first > 0)) {
              found = true ;
            }
          if(found)  {
            debugout << "chunks remaining in queue:" << endl ;
            for(int j=chunkQueueStart;j >=0;--j) 
              if((chunkQueue[j].first > 0)) {
                debugout << "chunk from p=" << chunkQueue[j].second.first
                         << ", time = " << chunkQueue[j].first << endl ;
              }
          }
        }

               
      }
#ifdef VERBOSE
      if(get_exec_rank() == 0) {
	debugout << "time xfers3 =";
	for(int i=0;i<get_exec_size();++i)
	  debugout << " " << time_xfers[i] ;
	debugout << endl ;
      }
#endif


	
      //------------------------------------------------------------------
      // convert sendto to execution and communication schedules
      //------------------------------------------------------------------
      
#ifdef VERBOSE
      Loci::debugout << "sendto:" << endl ;
      for(size_t i=0;i<sendto.size();++i)
        debugout << sendto[i].second.size() << "chunks to " <<
          sendto[i].first << endl ;
#endif
      int numChunks = chunkData.size() ;
      // Convert sendto to schedule
      vector<int> chunk_list(numChunks,-1) ;
      for(size_t i=0;i<sendto.size();++i)
	for(size_t j=0;j<sendto[i].second.size();++j)
	  chunk_list[sendto[i].second[j]] = sendto[i].first ;
      // Setup local execution list

      vector<int> local_list ;
      for(int i=0;i<numChunks;++i)
	if(chunk_list[i] == -1)
	  local_list.push_back(i) ;
      selfChunks.swap(local_list) ;

      int local_chk_info[2],global_chk_info[2] ;
      local_chk_info[0] = selfChunks.size() ;
      local_chk_info[1] = numChunks ;
      MPI_Allreduce(&local_chk_info[0], &global_chk_info[0], 2,
                    MPI_INT, MPI_SUM, facts1->get_comm()) ;
      if(get_exec_rank() == 0) {
        debugout << "IWS Schedule: Communicating "
                 << global_chk_info[1]-global_chk_info[0]
                 << " chunks, "
                 << 100.0*(1.-double(global_chk_info[0])/double(global_chk_info[1])) << "% of all chunks." << endl ;
      }
         
                           

      vector<chunkCommInfo> local_send_list ;
      for(size_t i=0;i<sendto.size();++i) {
	chunkCommInfo tmp ;
	tmp.proc = sendto[i].first ;
	tmp.chunkList = sendto[i].second ;
	entitySet chunk = interval(0,iwsSize-1) ;
	int szs = 0 ;
	for (variableSet::const_iterator vi = inputs.begin ();
	     vi != inputs.end (); ++vi) {
	  storeRepP s_ptr = facts1->get_variable (*vi);
	  szs += s_ptr->pack_size(chunk);
	}
	tmp.send_size = szs*tmp.chunkList.size() ;
	int szr = 0 ;
	for (variableSet::const_iterator vi = outputs.begin ();
	     vi != outputs.end (); ++vi) {
	  storeRepP s_ptr = facts1->get_variable (*vi);
	  szr += s_ptr->pack_size(chunk);
	}
	tmp.recv_size = szr*tmp.chunkList.size() ;
	local_send_list.push_back(tmp) ;
      }
      sendChunks.swap(local_send_list) ;

      // now invert sendChunks
      vector<int> sendSizes(get_exec_size(),0) ;
      for(size_t i=0;i<sendChunks.size();++i)
	sendSizes[sendChunks[i].proc] = sendChunks[i].chunkList.size() ;
      
      vector<int> recvSizes(get_exec_size(),0) ;
      MPI_Alltoall(&sendSizes[0], 1, MPI_INT,
		   &recvSizes[0], 1, MPI_INT,
		   facts1->get_comm()) ;
      numRemoteChunks = 0 ;

      vector<chunkCommInfo> local_recv_list ;
      for(int i=0;i<get_exec_size();++i) {
	numRemoteChunks += recvSizes[i] ;
	if(recvSizes[i]!=0) {
	  chunkCommInfo tmp ;
	  tmp.proc = i ;
	  for(int k=0;k<recvSizes[i];++k)
	    tmp.chunkList.push_back(k) ;
	  local_recv_list.push_back(tmp) ;
	}
      }
      if(numRemoteChunks != 0 && sendChunks.size() !=0) {
	cerr << "logic error in iterative weighted static LB method" << endl ;
	Loci::Abort() ;
      }
      for(size_t i=0;i<local_recv_list.size();++i) {
        MPI_Status tStatus;
	MPI_Recv(&local_recv_list[i].recv_size, 1, MPI_INT,
		 local_recv_list[i].proc, TAG_INFO,
		 facts1->get_comm(), &tStatus) ;
      }
      for(size_t i=0;i<sendChunks.size();++i) {
	MPI_Send(&sendChunks[i].send_size, 1, MPI_INT,
		 sendChunks[i].proc, TAG_INFO,
		 facts1->get_comm()) ;
      }
      for(size_t i=0;i<local_recv_list.size();++i) {
        MPI_Status tStatus;
	MPI_Recv(&local_recv_list[i].send_size, 1, MPI_INT,
		 local_recv_list[i].proc, TAG_INFO,
		 facts1->get_comm(), &tStatus) ;
      }
      for(size_t i=0;i<sendChunks.size();++i) {
	MPI_Send(&sendChunks[i].recv_size, 1, MPI_INT,
		 sendChunks[i].proc, TAG_INFO,
		 facts1->get_comm()) ;
      }
      recvChunks.swap(local_recv_list) ;
#ifdef VERBOSE
      if(sendChunks.size() != 0) {
      	debugout << "sendChunks: " << endl ;
	for(size_t i=0;i<sendChunks.size();++i)
	  debugout << sendChunks[i].proc << ' ' << sendChunks[i].chunkList.size() << ' ' << sendChunks[i].send_size << ' ' << sendChunks[i].recv_size << endl ;
      }
      if(recvChunks.size() != 0) {
	debugout << "recvChunks: " << endl ;
	for(size_t i=0;i<recvChunks.size();++i)
	  debugout << recvChunks[i].proc << ' ' << recvChunks[i].chunkList.size() << ' ' << recvChunks[i].send_size << ' ' << recvChunks[i].recv_size << endl ;
      }
#endif
      int nsendrecvl = recvChunks.size()+sendChunks.size() ;
      int nsendrecvg = 0 ;
      int nsendrecvs = 0 ;
      MPI_Allreduce(&nsendrecvl, &nsendrecvg, 1, MPI_INT, MPI_MAX,
                    facts1->get_comm()) ;
      MPI_Allreduce(&nsendrecvl, &nsendrecvs, 1, MPI_INT, MPI_SUM,
                    facts1->get_comm()) ;
      if(get_exec_rank() == 0) {
        debugout << "IWS Schedule: Each processor communicating with an average of "
                 << double(nsendrecvs)/double(get_exec_size())
                 << " processors." << endl
                 << "IWS Schedule: Maximum number of communicating partners: "
                 << nsendrecvg << " processors." << endl ;
      }
         
    }
    
#ifdef VERBOSE
    double sched_time = timer.stop() ;
    double timesin[2] = {comm_time,sched_time} ;
    double timesmx[2] = {0,0} ;

    MPI_Allreduce(&timesin[0], &timesmx[0], 2, MPI_DOUBLE, MPI_MAX,
                  facts1->get_comm()) ;
    if(get_exec_rank() == 0)
      Loci::debugout << "comm_time = " << timesmx[0] << ",sched_time=" << timesmx[1] << endl ;
#endif
    //=======================================================
    // If past trigger point, reset timing data
    //=======================================================
    if((numBalances&0x7) == 0) {
      for(size_t i=0;i<chunkData.size();++i) {
        chunkData[i].chunkTime[1] = chunkData[i].chunkTime[0] ;
        chunkData[i].chunkTime[0] = 0 ;
      }
    }

    numBalances++ ;

  }

  //=======================================================
  // routines called by DLS
  //=======================================================


  void dynamic_schedule_rule::SendInfo (int dest, int action, int chunkStart,
                                        int chunkSize, int chunkDest,
                                        double mu) {
    int bufInfo[3], pSize, tPos;
    int size ;
    MPI_Comm comm = facts1->get_comm() ;
    MPI_Pack_size(3, MPI_INT, comm, &size) ;
    pSize = size ;
    MPI_Pack_size(1, MPI_DOUBLE, comm, &size) ;
    pSize += size ;

    GetBuffer (pSize);
    bufInfo[0] = chunkStart;
    bufInfo[1] = chunkSize;
    bufInfo[2] = chunkDest;
    tPos = 0;
    MPI_Pack (bufInfo, 3, MPI_INT, buf, pSize, &tPos, comm) ;
    MPI_Pack (&mu, 1, MPI_DOUBLE, buf, pSize, &tPos, comm) ;
    MPI_Send (buf, pSize, MPI_PACKED, dest, action, comm) ;
#if DEBUGOUT
    //if (myRank == foreMan)
    Loci::debugout << "Info: action=" << action << ", " << myRank << "->"
                   << dest << "->" << chunkDest << " : " << chunkStart
                   << "," << chunkSize << endl;
#endif
  }
  
  void dynamic_schedule_rule:: RecvInfo (int src, int action, int *chunkStart,
                                         int *chunkSize, int *chunkDest,
                                         double *mu) {
    int bufInfo[3], pSize, tPos;
    MPI_Status tStatus;

    MPI_Comm comm = facts1->get_comm() ;
    int size ;
    MPI_Pack_size(3, MPI_INT, comm, &size) ;
    pSize = size ;
    MPI_Pack_size(1, MPI_DOUBLE, comm, &size) ;
    pSize += size ;

    GetBuffer (pSize);
    MPI_Recv (buf, pSize, MPI_PACKED, src, action, comm, &tStatus) ;
    tPos = 0;
    MPI_Unpack (buf, pSize, &tPos, bufInfo, 3, MPI_INT, comm) ;
    MPI_Unpack (buf, pSize, &tPos, mu, 1, MPI_DOUBLE, comm) ;
    *chunkStart = bufInfo[0];
    *chunkSize = bufInfo[1];
    *chunkDest = bufInfo[2];
  }

  void dynamic_schedule_rule::SendInput (int dest, int tStart, int tSize)
  {
    int position, pSize ;
    int bufInfo[3];

    MPI_Comm comm = facts1->get_comm() ;
    int size ;
    MPI_Pack_size(3, MPI_INT, comm, &size) ;
    pSize = size ;
    pSize += inputPackSize (tStart, tSize);

    GetBuffer (pSize);
    //Pack chunk info
    bufInfo[0] = tStart;
    bufInfo[1] = tSize;
    bufInfo[2] = pSize;
    position = 0;
    MPI_Pack (bufInfo, 3, MPI_INT, buf, pSize, &position, comm) ;
    entitySet myent = interval (tStart, tStart + tSize - 1);
    for (variableSet::const_iterator vi = inputs.begin ();
	 vi != inputs.end (); ++vi) {
      storeRepP s_ptr = facts1->get_variable (*vi);
      s_ptr->pack (&buf[0], position, pSize, myent);
    }
    MPI_Send (buf, pSize, MPI_PACKED, dest, WORK_REMOTE, comm) ;

#if DEBUGOUT
    Loci::debugout << "SendInput " << myRank << "->" << dest << ", chunk="
                   << tStart << "," << tSize << "," << pSize << endl;
#endif
  }

  void dynamic_schedule_rule::ReceiveInput (int src, int msgSize, int *tStart,
                                            int *tSize) {
    int position, size;
    MPI_Status tStatus;
    int bufInfo[3];

    MPI_Comm comm = facts1->get_comm() ;
    GetBuffer (msgSize);
    MPI_Recv (buf, msgSize, MPI_PACKED, src, WORK_REMOTE, comm,
	      &tStatus) ;
    //unpack chunk info
    position = 0;
    MPI_Unpack (buf, msgSize, &position, bufInfo, 3, MPI_INT, comm) ;
    *tStart = bufInfo[0];
    *tSize = bufInfo[1];
    size = bufInfo[2];
    //unpack inputs into local facts
    for (variableSet::const_iterator vi = inputs.begin ();
	 vi != inputs.end (); ++vi) {
      storeRepP s_ptr = local_facts1->get_variable (*vi);
      s_ptr->unpack (&buf[0], position, size,
                     sequence (interval (0, *tSize - 1)));
    }

#if DEBUGOUT
    Loci::debugout << "ReceiveInput " << myRank << "<-" << src
                   << ", chunk=" << *tStart << "," << *tSize << ","
                   << size << endl;
#endif
  }

  void dynamic_schedule_rule::SendOutput (int dest, int tStart, int tSize,
                                          double *tTime) {
    int position, pSize;
    int bufInfo[3];

    //compute buffer size
    //    size = 0;
    //    for (variableSet::const_iterator vi = outputs.begin ();
    //	 vi != outputs.end (); ++vi)
    //      {
    //	storeRepP s_ptr = local_facts1->get_variable (*vi);
    //	size += s_ptr->pack_size (interval (0, tSize - 1));
    //      }
    //    pSize = size + 3 * sizeof (int) + sizeof (double);
    MPI_Comm comm = facts1->get_comm() ;
    int size ;
    MPI_Pack_size(3, MPI_INT, comm, &size) ;
    pSize = size ;
    MPI_Pack_size(1, MPI_DOUBLE, comm, &size) ;
    pSize += size ;
    pSize += outputPackSize (0, tSize);

    GetBuffer (pSize);
    //Pack chunk info
    bufInfo[0] = tStart;
    bufInfo[1] = tSize;
    bufInfo[2] = pSize;
    position = 0;
    MPI_Pack (bufInfo, 3, MPI_INT, buf, pSize, &position, comm) ;
    MPI_Pack (tTime, 1, MPI_DOUBLE, buf, pSize, &position, comm) ;
    //Pack outputs
    entitySet myent2 = interval (0, tSize - 1);
    for (variableSet::const_iterator vi = outputs.begin ();
	 vi != outputs.end (); ++vi) {
      storeRepP s_ptr = local_facts1->get_variable (*vi);
      s_ptr->pack (&buf[0], position, pSize, myent2);
    }
    //Send outputs
    MPI_Send (buf, pSize, MPI_PACKED, dest, RECV_OUTPUT, comm) ;

#if DEBUGOUT
    Loci::debugout << "SendOutput " << myRank << "->" << dest << ", chunk="
                   << tStart << "," << tSize << "," << pSize << endl;
#endif
  }

  void dynamic_schedule_rule::ReceiveOutput (int src, int msgSize, int *iters,
                                             double *tTime) {
    MPI_Status tStatus;
    int position, tStart, tSize, size;
    int bufInfo[3];

    MPI_Comm comm = facts1->get_comm() ;
    GetBuffer (msgSize);
    MPI_Recv (buf, msgSize, MPI_PACKED, src, RECV_OUTPUT, comm,
	      &tStatus) ;
    //unpack chunk info
    position = 0;
    MPI_Unpack (buf, msgSize, &position, bufInfo, 3, MPI_INT, comm) ;
    MPI_Unpack (buf, msgSize, &position, tTime, 1, MPI_DOUBLE, comm) ;
    tStart = bufInfo[0];
    tSize = bufInfo[1];
    size = bufInfo[2];
    *iters = tSize;
    //unpack outputs into facts

    for (variableSet::const_iterator vi = outputs.begin ();
	 vi != outputs.end (); ++vi) {
      storeRepP s_ptr = facts1->get_variable (*vi);
      s_ptr->unpack (&buf[0], position, size,
                     sequence (interval (tStart, tStart + tSize - 1)));
    }

#if DEBUGOUT
    Loci::debugout << "ReceiveOutput " << myRank << "<-" << src << ", chunk="
                   << tStart << "," << tSize << "," << size << endl;
#endif
  }


  void dynamic_schedule_rule::GetChunkSize (int method, int minChunkSize,
                                            int source, int *yMap,
                                            int *chunkSize, int *batchSize,
                                            int *batchRem) {
    int i, tChunk, rem;

    rem = 0;
    for (i = 0; i < nProcs; i++)
      rem += yMap[2 * i + 1];

    switch (method) {

    case FSC:		//fixed size scheduling, P * log(N / P) chunks
      i = (allItems + nProcs - 1) / nProcs;
      tChunk = int ((0.55 + i * log (2.0e0) / log ((double) i)));
      tChunk = min (tChunk, rem);
      *batchSize = tChunk;
      *batchRem = min (*batchSize, rem);
      break;

    case GSS:		//guided self scheduling
      tChunk = max (minChunkSize, (rem + nProcs - 1) / nProcs);
      tChunk = min (rem, tChunk);
      *batchSize = tChunk;
      *batchRem = min (*batchSize, rem);
      break;

    case FAC:		//factoring
      if (*batchRem == 0) {
        tChunk = max (minChunkSize, (rem + 2 * nProcs - 1) / (2 * nProcs));
        *batchSize = nProcs * tChunk;
        *batchRem = min (*batchSize, rem);
      }
      // use current batchSize
      tChunk = max (minChunkSize, *batchSize / nProcs);
      tChunk = min (rem, tChunk);
      break;

    default:
      //no scheduling
      tChunk = yMap[2 * source + 1];
      *batchSize = tChunk;
      *batchRem = min (*batchSize, rem);
      break;

    }

    //adjust according to Remaining[source]
    if (yMap[2 * source + 1] > (tChunk + minChunkSize / 2))
      *chunkSize = tChunk;
    else
      *chunkSize = yMap[2 * source + 1];

    //adjust remaining in batch
    *batchRem -= *chunkSize;
    if (*batchRem < minChunkSize)
      *batchRem = 0;

    //adjust counters for source
    yMap[2 * source] += *chunkSize;
    yMap[2 * source + 1] -= *chunkSize;
  }


  void dynamic_schedule_rule::loop_scheduling ()
  {
    vector<int> yMap(2*nProcs) ;
    

    MPI_Status mStatus;
    int action = 0, lastAction = 0;
    int chunkSize = 0, chunkStart = 0, chunkDest = 0;
    int batchSize = 0, batchRem = 0, minChunkSize = 0;
    vector<int> inputSent(nProcs) ;

    int localStart, localSize;
    int remoteStart, remoteSize;
    int saveStart, saveSize, saveSrc;

    int numIdle;
    int returns;		//output data
    int incoming;		//input data
    int gotWork;

    int msgSrc = 0, msgLen = 0;
    int tStart = 0, tSize = 0, tSource = 0;
    int MsgInQueue = 0;
    int sendRequest, probeFreq;

    double tTime, currentChunkTime, latestFinish, execTime;
    double  timerDiff;
    //double idleTime;          // time in MPI_Probe

    int i;

    local1 = 0.0;
    remote1 = 0.0;
    currentChunkTime = 0.0;
    if(ewt.size() == 0) {
      vector<double> tmp(nProcs) ;
      ewt.swap(tmp) ;
    }
    for (i = 0; i < nProcs; i++) {
      workTime[i] = 0.0;
      ewt[i] = aveWorkTime[i];
    }

    // allocate space for LB
    tSize = (allItems + 4 * nProcs - 1) / (4 * nProcs);
    // add 25% saftey margin for allocation  (some problems found
    // with under-allocation)
    AllocateLBspace (tSize+tSize/4+1);

    returns = 0;
    incoming = 0;
    numIdle = 0;
    probeFreq = MIN_PROBE_FREQ;
    sendRequest = MIN_SEND_REQUEST;
    gotWork = 1;

    localSize = 0;
    localStart = 0;

    remoteSize = 0;
    remoteStart = 0;
    saveSrc = 0;
    saveSize = 0;
    saveStart = 0;


    // =======================================================
    // Initial loads
    // =======================================================
    // if (i_am_foreMan)
    //   Send initial local work info to workers
    //   Set own local work info, nextAction = WORK_LOCAL
    // else // i_am_worker
    //   Set nextAction = WAIT4_MSG

    if (myRank == foreMan) {
      //scheduler initializations
      batchSize = 0;
      batchRem = 0;
      minChunkSize = 2 * max (MIN_PROBE_FREQ, MIN_SEND_REQUEST);
      for (i = 0; i < (2 * nProcs); i++)
        yMap[i] = yMapSave[i];

      // use static scheduling during first call to execute()
      //   to initialize average work time
      if (nCalls)
        i = LBMethod;
      else
        i = -1;
      //send local work info to others
      for (tSource = 0; tSource < nProcs; tSource++)
        if (tSource != myRank) {
          inputSent[tSource] = 0;
          tStart = yMap[2 * tSource];
          GetChunkSize (i, minChunkSize, tSource,
                        &yMap[0], &tSize, &batchSize, &batchRem);
          SendInfo (tSource, WORK_LOCAL, tStart, tSize, -1, 0.0);
        }
      //foreman collects work for itself
      inputSent[myRank] = 0;
      chunkStart = yMap[2 * myRank];
      GetChunkSize (i, minChunkSize, myRank,
                    &yMap[0], &chunkSize, &batchSize, &batchRem);
      action = WORK_LOCAL;
    } else {
      //workers block for first message
      action = WAIT4_MSG;
    }
    //
    // =======================================================
    // Main loop
    // =======================================================
    stopWatch sc ;

   
    while (gotWork + localSize + remoteSize + incoming + returns) {
      sequence seq ;
      switch (action) {

        // =======================================================
        //     case WAIT4_MSG
        //       MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, ...)
        //       Set nextAction = RETRIEVE_MSG, lastAction = WAIT4_MSG
        // =======================================================

      case WAIT4_MSG:	// wait for a message
        //timerStart = MPI_Wtime ();
#if DEBUGOUT
        Loci::debugout << "MPI_Probe() ..." << endl;
#endif
        MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG,
                   facts1->get_comm(), &mStatus) ;
        //timerDiff = MPI_Wtime () - timerStart;
        //idleTime += timerDiff;
        action = RETRIEVE_MSG;
        lastAction = WAIT4_MSG;
        break;

        // =======================================================
        //     case TEST4_MSG
        //       MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, ...)
        //       if (there was a message) Set nextAction = RETRIEVE_MSG
        //       else Set nextAction = lastAction
        // =======================================================

      case TEST4_MSG:	// test for a message
        MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG,
                    facts1->get_comm(), &MsgInQueue, &mStatus) ;
        if (MsgInQueue)
          action = RETRIEVE_MSG;
        else
          action = lastAction;
        break;

        // =======================================================
        //     case RETRIEVE_MSG
        //       Set nextAction from message status.MPI_TAG (one of the #define)
        //       Retrieve content of message, which may be
        //         - signal to quit
        //         - info on local work to execute
        //         - info on input items to migrate
        //         - a request for load & timing data (for foreman)
        //         - input items
        //         - output items & timing data
        // =======================================================

      case RETRIEVE_MSG:	// retrieve message contents according to tag
        msgSrc = mStatus.MPI_SOURCE;
        action = mStatus.MPI_TAG;
        if ((action == QUIT) || (action == WORK_LOCAL) ||
            (action == SEND_INPUT) || (action == FILL_REQUEST)) {
          RecvInfo (msgSrc, action, &chunkStart, &chunkSize, &chunkDest,
                    &execTime);
        } else {			// message contains MPI_PACKED data

          MPI_Get_count (&mStatus, MPI_PACKED, &msgLen);
          if (action == WORK_REMOTE)
            ReceiveInput (msgSrc, msgLen, &chunkStart, &chunkSize);
          else	{	// if (action==RECV_OUTPUT)
            ReceiveOutput (msgSrc, msgLen, &i, &tTime);
            remote1 += tTime;
          }

        }
        break;

        // =======================================================
        //     case QUIT
        //       Set one of the termination flags, nextAction = WAIT4_MSG
        // =======================================================

      case QUIT:		// no more work, but there may be messages
        gotWork = 0;
        action = WAIT4_MSG;
        break;

        // =======================================================
        //     case WORK_LOCAL
        //       Update local work info
        //       Compute probing frequency, when to request next work info
        //       Set nextAction = EXEC_LOCAL_PART
        // =======================================================

      case WORK_LOCAL:	// update local work info
        localStart = chunkStart;
        localSize = chunkSize;
        currentChunkTime = 0.0;
        if (myRank != foreMan) {
          probeFreq = max (MIN_PROBE_FREQ, int (localSize * PCT_PROBE_FREQ));
          sendRequest = max (MIN_SEND_REQUEST,
                             int (localSize * PCT_SEND_REQUEST));
          if (returns)
            probeFreq = MIN_PROBE_FREQ;
        }
        action = EXEC_LOCAL_PART;	// execute local subchunk
        break;

        // =======================================================
        //     case WORK_REMOTE
        //       Update remote work info
        //       if (i_am_foreman) Set nextAction = EXEC_REMOTE_PART
        //       else Set nextAction = EXEC_REMOTE_WHOLE
        // =======================================================

      case WORK_REMOTE:	// receive remote input data
        saveStart = chunkStart;
        saveSize = chunkSize;
        saveSrc = msgSrc;
        remoteStart = 0;
        remoteSize = chunkSize;
        currentChunkTime = 0.0;
        if (myRank == foreMan)
          action = EXEC_REMOTE_PART;	// subchunk
        else
          action = EXEC_REMOTE_WHOLE;	// whole chunk
        break;

        // =======================================================
        //     case SEND_INPUT
        //       Send input according to info retrieved from message
        //       Set nextAction = lastAction
        // =======================================================

      case SEND_INPUT:	// send input data
        SendInput (chunkDest, chunkStart, chunkSize);
        returns++;
        probeFreq = min (probeFreq, MIN_PROBE_FREQ);
        action = lastAction;
        break;

        // =======================================================
        //     case RECV_OUTPUT
        //       Decrement counter for outputs to receive
        //       Set nextAction = lastAction
        // =======================================================

      case RECV_OUTPUT:	// receive output
        returns--;
        action = lastAction;
        break;

        // =======================================================
        //     case EXEC_LOCAL_PART
        //       If this is the last subchunk, then send request for next work info
        //       Execute a subchunk from current local work
        //       Update counters and accumulators
        //       if (current local chunk not finished)
        //         Set nextAction = TEST4_MSG, lastAction = EXEC_LOCAL_PART
        //       else if (i_am_foreman) // foreman finished current local chunk
        //         // pretend foreman sent request (to itself)
        //         Set lastAction = WAIT4_MSG, nextAction = FILL_REQUEST
        //       else Set nextAction = WAIT4_MSG
        // =======================================================

      case EXEC_LOCAL_PART:	// execute a local subchunk
        // Not yet the last subchunk?
        if (localSize > (sendRequest + probeFreq))
          tSize = probeFreq;	// no
        else
          tSize = localSize;	// yes
        //Send request before doing last subchunk
        if ((localSize == tSize) && (myRank != foreMan)) {
          SendInfo (foreMan, FILL_REQUEST, 0, 0, myRank, currentChunkTime);
        }

        sc.start() ;
	seq = sequence(interval (localStart, localStart + tSize - 1)) ;
        main_comp->prelude(seq) ;
        main_comp->compute(seq) ;
        timerDiff = sc.stop() ;
        comp_timer.addTime(timerDiff,1) ;
        workTime[myRank] += timerDiff;

        localStart += tSize;
        localSize -= tSize;
        currentChunkTime += timerDiff;
        if (localSize) {	// current local chunk not finished
          action = TEST4_MSG;	// test for a message
          lastAction = EXEC_LOCAL_PART;
        } else if (myRank == foreMan) {	// foreman finished current local chunk
          // pretend foreman sent request (to itself)
          lastAction = WAIT4_MSG;	// foreman was waiting for a message
          msgSrc = myRank;	// foreman received a message (from self)
          ewt[myRank] -= currentChunkTime;	// decrement remaining work time
          action = FILL_REQUEST;	// message was a request
        } else {
          action = WAIT4_MSG;	// wait for a message
        }
        break;

        // =======================================================
        //     case EXEC_REMOTE_PART
        //       Foreman executes subchunk of remote items
        //       if (chunk not finished)
        //         Set nextAction = TEST4_MSG, lastAction = EXEC_REMOTE_PART
        //       else // foreman finished current remote chunk
        //         Return outputs, Set lastAction = WAIT4_MSG
        //         // pretend foreman sent a request to self
        //         Set nextAction = FILL_REQUEST
        // =======================================================

      case EXEC_REMOTE_PART:	// foreman executes remote subchunk
        tSize = min (remoteSize, MIN_PROBE_FREQ);

        sc.start() ;

	seq = sequence(interval (remoteStart,remoteStart + tSize - 1)) ;

	local_comp1-> prelude (seq) ;
	local_comp1-> compute (seq) ;


        timerDiff = sc.stop() ;
        comp_timer.addTime(timerDiff,1) ;
        currentChunkTime += timerDiff;

        remoteStart += tSize;
        remoteSize -= tSize;
        if (remoteSize) {		// current remote chunk not finished
          action = TEST4_MSG;   	// test for a message
          lastAction = EXEC_REMOTE_PART;
        } else {		// foreman finished current remote chunk
          SendOutput (saveSrc, saveStart, saveSize, &currentChunkTime);
          workTime[saveSrc] += currentChunkTime;
          remote1 += currentChunkTime;
          incoming--;	// decrement counter for incoming data
          // pretend foreman sent request (to itself)
          lastAction = WAIT4_MSG;	// foreman waited for a message
          msgSrc = myRank;	// foreman received a message (from self)
          action = FILL_REQUEST;	// message was a request
          // RecvInfo (msgSrc, action, &chunkStart, &chunkSize, &chunkDest,
          // 	  &execTime);
          chunkStart = 0;
          chunkSize  = 0;
          chunkDest  = saveSrc;
          execTime = currentChunkTime;
        }
        break;

        // =======================================================
        //     case EXEC_REMOTE_WHOLE
        //       Worker executes remote chunk, return the outputs,
        //         requests for new work info
        //       Set nextAction = WAIT4_MSG
        // =======================================================

      case EXEC_REMOTE_WHOLE:	// execute remote chunk, send output, request
        sc.start() ;
	seq = sequence(interval (0, remoteSize - 1));
	local_comp1->prelude (seq) ;
	local_comp1->compute (seq) ;

        currentChunkTime = sc.stop() ;
        comp_timer.addTime(currentChunkTime,1) ;
        workTime[saveSrc] += currentChunkTime;
        remote1 += currentChunkTime;

        SendOutput (saveSrc, saveStart, saveSize, &currentChunkTime);
        SendInfo (foreMan, FILL_REQUEST, 0, 0, saveSrc, currentChunkTime);
        action = WAIT4_MSG;	// wait for a message
        remoteSize = 0;
        break;

        // =======================================================
        //     case FILL_REQUEST
        //       Foreman takes one of the following actions,
        //         - inform requesting proc that there are no more items
        //              to execute
        //         - send new local work info to requesting proc
        //         - inform a heavily loaded worker to send items to
        //              requesting proc
        // =======================================================

      case FILL_REQUEST:	// received a request (foreman only)
        ewt[chunkDest] -= execTime;	// decrement remaining work time of chunk owner
        // find a data source
        if (yMap[2 * msgSrc + 1] > 0)
          tSource = msgSrc;	// requester has unfinished data
        else if (inputSent[msgSrc])	// requester has none
          tSource = -1;	// prohibit requester from accepting work since it gave away data
        else {			// find worker expected to finish last
          tSource = -1;
          latestFinish = 0.0;
          for (i = 0; i < nProcs; i++)
            if (yMap[2*i+1] && (ewt[i]>latestFinish)) {
              tSource = i;
              latestFinish = ewt[i];
            }
          // Uncomment/comment next line to disable/enable load balancing
          //tSource = -1;
        }

        // who's requesting ?
        if (msgSrc == foreMan) {
          chunkStart = 0;
          chunkSize = 0;
          chunkDest = 0;
          if (tSource == -1) {		// no more unscheduled iterates
            action = WAIT4_MSG;	// wait for a message
          }		// end if (tSource == -1)
          else if (tSource == msgSrc) {		// msgSrc has some
            chunkStart = yMap[2 * msgSrc];
            GetChunkSize (LBMethod, minChunkSize, msgSrc,
                          &yMap[0], &chunkSize, &batchSize, &batchRem);
            action = WORK_LOCAL;	// update local work info
#if DEBUGOUT
            Loci::debugout << "Info0:" << myRank << "->" << tSource
                           << "->-1 : " << chunkStart << "," << chunkSize
                           << endl;
#endif
          }		//end if (tSource == msgSrc)
          else {		// tSource sends input data to msgSrc
            tStart = yMap[2 * tSource];
            GetChunkSize (LBMethod, minChunkSize, tSource,
                          &yMap[0], &tSize, &batchSize, &batchRem);
            SendInfo (tSource, SEND_INPUT, tStart, tSize, msgSrc,
                      0.0);
            inputSent[tSource]++;	// mark tSource as data giver
            action = WAIT4_MSG;	// wait for a message
            incoming++;	// increment counter for incoming work
          }
        }			// end if (msgSrc==foreMan)
        else {			// request was from a worker
          if (tSource == -1) {		// no more unscheduled iterates
            numIdle++;
            SendInfo (msgSrc, QUIT, 0, 0, 0, 0.0);
            gotWork = numIdle != (nProcs - 1);
            action = TEST4_MSG;	// test for a message
          }		// end if (tSource == -1)
          else if (tSource == msgSrc) {		// msgSrc has some
            tStart = yMap[2 * msgSrc];
            GetChunkSize (LBMethod, minChunkSize, msgSrc,
                          &yMap[0], &tSize, &batchSize, &batchRem);
            SendInfo (msgSrc, WORK_LOCAL, tStart, tSize, -1, 0.0);
            action = TEST4_MSG;	// test for a message
          }		//end if (tSource == msgSrc)
          else {		// tSource sends input data to msgSrc
            tStart = yMap[2 * tSource];
            GetChunkSize (LBMethod, minChunkSize, tSource,
                          &yMap[0], &tSize, &batchSize, &batchRem);
            inputSent[tSource]++;	// mark tSource as data giver
            if (tSource != foreMan) {
              SendInfo (tSource, SEND_INPUT, tStart, tSize, msgSrc,
                        0.0);
              action = TEST4_MSG;	// test for a message
            } else {
              // pretend foreman sent request (to itself)
              // pretend foreman waited for a message
              // pretend foreman received a message (from self)
              action = SEND_INPUT;	// message was to send input
              chunkStart = tStart;
              chunkSize = tSize;
              chunkDest = msgSrc;
#if DEBUGOUT
              Loci::debugout << "Info0:" << myRank << "->" << tSource
                             << "->" << chunkDest << " : " << chunkStart
                             << "," << chunkSize << endl;
#endif
            }
          }
        }			// end else (msgSrc!=foreMan)
        break;

      }			// end switch (action)
    }				// end while (gotWork + localSize + remoteSize + returns)

    //
    // MPI_Barrier (MPI_COMM_WORLD)
    //



#if DEBUGOUT
    Loci::debugout << "Exits: " << myRank << endl;
#endif
    FreeLBspace ();

#if SHOWTIMES
    wall1 = MPI_Wtime () - wall1;
    local1 = workTime[myRank];

    // BARRIER: gather work times on all procs
    for (i = 0; i < nProcs; i++)
      ewt[i] = 0.0;

  
    MPI_Allreduce (&workTime[0], &ewt[0], nProcs, MPI_DOUBLE, MPI_SUM,
		   facts1->get_comm()) ;
  

    // running average of work times
    if (nCalls)
      for (i = 0; i < nProcs; i++)
	aveWorkTime[i] = (1.0-CONTRIB )*aveWorkTime[i] + CONTRIB * ewt[i];
    else
      for (i = 0; i < nProcs; i++)
	aveWorkTime[i] = ewt[i];

    // perfectly balance work time
    ideal1 = 0.0;
    for (i = 0; i < nProcs; i++)
      ideal1 += ewt[i];
    ideal1 /= nProcs;
#else
    // simple BARRIER
    MPI_Barrier (facts1->get_comm()) ;
#endif
  }


  //=======================================================
  // The execute() routine
  //=======================================================

  void dynamic_schedule_rule::execute (fact_db & facts, sched_db &scheds)
  {

    stopWatch s ;
    s.start() ;

    // Hack to handle non-contiguos sets.
    if(compress_set) { 
#ifdef VERBOSE
      debugout << "allocating over set" << exec_set << endl ;
#endif
      // need to allocate sets and copy from main fact_db
      for (variableSet::const_iterator vi = inputs.begin ();
           vi != inputs.end (); ++vi) {
        storeRepP sp = backup_facts.get_variable (*vi);
        sp->allocate (exec_set);
      }
      for (variableSet::const_iterator vi = outputs.begin ();
           vi != outputs.end (); ++vi) {
        storeRepP sp = backup_facts.get_variable (*vi);
        sp->allocate (exec_set);
      }
      // now copy input data
      int bufsize = 0 ;
      for (variableSet::const_iterator vi = inputs.begin ();
           vi != inputs.end (); ++vi) {
        storeRepP s_ptr = facts.get_variable (*vi);
        bufsize += s_ptr->pack_size (given_exec_set);
      }
#ifdef VERBOSE
      debugout << "copying data " << bufsize << " bytes" << endl ;
#endif
      vector<unsigned char> data(bufsize) ;
      int position = 0 ;
      for (variableSet::const_iterator vi = inputs.begin ();
           vi != inputs.end (); ++vi) {
        storeRepP s_ptr = facts.get_variable(*vi) ;
        s_ptr->pack(&data[0],position,bufsize,given_exec_set) ;
      }
      position = 0 ;
      for (variableSet::const_iterator vi = inputs.begin ();
           vi != inputs.end (); ++vi) {
        storeRepP s_ptr = backup_facts.get_variable(*vi) ;
        s_ptr->unpack(&data[0],position,bufsize,exec_set) ;
      }
    }
    wall1 = MPI_Wtime ();
    wall2 = wall1;
    if (LBMethod == NLB) {
      stopWatch sc ;
      sc.start() ;
      sequence seq(exec_set);
      main_comp->prelude (seq) ;
      main_comp->compute (seq) ;
      comp_timer.addTime(sc.stop(),1) ;
      wall1 = MPI_Wtime () - wall1;
      local1 = wall1;
      ideal1 = wall1;
      wall2 = wall1;
    } else {
      // globals
      local_comp1 = local_compute1;
      
      facts1 = &facts;
      if(compress_set)
	facts1 = &backup_facts ;
      local_facts1 = &local_facts ;

      if (LBMethod == IWS)
        iterative_weighted_static ();
      else
        loop_scheduling ();

      // ideal1, local1, remote1 & wall1 are computed
      //   by the routines
      wall2 = MPI_Wtime () - wall2;

      facts1 = 0;
      local_facts1 = 0;
    }
#if SHOWTIMES
    Loci::debugout << LBMethod << "," << nProcs << "," << myRank << ","
                   << nCalls << ", ideal=" << ideal1 << ", wall1=" << wall1
                   << " (" << int (1000.0 * (ideal1 - wall1) / ideal1) / 10.0
                   << "%)" << ", wall2=" << wall2
                   << " (" << int (1000.0 * (ideal1 - wall2) / ideal1) / 10.0
                   << "%)" << ", local=" << local1 << ", remote=" << remote1
                   << ", ovrhd1="
                   << int (1000.0 * (wall1 - local1 - remote1) / wall1) / 10.0
                   << "%" << endl;

#endif

    if(compress_set) { // copy output back to main fact database, deallocate
      // memory

      // now copy input data
      int bufsize = 0 ;
      for (variableSet::const_iterator vi = outputs.begin ();
           vi != outputs.end (); ++vi) {
        storeRepP s_ptr = backup_facts.get_variable (*vi);
        bufsize += s_ptr->pack_size (exec_set);
      }
      vector<unsigned char> data(bufsize) ;
      int position = 0 ;
      for (variableSet::const_iterator vi = outputs.begin ();
           vi != outputs.end (); ++vi) {
        storeRepP s_ptr = backup_facts.get_variable(*vi) ;
        s_ptr->pack(&data[0],position,bufsize,exec_set) ;
      }
      position = 0 ;
      for (variableSet::const_iterator vi = outputs.begin ();
           vi != outputs.end (); ++vi) {
        storeRepP s_ptr = facts.get_variable(*vi) ;
        s_ptr->unpack(&data[0],position,bufsize,given_exec_set) ;
      }
      
      // release memory
      for (variableSet::const_iterator vi = inputs.begin ();
           vi != inputs.end (); ++vi) {
        storeRepP sp = backup_facts.get_variable (*vi);
        sp->allocate (EMPTY);
      }
      for (variableSet::const_iterator vi = outputs.begin ();
           vi != outputs.end (); ++vi) {
        storeRepP sp = backup_facts.get_variable (*vi);
        sp->allocate (EMPTY);
      }
    }      
    nCalls++;
    timer.addTime(s.stop(),1) ;
  }


  void dynamic_schedule_rule::Print (std::ostream & s) const {
    s << "dynamic schedule " << rule_tag << "  over set " << exec_set << endl;
  }

  void dynamic_schedule_rule::dataCollate(collectData &data_collector) const {
    ostringstream oss ;
    oss << "dynamic_schedule: " << rule_tag;
    int group = data_collector.openGroup("dynamic-schedule") ;
    data_collector.accumulateTime(comp_timer,EXEC_COMPUTATION,oss.str()) ;
    data_collector.closeGroup(group) ;

    group = data_collector.openGroup("dynamic-schedule-overhead") ;
    timeAccumulator ov ;
    ov.addTime(timer.getTime()-comp_timer.getTime(),timer.getEvents()) ;

    data_collector.accumulateTime(ov,EXEC_COMMUNICATION,oss.str())  ;
    data_collector.closeGroup(group) ;
  }

}
