//#############################################################################
//#
//# Copyright 2014-2025, Mississippi State University
//#
//# The GridMover is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The GridMover software is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the GridMover software. If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <vector>
#include <iostream>
#include <algorithm>
#include "mpi.h"
#include <math.h>
#include <unistd.h>

#define TAG_INFO     10002		// size of succeeding message
#define TAG_INPUT    10003		// input data
#define TAG_OUTPUT   10004		// output data
#define TAG_TIMES    10005          // timing data

//#define VERBOSE
//#define SILENT

namespace gridMotion {
  class stopWatch {
    double start_time ;
  public:
    void start() { // This method resets the clock
      start_time = MPI_Wtime() ;
    }
    double stop() { // This method returns time since last start call
      return MPI_Wtime()-start_time ;
    }
  } ;
  
  // Structure to hold iterative weighted static schedule (and other information
  // needed to generate new IWS schedules).
  struct IWSSchedule {
    MPI_Comm comm;
    int MPI_rank ;
    int MPI_processes ;
    struct chunkCommInfo {
      int proc ;
      std::vector<int> chunkList ;
      int send_size,recv_size ;
    } ;
    
    struct chunkInfo {
      int start;
      int end;
      float chunkTime[2];
    } ;

    std::vector<chunkInfo> chunkData;
    std::vector<int> selfChunks;
    int numBalances;
    int numRemoteChunks;
    
    std::vector<chunkCommInfo> sendChunks;
    std::vector<chunkCommInfo> recvChunks;
    
    void initialize(MPI_Comm C, int local_size) ;
    void generateSchedule() ;
  } ;

  // Perform IWS scheduling executing workFunc passing it data from input_data
  // and returning data to output_data
  template<typename TO, typename TI, typename F>
  void iterative_weighted_static (IWSSchedule &iwss,
                                  // Note, size of array defined by initialize
                                  TO * output_data,
                                  const TI * input_data,
                                  const F &workFunc) 
  {
    using std::vector ;
    using std::abs ;
    using std::max ;
    using std::cout ;
    using std::cerr ;
    using std::endl ;
    using std::pair ;
    
    vector<IWSSchedule::chunkCommInfo> &sendChunks = iwss.sendChunks ;
    vector<IWSSchedule::chunkCommInfo> &recvChunks = iwss.recvChunks ;
    int remote_size = 0 ;
    for(size_t i = 0; i < recvChunks.size(); ++i) {
      remote_size += recvChunks[i].recv_size;
    }

    vector<TI> remote_data(remote_size) ;
    vector<TO> remote_output(remote_size) ;

    int MPI_rank = iwss.MPI_rank ;
    int MPI_processes = iwss.MPI_processes;
  
    stopWatch timer;
    timer.start();
    //=======================================================
    // Send chunks to remote processor
    //=======================================================

    // Send work to remote processors
    if (sendChunks.size() != 0) {
      vector<MPI_Request> req_list(sendChunks.size());
      int tot_buf_size = 0 ;
      for (size_t i = 0; i < sendChunks.size(); ++i) 
        tot_buf_size += sendChunks[i].send_size ;
      
      vector<TI> buf(tot_buf_size);
      
      int position = 0;
      for (size_t i = 0; i < sendChunks.size(); ++i) {
        int buf_size = sendChunks[i].send_size ;
        for (size_t j = 0; j < sendChunks[i].chunkList.size(); ++j) {
          //ch is the index of a chunkInfo object in iwss.chunkData,
          //which contains a chunk range
          int ch = sendChunks[i].chunkList[j];
          for (size_t k = iwss.chunkData[ch].start;
               k <= iwss.chunkData[ch].end; ++k) {
            buf[k-iwss.chunkData[ch].start + position] = input_data[k];
          }
          position += iwss.chunkData[ch].end - iwss.chunkData[ch].start + 1;
        }
        
        int dest = sendChunks[i].proc;
        
        MPI_Isend (&buf[position-buf_size], buf_size * sizeof(TI),
                   MPI_BYTE, dest, TAG_INPUT, iwss.comm, &req_list[i]);
      }
      // Now we will wait for the sends to complete
      int nreqs = req_list.size();
      vector<MPI_Status> status_list(nreqs);
      MPI_Waitall(nreqs, &req_list[0], &status_list[0]);
      
    }
    
    //Receive Stage
    else if(recvChunks.size() != 0) {
      // Irecv first
      int tot_buf_size = 0;
      for (size_t i = 0; i < recvChunks.size(); ++i) {
        tot_buf_size += recvChunks[i].recv_size;
      }
      
      vector<MPI_Request> req_list(recvChunks.size());
      
      vector<TI> buf(tot_buf_size);
      int offset = 0;
      for (size_t i = 0; i < recvChunks.size(); ++i) {
        int buf_size = recvChunks[i].recv_size ;
        int src = recvChunks[i].proc;
        MPI_Irecv (&remote_data[offset], buf_size * sizeof(TI),
                   MPI_BYTE, src, TAG_INPUT, iwss.comm, &req_list[i]);
        offset += buf_size;
      }

      // Allocate space for data being received
      int nchunks = 0;
      for(size_t i = 0; i < recvChunks.size(); ++i) {
        nchunks += recvChunks[i].chunkList.size();
      }

      // Now we will wait for the receives.
      int nreqs = req_list.size();
      vector<MPI_Status> status_list(nreqs);
      MPI_Waitall(nreqs, &req_list[0], &status_list[0]);
      
      
    }

    //End communication timing.
    double comm_time = timer.stop();

    //Done with sends and receives; execution stage.
    // Begin execution timing
    vector<float> remote_times(iwss.numRemoteChunks, 0);
    stopWatch exec_clock;
    exec_clock.start();

    //=======================================================
    // Do computational Work
    //=======================================================
    //Compute local work
    int res = 0;

    for(size_t i = 0; i < iwss.selfChunks.size(); ++i) {
      
      int ch = iwss.selfChunks[i];
      stopWatch s;
      s.start();
      for(size_t j = iwss.chunkData[ch].start; j <= iwss.chunkData[ch].end; ++j) {
        workFunc(output_data[j],input_data[j]) ;
      }
      double elapsed_time = s.stop();
      iwss.chunkData[ch].chunkTime[0] += elapsed_time;
    }
  
    //Now compute remote work
    int offset = 0;
    for(size_t i = 0; i < recvChunks.size(); ++i) {
      stopWatch s;
      s.start();
      for(size_t j = 0; j < recvChunks[i].recv_size; ++j) {
        workFunc(remote_output[offset+j],remote_data[offset+j]) ;
      }
      double elapsed_time = s.stop();
      offset += recvChunks[i].recv_size;
      remote_times[i] = elapsed_time;
    }

    float exec_time = exec_clock.stop();
    //=======================================================
    // Return chunks to owning processor
    //=======================================================
    timer.start();

    if(sendChunks.size() != 0) {
      //Post Irecvs
      int tot_buf_size = 0;

      for(size_t i = 0; i < sendChunks.size(); ++i) {
        tot_buf_size += sendChunks[i].recv_size;
      }

      vector<TO> buf(tot_buf_size);
      int tot_chunks = 0;

      for(size_t i = 0; i < sendChunks.size(); ++i) {
        tot_chunks += sendChunks[i].chunkList.size();
      }

      vector<float>times(tot_chunks);
      vector<MPI_Request> req_list(sendChunks.size()*2);

      int offset1 = 0;
      int offset2 = 0;

      for(size_t i = 0; i < sendChunks.size(); ++i) {
        int buf_size = sendChunks[i].recv_size;
        int dest = sendChunks[i].proc;

#ifdef VERBOSE
        cout << "receiving computed data from " << dest << endl;
#endif
        MPI_Irecv (&buf[offset1], buf_size * sizeof(TO), MPI_BYTE, dest, TAG_OUTPUT,
                   iwss.comm, &req_list[i*2]);
        offset1 += buf_size;
        int nchunks = sendChunks[i].chunkList.size();

        MPI_Irecv(&times[offset2], nchunks, MPI_FLOAT, dest, TAG_TIMES,
                  iwss.comm, &req_list[i*2+1]);
        offset2 += nchunks;
      }

      // wait on recvs to complete
      int nreqs = req_list.size();
      vector<MPI_Status> status_list(nreqs);
      MPI_Waitall(nreqs, &req_list[0], &status_list[0]);
    
      offset1 = 0;
      offset2 = 0;
      int position = 0;
      for(size_t i = 0; i < sendChunks.size(); ++i) {
        int buf_size = sendChunks[i].recv_size;

        for(size_t j = 0; j < sendChunks[i].chunkList.size(); ++j) {
          int ch = sendChunks[i].chunkList[j];
          // Get values from chunkData[ch], and unpack
          for (size_t k = iwss.chunkData[ch].start;
               k <= iwss.chunkData[ch].end; ++k)
            {
              output_data[k] = buf[position++]  ;
            }
        }


        offset1 += buf_size;

        //Receive times
        int nchunks = sendChunks[i].chunkList.size();
        for(int j = 0; j < nchunks; ++j) {
          int ch = sendChunks[i].chunkList[j];
          iwss.chunkData[ch].chunkTime[0] += times[offset2+j];
        }
        offset2 += nchunks;
      }
    }
    else if(recvChunks.size() != 0) {
      vector<MPI_Request> req_list(recvChunks.size()*2);
      // Now send output data back
      int cnt = 0;
      int cnt2 = 0;
      for(size_t i = 0; i < recvChunks.size(); ++i) {
        int buf_size = recvChunks[i].send_size;
        int dest = recvChunks[i].proc;
#ifdef VERBOSE
        cout << "sending computed data to " << dest << endl ;
#endif
        MPI_Isend (&remote_output[cnt], buf_size * sizeof(TO), MPI_BYTE,
                  dest, TAG_OUTPUT,iwss.comm, &req_list[i*2]);
        
        cnt += recvChunks[i].recv_size ;
        int nchunks = recvChunks[i].chunkList.size() ;
        MPI_Isend(&remote_times[cnt2],nchunks,MPI_FLOAT,
                 dest,TAG_TIMES,iwss.comm, &req_list[i*2+1]) ;
        cnt2 += nchunks ;
      }
      // wait on sends to complete
      int nreqs = req_list.size();
      vector<MPI_Status> status_list(nreqs);
      MPI_Waitall(nreqs, &req_list[0], &status_list[0]);
    
    }
    comm_time += timer.stop();

#ifdef VERBOSE
    cout << MPI_rank << " done communication, collect times" << endl;
#endif
    //=======================================================
    // Estimate parallel efficiency
    //=======================================================
    float total_time = 0 ;
    MPI_Allreduce(&exec_time,&total_time,1,MPI_FLOAT,MPI_SUM,iwss.comm) ;

    float max_time = 0 ;
    MPI_Allreduce(&exec_time,&max_time,1,MPI_FLOAT,MPI_MAX,iwss.comm) ;
    double ef = total_time/(max_time*float(MPI_processes)) ;
    
#ifndef SILENT
    if(MPI_rank == 0) 
      cout << "parallel efficiency of IWS scheduled computation is "
           << ef*100.0 << "%, time =" << max_time << endl ;
#endif
    //=======================================================
    // If efficiency too low, regenerate load balance schedule
    //=======================================================
    timer.start() ;
    if((iwss.numBalances == 0 && ef < .5) ||
       (iwss.numBalances != 0 && ef < 0.89 && (iwss.numBalances&0x3) == 0) ) {
      iwss.generateSchedule() ;
    }
  
#ifdef VERBOSE
    double sched_time = timer.stop() ;
    double timesin[2] = {comm_time,sched_time} ;
    double timesmx[2] = {0,0} ;
  
    MPI_Allreduce(&timesin[0],&timesmx[0],2,MPI_DOUBLE,MPI_MAX,iwss.comm) ;
    if(MPI_rank == 0)
    cout << "comm_time = " << timesmx[0] << ",sched_time=" << timesmx[1] << endl ;
#endif
    //=======================================================
    // If past trigger point, reset timing data
    //=======================================================
    if((iwss.numBalances&0x7) == 0) {
      for(size_t i=0;i<iwss.chunkData.size();++i) {
        iwss.chunkData[i].chunkTime[1] = iwss.chunkData[i].chunkTime[0] ;
        iwss.chunkData[i].chunkTime[0] = 0 ;
      }
    }
  
    iwss.numBalances++ ;
  }
  
}
