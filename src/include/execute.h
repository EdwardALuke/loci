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
#ifndef EXECUTE_H
#define EXECUTE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <ostream>
#include <vector>
#include <string>
#include <list>
using std::string;

#include <Config/conf.h>
#include <Tools/cptr.h>
#include <Tools/intervalSet.h>
#include <mpi.h>

#define PROFILER
#ifdef USE_PAPI
#include <papi.h>
#include <papiStdEventDefs.h>
#endif

#ifdef PAPI_DEBUG
#include <papi.h>
#else
#define long_long long
#endif

namespace Loci {

  /** ************************************************************************
   * @brief This class is used for timing a period of execution to the
   * microsecond resolution in seconds
   ***************************************************************************/
  class stopWatch {
    double start_time ;
  public:
    /**
     * @brief The constructor initializes the timer to record that the timer
     * has not been started
     **/
    stopWatch() :start_time(-1.0) {}
    /**
     * @brief A static method that returns the current time with microsecond
     * resolution.  This could be used as a portable way to access time, but
     * generally this method would be used internally to the stopWatch class.
     * @return double
     **/
    static double timer() {
#ifdef PROFILER
#ifdef USE_PAPI
      return 1e-6*double(PAPI_get_real_usec());
#else
      return MPI_Wtime() ;
#endif
#else
      timeval tv ;
      struct timezone tz ;
      gettimeofday(&tv,&tz) ;
      return double(tv.tv_sec)+1e-6*double(tv.tv_usec) ;
#endif
    }
    /**
     * @brief This method starts the timer.
     **/
    void start() { 
      start_time = stopWatch::timer() ;
    }
    /**
     * @brief This method returns the time since the last call to start the
     * timer.  If start has not been called, it returns a negative value.
     * Method is declared a const as it does not change the state of the timer
     **/
    double stop() const { 
      return start_time>0?stopWatch::timer()-start_time:start_time ;
    }
  } ; // end of class stopWatch
  

  // This object is responsible for counting time and events that occur for
  // some entity over time.
  class timeAccumulator {
#ifdef PROFILER
    double elapsed_time ;
    size_t elapsed_events ;
#endif
  public:
    timeAccumulator() {
#ifdef PROFILER
      elapsed_time = 0 ;
      elapsed_events = 0 ;
#endif
    }
    void addTime(double time, size_t num_events) {
#ifdef PROFILER
      elapsed_time += time ;
      elapsed_events += num_events ;
#endif
    }
    double getTime() const {
#ifdef PROFILER
      return elapsed_time ;
#else
      return 0 ;
#endif
    }
    size_t getEvents() const {
#ifdef PROFILER
      return elapsed_events ;
#else
      return 0 ;
#endif
    }
  } ;

  enum allocEventType { ALLOC_CREATE, ALLOC_DELETE } ;
  enum executeEventType { EXEC_COMMUNICATION, EXEC_COMPUTATION, EXEC_CONTROL } ;

  // This object is the abstract base class for a data collator.  Currently
  // we are only collecting timing information via the accumulateTime method.
  // However, collecting other data (such as floating point operations or
  // cache misses) can be accomplished by adding additional pure virtual
  // methods to this base class (and all data collator objects.)
  // For an example of a data collator object, inspect the collectTiming
  // object defined in System/scheduler.cc
  class collectData {
  protected:
    std::vector<std::string> groups ;
  public:
    collectData() {}
    virtual ~collectData() {}
    // Create a group for categorizing objects.  Note, groups can form
    // a hierarchy, so we can open multiple levels of groups for hierarchical
    // classification of data
    int openGroup(std::string groupName) {
      groups.push_back(groupName) ;
      return groups.size() ;
    }
    // every group needs to be closed
    void closeGroup(int id) {
      warn(int(groups.size()) != id) ;
      if(int(groups.size()) != id) { // If open and close don't align, try to
        // make things right by poping off groups until this one
        while(id > int(groups.size()))
          groups.pop_back() ;
      }
      // Check for empty which may occur if openGroup does not match closeGroup
      if(!groups.empty()) 
        groups.pop_back() ; // remove group from list
    }
    // Pure virtual method which will be called by every entity that has
    // timing information.  Events are classified as either computation
    // communication, or control
    virtual void accumulateTime(const timeAccumulator &ta, executeEventType t, string eventName) = 0 ;
    virtual void accumulateMemory(const std::string &var,
                                  allocEventType t,
                                  double maxMallocMemory,
                                  double maxBeanMemory) = 0 ;
    virtual void accumulateSchedMemory(const std::string& eventName,
                                       double bytes) = 0;
    virtual void accumulateDCM(const std::string& eventName,
                               long_long l1_dcm, long_long l2_dcm) = 0;
  } ;
  
    
    
  class fact_db ;
  class sched_db ;
  class execute_modules : public CPTR_type {
  public:
    execute_modules() {}
    virtual ~execute_modules() {}
    virtual void execute(fact_db &facts, sched_db &scheds) = 0 ;
    virtual void execute_kernel(const sequence&) {}
    virtual void execute_prelude(const sequence&) {}
    virtual void execute_postlude(const sequence&) {}
    virtual void Print(std::ostream &s) const = 0 ;
    virtual string getName() = 0;
    virtual void dataCollate(collectData &data_collector) const = 0 ;
  } ;

  typedef CPTR<execute_modules> executeP ;

  class execute_list : public execute_modules {
    std::vector<executeP> elist ;
  public:
    execute_list() { }
    virtual ~execute_list() {}
    // Execute object code
    virtual void execute(fact_db &facts, sched_db &scheds) ;
    // Print schedule to stream s
    virtual void Print(std::ostream &s) const ;
    // Collect accumulated data
    virtual void dataCollate(collectData &data_collector) const ;
    // get module name
    virtual string getName() { return "execute_list";};
    // execute_list methods
    void append_list(const executeP &emodule) 
    { if(emodule != 0) elist.push_back(emodule) ; }
    int size() const { return elist.size() ; }
  };
    
  class execute_sequence : public execute_modules {
    std::vector<executeP> elist ;
  public:
    execute_sequence() { }
    virtual ~execute_sequence() {}
    // Execute object code
    virtual void execute(fact_db &facts, sched_db& scheds) ;
    // Print schedule to stream s
    virtual void Print(std::ostream &s) const ;
    // collect accumulated data
    virtual void dataCollate(collectData &data_collector) const ;
    // get name of module
    virtual string getName() { return "execute_sequence";};
    // execute_sequence methods:
    void append_list(const executeP &emodule) 
    { if(emodule != 0) elist.push_back(emodule) ; }
    int size() const { return elist.size() ; }
  };

}
#endif
