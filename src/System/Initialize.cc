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
#include <signal.h>

#ifdef USE_PETSC
#include "petsc.h"
#include <petscerror.h>
#include <petscvec.h>
#include <petscksp.h>

#if ((PETSC_VERSION_MAJOR > 3) || (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR > 2))
#define PETSC_33_API
#endif
#if ((PETSC_VERSION_MAJOR > 3) || (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR > 6))
#define PETSC_37_API
#endif

#if ((PETSC_VERSION_MAJOR > 3) || (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR > 16))
#define PETSC_317_API
#endif



// Force key petsc functions to load. (Helps when using static libraries)
void dummyFunctionDependencies(int i) {
  int localSize=0,globalSize=0 ;
  struct _p_Vec *v=0 ;
  struct _p_Mat *m=0 ;
  struct _p_KSP *ksp=0 ;
  int ierr = VecSetSizes(v,localSize,globalSize) ;
#ifdef PETSC_33_API
  ierr = MatCreateAIJ(MPI_COMM_WORLD,0,0,0,0,0,0,0,0,&m) ;
#else
  ierr = MatCreateMPIAIJ(MPI_COMM_WORLD,0,0,0,0,0,0,0,0,&m) ;
#endif
  ierr = KSPSetFromOptions(ksp) ;
  if(i==0) 
    dummyFunctionDependencies(ierr) ;
}
#endif

#ifdef USE_PAPI
#include <papi.h>
#include <papiStdEventDefs.h>
#endif

#include <rule.h>
#include <keyspace.h>
#include <mod_db.h>
#include "dist_tools.h"
#include "loci_globs.h"
#include "comp_tools.h"
#include "thread.h"
#include <Tools/debug.h>
#include <Tools/except.h>
#include "docreport.h"
#include <new>
using std::bad_alloc ;

#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>
#include <unistd.h>

#include "Tools/debugger.h"

#include <mpi.h>

#define SIMPLE_SPRNG
#define USE_MPI

#include <sprng.h>

#include <iostream>
#include <fstream>
#include <vector>

#ifdef PAPI_DEBUG
#include <papi.h>
#endif

using std::cout ;
using std::cerr ;
using std::endl ;
using std::ios ;
using std::ifstream ;
using std::ostream ;
using std::ofstream ;
using std::string ;
using std::ostringstream ;
using std::vector ;

namespace Loci {
  MPI_Datatype MPI_FADD ;
  MPI_Op MPI_FADD_SUM ;
  MPI_Op MPI_FADD_PROD ;
  MPI_Op MPI_FADD_MIN ;
  MPI_Op MPI_FADD_MAX ;

  MPI_Datatype MPI_FADD2 ;
  MPI_Op MPI_FADD2_SUM ;
  MPI_Op MPI_FADD2_PROD ;
  MPI_Op MPI_FADD2_MIN ;
  MPI_Op MPI_FADD2_MAX ;

  MPI_Datatype MPI_MFADD ;
  MPI_Op MPI_MFADD_SUM ;
  MPI_Op MPI_MFADD_PROD ;
  MPI_Op MPI_MFADD_MIN ;
  MPI_Op MPI_MFADD_MAX ;

  MPI_Info PHDF5_MPI_Info ;
  
  int MPI_processes = 1;
  int MPI_rank = 0 ; 
  int MPI_processes_per_host = 1 ;

  bool useDebugDir = true ;
  bool useDomainKeySpaces = false ;

  int method = 1 ; // Iterative Weighted Static
  bool verbose = false ;
  /////////////////////////////
  // flags to turn on/off the visualization feature
  bool show_graphs = false ;
  bool show_decoration = false ;
  // flag to enable/disable the dynamic memory management
  bool use_dynamic_memory = true ;
  // flag to enable/disable output of dynamic memory
  // and multilevel graph decoration information
  bool show_dmm_verbose = false ;
  // flag to enable/disable chomping
  bool use_chomp = true ;
  // flag to visualize the chomping graph
  bool show_chomp = false ;
  // flag to turn on the summary report on chomping
  bool chomp_verbose = false ;
  // flag to enable outputing schedule to file
  bool schedule_output = false ;
  // flag to enable dynamic scheduling feature
  bool use_dynamic_scheduling = true ;
  // chomping size(unit is KB), default is 128KB
  int chomping_size = 128 ;
  // flag to enable memory profiling
  bool profile_memory_usage = false ;
  // flag to enable memory information gathering.
  // this is different than the "profile_memory_usage" flag.
  // this flag reports the memory information using the
  // "data collector" framework.
  bool collect_memory_info = false ;
  // flag to use a different scheduler
  // (more memory conservative and more synchronization
  //  points will be generated)
  bool memory_greedy_schedule = true ;
  
  //use randomized greedy DAG scheduler. Firs the sizes of variables are queried.
  //Then the scheduler is rerun using the sizes as weight.
  bool randomized_memory_greedy_schedule = false;
  
  // flag to disable memory deallocation decoration
  // if use_dynamic_memory==true. this is more of a
  // debugging tool. in this mode, we only have
  // memory allocation and no memory deallocation
  bool dmm_no_deallocation = false ;

  // flag to enable work duplication techniques
  bool duplicate_work = false;
  bool multilevel_duplication = false;
  bool reduction_duplication = false;
  bool pointwise_duplication = false;
  bool extended_duplication = false;
  bool collect_timings = false;
  double time_duration_to_collect_data = 0 ;
  bool use_duplicate_model = false;
  bool use_simple_partition=false;

#ifdef H5_HAVE_PARALLEL
  bool use_parallel_io = false ; //true ;
#else
  bool use_parallel_io = false ;
#endif
  
  // space filling curve partitioner
#ifdef LOCI_USE_METIS
  bool use_sfc_partition = false ;
#else 
  // RSM COMMENT 20181108 setting this to true,
  // automatically disables calls to METIS decomposition
  // when we add support for ZOLTAN this will need to be 
  // REVISITED: METIS_DISABLE_NOTE
  bool use_sfc_partition=true;
#endif /* ifndef LOCI_USE_METIS */
  bool use_orb_partition = false ;

  int metis_cpp_threshold = 384 ;
  
  extern int factdb_allocated_base ;

  string PFS_Script ; // Parallel File System Striping Script

  char * model_file;
  // these flags are used to indicate additional weights
  // for mesh cells during the initial partitioning stage
  // Note: this is currently only designed for the ParMETIS
  // partitioning routine
  bool load_cell_weights = false ;
  string cell_weight_file ;
  
  // flag to indicate whether multithreading is used for each type of rules
  bool threading_pointwise = false;
  bool threading_global_reduction = false;
  bool threading_local_reduction = false;
  bool threading_chomping = false;
  bool threading_recursion = false;
  int num_threads = 0;
  int num_thread_blocks = 10; // default 10 blocks per threads


  ofstream debugout ;
  ofstream timeout;

  double total_memory_usage = 0 ;

  extern int current_rule_id ;

  size_t MPI_process_mem_avail() {
     size_t mem = 0 ;
     ifstream memfile("/proc/meminfo") ;
     if(memfile.fail())
       return mem ;
     string key="MemAvailable:" ;
     while(!memfile.fail() || !memfile.eof()) {
       string type ;
       int val ;
       string unit ;
       memfile >> type >> val ;
       std::getline(memfile,unit) ;
       if(type==key) {
	 mem = val ;
	 mem *= 1024 ;
	 mem /= MPI_processes_per_host ;
	 return mem ;
       }
     }
     return mem ;
  }

  void disableDebugDir() {useDebugDir = false ;}

  void debug_print_rule() {
    if(current_rule_id != 0) {
      rule r(current_rule_id) ;
      current_rule_id = 0 ;
      cerr << "crash occured in rule " << r << endl ;

      Loci::debugout << "crash occured in rule " << r << endl ;

      if(verbose && exec_current_fact_db != 0) {
        char buf[128] ;
        snprintf(buf,128,"crash_dump.%d",MPI_rank) ;
        ofstream cfile(buf,ios::out) ;
        cfile << "rule: " << r << endl ;

        variableSet v = r.sources() ;

        cfile << "facts = {" << endl ;
        for(variableSet::const_iterator vi=v.begin();vi!=v.end();++vi) {
          cfile << *vi << ':' ;
          storeRepP p = Loci::exec_current_fact_db->get_variable(*vi);
          p->Print(cfile) ;
        }
        cfile << "}" << endl ;
      }
    }
  }
}

extern "C" {
  void MPI_errors_reporter(MPI_Comm *comm, int *err,...) {
    int len = 1024 ;
    char error_string[1024] ;
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank) ;
    MPI_Error_string(*err,error_string ,&len) ;
    cerr << "MPI Error: Rank=" << rank
         << ", Error='"<< error_string << "'" << endl ;
    Loci::Abort() ;
  }

  void TerminateSignal(int sig) {
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank) ;
    if(rank == 0) {
      cerr << "Processor 0 recieved a terminate signal" << endl ;
    }
    if(Loci::current_rule_id != 0) {
      Loci::rule r(Loci::current_rule_id) ;
      Loci::current_rule_id = 0 ;
      Loci::debugout << "Terminate occured in rule " << r << endl ;
    }
    MPI_Abort(MPI_COMM_WORLD,-1) ;
  }
}

namespace Loci {
  extern   void register_closing_function(void (*fptr)(int code)) ;

#ifdef USE_CUDA_RT
  extern int setCudaDevice() ;
#endif

  void closeoutMPI(int code) {
    if(code == -1) {
      int MPI_rank ;
      MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank) ;
      cerr << "program failed on processor " << MPI_rank << endl ;
    }
    if(code == -2) {
      MPI_Abort(MPI_COMM_WORLD,-1) ;
    }
    
  }

  void sumFADd(FADd *rin, FADd *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] += rin[i] ;
  }
  void prodFADd(FADd *rin, FADd *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] *= rin[i] ;
  }
  void maxFADd(FADd *rin, FADd *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] = max(rinout[i],rin[i]) ;
  }
  void minFADd(FADd *rin, FADd *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] = min(rinout[i],rin[i]) ;
  }

  void sumFAD2d(FAD2d *rin, FAD2d *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] += rin[i] ;
  }
  void prodFAD2d(FAD2d *rin, FAD2d *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] *= rin[i] ;
  }
  void maxFAD2d(FAD2d *rin, FAD2d *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] = max(rinout[i],rin[i]) ;
  }
  void minFAD2d(FAD2d *rin, FAD2d *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] = min(rinout[i],rin[i]) ;
  }

  void sumMFADd(MFADd *rin, MFADd *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] += rin[i] ;
  }
  void prodMFADd(MFADd *rin, MFADd *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] *= rin[i] ;
  }
  void maxMFADd(MFADd *rin, MFADd *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] = max(rinout[i],rin[i]) ;
  }
  void minMFADd(MFADd *rin, MFADd *rinout, int *len, MPI_Datatype *dtype) {
    for(int i=0;i<*len;++i)
      rinout[i] = min(rinout[i],rin[i]) ;
  }

  MPI_Errhandler Loci_MPI_err_handler ;
  

  //This is the first call to be made for any Loci program be it
  //sequential or parallel.
  void Init(int* argc, char*** argv)  {
    char *execname = (*argv)[0] ;
    const char *hostname = "localhost" ;
    const char *debug = "gdb" ;
    //Setting up of the global variables for processor ID and the
    //total number of processes.
#ifdef USE_PETSC
    PetscInitialize(argc,argv,(char*)0,(char*)0) ;
#ifdef PETSC_37_API
    PetscOptionsSetValue(0,"-options_left","false") ;
#else
    PetscOptionsSetValue("-options_left","false") ;
#endif
    PetscPopErrorHandler() ;
#ifdef PETSC_317_API
    PetscPushErrorHandler(PetscIgnoreErrorHandler,PETSC_NULLPTR) ;
#else
    PetscPushErrorHandler(PetscIgnoreErrorHandler,PETSC_NULL) ;
#endif
#else
    MPI_Init(argc, argv) ;
#endif

    
    create_mpi_info(&PHDF5_MPI_Info) ;
#ifndef MPI_STUBB
    {
      // Create FADD type
      int count = 2 ;
      int blocklens[] = {1,1} ;
      MPI_Aint indices[2] ;
      FADd tmp ;
      indices[0] = (MPI_Aint)((char *) &(tmp.value) - (char *) &tmp) ;
      indices[1] = (MPI_Aint)((char *) &(tmp.grad) - (char *) &tmp) ;
      MPI_Datatype typelist[] = {MPI_DOUBLE,MPI_DOUBLE} ;

      MPI_Type_create_struct(count,blocklens,indices,typelist,&MPI_FADD) ;
      MPI_Type_commit(&MPI_FADD) ;

      MPI_Op_create((MPI_User_function *)sumFADd,1,&MPI_FADD_SUM) ;
      MPI_Op_create((MPI_User_function *)prodFADd,1,&MPI_FADD_PROD) ;
      MPI_Op_create((MPI_User_function *)maxFADd,1,&MPI_FADD_MAX) ;
      MPI_Op_create((MPI_User_function *)minFADd,1,&MPI_FADD_MIN) ;
    }

    {
      // Create FADD type
      int count = 3 ;
      int blocklens[] = {1,1,1} ;
      MPI_Aint indices[3] ;
      FAD2d tmp ;
      indices[0] = (MPI_Aint)((char *) &(tmp.value) - (char *) &tmp) ;
      indices[1] = (MPI_Aint)((char *) &(tmp.grad) - (char *) &tmp) ;
      indices[2] = (MPI_Aint)((char *) &(tmp.grad2) - (char *) &tmp) ;

      MPI_Datatype typelist[] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE} ;
      MPI_Type_create_struct(count,blocklens,indices,typelist,&MPI_FADD2) ;
      MPI_Type_commit(&MPI_FADD2) ;

      MPI_Op_create((MPI_User_function *)sumFAD2d,1,&MPI_FADD2_SUM) ;
      MPI_Op_create((MPI_User_function *)prodFAD2d,1,&MPI_FADD2_PROD) ;
      MPI_Op_create((MPI_User_function *)maxFAD2d,1,&MPI_FADD2_MAX) ;
      MPI_Op_create((MPI_User_function *)minFAD2d,1,&MPI_FADD2_MIN) ;

    }

    {
      // Create MFADD type
      //int blocklens[] = {1,1,1,1,1,1,1,1,1,1,1} ;
      int blocklens[MFAD_SIZE+1] ;
      for (int i=0;i<MFAD_SIZE+1;i++) blocklens[i] = 1 ;
      MPI_Aint indices[MFAD_SIZE+1] ;
      MFADd tmp ;
      indices[0] = (MPI_Aint)((char *) &(tmp.value) - (char *) &tmp) ;
      for (int i=1;i<MFAD_SIZE+1;i++) indices[i] = (MPI_Aint)((char *) &(tmp.grad[i-1]) - (char *) &tmp) ;

      MPI_Datatype typelist[MFAD_SIZE+1] ;
      for (int i=0;i<MFAD_SIZE+1;i++) typelist[i] = MPI_DOUBLE ;
      MPI_Type_create_struct(MFAD_SIZE+1,blocklens,indices,typelist,&MPI_MFADD) ;
      MPI_Type_commit(&MPI_MFADD) ;

      MPI_Op_create((MPI_User_function *)sumMFADd,1,&MPI_MFADD_SUM) ;
      MPI_Op_create((MPI_User_function *)prodMFADd,1,&MPI_MFADD_PROD) ;
      MPI_Op_create((MPI_User_function *)maxMFADd,1,&MPI_MFADD_MAX) ;
      MPI_Op_create((MPI_User_function *)minMFADd,1,&MPI_MFADD_MIN) ;

    }
#endif
    
    time_duration_to_collect_data = MPI_Wtick()*20;
    signal(SIGTERM,TerminateSignal) ;

    //    MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN) ;
    //    MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_ARE_FATAL) ;

    MPI_Comm_create_errhandler(MPI_errors_reporter,&Loci_MPI_err_handler) ;
    MPI_Comm_set_errhandler(MPI_COMM_WORLD,Loci_MPI_err_handler) ;
    //    MPI_Errhandler_create(&MPI_errors_reporter,&err_handler) ;
    //    MPI_Errhandler_set(MPI_COMM_WORLD,err_handler) ;

    MPI_Comm_size(MPI_COMM_WORLD, &MPI_processes) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank) ;

    int sprng_seed = 985456376 ;
    int sprng_gtype = SPRNG_LFG ; // sprng generator type
    // LFG   - 0 - Additive Lagged Fibonaci Generator
    // LCG   - 1 - Linear Congruential Generator
    // LCG64 - 2 - Linear Congruential Generator (64 bit)
    // CMRG  - 3 - Combined Multiple Recursive Generator
    // MLFG  - 4 - Multiplicative Lagged Fibonacci Generator

    char *p = 0 ;
    if((p=getenv("LOCI_PFS_SCRIPT")) != 0) {
      PFS_Script = string(p) ;
    } else {
      PFS_Script = string ("") ;
    }
    try {

      ostringstream oss ;
      char *p = 0 ;
      if((p = getenv("LOCI_MODULE_PATH")) == 0)
        p = getenv("LD_LIBRARY_PATH") ;
      if(p != 0) {
        string pathlist = string(p) ;
        string path ;
        for(string::const_iterator si = pathlist.begin();
            si!= pathlist.end();
            ++si) {
          if(*si != ':')
            path += *si ;
          else {
            if(path != "")
              AddModuleSearchDir(path) ;
            path = "" ;
          }
        }
        if(path != "")
          AddModuleSearchDir(path) ;
      }

#ifdef LOCI_RPATH
      { string rpath = LOCI_RPATH ; AddModuleSearchDir(rpath) ; }
#endif

      //    if(GLOBAL_OR(stat("output",&statbuf)!=0)) {
      //      if(MPI_rank == 0)
      //        mkdir("output",0755) ;
      //    }

      // All the rules in an unnamed namespace are first copied into the
      // global rule list. To add rules to the rule database we just
      // neeed to use the global_rule_list. Inititally when the rules
      // are registered using the register rule it gets pushed into the
      // register_rule_list which is a static rule list.
      if(!register_rule_list.empty()) {
        global_rule_list.copy_rule_list(register_rule_list) ;
        register_rule_list.clear() ;
      }
#ifdef DYNAMICSCHEDULING
      // do the same to get all the defined keyspace
      if(!register_key_space_list.empty()) {
        global_key_space_list.copy_space_list(register_key_space_list) ;
        register_key_space_list.clear() ;
      }
#endif
      bool debug_setup = false ;
      int i = 1 ;
      int k = 1 ; // copy cursor for removing processed arguments from argv
      while(i<*argc) {
        if(!strcmp((*argv)[i],"--display")) {
          debug_setup = true ;
          hostname = (*argv)[i+1] ;
          i+=2 ;
        } else if(!strcmp((*argv)[i],"--debug")) {
          debug = (*argv)[i+1] ;
          schedule_output = true ;
          i+=2 ;
        } else if(!strcmp((*argv)[i],"--sprng_gen")) {
          // set type of random number generator for sprng
          if(!strcmp((*argv)[i+1],"LFG")) {
            sprng_gtype = SPRNG_LFG ;
          } else if(!strcmp((*argv)[i+1],"LCG")) {
            sprng_gtype = SPRNG_LCG ;
          } else if(!strcmp((*argv)[i+1],"LCG64")) {
            sprng_gtype = SPRNG_LCG64 ;
          } else if(!strcmp((*argv)[i+1],"CMRG")) {
            sprng_gtype = SPRNG_CMRG ;
          } else if(!strcmp((*argv)[i+1],"MLFG")) {
            sprng_gtype = SPRNG_MLFG ;
          } else {
            cerr << "unknown generator type " << (*argv[i+1]) << endl ;
          }
          i+=2 ;
        } else if(!strcmp((*argv)[i],"--sprng_gen_seed")) {
          sprng_seed = make_sprng_seed() ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--verbose")) {
          verbose = true ;
          schedule_output = true ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--scheduleoutput")) {
          schedule_output = true ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--graphs")) {
          show_graphs = true ; // visualize the dependency graph &
          // the decomposed graph & every supernode
          i++ ;
        } else if(!strcmp((*argv)[i],"--decoration")) {
          show_decoration = true ; // visualize the decorated multilevel graph
          i++ ;
	} else if(!strcmp((*argv)[i],"--pio")) { // turn on parallel IO
#ifdef H5_HAVE_PARALLEL
	  if(MPI_processes > 1) 
	    use_parallel_io = true ;
	  else {
	    use_parallel_io = false ;
	    cerr << "Parallel I/O not used when running on a single CPU!" << endl ;
	  }
#else
	  cerr << "Cannot enable parallel I/O, Loci linked with serial HDF5 library!" << endl ;
	  use_parallel_io = false ;
#endif
	  i++ ;
	} else if(!strcmp((*argv)[i],"--nopio")) { // turn off parallel IO
#ifndef H5_HAVE_PARALLEL
	  cerr << "Cannot disable parallel I/O, Loci linked with serial HDF5 library!" << endl ;
#endif
	  i++ ;
	  use_parallel_io = false ;
	} else if(!strcmp((*argv)[i],"--max_gpu_stream")) {
	  int num = atoi((*argv)[i+1]) ;
	  i+=2 ;
	  num = max(1,min(8,num)) ;
#ifdef USE_CUDA_RT
	  MAXGPUStreamAlloc = 1<<num ;
#endif
	} else if(!strcmp((*argv)[i],"--metis_cpp_threshold")) {
            metis_cpp_threshold = atoi((*argv)[i+1]);
	    i+=2;
        } else if(!strcmp((*argv)[i],"--simple_partition")) {
          use_simple_partition = true ; //  partition domain using n/p cuts
          i++ ;
        } else if(!strcmp((*argv)[i],"--orb_partition")) {
          use_orb_partition = true ; // partition mesh using ORB method
          i++ ;
        } else if(!strcmp((*argv)[i],"--sfc_partition")) {
          use_sfc_partition = true ; // partition mesh using SFC method
          i++ ;
        } else if(!strcmp((*argv)[i],"--dmm")) {
          use_dynamic_memory = true ; // use dynamic memory management
          i++ ;
        } else if(!strcmp((*argv)[i],"--dmmnofree")) {
          dmm_no_deallocation = true ; // do not do deallocation
          i++ ;
        } else if(!strcmp((*argv)[i],"--nodmm")) {
          use_dynamic_memory = false ; // use static preallocation
          i++ ;
        } else if(!strcmp((*argv)[i],"--dmmverbose")) {
          show_dmm_verbose = true ; // output some info about dmm
          i++ ;
        } else if(!strcmp((*argv)[i],"--chomp")) {
          use_chomp = true ; // use the chomping scheme
          i++ ;
        } else if(!strcmp((*argv)[i],"--nochomp")) {
          use_chomp = false ; // don't use the chomping scheme
          i++ ;
        } else if(!strcmp((*argv)[i],"--showchomp")) {
          show_chomp = true ; // visualize the chomp graph
          i++ ;
        } else if(!strcmp((*argv)[i],"--chompverbose")) {
          chomp_verbose = true ; // enable chomping summary report
          i++ ;
        } else if(!strcmp((*argv)[i],"--chompingsize")) {
          std::stringstream ss ;
          ss << (*argv)[i+1] ;
          ss >> chomping_size ;
          if(chomping_size > 500*1024*1024/*500MB*/)
            cerr << "WARNING: setting chompingsize too large"
                 << " may crash the program" << endl ;
          i+=2 ;
        } else if(!strcmp((*argv)[i],"--memprofile")) {
          // profiling memory usage at run time
          profile_memory_usage = true ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--meminfo")) {
          collect_memory_info = true ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--memgreedy")) {
          memory_greedy_schedule = true ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--nomemgreedy")) {
          memory_greedy_schedule = false ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--memrandomized")) {
          randomized_memory_greedy_schedule = true ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--nomemrandomized")) {
          randomized_memory_greedy_schedule = false ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--method")) {
          if(!strcmp((*argv)[i+1],"none")) {
            method = 0 ;
          } else  if(!strcmp((*argv)[i+1],"IWS")) {
            method = 1 ;
          } else  if(!strcmp((*argv)[i+1],"FSC")) {
            method = 2 ;
          } else  if(!strcmp((*argv)[i+1],"FAC")) {
            method = 3 ;
          } else  if(!strcmp((*argv)[i+1],"GSS")) {
            method = 4 ;
          } else
            method = atoi((*argv)[i+1]);
          i+=2;
        } else if(!strcmp((*argv)[i],"--balance")) {
          use_dynamic_scheduling = true ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--nobalance")) {
          use_dynamic_scheduling = false ;
          i++ ;
        } else if(!strcmp((*argv)[i],"--duplicate_work")){
          duplicate_work = true;
          reduction_duplication = true;
          pointwise_duplication = true;
          extended_duplication = true;
          i++;
        } else if(!strcmp((*argv)[i],"--no_reduction_duplication")){
          reduction_duplication = false;
          i++;
        } else if(!strcmp((*argv)[i],"--no_pointwise_duplication")){
          pointwise_duplication = false;
          i++;
        } else if(!strcmp((*argv)[i],"--no_extended_duplication")){
          extended_duplication = false;
          i++;
        } else if(!strcmp((*argv)[i],"--multilevel_duplication")){
          multilevel_duplication = true;
          i++;
        } else if(!strcmp((*argv)[i],"--extreme_duplication")){
          duplicate_work = true;
          multilevel_duplication = true;
          reduction_duplication = true;
          pointwise_duplication = true;
          extended_duplication = true;
          i++;
        } else if(!strcmp((*argv)[i],"--collect_timings")){
          collect_timings = true;
          i++;
        } else if(!strcmp((*argv)[i],"--use_duplicate_model")){
          use_duplicate_model = true;
          model_file = (*argv)[i+1];
          i += 2;
          //      } else if(!strcmp((*argv)[i],"--measure_rule_timings")){
          //	measure_rule_timings = true;
          //	i++;
        } else if(!strcmp((*argv)[i],"--load_cell_weights")){
          load_cell_weights = true ;
          cell_weight_file = (*argv)[i+1] ;
          i+=2 ;
	} else if(!(strcmp((*argv)[i],"--test"))) {
	  useDomainKeySpaces = true ;
	  i++ ;
	} else if(!(strcmp((*argv)[i],"--set_4gig_entity_space")) ||
		  !(strcmp((*argv)[i],"--big"))) {
	  factdb_allocated_base = std::numeric_limits<int>::min() + 2048 ;
	  useDomainKeySpaces = true ;
	  i++ ;
        } else if(!strcmp((*argv)[i],"--threads")) {
          // determine the number of threads to use.
          // ideally we would like to detect the available
          // hardware threads on a given platform automatically,
          // this is a work in the future.
          
          // but right now, we will accept the user inputs
          // and check if it is reasonable
          int nt = atoi( (*argv)[i+1]);
          if(nt < 0 || nt > 128) 
            nt = 2;
          num_threads = nt;
          if(num_threads > 1) {
            threading_pointwise = true;
            threading_global_reduction = true;
            threading_local_reduction = true;
            threading_chomping = true;
            threading_recursion = true;
          }
          i+=2;
        } else if(!strcmp((*argv)[i],"--thread_blocks")) {
          // determine the number of thread blocks to use
          int nb = atoi( (*argv)[i+1]);
          if(nb <= 0 || nb > 100)
            nb = 10; // prevent unreasonable values
          num_thread_blocks = nb;
          i+=2;
        } else if(!strcmp((*argv)[i],"--no_threading_pointwise")) {
          threading_pointwise = false;
          i++;
        } else if(!strcmp((*argv)[i],"--no_threading_global_reduction")) {
          threading_global_reduction = false;
          i++;
        } else if(!strcmp((*argv)[i],"--no_threading_local_reduction")) {
          threading_local_reduction = false;
          i++;
        } else if(!strcmp((*argv)[i],"--no_threading_chomp")) {
          threading_chomping = false;
          i++;
        } else if(!strcmp((*argv)[i],"--no_threading_recursion")) {
          threading_recursion = false;
          i++;
        } else {
          // copy items to their new locations after removed items
          // we only need to copy once arguments are processed.
          if(k!=i)
            (*argv)[k] = (*argv)[i] ;
	  i++ ;
          k++ ;
	}
      }

      // reset the argument counter to account for processed flags
      *argc = k ;

      if(useDebugDir) {
        //Create a debug file for each process
        // if output directory doesn't exist, create one
        bool debug_is_directory = true ;
        struct stat statbuf ;
        if(GLOBAL_OR(stat("debug",&statbuf)!=0)) {
          if(MPI_rank == 0)
            mkdir("debug",0755) ;
          for(int i=0;i<1000;++i) {
            if(GLOBAL_AND(stat("debug",&statbuf)==0))
              break ;
          }
        } else {
          if(!S_ISDIR(statbuf.st_mode)) {
            cerr << "file 'debug' should be a directory!, rename 'output' and start again."
                 << endl ;
            debug_is_directory = false ;

          }
        }

        if(debug_is_directory) {
          if(MPI_processes == 1)
            oss << "debug/debug" ;
          else {
	    if((MPI_rank == 0) || verbose || schedule_output) {
	      oss << "debug/debug."<<MPI_rank ;
	    } else {
	      oss << "/dev/null" ;
	    }
	  }
        } else
          oss << "debug."<< MPI_rank ;

        string filename  = oss.str() ;
        debugout.open(filename.c_str(),ios::out) ;
      } else {
        debugout.open("/dev/null",ios::out) ;
      }
    

#ifdef USE_CUDA_RT
      setCudaDevice() ;
#endif
    
      // Find number of mpi processes per host
      {
	long hid = gethostid() ;
	vector<long> host_list(MPI_processes) ;
	MPI_Allgather(&hid,1,MPI_LONG,&host_list[0],1,MPI_LONG,MPI_COMM_WORLD) ;
	int cnt = 0 ;
	for(int i=0;i<MPI_processes;++i)
	  if(hid == host_list[i])
	    cnt++ ;
	MPI_processes_per_host = cnt ;
	debugout << "mpi processes per host = " << cnt << endl ;
      }

      init_sprng(sprng_gtype,sprng_seed,SPRNG_DEFAULT) ;


      if(collect_timings) {
        oss.str("");
        if(MPI_processes == 1)
          oss << "comp_timings";
        else
          oss << "comp_timings-"  << MPI_rank ;
        string filename = oss.str();
        timeout.open(filename.c_str(), ios::out);
      }

      /*
        if(measure_rule_timings) {
        oss.str("");
        if(MPI_processes == 1)
	oss << "rule_timings";
        else
	oss << "rule_timings-"  << MPI_rank ;
        filename = oss.str();
        ruleTimeOut.open(filename.c_str(), ios::out);
        ruleTimeOut << "Output Format:" << endl;
        ruleTimeOut << "===============================================" << endl;
        ruleTimeOut << endl;
        ruleTimeOut << "Rule Name " << endl;
        ruleTimeOut << "Context Size\tTime Taken\tAverage Time" << endl;
        ruleTimeOut << endl;
        ruleTimeOut << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        ruleTimeOut << endl;
        }
      */
      set_debug_callback(debug_print_rule) ;
      if(debug_setup) {
        setup_debugger(execname,debug,hostname) ;
      }
      chopsigs_() ;
    } catch(const BasicException &err) {
      cerr << "Caught exception in Loci::Initialize()"<<endl ;
      err.Print(cerr) ;
    } catch(const bad_alloc &x) {
      cerr << "Out of memory: " << x.what() <<endl ;
      Loci::Abort() ;
    }
    register_closing_function(closeoutMPI) ;

#ifdef PAPI_DEBUG
    // // initialize PAPI
    if( (PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT)
      cerr << "PAPI library init error!" << endl;
    // if(num_threads > 0) {
    //   if (PAPI_thread_init(pthread_self) != PAPI_OK) {
    //     cerr << "PAPI library thread init error!" << endl;
    //     Loci::Abort();
    //   }
    // }
#endif
  }
  //All Loci programs must end with this call.
  extern void call_closing_functions(int code) ;

  void Finalize() {
    call_closing_functions(0) ;
    register_rule_list.clear() ;
    global_rule_list.clear() ;
    rule::rdb_cleanup() ;
#ifdef DYNAMICSHCEDULING
    register_key_space_list.clear() ;
    global_key_space_list.clear() ;
#endif
    //    storeAllocateData.clear() ;
    //    GPUstoreAllocateData.clear() ;
    exec_current_fact_db = 0 ;
#ifdef USE_PETSC
    PetscFinalize() ;
#else
    MPI_Finalize() ;
#endif
  }

  void Abort() {
    debugger_() ;
    MPI_Abort(MPI_COMM_WORLD,-1) ;
  }

  double random() {
    return sprng() ;
  }

  int irandom() {
    return isprng() ;
  }

}
