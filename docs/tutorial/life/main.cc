#include <Loci.h>
#include <istream>
#include <iostream>

inline int positive_mod(int v,int m) {
    if(v>=0)
	return v%m ;
    if(v<0)
	return m+(v%m) ;
}

inline int array_compute(int i, int j, int ni, int nj, int base) 
{
  return base+positive_mod(j,nj)+positive_mod(i,ni)*nj ;
}

void create_torus(int ni, int nj, fact_db &facts) {
  int total_entities = ni*nj ;
  entitySet alloc = facts.get_allocation(total_entities) ;
  Map N,S,E,W ;
  N.allocate(alloc) ;
  S.allocate(alloc) ;
  E.allocate(alloc) ;
  W.allocate(alloc) ;
  int base = alloc.Min() ;
  for(int i=0;i<ni;++i) {
    for(int j=0;j<nj;++j) {
      int e = array_compute(i,j,ni,nj,base) ;
      N[e] = array_compute(i,j-1,ni,nj,base) ;
      S[e] = array_compute(i,j+1,ni,nj,base) ;
      E[e] = array_compute(i+1,j,ni,nj,base) ;
      W[e] = array_compute(i-1,j,ni,nj,base) ;
    }
  }

  facts.create_fact("N",N) ;
  facts.create_fact("S",S) ;
  facts.create_fact("E",E) ;
  facts.create_fact("W",W) ;
}

void setup_facts(fact_db &facts) {
  create_torus(3,3,facts) ;

  Map N,S,E,W ;
  N=facts.get_variable("N") ;
  S=facts.get_variable("S") ;
  E=facts.get_variable("E") ;
  W=facts.get_variable("W") ;
  entitySet life_dom = N.domain()&S.domain()&E.domain()&W.domain() ;

  store<int> initial_life ;
  initial_life.allocate(life_dom) ;

  for(entitySet::const_iterator ei=life_dom.begin();
      ei!=life_dom.end();
      ++ei) {
    if(drand48()>0.5)
      initial_life[*ei] = 1 ;
    else
      initial_life[*ei] = 0 ;
  }
  facts.create_fact("initial_life",initial_life) ;
}

int main(int argc, char *argv[])
{
  Loci::Init(&argc,&argv) ;

  set_fpe_abort() ;
  
  string query = "solution" ;

  while(argc>=2 && argv[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(argc >= 2 && !strcmp(argv[1],"-q")) {
      query = argv[2] ;
      argc -= 2 ;
      argv += 2 ;
    } else {
      cerr << "argument " << argv[1] << " is not understood." << endl ;
      argc-- ;
      argv++ ;
    }
  }      

  fact_db facts ;

  setup_facts(facts) ;

  param<int> max_iteration ;
  *max_iteration = 10 ;
  facts.create_fact("max_iteration",max_iteration) ;

  // Create a rule database called rdb
  rule_db rdb ;
  // Add all of the rules that were inserted into the global_rule_list
  // by register_rule<> types into the rule database rdb
  rdb.add_rules(global_rule_list) ;

  
  int num_procs = Loci::MPI_processes ;
  int myid = Loci::MPI_rank ;
  
  std::vector<entitySet> partition = Loci::generate_distribution(facts,rdb) ;
  Loci::distribute_facts(partition, facts, rdb) ;

  executeP schedule = create_execution_schedule(rdb,facts,query) ;

  if(schedule == 0) {
    cerr << "unable to produce execution schedule to satisfy query for "
         << query << endl ;
  } else {
    // Save the schedule in the file .schedule for reference
    ostringstream oss ;
    oss << ".schedule" ;

    if(num_procs > 1) {
      oss << "-" << myid ;
    }
    string sched_filename = oss.str() ;
    ofstream sched_file(sched_filename.c_str(),ios::out) ;
    schedule->Print(sched_file) ;
    sched_file.close() ;

    // execute schedule
    schedule->execute(facts) ;

    Loci::storeRepP query_var = facts.get_fact(query) ;
    cout << query << " = " << endl ;
    query_var->Print(cout) ;
  }
    
  Loci::Finalize() ;
  return 0 ;
}
