#include <Loci.h>
#include <iostream>

using namespace std ;

// Generate a 1d grid over the number line from [0,1] consiting of N segments
void generate_grid(fact_db &facts, int N) {
  // Allocate the nodes and cells of the grid
  entitySet nodes = facts.get_allocation(N+1) ;
  entitySet cells = facts.get_allocation(N) ;

  // setup x coordinates for nodes
  store<float> x ;
  x.allocate(nodes) ;
  float dx = 1./float(N) ; // Uniform delta x
  entitySet::const_iterator ni ;
  float xtmp = 0 ;
  for(ni=nodes.begin();ni!=nodes.end();++ni) {
    x[*ni] = xtmp ;
    xtmp += dx ;
  }
  // Add node positions to facts
  facts.create_fact("x",x) ;

  // Find the nodes that are on the left and right side of cells
  // by shifting the allocated numberings to the left or right by one
  entitySet left_nodes = (nodes >> 1) & nodes ;
  entitySet right_nodes = (nodes << 1) & nodes ;

  // Allocate maps for the left cell and right cell of a node
  Map cl,cr,il,ir ;
  cl.allocate(left_nodes) ;
  cr.allocate(right_nodes) ;
  il.allocate(cells) ;
  ir.allocate(cells) ;
  entitySet::const_iterator ci ;
  // Assign left nodes to cells in consecutive order
  ci = cells.begin() ;
  for(ni=left_nodes.begin();ni!=left_nodes.end();++ni,++ci) {
    cl[*ni] = *ci ;
    ir[*ci] = *ni ;
  }
  // Assign right nodes to cells in consecutive order
  ci = cells.begin() ;
  for(ni=right_nodes.begin();ni!=right_nodes.end();++ni,++ci) {
    cr[*ni] = *ci ;
    il[*ci] = *ni ;
  }

  facts.create_fact("cl",cl) ;
  facts.create_fact("cr",cr) ;
  facts.create_fact("il",il) ;
  facts.create_fact("ir",ir) ;
  
  constraint geom_cells ;
  *geom_cells = cells ;
  facts.create_fact("geom_cells",geom_cells) ;

  // Identify boundary conditions
  constraint left_boundary ;
  constraint right_boundary ;
  *right_boundary = cl.domain() - cr.domain() ;
  *left_boundary = cr.domain() - cl.domain() ;

  facts.create_fact("left_boundary",left_boundary) ;
  facts.create_fact("right_boundary",right_boundary) ;
}


int main(int argc, char *argv[])
{
  Loci::Init(&argc,&argv) ;
    
   
  set_fpe_abort() ;
  
  // query for this variable by default
  std::string query = "solution" ;
  // do not write out facts by default
  bool write_facts = false;
  // do not write out rules by default
  bool write_rules = false;

  // Parse out the LOCI specific arguments.
  //	-q queryvar	make alternative queries
  //	-fact		print out the fact database during program execution
  //	-rule		print out the rule database during program execution
  for (int i=1; i<argc; i++) {
    if (!strcmp(argv[i], "-q") && (i+1) < argc) {
      query = argv[i+1];
      for (int j=0; j<argc-i-2; j++) {
	argv[i+j] = argv[i+j+2];
      }
      argc -= 2;
      i--;
    }
    else if (!strcmp(argv[i], "-fact")) {
      write_facts = !write_facts;
      for (int j=0; j<argc-i-1; j++) {
	argv[i+j] = argv[i+j+1];
      }
      argc -= 1;
      i--;
    }
    else if (!strcmp(argv[i], "-rule")) {
      write_rules = !write_rules;
      for (int j=0; j<argc-i-1; j++) {
	argv[i+j] = argv[i+j+1];
      }
      argc -= 1;
      i--;
    } else {
      std::cerr << "argument " << argv[i] << " is not understood."
		<< std::endl;
    }
  }      

  //-----------------------------------------------------------------
  // Create Rule Database
  //-----------------------------------------------------------------
  
  rule_db rdb ;
  rdb.add_rules(global_rule_list) ;

  //-----------------------------------------------------------------
  // Create Fact Database
  //-----------------------------------------------------------------
  fact_db facts ;

  string varsFile = "heat.vars" ;
  facts.read_vars(varsFile,rdb) ;
  
  param<int> Nin ;
  Nin = facts.get_variable("N") ;
  
  const int N = *Nin ; // Number of points in grid.

  generate_grid(facts,N) ;

  if (Loci::MPI_processes == 1) {
    // Write out the initial fact database if -fact
    if (write_facts) {
      char fn_facts[11];
      strcpy(fn_facts, "facts_init");
      std::ofstream of_facts(fn_facts);

      facts.write(of_facts);

      of_facts.close();
    }

    // Write out the rule database if -rule
    if (write_rules) {
      char fn_rules[6];
      strcpy(fn_rules, "rules");
      std::ofstream of_rules(fn_rules);

      Loci::ruleSet rset = rdb.all_rules();
      Loci::ruleSet::const_iterator itr;
      for (itr = rset.begin(); itr!=rset.end(); itr++) {
	of_rules << *itr << std::endl;
      }

      of_rules.close();
    }
  }

  // Query Loci for fact derived fact 'solution'
//  cout<<"\nQuery Started\n";
  if(!Loci::makeQuery(rdb,facts,query)) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }
//  cout<<"\nQuery Ended\n";
  
  //-----------------------------------------------------------------
  // Final Step: Query the database for solution:
  //-----------------------------------------------------------------

  if(query == "solution") {
    store<float> usol ;
    usol = facts.get_variable("solution") ;

    cout << "The solution is : " <<endl;
    entitySet::const_iterator ei ;
    entitySet dom = usol.domain() ;
    for(ei=dom.begin();ei!=dom.end();++ei) 
      cout << ""<< *ei<<" "<<usol[*ei]<<endl ;
  } else {
    // If we queried for something else, print out the results
    using Loci::variableSet ;
    variableSet queries = variableSet(Loci::expression::create(query)) ;
    for(variableSet::const_iterator vi=queries.begin();vi!=queries.end();++vi){
      Loci::storeRepP rep = facts.get_variable(*vi) ;
      cout << "variable " << *vi << "= " << endl ;
      rep->Print(cout) ;
    }
  }
   
  //-----------------------------------------------------------------
  // End of computations::
  //-----------------------------------------------------------------

  Loci::Finalize() ;
  return 0;

}
//*******************************************************************
