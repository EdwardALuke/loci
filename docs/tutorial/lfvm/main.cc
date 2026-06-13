#include <Loci.h>
#include <istream>
#include <iostream>
#include <string>
#include <fstream>
#include "prototypes.h"

int main(int argc, char *argv[])
{
  using std::cerr ;
  using std::endl ;
  
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
    }
    else {
      std::cerr << "argument " << argv[i] << " is not understood."
		<< std::endl;
    }
  }      

  char *grid_filename = "example_grid" ;
  if(argc >= 2)
    grid_filename = argv[1] ;

  fact_db facts ;


  std::cout << "reading grid file " << grid_filename << std::endl ;
  std::ifstream gridf(grid_filename,std::ios::in) ;
  read_triangles(gridf,facts) ;

  setup_edges(facts) ;

  param<double> T_initial,heat_capacity,conductivity, density, bc_heat_density ;
  double thickness = 0.001 ;
  *T_initial = 250 ; // Kelvin
  *bc_heat_density = 10.0/0.0275 ; // 10 watts over 0.0275 meters 
  *heat_capacity = 24.435/63.483*1000.0 ; // J/Kg/K
  *conductivity = 401.0*thickness ; // Watts/m/K
  *density = 8.92*1000.0*thickness ; // kg/m^2

  facts.create_fact("T_initial",T_initial) ;
  facts.create_fact("bc_heat_density",bc_heat_density) ;
  facts.create_fact("heat_capacity",heat_capacity) ;
  facts.create_fact("conductivity",conductivity) ;
  facts.create_fact("density",density) ;

  param<int> max_iteration ;
  *max_iteration = 10000 ; 
  param<double> max_time ;
  *max_time = 1e30 ; // Simulate for 1 second
  
  facts.create_fact("max_iteration",max_iteration) ;
  facts.create_fact("max_time",max_time) ;

  std::ifstream infile("heat.vars",std::ios::in) ;
  if(!infile.fail()) {
    std::cout << "reading heat.vars"<< std::endl ;
    facts.read(infile) ;
  }
  std::cout << "max_iteration="<<*max_iteration << std::endl;

  // Setup boundary conditions

  constraint heat_flux_boundary, tmp ;
  tmp = facts.get_fact("BC_top") ;
  *heat_flux_boundary = *tmp ;
  
  facts.create_fact("heat_flux_boundary",heat_flux_boundary) ;

  entitySet Tset_BC ;
  tmp = facts.get_fact("BC_bottom") ;
  Tset_BC += *tmp ;
  tmp = facts.get_fact("BC_left") ;
  Tset_BC += *tmp ;
  tmp = facts.get_fact("BC_right") ;
  Tset_BC += *tmp ;

  constraint set_temperature_boundary ;
  set_temperature_boundary = Tset_BC ;

  facts.create_fact("set_temperature_boundary",set_temperature_boundary) ;

  std::cout << "heat flux boundary nodes = "<< *heat_flux_boundary<<std::endl ;
  std::cout << "temperature prescribed boundary nodes = " <<
    *set_temperature_boundary << std::endl ;
  

  // Create a rule database called rdb
  rule_db rdb ;
  // Add all of the rules that were inserted into the global_rule_list
  // by register_rule<> types into the rule database rdb
  rdb.add_rules(global_rule_list) ;

  int num_procs = Loci::MPI_processes ;
  int myid = Loci::MPI_rank ;

  // Do some output if processor 0
  if (myid == 0) {
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
  if(!Loci::makeQuery(rdb,facts,query)) {
    cerr << "query failed!" << endl ;
    Loci::Abort() ;
  }


  if(query == "solution") {
    store<double> sol(facts.get_fact("solution")) ;
    std::cout << "sol.domain() = " << sol.domain() << std::endl ;
    std::ofstream ofile("res.2dgv",std::ios::out) ;
    if(ofile.fail()) {
      std::cerr << "unable to open 'res.2dgv'" << std::endl ;
      exit(-1) ;
    }
    ofile << "general"<< std::endl ;
      
    store<vector2d<double> > pos(facts.get_fact("pos")) ;
    entitySet nodes = pos.domain() ;
    Map cl(facts.get_fact("cl")) ;
    Map cr(facts.get_fact("cr")) ;
    MapVec<2> edge_nodes(facts.get_fact("edge_nodes")) ;
    entitySet faces = edge_nodes.domain() ;
    constraint boundary_edges(facts.get_fact("boundary_edges")) ;
    entitySet boundaries = *boundary_edges ;
    entitySet interior_faces = faces - boundaries ;
    entitySet cells = cl.image(faces) + cr.image(interior_faces) ;
    
    
    ofile << nodes.size() << ' ' << nodes.Min() << std::endl ;
    
    entitySet::const_iterator ei ;
    for(ei=nodes.begin();ei!=nodes.end();++ei)
      ofile << pos[*ei] << std::endl ;
    ofile << faces.size() << ' ' << faces.Min() << ' ' ;
    ofile << cells.size() << ' ' << cells.Min() << std::endl ;
    for(ei=interior_faces.begin();ei!=interior_faces.end();++ei) 
      ofile << edge_nodes[*ei][0] << ' '
            << edge_nodes[*ei][1] << ' '
            << cl[*ei] <<' ' << cr[*ei]
            << std::endl ;
    boundaries = faces-interior_faces ;
    for(ei=boundaries.begin();ei!=boundaries.end();++ei)
      ofile << edge_nodes[*ei][0] << ' '
            << edge_nodes[*ei][1] << ' '
            << cl[*ei] << " -1" << std::endl ;
    
    for(ei=nodes.begin();ei!=nodes.end();++ei)
      ofile << sol[*ei] << std::endl ;
  } else {
    Loci::storeRepP query_var = facts.get_fact(query) ;
    std::cout << query << " = " << std::endl ;
    query_var->Print(std::cout) ;
  }
  
  Loci::Finalize() ;
  return 0 ;
}
