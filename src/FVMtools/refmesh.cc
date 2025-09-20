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

#include <iostream>
#include <string>
#include <Loci.h>
#include "./FVMAdapt/defines.h"
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ios;

namespace Loci{
  void parallelClassifyCell(fact_db &facts) ;
  void createEdgesPar(fact_db &facts);
}


int main(int argc, char ** argv) {
  // Let Loci initialize itself.
  // This also gives Loci first dibs on the command line arguments.
  Loci::Init(&argc, &argv);
  

  // This is the name of the mesh file that you want to read in.
  string meshfile;
  //This is the name of the refinement plan file
  string planFile;
  //This is the name of the output gridfile
  string outFile = "out.vog";
  string c2pFile;
  string parentPlanFile;//the plan file from former cycle, only need it when cell2parent_output is queried
  
  // Here's where we parse out the command line arguments that are
  // relevant to this program.
  int j=1;
  string pathname = "";
  bool cell2parent = false;
  bool restart = false;
  
  //print out help info
  if( (argc == 1)||(argc==2) ){
    if(Loci::MPI_rank == 0){
      cout<<"command line:" << endl;
      cout <<"refmesh <options> <filename> <options> <filename> ... "<< endl;
      cout << endl;
      cout << "options:" << endl;
      cout <<"-g <file> -- original grid file, refinement plans are based on this grid" << endl;
      cout <<"-r <file> -- input refinement plan file" <<endl;
      cout <<"-pr <file> -- input refinement plan file from former cycle" <<endl;
      cout <<"-o <file> -- output grid file" << endl;
      cout <<"-c2p <file> -- output cell2parent map file" << endl;
    }
    return 0;
  }
  
  for (int i=1; i<argc; i++) {
    // Let's look at the i'th argument
    string arg(argv[i]);
    
    if (arg == "-g" && (i+1) < argc) {
      // Replace the mesh filename with the next argument
      meshfile =  argv[++i];      
    }
    else if(arg == "-o" && (i+1) < argc){
      //replace the output filename with the next argument
      outFile =  argv[++i];
    }
    else if(arg == "-c2p" && (i+1) < argc){
      //replace the filename with the next argument
      c2pFile =  argv[++i];
      cell2parent = true;
    }
    else if(arg == "-r" && (i+1) < argc){
      //replace the input refinement plan filename with the next argument
      planFile =  argv[++i];
      
    }else if(arg == "-pr" && (i+1) < argc){
      //replace the input refinement plan filename with the next argument
      parentPlanFile =  argv[++i];
      restart = true;
      
    }
    else{
      // Anything that we don't recognize gets recycled back into argv
      argv[j++] = argv[i];
    }
  }
  // Set the new size of the argument list to however many arguments
  // we put back in argv.
  argc = j;
  meshfile = pathname + meshfile;
  outFile = pathname + outFile;
  planFile = pathname + planFile;
  if(cell2parent)c2pFile = pathname + c2pFile;
  if(restart)parentPlanFile = pathname + parentPlanFile;


  // Setup the fact database.
  fact_db facts;
  

  
  
  param<std::string> planfile_par ;
  *planfile_par = planFile;
  facts.create_fact("balanced_planfile_par",planfile_par) ;

  param<std::string> meshfile_par ;
  *meshfile_par = meshfile;
  facts.create_fact("meshfile_par",meshfile_par) ;

  param<std::string> outfile_par ;
  *outfile_par = outFile;
  facts.create_fact("outfile_par",outfile_par) ;
  
  if(cell2parent){
    param<std::string> c2pfile_par ;
    *c2pfile_par = c2pFile;
    facts.create_fact("cell2parent_file_par",c2pfile_par) ;
    if(restart){
      param<std::string> parent_planfile_par ;
      *parent_planfile_par = parentPlanFile;
      facts.create_fact("parent_planfile_par",parent_planfile_par) ;
    }
  }
  
  // Setup the rule database.
  // Add all registered rules.
  rule_db rules;
  rules.add_rules(global_rule_list);
  
 
  if(Loci::MPI_rank == 0) std::cout <<"reading in meshfile" << std::endl;
  // Read in the mesh file.  Setup Loci datastructures
  if(!Loci::setupFVMGrid(facts,meshfile)) {
    std::cerr << "unable to read grid file '" << meshfile << "'" << std::endl ;
    Loci::Abort() ;
  }
  
  Loci::createLowerUpper(facts) ;
  Loci::createEdgesPar(facts) ;
 
  Loci:: parallelClassifyCell(facts);
  
 Loci::load_module("fvmadapt", rules);
 // if(Loci::MPI_rank==0){
//     Loci::ruleSet all_rules = rules.all_rules();
//     for(Loci::ruleSet::const_iterator ri = all_rules.begin();
//         ri != all_rules.end(); ri++){
//       cout << *ri << endl;
//     }
//     cout<< endl;
//   }
  if(cell2parent){
    if(!Loci::makeQuery(rules, facts, "cell2parent_output")) {
      std::cerr << "query failed!" << std::endl;
      Loci::Abort();
    }
    
  }
 
    
  if(!Loci::makeQuery(rules, facts, "node_output")) {
    std::cerr << "query failed!" << std::endl;
    Loci::Abort();
  }
   
  if(!Loci::makeQuery(rules, facts, "face_output")) {
    std::cerr << "query failed!" << std::endl;
    Loci::Abort();
  }
  if(!Loci::makeQuery(rules, facts, "volTag_output")) {
    std::cerr << "query failed!" << std::endl;
    Loci::Abort();
  }

 
  
  
  // Tell Loci to cleanup after itself and exit gracefully.
  Loci::Finalize();
  
}
