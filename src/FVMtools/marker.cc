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
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                        test8.cc
//
//   this file read in a coarse grid file, a refinement plan file(based on the coarse grid),
//and a tag file(if input refinement file is provided, the tag file is for new grid generated using the input
//refinement plan, otherwise, the tag file is for the coarse grid file),
//and then write out a new  refinement plan file(the plan also based on the coarse grid).
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////



#include <iostream>
#include <string>

#include <Loci.h>
#include <stdlib.h>
#include <limits>
#include "./FVMAdapt/globals.h"
#include <Loci>

using std::cout ;


using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;
using Loci::MPI_processes;
namespace Loci{
  void parallelClassifyCell(fact_db &facts) ;
  void createEdgesPar(fact_db &facts) ;
}

  
int main(int argc, char ** argv) {
  // Let Loci initialize itself.
  // This also gives Loci first dibs on the command line arguments.
  Loci::Init(&argc, &argv);
  
  
  // This is the name of the mesh file that you want to read in.
  // This may be overridden by the command line argument "-g file.xdr"
  //  string meshfile = "testGrid.xdr";
  string meshfile ;
  //This is the name of the input refinement plan file
  string planFile = "out1.plan"; //default file can be empty file
  //This is the name of the output refinement file
  string outFile  = "out.plan";
  //this is the name of the tag file
  string tagFile ;
  //this is the name of the tolerance file, this option allows
  //each marked node has a unique tolerance value
  string parFile;

  //this is the name of xml file
  string xmlFile;
  // Here's where we parse out the command line arguments that are
  // relevant to this program.
  int j=1;
  string pathname = "";
  //if the input refinement filename is provided in command line, restart is true,
  //otherwise, restart is false
  bool restart = false;
  //if xml_input is true, cell plan is defines by a region and tol value. otherwise,
  //it's decided by tagfile
  bool xml_input = false;
  bool tol_input = false;
  bool tag_input = false;
  bool par_input = false;
  bool levels_input = false;
  bool ctag_input = false;
  int split_mode = 0;//default mode, hexcells and prisms split according to edge length
  //print out help info
  if( (argc == 1)||(argc==2) ){
  
    if(Loci::MPI_rank == 0){
      cout<<"Usuage:" << endl;
      cout <<"marker <options> <filename/input value> <options> <filename/input value> ... "<< endl;
      cout << endl;
      cout << "options:" << endl;
      cout <<"-g <file> -- original grid file(required)," << endl;
      cout <<"             refinement plans are based on this grid" << endl;
      cout <<"-r <file> -- input refinement plan file(optional)" <<endl;
      
      cout << "-xml <file> -- input xml file(optional), the file define a geometric region" << endl;
      cout << "               all cells inside the region will be refined to a tolerance value" << endl;
      cout << "               -xml option is uaually used with -tol option" << endl;
      
      cout <<"-tag <file> -- input node tag file(optional), " << endl;
      cout <<"               if there is an input refinement plan file,"<<endl;
      cout <<"               the tag file is for the refined grid" << endl;
      cout <<"               otherwise, the tag file is for the original grid" << endl;
      cout <<"               -tag option might be used with -levels option" << endl;

      cout <<"-ctag <file> -- input cell tag file(optional), " << endl;
      cout <<"               if there is an input refinement plan file,"<<endl;
      cout <<"               the tag file is for the refined grid" << endl;
      cout <<"               otherwise, the tag file is for the original grid" << endl;
      cout <<"               -tag option might be used with -levels option" << endl;
      
      cout <<"-par <file> -- input parameter file(optional), " << endl;
      cout <<"               The parameters will decide recursively if each cell need to split" <<endl;               
      cout << "              -xml option and -tag option and --par can not be selected at the same time" <<endl;
      cout<< "               and one of them must be selected"<<endl;
      
      cout <<"-o <file> -- output refinement plan file" << endl;
      cout <<"-tol <double> -- tolerance, minimum grid spacing allowed(default value: 1e-10), need to be specified for -xml option" << endl;
      cout <<"-fold <double> -- twist value, maximum face folding allowed(default value: 1.5708)" << endl;
      cout <<"-factor <double> -- anisotripic factor value, minimum allowed ratio between average edge length in different direction (default value: 2)" << endl;
      cout <<"-levels <int> --  levels of refinement(default value: 1), for anisotropic refinement, levels can only be 1" << endl;
      cout << "-mode <int: 0, 1, 2, 3> -- split mode 0: anisotropic split according to edge length, this is the default value " << endl;
      cout << "                        -- split mode 1: don't refine in z direction" << endl;
      cout << "                        -- split mode 2: fully isotropic split" << endl;
      cout <<"-balance<int:0, 1, 2> -- balance option 0: no edge's depth is greater than 1, this is the default value" << endl;
      cout <<"                      -- balance option 1: option 0 plus no cell has more than half of its faces split" << endl;
      cout <<"                      -- balance option 2:  option 0, option 1 plus no cell has two opposite faces split" << endl; 
      
    
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
    
    else if(arg == "-r" && (i+1) < argc){
      //replace the input filename with the next argument
      planFile =  argv[++i];
      restart = true;
    }
    
    else if(arg == "-tag" && (i+1) < argc){
      //replace the tag filename with the next argument
      tagFile =  argv[++i];
      tag_input = true;
    }
    else if(arg == "-ctag" && (i+1) < argc){
      //replace the tag filename with the next argument
      tagFile =  argv[++i];
      ctag_input = true;
    }
    else if(arg == "-par" && (i+1) < argc){
      //replace the tag filename with the next argument
      parFile =  argv[++i];
      par_input = true;
    }
    else if(arg == "-xml" && (i+1) < argc){
      //replace the xml filename with the next argument
      xmlFile =  argv[++i];
      xml_input = true;
    }

     
    else if(arg == "-tol" && (i+1) < argc){
      //replace the tolarence(minimum grid spacing) with the next argument
      char** endptr = 0;
      Globals::tolerance = strtod(argv[++i], endptr);
      tol_input = true;
    }
    else if(arg == "-balance" && (i+1) < argc){
      //replace the balance option with the next argument
      
      Globals::balance_option = atoi(argv[++i]);
      
    }
    else if(arg == "-fold" && (i+1) < argc){
      //replace the fold(maximum face fold value)with the next argument
      char** endptr = 0;
      Globals::fold = strtod(argv[++i], endptr);
      
    }
    else if(arg == "-factor" && (i+1) < argc){
      //replace the fold(maximum face fold value)with the next argument
      char** endptr = 0;
      Globals::factor = strtod(argv[++i], endptr);
      
    }

    else if(arg == "-mode" && (i+1) < argc){
      //replace the split mode with the next argument
      split_mode = atoi(argv[++i]);
      
    }
       
    else if(arg == "-levels" && (i+1) < argc){
      //replace the fold(maximum face fold value)with the next argument
      Globals::levels = atoi(argv[++i]);
      levels_input = true;
       
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
  tagFile = pathname + tagFile;
  xmlFile = pathname + xmlFile;
  parFile = pathname + parFile;
  
  if(Loci::MPI_rank == 0){
    cout <<"Marker running" <<endl;
    cout << "fold: " << Globals::fold << endl;
    cout << "tolerance: "<<Globals::tolerance <<endl;
    cout << "levels: " <<Globals::levels << endl;
    cout << "restart: " << restart << endl;
    cout << "split_mode: " << split_mode << endl;
    cout << "balance_option: " <<Globals::balance_option << endl;
        
  }

  if(split_mode == 0 && Globals::levels > 1){
    if(Loci::MPI_rank == 0) cerr << "WARNING: multi-level refinement is not allowed in anisotropic refinement"<< endl;
    Loci::Abort();
  }
  
  if(xml_input && !tol_input){
    if(Loci::MPI_rank == 0)  cerr <<"WARNING: The user need specify tolerance value at -xml option"<< endl;
    Loci::Abort();
  }
  
  if(xml_input &&  Globals::tolerance < 1e-37){
    if(Loci::MPI_rank == 0)  cerr <<"WARNING: small tolerance value can cause infinite loop at -xml option"<< endl;
    Loci::Abort();
  }
  
  if(xml_input &&  Globals::tolerance < 1e-20){
    if(Loci::MPI_rank == 0)   cerr <<"WARNING: small tolerance value can cause infinite loop at -xml option"<< endl;
  }

  if(xml_input &&  levels_input){
    if(Loci::MPI_rank == 0)  { cerr <<"WARNING: -xml options will split a cell until the tolerance value is met," <<endl;
      cerr<<"          the value of levels is not used"<< endl;
    }
  }
  if((!xml_input) && (!tag_input)&&(!par_input)&&(!ctag_input)){
    if(Loci::MPI_rank == 0)  cerr <<"WARNING: one option has to be specified, \neither  -tag option or -xml option, -ctag option or -par option"<< endl;
    Loci::Abort();
  }
  
  if(xml_input && tag_input){
    if(Loci::MPI_rank == 0)  cerr <<"WARNING: only one option need to be specified, either  -tag option or -xml option"<< endl;
    Loci::Abort();
  }
  
  if(xml_input && par_input){
    if(Loci::MPI_rank == 0)  cerr <<"WARNING: only one option need to be specified, either  -par option or -xml option"<< endl;
    Loci::Abort();
  }

  if(tag_input && par_input){
    if(Loci::MPI_rank == 0)  cerr <<"WARNING: only one option need to be specified, either  -par option or -tag option"<< endl;
    Loci::Abort();
  }
  if(ctag_input && par_input){
    if(Loci::MPI_rank == 0)  cerr <<"WARNING: only one option need to be specified, either  -par option or -ctag option"<< endl;
    Loci::Abort();
  }
  
  if(ctag_input && tag_input){
    if(Loci::MPI_rank == 0)  cerr <<"WARNING: only one option need to be specified, either  -tag option or -ctag option"<< endl;
    Loci::Abort();
  }
  if(ctag_input && xml_input){
    if(Loci::MPI_rank == 0)  cerr <<"WARNING: only one option need to be specified, either  -xml option or -ctag option"<< endl;
  }
  if(tol_input && par_input){
    if(Loci::MPI_rank == 0){
      cerr<<"WARNING: in -par option, the value of -tol option is not used"<< endl;
    }
  }

  if(levels_input && par_input){
    if(Loci::MPI_rank == 0){
      cerr <<"WARNING: in -par options , the value of -levels option is not used"<< endl;
    }
  } 
 

  // Setup the rule database.
  // Add all registered rules.  
  rule_db rules;
  rules.add_rules(global_rule_list);
  
  // Setup the fact database.
  fact_db facts;
  
  if(Loci::MPI_rank == 0) cout <<"reading in meshfile" << std::endl;

  // Read in the mesh file.  
  if(!Loci::setupFVMGrid(facts,meshfile)) {
    std::cerr << "unable to read grid file '" << meshfile << "'" << std::endl ;
    Loci::Abort() ;
  }
  
  //Setup Loci datastructures
  Loci::createLowerUpper(facts) ;
  Loci::createEdgesPar(facts) ;
  Loci:: parallelClassifyCell(facts);
  
  //this is a dummy parameter to trick Loci scheduler
  param<bool> beginWithMarker;
  *beginWithMarker = true;
  facts.create_fact("beginWithMarker",beginWithMarker) ; 
  
  param<std::string> outfile_par ;
  *outfile_par = outFile;
  facts.create_fact("plan_outfile_par",outfile_par) ; 
 
  if(tag_input){
    param<std::string> tagfile_par ;
    *tagfile_par = tagFile;
    facts.create_fact("tagfile_par",tagfile_par) ;
  }
 
  if(ctag_input){
    param<std::string> tagfile_par ;
    *tagfile_par = tagFile;
    facts.create_fact("cell_tagfile_par",tagfile_par) ;
  }
   
  if(par_input){
    param<std::string> parfile_par ;
    *parfile_par = parFile;
    facts.create_fact("parfile_par",parfile_par) ;
  }
  
  if(xml_input){
    param<std::string> xmlfile_par ;
    *xmlfile_par = xmlFile;
    facts.create_fact("xmlfile_par",xmlfile_par) ;
  }
  
  //parameters to identify different options
  if(restart){
    param<std::string> planfile_par ;
    *planfile_par = planFile;
    facts.create_fact("planfile_par",planfile_par) ;
    
    if(xml_input){
      param<int> restart_xml_par;
      *restart_xml_par = 1;
      facts.create_fact("restart_xml_par",restart_xml_par);
    } else if(tag_input||ctag_input){
      param<int> restart_tag_par;
      *restart_tag_par = 1;
      facts.create_fact("restart_tag_par",restart_tag_par);
    }else if(par_input){
      param<int> restart_par_par;
      *restart_par_par = 1;
      facts.create_fact("restart_par_par",restart_par_par); 
    }
  }else{
    if(xml_input){
      param<int> norestart_xml_par;
      *norestart_xml_par = 1;
      facts.create_fact("norestart_xml_par",norestart_xml_par);
    }else if(tag_input||ctag_input){
      param<int> norestart_tag_par;
      *norestart_tag_par = 1;
      facts.create_fact("norestart_tag_par",norestart_tag_par);
    }else if(par_input){
      param<int> norestart_par_par;
      *norestart_par_par = 1;
      facts.create_fact("norestart_par_par",norestart_par_par);
    }
  }
  
  param<int> split_mode_par;
  *split_mode_par = split_mode;
  facts.create_fact("split_mode_par", split_mode_par);

  Loci::load_module("fvmadapt", rules);
  //  if(Loci::MPI_rank==0){
  //     Loci::ruleSet all_rules = rules.all_rules();
  //     for(Loci::ruleSet::const_iterator ri = all_rules.begin();
  //         ri != all_rules.end(); ri++){
  //       cout << *ri << endl;
  //     }
  //     cout<< endl;
  //   }
  
  if(!Loci::makeQuery(rules, facts, "cellplan_output")) {
    std::cerr << "query failed!" << std::endl;
    Loci::Abort();
  }
  
  // Tell Loci to cleanup after itself and exit gracefully.
  Loci::Finalize();
}
  
