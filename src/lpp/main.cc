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
#include "lpp.h"
#include <Config/version_details.h>
#include <unistd.h>

using namespace std ;

list<string> include_dirs ;

bool prettyOutput = false ;
namespace {

  std::string version() {
    string rn = string(LOCI_BRANCH) + "-" + string(LOCI_VERSION) ;
    rn += " lpp compiled at " ;
    rn += __DATE__ ;
    rn += " " ;
    rn += __TIME__ ;
    return rn ;
  }
}


/// Storage for system variables
std::map<std::string,std::string> lpp_vars ;

// --------------------------------------------------------------------------
/// @brief setSystemVar() function to set system variable value in database
/// @param [var] - input variable to set
/// @param [val] - input value to assign to variable, replaces previous
/// definition of variable, or creates new variable if one doesn't exist
void setSystemVar(const std::string &var, const std::string &val) {
  lpp_vars[var] = val ;
}

// --------------------------------------------------------------------------
/// @brief checkSystemVar() returns true if provided variable exists
/// @param [var] - input variable
/// @returns boolean indicating existence of variable
bool checkSystemVar(const std::string &var) {
  return lpp_vars.find(var) != lpp_vars.end() ;
}

// --------------------------------------------------------------------------
/// @brief getSystemVar() returns value assigned to provided system variable
/// @param [var] - input variable
/// @returns string contained in the database for [var]
std::string getSystemVar(const std::string &var) {
  auto fp = lpp_vars.find(var) ;
  if(fp != lpp_vars.end())
    return fp->second ;
  return "" ;
}


// --------------------------------------------------------------------------
/// @brief evalSystemVar() evaluates system variables as integers
/// @param [varmap] - output map giving integer values to every variable
void evalSystemVars(std::map<string,int> &varmap) {
  for(auto ip=lpp_vars.begin();ip!=lpp_vars.end();++ip) {
    string name = ip->first ;
    int val = 0 ;
    if(ip->second.size() > 0 && ip->second[0]>='0' && ip->second[0]<='9') {
      val = atoi(ip->second.c_str()) ;
    }
    varmap[name] = val ;
  }
}
    

void Usage(int argc, char *argv[]) {
  cerr << "Loci Pre-Processor Usage:" << endl ;
  cerr << argv[0] <<" -I<dir> filename.loci -o filename.cc" << endl ;
  exit(-1) ;
}

bool no_cuda = false ;

// --------------------------------------------------------------------------
/// @brief intToString() is a utility routine for converting an integer to
/// its string representation
/// @param [val] input integer
/// @returns string representation of [val]

std::string intToString(int val) {
  std::ostringstream oss ;
  oss << val ;
  return std::string(oss.str()) ;
}

int main(int argc, char *argv[]) {

  bool file_given = false;
  string filename ;
  bool out_given = false ;
  string outfile ;
  bool test_parse = false ;
  int diag_level = 0 ;
  for(int i=1;i<argc;++i) {
    if(argv[i][0] == '-') {
      if(argv[i][1] == 'I') {
        string dir = &argv[i][2] ;
        include_dirs.push_back(dir) ;
      } else if(argv[i][1] == 'p') {
        prettyOutput = true ;
      } else  if(argv[i][1] == 'x') {
	no_cuda = true ;
      } else if(argv[i][1] == 't') {
        test_parse = true ;
      } else if(argv[i][1] == 'd') {
        if(argv[i][2] >= '0' && argv[i][2] <= '9') 
          diag_level = argv[i][2] - '0' ;
        if(argv[i][2] == '\0')
          diag_level = 10 ;
      } else if(argv[i][1] == 'o') {
        if(i+1>argc || out_given)
          Usage(argc,argv) ;
        outfile = argv[i+1] ;
        i++ ;
        out_given = true ;
      } else if(argv[i][1] == 'v') {
        cout << "Loci version: " << version() << endl ;
      } else if(argv[i][1] == 'V') {
        cout << "Loci version: " << version() << endl ;
      } else if(argv[i][1] == 'D') {
        string var ;
        string val ;
        for(int j=2;argv[i][j] != '\0';++j) {
          if(argv[i][j] == '=') {
            val = string(&argv[i][j+1]) ;
          }
          var += argv[i][j] ;
        }
        setSystemVar(var,val) ;
      } else {
	cerr << "Warning: Unknown option " << argv[i] << endl ;
	//	exit(-1) ;
      }
    } else {
      if(file_given == true) {
        cerr << "multiple filenames given" << endl ;
        Usage(argc,argv) ;
      }
      filename = argv[i] ;
      file_given = true ;
    }
  }
  // Add system version variables to system vars
  setSystemVar("LOCI_VERSION_MAJOR",intToString(LOCI_VERSION_MAJOR)) ;
  setSystemVar("LOCI_VERSION_MINOR",intToString(LOCI_VERSION_MINOR)) ;

  if(!file_given) {
    cerr << "no filename" << endl ;
    Usage(argc,argv) ;
  }
  parseSharedInfo parseInfo ;
#ifndef USE_CUDA_RT
  no_cuda = true ;
#endif
  parseInfo.no_cuda = no_cuda ;
  parseInfo.test_parse = test_parse ;
  parseInfo.diag_level = diag_level ;
  
  parseFile parser ;
  try {
    if(out_given) {
      ofstream file(outfile.c_str(),ios::out) ;
      if(file.fail()) {
        cerr << "unable to open file " << outfile << " for writing!" << endl ;
        exit(-1) ;
      }
      parser.processFile(filename,file,parseInfo) ;
    } else {
      parser.processFile(filename,cout,parseInfo) ;
    }
  } catch(parseError pe) {
    if(pe.error_type != "syntax error")
      cerr << pe.error_type << endl ;
    if(out_given)
      unlink(outfile.c_str()) ;
    exit(-1) ;
  } catch(...) {
    cerr << "Unknown exception caught!" << endl ;
    if(out_given)
      unlink(outfile.c_str()) ;
    exit(-1) ;
  }

  return 0 ;
}
