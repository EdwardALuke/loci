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
#include "lpp.h"
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

void Usage(int argc, char *argv[]) {
  cerr << "Loci Pre-Processor Usage:" << endl ;
  cerr << argv[0] <<" -I<dir> filename.loci -o filename.cc" << endl ;
  exit(-1) ;
}

bool no_cuda = false ;

int main(int argc, char *argv[]) {


  bool file_given = false;
  string filename ;
  bool out_given = false ;
  string outfile ;
  for(int i=1;i<argc;++i) {
    if(argv[i][0] == '-') {
      if(argv[i][1] == 'I') {
        string dir = &argv[i][2] ;
        include_dirs.push_back(dir) ;
      } else if(argv[i][1] == 'p') {
        prettyOutput = true ;
      } else  if(argv[i][1] == 'x') {
	no_cuda = true ;
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
	// ignore this option
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
  if(!file_given) {
    cerr << "no filename" << endl ;
    Usage(argc,argv) ;
  }
  parseSharedInfo parseInfo ;
#ifndef USE_CUDA_RT
  no_cuda = true ;
  parseInfo.no_cuda = true ;
#endif

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
