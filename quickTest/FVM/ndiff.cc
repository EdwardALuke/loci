//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#include <Config/conf.h>
#include <iostream>
#include <fstream>
#include <vector>
#ifdef NO_CMATH
#include <math.h>
#else
#include <cmath>
#endif
#include <stdlib.h>
using namespace std ;

int main(int ac, char *av[]) {
  double min_error = 0.0 ;
  if(ac <= 3 || ac >= 6) {
    cerr << "3 arguments required!" << endl ;
    for(int i=0;i<ac;++i)
      cerr << "'" << av[i] << "' " ;
    cerr << endl ;
    exit(EXIT_FAILURE) ;
  }
  char *file1 = av[1] ;
  char *file2 = av[2] ;
  double tol = atof(av[3]) ;
  ifstream in1(file1,ios::in) ;
  ifstream in2(file2,ios::in) ;  

  if(ac == 5)
    min_error = atof(av[4]) ;

  if( in1.fail() ){
    cerr << "unable to open " << file1 << endl ;
    exit(EXIT_FAILURE) ;
  }
  if( in2.fail() ) {
    cerr << "unable to open " << file2 << endl ;
    exit(EXIT_FAILURE) ;
  }

  vector<double> v1 ;
  vector<double> v2 ;
  while(!in1.eof() || !in2.eof()) {
    double d1=1e33,d2=1e33 ;
    in1 >> d1 ;
    in2 >> d2 ;
    v1.push_back(d1) ;
    v2.push_back(d2) ;
    
    while(in1.peek() == ' ' || in1.peek() == '\t')
      in1.get() ;
    while(in2.peek() == ' ' || in2.peek() == '\t')
      in2.get() ;
    if(in1.peek() == '\r' || in1.peek() == '\n')
      break ;
    if(v1.size() > 100)
      break ;
  }

  const int num_cols = v1.size() ;
  vector<double> RMSdelta(num_cols) ;
  vector<double> RMS(num_cols) ;

  for(int i=0;i<num_cols;++i) {
    RMSdelta[i] = (v1[i]-v2[i])*(v1[i]-v2[i]) ;
    RMS[i] = max(v1[i]*v1[i],v2[i]*v2[i]) ;
  }
  int count = 1 ;
  do {
    for(int i=0;i<num_cols;++i) {
      double d1=1e33,d2=1e33 ;
      in1 >> d1 ;
      in2 >> d2 ;
      if((in1.eof() && !in2.eof()) || (!in1.eof() && in2.eof())) {
        cerr << "files mismatched!" << endl ;
        exit(EXIT_FAILURE) ;
      }
      if(!in1.eof() || !in2.eof()) {
        if(0==i) 
          count++ ;
        RMSdelta[i] += (d1-d2)*(d1-d2) ;
        RMS[i] += max(d1*d1,d2*d2) ;
      }
    }
  }  while(!in1.eof() || !in2.eof()) ;

  bool fail = false  ;
  for(int i=0;i<num_cols;++i) {
    double delta_abs = sqrt(RMSdelta[i]) ; // rms of difference
    double rms_abs = sqrt(RMS[i]) ;        // rms of signal
    double delta_rel = delta_abs/(rms_abs+min_error) ; // scaled error
    if(delta_abs > min_error && delta_rel > tol) {
      cerr << "compare failed on column " << i
           << " with a relative rms error of " << delta_rel << endl;
      fail = true ;
    }
    // Output the relative error and the rms average value for each column
    cout << i << "-" << delta_rel << ' ' << rms_abs/sqrt(double(count)) << endl ;
  }
  if(fail) 
    exit(EXIT_FAILURE) ;
    
  exit(EXIT_SUCCESS) ;
}
