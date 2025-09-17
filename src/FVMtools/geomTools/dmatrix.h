//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
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
#ifndef DMATRIX_H
#define DMATRIX_H
#include <vector>

class dmatrix {
  int rstart, cstart ;
  int nr, nc ;
  std::vector<double> data ;
 public:
  dmatrix(int nrl, int nrh, int ncl, int nch) {
    nr = nrh-nrl+1 ;
    nc = nch-ncl+1 ;
    rstart = nrl ;
    cstart = ncl ;
    std::vector<double> tmp(nr*nc) ;
    data.swap(tmp) ;
  }
  double *operator[](int r) { double *col =  &data[(r-rstart)*nc] ; return col-rstart ; }
} ;

class dvector {
  int start ;
  std::vector<double> data ;
 public:
  dvector(int nl, int nh) {
    std::vector<double> tmp(nh-nl+1) ;
    start = nl ;
    data.swap(tmp) ;
  }
  double &operator[](int i) { return data[i-start]; }
} ;

void svdcmp(dmatrix &a, int m, int n, dvector &w, dmatrix &v);

#endif
