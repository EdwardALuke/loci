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

#ifndef READ_PAR_H
#define READ_PAR_H
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include "node_edge.h"

using std::ofstream;
using std::vector;
struct source_par{
  vect3d p1;
  vect3d p2;
  double r0;
  double s0;
  double r1;
  double s1;
  double a;
};

void readPar(string filename, vector<source_par>& source_pars); 
double get_spacing(const vect3d& p, const source_par& s);
double get_min_spacing(const vector<Node*>& nodes, const vector<source_par>& ss);

int tag_cell(const vector<Node*>& nodes, const vector<source_par>& source_pars, double min_edge_len);

#endif

