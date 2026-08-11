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

/**
 * @file read_par.h
 * @ingroup fvmadapt_input
 * @brief Source-parameter helpers used by FVMAdapt refinement tagging.
 */

/**
 * Source-based spacing request used by tag_cell().
 *
 * @ingroup fvmadapt_input
 *
 * The source is represented by a segment from `p1` to `p2`. The remaining
 * fields define requested spacing near and away from that segment; see
 * get_spacing() for the exact evaluation used by the current implementation.
 */
struct source_par {
  vect3d p1 ;
  vect3d p2 ;
  double r0 ;
  double s0 ;
  double r1 ;
  double s1 ;
  double a ;
};

/**
 * Reads source_par entries from a text file.
 */
void readPar(string filename, vector<source_par>& source_pars) ;

/**
 * Evaluates one source_par spacing request at a point.
 */
double get_spacing(const vect3d& p, const source_par& s) ;

/**
 * Computes the minimum requested spacing over a cell's vertices and center.
 */
double get_min_spacing(const vector<Node*>& nodes, const vector<source_par>& ss) ;

/**
 * Returns `1` when a cell should be refined for the supplied source requests.
 */
int tag_cell(const vector<Node*>& nodes, const vector<source_par>& source_pars, double min_edge_len) ;

#endif
