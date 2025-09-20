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
#ifndef INITIALCONDITIONS_HEADER
#define INITIALCONDITIONS_HEADER

#include <Loci.h>
#include "flowTypes.h"
#include <vector>

namespace flowPsi {
  class geomTest : public Loci::CPTR_type {
  public:
    virtual bool inGeomPt(vect3d pt) const = 0 ;
  } ;
  Loci::CPTR<geomTest> geomTestFactory(string name, const options_list ol) ;

  struct ICstate_info {
    string name ;
    Loci::options_list state_info ;
  } ;
  struct ICgeom_info {
    Loci::CPTR<geomTest> geomTestFunc;
    int id ; // State for this geometry
  } ;
  
}



#endif
