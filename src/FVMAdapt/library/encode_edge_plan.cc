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

#include <queue>
#include <utility>

#include "plan_operations.h"

void encode_edge_plan(const Loci::SetLong& splitCoordinates,
                      std::vector<char>& edgePlan) {
  edgePlan.clear() ;
  if(splitCoordinates.aset.empty()) { return ; }

  // Define the integer-coordinate range of the edge.
  int64 maxValue = int64(1) << MAXLEVEL ;

  // Build the refinement plan from the requested split coordinates.
  std::queue<std::pair<int64, int64> > ranges ;
  ranges.push(std::make_pair(int64(0), maxValue)) ;

  while(!ranges.empty()) {
    const std::pair<int64, int64> range = ranges.front() ;
    const int64 midpoint = (range.first + range.second) / 2 ;

    if(splitCoordinates.aset.find(midpoint) != splitCoordinates.aset.end()) {
      ranges.push(std::make_pair(range.first, midpoint)) ;
      ranges.push(std::make_pair(midpoint, range.second)) ;
      edgePlan.push_back(1) ;
    }else {
      edgePlan.push_back(0) ;
    }
    ranges.pop() ;
  }

  while(!edgePlan.empty() && edgePlan.back() == 0) {
    edgePlan.pop_back() ;
  }
  reduce_vector(edgePlan) ;
}
