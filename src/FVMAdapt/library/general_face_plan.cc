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

using std::make_pair;
using std::pair;

std::vector<char> merge_faceplan(std::vector<char>& facePlanL,
                                 std::vector<char>& facePlanR,
                                 int faceNodeCount) {
  if(facePlanL.size() == 0) { return facePlanR ; }
  if(facePlanR.size() == 0) { return facePlanL ; }
  std::vector<char> mergedFacePlan ;
  size_t leftPlanIndex = 0 ;
  size_t rightPlanIndex = 0 ;

  std::queue<pair<char, char> > pendingCodePairs ;
  char leftSplitCode ;
  char rightSplitCode ;

  // assume the first code of both facePlanL and facePlanR is 1
  leftPlanIndex++ ;
  rightPlanIndex++ ;

  mergedFacePlan.push_back(1) ;
  for(int i=0; i<faceNodeCount; i++) {
    if(leftPlanIndex < facePlanL.size()) {
      leftSplitCode = facePlanL[leftPlanIndex] ;
    }else {
      leftSplitCode = 0 ;
    }

    if(rightPlanIndex < facePlanR.size()) {
      rightSplitCode = facePlanR[rightPlanIndex] ;
    }else {
      rightSplitCode = 0 ;
    }
    pendingCodePairs.push(make_pair(leftSplitCode, rightSplitCode)) ;
    leftPlanIndex++ ;
    rightPlanIndex++ ;
  }

  while(!pendingCodePairs.empty()) {
    if(pendingCodePairs.front().first == 1 &&
       pendingCodePairs.front().second == 1) {
      mergedFacePlan.push_back(1) ;

      for(int i=0; i<4; i++) {
        if(leftPlanIndex < facePlanL.size()) {
          leftSplitCode = facePlanL[leftPlanIndex] ;
        }else {
          leftSplitCode = 0 ;
        }

        if(rightPlanIndex < facePlanR.size()) {
          rightSplitCode = facePlanR[rightPlanIndex] ;
        }else {
          rightSplitCode = 0 ;
        }
        pendingCodePairs.push(make_pair(leftSplitCode, rightSplitCode)) ;
        leftPlanIndex++ ;
        rightPlanIndex++ ;
      }
    }else if(pendingCodePairs.front().first == 1 &&
             pendingCodePairs.front().second == 0) {
      mergedFacePlan.push_back(1) ;
      for(int i=0; i<4; i++) {
        if(leftPlanIndex < facePlanL.size()) {
          leftSplitCode = facePlanL[leftPlanIndex] ;
        }else {
          leftSplitCode = 0 ;
        }
        pendingCodePairs.push(make_pair(leftSplitCode, 0)) ;
        leftPlanIndex++ ;
      }
    }else if(pendingCodePairs.front().first == 0 &&
             pendingCodePairs.front().second == 1) {
      mergedFacePlan.push_back(1) ;
      for(int i=0; i<4; i++) {
        if(rightPlanIndex < facePlanR.size()) {
          rightSplitCode = facePlanR[rightPlanIndex] ;
        }else {
          rightSplitCode = 0 ;
        }
        pendingCodePairs.push(make_pair(0, rightSplitCode)) ;
        rightPlanIndex++ ;
      }
    }else {
      mergedFacePlan.push_back(0) ;
    }
    pendingCodePairs.pop() ;
  }

  while(mergedFacePlan.size() != 0 && mergedFacePlan.back() == 0) {
    mergedFacePlan.pop_back() ;
  }
  reduce_vector(mergedFacePlan) ;
  return mergedFacePlan ;
}
