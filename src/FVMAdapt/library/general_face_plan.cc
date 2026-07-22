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

std::vector<char> merge_faceplan(std::vector<char>& planl, std::vector<char>& planr, int numNodes) {
  if(planl.size() == 0) { return planr ; }
  if(planr.size() == 0) { return planl ; }
  std::vector<char> fplan ;
  size_t ptl = 0 ;
  size_t ptr = 0 ;

  std::queue<pair<char, char> > Q ;
  char codel ;
  char coder ;

  // assume the first code of both planl and planr is 1
  ptl++ ;
  ptr++ ;

  fplan.push_back(1) ;
  for(int i=0; i<numNodes; i++) {
    if(ptl < planl.size()) {
      codel = planl[ptl] ;
    }else {
      codel = 0 ;
    }

    if(ptr < planr.size()) {
      coder = planr[ptr] ;
    }else {
      coder = 0 ;
    }
    Q.push(make_pair(codel, coder)) ;
    ptl++ ;
    ptr++ ;
  }

  while(!Q.empty()) {
    if(Q.front().first ==1 && Q.front().second == 1) {
      fplan.push_back(1) ;

      for(int i=0; i<4; i++) {
        if(ptl < planl.size()) {
          codel = planl[ptl] ;
        }else {
          codel = 0 ;
        }

        if(ptr < planr.size()) {
          coder = planr[ptr] ;
        }else {
          coder = 0 ;
        }
        Q.push(make_pair(codel, coder)) ;
        ptl++ ;
        ptr++ ;
      }
    }else if(Q.front().first == 1 && Q.front().second == 0) {
      fplan.push_back(1) ;
      for(int i=0; i<4; i++) {
        if(ptl < planl.size()) {
          codel = planl[ptl] ;
        }else {
          codel = 0 ;
        }
        Q.push(make_pair(codel, 0)) ;
        ptl++ ;
      }
    }else if(Q.front().first == 0 && Q.front().second == 1) {
      fplan.push_back(1) ;
      for(int i=0; i<4; i++) {
        if(ptr < planr.size()) {
          coder = planr[ptr] ;
        }else {
          coder = 0 ;
        }
        Q.push(make_pair(0, coder)) ;
        ptr++ ;
      }
    }else {
      fplan.push_back(0) ;
    }
    Q.pop() ;
  }

  while(fplan.size() != 0 && fplan.back() == 0) {
    fplan.pop_back() ;
  }
  reduce_vector(fplan) ;
  return fplan ;
}
