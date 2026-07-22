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

#include <cstddef>
#include <iostream>
#include <iostream>
#include <queue>
#include <utility>

#include "plan_operations.h"

std::vector<Loci::SetLong> project_face_plan_to_edge_splits(
  const std::vector<char>& facePlan,
  bool isQuadFace,
  const std::vector<bool>& edgeReversed) {
  const int64 maxValue = int64(1) << MAXLEVEL ;
  std::vector<Loci::SetLong> edgeSplits(edgeReversed.size()) ;

  if(facePlan.empty()) { return edgeSplits ; }

  if(isQuadFace) {
    for(std::size_t i=0; i<edgeReversed.size(); i++) {
      std::vector<char> plan ;
      extract_quad_edge(facePlan, plan, i) ;
      if(plan.size() == 0) { continue ; }

      std::queue<std::pair<int64, int64> > Q ;
      char code = plan.front() ;
      unsigned int index = 0 ;
      std::pair<int64, int64> pRange =
        std::make_pair(int64(0), maxValue) ;
      Q.push(pRange) ;

      int64 midP ;

      while(!Q.empty()) {
        // read in a code
        if(index >= plan.size()) {
          code = 0 ;
        }else {
          code = plan[index++] ;
        }

        // read in the range from the queue
        pRange = Q.front() ;

        // process the cell, push children ranges into Q
        // and put the middle point into points
        switch(code) {
          case 1: {
            midP = (pRange.first + pRange.second)/2 ;
            Q.push(std::make_pair(pRange.first, midP)) ;
            Q.push(std::make_pair(midP, pRange.second)) ;

            if(edgeReversed[i]) {
              edgeSplits[i].aset.insert(maxValue - midP) ;
            }else {
              edgeSplits[i].aset.insert(midP) ;
            }
            break ;
          }

          case 0:
            break ;

          case 8:
            Q.push(pRange) ;
            break ;

          default:
            std::cerr << "WARNING: illegal splitCode in rule edge_points_apply"
                      << std::endl ;
            break ;
        }
        Q.pop() ;
      }
    }
    return edgeSplits ;
  }

  const int nfnode = edgeReversed.size() ;

  // first time split face, put all the edgecenter into pointSet
  std::queue<std::pair<int, TwoEdge> > Q ;
  std::vector<TwoEdge> child(nfnode) ;
  for(int i=0; i<nfnode; i++) {
    if(edgeReversed[i]) {
      child[i].e0 = std::make_pair(maxValue, maxValue/2) ;
    }else {
      child[i].e0 = std::make_pair(int64(0),maxValue/2) ;
    }

    if(edgeReversed[i==0?nfnode-1:i-1]) {
      child[i].e3 = std::make_pair(maxValue/2, int64(0)) ;
    }else {
      child[i].e3 = std::make_pair(maxValue/2, maxValue) ;
    }

    Q.push(std::make_pair(i, child[i])) ;
    edgeSplits[i].aset.insert(maxValue/2) ;
  }

  TwoEdge current ;
  unsigned int index = 1 ;
  char mySplitCode ;
  int edgeID ;

  while(!Q.empty()) {
    current = Q.front().second ;
    edgeID = Q.front().first ;
    if(index >= facePlan.size()) {
      mySplitCode = 0 ;
    }else {
      // take a code from splitcode
      mySplitCode = facePlan[index] ;
      index++ ;
    }

    if(mySplitCode == 1) {
      // define new edgeID and push the children into Q
      // put edgecenter into pointSet if necessary
      int newID[4] = {-1,-1,-1,-1} ;

      // define the 4 children
      TwoEdge child[4] ;

      // if edgeID== -1, each newID is -1, no edgecenter need to be put
      // into pointSet
      if(edgeID == -1) {
        for(int i=0; i<3; i++) {
          newID[i] = -1 ;
        }
      }else if(edgeID >= 0 && edgeID < nfnode) {
        // if edgeID is in [0, nfnode), edge 0 is on edgeID, edge3 is on edgeID-1
        newID[2] = -1 ;
        newID[0] = edgeID ; //e0 on edgeID, e3 on edgeID-1
        child[0].e0 = std::make_pair(
          current.e0.first, (current.e0.first + current.e0.second)/2) ;
        child[0].e3 = std::make_pair(
          (current.e3.first + current.e3.second)/2, current.e3.second) ;

        newID[1] = 2*nfnode+edgeID ; //e3 on edgId
        child[1].e3 = std::make_pair(
          (current.e0.first + current.e0.second)/2, current.e0.second) ;

        newID[3] = nfnode +(edgeID==0?nfnode-1:edgeID-1) ; //e0 on edgeId
        child[3].e0 = std::make_pair(
          current.e3.first, (current.e3.first + current.e3.second)/2) ;

        edgeSplits[edgeID].aset.insert(
          (current.e0.first + current.e0.second)/2) ;
        edgeSplits[edgeID==0?nfnode-1:edgeID-1].aset.insert(
          (current.e3.first + current.e3.second)/2) ;
      }else if(edgeID >=nfnode && edgeID <2*nfnode) {
        // edge 0 on edgeID-nfnode
        newID[2] = newID[3] = -1 ;
        newID[0] = edgeID ; //e0 on edgeId
        child[0].e0 = std::make_pair(
          current.e0.first, (current.e0.first + current.e0.second)/2) ;

        newID[1] = nfnode + edgeID ; //e3 on edgeId
        child[1].e3 = std::make_pair(
          (current.e0.first + current.e0.second)/2, current.e0.second) ;

        edgeSplits[edgeID-nfnode].aset.insert(
          (current.e0.first + current.e0.second)/2) ;
      }else if(edgeID >= 2*nfnode && edgeID < 3*nfnode) {
        // edge 3 on edgeID-2*nfnode
        newID[1] = newID[2] = -1 ;
        newID[0] = edgeID ; //e3 on edgeId
        child[0].e3 = std::make_pair(
          (current.e3.first + current.e3.second)/2, current.e3.second) ;

        newID[3] = edgeID - nfnode ; //e0 on edgeId
        child[3].e0 = std::make_pair(
          current.e3.first, (current.e3.first + current.e3.second)/2) ;

        edgeSplits[edgeID-2*nfnode].aset.insert(
          (current.e3.first + current.e3.second)/2) ;
      }else {
        std::cerr << " WARNING: illegal edgeID" << std::endl ;
        Loci::Abort() ;
      }

      for(int i=0; i<4; i++) {
        Q.push(std::make_pair(newID[i],child[i])) ;
      }
    }
    Q.pop() ;
  }

  return edgeSplits ;
}
