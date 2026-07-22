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
#include <queue>
#include <utility>

#include "plan_operations.h"

namespace {

struct BoundaryEdgeRanges {
  std::pair<int64, int64> edge0Range ;
  std::pair<int64, int64> edge3Range ;
} ;

} // namespace

std::vector<Loci::SetLong> project_face_plan_to_edge_splits(
  const std::vector<char>& facePlan,
  bool isQuadFace,
  const std::vector<bool>& edgeIsReversed) {
  const int64 maxCoordinate = int64(1) << MAXLEVEL ;
  std::vector<Loci::SetLong> edgeSplits(edgeIsReversed.size()) ;

  if(facePlan.empty()) { return edgeSplits ; }

  if(isQuadFace) {
    for(std::size_t edgeIndex=0;
        edgeIndex<edgeIsReversed.size(); edgeIndex++) {
      std::vector<char> edgePlan ;
      extract_quad_edge(facePlan, edgePlan, edgeIndex) ;
      if(edgePlan.size() == 0) { continue ; }

      std::queue<std::pair<int64, int64> > pendingRanges ;
      char splitCode = edgePlan.front() ;
      unsigned int planIndex = 0 ;
      std::pair<int64, int64> edgeRange =
        std::make_pair(int64(0), maxCoordinate) ;
      pendingRanges.push(edgeRange) ;

      int64 midpoint ;

      while(!pendingRanges.empty()) {
        // read in a code
        if(planIndex >= edgePlan.size()) {
          splitCode = 0 ;
        }else {
          splitCode = edgePlan[planIndex++] ;
        }

        // read in the range from the queue
        edgeRange = pendingRanges.front() ;

        // process the edge, push child ranges into the queue,
        // and record the midpoint
        switch(splitCode) {
          case 1: {
            midpoint = (edgeRange.first + edgeRange.second)/2 ;
            pendingRanges.push(
              std::make_pair(edgeRange.first, midpoint)) ;
            pendingRanges.push(
              std::make_pair(midpoint, edgeRange.second)) ;

            if(edgeIsReversed[edgeIndex]) {
              edgeSplits[edgeIndex].aset.insert(maxCoordinate - midpoint) ;
            }else {
              edgeSplits[edgeIndex].aset.insert(midpoint) ;
            }
            break ;
          }

          case 0:
            break ;

          case 8:
            pendingRanges.push(edgeRange) ;
            break ;

          default:
            std::cerr << "WARNING: illegal split code in "
                         "project_face_plan_to_edge_splits"
                      << std::endl ;
            break ;
        }
        pendingRanges.pop() ;
      }
    }
    return edgeSplits ;
  }

  const int faceEdgeCount = edgeIsReversed.size() ;

  // The root split contributes the midpoint of every face-local edge.
  std::queue<std::pair<int, BoundaryEdgeRanges> > pendingBoundaryStates ;
  std::vector<BoundaryEdgeRanges> rootChildRanges(faceEdgeCount) ;
  for(int edgeIndex=0; edgeIndex<faceEdgeCount; edgeIndex++) {
    if(edgeIsReversed[edgeIndex]) {
      rootChildRanges[edgeIndex].edge0Range =
        std::make_pair(maxCoordinate, maxCoordinate/2) ;
    }else {
      rootChildRanges[edgeIndex].edge0Range =
        std::make_pair(int64(0),maxCoordinate/2) ;
    }

    const int previousEdgeIndex =
      edgeIndex==0 ? faceEdgeCount-1 : edgeIndex-1 ;
    if(edgeIsReversed[previousEdgeIndex]) {
      rootChildRanges[edgeIndex].edge3Range =
        std::make_pair(maxCoordinate/2, int64(0)) ;
    }else {
      rootChildRanges[edgeIndex].edge3Range =
        std::make_pair(maxCoordinate/2, maxCoordinate) ;
    }

    pendingBoundaryStates.push(
      std::make_pair(edgeIndex, rootChildRanges[edgeIndex])) ;
    edgeSplits[edgeIndex].aset.insert(maxCoordinate/2) ;
  }

  BoundaryEdgeRanges currentRanges ;
  unsigned int planIndex = 1 ;
  char splitCode ;
  int boundaryState ;

  while(!pendingBoundaryStates.empty()) {
    currentRanges = pendingBoundaryStates.front().second ;
    boundaryState = pendingBoundaryStates.front().first ;
    if(planIndex >= facePlan.size()) {
      splitCode = 0 ;
    }else {
      // Read the next face split code.
      splitCode = facePlan[planIndex] ;
      planIndex++ ;
    }

    if(splitCode == 1) {
      // Define the child boundary states and record boundary midpoints.
      int childBoundaryStates[4] = {-1,-1,-1,-1} ;

      // Define the four children.
      BoundaryEdgeRanges childRanges[4] ;

      // If boundaryState == -1, every child is interior.
      if(boundaryState == -1) {
        for(int childIndex=0; childIndex<3; childIndex++) {
          childBoundaryStates[childIndex] = -1 ;
        }
      }else if(boundaryState >= 0 && boundaryState < faceEdgeCount) {
        // In this band, local edge 0 is on the face edge indexed by
        // boundaryState, and local edge 3 is on the preceding face edge.
        childBoundaryStates[2] = -1 ;
        childBoundaryStates[0] = boundaryState ;
        childRanges[0].edge0Range = std::make_pair(
          currentRanges.edge0Range.first,
          (currentRanges.edge0Range.first +
           currentRanges.edge0Range.second)/2) ;
        childRanges[0].edge3Range = std::make_pair(
          (currentRanges.edge3Range.first +
           currentRanges.edge3Range.second)/2,
          currentRanges.edge3Range.second) ;

        childBoundaryStates[1] = 2*faceEdgeCount+boundaryState ;
        childRanges[1].edge3Range = std::make_pair(
          (currentRanges.edge0Range.first +
           currentRanges.edge0Range.second)/2,
          currentRanges.edge0Range.second) ;

        const int previousEdgeIndex = boundaryState==0 ?
          faceEdgeCount-1 : boundaryState-1 ;
        childBoundaryStates[3] = faceEdgeCount + previousEdgeIndex ;
        childRanges[3].edge0Range = std::make_pair(
          currentRanges.edge3Range.first,
          (currentRanges.edge3Range.first +
           currentRanges.edge3Range.second)/2) ;

        edgeSplits[boundaryState].aset.insert(
          (currentRanges.edge0Range.first +
           currentRanges.edge0Range.second)/2) ;
        edgeSplits[previousEdgeIndex].aset.insert(
          (currentRanges.edge3Range.first +
           currentRanges.edge3Range.second)/2) ;
      }else if(boundaryState >= faceEdgeCount &&
               boundaryState < 2*faceEdgeCount) {
        // In this band, local edge 0 lies on one face boundary edge.
        childBoundaryStates[2] = childBoundaryStates[3] = -1 ;
        childBoundaryStates[0] = boundaryState ;
        childRanges[0].edge0Range = std::make_pair(
          currentRanges.edge0Range.first,
          (currentRanges.edge0Range.first +
           currentRanges.edge0Range.second)/2) ;

        childBoundaryStates[1] = faceEdgeCount + boundaryState ;
        childRanges[1].edge3Range = std::make_pair(
          (currentRanges.edge0Range.first +
           currentRanges.edge0Range.second)/2,
          currentRanges.edge0Range.second) ;

        edgeSplits[boundaryState-faceEdgeCount].aset.insert(
          (currentRanges.edge0Range.first +
           currentRanges.edge0Range.second)/2) ;
      }else if(boundaryState >= 2*faceEdgeCount &&
               boundaryState < 3*faceEdgeCount) {
        // In this band, local edge 3 lies on one face boundary edge.
        childBoundaryStates[1] = childBoundaryStates[2] = -1 ;
        childBoundaryStates[0] = boundaryState ;
        childRanges[0].edge3Range = std::make_pair(
          (currentRanges.edge3Range.first +
           currentRanges.edge3Range.second)/2,
          currentRanges.edge3Range.second) ;

        childBoundaryStates[3] = boundaryState - faceEdgeCount ;
        childRanges[3].edge0Range = std::make_pair(
          currentRanges.edge3Range.first,
          (currentRanges.edge3Range.first +
           currentRanges.edge3Range.second)/2) ;

        edgeSplits[boundaryState-2*faceEdgeCount].aset.insert(
          (currentRanges.edge3Range.first +
           currentRanges.edge3Range.second)/2) ;
      }else {
        std::cerr << "WARNING: illegal boundary state in "
                     "project_face_plan_to_edge_splits" << std::endl ;
        Loci::Abort() ;
      }

      for(int childIndex=0; childIndex<4; childIndex++) {
        pendingBoundaryStates.push(
          std::make_pair(childBoundaryStates[childIndex],
                         childRanges[childIndex])) ;
      }
    }
    pendingBoundaryStates.pop() ;
  }

  return edgeSplits ;
}
