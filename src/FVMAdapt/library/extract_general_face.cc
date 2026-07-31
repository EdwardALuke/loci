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

#include <list>
#include <queue>
#include <utility>
#include <vector>

#include <Loci.h>

#include "diamondcell.h"
#include "plan_operations.h"

using std::cerr;
using std::endl;
using std::list;
using std::queue;

// This file extracts a facePlan from a cellPlan, and then merges two facePlan together
// the extracting algrithm uses the numbering system of DiamondCell and is faster than
// straight forward approach


std::vector<char> extract_general_face(const Entity* lower, int lower_size,
                                       const Entity* upper, int upper_size,
                                       const Entity* boundary_map, int boundary_map_size,
                                       const const_multiMap& face2node,
                                       const const_multiMap& face2edge,
                                       const const_MapVec<2>& edge2node,
                                       const std::vector<char>& cellPlan,
                                       Entity ff,
                                       const const_store<int>& node_remap
                                       ) {

  // if cellPlan empty, facePlan also empty
  std::vector<char> facePlan ;
  if(cellPlan.size() == 0) {
    reduce_vector(facePlan) ;
    return facePlan ;
  }

  // first build a Cell from cellPlan
  std::list<Node*> node_list ;
  std::list<Edge*> edge_list ;
  std::list<Face*> face_list ;
  Cell* aCell = build_general_cell(lower, lower_size,
                                   upper, upper_size,
                                   boundary_map, boundary_map_size,
                                   face2node,
                                   face2edge,
                                   edge2node,
                                   node_list,
                                   edge_list,
                                   face_list,
                                   node_remap);

  // split  a general cell isotropically according to cellPlan,
  // only child is defined
  aCell->empty_resplit(cellPlan) ;

  int findex = find_face_index(lower, lower_size, upper, upper_size,
                               boundary_map, boundary_map_size,
                               face2node, ff, node_remap) ;

 char faceOrient = aCell->faceOrient[findex] ;

  facePlan.push_back(1) ;
  std::queue<pair<DiamondCell*, int> > Q ;

  // each node on ff corresponds to a child cell, and each
  // cell corresponds to a tree in diamonds, find  which tree correspond to
  // the node , and put it in the Q

  // At the same time , find out faceID in the childCell
  Face* theFace = aCell->face[findex] ;

  for(int i=0; i<theFace->numEdge; i++) {
    Node* theNode ;
    if(theFace->needReverse[i]) {
      theNode = theFace->edge[i]->tail ;
    }else {
      theNode = theFace->edge[i]->head ;
    }

    int childID = 0 ;
    for(childID = 0; childID < aCell->numNode; childID++) {
      if(aCell->node[childID] == theNode) { break ; }
    }
    if(childID == aCell->numNode) {
      cerr << " WARNING: face nodes not in cell" << endl ;
      Loci::Abort() ;
    }

    // to find faceID, first set up n2f and n2e
    std::vector<Face*> n2f ;
    std::vector<Edge*> n2e ;
    std::vector<int> rot ;

    aCell->set_n2f_n2e(n2f, n2e, rot, childID) ;

    // find faceID
    int faceID ;
    int n2f_size = n2f.size() ;
    for(faceID=0; faceID<n2f_size; faceID++) {
      if(n2f[faceID] == aCell->face[findex]) { break ; }
    }
    if(faceID == n2f_size) {
      cerr << " WARNING: can not find faceID" << endl ;
      Loci::Abort() ;
    }

    // when split a general cell, center of face n2f[i] -> vertex[i+2],
    // vertex[1] and vertex[i+2] belongs to face nfold+i in the childCell

    faceID = faceID + aCell->child[childID]->getNfold() ;
    Q.push(make_pair(aCell->child[childID], faceID)) ;
  } //the general face is split

  DiamondCell* current ;
  int currentFace ;
  while(!Q.empty()) {
    current = Q.front().first ;
    currentFace = Q.front().second ;
    int currentNfold = current->getNfold() ;
    if(current->getChildCell() != 0) {
      facePlan.push_back(1) ;

      // push the selected four children in the Q
      if(faceOrient == 0) {
        Q.push(make_pair(current->getChildCell(1), currentFace)) ;
        Q.push(make_pair(current->getChildCell(currentFace == currentNfold? (2*currentNfold+1):(currentFace+1)),4)) ;
        Q.push(make_pair(current->getChildCell(currentFace-currentNfold+2), 4)) ;
        Q.push(make_pair(current->getChildCell(currentFace +2), 5)) ;
      }else {
        Q.push(make_pair(current->getChildCell(1), currentFace)) ;
        Q.push(make_pair(current->getChildCell(currentFace +2), 5)) ;
        Q.push(make_pair(current->getChildCell(currentFace-currentNfold+2), 4)) ;
        Q.push(make_pair(current->getChildCell(currentFace == currentNfold? (2*currentNfold+1):(currentFace+1)),4)) ;
      }
    }else {
      facePlan.push_back(0) ;
    }
    Q.pop() ;
  }
  while(facePlan.size() != 0 && facePlan.back() == 0 ) {
    facePlan.pop_back() ;
  }
  reduce_vector(facePlan) ;

  //clean up
  if(aCell != 0) {
    delete aCell ;
    aCell = 0 ;
  }
  cleanup_list(node_list, edge_list, face_list) ;
  return facePlan ;
}
