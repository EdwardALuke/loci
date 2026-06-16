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
#include <iostream>
#include <vector>
#include "node_edge.h"

/**
 * @file node_edge.cc
 * @brief Edge-tree operations used by FVMAdapt refinement plans.
 *
 * Edges are represented as binary trees. A split creates two child edges and a
 * midpoint node; resplitting replays a stored edge plan in breadth-first order.
 */


void Edge::resplit(const std::vector<char>& edgePlan, bool needReverse,
                   std::list<Node*>& node_list) {
     if(edgePlan.size() == 0) {
      return ;
    }

    std::queue<Edge*> Q;
    Q.push(this) ;

    Edge* current ;
    unsigned int index = 0 ;
    char currentCode ;

    while(!Q.empty()) {
      current = Q.front() ;
      if(index >= edgePlan.size()) {
        currentCode = 0 ;
      }else {
        // take a code from edgePlan
        currentCode = edgePlan[index] ;
        index++ ;
      }
      switch(currentCode) {
        // 0 no split, this is a leaf
        case 0:
          break ;
        case 1:
          current->split(node_list) ;
          if(needReverse) {
            Q.push(current->child[1]) ;
            Q.push(current->child[0]) ;
          }else {
            Q.push(current->child[0]) ;
            Q.push(current->child[1]) ;
          }
          break ;
        default:
          cerr << "WARNING: illegal splitcode in function Edge::reSplit()" << endl ;
          break ;
      }
      Q.pop() ;
    }
}

void Edge::resplit(const std::vector<char>& edgePlan,
                   std::list<Node*>& node_list) {

  if(edgePlan.size() == 0) {
    return ;
  }

  std::queue<Edge*> Q ;
  Q.push(this) ;

  Edge* current ;
  unsigned int index = 0 ;
  char currentCode ;

  while(!Q.empty()) {
    current = Q.front() ;

    if(index >= edgePlan.size()) {
      currentCode = 0 ;
    }else {
      // take a code from edgePlan
      currentCode = edgePlan[index] ;
      index++ ;
    }

    switch(currentCode) {
      // 0 no split, this is a leaf
      case 0:
        break ;

      case 1:
        current->split(node_list) ;
        Q.push(current->child[0]) ;
        Q.push(current->child[1]) ;
        break ;

      default:
        cerr << "WARNING: illegal splitcode in function Edge::reSplit()" << endl ;
        break ;
    }
    Q.pop() ;
  }
}

void Edge::sort_leaves(std::list<Edge*>& leaves) {
  if(child != 0) {
    child[0]->sort_leaves(leaves) ;
    child[1]->sort_leaves(leaves) ;
  }else {
    leaves.push_back(this) ;
  }
}

int Edge::get_depth() {
  std::list<Edge*> leaves ;
  int depth = 0 ;
  for(std::list<Edge*>::const_iterator p = leaves.begin(); p!=leaves.end(); p++) {
    depth = max(depth, (*p)->level) ;
  }
  return depth ;
}
