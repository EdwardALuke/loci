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

/**
 * Replays an edge refinement plan with optional reversed child traversal.
 *
 * Code `1` splits the current edge and queues its children. Code `0`, or a
 * missing trailing code after the plan vector is exhausted, leaves the current
 * edge as a leaf. When @p needReverse is true, children are queued tail-to-head
 * so the caller can match the owning face orientation.
 *
 * @param edgePlan    Breadth-first edge refinement plan.
 * @param needReverse Whether to queue children in reverse edge order.
 * @param node_list   Receives midpoint nodes created by split().
 */
void Edge::resplit(const std::vector<char>& edgePlan, bool needReverse,
                   std::list<Node*>& node_list){
     if(edgePlan.size() == 0) {
      return;
    }
    
    std::queue<Edge*> Q;
    Q.push(this);
    
    Edge* current;
    unsigned int index = 0;
    char currentCode;
    
    while(!Q.empty()){
      current = Q.front();
      if(index >= edgePlan.size()){
      currentCode = 0;
      }
      else{ 
        //take a code from facePlan
        currentCode = edgePlan[index];
        index++;  
      }
      switch(currentCode)
        {
          
          //0 no split,this is a leaf, output faces
        case 0:
          break;
        case 1:
          current->split(node_list);
          if(needReverse){
            Q.push(current->child[1]);
            Q.push(current->child[0]);
          }
          else{
          Q.push(current->child[0]);
          Q.push(current->child[1]);
          }
          break;
        default:
          cerr <<"WARNING: illegal splitcode in function Edge::reSplit()" << endl;
          break;
        }
      Q.pop();
    }
}

/**
 * Replays an edge refinement plan in head-to-tail child order.
 *
 * This overload is used when the edge orientation already matches the caller's
 * local ordering.
 *
 * @param edgePlan  Breadth-first edge refinement plan.
 * @param node_list Receives midpoint nodes created by split().
 */
void Edge::resplit(const std::vector<char>& edgePlan,
                   std::list<Node*>& node_list){
  
  if(edgePlan.size() == 0) {
    return;
  }
  
  
 
  std::queue<Edge*> Q;
  Q.push(this);
  
  Edge* current;
  unsigned int index = 0;
  char currentCode;
 
  while(!Q.empty()){
    current = Q.front();

    if(index >= edgePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from facePlan
      currentCode = edgePlan[index];
      index++;  
    }
    
    
    switch(currentCode)
      {
        
        //0 no split,this is a leaf, output faces
      case 0:
        
        break;
        
      case 1:
        current->split(node_list);
        Q.push(current->child[0]);
        Q.push(current->child[1]);
        break;
        
      default:
        cerr <<"WARNING: illegal splitcode in function Edge::reSplit()" << endl;
        break;
      }
    
    Q.pop();
  }
  
}
  
/**
 * Appends leaf edges in head-to-tail tree order.
 *
 * The existing contents of @p leaves are preserved; this function only appends
 * leaves from the receiver.
 *
 * @param leaves Output list receiving leaf edges.
 */
void Edge::sort_leaves(std::list<Edge*>& leaves){
  if(child != 0){
    child[0]->sort_leaves(leaves);
    child[1]->sort_leaves(leaves);
  }
  else{
    leaves.push_back(this);
  }
}
  
 

int Edge::get_depth(){
  std::list<Edge*> leaves;
  int depth = 0;
  for(std::list<Edge*>::const_iterator p = leaves.begin(); p!=leaves.end(); p++){
    depth = max(depth, (*p)->level);
  }
  return depth;
}
