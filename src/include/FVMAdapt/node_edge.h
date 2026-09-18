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

#ifndef NODE_EDGE_H
#define NODE_EDGE_H

#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include <queue>
#include "defines.h"
using std::queue;

/**
 * @file node_edge.h
 *
 * Nodes and binary edge refinement trees. An Edge points from its stored head
 * to tail.
 */


/**
 * Node coordinates, a caller-assigned node index, and a refinement tag. Tag
 * values are 0 for unchanged, 1 for refinement, and 2 for derefinement.
 */
class Node{
public:
  Node():index(0),tag(0){}
  Node(vect3d& p0):p(p0), index(0), tag(0){}
  Node(const Loci::vector3d<double>& p0):p(p0),index(0),tag(0){}
  Node(vect3d& p0, int32 n):p(p0), index(n), tag(0){}
  Node(const Loci::vector3d<double>& p0, int32 n):p(p0), index(n), tag(0){}
public:
  /// Coordinate of the node
  vect3d p ;

  /// Node index assigned by the caller, for example from node_l2f or
  /// node_offset. Constructors initialize it to 0.
  int32 index ;

  /// Node tag: 0 unchanged, 1 refine, 2 derefine.
  char tag ;
};

/**
 * Edge from head to tail, split at its midpoint into two child edges. The
 * edge vector is tail->p - head->p.
 *
 * An edge plan uses code 1 to split the current edge and code 0 to request no
 * split. resplit() can reverse child traversal to match the edge direction
 * needed by the face.
 */
class Edge{
public:
  Edge(Node* p0, Node* p1):head(p0), tail(p1), child(0), parent(0),level(0) {}
  Edge(Node* p0, Node* p1, int lev):head(p0), tail(p1),child(0),parent(0), level(lev){}
  Edge(Node* p0, Node* p1, int lev, Edge* p):head(p0), tail(p1),child(0),parent(p), level(lev){}
  Edge():child(0), parent(0),level(0){}

  ~Edge(){
    if(child!=0){
      if(child[1] != 0){
        delete child[1];
        child[1] = 0;
      }
      if(child[0] != 0){
        delete child[0];
        child[0] = 0;
      }
      delete[] child;
      child = 0;
    }
  }

  /// Return a new Node at the edge midpoint. The caller deletes it.
  Node* centroid(){
    return new Node(0.5*(head->p + tail->p));
  }

  /// Return the edge length.
  double length(){
    return norm(head->p - tail->p);
  }

  /// Return the edge length, as in length().
  double get_length(){
    return norm(tail->p - head->p);
  }

  /// Return 0 if theNode is head, 1 if it is tail, or -1 if neither endpoint
  /// matches.
  int containNode(Node* theNode){
    if(head == theNode) return 0 ;
    if(tail == theNode) return 1 ;
    return -1 ;
  }

  /**
  * Split this edge at its midpoint if it is still a leaf.
  *
  * The new midpoint node is appended to @p node_list. The two child edges are
  * stored as `child[0]` from `head` to midpoint and `child[1]` from midpoint to
  * `tail`.
  *
  * @param[in,out] node_list Receives the midpoint node created by this split.
  */
  void split(std::list<Node*>& node_list){
    if(child == 0) {
      Node *center = centroid() ;
      node_list.push_back(center) ;

      child = new Edge*[2] ;
      child[0] = new Edge(head, center, level+1, this) ;
      child[1] = new Edge(center, tail, level+1, this) ;
    }
  }

  /// Return true if either child edge has children.
  bool depth_greater_than_1(){
    if(child == 0) return false ;
    if(child[0]->child != 0) return true ;
    if(child[1]->child != 0) return true ;
    return false ;
  }

  /**
  * Append this edge tree's leaf edges in child[0], then child[1] order.
  *
  * @param leaves In/out list; leaf edges are appended in child[0], child[1] order.
  */
  void sort_leaves(std::list<Edge*>& leaves) ;

  /// Currently returns 0; the implementation does not collect the leaf edges
  /// before checking their levels.
  int get_depth() ;

  /**
   * Apply edgePlan in breadth-first order. Code 1 splits the current edge and
   * queues its children; code 0 or an omitted entry requests no split.
   * Existing children are not removed.
   *
   * If needReverse is true, queue child[1] before child[0]. Otherwise use
   * head-to-tail order. Append new midpoint nodes to node_list.
   */
  void resplit(const std::vector<char>& edgePlan, bool needReverse,
               std::list<Node*>& node_list) ;

  /**
   * Apply edgePlan in breadth-first, head-to-tail order and append new
   * midpoint nodes to node_list. Existing children are not removed.
   */
  void resplit(const std::vector<char>& edgePlan,
               std::list<Node*>& node_list) ;

public:
  Node* head ;  ///< Start node in the stored edge direction.
  Node* tail ; ///< End node in the stored edge direction.
  Edge** child ; ///< Two midpoint-split child edges, or null for a leaf edge.
  Edge* parent ; ///< Parent edge in the refinement tree, or null for the root.
  int level ; ///< Stored level; children have level+1, roots can start above 0.
};

/**
 * Delete the listed nodes and root edges, then clear the lists. Each object
 * must appear only once. An Edge deletes its children, so child edges must
 * not also appear in edge_list. Edge nodes are deleted through node_list.
 */
inline void cleanup_list(std::list<Node*>& node_list,
                         std::list<Edge*>& edge_list) {

  for(std::list<Node*>::iterator p = node_list.begin(); p != node_list.end(); p++){
    if((*p) != 0) {
      delete (*p) ;
      (*p) = 0 ;
    }
  }
  node_list.clear() ;

  for(std::list<Edge*>::iterator p = edge_list.begin(); p != edge_list.end(); p++){
    if((*p) != 0) {
      delete (*p) ;
      (*p) = 0 ;
    }
  }
  edge_list.clear() ;
}

/**
 * Delete the listed nodes and clear the list. Each node must appear only
 * once.
 */
inline void cleanup_list(std::list<Node*>& node_list) {

  for(std::list<Node*>::iterator p = node_list.begin(); p != node_list.end(); p++){
    if((*p) != 0) {
      delete (*p) ;
      (*p) = 0 ;
    }
  }
  node_list.clear() ;
}

#endif
