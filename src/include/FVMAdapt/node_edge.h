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
 * @ingroup fvmadapt_elements
 * @brief Node and binary-edge helper trees used by FVMAdapt refinement plans.
 *
 * For reference here:
 * head--------> tail, edges point from head to tail
 * edge2node[0]-------->edge2node[1]
 */


/**
 * @brief Node record used by the FVMAdapt helper-tree builders.
 * @ingroup fvmadapt_elements
 *
 * A Node stores the coordinates of a mesh vertex, along with an optional
 * one-based grid-file index and an adaptation tag. Nodes are used in edge,
 * face, and cell refinement trees.
 */
class Node{
public:
  /// Constructors
  Node():index(0),tag(0){}
  Node(vect3d& p0):p(p0), index(0), tag(0){}
  Node(const Loci::vector3d<double>& p0):p(p0),index(0),tag(0){}
  Node(vect3d& p0, int32 n):p(p0), index(n), tag(0){}
  Node(const Loci::vector3d<double>& p0, int32 n):p(p0), index(n), tag(0){}
public:
  /// Coordinate of the node
  vect3d p ;

  /// The index of node in input or output grid file, starts with 1
  int32 index ;

  /// Marker used by edge, face, and cell logic. 1 or 0, indicates if the
  /// node needs to be refined or not.
  char tag ;
};

/**
 * @brief Oriented mesh edge with a binary midpoint-refinement tree.
 * @ingroup fvmadapt_elements
 *
 * An Edge connects `head` to `tail` and may be split at its midpoint into two
 * child edges. Repeated splits form a binary tree that can be replayed from an
 * edge plan: code `1` splits the current edge, while code `0` leaves it as a
 * leaf.
 *
 * The reversed `resplit()` overload lets faces replay the same edge plan
 * in the orientation needed by their local edge ordering.
 *
 *
 * head-------->tail. The edge points from head to tail.
 * The vector is: tail->p - head->p
 *
 */
class Edge{
public:
  /// Constructors
  Edge(Node* p0, Node* p1):head(p0), tail(p1), child(0), parent(0),level(0) {}
  Edge(Node* p0, Node* p1, int lev):head(p0), tail(p1),child(0),parent(0), level(lev){}
  Edge(Node* p0, Node* p1, int lev, Edge* p):head(p0), tail(p1),child(0),parent(p), level(lev){}
  Edge():child(0), parent(0),level(0){}

  /// Destructor
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

  /// Calculate the middle point of the edge
  Node* centroid(){
    return new Node(0.5*(head->p + tail->p));
  }

  /// Euclidean edge length
  double length(){
    return norm(head->p - tail->p);
  }

  /// Duplicate of the length() method. TODO: One of these can be deprecated.
  double get_length(){
    return norm(tail->p - head->p);
  }

  /// Check if the edge contains the `theNode`.
  /// return 0 -> the node is the head
  /// return 1 -> the node is the tail
  /// return -1 -> the edge does not contain the node
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
  inline void split(std::list<Node*>& node_list){
    if(child == 0) {
      Node *center = centroid() ;
      node_list.push_back(center) ;

      child = new Edge*[2] ;
      child[0] = new Edge(head, center, level+1, this) ;
      child[1] = new Edge(center, tail, level+1, this) ;
    }
  }

  /// Return true if this edge has been split more than one level below itself.
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

  /// Return the maximum refinement-tree level below this edge.
  int get_depth() ;

  /**
  * Replay a breadth-first edge plan, optionally reversing child traversal.
  *
  * Each plan entry is applied to the next queued edge. Code `1` splits that
  * edge at its midpoint and queues its two children; code `0` leaves it as a
  * leaf. If `needReverse` is true, child edges are queued in reverse order.
  * New midpoint nodes are appended to `node_list`.
  *
  * @param edgePlan    Breadth-first edge refinement plan.
  * @param needReverse Whether to queue children in reverse edge order.
  * @param node_list   Receives midpoint nodes created by split().
  */
  void resplit(const std::vector<char>& edgePlan, bool needReverse,
               std::list<Node*>& node_list) ;

  /**
  * Replay a breadth-first edge plan in this edge's head-to-tail order.
  *
  * This overload is used when the edge orientation already matches the caller's
  * local ordering.
  *
  * @param edgePlan  Breadth-first edge refinement plan.
  * @param node_list Receives midpoint nodes created by split().
  */
  void resplit(const std::vector<char>& edgePlan,
               std::list<Node*>& node_list) ;

public:
  Node* head ;  ///< Start node in this edge's local orientation. edge2node[0]
  Node* tail ; ///< End node in this edge's local orientation. edge2node[1]
  Edge** child ; ///< Two midpoint-split child edges, or null for a leaf edge.
  Edge* parent ; ///< Parent edge in the refinement tree, or null for the root.
  int level ; ///< Refinement tree depth of this edge. Root edges a level 0
};

/**
 * Delete all Node and Edge pointers in the supplied lists and clear the lists.
 *
 * Each pointer is deleted at most once by this function, so callers must only
 * pass lists that own the objects they contain. Deleting an Edge also deletes
 * its child-edge subtree; Node objects referenced by edges must still be owned
 * and deleted through node_list.
 *
 * @param[in,out] node_list Owned node pointers to delete; cleared on return.
 * @param[in,out] edge_list Owned root edge pointers to delete; cleared on return.
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
 * Delete all Node pointers in the supplied list and clear the list.
 *
 * @param[in,out] node_list Owned node pointers to delete; cleared on return.
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


