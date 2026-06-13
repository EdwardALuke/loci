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
 * head--------> tail, edge point from head to tail
 * edge2node[0]-------->edge2node[1]
 */

//declaration of class Node, Edge and Face
/**
 * @brief Point record used by the FVMAdapt helper-tree builders.
 * @ingroup fvmadapt_elements
 *
 * `index` is the one-based node id used when reading or writing grid data.
 * `tag` is used by selected tagging/resplit paths; it is not a general plan
 * code.
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
  /// Coordinates
  vect3d p;

  /// The index of node in input or output grid file, starts with 1
  int32 index;

  /// 1 or 0, indicates if the node needs to be refined or not
  char tag;
};

/**
 * @brief Binary refinement tree for one mesh edge.
 * @ingroup fvmadapt_elements
 *
 * Edge plans use code `0` for a leaf and `1` for a midpoint split. Child order
 * can be replayed in the owning face direction through the `needReverse`
 * overload of resplit().
 *
 * head-------->tail. The edge points from head to tail.
 * The vector is: tail->p - head->p
 *
 * @warning The destructor owns and deletes child edges, but not the `head` and
 * `tail` nodes. Builder code usually tracks allocated nodes separately in
 * lists that are later passed to cleanup_list().
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
  inline Node* centroid(){
    return new Node(0.5*(head->p + tail->p));
  }

  inline double length(){
    return norm(head->p - tail->p);
  }

  inline double get_length(){
    return norm(tail->p - head->p);
  }

  /// Check if the edge contains the theNode.
  /// return 0 -> the node is the head
  /// return 1 -> the node is the tail
  /// return -1 -> the edge does not contain the node
  inline  int containNode(Node* theNode){
    if(head == theNode) return 0;
    if(tail == theNode) return 1;
    return -1;
  }

  /// Split function, split the edge in the middle
  inline void split(std::list<Node*>& node_list){
    if(child == 0){
      Node *center = centroid();
      node_list.push_back(center);

      child = new Edge*[2];
      child[0] = new Edge(head, center, level+1, this);
      child[1] = new Edge(center, tail, level+1, this);
    }
  }

  inline bool depth_greater_than_1(){
    if(child == 0) return false;
    if(child[0]->child != 0) return true;
    if(child[1]->child != 0) return true;
    return false;
  }

  void sort_leaves(std::list<Edge*>& leaves);

  int get_depth();

  // For QuadFace
  void resplit(const std::vector<char>& edgePlan, bool needReverse,
               std::list<Node*>& node_list);

  // For Face
  void resplit(const std::vector<char>& edgePlan,
               std::list<Node*>& node_list);

public:
  Node* head;  //edge2node[0]
  Node* tail; //edge2node[1]
  Edge** child;
  Edge* parent;

  /// The level of tree structure
  int level;
};


inline void cleanup_list(std::list<Node*>& node_list,
                         std::list<Edge*>& edge_list){

  for(std::list<Node*>::iterator p = node_list.begin(); p != node_list.end(); p++){
    if((*p) != 0){
      delete (*p);
      (*p) = 0;
    }
  }
  node_list.clear();

  for(std::list<Edge*>::iterator p = edge_list.begin(); p != edge_list.end(); p++){
    if((*p) != 0){
      delete (*p);
      (*p) = 0;
    }
  }
  edge_list.clear();

}

inline void cleanup_list(std::list<Node*>& node_list){

  for(std::list<Node*>::iterator p = node_list.begin(); p != node_list.end(); p++){
    if((*p) != 0){
      delete (*p);
      (*p) = 0;
    }
  }
  node_list.clear();
}

#endif


