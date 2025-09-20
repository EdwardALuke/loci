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
/////////////////////////////////////////////////////////////////////////////////////
//                                face.h
//
// this file include the declaration of class Face. Class Face is the abstraction
// of general face(polygon). It's defined as a collection of edges and the directions
// of edges. it supports only isotropic refinement.
// If the face is built as face2node, it can resplit without orientCode(as used in
// DiamondCell and  general Cell), A face can also be built as defined in cell, and
// resplit with orientCode(as used in Prism)

#ifndef FACE_H
#define FACE_H

#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include "node_edge.h"



//general face, i.e., polygon
class Face{
public:
  //constructors, it's important for child to be initialized to zero
  //for both memory management and for split

  Face(int n):numEdge(n),edge(new Edge*[n]),needReverse(new bool[n]),child(0){}
  
  Face(int n, Edge** e, bool* r):numEdge(n),edge(e), needReverse(r),child(0){}

  Face():numEdge(0),edge(0),needReverse(0),child(0){}
  
  //destructor
  ~Face(){
    if(child!= 0){
      for(int i = 0; i < numEdge; i++){
	if(child[i] != 0) {
	  delete child[i];
	  child[i] = 0;
	}
      }
      delete [] child;
      child = 0;
    }
    if(edge!=0){
      delete [] edge;
      edge = 0;
    }
    
    if(needReverse !=0){
      delete [] needReverse;
      needReverse = 0;
    }
  }

  inline Node* simple_center(){
    std::vector<vect3d> nodes(numEdge);
    for(int i = 0; i < numEdge; i++){
      nodes[i] =  edge[i]->head->p;
    }
    //calculate the mass center of the edge centers
    vect3d p = point_center(nodes);
    return new Node(p);
  }
    
  inline double area(){
    Node* c = simple_center();
    vect3d tmp_center = c->p;
    vect3d sum = vect3d(0.0, 0.0, 0.0);
    for(int i = 0; i < numEdge; i++){
      sum += cross((edge[i]->tail->p - tmp_center), (edge[i]->head->p - tmp_center));
    }
    if(c!=0)delete c;
    return 0.5*norm(sum);
  }
  
  //the center of the face, defined as the mass center of edge centers
  //precondition:: all its edges have been split
  inline Node* centroid(){
    switch(CENTROID){
    case 0:
      return simple_center();
    case 1:
      return wireframe();
    default:
      return wireframe();
    }
   
  }

  inline Node* wireframe(){
    
    //allocate edgecenter
    std::vector<vect3d> edgecenter(numEdge);
    std::vector<double> len(numEdge);
    
    //get edge centers
    for(int i = 0; i < numEdge; i++){
      edgecenter[i] = edge[i]->child[0]->tail->p;
      len[i] = edge[i]->length();
    }
   
    //calculate the mass center of the edge centers
    vect3d p = weighted_center(edgecenter, len);
    return new Node(p);
  }



  
  inline int  getLevel(){return edge[0]->level;};

  //precondition: all its edges have been split
  //condition: edgecenter must be allocated and deallocated by caller
  inline void getEdgeCenter(Node** edgecenter)const{
    for(int i = 0; i < numEdge; i++){
      edgecenter[i] = edge[i]->child[0]->tail;
    }
  }

  //check if one of the child is theFace
  //return  -1: no child is theFace
  //return i in range[0, numEdge): child[i] is theFace
  inline int containFace(Face* theFace)const{
    if(child !=0){
      for(int i = 0; i < numEdge; i++){
        if(child[i] == theFace){
          return i;
        }
      }
    }
    return -1;
  }
  
  //check if one of edge is theEdge
  //return  -1: no edge is theEdge
  //return i in range[0, numEdge): edge[i] is theEdge
  inline  int containEdge(Edge* theEdge)const{
    for(int i = 0; i < numEdge; i++){
      if(edge[i] == theEdge){
        return i;
      }
    }
    return -1;
  }

  //check if theNode is one of the face's vertexes
  //return -1: theNode is not one of the vertexes
  //return i in range[0, numEdge): theNode is f2n[i]
  inline int containNode(const Node* theNode)const{
    std::vector<Node*> f2n(numEdge);
    for(int i = 0; i < numEdge; i++){
      if(needReverse[i]) f2n[i] = edge[i]->tail;
      else f2n[i] = edge[i]->head;
    }
    int nodeID = -1;
    for(int i = 0; i < numEdge; i++){
      if(f2n[i] == theNode) {
        nodeID = i;
        break;
      }
    }
    return nodeID;
  }
  
  //get all the leaves of this
  void get_leaves(std::vector<Face*>& leaves) ;
 
  //define face2node
  void set_f2n(std::list<int32>& f2n);
  
  //split isotropically, all new nodes and edges are put into node_list and edge_list 
  void split(std::list<Node*>& node_list, std::list<Edge*>& edge_list);

  //for prism, when built prism cell, triangle face is built as node 0->1->2->0 and
  //node 3->4->->5->3, and it will split with orientCode
  void split(char orientCode, std::list<Node*>& node_list, std::list<Edge*>& edge_list);
 
  void empty_split();
  
  void resplit(const std::vector<char>& facePlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::vector<Face*>& fine_face);

  void resplit(const std::vector<char>& facePlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list);

  //for prism, when built prism cell, triangle face is built as node 0->1->2->0 and
  //node 3->4->->5->3, and it will split with orientCode
  void resplit(const std::vector<char>& facePlan,
               char orientCode,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list);

  //number of fine faces is returned 
  int empty_resplit(const std::vector<char>& facePlan);

  void empty_resplit(const std::vector<char>& facePlan, std::vector<Face*>& leaves);

  //this function is for merge_general_face_pp, 
  void empty_resplit(const std::vector<char>& facePlan, char orientCode);

  //compile the facePlan according the tree structure of aFace
  std::vector<char> make_faceplan();

  int get_num_leaves()const;//for mxfpc
public:
  int numEdge;
  Edge** edge;
  //if each edge is built as edge2node, and the face is built as face2node, needReverse is
  //true if face2node[i] == edge2node[face2edge[i]][1]

  //if each edge is defined as in prism, and the face is also defined as in prism,
  //needReverse is also false

  //during split, needReverse is decided by both the way edge is defined and the way the
  //face is defined
  bool* needReverse;
  
  Face** child;
};




//build a Face from Loci data structures, the locations of nodes are defined
//and edges are split according to edgePlan
Face* build_general_face( const Entity* face2node, int num_edge,
                          const Entity* face2edge,
                          const const_MapVec<2>& edge2node,
                          const const_store<vect3d>& pos,
                          const const_store<std::vector<char> >& edgePlan,
                          std::list<Node*>& bnode_list,
                          std::list<Edge*>& edge_list);

Face* build_general_face( const Entity* face2node, int num_edge,
                          const Entity* face2edge,
                          const const_MapVec<2>& edge2node,
                          const const_store<vect3d>& pos,
                          const const_store<int>& node_offset,
                          const const_store<std::vector<char> >& edgePlan,
                          std::list<Node*>& bnode_list,
                          std::list<Edge*>& edge_list,
                          const const_store<int>& node_remap);



//this function is used in build_general_cell with quadface
Face* build_tmp_general_face( const Entity* face2node, int num_edge,
                              const Entity* face2edge,
                              const const_MapVec<2>& edge2node,
                              const const_store<std::vector<char> >& edgePlan,
                              std::list<Node*>& bnode_list,
                              std::list<Edge*>& edge_list);


//if the intersection if leaves of f1 and the leaves of f2 is empty
bool is_overlapped( Face* f1,  Face* f2);

//for prism, used in merge_prism_face, the face is built as face2node
// and the plan is extracted from the cell.
int general_childID_orient_c2f(int childID_c, char orientCode, int numEdge);

//for prism, used in built_prism_cell and get_c1_prism,the face is built as in prism
// the facePlan is for face defined by face2node
int general_childID_orient_f2c(int childID_f, char orientCode, int numEdge);

//for prism, used in Face::split(), when face split with orientCode, to ensure the node_list
// generated is exact the same as split without orientCode, the edgeID need to be oriented 
int general_edgeID_orient_f2c(int i, char orientCode, int numEdge);


inline void cleanup_list(std::list<Node*>& node_list,
                         std::list<Edge*>& edge_list,
                         std::list<Face*>& face_list){
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

  for(std::list<Face*>::iterator p = face_list.begin();  p != face_list.end(); p++){
    if((*p) != 0){
      delete (*p);
      (*p) = 0;
    }
  }
  face_list.clear();
}

inline void cleanup_list( std::list<Face*>& face_list){
  for(std::list<Face*>::iterator p = face_list.begin();  p != face_list.end(); p++){
    if((*p) != 0){
      delete (*p);
      (*p) = 0;
    }
  }
  face_list.clear();
}


#endif



