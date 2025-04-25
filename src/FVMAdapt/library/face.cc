//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <set>
#include <stack>
#include "face.h"
#include "defines.h"


using std::queue;
using std::cerr;
using std::endl;
using std::cout;

//for prism, used in Face::split(), when face split with orientCode, to ensure the node_list
// generated is exact the same as split without orientCode, the edgeID need to be oriented 

int general_edgeID_orient_f2c(int i, char orientCode, int numEdge){
  if(numEdge == 4){
    if(orientCode >=4) return 3-i;
    else return i;
  }
  else{//numEdge == 3

    if(orientCode/4 == 0)return (i+3 -orientCode)%3;
    else {
      switch(orientCode){
      case 4:
        return 2-i;
        break;
      case 5:
        return (3-i)%3;
        break;
      case 6:
        return (4-i)%3;
        break;
      }
      
    }
  }
  cerr << "WARNING:reach dummy code" << endl;
  exit(0);
  return 0 ;
}
      

//for prism only
int general_childID_orient_c2f(int childID_c, char orientCode, int numEdge){
  if(numEdge == 3){
    if(orientCode/4 == 0)return (childID_c +orientCode)%3;
    else return (2-childID_c+ orientCode)%3;
  }  
  else{
    
    if(orientCode/4 == 0)return childID_c ;
    else return (4-childID_c)%4;
  }
}



int general_childID_orient_f2c(int childID_f, char orientCode, int numEdge){
  if(orientCode == 3)cerr << "WARNING: illegal orientCode" << endl;
  if(numEdge == 3){
    if(orientCode/4 == 0)return (childID_f+3 -orientCode)%3;
    else return (2-childID_f+ orientCode)%3;
  }  
  else{
    
    if(orientCode/4 == 0)return childID_f ;
    else return (4-childID_f)%4;
  }
}

  

void Face::split(std::list<Node*>& node_list, std::list<Edge*>& edge_list){
  if(child != 0)return;
  //split edges
  for(int i = 0; i < numEdge; i++){
    if(edge[i]->child == 0) edge[i]->split(node_list);
  }
  
  //define facecenter is the center of edgecenters
  //calculate center and put it in node_list
  Node* facecenter = centroid();
  node_list.push_back(facecenter);
  
  //get edgecenter
  std::vector<Node*> edgecenter(numEdge);
  getEdgeCenter(&edgecenter[0]);
  
  
  //define new edges(from facecenter to edgecenetr) and put them into edge_list
  std::vector<Edge*> newEdges(numEdge);
  for(int i = 0; i < numEdge; i++){
    newEdges[i] = new Edge(facecenter, edgecenter[i], this->getLevel()+1);
    edge_list.push_back(newEdges[i]);
  }
  
  //define child
  child = new Face*[numEdge];
  for(int i = 0; i < numEdge; i++){
    child[i] = new Face(4);
   
    if(needReverse[i]) {
      child[i]->edge[0] = edge[i]->child[1];
      child[i]->needReverse[0] = true;
    }
    else{
      child[i]->edge[0] = edge[i]->child[0];
      child[i]->needReverse[0] = false;
    }

    child[i]->edge[1] = newEdges[i];
    child[i]->needReverse[1] = true;
    
    child[i]->edge[2] = newEdges[i==0?(numEdge-1):(i-1)];
    child[i]->needReverse[2] = false;

    if(needReverse[i==0?(numEdge-1):(i-1)]){
      child[i]->edge[3] = edge[i==0?(numEdge-1):(i-1)]->child[0];
      child[i]->needReverse[3] = true;
    }
    else{
      child[i]->edge[3] = edge[i==0?(numEdge-1):(i-1)]->child[1];
      child[i]->needReverse[3] = false;
    }
  }
}

//for build prism, when built prism cell, triangle face is built as node 0->1->2->0 and
//node 3->4->->5->3, and it will split with orientCode
void Face::split(char orientCode, std::list<Node*>& node_list, std::list<Edge*>& edge_list){
  if(child!=0)return;
  //split edges
  for(int i = 0; i < numEdge; i++){
    //build as in cell, and split as in face, so f2c
    int edgeID = general_edgeID_orient_f2c(i, orientCode, numEdge);
    if(edge[edgeID]->child == 0) edge[edgeID]->split(node_list);
  }
  
  //define facecenter is the center of edgecenters
  //calculate center and put it in node_list
  Node* facecenter = centroid();
  node_list.push_back(facecenter);


  //get edgecenter
  std::vector<Node*> edgecenter(numEdge);
  getEdgeCenter(&edgecenter[0]);
  
  
  //define new edges(from facecenter to edgecenetr) and put them into edge_list
  std::vector<Edge*> newEdges(numEdge);
  for(int i = 0; i < numEdge; i++){
    newEdges[i] = new Edge(facecenter, edgecenter[i], this->getLevel()+1);
    edge_list.push_back(newEdges[i]);
  }
  
  //define child
  child = new Face*[numEdge];
  for(int i = 0; i < numEdge; i++){
    child[i] = new Face(4);
   
    if(needReverse[i]) {
      child[i]->edge[0] = edge[i]->child[1];
      child[i]->needReverse[0] = true;
    }
    else{
      child[i]->edge[0] = edge[i]->child[0];
      child[i]->needReverse[0] = false;
    }

    child[i]->edge[1] = newEdges[i];
    child[i]->needReverse[1] = true;
    
    child[i]->edge[2] = newEdges[i==0?(numEdge-1):(i-1)];
    child[i]->needReverse[2] = false;

    if(needReverse[i==0?(numEdge-1):(i-1)]){
      child[i]->edge[3] = edge[i==0?(numEdge-1):(i-1)]->child[0];
      child[i]->needReverse[3] = true;
    }
    else{
      child[i]->edge[3] = edge[i==0?(numEdge-1):(i-1)]->child[1];
      child[i]->needReverse[3] = false;
    }
  }
}



int Face::get_num_leaves()const{
  if(child == 0) return 1;
  else {
    int count = 0;
    for(int i = 0; i < numEdge; i++){
      count += child[i]->get_num_leaves();
    }
    return count;
  }
}


void Face::empty_split(){
  //define child
  if(child == 0){
    child = new Face*[numEdge];
    for(int i = 0; i < numEdge; i++){
      child[i] = new Face();
      child[i]->numEdge = 4;
    }
  }
}

//define face2node
void Face::set_f2n(std::list<int32>& f2n){
  f2n.clear();
  for(int i = 0; i < numEdge; i++){

    //each edge sort leaves
    std::list<Edge*> edge_leaves;
    edge[i]->sort_leaves(edge_leaves);
   
    //if the edge needReverse, take the tail index value  
    if(needReverse[i]){
      for(std::list<Edge*>::reverse_iterator ep = edge_leaves.rbegin();
          ep != edge_leaves.rend(); ep++){
        f2n.push_back((*ep)->tail->index);
      }
    }
    else{//otherwise take the head index value
      for(std::list<Edge*>::iterator ep = edge_leaves.begin();
          ep != edge_leaves.end(); ep++){
        f2n.push_back((*ep)->head->index);
      } 
    } 
  }
}

// get all the leaves of this
void Face::get_leaves(std::vector<Face*>& leaves){
  if(child == 0){
    leaves.push_back(this);
    return;
  }
  else{
    for(int i = 0; i < numEdge; i++){
      child[i]->get_leaves(leaves);
    }
  }
}

//if the intersection if leaves of f1 and the leaves of f2 is empty
bool is_overlapped( Face* f1,  Face* f2){
  if(f1 == f2) return true;
  
  std::vector<Face*> leaves1;
  f1->get_leaves(leaves1);
  std::vector<Face*> leaves2;
  f2->get_leaves(leaves2);

  std::vector<Face*> intersection;
  std::sort(leaves1.begin(), leaves1.end());
  std::sort(leaves2.begin(), leaves2.end());
  std::set_intersection(leaves1.begin(), leaves1.end(), leaves2.begin(), leaves2.end(), std::inserter(intersection, intersection.begin()));
  return !(intersection.empty());
}




//compile the facePlan according the tree structure of aFace
std::vector<char> Face::make_faceplan(){
  std::vector<char> facePlan;
  std::queue<Face*> Q;
  Q.push(this);
  Face* current;
  while(!Q.empty()){
    current = Q.front();
    if(current->child != 0){
      facePlan.push_back(1);
      for(int i = 0; i < current->numEdge; i++){
        Q.push(current->child[i]);
      }
    }
    else{
      facePlan.push_back(0);
    }
    Q.pop();
  }                 
  while(facePlan.size() != 0 && facePlan.back() == 0) facePlan.pop_back();
  reduce_vector(facePlan);
  return facePlan;
}

//build a Face from Loci data structures, the locations of nodes are defined
//and edges are split according to edgePlan
Face* build_general_face( const Entity* face2node, int num_edge,
                          const Entity* face2edge,
                          const const_MapVec<2>& edge2node,
                          const const_store<vect3d>& pos,
                          const const_store<std::vector<char> >& edgePlan,
                          std::list<Node*>& bnode_list,
                          std::list<Edge*>& edge_list){
  

  std::vector<Node*> node(num_edge);
  
  for(int nindex = 0; nindex < num_edge; nindex++){
    node[nindex] = new Node(pos[face2node[nindex]]);
    bnode_list.push_back(node[nindex]);
  }
  
  //define each edge and put it into edge_list
  
  Edge** edge = new Edge*[num_edge];
  bool* needReverse = new bool[num_edge];
  
  for(int eindex = 0; eindex < num_edge; eindex++){
    //define the edge
    edge[eindex] = new Edge();
    edge_list.push_back(edge[eindex]);
    
    if(edge2node[face2edge[eindex]][0] == face2node[eindex]  && edge2node[face2edge[eindex]][1] == face2node[eindex==(num_edge-1)?0:eindex+1])
      {
        edge[eindex]->head = node[eindex];
        edge[eindex]->tail = node[eindex==(num_edge-1)?0:(eindex+1)];
        needReverse[eindex] = false;
      }
    else{
      edge[eindex]->tail = node[eindex];
      edge[eindex]->head = node[eindex==(num_edge-1)?0:(eindex+1)];
      needReverse[eindex] = true;
    }
    
    //replit the edge
    edge[eindex]->resplit(edgePlan[face2edge[eindex]], bnode_list);
  }
  
  //define the face
  Face* aFace = new Face(num_edge, edge, needReverse);
 
  return aFace;
}

//parallel version, build a face and index all the boundary nodes
Face* build_general_face( const Entity* face2node, int num_edge,
                          const Entity* face2edge,
                          const const_MapVec<2>& edge2node,
                          const const_store<vect3d>& pos,
                          const const_store<int>& node_offset,
                          const const_store<std::vector<char> >& edgePlan,
                          std::list<Node*>& bnode_list,
                          std::list<Edge*>& edge_list,
                          const const_store<int>& node_l2f){

  std::vector<Node*> node(num_edge);
  for(int nindex = 0; nindex < num_edge; nindex++){
   
    node[nindex] = new Node(pos[face2node[nindex]], node_l2f[face2node[nindex]]);
    bnode_list.push_back(node[nindex]);
  }
  
  //define each edge and put it into edge_list
  
  Edge** edge = new Edge*[num_edge];
  bool* needReverse = new bool[num_edge];
  
  //define edges and index its inner nodes
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
  
  for(int eindex = 0; eindex < num_edge; eindex++){
    //define the edge
    edge[eindex] = new Edge();
    edge_list.push_back(edge[eindex]);
    
    if(edge2node[face2edge[eindex]][0] == face2node[eindex] && edge2node[face2edge[eindex]][1] == face2node[eindex==(num_edge-1)?0:eindex+1])
      {
        edge[eindex]->head = node[eindex];
        edge[eindex]->tail = node[eindex==(num_edge-1)?0:(eindex+1)];
        needReverse[eindex] = false;
      }
    else{
      edge[eindex]->tail = node[eindex];
      edge[eindex]->head = node[eindex==(num_edge-1)?0:(eindex+1)];
      needReverse[eindex] = true;
    }
    
    //replit the edge
    edge[eindex]->resplit(edgePlan[face2edge[eindex]], bnode_list);
    int nindex = node_offset[face2edge[eindex]];
    
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->index =  nindex++;
    }
    
    bnode_begin = --(bnode_list.end());
    
  }
  
  //define the face
  Face* aFace = new Face(num_edge, edge, needReverse);
  return aFace;
}

//this function is used in build_general_cell with quadface
Face* build_tmp_general_face( const Entity* face2node, int num_edge,
                              const Entity* face2edge,
                              const const_MapVec<2>& edge2node,
                              const const_store<std::vector<char> >& edgePlan,
                              std::list<Node*>& bnode_list,
                              std::list<Edge*>& edge_list){


  
  std::vector<Node*> node(num_edge);

  
 
  vect3d p[4];
  int64 maxX = int64(1) << MAXLEVEL;
  int64 maxY = int64(1) << MAXLEVEL;
  p[0] = vect3d(0.0, 0.0, 0.0);
  p[1] = vect3d(maxX, 0.0, 0.0);
  p[2] = vect3d(maxX, maxY, 0.0);
  p[3] = vect3d(0.0, maxY, 0.0);
    for(int nindex = 0; nindex < num_edge; nindex++){
   
    node[nindex] = new Node(p[nindex]);
    bnode_list.push_back(node[nindex]);
  }
  //define each edge and put it into edge_list
  
  Edge** edge = new Edge*[num_edge];
  bool* needReverse = new bool[num_edge];
  
    
  for(int eindex = 0; eindex < num_edge; eindex++){
    //define the edge
    edge[eindex] = new Edge();
    edge_list.push_back(edge[eindex]);
    
    if(edge2node[face2edge[eindex]][0] == face2node[eindex] && edge2node[face2edge[eindex]][1] == face2node[eindex==(num_edge-1)?0:eindex+1])
      {
        edge[eindex]->head = node[eindex];
        edge[eindex]->tail = node[eindex==(num_edge-1)?0:(eindex+1)];
        needReverse[eindex] = false;
      }
    else{
      edge[eindex]->tail = node[eindex];
      edge[eindex]->head = node[eindex==(num_edge-1)?0:(eindex+1)];
      needReverse[eindex] = true;
    }
    
    //replit the edge
    edge[eindex]->resplit(edgePlan[face2edge[eindex]], bnode_list);
  }
  
  //define the face
  Face* aFace = new Face(num_edge, edge, needReverse);
  return aFace;
}
















//this function split  a general face according to facePlan,
//all fine quadface in faces
void Face::resplit(const std::vector<char>& facePlan,
                   std::list<Node*>& node_list,
                   std::list<Edge*>& edge_list,
                    std::vector<Face*>& fine_face){
  
                          
  
  if(facePlan.size() == 0) {
    fine_face.push_back(this);
    reduce_vector(fine_face);
    return;
  }
  

  //assume the first code in facePlan is 1
  std::queue<Face*> Q;
  Q.push(this);
  
  Face* current;
  unsigned int index = 0;
  char currentCode;

  while(!Q.empty()){
    current = Q.front();
    if(index >= facePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from facePlan
      currentCode = facePlan[index];
      index++;  
    }
    
    
    switch(currentCode)
      {
        
        //0 no split,this is a leaf, output faces
      case 0:
        fine_face.push_back(current);
        break;
        
      case 1:
        current->split(node_list, edge_list);
        
        for(int i = 0; i < current->numEdge; i++){
          Q.push(current->child[i]); 
        }        
        break;
        
      default:
        cerr <<"WARNING: illegal splitcode in function Face::resplit()" << endl;
     
        break;
      }
    
    Q.pop();
  }
  reduce_vector(fine_face);
}
//this function split  a general face according to facePlan,

void Face::resplit(const std::vector<char>& facePlan,
                   std::list<Node*>& node_list,
                   std::list<Edge*>& edge_list){
  
                          
  
  if(facePlan.size() == 0) {
    return;
  }
  

  //assume the first code in facePlan is 1
  std::queue<Face*> Q;
  Q.push(this);
  
  Face* current;
  unsigned int index = 0;
  char currentCode;

  while(!Q.empty()){
    current = Q.front();
    if(index >= facePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from facePlan
      currentCode = facePlan[index];
      index++;  
    }
    
    
    switch(currentCode)
      {
        
        //0 no split,this is a leaf, output faces
      case 0:
        break;
        
      case 1:
        current->split(node_list, edge_list);
        
        for(int i = 0; i < current->numEdge; i++){
          Q.push(current->child[i]); 
        }        
        break;
        
      default:
        cerr <<"WARNING: illegal splitcode in function Face::resplit()" << endl;
     
        break;
      }
    
    Q.pop();
  }
 
}
//for build prism when built prism cell, triangle face is built as node 0->1->2->0 and
//node 3->4->->5->3, and it will split with orientCode
void Face::resplit(const std::vector<char>& facePlan,
                   char orientCode,
                   std::list<Node*>& node_list,
                   std::list<Edge*>& edge_list){
  
                          
  
  if(facePlan.size() == 0) {
    return;
  }
  
  std::queue<Face*> Q;
  Q.push(this);
  
  Face* current;
  unsigned int index = 0;
  char currentCode;
  
  while(!Q.empty()){
    current = Q.front();
    if(index >= facePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from facePlan
      currentCode = facePlan[index];
      index++;  
    }
    
    
    switch(currentCode)
      {
        
        //0 no split,this is a leaf, output faces
      case 0:
        
        break;
        
      case 1:
        current->split(orientCode, node_list, edge_list);
        //build as in cell, and  split as in face, so f2c
        for(int i = 0; i < current->numEdge; i++){
          Q.push(current->child[general_childID_orient_f2c(i, orientCode, current->numEdge)]); 
        }        
        break;
        
      default:
        cerr <<"WARNING: illegal splitcode in function Face::reSplit(char orientCode)" << endl;
        break;
      }
    
    Q.pop();
  }

}


void  Face::empty_resplit(const std::vector<char>& facePlan,
                          std::vector<Face*>& leaves){
  if(facePlan.size() == 0) {
    leaves.push_back(this);
    reduce_vector(leaves);
    return ;
  }
  
  
  //assume the first code in facePlan is 1
  std::queue<Face*> Q;
  Q.push(this);
  
  Face* current;
  unsigned int index = 0;
  char currentCode;

  while(!Q.empty()){
    current = Q.front();
    if(index >= facePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from facePlan
      currentCode = facePlan[index];
      index++;  
    }
    
    
    switch(currentCode)
      {
        
        //0 no split,this is a leaf, output faces
      case 0:
        leaves.push_back(current);
        break;
        
      case 1:
        current->empty_split();
        
        for(int i = 0; i < current->numEdge; i++){
          Q.push(current->child[i]); 
        }        
        break;
        
      default:
        cerr <<"WARNING: illegal splitcode in function Face::empty_resplit()" << endl;
        break;
      }
    
    Q.pop();
  }
  return;
}


//this function split  a general face according to facePlan,

int  Face::empty_resplit(const std::vector<char>& facePlan){
  
  int num_face = 0;
  
  if(facePlan.size() == 0) {
    return 1;
  }
  
  
  //assume the first code in facePlan is 1
  std::queue<Face*> Q;
  Q.push(this);
  
  Face* current;
  unsigned int index = 0;
  char currentCode;

  while(!Q.empty()){
    current = Q.front();
    if(index >= facePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from facePlan
      currentCode = facePlan[index];
      index++;  
    }
    
    
    switch(currentCode)
      {
        
        //0 no split,this is a leaf, output faces
      case 0:
        num_face++;
        break;
        
      case 1:
        current->empty_split();
        
        for(int i = 0; i < current->numEdge; i++){
          Q.push(current->child[i]); 
        }        
        break;
        
      default:
        cerr <<"WARNING: illegal splitcode in function Face::empty_resplit()" << endl;
      
        break;
      }
    
    Q.pop();
  }
  return num_face;
}


//this function is for merge_general_face_pp, 
void Face::empty_resplit(const std::vector<char>& facePlan, char orientCode){
  if(facePlan.size() == 0) {
    return;
  }
  
  std::queue<Face*> Q;
  Q.push(this);
  
  Face* current;
  unsigned int index = 0;
  char currentCode;

  while(!Q.empty()){
    current = Q.front();
    if(index >= facePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from facePlan
      currentCode = facePlan[index];
      index++;  
    }
    
    
    if(currentCode == 1){
      current->empty_split();
      for(int i = 0; i < current->numEdge; i++){
        Q.push(current->child[general_childID_orient_c2f(i, orientCode, current->numEdge)]); 
      }        
    }
    else if(currentCode == 8) Q.push(current);
    
    Q.pop();
  }
  return ;
}



