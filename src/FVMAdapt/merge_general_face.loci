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
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//  This file extract a facePlan from a cellPlan, and then merge two facePlan together
// the extracting algrithm uses the numbering system of DiamondCell and is faster than
//straight forward approach
//
///////////////////////////////////////////////////////////////////////////////////////




#include <queue>
#include <vector>
#include <utility>
#include <list>
#include "diamondcell.h"
#include <Loci.h>
#include <algorithm>
using std::list;
using std::queue;

using std::cerr;
using std::endl;
using std::cout;


std::vector<char>   extract_general_face(const Entity* lower, int lower_size,
                                         const Entity* upper, int upper_size,
                                         const Entity* boundary_map, int boundary_map_size,
                                         const const_multiMap& face2node,
                                         const const_multiMap& face2edge,
                                         const const_MapVec<2>& edge2node,
                                         const std::vector<char>& cellPlan,
                                         Entity ff,
                                         const const_store<int>& node_remap
                                         ){
 
  //if cellPlan empty, facePlan also empty
  std::vector<char> facePlan;
  if(cellPlan.size() == 0){
    reduce_vector(facePlan);
    return facePlan;
  }
  
  //first build a Cell from cellPlan
  std::list<Node*> node_list;
  std::list<Edge*> edge_list;
  std::list<Face*> face_list;
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
  //only child is defined
  aCell->empty_resplit(cellPlan);
  
 

  int findex = find_face_index(lower, lower_size,
                               upper, upper_size,
                               boundary_map, boundary_map_size,
                               face2node,
                               ff,
                               node_remap);


 char faceOrient = aCell->faceOrient[findex];
 
  facePlan.push_back(1);
  std::queue<pair<DiamondCell*, int> > Q;

  //each node on ff corresponds to a child cell, and each
  //cell corresponds to a tree in diamonds, find  which tree correspond to
  // the node , and put it in the Q
  
  //At the same time , find out faceID in the childCell
  Face* theFace = aCell->face[findex];
  
  for(int i = 0; i < theFace->numEdge; i++){

    Node* theNode;
    if(theFace->needReverse[i]) theNode = theFace->edge[i]->tail;
    else theNode = theFace->edge[i]->head;
    
    int childID = 0;
    
    for(childID = 0; childID < aCell->numNode; childID++){
      if(aCell->node[childID] == theNode) break;
    }
    if(childID == aCell->numNode){
      cerr << " WARNING: face nodes not in cell" << endl;
      Loci::Abort();
    }
    
    // to find faceID, first set up n2f and n2e
    std::vector<Face*> n2f;
    std::vector<Edge*> n2e;
    std::vector<int> rot;
    
    aCell->set_n2f_n2e(n2f, n2e, rot, childID); 

    //find faceID
    int faceID;
    int n2f_size = n2f.size();
    for( faceID = 0; faceID < n2f_size; faceID++){
      if(n2f[faceID] == aCell->face[findex]) break;
    }
    if(faceID == n2f_size){
      cerr<<" WARNING: can not find faceID" << endl;
      Loci::Abort();
    }
    
    //when split a general cell, center of face n2f[i] -> vertex[i+2],
    // vertex[1] and vertex[i+2] belongs to face nfold+i in the childCell
    
    faceID = faceID + aCell->child[childID]->getNfold();
    Q.push(make_pair(aCell->child[childID], faceID));
  }//the general face is splitted
  
  DiamondCell* current;
  int currentFace;
  while(!Q.empty()){
    current = Q.front().first;
    currentFace = Q.front().second;
    int currentNfold = current->getNfold();
    if(current->getChildCell() != 0){
      facePlan.push_back(1);
      //push the selected four children in the Q
      if(faceOrient == 0){
        
        Q.push(make_pair(current->getChildCell(1), currentFace));
        Q.push(make_pair(current->getChildCell(currentFace == currentNfold? (2*currentNfold+1):(currentFace+1)),4));
        Q.push(make_pair(current->getChildCell(currentFace-currentNfold+2), 4));
        Q.push(make_pair(current->getChildCell(currentFace +2), 5));
        
      }
      else{

        Q.push(make_pair(current->getChildCell(1), currentFace));
        Q.push(make_pair(current->getChildCell(currentFace +2), 5));
        Q.push(make_pair(current->getChildCell(currentFace-currentNfold+2), 4));
        Q.push(make_pair(current->getChildCell(currentFace == currentNfold? (2*currentNfold+1):(currentFace+1)),4));
        
      }
    }
    else{
      facePlan.push_back(0);
    }
    
    Q.pop();
  }
  while(facePlan.size() != 0 && facePlan.back() == 0 )facePlan.pop_back();
  reduce_vector(facePlan);

 //clean up
    if(aCell != 0){
      delete aCell;
      aCell = 0;
    }
    cleanup_list(node_list, edge_list, face_list);
   
    return facePlan;  
}




std::vector<char> merge_faceplan(std::vector<char>& planl, std::vector<char>& planr, int numNodes){
  if(planl.size() == 0) return planr;
  if(planr.size() == 0) return planl;
  std::vector<char> fplan;
  size_t ptl =0;
  size_t ptr =0;

  std::queue<pair<char, char> > Q;
  char codel;
  char coder;
 
  //assume the first code of both planl and planr is 1
  ptl++;
  ptr++;
  
  fplan.push_back(1);
  for(int  i = 0; i < numNodes; i++){
    if(ptl < planl.size()) {codel = planl[ptl];}
    else{ codel = 0;}
    if(ptr < planr.size()){ coder = planr[ptr];}
    else {coder = 0;}
    Q.push(make_pair(codel, coder));
    ptl++;
    ptr++;
  }

    
  while(!Q.empty()){
   
    if(Q.front().first ==1 && Q.front().second == 1) {
      fplan.push_back(1);
      
      for(int  i = 0; i < 4; i++){
        if(ptl < planl.size()){ codel = planl[ptl];}
        else {codel = 0;}
        if(ptr < planr.size()) {coder = planr[ptr];}
        else {coder = 0;}
        Q.push(make_pair(codel, coder));
        ptl++;
        ptr++;
      }
       
    }

    else if(Q.front().first == 1 && Q.front().second ==0){
      fplan.push_back(1);
      
      for(int  i = 0; i < 4; i++){
        if(ptl < planl.size()){ codel = planl[ptl];}
        else {codel = 0;}
        Q.push(make_pair(codel, 0));
        ptl++;
      }
        
    }

    else if(Q.front().first == 0 && Q.front().second ==1){
      fplan.push_back(1);

      for(int  i = 0; i < 4; i++){
        if(ptr < planr.size()) {coder = planr[ptr];}
        else {coder = 0;}
        Q.push(make_pair(0, coder));
        ptr++;
      }
         
    }
    else{
      fplan.push_back(0);
        
    }
    
    Q.pop();
    
  }

  while(fplan.size() != 0 && fplan.back() == 0 )fplan.pop_back();
  reduce_vector(fplan);
  

  return fplan;
}

class merge_general_interior_face:public pointwise_rule{
  const_Map cl;
  const_Map cr;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_store<vect3d> pos; //dummy
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;

  const_store<int> node_l2f;
public:
  merge_general_interior_face(){
    name_store("cl", cl);
    name_store("cr", cr);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("fileNumber(pos)", node_l2f);
   
    input("(cl,cr)->cellPlan");
    input("(cl, cr)->(lower, upper, boundary_map)->face2node->(pos,fileNumber(pos))");
   
    input("(cl, cr)->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    output("facePlan");
    constraint("(cl, cr)->gnrlcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
      
      do_loop(seq, this);
    }
  }
  void calculate(Entity f){
   
   
    
    facePlan[f].clear();
    std::vector<char> facePlanL = extract_general_face(lower[cl[f]].begin(), lower.num_elems(cl[f]),
                                                       upper[cl[f]].begin(), upper.num_elems(cl[f]),
                                                       boundary_map[cl[f]].begin(), boundary_map.num_elems(cl[f]),
                                                       face2node,
                                                       face2edge,
                                                       edge2node,
                                                       cellPlan[cl[f]],
                                                       f,
                                                       node_l2f);
    
    std::vector<char> facePlanR = extract_general_face(lower[cr[f]].begin(), lower.num_elems(cr[f]),
                                                       upper[cr[f]].begin(), upper.num_elems(cr[f]),
                                                       boundary_map[cr[f]].begin(), boundary_map.num_elems(cr[f]),
                                                       face2node,
                                                       face2edge,
                                                       edge2node,
                                                       cellPlan[cr[f]],
                                                       f, node_l2f);
    
    
    
    facePlan[f] = merge_faceplan(facePlanL, facePlanR, face2node.num_elems(f));
    reduce_vector(facePlan[f]);
  }
};

register_rule<merge_general_interior_face> register_merge_general_interior_face;



class merge_general_boundary_face:public pointwise_rule{
  const_Map cl;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_store<vect3d> pos;//dummy
  const_store<std::vector<char> > cellPlan;
  store<std::vector<char> > facePlan;

   const_store<int> node_l2f;
public:
  merge_general_boundary_face(){
    name_store("cl", cl);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("cellPlan", cellPlan);
    name_store("facePlan", facePlan);
    name_store("fileNumber(pos)", node_l2f);
    
    input("cl->cellPlan");
    input("cl->(lower, upper, boundary_map)->face2node->(pos,fileNumber(pos))");
   
    input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    output("facePlan");
    constraint("boundary_faces, cl->gnrlcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
  
    do_loop(seq, this);
    }
   
  }
  void calculate(Entity f){
   
  
    facePlan[f] = extract_general_face(lower[cl[f]].begin(), lower.num_elems(cl[f]),
                                       upper[cl[f]].begin(), upper.num_elems(cl[f]),
                                       boundary_map[cl[f]].begin(), boundary_map.num_elems(cl[f]),
                                       face2node,
                                       face2edge,
                                       edge2node,
                                       cellPlan[cl[f]],
                                       f,
                                       node_l2f);  
    reduce_vector(facePlan[f]);
  }
  
};

register_rule<merge_general_boundary_face> register_merge_general_boundary_face;
