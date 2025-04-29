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
#include <Loci.h>
#include "diamondcell.h"
#include "defines.h"
#include "hex_defines.h"
using std::cerr;
using std::endl;
using std::swap;
using std::cout;

std::vector<char> transfer_plan_q2g(const std::vector<char>& facePlan);
std::vector<int32> contain_2d(const std::vector<pair<Range2d, int32> >& faceMap,
                              const std::vector<Range2d>& leaves);

struct Cell_Face{
  DiamondCell* c;
  int fi;
  Face* f;
  Cell_Face(DiamondCell* c1, int i, Face* f1):c(c1), fi(i), f(f1){};
};


std::vector<int32> get_c1(const Entity* lower, int lower_size,
                          const Entity* upper, int upper_size,
                          const Entity* boundary_map, int boundary_map_size,
                          const const_multiMap& face2node, 
                          const const_multiMap& face2edge,
                          const const_MapVec<2>& edge2node,
                          const std::vector<char>& cellPlan,
                          const std::vector<char>& facePlan,
                          Entity ff,
                          const const_store<int>& node_remap
                          ){

    
  
  std::vector<int32> c1;
  
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

 
  Face* aFace = aCell->face[findex];

 int num_face =  aFace->empty_resplit(facePlan);

  if(cellPlan.size() != 0){
    
    std::queue<Cell_Face> Q;
    std::queue<pair<int, Face*> >face_Q;
    
    //each node on ff corresponds to a child cell, and each
    //cell corresponds to a tree in diamonds, find  which tree correspond to
    // the node , and put it in the Q
    
    //At the same time , find out faceID in the childCell
    for(int i = 0; i < aFace->numEdge; i++){
      
      Node* theNode;
      if(aFace->needReverse[i]) theNode = aFace->edge[i]->tail;
      else theNode = aFace->edge[i]->head;
      
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
      Q.push(Cell_Face(aCell->child[childID], faceID, aFace->child[i]));
    }//the general face is splitted
   
    DiamondCell* current;
    int currentFace;
    Face* currentF;

   
    while(!Q.empty()){
      current = Q.front().c;
      currentFace = Q.front().fi;
      currentF = Q.front().f;
      int currentNfold = current->getNfold();
      if(current->getChildCell() != 0){
        if(currentF->child == 0){
          cerr << "WARNING: facePlan not consistent with cellPlan" << endl;
          Loci::Abort();
        }
        //push the selected four children in the Q
        if(faceOrient == 0){
          
          Q.push(Cell_Face(current->getChildCell(1), currentFace, currentF->child[0]));
          Q.push(Cell_Face(current->getChildCell(currentFace == currentNfold? (2*currentNfold+1):(currentFace+1)),
                           4, currentF->child[1]));
          Q.push(Cell_Face(current->getChildCell(currentFace-currentNfold+2), 4, currentF->child[2]));
          Q.push(Cell_Face(current->getChildCell(currentFace +2), 5, currentF->child[3]));
          
        }
        else{
          
          Q.push(Cell_Face(current->getChildCell(1), currentFace, currentF->child[0]));
          Q.push(Cell_Face(current->getChildCell(currentFace +2), 5, currentF->child[1]));
          Q.push(Cell_Face(current->getChildCell(currentFace-currentNfold+2), 4, currentF->child[2]));
          Q.push(Cell_Face(current->getChildCell(currentFace == currentNfold? (2*currentNfold+1):(currentFace+1)),
                           4, currentF->child[3]));
          
        }
      }
      
      else if(currentF->child != 0){
        
        Q.push(Cell_Face(current,currentFace, currentF->child[0]));
        Q.push(Cell_Face(current,currentFace, currentF->child[1]));
        Q.push(Cell_Face(current,currentFace, currentF->child[2]));
        Q.push(Cell_Face(current,currentFace, currentF->child[3]));    
      }
      else{
        c1.push_back(current->getCellIndex());
      }
      
      Q.pop();
      
    }
    
  }//end of if(cellPlan.size() != 0)
  else{
    
    for(int i = 0 ; i <num_face; i++){
      c1.push_back(int32(1));
    }
  }
  //clean up
  if(aCell != 0){
    delete aCell;
    aCell = 0;
  }
  cleanup_list(node_list, edge_list, face_list);
  
  reduce_vector(c1);
  return c1;
  
}

//need write another function get_c1_general from quadface to general cell
std::vector<int32> get_c1_general(const Entity* lower, int lower_size,
                                  const Entity* upper, int upper_size,
                                  const Entity* boundary_map, int boundary_map_size,
                                  bool is_quadface,
                                  const const_multiMap& face2node, 
                                  const const_multiMap& face2edge,
                                  const const_MapVec<2>& edge2node,
                                  const std::vector<char>& cellPlan,
                                  const std::vector<char>& facePlan,
                                  Entity ff,
                                  const const_store<int>& node_remap
                                  ){

  std::vector<int32> c1;
  if(is_quadface){
    std::vector<char> newPlan = transfer_plan_q2g(facePlan);
    std::vector<int32> tmp_c1 = get_c1(lower, lower_size,
                                       upper, upper_size,
                                       boundary_map, boundary_map_size,
                                       face2node,
                                       face2edge,
                                       edge2node,
                                       cellPlan,
                                       newPlan,
                                       ff,
                                       node_remap);
    
 
    
    //find the relation of quadface leaves and face leaves
    //transfer tmp_c1 to c1

    //first create a faceMap
    int64 maxX = int64(1) << MAXLEVEL;
    int64 maxY = int64(1) << MAXLEVEL;
    
    Range2d rg = Range2d(Point2d(0, 0), Point2d(maxX, maxY));
    //resplit face
    std::vector<pair<Range2d, int32> > faceMap(tmp_c1.size());
    int mapPointer = 0;
    std::queue<pair<Range2d, int> > Q;
    Q.push(make_pair(rg, 0));
   

    Range2d currentF;
     int orientCode;
    unsigned int index =0;
    char mySplitCode;
     int quad_childID[4]= {0, 2, 3, 1};//from general childID to quad childID
   

    while(!Q.empty()){
      currentF = Q.front().first;
      orientCode = Q.front().second;
      
      if(index >= newPlan.size()){
        mySplitCode = 0;
      }
      else{ 
        //take a code from splitcode
        mySplitCode = newPlan[index];
        index++;  
      }
      
     
      //the middle points used for children range
      //       p4     pmax
      //
      //p1     p5     p3
      // 
      //pmin   p2
      
     
      Point2d p1, p2, p3, p4, p5, pmin, pmax;
      
      
      
      //define positions
      pmin = currentF.minP;
      pmax = currentF.maxP;
      p1.x = pmin.x;
      p1.y = (pmin.y + pmax.y)/2;
      p2.x = (pmin.x + pmax.x)/2;
      p2.y = pmin.y;
      p3.x = pmax.x;
      p3.y = p1.y;
      p4.x = p2.x;
      p4.y = pmax.y;
      p5.x = p2.x;
      p5.y = p1.y;

      Range2d child[4] = {Range2d(pmin, p5),Range2d(p1, p4),Range2d(p2, p3),Range2d(p5, pmax)};
     
     //push the selected  children in the Q
      int childID;
      switch(mySplitCode){
      case 1:
        for(int i = 0; i < 4; i++){//i is general childID
          //   childID = (quad_childID[i] - orientCode +4)%4; //quad childID
          childID = quad_childID[(i + orientCode)%4];
          Q.push(make_pair(child[childID], (orientCode + i)%4));
        } 
        break;
     case 0:
       
       faceMap[mapPointer] = make_pair(currentF, tmp_c1[mapPointer]);
       mapPointer++;
       
       break;
     default:
       cerr << "illegal face code in get_c1_general() " << endl; 
       break;
      }
      
      Q.pop();
     
    }//while(!Q.empty())
  
    if(mapPointer != int(tmp_c1.size())){
      cerr << "WARNING: size of leaves and size of tmp_c1 not equal" << endl;
      exit(0);
    }
    //resplit rg with facePlan

    std::vector<Range2d> quad_leaves;
    std::queue<Range2d> quad_Q;
    quad_Q.push(rg);
    index = 0;      
    while(!quad_Q.empty()){
      currentF = quad_Q.front();
            
      if(index >= facePlan.size()){
        mySplitCode = 0;
      }
      else{ 
        //take a code from splitcode
        mySplitCode = facePlan[index];
        index++;  
      }
      
     
      //the middle points used for children range
      //       p4     pmax
      //
      //p1     p5     p3
      // 
      //pmin   p2
      
     
      Point2d p1, p2, p3, p4, p5, pmin, pmax;
      
      
      
      //define positions
      pmin = currentF.minP;
      pmax = currentF.maxP;
      p1.x = pmin.x;
      p1.y = (pmin.y + pmax.y)/2;
      p2.x = (pmin.x + pmax.x)/2;
      p2.y = pmin.y;
      p3.x = pmax.x;
      p3.y = p1.y;
      p4.x = p2.x;
      p4.y = pmax.y;
      p5.x = p2.x;
      p5.y = p1.y;

   
     
     //push the selected  children in the Q
   
      switch(mySplitCode){
      case 1:
        quad_Q.push( Range2d(pmin,p3));
        quad_Q.push( Range2d(p1, pmax));
        break;
       
     case 2:
       quad_Q.push(Range2d(pmin, p4));
       quad_Q.push(Range2d(p2, pmax));
       break;
       
      case 3:
        quad_Q.push(Range2d(pmin, p5));
        quad_Q.push(Range2d(p1, p4));
        quad_Q.push(Range2d(p2, p3));
        quad_Q.push(Range2d(p5, pmax));
        break;
      case 0:
       
        quad_leaves.push_back(currentF);
       
       break;
      default:
        cerr << "illegal face code in resplit(): " << endl; 
       break;
      }
      
      quad_Q.pop();
      
    }//while(!Q.empty())
    
    
    c1 = contain_2d(faceMap, quad_leaves);

  }
  else{
  
    c1 =  get_c1(lower, lower_size,
                 upper, upper_size,
                 boundary_map, boundary_map_size,
                 face2node,
                 face2edge,
                 edge2node,
                 cellPlan,
                  facePlan,
                 ff,
                 node_remap);
  }
  return c1;
}
