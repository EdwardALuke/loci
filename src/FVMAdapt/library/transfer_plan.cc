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
#include <utility>
#include <Loci.h>
#include "hex_defines.h"
#include "defines.h"
#include "face.h"
#include "quadface.h"

using std::vector;
using std::queue;

using std::cerr;
using std::endl;
using std::cout;
std::vector<char> transfer_plan_g2q(std::vector<char>& facePlan){
  //first built a tree, empty split it
  Face* aFace = new Face();
  aFace->numEdge = 4;
  aFace->empty_resplit(facePlan);

  std::vector<char> newFacePlan;
  int general_childID[4]= {0, 3, 1, 2};//from quad childID to general childID
  std::queue<pair<Face*, int> > Q;
  Q.push(make_pair(aFace, 0));

  Face* current;
  int orientCode;
 
  
  while(!Q.empty()){
    current = Q.front().first;
    orientCode = Q.front().second;
    if(current->child == 0) newFacePlan.push_back(0);
    else{
      newFacePlan.push_back(3);
      int childID ; 
      for(int i = 0; i< 4; i++){
        childID = (general_childID[i] - orientCode +4)%4;
        // Q.push(make_pair(current->child[childID], (orientCode + general_childID[i])%4));
       Q.push(make_pair(current->child[childID], (orientCode + childID)%4));  
      }
    }
    Q.pop();
  }
  while(newFacePlan.size() != 0 && newFacePlan.back() == 0) newFacePlan.pop_back();
  reduce_vector(newFacePlan);
  delete aFace;
  return newFacePlan;
}



std::vector<char> transfer_plan_q2g(const std::vector<char>& facePlan){
  

 

  //built an empty tree
  QuadFace* aFace = new QuadFace(4);
  std::vector<QuadFace*> fine_faces;  
  aFace->empty_resplit(facePlan, 0, fine_faces);
  //rewrite into quadPlan


 
 std::vector<char> newFacePlan;
  int quadID[4]= {0, 2, 3, 1};//from general childID to quad childID
  std::queue<pair<QuadFace*, int> > Q;
  Q.push(make_pair(aFace, 0));

  QuadFace* current;
  int orientCode;
 
  
  while(!Q.empty()){
    current = Q.front().first;
    orientCode = Q.front().second;
    if(current->code != 3) newFacePlan.push_back(0);
    else{
      newFacePlan.push_back(1);
      int childID ; 
      for(int i = 0; i< 4; i++){
        childID = quadID[(i+orientCode)%4];
        Q.push(make_pair(current->child[childID], (orientCode + i)%4));
      }
    }
    Q.pop();
  }
  while(newFacePlan.size() != 0 && newFacePlan.back() == 0) newFacePlan.pop_back();
  reduce_vector(newFacePlan);
  delete aFace;


  return newFacePlan;
  
 
}
