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
//************************************************************************
// this file extract the edge refinement plan from the face refinement
//plan. A queue is used to simulate the tree-building process but no tree is
//actually built.
// Two tables are used. edgeIDTable indicates the IDs of the children
//whose codes need to be extracted. edgeCodeTable transfers from face code to
//edge code.
//***************************************************************************

#include <vector>
#include <queue>
#include <iostream>
#include "tables.h"

using std::queue;
using std::cerr;
using std::endl;
//using namespace std;

void  extract_quad_edge(const  std::vector<char>& facePlan, std::vector<char>& edgePlan, unsigned int dd){
  
  //output edgeCodeTable
  /* for(int i=0; i<12; i++){
    if(i%3==0) cout<< endl;
    cout << edgeCodeTable[i] <<',';
      }
  cout << endl << endl;
  */
  //output edgeIDTable
  /* for(int i=0; i<12; i++){
    if(i%3 == 0)cout <<endl; 
    cout <<'{';
      for(int j=0; j<edgeIDTable[i].size(); j++) cout<<edgeIDTable[i][j]<<',';
    cout <<'}';
   
  }
  cout << endl;
  */

  if(facePlan.size() == 0){
    edgePlan.clear();
    return;
  }
  unsigned int index = 0;  
  queue<bool> Q;
  Q.push(true);
  bool needExtract; 
  
  char  faceCode, edgeCode;
  while(!Q.empty()){
    needExtract = Q.front();
    if(index >= facePlan.size()){
      break;
    }
    else{
      faceCode = facePlan[index++];
    }
    //  cout << "faceCode: " << faceCode <<endl;
    if(needExtract){
      if(faceCode == 0){
        edgeCode = 0;
      }
      else{
        edgeCode = edgeCodeTable[dd*3+ faceCode-1];
      }
     
      edgePlan.push_back(edgeCode);
      
      //cout << "face: " << faceCode <<"  " <<"edge: " << edgeCode <<endl;
    }
    
    
    if(faceCode != 0){
      std::vector<bool> childrenID = edgeIDTable[dd*3+faceCode-1];
      if(!needExtract){
        for(unsigned int i = 0; i < childrenID.size(); i++)childrenID[i] = 0;
      }
      
      for(unsigned int i = 0; i < childrenID.size(); i++){
        Q.push(childrenID[i]);
      }
    }
    Q.pop();
  }
  //delete 0 and 8 in the beginning of edgeplan
  while((edgePlan.size() != 0) && ((edgePlan.front() == 0)||(edgePlan.front() == 8)))
    edgePlan.erase(edgePlan.begin());
  //delete 0 and 8 at the end of faceplan
  while((edgePlan.size() != 0) && ((edgePlan.back() == 0)||(edgePlan.back() == 8)))
    edgePlan.pop_back(); 
   
} 
  
