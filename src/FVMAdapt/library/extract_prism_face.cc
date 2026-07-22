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
//************************************************************************
// this file extract the face refinement plan from the cell refinement
//plan. A queue is used to simulate the tree-building process but no tree is
//actually built.
// Two tables are used. childrenID table indicates the IDs of the children
//whose codes need to be extracted. faceCodeTable transfers from cell code to
//face code.
//***************************************************************************

#include <vector>
#include <queue>
#include "plan_operations.h"
#include "prism.h"
#include "tables.h"
using std::queue;
using std::cerr;
using std::endl;
using std::cout;
//using namespace std;


std::vector<char> extract_prism_face(const  std::vector<char>& cellPlan,  int dd){
  //Since when extract the code from childcell, the order can be 3, 0 or 2, 0,
  //it's more convenient to have an empty tree split
  std::vector<char> facePlan;
  if(cellPlan.size() == 0) {
    return facePlan;
  }

  Prism *aCell = new Prism(); // the root, nfold = 3
  aCell->empty_resplit(cellPlan); //only nfold and childCell is defined
  


  queue<pair<Prism*, int> > Q;
  Q.push(make_pair(aCell, dd));
  
  
  char cellCode, faceCode;
  Prism *current;
  int faceID;
  int nfold;
  while(!Q.empty()){
    current = Q.front().first;
    cellCode = current->mySplitCode;
    nfold = current-> getNfold();
    faceID = Q.front().second;
  
    //    cout << "cellCode: " << cellCode <<endl;
  
    if(cellCode == 0){
      facePlan.push_back(char(0));
      Q.pop();
      continue;
    }
    else if(faceID >= 2)faceCode = cellCode;
    else if(cellCode == 1)faceCode = 8; //for triface, facecode is 8 for cellcode 1
    else faceCode =1; //for triface, facecode is 1 for cellcode 2, 3
    
    facePlan.push_back(faceCode);
    //define chidID and child's faceID
    switch(cellCode){
    case 1:
      if(faceID == 0)Q.push(make_pair(current->childCell[0], faceID)); 
      else if (faceID == 1)Q.push(make_pair(current->childCell[1], faceID));
      else{
        Q.push(make_pair(current->childCell[0], faceID));
        Q.push(make_pair(current->childCell[1], faceID));
      }
      break;
      
    case 2:
      if(faceID < 2){
        for(int i = 0; i < nfold; i++){
          Q.push(make_pair(current->childCell[i], faceID)); 
        }
      }
      else{
        Q.push(make_pair(current->childCell[faceID-2], 2));
        Q.push(make_pair(current->childCell[(faceID-2)== nfold-1? 0:(faceID-1)], 5));
      }
      break;
      
    case 3:
      if(faceID == 0){
        for(int i = 0; i < nfold; i++){
          Q.push(make_pair(current->childCell[i], faceID)); 
        }
        
      }
      else if(faceID == 1){
        for(int i = 0; i < nfold; i++){
          Q.push(make_pair(current->childCell[i+nfold], faceID)); 
        }
     
      }
      else {
        Q.push(make_pair(current->childCell[faceID-2], 2));
        Q.push(make_pair(current->childCell[faceID-2+nfold], 2));
        Q.push(make_pair(current->childCell[(faceID-2)== nfold-1? 0:(faceID-1)], 5));
        Q.push(make_pair(current->childCell[(faceID-2)== nfold-1? nfold:(faceID-1+nfold)], 5));
      }
      break;
    default:
      cerr<<"WARNING: illegal cellCode" << endl;
      break;
      
      // cout << "cell: " << cellCode <<"  " <<"face: " << faceCode <<endl;
    }
  
    Q.pop();
  }
  //delete 0 and 8 in the beginning of faceplan
  while(facePlan.size() != 0 && ((facePlan.front() == 0)||(facePlan.front() == 8))){
    facePlan.erase(facePlan.begin());}
  //delete 0 and 8 at the end of faceplan
  while(facePlan.size() != 0 && ((facePlan.back() == 0)||(facePlan.back() == 8))){
    facePlan.pop_back();
  }
  //clean up
  delete aCell;
  return facePlan;
}




std::vector<char> merge_quad_face_pp(const std::vector<char>& cellPlan1,
                                     int dd1, char orientCode1,
                                     const std::vector<char>& cellPlan2,
                                     int dd2, char orientCode2){
  std::vector<char> facePlan1 = extract_prism_face(cellPlan1, dd1);
  std::vector<char> facePlan2 = extract_prism_face(cellPlan2, dd2);

  return merge_quad_face(facePlan1, orientCode1, facePlan2, orientCode2);
}




std::vector<char> merge_quad_face_p(const std::vector<char>& cellPlan1,
                                    int dd1, char orientCode1){
  std::vector<char> facePlan = extract_prism_face(cellPlan1, dd1);

  return merge_quad_face(facePlan, orientCode1);
}




std::vector<char>  merge_tri_face_pp(const  std::vector<char>& cellPlan1, int dd1, char orientCode1,
                                     const  std::vector<char>& cellPlan2, int dd2, char orientCode2){

  
  std::vector<char> plan1;
  std::vector<char> plan2;
  if(dd1>=2 || dd2 >= 2) {
    cerr<< "WARNING: illegal  faceID in merge_tri_face_pp" << endl;
    return plan1;
  }

 
  plan1 = extract_prism_face(cellPlan1, dd1);
  plan2 = extract_prism_face(cellPlan2, dd2);

  if(plan1.size() == 0 && plan2.size() == 0) return plan1;
  
  Face* aFace = new Face();
  aFace->numEdge = 3;
  aFace->empty_resplit(plan1, orientCode1);
  aFace->empty_resplit(plan2, orientCode2);
  plan1.clear();
  plan1 = aFace->make_faceplan();

  delete aFace;
  return plan1;
}
    
std::vector<char>  merge_tri_face_p(const  std::vector<char>& cellPlan1, int dd1, char orientCode1){
  
  
  std::vector<char> plan1;
  
  if(dd1>=2 ) {
    cerr<< "WARNING: illegal  faceID in merge_tri_face_pp" << endl;
    return plan1;
  }
  
 
  plan1 = extract_prism_face(cellPlan1,dd1);
  if(plan1.size() == 0) return plan1;
  
  Face* aFace = new Face();
  aFace->numEdge = 3;
  aFace->empty_resplit(plan1, orientCode1);
 
  plan1.clear();
  plan1 =  aFace->make_faceplan();
  delete aFace;
  return plan1;
}
