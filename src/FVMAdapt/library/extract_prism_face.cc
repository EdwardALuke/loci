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


std::vector<char> extract_prism_face(const std::vector<char>& cellPlan,
                                     int faceID){
  //Since when extract the code from childcell, the order can be 3, 0 or 2, 0,
  //it's more convenient to have an empty tree split
  std::vector<char> facePlan;
  if(cellPlan.size() == 0) {
    return facePlan;
  }

  Prism *rootCell = new Prism(); // the root, nfold = 3
  rootCell->empty_resplit(cellPlan); //only nfold and childCell is defined
  


  queue<pair<Prism*, int> > pendingFaces;
  pendingFaces.push(make_pair(rootCell, faceID));
  
  
  char cellSplitCode, faceSplitCode;
  Prism *currentCell;
  int currentFaceID;
  int currentNfold;
  while(!pendingFaces.empty()){
    currentCell = pendingFaces.front().first;
    cellSplitCode = currentCell->mySplitCode;
    currentNfold = currentCell-> getNfold();
    currentFaceID = pendingFaces.front().second;
  
    //    cout << "cellSplitCode: " << cellSplitCode <<endl;
  
    if(cellSplitCode == 0){
      facePlan.push_back(char(0));
      pendingFaces.pop();
      continue;
    }
    else if(currentFaceID >= 2)faceSplitCode = cellSplitCode;
    else if(cellSplitCode == 1)faceSplitCode = 8; //for triface, facecode is 8 for cellcode 1
    else faceSplitCode =1; //for triface, facecode is 1 for cellcode 2, 3
    
    facePlan.push_back(faceSplitCode);
    //define chidID and child's faceID
    switch(cellSplitCode){
    case 1:
      if(currentFaceID == 0)
        pendingFaces.push(make_pair(currentCell->childCell[0], currentFaceID));
      else if (currentFaceID == 1)
        pendingFaces.push(make_pair(currentCell->childCell[1], currentFaceID));
      else{
        pendingFaces.push(make_pair(currentCell->childCell[0], currentFaceID));
        pendingFaces.push(make_pair(currentCell->childCell[1], currentFaceID));
      }
      break;
      
    case 2:
      if(currentFaceID < 2){
        for(int i = 0; i < currentNfold; i++){
          pendingFaces.push(make_pair(currentCell->childCell[i], currentFaceID));
        }
      }
      else{
        pendingFaces.push(make_pair(currentCell->childCell[currentFaceID-2], 2));
        pendingFaces.push(make_pair(
          currentCell->childCell[(currentFaceID-2)== currentNfold-1?
                                 0:(currentFaceID-1)], 5));
      }
      break;
      
    case 3:
      if(currentFaceID == 0){
        for(int i = 0; i < currentNfold; i++){
          pendingFaces.push(make_pair(currentCell->childCell[i], currentFaceID));
        }
        
      }
      else if(currentFaceID == 1){
        for(int i = 0; i < currentNfold; i++){
          pendingFaces.push(make_pair(currentCell->childCell[i+currentNfold], currentFaceID));
        }
     
      }
      else {
        pendingFaces.push(make_pair(currentCell->childCell[currentFaceID-2], 2));
        pendingFaces.push(make_pair(currentCell->childCell[currentFaceID-2+currentNfold], 2));
        pendingFaces.push(make_pair(
          currentCell->childCell[(currentFaceID-2)== currentNfold-1?
                                 0:(currentFaceID-1)], 5));
        pendingFaces.push(make_pair(
          currentCell->childCell[(currentFaceID-2)== currentNfold-1?
                                 currentNfold:(currentFaceID-1+currentNfold)],
          5));
      }
      break;
    default:
      cerr<<"WARNING: illegal cell split code" << endl;
      break;
      
      // cout << "cell: " << cellSplitCode <<"  " <<"face: " << faceSplitCode <<endl;
    }
  
    pendingFaces.pop();
  }
  //delete 0 and 8 in the beginning of faceplan
  while(facePlan.size() != 0 && ((facePlan.front() == 0)||(facePlan.front() == 8))){
    facePlan.erase(facePlan.begin());}
  //delete 0 and 8 at the end of faceplan
  while(facePlan.size() != 0 && ((facePlan.back() == 0)||(facePlan.back() == 8))){
    facePlan.pop_back();
  }
  //clean up
  delete rootCell;
  return facePlan;
}




std::vector<char> merge_quad_face_pp(const std::vector<char>& cellPlanL,
                                     int faceIDL, char orientCodeL,
                                     const std::vector<char>& cellPlanR,
                                     int faceIDR, char orientCodeR){
  std::vector<char> facePlanL = extract_prism_face(cellPlanL, faceIDL);
  std::vector<char> facePlanR = extract_prism_face(cellPlanR, faceIDR);

  return merge_quad_face(facePlanL, orientCodeL, facePlanR, orientCodeR);
}




std::vector<char> merge_quad_face_p(const std::vector<char>& cellPlan,
                                    int faceID, char orientCode){
  std::vector<char> facePlan = extract_prism_face(cellPlan, faceID);

  return merge_quad_face(facePlan, orientCode);
}




std::vector<char> merge_tri_face_pp(const std::vector<char>& cellPlanL,
                                    int faceIDL, char orientCodeL,
                                    const std::vector<char>& cellPlanR,
                                    int faceIDR, char orientCodeR){

  
  std::vector<char> facePlanL;
  std::vector<char> facePlanR;
  if(faceIDL>=2 || faceIDR >= 2) {
    cerr << "WARNING: illegal face ID in merge_tri_face_pp" << endl;
    return facePlanL;
  }

 
  facePlanL = extract_prism_face(cellPlanL, faceIDL);
  facePlanR = extract_prism_face(cellPlanR, faceIDR);

  if(facePlanL.size() == 0 && facePlanR.size() == 0) return facePlanL;
  
  Face* mergedFace = new Face();
  mergedFace->numEdge = 3;
  mergedFace->empty_resplit(facePlanL, orientCodeL);
  mergedFace->empty_resplit(facePlanR, orientCodeR);
  facePlanL.clear();
  facePlanL = mergedFace->make_faceplan();

  delete mergedFace;
  return facePlanL;
}
    
std::vector<char> merge_tri_face_p(const std::vector<char>& cellPlan,
                                   int faceID, char orientCode){
  
  
  std::vector<char> facePlan;
  
  if(faceID>=2 ) {
    cerr << "WARNING: illegal face ID in merge_tri_face_p" << endl;
    return facePlan;
  }
  
 
  facePlan = extract_prism_face(cellPlan,faceID);
  if(facePlan.size() == 0) return facePlan;
  
  Face* orientedFace = new Face();
  orientedFace->numEdge = 3;
  orientedFace->empty_resplit(facePlan, orientCode);
 
  facePlan.clear();
  facePlan =  orientedFace->make_faceplan();
  delete orientedFace;
  return facePlan;
}
