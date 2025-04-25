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

/////////////////////////////////////////////////////////////////////////////////////////////////////
//                                prism.cc
//                                by: Qiuhan Xue
//    This function include definition of class Prism, struct NeibIndex, IntegerPoint,
// IntegerFace, and related friend function setFaces()
//
/////////////////////////////////////////////////////////////////////////////////////////////////////
 
#include <stack>
#include <queue>
#include <vector>
#include <map>
#include <bitset>
#include "prism.h"
#include "globals.h"


//in each cell, only childCell, mySplitCode, and cellIndex is defined
//this function will not set mysplitcode unless the splitCode got from cellPlan
//is nonzero, and the the current cell from Q is never split
//
// this function will not set cellIndex unless the splitCode got from cellPlan is zero

int Prism::empty_resplit( const std::vector<char>& cellPlan){
  if(cellPlan.size() == 0){
    this->cellIndex = 1;
    return 1;
  }
  
  queue<Prism*> Q;
  Q.push(this);
  Prism* current;
  unsigned int index =0;
  int32 cIndex = 0;
  char currentCode=0;
  while(!Q.empty()){
    current = Q.front();
    
    if(index >= cellPlan.size()){
      currentCode = 0;
    }else{ 
      //take a code from splitcode
      currentCode = cellPlan[index];
      index++;  
    }
    
    if(currentCode ==0){
      current->cellIndex = ++cIndex;
    }else if((current->mySplitCode)!=0 && currentCode != (current->mySplitCode)){
      Loci::debugout << "WARNING: split code is not consistent" << endl;
      Loci::debugout << int(currentCode) << " intree " << int( current->mySplitCode)<< endl; 
      Loci::Abort();
      
    }else{
      if((current->childCell)==0){
        current->mySplitCode = currentCode;
        current->empty_split();
      }
      for(int i = 0; i <current->numChildren(); i++){
        Q.push(current->childCell[i]);
      }
    }
    
    Q.pop();
  }
  return cIndex;
}




void Prism::resplit( const std::vector<char>& cellPlan,
                     std::list<Node*>& node_list,
                     std::list<Edge*>& edge_list,
                     std::list<QuadFace*>& quadface_list,
                     std::list<Face*>& face_list,
                     std::vector<Prism*>& cells){
  if(cellPlan.size() == 0){
    this->cellIndex = 1;
    reduce_vector(cells);
    return;
  }
  
  queue<Prism*> Q;
  Q.push(this);
  Prism* current;
  unsigned int index =0;
  int32 cIndex = 0;
  char currentCode =0;
  
  while(!Q.empty()){
    current = Q.front();
    
    if(index >= cellPlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from splitcode
      currentCode = cellPlan[index];
      index++;  
    }

    if(currentCode ==0){
      current->cellIndex = ++cIndex;
      cells.push_back(current);
    }else if((current->mySplitCode)!=0 && currentCode != (current->mySplitCode)){
      Loci::debugout << "WARNING: split code is not consistent" << endl;
      Loci::debugout << int(currentCode) << " intree " << int( current->mySplitCode)<< endl; 
      Loci::Abort();
    }else{
      if((current->childCell)==0){
        current->mySplitCode = currentCode;
        current->split(node_list, edge_list, quadface_list, face_list);
      }
      for(int i = 0; i <current->numChildren(); i++){
        Q.push(current->childCell[i]);
      }
    }
    
    Q.pop();
  } 
}

//used in make_prism_cellplan.cc
void Prism::resplit( int level,
                     std::list<Node*>& node_list,
                     std::list<Edge*>& edge_list,
                     std::list<QuadFace*>& quadface_list,
                     std::list<Face*>& face_list){
 
  if(level <= 0) return;
  int currentLevel = level;
  
  queue<Prism*> Q;
  Q.push(this);
  Prism* current;
  
 
  while(!Q.empty()){
    current = Q.front();
    
    if(currentLevel > 0){
      current-> mySplitCode = 3;
      current->split(node_list, edge_list, quadface_list, face_list);
      
      for(int i = 0; i <current->numChildren(); i++){
        Q.push(current->childCell[i]);
      }
      currentLevel--;
    }
    else{
      current-> mySplitCode = 0;
    }
    Q.pop();
  } 
}




//This function check if aCell is my neighbor is direction dd
// edge neighbor is excluded. if two faces overlap area is not zero
// their cells are neighbors
bool  Prism::isNeighbor(const Prism* aCell, int dd, int& nf)const{
  //exclude myself
  if(aCell == this) return false;
  if(isSiblingNeighbor(aCell, dd, nf)) return true;//nf already modified here
  
  if(this->nfold == 3 || aCell->nfold == 3){
    
    switch(dd){
    case 0:
      nf = 1;
      return is_overlapped(gnrlface[0], aCell->gnrlface[1]);
      break;
     
    case 1:
      nf = 0;
      return is_overlapped(gnrlface[1], aCell->gnrlface[0]);
      break;
     
    case 2:
    case 3:
    case 4:
      return false;
      break;
    default:
      cerr<<"WARNING: illegal direction/face ID" << endl;
      break;
    }
  }
  
  else{//nfold ==4
    int j  =0;
    switch(dd){
    case 0:
      nf = 1;
      return is_overlapped(gnrlface[0], aCell->gnrlface[1]);
      break;
       
    case 1:
      nf = 0;
      return is_overlapped(gnrlface[1], aCell->gnrlface[0]);
      break;
     
    case 2:
    case 3:
    case 4:
    case 5:
     
      for( j = 2; j <=5; j++){
        if(is_overlapped(quadface[dd-2], aCell->quadface[j-2])){
          nf = j ;
          return true;
        }
      }
      if(j == 6) {
        nf = -1;
        return false;
      }
      break;
  
    default:
      cerr<<"WARNING: illegal direction/face ID" << endl;
      break;
    }
  }
  return false;
}  

 

//this function check if aCell is my sibling neighbor.
//sibling means same size face, not necessarily has same parent
bool  Prism::isSiblingNeighbor(const Prism* aCell, int dd, int &nf)const{
  if(aCell == this) return false;
  if(this->nfold == 3 || aCell->nfold ==3){
    switch(dd){
    case 0:
      nf = 1;
      return (gnrlface[0] == aCell->gnrlface[1]);
      break;
    case 1:
      nf = 0;
      return (gnrlface[1] == aCell->gnrlface[0]);
      break;
    case 2:
    case 3:
    case 4:
      nf = -1;
      return false;
      break;
    default:
      cerr<<"WARNING: illegal direction/face ID" << endl;
      break;
    }
  }
  
  else if(this->nfold == 4){
    switch(dd){
    case 0:
      nf = 1;
      return (gnrlface[0] == aCell->gnrlface[1]);
      break;
    case 1:
      nf = 0;
      return (gnrlface[1] == aCell->gnrlface[0]);
      break;

    case 2:
    case 3:
    case 4:
    case 5:
      
      for(int j = 2; j <=5; j++){
        if(quadface[dd-2] == aCell->quadface[j-2]){
          nf = j ;
          return true;
        }
      }
      nf = -1;
      return false;
      break;
    default:
      cerr<<"WARNING: illegal direction/face ID" << endl;
      break;
    }
  }
  cerr<<"WARNING: reach dummy code in Prism::isSiblingNeib()" << endl;
  return false;
}

Prism* Prism::getSiblingNeib(int dd, int& nf){
  
  nf = -1;
  if(parentCell == 0){
    return 0; //root has no neib
  }
  
  
  int intCode = parentCell->mySplitCode;
  int childID = whichChild();
  int pnfold = parentCell->nfold;
 
  //if sibling neib exist, return sibling neib
  if( intCode == 1){
    if(dd == 0 && childID == 1){
      nf = 1;
      return parentCell->childCell[0];
    }
    if(dd == 1 && childID == 0){
      nf = 0;
      return parentCell->childCell[1];
    }
  }
  else if(intCode == 2){
    //find sibling neighbor(real sibling)
    if(dd == 3){
      nf = 4;
      return parentCell->childCell[childID == (pnfold-1)?0:(childID+1)];
    }
    else if(dd == 4){
      nf = 3;
      return parentCell->childCell[childID ==0?(pnfold -1): (childID-1)];
    }
  }
  else if(intCode == 3){
    if(dd == 0 && childID >= pnfold){
      nf = 1;
      return parentCell->childCell[childID - pnfold];
    }
    if(dd == 1 && childID < pnfold){
      nf = 0;
      return parentCell->childCell[childID + pnfold];
    }
    if(dd == 3){
      nf = 4;
      if(childID < pnfold) return parentCell->childCell[childID == (pnfold-1)?0:(childID+1)];
      else return parentCell->childCell[childID == (2*pnfold-1)?pnfold:(childID+1)];
    }
    else if(dd == 4){
      nf = 3;
      if(childID < pnfold)return parentCell->childCell[childID ==0?(pnfold -1): (childID-1)];
      else return parentCell->childCell[childID ==pnfold?(2*pnfold -1): (childID-1)];
    }   
    
  }
  return 0;
}





//this function find neighbor in direction dd, the neighbor's face is greater or
// same as mine
Prism* Prism::findNeighbor(int dd, int& nf)
{

  if(parentCell==0) return 0; //root has no neib;

  Prism* tempNeighbor = getSiblingNeib(dd, nf);
  if(tempNeighbor != 0) return tempNeighbor;

  //otherwise find the nearest common ancestor 
 
  if(dd < 2){
    tempNeighbor = parentCell->findNeighbor(dd, nf);
    if(tempNeighbor == 0) return 0;  //no neib
    
    //then climb down to find my neighbor
    while(tempNeighbor != 0){
          
      //if tempNeighbor is a leaf and is my sibling neighbor, return it
      int tempNf = -1;
      if(isSiblingNeighbor(tempNeighbor, dd, tempNf)&& tempNeighbor->mySplitCode == 0 ) return tempNeighbor;
      
      
      //if tempNeighbor's face has been smaller than mine, return tempneighbor's parentCell 
      
      if(tempNeighbor->getLevel(0) > this->getLevel(0))return tempNeighbor->parentCell;
      
      //if tempNeighbor is a leaf(whose face is greater than mine) , return
      if(tempNeighbor->mySplitCode == 0)return tempNeighbor; 
      
      //if tempNeighbor is not a leaf, and her face still greater than mine, climb down
      //caution: here if no child is neighbor, a dead loop will happen
      int i;
      int nc = tempNeighbor->numChildren(); 
      for(i = 0; i < nc; i++){
        if( isNeighbor(tempNeighbor->childCell[i], dd, nf)){
        
          tempNeighbor = tempNeighbor->childCell[i];
          break;
          //only need find one neighbor child.Because if there are two or more neighbors,
          //they must be smaller, so the switch() statement will return the parent cell
        }
      }
      if( i == nc){
        cerr<<"WARNING: get lost when climbing down the tree" << endl;
        exit(0);
      }
    }//while
  }
  else if(dd >=2){//this has no sibling neib
    //problem: faceStack can not be used here because can get lost when climb down,
    //if faceStack.top() is overlap with two
    //childCells, and only select the first one found, when faceStack.pop(), it's found that the wrong one
    //is selected;

    //if the stack take one step and pNeighbor 
    // take two steps, the faceStack will be empty before the real neighbor
    // is found
    
    //solution: don't need a faceStack, when climb down, overlap with this->quadface[dd-2]    

    //first find NCA
    Prism* current = this;
    int   currentFace = dd;
    Prism* NCA;
    
    while(true){
      
      if(current->parentCell == 0) return 0; 
      if((currentFace==3 || currentFace == 4) && current->parentCell->mySplitCode != 1){
        NCA = current->getSiblingNeib(currentFace, nf);
        break;
      }
           
      int newFaceID = current->parentFace(currentFace);
      if(newFaceID < 0) {
        cerr << "WARNING: illegal faceID" << endl;
        Loci::Abort();
      }
      //climb up
      current = current->parentCell;
      currentFace = newFaceID;
    }//while(true)

    
    if(NCA == 0){
      return 0;  //no neib
    }
    
    Prism* pNeighbor = NCA;
  
    
    //if pNeighbor is a leaf, return it
    if(pNeighbor->childCell == 0) return pNeighbor;
    
    
    //if pNeighbor is not a leaf, and her face still greater than mine,
    //climb down to my level
    bool childFound;
    int cID = 0;
    //if allow ==, next step can be two small
    while( pNeighbor->getLevel(0) <= getLevel(0) &&
           pNeighbor->getLevel(1) <= getLevel(1) &&
           pNeighbor->childCell!= 0){

      if(pNeighbor->getLevel(0) == getLevel(0) && pNeighbor->getLevel(1)==getLevel(1)) return pNeighbor;
      //some check should be done here
     
      
      childFound = false;
      for(cID = 0 ; cID < pNeighbor->numChildren(); cID++){
        for(int j = 0; j <(pNeighbor->childCell[cID]->nfold); j++){
          if(is_overlapped(this->quadface[dd-2], pNeighbor->childCell[cID]->quadface[j])){
            nf = j+2;
            childFound = true;
            break;
          }
        }
        if(childFound) break;
        
      }
      if(!childFound){
        cerr <<"WARNING: get lost when climbing down the tree" << endl;
        Loci::Abort();
      } 
      
      pNeighbor = pNeighbor->childCell[cID];

      
      if(pNeighbor->getLevel(0) > getLevel(0) || pNeighbor->getLevel(1) > getLevel(1)) return pNeighbor->parentCell;
      
    }
    
    //sanity check
    if(!is_overlapped(this->quadface[dd-2], pNeighbor->quadface[nf-2])){
      cerr << " WARNING: get the wrong neighbor" << endl;
      Loci::Abort();
    }

    return pNeighbor;
  }//else if(dd >=2)
  return 0;
}



// this function check each cell in cells, find all the individual faces,put their connectivity info
// into faces. 

// cells: in, all the leaves cells in the 2N tree
// faces: in and out, initially empty, return all the indivial faces and its two cell inex

void set_prism_faces(const std::vector<Prism*>& cells,
                     std::map<QuadFace*, NeibIndex>& qfaces,
                     std::map<Face*, NeibIndex>& gfaces){
  
  
  Prism* tempNeib;
  int numCells = cells.size();
  
  if(numCells > 0){
    int nf = -1;   
    //for each cell
    for(int i = 0; i < numCells; i++){
      int numFaces = (cells[i]->nfold == 3?5:6);
      for(int dd = 0; dd <numFaces; dd++){
        if(!((cells[i]->faceMarked).test(dd))){
          cells[i]->faceMarked.set(dd);
          tempNeib = cells[i]->findNeighbor(dd, nf);
          //if no neib, don't output the face
          //there is a neib
          if(tempNeib != 0){
            int tempnf = -1;
            bool isSibling = cells[i]->isSiblingNeighbor(tempNeib,dd, tempnf);
            bool isALeaf =  (tempNeib->mySplitCode == 0);
                     
            //if neib is a leaf, marked its face
            if(isALeaf){
              // tempNeib ->faceMarked.set(nf);
             
              if(dd < 2){
                if(dd == 1) gfaces[cells[i]->gnrlface[dd]] =  NeibIndex(cells[i]->cellIndex, tempNeib ->cellIndex);
                else gfaces[cells[i]->gnrlface[dd]] =  NeibIndex( tempNeib ->cellIndex,cells[i]->cellIndex);
                tempNeib ->faceMarked.set(nf);
              }
              else{
                std::vector<QuadFace*> commonfaces = overlap(cells[i]->quadface[dd-2], tempNeib->quadface[nf-2]);
                if(commonfaces.size() != 1){
                  cerr << "WARNING: more than one commom face" << endl;
                }
                //faceOrient == 0, point outward
                if(!((cells[i]->faceOrient).test(dd-2))) qfaces[commonfaces[0]] = NeibIndex(cells[i]->cellIndex, tempNeib->cellIndex);
                else qfaces[commonfaces[0]] = NeibIndex( tempNeib->cellIndex,cells[i]->cellIndex);
              }
            }
            
            else if(!isSibling) {
              //this stack hold all the neighbor leaves from  tempNeib
              stack<pair<Prism*, int> > cellStack;
              Prism* current=0;
              cellStack.push(make_pair(tempNeib, nf));
              int currentNf;
              while(!cellStack.empty()){
                //current takes out a cell from stack
                current = cellStack.top().first;
                currentNf = cellStack.top().second;
                cellStack.pop();
                
                
                int  numNeibChildren = current->numChildren();
                if( numNeibChildren != 0){
                  for(int j = 0; j < numNeibChildren; j++){
                    if(cells[i]->isNeighbor(current->childCell[j], dd, nf)){
                      cellStack.push(make_pair(current->childCell[j], nf));
                    }
                  }
                }//if not leaf, put it children into stack
              
               
                
                //proocess leaves
                else if(dd >= 2){
                  //find the intersection of cells[i]->face[dd] and current->face[nf] leaves
                  std::vector<QuadFace*> commonfaces = overlap(cells[i]->quadface[dd-2], current->quadface[currentNf-2]);
                  
                  if(commonfaces.size() != 1){
                    cerr << "WARNING: more than one commom face" << endl;
                  }
                  //faceOrient == 0 point outward
                  if(!((cells[i]->faceOrient).test(dd-2)))  qfaces[commonfaces[0]] = NeibIndex(cells[i]->cellIndex, current->cellIndex);
                  else qfaces[commonfaces[0]] = NeibIndex( current->cellIndex,cells[i]->cellIndex);
                }
                else{//dd<2, leaves
                  //dd == 1 point outward                                                      
                  if(dd == 1) gfaces[current->gnrlface[currentNf]] = NeibIndex(cells[i]->cellIndex, current->cellIndex);
                  else  gfaces[current->gnrlface[currentNf]] = NeibIndex(current->cellIndex,cells[i]->cellIndex);
                  current->faceMarked.set(currentNf);
                }
                
               
              }//while
            }//elseif(!isSibling)
            
          }//if(tempNeib !=0)
        }//if(!faceMarked.test[dd])
      }//for(int dd = 0; dd < numFaces; dd++) 
    }//for(int i = 0; i < numcells; i++)
  }//if(numCells > 0)
}


void Prism::empty_split(){
  if(childCell!=0)return;
  switch(mySplitCode)
    {
      //000 no split,this is a leaf, compare and set the maxLevels, output cells
    case 0:
      break;
      //100  x direction is splitted 
    case 1:
      childCell = new Prism*[2];
      for(int i = 0; i <2; i++){
        childCell[i] = new Prism();
        childCell[i] ->nfold = nfold;
      }
      break;
      
      //010  y direction is splitted 
    case 2:
      childCell = new Prism*[nfold];
      for(int i = 0; i < nfold; i++){
        childCell[i] = new Prism();
        childCell[i]->nfold = 4;
      }
      break;
      //001  z direction is splitted 
    case 3:
      childCell = new Prism*[2*nfold];
      for(int i = 0; i <2*nfold; i++){
        childCell[i] = new Prism();
        childCell[i]->nfold = 4;
      }
      break;
    default:
      cerr <<"WARNING: illegal splitcode in function Prism::reSplit()" << endl;
      break;
    }
}



  
void Prism::split( std::list<Node*>& node_list,
                   std::list<Edge*>& edge_list,
                   std::list<QuadFace*>& quadface_list,
                   std::list<Face*>& face_list){
  if(childCell!=0)return;
  Face* newFace = 0;
  QuadFace* newface[8];
  Face* newgface[4];
  Edge* newedge[6];
  Node* ccenter = 0;
  
  switch(mySplitCode)
    {
      
      //000
    case 0:
      break;
      
      
    case 1:
      
      childCell = new Prism*[2];
      for(int i = 0; i <2; i++){
        childCell[i] = new Prism(nfold);
        childCell[i]->parentCell = this;
      }
      
      //quadface split in y direction
      for(int i = 0; i < nfold; i++){
        quadface[i]->split(char(1),char(0), node_list, edge_list);
      }
      
      
      //define new face
      newFace = new Face(nfold);
      face_list.push_back(newFace);
      
      for(int i = 0; i < nfold; i++){
        newFace->edge[i] = quadface[i]->childy[1]->edge[0];
        newFace->needReverse[i] = faceOrient.test(i); 
      }

      //define childCell
      childCell[0]->gnrlface[1]= childCell[1]->gnrlface[0] = newFace;
      childCell[0]->gnrlface[0] = gnrlface[0];
      childCell[1]->gnrlface[1] = gnrlface[1];
      for(int i = 0; i < nfold; i++){
        childCell[0]->quadface[i] = quadface[i]->childy[0];
        childCell[1]->quadface[i] = quadface[i]->childy[1];
      }
      
      childCell[0]->faceOrient = childCell[1]->faceOrient = faceOrient;
      break;
      
      
      //010  y direction is splitted 
    case 2:
      
      childCell = new Prism*[nfold];
     
      for(int i = 0; i < nfold; i++){
        childCell[i] = new Prism(4);
        //childCell[i]->whichChild = i;
        childCell[i]->parentCell = this;
      }
     
      //quadface split in x direction
      for(int i = 0; i < nfold; i++){
        quadface[i]->split(char(2),char(0), node_list, edge_list);
      }
      //gnrlface split
      for(int i = 0; i < 2; i++){
        if(gnrlface[i]->child == 0) gnrlface[i]->split( node_list, edge_list);
      } 

      //define new edge
      newedge[0] = new Edge(getFaceCenter(0), getFaceCenter(1), quadface[0]->edge[1]->level);
      edge_list.push_back(newedge[0]);
     
     
     
      //define new face
      for(int i = 0; i < nfold; i++){
        newface[i] = new QuadFace(4);
        quadface_list.push_back(newface[i]);
      }   
     
  
 
     
      //defines new face
      for(int i = 0; i < nfold; i++){
        newface[i]-> edge[0] = gnrlface[0]->child[i]->edge[1];
        newface[i]->edge[1] = quadface[i]->childx[0]->edge[1];
        newface[i]->edge[2] = gnrlface[1]->child[i]->edge[1];
        newface[i]->edge[3] = newedge[0];
        //default; no edge need reversed
      }
       
       
      //define childCell 
      for(int i = 0; i < nfold; i++){
        childCell[i]->gnrlface[0] = gnrlface[0]->child[i];
        childCell[i]->gnrlface[1] = gnrlface[1]->child[i];
        childCell[i]->quadface[0] = quadface[i]->childx[faceOrient.test(i)?1:0];
        childCell[i]->faceOrient[0] = faceOrient[i];
        childCell[i]->quadface[1] = newface[i];
        childCell[i]->faceOrient[1]= 1;
        childCell[i]->quadface[2] = newface[i==0? (nfold-1):(i-1)];
        childCell[i]->quadface[3] = quadface[i==0?(nfold-1):(i-1)]->childx[faceOrient.test(i==0?(nfold-1):(i-1))?0:1];
        childCell[i]->faceOrient[3] = faceOrient[i==0?(nfold-1):(i-1)];
      }
             
      break;
             
    case 3:
      childCell = new Prism*[2*nfold];
      
      for(int i = 0; i < 2*nfold; i++){
        childCell[i] = new Prism(4);
        //childCell[i]->whichChild = i;
        childCell[i]->parentCell = this;
      }
                  
      //quadface split in x and y direction
      for(int i = 0; i < nfold; i++){
        quadface[i]->split(char(3),char(0), node_list, edge_list);
      }
      //gnrlface split
      for(int i = 0; i < 2; i++){
        if(gnrlface[i]->child ==0) gnrlface[i]->split( node_list, edge_list);
      }

      
      //compute cell center
      ccenter = centroid();
      node_list.push_back(ccenter);
        
      //define newedge[0] from facecenter[0] to cellcenter
      newedge[0] = new Edge(getFaceCenter(0), ccenter, getLevel(1)+1);
      edge_list.push_back(newedge[0]);
        
      //define newedge[1] from cellcenter to facecenter[1] 
      newedge[1] = new Edge(ccenter, getFaceCenter(1), getLevel(1)+1);
      edge_list.push_back(newedge[1]);

      //define other new edges from cellcenter to quadface centers
      for(int i = 0; i< nfold; i++){
        newedge[i+2] = new Edge(ccenter, quadface[i]->getCenter(), getLevel(0)+1);
        edge_list.push_back(newedge[i+2]);
      }
        
      //define new face
      for(int i = 0; i < 2*nfold; i++){
        newface[i] = new QuadFace(4);
        quadface_list.push_back(newface[i]);
      }
      
      for(int i = 0; i < nfold; i++){
        newgface[i] = new Face(4);
        face_list.push_back(newgface[i]);
      }
      
  
     
        
      //build new quad face
      for(int i = 0; i < nfold; i++){
        newgface[i]->edge[0] = quadface[i]->child[faceOrient.test(i)?2:0]->edge[2];
        newgface[i]->needReverse[0] = faceOrient.test(i);
        newgface[i]->edge[1] = newedge[i+2];
        newgface[i]->needReverse[1] = true;
        newgface[i]->edge[2] = newedge[(i==0)?(nfold-1+2):(i-1+2)];
        newgface[i]->needReverse[2] = false;
        newgface[i]->edge[3] = quadface[i==0?(nfold-1):(i-1)]->child[faceOrient.test(i==0?(nfold-1):(i-1))?0:2]->edge[2];
        newgface[i]->needReverse[3] = faceOrient.test(i==0?(nfold-1):(i-1));
      }
        
      //build new quadface
      for(int i = 0; i < nfold; i++){
        newface[i]->edge[0] = gnrlface[0]->child[i]->edge[1]; 
        newface[i]->edge[1] = quadface[i]->child[0]->edge[1]; 
        newface[i]->edge[2] = newedge[i+2]; 
        newface[i]->edge[3] = newedge[0];
        //no edge need reverse
      }
      for(int i = 0; i < nfold; i++){
        newface[i+nfold]->edge[0] = newedge[i+2];
        newface[i+nfold]->edge[1] = quadface[i]->child[1]->edge[1];
        newface[i+nfold]->edge[2] = gnrlface[1]->child[i]->edge[1];
        newface[i+nfold]->edge[3] = newedge[1]; 
        //no edge need reverse
      }

       
    
      //build childcell
      for(int i = 0; i < nfold; i++){
        childCell[i]->gnrlface[0] = gnrlface[0]->child[i];
        childCell[i]->gnrlface[1] = newgface[i];
        childCell[i]->quadface[0] = quadface[i]->child[faceOrient.test(i)?2:0];
        childCell[i]->faceOrient[0] = faceOrient[i];
        childCell[i]->quadface[1] = newface[i];
        childCell[i]->faceOrient[1] = 1;
        childCell[i]->quadface[2] = newface[i==0? (nfold-1):(i-1)];
        childCell[i]->quadface[3] = quadface[i==0?(nfold-1):(i-1)]->child[faceOrient.test(i==0?(nfold-1):(i-1))?0:2];
        childCell[i]->faceOrient[3] = faceOrient[i==0?(nfold-1):(i-1)];
      }
      
      for(int i = 0; i < nfold; i++){
        childCell[i+nfold]->gnrlface[0] = newgface[i];
        childCell[i+nfold]->gnrlface[1] = gnrlface[1]->child[i];
        childCell[i+nfold]->quadface[0] = quadface[i]->child[faceOrient.test(i)?3:1];
        childCell[i+nfold]->faceOrient[0] = faceOrient[i];
        childCell[i+nfold]->quadface[1] = newface[i+nfold];
        childCell[i+nfold]->faceOrient[1] = 1;
        childCell[i+nfold]->quadface[2] = newface[i==0? (2*nfold-1):(i+nfold-1)];
        childCell[i+nfold]->quadface[3] = quadface[i==0?(nfold-1):(i-1)]->child[faceOrient.test(i==0?(nfold-1):(i-1))?1:3];
        childCell[i+nfold]->faceOrient[3] = faceOrient[i==0?(nfold-1):(i-1)];
      }
      
      break;
      
    default:
      cerr <<"WARNING: illegal splitcode in function split()" << endl;
      break;
    }
}

int Prism::get_num_fine_faces(){
  int num_faces = 0;
  
  for(int i = 0; i < 2; i++){
    num_faces += gnrlface[i]->get_num_leaves();
  }
  
  for(int i = 0; i < nfold; i++){
    num_faces += quadface[i]->get_num_leaves();
  }
  return num_faces;
}


int Prism::get_tagged(){
  std::vector<Node*> nodes(2*nfold);
  get_nodes(nodes); 
    
  //if all nodes get detagged, the cell is detagged
  bool detagged = true;
  for(std::vector<Node*>::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    if((*np)->tag != 2)detagged = false;
  }
  if(detagged) return 2;
  for(std::vector<Node*>::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    if((*np)->tag == 1)return 1;
    
  }
  //otherwise, the cell remains unchanged
  return 0;
} 
  

int Prism::get_tagged(const vector<source_par>& sources){
  std::vector<Node*> nodes(2*nfold);
  get_nodes(nodes); 
  
  double min_len = get_min_edge_length();
  if(tag_cell(nodes, sources, min_len)){
    mySplitCode = 3;
    return 1;
  }else{
    mySplitCode = 0;
    return 0;
  }
  return 0;
}
//find the minimum edge length in a cell(before split)
double Prism::get_min_edge_length(){
  std::vector<Edge*> edges = get_edges();
  std::vector<double> edge_length(3*nfold);
  for(int i = 0; i < 3*nfold; i++)edge_length[i] = edges[i]->get_length();
    
  double min_length = edge_length[0];
  for(int i = 1; i < 3*nfold; i++){
    min_length = min(min_length, edge_length[i]);
  }
  return min_length;
}


void Prism::setSplitCode(int split_mode, double tol){
  if(mySplitCode != 0){
    cerr<<"WARNING:mySplitCode is not zero in setSplitCode of prism: " << int(mySplitCode) <<endl;
    exit(0);
  }
  //find the edge_length of all edges
  std::vector<Edge*> edges = get_edges();
  std::vector<double> edge_length(3*nfold);
  for(int i = 0; i < 3*nfold; i++)edge_length[i] = edges[i]->get_length();
    
  //find the min_edge_length in XY direction
  double average_length0, average_length1;
    
  average_length0  = 0.0;
  for(int i = 0; i<2*nfold; i++) average_length0 +=  edge_length[i];
  average_length0 /= double(2*nfold);
    
  average_length1 = 0.0;
  for(int i = 2*nfold; i < 3*nfold; i++)average_length1 +=  edge_length[i];
  average_length1 /= double(nfold);
    
  bitset<2> tolerance_mask(3); //all 1s
  if(average_length0 < 2*tol)tolerance_mask.reset(1);
  if(average_length1 < 2*tol)tolerance_mask.reset(0);
   
  double minimum_length = min(average_length0, average_length1);
  if(split_mode == 0){
    
    if(average_length1/average_length0 > Globals::factor){
      bitset<2> oldCode(1);
      oldCode = oldCode & tolerance_mask;
      mySplitCode = char(oldCode.to_ulong());
      //mySplitCode = 1; //split in z direction
      return;
    }
    else if(average_length0/average_length1 > Globals::factor){
      bitset<2> oldCode(2);
      oldCode = oldCode & tolerance_mask;
      mySplitCode = char(oldCode.to_ulong());
      //mySplitCode = 2;//split in xy direction
      return;
    }
    else if(minimum_length > 2.0*tol){
      mySplitCode = 3;
      return;
    }
    else{
      mySplitCode = 0;
      return;  
    }
      
    cerr<< " WARNING: reach dummy code in setSplitCode()" << endl;
  }
  else if(split_mode == 1){
    if(average_length0/tol > 2.0){
      // bitset<2> oldCode(2);
      // oldCode = oldCode & tolerance_mask;
      // mySplitCode = char(oldCode.to_ulong());
        
      mySplitCode = 2;//split in xy direction
      return;
    }
    else{
      mySplitCode = 0;
      return;  
    }  
  }
  else if(split_mode == 2){
    if (minimum_length > 2.0*tol){
      mySplitCode = 3;
      return;
    }
    else{
      mySplitCode = 0;
      return;
    }
  }
}
  

//after a cell is split, compose the cell plan according to the tree structure      
std::vector<char> Prism::make_cellplan(){
  
  std::vector<char> cellPlan;
  
  if(childCell == 0) {
    reduce_vector(cellPlan);
    return cellPlan;
  }
 
  std::queue<Prism*> Q;
  Q.push(this);
  
  Prism* current;
  while(!Q.empty()){
    current = Q.front();
    if(current->getMySplitCode() != 0){
      cellPlan.push_back(current->getMySplitCode());
      for(int i = 0; i < current->numChildren(); i++){
        Q.push(current->getChildCell(i));
      }
    }
    else{
      cellPlan.push_back(0);
    }
    Q.pop();
  }                 
  while(cellPlan.size() != 0 && cellPlan.back() == 0) cellPlan.pop_back();
  reduce_vector(cellPlan);
  return cellPlan;
}

std::vector<char> Prism::make_cellplan(int level){
  
  std::vector<char> cellPlan;
 
  
  if(level <= 0) {
    reduce_vector(cellPlan);
    return cellPlan;
  }

  cellPlan.push_back(3);
  for(int i = 1; i < level; i++){
    int total_num_children = 6*(1 << (3*(i -1)));
    for(int j = 0; j < total_num_children; j++) cellPlan.push_back(3);
  }
  reduce_vector(cellPlan);
  return cellPlan;
}


bool Prism::balance_cell(int split_mode,
                         std::list<Node*>& node_list,
                         std::list<Edge*>& edge_list,
                         std::list<QuadFace*>& qface_list,
                         std::list<Face*>& gface_list){ 

  
  
  bool  needBalance = false;
  if(childCell!=0){
    std::list<Prism*> leaves;
    sort_leaves(leaves);
    needBalance = false;
    for(std::list<Prism*>::const_iterator p = leaves.begin(); p != leaves.end(); p++){
      bool tmp =  (*p)->balance_cell(split_mode, node_list, edge_list, qface_list, gface_list);
      needBalance = tmp||needBalance;
    }
    return needBalance;
  }
  
  else if(childCell == 0){
    std::vector<Edge*> edge= get_edges();

    needBalance = false;
   

    if(split_mode == 2){

      for(int i = 0; i < 3*nfold; i++){
        if(edge[i]->depth_greater_than_1()){
          mySplitCode = 3;
          break;
        }
      }

      if(Globals::balance_option >= 1 && mySplitCode == 0){
        if((gnrlface[0]->child != 0) && (gnrlface[1]->child != 0)){
          mySplitCode = 3;
        }else{
          int num_faces_split = 0;
          for(int i = 0; i < nfold; i++){
            if(quadface[i]-> code != 0 ) num_faces_split++;
          }
          if(num_faces_split > (nfold/2)){
            
            mySplitCode = 3;
          }
        }
      }

      if(Globals::balance_option >= 2 && mySplitCode == 0){
        if( nfold ==4){
          if(((quadface[0]-> code != 0) && (quadface[2]->code != 0))||
             ((quadface[1]-> code != 0) && (quadface[3]->code != 0))
             )mySplitCode = 3;
        }
        
      }
      
    }else{
      bitset<2> code(0);
      
      for(int i = 0; i < 2*nfold; i++){
        if( edge[i]->depth_greater_than_1()){
          code.set(1);
          break;
        }
      }
      for(int i = 2*nfold; i < 3*nfold; i++){
        if( edge[i]->depth_greater_than_1()){
          code.set(0);
          break;
        }
      }
    
      mySplitCode = char(code.to_ulong());

      
      if(Globals::balance_option >= 1 && mySplitCode == 0){
        int num_faces_split = 0;
        bitset<2> face_code;
        if(gnrlface[0]->child != 0 && gnrlface[1]->child!=0)face_code.set(1);
         
        for(int i = 0; i < nfold; i++){
          if(quadface[i]-> code != 0 ) num_faces_split++;
        }
         
        if(num_faces_split > (nfold/2)){
          for(int i = 0; i< nfold; i++) face_code |= bitset<2>(quadface[i]->code);
          
        }
        mySplitCode = char(face_code.to_ulong());
      }

      
      if(Globals::balance_option >= 2 && mySplitCode == 0){
        bitset<2> the_code;
        the_code.reset();
        if((gnrlface[0]->child != 0) && (gnrlface[1]->child != 0)){
          the_code.set(1);
        }
        if( nfold ==4){
          if(((quadface[0]-> code != 0) && (quadface[2]->code != 0))||
             ((quadface[1]-> code != 0) && (quadface[3]->code != 0))
             ){
            
            for(int i = 0; i< nfold; i++) the_code |= bitset<2>(quadface[i]->code);
            
          }
        }
          
        mySplitCode = char(the_code.to_ulong());  
      }
      
    }
    needBalance = (mySplitCode != 0);
    if(needBalance) {
      split(node_list, edge_list, qface_list, gface_list);
    
      for(int i = 0; i < numChildren(); i++){
        childCell[i]->balance_cell(split_mode,node_list, edge_list, qface_list,gface_list );
      }
    }
  }
  
 
  return needBalance;
}



void Prism::sort_leaves(std::list<Prism*>& leaves){
  if(childCell != 0){
    for(int i = 0; i < numChildren(); i++)childCell[i]->sort_leaves(leaves);
  }
  else{
    leaves.push_back(this);
  }
}


void Prism::rebalance_cells(int split_mode,
                            std::list<Node*>& node_list,
                            std::list<Edge*>& edge_list,
                            std::list<QuadFace*>& qface_list,
                            std::list<Face*>& gface_list){

  bool need_balance_more = true;

  while(need_balance_more){
    need_balance_more = balance_cell(split_mode, node_list, edge_list, qface_list, gface_list); 
  }
}

int32 Prism::traverse(const std::vector<char>& parentPlan,  vector<pair<int32, int32> >& indexMap){
  indexMap.clear();
  if(parentPlan.size() == 0){
    if(numChildren()!=0){
      list<Prism*> leaves;
      sort_leaves(leaves); 
      for(std::list<Prism*>::const_iterator p = leaves.begin(); p != leaves.end(); p++)
        indexMap.push_back(make_pair((*p)->cellIndex, 1));
      return 1;
    }else{
      indexMap.push_back(make_pair(1,1));
      return 1;
    }
  }
  std::queue<Prism*> Q;
  Q.push(this);
  Prism* current;
  unsigned int index =0;
  int32 cIndex = 0;
  char currentCode;
  while(!Q.empty()){
    current = Q.front();
    if(index >= parentPlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from splitcode
      currentCode = parentPlan[index];
      index++;  
    }
    list<Prism*> leaves;
    switch(currentCode){
      //0 no split,this is a leaf,this is a leaf for parentPlan
    case 0:
      ++cIndex;
      current->sort_leaves(leaves);
      if(leaves.front() != current){//current is not a leaf for cellPlan
        for(std::list<Prism*>::const_iterator p = leaves.begin(); p != leaves.end(); p++){
          indexMap.push_back(make_pair((*p)->cellIndex, cIndex));
        }
      }else{
        //current is a leaf for cellPlan and for parentPlan
        if(current->cellIndex != 0) indexMap.push_back(make_pair(current->cellIndex, cIndex));  
      }
      break;
    default:
      if(current->cellIndex !=0){//derefinement happen, current is a leaf for cellPlan
        current->sort_leaves(leaves);
        for(std::list<Prism*>::const_iterator p = leaves.begin(); p != leaves.end(); p++){
          indexMap.push_back(make_pair(current->cellIndex, (*p)->cellIndex));
          (*p)->cellIndex = 0; //so that *p will not appear in indexMap when it is popped out of Q
        }
      }
      for(int i = 0; i <current->numChildren(); i++){
        Q.push(current->childCell[i]);
      }
    }
    
    Q.pop();
  }
  return cIndex;
}

bool Prism::needDerefine(){
  if(childCell != 0){
    bool derefine = true;
    for(int i = 0; i < numChildren(); i++){
      if( (childCell[i] ->get_tagged()) != 2) return false;
      if (childCell[i]->childCell !=0)return false;
    }
    if(derefine){
      std::vector<Edge*> edge= get_edges();
      for(unsigned int i = 0; i < edge.size(); i++){
	if(edge[i]->depth_greater_than_1())return false;
      }
      return true;
    }
  }
  return false;
}
bool Prism::needDerefine_ctag(){
  if(childCell != 0){
    bool derefine = true;
    for(int i = 0; i < numChildren(); i++){
      if( (childCell[i] ->getTag()) != 2) return false;
      if (childCell[i]->childCell !=0)return false;
    }
    if(derefine){
      std::vector<Edge*> edge= get_edges();
      for(unsigned int i = 0; i < edge.size(); i++){
	if(edge[i]->depth_greater_than_1())return false;
      }
      return true;
    }
  }
  return false;
}


void Prism::derefine(){
  if(childCell != 0){
    int nc = numChildren();
    for(int i = 0; i < nc; i++){
      if(childCell[i] != 0){
	delete  childCell[i];
	childCell[i] = 0;
      }
    }
    delete [] childCell;
    childCell = 0;
    mySplitCode = 0;
  }
}

