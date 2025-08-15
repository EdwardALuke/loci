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
//                                hexcell.cc
//                                by: Qiuhan Xue
//    This function include definition of class HexCell
// 
//
/////////////////////////////////////////////////////////////////////////////////////////////////////
 
#include <stack>
#include <queue>
#include <vector>
#include <map>
#include <algorithm>
#include "hexcell.h"
#include "globals.h"
using std::stack;
using std::queue;
using std::cerr;
using std::endl;
using std::cout;


//splitCode: in, the refinement plan

//maxLevels: in and out, initial value is 0, return value the max levels in each dimension
// of the whole tree. will be used to set the integer coordinates


//cellID: in and out, initial value 0, each time a leaf is met, its value is sent to the leaf's
//cellIndex, and cellID++

//cells: in and out, return all the leaves cells in a std::vector


void HexCell::resplit( const std::vector<char>& cellPlan,
                       std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list,
                       std::list<QuadFace*>& face_list,
                       std::vector<HexCell*>& cells){


  
  if(cellPlan.size() == 0){
    this->cellIndex = 1;
    reduce_vector(cells);
    return;
  }
  
  queue<HexCell*> Q;
  Q.push(this);
  HexCell* current;
  unsigned int index =0;
  int32 cIndex = 0;
  char currentCode;
  
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
      cells.push_back(current);
    }else if( (current->mySplitCode != 0)  && (currentCode != current->mySplitCode)){//consistency check
      Loci::debugout << " nonconsistent split code in hexcell resplit" << endl;
      Loci::Abort();
    }else{
      current->mySplitCode = currentCode;
      current->split(node_list, edge_list, face_list);
      for(int i = 0; i <current->numChildren(); i++){
        Q.push(current->childCell[i]);
      }
    }
    
    Q.pop();
  }

}

void HexCell::resplit(int level,
                      std::list<Node*>& node_list,
                      std::list<Edge*>& edge_list,
                      std::list<QuadFace*>& face_list){
  if(level <= 0) return;
  int currentLevel = level;
  queue<HexCell*> Q;
  Q.push(this);
  HexCell* current;
  
  while(!Q.empty()){
    current = Q.front();
    
    if(currentLevel > 0){
      current-> mySplitCode = 7;
      current->split(node_list, edge_list, face_list);
      
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



//this function will return num_fine_cells
//this function will not set mysplitcode unless the splitCode got from cellPlan
//is nonzero, and the the current cell from Q is never split
//
// this function will not set cellIndex unless the splitCode got from cellPlan is zero
int HexCell::empty_resplit( const std::vector<char>& cellPlan){
  if(cellPlan.size() == 0){
    this->cellIndex = 1;
    return 1;
  }
  
  queue<HexCell*> Q;
  Q.push(this);
  HexCell* current;
  unsigned int index =0;
  int32 cIndex = 0;
  char currentCode=0;
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



int32 HexCell::num_fine_cells( const std::vector<char>& cellPlan)const{
  if(cellPlan.size() == 0){
    // this->cellIndex = 1;
    return 1;
  }
  
  queue<char> Q;
  unsigned int index =0;
  int32 cIndex = 0;
  Q.push(cellPlan[0]);
  char current;
  
  while(!Q.empty()){
    current = Q.front();
    if(current==0){
      cIndex++;
    }else{
      for(int i = 0; i <numChildren(current); i++){
        index++;
        if(index >= cellPlan.size()){
          Q.push(0);
        }else{
          Q.push(cellPlan[index]);
        }
      }
    }
    
    Q.pop();
  }
  return cIndex++;
}







//This function check if aCell is my neighbor is direction dd
// edge neighbor is excluded. if two faces overlap area is not zero
// their cells are neighbors
bool  HexCell::isNeighbor(const HexCell* aCell, DIRECTION dd)const{
  //exclude myself
  if(aCell == this) return false;
  return is_overlapped(face[dd], aCell->face[(dd%2)==0?(dd+1):(dd-1)]);
}
  

 

//this function check if aCell is my sibling neighbor.
//sibling means same size face, not necessarily has same parent
bool  HexCell::isSiblingNeighbor(const HexCell* aCell, DIRECTION dd)const{
  if(aCell == this) return false;
  return face[dd]==aCell->face[(dd%2)==0?(dd+1):(dd-1)];
}

//this function find neighbor in direction dd, the neighbor's face is greater or
// same as mine
HexCell* HexCell::findNeighbor(DIRECTION dd)
{
  
  HexCell* tempNeighbor = 0;
  if(parentCell == 0) return 0; //root has no neib

  int intCode = parentCell->mySplitCode;
  //find sibling neighbor(real sibling)
  switch(dd){
  case RIGHT:
  case LEFT:
    
    //check if parent cell is splitted in x direction 
    if(intCode & 4){
      for(int i = 0; i<parentCell->numChildren(); i++){
        if(isSiblingNeighbor(parentCell->childCell[i], dd))
          return parentCell->childCell[i];
      }
    }
    
    break;
    
  case FRONT:
  case BACK:

    //check if parent cell is splitted in y direction 
    if(intCode & 2)
      {
        for(int i = 0; i<parentCell->numChildren(); i++){
          if(isSiblingNeighbor(parentCell->childCell[i], dd))
            return parentCell->childCell[i];
        }
      }
    break;


  case UP:
  case DOWN:
    
    //check if parent cell is splitted in z direction 
    if(intCode & 1)
      {
        for(int i = 0; i<parentCell->numChildren(); i++){
          if( isSiblingNeighbor(parentCell->childCell[i], dd))
            return parentCell->childCell[i];
        }
      }
    break;
  default:
    cerr <<"WARNING: illegal DIRECTION in function finfNeighbor()" << endl;
    break;
  }
  


  //otherwise find the nearest common ancestor
  
  tempNeighbor = parentCell->findNeighbor(dd);
  if(tempNeighbor == 0) return 0;  //no neib
  

  //then climb down to find my neighbor
  while(tempNeighbor != 0){
    
    
    //if tempNeighbor is a leaf and is my sibling neighbor, return it
    if(isSiblingNeighbor(tempNeighbor, dd)&& tempNeighbor->mySplitCode == 0 ) return tempNeighbor;
    
    //if tempNeighbor's face has been smaller than mine, return tempneighbor's parentCell 
    switch(dd){
    case RIGHT:
    case LEFT:
      if((tempNeighbor->getLevel(YY) > this->getLevel(YY)) ||
         (tempNeighbor->getLevel(ZZ) > this->getLevel(ZZ)))return tempNeighbor->parentCell;
      break;
      
    case FRONT:
    case BACK:
      if((tempNeighbor->getLevel(XX) > this->getLevel(XX)) ||
         (tempNeighbor->getLevel(ZZ) > this->getLevel(ZZ)))return tempNeighbor->parentCell;
      break;
      
    case UP:
    case DOWN:
      if((tempNeighbor->getLevel(XX) > this->getLevel(XX)) ||
         (tempNeighbor->getLevel(YY) > this->getLevel(YY)))return tempNeighbor->parentCell;
      break;
      
    default:
      cerr << "WARNING: illegal DIRECTION in function findNeighbor()" << endl;
      break;
    }
    
    //if tempNeighbor is a leaf(whose face is greater than mine) , return
    if(tempNeighbor->mySplitCode == 0)return tempNeighbor; 
    
    //if tempNeighbor is not a leaf, and her face still greater than mine, climb down
    //caution: here if no child is neighbor, a dead loop will happen
    int i;
    int nc = tempNeighbor->numChildren(); 
    for( i = 0; i < tempNeighbor->numChildren(); i++){
      if( isNeighbor(tempNeighbor->childCell[i], dd)){
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
  
  return 0;
  
}


//this function check each cell in cells, find all the individual faces,put their connectivity info
// into faces. if necessary, split the face and put new points in nodes

// cells: in, all the leaves cells in the 2N tree
//nodes: in and out, initially the nodes after resplit(), when return, all the new points
//from split the faces are added
// faces: in and out, initially empty, return all the indivial faces and its two cell inex
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wchar-subscripts"
#endif

void set_hex_faces(const std::vector<HexCell*>& cells,
                   std::map<QuadFace*, NeibIndex>& faces){
  
  
  HexCell* tempNeib;
  int numCells = cells.size();
  
  if(numCells > 0){
    
    //for each cell
    for(int i = 0; i < numCells; i++){
      
      for(DIRECTION dd = RIGHT; dd <= DOWN; dd = static_cast<DIRECTION>(dd+1)){
        
        
        if(!((cells[i]->faceMarked).test(dd))){
          cells[i]->faceMarked.set(dd);
          
          tempNeib = cells[i]->findNeighbor(dd);
         
          //if no neib, don't output the face
          //there is a neib
          if(tempNeib != 0){
            bool isSibling = cells[i]->isSiblingNeighbor(tempNeib,dd);
            bool isALeaf =  (tempNeib->mySplitCode == 0);
            
          
            //if neib is a leaf, marked its face
            if(isALeaf){
              char nf = (dd%2==0?(dd+1):(dd-1));
              tempNeib ->faceMarked.set(nf);
              if(cells[i]->face[dd]->code != 0){
                cerr << "WARNING: face not a leaf" << endl;
              }
              
              if(dd==0 || dd==3 || dd==4) faces[cells[i]->face[dd]] =  NeibIndex(cells[i]->cellIndex, tempNeib ->cellIndex);
              else  faces[cells[i]->face[dd]] =  NeibIndex( tempNeib ->cellIndex,cells[i]->cellIndex);
            }
            
            else if(!isSibling) {
              //this stack hold all the neighbor leaves from  tempNeib
              stack<HexCell*> cellStack;
              
              HexCell* current=0;
              cellStack.push(tempNeib);
              
              while(!cellStack.empty()){
                //current takes out a cell from stack
                current = cellStack.top();
                cellStack.pop();
                
                int  numNeibChildren = current->numChildren();
                if( numNeibChildren != 0){
                  for(int j = 0; j < numNeibChildren; j++){
                    if(cells[i]->isNeighbor(current->childCell[j], dd)){
                      cellStack.push(current->childCell[j]);
                    }
                  }
                }
                else{
                  char nf = ((dd%2)==0?(dd+1):(dd-1));
                  //find the intersection of cells[i]->face[dd] and current->face[nf] leaves
                  std::vector<QuadFace*> commonfaces = overlap(cells[i]->face[dd], current->face[nf]);
                  if(commonfaces.size() != 1){
                    cerr << "WARNING: more than one commom face" << endl;
                  }
                  if(dd==0 || dd ==3 || dd ==4) faces[commonfaces[0]] = NeibIndex(cells[i]->cellIndex, current->cellIndex);
                  else faces[commonfaces[0]] = NeibIndex(current->cellIndex,cells[i]->cellIndex);
                }
                
                //if myself > neighbor do nothing
              }
            }//elseif(!isSibling)
            
          }//if(tempNeib !=0)
        }
      }
    }
  }
}
       

void HexCell::empty_split(){
  if(childCell!=0) return;
  switch(mySplitCode)
    {
      
      //000 no split,this is a leaf
    case 0:
      break;
      
      //100  x direction is splitted 
    case 4:
      childCell = new HexCell*[2];
      for(int i = 0; i <2; i++){
        childCell[i] = new HexCell();
      }
      
      break;
        
      //010  y direction is splitted 
    case 2:
      childCell = new HexCell*[2];
      for(int i = 0; i < 2; i++){
        childCell[i] = new HexCell();
      }
      break;
        
        
      //001  z direction is splitted 
    case 1:
      childCell = new HexCell*[2];
      for(int i = 0; i <2; i++){
        childCell[i] = new HexCell();
      }
      break;
              
      //011  y and z directions are splitted 
    case 3:
      childCell = new HexCell*[4];
      for(int i = 0; i < 4; i++){
        childCell[i] = new HexCell();
      }
      break;
      
      //101  x and z directions are splitted 
    case 5:
      childCell = new HexCell*[4];
      for(int i = 0; i < 4; i++){
        childCell[i] = new HexCell();
      }
      break;
        
      //110  y and z directions are splitted 
    case 6:
      childCell = new HexCell*[4];
      for(int i = 0; i < 4; i++){
        childCell[i] = new HexCell();
      }
      break;
      
      //111  x, y and z directions are splitted 
    case 7:
      
      childCell = new HexCell*[8];
      
      for(int i = 0; i < 8; i++){
	childCell[i] = new HexCell();
      }
      break;
      
    default:
      cerr <<"WARNING: illegal splitcode in function HexCell::reSplit()" << endl;
      break;
    }
}


  
void HexCell::split( std::list<Node*>& node_list,
                     std::list<Edge*>& edge_list, std::list<QuadFace*>& face_list){
  
  if(childCell != 0)return;
  
  QuadFace* newFace = 0;
  QuadFace* newface[12];
  Edge* newedge[6];
  Node* ccenter = 0;
  
  switch(mySplitCode)
    {
     
      //000
    case 0:
      break;
      
      //100  x direction is splitted 
    case 4:
    
      childCell = new HexCell*[2];
      for(int i = 0; i <2; i++){
        childCell[i] = new HexCell();
        childCell[i]->face = new QuadFace*[6];
          
        //	childCell[i]->whichChild = i;
        childCell[i]->parentCell = this;
      }
        
      //face 2-5 split in x direction
      for(int i = 2; i < 6; i++){
        face[i]->split(char(2),char(0), node_list, edge_list);
      }
        
        
      //define new face
      newFace = new QuadFace(4);
      face_list.push_back(newFace);
        
      newFace->edge[0] = face[5]->childx[0]->edge[1];
      childCell[0]->face[5] = face[5]->childx[0];
      childCell[1]->face[5] = face[5]->childx[1];
        
      newFace->edge[1] = face[2]->childx[0]->edge[1];
      childCell[0]->face[2] = face[2]->childx[0];
      childCell[1]->face[2] = face[2]->childx[1];
        
      newFace->edge[2] = face[4]->childx[0]->edge[1];
      childCell[0]->face[4] = face[4]->childx[0];
      childCell[1]->face[4] = face[4]->childx[1];
        
      newFace->edge[3] = face[3]->childx[0]->edge[1];
      childCell[0]->face[3] = face[3]->childx[0];
      childCell[1]->face[3] = face[3]->childx[1];
        
      childCell[0]->face[0] = newFace;
      childCell[0]->face[1] = face[1];
      childCell[1]->face[0] = face[0];
      childCell[1]->face[1] = newFace;
     
      break;
     
        
      //010  y direction is splitted 
    case 2:
      
      childCell = new HexCell*[2];
     
      for(int i = 0; i < 2; i++){
        childCell[i] = new HexCell();
        childCell[i]->face = new QuadFace*[6];
       
        //childCell[i]->whichChild = i;
        childCell[i]->parentCell = this;
      }
     
      //face 2-5 split 
     
     
      face[0]->split(char(2),char(0), node_list, edge_list);
      face[1]->split(char(2),char(0), node_list, edge_list);
      face[4]->split(char(1),char(0), node_list, edge_list);
      face[5]->split(char(1),char(0), node_list, edge_list);
      
      
      
      //define new face
      newFace = new QuadFace(4);
      face_list.push_back(newFace);
      
      newFace->edge[0] = face[5]->childy[1]->edge[0];
      childCell[0]->face[5] = face[5]->childy[0];
      childCell[1]->face[5] = face[5]->childy[1];
      
      newFace->edge[1] = face[0]->childx[0]->edge[1];
      childCell[0]->face[0] = face[0]->childx[0];
      childCell[1]->face[0] = face[0]->childx[1];
      
      newFace->edge[2] = face[4]->childy[1]->edge[0];
      childCell[0]->face[4] = face[4]->childy[0];
      childCell[1]->face[4] = face[4]->childy[1];
      
      newFace->edge[3] = face[1]->childx[0]->edge[1];
      childCell[0]->face[1] = face[1]->childx[0];
      childCell[1]->face[1] = face[1]->childx[1];
      
      childCell[0]->face[2] = newFace;
      childCell[0]->face[3] = face[3];
      childCell[1]->face[2] = face[2];
      childCell[1]->face[3] = newFace;
      
      break;
        
        
      //001  z direction is splitted 
    case 1:
      
      childCell = new HexCell*[2];
      
      for(int i = 0; i <2; i++){
        childCell[i] = new HexCell();
        childCell[i]->face = new QuadFace*[6];
   
        //	childCell[i]->whichChild = i;
        childCell[i]->parentCell = this;
      }

      //face 3,0, 2,1 split 
      
        
     
      face[0]->split(char(1),char(0), node_list, edge_list);
      face[1]->split(char(1),char(0), node_list, edge_list);
      face[2]->split(char(1),char(0), node_list, edge_list);
      face[3]->split(char(1),char(0), node_list, edge_list);
      
      //define new face
      newFace = new QuadFace(4);
      face_list.push_back(newFace);
      
     
      newFace->edge[0] = face[3]->childy[1]->edge[0];
      childCell[0]->face[3] = face[3]->childy[0];
      childCell[1]->face[3] = face[3]->childy[1];
      
      newFace->edge[1] = face[0]->childy[1]->edge[0];
      childCell[0]->face[0] = face[0]->childy[0];
      childCell[1]->face[0] = face[0]->childy[1];
      
      newFace->edge[2] = face[2]->childy[1]->edge[0];
      childCell[0]->face[2] = face[2]->childy[0];
      childCell[1]->face[2] = face[2]->childy[1];
      
      newFace->edge[3] = face[1]->childy[1]->edge[0];
      childCell[0]->face[1] = face[1]->childy[0];
      childCell[1]->face[1] = face[1]->childy[1];
      
      childCell[0]->face[5] = face[5];
      childCell[0]->face[4] = newFace;
      childCell[1]->face[4] = face[4];
      childCell[1]->face[5] = newFace;
      break;
      
      
      
      //011  y and z directions are splitted 
    case 3:
     
      childCell = new HexCell*[4];
     
      for(int i = 0; i < 4; i++){
        childCell[i] = new HexCell();
        childCell[i]->face = new QuadFace*[6];
        
        //	childCell[i]->whichChild = i;
        childCell[i]->parentCell = this;
      }

      //face split 
     
      face[0]->split(char(3),char(0), node_list, edge_list);
      face[1]->split(char(3),char(0), node_list, edge_list);
      face[2]->split(char(1),char(0), node_list, edge_list);
      face[3]->split(char(1),char(0), node_list, edge_list);
      face[4]->split(char(1),char(0), node_list, edge_list);
      face[5]->split(char(1),char(0), node_list, edge_list);
      
      
      //define new face
      
      for(int i = 0; i < 4; i++){
        newface[i] = new QuadFace(4);
        face_list.push_back(newface[i]);
      }
      
      //define new edge
      newedge[0] = new Edge(face[1]->getCenter(), face[0]->getCenter(), face[2]->edge[0]->level);
      edge_list.push_back(newedge[0]);
      
      //connect new edge
      newface[0]->edge[2] = newface[1]->edge[0] = newface[2]->edge[0] = newface[3]->edge[2]  = newedge[0];

      //connect short edges
      newface[0]->edge[1] = face[0]->child[0]->edge[1];
      newface[1]->edge[1] = face[0]->child[3]->edge[0];
      newface[2]->edge[1] = face[0]->child[3]->edge[3];
      newface[3]->edge[1] = face[0]->child[0]->edge[2];

      newface[0]->edge[3] = face[1]->child[0]->edge[1];
      newface[1]->edge[3] = face[1]->child[3]->edge[0];
      newface[2]->edge[3] = face[1]->child[3]->edge[3];
      newface[3]->edge[3] = face[1]->child[0]->edge[2];

      //connect square face
      for(int i = 0; i < 4; i++){
        childCell[i]->face[0] = face[0]->child[i];
        childCell[i]->face[1] = face[1]->child[i];
      }
      
      //connect newface
      childCell[0]->face[2] = childCell[2]->face[3] = newface[0];
      childCell[2]->face[4] = childCell[3]->face[5] = newface[1];
      childCell[1]->face[2] = childCell[3]->face[3] = newface[2];
      childCell[0]->face[4] = childCell[1]->face[5] = newface[3];

      //create and connect rectangle faces
      
      newface[0]->edge[0] = face[5]->childy[1]->edge[0];
      childCell[0]->face[5] = face[5]->childy[0];
      childCell[2]->face[5] = face[5]->childy[1];
     

      
      newface[1]->edge[2] = face[2]->childy[1]->edge[0];
      childCell[2]->face[2] = face[2]->childy[0];
      childCell[3]->face[2] = face[2]->childy[1];
     
      
      newface[2]->edge[2] = face[4]->childy[1]->edge[0];
      childCell[1]->face[4] = face[4]->childy[0];
      childCell[3]->face[4] = face[4]->childy[1];
      
      newface[3]->edge[0] = face[3]->childy[1]->edge[0];
      childCell[0]->face[3] = face[3]->childy[0];
      childCell[1]->face[3] = face[3]->childy[1];
      
      break;
      
      
      
      //101  x and z directions are splitted 
    case 5:
        
      childCell = new HexCell*[4];
        
      for(int i = 0; i < 4; i++){
        childCell[i] = new HexCell();
        childCell[i]->face = new QuadFace*[6];
          
        //	childCell[i]->whichChild = i;
        childCell[i]->parentCell = this;
      }

      //face split 
      face[2]->split(char(3),char(0), node_list, edge_list);
      face[3]->split(char(3),char(0), node_list, edge_list);
      face[0]->split(char(1),char(0), node_list, edge_list);
      face[1]->split(char(1),char(0), node_list, edge_list);
      face[4]->split(char(2),char(0), node_list, edge_list);
      face[5]->split(char(2),char(0), node_list, edge_list);
      
      
      //define new face
     
      for(int i = 0; i < 4; i++){
        newface[i] = new QuadFace(4);
        face_list.push_back(newface[i]);
      }
        
      //define new edge
      newedge[0] = new Edge(face[3]->getCenter(), face[2]->getCenter(), face[0]->edge[0]->level);
      edge_list.push_back(newedge[0]);

      //connect newedge
      newface[0]->edge[2] = newface[1]->edge[3] = newface[2]->edge[0] = newface[3]->edge[1]  = newedge[0];
      //connect short edges
      newface[0]->edge[1] = face[2]->child[0]->edge[1];
      newface[1]->edge[2] = face[2]->child[3]->edge[0];
      newface[2]->edge[1] = face[2]->child[3]->edge[3];
      newface[3]->edge[2] = face[2]->child[0]->edge[2];

      newface[0]->edge[3] = face[3]->child[0]->edge[1];
      newface[1]->edge[0] = face[3]->child[3]->edge[0];
      newface[2]->edge[3] = face[3]->child[3]->edge[3];
      newface[3]->edge[0] = face[3]->child[0]->edge[2];

      //connect square face
      for(int i = 0; i < 4; i++){
        childCell[i]->face[2] = face[2]->child[i];
        childCell[i]->face[3] = face[3]->child[i];
      }

      //connect newface
      childCell[0]->face[0] = childCell[2]->face[1] = newface[0];
      childCell[2]->face[4] = childCell[3]->face[5] = newface[1];
      childCell[1]->face[0] = childCell[3]->face[1] = newface[2];
      childCell[0]->face[4] = childCell[1]->face[5] = newface[3];

      //create and connect rectangle faces
       
      newface[0]->edge[0] = face[5]->childx[0]->edge[1];
      childCell[0]->face[5] = face[5]->childx[0];
      childCell[2]->face[5] = face[5]->childx[1];
       
      newface[1]->edge[1] = face[0]->childy[1]->edge[0];
      childCell[2]->face[0] = face[0]->childy[0];
      childCell[3]->face[0] = face[0]->childy[1];
     
      newface[2]->edge[2] = face[4]->childx[0]->edge[1];
      childCell[1]->face[4] = face[4]->childx[0];
      childCell[3]->face[4] = face[4]->childx[1];
              
      newface[3]->edge[3] = face[1]->childy[1]->edge[0];
      childCell[0]->face[1] = face[1]->childy[0];
      childCell[1]->face[1] = face[1]->childy[1];
       
        
      break;
        
      //110  y and z directions are splitted 
    case 6:

      childCell = new HexCell*[4];
        
      for(int i = 0; i < 4; i++){
        childCell[i] = new HexCell();
        childCell[i]->face = new QuadFace*[6];
   
        //	childCell[i]->whichChild = i;
        childCell[i]->parentCell = this;
      }
      //face split 
      face[4]->split(char(3),char(0), node_list, edge_list);
      face[5]->split(char(3),char(0), node_list, edge_list);
      face[0]->split(char(2),char(0), node_list, edge_list);
      face[1]->split(char(2),char(0), node_list, edge_list);
      face[2]->split(char(2),char(0), node_list, edge_list);
      face[3]->split(char(2),char(0), node_list, edge_list);
      
      
      //define new face
     
      for(int i = 0; i < 4; i++){
        newface[i] = new QuadFace(4);
        face_list.push_back(newface[i]);
      }
        
      //define new edge
      newedge[0] = new Edge(face[5]->getCenter(), face[4]->getCenter(), face[0]->edge[1]->level);
      edge_list.push_back(newedge[0]);
        
      //connect newedge
      newface[0]->edge[1] = newface[1]->edge[3] = newface[2]->edge[3] = newface[3]->edge[1]  = newedge[0];
      //connect short edges
      newface[0]->edge[2] = face[4]->child[0]->edge[1];
      newface[1]->edge[2] = face[4]->child[3]->edge[0];
      newface[2]->edge[2] = face[4]->child[3]->edge[3];
      newface[3]->edge[2] = face[4]->child[0]->edge[2];

      newface[0]->edge[0] = face[5]->child[0]->edge[1];
      newface[1]->edge[0] = face[5]->child[3]->edge[0];
      newface[2]->edge[0] = face[5]->child[3]->edge[3];
      newface[3]->edge[0] = face[5]->child[0]->edge[2];

      //connect square face
      for(int i = 0; i < 4; i++){
        childCell[i]->face[4] = face[4]->child[i];
        childCell[i]->face[5] = face[5]->child[i];
      }

      //connect newface
      childCell[0]->face[0] = childCell[2]->face[1] = newface[0];
      childCell[2]->face[2] = childCell[3]->face[3] = newface[1];
      childCell[1]->face[0] = childCell[3]->face[1] = newface[2];
      childCell[0]->face[2] = childCell[1]->face[3] = newface[3];

      //create and connect rectangle faces
    
      newface[0]->edge[3] = face[3]->childx[0]->edge[1];
      childCell[0]->face[3] = face[3]->childx[0];
      childCell[2]->face[3] = face[3]->childx[1];
        
      newface[1]->edge[1] = face[0]->childx[0]->edge[1];
      childCell[2]->face[0] = face[0]->childx[0];
      childCell[3]->face[0] = face[0]->childx[1];
        
        
      newface[2]->edge[1] = face[2]->childx[0]->edge[1];
      childCell[1]->face[2] = face[2]->childx[0];
      childCell[3]->face[2] = face[2]->childx[1];
        
      newface[3]->edge[3] = face[1]->childx[0]->edge[1];
      childCell[0]->face[1] = face[1]->childx[0];
      childCell[1]->face[1] = face[1]->childx[1];
        
        
      break;
        
      //111  x, y and z directions are splitted 
    case 7:

      childCell = new HexCell*[8];
        
      for(int i = 0; i < 8; i++){
        childCell[i] = new HexCell();
        childCell[i]->face = new QuadFace*[6];
   
        //	childCell[i]->whichChild = i;
        childCell[i]->parentCell = this;
      }

      //face split 
      for(int i = 0; i< 6; i++)face[i]->split(char(3), char(0), node_list, edge_list);
      
      //define new node
      ccenter = centroid();
      node_list.push_back(ccenter);
        
      //define new face
        
      for(int i = 0; i < 12; i++){
        newface[i] = new QuadFace(4);
        face_list.push_back(newface[i]);
      }
        
      //define new edge
      for(int i =0; i<6; i++){
        if(i%2==1){
          newedge[i] = new Edge(face[i]->getCenter(), ccenter);
        }
        else{
          newedge[i] = new Edge(ccenter, face[i]->getCenter());
        }
        edge_list.push_back(newedge[i]);
      }
      newedge[0]->level = newedge[1]->level = face[2]->edge[0]->level+1;
      newedge[2]->level = newedge[3]->level = face[0]->edge[0]->level+1; 
      newedge[4]->level = newedge[5]->level = face[0]->edge[1]->level+1;
        


      newface[0]->edge[0] = face[5]->child[0]->edge[1];
      newface[0]->edge[1] = newedge[5];
      newface[0]->edge[2] = newedge[3];
      newface[0]->edge[3] = face[3]->child[0]->edge[1];

      newface[1]->edge[0] = newedge[3];
      newface[1]->edge[1] = newedge[4];
      newface[1]->edge[2] = face[4]->child[0]->edge[1];
      newface[1]->edge[3] = face[3]->child[1]->edge[1];

      newface[2]->edge[0] = face[5]->child[1]->edge[1];
      newface[2]->edge[1] = face[2]->child[0]->edge[1];
      newface[2]->edge[2] = newedge[2];
      newface[2]->edge[3] = newedge[5];

      newface[3]->edge[0] = newedge[2];
      newface[3]->edge[1] = face[2]->child[1]->edge[1];
      newface[3]->edge[2] = face[4]->child[1]->edge[1];
      newface[3]->edge[3] = newedge[4];

      newface[4]->edge[0] = face[5]->child[0]->edge[2];
      newface[4]->edge[1] = newedge[5];
      newface[4]->edge[2] = newedge[1];
      newface[4]->edge[3] = face[1]->child[0]->edge[1];

      newface[5]->edge[0] = newedge[1];
      newface[5]->edge[1] = newedge[4];
      newface[5]->edge[2] = face[4]->child[0]->edge[2];
      newface[5]->edge[3] = face[1]->child[1]->edge[1];

      newface[6]->edge[0] = face[5]->child[2]->edge[2];
      newface[6]->edge[1] = face[0]->child[0]->edge[1];
      newface[6]->edge[2] = newedge[0];
      newface[6]->edge[3] = newedge[5];

      newface[7]->edge[0] = newedge[0];
      newface[7]->edge[1] = face[0]->child[1]->edge[1];
      newface[7]->edge[2] = face[4]->child[2]->edge[2];
      newface[7]->edge[3] = newedge[4];

      newface[8]->edge[0] = face[3]->child[0]->edge[2];
      newface[8]->edge[1] = newedge[3];
      newface[8]->edge[2] = newedge[1];
      newface[8]->edge[3] = face[1]->child[0]->edge[2];

      newface[9]->edge[0] = newedge[1];
      newface[9]->edge[1] = newedge[2];
      newface[9]->edge[2] = face[2]->child[0]->edge[2];
      newface[9]->edge[3] = face[1]->child[2]->edge[2];

      newface[10]->edge[0] = face[3]->child[2]->edge[2];
      newface[10]->edge[1] = face[0]->child[0]->edge[2];
      newface[10]->edge[2] = newedge[0];
      newface[10]->edge[3] = newedge[3];

      newface[11]->edge[0] = newedge[0];
      newface[11]->edge[1] = face[0]->child[2]->edge[2];
      newface[11]->edge[2] = face[2]->child[2]->edge[2];
      newface[11]->edge[3] = newedge[2];
        
      //connect square face
      for(int i = 0; i < 4; i++){
        childCell[i]->face[1] = face[1]->child[i];
      }
      for(int i = 4; i < 8; i++){
        childCell[i]->face[0] = face[0]->child[i-4];
      }

      childCell[2]->face[2] = face[2]->child[0];
      childCell[3]->face[2] = face[2]->child[1];
      childCell[6]->face[2] = face[2]->child[2];
      childCell[7]->face[2] = face[2]->child[3];

      childCell[0]->face[3] = face[3]->child[0];
      childCell[1]->face[3] = face[3]->child[1];
      childCell[4]->face[3] = face[3]->child[2];
      childCell[5]->face[3] = face[3]->child[3];

      childCell[1]->face[4] = face[4]->child[0];
      childCell[3]->face[4] = face[4]->child[1];
      childCell[5]->face[4] = face[4]->child[2];
      childCell[7]->face[4] = face[4]->child[3];
        
      childCell[0]->face[5] = face[5]->child[0];
      childCell[2]->face[5] = face[5]->child[1];
      childCell[4]->face[5] = face[5]->child[2];
      childCell[6]->face[5] = face[5]->child[3];

        
      //connect newface
      childCell[0]->face[0] = childCell[4]->face[1] = newface[0];
      childCell[1]->face[0] = childCell[5]->face[1] = newface[1];
      childCell[2]->face[0] = childCell[6]->face[1] = newface[2];
      childCell[3]->face[0] = childCell[7]->face[1] = newface[3];

      childCell[0]->face[2] = childCell[2]->face[3] = newface[4];
      childCell[1]->face[2] = childCell[3]->face[3] = newface[5];
      childCell[4]->face[2] = childCell[6]->face[3] = newface[6];
      childCell[5]->face[2] = childCell[7]->face[3] = newface[7];

      childCell[0]->face[4] = childCell[1]->face[5] = newface[8];
      childCell[2]->face[4] = childCell[3]->face[5] = newface[9];
      childCell[4]->face[4] = childCell[5]->face[5] = newface[10];
      childCell[6]->face[4] = childCell[7]->face[5] = newface[11];

        
      break;
        
    default:
      cerr <<"WARNING: illegal splitcode in function split()" << endl;
      break;
    }
}
  
int HexCell::get_num_fine_faces()const{
  int num_faces = 0;
 
  for(int i = 0; i < 6; i++){
    num_faces += face[i]->get_num_leaves();
  }
  return num_faces;
}

int HexCell::get_tagged(){
  std::vector<Node*> nodes(8);
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
      
int HexCell::get_tagged(const vector<source_par>& sources){
  std::vector<Node*> nodes(8);
  get_nodes(nodes);
  
  double min_len = get_min_edge_length();
  if(tag_cell(nodes, sources, min_len)){
    mySplitCode = 7;
    return 1;
  }else{
    mySplitCode = 0;
    return 0;
  }
  return 0;
} 
//find the minimum edge length in a cell(before split)
double HexCell::get_min_edge_length(){
  std::vector<Edge*> edges = get_edges();
  std::vector<double> edge_length(12);
  for(int i = 0; i < 12; i++)edge_length[i] = edges[i]->get_length();
    
  double min_length = edge_length[0];
  for(int i = 1; i < 12; i++){
    min_length = min(min_length, edge_length[i]);
  }
  return min_length;
} 


//return a splitCode
//find   min_edge_length in XX , YY and ZZ directions
//find minimun_edge_length in all directions
  
//if 
//if max_edge_length/min_edge_length > Globals::factor1 and they are in different direction
//split the max_length edge
void HexCell::setSplitCode(int split_mode, double tol){

   
  //face2edge in build_hexcell
  // int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
  //              {1, 7, 3, 5}, {0, 6, 2, 4}};
    
  //find the edge_length of all edges
  std::vector<double> edge_length(12);
  edge_length[0] = face[3]->edge[0]->get_length();
  edge_length[1] = face[3]->edge[2]->get_length();
  edge_length[2] = face[2]->edge[0]->get_length();
  edge_length[3] = face[2]->edge[2]->get_length();
  edge_length[4] = face[1]->edge[0]->get_length();
  edge_length[5] = face[1]->edge[2]->get_length();
  edge_length[6] = face[0]->edge[0]->get_length();
  edge_length[7] = face[0]->edge[2]->get_length();
  edge_length[8] = face[1]->edge[3]->get_length();
  edge_length[9] = face[1]->edge[1]->get_length();
  edge_length[10] = face[0]->edge[3]->get_length();
  edge_length[11] = face[0]->edge[1]->get_length();
    

    
    
  //find the min_edge_length in XX direction
  std::vector<double> average_length(3);
  std::vector<double> min_length(3);
    
  average_length[0]  = 0.0;
  min_length[0] = 1e10;
  for(int i = 0; i<=3; i++){
    average_length[0] +=  edge_length[i];
    min_length[0] = min(min_length[0], edge_length[i]);
  }
  average_length[0] /= 4.0;
    
  average_length[1] = 0.0;
  min_length[1] = 1e10;
  for(int i = 4; i <= 7; i++){
    average_length[1] +=  edge_length[i];
    min_length[1] = min(min_length[1], edge_length[i]);
  }
  average_length[1] /= 4.0;
    
  average_length[2] = 0.0;
  min_length[2] = 1e10;
  for(int i = 8; i <= 11; i++){
    average_length[2] += edge_length[i];
    min_length[2] = min(min_length[2], edge_length[i]);
  }
  average_length[2] /= 4.0;

  bitset<3> tolerance_mask(7); //all 1s
  if(min_length[0] < 2*tol)tolerance_mask.reset(2);
  if(min_length[1] < 2*tol)tolerance_mask.reset(1);
  if(min_length[2] < 2*tol) tolerance_mask.reset(0);
  //    cout << "tolerance: " << tolerance_mask.to_ulong() << endl;
    
  double minimum_length = *std::min_element(min_length.begin(), min_length.end());
    
  if(split_mode == 0){
    if(mySplitCode != 0){
      cerr<<"WARNING:mySplitCode is not zero in setSplitCode of HexCell" << endl;
      exit(0);
    }
    
    if((average_length[2]/average_length[0] > Globals::factor)
       && (average_length[1]/average_length[0] >Globals::factor)){

      bitset<3> oldCode(3);
      oldCode = oldCode & tolerance_mask;
      mySplitCode = char(oldCode.to_ulong());
      //yz direction get split
        
      return;
    }
      
    else if((average_length[2]/average_length[1] > Globals::factor)
            && (average_length[0]/average_length[1] >Globals::factor) ){

      bitset<3> oldCode(5);
      oldCode = oldCode & tolerance_mask;
      mySplitCode = char(oldCode.to_ulong());
      //mySplitCode = 5;//xz direction get split
      return;
    }
      
    else if((average_length[1]/average_length[2] > Globals::factor)
            && (average_length[0]/average_length[2] >Globals::factor) ){
      bitset<3> oldCode(6);
      oldCode = oldCode & tolerance_mask;
      mySplitCode = char(oldCode.to_ulong());
      // mySplitCode = 6;//xy direction get split
      return;
    }
      
    else  if(((average_length[0]/average_length[1] > Globals::factor)
              && (average_length[2] > average_length[1])) ||
             ((average_length[0]/average_length[2] > Globals::factor)
              && (average_length[1] > average_length[2]))){
      bitset<3> oldCode(4);
      oldCode = oldCode & tolerance_mask;
      mySplitCode = char(oldCode.to_ulong());
      // mySplitCode = 4;//x direction get split
      return;
    }
      
    else  if(((average_length[1]/average_length[2] > Globals::factor)
              && (average_length[0] > average_length[2])) ||
             ((average_length[1]/average_length[0] > Globals::factor)
              && (average_length[2] > average_length[0]))){
      bitset<3> oldCode(2);
      oldCode = oldCode & tolerance_mask;
      mySplitCode = char(oldCode.to_ulong());
      // mySplitCode = 2;//y direction get split
      return;
    }
      
    else  if(((average_length[2]/average_length[0] > Globals::factor)
              && (average_length[1] > average_length[0])) ||
             ((average_length[2]/average_length[1] > Globals::factor)
              && (average_length[0] > average_length[1]))){
      bitset<3> oldCode(1);
      oldCode = oldCode & tolerance_mask;
      mySplitCode = char(oldCode.to_ulong());
      // mySplitCode = 1;//z direction get split
      return;
    }
    else if(minimum_length > 2.0*tol){
        
      mySplitCode = 7;
      return;
    }
    else{
      mySplitCode = 0;
      return;
    }
    cerr<< " WARNING: reach dummy code in setSplitCode()" << endl;
  }//end of if(split_mode==0)
  //when split_specified
  else if(split_mode == 1){
    if(mySplitCode != 0){
      cerr<<"WARNING:mySplitCode is not zero in setSplitCode of hexcell" << endl;
      exit(0);
    }
    //   setSplitCode(0);

      
    //face2edge in build_hexcell
    // int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
    //              {1, 7, 3, 5}, {0, 6, 2, 4}};
      
    //find  all edges
    vect3d edges[12];
      
    edges[0] = face[3]->edge[0]->head->p -  face[3]->edge[0]->tail->p;
    edges[1] = face[3]->edge[2]->head->p - face[3]->edge[2]->tail->p; 
    edges[2] = face[2]->edge[0]->head->p - face[2]->edge[0]->tail->p;
    edges[3] = face[2]->edge[2]->head->p - face[2]->edge[2]->tail->p;
    edges[4] = face[1]->edge[0]->head->p - face[1]->edge[0]->tail->p ;
    edges[5] = face[1]->edge[2]->head->p - face[1]->edge[2]->tail->p ;
    edges[6] = face[0]->edge[0]->head->p - face[0]->edge[0]->tail->p ;
    edges[7] = face[0]->edge[2]->head->p - face[0]->edge[2]->tail->p;
    edges[8] = face[1]->edge[3]->head->p -  face[1]->edge[3]->tail->p ;
    edges[9] = face[1]->edge[1]->head->p - face[1]->edge[1]->tail->p;
    edges[10] = face[0]->edge[3]->head->p - face[0]->edge[3]->tail->p ;
    edges[11] = face[0]->edge[1]->head->p -  face[0]->edge[1]->tail->p ;
    //  for(int i = 0 ; i < 12; i++) normalize(edges[i]);
    //since edges are built in the same direction,no edge need reverse
    //it's OK to average the directions
    vect3d average_direction[3];
    average_direction[0]  = vect3d(0.0, 0.0, 0.0);
    for(int i = 0; i<=3; i++) average_direction[0] = average_direction[0] +  edges[i];
    average_direction[0] = average_direction[0] / 4.0;
         
    average_direction[1] = vect3d(0.0, 0.0, 0.0);
    for(int i = 4; i <= 7; i++)average_direction[1] = average_direction[1] +  edges[i];
    average_direction[1] = average_direction[1] / 4.0;
      
    average_direction[2] = vect3d(0.0, 0.0, 0.0);
    for(int i = 8; i <= 11; i++)average_direction[2] = average_direction[2] + edges[i];
    average_direction[2] = average_direction[2] / 4.0;
      
    for(int i = 0; i < 3; i++)normalize(average_direction[i]);
      
    // bitset<3> oldCode(mySplitCode);
    //    bitset<3> mask(7); //all bits are 1s;
    int z_direction =     -1;
    for(int i = 0; i < 3; i++){
      if(abs(average_direction[i].x) < 0.01 && abs(average_direction[i].y) < 0.01) {
        z_direction = i;
        break;
      }
    }
    if(z_direction == -1) {
      cerr<< "WARNING: can not find z direction in hexcell.cc" << endl;
      exit(0);
    }

    bitset<3> oldCode(0);
    switch(z_direction){
    case 0:
      if(average_length[1]/average_length[2] > Globals::factor) oldCode.set(1);
      else if(average_length[2]/average_length[1] > Globals::factor) oldCode.set(0);
      else{
        oldCode.set(1);
        oldCode.set(0);
      }
      break;
    case 1:
      if(average_length[0]/average_length[2] > Globals::factor) oldCode.set(2);
      else if(average_length[2]/average_length[0] > Globals::factor) oldCode.set(0);
      else{
        oldCode.set(2);
        oldCode.set(0);
      }


        
      break;
    case 2:

      if(average_length[0]/average_length[1] > Globals::factor) oldCode.set(2);
      else if(average_length[1]/average_length[0] > Globals::factor) oldCode.set(1);
      else{
        oldCode.set(2);
        oldCode.set(1);
      }

      break;
    default:
      cerr<<"WARNING: invalid z_direction" << endl;
      break;
    }
    oldCode = oldCode  & tolerance_mask;
    mySplitCode = char(oldCode.to_ulong());
     
    return;
  }//end if(split_mode == 1)
    
  else if(split_mode == 2){
    if(minimum_length/tol  > 2.0){
      mySplitCode = 7;
      return;
    }
    else{
      mySplitCode = 0;
      return;
    }
  }
    
 
 
}

//after a cell is split, compose the cell plan according to the tree structure      
std::vector<char> HexCell::make_cellplan(){
  
  std::vector<char> cellPlan;
  
  if(childCell == 0) {
    reduce_vector(cellPlan);
    return cellPlan;
  }
 
  std::queue<HexCell*> Q;
  Q.push(this);
  
  HexCell* current;
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
std::vector<char> HexCell::make_cellplan(int level){
  
  std::vector<char> cellPlan;
 
  
  if(level <= 0) {
    reduce_vector(cellPlan);
    return cellPlan;
  }

  cellPlan.push_back(7);
  for(int i = 1; i < level; i++){
    int total_num_children = 1 << (3*i);
    for(int j = 0; j < total_num_children; j++) cellPlan.push_back(7);
  }
  reduce_vector(cellPlan);
  return cellPlan;
}
bool HexCell::balance_cell(int split_mode,
                           std::list<Node*>& node_list,
                           std::list<Edge*>& edge_list,
                           std::list<QuadFace*>& face_list){ 

 
  bool  needBalance = false;
  if(childCell != 0){

    std::list<HexCell*> leaves;
    sort_leaves(leaves);
    needBalance = false;
    for(std::list<HexCell*>::const_iterator p = leaves.begin(); p != leaves.end(); p++){
      bool tmp =  (*p)->balance_cell(split_mode, node_list, edge_list, face_list);
      needBalance = tmp||needBalance;
    }
    return needBalance;
  }
  if(childCell == 0){
    std::vector<Edge*> edge = get_edges();
    needBalance = false;
    if(split_mode == 2){
      for(int i = 0; i < 12; i++){
        if( edge[i]->depth_greater_than_1()){
          mySplitCode = 7;
          break;
        }
      }
      if(Globals::balance_option >= 1 && mySplitCode == 0){
        int num_faces_split = 0;
        for(int i = 0; i < 6; i++){
          if(face[i]-> code != 0 ) num_faces_split++;
        }
        if(num_faces_split > 3){
          
          mySplitCode = 7;
        }
      }

      if(Globals::balance_option >= 2 && mySplitCode == 0){
        if((face[0]->code != 0 && face[1]->code != 0)||
           (face[2]->code != 0 && face[3]->code != 0)||
           (face[4]->code != 0 && face[5]->code != 0))                                               
          
          mySplitCode = 7;
      }
      
    }
    else{
      bitset<3> code;
      for(int i = 0; i < 4; i++){
        if( edge[i]->depth_greater_than_1()){
          code.set(2);
        }
      }
      for(int i = 4; i < 8; i++){
        if( edge[i]->depth_greater_than_1()){
          code.set(1);
        }
      }
      for(int i = 8; i < 12; i++){
        if( edge[i]->depth_greater_than_1()){
          code.set(0);
        }
      }
      mySplitCode = char(code.to_ulong());

      if(Globals::balance_option >= 1 && mySplitCode == 0){
        int num_faces_split = 0;
        for(int i = 0; i < 6; i++){
          if(face[i]-> child != 0 ) num_faces_split++;
        }
        if(num_faces_split > 3){
          bitset<3> the_code;
          bitset<2> facexy = bitset<2>(face[4]->code) | bitset<2>(face[5]->code);
          bitset<2> faceyz = bitset<2>(face[0]->code) | bitset<2>(face[1]->code);
          bitset<2> facexz = bitset<2>(face[2]->code) | bitset<2>(face[3]->code);
          the_code.reset();
          the_code[0] = faceyz[0] | facexz[0];//x bit
          the_code[1] = facexy[0] | faceyz[1];//y bit
          the_code[2] = facexy[1] | facexz[1];//z bit
          
          mySplitCode = char(the_code.to_ulong());
          
        }
      }
      
      if(Globals::balance_option >= 2 && mySplitCode == 0){

        bitset<3> the_code;
        bitset<2> facexy = bitset<2>(face[4]->code) | bitset<2>(face[5]->code);
        bitset<2> faceyz = bitset<2>(face[0]->code) | bitset<2>(face[1]->code);
        bitset<2> facexz = bitset<2>(face[2]->code) | bitset<2>(face[3]->code);
        the_code.reset();
        if(face[0]->child != 0 && face[1]->child != 0){//faceyz
          the_code[1] = the_code[1] | faceyz[1];
          the_code[0] = the_code[0] | faceyz[0];
        }
        if (face[2]->child != 0 && face[3]->child != 0){//facexz
          the_code[2] = the_code[2] | facexz[1];
          the_code[0] = the_code[0] | facexz[0];
        }
        if(face[4]->child != 0 && face[5]->child != 0){//facexy

          the_code[2] = the_code[2] | facexy[1];
          the_code[1] = the_code[1] | facexy[0];
        }
          
        mySplitCode = char(the_code.to_ulong());  
      }
      
    }
    if(mySplitCode !=0){
      needBalance = true;
      
      split(node_list, edge_list, face_list);
    
      for(int i = 0; i < numChildren(); i++){
        childCell[i]->balance_cell(split_mode, node_list, edge_list, face_list);
      }
    }
  }

 
  return needBalance;
}

void HexCell::sort_leaves(std::list<HexCell*>& leaves){
  if(childCell != 0){
    for(int i = 0; i < numChildren(); i++)childCell[i]->sort_leaves(leaves);
  }
  else{
    leaves.push_back(this);
  }
}

void HexCell::rebalance_cells(int split_mode,
                              std::list<Node*>& node_list,
                              std::list<Edge*>& edge_list,
                              std::list<QuadFace*>& face_list){ 
  bool need_balance_more = true;

  while(need_balance_more){
    need_balance_more = balance_cell(split_mode, node_list, edge_list, face_list); 
  }
}

int32 HexCell::traverse(const std::vector<char>& parentPlan,  vector<pair<int32, int32> >& indexMap){
  indexMap.clear();
  if(parentPlan.size() == 0){
    if(numChildren()!=0){
      list<HexCell*> leaves;
      sort_leaves(leaves); 
      for(std::list<HexCell*>::const_iterator p = leaves.begin(); p != leaves.end(); p++)
        indexMap.push_back(make_pair((*p)->cellIndex, 1));
      return 1;
    }else{
      indexMap.push_back(make_pair(1,1));
      return 1;
    }
  }
  std::queue<HexCell*> Q;
  Q.push(this);
  HexCell* current;
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
    list<HexCell*> leaves;
    switch(currentCode){
      //0 no split,this is a leaf, this is a leaf for parentPlan
    case 0:
      ++cIndex;
      current->sort_leaves(leaves);
      if(leaves.front() != current){//current is not a leaf for cellPlan
        for(std::list<HexCell*>::const_iterator p = leaves.begin(); p != leaves.end(); p++){
          indexMap.push_back(make_pair((*p)->cellIndex, cIndex));
        }
      }else{
        //current is a leaf for cellPlan and for parentPlan
        if(current->cellIndex != 0) indexMap.push_back(make_pair(current->cellIndex, cIndex));  
      }
      break;
    default:
      if(current->cellIndex !=0){//derefinement happen,current is a leaf for cellPlan
        current->sort_leaves(leaves);
        for(std::list<HexCell*>::const_iterator p = leaves.begin(); p != leaves.end(); p++){
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







  
//assume with derefinement, the balance option is always no edge has levels greater than 1 
bool HexCell::needDerefine(){
  if(childCell != 0){
    bool derefine = true;
    for(int i = 0; i < numChildren(); i++){
      if( (childCell[i] ->get_tagged()) != 2)return false;
      if(childCell[i]->childCell!=0) return false;
    }
    if(derefine){
      std::vector<Edge*> edge = get_edges();
      for(int i = 0; i < 12; i++){
	if( edge[i]->depth_greater_than_1())return false;
      }
      return true;
    }
  }
  return false;
}
bool HexCell::needDerefine_ctag(){
  if(childCell != 0){
    bool derefine = true;
    for(int i = 0; i < numChildren(); i++){
      if( (childCell[i] ->getTag()) != 2)return false;
      if(childCell[i]->childCell!=0) return false;
    }
    if(derefine){
      std::vector<Edge*> edge = get_edges();
      for(int i = 0; i < 12; i++){
	if( edge[i]->depth_greater_than_1())return false;
      }
      return true;
    }
  }
  return false;
}


      
void HexCell::derefine(){
  if(childCell != 0){
    for(int i = 0; i < numChildren(); i++){
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
void reorder_faces(const const_store<int>& node_remap, std::vector<Entity>& lower){
  
  //reverse the map 
  std::vector<pair<int, Entity> > node_f2l(lower.size());
  for(unsigned int  index  = 0; index < lower.size(); index++){
    node_f2l[index] = pair<int, Entity>(node_remap[lower[index]],lower[index]);
  }
  std::sort(node_f2l.begin(), node_f2l.end());
  
  for( unsigned int i= 0; i < lower.size(); i++){
    lower[i] = node_f2l[i].second;
  }
  
}


void reorder_faces(const const_store<int>& node_remap, std::vector<Entity>& lower,
                   std::vector<Entity>& upper,
                   std::vector<Entity>& boundary_map){
  
  reorder_faces(node_remap, lower);
  reorder_faces(node_remap, upper);
  reorder_faces(node_remap, boundary_map);
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
