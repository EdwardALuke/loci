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
///////////////////////////////////////////////////////////////////////////////
//                                diamondcell.cc
//                                by: Qiuhan Xue
//    This function include definition of class DiamondCell,
//    and class Face
//
///////////////////////////////////////////////////////////////////////////////
 
#include <queue>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <set>
#include <stack>
#include "diamondcell.h"
#include "defines.h"
#include "globals.h"

using std::queue;
using std::cerr;
using std::endl;
using std::cout;



//this function splits diamondcell isotropially once
// newly created nodes, edges and faces are put into the lists
void DiamondCell::split(std::list<Node*>& node_list,
                        std::list<Edge*>& edge_list,
                        std::list<Face*>& face_list){

  if(childCell!=0)return;
  //split each face
  for(int i = 0; i < 2*nfold; i++){
    if(face[i]->child ==0) face[i]->split(node_list,edge_list);
  }
  
  //calculate the cellcenter and put it in node_list
  Node*  cellcenter = centroid();
  node_list.push_back(cellcenter);
  
  //get facecenter
  vector<Node*> facecenter(2*nfold);
  getFaceCenter(&facecenter[0]);
  
  //get edgecenter
  vector<Node*> edgecenter(4*nfold);
  getEdgeCenter(&edgecenter[0]);
  
  //defines new edges
  vector<Edge*> newEdge(2*nfold);
  for(int i = 0; i < 2*nfold; i++){
    newEdge[i] = new Edge(cellcenter, facecenter[i], getLevel()+1);
    edge_list.push_back(newEdge[i]);
  }

  //define new face
  vector<Face*> newFace(4*nfold);
  for(int i = 0; i< nfold; i++){
    newFace[i] = new Face(4);
    
    face_list.push_back(newFace[i]);
    
    //cellcenter->facecenter[i]->edgecenter[i]->facecenter[i-1]    
    newFace[i]->edge[0] = newEdge[i];
    newFace[i]->needReverse[0] = false;
    if(faceOrient[i] == 1){
      newFace[i]->edge[1] = face[i]->child[0]->edge[1];
      newFace[i]->needReverse[1] = false; //from facecenter to edgecenter
    }
    else{
      newFace[i]->edge[1] = face[i]->child[0]->edge[2];
      newFace[i]->needReverse[1] = false;
    }
    
    if(faceOrient[i==0?(nfold-1):(i-1)]==1){
      newFace[i]->edge[2] = face[i==0?(nfold-1):(i-1)]->child[0]->edge[2];
      newFace[i]->needReverse[2] = true;
    }
    else{
      newFace[i]->edge[2] = face[i==0?(nfold-1):(i-1)]->child[0]->edge[1];
      newFace[i]->needReverse[2] = true;
    }
    newFace[i]->edge[3] = newEdge[i==0?(nfold-1):(i-1)];
    newFace[i]->needReverse[3] = true;
    //face i+2*nfold
    newFace[i+2*nfold] = new Face(4);
    face_list.push_back(newFace[i+2*nfold]);
    
    //cellcenter->facecenter[i]->edgecenter[i+2*nfold]->facecenter[i+nfold]
    newFace[i+2*nfold]->edge[0] = newEdge[i];
    newFace[i+2*nfold]->needReverse[0] = false;
    if(faceOrient[i] == 1){
      newFace[i+2*nfold]->edge[1] = face[i]->child[2]->edge[2];
      newFace[i+2*nfold]->needReverse[1] = false; //from facecenter to edgecenter
    }
    else{
      newFace[i+2*nfold]->edge[1] = face[i]->child[2]->edge[1];
      newFace[i+2*nfold]->needReverse[1] = false;
    }
    if(faceOrient[i+nfold]==1){
      newFace[i+2*nfold]->edge[2] = face[i+nfold]->child[2]->edge[2];
      newFace[i+2*nfold]->needReverse[2] = true;
    }
    else{
      newFace[i+2*nfold]->edge[2] = face[i+nfold]->child[2]->edge[1];
      newFace[i+2*nfold]->needReverse[2] = true;
    }
    newFace[i+2*nfold]->edge[3] = newEdge[i+nfold];
    newFace[i+2*nfold]->needReverse[3] = true;
  }

  for(int i = nfold; i<2*nfold; i++){
    newFace[i] = new Face(4);
    face_list.push_back(newFace[i]);
    
    //cellcenter->facecenter[i]->edgecenter[i]->facecenter[i+1]    
    newFace[i]->edge[0] = newEdge[i];
    newFace[i]->needReverse[0] = false;
    if(faceOrient[i] == 1){
      newFace[i]->edge[1] = face[i]->child[0]->edge[1];
      newFace[i]->needReverse[1] = false; //from facecenter to edgecenter
    }
    else{
      newFace[i]->edge[1] = face[i]->child[0]->edge[2];
      newFace[i]->needReverse[1] = false;
    }
    if(faceOrient[i==(2*nfold-1)?nfold:(i+1)]==1){
      newFace[i]->edge[2] = face[i==(2*nfold-1)?nfold:(i+1)]->child[0]->edge[2];
      newFace[i]->needReverse[2] = true;
    }
    else{
      newFace[i]->edge[2] = face[i==(2*nfold-1)?nfold:(i+1)]->child[0]->edge[1];
      newFace[i]->needReverse[2] = true;
    }
    newFace[i]->edge[3] = newEdge[i==(2*nfold-1)? nfold:(i+1)];
    newFace[i]->needReverse[3] = true;
  
  
    //face i+2*nfold
  
    newFace[i+2*nfold] = new Face(4);
    face_list.push_back(newFace[i+2*nfold]);
    
    //cellcenter->facecenter[i]->edgecenter[i+2*nfold]->facecenter[i-nfold-1]
    
    newFace[i+2*nfold]->edge[0] = newEdge[i];
    newFace[i+2*nfold]->needReverse[0] = false;
    if(faceOrient[i] == 1){
      newFace[i+2*nfold]->edge[1] = face[i]->child[2]->edge[1];
      newFace[i+2*nfold]->needReverse[1] = false; //from facecenter to edgecenter
    }
    else{
      newFace[i+2*nfold]->edge[1] = face[i]->child[2]->edge[2];
      newFace[i+2*nfold]->needReverse[1] = false;
    }
    if(faceOrient[i==nfold?(nfold-1):(i-nfold-1)]==1){
      newFace[i+2*nfold]->edge[2] = face[i==nfold?(nfold-1):(i-nfold-1)]->child[2]->edge[1];
      newFace[i+2*nfold]->needReverse[2] = true;
    }
    else{
      newFace[i+2*nfold]->edge[2] = face[i==nfold?(nfold-1):(i-nfold-1)]->child[2]->edge[2];
      newFace[i+2*nfold]->needReverse[2] = true;
    }
    newFace[i+2*nfold]->edge[3] = newEdge[i==nfold? (nfold-1):(i-nfold-1)];
    newFace[i+2*nfold]->needReverse[3] = true;
  }
  //finish define newFace
  
  //allocate childCell
  childCell = new DiamondCell*[2*nfold +2];
  
  //child 0 and 1 is n-fold
  childCell[0]  = new DiamondCell(nfold);
  childCell[1] = new DiamondCell(nfold);
  //other children are 3-fold
  for(int i = 2; i < 2*nfold+2; i++){
    childCell[i] = new DiamondCell(3);
  }
  
  
  //set up level and parentCell 
  for(int i = 0; i < 2*nfold+2; i++){
    childCell[i] -> parentCell = this;
    childCell[i] ->whichChild  = i;
  }
  
  //define face and faceOrient

  for(int i = 0; i < nfold; i++){
    childCell[0]->face[i] = newFace[nfold-i-1];
    childCell[0]->faceOrient[i] = 1;//inward
  }
  for(int i = nfold; i< 2*nfold; i++){
    childCell[0]->face[i] = face[2*nfold-i-1]->child[0];
    childCell[0]->faceOrient[i] = faceOrient[2*nfold-i-1];
  }//finish childCell 0

 
  for(int i = 0; i < nfold; i++){
    childCell[1]->face[i] = newFace[i+nfold];
    childCell[1]->faceOrient[i] = 1; //inward
  }
  for(int i = nfold; i < 2*nfold; i++){
    childCell[1]->face[i] = face[i]->child[0];
    childCell[1]->faceOrient[i] = faceOrient[i];
  }//finish childCell 1

  for(int childID = 2; childID < 2+nfold; childID++){
   
    childCell[childID]->face[0] = newFace[2*nfold+childID-2];
    childCell[childID]->faceOrient[0] = 1;
    
    childCell[childID]->face[1] = newFace[3*nfold+childID-2];
    childCell[childID]->faceOrient[1] = 1;
    
    childCell[childID]->face[2] = newFace[childID-2];
    childCell[childID]->faceOrient[2] = 0;
    
    childCell[childID]->face[3] = face[childID-2]->child[faceOrient[childID-2]==0?3:1];
    childCell[childID]->faceOrient[3] = faceOrient[childID-2];
    
    childCell[childID]->face[4] = face[nfold+childID-2]->child[2];
    childCell[childID]->faceOrient[4] = faceOrient[nfold+childID-2];
    
    childCell[childID]->face[5] = face[childID==2?(nfold-1):(childID-3)]->
      child[faceOrient[childID==2?(nfold-1):(childID-3)]==0?1:3];
    childCell[childID]->faceOrient[5] = faceOrient[childID==2?(nfold-1):(childID-3)];
  }//finish childCell [2, nfold+2)

  for(int childID = nfold+2; childID < 2*nfold +2; childID++){
   
    childCell[childID]->face[0] = newFace[childID==(2*nfold+1)?(3*nfold):(2*nfold+childID-1)];
    childCell[childID]->faceOrient[0] = 0;
    
    childCell[childID]->face[1] = newFace[childID-2];
    childCell[childID]->faceOrient[1] = 0;
    
    childCell[childID]->face[2] = newFace[nfold+childID-2];
    childCell[childID]->faceOrient[2] = 0;
    
    childCell[childID]->face[3] = face[childID-2-nfold]->child[2];
    childCell[childID]->faceOrient[3] = faceOrient[childID-2-nfold];
    
    childCell[childID]->face[4] = face[childID==(2*nfold+1)?nfold:(childID-1)]->
      child[faceOrient[childID==(2*nfold+1)?nfold:(childID-1)]==0?1:3];
    childCell[childID]->faceOrient[4] = faceOrient[childID==(2*nfold+1)?nfold:(childID-1)];
    
    childCell[childID]->face[5] = face[childID-2]->child[faceOrient[childID-2]==0?3:1];
    childCell[childID]->faceOrient[5] = faceOrient[childID-2]; 
  }//finish all faces
  
}
  
void DiamondCell::empty_split(){
  if(childCell!=0)return;
  //allocate childCell
  childCell = new DiamondCell*[2*nfold +2];
  
  //child 0 and 1 is n-fold
  childCell[0]  = new DiamondCell(nfold);
  childCell[1] = new DiamondCell(nfold);
  //other children are 3-fold
  for(int i = 2; i < 2*nfold+2; i++){
    childCell[i] = new DiamondCell(3);
  }
};

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wchar-subscripts"
#endif

//this function return my sibling neighbor, mf is my faceID,
//nf is neighbor's faceId 
DiamondCell* DiamondCell::getSiblingNeib(const Cell* aCell,
                                         const std::vector<std::vector<Edge*> >& n2e,
                                         int mf, int& nf)const{
  
  //only face [0, nfold) has sibling neib 
  if(mf<0 || mf >= nfold){
    cerr <<"WARNING: faceID out of range" << endl;
    Loci::Abort();
  }
  
  char edgeID = -1;
  char neibID = -1;
  //if the sibling is inside a diamondcell,
  //define the sibling according to the numbering system  
  if(parentCell != 0){
    char n = parentCell->nfold;
    
    if(whichChild == 0){
      edgeID = n - mf -1;
      neibID = edgeID + 2;
      nf = 2;
    }
    else if(whichChild == 1){
      edgeID = n + mf;
      neibID = edgeID + 2;
      nf = 1;
    }
    
    else if (whichChild >= 2 && whichChild < n+2){
      switch(mf){
      case 0:
        // edgeID = 2*n +whichChild -2;
        neibID = whichChild + n;
        nf = 2;
        break;
      case 1:
        // edgeID = 3*n +whichChild -2;
        neibID = whichChild==2?(2*n+1):(whichChild +n -1);
        nf = 0;
        break;
      case 2:
        // edgeID = whichChild -2;
        neibID = 0;
        nf = n+1-whichChild;
        break;
      default:
        cerr<<" WARNING: illegal mf" << endl;
        Loci::Abort();
        break;
      }
    }
    else if( whichChild >= (n+2) && whichChild < (2*n+2)){
      switch(mf){
      case 0:
        
        neibID = (whichChild==(2*n+1))?2:(whichChild-n+1);
        nf = 1;
        break;
      case 1:
        
        neibID = 1;
        nf = whichChild-n-2;;
        break;
      case 2:
        
        neibID = whichChild -n;
        nf = 0;
        break;
      default:
        cerr<<" WARNING: illegal mf" << endl;
        Loci::Abort();
        break;
      }
    }  
    
    else{
      cerr<< " WARNING: illegal faceID" << endl;
      Loci::Abort();
    }
    //sanity check
    if(this->face[mf] != parentCell->childCell[neibID]->face[nf]){
      cerr<<"WARNING: find the wrong sibling Neib" << endl;
      Loci::Abort();
    }
    return parentCell->childCell[neibID];
  }
  //else find sibling in aCell
  else if(parentCell == 0){
    Edge* theEdge =  n2e[whichChild][mf];
    int childID;
    
    for( childID = 0; childID < aCell->numNode; childID++){
      if(childID != whichChild){
        for(nf = 0; nf < int(n2e[childID].size()); nf++){
          if(n2e[childID][nf] == theEdge){
            break;
          }
        }
        if(nf != int(n2e[childID].size()))break;
      }
    }
    
    if(childID >= aCell->numNode){
      cerr<<"WARNING: can not find sibling neib in diamonds" << endl;
      return 0;
    }
    //sanity check
    if(this->face[mf] != aCell->child[childID]->face[nf]){
      cerr << " WARNING: find the wrong sibling neib" << endl;
      Loci::Abort();
    }
    return aCell->child[childID];
    
  }
  return 0;
}




int DiamondCell::parentFace( int faceID)const{
  if(faceID < nfold || faceID >= 2*nfold){
    cerr<<"WARNING: no parent faceID" << endl;
    return -1;
  }
  
  for(int i = 0; i <2*(parentCell->nfold); i++){
    int childID  = parentCell->face[i]->containFace(face[faceID]);
    if(childID >=0){
      return i;
    }
  }
  cerr<<"WARNING: can not find parentFace"<< endl;
  Loci::Abort();
  return -1; //dummy code
}

int DiamondCell::get_num_fine_faces(){
  int count = 0;
  for(int i = 0; i < 2*nfold; i++){
    count += face[i]->get_num_leaves();
  }
  return count;
}


//this function find face neighbor,

//condition: the neighbor's face is greater or same as mine
//so only one cell is returned

//condition on ff: in, my faceID, ff >=0 && ff < 2*nfold
DiamondCell*  DiamondCell::findNeighbor(const Cell* aCell,
                                        const std::vector<std::vector<Edge*> >& n2e,
                                        int ff, int& nf)const{
  //security check
  if( ff < 0 || ff >=2*nfold){
    cerr << "WARNING: illegal faceID in findNeighbor()" << endl;
    exit(1);
  }
  
  //find sibling neighbor(real sibling)
  if(ff<nfold){
    return getSiblingNeib(aCell, n2e,ff, nf);    
  }
  //subroot has no neib
  if(parentCell == 0 && ff >= nfold){
    //if no face neib inside a genaral cell, nf is meaningless 
    return 0;
  }
  
   
  //otherwise find the nearest common ancestor,
  else{//ff>=nfold

    std::stack<Face*> faceStack;
    const  DiamondCell* current = this;
    int   currentFace = ff;
    DiamondCell* NCA;
    
    while(true){
      
      if(currentFace < current->nfold){
        NCA = current->getSiblingNeib(aCell, n2e, currentFace, nf);
        break;
      }
    
      if((current->parentCell == 0) && (currentFace >= current->nfold)) return 0; 
      
      //current->parentCell != 0 && currentFace >= current->nfold
      faceStack.push(current->face[currentFace]);
      int newFaceID = current->parentFace(currentFace);
      if(newFaceID < 0) {
        cerr << "WARNING: illegal faceID" << endl;
        Loci::Abort();
      }
      
      current = current->parentCell;
      currentFace = newFaceID;
    }
    if(NCA == 0){
      return 0;  //no neib
    }
    
    DiamondCell* pNeighbor = NCA;
    
    //if pNeighbor is a leaf, or it's not bigger than me, return it
    if(pNeighbor->childCell == 0 ||pNeighbor->getLevel() >= getLevel() ) return pNeighbor;
    
    
    //if pNeighbor is not a leaf, and her face still greater than mine,
    //climb down to my level
    bool childFound;
    int cID = 0;
    while(pNeighbor->getLevel() < getLevel() && pNeighbor->childCell != 0){
      childFound = false;
      for(cID = 0 ; cID < 2*(pNeighbor->nfold)+2; cID++){
        nf = pNeighbor->childCell[cID]->containFace(faceStack.top());
        if(nf>=0){
          childFound = true;
          break;
        }
      }
      if(!childFound){
        cerr <<"WARNING: get lost when climbing down the tree" << endl;
        Loci::Abort();
      } 
      
      pNeighbor = pNeighbor->childCell[cID];
      faceStack.pop();
    }
    
    //sanity check
    if(this->getLevel() == pNeighbor->getLevel() && this->face[ff] != pNeighbor->face[nf]){
      cerr << " WARNING: get the wrong neighbor" << endl;
      Loci::Abort();
    }
    return pNeighbor;
  }
  // return 0;
}


//this function splits DiamondCell isotropically level times
// newly created nodes, edges and faces are put into the lists
void DiamondCell::resplit(int level, std::list<Node*>& node_list,
                          std::list<Edge*>& edge_list, std::list<Face*>& face_list ){
 
  if(level <= 0) return;

  int initialLevel  = this->getLevel();
  int finalLevel = initialLevel + level;
  std::queue<DiamondCell*> Q;
  Q.push(this);
  
  DiamondCell* current;
  
  
  while(!Q.empty()){
    current = Q.front();
    current->split(node_list, edge_list, face_list);
    if((current->getLevel()+1)< finalLevel){
      for(int i = 0; i <(2*current->nfold +2); i++){
        Q.push(current->childCell[i]);
      }
    }
  
    Q.pop();
  }
}

//get the 4*nfold edges of DiamondCell
void DiamondCell::get_edges(std::set<Edge*>& edges){
  for(int i = 0; i < 2*nfold; i++){
    for(int j = 0; j < face[i]->numEdge; j++){
      edges.insert(face[i]->edge[j]);
    }
  }
  if(int(edges.size()) != 4*nfold){
    cerr<< "WARNING: incorrect num of edges" << endl;
    cerr << "edges_size " <<edges.size() << " 4*nfold " << 4*nfold << endl;
    Loci::Abort();
  }
}

//get 2*nfold+2 nodes
void DiamondCell::get_nodes(std::set<Node*>& nodes){
  std::set<Edge*> edges;
  get_edges(edges);
  for(std::set<Edge*>::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    nodes.insert((*ep)->head);
    nodes.insert((*ep)->tail);
  }
  if(int(nodes.size()) != (2*nfold+2)){
    cerr<< " WARNING: incorrect num of nodes" << endl;
    Loci::Abort();
  }
}


//if all nodes get detagged, the cell get detagged, return 2
//if any node get tagged, the cell get tagged, return 1
//other wise return 0
int DiamondCell::get_tagged(){
  std::set<Node*> nodes;
  get_nodes(nodes);
  //if all nodes get detagged, the cell is detagged
  bool detagged = true;
  for(std::set<Node*>::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    if((*np)->tag != 2)detagged = false;
  }
  if(detagged) return 2;
  //if any node gets tagged, the cell is tagged
  for(std::set<Node*>::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    if((*np)->tag == 1)return 1;
  }
  //otherwise, the cell remains unchanged
  return 0;
}

int DiamondCell::get_tagged(const vector<source_par>& sources){
  std::set<Node*> node_set;
  get_nodes(node_set);
  
  std::vector<Node*> nodes(node_set.size());
  int vi = 0;
  for(std::set<Node*>::const_iterator ni = node_set.begin(); ni != node_set.end();
      ni++, vi++){
    nodes[vi] = *ni;
  }
  double min_len = get_min_edge_length();
  return tag_cell(nodes, sources, min_len);
}
void Cell::resplit( const std::vector<char>& cellPlan,
                    std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<Face*>& face_list,
                    std::vector<DiamondCell*>& cells){

  int cellID =0;  
  if(cellPlan.size() == 0){
    reduce_vector(cells);
    return;
  }
  char currentCode = cellPlan[0];
  if(currentCode == 0){
    reduce_vector(cells);
    return;
  }
  
  std::queue<DiamondCell*> Q;
  //assume if cellPlan not empty, it always begin with 1
  if(child ==0)split(node_list, edge_list, face_list);
  
  for(int i = 0; i < numNode; i++){ 
    Q.push(child[i]);
  }
  
  //process the DiamondCells in the Q until Q is empty
  
  DiamondCell* current;
  unsigned int index =1;
 
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
             
    switch(currentCode)
      {
        
        //0 no split,this is a leaf, output cells
      case 0:
        current->cellIndex = ++cellID;
        cells.push_back(current);
        break;
         
      case 1:
         
        if(current->childCell==0)current->split(node_list, edge_list, face_list);
         
        for(int i = 0; i < 2*(current->nfold) +2; i++){
          Q.push(current->childCell[i]);
        }
         
        break;
        
      default:
        cerr <<"WARNING: illegal splitcode in function Cell::reSplit()" << endl;
        break;
      }
    
    Q.pop();
  }
  reduce_vector(cells);
}






//this function check each cell in cells, find all the individual faces,
//put their connectivity info into faces. 
// cells: in, all the leaves cells in the  tree
// faces: out, initially empty, return all the indivial faces and its two cell inex

void set_general_faces(const Cell* aCell,
                       const std::vector<DiamondCell*>& cells,
                       const std::vector<std::vector<Edge*> >& n2e,
                       std::list<pair<Face*, NeibIndex> >& faces){
  
  DiamondCell* tempNeib = 0;
  int numCells = cells.size();
  
  if(numCells > 0){
    //first initialize faceMarked
    for(int i = 0; i < numCells; i++){
      cells[i] -> faceMarked = new bool[2*cells[i]->nfold];
      for(int j= 0; j < 2*cells[i]->nfold; j++){
        cells[i]->faceMarked[j] = false;
      }
    }

    
    //for each cell
    for(int i = 0; i < numCells; i++){
      
      int numFaces = 2* cells[i]->nfold;
      //for each face
      for(int ff = 0; ff <numFaces; ff++){
        //if not marked
        if(!(cells[i]->faceMarked[ff])){
          //mark it
          cells[i]->faceMarked[ff] = true;
          
          
          //find face neigbor
          int nf = -1;
          tempNeib = cells[i]->findNeighbor(aCell,n2e,ff, nf);
          
          //if no neib, don't output the face
          //there is a neib
          if(tempNeib != 0){
            bool isALeaf =  tempNeib->childCell == 0;
            //if neib is a leaf, marked its face
            if(isALeaf){
              tempNeib ->faceMarked[nf] = true;
              //if face ff point outward cells[i]
              if(cells[i]->faceOrient[ff] == 0)faces.push_back(make_pair(cells[i]->face[ff], NeibIndex(cells[i]->cellIndex, tempNeib ->cellIndex)));
              else faces.push_back(make_pair(cells[i]->face[ff], NeibIndex( tempNeib ->cellIndex,cells[i]->cellIndex)));
            }
            //if myself > neighbor do nothing
            
          }
        }
      }
      
    }
    
  }
}
 
 


//this function balance a leaf cell current_cell

bool DiamondCell::balance_cell(std::list<Node*>& node_list,
                               std::list<Edge*>& edge_list,
                               std::list<Face*>& face_list){ 
  bool  needBalance = false;
  
  if(childCell!=0){
    for(int i = 0; i < 2*nfold+2; i++){
      bool tmp = childCell[i]->balance_cell(node_list, edge_list, face_list);
      needBalance = needBalance ||tmp; 
    }
    return needBalance;
  }
  else{
 

    if(childCell ==0){  
      //gather all the edges of this
      std::set<Edge*> edges;
      get_edges(edges);
      needBalance = false;
      for(std::set<Edge*>::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
       
        if( (*ep)->depth_greater_than_1()){
          needBalance = true;
          split(node_list, edge_list, face_list);
         
          for(int i = 0; i < (2*nfold+2); i++){
            childCell[i]->balance_cell(node_list, edge_list, face_list);
          }
          break;
        }
       
      }
    }
    //if all the edges are balanced, check the faces
    if(Globals::balance_option >=1 && childCell == 0){
      int num_faces_split = 0;
      for(int i = 0; i < 2*nfold; i++){
        if(face[i]-> child != 0) num_faces_split++;
      }
      if(num_faces_split > nfold){
        needBalance = true;
        split(node_list, edge_list, face_list);
        for(int i = 0; i < (2*nfold+2); i++){
          childCell[i]->balance_cell(node_list, edge_list, face_list);
        } 
      }
    }

    if(Globals::balance_option >=2 && childCell == 0){
      for(int i = 0; i < nfold; i++){
        if(face[i]-> child != 0 && face[(i+2)%nfold+nfold]->child != 0){
          needBalance = true;
          split(node_list, edge_list, face_list);
          for(int i = 0; i < (2*nfold+2); i++){
            childCell[i]->balance_cell(node_list, edge_list, face_list);
          }
          break;
        }
      }
    }
   
    return needBalance;
  
  }
}

bool Cell::balance_cell(std::list<Node*>& node_list,
                        std::list<Edge*>& edge_list,
                        std::list<Face*>& face_list){ 
  bool  needBalance = false;
  if(child != 0){
    needBalance = false;
    std::list<DiamondCell*> leaves;
    for(int i = 0; i < numNode; i++){
      child[i]->sort_leaves(leaves);
    }
    for(std::list<DiamondCell*>::const_iterator p = leaves.begin(); p != leaves.end(); p++){
      bool tmp =  (*p)->balance_cell( node_list, edge_list, face_list);
      needBalance = tmp||needBalance;
    }
    return needBalance;  
  }
 
  if(child == 0){
   
    needBalance = false;
    for(int i = 0; i < numEdge; i++){
      if( edge[i]->depth_greater_than_1()){
        needBalance = true;
        split(node_list, edge_list, face_list);
        
        for(int i = 0; i < numNode; i++){
          child[i]->balance_cell(node_list, edge_list, face_list);
        }
        break;
      }
    }
   
  }

  //if all the edges are balanced, check the faces
  if(Globals::balance_option >=1){
    if(child ==0){
      int num_faces_split = 0;
      for(int i = 0; i < numFace; i++){
        if(face[i]-> child != 0) num_faces_split++;
      }
      if(num_faces_split > numFace/2){
        needBalance = true;
        split(node_list, edge_list, face_list);
       
        for(int i = 0; i < numNode; i++){
          child[i]->balance_cell(node_list, edge_list, face_list);
        } 
      }

    }
  }

  

 
  return needBalance;
}
//depth first sort_leaves
void DiamondCell::sort_leaves(std::list<DiamondCell*>& leaves){
  if(childCell != 0){
    for(int i = 0; i <(2*nfold+2); i++)childCell[i]->sort_leaves(leaves);
  }
  else{
    leaves.push_back(this);
  }
}

// //breadth first get leaves
// void Cell::get_leaves(std::vector<DiamondCell*>& leaves){
//   std::queue<DiamondCell*> Q;

//   //no leaves
//   if(child ==0)return;
  
//   for(int i = 0; i < numNode; i++){ 
//     Q.push(child[i]);
//   }
  
//   //process the DiamondCells in the Q until Q is empty
  
//   DiamondCell* current;
   
//   while(!Q.empty()){
//     current = Q.front();
//     if(current->childCell==0)leaves.push_back(current);
//     else{
//       for(int i = 0; i < 2*(current->nfold) +2; i++){
//         Q.push(current->childCell[i]);
//       }
//     }
    
//     Q.pop();
//   }
// }



void Cell::rebalance_cells(std::list<Node*>& node_list,
                           std::list<Edge*>& edge_list,
                           std::list<Face*>& face_list){ 
  bool need_balance_more = true;

  while(need_balance_more){


    need_balance_more = balance_cell(node_list, edge_list, face_list); 
  }
}

void Cell::sort_leaves(std::list<DiamondCell*>& leaves){
  if(child !=0){
    for(int i = 0; i <numNode; i++)child[i]->sort_leaves(leaves);
  }
}
int Cell::get_tagged(){
  //if all nodes get detagged, the cell is detagged
  bool detagged = true;
  for(int i = 0; i<numNode; i++){
    if(node[i]->tag != 2)detagged = false;
  }
  if(detagged) return 2;
  //if any node gets tagged, the cell is tagged
  for(int i = 0; i<numNode; i++){
    if(node[i]->tag ==1) return 1;
  }
  //otherwise remain unchanged
  return 0;
}

int Cell::get_tagged(const vector<source_par>& sources){
  vector<Node*> nodes(numNode);
  for(int i = 0; i<numNode; i++){
    nodes[i] = node[i];
  }
  double min_len = get_min_edge_length();
  return tag_cell(nodes, sources, min_len);  
}


int Cell::get_num_fine_faces(){
  int count = 0;
  for(int i = 0; i < numFace; i++){
    count += face[i]->get_num_leaves();
  }
  return count;
}
void Cell::set_n2f_n2e(std::vector<Face*>& n2f, std::vector<Edge*>& n2e, std::vector<int>& rot, int nindex){

  n2e.clear();
  n2f.clear();
  rot.clear();
  
  for( int i=0; i< numFace; i++){
    
    int j = face[i]->containNode(node[nindex]);
    if(j>=0){
      n2f.push_back(face[i]);
      if(faceOrient[i] == 0) rot.push_back(j);
      else rot.push_back(-(j+1));
    }
  }
 
  //rearrange n2f, so that n2f[i] and n2f[i+1] will share n2e[i]
  //procedure is in zig-zag way, n2f[i]-> n2e[i]->n2f[i+1],  
  for(unsigned int i = 0; i < n2f.size(); i++){
    if(rot[i]>=0) n2e.push_back(n2f[i]->edge[rot[i]==0?(n2f[i]->numEdge - 1):(rot[i]-1)]);
    else n2e.push_back(n2f[i]->edge[-rot[i]-1]);
    if(i != n2f.size() -1){
      for(unsigned int f =i+1; f< n2f.size(); f++){
        int e;
        for( e=0; e < n2f[f]->numEdge; e++){
          if(n2f[f]->edge[e] == n2e.back()){
            if(f != i+1){
              std::swap(n2f[f], n2f[i+1]);
              std::swap(rot[f], rot[i+1]);
            }
            break;
          }
          
        }
        if(e < n2f[f]->numEdge) break; 
      }
    }
  }
  
  if(n2f.size() != n2e.size()){
    cerr << " WARNING: numNodes and numEdges not equal in a face" << endl;
    Loci::Abort();
  }
 

  
  reduce_vector(n2f);
  reduce_vector(n2e);
  reduce_vector(rot);
 
}


std::vector<std::vector<Edge*> > Cell::set_n2e(){
  std::vector<std::vector<Edge*> >  n2e(numNode);
  std::vector<Face*> n2f;
  std::vector<int> rot;
  
  for(int nindex = 0; nindex < numNode; nindex++){
    n2f.clear();
    rot.clear();
    
    for( int i=0; i< numFace; i++){
      int j = face[i]->containNode(node[nindex]);
      if(j>=0){
        n2f.push_back(face[i]);
        if(faceOrient[i] == 0) rot.push_back(j);
        else rot.push_back(-(j+1));
      }
    }
    //rearrange n2f, so that n2f[i] and n2f[i+1] will share n2e[i]
    //procedure is in zig-zag way, n2f[i]-> n2e[i]->n2f[i+1],  
    for(unsigned int i = 0; i < n2f.size(); i++){
      if(rot[i]>=0) n2e[nindex].push_back((n2f[i])->edge[(rot[i]==0)?(n2f[i]->numEdge - 1):(rot[i]-1)]);
      else (n2e[nindex]).push_back(n2f[i]->edge[-rot[i]-1]);
      if(i != n2f.size() -1){
        for(unsigned int f =i+1; f< n2f.size(); f++){
          int e;
          for( e=0; e < n2f[f]->numEdge; e++){
            if(n2f[f]->edge[e] == (n2e[nindex]).back()){
              if(f != i+1){
                std::swap(n2f[f], n2f[i+1]);
                std::swap(rot[f], rot[i+1]);
              }
              break;
            }
            
          }
          if(e < n2f[f]->numEdge) break; 
        }
      }
    }
    
    if(n2f.size() != n2e[nindex].size()){
      cerr << " WARNING: numFaces and numEdges not equal in a cell" << endl;
      Loci::Abort();
    }
   
    reduce_vector(n2e[nindex]);
  
  }
  return n2e; 
}

  

//this function split  a general cel,
void Cell::split(std::list<Node*>& node_list, std::list<Edge*>& edge_list, std::list<Face*>& face_list){
  if(child!=0)return;
  //split each face
  for(int i = 0; i <numFace; i++){
    if(face[i]->child==0)face[i]->split(node_list, edge_list);
  }
  
  
  //calculate cellcenter and put it into node_list
  Node* cellcenter = centroid();
  node_list.push_back(cellcenter);

  
  //get facecenter
  std::vector<Node*> facecenter(numFace);
  getFaceCenter(&facecenter[0]);
  
  //define new edges , each edge points from cellcenter to facecenter
  std::vector<Edge*> newEdges(numFace);
  for(int i = 0; i< numFace; i++){
    newEdges[i] = new Edge(cellcenter, facecenter[i], 1);
    edge_list.push_back(newEdges[i]);
  }
  
  //define new faces 
  //each edge generates a new quad face, 
  //for each new quad face, 
  std::vector<Face*> newFaces(numEdge);
  for(int i = 0; i < numEdge; i++){
    newFaces[i] = new Face(4);
   
   
    
    //find two faces that share the edge ,
    int findex1=0, findex2 = 0;
    int edgeID1 = 0, edgeID2 = 0;

    for(findex1 = 0; findex1 < numFace; findex1++){
      edgeID1= face[findex1]->containEdge(edge[i]);
      if(edgeID1 >= 0)break;
    }
   
    for(findex2 = (findex1+1); findex2 < numFace; findex2++){
      edgeID2 = face[findex2]->containEdge(edge[i]);
      if(edgeID2 >= 0)break;
    }
    if(findex2 == numFace){
      cerr<< "WARNING: findex out of range" << endl;
      Loci::Abort();
    }
    //findex1 < findex2
    
    //edge[0] points from cellcenter to facecenter[findex1]
    newFaces[i]->edge[0] = newEdges[findex1];
    newFaces[i]->needReverse[0] = false;
    
    //edge[3] points from  cellcenter to facecenter[findex2]
    newFaces[i]->edge[3] = newEdges[findex2];
    newFaces[i]->needReverse[3] = true;
    
    newFaces[i]->edge[1] = face[findex1]->child[edgeID1]->edge[1];
    newFaces[i]->needReverse[1] = false;
    newFaces[i]->edge[2] = face[findex2]->child[edgeID2]->edge[1];
    newFaces[i]->needReverse[2] = true;
    face_list.push_back(newFaces[i]); 
  }
    
  //for each node, define the child cell containing it
  child = new DiamondCell*[numNode];
  for(int  nindex = 0; nindex < numNode; nindex++){
    //find all the faces that connected to the node, put them into n2f
    std::vector<Face*> n2f;
    std::vector<Edge*> n2e;
    std::vector<int> rot;
    set_n2f_n2e(n2f, n2e,rot, nindex);
    //define the child cell
    child[nindex] = new DiamondCell(n2e.size());
    child[nindex]->parentCell = 0;
    child[nindex]->whichChild = nindex;
   

    //define face [0, nfold],
    //each edge in n2e generates a new face
    for(unsigned int i = 0; i < n2e.size(); i++){
      
      //for each edge in n2e, find eindex
      int eindex = 0;
      for(eindex = 0; eindex < numEdge; eindex++){
        if(edge[eindex] == n2e[i])break;
      }
      
      child[nindex]->face[i] = newFaces[eindex];

      //if findex1(n2f[i]) < findex2(n2f[i+1] point inward, otherwise point outward
      int findex1, findex2;
      for(findex1 = 0; findex1 < numFace; findex1++){
        if(face[findex1] == n2f[i]) break;
      }

      for(findex2 = 0; findex2 < numFace; findex2++){
        if(face[findex2] == n2f[i==(n2f.size()-1)?0:i+1])break;
      }
      
      if(findex1 < findex2) child[nindex]->faceOrient[i] = 1;
      else child[nindex]->faceOrient[i] = 0;
    }
    //define face [nfold, 2*nfold)
    for(unsigned int i = n2e.size(); i < 2*n2e.size(); i++){
      
      child[nindex]->face[i] = n2f[i-n2e.size()]->child[rot[i-n2e.size()]>=0?
                                                        rot[i-n2e.size()]:(-rot[i-n2e.size()]-1)];
      child[nindex]->faceOrient[i] = rot[i-n2e.size()]>=0?0:1;
    }
  }
 
}

std::vector<char> Cell::make_cellplan(){
  
  std::vector<char> cellPlan;
  if(child == 0) {
    reduce_vector(cellPlan);
    return cellPlan;
  }
  cellPlan.push_back(1);
  std::queue<DiamondCell*> Q;
  for(int i = 0; i < numNode; i++){
    Q.push(child[i]);
  }
  DiamondCell* current;
  while(!Q.empty()){
    current = Q.front();
    if(current->getChildCell() != 0){
      cellPlan.push_back(1);
      for(int i = 0; i < 2*current->getNfold()+2; i++){
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

void Cell::empty_split(){
  if(child!=0)return;
  child = new DiamondCell*[numNode];
  for(int  nindex = 0; nindex < numNode; nindex++){
    //find all the faces that connected to the node, put them into n2f
    std::vector<Face*> n2f;
    std::vector<Edge*> n2e;
    std::vector<int> rot;
    set_n2f_n2e(n2f, n2e,rot, nindex);
    //define the child cell
    child[nindex] = new DiamondCell(n2e.size());
  }
 
}
int32 Cell::empty_resplit(const std::vector<char>& cellPlan){
  
                     
  if(cellPlan.size() == 0){
    
    return 1;
  }
  std::queue<DiamondCell*> Q;
  //assume if cellPlan not empty, it always begin with 1
  if(child==0)empty_split();
  
  for(int i = 0; i < numNode; i++){ 
    Q.push(child[i]);
  }
  
  //process the DiamondCells in the Q until Q is empty
  int cIndex = 0; //the local index starts with 1 
  DiamondCell* current;
  unsigned int index =1;
  char currentCode;
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
      
    switch(currentCode)
      {
          
        //0 no split,this is a leaf, output cells
      case 0:
        current->cellIndex = ++cIndex;
        break;
         
      case 1:
          
        if((current->childCell)==0)current->empty_split();
         
        for(int i = 0; i < 2*(current->nfold) +2; i++){
          Q.push(current->childCell[i]);
        }
         
        break;
          
      default:
        cerr <<"WARNING: illegal splitcode in function Cell::empty_resplit(): " << currentCode+'0'<< endl;
        break;
      }
      
    Q.pop();
  }
  return cIndex;  
}



int32 Cell::traverse(const std::vector<char>& parentPlan,  vector<pair<int32, int32> >& indexMap){
  indexMap.clear();
  if(parentPlan.size() == 0){
    
    if(child==0){
      indexMap.push_back(make_pair(1,1));
      return 1;
    }else{
      
      list<DiamondCell*> leaves;
      for(int i = 0; i < numNode; i++){
        child[i]->sort_leaves(leaves);
      }
      for(std::list<DiamondCell*>::const_iterator p = leaves.begin(); p != leaves.end(); p++){
        indexMap.push_back(make_pair((*p)->cellIndex, 1));
       
      }
      return 1;
    }
  }
  
  std::queue<DiamondCell*> Q;
  //assume if parentPlan not empty, it always begin with 1
  for(int i = 0; i < numNode; i++){ 
    Q.push(child[i]);
  }
  //process the DiamondCells in the Q until Q is empty
  int32 cIndex = 0;  
  DiamondCell* current;
  unsigned int index =1;
  char currentCode;
  while(!Q.empty()){
    current = Q.front();
    //take a code from splitcode

    if(index >= parentPlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from splitcode
      currentCode = parentPlan[index];
      index++;  
    }
  
    list<DiamondCell*> leaves;
    switch(currentCode)
      {
        
        //0 no split,this is a leaf for parentPlan
      case 0:
        ++cIndex;
        current->sort_leaves(leaves);
        if(leaves.front() != current){//current is not a leaf for cellPlan
          for(std::list<DiamondCell*>::const_iterator p = leaves.begin(); p != leaves.end(); p++){
            indexMap.push_back(make_pair((*p)->cellIndex, cIndex));
          }
        }else{
          //current is a leaf for cellPlan and for parentPlan
          if(current->cellIndex != 0) indexMap.push_back(make_pair(current->cellIndex, cIndex));  
        }
        break;
         
      case 1:
        if(current->cellIndex !=0){//derefinement happen, current is a leaf for cellPlan
          current->sort_leaves(leaves);
          for(std::list<DiamondCell*>::const_iterator p = leaves.begin(); p != leaves.end(); p++){
            indexMap.push_back(make_pair(current->cellIndex, (*p)->cellIndex));
            (*p)->cellIndex = 0; //so that *p will not appear in indexMap when it is popped out of Q
          }
        }
        for(int i = 0; i < 2*(current->nfold) +2; i++){
          Q.push(current->childCell[i]);
        }
         
        break;
          
      default:
        cerr <<"WARNING: illegal splitcode in function Cell:traverse(): " <<currentCode+'0'<< endl;
        break;
      }
      
    Q.pop();
  }
  
  return cIndex;  
}

bool Cell::needDerefine(){
  if(child!= 0){
    bool derefine = true;
    for(int i = 0; i < numNode; i++){
      if(child[i] != 0){
	if((child[i]->get_tagged())!=2) return false;
	if(child[i]->childCell !=0) return false;
      }
    }
    if(derefine){
      for(int i = 0; i < numEdge; i++){
	if( edge[i]->depth_greater_than_1())return false;
      }
      return true;
    }
  }
  return false;
}

bool Cell::needDerefine_ctag(){
  if(child!= 0){
    bool derefine = true;
    for(int i = 0; i < numNode; i++){
      if(child[i] != 0){
	if((child[i]->getTag())!=2) return false;
	if(child[i]->childCell !=0) return false;
      }
    }
    if(derefine){
      for(int i = 0; i < numEdge; i++){
	if( edge[i]->depth_greater_than_1())return false;
      }
      return true;
    }
  }
  return false;
}


void Cell::derefine(){
  if(child!= 0){
    for(int i = 0; i < numNode; i++){
      if(child[i] != 0){
	delete child[i];
	child[i] = 0;
      }
    }
    delete [] child;
    child = 0;
    
  }
}

bool DiamondCell::needDerefine(){
  if(childCell != 0){
    bool derefine = true;
    for(int i = 0; i < 2*nfold +2; i++){
      if(childCell[i]!=0){
	if( (childCell[i] ->get_tagged()) != 2)return false;
	if(childCell[i]->childCell !=0) return false;
      }
    }
    if(derefine){
      std::set<Edge*> edges;
      get_edges(edges);
      for(std::set<Edge*>::const_iterator ei = edges.begin(); ei != edges.end(); ei++){
	if((*ei)->depth_greater_than_1())return false;
      }
      return true;
    }
  }
  return false;
}

bool DiamondCell::needDerefine_ctag(){
  if(childCell != 0){
    bool derefine = true;
    for(int i = 0; i < 2*nfold +2; i++){
      if(childCell[i]!=0){
	if( (childCell[i] ->getTag()) != 2)return false;
	if(childCell[i]->childCell !=0) return false;
      }
    }
    if(derefine){
      std::set<Edge*> edges;
      get_edges(edges);
      for(std::set<Edge*>::const_iterator ei = edges.begin(); ei != edges.end(); ei++){
	if((*ei)->depth_greater_than_1())return false;
      }
      return true;
    }
  }
  return false;
}

void DiamondCell::derefine(){
  if(childCell != 0){                   
    for(int i = 0; i < 2*nfold +2; i++){
      if(childCell[i] != 0){
	delete  childCell[i];
	childCell[i] = 0;
      }
    }
    delete [] childCell;
    childCell = 0;
    
  }
}
 
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
