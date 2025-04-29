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
#include "quadface.h"
#include "hex_defines.h"
using std::cout;
using std::endl;


//for get_c1_hex and get_c1_prism,when a cell is empty_split, faceMap is created
// and face is split to generate leaves, by look-up the faceMap, c1 can be computed
//use Range2d to avoid full split cell and face.  
std::vector<int32> contain_2d(const std::vector<pair<Range2d, int32> >& faceMap,
                              const std::vector<Range2d>& leaves){
  std::vector<int32> c1(leaves.size());
  for(unsigned int i = 0; i < leaves.size(); i++){
    unsigned int j = 0;
    for( j = 0; j < faceMap.size(); j++){
      Range2d bigface = faceMap[j].first;
      if(leaves[i].minP.x >= bigface.minP.x &&
         leaves[i].maxP.x <= bigface.maxP.x &&
         leaves[i].minP.y >= bigface.minP.y &&
         leaves[i].maxP.y <= bigface.maxP.y){
        c1[i] = faceMap[j].second;
        break;
      }
    }
    if(j == faceMap.size()){
      cerr<<"WARNING: can not find a Range that contain the face" << endl;
      Loci::Abort();
    }
  }
  return c1;
}

//f2c orient functions are uesd when a quadface is built as in cell
// and the facePlan is for the face defined by face2node 
char orient_splitCode_f2c(char splitCode, char orientCode){
  switch(orientCode){
  case 0:
  case 2:
  case 5:
  case 7:  
    return splitCode;
    break;
  case 1:
  case 3:
  case 4:
  case 6:
    switch(splitCode){
    case 1:
      return 2;
      break;
    case 2:
      return 1;
      break;
    case 0:
    case 3:
      return splitCode;
      break;
    default:
      cerr<<"WARNING: reach dummy code 1 in orient_splitCode_f2c" << endl;
      break;
    }
    break;
  default:
    cerr<<"WARNING: reach dummy code 2 in orient_splitCode_f2c" << endl;
    break;
  }
  cerr<<"WARNING: reach dummy code 3 in orient_splitCode_f2c" << endl;
  return 0;//dummy code
}

//f2c orient functions are uesd when a quadface is built as in cell
// and the facePlan is for the face defined by face2node 
//splitcode is for face
char orient_childID_f2c(char childID, char orientCode, char splitCode){
  Point2d p;
  Point2d pf;
  switch(splitCode){
  case 3:
    
    p.x = childID/2;
    p.y = childID%2;
    pf = transfer_f2c(p, Point2d(1,1),orientCode);
   
    return char(pf.x *2 + pf.y);
    break;
  case 2:
    switch(orientCode){
    case 0:
    case 3:
    case 4:
    case 7:
      return childID;
      break;
    case 1:
    case 2:
    case 5:
    case 6:
      return (1-childID);
      break;
    }
    break;
  case 1:
    switch(orientCode){
    case 0:
    case 1:
    case 4:
    case 5:
      return childID;
      break;
    case 2:
    case 3:
    case 6:
    case 7:
      return 1-childID;
      break;
    }
    break;
  default:
    cerr<< "WARNING: illegal splitcode, " << int(splitCode)<< " reach dummy code in orient_childID_f2c" <<endl;
    break;
  }
  cerr<<"WARNING: reach dummy code in orient_childID_f2c" << endl;
  return -1;
}
      
//c2f orient functions are used when the face is defined by face2node
//and facePlan is extracted from cell 
char orient_edgeID_c2f(char edgeID, char orientCode){
  if(orientCode/4 == 0){
    return char((edgeID +orientCode)%4);
  }
  else{
   
    return char((3-edgeID+orientCode)%4);
  }
}

//f2c orient functions are uesd when a quadface is built as in cell
// and the facePlan is for the face defined by face2node 
char orient_edgeID_f2c(char edgeID, char orientCode){
  if(orientCode/4 ==0){
    return char((edgeID + 4-orientCode)%4);
  }
  else{
   
    return char((3-edgeID+orientCode)%4);
  }
}

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wchar-subscripts"
#endif

//used in building cells, quadface is built as defined in cell, and split with orientCode
//all new nodes and edges are put into node_list and edge_list
void QuadFace::split(char splitCode, char orientCode, 
                     std::list<Node*>& node_list,
                     std::list<Edge*>& edge_list
                     ){
 
  Edge* newEdge = 0;
  Edge* tempEdge = 0;
  Node* facecenter = 0;
  Edge* newedge[4];
  Node* edgecenter[4];
  char formercode = this->code;

  
  switch(splitCode){
  case 0:
    break;
  case 1:
    switch(formercode){
    case 0:
      
      //split edge[1] and edge[3]
      if(orient_edgeID_c2f(1,orientCode) < orient_edgeID_c2f(3, orientCode)){
        edge[1]->split(node_list);
        edge[3]->split(node_list);
        
      }
      else{
        edge[3]->split(node_list);
        edge[1]->split(node_list);
      }
        
      //define newEdge
      newEdge = new Edge(edge[3]->child[0]->tail,
                         edge[1]->child[0]->tail,
                         edge[0]->level);
      edge_list.push_back(newEdge);
      //define child
      childy = new QuadFace*[2];
       
      childy[0] = new QuadFace(4);
       
      
      childy[0]->edge[0] = edge[0];
      childy[0]->edge[1] = edge[1]->child[0];
      childy[0]->edge[2] = newEdge;
      childy[0]->edge[3] = edge[3]->child[0];
      
       
      childy[1] = new QuadFace(4);
       
      
      childy[1]->edge[0] = newEdge;
      childy[1]->edge[1] = edge[1]->child[1];
      childy[1]->edge[2] = edge[2];
      childy[1]->edge[3] = edge[3]->child[1];
      
      //define code
      code = 1;
      break;
    case 1:
      
      break;
    case 2:
     //  if(childx[0]->childx != 0) cerr<<"WARNING: 1->2 quadface split with deeper childx"<< endl;
//       if(childx[0]->child != 0) cerr<<"WARNING: 1->2 quadface split with deeper child"<< endl;
//       if(childx[1]->childx != 0) cerr<<"WARNING: 1->2 quadface split with deeper childx"<< endl;
//       if(childx[1]->child != 0) cerr<<"WARNING: 1->2 quadface split with deeper child"<< endl<< endl;
      
      //split formerChild
      childx[0]->split(char(1),orientCode, node_list, edge_list);
      childx[1]->split(char(1),orientCode, node_list, edge_list);
     
    
      
      
      //define child
      if(child == 0){
        child = new QuadFace*[4];
        child[0] = childx[0]->childy[0];
        child[1] = childx[0]->childy[1];
        child[2] = childx[1]->childy[0];
        child[3] = childx[1]->childy[1];
      }
      //define childy
      if(childy == 0){
        if(child[1]->edge[0]->parent == 0){
          
          tempEdge = new Edge(child[1]->edge[0]->head, child[3]->edge[0]->tail, edge[0]->level);
          edge_list.push_back(tempEdge);
         
          tempEdge->child  = new Edge*[2];
          tempEdge->child[0] = child[1]->edge[0];
          tempEdge->child[1] = child[3]->edge[0];

          edge_list.remove(tempEdge->child[0]);
          edge_list.remove(tempEdge->child[1]);
          
          child[1]->edge[0]->parent = child[3]->edge[0]->parent = tempEdge;
        }
          
        childy = new QuadFace*[2];
        childy[0] = new QuadFace(4);
        childy[1] = new QuadFace(4);
        
        
        childy[0]->edge[0] = childx[0]->edge[0]->parent;
        childy[0]->edge[1] = child[2]->edge[1];
        childy[0]->edge[2] = child[1]->edge[0]->parent;
        childy[0]->edge[3] = child[0]->edge[3];
        

        
        childy[1]->edge[0] = child[1]->edge[0]->parent;
        childy[1]->edge[1] = child[3]->edge[1];
        childy[1]->edge[2] = child[1]->edge[2]->parent;
        childy[1]->edge[3] = child[1]->edge[3];
        
        
        childy[0]->code = childy[1]->code = 2;
        childy[0]->childx = new QuadFace*[2];
        childy[1]->childx = new QuadFace*[2];
        childy[0]->childx[0] = child[0];
        childy[0]->childx[1] = child[2];
        childy[1]->childx[0] = child[1];
        childy[1]->childx[1] = child[3];
      }
      
      //define code
      code = 3;
      break;
    case 3:
      
      break;
    default:
      cerr<<"WARNING: illegal face code " <<char(formercode + '0') << endl;
      break;
    }
    break;
  case 2:
    switch(formercode){
    case 0:
      
      //split edge[0] and edge[2]
      if(orient_edgeID_c2f(0, orientCode)< orient_edgeID_c2f(2,orientCode)){
        edge[0]->split(node_list);
        edge[2]->split(node_list);
       
      }
      else{
        edge[2]->split(node_list);
        edge[0]->split(node_list);
       
      }
      //define newEdge
      newEdge = new Edge(edge[0]->child[0]->tail,
                         edge[2]->child[0]->tail,
                         edge[1]->level);
      edge_list.push_back(newEdge);

      //define child
      childx = new QuadFace*[2];
      
      childx[0] = new QuadFace(4);
      
    
      childx[0]->edge[0] = edge[0]->child[0];
      childx[0]->edge[1] = newEdge;
      childx[0]->edge[2] = edge[2]->child[0];
      childx[0]->edge[3] = edge[3];
     
      
      childx[1] = new QuadFace(4);
     
      childx[1]->edge[0] = edge[0]->child[1];
      childx[1]->edge[1] = edge[1];
      childx[1]->edge[2] = edge[2]->child[1];
      childx[1]->edge[3] = newEdge;
     
      
      //define code
      code = 2;
      break;
    case 1:
      // if(childy[0]->childy != 0) cerr<<"WARNING: 2->1 quadface split with deeper childy"<< endl;
//       if(childy[0]->child != 0) cerr<<"WARNING: 2->1 quadface split with deeper child"<< endl;
//       if(childy[1]->childy != 0) cerr<<"WARNING: 2->1 quadface split with deeper childy"<< endl;
//       if(childy[1]->child != 0) cerr<<"WARNING: 2->1 quadface split with deeper child"<< endl<< endl;

      //split formerChild
      childy[0]->split(char(2),orientCode, node_list, edge_list);
      childy[1]->split(char(2),orientCode, node_list, edge_list);
      
      //define child
      child = new QuadFace*[4];
      
      child[0] = childy[0]->childx[0];
      child[2] = childy[0]->childx[1];
      child[1] = childy[1]->childx[0];
      child[3] = childy[1]->childx[1];

      //define childy
      if(childx == 0){
        if(child[0]->edge[1]->parent == 0){
          tempEdge = new Edge(child[0]->edge[1]->head, child[1]->edge[1]->tail, edge[1]->level);
          edge_list.push_back(tempEdge);
          tempEdge->child  = new Edge*[2];
          tempEdge->child[0] = child[0]->edge[1];
          tempEdge->child[1] = child[1]->edge[1];

          edge_list.remove(tempEdge->child[0]);
          edge_list.remove(tempEdge->child[1]);
          child[0]->edge[1]->parent = child[1]->edge[1]->parent = tempEdge;
        }
        
        childx = new QuadFace*[2];
        childx[0] = new QuadFace(4);
        childx[1] = new QuadFace(4);

       
        childx[0]->edge[0] = child[0]->edge[0];
        childx[0]->edge[1] = child[0]->edge[1]->parent;
        childx[0]->edge[2] = child[1]->edge[2];
        childx[0]->edge[3] = child[0]->edge[3]->parent;
       

       
        childx[1]->edge[0] = child[2]->edge[0];
        childx[1]->edge[1] = child[2]->edge[1]->parent;
        childx[1]->edge[2] = child[3]->edge[2];
        childx[1]->edge[3] = child[0]->edge[1]->parent;
       
        
        childx[0]->code = childx[1]->code = 1;
        childx[0]->childy = new QuadFace*[2];
        childx[1]->childy = new QuadFace*[2];
        childx[0]->childy[0] = child[0];
        childx[0]->childy[1] = child[1];
        childx[1]->childy[0] = child[2];
        childx[1]->childy[1] = child[3];
      }
      //define code
      code = 3;
      break;
    case 2:
      
      break;
    case 3:
      break;
    default:
      cerr <<"WARNING: illegal face code  " <<char(formercode + '0') << endl; 
      break;
    }
    break;
  case 3:
    switch(formercode){
    case 0:
      //split edges
      for(char i = 0; i < 4; i++){
        char edgeID = orient_edgeID_f2c(i, orientCode);
        if(edge[edgeID]->child == 0) edge[edgeID]->split(node_list);
      }
      
      //define new Node
      facecenter = centroid();
      node_list.push_back(facecenter);

      //get edgecenter
      getEdgeCenter(edgecenter);
      
      //define new edges
      newedge[0] = new Edge(edgecenter[0],facecenter,(edge[1]->level)+1);
      newedge[1] = new Edge(facecenter, edgecenter[1], (edge[0]->level)+1);
      newedge[2] = new Edge(facecenter, edgecenter[2], (edge[1]->level)+1);
      newedge[3] = new Edge(edgecenter[3],facecenter,(edge[0]->level)+1);
      
      tempEdge = new Edge(newedge[0]->head, newedge[2]->tail, edge[1]->level);
      edge_list.push_back(tempEdge);
      tempEdge->child = new Edge*[2];
      tempEdge->child[0] = newedge[0];
      tempEdge->child[1] = newedge[2];
      newedge[0]->parent = newedge[2]->parent = tempEdge;
     
       
      tempEdge = new Edge(newedge[3]->head, newedge[1]->tail, edge[0]->level);
      edge_list.push_back(tempEdge);
      tempEdge->child = new Edge*[2];
      tempEdge->child[0] = newedge[3];
      tempEdge->child[1] = newedge[1];
      newedge[3]->parent = newedge[1]->parent = tempEdge;
     
      
      
      //define child
      if(child !=0)cerr<<"WARNING: error in quadface split" << endl;
      child = new QuadFace*[4];
      for(int i = 0; i<4; i++){
        child[i] = new QuadFace(4);
       
      }


      
      child[0]->edge[0] = edge[0]->child[0];
      child[0]->edge[1] = newedge[0];
      child[0]->edge[2] = newedge[3];
      child[0]->edge[3] = edge[3]->child[0];
      
      
      
      
      child[1]->edge[0] = newedge[3];
      child[1]->edge[1] = newedge[2];
      child[1]->edge[2] = edge[2]->child[0];
      child[1]->edge[3] = edge[3]->child[1];
      
      
      child[2]->edge[0] = edge[0]->child[1];
      child[2]->edge[1] = edge[1]->child[0];
      child[2]->edge[2] = newedge[1];
      child[2]->edge[3] = newedge[0];
      
      
      
      child[3]->edge[0] = newedge[1];
      child[3]->edge[1] = edge[1]->child[1];
      child[3]->edge[2] = edge[2]->child[1];
      child[3]->edge[3] = newedge[2];
      
      
      //define childx
      childx = new QuadFace*[2];
      
      for(int i = 0; i<2; i++){
        childx[i] = new QuadFace(4);
      }

      childx[0]->edge[0] = child[0]->edge[0];
      childx[0]->edge[1] = child[0]->edge[1]->parent;
      childx[0]->edge[2] = child[1]->edge[2];
      childx[0]->edge[3] = child[0]->edge[3]->parent;
       
      
      childx[1]->edge[0] = child[2]->edge[0];
      childx[1]->edge[1] = child[2]->edge[1]->parent;
      childx[1]->edge[2] = child[3]->edge[2];
      childx[1]->edge[3] = child[0]->edge[1]->parent;
       
      
      childx[0]->code = childx[1]->code = 1;
      childx[0]->childy = new QuadFace*[2];
      childx[1]->childy = new QuadFace*[2];
      childx[0]->childy[0] = child[0];
      childx[0]->childy[1] = child[1];
      childx[1]->childy[0] = child[2];
      childx[1]->childy[1] = child[3];

      childy = new QuadFace*[2];
      childy[0] = new QuadFace(4);
      childy[1] = new QuadFace(4);
       
       
      childy[0]->edge[0] = childx[0]->edge[0]->parent;
      childy[0]->edge[1] = child[2]->edge[1];
      childy[0]->edge[2] = child[1]->edge[0]->parent;
      childy[0]->edge[3] = child[0]->edge[3];
       
       
      childy[1]->edge[0] = child[1]->edge[0]->parent;
      childy[1]->edge[1] = child[3]->edge[1];
      childy[1]->edge[2] = child[1]->edge[2]->parent;
      childy[1]->edge[3] = child[1]->edge[3];
       
       
      childy[0]->code = childy[1]->code = 2;
      childy[0]->childx = new QuadFace*[2];
      childy[1]->childx = new QuadFace*[2];
      childy[0]->childx[0] = child[0];
      childy[0]->childx[1] = child[2];
      childy[1]->childx[0] = child[1];
      childy[1]->childx[1] = child[3];
      //define code
      code = 3;
      break;
    case 1:
      // if(childy[0]->childy != 0) cerr<<"WARNING: 3->1 quadface split with deeper childy"<< endl;
//       if(childy[0]->child != 0) cerr<<"WARNING: 3->1 quadface split with deeper child"<< endl;
//       if(childy[1]->childy != 0) cerr<<"WARNING: 3->1 quadface split with deeper childy"<< endl;
//       if(childy[1]->child != 0) cerr<<"WARNING: 3->1 quadface split with deeper child"<<endl <<endl;

      //split formerChild
      childy[0]->split(char(2),orientCode, node_list, edge_list);
      childy[1]->split(char(2),orientCode, node_list, edge_list);
      
      //define child
      child = new QuadFace*[4];
      
      child[0] = childy[0]->childx[0];
      child[2] = childy[0]->childx[1];
      child[1] = childy[1]->childx[0];
      child[3] = childy[1]->childx[1];

      //define childx
      if(childx == 0){
        if(child[0]->edge[1]->parent == 0){
          tempEdge = new Edge(child[0]->edge[1]->head, child[1]->edge[1]->tail, edge[1]->level);
          edge_list.push_back(tempEdge);
          tempEdge->child  = new Edge*[2];
          tempEdge->child[0] = child[0]->edge[1];
          tempEdge->child[1] = child[1]->edge[1];
          edge_list.remove(tempEdge->child[0]);
          edge_list.remove(tempEdge->child[1]);
          child[0]->edge[1]->parent = child[1]->edge[1]->parent = tempEdge;
        }
        
        childx = new QuadFace*[2];
        childx[0] = new QuadFace(4);
        childx[1] = new QuadFace(4);
        
        
        childx[0]->edge[0] = child[0]->edge[0];
        childx[0]->edge[1] = child[0]->edge[1]->parent;
        childx[0]->edge[2] = child[1]->edge[2];
        childx[0]->edge[3] = child[0]->edge[3]->parent;
        
        
        childx[1]->edge[0] = child[2]->edge[0];
        childx[1]->edge[1] = child[2]->edge[1]->parent;
        childx[1]->edge[2] = child[3]->edge[2];
        childx[1]->edge[3] = child[0]->edge[1]->parent;
        
        
        childx[0]->code = childx[1]->code = 1;
        childx[0]->childy = new QuadFace*[2];
        childx[1]->childy = new QuadFace*[2];
        childx[0]->childy[0] = child[0];
        childx[0]->childy[1] = child[1];
        childx[1]->childy[0] = child[2];
        childx[1]->childy[1] = child[3];
      }

      
      code = 3;
      break;
    case 2:
     //  if(childx[0]->childx != 0) cerr<<"WARNING: 3->2 quadface split with deeper childx"<< endl;
//       if(childx[0]->child != 0) cerr<<"WARNING: 3->2 quadface split with deeper child"<< endl;
//       if(childx[1]->childx != 0) cerr<<"WARNING: 3->2 quadface split with deeper childx"<< endl;
//       if(childx[1]->child != 0) cerr<<"WARNING: 3->2 quadface split with deeper child"<< endl<<endl;;

      //split formerChild
      childx[0]->split(char(1),orientCode, node_list, edge_list);
      childx[1]->split(char(1),orientCode, node_list, edge_list);
      //define child
      if(child == 0){
        child = new QuadFace*[4];
        child[0] = childx[0]->childy[0];
        child[1] = childx[0]->childy[1];
        child[2] = childx[1]->childy[0];
        child[3] = childx[1]->childy[1];
      }
      //define childy
      if(childy == 0){
        if(child[1]->edge[0]->parent == 0){
          tempEdge = new Edge(child[1]->edge[0]->head, child[3]->edge[0]->tail, edge[0]->level);
          edge_list.push_back(tempEdge);
          tempEdge->child  = new Edge*[2];
          tempEdge->child[0] = child[1]->edge[0];
          tempEdge->child[1] = child[3]->edge[0];
          edge_list.remove(tempEdge->child[0]);
          edge_list.remove(tempEdge->child[1]);
          child[1]->edge[0]->parent = child[3]->edge[0]->parent = tempEdge;
        }
        
        childy = new QuadFace*[2];
        childy[0] = new QuadFace(4);
        childy[1] = new QuadFace(4);
       
        
        childy[0]->edge[0] = childx[0]->edge[0]->parent;
        childy[0]->edge[1] = child[2]->edge[1];
        childy[0]->edge[2] = child[1]->edge[0]->parent;
        childy[0]->edge[3] = child[0]->edge[3];
       
        
        childy[1]->edge[0] = child[1]->edge[0]->parent;
        childy[1]->edge[1] = child[3]->edge[1];
        childy[1]->edge[2] = child[1]->edge[2]->parent;
        childy[1]->edge[3] = child[1]->edge[3];
       
        
        childy[0]->code = childy[1]->code = 2;
        childy[0]->childx = new QuadFace*[2];
        childy[1]->childx = new QuadFace*[2];
        childy[0]->childx[0] = child[0];
        childy[0]->childx[1] = child[2];
        childy[1]->childx[0] = child[1];
        childy[1]->childx[1] = child[3];

      }
      
      code = 3;
      break;
    case 3:
      break;
    default:
      cerr <<"WARNING: illegal face code     " <<char(formercode + '0')<< endl;
      break;
    }
    break;
  default:
    cerr <<"WARNING: illegal split code " << char(splitCode +'0') <<endl;
    break;
  }
   
}

//this function cannot apply to the same tree more than once
void QuadFace::empty_split(char splitCode ){
  if(code != 0){
    Loci::debugout<< "WARNING: QuadFace::empty_split() is called on a face that code is nonzero." <<endl;
    return;
  }
  code = splitCode;
  switch(code){
  case 0:
    break;
  case 1:
    //define child
    childy = new QuadFace*[2];
    childy[0] = new QuadFace();
    childy[1] = new QuadFace();
    break;
  case 2:
    childx = new QuadFace*[2];
    childx[0] = new QuadFace();
    childx[1] = new QuadFace();
    break;
  case 3:
    if(child!=0) cerr<<"WARNING: errro in empty_split " << endl;
    child = new QuadFace*[4];
    for(int i = 0; i<4; i++){
      child[i] = new QuadFace();
    }
    //define childx
    childx = new QuadFace*[2];
    
    for(int i = 0; i<2; i++){
      childx[i] = new QuadFace();
    }
    //define childx
    childy = new QuadFace*[2];
    
    for(int i = 0; i<2; i++){
      childy[i] = new QuadFace();
    }
    childx[0]->code = childx[1]->code = 1;
    childx[0]->childy = new QuadFace*[2];
    childx[1]->childy = new QuadFace*[2];
    childx[0]->childy[0] = child[0];
    childx[0]->childy[1] = child[1];
    childx[1]->childy[0] = child[2];
    childx[1]->childy[1] = child[3]; 

    childy[0]->code = childy[1]->code = 2;
    childy[0]->childx = new QuadFace*[2];
    childy[1]->childx = new QuadFace*[2];
    childy[0]->childx[0] = child[0];
    childy[0]->childx[1] = child[2];
    childy[1]->childx[0] = child[1];
    childy[1]->childx[1] = child[3];

    break;
  default:
    cerr <<"WARNING: illegal split code " << char(splitCode +'0') <<endl;
    break;
  }
}

//used in building cells, quadface is built as defined in cell, and split with orientCode
//all new nodes and edges are put into node_list and edge_list 
void QuadFace::resplit(const std::vector<char>& facePlan, char orientCode,
                       std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list ){
 
  if(facePlan.size() == 0) return;
  std::queue<QuadFace*> Q;
  Q.push(this);
  QuadFace* current;
  unsigned int index =0;
  char currentCode;

  while(!Q.empty()){
    current = Q.front();
    
    if(index >= facePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from splitcode
      currentCode = facePlan[index];
      index++;  
    }
    
    if(currentCode !=0){
      char tmpCode =orient_splitCode_f2c(currentCode, orientCode); 
      current->split(orient_splitCode_f2c(currentCode, orientCode),orientCode, node_list, edge_list);
     
      switch(tmpCode){
      case 1:
        for(int i = 0; i < 2; i++){
          Q.push(current->childy[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      case 2:
        for(int i = 0; i < 2; i++){
          Q.push(current->childx[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      case 3:
        for(int i = 0; i < 4; i++){
          Q.push(current->child[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      case 0:
        break;
      default:
        cerr << "illegal face code in resplit(): " <<char(current->code + '0') << endl; 
        break;
      }
    }
    Q.pop();
  }
}



//this function split the local face of a cell according to global face(defined by face2node)'s plan
void QuadFace::resplit(const std::vector<char>& facePlan, char orientCode,
                       std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list,
                       std::vector<QuadFace*>& fine_faces ){
 
  if(facePlan.size() == 0) {
    fine_faces.push_back(this);
    reduce_vector(fine_faces);
    return;
  }
  
  
  std::queue<QuadFace*> Q;
  Q.push(this);
  QuadFace* current;
  unsigned int index =0;
  char currentCode;
  
  while(!Q.empty()){
    current = Q.front();
    
    if(index >= facePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from splitcode
      currentCode = facePlan[index];
      index++;  
    }

    if(currentCode == 0){
      fine_faces.push_back(current);
    }else{
      char tmpCode =orient_splitCode_f2c(currentCode, orientCode); 
      current->split(orient_splitCode_f2c(currentCode, orientCode),orientCode, node_list, edge_list);
      
      switch(tmpCode){
      case 1:
        for(int i = 0; i < 2; i++){
          Q.push(current->childy[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      case 2:
        for(int i = 0; i < 2; i++){
          Q.push(current->childx[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      case 3:
        for(int i = 0; i < 4; i++){
          Q.push(current->child[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      default:
        break;
      }
    }
    Q.pop();
  }
}

// this has been built as defined in cell, and has been resplit with orientCode
//get leaves that is in the same order as the face is bulit as defined by face2node and resplit without orientCode
void QuadFace::get_leaves(const std::vector<char>& facePlan, char orientCode,
                          std::vector<QuadFace*>& fine_faces ){
  if(facePlan.size() == 0) {
    fine_faces.push_back(this);
    reduce_vector(fine_faces);
    return;
  }
  
  std::queue<QuadFace*> Q;
  Q.push(this);
  QuadFace* current;
  unsigned int index =0;
  char currentCode;
  
  while(!Q.empty()){
    current = Q.front();
    
    if(index >= facePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from splitcode
      currentCode = facePlan[index];
      index++;  
    }
    
    if(currentCode == 0){
      fine_faces.push_back(current);
    }
    else {
      switch(current->code){
      case 1:
        for(int i = 0; i < 2; i++){
          Q.push(current->childy[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      case 2:
        for(int i = 0; i < 2; i++){
          Q.push(current->childx[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      case 3:
        for(int i = 0; i < 4; i++){
          Q.push(current->child[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      default:
        break;
      }
    }
    Q.pop();
  }
}

//this function split the local face of a cell according to global face(define by face2node)'s plan
//fine_faces is in the same order as global face's leaves
//this function can not apply to the same tree more than once
void QuadFace::empty_resplit(const std::vector<char>& facePlan, char orientCode,
                             std::vector<QuadFace*>& fine_faces ){
  
  if(facePlan.size() == 0) {
    fine_faces.push_back(this);
    reduce_vector(fine_faces);
    return;
  }
  
  std::queue<QuadFace*> Q;
  Q.push(this);
  QuadFace* current;
  unsigned int index =0;
  char currentCode;
  
  while(!Q.empty()){
    current = Q.front();
    
    if(index >= facePlan.size()){
      currentCode = 0;
    }
    else{ 
      //take a code from splitcode
      currentCode = facePlan[index];
      index++;  
    }

    if(currentCode == 0){
      fine_faces.push_back(current);
    }else if((current->code)!=0 && currentCode != (current->code)){
      Loci::debugout << "WARNING: split code is not consistent" << endl;
      Loci::debugout << int(currentCode) << " intree " << int( current->code)<< endl; 
      Loci::Abort();
    }else{
      if((current->code)==0){
        current->empty_split(orient_splitCode_f2c(currentCode, orientCode));
      }
      switch(current->code){
      case 1:
        for(int i = 0; i < 2; i++){
          Q.push(current->childy[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      case 2:
        for(int i = 0; i < 2; i++){
          Q.push(current->childx[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      case 3:
        for(int i = 0; i < 4; i++){
          Q.push(current->child[orient_childID_f2c(i, orientCode, currentCode)]);
        }
        break;
      default:
        break;
      }
    }
    Q.pop();
  }
}



// //get all the leaves of this
void QuadFace::get_leaves(std::vector<QuadFace*>& leaves){
  if(code == 0){
    leaves.push_back(this);
    return;
  }
  switch(code){
  case 1:
    for(int i = 0; i <2; i++){
      childy[i]->get_leaves(leaves);
    }
    break;
  case 2:
    for(int i = 0; i <2; i++){
      childx[i]->get_leaves(leaves);
    }
    break;
  case 3:
    for(int i = 0; i <4; i++){
      child[i]->get_leaves(leaves);
    }
    break;
  default:
    cerr<< " WARNING: illegal code in QuadFace::get_leaves() " << char(code+'0')<< endl;
    break;
  }
}

//if the intersection of the leaves of f1 and the leaves of f2 is empty 
bool is_overlapped( QuadFace* f1,  QuadFace* f2){
  if(f1 == f2) return true;
  
  std::vector<QuadFace*> leaves1;
  f1->get_leaves(leaves1);
  std::vector<QuadFace*> leaves2;
  f2->get_leaves(leaves2);

  std::vector<QuadFace*> intersection;
  std::sort(leaves1.begin(), leaves1.end());
  std::sort(leaves2.begin(), leaves2.end());
  std::set_intersection(leaves1.begin(), leaves1.end(), leaves2.begin(), leaves2.end(), std::inserter(intersection, intersection.begin()));
  return !(intersection.empty());
}

//return the intersection of leaves of f1 and leaves of f2
std::vector<QuadFace*> overlap( QuadFace* f1,   QuadFace* f2){
  std::vector<QuadFace*> intersection;
    
  std::vector<QuadFace*> leaves1;
  f1->get_leaves(leaves1);
  std::vector<QuadFace*> leaves2;
  f2->get_leaves(leaves2);
 
  
  
  std::sort(leaves1.begin(), leaves1.end());
  std::sort(leaves2.begin(), leaves2.end());
  std::set_intersection(leaves1.begin(), leaves1.end(), leaves2.begin(), leaves2.end(), std::inserter(intersection, intersection.begin()));
 
  return intersection;
}

int QuadFace::get_num_leaves()const{
  if(code ==0) return 1;
  
  int count = 0;
  switch(code){
  case 1:
    for(int i = 0; i < 2; i++){
      count += childy[i]->get_num_leaves();
    }
    break;
  case 2:
    for(int i = 0; i < 2; i++){
      count += childx[i]->get_num_leaves();
    }
    break;
  case 3:
    for(int i = 0; i < 4; i++){
      count += child[i]->get_num_leaves();
    }
    break;
  default:
    cerr << "WARNING: illegal code in QuadFace::get_num_leaves()" << endl;
    break;
  }
  return count;
}

//define face2node
void QuadFace::set_f2n(std::list<int32>& f2n){
  f2n.clear();
  for(int i = 0; i < 4; i++){
    
    //each edge sort leaves
    std::list<Edge*> edge_leaves;
    edge[i]->sort_leaves(edge_leaves);

    
    if(i >=2 ){
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
QuadFace* build_quad_face( const Entity* face2node, 
                           const Entity* face2edge,
                           const const_MapVec<2>& edge2node,
                           const const_store<vect3d>& pos,
                           const const_store<std::vector<char> >& edgePlan,
                           std::list<Node*>& bnode_list,
                           std::list<Edge*>& edge_list){
  

  Node* node[4];
  
  for(int nindex = 0; nindex < 4; nindex++){
    node[nindex] = new Node(pos[face2node[nindex]]);
    bnode_list.push_back(node[nindex]);
  }
  
  //define each edge and put it into edge_list
  
  Edge** edge = new Edge*[4];
  bool needReverse[4];
          
  for(int eindex = 0; eindex < 4; eindex++){
    //define the edge
    edge[eindex] = new Edge();
    edge_list.push_back(edge[eindex]);
  }
   
  edge[0]->head = node[0];
  edge[0]->tail = node[1];
  needReverse[0] = (edge2node[face2edge[0]][0] == face2node[1]);  
   

  edge[1]->head = node[1];
  edge[1]->tail = node[2];
  needReverse[1] = (edge2node[face2edge[1]][0] == face2node[2]);

  
  edge[2]->head = node[3];
  edge[2]->tail = node[2];
  needReverse[2] = (edge2node[face2edge[2]][0] == face2node[2]);
  
  
  edge[3]->head = node[0];
  edge[3]->tail = node[3];
  needReverse[3] = (edge2node[face2edge[3]][0] == face2node[3]);
  
  for(int eindex = 0; eindex < 4; eindex++){
    //resplit the edge
    edge[eindex]->resplit(edgePlan[face2edge[eindex]],needReverse[eindex],bnode_list);
    
  }
  
  //define the face
  QuadFace* aFace = new QuadFace(edge);
 
  return aFace;
}


//parallel version
QuadFace* build_quad_face( const Entity* face2node, 
                           const Entity* face2edge,
                           const const_MapVec<2>& edge2node,
                           const const_store<vect3d>& pos,
                           const const_store<std::vector<char> >& edgePlan,
                           const const_store<int>& node_offset,
                           const const_store<int>&  node_l2f,
                           std::list<Node*>& bnode_list,
                           std::list<Edge*>& edge_list){
  
  
  Node* node[4];
  
  for(int nindex = 0; nindex < 4; nindex++){
    node[nindex] = new Node(pos[face2node[nindex]], node_l2f[face2node[nindex]]);
    bnode_list.push_back(node[nindex]);
  }
  
  //define each edge and put it into edge_list
  Edge** edge = new Edge*[4];
  bool needReverse[4];

  //define edges and index its inner nodes
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());        

  for(int eindex = 0; eindex < 4; eindex++){
    //define the edge
    edge[eindex] = new Edge();
    edge_list.push_back(edge[eindex]);
  }
   
  edge[0]->head = node[0];
  edge[0]->tail = node[1];
  needReverse[0] = (edge2node[face2edge[0]][0] == face2node[1]);  
   

  edge[1]->head = node[1];
  edge[1]->tail = node[2];
  needReverse[1] = (edge2node[face2edge[1]][0] == face2node[2]);

  edge[2]->head = node[3];
  edge[2]->tail = node[2];
  needReverse[2] = (edge2node[face2edge[2]][0] == face2node[2]);

  edge[3]->head = node[0];
  edge[3]->tail = node[3];
  needReverse[3] = (edge2node[face2edge[3]][0] == face2node[3]);
    
  for(int eindex = 0; eindex < 4; eindex++){
   
    edge[eindex]->resplit(edgePlan[face2edge[eindex]],needReverse[eindex],bnode_list);
    
    
    int nindex = node_offset[face2edge[eindex]];
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->index =  nindex++;
    }
        
    bnode_begin = --(bnode_list.end());
  }
  
  //define the face
  QuadFace* aFace = new QuadFace(edge);
 
  return aFace;
}
//parallel version, this function is used in build_general_cell with quadface
QuadFace* build_tmp_quad_face( const Entity* face2node, 
                               const Entity* face2edge,
                               const const_MapVec<2>& edge2node,
                               const const_store<std::vector<char> >& edgePlan,
                               std::list<Node*>& bnode_list,
                               std::list<Edge*>& edge_list){
  
  
  Node* node[4];
  vect3d p[4];
  int64 maxX = int64(1) << MAXLEVEL;
  int64 maxY = int64(1) << MAXLEVEL;
  p[0] = vect3d(0.0, 0.0, 0.0);
  p[1] = vect3d(maxX, 0.0, 0.0);
  p[2] = vect3d(maxX, maxY, 0.0);
  p[3] = vect3d(0.0, maxY, 0.0);
  for(int nindex = 0; nindex < 4; nindex++){
    node[nindex] = new Node(p[nindex]);
    bnode_list.push_back(node[nindex]);
  }
  
  //define each edge and put it into edge_list
  Edge** edge = new Edge*[4];
  bool needReverse[4];

  //define edges
  for(int eindex = 0; eindex < 4; eindex++){
    //define the edge
    edge[eindex] = new Edge();
    edge_list.push_back(edge[eindex]);
  }
   
  edge[0]->head = node[0];
  edge[0]->tail = node[1];
  needReverse[0] = (edge2node[face2edge[0]][0] == face2node[1]);  
   

  edge[1]->head = node[1];
  edge[1]->tail = node[2];
  needReverse[1] = (edge2node[face2edge[1]][0] == face2node[2]);

  edge[2]->head = node[3];
  edge[2]->tail = node[2];
  needReverse[2] = (edge2node[face2edge[2]][0] == face2node[2]);

  edge[3]->head = node[0];
  edge[3]->tail = node[3];
  needReverse[3] = (edge2node[face2edge[3]][0] == face2node[3]);
    
  for(int eindex = 0; eindex < 4; eindex++){
    edge[eindex]->resplit(edgePlan[face2edge[eindex]],needReverse[eindex],bnode_list);
  }
  
  //define the face
  QuadFace* aFace = new QuadFace(edge);
  return aFace;
}

void tag_quad_face( const Entity* face2node, 
                    const Entity* face2edge,
                    const const_MapVec<2>& edge2node,
                    const const_store<std::vector<char> >& edgePlan,
                    const std::vector<char>& facePlan, char orientCode,
                    const std::vector<char>& nodeTag,//the tag for facePlan 
                    const std::vector<char>& facePlan1,
                    std::list<Node*>& bnode_list,//node list from facePlan1
                    std::list<Node*>::const_iterator bnode_begin){

  std::list<Node*> tmp_bnode_list, tmp_node_list, tmp_node_list2;
  std::list<Edge*> edge_list;
  QuadFace* tmp_qface = build_tmp_quad_face(face2node,
                                            face2edge,
                                            edge2node,
                                            edgePlan,
                                            tmp_bnode_list,
                                            edge_list);
   //resplit the quadface to get node index
  tmp_qface->resplit(facePlan,
                     char(0),
                     tmp_node_list,
                     edge_list);


  
  int   nindex = 0;
  for(std::list<Node*>::const_iterator np = tmp_node_list.begin(); np!= tmp_node_list.end(); np++){
    (*np)->tag = nodeTag[nindex++];
  }

  
  QuadFace* tmp_qface1 = build_tmp_quad_face(face2node,
                                             face2edge,
                                             edge2node,
                                             edgePlan,
                                             tmp_bnode_list,
                                             edge_list);
  tmp_qface1->resplit(facePlan1,
                      char(0),
                      tmp_node_list2,
                      edge_list);

  std::list<Node*>::const_iterator tmp_np2 = tmp_node_list2.begin();
  for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++, tmp_np2++){
    for(std::list<Node*>::const_iterator tmp_np = tmp_node_list.begin(); tmp_np!= tmp_node_list.end(); tmp_np++){
      if(int_equal((*tmp_np)->p, (*tmp_np2)->p)){
        (*np)->tag = (*tmp_np)->tag;
        break;
      }
    }
    //it's OK if node_found is false
  }
  //cleanup
  delete tmp_qface;
  delete tmp_qface1;
  cleanup_list(tmp_node_list, edge_list);
  cleanup_list(tmp_node_list2);
  cleanup_list(tmp_bnode_list);
}
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
