//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
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
// This file build a Cell according to Loci data structures and refinementplans
#include <vector>
#include <list>
#include <Loci.h>
#include "diamondcell.h"
#include "defines.h"
#include "quadface.h"
using std::cerr;
using std::endl;
using std::swap;
using std::cout;
using std::vector;
// void checkNode(Node* nd){
//   vect3d p = nd->p;
//   if(abs(p.x)>0.2 && abs(p.x)<0.5){
//     if(nd->tag != 1) cerr << " p: " << p << " ," << char(nd->tag + '0')<< endl;
//   }
//   else if(nd->tag == 1) cerr << " p: " << p << " ," << char(nd->tag + '0')<< endl;
// }
          

std::vector<char> transfer_plan_q2g(const std::vector<char>& facePlan);




std::vector<Entity> reorder_nodes(const const_store<int>& node_remap, const entitySet& localSet){
  
  //reverse the map 
  std::vector<pair<int, Entity> > node_f2l(localSet.size());
  int index = 0;
  for(entitySet::const_iterator ei = localSet.begin(); ei != localSet.end(); ei++, index++){
    node_f2l[index] = pair<int, Entity>(node_remap[*ei], *ei);
  }
  std::sort(node_f2l.begin(), node_f2l.end());
  std::vector<Entity> orderedSet(localSet.size());
  for( size_t i= 0; i < localSet.size(); i++){
    orderedSet[i] = node_f2l[i].second;
  }
  return orderedSet;
}

std::vector<Entity> reorder_edges(const const_store<int>& node_remap,const const_MapVec<2>& edge2node, const entitySet& localSet){
  
  //reverse the map 
  std::vector<pair<std::vector<int>, Entity> > node_f2l(localSet.size());
  int index = 0;
  for(entitySet::const_iterator ei = localSet.begin(); ei != localSet.end(); ei++, index++){
    vector<int> e2n(2);
    //  e2n[0] = min(node_remap[edge2node[*ei][0]],node_remap[edge2node[*ei][1]]) ;
    //     e2n[1] = max(node_remap[edge2node[*ei][0]],node_remap[edge2node[*ei][1]]) ;

    e2n[0] = node_remap[edge2node[*ei][0]] ;
    e2n[1] = node_remap[edge2node[*ei][0]] ;
    node_f2l[index] = pair<vector<int>, Entity>(e2n,*ei);
  }
  std::sort(node_f2l.begin(), node_f2l.end());
  std::vector<Entity> orderedSet(localSet.size());
  for( size_t i= 0; i < localSet.size(); i++){
    orderedSet[i] = node_f2l[i].second;
  }
  return orderedSet;
}


void reorder_faces(const const_store<int>& node_remap, const const_multiMap& face2node, std::vector<Entity>& localSet, char* orient){
  
  //reverse the map 
  std::vector<pair<vector<int>, pair<Entity, char> > > node_f2l(localSet.size());
  for(unsigned int  index  = 0; index < localSet.size(); index++){
    vector<int> f2n(face2node.num_elems(localSet[index]));
    for(int i = 0; i< face2node.num_elems(localSet[index]); i++){
      f2n[i] = node_remap[face2node[localSet[index]][i]];
    }
    
    node_f2l[index] = pair<vector<int>, pair<Entity, char> >(f2n,
                                                             pair<Entity, char>(localSet[index], orient[index]));
  }
  std::sort(node_f2l.begin(), node_f2l.end());
  
  for( unsigned int i= 0; i < localSet.size(); i++){
    localSet[i] = (node_f2l[i].second).first;
    orient[i] = (node_f2l[i].second).second;
  }
  
}




//parallel version in set_general_nums.cc
Cell* build_general_cell(const Entity* lower, int lower_size,
                         const Entity* upper, int upper_size,
                         const Entity* boundary_map, int boundary_map_size,
                         const const_store<bool>& is_quadface,
                         const const_multiMap& face2node,
                         const const_multiMap& face2edge,
                         const const_MapVec<2>& edge2node,
                         const const_store<vect3d>& pos,
                         const const_store<std::vector<char> >& edgePlan,
                         const const_store<std::vector<char> >& facePlan,
                         std::list<Node*>& bnode_list,
                         std::list<Edge*>& edge_list,
                         std::list<Face*>& face_list,
                         const const_store<int>& node_remap){
  int numFaces = lower_size + upper_size + boundary_map_size;
  std::vector<Entity> faces(numFaces);
  
  //orient value: upper and boundary_map: 1
  //              lower: -1
  char* orient = new char[numFaces];
  
  int findex = 0;
  for(int i=0; i<lower_size; i++){
    faces[findex] = lower[i];
    orient[findex++] = 1;
  }
  
  for(int i=0; i<upper_size; i++){
    faces[findex] = upper[i];
    orient[findex++] = 0;
    
  }
  for(int i=0; i<boundary_map_size; i++){
    faces[findex] = boundary_map[i];
    orient[findex++] = 0;
    
  } 
  
  
  //collect all the nodes and edges of the cell
  entitySet edges, nodes;
  for(unsigned int i=0; i< faces.size(); i++){
    for( int j= 0; j< face2node.num_elems(faces[i]); j++){
      nodes += face2node[faces[i]][j];
    }
    
    for( int j = 0; j < face2edge.num_elems(faces[i]); j++){
      edges += face2edge[faces[i]][j];
    }
  }

  //reorder entities
  std::vector<Entity> orderedNodes = reorder_nodes(node_remap, nodes);
  std::vector<Entity> orderedEdges = reorder_edges(node_remap, edge2node,edges);
  reorder_faces(node_remap, face2node, faces, orient);
  
  int numEdges = edges.size();
  int numNodes = nodes.size();
  
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(entitySet::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    Node* aNode = new Node(pos[*np]);
    bnode_list.push_back(aNode);
    n2n[*np] = aNode;
  }
  
  //define each edge and put it into edge_list
  std::map<Entity, Edge*> e2e;
  for(entitySet::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    Edge* anEdge = new Edge(n2n[edge2node[*ep][0]], n2n[edge2node[*ep][1]]);
    edge_list.push_back(anEdge);
    e2e[*ep] = anEdge;
    
    //resplit the edge
    anEdge->resplit(edgePlan[*ep], bnode_list);
  }
        
  
  
  //defines each face and put it into face_list
  Face** face = new Face*[numFaces];
  for(int i  = 0; i < numFaces; i++){
    face[i] = new Face(face2edge.num_elems(faces[i]));
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < face[i]->numEdge; j++){
      face[i]->edge[j] = e2e[face2edge[faces[i]][j]];
      if(edge2node[face2edge[faces[i]][j]][0] == face2node[faces[i]][j] &&
         edge2node[face2edge[faces[i]][j]][1] == face2node[faces[i]][j==face2node.num_elems(faces[i])-1? 0:j+1]) face[i]->needReverse[j] = false;
      else face[i]->needReverse[j] = true;
    }
    
    //resplit each face
 
    if(is_quadface[faces[i]]){
      std::vector<char> newPlan = transfer_plan_q2g(facePlan[faces[i]]);
     
      face[i]->resplit(newPlan, bnode_list, edge_list);
    }
    else{
      face[i]->resplit(facePlan[faces[i]], bnode_list, edge_list);
    }
  }
  
  Node** node = new Node*[numNodes];
  
  for(unsigned int nindex = 0; nindex < orderedNodes.size(); nindex++){
    node[nindex] = n2n[orderedNodes[nindex]];
  }
  
  Edge** edge = new Edge*[numEdges];
  
  for(unsigned int eindex = 0; eindex < orderedEdges.size(); eindex++){
    edge[eindex] = e2e[orderedEdges[eindex]];
  }
  Cell* aCell = new Cell(numNodes, numEdges, numFaces, node, edge, face,orient);
  return aCell;
}


int  find_face_index(const Entity* lower, int lower_size,
                     const Entity* upper, int upper_size,
                     const Entity* boundary_map, int boundary_map_size,
                     const const_multiMap& face2node,
                     Entity f, const const_store<int>& node_remap){
  
  int numFaces = lower_size + upper_size + boundary_map_size;
  std::vector<Entity> faces(numFaces);
  
  //orient value: upper and boundary_map: 1
  //              lower: -1
  std::vector<char> orient(numFaces);
  
  int findex = 0;
  for(int i=0; i<lower_size; i++){
    faces[findex] = lower[i];
    orient[findex++] = 1;
  }
  
  for(int i=0; i<upper_size; i++){
    faces[findex] = upper[i];
    orient[findex++] = 0;
    
  }
  for(int i=0; i<boundary_map_size; i++){
    faces[findex] = boundary_map[i];
    orient[findex++] = 0;
    
  } 
  
  
  reorder_faces(node_remap, face2node, faces, &orient[0]);

  findex = -1;
  for( int i = 0; i < numFaces; i++){
    if(faces[i] == f){
      findex = i;
      break;
    }
  }
    
  
  if(findex == -1){
    cerr << "WARNING: can not find the face" << endl;
    Loci::Abort();
  }
  return findex;
}



//parallel version in merge_general_face.cc, set_general_face.cc
Cell* build_general_cell(const Entity* lower, int lower_size,
                         const Entity* upper, int upper_size,
                         const Entity* boundary_map, int boundary_map_size,
                         const const_multiMap& face2node,
                         const const_multiMap& face2edge,
                         const const_MapVec<2>& edge2node,
                         std::list<Node*>& bnode_list,
                         std::list<Edge*>& edge_list,
                         std::list<Face*>& face_list,
                         const const_store<int>& node_remap){
  int numFaces = lower_size + upper_size + boundary_map_size;
  std::vector<Entity> faces(numFaces);
  
  //orient value: upper and boundary_map: 1
  //              lower: -1
  char* orient = new char[numFaces];
  
  int findex = 0;
  for(int i=0; i<lower_size; i++){
    faces[findex] = lower[i];
    orient[findex++] = 1;
  }
  
  for(int i=0; i<upper_size; i++){
    faces[findex] = upper[i];
    orient[findex++] = 0;
    
  }
  for(int i=0; i<boundary_map_size; i++){
    faces[findex] = boundary_map[i];
    orient[findex++] = 0;
    
  } 
  
  
  //collect all the nodes and edges of the cell
  entitySet edges, nodes;
  for(unsigned int i=0; i< faces.size(); i++){
    for( int j= 0; j< face2node.num_elems(faces[i]); j++){
      nodes += face2node[faces[i]][j];
    }
    
    for( int j = 0; j < face2edge.num_elems(faces[i]); j++){
      edges += face2edge[faces[i]][j];
    }
  }
  //reorder entities

  std::vector<Entity> orderedNodes = reorder_nodes(node_remap, nodes);
  std::vector<Entity> orderedEdges = reorder_edges(node_remap, edge2node, edges);
  reorder_faces(node_remap, face2node, faces, orient);
  
  int numEdges = edges.size();
  int numNodes = nodes.size();
  
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(entitySet::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    Node* aNode = new Node();
    bnode_list.push_back(aNode);
    n2n[*np] = aNode;
  }
  
  //define each edge and put it into edge_list
  std::map<Entity, Edge*> e2e;
  for(entitySet::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    Edge* anEdge = new Edge(n2n[edge2node[*ep][0]], n2n[edge2node[*ep][1]]);
    edge_list.push_back(anEdge);
    e2e[*ep] = anEdge;
  }
          
  //defines each face and put it into face_list
  Face** face = new Face*[numFaces];
  for(int i  = 0; i < numFaces; i++){
    face[i] = new Face(face2edge.num_elems(faces[i]));
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < face[i]->numEdge; j++){
      face[i]->edge[j] = e2e[face2edge[faces[i]][j]];
      if(edge2node[face2edge[faces[i]][j]][0] == face2node[faces[i]][j] &&
         edge2node[face2edge[faces[i]][j]][1] == face2node[faces[i]][j==face2node.num_elems(faces[i])-1? 0:j+1]) face[i]->needReverse[j] = false;
            
      else face[i]->needReverse[j] = true;
    }
  }
        
  Node** node = new Node*[numNodes];

  for(unsigned int nindex = 0; nindex < orderedNodes.size(); nindex++){
    node[nindex] = n2n[orderedNodes[nindex]];
  }
        
  Edge** edge = new Edge*[numEdges];
        
  for(unsigned int eindex = 0; eindex < orderedEdges.size(); eindex++){
    edge[eindex] = e2e[orderedEdges[eindex]];
  }
  Cell* aCell = new Cell(numNodes, numEdges, numFaces, node, edge, face,orient);
       
  return aCell;
}

//parallel version in make_general_cellplan.cc
Cell* build_general_cell(const Entity* lower, int lower_size,
                         const Entity* upper, int upper_size,
                         const Entity* boundary_map, int boundary_map_size,
                         const const_store<bool>& is_quadface,
                         const const_multiMap& face2node,
                         const const_multiMap& face2edge,
                         const const_MapVec<2>& edge2node,
                         const const_store<vect3d>& pos,
                         const const_store<std::vector<char> >& edgePlan,
                         const const_store<std::vector<char> >& facePlan,
                         const const_store<char>& posTag,
                         const const_store<std::vector<char> >& nodeTag,
                         std::list<Node*>& bnode_list,
                         std::list<Edge*>& edge_list,
                         std::list<Face*>& face_list,
                         const const_store<int>& node_remap){
  int numFaces = lower_size + upper_size + boundary_map_size;
  std::vector<Entity> faces(numFaces);
  
  //orient value: upper and boundary_map: 1
  //              lower: -1
  char* orient = new char[numFaces];
  
  int findex = 0;
  for(int i=0; i<lower_size; i++){
    faces[findex] = lower[i];
    orient[findex++] = 1;
  }
    
  for(int i=0; i<upper_size; i++){
    faces[findex] = upper[i];
    orient[findex++] = 0;
      
  }
  for(int i=0; i<boundary_map_size; i++){
    faces[findex] = boundary_map[i];
    orient[findex++] = 0;
      
  } 
    
    
  //collect all the nodes and edges of the cell
  entitySet edges, nodes;
  for(unsigned int i=0; i< faces.size(); i++){
    for( int j= 0; j< face2node.num_elems(faces[i]); j++){
      nodes += face2node[faces[i]][j];
    }
      
    for( int j = 0; j < face2edge.num_elems(faces[i]); j++){
      edges += face2edge[faces[i]][j];
    }
  }
  //reorder entities
  std::vector<Entity> orderedNodes = reorder_nodes(node_remap, nodes);
  std::vector<Entity> orderedEdges = reorder_edges(node_remap, edge2node, edges);
  reorder_faces(node_remap, face2node, faces, orient);
  
  int numEdges = edges.size();
  int numNodes = nodes.size();
          
  //define each node
  std::map<Entity, Node*> n2n;
  for(entitySet::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    Node* aNode = new Node(pos[*np]);
    aNode->tag = posTag[*np];
    bnode_list.push_back(aNode);
    n2n[*np] = aNode;
  }
          
  //define each edge 
    
  std::map<Entity, Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
         
          
  for(entitySet::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    Edge* anEdge = new Edge(n2n[edge2node[*ep][0]], n2n[edge2node[*ep][1]]);
    edge_list.push_back(anEdge);
    e2e[*ep] = anEdge;
            
    //replit the edge
      
    anEdge->resplit(edgePlan[*ep], bnode_list);
      
    //index the node
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->tag =  nodeTag[*ep][nindex++];
      // checkNode(*np);
    }
    bnode_begin = --(bnode_list.end());
  }
    
    
    
  //defines each face 
  Face** face = new Face*[numFaces];
  for(int i  = 0; i < numFaces; i++){
    face[i] = new Face(face2edge.num_elems(faces[i]));
    face_list.push_back(face[i]);
    //find each edge
    for(int j = 0; j < face[i]->numEdge; j++){
      face[i]->edge[j] = e2e[face2edge[faces[i]][j]];
      if(edge2node[face2edge[faces[i]][j]][0] == face2node[faces[i]][j] &&
         edge2node[face2edge[faces[i]][j]][1] == face2node[faces[i]][j==face2node.num_elems(faces[i])-1? 0:j+1]) face[i]->needReverse[j] = false;
        
      else face[i]->needReverse[j] = true;
    }
      
    //replit each face
      
     
    if(is_quadface[faces[i]]){
      std::vector<char> newPlan = transfer_plan_q2g(facePlan[faces[i]]);
        
      face[i]->resplit(newPlan, bnode_list, edge_list);
      //to tag the node of inner nodes of face[i], first build a quadface with integer coordinates,
      //resplit it
      // tag the inner nodes of the quadface
      //then build a general face with integer coordinates
      //resplit it
      //find each inner node of quadface in the general face's inner node list, get its tag,
      //assign it to face[i]'s inner nodes

      std::list<Node*> tmp_bnode_list, tmp_node_list, tmp_node_list2;
      //build a quadface
      QuadFace* tmp_qface = build_tmp_quad_face(face2node[faces[i]].begin(),
                                                face2edge[faces[i]].begin(),
                                                edge2node,
                                                edgePlan,
                                                tmp_bnode_list,
                                                edge_list);
    
      //resplit the quadface to get node index
      tmp_qface->resplit(facePlan[faces[i]],
                         char(0),
                         tmp_node_list,
                         edge_list);
      
      int   nindex = 0;
      for(std::list<Node*>::const_iterator np = tmp_node_list.begin(); np!= tmp_node_list.end(); np++){
        (*np)->tag = nodeTag[faces[i]][nindex++];
      }

      Face* tmp_face =  build_tmp_general_face(face2node[faces[i]].begin(), face2node.num_elems(faces[i]),
                                               face2edge[faces[i]].begin(),
                                               edge2node,
                                               edgePlan,
                                               tmp_bnode_list,
                                               edge_list);
      


      
      tmp_face->resplit(newPlan, tmp_node_list2, edge_list);

      std::list<Node*>::const_iterator tmp_np2 = tmp_node_list2.begin();
      for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++, tmp_np2++){
        bool node_found = false;
        for(std::list<Node*>::const_iterator tmp_np = tmp_node_list.begin(); tmp_np!= tmp_node_list.end(); tmp_np++){
          if(int_equal((*tmp_np)->p, (*tmp_np2)->p)){
            (*np)->tag = (*tmp_np)->tag;
            node_found = true;
            break;
          }
        }
        if(!node_found){
          cerr <<"WARNING: can not find the node in build_general_cell"<< endl;
          Loci::Abort();
        }
      
      }
      //cleanup
      delete tmp_qface;
      delete tmp_face;
      cleanup_list(tmp_node_list);
      cleanup_list(tmp_node_list2);
      cleanup_list(tmp_bnode_list);

        
    }
    else{
      face[i]->resplit(facePlan[faces[i]], bnode_list, edge_list);
      //index the node
      int   nindex = 0;
      for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
        (*np)->tag = nodeTag[faces[i]][nindex++];
        //  checkNode(*np);
      } 
    }
     
      
     
    bnode_begin= --(bnode_list.end());
  }
    
  Node** node = new Node*[numNodes];
  for(unsigned int nindex = 0; nindex < orderedNodes.size(); nindex++){
    node[nindex] = n2n[orderedNodes[nindex]];
  }
        
  Edge** edge = new Edge*[numEdges];
        
  for(unsigned int eindex = 0; eindex < orderedEdges.size(); eindex++){
    edge[eindex] = e2e[orderedEdges[eindex]];
  }   
     
  return  new Cell(numNodes, numEdges, numFaces, node, edge, face,orient);
        
        
}
//build a cell with edgePlan and facePlan, tag the nodes
//then resplit the edges and faces with edgePlan1 and facePlan1
Cell* build_resplit_general_cell(const Entity* lower, int lower_size,
                                 const Entity* upper, int upper_size,
                                 const Entity* boundary_map, int boundary_map_size,
                                 const const_store<bool>& is_quadface,
                                 const const_multiMap& face2node,
                                 const const_multiMap& face2edge,
                                 const const_MapVec<2>& edge2node,
                                 const const_store<vect3d>& pos,
                                 const const_store<std::vector<char> >& edgePlan,
                                 const const_store<std::vector<char> >& facePlan,
                                 const const_store<std::vector<char> >& edgePlan1,
                                 const const_store<std::vector<char> >& facePlan1,    
                                 const const_store<char>& posTag,
                                 const const_store<std::vector<char> >& nodeTag,
                                 std::list<Node*>& bnode_list,
                                 std::list<Node*>& node_list,
                                 std::list<Edge*>& edge_list,
                                 std::list<Face*>& face_list,
                                 const const_store<int>& node_remap,
                                 const std::vector<char>& cellPlan,
                                 const  std::vector<char>& cellNodeTag ){
  int numFaces = lower_size + upper_size + boundary_map_size;
  std::vector<Entity> faces(numFaces);
  
  //orient value: upper and boundary_map: 1
  //              lower: -1
  char* orient = new char[numFaces];
  
  int findex = 0;
  for(int i=0; i<lower_size; i++){
    faces[findex] = lower[i];
    orient[findex++] = 1;
  }
    
  for(int i=0; i<upper_size; i++){
    faces[findex] = upper[i];
    orient[findex++] = 0;
      
  }
  for(int i=0; i<boundary_map_size; i++){
    faces[findex] = boundary_map[i];
    orient[findex++] = 0;
      
  } 
    
    
  //collect all the nodes and edges of the cell
  entitySet edges, nodes;
  for(unsigned int i=0; i< faces.size(); i++){
    for( int j= 0; j< face2node.num_elems(faces[i]); j++){
      nodes += face2node[faces[i]][j];
    }
      
    for( int j = 0; j < face2edge.num_elems(faces[i]); j++){
      edges += face2edge[faces[i]][j];
    }
  }
  //reorder entities
  std::vector<Entity> orderedNodes = reorder_nodes(node_remap, nodes);
  std::vector<Entity> orderedEdges = reorder_edges(node_remap, edge2node, edges);
  reorder_faces(node_remap, face2node, faces, orient);
  
  int numEdges = edges.size();
  int numNodes = nodes.size();
          
  //define each node
  std::map<Entity, Node*> n2n;
  for(entitySet::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    Node* aNode = new Node(pos[*np]);
    aNode->tag = posTag[*np];
    bnode_list.push_back(aNode);
    n2n[*np] = aNode;
  }
          
  //define each edge 
    
  std::map<Entity, Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
         
          
  for(entitySet::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    Edge* anEdge = new Edge(n2n[edge2node[*ep][0]], n2n[edge2node[*ep][1]]);
    edge_list.push_back(anEdge);
    e2e[*ep] = anEdge;
            
    //replit the edge
      
    anEdge->resplit(edgePlan[*ep], bnode_list);
      
    //index the node
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->tag =  nodeTag[*ep][nindex++];
      // checkNode(*np);
    }
    bnode_begin = --(bnode_list.end());
  }
    
    
    
  //defines each face 
  Face** face = new Face*[numFaces];
  for(int i  = 0; i < numFaces; i++){
    face[i] = new Face(face2edge.num_elems(faces[i]));
    face_list.push_back(face[i]);
    //find each edge
    for(int j = 0; j < face[i]->numEdge; j++){
      face[i]->edge[j] = e2e[face2edge[faces[i]][j]];
      if(edge2node[face2edge[faces[i]][j]][0] == face2node[faces[i]][j] &&
         edge2node[face2edge[faces[i]][j]][1] == face2node[faces[i]][j==face2node.num_elems(faces[i])-1? 0:j+1]) face[i]->needReverse[j] = false;
        
      else face[i]->needReverse[j] = true;
    }
      
    //replit each face
      
     
    if(is_quadface[faces[i]]){
      std::vector<char> newPlan = transfer_plan_q2g(facePlan[faces[i]]);
        
      face[i]->resplit(newPlan, bnode_list, edge_list);
      //to tag the node of inner nodes of face[i], first build a quadface with integer coordinates,
      //resplit it
      // tag the inner nodes of the quadface
      //then build a general face with integer coordinates
      //resplit it
      //find each inner node of quadface in the general face's inner node list, get its tag,
      //assign it to face[i]'s inner nodes

      std::list<Node*> tmp_bnode_list, tmp_node_list, tmp_node_list2;
      //build a quadface
      QuadFace* tmp_qface = build_tmp_quad_face(face2node[faces[i]].begin(),
                                                face2edge[faces[i]].begin(),
                                                edge2node,
                                                edgePlan,
                                                tmp_bnode_list,
                                                edge_list);
    
      //resplit the quadface to get node index
      tmp_qface->resplit(facePlan[faces[i]],
                         char(0),
                         tmp_node_list,
                         edge_list);
      
      int   nindex = 0;
      for(std::list<Node*>::const_iterator np = tmp_node_list.begin(); np!= tmp_node_list.end(); np++){
        (*np)->tag = nodeTag[faces[i]][nindex++];
      }

      Face* tmp_face =  build_tmp_general_face(face2node[faces[i]].begin(), face2node.num_elems(faces[i]),
                                               face2edge[faces[i]].begin(),
                                               edge2node,
                                               edgePlan,
                                               tmp_bnode_list,
                                               edge_list);
      


      
      tmp_face->resplit(newPlan, tmp_node_list2, edge_list);

      std::list<Node*>::const_iterator tmp_np2 = tmp_node_list2.begin();
      for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++, tmp_np2++){
        bool node_found = false;
        for(std::list<Node*>::const_iterator tmp_np = tmp_node_list.begin(); tmp_np!= tmp_node_list.end(); tmp_np++){
          if(int_equal((*tmp_np)->p, (*tmp_np2)->p)){
            (*np)->tag = (*tmp_np)->tag;
            node_found = true;
            break;
          }
        }
        if(!node_found){
          cerr <<"WARNING: can not find the node in build_general_cell"<< endl;
          Loci::Abort();
        }
      
      }
      //cleanup
      delete tmp_qface;
      delete tmp_face;
      cleanup_list(tmp_node_list);
      cleanup_list(tmp_node_list2);
      cleanup_list(tmp_bnode_list);

        
    }
    else{
      face[i]->resplit(facePlan[faces[i]], bnode_list, edge_list);
      //index the node
      int   nindex = 0;
      for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
        (*np)->tag = nodeTag[faces[i]][nindex++];
        //  checkNode(*np);
      } 
    }
     
      
     
    bnode_begin= --(bnode_list.end());
  }
    
  Node** node = new Node*[numNodes];
  for(unsigned int nindex = 0; nindex < orderedNodes.size(); nindex++){
    node[nindex] = n2n[orderedNodes[nindex]];
  }
        
  Edge** edge = new Edge*[numEdges];
        
  for(unsigned int eindex = 0; eindex < orderedEdges.size(); eindex++){
    edge[eindex] = e2e[orderedEdges[eindex]];
  }
  Cell* aCell = new Cell(numNodes, numEdges, numFaces, node, edge, face,orient);
  
  //finish build
  
  //resplit cell
  std::vector<DiamondCell*> cells;
  aCell->resplit( cellPlan, 
                  node_list,
                  edge_list,
                  face_list,
                  cells);
#ifdef SIZE_DEBUG     
  if(node_list.size()!= cellNodeTag.size()){
    cerr<< " nodeTag size and node_list size mismatch(), nodeTag: " << cellNodeTag.size() << " node_list " << node_list.size() <<endl;
    Loci::Abort();
  }
#endif
  //tag nodes
  int nindex = 0;
  for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++){
    (*np)->tag = cellNodeTag[nindex++];
  }
  cells.clear();
  
  
  
  //resplit the edges again without tagging the node
  for(entitySet::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    e2e[*ep]->resplit(edgePlan1[*ep], bnode_list);
  }
  //resplit the faces again without tagging the node
  for(int i  = 0; i < numFaces; i++){
    if(is_quadface[faces[i]]){
      std::vector<char> newPlan = transfer_plan_q2g(facePlan1[faces[i]]);
      face[i]->resplit(newPlan, bnode_list, edge_list);
    }else{
      face[i]->resplit(facePlan1[faces[i]], bnode_list, edge_list);
    }
  }

  
  return aCell;  
        
        
}





//build a cell with edgePlan,facePlan and cellPlan, tag the cells
//then resplit the edges and faces with edgePlan1 and facePlan1
Cell* build_resplit_general_cell_ctag(const Entity* lower, int lower_size,
                                      const Entity* upper, int upper_size,
                                      const Entity* boundary_map, int boundary_map_size,
                                      const const_store<bool>& is_quadface,
                                      const const_multiMap& face2node,
                                      const const_multiMap& face2edge,
                                      const const_MapVec<2>& edge2node,
                                      const const_store<vect3d>& pos,
                                      const const_store<std::vector<char> >& edgePlan,
                                      const const_store<std::vector<char> >& facePlan,
                                      const const_store<std::vector<char> >& edgePlan1,
                                      const const_store<std::vector<char> >& facePlan1,    
                                      std::list<Node*>& bnode_list,
                                      std::list<Node*>& node_list,
                                      std::list<Edge*>& edge_list,
                                      std::list<Face*>& face_list,
                                      const const_store<int>& node_remap,
                                      const std::vector<char>& cellPlan,
                                      const  std::vector<char>& fineCellTag){



  int numFaces = lower_size + upper_size + boundary_map_size;
  std::vector<Entity> faces(numFaces);
  
  //orient value: upper and boundary_map: 1
  //              lower: -1
  char* orient = new char[numFaces];
  
  int findex = 0;
  for(int i=0; i<lower_size; i++){
    faces[findex] = lower[i];
    orient[findex++] = 1;
  }
    
  for(int i=0; i<upper_size; i++){
    faces[findex] = upper[i];
    orient[findex++] = 0;
      
  }
  for(int i=0; i<boundary_map_size; i++){
    faces[findex] = boundary_map[i];
    orient[findex++] = 0;
      
  } 
    
    
  //collect all the nodes and edges of the cell
  entitySet edges, nodes;
  for(unsigned int i=0; i< faces.size(); i++){
    for( int j= 0; j< face2node.num_elems(faces[i]); j++){
      nodes += face2node[faces[i]][j];
    }
      
    for( int j = 0; j < face2edge.num_elems(faces[i]); j++){
      edges += face2edge[faces[i]][j];
    }
  }
  //reorder entities
  std::vector<Entity> orderedNodes = reorder_nodes(node_remap, nodes);
  std::vector<Entity> orderedEdges = reorder_edges(node_remap, edge2node, edges);
  reorder_faces(node_remap, face2node, faces, orient);
  
  int numEdges = edges.size();
  int numNodes = nodes.size();
          
  //define each node
  std::map<Entity, Node*> n2n;
  for(entitySet::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    Node* aNode = new Node(pos[*np]);
    bnode_list.push_back(aNode);
    n2n[*np] = aNode;
  }
          
  //define each edge 
    
  std::map<Entity, Edge*> e2e;
           
  for(entitySet::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    Edge* anEdge = new Edge(n2n[edge2node[*ep][0]], n2n[edge2node[*ep][1]]);
    edge_list.push_back(anEdge);
    e2e[*ep] = anEdge;
            
    //replit the edge
      
    anEdge->resplit(edgePlan[*ep], bnode_list);
  }
    
    
    
  //defines each face 
  Face** face = new Face*[numFaces];
  for(int i  = 0; i < numFaces; i++){
    face[i] = new Face(face2edge.num_elems(faces[i]));
    face_list.push_back(face[i]);
    //find each edge
    for(int j = 0; j < face[i]->numEdge; j++){
      face[i]->edge[j] = e2e[face2edge[faces[i]][j]];
      if(edge2node[face2edge[faces[i]][j]][0] == face2node[faces[i]][j] &&
         edge2node[face2edge[faces[i]][j]][1] == face2node[faces[i]][j==face2node.num_elems(faces[i])-1? 0:j+1]) face[i]->needReverse[j] = false;
        
      else face[i]->needReverse[j] = true;
    }
      
    //replit each face
      
     
    if(is_quadface[faces[i]]){
      std::vector<char> newPlan = transfer_plan_q2g(facePlan[faces[i]]);
        
      face[i]->resplit(newPlan, bnode_list, edge_list);
      
    }
    else{
      face[i]->resplit(facePlan[faces[i]], bnode_list, edge_list);
     
    }
     
  }
    
  Node** node = new Node*[numNodes];
  for(unsigned int nindex = 0; nindex < orderedNodes.size(); nindex++){
    node[nindex] = n2n[orderedNodes[nindex]];
  }
        
  Edge** edge = new Edge*[numEdges];
        
  for(unsigned int eindex = 0; eindex < orderedEdges.size(); eindex++){
    edge[eindex] = e2e[orderedEdges[eindex]];
  }
  Cell* aCell = new Cell(numNodes, numEdges, numFaces, node, edge, face,orient);
  
  //finish build
  
  //resplit cell
  std::vector<DiamondCell*> cells;
  aCell->resplit( cellPlan, 
                  node_list,
                  edge_list,
                  face_list,
                  cells);
#ifdef SIZE_DEBUG     
  if(cells.size()!= fineCellTag.size()){
    cerr<< " fineCellTag size and cells size mismatch(), fineCellTag: " << fineCellTag.size() << " cells " << cells.size() <<endl;
    Loci::Abort();
  }
#endif
  //tag nodes
  int cindex = 0;
  for(std::vector<DiamondCell*>::const_iterator np = cells.begin(); np!= cells.end(); np++){
    (*np)->setTag(fineCellTag[cindex++]);
  }
  cells.clear();
  
  
  
  //resplit the edges again without tagging the node
  for(entitySet::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    e2e[*ep]->resplit(edgePlan1[*ep], bnode_list);
  }
  //resplit the faces again without tagging the node
  for(int i  = 0; i < numFaces; i++){
    if(is_quadface[faces[i]]){
      std::vector<char> newPlan = transfer_plan_q2g(facePlan1[faces[i]]);
      face[i]->resplit(newPlan, bnode_list, edge_list);
    }else{
      face[i]->resplit(facePlan1[faces[i]], bnode_list, edge_list);
    }
  }

  
  return aCell;  
        
        
}



//parallel version in make_general_cellplan.cc
Cell* build_general_cell(const Entity* lower, int lower_size,
                         const Entity* upper, int upper_size,
                         const Entity* boundary_map, int boundary_map_size,
                         const const_multiMap& face2node,
                         const const_multiMap& face2edge,
                         const const_MapVec<2>& edge2node,
                         const const_store<vect3d>& pos,
                         const const_store<char>& posTag,
                         std::list<Node*>& bnode_list,
                         std::list<Edge*>& edge_list,
                         std::list<Face*>& face_list,
                         const const_store<int>& node_remap){
  int numFaces = lower_size + upper_size + boundary_map_size;
  std::vector<Entity> faces(numFaces);
  
  //orient value: upper and boundary_map: 1
  //              lower: -1
  char* orient = new char[numFaces];
  
  int findex = 0;
  for(int i=0; i<lower_size; i++){
    faces[findex] = lower[i];
    orient[findex++] = 1;
  }
    
  for(int i=0; i<upper_size; i++){
    faces[findex] = upper[i];
    orient[findex++] = 0;
      
  }
  for(int i=0; i<boundary_map_size; i++){
    faces[findex] = boundary_map[i];
    orient[findex++] = 0;
      
  } 
    
    
  //collect all the nodes and edges of the cell
  entitySet edges, nodes;
  for(unsigned int i=0; i< faces.size(); i++){
    for( int j= 0; j< face2node.num_elems(faces[i]); j++){
      nodes += face2node[faces[i]][j];
    }
      
    for( int j = 0; j < face2edge.num_elems(faces[i]); j++){
      edges += face2edge[faces[i]][j];
    }
  }
  //reorder entities
  std::vector<Entity> orderedNodes = reorder_nodes(node_remap, nodes);
  std::vector<Entity> orderedEdges = reorder_edges(node_remap, edge2node, edges);
  reorder_faces(node_remap, face2node, faces, orient);
  
  int numEdges = edges.size();
  int numNodes = nodes.size();
          
  //define each node
  std::map<Entity, Node*> n2n;
  for(entitySet::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    Node* aNode = new Node(pos[*np]);
    aNode->tag = posTag[*np];
    bnode_list.push_back(aNode);
    n2n[*np] = aNode;
  }
          
  //define each edge 
    
  std::map<Entity, Edge*> e2e;

         
          
  for(entitySet::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    Edge* anEdge = new Edge(n2n[edge2node[*ep][0]], n2n[edge2node[*ep][1]]);
    edge_list.push_back(anEdge);
    e2e[*ep] = anEdge;
  }

    
    
    
  //defines each face 
  Face** face = new Face*[numFaces];
  for(int i  = 0; i < numFaces; i++){
    face[i] = new Face(face2edge.num_elems(faces[i]));
    face_list.push_back(face[i]);
    //find each edge
    for(int j = 0; j < face[i]->numEdge; j++){
      face[i]->edge[j] = e2e[face2edge[faces[i]][j]];
      if(edge2node[face2edge[faces[i]][j]][0] == face2node[faces[i]][j] &&
         edge2node[face2edge[faces[i]][j]][1] == face2node[faces[i]][j==face2node.num_elems(faces[i])-1? 0:j+1]) face[i]->needReverse[j] = false;
        
      else face[i]->needReverse[j] = true;
    }
  }
    
  Node** node = new Node*[numNodes];
  for(unsigned int nindex = 0; nindex < orderedNodes.size(); nindex++){
    node[nindex] = n2n[orderedNodes[nindex]];
  }
    
  Edge** edge = new Edge*[numEdges];
    
  for(unsigned int eindex = 0; eindex < orderedEdges.size(); eindex++){
    edge[eindex] = e2e[orderedEdges[eindex]];
  }   
    
  return  new Cell(numNodes, numEdges, numFaces, node, edge, face,orient);
    
    
}
//parallel version in make_general_cellplan.cc
Cell* build_general_cell(const Entity* lower, int lower_size,
                         const Entity* upper, int upper_size,
                         const Entity* boundary_map, int boundary_map_size,
                         const const_multiMap& face2node,
                         const const_multiMap& face2edge,
                         const const_MapVec<2>& edge2node,
                         const const_store<vect3d>& pos,
                         std::list<Node*>& bnode_list,
                         std::list<Edge*>& edge_list,
                         std::list<Face*>& face_list,
                         const const_store<int>& node_remap){
  int numFaces = lower_size + upper_size + boundary_map_size;
  std::vector<Entity> faces(numFaces);
  
  //orient value: upper and boundary_map: 1
  //              lower: -1
  char* orient = new char[numFaces];
  
  int findex = 0;
  for(int i=0; i<lower_size; i++){
    faces[findex] = lower[i];
    orient[findex++] = 1;
  }
    
  for(int i=0; i<upper_size; i++){
    faces[findex] = upper[i];
    orient[findex++] = 0;
      
  }
  for(int i=0; i<boundary_map_size; i++){
    faces[findex] = boundary_map[i];
    orient[findex++] = 0;
      
  } 
    
    
  //collect all the nodes and edges of the cell
  entitySet edges, nodes;
  for(unsigned int i=0; i< faces.size(); i++){
    for( int j= 0; j< face2node.num_elems(faces[i]); j++){
      nodes += face2node[faces[i]][j];
    }
      
    for( int j = 0; j < face2edge.num_elems(faces[i]); j++){
      edges += face2edge[faces[i]][j];
    }
  }
  //reorder entities
  std::vector<Entity> orderedNodes = reorder_nodes(node_remap, nodes);
    
  std::vector<Entity> orderedEdges = reorder_edges(node_remap, edge2node, edges);
  reorder_faces(node_remap, face2node, faces, orient);
  
  int numEdges = edges.size();
  int numNodes = nodes.size();
          
  //define each node
  std::map<Entity, Node*> n2n;
  for(entitySet::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    Node* aNode = new Node(pos[*np]);
    bnode_list.push_back(aNode);
    n2n[*np] = aNode;
  }
          
  //define each edge 
    
  std::map<Entity, Edge*> e2e;

         
          
  for(entitySet::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    Edge* anEdge = new Edge(n2n[edge2node[*ep][0]], n2n[edge2node[*ep][1]]);
    edge_list.push_back(anEdge);
    e2e[*ep] = anEdge;
  }

    
    
    
  //defines each face 
  Face** face = new Face*[numFaces];
  for(int i  = 0; i < numFaces; i++){
    face[i] = new Face(face2edge.num_elems(faces[i]));
    face_list.push_back(face[i]);
    //find each edge
    for(int j = 0; j < face[i]->numEdge; j++){
      face[i]->edge[j] = e2e[face2edge[faces[i]][j]];
      if(edge2node[face2edge[faces[i]][j]][0] == face2node[faces[i]][j] &&
         edge2node[face2edge[faces[i]][j]][1] == face2node[faces[i]][j==face2node.num_elems(faces[i])-1? 0:j+1]) face[i]->needReverse[j] = false;
        
      else face[i]->needReverse[j] = true;
    }
  }
    
  Node** node = new Node*[numNodes];
  for(unsigned int nindex = 0; nindex < orderedNodes.size(); nindex++){
    node[nindex] = n2n[orderedNodes[nindex]];
  }
    
  Edge** edge = new Edge*[numEdges];
    
  for(unsigned int eindex = 0; eindex < orderedEdges.size(); eindex++){
    edge[eindex] = e2e[orderedEdges[eindex]];
  }   
    
  return  new Cell(numNodes, numEdges, numFaces, node, edge, face,orient);
    
    
}

//parallel version, build  a cell and index all boundary nodes
Cell* build_general_cell(const Entity* lower, int lower_size,
                         const Entity* upper, int upper_size,
                         const Entity* boundary_map, int boundary_map_size,
                         const const_store<bool>& is_quadface,
                         const_multiMap& face2node,
                         const_multiMap& face2edge,
                         const_MapVec<2>& edge2node,
                         const_store<vect3d>& pos,
                         const_store<std::vector<char> >& edgePlan,
                         const_store<std::vector<char> >& facePlan,
                         const_store<int>& node_offset,
                         std::list<Node*>& bnode_list,
                         std::list<Edge*>& edge_list,
                         std::list<Face*>& face_list,
                         const const_store<int>& node_remap
                         ){

  int numFaces = lower_size + upper_size + boundary_map_size;
  std::vector<Entity> faces(numFaces);
  
  //orient value: upper and boundary_map: 1
  //              lower: -1
  char* orient = new char[numFaces];
  
  int findex = 0;
  for(int i=0; i<lower_size; i++){
    faces[findex] = lower[i];
    orient[findex++] = 1;
  }
  
  for(int i=0; i<upper_size; i++){
    faces[findex] = upper[i];
    orient[findex++] = 0;
  }
  for(int i=0; i<boundary_map_size; i++){
    faces[findex] = boundary_map[i];
    orient[findex++] = 0;
    
  } 
  
          
  //collect all the nodes and edges of the cell
  entitySet edges, nodes;
  for(unsigned int i=0; i< faces.size(); i++){
    for( int j= 0; j< face2node.num_elems(faces[i]); j++){
      nodes += face2node[faces[i]][j];
    }
    
    for( int j = 0; j < face2edge.num_elems(faces[i]); j++){
      edges += face2edge[faces[i]][j];
    }
  }
  //reorder entities
  std::vector<Entity> orderedNodes = reorder_nodes(node_remap, nodes);
  std::vector<Entity> orderedEdges = reorder_edges(node_remap, edge2node, edges);
  reorder_faces(node_remap, face2node, faces, orient);
  
  int numEdges = edges.size();
  int numNodes = nodes.size();


  
  
  
  
  //define each node
  std::map<Entity, Node*> n2n;
  for(entitySet::const_iterator np = nodes.begin(); np != nodes.end(); np++){
    // if(node_remap[*np] < 1) cout << node_remap[*np] << endl;
    Node* aNode = new Node(pos[*np], (node_remap[*np]));
    bnode_list.push_back(aNode);
    n2n[*np] = aNode;
  }
  
  //define each edge 
  std::map<Entity, Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
  
  
  for(entitySet::const_iterator ep = edges.begin(); ep != edges.end(); ep++){
    Edge* anEdge = new Edge(n2n[edge2node[*ep][0]], n2n[edge2node[*ep][1]]);
    edge_list.push_back(anEdge);
    e2e[*ep] = anEdge;
    
    //replit the edge
    anEdge->resplit(edgePlan[*ep],bnode_list);
    
    //index the node
    int nindex = node_offset[*ep];
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->index =  nindex++;
     
    }
    bnode_begin = --(bnode_list.end());
  }
  
  
  
  //defines each face 
  Face** face = new Face*[numFaces];
  for(int i  = 0; i < numFaces; i++){
    face[i] = new Face(face2edge.num_elems(faces[i]));
    face_list.push_back(face[i]);
    //find each edge
    for(int j = 0; j < face[i]->numEdge; j++){
      face[i]->edge[j] = e2e[face2edge[faces[i]][j]];
      if(edge2node[face2edge[faces[i]][j]][0] == face2node[faces[i]][j] &&
         edge2node[face2edge[faces[i]][j]][1] == face2node[faces[i]]
         [j==face2node.num_elems(faces[i])-1? 0:j+1]) face[i]->needReverse[j] = false;
      
      else face[i]->needReverse[j] = true;
    }
    
    //replit each face
    
   
    if(is_quadface[faces[i]]){
      //first transfer quadface plan to face plan
      std::vector<char> newPlan = transfer_plan_q2g(facePlan[faces[i]]);
      //resplit face[i]
      face[i]->resplit(newPlan, bnode_list, edge_list);
      //to index the node of inner nodes of face[i], first build a quadface with integer coordinates,
      //resplit it
      // index the inner nodes of the quadface
      //then build a general face with integer coordinates
      //resplit it
      //find each inner node of quadface in the general face's inner node list, get its index,
      //assign it to face[i]'s inner nodes

      std::list<Node*> tmp_bnode_list, tmp_node_list, tmp_node_list2;
      //build a quadface
      QuadFace* tmp_qface = build_tmp_quad_face(face2node[faces[i]].begin(),
                                                face2edge[faces[i]].begin(),
                                                edge2node,
                                                edgePlan,
                                                tmp_bnode_list,
                                                edge_list);
   
      //resplit the quadface to get node index
      tmp_qface->resplit(facePlan[faces[i]],
                         char(0),
                         tmp_node_list,
                         edge_list);
      
      int   nindex = node_offset[faces[i]];
      for(std::list<Node*>::const_iterator np = tmp_node_list.begin(); np!= tmp_node_list.end(); np++){
        (*np)->index = nindex++;
      }

      Face* tmp_face =  build_tmp_general_face(face2node[faces[i]].begin(), face2node.num_elems(faces[i]),
                                               face2edge[faces[i]].begin(),
                                               edge2node,
                                               edgePlan,
                                               tmp_bnode_list,
                                               edge_list);
      


      
      tmp_face->resplit(newPlan, tmp_node_list2, edge_list);

      std::list<Node*>::const_iterator tmp_np2 = tmp_node_list2.begin();
      for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++, tmp_np2++){
        bool node_found = false;
        for(std::list<Node*>::const_iterator tmp_np = tmp_node_list.begin(); tmp_np!= tmp_node_list.end(); tmp_np++){
          if(int_equal((*tmp_np)->p, (*tmp_np2)->p)){
            (*np)->index = (*tmp_np)->index;
            node_found = true;
            break;
          }
        }
        if(!node_found){
          cerr <<"WARNING: can not find the node in build_general_cell"<< endl;
          Loci::Abort();
        }
      
      }
      //cleanup
      delete tmp_qface;
      delete tmp_face;
      cleanup_list(tmp_node_list);
      cleanup_list(tmp_node_list2);
      cleanup_list(tmp_bnode_list);
    }
    else{
      face[i]->resplit(facePlan[faces[i]], bnode_list, edge_list);
      //index the node
      int   nindex = node_offset[faces[i]];
      for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
        (*np)->index = nindex++;
     
      }
    }
   
    
   
    bnode_begin= --(bnode_list.end());
  }
  
  Node** node = new Node*[numNodes];
  
  for(unsigned int nindex = 0; nindex < orderedNodes.size(); nindex++){
    node[nindex] = n2n[orderedNodes[nindex]];
  }
        
  Edge** edge = new Edge*[numEdges];
        
  for(unsigned int eindex = 0; eindex < orderedEdges.size(); eindex++){
    edge[eindex] = e2e[orderedEdges[eindex]];
  }
  return new Cell(numNodes, numEdges, numFaces, node, edge, face,orient);
}

