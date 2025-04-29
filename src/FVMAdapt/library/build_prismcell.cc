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
// This file contain funstions that build a prism according to Loci data structures
// and refinement plans. When build a prism, edges are built as local edges in cell
// and resplit with needReverse. both general faces and quadfaces faces are built as
// local faces in cell and resplit with orientCode. Because directions of local edges
// in cell are consistent with faces, the needReverse of general faces are always false.
// also the faceOrient of quadfaces are always reset(point outward)

#include <vector>
#include <list>
#include <Loci.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <Tools/tools.h>
#include "prism.h"


using std::cerr;
using std::endl;
using std::swap;
using std::cout;
using Loci::storeRepP;
using std::vector;
void reorder_faces(const const_store<int>& node_remap, std::vector<Entity>& lower,
                   std::vector<Entity>& upper,
                   std::vector<Entity>& boundary_map);
//this function defines the 5 faces of prism from loci data structures
Array<Entity, 5> collect_prism_faces( const Entity* lower, int lower_size,
                                      const Entity* upper, int upper_size,
                                      const Entity* boundary_map,int boundary_map_size,
                                      const Array<char, 5>& prism2face, const const_store<int>& node_remap ){
  //first create vectors for reordering
  vector<Entity> vlower(lower_size);
  vector<Entity> vupper(upper_size);
  vector<Entity> vboundary_map(boundary_map_size);
  int nf =0;
  for (int f=0; f<lower_size; f++) vlower[nf++] =lower[f]; 
  nf =0;
  for (int f=0; f<upper_size; f++) vupper[nf++] =upper[f]; 
  nf =0;
  for (int f=0; f<boundary_map_size; f++) vboundary_map[nf++] =boundary_map[f]; 
  reorder_faces(node_remap, vlower, vupper, vboundary_map);

  
  // Collect the entity designation of all faces of the cell cc in the all_faces array
  Array<Entity, 5> faces;
  for (int f=0; f<5; f++) {
    switch (prism2face[f]/6) {
    case 0:
      faces[f] = vlower[prism2face[f]%6];
      break;
    case 1:
      faces[f] = vupper[prism2face[f]%6];
      break;
    case 2:
      faces[f] = vboundary_map[prism2face[f]%6];
      break;
    default:
      cerr << " WARNING: illegal prism2face value" << endl;
      break;
    }
  }
   
  return faces;
}

//this function defines the 6 vertices of prism from loci data structures
Array<Entity, 6> collect_prism_vertices(const const_multiMap& face2node,
                                        const Array<Entity, 5>& faces,
                                        const Array<char, 6>& prism2node){
  Array<Entity, 6>  prism_vertices;
  for(int i = 0; i < 3; i++){
    prism_vertices[i] = face2node[faces[0]][prism2node[i]];
  }
  for(int i = 0; i < 3; i++){
    prism_vertices[i+3] = face2node[faces[1]][prism2node[i+3]];
  }
  return prism_vertices;
}




//collect the entity designation of all edges of the prism,
//the directions of edges are stored in needReverse 
Array<Entity, 9> collect_prism_edges(const Array<Entity, 5>& faces, const Array<Entity,6>& prism_vertices,
                                     const const_multiMap& face2edge, const const_MapVec<2>& edge2node,
                                     Array<bool, 9>& needReverse){
  entitySet all_edges;
  Array<Entity, 9> prism_edges;
  for(int fID = 0; fID < 2; fID++){
    for(int eID = 0; eID < 3; eID++){
      all_edges += face2edge[faces[fID]][eID];
    }
  }
  for(int fID = 2; fID < 5; fID++){
    for(int eID = 0; eID < 4; eID++){
      all_edges += face2edge[faces[fID]][eID];
    }
  }
  
  
  if(all_edges.size() != 9){
    cerr << "WARNING: the number of edges is not 12" << endl;
    Loci:: Abort();
  }
  Entity head, tail;
  int node0[9] = {0, 1, 2, 3, 4, 5, 0, 1, 2};
  int node1[9] = {1, 2, 0, 4, 5, 3, 3, 4, 5};
  
  for( int e = 0; e < 9; e++ ){
       
    head = prism_vertices[node0[e]];
    tail = prism_vertices[node1[e]];
    entitySet::const_iterator ee;
    for(ee = all_edges.begin(); ee != all_edges.end(); ee++){
      if(edge2node[*ee][0] == head && edge2node[*ee][1] == tail){
        prism_edges[e] = *ee;
        needReverse[e] = false;
        break;
      }
      else if(edge2node[*ee][1] == head && edge2node[*ee][0] == tail){
        prism_edges[e] = *ee;
        needReverse[e] = true;
        break;
      }
      
    }
    if(ee ==  all_edges.end()){
      cerr << "WARNING: can not find edge" << endl;
      Loci::Abort();
    }
  }
  
  return prism_edges;
}


Prism* build_prism_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,5>& prism2face,
                        const Array<char,6>& prism2node,
                        const Array<char,5>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        const const_store<std::vector<char> >& edgePlan,
                        const const_store<std::vector<char> >& facePlan,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& qface_list,
                        std::list<Face*>& gface_list,
                        const const_store<int>& node_remap){
  
  Array<Entity, 5> face_entity = collect_prism_faces(lower, lower_size,
                                                     upper, upper_size,
                                                     boundary_map,boundary_map_size,
                                                     prism2face, node_remap);
  
  Array<Entity, 6> node_entity = collect_prism_vertices(face2node,
                                                        face_entity,
                                                        prism2node);



  Array<bool, 9> edge_reverse;
  Array<Entity, 9> edge_entity = collect_prism_edges( face_entity,
                                                      node_entity,
                                                      face2edge,
                                                      edge2node,
                                                      edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 6; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  //edge is built according to the direction inside the prism
  //and resplit with needReverse
  std::map<Entity,Edge*> e2e;
  for(int i = 0; i < 9; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;
    
    
    //resplit the edge
    anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);
  }
  
  int f2e[3][4]= {{0, 7, 3, 6}, {1, 8, 4, 7}, {2, 6, 5, 8}};
  int gf2e[2][3] = {{0, 1, 2}, {3, 4, 5}};
  
  Face* gface;
  QuadFace* qface;
  //defines each face and put it into face_list
  Prism* aCell = new Prism(3);
  
  //gnrlface[0]
  for(int i = 0; i < 2; i++){
    gface = new Face(3);
    gface_list.push_back(gface);

    //define each edge
    for(int j = 0; j < 3; j++){
      gface->edge[j] = e2e[edge_entity[gf2e[i][j]]];
      gface->needReverse[j] = false; 
    }
    //resplit each face
    gface->resplit(facePlan[face_entity[i]], orientCode[i], bnode_list, edge_list);//resplit without orientCode
    aCell->setFace(i, gface);
  }
  
  //quadface
  for(int i = 0; i < 3; i++){
    qface = new QuadFace(4);
    qface_list.push_back(qface);
      
    //define each edge
    for(int j = 0; j < 4; j++){
      qface->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    qface->resplit(facePlan[face_entity[i+2]],orientCode[i+2], bnode_list, edge_list);
    aCell->setFace(i, qface);
  }
  
  return aCell;
}




  
Prism* build_prism_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,5>& prism2face,
                        const Array<char,6>& prism2node,
                        const Array<char,5>& orientCode,
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
                        std::list<QuadFace*>& qface_list,
                        std::list<Face*>& gface_list,
                        const const_store<int>& node_remap){
  
  Array<Entity, 5> face_entity = collect_prism_faces(lower, lower_size,
                                                     upper,upper_size,
                                                     boundary_map,boundary_map_size,
                                                     prism2face, node_remap);
  
  Array<Entity, 6> node_entity = collect_prism_vertices(face2node,
                                                        face_entity,
                                                        prism2node);



  Array<bool, 9> edge_reverse;
  Array<Entity, 9> edge_entity = collect_prism_edges( face_entity,
                                                      node_entity,
                                                      face2edge,
                                                      edge2node,
                                                      edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 6; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    aNode->tag = posTag[node_entity[i]];
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  //edge is built according to the direction inside the prism
  //and resplit with needReverse
  std::map<Entity,Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
  
  for(int i = 0; i < 9; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;
    
    
    //resplit the edge
    anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->tag =  nodeTag[edge_entity[i]][nindex++];
      // checkNode(*np);
    }
    bnode_begin = --(bnode_list.end());
  }
  
  int f2e[3][4]= {{0, 7, 3, 6}, {1, 8, 4, 7}, {2, 6, 5, 8}};
  int gf2e[2][3] = {{0, 1, 2}, {3, 4, 5}};
  
  Face* gface;
  QuadFace* qface;
  //defines each face and put it into face_list
  Prism* aCell = new Prism(3);
  
  //gnrlface[0]
  for(int i = 0; i < 2; i++){
    gface = new Face(3);
    gface_list.push_back(gface);

    //define each edge
    for(int j = 0; j < 3; j++){
      gface->edge[j] = e2e[edge_entity[gf2e[i][j]]];
      gface->needReverse[j] = false; 
    }
    //resplit each face
    gface->resplit(facePlan[face_entity[i]], orientCode[i], bnode_list, edge_list);//resplit without orientCode
  
    //index the node
    int   nindex = 0;
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->tag = nodeTag[face_entity[i]][nindex++];
      //  checkNode(*np);
    } 
      
     
      
     
    bnode_begin= --(bnode_list.end());
    aCell->setFace(i, gface);
  }
  
  //quadface
  for(int i = 0; i < 3; i++){
    qface = new QuadFace(4);
    qface_list.push_back(qface);
      
    //define each edge
    for(int j = 0; j < 4; j++){
      qface->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    qface->resplit(facePlan[face_entity[i+2]],orientCode[i+2], bnode_list, edge_list);

    //index the node
    int   nindex = 0;
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->tag = nodeTag[face_entity[i+2]][nindex++];
      //  checkNode(*np);
    } 
      
    bnode_begin= --(bnode_list.end());
    aCell->setFace(i, qface);
  }
  
  return aCell;
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

//build a cell with edgePlan and facePlan, tag the nodes
//then resplit the edges and faces with edgePlan1 and facePlan1
Prism* build_resplit_prism_cell(const Entity* lower, int lower_size,
                                const Entity* upper, int upper_size,
                                const Entity* boundary_map, int boundary_map_size,
                                const Array<char,5>& prism2face,
                                const Array<char,6>& prism2node,
                                const Array<char,5>& orientCode,
                                const const_multiMap& face2node,
                                const const_multiMap& face2edge,
                                const const_MapVec<2>& edge2node,
                                const const_store<vect3d>& pos,
                                const const_store<std::vector<char> >& edgePlan,//the plan from last cycle
                                const const_store<std::vector<char> >& facePlan,//the plan from last cycle
                                const const_store<std::vector<char> >& edgePlan1,//the plan from this cycle
                                const const_store<std::vector<char> >& facePlan1, //the plan from this cycle  
                                const const_store<char>& posTag,
                                const const_store<std::vector<char> >& nodeTag,//tag according to the plan from last cycle
                                std::list<Node*>& bnode_list,
                                std::list<Node*>& node_list,
                                std::list<Edge*>& edge_list,
                                std::list<QuadFace*>& qface_list,
                                std::list<Face*>& gface_list,
                                const const_store<int>& node_remap,
                                const std::vector<char>& cellPlan,
                                const  std::vector<char>& cellNodeTag
                                ){
  
  Array<Entity, 5> face_entity = collect_prism_faces(lower, lower_size,
                                                     upper,upper_size,
                                                     boundary_map,boundary_map_size,
                                                     prism2face, node_remap);
  
  Array<Entity, 6> node_entity = collect_prism_vertices(face2node,
                                                        face_entity,
                                                        prism2node);



  Array<bool, 9> edge_reverse;
  Array<Entity, 9> edge_entity = collect_prism_edges( face_entity,
                                                      node_entity,
                                                      face2edge,
                                                      edge2node,
                                                      edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 6; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    aNode->tag = posTag[node_entity[i]];
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  //edge is built according to the direction inside the prism
  //and resplit with needReverse
  std::map<Entity,Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
  
  for(int i = 0; i < 9; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;
    
    
    //resplit the edge
    anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->tag =  nodeTag[edge_entity[i]][nindex++];
      // checkNode(*np);
    }
    bnode_begin = --(bnode_list.end());
  }
  //resplit the edges again without tagging the node 
  for(int i = 0; i < 9; i++){
    e2e[edge_entity[i]]->resplit(edgePlan1[edge_entity[i]],edge_reverse[i], bnode_list);
  }
  bnode_begin = --(bnode_list.end());//set the tag begin point

  
  int f2e[3][4]= {{0, 7, 3, 6}, {1, 8, 4, 7}, {2, 6, 5, 8}};
  int gf2e[2][3] = {{0, 1, 2}, {3, 4, 5}};
  
  Face* gface;
  QuadFace* qface;
  //defines each face and put it into face_list
  Prism* aCell = new Prism(3);
  
  //gnrlface[0]
  for(int i = 0; i < 2; i++){
    gface = new Face(3);
    gface_list.push_back(gface);

    //define each edge
    for(int j = 0; j < 3; j++){
      gface->edge[j] = e2e[edge_entity[gf2e[i][j]]];
      gface->needReverse[j] = false; 
    }
    //resplit each face
    gface->resplit(facePlan[face_entity[i]], orientCode[i], bnode_list, edge_list);//resplit without orientCode
  
    //index the node
    int   nindex = 0;
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->tag = nodeTag[face_entity[i]][nindex++];
      //  checkNode(*np);
    } 
      
     
      
     
    bnode_begin= --(bnode_list.end());
    aCell->setFace(i, gface);
  }
  
  //quadface
  for(int i = 0; i < 3; i++){
    qface = new QuadFace(4);
    qface_list.push_back(qface);
      
    //define each edge
    for(int j = 0; j < 4; j++){
      qface->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    qface->resplit(facePlan1[face_entity[i+2]],orientCode[i+2], bnode_list, edge_list);

  
    
    tag_quad_face(  face2node[face_entity[i+2]].begin(), 
                    face2edge[face_entity[i+2]].begin(),
                    edge2node,
                    edgePlan,
                    facePlan[face_entity[i+2]], orientCode[i+2],
                    nodeTag[face_entity[i+2]],//the tag for facePlan 
                    facePlan1[face_entity[i+2]],
                    bnode_list,//node list from facePlan1
                    bnode_begin);


    
    bnode_begin= --(bnode_list.end());
    aCell->setFace(i, qface);
  }
  //finish build
  
  //resplit cell
  std::vector<Prism*> cells;
  aCell->resplit( cellPlan, 
                  node_list,
                  edge_list,
                  qface_list,
                  gface_list,
                  cells);

  //tag nodes
  int nindex = 0;
  for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++){
    (*np)->tag = cellNodeTag[nindex++];
  }
  cells.clear();
  




  
  
  //resplit the faces again without tagging the node
  //  std::list<QuadFace*>::const_iterator qiter = qface_list.begin();
  std::list<Face*>::const_iterator giter = gface_list.begin();
  for(int i = 0; i < 2; i++){
    (*giter)->resplit(facePlan1[face_entity[i]], orientCode[i], bnode_list, edge_list);//resplit without orientCode
  }
 
  
  return aCell;
}

//build a cell with edgePlan, facePlan and cellPlan, tag the cells
//then resplit the edges and faces with edgePlan1 and facePlan1
Prism* build_resplit_prism_cell_ctag(const Entity* lower, int lower_size,
                                     const Entity* upper, int upper_size,
                                     const Entity* boundary_map, int boundary_map_size,
                                     const Array<char,5>& prism2face,
                                     const Array<char,6>& prism2node,
                                     const Array<char,5>& orientCode,
                                     const const_multiMap& face2node,
                                     const const_multiMap& face2edge,
                                     const const_MapVec<2>& edge2node,
                                     const const_store<vect3d>& pos,
                                     const const_store<std::vector<char> >& edgePlan,//the plan from last cycle
                                     const const_store<std::vector<char> >& facePlan,//the plan from last cycle
                                     const const_store<std::vector<char> >& edgePlan1,//the plan from this cycle
                                     const const_store<std::vector<char> >& facePlan1, //the plan from this cycle  
                                     std::list<Node*>& bnode_list,
                                     std::list<Node*>& node_list,
                                     std::list<Edge*>& edge_list,
                                     std::list<QuadFace*>& qface_list,
                                     std::list<Face*>& gface_list,
                                     const const_store<int>& node_remap,
                                     const std::vector<char>& cellPlan,
                                     const  std::vector<char>& fineCellTag
                                     ){
  
  Array<Entity, 5> face_entity = collect_prism_faces(lower, lower_size,
                                                     upper,upper_size,
                                                     boundary_map,boundary_map_size,
                                                     prism2face, node_remap);
  
  Array<Entity, 6> node_entity = collect_prism_vertices(face2node,
                                                        face_entity,
                                                        prism2node);



  Array<bool, 9> edge_reverse;
  Array<Entity, 9> edge_entity = collect_prism_edges( face_entity,
                                                      node_entity,
                                                      face2edge,
                                                      edge2node,
                                                      edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 6; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  //edge is built according to the direction inside the prism
  //and resplit with needReverse
  std::map<Entity,Edge*> e2e;
  
  for(int i = 0; i < 9; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;
    
    
    //resplit the edge
    anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);

  }
  //resplit the edges again without tagging the node 
  for(int i = 0; i < 9; i++){
    e2e[edge_entity[i]]->resplit(edgePlan1[edge_entity[i]],edge_reverse[i], bnode_list);
  }
  
  int f2e[3][4]= {{0, 7, 3, 6}, {1, 8, 4, 7}, {2, 6, 5, 8}};
  int gf2e[2][3] = {{0, 1, 2}, {3, 4, 5}};
  
  Face* gface;
  QuadFace* qface;
  //defines each face and put it into face_list
  Prism* aCell = new Prism(3);
  
  //gnrlface[0]
  for(int i = 0; i < 2; i++){
    gface = new Face(3);
    gface_list.push_back(gface);

    //define each edge
    for(int j = 0; j < 3; j++){
      gface->edge[j] = e2e[edge_entity[gf2e[i][j]]];
      gface->needReverse[j] = false; 
    }
    //resplit each face
    gface->resplit(facePlan[face_entity[i]], orientCode[i], bnode_list, edge_list);//resplit without orientCode
    aCell->setFace(i, gface);
  }
  
  //quadface
  for(int i = 0; i < 3; i++){
    qface = new QuadFace(4);
    qface_list.push_back(qface);
      
    //define each edge
    for(int j = 0; j < 4; j++){
      qface->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    qface->resplit(facePlan1[face_entity[i+2]],orientCode[i+2], bnode_list, edge_list);
    aCell->setFace(i, qface);
  }
  //finish build
  
  //resplit cell
  std::vector<Prism*> cells;
  aCell->resplit( cellPlan, 
                  node_list,
                  edge_list,
                  qface_list,
                  gface_list,
                  cells);

  //tag nodes
  int nindex = 0;
  for(std::vector<Prism*>::const_iterator np = cells.begin(); np!= cells.end(); np++){
    (*np)->setTag(fineCellTag[nindex++]);
  }
  cells.clear();
    
  //resplit the faces again without tagging the node
  //  std::list<QuadFace*>::const_iterator qiter = qface_list.begin();
  std::list<Face*>::const_iterator giter = gface_list.begin();
  for(int i = 0; i < 2; i++){
    (*giter)->resplit(facePlan1[face_entity[i]], orientCode[i], bnode_list, edge_list);//resplit without orientCode
  }
   
  return aCell;
}
Prism* build_prism_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,5>& prism2face,
                        const Array<char,6>& prism2node,
                        const Array<char,5>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        const const_store<char>& posTag,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& qface_list,
                        std::list<Face*>& gface_list,
                        const const_store<int>& node_remap){
  
  Array<Entity, 5> face_entity = collect_prism_faces(lower,lower_size,
                                                     upper,upper_size,
                                                     boundary_map,boundary_map_size,
                                                     prism2face, node_remap);
  
  Array<Entity, 6> node_entity = collect_prism_vertices(face2node,
                                                        face_entity,
                                                        prism2node);



  Array<bool, 9> edge_reverse;
  Array<Entity, 9> edge_entity = collect_prism_edges( face_entity,
                                                      node_entity,
                                                      face2edge,
                                                      edge2node,
                                                      edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 6; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    aNode->tag = posTag[node_entity[i]];
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  //edge is built according to the direction inside the prism
  //and resplit with needReverse
  std::map<Entity,Edge*> e2e;
  for(int i = 0; i < 9; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;
    
    
  }
  
  int f2e[3][4]= {{0, 7, 3, 6}, {1, 8, 4, 7}, {2, 6, 5, 8}};
  int gf2e[2][3] = {{0, 1, 2}, {3, 4, 5}};
  
  Face* gface;
  QuadFace* qface;
  //defines each face and put it into face_list
  Prism* aCell = new Prism(3);
  
  //gnrlface[0]
  for(int i = 0; i < 2; i++){
    gface = new Face(3);
    gface_list.push_back(gface);

    //define each edge
    for(int j = 0; j < 3; j++){
      gface->edge[j] = e2e[edge_entity[gf2e[i][j]]];
      gface->needReverse[j] = false; 
    }
    
    aCell->setFace(i, gface);
  }
  
  //quadface
  for(int i = 0; i < 3; i++){
    qface = new QuadFace(4);
    qface_list.push_back(qface);
      
    //define each edge
    for(int j = 0; j < 4; j++){
      qface->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
  
    aCell->setFace(i, qface);
  }
  
  return aCell;
}

Prism* build_prism_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,5>& prism2face,
                        const Array<char,6>& prism2node,
                        const Array<char,5>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& qface_list,
                        std::list<Face*>& gface_list,
                        const const_store<int>& node_remap){
  
  Array<Entity, 5> face_entity = collect_prism_faces(lower,lower_size,
                                                     upper,upper_size,
                                                     boundary_map,boundary_map_size,
                                                     prism2face, node_remap);
  
  Array<Entity, 6> node_entity = collect_prism_vertices(face2node,
                                                        face_entity,
                                                        prism2node);



  Array<bool, 9> edge_reverse;
  Array<Entity, 9> edge_entity = collect_prism_edges( face_entity,
                                                      node_entity,
                                                      face2edge,
                                                      edge2node,
                                                      edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 6; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  //edge is built according to the direction inside the prism
  //and resplit with needReverse
  std::map<Entity,Edge*> e2e;
  for(int i = 0; i < 9; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;
    
    
  }
  
  int f2e[3][4]= {{0, 7, 3, 6}, {1, 8, 4, 7}, {2, 6, 5, 8}};
  int gf2e[2][3] = {{0, 1, 2}, {3, 4, 5}};
  
  Face* gface;
  QuadFace* qface;
  //defines each face and put it into face_list
  Prism* aCell = new Prism(3);
  
  //gnrlface[0]
  for(int i = 0; i < 2; i++){
    gface = new Face(3);
    gface_list.push_back(gface);

    //define each edge
    for(int j = 0; j < 3; j++){
      gface->edge[j] = e2e[edge_entity[gf2e[i][j]]];
      gface->needReverse[j] = false; 
    }
    
    aCell->setFace(i, gface);
  }
  
  //quadface
  for(int i = 0; i < 3; i++){
    qface = new QuadFace(4);
    qface_list.push_back(qface);
      
    //define each edge
    for(int j = 0; j < 4; j++){
      qface->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
  
    aCell->setFace(i, qface);
  }
  
  return aCell;
}


















//parallel version
Prism* build_prism_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,5>& prism2face,
                        const Array<char,6>& prism2node,
                        const Array<char,5>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        const const_store<std::vector<char> >& edgePlan,
                        const const_store<std::vector<char> >& facePlan,
                        const const_store<int>& node_offset,
                        const const_store<int>&  face_l2f,
                        const const_store<int>&  node_l2f,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& qface_list,
                        std::list<Face*>& gface_list){

  Array<Entity, 5> face_entity = collect_prism_faces(lower,lower_size,
                                                     upper,upper_size,
                                                     boundary_map,boundary_map_size,
                                                     prism2face, face_l2f);
  
  Array<Entity, 6> node_entity = collect_prism_vertices(face2node,
                                                        face_entity,
                                                        prism2node);



  Array<bool, 9> edge_reverse;
  Array<Entity, 9> edge_entity = collect_prism_edges( face_entity,
                                                      node_entity,
                                                      face2edge,
                                                      edge2node,
                                                      edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 6; i++){
    Node* aNode = new Node(pos[node_entity[i]], node_l2f[node_entity[i]]);
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity,Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
  for(int i = 0; i < 9; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;


    //resplit the edge
    anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);

    //index the node
    int nindex = node_offset[edge_entity[i]];
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->index =  nindex++;
    }
    bnode_begin = --(bnode_list.end());
  }
  
  int f2e[3][4]= {{0, 7, 3, 6}, {1, 8, 4, 7}, {2, 6, 5, 8}};
  int gf2e[2][3] = {{0, 1, 2}, {3, 4, 5}};
  Face* gface;
  QuadFace* qface;
  //defines each face and put it into face_list
  Prism* aCell = new Prism(3);
  
 
  for(int i = 0; i < 2; i++){


    gface = new Face(3);
    gface_list.push_back(gface);
    
    //define each edge
    for(int j = 0; j < 3; j++){
     
      gface->edge[j] =  e2e[edge_entity[gf2e[i][j]]];
      gface->needReverse[j] = false; 
    }
    //resplit each face
   
    gface->resplit(facePlan[face_entity[i]],orientCode[i], bnode_list, edge_list);//resplit without orientCode

    //index the node
    int nindex = node_offset[face_entity[i]];
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->index =  nindex++;
    }
    bnode_begin = --(bnode_list.end());
    aCell->setFace(i, gface);
  }
  
  //quadface
  for(int i = 0; i < 3; i++){
    qface = new QuadFace(4);
    qface_list.push_back(qface);

    //define each edge
    for(int j = 0; j < 4; j++){
      qface->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    qface->resplit(facePlan[face_entity[i+2]],orientCode[i+2], bnode_list, edge_list);
    //index the node
    int nindex = node_offset[face_entity[i+2]];
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->index =  nindex++;
    }
    bnode_begin = --(bnode_list.end());
    aCell->setFace(i, qface);
  }
  
  
 
  return aCell;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
