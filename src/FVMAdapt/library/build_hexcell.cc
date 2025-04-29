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
// This file containing the functions that build a HexCell according to Loci data structures
//and refinementplans, when build a hexcell, each edge is defined as local edge in cell and resplit with
//needReverse. Each face is defined as local face in hexcell and resplit with orientcode


#include <vector>
#include <list>
#include <Loci.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <Tools/tools.h>
#include "hexcell.h"


using std::cerr;
using std::endl;
using std::swap;
using std::cout;
using std::vector;

void reorder_faces(const const_store<int>& node_remap, std::vector<Entity>& lower,
                   std::vector<Entity>& upper,
                   std::vector<Entity>& boundary_map);
  
//this function define the 6 faces of hexcell from loci data structures
Array<Entity, 6> collect_hex_faces( const Entity* lower, int lower_size,
                                    const Entity* upper, int upper_size,
                                    const Entity* boundary_map, int boundary_map_size,
                                    const Array<char, 6>& hex2face,
                                    const const_store<int>& node_remap){
  
  // Collect the entity designation of all faces of the cell cc in the all_faces array

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

  Array<Entity, 6> faces;
  for (int f=0; f<6; f++) {
    switch (hex2face[f]/6) {
    case 0:
      faces[f] = vlower[hex2face[f]%6];
      break;
    case 1:
      faces[f] = vupper[hex2face[f]%6];
      break;
    case 2:
      faces[f] = vboundary_map[hex2face[f]%6];
      break;
    default:
      cerr << " WARNING: illegal hex2face value" << endl;
      break;
    } 
  }
  return faces;
}


//this function define the 8 vertices of hexcell from loci data structures
Array<Entity, 8> collect_hex_vertices(const const_multiMap& face2node, const Array<Entity, 6>& faces,
                                      const Array<char, 8>& hex2node){
  Array<Entity, 8>  hex_vertices;
  for(int i = 0; i < 4; i++){
    hex_vertices[i] = face2node[faces[DOWN]][hex2node[i]];
  }
  for(int i = 4; i < 8; i++){
    hex_vertices[i] = face2node[faces[UP]][hex2node[i]];
  }
  
  std::swap(hex_vertices[1], hex_vertices[4]);
  std::swap(hex_vertices[3], hex_vertices[6]);
  return hex_vertices;
}




//collect the entity designation of all edges of the cell , and the direction of edges are stored in needReverse
Array<Entity, 12> collect_hex_edges(const Array<Entity, 6>& faces, const Array<Entity,8>& hex_vertices,
                                    const const_multiMap& face2edge, const const_MapVec<2>& edge2node,
                                    Array<bool, 12>& needReverse){
  entitySet all_edges;
  Array<Entity, 12> hex_edges;
  for(int fID = 0; fID < 6; fID++){
    for(int eID = 0; eID < 4; eID++){
      all_edges += face2edge[faces[fID]][eID];
    }
  }      
  if(all_edges.size() != 12){
    cerr << "WARNING: the number of edges is not 12" << endl;
    Loci:: Abort();
  }
  Entity head, tail;
  int node0[12] = {0, 1, 2, 3, 0, 1, 4, 5, 0, 2, 4, 6};
  int node1[12] = {4, 5, 6, 7, 2, 3, 6, 7, 1, 3, 5, 7};
  
  for( int e = 0; e < 12; e++ ){
    head = hex_vertices[node0[e]];
    tail =  hex_vertices[node1[e]];
    entitySet::const_iterator ee;
    for(ee = all_edges.begin(); ee != all_edges.end(); ee++){
      if(edge2node[*ee][0] == head && edge2node[*ee][1] == tail){
        hex_edges[e] = *ee;
        needReverse[e] = false;
        break;
      }
      else if(edge2node[*ee][1] == head && edge2node[*ee][0] == tail){
        hex_edges[e] = *ee;
        needReverse[e] = true;
        break;
      }
      
    }
    if(ee ==  all_edges.end()){
      cerr << "WARNING: can not find edge" << endl;
      Loci::Abort();
    }
  }
  
  return hex_edges;
}






  
HexCell* build_hex_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,6>& hex2face,
                        const Array<char,8>& hex2node,
                        const Array<char,6>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        const const_store<std::vector<char> >& edgePlan,
                        const const_store<std::vector<char> >& facePlan,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& face_list,
                        const const_store<int>& node_remap ){

  Array<Entity, 6> face_entity = collect_hex_faces(lower,lower_size,
                                                   upper,upper_size,
                                                   boundary_map,boundary_map_size,
                                                   hex2face,
                                                   node_remap);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                      face_entity,
                                                      hex2node);



  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                     node_entity,
                                                     face2edge,
                                                     edge2node,
                                                     edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
  for(int i = 0; i < 12; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;


    //resplit the edge
    anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);
  }
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    face[i]->resplit(facePlan[face_entity[i]],orientCode[i], bnode_list, edge_list);
  }
  
  HexCell* aCell = new HexCell(face);
  return aCell;
}



HexCell* build_hex_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,6>& hex2face,
                        const Array<char,8>& hex2node,
                        const Array<char,6>& orientCode,
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
                        std::list<QuadFace*>& face_list,
                        const const_store<int>& node_remap){

  Array<Entity, 6> face_entity = collect_hex_faces(lower,lower_size,
                                                   upper,upper_size,
                                                   boundary_map,boundary_map_size,
                                                   hex2face, node_remap);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                      face_entity,
                                                      hex2node);



  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                     node_entity,
                                                     face2edge,
                                                     edge2node,
                                                     edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    aNode->tag = posTag[node_entity[i]];
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
  
  for(int i = 0; i < 12; i++){
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
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    face[i]->resplit(facePlan[face_entity[i]],orientCode[i], bnode_list, edge_list);

    int   nindex = 0;
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->tag = nodeTag[face_entity[i]][nindex++];
      //  checkNode(*np);
    } 
    
    bnode_begin= --(bnode_list.end());
    
  }
  
  HexCell* aCell = new HexCell(face);
  return aCell;
}
//build a cell with edgePlan and facePlan, tag the nodes
//then resplit the edges and faces with edgePlan1 and facePlan1
HexCell* build_resplit_hex_cell(const Entity* lower, int lower_size,
                                const Entity* upper, int upper_size,
                                const Entity* boundary_map, int boundary_map_size,
                                const Array<char,6>& hex2face,
                                const Array<char,8>& hex2node,
                                const Array<char,6>& orientCode,
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
                                std::list<QuadFace*>& face_list,
                                const const_store<int>& node_remap,
                                const std::vector<char>& cellPlan,
                                const  std::vector<char>& cellNodeTag ){
  
  Array<Entity, 6> face_entity = collect_hex_faces(lower,lower_size,
                                                   upper,upper_size,
                                                   boundary_map,boundary_map_size,
                                                   hex2face, node_remap);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                      face_entity,
                                                      hex2node);



  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                     node_entity,
                                                     face2edge,
                                                     edge2node,
                                                     edge_reverse);

 
  //define each node and put it into bnode_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    aNode->tag = posTag[node_entity[i]];
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
  
  for(int i = 0; i < 12; i++){
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
  for(int i = 0; i < 12; i++){
    e2e[edge_entity[i]]->resplit(edgePlan1[edge_entity[i]],edge_reverse[i], bnode_list);
  }
  bnode_begin = --(bnode_list.end()); //set the start point for tag
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    face[i]->resplit(facePlan1[face_entity[i]],orientCode[i], bnode_list, edge_list);
    tag_quad_face(  face2node[face_entity[i]].begin(), 
                    face2edge[face_entity[i]].begin(),
                    edge2node,
                    edgePlan,
                    facePlan[face_entity[i]], orientCode[i],
                    nodeTag[face_entity[i]],//the tag for facePlan 
                    facePlan1[face_entity[i]],
                    bnode_list,//node list from facePlan1
                    bnode_begin);
  
    
    bnode_begin= --(bnode_list.end());
  }
  HexCell* aCell = new HexCell(face);
  
  //finish build
  
  //resplit cell
  std::vector<HexCell*> cells;
  aCell->resplit( cellPlan, 
                  node_list,
                  edge_list,
                  face_list,
                  cells);

  //tag nodes
  int nindex = 0;
  for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++){
    (*np)->tag = cellNodeTag[nindex++];
  }
  cells.clear();

  return aCell;
}
//build a cell with edgePlan, facePlan and cellPlan, tag the cells
//then resplit the edges and faces with edgePlan1 and facePlan1
HexCell* build_resplit_hex_cell_ctag(const Entity* lower, int lower_size,
                                     const Entity* upper, int upper_size,
                                     const Entity* boundary_map, int boundary_map_size,
                                     const Array<char,6>& hex2face,
                                     const Array<char,8>& hex2node,
                                     const Array<char,6>& orientCode,
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
                                     std::list<QuadFace*>& face_list,
                                     const const_store<int>& node_remap,
                                     const std::vector<char>& cellPlan,
                                     const  std::vector<char>& fineCellTag ){
  
  Array<Entity, 6> face_entity = collect_hex_faces(lower,lower_size,
                                                   upper,upper_size,
                                                   boundary_map,boundary_map_size,
                                                   hex2face, node_remap);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                      face_entity,
                                                      hex2node);



  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                     node_entity,
                                                     face2edge,
                                                     edge2node,
                                                     edge_reverse);

 
  //define each node and put it into bnode_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
    
  for(int i = 0; i < 12; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;

    //resplit the edge
    anEdge->resplit(edgePlan[edge_entity[i]],edge_reverse[i], bnode_list);
  }
  //resplit the edges again without tagging the node
  for(int i = 0; i < 12; i++){
    e2e[edge_entity[i]]->resplit(edgePlan1[edge_entity[i]],edge_reverse[i], bnode_list);
  }
 
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    face[i]->resplit(facePlan1[face_entity[i]],orientCode[i], bnode_list, edge_list);
    
  }
  HexCell* aCell = new HexCell(face);
  
  //finish build
  
  //resplit cell
  std::vector<HexCell*> cells;
  aCell->resplit( cellPlan, 
                  node_list,
                  edge_list,
                  face_list,
                  cells);

  //tag nodes
  int nindex = 0;
  for(std::vector<HexCell*>::const_iterator np = cells.begin(); np!= cells.end(); np++){
    (*np)->setTag(fineCellTag[nindex++]);
  }
  cells.clear();
  
  return aCell;
}
//parallel version
HexCell* build_hex_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,6>& hex2face,
                        const Array<char,8>& hex2node,
                        const Array<char,6>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        const const_store<std::vector<char> >& edgePlan,
                        const const_store<std::vector<char> >& facePlan,
                        const const_store<int>& node_offset,
                        const const_store<int>& face_l2f,
                        const const_store<int>& node_l2f,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& face_list){

  Array<Entity, 6> face_entity = collect_hex_faces(lower,lower_size,
                                                   upper,upper_size,
                                                   boundary_map,boundary_map_size,
                                                   hex2face, face_l2f);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                      face_entity,
                                                      hex2node);


  
  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                     node_entity,
                                                     face2edge,
                                                     edge2node,
                                                     edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]], node_l2f[node_entity[i]]);
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
  std::list<Node*>::const_iterator bnode_begin = --(bnode_list.end());
  for(int i = 0; i < 12; i++){
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
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
    //resplit each face
    face[i]->resplit(facePlan[face_entity[i]],orientCode[i], bnode_list, edge_list);

    //index the node
    int nindex = node_offset[face_entity[i]];
    for(std::list<Node*>::const_iterator np = ++bnode_begin; np!= bnode_list.end(); np++){
      (*np)->index =  nindex++;
    }
    bnode_begin = --(bnode_list.end());
  }
  
  HexCell* aCell = new HexCell(face);
  return aCell;
}


HexCell* build_hex_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,6>& hex2face,
                        const Array<char,8>& hex2node,
                        const Array<char,6>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        const const_store<char>& posTag,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& face_list,
                        const const_store<int>& node_remap){

  Array<Entity, 6> face_entity = collect_hex_faces(lower,lower_size,
                                                   upper,upper_size,
                                                   boundary_map,boundary_map_size,
                                                   hex2face, node_remap);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                      face_entity,
                                                      hex2node);



  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                     node_entity,
                                                     face2edge,
                                                     edge2node,
                                                     edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    aNode->tag = posTag[node_entity[i]];
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
  for(int i = 0; i < 12; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;
    
  }
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
  }
  
  HexCell* aCell = new HexCell(face);
  return aCell;
}

//build without restart(no edgeplan or faceplan)
HexCell* build_hex_cell(const Entity* lower, int lower_size,
                        const Entity* upper, int upper_size,
                        const Entity* boundary_map, int boundary_map_size,
                        const Array<char,6>& hex2face,
                        const Array<char,8>& hex2node,
                        const Array<char,6>& orientCode,
                        const const_multiMap& face2node,
                        const const_multiMap& face2edge,
                        const const_MapVec<2>& edge2node,
                        const const_store<vect3d>& pos,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& face_list,
                        const const_store<int>& node_remap){

  Array<Entity, 6> face_entity = collect_hex_faces(lower, lower_size,
                                                   upper,upper_size,
                                                   boundary_map,boundary_map_size,
                                                   hex2face, node_remap);
  
  Array<Entity, 8> node_entity = collect_hex_vertices(face2node,
                                                      face_entity,
                                                      hex2node);



  Array<bool, 12> edge_reverse;
  Array<Entity, 12> edge_entity = collect_hex_edges( face_entity,
                                                     node_entity,
                                                     face2edge,
                                                     edge2node,
                                                     edge_reverse);

 
  //define each node and put it into node_list
  std::map<Entity, Node*> n2n;
  for(int i = 0; i < 8; i++){
    Node* aNode = new Node(pos[node_entity[i]]);
    bnode_list.push_back(aNode);
    n2n[node_entity[i]] = aNode;
  }
  
  
  std::map<Entity, Edge*> e2e;
  for(int i = 0; i < 12; i++){
    Edge* anEdge = new Edge(n2n[edge2node[edge_entity[i]][edge_reverse[i]?1:0]],
                            n2n[edge2node[edge_entity[i]][edge_reverse[i]?0:1]]);
    edge_list.push_back(anEdge);
    e2e[edge_entity[i]] = anEdge;
    
  }
  
  int f2e[6][4]= {{6, 11, 7, 10}, {4, 9, 5, 8}, {2, 11, 3, 9}, {0, 10, 1, 8},
                  {1, 7, 3, 5}, {0, 6, 2, 4}};
  
  //defines each face and put it into face_list
  QuadFace** face = new QuadFace*[6];
  for(int i  = 0; i < 6; i++){
    face[i] = new QuadFace(4);
    face_list.push_back(face[i]);
    //define each edge
    for(int j = 0; j < 4; j++){
      face[i]->edge[j] = e2e[edge_entity[f2e[i][j]]];
    }
  }
  
  HexCell* aCell = new HexCell(face);
  return aCell;
}


