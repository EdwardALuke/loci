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
//////////////////////////////////////////////////////////////////////////////
//                          hexcell.h
//  This file includes the declaration of class HexCell, it's designed for
//  anisotropic refinement of hexahedra.
//  In a HexCell, all edges and all faces  point to the positive 
//  x, y, or z direction  
//   
//////////////////////////////////////////////////////////////////////////////

#ifndef HEXCELL_H
#define HEXCELL_H
#include <Loci.h>
#include <vector>
#include <bitset>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <utility> 
#include "hex_defines.h"
#include "quadface.h"
#include "read_par.h"
using std::cerr;
using std::endl;
using std::list;
std::vector<int32> get_c1_hex(const std::vector<char>& cellPlan,
                              const std::vector<char>& facePlan,
                              char orientCode,
                              char findex);


class HexCell
{
public:
  
  //constructor; 
  HexCell():cellIndex(0), mySplitCode(0), face(0), parentCell(0),
            childCell(0),tag(0){}
  HexCell(QuadFace** f):cellIndex(0), mySplitCode(0), face(f), parentCell(0),
                        childCell(0),tag(0){}  
  //destructor
  ~HexCell(){
    if(childCell != 0){
      for(int i = 0; i < numChildren(); i++){
	if(childCell[i] != 0){
	  delete  childCell[i];
	  childCell[i] = 0;
	}
      }
      delete[] childCell;
      childCell = 0;
    }
    parentCell = 0;
    
    if(face != 0){
      delete[] face;
      face = 0;
    }
  }

  //if all children are tagged as 2, remove all children
  void derefine();
  bool needDerefine();
  bool needDerefine_ctag();
  inline char getTag() const {return tag;}
  inline void setTag(char c){tag=c;}
  
  inline int32 getCellIndex() const {return cellIndex;}
  inline char getMySplitCode() const{return mySplitCode;}
  inline HexCell* getChildCell(int i)const{return childCell[i];}
  inline HexCell* getParentCell(){return parentCell;}

  //define if this is tagged for refinement, derefinement or unchanged
  int get_tagged();
  int get_tagged(const vector<source_par>& s);
  //return a splitCode
  //find   average_edge_length in XX , YY and ZZ directions
  //find min_edge_length in all directions
  
  //if max_edge_length/min_edge_length > Globals::factor1 and there are in 
  //if max_edge_length/min_edge_length > Globals::factor1 and they are in different direction
  //split the max_length edge
  void setSplitCode(int split_mode, double tol);

  
  inline int getLevel(NORMAL_DIRECTION d)const{
    switch(d){
    case XX: //x direction
      return face[FRONT]->edge[0]->level;
      break;
    case YY: //y direction
      return face[RIGHT]->edge[0]->level;
      break;
    case ZZ: // z direction
      return face[RIGHT]->edge[1]->level;
      break;
    default:
      cerr << "WARNING: illegal levelID" << endl;
      break;
    }
    return 0;
  }

  inline std::vector<Edge*> get_edges(){
    std::vector<Edge*> edges(12);
    edges[0] = face[3]->edge[0] ;
    edges[1] = face[3]->edge[2] ;
    edges[2] = face[2]->edge[0] ;
    edges[3] = face[2]->edge[2] ;
    edges[4] = face[1]->edge[0] ;
    edges[5] = face[1]->edge[2] ;
    edges[6] = face[0]->edge[0] ;
    edges[7] = face[0]->edge[2] ;
    edges[8] = face[1]->edge[3] ;
    edges[9] = face[1]->edge[1] ;
    edges[10] = face[0]->edge[3] ;
    edges[11] = face[0]->edge[1] ;
    return edges;
  }
 
  inline int numChildren()const{
    switch(mySplitCode){
    case 0:
      return 0;
    case 1:
    case 2:
    case 4:
      return 2;
    case 3:
    case 5:
    case 6:
      return 4;
    case 7:
      return 8;
    default:
      cerr<< "WARNING: illegal split code" << endl;
      break;
    }
    return -1;
  }
  //find numChildren without actually building the cell   
  inline int numChildren(char splitCode)const{
    switch(splitCode){
    case 0:
      return 0;
    case 1:
    case 2:
    case 4:
      return 2;
    case 3:
    case 5:
    case 6:
      return 4;
    case 7:
      return 8;
    default:
      cerr<< "WARNING: illegal split code" << endl;
      break;
    }
    return -1;
  }
  //find num_fine_cells without  actually building the tree
  int32 num_fine_cells( const std::vector<char>& cellPlan)const;
  
  int get_num_fine_faces()const;//for mxfpc
  double get_min_edge_length();   
 
  void split(std::list<Node*>& node_list,
             std::list<Edge*>& edge_list,
             std::list<QuadFace*>& face_list);
  
 
  //only define childCell
  void empty_split();
  //return num_fine_cells
  int empty_resplit(const std::vector<char>& cellPlan);
  
  //after the cell is split into a tree, get the indexMap from current index to parent index
  int32 traverse(const std::vector<char>& parentPlan,  vector<pair<int32, int32> >& indexMap);

  void resplit(const std::vector<char>& cellPlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::list<QuadFace*>& face_list,
               std::vector<HexCell*>& cells);

  //used in make_hex_cellplan.cc
  void resplit(int level,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::list<QuadFace*>& face_list);

  //this function check if aFace is one of the face
  //return -1: No
  //return i in [0, 6), Yes, face[i] == aFace
  inline int containFace(QuadFace* aFace){
    for(int i = 0; i <6; i++){
      if(face[i] == aFace) return i;
    }
    return -1;
  }
  
  //this function check if aCell is my sibling neighbor.
  //sibling means same size face, not necessarily has same parent
  bool isSiblingNeighbor(const HexCell* aCell, DIRECTION dd)const;


  //This function check if aCell is my neighbor is direction dd
  // edge neighbor is excluded. if two faces overlap area is not zero
  // their cells are neighbors
  bool  isNeighbor(const HexCell* aCell, DIRECTION dd)const;
  
  //this function find  my face neighbor,
  //ff: in,  my faceID,
  //nf: out, the faceID of neibCell 
  HexCell*  findNeighbor(DIRECTION d);

  //after a cell is split, compose the cell plan according to the tree structure 
  std::vector<char> make_cellplan();
  //make a plan for isotropical split an original hex cell level levels 
  std::vector<char> make_cellplan(int level);






  //if any edge is more than 1 levels down than my level, split myself, then balance each child 
  bool balance_cell(int split_mode,
                    std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<QuadFace*>& face_list);
  void sort_leaves(std::list<HexCell*>& v1);
  void rebalance_cells(int split_mode,
                       std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list,
                       std::list<QuadFace*>& face_list);

  //after a hexcell is resplit, find all it's inner faces and their c1, c2
  friend void set_hex_faces(const std::vector<HexCell*>& cells,
                            std::map<QuadFace*, NeibIndex>& faces);
  

                         
  //find c1 for all the leaves of a face
  friend std::vector<int32> get_c1_hex(const std::vector<char>& cellPlan,
                                       const std::vector<char>& facePlan,
                                       char orientCode,
                                       char findex);
private:
 
  int32 cellIndex; //the index of the cell, start with 1

  //describe how the cell will be splitted. the order of dimensions is xyz, and 0 means no split;
  // 1 means split; for example, 011(3) means  y and z coordinations are splitted.
  char mySplitCode;

  //6 faces,  pointing to positive xi, eta, or zeta direction
  // the numbering of face:
  // 0(RIGHT): xi = 1
  // 1(LEFT):  xi = 0
  // 2(FRONT): eta = 1
  // 3(BACK):  eta  = 0
  // 4(UP):    zeta =1
  // 5(DOWN):  zeta = 0
  QuadFace** face;
  

  HexCell *parentCell; // the parent of the cell
  HexCell **childCell;// an dynamic array of pointers to children cells
  

  
  std::bitset<6> faceMarked; //if the face in direction RIGHT, LEFT... has been checked
  char tag;
  //  char whichChild;
  //assignment and copying are prohibited
  void operator=(const HexCell&);
  HexCell(const HexCell&);
  
private:
  
  //get 8 nodes
  inline void get_nodes(std::vector<Node*>& node){
    for(int i = 0; i < 4; i++){
      node[i] = face[0]->getNode(i);
      node[i+4] = face[1]->getNode(i);
    }
  }

 
  
  //calculate the centroid of the HexCell, it's defined as
  //the mean value of nodes
  inline Node* simple_center(){
    Node* cellcenter = new Node();
    std::vector<Node*> vertices(8);
    get_nodes(vertices);
    std::vector<vect3d> nodes(8);
    for(int i = 0; i<8; i++){
      nodes[i] = vertices[i]->p;
    }
    cellcenter->p = point_center(nodes);
    return cellcenter;
  }

  inline Node* wireframe(){
    
    //allocate edgecenter
    std::vector<vect3d> facecenter(6);
    std::vector<double> areas(6);
    
    //get edge centers
    for(int i = 0; i < 6; i++){
      facecenter[i]= getFaceCenter(i)->p;
      areas[i] = face[i]->area();
    }
   
    //calculate the mass center of the edge centers
    vect3d p = weighted_center(facecenter, areas);
    return new Node(p);
  }
  //the center of the face, defined as the mass center of edge centers
  //precondition:: all its edges have been split
  inline Node* centroid(){
    switch(CENTROID){
    case 0:
      return simple_center();
    case 1:
      return wireframe();
    default:
      return wireframe();
    }
   
  }

  
  //get the facecenter
  //condition: facecenter will be allocated and deallocated by the caller
  //precondition: the faces have been splitted
  inline void getFaceCenter(Node** facecenter){
    for(int i = 0; i <6; i++){
      facecenter[i] = face[i]->getCenter();
    }
  }

  inline Node* getFaceCenter(int faceID){
    return face[faceID]->getCenter();
  }
};

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
                        const const_store<int>& node_remap);



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
                        const const_store<int>& node_remap);

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
                                const  std::vector<char>& cellNodeTag);
 

//build a cell with edgePlan and facePlan, tag the nodes
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
                                     const  std::vector<char>& fineCellTag);
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
                        const const_store<int>&  face_l2f,
                        const const_store<int>&  node_l2f, 
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& face_list);




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
                        const const_store<int>& node_remap);



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
                        const const_store<int>& node_remap);
//this function collect 6 faces of a cell 
Array<Entity, 6> collect_hex_faces( const Entity*  lower,int lower_size,
                                    const Entity* upper,int upper_size,
                                    const Entity* boundary_map,int boundary_map_size,
                                    const Array<char, 6>& hex2face, const const_store<int>& node_remap);

//this function collects 8 vertices of a cell
Array<Entity, 8> collect_hex_vertices(const const_multiMap& face2node, const Array<Entity, 6>& faces,
                                      const Array<char, 8>& hex2node);


//collect the entity designation of all edges of the cell in the all_edges entitySet
Array<Entity, 12> collect_hex_edges(const Array<Entity, 6>& faces, const Array<Entity, 8>& hex_vertices,
                                    const const_multiMap& face2edge, const const_MapVec<2>& edge2node,
                                    Array<bool,12>& needReverse);


//this function will define face2node for each fine faces and write them out,
//at the same time, mxppf(max num of points per face) will be updated
// void  write_hex_inner_faces(const std::list<pair<QuadFace*, NeibIndex> >& faces,
//                             int cell_offset, int& mxppf,std::ofstream& ofile);

std::vector<char>  extract_hex_face(const  std::vector<char>& cellPlan,  DIRECTION dd);
#endif
  
