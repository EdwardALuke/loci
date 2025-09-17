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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                             diamondcell.h                                                            //
//    This file includes the declaration of class DiamondCell .
//                               
//
//   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef DIAMONDCELL_H
#define DIAMONDCELL_H
#include <Loci.h>
#include <vector>
#include <bitset>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <utility> 
//#include <cassert>
#include "defines.h"
#include "face.h"
#include "read_par.h"

using std::cerr;
using std::endl;
using std::list;

class Cell;
//when a general cell is splitted isotropically, it will become m DiamondCell
//where m is number of node of the general cell
//a DiamondCell with fold n will have 2n+2 nodes, 2n face, 4n edges
//n faces share node 0, and n faces share node 1, 
class DiamondCell
{
public:
  
  //constructor; 
  DiamondCell(char m):nfold(m),cellIndex(0),parentCell(0),childCell(0),
                      face(new Face*[2*m]), faceOrient(new char[2*m]), faceMarked(0),tag(0){}
  
  
  //destructor
  ~DiamondCell(){
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
    
    if(faceMarked != 0){
      delete[] faceMarked;
      faceMarked = 0;
    }
    
    if(faceOrient!=0){
      delete [] faceOrient;
      faceOrient =0;
    }
    
    if(face != 0){
      delete [] face;
      face = 0;
    }
    
    parentCell = 0;
  }
  //if all children are tagged as 2, remove all children
  void derefine();
  bool needDerefine();
   bool needDerefine_ctag();
  inline char getTag() const {return tag;}
  inline void setTag(char c){tag=c;}
    
  inline void setCellIndex(int32 cellID){cellIndex = cellID;}
  
  inline int32 getCellIndex() const {return cellIndex;}
  
  inline char getNfold() const{return nfold;}
  
  inline int getLevel()const{return face[0]->edge[0]->level;} 
 
  //define if this is tagged for refinement, derefinement, or unchanged
  int get_tagged();
  int get_tagged(const vector<source_par>& s);
  int get_num_fine_faces();//for mxfpc
  inline void setParentCell( DiamondCell* parent){parentCell = parent;}
  inline DiamondCell* getParentCell(){return parentCell;}
    
  inline DiamondCell* getChildCell(int i)const {return childCell[i];}
  
  inline DiamondCell** getChildCell()const {return childCell;}
  
  //when faceID >= nfold, a cell will share the same face with its parentCell
  //this function find the faceID in parentCell
  int parentFace(int faceID)const;
  
  //this function splits DiamondCell isotropically level times
  // newly created nodes, edges and faces are put into the lists
  void resplit(int level,std::list<Node*>& node_list,
               std::list<Edge*>& edge_list, std::list<Face*>& face_list);

  //this function splits diamondcell isotropially once
  // newly created nodes, edges and faces are put into the lists
  void split(std::list<Node*>& node_list,
             std::list<Edge*>& edge_list,
             std::list<Face*>& face_list);

  //this function splits diamondcell isotropically once
  //only define childCell
  void empty_split();

  //this function check if aFace is one of the face
  //return -1: No
  //return i in [0, 2*nfold), Yes, face[i] == aFace
  inline int containFace(Face* aFace){
    for(int i = 0; i <2*nfold; i++){
      if(face[i] == aFace) return i;
    }
    return -1;
  }
  
  
  //this function find  my face neighbor,
  //ff: in,  my faceID,
  //nf: out, the faceID of neibCell 
  DiamondCell*  findNeighbor(const Cell* aCell,
                             const std::vector<std::vector<Edge*> >& n2e,
                             int ff, int& nf)const;
  

  //this function find  my sibling face neighbor,
  //ff: in,  my faceID,
  //nf: out, the faceID of neibCell 
  DiamondCell* getSiblingNeib(const Cell* aCell,
                              const std::vector<std::vector<Edge*> >& n2e,
                              int mf, int& nf)const;
  
  
  //if any edge is more than 1 levels down than my level, split myself, then balance each child 
  bool balance_cell(std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<Face*>& face_list);
  
  void sort_leaves(std::list<DiamondCell*>& v1);  
  
  
  //after aCell is splitted into cells, for each inner fine face, define its c1, c2
  friend  void set_general_faces(const Cell* aCell,
                                 const std::vector<DiamondCell*>& cells,
                                 const std::vector<std::vector<Edge*> >& n2e,
                                 std::list<pair<Face*, NeibIndex> >& fine_face);
  
 

  //find the minimum edge length in a cell(before split)
  inline double get_min_edge_length(){
    std::set<Edge*> edge;
    get_edges(edge);
    std::set<Edge*>::const_iterator cur_edge = edge.begin();
    
    double min_length = norm((*cur_edge)->head->p - (*cur_edge)->tail->p);
    for(cur_edge= edge.begin();cur_edge != edge.end(); cur_edge++){
      min_length = min(min_length, norm((*cur_edge)->head->p - (*cur_edge)->tail->p));
    }
    return min_length;
  }

private:
  //a nfold diamondcell will have 2*nfold faces and 2*nfold+2 vertex, node 0 and node 1 have
  //nfold edges, other vertices has 3 edges 
  char nfold;
  
  int32 cellIndex; //the index of the cell, start with 1
    
  DiamondCell *parentCell; // the parent of the cell

  DiamondCell **childCell;// an dynamic array of pointers to children cells

  Face** face;
  // points inward, faceOrient[i] = 1; else faceOrient[i] = 0;
  char* faceOrient; //the face points inward or outward 
  
  bool* faceMarked;//if the face has been checked. size: 2*nfold

  char whichChild; //used in findNeighbor() function
  char tag;//cell tag
  //assignment and copying are prohibited
  void operator=(const DiamondCell&);
  DiamondCell(const DiamondCell&);
  friend class Cell;
private:
  //  //get all the leaves
  //   void get_leaves(std::vector<DiamondCell*>& leaf_cell);

  //get all the 2*nfold+2 nodes
  void get_nodes(std::set<Node*>& node);

  //get all the 4*nfold edges
  void get_edges(std::set<Edge*>& edge);
  
  //calculate the centroid of the diamondcell, it's defined as
  //the mean value of facecenters
  //precondition: all the faces have been splitted
  inline Node* simple_center(){
    Node* cellcenter = new Node();
    std::vector<vect3d> facecenter(2*nfold);
    for(int i = 0; i < 2*nfold ; i++){
      facecenter[i] = face[i]->child[0]->edge[2]->head->p;
    }
    cellcenter->p = point_center(facecenter);
   
    return cellcenter;
  }
  inline Node* wireframe(){
    
    //allocate edgecenter
    std::vector<vect3d> facecenter(2*nfold);
    std::vector<double> areas(2*nfold);
    
    //get edge centers
    for(int i = 0; i < 2*nfold; i++){
      facecenter[i]= face[i]->child[0]->edge[2]->head->p;
      areas[i] = face[i]->area();
    }
   
    //calculate the mass center of the edge centers
    vect3d p = weighted_center(facecenter, areas);

    
    return new Node(p);
  }
  
  
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
    for(int i = 0; i < 2*nfold; i++){
      facecenter[i] = face[i]->child[0]->edge[2]->head;
    }
  }

  //get the edgecenter
  //precondition: all the faces and edges has been splitted
  //condition: edgecenter need be alloctaed and deallocated by the caller
  //edge center is defined from the edges of faces according to the number system and faceOrient
  inline void getEdgeCenter(Node** edgecenter){
    for(int i = 0; i <nfold; i++){
      Node* ecenter[4];//edgecenter of a face
      face[i]->getEdgeCenter(ecenter);
      if(faceOrient[i]== 1){//inward
        edgecenter[i] = ecenter[0];
        edgecenter[i+2*nfold] = ecenter[1];
      }
      else{//outward
        edgecenter[i] = ecenter[3];
        edgecenter[i+2*nfold] = ecenter[2];
      }
      
    }
    for(int i = nfold; i <2*nfold; i++){
      Node* ecenter[4];//edgecenter of a face
      face[i]->getEdgeCenter(ecenter);
      if(faceOrient[i] == 0){//outward
        edgecenter[i] = ecenter[3];
        edgecenter[i+2*nfold] = ecenter[1];
      }
      else{//inward
        edgecenter[i] = ecenter[0];
        edgecenter[i+2*nfold] = ecenter[2];
      }
    }
  }
};

//this class defines the general cell
class Cell{
public:
  // constructors
  Cell(int nd, int ne, int nf, Node** n, Edge** e, Face** f, char* fo):
    numNode(nd), numEdge(ne), numFace(nf), node(n), edge(e), face(f),
    faceOrient(fo),child(0){}

  Cell():node(0), edge(0), face(0), faceOrient(0),child(0){}

  //destructor
  ~Cell(){
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
    if(faceOrient !=0){
      delete [] faceOrient;
      faceOrient = 0;
    }
    if(face!=0){
      delete [] face;
      face = 0;
    }
    
    if(edge !=0){
      delete [] edge;
      edge = 0;
    }
    if(node !=0){
      delete [] node;
      node = 0;
    }
  }

  //if all children are tagged as 2, remove all children
  bool needDerefine();
  bool needDerefine_ctag();
  void derefine();
  

  //center of the cell, defined as the mean value of the facecenter
  //precondition: all the face and edge have been splitted
  inline Node* simple_center(){
    Node* center = new Node();
    std::vector<vect3d> facecenter(numFace);
    for(int i = 0; i < numFace; i++){
      facecenter[i] = face[i]->child[0]->edge[2]->head->p;
    }
    center->p = point_center(facecenter);
    return center;
  }

  inline Node* wireframe(){
    
    //allocate edgecenter
    std::vector<vect3d> facecenter(numFace);
    std::vector<double> areas(numFace);
    
    //get edge centers
    for(int i = 0; i < numFace; i++){
      facecenter[i]=face[i]->child[0]->edge[2]->head->p;
      areas[i] = face[i]->area();
    }
   
    //calculate the mass center of the edge centers
    vect3d p = weighted_center(facecenter, areas);
    return new Node(p);
  }

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

  
  //precondition: all the faces have been splitted
  //condition: facecenter is allocated and deallocated by the caller
  inline void getFaceCenter(Node** facecenter){
    for(int i = 0; i < numFace; i++){
      facecenter[i] = face[i]->child[0]->edge[2]->head;
    }
  }

  void sort_leaves(std::list<DiamondCell*>& v1);//depth first search
 
  //find the minimum edge length in a cell(before split)
  inline double get_min_edge_length(){
    double min_length = norm(edge[0]->head->p - edge[0]->tail->p);
    for(int i = 1; i < numEdge; i++){
      min_length = min(min_length, norm(edge[i]->head->p - edge[i]->tail->p));
    }
    return min_length;
  }

  //split a general cell isotropically into DiamondCells 
  void split( std::list<Node*>& node_list,
              std::list<Edge*>& edge_list,
              std::list<Face*>& face_list);
  void resplit( const std::vector<char>& cellPlan,
                std::list<Node*>& node_list,
                std::list<Edge*>& edge_list,
                std::list<Face*>& face_list,
                std::vector<DiamondCell*>& cells);
  //split a general cell isotropically into DiamondCells
  //only define child
  void empty_split();
  //return num_fine_cell
  int32 empty_resplit(const std::vector<char>& cellPlan);
  //  void get_leaves(std::vector<DiamondCell*>& leaf_cell);
  //after the cell is split into a tree, get the indexMap from current index to parent index
  int32 traverse(const std::vector<char>& parentPlan,  vector<pair<int32, int32> >& indexMap);
  int get_num_fine_faces();//for calculating mxfpc
  
  //find the node2edge of for all nodes, for each node, n2e[i] and n2e[i+1] share a face
  std::vector<std::vector<Edge*> > set_n2e();

  //find the info for node nindex. rot.size() == n2f.size()
  //if j=rot[i] >= 0,the orient of n2f[i] is 1, and node[nindex] is the jth node of n2f[i]
  //if j=rot[i] < 0, the orient of n2f[i] is -1, and  node[nindex] is (-j-1)th node of n2f[i] 
  void set_n2f_n2e(std::vector<Face*>& n2f, std::vector<Edge*>& n2e, std::vector<int>& rot, int nindex);
  
  //check if the cell is tagged for refinement, derefinement, or unchanged
  int get_tagged();
  //check if the cell is tagged for refinement
  int get_tagged(const vector<source_par>& s);
  //if any edge is more than 1 levels down than my level, split myself, then balance each child 
  bool balance_cell(std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<Face*>& face_list);

  void rebalance_cells(std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list,
                       std::list<Face*>& face_list);
  std::vector<char> make_cellplan();

public:
  int numNode;
  int numEdge;
  int numFace;
  Node** node; 
  Edge** edge;
  Face** face;
  
  char* faceOrient; //the face points inward or outward, points inward, faceOrient[i] = 1; else faceOrient[i] = 0;
  DiamondCell **child;// an dynamic array of pointers to children cells
};

int  find_face_index(const Entity* lower, int lower_size,
                     const Entity* upper, int upper_size,
                     const Entity* boundary_map, int boundary_map_size,
                     const const_multiMap& face2node,
                     Entity f,
                     const const_store<int>& node_remap);





//build a Cell from Loci data structures, the locations of nodes are defined
//edges and faces are split according to edgePlan and facePlan
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
                         const const_store<int>& node_remap);

//build a Cell from Loci data structures, the locations of nodes are defined
//edges and faces are split according to edgePlan and facePlan
//and all boundary nodes are tagged
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
                         const const_store<int>& node_remap);


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
                                 const  std::vector<char>& cellNodeTag);





//build a cell with edgePlan and facePlan, tag the cells
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
                                      const  std::vector<char>& fineCellTag);


//build a Cell from Loci data structures, the locations of nodes are defined
//and all boundary nodes are tagged
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
                         const const_store<int>& node_remap);




//build a Cell from Loci data structures, the locations of nodes are defined

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
                         const const_store<int>& node_remap);


//build a Cell from Loci data structures, the locations of nodes are not defined
Cell* build_general_cell(const Entity* lower, int lower_size,
                         const Entity* upper, int upper_size,
                         const Entity* boundary_map, int boundary_map_size,
                         const const_multiMap& face2node,
                         const const_multiMap& face2edge,
                         const const_MapVec<2>& edge2node,
                         std::list<Node*>& bnode_list,
                         std::list<Edge*>& edge_list,
                         std::list<Face*>& face_list,
                         const const_store<int>& node_remap);



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
                         const const_store<int>& node_remap);










//this function first build a Face f and a Cell cl[f]/cr[f], then split the Face according to facePlan
// and split the Cell according to cellPlan. for each fine face, find the index of fine cell that it belongs to
std::vector<int32> get_c1(const Entity* lower, int lower_size,
                          const Entity* upper, int upper_size,
                          const Entity* boundary_map, int boundary_map_size,
                          const const_multiMap& face2node, 
                          const const_multiMap& face2edge,
                          const const_MapVec<2>& edge2node,
                          const std::vector<char>& cellPlan,
                          const std::vector<char>& facePlan,
                          Entity f,
                          const const_store<int>& node_remap);

std::vector<int32> get_c1_general(const Entity* lower, int lower_size,
                                  const Entity* upper, int upper_size,
                                  const Entity* boundary_map, int boundary_map_size,
                                  bool is_quadface,
                                  const const_multiMap& face2node, 
                                  const const_multiMap& face2edge,
                                  const const_MapVec<2>& edge2node,
                                  const std::vector<char>& cellPlan,
                                  const std::vector<char>& facePlan,
                                  Entity f,
                                  const const_store<int>& node_remap);


//this function merge two isotropical facePlan
std::vector<char> merge_faceplan(std::vector<char>& planl, std::vector<char>& planr, int numNodes);

//this function extract facePlan from a cellPlan
std::vector<char>   extract_general_face(const Entity* lower, int lower_size,
                                         const Entity* upper, int upper_size,
                                         const Entity* boundary_map, int boundary_map_size,
                                         const const_multiMap& face2node,
                                         const const_multiMap& face2edge,
                                         const const_MapVec<2>& edge2node,
                                         const std::vector<char>& cellPlan,
                                         Entity ff, const const_store<int>& node_remap);


#endif




  
  
  


