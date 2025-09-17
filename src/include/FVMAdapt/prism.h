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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
//                               
//
//   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PRISM_H
#define PRISM_H
#include <Loci.h>
#include <vector>
#include <bitset>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <utility>
#include "read_par.h"
#include "hex_defines.h"
#include "quadface.h"
#include "face.h"


using std::cerr;
using std::endl;
using std::vector;
using std::stack;
using std::queue;
using std::cout;
using std::list;
std::vector<int32> get_c1_prism(const std::vector<char>& cellPlan,
                                const std::vector<char>& facePlan,
                                char orientCode,
                                int faceID);


class Prism
{
public:
  
  //constructor; 

  //this constructor is used when faces and nodes are not actually built
  Prism():cellIndex(0), nfold(3),mySplitCode(0),gnrlface(0),quadface(0), parentCell(0),
          childCell(0),tag(0){faceOrient.reset();}
  
  Prism(int n):cellIndex(0), nfold(n),mySplitCode(0),gnrlface(new Face*[2]),
               quadface(new QuadFace*[n]), parentCell(0), childCell(0),tag(0){
    faceOrient.reset();
  }
  
  //destructor
  ~Prism(){
    if(childCell != 0){
      int nc =  numChildren();
      for(int i = 0; i < nc; i++){
	if(childCell[i] != 0){
	  delete  childCell[i];
	  childCell[i] = 0;
	}
      }
      delete[] childCell;
      childCell = 0;
    }
    parentCell = 0;
    
    if(gnrlface != 0){
      delete[] gnrlface;
      gnrlface = 0;
    }
    
    if(quadface != 0){
      delete[] quadface;
      quadface = 0;
    }
  }
  
  //if all children are tagged as 2, remove all children
  bool needDerefine();
  bool needDerefine_ctag();
  void derefine();
  inline char getTag() const {return tag;}
  inline void setTag(char c){tag=c;}
  
  inline int32 getCellIndex() const {return cellIndex;}
  
  inline int getLevel(int d){
    switch(d){
    case 0: //xy direction
      return gnrlface[0]->edge[0]->level;
      break;
    case 1: //z direction
      return quadface[0]->edge[1]->level;
      break;
    default:
      cerr << "WARNING: illegal levelID" << endl;
      break;
    }
    return 0;
  }
  
  
  
  inline int getNfold() const{
    return nfold;
  }
  inline char getMySplitCode() const{
    return mySplitCode;
  }
  inline Prism* getChildCell(int i) const{
    return childCell[i];
  }
  inline Prism* getParentCell(){return parentCell;}
  
  inline int numChildren()const{
    switch(mySplitCode){
    case 0:
      return 0;
    case 1:
      return 2;
    case 2:
      return nfold;
    case 3:
      return 2*nfold;
    default:
      cerr<< "WARNING: illegal split code" << endl;
      break;
    }
    return -1;
  }
   
  //used in build_prismcell.cc
  inline void setFace(int faceID, Face* aFace){
    gnrlface[faceID] = aFace;
  }
  
  inline void setFace(int faceID, QuadFace* aFace){
    quadface[faceID] = aFace;
  }
  
  double get_min_edge_length();
 
  int get_num_fine_faces();//for mxfpc
  

  

  inline int whichChild(){
    if(parentCell == 0) return -1;
    for(int i = 0; i < parentCell->numChildren(); i++){
      if(this == parentCell->childCell[i]) return i;
    }
    return -1;
  }
  
    
  //this function splits diamondcell isotropially once
  // newly created nodes, edges and faces are put into the lists
  void split(std::list<Node*>& node_list,
             std::list<Edge*>& edge_list,
             std::list<QuadFace*>& quadface_list,
             std::list<Face*>& face_list);
  
  

  void resplit(const std::vector<char>& cellPlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::list<QuadFace*>& quadface_list,
               std::list<Face*>& face_list,
               std::vector<Prism*>& prism_cells);

  //used in make_prism_cellplan.cc
  void resplit( int level,
                std::list<Node*>& node_list,
                std::list<Edge*>& edge_list,
                std::list<QuadFace*>& quadface_list,
                std::list<Face*>& face_list);
 
  void empty_split();
  int empty_resplit(const std::vector<char>& cellPlan);
  //after the cell is split into a tree, get the indexMap from current index to parent index
  int32 traverse(const std::vector<char>& parentPlan,  vector<pair<int32, int32> >& indexMap);
  
 
  
  //this function check if aCell is my sibling neighbor.
  //sibling means same size face, not necessarily has same parent
  bool  isSiblingNeighbor(const Prism* aCell, int dd, int &nf)const;

  //get real sibling neib
  Prism* getSiblingNeib(int dd, int& nf);

  //only use it when dd >=2
  inline int parentFace(int dd){
    if(parentCell->mySplitCode == 1) return dd;
    
    else if(dd == 2){
      return (whichChild()%(parentCell->nfold))+2;
    }
    else if(dd == 5){
      int childID = (whichChild()%parentCell->nfold);
     
      return (childID== 0)?
        (parentCell->nfold +1):(childID+1);
    }
    return -1;
  }
  
  //This function check if aCell is my neighbor is direction dd
  // edge neighbor is excluded. if two faces overlap area is not zero
  // their cells are neighbors
  bool  isNeighbor(const Prism* aCell, int dd, int& nf)const;
  
  //this function find  my face neighbor,
  //ff: in,  my faceID,
  //nf: out, the faceID of neibCell 
  Prism*  findNeighbor(int d, int& nf);

  //define if this is tagged for refinement, derefinement or unchanged
  int get_tagged();
  int get_tagged(const vector<source_par>& s);
  void setSplitCode(int split_mode, double tol);
  //after a cell is split, compose the cell plan according to the tree structure 
  std::vector<char> make_cellplan();
  
  //make a plan for isotropical split an original prism cell level levels 
  std::vector<char> make_cellplan(int level);
  
  inline std::vector<Edge*> get_edges(){
  
    std::vector<Edge*> edges(3*nfold);
    for(int i = 0; i < nfold; i++) edges[i] = gnrlface[0]->edge[i];
    for(int i = nfold; i< 2*nfold; i++) edges[i] = gnrlface[1]->edge[i-nfold];
    for(int i = 2*nfold; i < 3*nfold; i++){
      int j = i%nfold;
      edges[i] = quadface[j]->edge[faceOrient.test(j)?1:3];
    }

    return edges;
  }
  //if any edge is more than 1 levels down than my level, split myself, then balance each child 
  bool balance_cell(int split_mode,
                    std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<QuadFace*>& qface_list,
                    std::list<Face*>& gface_list);
  void sort_leaves(std::list<Prism*>& v1);  
  void rebalance_cells(int split_mode,
                       std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list,
                       std::list<QuadFace*>& qface_list,
                       std::list<Face*>& gface_list);
  
  friend void set_prism_faces(const std::vector<Prism*>& cells,
                              std::map<QuadFace*, NeibIndex>& quadfaces,
                              std::map<Face*, NeibIndex>& faces);
  
  friend std::vector<int32> get_c1_prism(const std::vector<char>& cellPlan,
                                         const std::vector<char>& facePlan,
                                         char orientCode,
                                         int faceID);
  friend std::vector<char>  extract_prism_face(const  std::vector<char>& cellPlan, int dd);
  
  void print();
private:
 
  int32 cellIndex; //the index of the cell, start with 1
  char nfold; //3 for normal prism, 4 for the children of prism when quadface is split
  
  //describe how the cell will be splitted. the order of dimensions is xyz, and 0 means no split;
  // 1 means split; for example, 011 means  y and z coordinations are splitted.
  char mySplitCode;
  Face** gnrlface;
  QuadFace** quadface;
  //orientCode is not necessary here
  
  Prism *parentCell; // the parent of the cell
  Prism **childCell;// an dynamic array of pointers to children cells

  
  
  std::bitset<6> faceMarked; //if the face in direction RIGHT, LEFT... has been checked
  //when a Prism is first created, gnrlface[0] points inward, gnrlface[1] points outward
  // thress quadfaces point outward. When it's split, gnrlface orientation won't change
  // quadface orientation of the children need to be defined because when a new face is created,
  //it can  not make both of two children sharing it keep their face orientation
  // so when nfold==4, if quadface[i] points inward, faceOrient[i] = 1; else faceOrient[i] = 0;
  std::bitset<4> faceOrient;
  char tag;
  //  char whichChild;
  //assignment and copying are prohibited
  void operator=(const Prism&);
  Prism(const Prism&);
  
private:
  //get all the leaves
  void get_leaves(std::vector<Prism*>& leaf_cell);

  //get 6 nodes
  inline void get_nodes(std::vector<Node*>& node){
    node.resize(2*nfold);
    for(int i = 0; i < nfold; i++){
      node[i] = (gnrlface[0]->needReverse[i])?(gnrlface[0]->edge[i]->tail):gnrlface[0]->edge[i]->head;
      node[i+nfold] = (gnrlface[1]->needReverse[i])?(gnrlface[1]->edge[i]->tail):gnrlface[1]->edge[i]->head;
    }
  }
  
  //get all the 4*3 edges
  // inline void get_edges(Edge** edge){}
  
  
  //calculate the centroid of the Prism, it's defined as
  //the mean value of nodes
  
  inline Node* simple_center(){
   
    Node* cellcenter = new Node();
    std::vector<Node*> vertices(2*nfold);
    get_nodes(vertices);
    std::vector<vect3d> nodes(2*nfold);
    for(int i = 0; i<2*nfold; i++){
      nodes[i] = vertices[i]->p;
    }
    cellcenter->p = point_center(nodes);
    return cellcenter;
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

  inline Node* wireframe(){
    
    //allocate edgecenter
    std::vector<vect3d> facecenter(nfold+2);
    std::vector<double> areas(nfold+2);
    
    //get edge centers
    for(int i = 0; i < nfold+2; i++){
      facecenter[i]= getFaceCenter(i)->p;
      if(i<2)areas[i] = gnrlface[i]->area();
      else areas[i] = quadface[i-2]->area();
     
    }
   
    //calculate the mass center of the edge centers
    vect3d p = weighted_center(facecenter, areas);
    return new Node(p);
  }

  
  
  

  inline Node* getFaceCenter(int faceID){
    if(faceID < 2) return  gnrlface[faceID]->child[0]->edge[2]->head;
    return quadface[faceID-2]->getCenter();
  }

 
};

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
                        const const_store<int>& node_remap);
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
                        const const_store<int>& node_remap);
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
                                const const_store<std::vector<char> >& edgePlan,
                                const const_store<std::vector<char> >& facePlan,
                                const const_store<std::vector<char> >& edgePlan1,
                                const const_store<std::vector<char> >& facePlan1,   
                                const const_store<char>& posTag,
                                const const_store<std::vector<char> >& nodeTag,
                                std::list<Node*>& bnode_list,
                                std::list<Node*>& node_list,
                                std::list<Edge*>& edge_list,
                                std::list<QuadFace*>& qface_list,
                                std::list<Face*>& gface_list,
                                const const_store<int>& node_remap,
                                const std::vector<char>& cellPlan,
                                const  std::vector<char>& cellNodeTag);

//build a cell with edgePlan and facePlan, tag the nodes
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
                                     const const_store<std::vector<char> >& edgePlan,
                                     const const_store<std::vector<char> >& facePlan,
                                     const const_store<std::vector<char> >& edgePlan1,
                                     const const_store<std::vector<char> >& facePlan1,   
                                     std::list<Node*>& bnode_list,
                                     std::list<Node*>& node_list,
                                     std::list<Edge*>& edge_list,
                                     std::list<QuadFace*>& qface_list,
                                     std::list<Face*>& gface_list,
                                     const const_store<int>& node_remap,
                                     const std::vector<char>& cellPlan,
                                     const  std::vector<char>& fineCellTag);
// for no restart
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
                        const const_store<int>& node_remap);



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
                        const const_store<int>& node_remap);

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
                        const const_store<int>& face_l2f,
                        const const_store<int>& node_l2f,
                        std::list<Node*>& bnode_list,
                        std::list<Edge*>& edge_list,
                        std::list<QuadFace*>& qface_list,
                        std::list<Face*>& gface_list);



int general_childID_orient_c2f(int childID_c, char orientCode, int numEdge);
int general_childID_orient_f2c(int childID_f, char orientCode, int numEdge);
#endif
  
