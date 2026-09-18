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

using std::cerr ;
using std::endl ;
using std::vector ;
using std::stack ;
using std::queue ;
using std::cout ;
using std::list ;

/**
 * @file prism.h
 *
 * Prism refinement with Face end faces and QuadFace side faces.
 */
std::vector<int32> get_c1_prism(const std::vector<char>& cellPlan,
                                const std::vector<char>& facePlan,
                                char orientCode,
                                int faceID) ;

/**
 * Cell with two nfold-sided Face end faces and nfold QuadFace side faces. An
 * original triangular prism has nfold == 3. Splitting the end faces creates
 * children with nfold == 4.
 *
 * mySplitCode selects an axial split (1), an end-face split (2), or both (3).
 * numChildren() returns 2, nfold, or 2*nfold for those codes.
 */
class Prism{
public:

  /// This constructor is used when faces and nodes are not actually built
  Prism():cellIndex(0), nfold(3),mySplitCode(0),gnrlface(0),quadface(0), parentCell(0),
          childCell(0),tag(0) { faceOrient.reset() ; }

  Prism(int n):cellIndex(0), nfold(n),mySplitCode(0),gnrlface(new Face*[2]),
               quadface(new QuadFace*[n]), parentCell(0), childCell(0),tag(0) {
    faceOrient.reset() ;
  }

  ~Prism(){
    if(childCell != 0) {
      int nc = numChildren() ;
      for(int i=0; i<nc; i++) {
        if(childCell[i] != 0) {
          delete childCell[i] ;
          childCell[i] = 0 ;
        }
      }
      delete[] childCell ;
      childCell = 0 ;
    }
    parentCell = 0 ;

    if(gnrlface != 0) {
      delete[] gnrlface ;
      gnrlface = 0 ;
    }

    if(quadface != 0) {
      delete[] quadface ;
      quadface = 0 ;
    }
  }

  /// Return true if this cell has children, every child is a leaf requesting
  /// derefinement through its node tags, and no boundary edge has
  /// grandchildren. This only tests eligibility; it does not remove children.
  bool needDerefine() ;
  /// Test the same leaf and edge conditions as needDerefine(), using each
  /// child's cell tag (getTag() == 2) instead of its node tags.
  bool needDerefine_ctag() ;
  /// Delete the children and reset mySplitCode to 0. This does not check tags
  /// or derefinement eligibility.
  void derefine() ;
  char getTag() const { return tag ; }
  void setTag(char c) { tag=c ; }

  int32 getCellIndex() const { return cellIndex ; }

  int getLevel(int d) {
    switch(d) {
    case 0: // xy direction
      return gnrlface[0]->edge[0]->level ;
      break ;
    case 1: //z direction
      return quadface[0]->edge[1]->level ;
      break ;
    default:
      cerr << "WARNING: illegal levelID" << endl ;
      break ;
    }
    return 0 ;
  }

  int getNfold() const {
    return nfold ;
  }

  char getMySplitCode() const {
    return mySplitCode ;
  }

  Prism* getChildCell(int i) const {
    return childCell[i] ;
  }

  Prism* getParentCell() { return parentCell ; }

  int numChildren() const {
    switch(mySplitCode) {
    case 0:
      return 0 ;
    case 1:
      return 2 ;
    case 2:
      return nfold ;
    case 3:
      return 2*nfold ;
    default:
      cerr << "WARNING: illegal split code" << endl ;
      break ;
    }
    return -1 ;
  }

  /// Set gnrlface[faceID] to aFace.
  void setFace(int faceID, Face* aFace) {
    gnrlface[faceID] = aFace ;
  }

  void setFace(int faceID, QuadFace* aFace) {
    quadface[faceID] = aFace ;
  }

  double get_min_edge_length() ;

  /// Return the number of leaf faces on the cell boundary, used to compute the
  /// maximum faces per cell (mxfpc).
  int get_num_fine_faces() ;

  int whichChild() {
    if(parentCell == 0) { return -1 ; }
    for(int i=0; i<parentCell->numChildren(); i++) {
      if(this == parentCell->childCell[i]) { return i ; }
    }
    return -1 ;
  }

  /// Split this Prism according to mySplitCode. Append new nodes, root edges,
  /// and root QuadFace and Face objects to the supplied lists.
  void split(std::list<Node*>& node_list,
             std::list<Edge*>& edge_list,
             std::list<QuadFace*>& quadface_list,
             std::list<Face*>& face_list) ;

  void resplit(const std::vector<char>& cellPlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::list<QuadFace*>& quadface_list,
               std::list<Face*>& face_list,
               std::vector<Prism*>& prism_cells) ;

  /// Apply level isotropic cell splits in breadth-first order. The counter
  /// decreases for each split cell, rather than for each tree level. Append
  /// new nodes, root edges, and root faces to the supplied lists.
  void resplit(int level,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::list<QuadFace*>& quadface_list,
               std::list<Face*>& face_list) ;

  void empty_split() ;
  int empty_resplit(const std::vector<char>& cellPlan) ;

  /// Replace indexMap with pairs of local fine-cell indices from the current
  /// tree and parentPlan. Refinement or derefinement can produce several pairs
  /// for one cell. Return the number of leaves in parentPlan.
  int32 traverse(const std::vector<char>& parentPlan,
                 vector<pair<int32, int32> >& indexMap) ;

  /// Return true if aCell is a neighbor across face dd with the same face
  /// size, and set nf to its face index. The cells need not have the same
  /// parent.
  bool isSiblingNeighbor(const Prism* aCell, int dd, int &nf) const ;

  /// Find the neighbor across face dd among children of parentCell, and set nf
  /// to its face index. Return 0 if there is no such sibling.
  Prism* getSiblingNeib(int dd, int& nf) ;

  /// Return the parent-face index for side face dd. Requires dd >= 2 and a
  /// nonnull parentCell.
  int parentFace(int dd) {
    if(parentCell->mySplitCode == 1) {
      return dd ;
    }else if(dd == 2) {
      return (whichChild()%(parentCell->nfold)) + 2 ;
    }else if(dd == 5) {
      int childID = (whichChild()%parentCell->nfold) ;
      return (childID== 0)?(parentCell->nfold +1):(childID+1) ;
    }
    return -1 ;
  }

  /// Return true if aCell shares a nonzero area with face dd, and set nf to
  /// its face index. Contact along an edge alone does not count.
  bool isNeighbor(const Prism* aCell, int dd, int& nf) const ;

  /// Find the neighbor across face d within the original Prism, setting nf to
  /// its face index. Return 0 at the original cell boundary. The returned cell
  /// may have children.
  Prism* findNeighbor(int d, int& nf) ;

  /// Return 2 if all cell nodes request derefinement, 1 if any requests
  /// refinement, and 0 otherwise.
  int get_tagged() ;
  /// Return the refinement request from the supplied spacing sources.
  int get_tagged(const vector<source_par>& s) ;
  void setSplitCode(int split_mode, double tol) ;

  /// Make a breadth-first cell refinement plan from this tree.
  std::vector<char> make_cellplan() ;

  /// Make a plan for level levels of isotropic refinement of an original
  /// triangular Prism, using split code 3.
  std::vector<char> make_cellplan(int level) ;

  std::vector<Edge*> get_edges() {

    std::vector<Edge*> edges(3*nfold) ;
    for(int i=0; i<nfold; i++) { edges[i] = gnrlface[0]->edge[i] ; }
    for(int i=nfold; i<2*nfold; i++) { edges[i] = gnrlface[1]->edge[i-nfold] ; }
    for(int i=2*nfold; i<3*nfold; i++) {
      int j = i%nfold ;
      edges[i] = quadface[j]->edge[faceOrient.test(j)?1:3] ;
    }
    return edges ;
  }

  /// Add the splits required by boundary-edge refinement and
  /// Globals::balance_option, using split_mode, then balance the children.
  /// Append new nodes, root edges, and root faces to the supplied lists.
  /// Return true if a split was added.
  bool balance_cell(int split_mode,
                    std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<QuadFace*>& qface_list,
                    std::list<Face*>& gface_list) ;

  void sort_leaves(std::list<Prism*>& v1) ;

  void rebalance_cells(int split_mode,
                       std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list,
                       std::list<QuadFace*>& qface_list,
                       std::list<Face*>& gface_list) ;

  friend void set_prism_faces(const std::vector<Prism*>& cells,
                              std::map<QuadFace*, NeibIndex>& quadfaces,
                              std::map<Face*, NeibIndex>& faces) ;

  friend std::vector<int32> get_c1_prism(const std::vector<char>& cellPlan,
                                         const std::vector<char>& facePlan,
                                         char orientCode,
                                         int faceID) ;

  friend std::vector<char>  extract_prism_face(const  std::vector<char>& cellPlan, int dd) ;

  void print() ;
private:

  // The index of the cell, start with 1
  int32 cellIndex ;

  // 3 for normal prism, 4 for the children of prism when quadface is split.
  char nfold ;

  /// Split code used by split() and empty_split(). Codes 0, 1, 2, and 3 give
  /// 0, 2, nfold, and 2*nfold children.
  char mySplitCode ;

  Face** gnrlface ;

  QuadFace** quadface ;
  //orientCode is not necessary here

  /// The parent of the cell
  Prism *parentCell ;

  /// Array of child-cell pointers.
  Prism **childCell ;

  // If the face in direction RIGHT, LEFT... has been checked
  std::bitset<6> faceMarked ;

  /// Side-face orientation flags: 1 points inward and 0 points outward. The
  /// three side faces of an original Prism point outward; children with nfold
  /// == 4 use these flags for shared faces. gnrlface[0] points inward and
  /// gnrlface[1] points outward throughout refinement.
  std::bitset<4> faceOrient ;

  char tag ;

  //  char whichChild;

  /// Assignment and copying are prohibited
  void operator=(const Prism&) ;
  Prism(const Prism&) ;

private:
  /// Get all the leaves
  void get_leaves(std::vector<Prism*>& leaf_cell) ;

  /// Resize node to 2*nfold entries and fill it with existing corner-node
  /// pointers, first from gnrlface[0], then gnrlface[1].
  void get_nodes(std::vector<Node*>& node) {
    node.resize(2*nfold) ;
    for(int i=0; i<nfold; i++) {
      node[i] = (gnrlface[0]->needReverse[i])?(gnrlface[0]->edge[i]->tail):gnrlface[0]->edge[i]->head ;
      node[i+nfold] = (gnrlface[1]->needReverse[i])?(gnrlface[1]->edge[i]->tail):gnrlface[1]->edge[i]->head ;
    }
  }

  // Get all the 4*3 edges
  // inline void get_edges(Edge** edge){}

  /// Return a new Node at the mean position of the 2*nfold corner nodes. The
  /// caller deletes it.
  Node* simple_center() {

    Node* cellcenter = new Node() ;
    std::vector<Node*> vertices(2*nfold) ;
    get_nodes(vertices) ;
    std::vector<vect3d> nodes(2*nfold) ;
    for(int i=0; i<2*nfold; i++) {
      nodes[i] = vertices[i]->p ;
    }
    cellcenter->p = point_center(nodes) ;
    return cellcenter ;
  }

  /// Return a new cell-center Node using the formula selected by CENTROID. The
  /// wireframe() formula requires existing face-center nodes. The caller
  /// deletes the returned Node.
  Node* centroid() {
    switch(CENTROID) {
    case 0:
      return simple_center() ;
    case 1:
      return wireframe() ;
    default:
      return wireframe() ;
    }
  }

  Node* wireframe() {

    // allocate edgecenter
    std::vector<vect3d> facecenter(nfold+2) ;
    std::vector<double> areas(nfold+2) ;

    // get edge centers
    for(int i=0; i<nfold+2; i++) {
      facecenter[i]= getFaceCenter(i)->p ;
      if(i<2) {
        areas[i] = gnrlface[i]->area() ;
      }else {
        areas[i] = quadface[i-2]->area() ;
      }
    }

    // calculate the mass center of the edge centers
    vect3d p = weighted_center(facecenter, areas) ;
    return new Node(p) ;
  }

  Node* getFaceCenter(int faceID) {
    if(faceID < 2) { return gnrlface[faceID]->child[0]->edge[2]->head ; }
    return quadface[faceID-2]->getCenter() ;
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
                        const const_store<int>& node_remap) ;


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
                        const const_store<int>& node_remap) ;


/// Build the cell using edgePlan, facePlan, and cellPlan, copy the node tags,
/// then apply edgePlan1 and facePlan1 to its boundary. Keep allocated objects
/// in the supplied lists for cleanup.
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
                                const  std::vector<char>& cellNodeTag) ;


/// Build the cell using edgePlan, facePlan, and cellPlan, copy fineCellTag to
/// the fine cells, then apply edgePlan1 and facePlan1 to its boundary. Keep
/// allocated objects in the supplied lists for cleanup.
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
                                     const  std::vector<char>& fineCellTag) ;

// For no restart
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
                        const const_store<int>& node_remap) ;


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
                        const const_store<int>& node_remap) ;


/// Build a Prism with boundary refinement and assign node indices using
/// node_l2f and node_offset. Use face_l2f to order the mesh faces. Append
/// allocated nodes, root edges, and root faces to the supplied lists; the
/// caller deletes the returned Prism.
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
                        std::list<Face*>& gface_list) ;


int general_childID_orient_c2f(int childID_c, char orientCode, int numEdge) ;
int general_childID_orient_f2c(int childID_f, char orientCode, int numEdge) ;

#endif
