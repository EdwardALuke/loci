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
#include "defines.h"
#include "face.h"
#include "read_par.h"

/**
 * @file diamondcell.h
 *
 * General Cell and DiamondCell refinement operations.
 */

using std::cerr ;
using std::endl ;
using std::list ;

class Cell ;

/**
 * Cell created when a general Cell is split at its vertices. Each original
 * vertex produces a DiamondCell whose nfold is the number of faces meeting at
 * that vertex.
 *
 * An nfold DiamondCell has 2*nfold quadrilateral faces, 2*nfold+2 nodes, and
 * 4*nfold edges. Nodes 0 and 1 each meet nfold edges; each other node meets
 * three edges.
 */
class DiamondCell {
public:

  DiamondCell(char m):nfold(m),cellIndex(0),parentCell(0),childCell(0),
              face(new Face*[2*m]), faceOrient(new char[2*m]), faceMarked(0),
              tag(0){}


  ~DiamondCell() {
    if(childCell != 0) {
      for(int i=0; i<2*nfold+2; i++) {
        if(childCell[i] != 0) {
          delete  childCell[i] ;
          childCell[i] = 0 ;
        }
      }
      delete [] childCell ;
      childCell = 0 ;
    }

    if(faceMarked != 0) {
      delete[] faceMarked ;
      faceMarked = 0 ;
    }

    if(faceOrient!=0) {
      delete [] faceOrient ;
      faceOrient =0 ;
    }

    if(face != 0) {
      delete [] face ;
      face = 0 ;
    }

    parentCell = 0 ;
  }

  /**
   * Delete this cell's children and make it a leaf. This does not check
   * whether derefinement is allowed.
   */
  void derefine() ;

  /**
   * Return true if this cell has children and they can be removed. Every
   * existing child must be a leaf with get_tagged() == 2, and no boundary
   * edge may have grandchildren. get_tagged() uses node tags.
   */
  bool needDerefine() ;

  /**
   * Test the same leaf and edge conditions as needDerefine(), using each
   * child's cell tag (getTag() == 2) instead of its node tags.
   */
  bool needDerefine_ctag() ;

  char getTag() const { return tag ; }

  void setTag(char c){ tag=c ; }

  void setCellIndex(int32 cellID){ cellIndex = cellID ; }

  int32 getCellIndex() const { return cellIndex ; }

  char getNfold() const{ return nfold ; }

  int getLevel() const{ return face[0]->edge[0]->level ; }

  /**
   * Return 2 if all nodes request derefinement, 1 if any node requests
   * refinement, and 0 otherwise. This reads node tags, not the cell's tag
   * member.
   */
  int get_tagged() ;

  /**
   * Return the refinement request from tag_cell(), using the cell vertices,
   * the shortest boundary edge, and the supplied source_par entries.
   */
  int get_tagged(const vector<source_par>& s) ;

  /**
   * Return the number of leaf faces on the cell boundary, used to compute the
   * maximum faces per cell (mxfpc).
   */
  int get_num_fine_faces() ;

  inline void setParentCell( DiamondCell* parent) { parentCell = parent ; }
  inline DiamondCell* getParentCell() { return parentCell ; }

  inline DiamondCell* getChildCell(int i) const { return childCell[i] ; }

  inline DiamondCell** getChildCell() const { return childCell ; }

  /**
   * Return the index of the parent face containing face[faceID]. parentCell
   * must be nonnull.
   *
   * Indices outside [nfold, 2*nfold) return -1 because they do not name
   * parent-boundary faces. An index in that range must have a matching parent
   * face; otherwise the routine aborts.
   */
  int parentFace(int faceID) const ;

  /**
   * Split toward the level getLevel() + level, appending new nodes, root
   * edges, and root faces to the supplied lists. Existing splits are
   * retained; level <= 0 does nothing.
   */
  void resplit(int level, std::list<Node*>& node_list,
               std::list<Edge*>& edge_list, std::list<Face*>& face_list) ;

  /**
   * Split this DiamondCell once, appending new nodes, root edges, and root
   * faces to the supplied lists. If childCell already exists, do nothing.
   */
  void split(std::list<Node*>& node_list, std::list<Edge*>& edge_list,
             std::list<Face*>& face_list);

  /// Create the child DiamondCells for one isotropic split, without nodes,
  /// edges, or faces.
  void empty_split() ;

  /// Return the index i for which face[i] == aFace, or -1 if no face matches.
  inline int containFace(Face* aFace) {
    for(int i=0; i<2*nfold; i++) {
      if(face[i] == aFace) { return i ; }
    }
    return -1 ;
  }

  /**
   * Find the neighbor across face mf within aCell, where mf is in [0,
   * 2*nfold). Return a coarser leaf or a neighbor at this cell's level; the
   * latter may have children.
   *
   * Set nf to the neighbor's face index. Return 0 if no neighbor is found
   * inside the original general Cell; nf is then unspecified.
   */
  DiamondCell* findNeighbor(const Cell* aCell,
                            const std::vector<std::vector<Edge*> >& n2e,
                            int mf, int& nf)const;

  /**
   * Find the sibling across internal face mf, where mf is in [0, nfold), and
   * set nf to its face index.
   *
   * For a cell directly under aCell, use the original cell's node-to-edge
   * table n2e. Otherwise use parentCell and whichChild. Return 0 if a sibling
   * in aCell cannot be found.
   */
  DiamondCell* getSiblingNeib(const Cell* aCell,
                              const std::vector<std::vector<Edge*> >& n2e,
                              int mf, int& nf) const ;

  /**
   * Add the splits required by boundary-edge refinement and the enabled face
   * checks in Globals::balance_option. Balance children recursively and
   * append new nodes, root edges, and root faces to the supplied lists.
   * Return true if this cell or a descendant was split.
   */
  bool balance_cell(std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<Face*>& face_list) ;

  void sort_leaves(std::list<DiamondCell*>& v1) ;

  /**
   * Append interior fine faces and their NeibIndex values to the output list.
   * c1 and c2 are local fine-cell indices ordered by face orientation. Use
   * aCell and its node-to-edge table n2e to find neighbors of the leaf cells.
   */
  friend void set_general_faces(const Cell* aCell,
                                const std::vector<DiamondCell*>& cells,
                                const std::vector<std::vector<Edge*> >& n2e,
                                std::list<pair<Face*, NeibIndex> >& fine_face);

  /// Return the minimum length of the cell boundary edges.
  inline double get_min_edge_length() {
    std::set<Edge*> edge ;
    get_edges(edge) ;
    std::set<Edge*>::const_iterator cur_edge = edge.begin() ;

    double min_length = norm((*cur_edge)->head->p - (*cur_edge)->tail->p) ;
    for(cur_edge= edge.begin(); cur_edge != edge.end(); cur_edge++) {
      min_length = min(min_length, norm((*cur_edge)->head->p - (*cur_edge)->tail->p)) ;
    }
    return min_length ;
  }

private:
  /// An nfold diamondcell will have 2*nfold faces and 2*nfold+2 vertices,
  /// node 0 and node 1 have nfold edges, other vertices have 3 edges
  char nfold;

  /// Local fine-cell index, starting at 1. Zero is used before numbering or
  /// for non-leaf state.
  int32 cellIndex ;

  /// Parent DiamondCell, or null for a child directly under the original Cell.
  DiamondCell *parentCell ;

  // A dynamic array of pointers to children cells
  DiamondCell **childCell ;

  Face** face ;

  /// Face orientation: 0 points outward from this cell; 1 points inward.
  char* faceOrient ; //the face points inward or outward

  /// Faces visited by set_general_faces(); 2*nfold entries.
  bool* faceMarked ;

  /// Index of this child in its parent's child array, used for neighbor
  /// lookup.
  char whichChild ;

  /// Cell tag: 0 unchanged, 1 refine, 2 derefine.
  char tag ;

  /// Assignment and copying are prohibited
  void operator=(const DiamondCell&) ;

  DiamondCell(const DiamondCell&) ;

  friend class Cell ;
private:
  //   void get_leaves(std::vector<DiamondCell*>& leaf_cell);

  /// Insert the existing 2*nfold+2 node pointers into node. The output set
  /// must initially be empty.
  void get_nodes(std::set<Node*>& node);

  /// Insert the existing 4*nfold edge pointers into edge. The output set must
  /// initially be empty.
  void get_edges(std::set<Edge*>& edge);

  /// Calculate the centroid of the DiamondCell, it's defined as the mean value
  /// of facecenters.
  /// Precondition: all the faces have been split
  inline Node* simple_center() {
    Node* cellcenter = new Node() ;
    std::vector<vect3d> facecenter(2*nfold) ;
    for(int i = 0; i < 2*nfold ; i++) {
      facecenter[i] = face[i]->child[0]->edge[2]->head->p ;
    }
    cellcenter->p = point_center(facecenter) ;
    return cellcenter ;
  }

  /**
  * Create a center node from area-weighted split face centers.
  *
  * Precondition: each of the `2*nfold` faces has already been split so its
  * center point is available from `face[i]->child[0]->edge[2]->head`.
  * The caller owns the returned node.
  */
  inline Node* wireframe() {

    // allocate facecenter
    std::vector<vect3d> facecenter(2*nfold) ;
    std::vector<double> areas(2*nfold) ;

    // get face centers
    for(int i = 0; i < 2*nfold; i++) {
      facecenter[i]= face[i]->child[0]->edge[2]->head->p ;
      areas[i] = face[i]->area() ;
    }

    // calculate the mass center of the face centers
    vect3d p = weighted_center(facecenter, areas) ;

    return new Node(p) ;
  }

  inline Node* centroid() {
    switch(CENTROID) {
    case 0:
      return simple_center() ;
    case 1:
      return wireframe() ;
    default:
      return wireframe() ;
    }
  }

  /// Fill the caller's array of 2*nfold pointers with the existing face-center
  /// nodes. All faces must already be split; no nodes are allocated or
  /// transferred.
  inline void getFaceCenter(Node** facecenter) {
    for(int i = 0; i < 2*nfold; i++) {
      facecenter[i] = face[i]->child[0]->edge[2]->head ;
    }
  }

  /// Fill the caller's array of 4*nfold pointers with the existing edge
  /// midpoints, in the DiamondCell edge order. All edges must already be
  /// split; no nodes are allocated or transferred.
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


/// General polyhedral cell represented by nodes, edges, and Face objects.
class Cell{
public:
  // constructors
  Cell(int nd, int ne, int nf, Node** n, Edge** e, Face** f, char* fo):
    numNode(nd), numEdge(ne), numFace(nf), node(n), edge(e), face(f),
    faceOrient(fo),child(0){}

  Cell():node(0), edge(0), face(0), faceOrient(0),child(0){}

  // destructor
  ~Cell() {
    if(child!= 0) {
      for(int i = 0; i < numNode; i++) {
        if(child[i] != 0) {
          delete child[i] ;
          child[i] = 0 ;
        }
      }
      delete [] child ;
      child = 0 ;
    }
    if(faceOrient !=0) {
      delete [] faceOrient ;
      faceOrient = 0 ;
    }
    if(face!=0) {
      delete [] face ;
      face = 0 ;
    }

    if(edge !=0) {
      delete [] edge ;
      edge = 0 ;
    }
    if(node !=0) {
      delete [] node ;
      node = 0 ;
    }
  }

  /// Return true if this cell has children, every existing child is a leaf
  /// requesting derefinement through its node tags, and no boundary edge has
  /// grandchildren. This only tests eligibility; it does not remove children.
  bool needDerefine() ;
  /// Test the same leaf and edge conditions as needDerefine(), using each
  /// child's cell tag (getTag() == 2) instead of its node tags.
  bool needDerefine_ctag() ;
  /// Delete the children. This does not check tags or derefinement
  /// eligibility.
  void derefine() ;


  /// Center of the cell, defined as the mean value of the facecenter.
  /// Precondition: all the face and edge have been splitted
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
    // allocate edgecenter
    std::vector<vect3d> facecenter(numFace);
    std::vector<double> areas(numFace);

    // get edge centers
    for(int i = 0; i < numFace; i++){
      facecenter[i]=face[i]->child[0]->edge[2]->head->p;
      areas[i] = face[i]->area();
    }

    // calculate the mass center of the edge centers
    vect3d p = weighted_center(facecenter, areas);
    return new Node(p);
  }

  inline Node* centroid() {
    switch(CENTROID) {
    case 0:
      return simple_center() ;
    case 1:
      return wireframe() ;
    default:
      return wireframe() ;
    }
  }

  /// Fill the caller's array of numFace pointers with existing face-center
  /// nodes. All faces must already be split; no nodes are allocated or
  /// transferred.
  inline void getFaceCenter(Node** facecenter) {
    for(int i = 0; i < numFace; i++) {
      facecenter[i] = face[i]->child[0]->edge[2]->head ;
    }
  }

  /// Append DiamondCell leaves in depth-first child order. An unsplit Cell
  /// contributes no DiamondCell leaves.
  void sort_leaves(std::list<DiamondCell*>& v1) ;

  /// Return the minimum length of the cell boundary edges.
  inline double get_min_edge_length() {
    double min_length = norm(edge[0]->head->p - edge[0]->tail->p) ;
    for(int i = 1; i < numEdge; i++) {
      min_length = min(min_length, norm(edge[i]->head->p - edge[i]->tail->p)) ;
    }
    return min_length ;
  }

  /// Split this Cell into one DiamondCell per vertex. Append new nodes, root
  /// edges, and root faces to the supplied lists.
  void split(std::list<Node*>& node_list,
             std::list<Edge*>& edge_list,
             std::list<Face*>& face_list);

  void resplit(const std::vector<char>& cellPlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::list<Face*>& face_list,
               std::vector<DiamondCell*>& cells);

  /// Create one child DiamondCell per vertex, without constructing split
  /// geometry.
  void empty_split();

  /// Apply cellPlan to the child-cell structure, assign local leaf indices
  /// starting at 1, and return the leaf count. No split geometry is created.
  int32 empty_resplit(const std::vector<char>& cellPlan);

  //  void get_leaves(std::vector<DiamondCell*>& leaf_cell);

  /// Replace indexMap with pairs of local fine-cell indices from the current
  /// tree and parentPlan. Refinement or derefinement can produce several pairs
  /// for one cell. Return the number of leaves in parentPlan.
  int32 traverse(const std::vector<char>& parentPlan,
                 vector<pair<int32, int32> >& indexMap);

  /// Return the number of leaf faces on the cell boundary, used to compute the
  /// maximum faces per cell (mxfpc).
  int get_num_fine_faces();

  /// For each node, return its incident edges in an order where consecutive
  /// edges share a face.
  std::vector<std::vector<Edge*> > set_n2e();

  /// Collect and order the faces and edges meeting node[nindex]. Consecutive
  /// faces n2f[i] and n2f[i+1] share n2e[i], with the last face followed by
  /// the first.
  ///
  /// For an outward face (faceOrient == 0), rot stores the node's position in
  /// that face. For an inward face (faceOrient == 1), rot stores -position-1.
  void set_n2f_n2e(std::vector<Face*>& n2f, std::vector<Edge*>& n2e, std::vector<int>& rot, int nindex);

  /// Return 2 if all cell nodes request derefinement, 1 if any requests
  /// refinement, and 0 otherwise.
  int get_tagged();

  /// Return the refinement request from the supplied spacing sources.
  int get_tagged(const vector<source_par>& s) ;

  /// Add the splits required by boundary-edge refinement and
  /// Globals::balance_option, then balance the children. Append new nodes,
  /// root edges, and root faces to the supplied lists. Return true if a split
  /// was added.
  bool balance_cell(std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<Face*>& face_list) ;

  void rebalance_cells(std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list,
                       std::list<Face*>& face_list) ;

  std::vector<char> make_cellplan() ;

public:
  int numNode ;
  int numEdge ;
  int numFace ;
  Node** node ;
  Edge** edge ;
  Face** face ;

  /// Face orientation: 0 points outward from this cell; 1 points inward.
  char* faceOrient ;

  // A dynamic array of pointers to children cells
  DiamondCell **child ;
};

int find_face_index(const Entity* lower, int lower_size,
                    const Entity* upper, int upper_size,
                    const Entity* boundary_map, int boundary_map_size,
                    const const_multiMap& face2node,
                    Entity f,
                    const const_store<int>& node_remap) ;

/// Build a Cell from the faces in lower, upper, and boundary_map. Order its
/// nodes, edges, and faces using node_remap, set node positions from pos, and
/// split boundary edges and faces according to edgePlan and facePlan.
///
/// The caller deletes the returned Cell. Allocated nodes, root edges, and root
/// faces are appended to bnode_list, edge_list, and face_list for cleanup.
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

/// Build a Cell with boundary refinement and copy posTag and nodeTag to its
/// boundary nodes. Append allocated nodes, root edges, and root faces to the
/// supplied lists; the caller deletes the returned Cell.
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

/// Build the cell using edgePlan, facePlan, and cellPlan, copy the node tags,
/// then apply edgePlan1 and facePlan1 to its boundary. Keep allocated objects
/// in the supplied lists for cleanup.
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

/// Build the cell using edgePlan, facePlan, and cellPlan, copy fineCellTag to
/// the fine cells, then apply edgePlan1 and facePlan1 to its boundary. Keep
/// allocated objects in the supplied lists for cleanup.
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

/// Build a Cell with node coordinates from pos and copy posTag to its boundary
/// nodes. Append allocated nodes, root edges, and root faces to the supplied
/// lists; the caller deletes the returned Cell.
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

/// Build an unsplit Cell with node coordinates from pos. Append allocated
/// nodes, root edges, and root faces to the supplied lists; the caller deletes
/// the returned Cell.
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

/// Build the Cell topology without setting node coordinates. Append allocated
/// nodes, root edges, and root faces to the supplied lists; the caller deletes
/// the returned Cell.
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


/// Return the local fine-cell index adjacent to each fine face on face ff.
/// Indices start at 1 within the original Cell and follow the leaf order of
/// facePlan.
///
/// Build the Cell from the mesh maps, apply cellPlan and facePlan without
/// creating split geometry, and follow the selected face through the cell
/// tree. The plans must describe compatible cell and face subdivisions.
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

/// Merge two isotropic Face refinement plans.
std::vector<char> merge_faceplan(std::vector<char>& planl, std::vector<char>& planr, int numNodes);

/// Extract the refinement plan for one face of a general Cell from cellPlan.
std::vector<char> extract_general_face(const Entity* lower, int lower_size,
                                       const Entity* upper, int upper_size,
                                       const Entity* boundary_map, int boundary_map_size,
                                       const const_multiMap& face2node,
                                       const const_multiMap& face2edge,
                                       const const_MapVec<2>& edge2node,
                                       const std::vector<char>& cellPlan,
                                       Entity ff, const const_store<int>& node_remap);

#endif
