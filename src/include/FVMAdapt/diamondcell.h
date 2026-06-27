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
 * @ingroup fvmadapt_elements
 * @brief General-cell refinement helpers based on diamond-cell decomposition.
 */

using std::cerr ;
using std::endl ;
using std::list ;

class Cell ;

/**
 * @brief Diamond-shaped subcell used in general-cell refinement trees.
 * @ingroup fvmadapt_elements
 *
 * Isotropically splitting a general cell creates one DiamondCell for each
 * original cell node. A DiamondCell is formed from each of the original
 * general cell's nodes.
 *
 * `nfold` is the number of faces that meet at the original node that was used
 * to form the DiamondCell. For a DiamondCell, there are two nodes that have
 * the highest number of faces that touch them, all of the other nodes of the
 * DiamondCell have fewer faces.
 *
 * An n-fold DiamondCell has `2n + 2` nodes, `2n` quadrilateral faces, and
 * `4n` edges. Of those faces, `n` share node `0` and `n` share node `1`.
 */
class DiamondCell {
public:

  /// Constructor
  DiamondCell(char m):nfold(m),cellIndex(0),parentCell(0),childCell(0),
              face(new Face*[2*m]), faceOrient(new char[2*m]), faceMarked(0),
              tag(0){}


  /// Destructor
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
  * Delete this cell's children and mark it as a leaf.
  */
  void derefine() ;

  /**
  * Return true if this DiamondCell can collapse its immediate children.
  *
  * Derefinement is allowed when this cell has children, each existing child
  * is tagged `2`, no child is itself refined, and no edge in this cell has
  * been split more than one level below the cell.
  */
  bool needDerefine() ;

  /**
  * Return true if this DiamondCell can collapse its immediate children using
  * cell tags.
  *
  * This is the cell-tag variant of needDerefine(): each existing child must have
  * `getTag() == 2`, no child may already be refined, and no edge in this cell may
  * be split more than one level below the cell.
  */
  bool needDerefine_ctag() ;

  char getTag() const { return tag ; }

  void setTag(char c){ tag=c ; }

  void setCellIndex(int32 cellID){ cellIndex = cellID ; }

  int32 getCellIndex() const { return cellIndex ; }

  char getNfold() const{ return nfold ; }

  int getLevel() const{ return face[0]->edge[0]->level ; }

  /**
  * Return this cell's adaptation state from its node tags.
  *
  * Returns `2` when every node is tagged for derefinement, `1` when at least
  * one node is tagged for refinement, and `0` otherwise. Partial derefinement
  * tags do not derefine the cell. This does not read the DiamondCell's own
  * `tag` field.
  */
  int get_tagged() ;

  /**
  * Return the refinement state implied by PAR-file spacing sources.
  *
  * The cell nodes and minimum edge length are passed to tag_cell(). The
  * `source_par` entries describe geometric spacing sources read from
  * `parfile_par`; this returns `1` when those sources indicate that the cell
  * should be refined and `0` otherwise.
  */
  int get_tagged(const vector<source_par>& s) ;

  /**
  * Return the total number of fine face pieces on this DiamondCell.
  *
  * This sums get_num_leaves() over the cell's `2*nfold` face trees.
  * (used for mxfpc...max faces per cell)
  */
  int get_num_fine_faces() ;

  inline void setParentCell( DiamondCell* parent) { parentCell = parent ; }
  inline DiamondCell* getParentCell() { return parentCell ; }

  inline DiamondCell* getChildCell(int i) const { return childCell[i] ; }

  inline DiamondCell** getChildCell() const { return childCell ; }

  /**
  * Return the parentCell face index that contains this cell's face[faceID].
  *
  * Only face indices `[nfold, 2*nfold)` map to parentCell faces. Indices
  * `[0, nfold)` are internal faces between children created by the split.
  *
  * @param faceID face index in `[nfold, 2*nfold)`.
  * @return Face index in `parentCell`, or `-1` if @p faceID has no parent-face
  *         mapping.
  */
  int parentFace(int faceID) const ;

  /**
  * Split this DiamondCell by `level` additional isotropic refinement levels.
  *
  * If `level <= 0`, this does nothing. Newly created nodes, edges, and
  * faces are appended to the supplied lists.
  */
  void resplit(int level, std::list<Node*>& node_list,
               std::list<Edge*>& edge_list, std::list<Face*>& face_list) ;

  /**
  * Split this DiamondCell once.
  *
  * If the cell is already split, this returns without changing the lists.
  * Otherwise it creates the child-cell geometry and appends newly allocated
  * nodes, edges, and faces to the supplied lists.
  */
  void split(std::list<Node*>& node_list, std::list<Edge*>& edge_list,
             std::list<Face*>& face_list);

  /// This function splits diamondcell isotropically once only define childCell
  void empty_split() ;

  /// This function checks if `aFace` is one of the faces of this cell.
  /// return -1: No
  /// return i in [0, 2*nfold), Yes, face[i] == aFace
  inline int containFace(Face* aFace) {
    for(int i=0; i<2*nfold; i++) {
      if(face[i] == aFace) { return i ; }
    }
    return -1 ;
  }

  /**
  * Return the single DiamondCell adjacent to this cell across local face `mf`.
  *
  * `mf` must be in `[0, 2*nfold)`. Internal faces `[0, nfold)` are resolved by
  * getSiblingNeib(). Faces `[nfold, 2*nfold)` are mapped up through parent
  * DiamondCells until a neighboring branch is found, then mapped back down when
  * needed. This routine is used where the matching neighbor face is represented
  * by one DiamondCell, not by a set of finer neighbors.
  *
  * `nf` is set to the matching local face index on the returned neighbor.
  * Returns `0` when the face has no neighbor inside the parent general cell.
  */
  DiamondCell* findNeighbor(const Cell* aCell,
                            const std::vector<std::vector<Edge*> >& n2e,
                            int mf, int& nf)const;

  /**
  * Return the DiamondCell that shares this cell's internal face.
  *
  * @param aCell General-cell parent used to resolve top-level DiamondCell
  *        neighbors.
  * @param n2e Node-to-edge table used when this cell has no DiamondCell parent.
  * @param mf Local face index on this cell; must be in `[0, nfold)`.
  * @param nf Output local face index on the returned neighboring DiamondCell.
  * @return Neighboring DiamondCell that shares `face[mf]`, or `0` if a
  *         top-level sibling cannot be found.
  */
  DiamondCell* getSiblingNeib(const Cell* aCell,
                              const std::vector<std::vector<Edge*> >& n2e,
                              int mf, int& nf) const ;

  /**
  * Balance this DiamondCell subtree by adding required splits.
  *
  * Existing children are balanced recursively. A leaf is split when any edge is
  * more than one level deeper than the cell, or when the enabled face-balance
  * checks require a split. Newly created nodes, edges, and faces are appended to
  * the supplied lists. Returns true if this cell or a descendant was split.
  */
  bool balance_cell(std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<Face*>& face_list) ;

  void sort_leaves(std::list<DiamondCell*>& v1) ;

  /**
  * Build oriented connectivity for interior fine faces of a split general cell.
  *
  * @param aCell Parent cell used to resolve DiamondCell neighbors.
  * @param cells Leaf DiamondCells in the split tree.
  * @param n2e Parent-cell node-to-edge table used for neighbor lookup.
  * @param fine_face Output `(face, NeibIndex)` pairs. `NeibIndex::c1` and
  *        `NeibIndex::c2` are the oriented local indices of the two leaf cells
  *        adjacent to that face.
  */
  friend void set_general_faces(const Cell* aCell,
                                const std::vector<DiamondCell*>& cells,
                                const std::vector<std::vector<Edge*> >& n2e,
                                std::list<pair<Face*, NeibIndex> >& fine_face);

  /// Find the minimum edge length in a cell(before split)
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

  /// The index of the cell. Cell indices start with 1
  int32 cellIndex ;

  /// The parent of the cell
  DiamondCell *parentCell ;

  // A dynamic array of pointers to children cells
  DiamondCell **childCell ;

  Face** face ;

  /// Orientation of a face.
  /// 0 - Point outward
  /// 1 - Points inward
  /// i.e. faceOrient[i] = 1; else faceOrient[i] = 0;
  char* faceOrient ; //the face points inward or outward

  /// If the face has been checked. size: 2*nfold
  bool* faceMarked ;

  /// Used in findNeighbor() function
  char whichChild ;

  /// Cell tag
  char tag ;

  /// Assignment and copying are prohibited
  void operator=(const DiamondCell&) ;

  DiamondCell(const DiamondCell&) ;

  friend class Cell ;
private:
  /// Get all the leaves
  //   void get_leaves(std::vector<DiamondCell*>& leaf_cell);

  /// Get the 2*nfold+2 nodes
  void get_nodes(std::set<Node*>& node);

  /// Get the 4*nfold edges of the DiamondCell
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

  /// Get the facecenter
  /// Condition: facecenter will be allocated and deallocated by the caller
  /// Precondition: the faces have been splitted
  inline void getFaceCenter(Node** facecenter) {
    for(int i = 0; i < 2*nfold; i++) {
      facecenter[i] = face[i]->child[0]->edge[2]->head ;
    }
  }

  /// Get the edgecenter
  /// Precondition: all the faces and edges has been split
  /// Condition: edgecenter needs to be allocated and deallocated by the caller
  /// edge center is defined from the edges of faces according to the number system and faceOrient
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


/// This class defines the general cell
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

  /// If all children are tagged as 2, remove all children
  bool needDerefine() ;
  bool needDerefine_ctag() ;
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

  /// Precondition: all the faces have been split.
  /// Condition: facecenter is allocated and deallocated by the caller
  inline void getFaceCenter(Node** facecenter) {
    for(int i = 0; i < numFace; i++) {
      facecenter[i] = face[i]->child[0]->edge[2]->head ;
    }
  }

  /// Depth first search
  void sort_leaves(std::list<DiamondCell*>& v1) ;

  /// Find the minimum edge length in a cell(before split)
  inline double get_min_edge_length() {
    double min_length = norm(edge[0]->head->p - edge[0]->tail->p) ;
    for(int i = 1; i < numEdge; i++) {
      min_length = min(min_length, norm(edge[i]->head->p - edge[i]->tail->p)) ;
    }
    return min_length ;
  }

  /// Split a general cell isotropically into DiamondCells
  void split(std::list<Node*>& node_list,
             std::list<Edge*>& edge_list,
             std::list<Face*>& face_list);

  void resplit(const std::vector<char>& cellPlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::list<Face*>& face_list,
               std::vector<DiamondCell*>& cells);

  /// Split a general cell isotropically into DiamondCells only define child
  void empty_split();

  /// Return num_fine_cell
  int32 empty_resplit(const std::vector<char>& cellPlan);

  //  void get_leaves(std::vector<DiamondCell*>& leaf_cell);

  /// After the cell is split into a tree, get the indexMap from current index
  /// to parent index
  int32 traverse(const std::vector<char>& parentPlan,
                 vector<pair<int32, int32> >& indexMap);

  /// For calculating mxfpc
  int get_num_fine_faces();

  /// Find the node2edge of for all nodes, for each node, n2e[i] and n2e[i+1]
  /// share a face.
  std::vector<std::vector<Edge*> > set_n2e();

  /// Find the info for node nindex. rot.size() == n2f.size()
  /// if j=rot[i] >= 0,the orient of n2f[i] is 1, and node[nindex] is the jth node of n2f[i]
  /// if j=rot[i] < 0, the orient of n2f[i] is -1, and  node[nindex] is (-j-1)th node of n2f[i]
  void set_n2f_n2e(std::vector<Face*>& n2f, std::vector<Edge*>& n2e, std::vector<int>& rot, int nindex);

  /// Check if the cell is tagged for refinement, derefinement, or unchanged
  int get_tagged();

  /// Check if the cell is tagged for refinement
  int get_tagged(const vector<source_par>& s) ;

  /// If any edge is more than 1 levels down than my level, split myself, then
  /// balance each child
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

  /// The face points inward or outward.
  /// points inward: faceOrient[i] = 1
  /// points outward: faceOrient[i] = 0
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

/// Build a Cell from Loci data structures, the locations of nodes are defined
/// edges and faces are split according to edgePlan and facePlan
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

/// Build a Cell from Loci data structures, the locations of nodes are defined
/// edges and faces are split according to edgePlan and facePlan
/// and all boundary nodes are tagged
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

/// Build a cell with edgePlan and facePlan, tag the nodes then resplit the
/// edges and faces with edgePlan1 and facePlan1
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

/// Build a cell with edgePlan and facePlan, tag the cells then resplit the
/// edges and faces with edgePlan1 and facePlan1
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

/// Build a Cell from Loci data structures, the locations of nodes are defined
/// and all boundary nodes are tagged
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

/// Build a Cell from Loci data structures, the locations of nodes are defined
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

/// Build a Cell from Loci data structures, the locations of nodes are not defined
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


/// This function first builds a Face f and a Cell cl[f]/cr[f], then splits
/// the Face according to facePlan and splits the Cell according to cellPlan.
/// For each fine face, find the index of fine cell that it belongs to.
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

/// This function merges two isotropical facePlan
std::vector<char> merge_faceplan(std::vector<char>& planl, std::vector<char>& planr, int numNodes);

/// This function extracts facePlan from a cellPlan
std::vector<char> extract_general_face(const Entity* lower, int lower_size,
                                       const Entity* upper, int upper_size,
                                       const Entity* boundary_map, int boundary_map_size,
                                       const const_multiMap& face2node,
                                       const const_multiMap& face2edge,
                                       const const_MapVec<2>& edge2node,
                                       const std::vector<char>& cellPlan,
                                       Entity ff, const const_store<int>& node_remap);

#endif
