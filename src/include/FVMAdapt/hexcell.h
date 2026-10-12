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

/**
 * @file hexcell.h
 *
 * HexCell refinement and fine-face connectivity.
 */

/**
 * Return the local fine-cell index adjacent to each fine face on face findex.
 * Indices start at 1 within the original HexCell and follow the leaf order of
 * facePlan.
 *
 * findex is in [0, 6). orientCode maps face2node order to the face order in
 * the HexCell. Use integer face ranges to match fine faces to the cells
 * produced by cellPlan.
 */
std::vector<int32> get_c1_hex(const std::vector<char>& cellPlan,
                              const std::vector<char>& facePlan,
                              char orientCode,
                              char findex) ;

/**
 * Hexahedral cell with directional refinement. mySplitCode uses bits 4, 2,
 * and 1 for the local xi, eta, and zeta directions. Splitting one, two, or
 * three directions creates two, four, or eight children.
 *
 * The face and child numbering is used by the extraction tables in tables.h.
 */
class HexCell
{
public:

  HexCell():cellIndex(0), mySplitCode(0), face(0), parentCell(0),
            childCell(0),tag(0){}

  HexCell(QuadFace** f):cellIndex(0), mySplitCode(0), face(f), parentCell(0),
                        childCell(0),tag(0){}
  ~HexCell(){
    if(childCell != 0) {
      for(int i = 0; i < numChildren(); i++) {
        if(childCell[i] != 0) {
          delete  childCell[i] ;
          childCell[i] = 0 ;
        }
      }
      delete[] childCell ;
      childCell = 0 ;
    }
    parentCell = 0 ;

    if(face != 0) {
      delete[] face ;
      face = 0 ;
    }
  }

  /// Delete the children and reset mySplitCode to 0. This does not check tags
  /// or derefinement eligibility.
  void derefine() ;
  /// Return true if this cell has children, every child is a leaf requesting
  /// derefinement through its node tags, and no boundary edge has
  /// grandchildren. This only tests eligibility.
  bool needDerefine() ;
  /// Test the same leaf and edge conditions as needDerefine(), using each
  /// child's cell tag (getTag() == 2) instead of its node tags.
  bool needDerefine_ctag() ;
  char getTag() const { return tag ; }
  void setTag(char c) { tag=c ; }

  int32 getCellIndex() const { return cellIndex ; }
  char getMySplitCode() const { return mySplitCode ; }
  HexCell* getChildCell(int i) const { return childCell[i] ; }
  HexCell* getParentCell() { return parentCell ; }

  /// Return 2 if all cell nodes request derefinement, 1 if any requests
  /// refinement, and 0 otherwise.
  int get_tagged() ;
  /// Return the refinement request from the supplied spacing sources.
  int get_tagged(const vector<source_par>& s) ;

  // return a splitCode
  // find average_edge_length in XX , YY and ZZ directions
  // find min_edge_length in all directions

  /// Choose mySplitCode using split_mode and tol. Automatic directional
  /// splitting compares ratios of average edge lengths in the local directions
  /// with Globals::factor.
  void setSplitCode(int split_mode, double tol) ;

  int getLevel(NORMAL_DIRECTION d) const {
    switch(d) {
    case XX: // x direction
      return face[FRONT]->edge[0]->level ;
      break ;
    case YY: // y direction
      return face[RIGHT]->edge[0]->level ;
      break ;
    case ZZ: // z direction
      return face[RIGHT]->edge[1]->level ;
      break ;
    default:
      cerr << "WARNING: illegal levelID" << endl ;
      break ;
    }
    return 0 ;
  }

  std::vector<Edge*> get_edges() {
    std::vector<Edge*> edges(12) ;
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
    return edges ;
  }

  int numChildren() const {
    switch(mySplitCode) {
    case 0:
      return 0 ;
    case 1:
    case 2:
    case 4:
      return 2 ;
    case 3:
    case 5:
    case 6:
      return 4 ;
    case 7:
      return 8 ;
    default:
      cerr << "WARNING: illegal split code" << endl ;
      break ;
    }
    return -1 ;
  }

  /// Find numChildren without actually building the cell
  int numChildren(char splitCode) const {
    switch(splitCode) {
    case 0:
      return 0 ;
    case 1:
    case 2:
    case 4:
      return 2 ;
    case 3:
    case 5:
    case 6:
      return 4 ;
    case 7:
      return 8 ;
    default:
      cerr << "WARNING: illegal split code" << endl ;
      break ;
    }
    return -1 ;
  }

  /// Find num_fine_cells without actually building the tree
  int32 num_fine_cells( const std::vector<char>& cellPlan) const ;

  /// Return the number of leaf faces on the cell boundary, used to compute the
  /// maximum faces per cell (mxfpc).
  int get_num_fine_faces() const ;
  double get_min_edge_length() ;

  void split(std::list<Node*>& node_list,
             std::list<Edge*>& edge_list,
             std::list<QuadFace*>& face_list) ;

  /// Create the child HexCells selected by mySplitCode, without nodes, edges,
  /// or faces.
  void empty_split() ;

  /// Apply cellPlan to the child-cell structure, assign local leaf indices
  /// starting at 1, and return the leaf count. No split geometry is created.
  int empty_resplit(const std::vector<char>& cellPlan) ;

  /// Replace indexMap with pairs of local fine-cell indices from the current
  /// tree and parentPlan. Refinement or derefinement can produce several pairs
  /// for one cell. Return the number of leaves in parentPlan.
  int32 traverse(const std::vector<char>& parentPlan,
                 vector<pair<int32, int32> >& indexMap) ;

  void resplit(const std::vector<char>& cellPlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::list<QuadFace*>& face_list,
               std::vector<HexCell*>& cells) ;

  //used in make_hex_cellplan.cc
  void resplit(int level,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::list<QuadFace*>& face_list) ;

  /// Return the index i for which face[i] == aFace, or -1 if no face matches.
  int containFace(QuadFace* aFace) {
    for(int i=0; i<6; i++) {
      if(face[i] == aFace) { return i ; }
    }
    return -1 ;
  }

  /// Return true if aCell is a neighbor in direction dd with the same face
  /// size. The cells need not have the same parent.
  bool isSiblingNeighbor(const HexCell* aCell, DIRECTION dd) const ;

  /// Return true if aCell shares a nonzero face area in direction dd. Contact
  /// along an edge alone does not count.
  bool isNeighbor(const HexCell* aCell, DIRECTION dd) const ;

  /// Find the face neighbor in direction d within the original HexCell. Return
  /// 0 at its boundary. The returned cell may have children.
  HexCell* findNeighbor(DIRECTION d) ;

  /// Make a breadth-first cell refinement plan from this tree.
  std::vector<char> make_cellplan() ;

  /// Make a plan for level levels of isotropic refinement of an original
  /// HexCell, using split code 7.
  std::vector<char> make_cellplan(int level) ;

  /// Add the splits required by boundary-edge refinement and
  /// Globals::balance_option, using split_mode, then balance the children.
  /// Append new nodes, root edges, and root faces to the supplied lists.
  /// Return true if a split was added.
  bool balance_cell(int split_mode,
                    std::list<Node*>& node_list,
                    std::list<Edge*>& edge_list,
                    std::list<QuadFace*>& face_list) ;

  void sort_leaves(std::list<HexCell*>& v1) ;

  void rebalance_cells(int split_mode,
                       std::list<Node*>& node_list,
                       std::list<Edge*>& edge_list,
                       std::list<QuadFace*>& face_list) ;

  /// Find interior fine faces and store their NeibIndex values in faces. c1
  /// and c2 are local fine-cell indices ordered by face orientation.
  friend void set_hex_faces(const std::vector<HexCell*>& cells,
                            std::map<QuadFace*, NeibIndex>& faces) ;

  friend std::vector<int32> get_c1_hex(const std::vector<char>& cellPlan,
                                       const std::vector<char>& facePlan,
                                       char orientCode,
                                       char findex) ;
private:

  /// Local fine-cell index, starting at 1. Zero is used before numbering or
  /// for non-leaf state.
  int32 cellIndex ;

  /// Local-direction split mask: bit 4 splits xi, bit 2 splits eta, and bit 1
  /// splits zeta. Code 0 requests no split; code 3 splits eta and zeta, for
  /// example.
  char mySplitCode ;

  /// 6 faces,  pointing to positive xi, eta, or zeta direction
  /// the numbering of face:
  /// 0(RIGHT): xi = 1
  /// 1(LEFT):  xi = 0
  /// 2(FRONT): eta = 1
  /// 3(BACK):  eta  = 0
  /// 4(UP):    zeta =1
  /// 5(DOWN):  zeta = 0
  QuadFace** face ;

  /// The parent of the cell
  HexCell *parentCell ;

  // A dynamic array of pointers to children cells
  HexCell **childCell ;

  /// Faces visited by set_hex_faces().
  std::bitset<6> faceMarked ;

  char tag ;

  //char whichChild ;

  /// Assignment and copying are prohibited
  void operator=(const HexCell&) ;
  HexCell(const HexCell&) ;

private:

  /// Fill eight entries of node with existing corner-node pointers. The caller
  /// must size the vector to at least eight entries.
  void get_nodes(std::vector<Node*>& node) {
    for(int i=0; i<4; i++) {
      node[i] = face[0]->getNode(i) ;
      node[i+4] = face[1]->getNode(i) ;
    }
  }

  /// Return a new Node at the mean position of the eight corner nodes. The
  /// caller deletes it.
  Node* simple_center() {
    Node* cellcenter = new Node() ;
    std::vector<Node*> vertices(8) ;
    get_nodes(vertices) ;
    std::vector<vect3d> nodes(8) ;
    for(int i=0; i<8; i++) {
      nodes[i] = vertices[i]->p ;
    }
    cellcenter->p = point_center(nodes) ;
    return cellcenter ;
  }

  Node* wireframe() {

    // allocate edgecenter
    std::vector<vect3d> facecenter(6) ;
    std::vector<double> areas(6) ;

    // get edge centers
    for(int i=0; i<6; i++){
      facecenter[i]= getFaceCenter(i)->p ;
      areas[i] = face[i]->area() ;
    }

    // calculate the mass center of the edge centers
    vect3d p = weighted_center(facecenter, areas) ;
    return new Node(p) ;
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

  /// Fill the caller's array of six pointers with existing face-center nodes.
  /// Each face must have split code 3; no nodes are allocated or transferred.
  void getFaceCenter(Node** facecenter) {
    for(int i=0; i<6; i++) {
      facecenter[i] = face[i]->getCenter() ;
    }
  }

  Node* getFaceCenter(int faceID) {
    return face[faceID]->getCenter() ;
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

/// Build the cell using edgePlan, facePlan, and cellPlan, copy the node tags,
/// then apply edgePlan1 and facePlan1 to its boundary. Keep allocated objects
/// in the supplied lists for cleanup.
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

/// Build the cell using edgePlan, facePlan, and cellPlan, copy fineCellTag to
/// the fine cells, then apply edgePlan1 and facePlan1 to its boundary. Keep
/// allocated objects in the supplied lists for cleanup.
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

/// Build a HexCell with boundary refinement and assign node indices using
/// node_l2f and node_offset. Use face_l2f to order the mesh faces. Append
/// allocated nodes, root edges, and root faces to the supplied lists; the
/// caller deletes the returned HexCell.
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

/// Collects 6 faces of a hexcell from loci data structures
Array<Entity, 6> collect_hex_faces(const Entity*  lower,int lower_size,
                                   const Entity* upper,int upper_size,
                                   const Entity* boundary_map,int boundary_map_size,
                                   const Array<char, 6>& hex2face, const const_store<int>& node_remap);

/// Collects 8 vertices of a hexcell from loci data structures
Array<Entity, 8> collect_hex_vertices(const const_multiMap& face2node, const Array<Entity, 6>& faces,
                                      const Array<char, 8>& hex2node);

/// Return the 12 edge entities in HexCell order and fill needReverse with
/// their directions relative to edge2node.
Array<Entity, 12> collect_hex_edges(const Array<Entity, 6>& faces, const Array<Entity, 8>& hex_vertices,
                                    const const_multiMap& face2edge, const const_MapVec<2>& edge2node,
                                    Array<bool,12>& needReverse);

//this function will define face2node for each fine faces and write them out,
//at the same time, mxppf(max num of points per face) will be updated
// void  write_hex_inner_faces(const std::list<pair<QuadFace*, NeibIndex> >& faces,
//                             int cell_offset, int& mxppf,std::ofstream& ofile);

std::vector<char>  extract_hex_face(const  std::vector<char>& cellPlan,  DIRECTION dd);
#endif
