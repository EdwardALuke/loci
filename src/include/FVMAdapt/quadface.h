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
#ifndef QUADFACE_H
#define QUADFACE_H
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include <bitset>
#include "node_edge.h"

using std::ofstream ;
using std::bitset ;
using std::cout ;
using std::endl ;

struct Range2d ;

/**
 * @file quadface.h
 * @ingroup fvmadapt_elements
 * @brief Directional quadrilateral face tree and orientation helpers.
 */

// f2c orient functions are used when a quadface is built as in cell and the
// facePlan is for the face defined by face2node
/**
 * @brief Convert a face-plan split code to the cell-local quadface axes.
 *
 * Leaf (`0`) and two-axis (`3`) split codes are unchanged. For orientations
 * that swap the local face axes, face-local y (`1`) and x (`2`) splits are
 * exchanged.
 */
char orient_splitCode_f2c(char splitCode, char orientCode) ;
char orient_childID_f2c(char childID, char orientCode, char splitCode) ;
char orient_edgeID_f2c(char edgeID, char orientCode) ;

/// c2f orient functions are used when the face is defined by face2node
/// and facePlan is extracted from cell
char orient_edgeID_c2f(char edgeID, char orientCode);

/// For get_c1_hex and get_c1_prism, when a cell is empty_split, faceMap is
/// created and face is split to generate leaves. By looking up faceMap, c1 can
/// be computed. Use Range2d to avoid fully splitting the cell and face.
std::vector<int32> contain_2d(const std::vector<pair<Range2d, int32> >& faceMap,
                              const std::vector<Range2d>& leaves);

/**
 * @brief Directional quadrilateral-face refinement tree.
 * @ingroup fvmadapt_elements
 *
 * QuadFace is the directional face tree used by @ref HexCell and @ref Prism
 * for anisotropic refinement of quadrilateral faces. It stores each face with
 * a fixed local corner ordering and zero-based `edge[]` ordering:
 *
 * <pre>
 *                         edge[2] (3 -> 2)
 *          node 3  ------------------------------>  node 2
 *            ^                                      ^
 *            |                                      |
 *            | edge[3] (0 -> 3)                    | edge[1] (1 -> 2)
 *            |                                      |
 *            |                                      |
 *          node 0  ------------------------------>  node 1
 *                         edge[0] (0 -> 1)
 * </pre>
 *
 * The corner labels are face-local node indices, not global mesh node IDs. The
 * edge arrows show each stored boundary-edge direction.
 *
 * For split code `3`, the face is split in both local directions. The four
 * `child[]` entries are the quadrants of the face, while `childx[]` and
 * `childy[]` provide the column and row views of the same split:
 *
 * <pre>
 *                              local x direction
 *                    childx[0]                  childx[1]
 *              +-----------------------+-----------------------+
 * childy[1]    |       child[1]        |       child[3]        |
 *              +-----------------------+-----------------------+
 * childy[0]    |       child[0]        |       child[2]        |
 *              +-----------------------+-----------------------+
 *                    local y increases upward
 * </pre>
 *
 * In other words, `childx[0]` is the left half, `childx[1]` is the right half,
 * `childy[0]` is the lower half, and `childy[1]` is the upper half.
 *
 * Quad face plans use code `0` for a leaf, `1` for a face-local y split, `2`
 * for a face-local x split, and `3` for both directions. The `childx`,
 * `childy`, and `child` arrays intentionally encode overlapping views of the
 * same directional split tree.
 *
 * @warning Child and orientation ordering are coupled to hex/prism face
 * extraction and merge helpers. Do not change child layout without updating the
 * orientation helpers and plan-extraction tables together.
 */
class QuadFace{
public:
  /// Constructors
  QuadFace(int numEdge):edge(new Edge*[numEdge]),child(0),childx(0), childy(0),code(char(0)){}

  //Constructor used for empty_split
  QuadFace():edge(0), child(0),childx(0), childy(0),code(char(0)){}
  QuadFace( Edge** e):edge(e),child(0), childx(0), childy(0),code(char(0)){}

 //  //destructor, it works this way without memory leakage
//   ~QuadFace(){

//     if(this != 0){
//       switch(code){
//       case 3:

//         if(childx != 0){
//           for(int i = 0; i < 2; i++){
//             if(childx[i] != 0){
//               //first detangle the pointer, cut offset all the children,  prevent address alias
//               if(childx[i]->childy !=0){
//                 childx[i]->childy[0] = 0;
//                 childx[i]->childy[1] = 0;
//               }

//                delete childx[i];
//               childx[i] = 0;
//             }
//           }
//           delete[] childx;
//           childx = 0;
//         }
//         if(childy != 0){
//           for(int i = 0; i < 2; i++){
//             if(childy[i] != 0) {
//               if(childy[i]->childx != 0){
//                 childy[i]->childx[0] = 0;
//                 childy[i]->childx[1] = 0;
//               }

//               delete childy[i];
//               childy[i] = 0;
//             }
//           }
//           delete[] childy;
//           childy = 0;
//         }


//         if(child!= 0){
//           for(int i = 0; i < 4; i++){
//             if(child[i] !=0)delete child[i];
//             child[i] = 0;
//           }
//           delete[] child;
//           child = 0;
//         }


//         break;

//       case 2:
//         if(childx != 0){
//           for(int i = 0; i < 2; i++){
//             if(childx[i] != 0){
//               delete childx[i];
//               childx[i] = 0;
//             }
//           }
//           delete[] childx;
//           childx = 0;
//         }
//         break;
//       case 1:
//         if(childy != 0){
//           for(int i = 0; i < 2; i++){
//             if(childy[i] != 0){
//               delete childy[i];
//               childy[i] = 0;
//             }
//           }
//           delete[] childy;
//           childy = 0;
//         }
//         break;
//       default:
//         break;
//       }

//       if(edge != 0){
//         delete [] edge;
//         edge = 0;
//       }
//     }
//   }
  /// Destructor, it works this way without memory leakage
  ~QuadFace() {
    switch(code) {
    case 3:
      if(childx != 0) {
        for(int i=0; i<2; i++) {
          if(childx[i] != 0) {
            delete childx[i] ;
            childx[i] = 0 ;
          }
        }
        delete[] childx ;
        childx = 0 ;
      }
      if(childy != 0) {
        for(int i=0; i<2; i++) {
          if(childy[i] != 0) {
            delete[] childy[i]->childx ;
            childy[i]->childx = 0 ;
            delete childy[i] ;
            childy[i] = 0 ;
          }
        }
        delete[] childy ;
        childy = 0 ;
      }

      if(child!= 0) {
        delete[] child ;
        child = 0 ;
      }
      break ;
    case 2:
      if(childx != 0) {
        for(int i=0; i<2; i++) {
          if(childx[i] != 0) {
            delete childx[i] ;
            childx[i] = 0 ;
          }
        }
        delete[] childx ;
        childx = 0 ;
      }
      break ;
    case 1:
      if(childy != 0) {
        for(int i=0; i<2; i++) {
          if(childy[i] != 0) {
            delete childy[i] ;
            childy[i] = 0 ;
          }
        }
        delete[] childy ;
        childy = 0 ;
      }
      break ;
    default:
      break ;
    }
    if(edge != 0) {
      delete [] edge ;
      edge = 0 ;
    }
  }

  /**
   * Estimates the unsigned quadrilateral area by triangulating about
   * simple_center().
   *
   * Returns half the norm of the summed triangle cross products.
   */
  double area() {
    Node* c = simple_center() ;
    vect3d tmp_center = c->p ;
    vect3d sum = vect3d(0.0, 0.0, 0.0) ;
    for(int i=0; i<2; i++) {
      sum += cross((edge[i]->tail->p - tmp_center), (edge[i]->head->p - tmp_center)) ;
    }
    for(int i=2; i<4; i++){
      sum += cross((edge[i]->head->p - tmp_center), (edge[i]->tail->p - tmp_center)) ;
    }
    if(c!=0) { delete c ; }
    return 0.5*norm(sum) ;
  }

  /**
   * Computes a length-weighted center of the boundary edge midpoints.
   *
   * Requires all boundary edges to have been split because it reads each
   * midpoint from edge[i]->child[0]->tail. Returns a newly allocated Node;
   * the caller owns it.
   */
  Node* wireframe() {

    // allocate edgecenter
    std::vector<vect3d> edgecenter(4) ;
    std::vector<double> len(4) ;

    // get edge centers
    for(int i=0; i<4; i++) {
      edgecenter[i] = edge[i]->child[0]->tail->p ;
      len[i] = edge[i]->length() ;
    }

    // calculate the weighted center of the edge centers
    vect3d p = weighted_center(edgecenter, len) ;
    return new Node(p) ;
  }

  /**
   * Computes the face-center node selected by CENTROID.
   *
   * If CENTROID selects wireframe(), all boundary edges must already be split.
   * Returns a newly allocated Node; the caller owns it.
   */
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

  /**
   * Computes the unweighted center of this face's four corner nodes.
   *
   * Returns a newly allocated Node; the caller owns it.
   */
  Node* simple_center() {
    std::vector<vect3d> nodes(4) ;
    // get nodes
    for(int i=0; i<2; i++) {
      nodes[i] = edge[i]->head->p ;
    }
    for(int i=2; i<4; i++) {
      nodes[i] = edge[i]->tail->p ;
    }

    // calculate the weighted center of nodes
    vect3d p = point_center(nodes) ;
    return new Node(p) ;
  }

  /**
   * Returns the midpoint node for each boundary edge.
   *
   * Requires all boundary edges to have been split. The caller provides an
   * array of at least four Node* entries; this function fills it with borrowed
   * midpoint-node pointers.
   */
  void getEdgeCenter(Node** edgecenter) const {
    for(int i=0; i<4; i++) {
      edgecenter[i] = edge[i]->child[0]->tail ;
    }
  }

  /**
   * Returns the midpoint node for edge[edgeID].
   *
   * Requires edge[edgeID] to have been split. The returned pointer is borrowed.
   */
  Node* getEdgeCenter(int edgeID) const {
    return edge[edgeID]->child[0]->tail ;
  }

  /**
   * Returns the QuadFace corner node numbered nodeID.
   *
   * Valid nodeID values are 0 through 3, matching the corner ordering in the
   * class diagram. The returned pointer is borrowed.
   */
  Node* getNode(int nodeID) const {
    if(nodeID==0 || nodeID == 1) {
      return edge[nodeID]->head ;
    }else {
      return edge[nodeID]->tail ;
    }
  }

  /// Only if code is 3
  Node* getCenter() const {
    return child[0]->edge[1]->tail ; //unsafe version
  }

  /**
   * Appends the terminal unsplit subfaces below this QuadFace to leaves.
   *
   * Existing contents of leaves are preserved. Any non-leaf node must already
   * have the child array for its split code populated.
   */
  void get_leaves(std::vector<QuadFace*>& leaves) ;

  /// This has been built as defined in cell, and has been resplit with
  /// orientCode. Get leaves that is in the same order as the face is bulit as
  /// defined by face2node and resplit without orientCode.
  void get_leaves(const std::vector<char>& facePlan, char orientCode,
                  std::vector<QuadFace*>& fine_faces) ;

  /// Define face2node
  void set_f2n(std::list<int32>& f2n) ;

  int get_num_leaves( )const ;

  /// Used in building cells, quadface is built as defined in cell, and split with orientCode.
  /// All new nodes and edges are put into node_list and edge_list.
  void split(char splitCode, char orientCode, std::list<Node*>& node_list,
             std::list<Edge*>& edge_list) ;

  /// Only used in transfer_plan_q2g
  void empty_split(char splitCode);

  /// Used in building cells, quadface is built as defined in cell, and split with orientCode.
  /// All new nodes and edges are put into node_list and edge_list.
  void resplit(const std::vector<char>& facePlan, char orientCode,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list);

  // Uused in building cells, quadface is built as defined in cell, and split with orientCode.
  // All new nodes and edges are put into node_list and edge_list.
  void resplit(const std::vector<char>& facePlan, char orientCode,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::vector<QuadFace*>& fine_faces);

  /// Only used in transfer_plan_q2g
  void empty_resplit(const std::vector<char>& facePlan, char orientCode,
                     std::vector<QuadFace*>& fine_faces) ;

public:
  Edge** edge ;
  QuadFace** child ;
  QuadFace** childx ;
  QuadFace** childy ;

  /// Split code, can change value during splitting.
  /// code= 1: split only in y direction, only childy is defined. childx = child = 0
  /// code= 2: split only in x direction, 2 childx, childy= child = 0
  /// code= 3: split in both x and y direction, 2 childx, 2 childy, 4 child
  char code ;
};

QuadFace* build_quad_face(const Entity* face2node,
                          const Entity* face2edge,
                          const const_MapVec<2>& edge2node,
                          const const_store<vect3d>& pos,
                          const const_store<std::vector<char> >& edgePlan,
                          std::list<Node*>& bnode_list,
                          std::list<Edge*>& edge_list) ;

/// Parallel version
QuadFace* build_quad_face(const Entity* face2node,
                          const Entity* face2edge,
                          const const_MapVec<2>& edge2node,
                          const const_store<vect3d>& pos,
                          const const_store<std::vector<char> >& edgePlan,
                          const const_store<int>& node_offset,
                          const const_store<int>& node_l2f,
                          std::list<Node*>& bnode_list,
                          std::list<Edge*>& edge_list) ;

/// This function is used in build_general_cell with quadface
QuadFace* build_tmp_quad_face(const Entity* face2node,
                              const Entity* face2edge,
                              const const_MapVec<2>& edge2node,
                              const const_store<std::vector<char> >& edgePlan,
                              std::list<Node*>& bnode_list,
                              std::list<Edge*>& edge_list) ;

/// If the intersection of the leaves of f1 and the leaves of f2 is empty
bool is_overlapped(QuadFace* f1, QuadFace* f2);

/// Return the intersection of leaves of f1 and leaves of f2
std::vector<QuadFace*> overlap(QuadFace* f1, QuadFace* f2) ;

/// For serial version, write out .vog file
void write_quad_inner_faces(const std::map<QuadFace*, NeibIndex>& faces,
                             int cell_offset, int& mxppf, ofstream& ofile) ;

/// When a quadface is resplit according to faceplan and its node need to be
/// tagged according to faceplan1, assume the node is stored in bnode_list
/// that start at bnode_begin++ until the end of list, build 2 temp quadface,
/// resplit them according to faceplan and faceplan1, find the node
/// correspondence and tag the node in bnode_list.
void tag_quad_face(const Entity* face2node,
                   const Entity* face2edge,
                   const const_MapVec<2>& edge2node,
                   const const_store<std::vector<char> >& edgePlan,
                   const std::vector<char>& facePlan, char orientCode,
                   const std::vector<char>& nodeTag, //the tag for facePlan
                   const std::vector<char>& facePlan1,
                   std::list<Node*>& bnode_list, //node list from facePlan1
                   std::list<Node*>::const_iterator bnode_begin) ; //the ++bnode_begin is the start point

/// Compile the facePlan according the tree structure of aQuadFace
/// std::vector<char> make_faceplan( QuadFace* aFace);
void cleanup_list(std::list<Node*>& node_list,
                  std::list<Edge*>& edge_list,
                  std::list<QuadFace*>& face_list) {

  for(std::list<Node*>::iterator p = node_list.begin(); p != node_list.end(); p++) {
    if((*p) != 0) {
      delete (*p) ;
      (*p) = 0 ;
    }
  }
  node_list.clear() ;

  for(std::list<Edge*>::iterator p = edge_list.begin(); p != edge_list.end(); p++) {
    if((*p) != 0) {
      delete (*p) ;
      (*p) = 0 ;
    }
  }
  edge_list.clear() ;

  for(std::list<QuadFace*>::iterator p = face_list.begin();  p != face_list.end(); p++) {
    if((*p) != 0) {
      delete (*p) ;
      (*p) = 0 ;
    }
  }
  face_list.clear() ;
}

void cleanup_list(std::list<QuadFace*>& face_list) {
  for(std::list<QuadFace*>::iterator p = face_list.begin();  p != face_list.end(); p++) {
    if((*p) != 0) {
      delete (*p) ;
      (*p) = 0 ;
    }
  }
  face_list.clear() ;
}

void extract_quad_edge(const std::vector<char>&, std::vector<char>&, unsigned int) ;
std::vector<char> merge_quad_face(std::vector<char>& facePlanL, char orientCodeL) ;
std::vector<char> merge_quad_face(std::vector<char>& facePlanL, char orientCodeL,
                                  std::vector<char>& facePlanR, char orientCodeR) ;
#endif
