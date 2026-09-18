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
 *
 * QuadFace directional refinement and face-orientation helpers.
 */

// f2c orient functions are used when a quadface is built as in cell and the
// facePlan is for the face defined by face2node
/**
 * Convert a QuadFace split code from face2node order to the face order in the
 * cell. orientCode is in [0, 8). Orientations that exchange x and y exchange
 * split codes 1 and 2; codes 0 and 3 are unchanged.
 */
char orient_splitCode_f2c(char splitCode, char orientCode) ;
char orient_childID_f2c(char childID, char orientCode, char splitCode) ;
char orient_edgeID_f2c(char edgeID, char orientCode) ;

/// Map a QuadFace edge index from the face order in the cell to face2node
/// order. edgeID is in [0, 4) and orientCode is in [0, 8).
char orient_edgeID_c2f(char edgeID, char orientCode);

/// For each range in leaves, return the index paired with the first faceMap
/// range that contains it. Every leaf range must have a match. get_c1_hex()
/// and get_c1_prism() use these ranges to find local fine-cell indices without
/// creating split geometry.
std::vector<int32> contain_2d(const std::vector<pair<Range2d, int32> >& faceMap,
                              const std::vector<Range2d>& leaves);

/**
 * Quadrilateral face with directional refinement, used by HexCell and
 * Prism. Corner indices and stored edge directions are:
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
 * For split code 3, child[] holds the four quarters. childx[] holds the
 * left and right halves; childy[] holds the lower and upper halves:
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
 * Plans use code 0 for no split, 1 for a y split, 2 for an x split, and
 * 3 for both directions. With code 3, childx and childy share the four
 * child faces shown above.
 */
class QuadFace{
public:
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
  /// Delete child faces and pointer arrays, accounting for children shared by
  /// childx and childy. Edge and Node objects are deleted separately.
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
   * Fill the caller's array of at least four pointers with existing edge
   * midpoints. All boundary edges must already be split; no nodes are
   * allocated or transferred.
   */
  void getEdgeCenter(Node** edgecenter) const {
    for(int i=0; i<4; i++) {
      edgecenter[i] = edge[i]->child[0]->tail ;
    }
  }

  /**
   * Return the existing midpoint node of edge[edgeID]. The edge must already
   * be split.
   */
  Node* getEdgeCenter(int edgeID) const {
    return edge[edgeID]->child[0]->tail ;
  }

  /**
   * Return the existing corner-node pointer for nodeID in [0, 4), using the
   * numbering in the class diagram.
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

  /// Append the leaf faces in depth-first child order to leaves.
  void get_leaves(std::vector<QuadFace*>& leaves) ;

  /// Append the faces selected by facePlan in the breadth-first order of the
  /// plan. The tree must already be split consistently with facePlan and
  /// orientCode.
  void get_leaves(const std::vector<char>& facePlan, char orientCode,
                  std::vector<QuadFace*>& fine_faces) ;

  /// Fill f2n with the boundary node indices in face order.
  void set_f2n(std::list<int32>& f2n) ;

  int get_num_leaves( )const ;

  /// Split this QuadFace according to splitCode in face2node order, using
  /// orientCode to map to its stored cell order. Append new nodes and root
  /// edges to the supplied lists.
  ///
  /// If the face is already split in one direction, a split in the other
  /// direction adds the missing children to form a four-way split.
  void split(char splitCode, char orientCode, std::list<Node*>& node_list,
             std::list<Edge*>& edge_list) ;

  /// Create child QuadFace objects for splitCode without nodes or edges. An
  /// already-split face is unchanged.
  void empty_split(char splitCode);

  /// Apply facePlan in breadth-first order, using orientCode to map face2node
  /// order to the face order in the cell. Append new nodes and root edges to
  /// the supplied lists. Existing splits are not removed.
  void resplit(const std::vector<char>& facePlan, char orientCode,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list);

  // Uused in building cells, quadface is built as defined in cell, and split with orientCode.
  // All new nodes and edges are put into node_list and edge_list.
  void resplit(const std::vector<char>& facePlan, char orientCode,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::vector<QuadFace*>& fine_faces);

  /// Apply facePlan to a QuadFace tree without creating nodes or edges, and
  /// append the selected faces in plan order. orientCode maps face2node order
  /// to the stored cell order.
  ///
  /// The tree must be consistent with the plan. A nonzero plan code is checked
  /// against an existing nonzero code; a zero or omitted entry appends the
  /// current face without checking or removing its children.
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

/// Build a QuadFace with split boundary edges and assign node indices from
/// node_l2f and node_offset. Append allocated nodes and root edges to the
/// supplied lists; the caller deletes the returned QuadFace.
QuadFace* build_quad_face(const Entity* face2node,
                          const Entity* face2edge,
                          const const_MapVec<2>& edge2node,
                          const const_store<vect3d>& pos,
                          const const_store<std::vector<char> >& edgePlan,
                          const const_store<int>& node_offset,
                          const const_store<int>& node_l2f,
                          std::list<Node*>& bnode_list,
                          std::list<Edge*>& edge_list) ;

/// Build a QuadFace on an integer reference square with edges split according
/// to edgePlan. Used to match nodes between face plans. Append allocated nodes
/// and root edges to the supplied lists; the caller deletes the returned
/// QuadFace.
QuadFace* build_tmp_quad_face(const Entity* face2node,
                              const Entity* face2edge,
                              const const_MapVec<2>& edge2node,
                              const const_store<std::vector<char> >& edgePlan,
                              std::list<Node*>& bnode_list,
                              std::list<Edge*>& edge_list) ;

/// Return true if the two face trees share a leaf face pointer.
bool is_overlapped(QuadFace* f1, QuadFace* f2);

/// Return the intersection of leaves of f1 and leaves of f2
std::vector<QuadFace*> overlap(QuadFace* f1, QuadFace* f2) ;

/// For serial version, write out .vog file
void write_quad_inner_faces(const std::map<QuadFace*, NeibIndex>& faces,
                             int cell_offset, int& mxppf, ofstream& ofile) ;

/// Transfer nodeTag from nodes created by facePlan to matching nodes created
/// by facePlan1, using temporary faces with integer coordinates to find
/// matches.
///
/// The destination nodes must occupy bnode_list from the entry after
/// bnode_begin to the end, in the order created by facePlan1. Nodes without a
/// match retain their tags.
void tag_quad_face(const Entity* face2node,
                   const Entity* face2edge,
                   const const_MapVec<2>& edge2node,
                   const const_store<std::vector<char> >& edgePlan,
                   const std::vector<char>& facePlan, char orientCode,
                   const std::vector<char>& nodeTag, //the tag for facePlan
                   const std::vector<char>& facePlan1,
                   std::list<Node*>& bnode_list, //node list from facePlan1
                   std::list<Node*>::const_iterator bnode_begin) ; //the ++bnode_begin is the start point

/// Delete the listed nodes, root edges, and root QuadFaces, then clear the
/// lists. Each object must appear only once. Child edges and faces are deleted
/// by their parents and must not also appear in the lists.
inline void cleanup_list(std::list<Node*>& node_list,
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

inline void cleanup_list(std::list<QuadFace*>& face_list) {
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
