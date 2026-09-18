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

#ifndef FACE_H
#define FACE_H

#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include "node_edge.h"

/**
 * @file face.h
 *
 * Polygonal Face objects and isotropic face refinement.
 */


/**
 * Polygonal face stored as an ordered array of Edge pointers. needReverse
 * records each edge's direction relative to the face node order. A split
 * creates one quadrilateral child per boundary edge.
 *
 * A Face deletes its edge and needReverse arrays and its child faces. The
 * Edge objects are shared and are deleted separately through the builders'
 * cleanup lists.
 */
class Face{
public:
  /**
   * Allocate arrays for n boundary edges and their direction flags. The
   * entries are uninitialized; callers must fill them before using the face
   * geometry.
   */
  Face(int n):numEdge(n),edge(new Edge*[n]),needReverse(new bool[n]),child(0){}

  /**
   * Use the supplied edge and direction arrays. The Face takes ownership of
   * arrays e and r, but does not delete the Edge objects referenced by e.
   */
  Face(int n, Edge** e, bool* r):numEdge(n),edge(e), needReverse(r),child(0){}

  /**
   * Create an empty Face without boundary arrays. Set numEdge before building
   * its child tree; geometric operations also require the boundary edges.
   */
  Face():numEdge(0),edge(0),needReverse(0),child(0){}

  /**
   * Delete child faces and the child, edge, and needReverse arrays. The Edge
   * objects are deleted separately.
   */
  ~Face(){
    if(child!= 0) {
      for(int i=0; i<numEdge; i++) {
        if(child[i] != 0) {
          delete child[i] ;
          child[i] = 0 ;
        }
      }
      delete [] child ;
      child = 0 ;
    }

    if(edge!=0) {
      delete [] edge ;
      edge = 0 ;
    }

    if(needReverse !=0) {
      delete [] needReverse ;
      needReverse = 0 ;
    }
  }


  /**
  * Computes an unweighted center from the stored edge-head vertices.
  *
  * @return Newly allocated Node at the arithmetic mean of edge[i]->head
  * coordinates; the caller owns the returned Node.
  */
  Node* simple_center() {
    std::vector<vect3d> nodes(numEdge) ;
    for(int i=0; i<numEdge; i++) {
      nodes[i] = edge[i]->head->p ;
    }
    // calculate the mass center of the edge centers
    vect3d p = point_center(nodes) ;
    return new Node(p) ;
  }


  /**
   * Estimates the unsigned polygon area by triangulating about simple_center().
   *
   * @return Half the norm of the summed triangle cross products.
   */
  double area() {
    Node* c = simple_center() ;
    vect3d tmp_center = c->p ;
    vect3d sum = vect3d(0.0, 0.0, 0.0) ;
    for(int i=0; i<numEdge; i++) {
      sum += cross((edge[i]->tail->p - tmp_center), (edge[i]->head->p - tmp_center)) ;
    }
    if(c!=0) { delete c ; }
    return 0.5*norm(sum) ;
  }


  /**
   * Computes the face-center node used by split().
   *
   * The selected formula is controlled by the global CENTROID setting. With the
   * current wireframe setting, all boundary edges must already have midpoint
   * children because wireframe() reads those edge-center nodes.
   *
   * @return Newly allocated center Node; the caller owns the returned Node.
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
   * Computes a length-weighted center of the boundary edge midpoints.
   *
   * @pre Each boundary edge has already been split, so edge[i]->child[0]->tail
   * is the midpoint node for that edge. The total boundary-edge length must be
   * nonzero because weighted_center() divides by the weight sum.
   *
   * @return Newly allocated Node at the weighted center; the caller owns the
   * returned Node.
   */
  inline Node* wireframe() {

    // allocate edgecenter
    std::vector<vect3d> edgecenter(numEdge) ;
    std::vector<double> len(numEdge) ;

    // get edge centers
    for(int i=0; i<numEdge; i++) {
      edgecenter[i] = edge[i]->child[0]->tail->p ;
      len[i] = edge[i]->length() ;
    }

    // calculate the mass center of the edge centers
    vect3d p = weighted_center(edgecenter, len) ;
    return new Node(p) ;
  }


  /**
   * Returns the refinement level recorded on the first boundary edge.
   *
   * @pre The face has at least one boundary edge.
   */
  int getLevel() { return edge[0]->level ; } ;


  /**
   * Fill the caller's array of at least numEdge pointers with existing edge
   * midpoints. All boundary edges must already be split; no nodes are
   * allocated or transferred.
   */
  void getEdgeCenter(Node** edgecenter) const {
    for(int i=0; i<numEdge; i++) {
      edgecenter[i] = edge[i]->child[0]->tail ;
    }
  }


  /**
   * Return the child index for theFace, or -1 if it is not an immediate
   * child.
   */
  int containFace(Face* theFace) const {
    if(child !=0) {
      for(int i=0; i<numEdge; i++) {
        if(child[i] == theFace) {
          return i ;
        }
      }
    }
    return -1 ;
  }


  /// Return the boundary-edge index for theEdge, or -1 if it is not found.
  int containEdge(Edge* theEdge) const {
    for(int i=0; i<numEdge; i++) {
      if(edge[i] == theEdge) {
        return i ;
      }
    }
    return -1 ;
  }


  /**
   * Return the vertex index for theNode in face node order, using needReverse
   * to select edge endpoints. Compare node pointers; return -1 if none
   * matches.
   */
  int containNode(const Node* theNode) const {
    std::vector<Node*> f2n(numEdge) ;
    for(int i=0; i<numEdge; i++) {
      if(needReverse[i]) f2n[i] = edge[i]->tail ;
      else f2n[i] = edge[i]->head ;
    }
    int nodeID = -1 ;
    for(int i=0; i<numEdge; i++) {
      if(f2n[i] == theNode) {
        nodeID = i ;
        break ;
      }
    }
    return nodeID ;
  }


  /// Append the leaf faces in depth-first child order to leaves.
  void get_leaves(std::vector<Face*>& leaves) ;


  /**
   * Writes the face-to-node ordering for the current refined face boundary.
   *
   * Each boundary edge contributes leaf-edge endpoint indices. Reversed edges
   * are traversed from tail to head to preserve the face-local orientation.
   *
   * @param f2n Output list replaced with node indices in face order.
   */
  void set_f2n(std::list<int32>& f2n) ;


  /**
   * Split each unsplit boundary edge, create the face-center node and edges
   * to the midpoints, and create one quadrilateral child per boundary edge.
   * Append new nodes and root edges to the supplied lists. An already-split
   * face is unchanged.
   */
  void split(std::list<Node*>& node_list, std::list<Edge*>& edge_list) ;


  /**
   * Split this face using the prism's orientCode to order boundary-edge
   * splits. The child layout is the same as split() without orientation, but
   * new nodes are appended in face2node order. An already-split face is
   * unchanged.
   */
  void split(char orientCode, std::list<Node*>& node_list, std::list<Edge*>& edge_list) ;


  /**
   * Create child Face objects without nodes or edges. Used when only the
   * refinement tree is needed.
   */
  void empty_split() ;


  /**
   * Apply facePlan in breadth-first order and append the selected faces to
   * fine_face. Code 1 splits a face and queues its children; code 0 or an
   * omitted entry selects the current face without removing existing
   * children. The tree must be consistent with the plan.
   *
   * Append new nodes and root edges to the supplied lists. See @ref
   * fvmadapt_plans_and_balancing for plan encoding.
   */
  void resplit(const std::vector<char>& facePlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::vector<Face*>& fine_face) ;


  /**
   * Apply the splits in facePlan in breadth-first order, without collecting
   * faces. Append new nodes and root edges to the supplied lists. Existing
   * splits are not removed.
   */
  void resplit(const std::vector<char>& facePlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list) ;


  /**
   * Apply facePlan to a prism face stored in cell order. Use orientCode to
   * follow the plan in face2node order, appending new nodes and root edges to
   * the supplied lists. Existing splits are not removed.
   */
  void resplit(const std::vector<char>& facePlan,
               char orientCode,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list) ;


  /**
   * Apply facePlan without creating nodes or edges and return the number of
   * faces selected by zero or omitted entries. Existing children are not
   * removed, so the tree must be consistent with the plan.
   */
  int empty_resplit(const std::vector<char>& facePlan) ;


  /**
   * Apply facePlan without creating nodes or edges, appending faces selected
   * by zero or omitted entries to leaves in plan order. Existing children are
   * not removed, so the tree must be consistent with the plan.
   */
  void empty_resplit(const std::vector<char>& facePlan, std::vector<Face*>& leaves) ;


  /**
   * Apply an extracted prism face plan to the child tree in face2node order.
   * orientCode maps from the prism face order. Code 1 splits the current
   * face; code 8 queues the same face again without splitting it. No nodes or
   * edges are created.
   */
  void empty_resplit(const std::vector<char>& facePlan, char orientCode) ;


  /**
   * Encodes the current face tree as a breadth-first refinement plan.
   *
   * Split nodes are written as code `1`; leaves are written as code `0`.
   * Trailing no-split entries are removed from the returned vector.
   *
   * @return Face refinement plan for the current tree shape.
   */
  std::vector<char> make_faceplan() ;


  /// Return the number of leaf faces in this tree.
  int get_num_leaves() const ; //for mxfpc
public:
  int numEdge ;
  Edge** edge ;

  // if each edge is built as edge2node, and the face is built as face2node,
  // needReverse is true if face2node[i] == edge2node[face2edge[i]][1]

  // if each edge is defined as in prism, and the face is also defined as in
  // prism, needReverse is also false.

  // during split, needReverse is decided by both the way edge is defined and
  // the way the face is defined.
  bool* needReverse ;
  Face** child ;
};


/**
 * Build a Face in face2node order, with node positions from pos and boundary
 * edges split according to edgePlan. Append allocated nodes and root edges to
 * bnode_list and edge_list. The caller deletes the returned Face.
 */
Face* build_general_face( const Entity* face2node, int num_edge,
                          const Entity* face2edge,
                          const const_MapVec<2>& edge2node,
                          const const_store<vect3d>& pos,
                          const const_store<std::vector<char> >& edgePlan,
                          std::list<Node*>& bnode_list,
                          std::list<Edge*>& edge_list) ;


/**
 * Build a Face with split boundary edges. Assign original node indices from
 * node_l2f and edge-interior indices from node_offset. Append allocated nodes
 * and root edges to bnode_list and edge_list. The caller deletes the returned
 * Face.
 */
Face* build_general_face( const Entity* face2node, int num_edge,
                          const Entity* face2edge,
                          const const_MapVec<2>& edge2node,
                          const const_store<vect3d>& pos,
                          const const_store<int>& node_offset,
                          const const_store<std::vector<char> >& edgePlan,
                          std::list<Node*>& bnode_list,
                          std::list<Edge*>& edge_list,
                          const const_store<int>& node_l2f) ;


/**
 * Build a four-edge Face on an integer reference square, applying edgePlan to
 * its boundary. num_edge must be 4. Used to match nodes with a temporary
 * QuadFace. Append allocated nodes and root edges to the supplied lists; the
 * caller deletes the returned Face.
 */
Face* build_tmp_general_face( const Entity* face2node, int num_edge,
                              const Entity* face2edge,
                              const const_MapVec<2>& edge2node,
                              const const_store<std::vector<char> >& edgePlan,
                              std::list<Node*>& bnode_list,
                              std::list<Edge*>& edge_list) ;


/// Return true if the two face trees share a leaf face pointer.
bool is_overlapped(Face* f1, Face* f2) ;


/**
 * Map a child index from the prism face order to face2node order. numEdge is
 * 3 or 4; orientCode uses the prism face convention.
 */
int general_childID_orient_c2f(int childID_c, char orientCode, int numEdge) ;


/**
 * Map a child index from face2node order to the prism face order. numEdge is
 * 3 or 4; orientCode uses the prism face convention.
 */
int general_childID_orient_f2c(int childID_f, char orientCode, int numEdge) ;


/**
 * Map an edge index from face2node order to the prism face order. numEdge is
 * 3 or 4, and i is in [0, numEdge).
 */
int general_edgeID_orient_f2c(int i, char orientCode, int numEdge) ;


/**
 * Delete the listed nodes, root edges, and root faces, then clear the lists.
 * Each object must appear only once. Edge and Face destructors delete their
 * children, so those children must not also appear in the lists.
 */
inline void cleanup_list(std::list<Node*>& node_list,
                         std::list<Edge*>& edge_list,
                         std::list<Face*>& face_list) {
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

  for(std::list<Face*>::iterator p = face_list.begin();  p != face_list.end(); p++) {
    if((*p) != 0) {
      delete (*p) ;
      (*p) = 0 ;
    }
  }
  face_list.clear() ;
}


/**
 * Delete the listed root faces and clear the list. Each face must appear only
 * once; child faces are deleted by their parent.
 */
inline void cleanup_list( std::list<Face*>& face_list) {
  for(std::list<Face*>::iterator p = face_list.begin();  p != face_list.end(); p++) {
    if((*p) != 0) {
      delete (*p) ;
      (*p) = 0 ;
    }
  }
  face_list.clear() ;
}

#endif
