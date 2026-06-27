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
 * @ingroup fvmadapt_elements
 * @brief Polygonal face tree used for general faces and prism end faces.
 *
 *  this file include the declaration of class Face. Class Face is the abstraction
 * of general face(polygon). It's defined as a collection of edges and the directions
 *  of edges. it supports only isotropic refinement.
 * If the face is built as face2node, it can resplit without orientCode(as used in
 * DiamondCell and  general Cell), A face can also be built as defined in cell, and
 * resplit with orientCode(as used in Prism)
 *
 */


/**
 * @brief Isotropic polygon-face refinement tree.
 * @ingroup fvmadapt_elements
 *
 * A `Face` is stored as an ordered set of `Edge*` values plus a `needReverse`
 * array that records how each edge direction relates to face-local node order.
 * Splitting creates one quadrilateral child per parent edge by connecting edge
 * midpoints through a face-center node.
 *
 * @warning The face owns the `edge` pointer array and any child faces, but it
 * does not own the pointed-to edges. Shared edge ownership is coordinated by
 * the builder cleanup lists.
 */
class Face{
public:
  /**
   * Constructs a face with storage for @p n boundary edges.
   *
   * The edge and needReverse arrays are allocated but not initialized; callers
   * fill them before using the face geometry. The child pointer is initialized
   * to null so split/replay code can test whether the face is a leaf.
   *
   * @param n Number of boundary edges in the polygon.
   */
  Face(int n):numEdge(n),edge(new Edge*[n]),needReverse(new bool[n]),child(0){}

  /**
   * Constructs a face around existing boundary-edge arrays.
   *
   * Ownership of the @p e and @p r arrays is transferred to the Face. The Face
   * destructor deletes those arrays, but it does not delete the Edge objects
   * referenced by @p e.
   *
   * @param n Number of boundary edges in the polygon.
   * @param e Boundary-edge pointer array.
   * @param r Per-edge orientation flags.
   */
  Face(int n, Edge** e, bool* r):numEdge(n),edge(e), needReverse(r),child(0){}

  /**
   * Constructs an empty face placeholder.
   *
   * Empty faces are used when replaying tree shape without constructing the
   * geometric edge/node data. Callers that use the placeholder as a child set
   * numEdge and any needed arrays afterward. Geometry helpers are not valid
   * until that topology has been populated.
   */
  Face():numEdge(0),edge(0),needReverse(0),child(0){}

  /**
   * Deletes child face subtrees and owned pointer arrays.
   *
   * The destructor owns the child array, the edge pointer array, and the
   * needReverse array. It does not own or delete the Edge objects referenced by
   * edge; those are tracked by builder cleanup lists.
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
   * Returns the midpoint node for each boundary edge.
   *
   * @pre All boundary edges have already been split.
   * @param edgecenter Caller-provided array of at least numEdge Node*
   * entries. The function fills it with borrowed midpoint-node pointers.
   */
  void getEdgeCenter(Node** edgecenter) const {
    for(int i=0; i<numEdge; i++) {
      edgecenter[i] = edge[i]->child[0]->tail ;
    }
  }


  /**
   * Finds an immediate child face by pointer identity.
   *
   * @param theFace Face pointer to search for.
   * @return Child index in [0, numEdge), or -1 if no immediate child matches.
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


  /**
   * Finds a boundary edge by pointer identity.
   *
   * @param theEdge Edge pointer to search for.
   * @return Edge index in [0, numEdge), or -1 if the edge is not on this face.
   */
  int containEdge(Edge* theEdge) const {
    for(int i=0; i<numEdge; i++) {
      if(edge[i] == theEdge) {
        return i ;
      }
    }
    return -1 ;
  }


  /**
   * Finds a face vertex by pointer identity in face-to-node order.
   *
   * The needReverse flags are used to choose the endpoint that represents each
   * face-local vertex.
   *
   * @param theNode Node pointer to search for.
   * @return Vertex index in [0, numEdge), or -1 if the node is not a face
   * vertex.
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


  /**
   * Collects leaf faces in tree traversal order.
   *
   * Existing contents of @p leaves are preserved; this function only appends
   * leaves from the receiver.
   *
   * @param leaves Output vector receiving leaf faces.
   */
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
   * Splits a general face once.
   *
   * The method first ensures every boundary edge has a midpoint split, then
   * creates a face-center node, radial edges from the center to each edge
   * midpoint, and one quadrilateral child face for each parent edge.
   *
   * @param node_list Receives newly created midpoint and face-center nodes.
   * @param edge_list Receives newly created edge objects.
   */
  void split(std::list<Node*>& node_list, std::list<Edge*>& edge_list) ;


  /**
   * Splits a general face once using a prism face orientation code.
   *
   * Boundary edges are split in the orientation-adjusted order returned by
   * general_edgeID_orient_f2c(). The created child topology is otherwise the
   * same as split(std::list<Node*>&, std::list<Edge*>&).
   *
   * @param orientCode Orientation code associated with the face in the prism.
   * @param node_list  Receives newly created midpoint and face-center nodes.
   * @param edge_list  Receives newly created edge objects.
   */
  void split(char orientCode, std::list<Node*>& node_list, std::list<Edge*>& edge_list) ;


  /**
   * Creates empty child Face objects without creating nodes or edges.
   *
   * This is used when replaying or converting plan vectors where only the tree
   * shape is needed.
   */
  void empty_split() ;


  /**
   * Replays a face refinement plan and returns the resulting leaf faces.
   *
   * Code `1` splits the current face and queues its children. Code `0`, or a
   * missing trailing code after the plan vector is exhausted, leaves the current
   * face as a leaf. See @ref fvmadapt_refinement_plans for the shared plan
   * encoding convention.
   *
   * @param facePlan Breadth-first general-face refinement plan.
   * @param node_list Receives nodes created during splitting.
   * @param edge_list Receives edges created during splitting.
   * @param fine_face Receives leaf faces in traversal order.
   */
  void resplit(const std::vector<char>& facePlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list,
               std::vector<Face*>& fine_face) ;


  /**
   * Replays a face refinement plan without collecting the leaf faces.
   *
   * This overload creates the same child geometry as resplit(facePlan,
   * node_list, edge_list, fine_face), but leaves callers to traverse the tree
   * later if they need the final leaves.
   *
   * @param facePlan Breadth-first general-face refinement plan.
   * @param node_list Receives nodes created during splitting.
   * @param edge_list Receives edges created during splitting.
   */
  void resplit(const std::vector<char>& facePlan,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list) ;


  /**
   * Replays a face refinement plan using a prism-local orientation code.
   *
   * This overload is used for triangular prism faces whose stored edge order
   * follows the prism cell while the refinement plan is interpreted in
   * face-to-node order.
   *
   * @param facePlan Breadth-first general-face refinement plan.
   * @param orientCode Orientation code associated with the face in the prism.
   * @param node_list Receives nodes created during splitting.
   * @param edge_list Receives edges created during splitting.
   */
  void resplit(const std::vector<char>& facePlan,
               char orientCode,
               std::list<Node*>& node_list,
               std::list<Edge*>& edge_list) ;


  /**
   * Replays a face plan into empty child faces and returns the leaf count.
   *
   * This builds tree shape only; it does not create geometric nodes or edges.
   *
   * @param facePlan Breadth-first general-face refinement plan.
   * @return Number of leaf faces represented by the plan.
   */
  int empty_resplit(const std::vector<char>& facePlan) ;


  /**
   * Replays a face plan into empty child faces and returns the leaves.
   *
   * Existing contents of @p leaves are preserved; this function only appends
   * leaves from the replayed tree.
   *
   * @param facePlan Breadth-first general-face refinement plan.
   * @param leaves Receives leaf faces in traversal order.
   */
  void empty_resplit(const std::vector<char>& facePlan, std::vector<Face*>& leaves) ;


  /**
   * Replays a face plan into empty child faces with prism orientation.
   *
   * This merge helper uses general_childID_orient_c2f() when queueing children.
   * Code `8` is treated as a propagation case: the current face is queued
   * again instead of being split or finalized.
   *
   * @param facePlan Breadth-first general-face refinement plan.
   * @param orientCode Orientation code associated with the face in the prism.
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


  /**
   * Counts leaf faces under this face tree.
   *
   * @return Number of leaves reachable from this face.
   */
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
 * Builds a polygon face from Loci face maps and edge plans.
 *
 * The face boundary nodes are created from @p pos. Each boundary edge is
 * created from @p edge2node and replayed with the corresponding entry in
 * @p edgePlan. Created nodes and edges are appended to the caller-owned cleanup
 * lists.
 *
 * @param face2node Face-to-node connectivity in face-local order.
 * @param num_edge Number of edges in the polygonal face.
 * @param face2edge Face-to-edge connectivity in face-local order.
 * @param edge2node Global edge-to-node connectivity.
 * @param pos Node coordinates.
 * @param edgePlan Refinement plans for boundary edges.
 * @param bnode_list Receives boundary and edge-interior nodes.
 * @param edge_list Receives created edge objects.
 * @return Newly allocated Face tree root.
 */
Face* build_general_face( const Entity* face2node, int num_edge,
                          const Entity* face2edge,
                          const const_MapVec<2>& edge2node,
                          const const_store<vect3d>& pos,
                          const const_store<std::vector<char> >& edgePlan,
                          std::list<Node*>& bnode_list,
                          std::list<Edge*>& edge_list) ;


/**
 * Builds a polygon face and assigns parallel node indices.
 *
 * This overload is used when boundary nodes need local-to-file indices while
 * edge-interior nodes are assigned from @p node_offset for each boundary edge.
 *
 * @param face2node Face-to-node connectivity in face-local order.
 * @param num_edge Number of edges in the polygonal face.
 * @param face2edge Face-to-edge connectivity in face-local order.
 * @param edge2node Global edge-to-node connectivity.
 * @param pos Node coordinates.
 * @param node_offset First index assigned to each edge's interior nodes.
 * @param edgePlan Refinement plans for boundary edges.
 * @param bnode_list Receives boundary and edge-interior nodes.
 * @param edge_list Receives created edge objects.
 * @param node_l2f Local-to-file node index map for boundary nodes.
 * @return Newly allocated Face tree root.
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
 * Builds a temporary polygon face for general-cell/quad-face interaction.
 *
 * Temporary nodes are placed on a reference square rather than using physical
 * coordinates. Boundary edges are replayed from @p edgePlan so the resulting
 * tree can be compared or merged with quad-face refinement state.
 *
 * @param face2node Face-to-node connectivity in face-local order.
 * @param num_edge Number of edges in the polygonal face.
 * @param face2edge Face-to-edge connectivity in face-local order.
 * @param edge2node Global edge-to-node connectivity.
 * @param edgePlan Refinement plans for boundary edges.
 * @param bnode_list Receives temporary boundary and edge-interior nodes.
 * @param edge_list Receives created edge objects.
 * @return Newly allocated temporary Face tree root.
 */
Face* build_tmp_general_face( const Entity* face2node, int num_edge,
                              const Entity* face2edge,
                              const const_MapVec<2>& edge2node,
                              const const_store<std::vector<char> >& edgePlan,
                              std::list<Node*>& bnode_list,
                              std::list<Edge*>& edge_list) ;


/**
 * Checks whether two face trees share any leaf Face pointer.
 *
 * @param f1 First face tree.
 * @param f2 Second face tree.
 * @return true when the two trees have at least one identical leaf pointer.
 */
bool is_overlapped(Face* f1, Face* f2) ;


/**
 * Maps a child index from cell-local order to face-local order.
 *
 * This conversion is used by prism face handling when child faces need to be
 * interpreted in the face-to-node ordering.
 *
 * @param childID_c Child index in cell-local ordering.
 * @param orientCode Orientation code associated with the face in the cell.
 * @param numEdge Number of parent face edges.
 * @return Child index in face-local ordering.
 */
int general_childID_orient_c2f(int childID_c, char orientCode, int numEdge) ;


/**
 * Maps a child index from face-local order to cell-local order.
 *
 * This is the opposite conversion from general_childID_orient_c2f() for the
 * same orientation-code convention.
 *
 * @param childID_f Child index in face-local ordering.
 * @param orientCode Orientation code associated with the face in the cell.
 * @param numEdge Number of parent face edges.
 * @return Child index in cell-local ordering.
 */
int general_childID_orient_f2c(int childID_f, char orientCode, int numEdge) ;


/**
 * Maps a face-local edge index to the cell-local edge index for an orientation.
 *
 * This helper is used when a general face is split through a prism-local
 * orientation. The implementation has explicit cases for triangular and
 * quadrilateral faces.
 *
 * @param i Face-local edge index.
 * @param orientCode Orientation code associated with the face in the cell.
 * @param numEdge Number of parent face edges.
 * @return Edge index in the oriented cell-local ordering.
 */
int general_edgeID_orient_f2c(int i, char orientCode, int numEdge) ;


/**
 * Deletes owned temporary nodes, edges, and faces, then clears the lists.
 *
 * These lists are used by builder and replay routines to track allocations
 * that are not owned directly by a single Face. Pass only lists that contain
 * unique owning pointers; child subtrees are handled by the destructors of the
 * listed Edge and Face roots.
 *
 * @param node_list Owning list of Node pointers to delete.
 * @param edge_list Owning list of Edge pointers to delete.
 * @param face_list Owning list of Face pointers to delete.
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
 * Deletes owned temporary faces, then clears the list.
 *
 * @param face_list Owning list of Face pointers to delete.
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
