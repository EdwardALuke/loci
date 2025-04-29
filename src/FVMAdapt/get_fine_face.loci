//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
///////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
// This file generate  fine_faces according refinement plans. The
// original grid is general element. Isotropic refinement is used. the rules should be
// used after node_offset and cell_offset has been set up.
// FineFaces is defined as c1, c2,followed by face2node.
//
////////////////////////////////////////////////////////////////////////////////////////                             


#include <iostream>
#include <vector>
#include <Loci.h>
#include <Tools/tools.h>
#include "diamondcell.h"

using std::cerr;
using std::endl;
using std::vector;
//get  fine_faces of general cells
class get_general_cell_faces : public pointwise_rule{
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos;
  const_store<int> node_offset;
  const_store<int> cell_offset;
  const_store<bool> is_quadface;
  store<Loci::FineFaces> fine_faces;

  const_store<int>  node_l2f;
public:
  get_general_cell_faces(){
    name_store("balancedCellPlan", cellPlan);
    name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("balanced_node_offset", node_offset);
    name_store("balanced_cell_offset", cell_offset);
    name_store("fine_faces", fine_faces);
    name_store("fileNumber(pos)", node_l2f);
    name_store("is_quadface", is_quadface);
    
    input("balancedCellPlan,balanced_node_offset,balanced_cell_offset");
    input("(lower, upper, boundary_map)->(balancedFacePlan,is_quadface, balanced_node_offset)");
    input("(lower, upper, boundary_map)->face2edge->(balancedEdgePlan,balanced_node_offset)");
    input("(lower, upper, boundary_map)->face2node->(pos, fileNumber(pos))");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
  
    output("fine_faces");
    constraint("gnrlcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
      
      do_loop(seq, this);
    }
  }
  void calculate(Entity cc){
    
    
    
    if(cellPlan[cc].size() == 0){
      vector<vector<int> >().swap(fine_faces[cc]);
      return;
    }
    
    std::list<Node*> bnode_list; //boundary node
    std::list<Edge*> edge_list;
    std::list<Face*> face_list;
    std::list<Node*> node_list; //inner node
    int c1, c2;

    //build a cell
    Cell* aCell =  build_general_cell(lower[cc].begin(), lower.num_elems(cc),
                                      upper[cc].begin(), upper.num_elems(cc),
                                      boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                      is_quadface,
                                      face2node,
                                      face2edge,
                                      edge2node,
                                      pos,
                                      edgePlan,
                                      facePlan,
                                      node_offset,
                                      bnode_list,
                                      edge_list,
                                      face_list,
                                      node_l2f);
      
          
          
          
   
    std::vector<DiamondCell*> cells;
    //split the cell
    aCell->resplit( cellPlan[cc], 
                    node_list,
                    edge_list,
                    face_list,
                    cells);
      
    // index the nodes
    
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++, nindex++){
    
      (*np)->index = node_offset[cc] + nindex;
    }
    //set up inner faces
    std::vector<std::vector<Edge*> > n2e = aCell->set_n2e();
    std::list<pair<Face*, NeibIndex> > inner_faces;
    set_general_faces(aCell, cells, n2e, inner_faces);
    
    //put inner faces into fine_faces[cc]
    std::vector<std::vector<int> >(inner_faces.size()).swap(fine_faces[cc]);
    
    int fptr = 0;                    
    for(std::list<pair<Face*, NeibIndex> >::const_iterator pf = inner_faces.begin();
        pf != inner_faces.end(); pf++,fptr++){

      // face2node is store in faceVertex 
      std::list<int32> faceVertex;
      //for each face, set up face2node 
      pf->first->set_f2n(faceVertex);
      //transfer c1 , c2 from local(inside a original cell) numbering  to file numbering
      c2 = (pf->second).c2+cell_offset[cc];
      c1 = (pf->second).c1+cell_offset[cc];
      
      //put c1, c2 and faceVertex into one vector
      int vec_size = faceVertex.size() +2;
      std::vector<int>(vec_size).swap(fine_faces[cc][fptr]);
      fine_faces[cc][fptr][0] = c1;
      fine_faces[cc][fptr][1] = c2;

      std::list<int32>::const_iterator nptr = faceVertex.begin();
      for(int i = 2; i < vec_size; i++,nptr++){
        fine_faces[cc][fptr][i] = *nptr;
      }
    } 
    
    //clean up
    if(aCell != 0){
      delete aCell;
      aCell = 0;
    }
    
    //aCell will clean up these
    cleanup_list(bnode_list, edge_list, face_list);
    cleanup_list(node_list);
  }
};
register_rule<get_general_cell_faces> register_get_general_cell_faces;  


//get fine_faces of interior faces  
class get_interior_general_face_faces : public pointwise_rule{
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<std::vector<char> > cellPlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos; 
  const_Map cl;
  const_Map cr;
  const_store<int> node_offset;
  const_store<int> cell_offset;
   
  store<Loci::FineFaces> fine_faces;

  const_store<int> node_l2f;
public:
  get_interior_general_face_faces(){
    name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("balancedCellPlan", cellPlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("cl", cl);
    name_store("cr", cr);
    name_store("balanced_node_offset", node_offset);
    name_store("balanced_cell_offset", cell_offset);
    name_store("fileNumber(pos)", node_l2f);
  
    name_store("fine_faces", fine_faces);
    input("balancedFacePlan, balanced_node_offset");
    input("face2edge->(balancedEdgePlan, balanced_node_offset)");
    input("(cl, cr)->(balancedCellPlan, balanced_cell_offset)");
    input("(cl, cr)->(upper, lower, boundary_map)->face2node->(pos, fileNumber(pos))");
   
    input("(cl, cr)->(upper, lower, boundary_map)-> face2edge->edge2node->pos");
    input("face2node->(pos, fileNumber(pos))");
   
    output("fine_faces");
    constraint("(cl, cr)->gnrlcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
    
      do_loop(seq, this);
    }
   
  }
  void calculate(Entity f){
   
   
    
    std::list<Node*> bnode_list, node_list;
    std::list<Edge*> edge_list;
    std::list<int32> faceVertex;
    int c1, c2;
    Face* aFace = build_general_face(face2node[f].begin(), face2node.num_elems(f),
                                     face2edge[f].begin(),
                                     edge2node,
                                     pos,
                                     node_offset,
                                     edgePlan,
                                     bnode_list,
                                     edge_list,
                                     node_l2f);
                                          

    if(facePlan[f].size() == 0){
    
      //num of fine faces is 1. 
      aFace->set_f2n(faceVertex);
      vector<vector<int> >(1).swap(fine_faces[f]); 
      int vec_size = faceVertex.size()+2;
      vector<int>(vec_size).swap(fine_faces[f][0]);

      //get c1, c2(local numbering start with 1)
      c1 = cell_offset[cl[f]] +1;
      c2 = cell_offset[cr[f]]+1;
      
      fine_faces[f][0][0] = c1;
      fine_faces[f][0][1] = c2;
      
      std::list<int32>::const_iterator nptr = faceVertex.begin();
      for(int i = 2; i < vec_size; i++, nptr++){
        fine_faces[f][0][i] = *nptr;
      }
      //clean up
      if(aFace != 0) delete aFace;
      cleanup_list(bnode_list, edge_list);  
      return;
    }

    //split the face
    std::vector<Face*> fine_face;
    aFace->resplit(facePlan[f],
                   node_list,
                   edge_list,
                   fine_face);
    
    // index the node_list
   
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++,nindex++){
   
      (*np)->index = node_offset[f] + nindex;
    }
   
  
    //get CL and CR
    std::vector<int32> CL = get_c1(lower[cl[f]].begin(), lower.num_elems(cl[f]),
                                   upper[cl[f]].begin(), upper.num_elems(cl[f]),
                                   boundary_map[cl[f]].begin(), boundary_map.num_elems(cl[f]),
                                   face2node,
                                   face2edge,
                                   edge2node,
                                   cellPlan[cl[f]],
                                   facePlan[f],
                                   f,
                                   node_l2f);
  
    if(fine_face.size() != CL.size()){
      cerr<<"WARNING: CL has the incorrect size" << endl;
      Loci::Abort();
    }
    
    std::vector<int32>    CR = get_c1(lower[cr[f]].begin(), lower.num_elems(cr[f]),
                                      upper[cr[f]].begin(), upper.num_elems(cr[f]),
                                      boundary_map[cr[f]].begin(), boundary_map.num_elems(cr[f]),
                                      face2node,
                                      face2edge,
                                      edge2node,
                                      cellPlan[cr[f]],
                                      facePlan[f],
                                      f,
                                      node_l2f);
    
    
    if(fine_face.size() != CR.size()){
      cerr<<"WARNING: CR has the incorrect size" << endl;
      Loci::Abort();
    }

    //set the size of fine_faces[f]
    vector<std::vector<int> >(fine_face.size()).swap(fine_faces[f]);
    
    //set fine_faces[f]
    for( unsigned int  fptr = 0; fptr < fine_face.size(); fptr++){
      
      //set c1
      if(cellPlan[cl[f]].size() == 0){
        c1= cell_offset[cl[f]] + 1;
      }
      else{
        c1= cell_offset[cl[f]] + CL[fptr];              
      }
      //set c2
      if(cellPlan[cr[f]].size() == 0){
        c2 = cell_offset[cr[f]] + 1;
      }
      else{
        c2 = cell_offset[cr[f]] + CR[fptr]; 
      }
      
      //set faceVertex
      faceVertex.clear();
      fine_face[fptr]->set_f2n(faceVertex);
      
      int vec_size = faceVertex.size() +2;
      vector<int>(vec_size).swap(fine_faces[f][fptr]);
      fine_faces[f][fptr][0] = c1;
      fine_faces[f][fptr][1] = c2;
      std::list<int32>::const_iterator nptr = faceVertex.begin();
      
      for(int i = 2; i < vec_size; i++,nptr++){
        fine_faces[f][fptr][i] = *nptr;
      }
          
    } 
    
    //clean up
    if(aFace != 0){
      delete aFace;
      aFace = 0;
    }
    cleanup_list(node_list);
    cleanup_list(bnode_list, edge_list);
  }
};
register_rule<get_interior_general_face_faces> register_get_interior_general_face_faces;

//get inner_nodes and fine_faces of boundary faces
class get_boundary_general_face_faces : public pointwise_rule{
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<std::vector<char> > cellPlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos; 
  const_Map cl;
  const_store<int> node_offset;
  const_store<int> cell_offset;
  const_store<std::string> boundary_names;
  const_Map ref;
  store<Loci::FineFaces> fine_faces;

  const_store<int> node_l2f;
public:
  get_boundary_general_face_faces(){
    name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("balancedCellPlan", cellPlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("cl", cl);
    name_store("balanced_node_offset", node_offset);
    name_store("balanced_cell_offset", cell_offset);
    name_store("boundary_tags", boundary_names);
    name_store("ref", ref);
  
    name_store("fine_faces", fine_faces);
    name_store("fileNumber(pos)", node_l2f);
    input("balancedFacePlan, balanced_node_offset");
    input("face2edge->(balancedEdgePlan, balanced_node_offset)");
    input("cl->(balancedCellPlan, balanced_cell_offset)");
    input("cl->(upper, lower, boundary_map)->face2node->(pos, fileNumber(pos))");
    input("cl->(upper, lower, boundary_map)-> face2edge->edge2node->pos");
   
    input("face2node->(pos, fileNumber(pos))");
    input("ref->boundary_tags");
    
    output("fine_faces");
    constraint("cl->gnrlcells, boundary_faces");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
     
      do_loop(seq, this);
    }
  }
  void calculate(Entity f){
   
   
    std::list<Node*> bnode_list, node_list;
    std::list<Edge*> edge_list;
    std::list<int32> faceVertex;
    int c1;

    int c2 = -atoi(&(boundary_names[ref[f]][3]));
    //build a face
    Face* aFace = build_general_face(face2node[f].begin(), face2node.num_elems(f),
                                     face2edge[f].begin(),
                                     edge2node,
                                     pos,
                                     node_offset,
                                     edgePlan,
                                     bnode_list,
                                     edge_list,
                                     node_l2f);
    
    
    if(facePlan[f].size() == 0){
      //num of inner nodes is 0
     
      
     
      //num of fine faces is 1. 
      vector< std::vector<int> >(1).swap(fine_faces[f]);

      aFace->set_f2n(faceVertex);
      int vec_size = faceVertex.size()+2;
      std::vector<int>(vec_size).swap(fine_faces[f][0]);

      //get c1
      c1 = cell_offset[cl[f]] +1;
     
      fine_faces[f][0][0] = c1;
      fine_faces[f][0][1] = c2;
      
      std::list<int32>::const_iterator nptr = faceVertex.begin();
      for(int i = 2; i < vec_size; i++,nptr++){
        fine_faces[f][0][i] = *nptr;
      }
      

      //clean up
      if(aFace != 0){
        delete aFace;
        aFace = 0;
      }
      cleanup_list(bnode_list, edge_list);  
      return;
    }

    //split the face
    std::vector<Face*> fine_face;
    aFace->resplit(facePlan[f],
                   node_list,
                   edge_list,
                   fine_face);


    // index the node_list
   
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++,nindex++){
   
      (*np)->index = node_offset[f] + nindex;
    }
   
        
   
    //get CL    
    std::vector<int32> CL = get_c1(lower[cl[f]].begin(), lower.num_elems(cl[f]),
                                   upper[cl[f]].begin(), upper.num_elems(cl[f]),
                                   boundary_map[cl[f]].begin(),
                                   boundary_map.num_elems(cl[f]),
                                   face2node,
                                   face2edge,
                                   edge2node,
                                   cellPlan[cl[f]],
                                   facePlan[f],
                                   f,
                                   node_l2f);
    
    if(fine_face.size() != CL.size()){
      cerr<<"WARNING: CL has the incorrect size" << endl;
      Loci::Abort();
    }
    
    
    vector<vector<int> >(fine_face.size()).swap(fine_faces[f]);                    
    //set face_vertices
    for( unsigned int  fptr = 0; fptr < fine_face.size(); fptr++){
      faceVertex.clear();
      fine_face[fptr]->set_f2n(faceVertex);
      
      //set c1
      if(cellPlan[cl[f]].size() == 0){
        c1= cell_offset[cl[f]] + 1;
      }
      else{
        c1= cell_offset[cl[f]] + CL[fptr];              
      }
       
    
      int vec_size = faceVertex.size() +2;
      vector<int>(vec_size).swap(fine_faces[f][fptr]);
      fine_faces[f][fptr][0] = c1;
      fine_faces[f][fptr][1] = c2;

      std::list<int32>::const_iterator nptr = faceVertex.begin();
      for(int i = 2; i < vec_size; i++, nptr++){
        fine_faces[f][fptr][i] = *nptr;
      }
    }
    
  
    //clean up
    if(aFace != 0){
      delete aFace;
      aFace = 0;
    }
    cleanup_list(node_list);
    cleanup_list(bnode_list, edge_list);
  }
};
register_rule<get_boundary_general_face_faces> register_get_boundary_general_face_faces;



