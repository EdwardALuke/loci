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

#include <iostream>
#include <vector>
#include <Loci.h>
#include <Tools/tools.h>
#include "hexcell.h"
#include "diamondcell.h"
#include "prism.h"
using std::cerr;
using std::endl;
using std::vector;




class get_interior_face_faces_gh:public pointwise_rule{

  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<std::vector<char> > cellPlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,6> > orientCode;
  const_store<Array<char,6> > hex2face;
  const_store<Array<char,8> > hex2node;
 
  const_store<char> fr;
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
  get_interior_face_faces_gh(){
   name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("balancedCellPlan", cellPlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
     name_store("hexOrientCode", orientCode);
    name_store("hex2face", hex2face);
    name_store("hex2node", hex2node);
  
    name_store("fr", fr);
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
    input("cr->(hex2node, hex2face, hexOrientCode)");
    input("fr");
    input("(cl, cr)->(upper, lower, boundary_map)->face2node->(pos, fileNumber(pos))");
    input("(cl, cr)->(upper, lower, boundary_map)-> face2edge->edge2node->pos");

    input("face2node->(pos, fileNumber(pos))");
    
    output("fine_faces");
    constraint("cl->gnrlcells, cr->hexcells");
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
      QuadFace* aFace = build_quad_face(face2node[f].begin(),
                                        face2edge[f].begin(),
                                        edge2node,
                                        pos,
                                        edgePlan,
                                        node_offset,
                                        node_l2f,
                                        bnode_list,
                                        edge_list);
      

      if(facePlan[f].size() == 0){
       
        //num of fine faces is 1. 
        aFace->set_f2n(faceVertex);
        vector<vector<int> >(1).swap(fine_faces[f]); 
        int vec_size = faceVertex.size()+2;
        vector<int>(vec_size).swap(fine_faces[f][0]);

      //get c1, c2(local numbering start with 1)
        c1 = cell_offset[cl[f]] +1;

        std::vector<int32>    CR = get_c1_hex(cellPlan[cr[f]],
                                              facePlan[f],
                                              orientCode[cr[f]][fr[f]],
                                              fr[f]);
    
    
        if( CR.size()!= 1){
          cerr<<"WARNING: CR has the incorrect size" << endl;
          Loci::Abort();
        }


        c2 = cell_offset[cr[f]]+CR[0];
        
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
      std::vector<QuadFace*> leaves;
      
      aFace->resplit(facePlan[f],
                     char(0),
                     node_list,
                     edge_list,
                     leaves);
      // and index the node_list
   
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++,nindex++){
   
      (*np)->index = node_offset[f] + nindex;
    }
   
  
    //get CL and CR
    std::vector<int32> CL = get_c1_general(lower[cl[f]].begin(), lower.num_elems(cl[f]),
                                           upper[cl[f]].begin(), upper.num_elems(cl[f]),
                                           boundary_map[cl[f]].begin(), boundary_map.num_elems(cl[f]),
                                           true,
                                           face2node,
                                           face2edge,
                                           edge2node,
                                           cellPlan[cl[f]],
                                           facePlan[f],
                                           f,
                                           node_l2f);
    if(leaves.size() != CL.size()){
      cerr<<"WARNING: CL has the incorrect size" << endl;
      Loci::Abort();
    }
    
    std::vector<int32>    CR = get_c1_hex(cellPlan[cr[f]],
                                          facePlan[f],
                                          orientCode[cr[f]][fr[f]],
                                          fr[f]);
    
    
    if(leaves.size() != CR.size()){
      cerr<<"WARNING: CR has the incorrect size" << endl;
      Loci::Abort();
    }

    //set the size of fine_faces[f]
    vector<std::vector<int> >(leaves.size()).swap(fine_faces[f]);
    
    //set fine_faces[f]
    for( unsigned int  fptr = 0; fptr < leaves.size(); fptr++){
      
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
      leaves[fptr]->set_f2n(faceVertex);
      
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

register_rule<get_interior_face_faces_gh> register_get_interior_face_faces_gh;


class get_interior_face_faces_hg:public pointwise_rule{
   const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<std::vector<char> > cellPlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,6> > orientCode;
  const_store<Array<char,6> > hex2face;
  const_store<Array<char,8> > hex2node;
 
  const_store<char> fl;
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
  get_interior_face_faces_hg(){
    name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("balancedCellPlan", cellPlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
     name_store("hexOrientCode", orientCode);
    name_store("hex2face", hex2face);
    name_store("hex2node", hex2node);
  
    name_store("fl", fl);
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
    input("cl->(hex2node, hex2face, hexOrientCode)");
    input("fl");
    input("(cl, cr)->(upper, lower, boundary_map)->face2node->(pos, fileNumber(pos))");
    input("(cl, cr)->(upper, lower, boundary_map)-> face2edge->edge2node->pos");

    input("face2node->(pos, fileNumber(pos))");
   
    output("fine_faces");
    constraint("cr->gnrlcells, cl->hexcells");
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
      QuadFace* aFace = build_quad_face(face2node[f].begin(),
                                        face2edge[f].begin(),
                                        edge2node,
                                        pos,
                                        edgePlan,
                                        node_offset,
                                        node_l2f,
                                        bnode_list,
                                        edge_list);
      

      if(facePlan[f].size() == 0){
       
        //num of fine faces is 1. 
        aFace->set_f2n(faceVertex);
        vector<vector<int> >(1).swap(fine_faces[f]); 
        int vec_size = faceVertex.size()+2;
        vector<int>(vec_size).swap(fine_faces[f][0]);

      //get c1, c2(local numbering start with 1)


        std::vector<int32> CL = get_c1_hex(cellPlan[cl[f]],
                                           facePlan[f],
                                           orientCode[cl[f]][fl[f]],
                                           fl[f]);
  
        if( CL.size()!= 1){
          cerr<<"WARNING: CL has the incorrect size" << endl;
          Loci::Abort();
        }

        c1 = cell_offset[cl[f]] +CL[0];
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
      std::vector<QuadFace*> leaves;
      
      aFace->resplit(facePlan[f],
                     char(0),
                     node_list,
                     edge_list,
                     leaves);
      // and index the node_list
   
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++,nindex++){
   
      (*np)->index = node_offset[f] + nindex;
    }
   
  
    //get CL and CR
    std::vector<int32> CL = get_c1_hex(cellPlan[cl[f]],
                                       facePlan[f],
                                       orientCode[cl[f]][fl[f]],
                                       fl[f]);
  
    if(leaves.size() != CL.size()){
      cerr<<"WARNING: CL has the incorrect size" << endl;
      Loci::Abort();
    }
    
    std::vector<int32>    CR = get_c1_general(lower[cr[f]].begin(), lower.num_elems(cr[f]),
                                              upper[cr[f]].begin(), upper.num_elems(cr[f]),
                                              boundary_map[cr[f]].begin(), boundary_map.num_elems(cr[f]),
                                              true,
                                              face2node,
                                              face2edge,
                                              edge2node,
                                              cellPlan[cr[f]],
                                              facePlan[f],
                                              f,
                                              node_l2f);
    
    if(leaves.size() != CR.size()){
      cerr<<"WARNING: CR has the incorrect size" << endl;
      Loci::Abort();
    }

    //set the size of fine_faces[f]
    vector<std::vector<int> >(leaves.size()).swap(fine_faces[f]);
    
    //set fine_faces[f]
    for( unsigned int  fptr = 0; fptr < leaves.size(); fptr++){
      
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
      leaves[fptr]->set_f2n(faceVertex);
      
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

register_rule<get_interior_face_faces_hg> register_get_interior_face_faces_hg;

class get_interior_face_faces_gp:public pointwise_rule{
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<std::vector<char> > cellPlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,5> > orientCode;
  const_store<Array<char,5> > prism2face;
  const_store<Array<char,6> > prism2node;

  const_store<char> fr;
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
  get_interior_face_faces_gp(){
     name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("balancedCellPlan", cellPlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prismOrientCode", orientCode);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    
    name_store("fr", fr);
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
    input("cr->(prism2node, prism2face, prismOrientCode)");
    input("fr");
    input("(cl, cr)->(upper, lower, boundary_map)->face2node->(pos, fileNumber(pos))");
    input("(cl, cr)->(upper, lower, boundary_map)-> face2edge->edge2node->pos");
    input("face2node->(pos, fileNumber(pos))");

    output("fine_faces");
    constraint("cl->gnrlcells, cr->prisms");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
     
      do_loop(seq, this);
    }
  }
  void calculate(Entity f){
    std::vector<QuadFace*> fine_qfaces;
    std::vector<Face*> fine_gfaces;
    QuadFace* aqFace = 0;
    Face* agFace = 0;
    unsigned int num_leaves = 0;
    
    std::list<Node*> bnode_list, node_list;
    std::list<Edge*> edge_list;
    std::list<int32> faceVertex;
    int c1, c2;
    
      if(face2node.num_elems(f) == 4){
        aqFace = build_quad_face(face2node[f].begin(),
                                 face2edge[f].begin(),
                                 edge2node,
                                 pos,
                                 edgePlan,
                                 node_offset,
                                 node_l2f,
                                 bnode_list,
                                 edge_list);
      

        if(facePlan[f].size() == 0){
          
          
          //num of fine faces is 1. 
          aqFace->set_f2n(faceVertex);
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
        if(aqFace != 0) delete aqFace;
        cleanup_list(bnode_list, edge_list);  
        return;
        }
        
        //split the face
        
        
        aqFace->resplit(facePlan[f],
                        char(0),
                        node_list,
                        edge_list,
                        fine_qfaces);
        num_leaves = fine_qfaces.size();
      }

      else{
        
        agFace = build_general_face(face2node[f].begin(),face2node.num_elems(f),
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
          agFace->set_f2n(faceVertex);
          vector<vector<int> >(1).swap(fine_faces[f]); 
          int vec_size = faceVertex.size()+2;
          vector<int>(vec_size).swap(fine_faces[f][0]);

          //get c1, c2(local numbering start with 1)
          c1 = cell_offset[cl[f]] +1;

          std::vector<int32>    CR = get_c1_prism(cellPlan[cr[f]],
                                                  facePlan[f],
                                                  orientCode[cr[f]][fr[f]],
                                                  fr[f]);
    
    
          if( CR.size()!= 1){
            cerr<<"WARNING: CR has the incorrect size" << endl;
            Loci::Abort();
          }

          c2 = cell_offset[cr[f]]+CR[0];
          
          fine_faces[f][0][0] = c1;
          fine_faces[f][0][1] = c2;
          
          std::list<int32>::const_iterator nptr = faceVertex.begin();
          for(int i = 2; i < vec_size; i++, nptr++){
            fine_faces[f][0][i] = *nptr;
          }
          //clean up
          if(agFace != 0) delete agFace;
          cleanup_list(bnode_list, edge_list);  
          return;
       }

 //split the face
       agFace->resplit(facePlan[f],
                       node_list,
                       edge_list,
                       fine_gfaces);
       num_leaves = fine_gfaces.size();
      }

      //and index the node_list
     
      int nindex = 0;
      for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++,nindex++){
     
        (*np)->index = node_offset[f] + nindex;
      }
   
  
    //get CL and CR
      std::vector<int32> CL = get_c1_general(lower[cl[f]].begin(), lower.num_elems(cl[f]),
                                           upper[cl[f]].begin(), upper.num_elems(cl[f]),
                                           boundary_map[cl[f]].begin(), boundary_map.num_elems(cl[f]),
                                           (face2node.num_elems(f) == 4),
                                           face2node,
                                           face2edge,
                                           edge2node,
                                           cellPlan[cl[f]],
                                           facePlan[f],
                                           f,
                                           node_l2f);
  
      if(num_leaves != CL.size()){
        cerr<<"WARNING: CL has the incorrect size" << endl;
        Loci::Abort();
      }
    
      std::vector<int32>    CR = get_c1_prism(cellPlan[cr[f]],
                                              facePlan[f],
                                              orientCode[cr[f]][fr[f]],
                                              fr[f]);
    
    
      if(num_leaves != CR.size()){
        cerr<<"WARNING: CR has the incorrect size" << endl;
        Loci::Abort();
      }
      
    //set the size of fine_faces[f]
      vector<std::vector<int> >(num_leaves).swap(fine_faces[f]);
    
      //set fine_faces[f]
      for( unsigned int  fptr = 0; fptr < num_leaves; fptr++){
        
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
      if(face2node.num_elems(f) == 4){
        fine_qfaces[fptr]->set_f2n(faceVertex);
      }
      else{fine_gfaces[fptr]->set_f2n(faceVertex);
      }

      
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
    if(aqFace != 0){
      delete aqFace;
      aqFace = 0;
    }

     if(agFace != 0){
      delete agFace;
      agFace = 0;
     }
     cleanup_list(node_list);
     cleanup_list(bnode_list, edge_list);
  }
};

register_rule<get_interior_face_faces_gp> register_get_interior_face_faces_gp;

class get_interior_face_faces_pg:public pointwise_rule{
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<std::vector<char> > cellPlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,5> > orientCode;
  const_store<Array<char,5> > prism2face;
  const_store<Array<char,6> > prism2node;

  const_store<char> fl;
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
  get_interior_face_faces_pg(){
      name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("balancedCellPlan", cellPlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prismOrientCode", orientCode);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    
    name_store("fl", fl);
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
    input("cl->(prism2node, prism2face, prismOrientCode)");
    input("fl");
    input("(cl, cr)->(upper, lower, boundary_map)->face2node->(pos, fileNumber(pos))");
    input("(cl, cr)->(upper, lower, boundary_map)-> face2edge->edge2node->pos");
    input("face2node->(pos, fileNumber(pos))");
  
    output("fine_faces");
    constraint("cr->gnrlcells, cl->prisms");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
     
      do_loop(seq, this);
    }
  }
  void calculate(Entity f){
   
  std::vector<QuadFace*> fine_qfaces;
    std::vector<Face*> fine_gfaces;
    QuadFace* aqFace = 0;
    Face* agFace = 0;
    unsigned int num_leaves = 0;
    
    std::list<Node*> bnode_list, node_list;
    std::list<Edge*> edge_list;
    std::list<int32> faceVertex;
    int c1, c2;
    
      if(face2node.num_elems(f) == 4){
        aqFace = build_quad_face(face2node[f].begin(),
                                 face2edge[f].begin(),
                                 edge2node,
                                 pos,
                                 edgePlan,
                                 node_offset,
                                 node_l2f,
                                 bnode_list,
                                 edge_list);
      

        if(facePlan[f].size() == 0){
                   
          //num of fine faces is 1. 
          aqFace->set_f2n(faceVertex);
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
        if(aqFace != 0) delete aqFace;
        cleanup_list(bnode_list, edge_list);  
        return;
        }
        
        //split the face
        
        
        aqFace->resplit(facePlan[f],
                        char(0),
                        node_list,
                        edge_list,
                        fine_qfaces);
        num_leaves = fine_qfaces.size();
      }

      else{
        
        agFace = build_general_face(face2node[f].begin(),face2node.num_elems(f),
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
          agFace->set_f2n(faceVertex);
          vector<vector<int> >(1).swap(fine_faces[f]); 
          int vec_size = faceVertex.size()+2;
          vector<int>(vec_size).swap(fine_faces[f][0]);

          //get c1, c2(local numbering start with 1)



          //get CL
          std::vector<int32> CL =  get_c1_prism(cellPlan[cl[f]],
                                                facePlan[f],
                                                orientCode[cl[f]][fl[f]],
                                                fl[f]);
  
  
          if( CL.size()!=1){
            cerr<<"WARNING: CL has the incorrect size" << endl;
            Loci::Abort();
          }




          c1 = cell_offset[cl[f]] +CL[0];
          c2 = cell_offset[cr[f]]+1;

        
          
          fine_faces[f][0][0] = c1;
          fine_faces[f][0][1] = c2;
          
          std::list<int32>::const_iterator nptr = faceVertex.begin();
          for(int i = 2; i < vec_size; i++, nptr++){
            fine_faces[f][0][i] = *nptr;
          }
          //clean up
          if(agFace != 0) delete agFace;
          cleanup_list(bnode_list, edge_list);  
          return;
       }

 //split the face
       agFace->resplit(facePlan[f],
                       node_list,
                       edge_list,
                       fine_gfaces);
       num_leaves = fine_gfaces.size();
      }

      // and index the node_list
      
      int nindex = 0;
      for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++,nindex++){
      
        (*np)->index = node_offset[f] + nindex;
      }
   
  
    //get CL and CR
      std::vector<int32> CL =  get_c1_prism(cellPlan[cl[f]],
                                           facePlan[f],
                                           orientCode[cl[f]][fl[f]],
                                           fl[f]);
  
  
      if(num_leaves != CL.size()){
        cerr<<"WARNING: CL has the incorrect size" << endl;
        Loci::Abort();
      }
    
      std::vector<int32>    CR = get_c1_general(lower[cr[f]].begin(), lower.num_elems(cr[f]),
                                           upper[cr[f]].begin(), upper.num_elems(cr[f]),
                                           boundary_map[cr[f]].begin(), boundary_map.num_elems(cr[f]),
                                           (face2node.num_elems(f)==4),
                                           face2node,
                                           face2edge,
                                           edge2node,
                                           cellPlan[cr[f]],
                                           facePlan[f],
                                           f,
                                           node_l2f);
      if(num_leaves != CR.size()){
        cerr<<"WARNING: CR has the incorrect size" << endl;
        Loci::Abort();
      }
      
    //set the size of fine_faces[f]
      vector<std::vector<int> >(num_leaves).swap(fine_faces[f]);
    
      //set fine_faces[f]
      for( unsigned int  fptr = 0; fptr < num_leaves; fptr++){
        
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
      if(face2node.num_elems(f) == 4){
        fine_qfaces[fptr]->set_f2n(faceVertex);
      }
      else{fine_gfaces[fptr]->set_f2n(faceVertex);
      }

      
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
    if(aqFace != 0){
      delete aqFace;
      aqFace = 0;
    }

     if(agFace != 0){
      delete agFace;
      agFace = 0;
     }
     cleanup_list(node_list);
     cleanup_list(bnode_list, edge_list);  
  }
};

register_rule<get_interior_face_faces_pg> register_get_interior_face_faces_pg;


class get_interior_face_faces_hp:public pointwise_rule{
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<std::vector<char> > cellPlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,6> > orientCodeH;
   const_store<Array<char,5> > orientCodeP;
  const_store<Array<char,6> > hex2face;
  const_store<Array<char,8> > hex2node;
   const_store<Array<char,5> > prism2face;
  const_store<Array<char,6> > prism2node;
  const_store<char> fl;
  const_store<char> fr;
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
  get_interior_face_faces_hp(){
    name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("balancedCellPlan", cellPlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("hexOrientCode", orientCodeH);
     name_store("prismOrientCode", orientCodeP);
    name_store("hex2face", hex2face);
    name_store("hex2node", hex2node);
      name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    name_store("fl", fl);
    name_store("fr", fr);
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
    input("cl->(hex2node, hex2face, hexOrientCode)");
    input("cr->(prism2node, prism2face, prismOrientCode)");
    input("(fl, fr)");
    input("(cl, cr)->(upper, lower, boundary_map)->face2node->pos");
    input("(cl, cr)->(upper, lower, boundary_map)-> face2edge->edge2node->pos");
    input("face2node->(pos, fileNumber(pos))");
   
    output("fine_faces");
    constraint("cl->hexcells, cr->prisms");
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
      QuadFace* aFace = build_quad_face(face2node[f].begin(),
                                        face2edge[f].begin(),
                                        edge2node,
                                        pos,
                                        edgePlan,
                                        node_offset,
                                        node_l2f,
                                        bnode_list,
                                        edge_list);
      

      if(facePlan[f].size() == 0){
        

        //num of fine faces is 1. 
        aFace->set_f2n(faceVertex);
        vector<vector<int> >(1).swap(fine_faces[f]); 
        int vec_size = faceVertex.size()+2;
        vector<int>(vec_size).swap(fine_faces[f][0]);

      //get c1, c2(local numbering start with 1)

        //get CL and CR
        std::vector<int32> CL = get_c1_hex(cellPlan[cl[f]],
                                           facePlan[f],
                                           orientCodeH[cl[f]][fl[f]],
                                           fl[f]);
  
        if(CL.size()!= 1){
          cerr<<"WARNING: CL has the incorrect size" << endl;
          Loci::Abort();
        }
        
      
        c1 = cell_offset[cl[f]] +CL[0];
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
      std::vector<QuadFace*> leaves;
      
      aFace->resplit(facePlan[f],
                     char(0),
                     node_list,
                     edge_list,
                     leaves);
      // and index the node_list
    
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++,nindex++){
    
      (*np)->index = node_offset[f] + nindex;
    }
   
  
    //get CL and CR
    std::vector<int32> CL = get_c1_hex(cellPlan[cl[f]],
                                           facePlan[f],
                                           orientCodeH[cl[f]][fl[f]],
                                           fl[f]);
  
    if(leaves.size() != CL.size()){
      cerr<<"WARNING: CL has the incorrect size" << endl;
      Loci::Abort();
    }
    
    std::vector<int32>    CR = get_c1_prism(cellPlan[cr[f]],
                                            facePlan[f],
                                            orientCodeP[cr[f]][fr[f]],
                                            fr[f]);
    
    
    if(leaves.size() != CR.size()){
      cerr<<"WARNING: CR has the incorrect size" << endl;
      Loci::Abort();
    }

    //set the size of fine_faces[f]
    vector<std::vector<int> >(leaves.size()).swap(fine_faces[f]);
    
    //set fine_faces[f]
    for( unsigned int  fptr = 0; fptr < leaves.size(); fptr++){
      
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
      leaves[fptr]->set_f2n(faceVertex);
      
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

register_rule<get_interior_face_faces_hp> register_get_interior_face_faces_hp;

class get_interior_face_faces_ph:public pointwise_rule{
 const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<std::vector<char> > cellPlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,6> > orientCodeH;
   const_store<Array<char,5> > orientCodeP;
  const_store<Array<char,6> > hex2face;
  const_store<Array<char,8> > hex2node;
   const_store<Array<char,5> > prism2face;
  const_store<Array<char,6> > prism2node;
  const_store<char> fl;
  const_store<char> fr;
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
  get_interior_face_faces_ph(){
     name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("balancedCellPlan", cellPlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("hexOrientCode", orientCodeH);
     name_store("prismOrientCode", orientCodeP);
    name_store("hex2face", hex2face);
    name_store("hex2node", hex2node);
      name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    name_store("fl", fl);
    name_store("fr", fr);
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
    input("cr->(hex2node, hex2face, hexOrientCode)");
    input("cl->(prism2node, prism2face, prismOrientCode)");
    input("(fl, fr)");
    input("(cl, cr)->(upper, lower, boundary_map)->face2node->pos");
    input("(cl, cr)->(upper, lower, boundary_map)-> face2edge->edge2node->pos");
    input("face2node->(pos, fileNumber(pos))");
   
    output("fine_faces");
   
    constraint("cr->hexcells, cl->prisms");
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
      QuadFace* aFace = build_quad_face(face2node[f].begin(),
                                        face2edge[f].begin(),
                                        edge2node,
                                        pos,
                                        edgePlan,
                                        node_offset,
                                        node_l2f,
                                        bnode_list,
                                        edge_list);
      

      if(facePlan[f].size() == 0){
       

        //num of fine faces is 1. 
        aFace->set_f2n(faceVertex);
        vector<vector<int> >(1).swap(fine_faces[f]); 
        int vec_size = faceVertex.size()+2;
        vector<int>(vec_size).swap(fine_faces[f][0]);

    
    
        std::vector<int32>    CR = get_c1_hex(cellPlan[cr[f]],
                                              facePlan[f],
                                              orientCodeH[cr[f]][fr[f]],
                                              fr[f]);
        
    
        if(CR.size()!= 1){
          cerr<<"WARNING: CR has the incorrect size" << endl;
          Loci::Abort();
        }




        c1 = cell_offset[cl[f]] + 1;
        c2 = cell_offset[cr[f]]+ CR[0];
        
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
      std::vector<QuadFace*> leaves;
      
      aFace->resplit(facePlan[f],
                     char(0),
                     node_list,
                     edge_list,
                     leaves);
      //and index the node_list
   
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++,nindex++){
   
      (*np)->index = node_offset[f] + nindex;
    }
   
  
    //get CL and CR
    std::vector<int32> CL = get_c1_prism(cellPlan[cl[f]],
                                           facePlan[f],
                                           orientCodeP[cl[f]][fl[f]],
                                           fl[f]);
  
    if(leaves.size() != CL.size()){
      cerr<<"WARNING: CL has the incorrect size" << endl;
      Loci::Abort();
    }
    
    std::vector<int32>    CR = get_c1_hex(cellPlan[cr[f]],
                                            facePlan[f],
                                            orientCodeH[cr[f]][fr[f]],
                                            fr[f]);
    
    
    if(leaves.size() != CR.size()){
      cerr<<"WARNING: CR has the incorrect size" << endl;
      Loci::Abort();
    }

    //set the size of fine_faces[f]
    vector<std::vector<int> >(leaves.size()).swap(fine_faces[f]);
    
    //set fine_faces[f]
    for( unsigned int  fptr = 0; fptr < leaves.size(); fptr++){
      
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
      leaves[fptr]->set_f2n(faceVertex);
      
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

register_rule<get_interior_face_faces_ph> register_get_interior_face_faces_ph;

