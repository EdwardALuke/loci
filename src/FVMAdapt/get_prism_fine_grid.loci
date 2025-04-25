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
////////////////////////////////////////////////////////////////////////////////////////
//                                get_prism_fine_grid.cc                                    //
//                                by: Qiuhan Xue                                      //
//                                                                                    //
// This file generate inner_nodes and fine_faces according refinement plans. The
//original grid is prism element. Anisotropic refinement is used. the rules should be
//used After node_offset and cell_offset has been set up.
// FineFaces is defined as c1, c2,followed by face2node.
//
////////////////////////////////////////////////////////////////////////////////////////                             


#include <iostream>
#include <vector>
#include <Loci.h>
#include <Tools/tools.h>
#include "prism.h"

using std::cerr;
using std::endl;
using std::vector;
//get inner_nodes  of prism cells
class get_prism_cell_nodes : public pointwise_rule{
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
   const_store<Array<char,5> > orientCode;
  const_store<Array<char,5> > prism2face;
  const_store<Array<char,6> > prism2node;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos;

  const_store<int> node_l2f;
  store<Loci::FineNodes> inner_nodes;
  
 
public:
  get_prism_cell_nodes(){
    name_store("balancedCellPlan", cellPlan);
    name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("prismOrientCode", orientCode);
    name_store("prism2face", prism2face);
    name_store("prism2node", prism2node);
    
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);

    
    name_store("inner_nodes", inner_nodes);
     name_store("fileNumber(face2node)", node_l2f);
  
    
    input("balancedCellPlan");
    input("(prism2face, prism2node, prismOrientCode)");
    input("(lower, upper, boundary_map)->balancedFacePlan");
    input("(lower, upper, boundary_map)->face2edge->balancedEdgePlan");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("(lower, upper, boundary_map)->fileNumber(face2node)");
    output("inner_nodes");
    constraint("prisms");
  }
  virtual void compute(const sequence &seq){
   
     if(seq.size()!=0){
      
        do_loop(seq, this);
      }
  
  
   
   
  }
  void calculate(Entity cc){
    
   
    
    if(cellPlan[cc].size() == 0){
      vector<vect3d>().swap(inner_nodes[cc]);

      return;
    }
    
    std::list<Node*> bnode_list; //boundary node
    std::list<Edge*> edge_list;
    std::list<QuadFace*> qface_list;
    std::list<Face*> gface_list;
    std::list<Node*> node_list; //inner node
   

   //build a Cell
      Prism* aCell =  build_prism_cell(lower[cc].begin(), lower.num_elems(cc),
                                       upper[cc].begin(), upper.num_elems(cc),
                                       boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                       prism2face[cc],
                                       prism2node[cc],
                                       orientCode[cc],
                                       face2node,
                                       face2edge,
                                       edge2node,
                                       pos,
                                       edgePlan,
                                       facePlan,
                                       bnode_list,
                                       edge_list,
                                       qface_list,
                                       gface_list,
                                       node_l2f);
          
      std::vector<Prism*> cells;
      
      aCell->resplit(cellPlan[cc], 
                     node_list,
                     edge_list,
                     qface_list,
                     gface_list,
                     cells);
          
 
    //put the node_list into inner_nodes, and index the nodes
    vector<vect3d>(node_list.size()).swap(inner_nodes[cc]); 
    int nindex = 0;
    for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++, nindex++){
      inner_nodes[cc][nindex]=(*np)->p;
     
    }
  
    //clean up
    if(aCell != 0){
      delete aCell;
      aCell = 0;
    }
    
    //aCell will clean up these
    cleanup_list(bnode_list, edge_list, qface_list);
    cleanup_list(gface_list);
    cleanup_list(node_list);
  }
};
register_rule<get_prism_cell_nodes> register_get_prism_cell_nodes;  


//get inner_nodes of  faces  
class get_face_nodes : public pointwise_rule{
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos; 
 
  const_store<bool> is_quadface;

  store<Loci::FineNodes> inner_nodes;
 
 
public:
  get_face_nodes(){
    name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("is_quadface", is_quadface);
    name_store("inner_nodes", inner_nodes);
    
    
    input("(balancedFacePlan,is_quadface)");
    input("face2edge->(balancedEdgePlan, edge2node)");
    input("face2node->pos");
    output("inner_nodes");
    constraint("faces");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
      do_loop(seq, this);
    }
    
  }
  void calculate(Entity f){

    if(facePlan[f].size() == 0){
      //num of inner nodes is 0
      vector<vect3d>().swap(inner_nodes[f]);
      return;
    }

    
  

    
    std::list<Node*> bnode_list, node_list;
    std::list<Edge*> edge_list;
    
   

    if(is_quadface[f]){

      QuadFace* aqFace  = build_quad_face(face2node[f].begin(),
                                          face2edge[f].begin(),
                                          edge2node,
                                          pos,
                                          edgePlan,
                                          bnode_list,
                                          edge_list);
             
      aqFace->resplit(facePlan[f],
                      char(0),
                      node_list,
                      edge_list);
      
    
      //put node_list into inner_nodes[f]
      vector<vect3d>(node_list.size()).swap(inner_nodes[f]);
      int nindex = 0;
      for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++,nindex++){
        inner_nodes[f][nindex] = (*np)->p;
      }
      
      
   
      //clean up
      if(aqFace != 0){
        delete aqFace;
        aqFace = 0;
      }
      
    }else{
      
      Face * agFace = build_general_face(face2node[f].begin(),face2node.num_elems(f),
                                         face2edge[f].begin(),
                                         edge2node,
                                         pos,
                                         edgePlan,
                                         bnode_list,
                                         edge_list);
      
      
      agFace->resplit(facePlan[f],
                      node_list,
                      edge_list);
      
    

      //put node_list into inner_nodes[f] and index the node_list
      vector<vect3d>(node_list.size()).swap(inner_nodes[f]);
      int nindex = 0;
      for(std::list<Node*>::const_iterator np = node_list.begin(); np!= node_list.end(); np++,nindex++){
        inner_nodes[f][nindex] = (*np)->p;
      }
   
  
     if(agFace != 0){
      delete agFace;
      agFace = 0;
     }
    }
    cleanup_list(node_list);
    cleanup_list(bnode_list, edge_list);
  }
};
register_rule<get_face_nodes> register_get_face_nodes;

