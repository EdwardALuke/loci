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
#include <queue>
#include <vector>
#include <map>
#include <set>
#include <Loci.h>
#include <Tools/tools.h>
#include "hexcell.h"

using std::queue;
using std::cerr;
using std::cout;
using std::endl;
using std::set;
using std::map;
using Loci::storeRepP;
class reset_hexcell_num_nodes : public pointwise_rule{
  const_store<std::vector<char> > cellPlan;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,6> > orientCode;
  const_store<Array<char,6> > hex2face;
  const_store<Array<char,8> > hex2node;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos;
  store<int> num_inner_nodes;

  const_store<int> node_l2f;
public:
  reset_hexcell_num_nodes(){
    name_store("balancedCellPlan", cellPlan);
    name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("hexOrientCode", orientCode);
    name_store("hex2face", hex2face);
    name_store("hex2node", hex2node);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("fileNumber(face2node)", node_l2f);
    name_store("balanced_num_inner_nodes", num_inner_nodes);
 

    input("balancedCellPlan");
    input("(hex2face, hex2node, hexOrientCode)");
    input("(lower, upper, boundary_map)->balancedFacePlan");
    input("(lower, upper, boundary_map)->face2edge->balancedEdgePlan");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node");
    input("(lower, upper, boundary_map)->(fileNumber(face2node),face2edge)");
    
    output("balanced_num_inner_nodes");
 
    constraint("hexcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
       
      do_loop(seq, this);
    }
  
  }
  void calculate(Entity cc){
       
    if(cellPlan[cc].size() == 0){
      num_inner_nodes[cc] = 0;
 
      return;
    }

    std::list<Node*> bnode_list; //boundary node
    std::list<Node*> node_list; //inner node
    std::list<Edge*> edge_list;
    std::list<QuadFace*> face_list;
    
    HexCell* aCell = build_hex_cell(lower[cc].begin(), lower.num_elems(cc),
                                    upper[cc].begin(), upper.num_elems(cc),
                                    boundary_map[cc].begin(), boundary_map.num_elems(cc),
                                    hex2face[cc],
                                    hex2node[cc],
                                    orientCode[cc],
                                    face2node,
                                    face2edge,
                                    edge2node,
                                    pos,
                                    edgePlan,
                                    facePlan,
                                    bnode_list,
                                    edge_list,
                                    face_list,
                                    node_l2f);
    
    std::vector<HexCell*> cells;
    aCell->resplit( cellPlan[cc],
                    node_list,
                    edge_list,
                    face_list,
                    cells);
  
  
  
    num_inner_nodes[cc] = node_list.size();

    
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
register_rule<reset_hexcell_num_nodes> register_reset_hexcell_num_nodes;  

class reset_hexcell_num_cells : public pointwise_rule{
  const_store<std::vector<char> > cellPlan;
  store<int> num_fine_cells;
public:
  reset_hexcell_num_cells(){
    name_store("balancedCellPlan", cellPlan);
    name_store("balanced_num_fine_cells", num_fine_cells);
    input("balancedCellPlan");
    output("balanced_num_fine_cells");
    constraint("hexcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){
      do_loop(seq, this);
    }
  }
  void calculate(Entity cc){
    HexCell* aCell = new HexCell;
    num_fine_cells[cc] = aCell->num_fine_cells(cellPlan[cc]);
    delete aCell;    
  }
};
register_rule<reset_hexcell_num_cells> register_reset_hexcell_num_cells;  





class reset_quadface_nums : public pointwise_rule{
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_multiMap face2node;
  const_multiMap face2edge;
  const_MapVec<2> edge2node;
  const_store<vect3d> pos;
  store<int> num_inner_nodes;
  
public:
  reset_quadface_nums(){
    name_store("balancedFacePlan", facePlan);
    name_store("balancedEdgePlan", edgePlan);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("priority::quadface::balanced_num_inner_nodes", num_inner_nodes);
    input("balancedFacePlan");
    input("face2edge->(balancedEdgePlan, edge2node)");
    input("face2node->pos");
    output("priority::quadface::balanced_num_inner_nodes");
    constraint("quadfaces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity f){
    if(facePlan[f].size() == 0){
      num_inner_nodes[f] = 0;
      return;
    }
    std::list<Node*> bnode_list;
    std::list<Node*> node_list;
    std::list<Edge*> edge_list;
   
    
    QuadFace* aFace = build_quad_face(face2node[f].begin(),
                                      face2edge[f].begin(),
                                      edge2node,
                                      pos,
                                      edgePlan,
                                      bnode_list,
                                      edge_list);
           
    aFace->resplit(facePlan[f],
                   char(0),
                   node_list,
                   edge_list);
    
    
    num_inner_nodes[f] = node_list.size();
    
    //clean up
    if(aFace != 0){
      delete aFace;
      aFace = 0;
    }
    
    cleanup_list(bnode_list, edge_list);
    cleanup_list(node_list);
    
  }
};
register_rule<reset_quadface_nums> register_reset_quadface_nums;



