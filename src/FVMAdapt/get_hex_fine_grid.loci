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
//                                get_hex_fine_grid.cc                                    //
//                                by: Qiuhan Xue                                      //
//                                                                                    //
// This file generate inner_nodes  according refinement plans. The
//original grid is hexcell element. Anisotropic refinement is used. the rules should be
//used After node_offset and cell_offset has been set up.
// FineFaces is defined as c1, c2,followed by face2node.
//
////////////////////////////////////////////////////////////////////////////////////////                             


#include <iostream>
#include <vector>
#include <Loci.h>
#include <Tools/tools.h>
#include "hexcell.h"

using std::cerr;
using std::endl;
using std::vector;
using Loci::storeRepP;
//get inner_nodes  of hex cells
class get_hex_cell_nodes : public pointwise_rule{
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
  

  store<Loci::FineNodes> inner_nodes;
 

  const_store<int> node_l2f;

public:
  get_hex_cell_nodes(){
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
  

    name_store("inner_nodes", inner_nodes);
    name_store("fileNumber(face2node)", node_l2f);

    
    input("balancedCellPlan");
    input("(hex2face, hex2node, hexOrientCode)");
    input("(lower, upper, boundary_map)->(fileNumber(face2node),balancedFacePlan)");
    input("(lower, upper, boundary_map)->face2edge->(balancedEdgePlan)");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    
    output("inner_nodes");
    constraint("hexcells");
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
    std::list<QuadFace*> face_list;
    std::list<Node*> node_list; //inner node
    

   //build a Cell
      HexCell* aCell =  build_hex_cell(lower[cc].begin(), lower.num_elems(cc),
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
      
      aCell->resplit(cellPlan[cc], 
                     node_list,
                     edge_list,
                     face_list,
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
    cleanup_list(bnode_list, edge_list, face_list);
    cleanup_list(node_list);
  }
};
register_rule<get_hex_cell_nodes> register_get_hex_cell_nodes;  


