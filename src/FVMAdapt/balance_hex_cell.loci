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
///////////////////////////////////////////////////////////////////////////////
//
// This file balance hex cellPlan between cells  
//
///////////////////////////////////////////////////////////////////////////////


#include <queue>
#include <vector>
#include <utility>
#include "hexcell.h"
#include <Loci.h>
#include <algorithm>
using std::queue;

using std::cerr;
using std::endl;

using Loci::storeRepP;


class init_hex_cell_updated:public pointwise_rule{
  const_store<std::vector<char> > newCellPlan;
  store<bool> cellUnchanged;
  store<std::vector<char> > tmpCellPlan;
  
public:
  init_hex_cell_updated(){
    name_store("newCellPlan", newCellPlan);
    name_store("tmpCellPlan{n=0}", tmpCellPlan);
    name_store("cellUnchanged{n=0}", cellUnchanged);
    output("(cellUnchanged{n=0},tmpCellPlan{n=0})");
    input("newCellPlan");
    constraint("hexcells");
  }
  virtual void compute(const sequence &seq){
    
    do_loop(seq, this);
  }
  void calculate(Entity cc){
    
    cellUnchanged[cc] = (newCellPlan[cc].size() == 0);
   
    tmpCellPlan[cc] = newCellPlan[cc];
    reduce_vector(tmpCellPlan[cc]);
    
    
  }
};

register_rule<init_hex_cell_updated> register_init_hex_cell_updated;



//this rule balance each cell with all its interior faces(lower, upper)
class advance_hex_cell_updated : public pointwise_rule{
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
   const_param<int> split_mode_par;
   const_store<Array<char, 6> > hex2face;
  const_store<Array<char, 8> > hex2node;
  const_store<Array<char, 6> > orientCode;
  const_multiMap face2node;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_store<vect3d> pos;
  const_store<bool> cellUnchangedn;
  const_store<std::vector<char> > tmpCellPlann;
  const_store<std::vector<char> > facePlan;
  const_store<std::vector<char> > edgePlan;
  const_store<bool> isIndivisible;
  store<bool> cellUnchangedn1;
  store<std::vector<char> > tmpCellPlann1;
  // const_blackbox<storeRepP>  node_remap;
  //Map node_l2f;
  const_store<int> node_l2f;
public:
  advance_hex_cell_updated(){
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("boundary_map", boundary_map);
     name_store("split_mode_par", split_mode_par);
    name_store("hex2face", hex2face);
    name_store("hex2node", hex2node);
    name_store("hexOrientCode", orientCode);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("pos", pos);
    name_store("cellUnchanged{n}", cellUnchangedn);
    name_store("cellUnchanged{n+1}", cellUnchangedn1);
    name_store("tmpCellPlan{n}", tmpCellPlann);
    name_store("tmpCellPlan{n+1}", tmpCellPlann1);
    name_store("tmpFacePlan{n}", facePlan);
    name_store("tmpEdgePlan{n}", edgePlan);
    name_store("isIndivisible", isIndivisible);
    // name_store("node_remap", node_remap);
    name_store("fileNumber(face2node)", node_l2f);
    input("split_mode_par");
    input("cellUnchanged{n},tmpCellPlan{n}, hex2face, hex2node, hexOrientCode, isIndivisible");
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge-> tmpEdgePlan{n}");
    input("(lower, upper, boundary_map)->face2edge-> edge2node->pos");
    input("(lower, upper, boundary_map)->tmpFacePlan{n}" );
    // input("node_remap");
    input("(lower, upper, boundary_map)->fileNumber(face2node)");
    output("cellUnchanged{n+1}, tmpCellPlan{n+1}");
    constraint("hexcells");
  }
  virtual void compute(const sequence &seq){
    if(seq.size()!=0){

      do_loop(seq, this);

    }
  }
  void calculate(Entity cc){
    //first assume the cell will not be updated
   

    cellUnchangedn1[cc] = true;

    std::list<Node*> node_list;
    std::list<Edge*> edge_list;
    std::list<QuadFace*> face_list;
    std::list<Node*> bnode_list;


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
    
    aCell->resplit( tmpCellPlann[cc], 
                    node_list,
                    edge_list,
                    face_list,
                    cells);
    
    
       

    //balance cells
    if(!isIndivisible[cc]){
      aCell->rebalance_cells(*split_mode_par, node_list, edge_list, face_list);
    }
        

    //write new cellPlan
    std::vector<char> newCellPlan (aCell->make_cellplan());

    
    cellUnchangedn1[cc] =  (newCellPlan == tmpCellPlann[cc]);

    
    tmpCellPlann1[cc] = newCellPlan;
    reduce_vector(tmpCellPlann1[cc]);
    
    //clean up
    
    if(aCell != 0){
        delete aCell;
        aCell = 0;
    }
    cleanup_list(node_list, edge_list, face_list);
    cleanup_list(bnode_list);
    
  }
  
};

register_rule<advance_hex_cell_updated> register_advance_hex_cell_updated;


