//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
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
#include <queue>
#include <vector>
#include <utility>
#include <list>
#include <bitset>
#include "hexcell.h"
#include <Loci.h>
#include <algorithm>
#include "globals.h"

using std::list;
using std::queue;

using std::cerr;
using std::endl;
using std::cout;



class check_folded_face:public pointwise_rule{
  const_multiMap face2node;
  const_store<vect3d> pos;
  store<bool> isFolded;
public:
  check_folded_face(){
    name_store("face2node", face2node);
    name_store("pos", pos);
    name_store("isFolded", isFolded);
    input("face2node->pos");
    output("isFolded");
    constraint("faces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity f){
    isFolded[f] = 0;
    int numNode = face2node.num_elems(f);
    if(numNode <= 3) return;
    
   
    std::vector<vect3d>  edge(numNode -1);
    std::vector<double> norms(numNode-2);


    //defines edges from node i+1 to node 0
    for(int i = 0; i < numNode-1; i++){
      edge[i] = pos[face2node[f][i+1]] - pos[face2node[f][0]];
    }

    //defines area of triangle
    for(int i = 0; i< numNode-2; i++){
      edge[i] = cross(edge[i], edge[i+1]);
    }
    
    
    double maxAngle = 0.0;
    double total_area  = 0.0;
   
   
    for(int i = 0; i< numNode-2; i++){
      norms[i] = norm(edge[i]);
      total_area += norms[i];
    }
    
    for(int i = 0; i< numNode-3; i++){
      if(norms[i] <= NORMALIZE_ZERO_THRESHOLD || norms[i+1] <= NORMALIZE_ZERO_THRESHOLD || norms[i]/total_area < 0.001)continue;
      else{
        double dot_value = dot(edge[i]/norms[i], edge[i+1]/norms[i+1]);
        // if(abs(dot_value) <1.0e-8 ) dot_value = abs(dot_value);//prevent error cause folded face
        maxAngle = max(acos(min(dot_value,1.0)), maxAngle);
      }
    }
    isFolded[f] = (maxAngle > Globals::fold);
    //if not folded, check again
    if(isFolded[f] == 0){
      //defines edges from node i+1 to node 0
      for(int i = 0; i < numNode-1; i++){
        edge[i] = pos[face2node[f][i]] - pos[face2node[f][numNode-1]];
      }

    //defines area of triangle
      for(int i = 0; i< numNode-2; i++){
        edge[i] = cross(edge[i], edge[i+1]);
      }
    
    
    double maxAngle = 0.0;
    double total_area  = 0.0;
   
   
    for(int i = 0; i< numNode-2; i++){
      norms[i] = norm(edge[i]);
      total_area += norms[i];
    }
    
    for(int i = 0; i< numNode-3; i++){
      if(norms[i] <= NORMALIZE_ZERO_THRESHOLD || norms[i+1] <= NORMALIZE_ZERO_THRESHOLD || norms[i]/total_area < 0.001)continue;
      else{
        double dot_value = dot(edge[i]/norms[i], edge[i+1]/norms[i+1]);
        // if(abs(dot_value) <1.0e-8 ) dot_value = abs(dot_value);//prevent error cause folded face
        maxAngle = max(acos(min(dot_value,1.0)), maxAngle);
      }
    }
    isFolded[f] = (maxAngle > Globals::fold);


    }
    //  if(isFolded[f]) cerr << " folded face exist " << f << endl;
   
  }
};
register_rule<check_folded_face> register_check_folded_face;

class check_folded_cell_unit:public unit_rule{
    store<bool> isIndivisible;
public:
  check_folded_cell_unit(){
    name_store("isIndivisible", isIndivisible);
    output("isIndivisible");
    constraint("geom_cells");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity cc){
    isIndivisible[cc] = false;
  }
};

register_rule<check_folded_cell_unit> register_check_folded_cell_unit;

struct logicalOr{
  void  operator()(bool& f1, const bool& f2){
    f1 = f1 || f2;
  }
};

class check_folded_cell_apply_interior:public apply_rule<store<bool>, logicalOr> {
  const_Map cl;
  const_Map cr;
  const_store<bool> isFolded;
  store<bool> isIndivisible;
public:
  check_folded_cell_apply_interior(){
    name_store("cl", cl);
    name_store("cr", cr);
    name_store("isFolded", isFolded);
    name_store("isIndivisible", isIndivisible);
    input("isFolded");
    input("(cl, cr)->isIndivisible");
    output("(cl, cr)->isIndivisible");
    constraint("interior_faces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity f){
    //  join(isIndivisible[cl[f]], isFolded[f]);
    // join(isIndivisible[cr[f]], isFolded[f]);
   
    isIndivisible[cl[f]]  = isIndivisible[cl[f]] || isFolded[f];
    isIndivisible[cr[f]]  = isIndivisible[cr[f]] || isFolded[f];
  }
};

register_rule<check_folded_cell_apply_interior> register_check_folded_cell_apply_interior;
  
class check_folded_cell_apply_boundary:public apply_rule<store<bool>, logicalOr> {
  const_Map cl;
  const_store<bool> isFolded;
  store<bool> isIndivisible;
public:
  check_folded_cell_apply_boundary(){
    name_store("cl", cl);
    name_store("isFolded", isFolded);
    name_store("isIndivisible", isIndivisible);
    input("isFolded");
    input("cl->isIndivisible");
    output("cl->isIndivisible");
    constraint("boundary_faces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
  }
  void calculate(Entity f){
    //join(isIndivisible[cl[f]], isFolded[f]);
    
    isIndivisible[cl[f]]  = isIndivisible[cl[f]] || isFolded[f];
  }
};

register_rule<check_folded_cell_apply_boundary> register_check_folded_cell_apply_boundary;
  
