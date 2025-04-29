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
#include <queue>
#include "quadface.h"
#include "hex_defines.h"
using std::cout;
using std::endl;
//dummy rules to create constraint quadface

class quad_face1:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face1(){
   
    name_store("quadfaces", quadfaces);
    output("quadfaces");
    constraint("boundary_faces, cl->hexcells");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face1> register_quad_face1;

class quad_face2:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face2(){
   
    name_store("quadfaces", quadfaces);
    output("quadfaces");
    constraint("boundary_faces, quadrangles, cl->prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face2> register_quad_face2;

class quad_face3:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face3(){
   
    name_store("priority::cr::quadfaces", quadfaces);
    output("priority::cr::quadfaces");
    constraint("interior_faces, quadrangles,  cl->prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face3> register_quad_face3;

class quad_face4:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face4(){
   
    name_store("priority::cr::quadfaces", quadfaces);
    output("priority::cr::quadfaces");
    constraint("interior_faces, cl->hexcells");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face4> register_quad_face4;

class quad_face5:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face5(){
   
    name_store("cr::quadfaces", quadfaces);
    output("cr::quadfaces");
    constraint("interior_faces, quadrangles, cr->prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face5> register_quad_face5;

class quad_face6:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face6(){
   
    name_store("cr::quadfaces", quadfaces);
    output("cr::quadfaces");
    constraint("interior_faces, cr->hexcells");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face6> register_quad_face6;


//class to check if a face is quadface
class identify_face:public pointwise_rule{
 
  store<bool > is_quadface;

public:
  identify_face(){
   
    name_store("quadfaces::is_quadface", is_quadface);
    output("quadfaces::is_quadface");
    constraint("faces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){
    is_quadface[f] = false;
  }
};

register_rule<identify_face> register_identify_face;


class identify_quadface:public pointwise_rule{
 
  store<bool > is_quadface;

public:
  identify_quadface(){
   
    name_store("priority::quadfaces::is_quadface", is_quadface);
    output("priority::quadfaces::is_quadface");
    constraint("quadfaces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
    //std::cout << "quadface in a rule: " << seq.size() << " myID " << Loci::MPI_rank << std::endl;
   
  }
  void calculate(Entity f){
    is_quadface[f] = true;
  }
};

register_rule<identify_quadface> register_identify_quadface;
