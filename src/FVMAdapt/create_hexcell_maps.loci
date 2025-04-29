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
#include <algorithm>
#include <iostream>
#include <Loci.h>
#include <Tools/tools.h>
#include "hexcell.h"
#include "defines.h"

using std::cerr;
using std::cout;
using std::endl;
using Loci::storeRepP;
using std::vector;

void reorder_faces(const const_store<int>& node_remap, std::vector<Entity>& lower) ;
void reorder_faces(const const_store<int>& node_remap, std::vector<Entity>& lower,
                   std::vector<Entity>& upper,
                   std::vector<Entity>& boundary_map);











/* create_hexcellmaps
   * Create an ordering for the faces and nodes of a hex cell.
   *
   * hex2face	This is an array of 6 integers which give an ordering
   *		to the entities which appear in the cell2face map.
   *		This ordering explicity maps the faces to the local
   *		coordinates of the hex cell as follows.
   *			0	xi = 1
   *			1	xi  = 0
   *			2	eta   = 1
   *			3	eta   =  0
   *			4	zeta  =  1
   *			5	zeta =  0
   *
   * hex2node	This is an array of 8 integers which give an ordering
   *		to the entities which appear in the face2node maps
   *		of hex2face[4] and hex2face[5].  The  entries 0, 2, 4, 6
   *		refer to nodes of hex2face[5] (zeta=0), while the entries 
   *		1, 3, ,5 7 refer to nodes of hex2face[4] (zeta=1).
   *		This ordering explicitly maps the nodes to the local
   *		coordinates of the hex cell as follows.
   *			0	(0, 0, 0)
   *			1	( 0, 0, 1)
   *			2	(0,  1, 0)
   *			3	( 0,  1, 1)
   *			4	(1, 0,  0)
   *			5	(1, 0,  1)
   *			6	(1,  1,  0)
   *			7	( 1,  1,  1)
   */
  class create_hexcell_maps : public pointwise_rule {
    /** lower:	Mapping from cells to faces */
    const_multiMap lower;
    /** upper:	Mapping from cells to faces */
    const_multiMap upper;
    /** boundary_map: Mapping from cells to faces */
    const_multiMap boundary_map;
    /** f2n:	Mapping from faces to nodes */
    const_multiMap face2node;
   
    store<Array<char,6> > rot;
    store<Array<char,6> > hex2face;
    store<Array<char,8> > hex2node;
    const_store<int>  node_l2f;
  public:
    create_hexcell_maps() {
      name_store("lower", lower);
      name_store("upper", upper);
      name_store("boundary_map", boundary_map);
      name_store("face2node", face2node);
      name_store("hexOrientCode", rot);
      name_store("hex2face", hex2face);
      name_store("hex2node", hex2node);
      name_store("fileNumber(face2node)", node_l2f);
      input("(upper,lower,boundary_map)->(face2node, fileNumber(face2node))");
   
      output("hexOrientCode");
      output("hex2face");
      output("hex2node");
      
      constraint("hexcells");
    }
    virtual void compute(const sequence & seq) {
      if(seq.size()!=0){
   
        do_loop(seq, this);
      }
    }
    void calculate(Entity cc) {
      // Collect the entity designations of all nodes in the nodes array
      Entity nodes[8];
      for(size_t i=0;i<8;++i)
        nodes[i] = 0 ;
      
      // Collect the entity designations of all faces in the faces array
        Entity faces[6];
      // Determine the relative orientations of the faces
      char orient[6];
      //first create vectors for reordering
      vector<Entity> vlower(lower.num_elems(cc));
      vector<Entity> vupper(upper.num_elems(cc));
      vector<Entity> vboundary_map(boundary_map.num_elems(cc));
      int nf =0;
      for (int f=0; f<lower.num_elems(cc); f++) vlower[nf++] =lower[cc][f]; 
      nf =0;
      for (int f=0; f<upper.num_elems(cc); f++) vupper[nf++] =upper[cc][f]; 
      nf =0;
      for (int f=0; f<boundary_map.num_elems(cc); f++) vboundary_map[nf++] =boundary_map[cc][f]; 
      reorder_faces(node_l2f, vlower, vupper, vboundary_map);
      
      // Initialize topology deduction data structures
      //this piece of code should work for parallel sine entityset not used
      nf=0;
      for (int f=0; f<lower.num_elems(cc); f++) {
	orient[nf] = 'r';
	hex2face[cc][nf] = f;
	rot[cc][nf] = 0;
	faces[nf++] = vlower[f];
      }
      for (int f=0; f<upper.num_elems(cc); f++) {
	orient[nf] = 'l';
	hex2face[cc][nf] = f + 6;
	rot[cc][nf] = 0;
	faces[nf++] = vupper[f];
      }
      for (int f=0; f<boundary_map.num_elems(cc); f++) {
	orient[nf] = 'l';
	hex2face[cc][nf] = f + 12;
	rot[cc][nf] = 0;
	faces[nf++] = vboundary_map[f];
      }


      
      // The first node of faces[0] will be nodes[0]
      nodes[0] = face2node[faces[0]][0];
      hex2node[cc][0] = 0;
      // The third node of faces[0] will be nodes[3]
      nodes[3] = face2node[faces[0]][2];
      hex2node[cc][3] = 2;
      // Determine the orientation of the first face.
      // If orient[0]=='r', then faces[0] has an inward facing normal.
      if (orient[0] == 'r') {
	// The second node of faces[0] will be nodes[1]
	nodes[1] = face2node[faces[0]][1];
	hex2node[cc][1] = 1;
	// The fourth node of faces[0] will be nodes[2]
	nodes[2] = face2node[faces[0]][3];
	hex2node[cc][2] = 3;
	// Set the rotation parameter to zero (no rotations, no flip)
	rot[cc][0] = 0; //qh: this is right for DOWN
      }
      // Otherwise, it's an outward facing normal, and we need to
      // reverse the direction of the first four nodes.
      else {
	// The fourth node of faces[0] will be nodes[1]
	nodes[1] = face2node[faces[0]][3];
	hex2node[cc][1] = 3;
	// The second node of faces[0] will be nodes[2]
	nodes[2] = face2node[faces[0]][1];
	hex2node[cc][2] = 1;
	// Set the rotation parameter to four (no rotations, one flip)
	rot[cc][0] = 4; //qh: this is right for DOWN
      }
      // Next, let's flag the other two faces which share nodes[0].
      // First up, faces[1] (eta=-1) will share nodes[0] and nodes[1]
      for (int f=1; f<6; f++) {
	for (int n=0; n<4; n++) {
	  if (face2node[faces[f]][n] == nodes[0]) {
	    // If the next node in the face2node map is nodes[1]
	    // then this face has an outward facing normal.
	    if (face2node[faces[f]][(n+1)%4] == nodes[1]) {
	      // Double check our deduction
	      if (orient[f] != 'l') {
		cerr << "WARNING[1a]: Face orientation test failed!" << endl;
		throw -1;
	      }
	      // (n rotations, one flip)
              // rot[cc][1] = n + 4;//qh: need change here
              rot[cc][1] = n ; //qh: if 'l', no flip 
	      // Make this face, faces[1]
	      std::swap(hex2face[cc][1], hex2face[cc][f]);
	      std::swap(faces[1], faces[f]);
	      std::swap(orient[1], orient[f]);
	      // Store the other corners of the face in the node list
	      nodes[5] = face2node[faces[1]][(n+2)%4];
	      nodes[4] = face2node[faces[1]][(n+3)%4];
	      // Expedite the exit from these loops
	      n = 4;
	      f = 6;
	    }
	    // If the previous node in the face2node map is nodes[1]
	    // then this face has an inward facing normal.
	    else if (face2node[faces[f]][(n+3)%4] == nodes[1]) {
	      // Double check our deduction
	      if (orient[f] != 'r') {
		cerr << "WARNING[1A]: Face orientation test failed!" << endl;
		throw -1;
	      }
	      // (n rotations, no flip)
              //  rot[cc][1] = n;//qh: need change here
              rot[cc][1] = n + 4; //qh: if 'r', flip
	      // Make this face, faces[1]
	      std::swap(hex2face[cc][1], hex2face[cc][f]);
	      std::swap(faces[1], faces[f]);
	      std::swap(orient[1], orient[f]);
	      // Store the other corners of the face in the node list
	      nodes[5] = face2node[faces[1]][(n+2)%4];
	      nodes[4] = face2node[faces[1]][(n+1)%4];
	      // Expedite the exit from these loops
	      n = 4;
	      f = 6;
	    }
	  }
	}
      }
      // Next, faces[2] (xi=-1) will share nodes[0] and nodes[2]
      for (int f=2; f<6; f++) {
	for (int n=0; n<4; n++) {
	  if (face2node[faces[f]][n] == nodes[0]) {
	    // If the next node in the face2node map is nodes[2]
	    // then this face has an inward facing normal.
	    if (face2node[faces[f]][(n+1)%4] == nodes[2]) {
	      // Double check our deduction
	      if (orient[f] != 'r') {
		cerr << "WARNING[2a]: Face orientation test failed!" << endl;
		throw -1;
	      }
	      // nodes[4] is also shared with faces[1]
	      if (face2node[faces[f]][(n+3)%4] != nodes[4]) {
		cerr << "WARNING[2b]: Coincident node test failed!" << endl;
		throw -1;
	      }
	      // (n rotations, no flip)
	      rot[cc][2] = n;
	      // Make this face, faces[2]
	      std::swap(hex2face[cc][2], hex2face[cc][f]);
	      std::swap(faces[2], faces[f]);
	      std::swap(orient[2], orient[f]);
	      // Store the remaining corner of this face in the node list
	      nodes[6] = face2node[faces[2]][(n+2)%4];
	      // Expedite the exit from these loops
	      n = 4;
	      f = 6;
	    }
	    // If the previous node in the face2node map is nodes[2]
	    // then this face has an outward facing normal.
	    else if (face2node[faces[f]][(n+3)%4] == nodes[2]) {
	      // Double check our deduction
	      if (orient[f] != 'l') {
		cerr << "WARNING[2A]: Face orientation test failed!" << endl;
		throw -1;
	      }
	      // nodes[4] is also shared with faces[1]
	      if (face2node[faces[f]][(n+1)%4] != nodes[4]) {
		cerr << "WARNING[2B]: Coincident node test failed!" << endl;
		throw -1;
	      }
	      // (n rotations, one flip)
	      rot[cc][2] = n + 4;
	      // Make this face, faces[2]
	      std::swap(hex2face[cc][2], hex2face[cc][f]);
	      std::swap(faces[2], faces[f]);
	      std::swap(orient[2], orient[f]);
	      // Store the remaining corner of this face in the node list
	      nodes[6] = face2node[faces[2]][(n+2)%4];
	      // Expedite the exit from these loops
	      n = 4;
	      f = 6;
	    }
	  }
	}
      }
      // Now, faces[3] (xi=1) which will share nodes[1] and nodes[3]
      for (int f=3; f<6; f++) {
	for (int n=0; n<4; n++) {
	  if (face2node[faces[f]][n] == nodes[1]) {
	    // If the next node in the face2node map is nodes[3]
	    // then this face has an outward facing normal
	    if (face2node[faces[f]][(n+1)%4] == nodes[3]) {
	      // Double check our deduction
	      if (orient[f] != 'l') {
		cerr << "WARNING[3a]: Face orientation test failed!" << endl;
		throw -1;
	      }
	      // nodes[5] is also shared with faces[1]
	      if (face2node[faces[f]][(n+3)%4] != nodes[5]) {
		cerr << "WARNING[3b]: Coincident node test failed!" << endl;
		throw -1;
	      }
	      // (n rotations, no flip)
	      rot[cc][3] = n;
	      // Make this face, faces[3]
	      std::swap(hex2face[cc][3], hex2face[cc][f]);
	      std::swap(faces[3], faces[f]);
	      std::swap(orient[3], orient[f]);
	      // Store the opposite corner in nodes[7]
	      nodes[7] = face2node[faces[3]][(n+2)%4];
	      // Expedite the exit from these loops
	      n = 4;
	      f = 6;
	    }
	    // Otherwise the normal is inward facing
	    else if (face2node[faces[f]][(n+3)%4] == nodes[3]) {
	      // Double check our deduction
	      if (orient[f] != 'r') {
		cerr << "WARNING[3A]: Face orientation test failed!" << endl;
		throw -1;
	      }
	      // nodes[5] is also shared with faces[1]
	      if (face2node[faces[f]][(n+1)%4] != nodes[5]) {
		cerr << "WARNING[3B]: Coincident node test failed!" << endl;
		throw -1;
	      }
	      // (n rotations, one flip)
	      rot[cc][3] = n + 4;
	      // Make this face, faces[3]
	      std::swap(hex2face[cc][3], hex2face[cc][f]);
	      std::swap(faces[3], faces[f]);
	      std::swap(orient[3], orient[f]);
	      // Store the opposite corner in nodes[7]
	      nodes[7] = face2node[faces[3]][(n+2)%4];
	      // Expedite the exit from these loops
	      n = 4;
	      f = 6;
	    }
	  }
	}
      }
      // And faces[4] (eta=1) which will share nodes[2] and nodes[3]
      for (int f=4; f<6; f++) {
	for (int n=0; n<4; n++) {
	  if (face2node[faces[f]][n] == nodes[2]) {
	    // If the next node in the face2node map is nodes[3]
	    // then this face has an inward facing normal
	    if (face2node[faces[f]][(n+1)%4] == nodes[3]) {
	      // Double check our deduction
	      if (orient[f] != 'r') {
		cerr << "WARNING[4a]: Face orientation test failed!" << endl;
		throw -1;
	      }
	      // nodes[6] is also shared with faces[2]
	      if (face2node[faces[f]][(n+3)%4] != nodes[6]) {
		cerr << "WARNING[4b]: Coincident node test failed!" << endl;
		throw -1;
	      }
	      // nodes[7] is also shared with faces[3]
	      if (face2node[faces[f]][(n+2)%4] != nodes[7]) {
		cerr << "WARNING[4c]: Coincident node test failed!" << endl;
		throw -1;
	      }
	      // (n rotations, one flip)
              // rot[cc][4] = n + 4; //qh:change here
              rot[cc][4] = n;//qh: for FRONT face, if'r', no flip
	      // Make this face, faces[4]
	      std::swap(hex2face[cc][4], hex2face[cc][f]);
	      std::swap(faces[4], faces[f]);
	      std::swap(orient[4], orient[f]);
	      // Expedite the exit from these loops
	      n = 4;
	      f = 6;
	    }
	    // Otherwise the normal is outward facing
	    else if (face2node[faces[f]][(n+3)%4] == nodes[3]) {
	      // Double check our deduction
	      if (orient[f] != 'l') {
		cerr << "WARNING[4A]: Face orientation test failed!" << endl;
		throw -1;
	      }
	      // nodes[6] is also shared with faces[2]
	      if (face2node[faces[f]][(n+1)%4] != nodes[6]) {
		cerr << "WARNING[4B]: Coincident node test failed!" << endl;
		throw -1;
	      }
	      // nodes[7] is also shared with faces[3]
	      if (face2node[faces[f]][(n+2)%4] != nodes[7]) {
		cerr << "WARNING[4C]: Coincident node test failed!" << endl;
		throw -1;
	      }
	      // (n rotations, no flip)
	      //rot[cc][4] = n; //qh:need change here
              rot[cc][4] = n + 4; //qh: for FRONT face, if 'l', flip
	      // Make this face, faces[4]
	      std::swap(hex2face[cc][4], hex2face[cc][f]);
	      std::swap(faces[4], faces[f]);
	      std::swap(orient[4], orient[f]);
	      // Expedite the exit from these loops
	      n = 4;
	      f = 6;
	    }
	  }
	}
      }
      // And finally faces[5] (zeta=1), should be the only one remaining
      for (int n=0; n<4; n++) {
	// faces[5] should share nodes[4], nodes[5], nodes[6], and nodes[7]
	if (face2node[faces[5]][n] == nodes[4]) {
	  // If the next node in the face2node map is nodes[5]
	  // then this face has an outward facing normal
	  if (face2node[faces[5]][(n+1)%4] == nodes[5]) {
	    // Double check our deduction
	    if (orient[5] != 'l') {
	      cerr << "WARNING[5a]: Face orientation test failed!" << endl;
	      throw -1;
	    }
	    // nodes[6] is also shared with faces[2]
	    if (face2node[faces[5]][(n+3)%4] != nodes[6]) {
	      cerr << "WARNING[5b]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // nodes[7] is also shared with faces[3]
	    if (face2node[faces[5]][(n+2)%4] != nodes[7]) {
	      cerr << "WARNING[5c]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // (n rotations, no flip)
	    rot[cc][5] = n;
	    // Node n of faces[5] is nodes[4]
            // hex2node[cc][4] = 5*8 + n;//qh:
            hex2node[cc][4] =  n;//qh: changed next 4 lines
	    // Node (n+1)%4 of faces[5] is nodes[5]
	    hex2node[cc][5] =  (n+1)%4;
	    // Node (n+2)%4 of faces[5] is nodes[7]
	    hex2node[cc][7] = (n+2)%4;
	    // Node (n+3)%4 of faces[5] is nodes[6]
	    hex2node[cc][6] =  (n+3)%4;
	  }
	  // Otherwise the normal is inward facing
	  else if (face2node[faces[5]][(n+3)%4] == nodes[5]) {
	    // Double check our deduction
	    if (orient[5] != 'r') {
	      cerr << "WARNING[5A]: Face orientation test failed!" << endl;
	      throw -1;
	    }
	    // nodes[6] is also shared with faces[2]
	    if (face2node[faces[5]][(n+1)%4] != nodes[6]) {
	      cerr << "WARNING[5B]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // nodes[7] is also shared with faces[3]
	    if (face2node[faces[5]][(n+2)%4] != nodes[7]) {
	      cerr << "WARNING[5C]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // (n rotations, one flip)
	    rot[cc][5] = n + 4;
	    // Node n of faces[5] is nodes[4]
            // hex2node[cc][4] = 5*8 + n;//qh:
            hex2node[cc][4] =  n; //qh: changed next 4 lines
	    // Node (n+1)%4 of faces[5] is nodes[6]
	    hex2node[cc][6] =  (n+1)%4;
	    // Node (n+2)%4 of faces[5] is nodes[7]
	    hex2node[cc][7] =  (n+2)%4;
	    // Node (n+3)%4 of faces[5] is nodes[5]
	    hex2node[cc][5] =  (n+3)%4;
	  }
	}
      
      }
      
            //qh:: switch to my numbering system
      Array<char, 6>h2f_temp;
      h2f_temp = hex2face[cc];
      Array<char, 6> rot_temp;
      rot_temp = rot[cc];
      
      hex2face[cc][0] = h2f_temp[3];
      rot[cc][0] = rot_temp[3];
      
      hex2face[cc][1] = h2f_temp[2];
      rot[cc][1] = rot_temp[2];
      
      hex2face[cc][2] = h2f_temp[4];
      rot[cc][2] = rot_temp[4];
      
      hex2face[cc][3] = h2f_temp[1];
      rot[cc][3] = rot_temp[1];
      
      hex2face[cc][4] = h2f_temp[5];
      rot[cc][4] = rot_temp[5];
      
      hex2face[cc][5] = h2f_temp[0];
      rot[cc][5] = rot_temp[0];
      
      //for print out
      /*
      if(true){
        for (int i = 0; i < 8; i++){
          hex2node[cc][i] += '0';
        }
        for(int i = 0; i < 6; i++){
          hex2face[cc][i] = hex2face[cc][i]%6 + '0';
          rot[cc][i] += '0';
        }
      }
      */
    }
  
};

register_rule<create_hexcell_maps> register_create_hexcell_maps;


  class determine_fr : public pointwise_rule {
    const_Map cr;
    const_multiMap lower;
    const_multiMap upper;
    const_multiMap boundary_map;
    const_store<Array<char,6> > hex2face;
    const_store<int> node_l2f;
    store<char> fr;
  public:
    determine_fr() {
      name_store("cr", cr);
      name_store("lower", lower);
      name_store("upper", upper);
      name_store("boundary_map", boundary_map);
      name_store("hex2face", hex2face);
      name_store("fr", fr);
      name_store("fileNumber(face2node)", node_l2f);
      input("cr->hex2face");
      input("cr->(lower,upper,boundary_map)->fileNumber(face2node)");
      input("cr->(lower,upper,boundary_map)");
      output("fr");
      constraint("cr->hexcells");
    }
    virtual void compute(const sequence & seq) {
      if(seq.size()!=0){
     
      do_loop(seq, this);
      }
    }
    void calculate(Entity ff) {
      // Collect the entity designations of all faces in the faces array
      Array<Entity, 6> faces = collect_hex_faces(lower[cr[ff]].begin(),lower.num_elems(cr[ff]),
                                                 upper[cr[ff]].begin(),upper.num_elems(cr[ff]),
                                                 boundary_map[cr[ff]].begin(), boundary_map.num_elems(cr[ff]),
                                                 hex2face[cr[ff]], node_l2f);
    
      // Store the original location of each face entity in faces
      int f;
      for ( f=0; f<6; f++) {
        if (node_l2f[ff] == node_l2f[faces[f]]) {
          fr[ff] = f;
          break;
        }
      }
      if(f ==6){
        cerr << "WARNING: can not find fr" << endl;
        Loci::Abort();
      }      
      //for print out
      /*      if(true){
        fr[ff] += '0';
        }*/
    }
  };

  register_rule<determine_fr> register_determine_fr;

 class determine_fl : public pointwise_rule {
   const_Map cl;
   const_multiMap lower;
   const_multiMap upper;
   const_multiMap boundary_map;
   
   const_store<int> node_l2f;
   const_store<Array<char,6> > hex2face;
   store<char> fl;
  public:
    determine_fl() {
      name_store("cl", cl);
      name_store("lower", lower);
      name_store("upper", upper);
      name_store("boundary_map", boundary_map);
      name_store("hex2face", hex2face);
      name_store("fl", fl);
      name_store("fileNumber(face2node)", node_l2f);
      
    
      input("cl->hex2face");
      input("cl->(lower,upper,boundary_map)");
      input("cl->(lower,upper,boundary_map)->fileNumber(face2node)");
      output("fl");
      constraint("cl->hexcells");
    }
   virtual void compute(const sequence & seq) {
     if(seq.size()!=0){
  
     do_loop(seq, this);
     }
   }
   void calculate(Entity ff) {
      // Collect the entity designations of all faces in the faces array
     Array<Entity, 6> faces = collect_hex_faces(lower[cl[ff]].begin(), lower.num_elems(cl[ff]),
                                                upper[cl[ff]].begin(), upper.num_elems(cl[ff]),
                                                boundary_map[cl[ff]].begin(),boundary_map.num_elems(cl[ff]),
                                                hex2face[cl[ff]],node_l2f);
      
      // Store the original location of each face entity in faces
      int f;
      for ( f=0; f<6; f++) {
        if (node_l2f[ff] == node_l2f[faces[f]]) {
          fl[ff] = f;
          break;
        }
      }
      if(f ==6){
        cerr << "WARNING: can not find fl" << endl;
        Loci::Abort();
      }
      //for print out
      /* if(true){
        fl[ff] += '0';
        }*/
    }
  };
 
  register_rule<determine_fl> register_determine_fl;

