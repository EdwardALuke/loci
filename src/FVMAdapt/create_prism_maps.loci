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
#include "hex_defines.h"

using std::cerr;
using std::cout;
using std::endl;
using Loci::storeRepP;
using std::vector;

Array<Entity, 5> collect_prism_faces( const Entity* lower, int lower_size,
                                      const Entity* upper,int upper_size,
                                      const Entity* boundary_map,int boundary_map_size,
                                      const Array<char, 5>& prism2face, const const_store<int>& node_remap);
                                     
void reorder_faces(const const_store<int>& node_remap, std::vector<Entity>& lower,
                   std::vector<Entity>& upper,
                   std::vector<Entity>& boundary_map);
/* create_prismmaps
   * Create an ordering for the faces and nodes of a prism cell.
   *
   * prism2face	This is an array of 5 integers which give an ordering
   *		to the entities which appear in the cell2face map.
   * prism2node	This is an array of 6 integers which give an ordering
   *		to the entities which appear in the face2node maps
   *		of prism2triface[0] and prism2triface[1].  The first 3 entries
   *		refer to nodes of prism2face[0]  while the last
   *		3 entries refer to nodes of prism2triface[1] (zeta=1).
   *		coordinates of the prism cell as follows.
   */
  class create_prism_maps : public pointwise_rule {
    /** lower:	Mapping from cells to faces */
    const_multiMap lower;
    /** upper:	Mapping from cells to faces */
    const_multiMap upper;
    /** boundary_map: Mapping from cells to faces */
    const_multiMap boundary_map;
    /** face2node	Mapping from faces to nodes */
    const_multiMap face2node;
    store<Array<char,5> > orientCode;
    store<Array<char,5> > prism2face;
    store<Array<char,6> > prism2node;

    
  const_store<int> node_l2f;
  public:
    create_prism_maps() { 
      name_store("lower", lower);
      name_store("upper", upper);
      name_store("boundary_map", boundary_map);
      name_store("face2node", face2node);
      name_store("prismOrientCode", orientCode);
      name_store("prism2face", prism2face);
      name_store("prism2node", prism2node);
       name_store("fileNumber(face2node)", node_l2f);
      input("(upper,lower,boundary_map)->(face2node, fileNumber(face2node))");

      output("prismOrientCode");
      output("prism2face");
      output("prism2node");
      
      constraint("prisms");
    }
    virtual void compute(const sequence & seq) {
      if(seq.size()!=0){
       
        do_loop(seq, this);
      }
      
    }
    void calculate(Entity cc) {

 //first create vectors and reorder faces
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











      // Collect the entity designations of all nodes in the nodes array
      Entity nodes[6];
      for(size_t i=0;i<6;++i)
        nodes[i] = 0 ;
      // Collect the entity designations of all faces in the faces array
      Entity trifaces[2];
      Entity quadfaces[3];
      // Determine the relative orientations of the faces
      char orient_tri[2]; // value 'l' or 'r'
      char orient_quad[3]; //value 'l', or 'r'
      Array<char, 2> trirot;
      Array<char, 3> quadrot;
      Array<char, 2> prism2triface;
      Array<char, 3> prism2quadface;
      
      // Initialize topology deduction data structures
      int ntf=0;
      int nqf = 0;
      for (int f=0; f<lower.num_elems(cc); f++) {

        if(face2node.num_elems(vlower[f]) == 3){
          prism2triface[ntf] = f;
          orient_tri[ntf] = 'r';
          trirot[ntf] = 0;
          trifaces[ntf++] = vlower[f];
        }
        else if (face2node.num_elems(vlower[f]) == 4){
          prism2quadface[nqf] = f;
          orient_quad[nqf] = 'r';
          quadrot[nqf] = 0;
          quadfaces[nqf++] = vlower[f];
        }
        
          
        else{
          cerr << "WARNING: not a prism " << endl;
          Loci::Abort();
        }
      }
     
      
      for (int f=0; f<upper.num_elems(cc); f++) {
        if(face2node.num_elems(vupper[f]) == 3){
          prism2triface[ntf] = f+6;
          orient_tri[ntf] = 'l';
          trirot[ntf] = 0;
          trifaces[ntf++] = vupper[f];
        }
        else if (face2node.num_elems(vupper[f]) == 4){
          prism2quadface[nqf] = f+6;
          orient_quad[nqf] = 'l';
          quadrot[nqf] = 0;
          quadfaces[nqf++] = vupper[f];
        }
        
        
        else{
          cerr << "WARNING: not a prism " << endl;
          Loci::Abort();
        }
      }
      for (int f=0; f<boundary_map.num_elems(cc); f++) {

        if(face2node.num_elems(vboundary_map[f]) == 3){
          prism2triface[ntf] = f+12;
          orient_tri[ntf] = 'l';
          trirot[ntf] = 0;
          trifaces[ntf++] = vboundary_map[f];
        }
        else if (face2node.num_elems(vboundary_map[f]) == 4){
          prism2quadface[nqf] = f+12;
          orient_quad[nqf] = 'l';
          quadrot[nqf] = 0;
          quadfaces[nqf++] = vboundary_map[f];
        }
        
        
        else{
          cerr << "WARNING: not a prism " << endl;
          Loci::Abort();
        }
      }

      if(ntf != 2 || nqf != 3){
        cerr<<"WARNING: error in clollecting prism faces" << endl;
        Loci::Abort();
      }
      //order the faces and everything according node_remap here

      
      // The first node of trifaces[0] will be nodes[0]
      nodes[0] = face2node[trifaces[0]][0];
      prism2node[cc][0] = 0;
      
      
      
      // Determine the orientation of the first face.
      // If orient[0]=='r', then faces[0] has an inward facing normal.
      if (orient_tri[0] == 'r') {
	// The second node of trifaces[0] will be nodes[1]
	nodes[1] = face2node[trifaces[0]][1];
	prism2node[cc][1] = 1;
	// The third node of faces[0] will be nodes[2]
	nodes[2] = face2node[trifaces[0]][2];
	prism2node[cc][2] = 2;
	// Set the rotation parameter to zero (no rotations, no flip)
	trirot[0] = 0;
      }
      // Otherwise, it's an outward facing normal, and we need to
      // reverse the direction of the first 3 nodes.
      else {
	// The 3rd node of faces[0] will be nodes[1]
	nodes[1] = face2node[trifaces[0]][2];
	prism2node[cc][1] = 2;
	// The 2nd node of faces[0] will be nodes[2]
	nodes[2] = face2node[trifaces[0]][1];
	prism2node[cc][2] = 1;
	// Set the rotation parameter to four (no rotations, one flip)
	trirot[0] = 4;
      }



      // Next, let's flag the other two faces which share nodes[0].
      // First, quadfaces[0]  will share nodes[0] and nodes[1]
      for (int f=0; f<3; f++) {
	for (int n=0; n<4; n++) {
	  if (face2node[quadfaces[f]][n] == nodes[0]) {
	    // If the next node in the face2node map is nodes[1]
	    // then this face has an outward facing normal.
	    if (face2node[quadfaces[f]][(n+1)%4] == nodes[1]) {
	      // Double check our deduction
	      if (orient_quad[f] != 'l') {
		cerr << "WARNING[1a]: Face orientation test failed!" << endl;
		throw -1;
	      }
	      // (n rotations)
             quadrot[0] = n ; // no flip 
	      // Make this face, quadfaces[0]
	      std::swap(prism2quadface[0], prism2quadface[f]);
	      std::swap(quadfaces[0], quadfaces[f]);
	      std::swap(orient_quad[0], orient_quad[f]);
	      // Store the other corners of the face in the node list
	      nodes[4] = face2node[quadfaces[0]][(n+2)%4];
	      nodes[3] = face2node[quadfaces[0]][(n+3)%4];
	      // Expedite the exit from these loops
	      n = 4;
	      f = 3;
	    }
	    // If the previous node in the face2node map is nodes[1]
	    // then this face has an inward facing normal.
	    else if (face2node[quadfaces[f]][(n+3)%4] == nodes[1]) {
	      // Double check our deduction
	      if (orient_quad[f] != 'r') {
		cerr << "WARNING[1A]: Face orientation test failed!" << endl;
		throw -1;
	      }
              quadrot[0] = n + 4; //qh: if 'r', flip
	      // Make this face, faces[1]
	      std::swap(prism2quadface[0], prism2quadface[f]);
	      std::swap(quadfaces[0], quadfaces[f]);
	      std::swap(orient_quad[0], orient_quad[f]);
	      // Store the other corners of the face in the node list
	      nodes[4] = face2node[quadfaces[0]][(n+2)%4];
	      nodes[3] = face2node[quadfaces[0]][(n+1)%4];
	      // Expedite the exit from these loops
	      n = 4;
	      f = 3;
	    }
	  }
	}
      }
      // Next, quadfaces[2](2, 0, 3, 5)  will share nodes[2] and nodes[0]
      for (int f=1; f<3; f++) {
	for (int n=0; n<4; n++) {
	  if (face2node[quadfaces[f]][n] == nodes[2]) {
	    // If the next node in the face2node map is nodes[0]
	    // then this face has an outward facing normal.
	    if (face2node[quadfaces[f]][(n+1)%4] == nodes[0]) {
	      // Double check our deduction
	      if (orient_quad[f] != 'l') {
		cerr << "WARNING[2a]: Face orientation test failed!" << endl;
		throw -1;
	      }
              // (n rotations,   no flip)
	      quadrot[2] = n;
	      // Make this face, faces[2]
	      std::swap(prism2quadface[2], prism2quadface[f]);
	      std::swap(quadfaces[2], quadfaces[f]);
	      std::swap(orient_quad[2], orient_quad[f]);
	      // Store the remaining corner of this face in the node list
              nodes[3] = face2node[quadfaces[2]][(n+2)%4];
              nodes[5] = face2node[quadfaces[2]][(n+3)%4];
              // Expedite the exit from these loops
              n = 4;
	      f = 3;
	    }
            // If the previous node in the face2node map is nodes[0]
	    // then this face has an inward facing normal.
            else if (face2node[quadfaces[f]][(n+3)%4] == nodes[0]) {
	      // Double check our deduction
	      if (orient_quad[f] != 'r') {
		cerr << "WARNING[2A]: Face orientation test failed!" << endl;
		throw -1;
	      }
              
	      // (n rotations, one flip)
              quadrot[2] = n + 4;
	      // Make this face, faces[2]
	      std::swap(prism2quadface[2], prism2quadface[f]);
	      std::swap(quadfaces[2], quadfaces[f]);
	      std::swap(orient_quad[2], orient_quad[f]);
	      // Store the remaining corner of this face in the node list
              nodes[3] = face2node[quadfaces[2]][(n+2)%4];
              nodes[5] = face2node[quadfaces[2]][(n+1)%4];
              
	      // Expedite the exit from these loops
	      n = 4;
	      f = 3;
	    }
	  }
        }
      }
    
      // And trifaces[1]  will share nodes[3] and nodes[4] and nodes[5]
      for (int n=0; n<4; n++) {
        if (face2node[trifaces[1]][n] == nodes[3]) {
          // If the next node in the face2node map is nodes[4]
          // then this face has an outward facing normal
          if (face2node[trifaces[1]][(n+1)%3] == nodes[4]) {
            // Double check our deduction
            if (orient_tri[1] != 'l') {
              cerr << "WARNING[4a]: Face orientation test failed!" << endl;
              throw -1;
            }
            // nodes[5] is also shared with trifaces[1]
            if (face2node[trifaces[1]][(n+2)%3] != nodes[5]) {
              cerr << "WARNING[4b]: Coincident node test failed!" << endl;
              throw -1;
            }
              
            // (n rotations, no flip)
            trirot[1] = n;
            prism2node[cc][3] =  n;//qh: changed next 4 lines
            prism2node[cc][4] =  (n+1)%3;
            prism2node[cc][5] = (n+2)%3;
            n = 4;
          }
          // Otherwise the normal is inward facing
          else if (face2node[trifaces[1]][(n+1)%3] == nodes[5]) {
            // Double check our deduction
            if (orient_tri[1] != 'r') {
              cerr << "WARNING[4A]: Face orientation test failed!" << endl;
              throw -1;
            }
            // nodes[4] is also shared with trifaces[1]
            if (face2node[trifaces[1]][(n+2)%3] != nodes[4]) {
              cerr << "WARNING[4B]: Coincident node test failed!" << endl;
              throw -1;
            }
            
            // (n rotations, one flip)
	     
            trirot[1] = n + 4; 
            prism2node[cc][3] =  n;//qh: changed next 4 lines
            prism2node[cc][4] =  (n+2)%3;
            prism2node[cc][5] = (n+1)%3;
            n = 4;
            
          }
        }
      }
      // And finally quadfaces[1] should be the only one remaining
      for (int n=0; n<4; n++) {
	// quadfaces[1] should share nodes[1], nodes[2], nodes[5], and nodes[4]
	if (face2node[quadfaces[1]][n] == nodes[1]) {
	  // If the next node in the face2node map is nodes[2]
	  // then this face has an outward facing normal
	  if (face2node[quadfaces[1]][(n+1)%4] == nodes[2]) {
	    // Double check our deduction
	    if (orient_quad[1] != 'l') {
	      cerr << "WARNING[5a]: Face orientation test failed!" << endl;
	      throw -1;
	    }
	    // nodes[5] is also shared with quadfaces[1]
	    if (face2node[quadfaces[1]][(n+2)%4] != nodes[5]) {
	      cerr << "WARNING[5b]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // nodes[4] is also shared with faces[3]
	    if (face2node[quadfaces[1]][(n+3)%4] != nodes[4]) {
	      cerr << "WARNING[5c]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // (n rotations, no flip)
	    quadrot[1] = n;
	    
          
	  }
	  // Otherwise the normal is inward facing
	  else if (face2node[quadfaces[1]][(n+1)%4] == nodes[4]) {
	    // Double check our deduction
	    if (orient_quad[1] != 'r') {
	      cerr << "WARNING[5A]: Face orientation test failed!" << endl;
	      throw -1;
	    }
	    // nodes[6] is also shared with faces[2]
	    if (face2node[quadfaces[1]][(n+2)%4] != nodes[5]) {
	      cerr << "WARNING[5B]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // nodes[7] is also shared with faces[3]
	    if (face2node[quadfaces[1]][(n+3)%4] != nodes[2]) {
	      cerr << "WARNING[5C]: Coincident node test failed!" << endl;
	      throw -1;
	    }
	    // (n rotations, one flip)
	   quadrot[1] = n + 4;
          }
	}
        
      }
      
      for(int i = 0; i < 2; i++){
        prism2face[cc][i] = prism2triface[i];
        orientCode[cc][i] = trirot[i];
      }
      for(int i = 0; i < 3; i++){
        prism2face[cc][i+2] = prism2quadface[i];
        orientCode[cc][i+2] = quadrot[i];
      }  


      //for print out
      
      if(false){
        cout << "cc " << cc << endl;
        cout << " prism2node: " ;
        for (int i = 0; i < 6; i++){
          cout << char(prism2node[cc][i] + '0') << " ";
        }
        cout << endl;
        
        for(int i = 0; i < 2; i++){
          cout <<  int(prism2triface[i]+ int(0)) << "; " << int(trirot[i]+int(0)) <<  "           ";
          
        }
        cout << endl;
        for(int i = 0; i < 3; i++){
          cout <<  int(prism2quadface[i]) << "; " << int(quadrot[i]) <<  "           ";
          
        }
        cout << endl;
        
        
      }
      
    }//calculate
    
  };

register_rule<create_prism_maps> register_create_prism_maps;


class determine_prism_fr : public pointwise_rule {
  const_Map cr;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,5> > prism2face;

  const_store<int> node_l2f;
  store<char> fr;
public:
  determine_prism_fr() {
      name_store("cr", cr);
      name_store("lower", lower);
      name_store("upper", upper);
      name_store("boundary_map", boundary_map);
      name_store("prism2face", prism2face);
      name_store("fileNumber(face2node)", node_l2f);
      name_store("fr", fr);
      
    
      input("cr->prism2face");
      input("cr->(lower,upper,boundary_map)");
      input("cr->(lower,upper,boundary_map)->fileNumber(face2node)");
      output("fr");
      constraint("cr->prisms");
    }
    virtual void compute(const sequence & seq) {
      if(seq.size()!=0){
       
        do_loop(seq, this);
      }
    }
  void calculate(Entity ff) {
    // Collect the entity designations of all faces in the faces array
    Array<Entity, 5> faces = collect_prism_faces(lower[cr[ff]].begin(), lower.num_elems(cr[ff]),
                                                 upper[cr[ff]].begin(), upper.num_elems(cr[ff]),
                                                 boundary_map[cr[ff]].begin(), boundary_map.num_elems(cr[ff]),
                                                 prism2face[cr[ff]], node_l2f);
    
    // Store the original location of each face entity in faces
    int f;
    for ( f=0; f<5; f++) {
      if (node_l2f[ff] == node_l2f[faces[f]]) {
        fr[ff] = f;
        break;
      }
    }
      if(f ==5){
        cerr << "WARNING: can not find fr" << endl;
        Loci::Abort();
      }      
      //for print out
     //  if(false){
//         cout << "ff: " << ff << "fr[ff]: " <<char(fr[ff] + '0') <<endl ;
//       }
  }
};

register_rule<determine_prism_fr> register_determine_prism_fr;

class determine_prism_fl : public pointwise_rule {
  const_Map cl;
  const_multiMap lower;
  const_multiMap upper;
  const_multiMap boundary_map;
  const_store<Array<char,5> > prism2face;
  const_store<int> node_l2f;
 
  store<char> fl;
public:
  determine_prism_fl() {
      name_store("cl", cl);
      name_store("lower", lower);
      name_store("upper", upper);
      name_store("boundary_map", boundary_map);
      name_store("prism2face", prism2face);
      name_store("fileNumber(face2node)", node_l2f);
      name_store("fl", fl);
      
    
      input("cl->prism2face");
      input("cl->(lower,upper,boundary_map)");
       input("cl->(lower,upper,boundary_map)->fileNumber(face2node)");
      output("fl");
      constraint("cl->prisms");
  }
  virtual void compute(const sequence & seq) {
    if(seq.size()!=0){
     
      do_loop(seq, this);
    }
  }
  void calculate(Entity ff) {
      // Collect the entity designations of all faces in the faces array
    Array<Entity, 5> faces = collect_prism_faces(lower[cl[ff]].begin(), lower.num_elems(cl[ff]),
                                                 upper[cl[ff]].begin(), upper.num_elems(cl[ff]),
                                                 boundary_map[cl[ff]].begin(), boundary_map.num_elems(cl[ff]),
                                                 prism2face[cl[ff]], node_l2f);
    
    // Store the original location of each face entity in faces
    int f;
    for ( f=0; f<5; f++) {
      if (node_l2f[ff] == node_l2f[faces[f]]) {
        fl[ff] = f;
          break;
      }
    }
    if(f ==5){
      cerr << "WARNING: can not find fl" << endl;
        Loci::Abort();
    }
    // if(false){
//       cout << "ff: " << ff << "fr[ff]: " <<char(fl[ff] + '0') <<endl;
//     } 
  }
};

register_rule<determine_prism_fl> register_determine_prism_fl;

