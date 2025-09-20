//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
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
////////////////////////////////////////////////////////////////////////////////////////////
// This file reads in the .surf and .surfn files, compute the geometry from exact normals,
// and output the geometry into .coeff file
//
///////////////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <map>
#include <sstream>
#include <cmath>
#include <stddef.h>
#include <stdlib.h>
using namespace std;

//the definition of vector3d is copied from Loci Loci_Datatypes.h
  //---------------------vector3d------------------//
template <class T> 
struct vector3d {
  T x,y,z ;
  vector3d() {} 
  vector3d(T xx,T yy, T zz) : x(xx),y(yy),z(zz) {}
  vector3d(const vector3d &v) {x=v.x;y=v.y;z=v.z;}
} ;

template <class T> inline std::ostream & operator<<(std::ostream &s, const vector3d<T> &v)
{
  s << v.x << ' ' << v.y << ' ' << v.z << ' ' ;
  return s ;
}

template <class T> inline std::istream &operator>>(std::istream &s, vector3d<T> &v)
{
  s >> v.x >> v.y >> v.z ;
  return s ;
}

template <class T> inline T dot(const vector3d<T> &v1, const vector3d<T> &v2) {
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z ;
}

template <class T> inline T norm(const vector3d<T> &v) {
  return sqrt(v.x*v.x+v.y*v.y+v.z*v.z) ;
}

template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const vector3d<T> &v2) {
  return vector3d<T>(v1.y*v2.z-v1.z*v2.y,
                     v1.z*v2.x-v1.x*v2.z,
                     v1.x*v2.y-v1.y*v2.x) ;
}

template<class T> inline vector3d<T> cross(const vector3d<T> &v1, const T ra2[]) {
  return vector3d<T>(v1.y*ra2[2]-v1.z*ra2[1],
                     v1.z*ra2[0]-v1.x*ra2[2],
                     v1.x*ra2[1]-v1.y*ra2[0]) ;
}
template<class T> inline vector3d<T> cross(const T ra1[], const vector3d<T> &v2) {
  return vector3d<T>(ra1[1]*v2.z-ra1[2]*v2.y,
                     ra1[2]*v2.x-ra1[0]*v2.z,
                     ra1[0]*v2.y-ra1[1]*v2.x) ;
}

template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, float val) {
    target.x *= val ;
    target.y *= val ;
    target.z *= val ;
    return target ;
}

template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, float val) {
  target.x /= val ;
  target.y /= val ;
  target.z /= val ;
  return target ;
}

template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, double val) {
  target.x *= val ;
  target.y *= val ;
  target.z *= val ;
  return target ;
  }

template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, double val) {
  target.x /= val ;
  target.y /= val ;
  target.z /= val ;
  return target ;
}

template<class T> inline vector3d<T> &operator*=(vector3d<T> &target, long double val) {
  target.x *= val ;
  target.y *= val ;
  target.z *= val ;
  return target ;
}

template<class T> inline vector3d<T> &operator/=(vector3d<T> &target, long double val) {
  target.x /= val ;
  target.y /= val ;
  target.z /= val ;
  return target ;
}

template<class T> inline vector3d<T> operator+=(vector3d<T> &target, const vector3d<T> &val) {
  target.x += val.x ;
  target.y += val.y ;
  target.z += val.z ;
  return target ;
}

template<class T> inline vector3d<T> operator-=(vector3d<T> &target, const vector3d<T> &val) {
  target.x -= val.x ;
  target.y -= val.y ;
  target.z -= val.z ;
  return target ;
}

template<class T> inline vector3d<T> operator+(const vector3d<T> &v1, const vector3d<T> &v2) {
  return vector3d<T>(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z) ;
}

template<class T> inline vector3d<T> operator-(const vector3d<T> &v1, const vector3d<T> &v2) {
  return vector3d<T>(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z) ;
}

template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, float r2) {
  return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
}

template<class T> inline vector3d<T> operator*(float r1, const vector3d<T> &v2) {
  return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
}

template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, float r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
}

  template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, double r2) {
    return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
  }

template<class T> inline vector3d<T> operator*(double r1, const vector3d<T> &v2) {
    return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
}

template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, double r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
}

template<class T> inline vector3d<T> operator*(const vector3d<T> &v1, long double r2) {
  return vector3d<T>(v1.x*r2,v1.y*r2,v1.z*r2) ;
}

template<class T> inline vector3d<T> operator*(long double r1, const vector3d<T> &v2) {
  return vector3d<T>(v2.x*r1,v2.y*r1,v2.z*r1) ;
}

template<class T> inline vector3d<T> operator/(const vector3d<T> &v1, long double r2) {
    return vector3d<T>(v1.x/r2,v1.y/r2,v1.z/r2) ;
}


typedef vector3d<double>  vect3d;

//data structure for geometry
 struct GeomCoeff {
   
   vect3d b300;
   vect3d b030;
   vect3d b003;
   vect3d b210;
   vect3d b120;
   vect3d b021;
   vect3d b012;
   vect3d b102;
   vect3d b201;
   vect3d b111;
   vect3d N1,N2,N3 ;
   int pt1,pt2,pt3 ;
 } ;

struct GeomInfo {
  vect3d P1,P2,P3 ;
  vect3d N1,N2,N3 ;
} ;


typedef struct my_type2{
  //indexes of nodes 
  int p1; 
  int p2;
  int p3;
  //normals of nodes
  vect3d N1;
  vect3d N2;
  vect3d N3;
} triangle;

typedef struct my_type3{
  //indexes of nodes  
  int p1;
  int p2;
  int p3;
  int p4;
  //normals of nodes
  vect3d N1;
  vect3d N2;
  vect3d N3;
  vect3d N4;
} quad;


using std::vector ;
using std::string ;
using std::cerr ;
using std::endl ;
using std::ifstream ;
using std::ios ;
using std::cout;

int main(int argc, char *argv[]) {
  
  
  if(argc != 2) {
    cerr << "must pass in argument for grid to process"<< endl ;
    exit(-1) ;
  }
  string base = argv[1] ;
  
  string outfile = base +".geom";
  string surffile = base +".surf";
  string surfnfile = base +".surfn";

  //open .surf file and read in the positions of nodes
  //and face2node information
  ifstream inFile(surffile.c_str(),ios::in) ;
  if(inFile.fail()) {
    cerr<<"can't open " << surffile << endl ;
    exit(1);
  }
  int num_tris = 0;
  int num_quads = 0;
  int num_nodes = 0;
  inFile >> num_tris >> num_quads >> num_nodes;
  

  string ins; //a variable used to skip unwanted part of the input file 
  
  vect3d* pos = 0;
  triangle* tris = 0;
  quad* quads = 0;
  
  if(num_nodes <= 0){
    cerr << "num_nodes is zero" << endl;
    exit(1);
  }
  
  //read in positions of nodes
  pos = new vect3d[num_nodes];
  for(int i = 0; i < num_nodes; i++){
    inFile >> pos[i].x >> pos[i].y >> pos[i].z;
    getline(inFile, ins); // for optional items
  }
  
  //read in face2node of triangles
  if(num_tris > 0){
    tris = new triangle[num_tris];
    for(int i = 0; i < num_tris; i++){
      inFile >> tris[i].p1 >> tris[i].p2 >> tris[i].p3;
      getline(inFile, ins);
    }
  }

  
  //read in face2node information of quads
  if(num_quads > 0){
    quads = new quad[num_quads];
    for(int i = 0; i < num_quads; i++){
      inFile >> quads[i].p1 >> quads[i].p2 >> quads[i].p3>> quads[i].p4;
      getline(inFile, ins);
    }
  }
  inFile.close();

  //open .surfn file and read in normals for each boundary face
  inFile.open(surfnfile.c_str()) ;
  if(inFile.fail()) {
    cerr<<"can't open " << surfnfile << endl ;
    exit(1);
  }
  int n_tris, n_quads;
  inFile >> n_tris >> n_quads;
  if(num_tris != n_tris || num_quads != n_quads){
    cerr<<"num_tris or num_quads in .surf file and .surfn file are not consistent"<< endl;
    exit(1);
  }
  
  double tol = 1e-5;
  char c;
  //read in tri normals
  if(num_tris > 0){
    getline(inFile, ins);
    for(int i = 0; i < num_tris; i++){
      c = inFile.peek();
      if((c == '-') ||((c >= '0') && (c <='9')) ){
        inFile >> tris[i].N1;
        getline(inFile, ins);
      }
      else{
        tris[i].N1 = vect3d(0.0, 0.0, 0.0);
        getline(inFile, ins);
      }

      c = inFile.peek();
      if((c == '-') ||((c >= '0') && (c <='9')) ){
        inFile >> tris[i].N2;
        getline(inFile, ins);
      }
      else{
        tris[i].N2 = vect3d(0.0, 0.0, 0.0);
        getline(inFile, ins);
      }

      
      c = inFile.peek();
      if((c == '-') ||((c >= '0') && (c <='9')) ){
        inFile >> tris[i].N3;
        getline(inFile, ins);
      }
      else{
        tris[i].N3 = vect3d(0.0, 0.0, 0.0);
        getline(inFile, ins);
      }
    }
  }

  
  //read in quad normals
  if(num_quads > 0){
    for(int i = 0; i < num_quads; i++){
      c = inFile.peek();
      if((c == '-') ||((c >= '0') && (c <='9')) ){
        inFile >> quads[i].N1;
        getline(inFile, ins);
      }
      else{
        quads[i].N1 = vect3d(0.0, 0.0, 0.0);
        getline(inFile, ins);
      }

      c = inFile.peek();
      if((c == '-') ||((c >= '0') && (c <='9')) ){
        inFile >> quads[i].N2;
        getline(inFile, ins);
      }
      else{
        quads[i].N2 = vect3d(0.0, 0.0, 0.0);
        getline(inFile, ins);
      }

      
      c = inFile.peek();
      if((c == '-') ||((c >= '0') && (c <='9')) ){
        inFile >> quads[i].N3;
        getline(inFile, ins);
      }
      else{
        quads[i].N3 = vect3d(0.0, 0.0, 0.0);
        getline(inFile, ins);
      }

      c = inFile.peek();
      if((c == '-') ||((c >= '0') && (c <='9')) ){
        inFile >> quads[i].N4;
        getline(inFile, ins);
      }
      else{
        quads[i].N4 = vect3d(0.0, 0.0, 0.0);
        getline(inFile, ins);
        
      }
      
    }
  }
  

  //compute geomerty coefficient
  
  int num_geom = num_tris + 2*num_quads;
  
  
  GeomCoeff* geom = new GeomCoeff[num_geom];
  for(int i = 0; i < num_tris; i++){
    vect3d P1, P2, P3, N1, N2, N3;

    int pt1 = tris[i].p1-1 ;
    int pt2 = tris[i].p2-1 ;
    int pt3 = tris[i].p3-1 ;
    
      P1 = pos[pt1];
      P2 = pos[pt2];
      P3 = pos[pt3];

      
      N1 = tris[i].N1;
      N2 = tris[i].N2;
      N3 = tris[i].N3;
      if(norm(N1) < tol ||norm(N2) < tol ||norm(N3) < tol){
        vect3d total_vec = cross(P2-P1, P3-P1);
        double a = norm(total_vec);
        if(a < tol){
          cerr << "folded triangle" << endl;
          exit(1);
        }
        total_vec = total_vec / a;
        if(norm(N1) < tol) N1 = total_vec;
        if(norm(N2) < tol) N2 = total_vec;
        if(norm(N3) < tol) N3 = total_vec;
      }
      
      
      geom[i].b300 = P1;
      geom[i].b030 = P2;
      geom[i].b003 = P3;
      geom[i].N1=N1 ;
      geom[i].N2=N2 ;
      geom[i].N3=N3 ;
      geom[i].pt1 = pt1 ;
      geom[i].pt2 = pt2 ;
      geom[i].pt3 = pt3 ;
      
      double w12 = dot((P2-P1), N1);
      double w21 = dot((P1-P2), N2);
      double w23 = dot((P3-P2), N2);
      double w32 = dot((P2-P3), N3);
      double w31 = dot((P1-P3), N3);
      double w13 = dot((P3-P1), N1);
      geom[i].b210 = (2.0*P1 + P2 - w12*N1)/3.0;
      geom[i].b120 = (2.0*P2 + P1 - w21*N2)/3.0;
      geom[i].b021 = (2.0*P2 + P3 - w23*N2)/3.0;
      geom[i].b012 = (2.0*P3 + P2 - w32*N3)/3.0;
      geom[i].b102 = (2.0*P3 + P1 - w31*N3)/3.0;
      geom[i].b201 = (2.0*P1 + P3 - w13*N1)/3.0;
      vect3d E = (geom[i].b210 + geom[i].b120 + geom[i].b021 + geom[i].b012
                  + geom[i].b102 + geom[i].b201)/6.0;
      vect3d V = (P1 + P2 + P3)/3.0;
      geom[i].b111 = E + 0.5 *(E-V);
  }

  int geomID = num_tris;
  for(int i = 0; i < num_quads; i++){
    vect3d P1, P2, P3, N1, N2, N3;
    int pt1 = quads[i].p1-1 ;
    int pt2 = quads[i].p2-1 ;
    int pt3 = quads[i].p3-1 ;
  
      P1 = pos[pt1];
      P2 = pos[pt2];
      P3 = pos[pt3];
      N1 = quads[i].N1;
      N2 = quads[i].N2;
      N3 = quads[i].N3;

      if(norm(N1) < tol ||norm(N2) < tol ||norm(N3) < tol){
        vect3d total_vec = cross(P2-P1, P3-P1);
        double a = norm(total_vec);
        if(a < tol){
          cerr << "folded triangle" << endl;
          exit(1);
        }
        total_vec = total_vec / a;
        if(norm(N1) < tol) N1 = total_vec;
        if(norm(N2) < tol) N2 = total_vec;
        if(norm(N3) < tol) N3 = total_vec;
      }

      
      geom[geomID].b300 = P1;
      geom[geomID].b030 = P2;
      geom[geomID].b003 = P3;
      
      geom[geomID].N1 = N1;
      geom[geomID].N2 = N2;
      geom[geomID].N3 = N3;

      geom[geomID].pt1 = pt1;
      geom[geomID].pt2 = pt2;
      geom[geomID].pt3 = pt3;

      
      double w12 = dot((P2-P1), N1);
      double w21 = dot((P1-P2), N2);
      double w23 = dot((P3-P2), N2);
      double w32 = dot((P2-P3), N3);
      double w31 = dot((P1-P3), N3);
      double w13 = dot((P3-P1), N1);
      geom[geomID].b210 = (2.0*P1 + P2 - w12*N1)/3.0;
      geom[geomID].b120 = (2.0*P2 + P1 - w21*N2)/3.0;
      geom[geomID].b021 = (2.0*P2 + P3 - w23*N2)/3.0;
      geom[geomID].b012 = (2.0*P3 + P2 - w32*N3)/3.0;
      geom[geomID].b102 = (2.0*P3 + P1 - w31*N3)/3.0;
      geom[geomID].b201 = (2.0*P1 + P3 - w13*N1)/3.0;
      vect3d E = (geom[geomID].b210 + geom[geomID].b120 + geom[geomID].b021 + geom[geomID].b012
                  + geom[geomID].b102 + geom[geomID].b201)/6.0;
      vect3d V = (P1 + P2 + P3)/3.0;
      geom[geomID].b111 = E + 0.5 *(E-V);

      //compute next triangle P1 P3 P4
      geomID++;
      pt1 = quads[i].p1-1 ;
      pt2 = quads[i].p3-1 ;
      pt3 = quads[i].p4-1 ;
      P1 = pos[pt1];
      P2 = pos[pt2];
      P3 = pos[pt3];
      N1 = quads[i].N1;
      N2 = quads[i].N3;
      N3 = quads[i].N4; 


      if(norm(N1) < tol ||norm(N2) < tol ||norm(N3) < tol){
        vect3d total_vec = cross(P2-P1, P3-P1);
        double a = norm(total_vec);
        if(a < tol){
          cerr << "folded triangle" << endl;
          exit(1);
        }
        total_vec = total_vec / a;
        if(norm(N1) < tol) N1 = total_vec;
        if(norm(N2) < tol) N2 = total_vec;
        if(norm(N3) < tol) N3 = total_vec;
      }
      
      geom[geomID].b300 = P1;
      geom[geomID].b030 = P2;
      geom[geomID].b003 = P3;
      geom[geomID].N1 = N1;
      geom[geomID].N2 = N2;
      geom[geomID].N3 = N3;
      geom[geomID].pt1 = pt1;
      geom[geomID].pt2 = pt2;
      geom[geomID].pt3 = pt3;
      
       w12 = dot((P2-P1), N1);
       w21 = dot((P1-P2), N2);
       w23 = dot((P3-P2), N2);
       w32 = dot((P2-P3), N3);
       w31 = dot((P1-P3), N3);
       w13 = dot((P3-P1), N1);
      geom[geomID].b210 = (2.0*P1 + P2 - w12*N1)/3.0;
      geom[geomID].b120 = (2.0*P2 + P1 - w21*N2)/3.0;
      geom[geomID].b021 = (2.0*P2 + P3 - w23*N2)/3.0;
      geom[geomID].b012 = (2.0*P3 + P2 - w32*N3)/3.0;
      geom[geomID].b102 = (2.0*P3 + P1 - w31*N3)/3.0;
      geom[geomID].b201 = (2.0*P1 + P3 - w13*N1)/3.0;
       E = (geom[geomID].b210 + geom[geomID].b120 + geom[geomID].b021 + geom[geomID].b012
                  + geom[geomID].b102 + geom[geomID].b201)/6.0;
       V = (P1 + P2 + P3)/3.0;
       geom[geomID].b111 = E + 0.5 *(E-V);

       geomID++;
      
  }

  //smooth up the edges, 
  map<pair<int, int>, vector<pair<vect3d, vect3d> > > edge_map;
  map<pair<int,int>,vector<pair<int,int> > > edge_info ;
  //for each triangle, put the control points of each edge into a map
  for(int i = 0; i < num_geom; i++){
    int dir = (geom[i].pt1 < geom[i].pt2)?0:1 ;
    pair<int,int> e1 = make_pair(min(geom[i].pt1,geom[i].pt2),
                               max(geom[i].pt1,geom[i].pt2)) ;
      
    edge_map[e1].push_back(make_pair(geom[i].b210, geom[i].b120));
    edge_info[e1].push_back(make_pair(i,dir+0)) ;

    dir = (geom[i].pt2 < geom[i].pt3)?0:1 ;
    pair<int,int> e2 = make_pair(min(geom[i].pt2,geom[i].pt3),
                                 max(geom[i].pt2,geom[i].pt3)) ;

    edge_map[e2].push_back(make_pair(geom[i].b021, geom[i].b012));
    edge_info[e2].push_back(make_pair(i,dir+2)) ;
    
    dir = (geom[i].pt1 < geom[i].pt3)?0:1 ;
    pair<int,int> e3 = make_pair(min(geom[i].pt1,geom[i].pt3),
                                 max(geom[i].pt1,geom[i].pt3)) ;

    edge_map[e3].push_back(make_pair(geom[i].b102, geom[i].b201));
    edge_info[e3].push_back(make_pair(i,dir+4)) ;
  }

  //compute the mean value of control points for each edge in the map
  map<pair<int, int>, vector<pair<vect3d, vect3d> > >::iterator empt;
  
  for(empt= edge_map.begin(); empt != edge_map.end(); empt++){
    if(empt->second.size() < 1 || empt->second.size() > 2){
      cerr << "something is wrong in edge smoothing "<< empt->second.size()
           <<" " << empt->first.first <<" " << empt->first.second << endl;
      for(unsigned int j = 0; j < empt->second.size(); j++){
        cerr << empt->second[j].first << "   " <<empt->second[j].second << endl;
      }
      exit(0);
    }
    if(empt->second.size() == 2){
      vect3d N0,N1,X0,X1,P0,P1;
      int x0id ;

      pair<int,int> pdat = edge_info[empt->first][0] ;
      int gid0 = pdat.first ;
      int mode0 = pdat.second ;
      pdat = edge_info[empt->first][1] ;
      int gid1 = pdat.first ;
      int mode1 = pdat.second ;
      if((mode0&6) == 0) {
        N0 = geom[gid0].N1 ;
        X0 = geom[gid0].b300 ;
        x0id = geom[gid0].pt1 ;
        P0 = geom[gid0].b210 ;
      } else if((mode0&6) == 2) {
        N0 = geom[gid0].N2 ;
        X0 = geom[gid0].b030 ;
        x0id = geom[gid0].pt2 ;
        P0 = geom[gid0].b021 ;
      } else {
        N0 = geom[gid0].N3 ;
        X0 = geom[gid0].b003 ;
        x0id = geom[gid0].pt3 ;
        P0 = geom[gid0].b102 ;
      }
      if((mode1&6) == 0) {
        if(x0id == geom[gid1].pt1) {
          N1 = geom[gid1].N1 ;
          X1 = geom[gid1].b300 ;
          P1 = geom[gid1].b210 ;
        } else {
          N1 = geom[gid1].N2 ;
          X1 = geom[gid1].b030 ;
          P1 = geom[gid1].b120 ;
        }
      } else if((mode1&6) == 2) {
        if(x0id == geom[gid1].pt2) {
          N1 = geom[gid1].N2 ;
          X1 = geom[gid1].b030 ;
          P1 = geom[gid1].b021 ;
        } else {
          N1 = geom[gid1].N3 ;
          X1 = geom[gid1].b003 ;
          P1 = geom[gid1].b012 ;
        }
      } else {
        if(x0id == geom[gid1].pt3) {
          N1 = geom[gid1].N3 ;
          X1 = geom[gid1].b003 ;
          P1 = geom[gid1].b102 ;
        } else {
          N1 = geom[gid1].N1 ;
          X1 = geom[gid1].b300 ;
          P1 = geom[gid1].b201 ;
        }
      }
      if(X0.x != X1.x || X0.y != X1.y || X0.z != X1.z) {
        cerr << "coordinates not equal1!" << endl ;
        cerr << "X0=" << X0 << " X1=" << X1 << endl ;
        cerr << "mode0 = " << mode0 << ", mode1=" << mode1 << endl ;
        cerr << "g0=" << geom[gid0].pt1 <<' '<< geom[gid0].pt2 << ' '
             << geom[gid0].pt3
             << "g1=" << geom[gid1].pt1 <<' '<< geom[gid1].pt2 << ' '
             << geom[gid1].pt3  <<endl ;
      }
      
      vect3d NC = cross(N0,N1) ;
      double nNc = norm(NC) ;
      vect3d PM = .5*(P0+P1) ;
      if(nNc > 4e-2) {
        vect3d dV = PM-X0 ;
        double ldV = norm(dV) ;
        NC *= ldV/nNc ; // scale Nc to have length of dV
        if(dot(dV,NC) < 0) // make Nc point in the same direction as DV 
          NC *= -1. ;
        PM = X0+NC ;
      }
      // put PM back
      if((mode0&6) == 0) {
        geom[gid0].b210 = PM ;
      } else if((mode0&6) == 2) {
        geom[gid0].b021 = PM ;
      } else {
        geom[gid0].b102 = PM ;
      }
      if((mode1&6) == 0) {
        if(x0id == geom[gid1].pt1) {
          geom[gid1].b210 = PM ;
        } else {
          geom[gid1].b120 = PM ;
        }
      } else if((mode1&6) == 2) {
        if(x0id == geom[gid1].pt2) {
          geom[gid1].b021 = PM ;
        } else {
          geom[gid1].b012 = PM ;
        }
      } else {
        if(x0id == geom[gid1].pt3) {
          geom[gid1].b102  = PM ;
        } else {
          geom[gid1].b201 = PM ;
        }
      }

      // Now work on other side of edge
      
      if((mode0&6) == 0) {
        N0 = geom[gid0].N2 ;
        X0 = geom[gid0].b030 ;
        P0 = geom[gid0].b120 ;
      } else if((mode0&6) == 2) {
        N0 = geom[gid0].N3 ;
        X0 = geom[gid0].b003 ;
        P0 = geom[gid0].b012 ;
      } else {
        N0 = geom[gid0].N1 ;
        X0 = geom[gid0].b300 ;
        P0 = geom[gid0].b201 ;
      }
      if((mode1&6) == 0) {
        if(x0id != geom[gid1].pt1) {
          N1 = geom[gid1].N1 ;
          X1 = geom[gid1].b300 ;
          P1 = geom[gid1].b210 ;
        } else {
          N1 = geom[gid1].N2 ;
          X1 = geom[gid1].b030 ;
          P1 = geom[gid1].b120 ;
        }
      } else if((mode1&6) == 2) {
        if(x0id != geom[gid1].pt2) {
          N1 = geom[gid1].N2 ;
          X1 = geom[gid1].b030 ;
          P1 = geom[gid1].b021 ;
        } else {
          N1 = geom[gid1].N3 ;
          X1 = geom[gid1].b003 ;
          P1 = geom[gid1].b012 ;
        }
      } else {
        if(x0id != geom[gid1].pt3) {
          N1 = geom[gid1].N3 ;
          X1 = geom[gid1].b003 ;
          P1 = geom[gid1].b102 ;
        } else {
          N1 = geom[gid1].N1 ;
          X1 = geom[gid1].b300 ;
          P1 = geom[gid1].b201 ;
        }
      }
      if(X0.x != X1.x || X0.y != X1.y || X0.z != X1.z) {
        cerr << "coordinates not equal2!" << endl ;
        cerr << "X0=" << X0 << " X1=" << X1 << endl ;
        cerr << "mode0 = " << mode0 << ", mode1=" << mode1 << endl ;
      }
      NC = cross(N0,N1) ;
      nNc = norm(NC) ;
      PM = .5*(P0+P1) ;
      if(nNc > 4e-2) {
        vect3d dV = PM-X0 ;
        double ldV = norm(dV) ;
        NC *= ldV/nNc ; // scale Nc to have length of dV
        if(dot(dV,NC) < 0) // make Nc point in the same direction as DV 
          NC *= -1. ;
        PM = X0+NC ;
      }
      // put PM back
      if((mode0&6) == 0) {
        geom[gid0].b120 = PM ;
      } else if((mode0&6) == 2) {
        geom[gid0].b012 = PM;
      } else {
        geom[gid0].b201 = PM ;
      }
      if((mode1&6) == 0) {
        if(x0id != geom[gid1].pt1) {
          geom[gid1].b210 = PM ;
        } else {
          geom[gid1].b120 = PM ;
        }
      } else if((mode1&6) == 2) {
        if(x0id != geom[gid1].pt2) {
          geom[gid1].b021 = PM ;
        } else {
          geom[gid1].b012 = PM ;
        }
      } else {
        if(x0id != geom[gid1].pt3) {
          geom[gid1].b102  = PM ;
        } else {
          geom[gid1].b201 = PM ;
        }
      }
    }
  }
  
  //output the geometry into .coeff file 
  ofstream ofile(outfile.c_str()); 

  ofile.precision(14) ;
  //write out positions  of nodes
  ofile << num_nodes << endl;
  for(int i = 0; i < num_nodes; i++){
    ofile << pos[i].x<< " "  << pos[i].y<< " " <<pos[i].z<< endl;
  }
  
  //write out face2node
  ofile << num_geom << endl;
  for(int i = 0; i < num_tris; i++){
    ofile << tris[i].p1-1 << ' '
          << tris[i].p2-1 << ' '
          << tris[i].p3-1 << endl;
  }
  for(int i = 0; i < num_quads; i++){
    ofile << quads[i].p1-1 << ' '
          << quads[i].p2-1 << ' '
          << quads[i].p3-1 <<endl;
    ofile  << quads[i].p1-1 << ' '
          << quads[i].p3-1 << ' '
          << quads[i].p4-1 <<endl;
  }

  //write out GeomCoeff
  ofile.precision(10);
  for(int i = 0; i < num_geom; i++){
    GeomCoeff g = geom[i];
      ofile << g.b300<<' '
            << g.b030<<' '
            << g.b003<<' '
            << g.b210<<' '
            << g.b120<<' '
            << g.b021<<' '
            << g.b012<<' '
            << g.b102<<' '
            << g.b201<<' '
            << g.b111<<endl;
  }
  ofile.close();
  delete [] pos;
  if(tris != 0)delete [] tris;
  if(quads != 0) delete [] quads;
  if(geom != 0) delete [] geom;
}
