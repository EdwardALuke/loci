//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
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
#ifndef HEX_DEFINES_H
#define HEX_DEFINES_H
#include <Loci.h>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include "defines.h"

/// Some machine, long int and int have the same size, both are 32 bits
typedef Loci::vector3d<int64> IntegerPoint;
typedef Loci::vector2d<int64>  Point2d;

class HexCell;
class Node;
class Edge;
class QuadFace;

/// The direction of neighbor cells, also as faceID
/// RIGHT: x = 1; LEFT: x = -1
/// FRONT: y = 1; BACK: y = -1
/// UP:    z = 1; DOWN: z = -1
enum DIRECTION{RIGHT, LEFT, FRONT, BACK, UP, DOWN};

/// The direction of normal of faces
enum NORMAL_DIRECTION{XX, YY, ZZ};


namespace Loci {

  template<class T> inline bool operator < (const vector3d<T> &p1,const vector3d<T> &p2){
    if(p1.x < p2.x) return true;

    else if(p1.x == p2.x){
      if(p1.y < p2.y) return true;
      else if (p1.y == p2.y){
        if(p1.z < p2.z) return true;
      }
    }
    return false;
  }

  template<class T> inline bool operator == (const  vector3d<T>& p1,const vector3d<T>& p2){
    return (p1.x == p2.x) && (p1.y == p2.y) && (p1.z == p2.z);
  }


  /// Function object for building a map such as yNodes
  class yless{
  public:
    template<class T> inline  bool operator()(const vector3d<T>& p1, const vector3d<T>& p2)const{
      if(p1.x < p2.x) return true;
      else if(p1.x == p2.x){
        if(p1.z < p2.z) return true;
        else if (p1.z == p2.z){
          if(p1.y < p2.y) return true;
        }
      }
      return false;
    }
  };

  /// Function object for building a map such as xNodes
  class xless{
  public:
    template<class T>  inline  bool operator()(const vector3d<T>& p1, const vector3d<T>& p2)const{
      if(p1.y < p2.y) return true;
      else if(p1.y == p2.y){
        if(p1.z < p2.z) return true;
        else if (p1.z == p2.z){
          if(p1.x < p2.x) return true;
        }
      }
      return false;
    }
  };
}

struct Edge2d{
  int64 pos;
  int64 head;
  int64 tail;
  Edge2d():pos(0),head(0),tail(0){}
  Edge2d(int64 p,int64 h, int64 t):pos(p), head(h),tail(t){}
  bool operator==(const Edge2d & ref) { return pos==ref.pos && head==ref.head && tail == ref.tail ; }
  bool operator!=(const Edge2d & ref) { return !operator==(ref) ; }
};

struct Range2d{
  Point2d minP;
  Point2d maxP;
  Range2d():minP(Point2d(0, 0)), maxP(Point2d(0, 0)){}
  Range2d(Point2d p0, Point2d p1):minP(p0), maxP(p1){}
};

inline bool operator < (const Edge2d& e1,const Edge2d& e2){
  if(e1.pos < e2.pos) return true;

  else if(e1.pos == e2.pos){
    if(e1.head < e2.head) return true;
    else if (e1.head == e2.head){
      if(e1.tail < e2.tail)return true;
    }

  }
  return false;
}

inline bool operator == (const  Edge2d& e1,const Edge2d& e2){
  return ((e1.pos == e2.pos) && (e1.head == e2.head)&&
          (e1.tail == e2.tail));
}

//transfer from face's local coordinates to cell's local coordinates
Point2d transfer_f2c(Point2d p,Point2d maxPc, char orientCode);
Point2d transfer_c2f(Point2d p, Point2d maxPf, char orientCode);
Range2d transfer_f2c(Range2d c, Point2d maxPc, char orientCode);

#endif
