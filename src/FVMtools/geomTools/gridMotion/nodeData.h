//#############################################################################
//#
//# Copyright 2014-2025, Mississippi State University
//#
//# The GridMover is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The GridMover software is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the GridMover software. If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef NODEDATA_H
#define NODEDATA_H

#include "rotor.h"

namespace gridMotion {

  struct nodeDataBase {
    vector3d<float> pos ;
    float weight ;
    int order ;
  } ;

  struct NodeData {
    vector3d<double> pos;
    vector3d<float> disp;
    Rotor            rot;
    float           weight ;
  };

  inline std::ostream & operator<<(std::ostream & lhs, const NodeData & rhs) {
    lhs << rhs.pos       << std::endl;
    lhs << rhs.disp      << std::endl;
    lhs << rhs.rot.alpha << std::endl;
    lhs << rhs.rot.beta  << std::endl;
    lhs << rhs.weight    << std::endl;
    return lhs;
  }

  inline std::istream & operator>>(std::istream & lhs, NodeData & rhs) {
    lhs >> rhs.pos >> rhs.disp >> rhs.rot.alpha >> rhs.rot.beta >> rhs.weight;
    return lhs;
  }

  template<class T>
  struct Append {
    void operator()(T & lhs, const T & rhs) {
      lhs.insert(lhs.end(),rhs.begin(),rhs.end());
    }
  };
  
  
#define C_PLANE    (1)
#define C_DISC     (2)
#define C_CYLINDER (3)
  
  struct constraintData {
    int constraintType ;
    vector3d<double> v1,v2 ;
    double r,weight ;
    int id ;
  } ;

  inline std::ostream & operator<<(std::ostream & lhs, const constraintData & rhs) {
    lhs << rhs.constraintType      << std::endl;
    lhs << rhs.v1      << std::endl;
    lhs << rhs.v2 << std::endl;
    lhs << rhs.r  << std::endl;
    lhs << rhs.weight  << std::endl;
    lhs << rhs.id  << std::endl;
    return lhs;
  }

  inline std::istream & operator>>(std::istream & lhs, constraintData & rhs) {
    lhs >> rhs.constraintType >> rhs.v1 >> rhs.v2 >> rhs.r >> rhs.weight >> rhs.id;
    return lhs;
  }
}

#endif // NODEDATA_H
