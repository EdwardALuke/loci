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
#ifndef NODETREE_H
#define NODETREE_H

#include "gridTypes.h"

namespace gridMotion {

  struct treeInfoBase {
    int left, right, start, end, coord ;
    float split ;
  } ;
  
  struct treeInfoDisplacement {
    Array<vector3d<float>,4 > displacement ;
  } ;
  struct treeInfoRotation {
    Array<tensor3d<float>,4 > rotation ;
    float drot ;
  } ;
      
  struct tree_info {
    // approximate displacement field values for 4 points
    Array<vector3d<double>,4> displacement ;
    Array<tensor3d<double>,4> rotation ;
    int left,right ; // left and right tree branches
    int start,end ; // start and end pos in split_list
    int coord ; // Coordinate to split
    float split ; // value to split
    Array<vector3d<float>,4> q ; // quad point approximation
    vector3d<float> centroid ; // centroid of collection
    float weight ;
    // error estimate parameters
    float radius ;
    float err ;
    float drot ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, const tree_info &ti) {
    s << ti.left << " ";
    s << ti.right << " ";
    s << ti.start << " ";
    s << ti.end << " ";
    s << ti.coord << " ";
    s << ti.split << " ";
    s << ti.q[0] << " ";
    s << ti.q[1] << " ";
    s << ti.q[2] << " ";
    s << ti.q[3] << " ";
    s << ti.centroid << " ";
    s << ti.weight << " ";
    s << ti.displacement[0] << " ";
    s << ti.displacement[1] << " ";
    s << ti.displacement[2] << " ";
    s << ti.displacement[3] << " ";
    s << ti.rotation[0] << " ";
    s << ti.rotation[1] << " ";
    s << ti.rotation[2] << " ";
    s << ti.rotation[3] << " ";
    s << ti.radius << " ";
    s << ti.err << " ";
    s << ti.drot;
    return s ;
  }
  inline std::istream & operator>>(std::istream &s, tree_info &ti) {
    s >> ti.left ;
    s >> ti.right ;
    s >> ti.start ;
    s >> ti.end ;
    s >> ti.coord ;
    s >> ti.split ;
    s >> ti.q[0] ;
    s >> ti.q[1] ;
    s >> ti.q[2] ;
    s >> ti.q[3] ;
    s >> ti.centroid ;
    s >> ti.weight ;
    s >> ti.displacement[0] ;
    s >> ti.displacement[1] ;
    s >> ti.displacement[2] ;
    s >> ti.displacement[3] ;
    s >> ti.rotation[0] ;
    s >> ti.rotation[1] ;
    s >> ti.rotation[2] ;
    s >> ti.rotation[3] ;
    s >> ti.radius ;
    s >> ti.err ;
    s >> ti.drot;
    return s ;
  }
}


#endif
