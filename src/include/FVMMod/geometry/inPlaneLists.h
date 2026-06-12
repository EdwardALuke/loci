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
#ifndef INPLANELIST_H
#define INPLANELIST_H
#include "geomTest.h"

/// @brief Class that creates a bounding cone for two given points and two radii.
class inPlaneList : public geomTest {
    vector<vector3d<real_t> > pts ;
    vector<vector3d<real_t> > normals ;
  public:
    planelist_type() {
    }

    planelist_type(std::vector<vector3d<real_t> > &ps,
                   std::vector<vector3d<real_t> > &ns) : pts(ps),normals(ns) {}


    bool inGeomPt(vector3d<real_t> pt) const {
      bool ingeometry = true ;
      for(size_t i=0;i<pts.size();++i) {
      if(dot(pt-pts[i],normals[i]) > 0)
        ingeometry = false ;
      }
      return ingeometry ;
    }
        
    real_t distToSurface(vector3d<real_t> pt) const {
      real_t maxdist = 0 ;
      for(size_t i=0;i<pts.size();++i) {
        maxdist = max(maxdist,dot(pt-pts[i],normals[i])) ;
      }
      return max<real_t>(maxdist,0.0) ;
    }
};

#endif // INPLANELISTS_H
