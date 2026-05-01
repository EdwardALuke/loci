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
#ifndef INBOX_H
#define INBOX_H
#include "geomTest.h"

/// @brief Class that creates a bounding box for two given points.
class inBox : public geomTest
{
  real_t xmax; //!< [m]
  real_t ymax; //!< [m]
  real_t zmax; //!< [m]

  real_t xmin; //!< [m]
  real_t ymin; //!< [m]
  real_t zmin; //!< [m]
public:

  inBox(){
    xmax = 1.0;
    xmin = 0.0;
    ymax = 1.0;
    ymin = 0.0;
    zmax = 1.0;
    zmin = 0.0;
  }

  /** **************************************************************************
   * @brief Construct the inBox object given two points.
   * @param p1 Position vector of first point
   * @param p2 Positoin vector of second point
   ****************************************************************************/
  inBox(vector3d<real_t> p1, vector3d<real_t> p2)
  {
    xmax = max(p1.x, p2.x);
    xmin = min(p1.x, p2.x);

    ymax = max(p1.y, p2.y);
    ymin = min(p1.y, p2.y);

    zmax = max(p1.z, p2.z);
    zmin = min(p1.z, p2.z);
  } // End of inBox()

  /** **************************************************************************
   * @brief Construct the inBox object based on user input options_list.
   * @param ol Referenced options_list pointer
   ****************************************************************************/
  inBox(const options_list &ol)
  {
    if(!ol.optionExists("p1") || !ol.optionExists("p2"))
    {
      printf("ERROR::inBox: Need two points, 'p1' and 'p2'\n");
      Loci::Abort();
    }

    vector3d<real_t> p1;
    vector3d<real_t> p2;
    ol.getOptionUnits("p1", "m", p1);
    ol.getOptionUnits("p2", "m", p2);
    xmax = max(p1.x, p2.x);
    xmin = min(p1.x, p2.x);

    ymax = max(p1.y, p2.y);
    ymin = min(p1.y, p2.y);

    zmax = max(p1.z, p2.z);
    zmin = min(p1.z, p2.z);
  } // End of inBox()

  /** **************************************************************************
   * @brief  Determine if a position vector 'p' is within the bounding box.
   * @param  p  [m] Position vector of a point.
   * @return true  - If 'p' is inside bounding box.
   * @return false - If 'p' is outside bounding box.
   ****************************************************************************/
  bool inGeomPt(vector3d<real_t> p) const
  {
    if(p.x >= xmin && p.x <= xmax &&
       p.y >= ymin && p.y <= ymax &&
       p.z >= zmin && p.z <= zmax)
      return true;
    return false;
  } // End of inGeomPt()

  real_t distToSurface(vector3d<real_t> p) const
  {
    vector3d<real_t> c = 0.5*vector3d<real_t>(xmin+xmax,ymin+ymax,zmin+zmax);
    vector3d<real_t> h = 0.5*vector3d<real_t>(xmax-xmin,ymax-ymin,zmax-zmin);
    real_t dx = abs(p.x-c.x) - h.x;
    real_t dy = abs(p.y-c.y) - h.y;
    real_t dz = abs(p.z-c.z) - h.z;

    real_t ex = max(dx,0.0);
    real_t ey = max(dy,0.0);
    real_t ez = max(dz,0.0);
    // external distance
    real_t ed = sqrt(ex*ex + ey*ey + ez*ez);

    // internal distance
    real_t id = min(max({dx,dy,dz}),0.0);

    return ed + id; // distance to inside/outside
  }

}; // End Class inBox

#endif // INBOX_H
