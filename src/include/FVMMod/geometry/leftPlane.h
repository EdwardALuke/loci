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
#ifndef LEFTPLANE_H
#define LEFTPLANE_H
#include "geomTest.h"

/// @brief Class that creates a bounding plane for a given point and a normal
///        vector. All cells that are opposite of the normal vector are
///        considered inside the plane.
class leftPlane : public geomTest
{
  vector3d<real_t> pt; //!< [m] Reference point to begin plane.
  vector3d<real_t> n;  //!< [-] Normal vector to the reference point.
public:

  leftPlane()
  {
    pt = vector3d<real_t>(0.0,0.0,0.0);
    n = vector3d<real_t>(1.0,0.0,0.0);
  }

  /** **************************************************************************
   * @brief Construct the leftPlane object with a 'point' and 'normal' vector.
   * @param p1 [m] Reference point to begin plane.
   * @param n1 Normal vector to the reference point.
   ****************************************************************************/
  leftPlane(vector3d<real_t> p1, vector3d<real_t> n1) : pt(p1), n(n1)
  {
    n *= 1./norm(n);
  } // End of leftPlane()

  /** **************************************************************************
   * @brief Construct the leftPlane object based on user input options_list.
   * @param ol Referenced options_list pointer
   ****************************************************************************/
  leftPlane(const options_list &ol)
  {
    if(!ol.optionExists("point") || !ol.optionExists("normal"))
    {
      printf("ERROR::leftPlane: Need 'point' and 'normal'\n");
      Loci::Abort();
    }
    ol.getOptionUnits("point",  "m", pt);
    ol.getOptionUnits("normal", "",  n);
    n *= 1./max<real_t>(norm(n),  1e-33);
  } // End of leftPlane()

  /** **************************************************************************
   * @brief  Check if point 'p' is to the left (negative direction) of the plane
   *         described by the point 'pt' and normal vector 'n'.
   * @param  p  [m] Position vector of a point.
   * @return true  - If 'p' is inside plane.
   * @return false - If 'p' is outside plane.
   ****************************************************************************/
  bool inGeomPt(vector3d<real_t> p) const
  {
    if(dot(p-pt,n) < 0.0)
      return true;
    return false;
  } // End of inGeomPt()

  real_t distToSurface(vector3d<real_t> p) const
  {
    return dot(p-pt,n);
  }
}; // End Class leftPlane

#endif // LEFTPLANE_H
