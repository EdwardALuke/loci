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
#ifndef INSPHERE_H
#define INSPHERE_H
#include "geomTest.h"

/// @brief Class that creates a bounding sphere for a given center point and a
///        radius.
class inSphere : public geomTest
{
  real_t   r;      //!< [m]
  real_t   r2;     //!< [m^2]
  vector3d<real_t> center; //!< [m]
public:

  inSphere()
  {
    r1 = 1.0; 
    r2 = 1.0;
    center = vector3d<real_t>(0.0,0.0,0.0);
  }


  /** **************************************************************************
   * @brief
   * @param ri
   * @param c
   ****************************************************************************/
  inSphere(real_t ri, vector3d<real_t> c)
  {
    r      = ri;
    center = c;
    r2     = r*r;
  } // End of inSphere()

  /** **************************************************************************
   * @brief Construct the inSphere object based on user input options_list.
   * @param ol Referenced options_list pointer
   ****************************************************************************/
  inSphere(const options_list &ol)
  {
    // Here get sphere information
    if(!ol.optionExists("radius") || !ol.optionExists("center"))
    {
      printf("ERROR::inSphere: Need a radius and center\n");
      Loci::Abort();
    }
    ol.getOptionUnits("radius","m",r);
    ol.getOptionUnits("center","m",center);
    r2 = r*r;
  } // End of inSphere()

  /** **************************************************************************
   * @brief  Determine if a position vector 'p' is within the bounding sphere.
   * @param  p  [m] Position vector of a point.
   * @return true  - If 'p' is inside bounding sphere.
   * @return false - If 'p' is outside bounding sphere.
   ****************************************************************************/
  bool inGeomPt(vector3d<real_t> p) const
  {
    if(dot(p-center, p-center) < r2)
      return true;
    return false;
  } // End of inGeomPt()

  real_t distToSurface(vector3d<real_t> pt) const {
    vector3d<real_t> p = pt-center ;
    real_t r2 = dot(p,p) ;
    if(r2 < radius2) {
      return 0 ;
    }
    return sqrt(r2)-radius ;
  }

}; // End Class inSphere
#endif // INSPHERE_H
