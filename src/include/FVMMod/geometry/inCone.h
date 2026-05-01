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
#ifndef INCONE_H
#define INCONE_H
#include "geomTest.h"

/// @brief Class that creates a bounding cone for two given points and two radii.
class inCone : public geomTest
{
  real_t   mag; //!< [m]
  real_t   r1;  //!< [m]
  real_t   r2;  //!< [m]
  vector3d<real_t> p1;  //!< [m]
  vector3d<real_t> p2;  //!< [m]
  vector3d<real_t> n;   //!< [-]
public:

  inCone()
  {
    p1 = vector3d<real>(0.0,0.0,0.0);
    p2 = p2;
    r1 = 1.0; 
    r2 = 1.0;
    mag = 1.0;
    n = p1;
  }

  /** **************************************************************************
   * @brief
   * @param r1i
   * @param r2i
   * @param p1i
   * @param p2i
   ****************************************************************************/
  inCone(real_t r1i, real_t r2i, vector3d<real_t> p1i, vector3d<real_t> p2i)
  {
    r1   = r1i;
    r2   = r2i;
    p1   = p1i;
    p2   = p2i;
    n    = p2-p1;
    mag  = norm(n);
    n   *= 1./mag;
  } // End of inCone()

  /** **************************************************************************
   * @brief Construct the inCone object based on user input options_list.
   * @param ol Referenced options_list pointer
   ****************************************************************************/
  inCone(const options_list &ol)
  {
    // Here get sphere information
    if(!ol.optionExists("r1") || !ol.optionExists("r2") ||
       !ol.optionExists("p1") || !ol.optionExists("p2"))
    {
      printf("ERROR::inCone: Needs two radii (r1 and r2) and two axis points "
             "(p1 and p2)\n");
      Loci::Abort();
    }

    ol.getOptionUnits("r1","m",r1);
    ol.getOptionUnits("r2","m",r2);
    ol.getOptionUnits("p1","m",p1);
    ol.getOptionUnits("p2","m",p2);
    n    = p2-p1;
    mag  = norm(n);
    n   *= 1./mag;
  } // End of inCone()

  /** **************************************************************************
   * @brief  Determine if a position vector 'p' is within the bounding cone.
   * @param  p  [m] Position vector of a point.
   * @return true  - If 'p' is inside bounding cone.
   * @return false - If 'p' is outside bounding cone.
   ****************************************************************************/
  bool inGeomPt(vector3d<real_t> p) const
  {
    real_t v = dot(p-p1,n);
    if(v < 0.0 || v > mag)
      return false;

    // project point onto axis
    real_t   t     = v/mag;
    real_t   rx    = r1*(1.-t) + r2*t;
    real_t   rx2   = rx*rx;
    vector3d<real_t> paxis = p1 + n*v;
    if(dot(p-paxis, p-paxis) > rx2) // Check radius
      return false; // outside
    return true;
  } // End of inGeomPt()

  real_t distToSurface(vector3d<real_t> pt) const {
    vector3d<real_t> v = pt-pt2 ;
    vector3d<real_t> h = pt2-pt1;
    real_t t = dot(p,h)*rlen2 ;
    real_t hn = norm(h);
    real_t rad = sqrt(radius*radius + hn*hn);
    real_t orad = rad;
    real_t sint = radius*orad;
    real_t cost = hn*orad;

    real_t dperp = norm(v - (t * h));
    real_t d2side = dperp * cost - t * sint;
    real_t daside = dperp * sint + t * cost;
    real_t dbase = hn - t;
    // outside the body
    if (t < 0) { // below the base
        return sqrt(dperp*dperp + t*t);
    }  else if (dperp*hn > t*radius) { // outside
      if (d2side < 0) return sqrt(dperp*dperp + t*t);
      if (daside > rad) return sqrt(pow(dperp-radius,2) + pow(t-hn,2));
    } else { // inside
      return max(d2side,t-h);
    }
  }

}; // End Class inCone

#endif // INCONE_H
