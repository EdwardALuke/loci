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
#ifndef NEAR_RECTANGLE_H
#define NEAR_RECTANGLE_H
#include "geomTest.h"
#include "mathUtils.h"
using namespace std;

/// @brief Class that 
class nearRectangle : public geomTest
{
  vector3d<real_t>  p1; //!< [m] Position vector of 1st corner of rectangle
  vector3d<real_t>  p2; //!< [m] Position vecotr of 2nd, opposite corner of rectangle
  vector3d<real_t>  pC; //!< [m] Position vecotr of center of the rectangle
  vector3d<real_t>  p3; //!< [m] Position vecotr of 2nd, opposite corner of rectangle
  vector3d<real_t>  p4; //!< [m] Position vecotr of 2nd, opposite corner of rectangle
  vector3d<real_t>  n;  //!< [-] Unit normal vector between p1 and p2

  real_t xmax; //!< [m]
  real_t ymax; //!< [m]
  real_t zmax; //!< [m]

  real_t xmin; //!< [m]
  real_t ymin; //!< [m]
  real_t zmin; //!< [m]
  vector<vector3d<real_t>> vector1 = vector<vector3d<real_t>>(9);


public:
  real_t    totalArea; //!< [m^2]

  /** **************************************************************************
   * @brief 
   * @param vecP1 
   * @param vecP2 
   * @param vecN 
   ****************************************************************************/
  nearRectangle(vector3d<real_t> vecP1, vector3d<real_t> vecP2, vector3d<real_t> vecN)
  {
    p1 = vecP1;
    p2 = vecP2;
    n  = vecN;
    vUpdatePoints();
  } // End of nearRectangle()

  /** **************************************************************************
   * @brief Construct the nearRectangle object based on user input options_list.
   * @param ol Referenced options_list pointer
   ****************************************************************************/
  nearRectangle(const options_list &ol)
  {
    if(!ol.optionExists("p1") || !ol.optionExists("p2") || !ol.optionExists("n"))
    {
      printf("ERROR::nearRectangle: Need two points and a normal vector "
             "('p1', 'p2', and 'n').\n");
      Loci::Abort();
    }

    ol.getOptionUnits("p1", "m", p1);
    ol.getOptionUnits("p2", "m", p2);
    ol.getOptionUnits("n",   "", n);
    n *= 1./max<real_t>(norm(n),  1e-33);
    vUpdatePoints();

  } // End of nearRectangle()

  /** **************************************************************************
   * @brief 
   * @param p 
   * @return true 
   * @return false 
   ****************************************************************************/
  bool inGeomPt(vector3d<real_t> p) const
  {
    if(p.x >= xmin && p.x <= xmax &&
       p.y >= ymin && p.y <= ymax &&
       p.z >= zmin && p.z <= zmax)
    {
      return true;
    }

    // Check if the x-axis is close
    if(bEssentiallyEqual(p.x, xmax, 1e-3) &&
       p.y >= ymin && p.y <= ymax &&
       p.z >= zmin && p.z <= zmax)
    {
      return true;
    }

    return false;
  } // End of inGeomPt()

  /** **************************************************************************
   * @brief 
   * @return vector<vector3d<real_t>> 
   ****************************************************************************/
  vector<vector3d<real_t>> getPoints() const
  {
    vector<vector3d<real_t>> vec;
    vec.push_back(n);
    vec.push_back(pC);
    vec.push_back(p1);
    vec.push_back(p2);
    vec.push_back(p3);
    vec.push_back(p4);

    return vec;
  } // End of getPoints()

  real_t levelSet(vector3d<real_t> p) const
  {
    return 0.;
  }
private:
  /** **************************************************************************
   * @brief
   ****************************************************************************/
  void vUpdatePoints()
  {
    real_t dX = p1.x - p2.x;
    real_t dY = p1.y - p2.y;
    real_t dZ = p1.z - p2.z;

    bool bDX = bEssentiallyEqual(dX, 0, 1e-3);
    bool bDY = bEssentiallyEqual(dY, 0, 1e-3);
    bool bDZ = bEssentiallyEqual(dZ, 0, 1e-3);

    if((bDX && bDY) || (bDX && bDZ) || (bDY && bDZ))
    {
      printf("ERROR::nearRectangle: Width between the two points 'p1' and 'p2' "
             "can't be zero. The two points must create the opposite corners "
             "of a rectangle.\n");
      Loci::Abort();
    } // End If(Delta Check)

    // real_t   rDelta     = -1e-8;
    real_t   rDelta     = 0;
    vector3d<real_t> vecNegNorm = rDelta*n;
    vector3d<real_t> vecDeltaP2 = vecNegNorm + p2; // [m] Opposite corner of very small box

    if(bDX)
    {
      p3.x = p1.x;
      p3.y = p1.y;
      p3.z = p2.z;

      p4.x = p1.x;
      p4.y = p2.y;
      p4.z = p1.z;

      pC.x = p1.x;
      pC.y = p1.y+abs(dY)/2.;
      pC.z = p1.z+abs(dZ)/2.;
      totalArea = abs(dY*dZ);
    }else if(bDY)
    {
      pC.x = p1.x+abs(dX)/2.;
      pC.y = p1.y;
      pC.z = p1.z+abs(dZ)/2.;
      p3.x = p1.x;
      p3.y = p1.y;
      p3.z = p2.z;

      p4.x = p2.x;
      p4.y = p1.y;
      p4.z = p1.z;
      totalArea = abs(dX*dZ);
    }else
    {
      pC.x = p1.x+abs(dX)/2.;
      pC.y = p1.y+abs(dY)/2.;
      pC.z = p1.z;
      p3.x = p1.x;
      p3.y = p2.y;
      p3.z = p1.z;

      p4.x = p2.x;
      p4.y = p1.y;
      p4.z = p1.z;
      totalArea = abs(dX*dY);
    }

    xmax = max(p1.x, vecDeltaP2.x);
    xmin = min(p1.x, vecDeltaP2.x);

    ymax = max(p1.y, vecDeltaP2.y);
    ymin = min(p1.y, vecDeltaP2.y);

    zmax = max(p1.z, vecDeltaP2.z);
    zmin = min(p1.z, vecDeltaP2.z);

    return;
  } // End of vUpdatePoints()
}; // End Class nearRectangle
#endif // NEAR_RECTANGLE_H
