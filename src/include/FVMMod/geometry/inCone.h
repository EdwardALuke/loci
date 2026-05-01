/** ****************************************************************************
 * @file      inCone.h
 * @author    Mark A. Hunt (CFDRC)
 * @brief     Contains the inCone class. Mostly used for defining regions in the
 *            vars file.
 * @date      2022-11-07
 * @copyright CFDRC Copyright (c) 2022
 * @version   0.1
 * @details
 ******************************************************************************/
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

}; // End Class inCone

#endif // INCONE_H
