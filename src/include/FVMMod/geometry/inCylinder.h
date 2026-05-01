/** ****************************************************************************
 * @file      inCylinder.h
 * @author    Mark A. Hunt (CFDRC)
 * @brief     Contains the inCylinder class. Mostly used for defining regions in
 *            the vars file.
 * @date      2022-11-07
 * @copyright CFDRC Copyright (c) 2022
 * @version   0.1
 * @details
 ******************************************************************************/
#ifndef INCYLINDER_H
#define INCYLINDER_H
#include "geomTest.h"

/// @brief Class that creates a bounding cylinder for two given points and a
///        radius.
class inCylinder : public geomTest
{
  real_t  mag; //!< [m]
  real_t   r;   //!< [m]
  real_t   r2;  //!< [m^2]
  vector3d<real_t> p1;  //!< [m]
  vector3d<real_t> p2;  //!< [m]
  vector3d<real_t> n;   //!< [-]
public:

/**
 * @brief class constructor
 */
 inCylinder() {
  p1 = vector3d<real_t>(0.,0.,0.);
  p2 = vector3d<real_t>(0.,0.,1.);
  r = 0.;
  r2 = 0.;
  mag = 1.;
 }

  /** **************************************************************************
   * @brief
   * @param ri
   * @param p1i
   * @param p2i
   ****************************************************************************/
  inCylinder(real_t ri,  vector3d<real_t> p1i,vector3d<real_t> p2i)
  {
    r    = ri;
    r2   = r*r;
    p1   = p1i;
    p2   = p2i;
    n    = p2-p1;
    mag  = norm(n);
    n   *= 1./mag;
  } // End of inCylinder()

  /** **************************************************************************
   * @brief Construct the inCylinder object based on user input options_list.
   * @param ol Referenced options_list pointer
   ****************************************************************************/
  inCylinder(const options_list &ol)
  {
    // Here get sphere information
    if(!ol.optionExists("radius") ||
       !ol.optionExists("p1")     ||
       !ol.optionExists("p2"))
    {
      printf("ERROR::inCylinder: Need a radius and two axis points (p1 and p2)\n");
      Loci::Abort();
    }

    ol.getOptionUnits("radius", "m", r);
    ol.getOptionUnits("p1",     "m", p1);
    ol.getOptionUnits("p2",     "m", p2);
    n    = p2-p1;
    mag  = norm(n);
    n   *= 1./mag;
    r2   = r*r;
  } // End of inCylinder()

  /** **************************************************************************
   * @brief  Determine if a position vector 'p' is within the bounding cylinder.
   * @param  p  [m] Position vector of a point.
   * @return true  - If 'p' is inside bounding cylinder.
   * @return false - If 'p' is outside bounding cylinder.
   ****************************************************************************/
  bool inGeomPt(vector3d<real_t> p) const
  {
    real_t v = dot(p-p1, n);
    if(v < 0.0 || v > mag)
      return false;
    // project point onto axis
    vector3d<real_t> paxis = p1 + n*v;
    if(dot(p-paxis, p-paxis) > r2) // Check radius
      return false; // outside
    return true;
  } // End of inGeomPt()

   real_t distToSurface(vector3d<real_t> pt) const {
    vector3d<real_t> p = pt-pt1 ;
    real_t t = dot(p,pt2-pt1)*rlen2 ;
    if(t < 0 || t > 1) {
      vector3d<real_t> v(0,0,0) ;
      vector3d<real_t> n = pt2-pt1 ;
      n*= 1./norm(n) ;
      if(t > 0.5) {
	      p = pt-pt2 ;
      }
      vector3d<real_t> ps = p-dot(p,n)*n ;
      real_t rt2 = dot(ps,ps) ;
      if(rt2 > radius2)
	      ps *= sqrt(radius2/rt2) ;
      return norm(ps-p) ;
    }

    vector3d<real_t> v = t*(pt2-pt1) ;
    real_t r2axis = dot(v-p,v-p) ;
    if(r2axis < radius2)
      return 0 ;
    return sqrt(r2axis)-radius ;
  }
}
#endif // INCYLINDER_H
