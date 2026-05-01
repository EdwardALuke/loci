/** ****************************************************************************
 * @file      inSphere.h
 * @author    Mark A. Hunt (CFDRC)
 * @brief     Contains the inSphere class. Mostly used for defining regions in
 *            the vars file.
 * @date      2022-11-07
 * @copyright CFDRC Copyright (c) 2022
 * @version   0.1
 * @details
 ******************************************************************************/
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
