/** ****************************************************************************
 * @file      leftPlane.h
 * @author    Mark A. Hunt (CFDRC)
 * @brief     Contains the leftPlane class. Mostly used for defining regions in
 *            the vars file.
 * @date      2022-11-07
 * @copyright CFDRC Copyright (c) 2022
 * @version   0.1
 * @details
 ******************************************************************************/
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
