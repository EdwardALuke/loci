/** ****************************************************************************
 * @file      geomTest.h
 * @author    Mark A. Hunt (CFDRC)
 * @brief
 * @date      2022-11-07
 * @copyright CFDRC Copyright (c) 2022
 * @version   0.1
 * @details
 ******************************************************************************/
#ifndef GEOMTEST_H
#define GEOMTEST_H

#include <Loci.h>
#include "sgn.h"

/// @brief Class that
class geomTest : public Loci::CPTR_type
{
public:
  virtual bool inGeomPt(vector3d<real_t> pt) const = 0;
  virtual real_t levelSet(vector3d<real_t> pt) const = 0; // +inside, -outside
  virtual real_t distToSurface(vector3d<real_t> pt) const = 0;
};

// level-set operations for geometry
inline real_t LS_union(real_t pbase, real_t pnew)
{
  real_t sol = (pnew<0?-1.:1.)*min<real_t>(abs(pbase), abs(pnew));
  if(pbase>0. || pbase*pnew<=0.)
    sol = max<real_t>(pbase, pnew);
  return sol;
}

inline real_t LS_minus(real_t pbase, real_t pnew)
{
  real_t pbase0 = -pbase;
  real_t sol = (pnew<0?-1.:1.)*min<real_t>(abs(pbase0), abs(pnew));
  if(pbase0>0. || pbase0*pnew<=0.)
    sol = max<real_t>(pbase0, pnew);
  return -sol;
}

#endif // GEOMTEST_H
