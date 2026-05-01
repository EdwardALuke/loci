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
