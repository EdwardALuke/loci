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
#ifndef INTERVALSET_H
#define INTERVALSET_H 
#define ENTITY64

#include <Tools/intervalSet_def.h>
#include <Tools/intervalSet_impl.h>
namespace Loci {
  typedef  genIntervalSet<int_type> intervalSet;
  typedef genSequence<int_type> sequence;
  typedef genIntervalSet<gEntity> gEntitySet;
  typedef genSequence<gEntity> gSequence;
  extern const genIntervalSet<int_type> EMPTY ;
  extern const genIntervalSet<gEntity> GEMPTY ;
  extern const int_type UNIVERSE_MIN ;
  extern const int_type UNIVERSE_MAX ;
  extern const gEntity GUNIVERSE_MIN ;
  extern const gEntity GUNIVERSE_MAX ;
  typedef std::pair<int_type, int_type> interval;
  typedef std::pair<gEntity, gEntity> gInterval;
  inline intervalSet gconvert(gEntitySet g) {
    int nivals = g.num_intervals() ;
    if(g==~GEMPTY)
      return ~EMPTY ;
    intervalSet tmp ;
    for(int i=0;i<nivals;++i)
      tmp += interval(g[i].first,g[i].second) ;
    return tmp ;
  }
}
#endif
