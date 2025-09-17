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
#ifndef MAPVEC_H
#define MAPVEC_H

#include <DMapVec.h>
#include <MapVec_def.h>
#include <MapVec_impl.h>

namespace Loci {
  template <int M>
  inline void remapMap(MapVec<M> &m, const Map &remap, Map &tmp) {
    for(int i=0;i<M;++i) {
      FORALL(m.domain(),fc) {
        tmp[remap[fc]] = m[fc][i] ;
      } ENDFORALL ;
      FORALL(m.domain(),fc) {
        m[fc][i] = tmp[fc] ;
      } ENDFORALL ;
    }
  }
}
    

#endif
