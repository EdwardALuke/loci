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
#include <Loci>

namespace Loci {
  inline void computeLSWeights(vector3d<real_t> *wvec,
                               const vector3d<real_t> *deltas,
                               const real_t *areas,
                               int ndeltas) {
    tmp_array<vector3d<real_t> > weights(ndeltas) ;

    real_t r[3][3] ;
    r[0][0] = 0 ;
    for(int i=0;i<ndeltas;++i) {
      r[0][0] += deltas[i].x*deltas[i].x ;
      weights[i] = deltas[i] ;
    }
    r[0][0] = sqrt(r[0][0]) ;
    r[0][0] = max(r[0][0],real_t(1e-30)) ;

    const real_t rr00 = 1./r[0][0] ;
    for(int i=0;i<ndeltas;++i)
      weights[i].x *= rr00 ;
    r[0][1] = 0 ;
    for(int i=0;i<ndeltas;++i)
      r[0][1] += weights[i].x*weights[i].y ;
    for(int i=0;i<ndeltas;++i)
      weights[i].y -= r[0][1]*weights[i].x ;
    r[0][2] = 0 ;
    for(int i=0;i<ndeltas;++i)
      r[0][2] += weights[i].x*weights[i].z ;
    for(int i=0;i<ndeltas;++i)
      weights[i].z -= r[0][2]*weights[i].x ;
    r[1][1] = 0 ;
    for(int i=0;i<ndeltas;++i)
      r[1][1] += weights[i].y*weights[i].y ;
    r[1][1] = sqrt(r[1][1]) ;
    r[1][1] = max(r[1][1],real_t(1e-30)) ;
    const real_t rr11 = 1./r[1][1] ;
    for(int i=0;i<ndeltas;++i)
      weights[i].y *= rr11 ;

    r[1][2] = 0 ;
    for(int i=0;i<ndeltas;++i)
      r[1][2] += weights[i].y*weights[i].z ;
    for(int i=0;i<ndeltas;++i)
      weights[i].z -= r[1][2]*weights[i].y ;
    r[2][2] = 0 ;
    for(int i=0;i<ndeltas;++i)
      r[2][2] += weights[i].z*weights[i].z ;
    r[2][2] = sqrt(r[2][2]) ;
    r[2][2] = max(r[2][2],real_t(1e-30)) ;
    const real_t rr22 = 1./r[2][2] ;
    for(int i=0;i<ndeltas;++i)
      weights[i].z *= rr22 ;
    // Compute R inverse
    //
    // Rinv[0] = R^-1[0][0]
    // Rinv[1] = R^-1[0][1]
    // Rinv[2] = R^-1[0][2]
    // Rinv[3] = R^-1[1][1]
    // Rinv[4] = R^-1[1][2]
    // Rinv[5] = R^-1[2][2]
    real_t Rinv[6] ;
    Rinv[0] = rr00 ;
    Rinv[1] = -r[0][1]*rr00*rr11 ;
    Rinv[2] = (r[0][1]*r[1][2]-r[0][2]*r[1][1])*rr00*rr11*rr22 ;
    Rinv[3] = rr11 ;
    Rinv[4] = -r[1][2]*rr11*rr22 ;
    Rinv[5] = rr22 ;

    for(int i=0;i<ndeltas;++i) {
      real_t weightx = Rinv[0]*weights[i].x + Rinv[1]*weights[i].y +
        Rinv[2]*weights[i].z ;
      real_t weighty = Rinv[3]*weights[i].y + Rinv[4]*weights[i].z ;
      real_t weightz = Rinv[5]*weights[i].z ;
      wvec[i] = vector3d<real_t>(weightx,weighty,weightz)*areas[i] ;
    }
  }
}
