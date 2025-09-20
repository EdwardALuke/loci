//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <Loci.h>
#include <flowTypes.h>

namespace flowPsi {
  inline void inviscidFlux(real iflux[],
			   real pg, real T, vect3d U,
			   vect3d an, real area,
			   real pambient, real Rt, real gamma) {
    real gm1 = gamma-1 ;
    real rho = (pg+pambient)/(Rt*T) ;
    
    real h0 = Rt*T*gamma/gm1 + 0.5*dot(U,U) ;
    const real e0 = h0-Rt*T ;
    const real uta = dot(U,an) ;
    const real ut = uta ; 
    real mdot = area*rho*ut ;
    real p = pg+pambient ;

    iflux[0] = mdot ;
    iflux[1] = mdot*U.x+area*pg*an.x ;
    iflux[2] = mdot*U.y+area*pg*an.y ;
    iflux[3] = mdot*U.z+area*pg*an.z ;
    iflux[4] = mdot*e0+area*p*uta ;
  }
    


  void hllc_flux(real iflux[],
		 real pgl, real Tl, vect3d Ul,
		 real pgr, real Tr, vect3d Ur,
		 vect3d an, real area,
		 real pambient, real Rt,real gamma) ;

  inline void inviscidRiemannFlux(real iflux[],
				  real pgl, real Tl, vect3d Ul,
				  real pgr, real Tr, vect3d Ur,
				  vect3d an, real area,
				  real pambient, real Rt,real gamma) {
    hllc_flux(iflux,pgl,Tl,Ul,pgr,Tr,Ur,an,area,pambient,Rt,gamma) ;
  }

}
