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
#include "flowTypes.h"

namespace flowPsi {
  inline void cv_inflow(real &pgb, real &Tb, vect3d &ub,
			real pgi, real Ti, vect3d ui,
			real pg0, real T0, vect3d u0,
			vect3d n, real us_n, 
			real Pambient, real gamma, real Rtilde,
			real Beta) {
    const real uit = dot(ui,n)-us_n ; //dot(uin-us,n) ;
    const real u0t = dot(u0,n)-us_n ; //dot(uref-us,n) ;

    const real a02 = gamma*Rtilde*T0 ;
    const real a0 = sqrt(a02) ;
    const real r0 = (pg0+Pambient)/(Rtilde*T0) ;
    const real MinPressure = 0.1*(pg0+Pambient)-Pambient ;
    const real BetaM = 0.5*(1.0-Beta) ;
    const real uitBeta  = uit*(1.-Beta) ;
    const real sigma2 = uitBeta*uitBeta + 4.*Beta*a0*a0 ;
    const real sigma = sqrt(sigma2) ;

    pgb = 0.5*(pg0 + pgi - 0.5*r0*sigma*(u0t-uit)) 
      -((pgi-pg0) + r0*(uit-u0t)*uit*BetaM)*uit*BetaM/(2.*sigma) ;
    pgb = max(pgb,MinPressure) ;

    ub = u0 + 2.*n*(pgb-pg0)/(r0*(sigma+uit*BetaM)) ;
    real rhob = r0 + (pgb-pg0)*(sigma-uit*BetaM)/
      (Beta*a0*a0*(sigma+uit*BetaM)) ;

    Tb = (pgb+Pambient)/(Rtilde*rhob) ;
  }
  inline void cv_outflow(real &pgb, real &Tb, vect3d &ub,
			 real pgi, real Ti, vect3d ui,
			 real pg0, real T0, vect3d u0,
			 vect3d n, real us_n, 
			 real Pambient, real gamma, real Rtilde,
			 real Beta) {
    const real uit = dot(ui,n)-us_n ; //dot(uin-us,n) ;
    const real u0t = dot(u0,n)-us_n ; //dot(uref-us,n) ;

    const real a02 = gamma*Rtilde*T0 ;
    const real a0 = sqrt(a02) ;
    const real r0 = (pg0+Pambient)/(Rtilde*T0) ;
    const real MinPressure = 0.1*(pgi+Pambient)-Pambient ;
    const real BetaM = 0.5*(1.0-Beta) ;

    const real uitBeta  = uit*(1.-Beta) ;
    const real sigma2 = uitBeta*uitBeta + 4.*Beta*a0*a0 ;
    const real sigma = sqrt(sigma2) ;

    pgb = .5*(pgi+pg0-0.5*r0*sigma*(u0t-uit)) 
      +((pgi-pg0) + r0*(uit-u0t)*uit*BetaM)*uit*BetaM/(2.*sigma) ;
    pgb = max(pgb,MinPressure) ;
    ub = ui + 2.*n*(pgi-pgb)/(r0*(sigma+uit*BetaM)) ;

    const real rhob = r0 + (pgb-pgi)*(sigma-uit*BetaM)/
      (Beta*a0*a0*(sigma+uit*BetaM)) ;

    Tb = (pgb+Pambient)/(Rtilde*rhob) ;
  }      
}
