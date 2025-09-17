//#############################################################################
//#
//# Copyright 2014-2025, Mississippi State University
//#
//# The GridMover is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The GridMover software is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the gridMover software. If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include "gridMotion/rotor.h"
#include "gridMotion/gridTypes.h"

namespace gridMotion {

  using std::pair;
  using std::vector;
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::setw;
  using std::min ;
  using std::max ;
  /**
   * Input: 3x3 rotation matrix
   *
   * Output: A unit rotor which approximates the rotation component
   * of the given matrix.
   *
   * NOTE: Least squares approximation is utilized to help minimize
   * the error associated with this change in representation.  There
   * are nine equations and four unknowns.
   */
  Rotor::Rotor(const tens3d & M) {
    const vect3d X1(1,0,0);
    const vect3d X2 = dot(M,X1);
    const vect3d Y1(0,1,0);
    const vect3d Y2 = dot(M,Y1);
    const vect3d Z1(0,0,1);
    const vect3d Z2 = dot(M,Z1);
    vector< pair<vect3d, vect3d> > edgeData;
    edgeData.push_back(make_pair(X1,X2-X1));
    edgeData.push_back(make_pair(Y1,Y2-Y1));
    edgeData.push_back(make_pair(Z1,Z2-Z1));
    Rotor R(edgeData);
    alpha = R.alpha;
    beta = R.beta;
  }

  Rotor::Rotor(const vect3d & u, const vect3d & v)
    : alpha(1), beta(0,0,0) {

    const double uSq = dot(u,u);
    const double vSq = dot(v,v);
    if (uSq == 0.0 || vSq == 0.0) {
      // Degenerate vectors so rotation is zero
      return ;
    }

    const vect3d U = u/sqrt(uSq);
    const vect3d V = v/sqrt(vSq);
    if (1 - dot(U,V) == 0.0) {
      // Both vectors are oriented in the same direction.
      // This implies no rotation.  Return the default - NULL rotor.
      return;
    }
    else if (1 + dot(U,V) == 0.0) {
      // Both vectors are oriented in opposite directions.
      // This implies a Pi-radian rotation about an axis perpendicular
      // to the two given vectors.
      vect3d n(randf(), randf(), randf());
      // Orthogonalize this vector with respect to U.
      n = n - dot(n,U)*U;
      alpha = 0;
      beta = n/norm(n);
    }
    else {
      // Create a bisecting vector (oriented halfway between U and V).
      const vect3d w = U+V;
      const vect3d W = w/norm(w);
      // Create a rotor that will rotate U into V
      alpha = dot(U,W);
      beta = cross(U,W);
    }
  }

  /**
   * Input: Two edge vectors and two displacement vectors.  Here we
   * assume that the two edge vectors have a common origin and when the
   * displacements are added to the edge vectors, we get the two
   * displaced edges, also sharing a common origin.
   *
   * Algorithm: Generate a two rotors.  The first rotor rotates the
   * first edge vector into the same orientation as its displaced
   * counter-part. (If the length of the edge changes, it will have no
   * effect on the rotor.)  The first rotor is also applied to the
   * second edge to obtain an intermediate edge vector.  A second rotor
   * is then created which rotates this intermediate edge vector into
   * the same orientation as the second displaced edge vector.  (Again,
   * change in edge length will not matter.)  The second rotor is also
   * constructed in such a way that its axis of rotation is parallel to
   * the displaced first edge vector.  In this way, the second rotor
   * should leave the displaced first edge invariant.  Finally, the two
   * rotors are composed together to get a single rotation which will
   * rotate both edge vectors into the same orientation as their
   * displaced counter-parts.
   *
   * NOTE: If the angle between the edges should change, then there is
   * no guarantee that the obtained rotor will be able to reproduce both
   * of the desired orientation changes exactly.  This is to be expected
   * since the nodes are obviously not undergoing rigid body motion.  It
   * should be noted, however, that this algorithm will produce a rotor
   * which rotates vectors aligned with the first edge into the
   * orientation of its displaced counter-part.
   */
  Rotor::Rotor( const pair<vect3d,vect3d> & u,
                const pair<vect3d,vect3d> & v )
    : alpha(1.0), beta(0,0,0) {

    // First edge vector
    const vect3d u1 = u.first;
    // Second edge vector
    const vect3d v1 = v.first;
    // First edge displacement
    const vect3d du = u.second;
    // Second edge displacement
    const vect3d dv = v.second;

    const bool nullU1 = (dot(u1,u1) == 0.0);
    const bool nullV1 = (dot(v1,v1) == 0.0);

    if (nullU1 || nullV1) {
      // Degenerate edges, so no rotation
      return ;
    }
    // Let U1 be u1 normalized
    const vect3d U1 = u1/norm(u1);
    // Let V1 be v1 normalized
    const vect3d V1 = v1/norm(v1);

    const vect3d u2 = u1 + du;
    const vect3d v2 = v1 + dv;

    const bool nullDU = (dot(du,du) == 0.0);
    const bool nullDV = (dot(dv,dv) == 0.0);
    const bool nullRU = ( nullDU || ( (norm(du) - fabs(dot(du,U1)) == 0.0) ) );
    const bool nullRV = ( nullDV || ( (norm(dv) - fabs(dot(dv,V1)) == 0.0) ) );
    if (nullRU && nullRV) {
      // Orientation of both edges are unchanged. This implies that
      // the mesh is either unchanged, or undergoing translation
      // and/or scaling only.  Either way, there are no rotations
      // present. Return default NULL rotor.
      return;
    }
    else if (nullRU) {
      // If only the orientation of the first edge is unchanged, then
      // generate a rotor which has its axis of rotation aligned with
      // the first vector.
      Rotor RV = findOrthogonalRotor(u1, v1, v2);
      alpha = RV.alpha;
      beta = RV.beta;
    }
    else if (nullRV) {
      // If only the orientation of the second edge is unchanged, then
      // generate a rotor which has its axis of rotation aligned with
      // the second vector.
      Rotor RU = findOrthogonalRotor(v1, u1, u2);
      alpha = RU.alpha;
      beta = RU.beta;
    }
    else {
      // Find the rotor for u1 -> u2
      Rotor R1(u1, u2);
      // Then find an orthogonal rotor for R*v1 -> v2
      Rotor R2 = findOrthogonalRotor(u2, R1*v1, v2);

      // Compose the two rotations into one
      Rotor R = R2*R1;
      alpha = R.alpha ;
      beta = R.beta ;
    }
  }

  /**
   * Input: position and displacement of the current node under
   * consideration, and a list of positions and displacments of the
   * adjacent nodes relative to the current node and its displacement.
   * This list consists of all edges adjacent to a given node.  The
   * positions are the locations of the adjacent nodes relative to the
   * current node, and the displacements are the nodal displacements
   * minus the displacement of the current node.
   *
   * Algorithm: Use nonlinear least-squares fit to find the four
   * quaternion parameters.  (Actually, since the quaternion is
   * normalized, there are only three free parameters.)
   */
  Rotor::Rotor( const vector< pair<vect3d, vect3d> > & edgeDatai )
    : alpha(1), beta(0,0,0) {
    const int N = edgeDatai.size();
    if (N<1) {
      // no edges, so no rotation
      return ;
    }

    vector<pair<vect3d, vect3d> > edgeData(N) ;
    // Normalize deformations so that angles preserved, length of all vectors
    // before and after deformation are unity
    for(int i=0;i<N;++i) {
      vect3d v1 = edgeDatai[i].first ; // Initial location of edge
      vect3d v2 = v1+edgeDatai[i].second ; // deformed location of edge
      v1 *= 1./norm(v1) ; // normalize vectors
      v2 *= 1./norm(v2) ;
      edgeData[i].first = v1 ; // normalized edge vector
      edgeData[i].second = v2-v1 ; // displacement to normalized deformed edge
    }
    const int maxNewtonIter = 50;

    // Construct an initial guess for alpha and beta
    vector< pair<vect3d,vect3d> >::const_iterator itr;
    for (itr=edgeData.begin(); itr!=edgeData.end(); ++itr) {
      // Search until valid data is found
      vect3d u0  = itr->first;
      vect3d du0 = itr->second;
      if ( pow(dot(u0,du0),2) != dot(u0,u0)*dot(du0,du0) ) {
        // Compute the displaced edge
        vect3d v0  = u0 + du0;
        // Normalize the edges
        u0 = u0/(norm(u0)+1e-30);
        v0 = v0/(norm(v0)+1e-30);
        vect3d w0 = (u0+v0)/norm(u0+v0);
        // Initialize alpha and beta
        alpha = dot(u0, w0);
        beta = cross(u0, w0);
        break;
      }
    }
    if ( itr==edgeData.end() ) {
      // No valid displacements found - return identity rotor
      return;
    }

    for (int n=0; n<maxNewtonIter; ++n) {

      tens3d R(vect3d(pow(alpha,2)+pow(beta.x,2)-pow(beta.y,2)-pow(beta.z,2),
                      2*(beta.x*beta.y - alpha*beta.z),
                      2*(beta.x*beta.z + alpha*beta.y)),
               vect3d(2*(beta.y*beta.x + alpha*beta.z),
                      pow(alpha,2)-pow(beta.x,2)+pow(beta.y,2)-pow(beta.z,2),
                      2*(beta.y*beta.z - alpha*beta.x)),
               vect3d(2*(beta.z*beta.x - alpha*beta.y),
                      2*(beta.z*beta.y + alpha*beta.x),
                      pow(alpha,2)-pow(beta.x,2)-pow(beta.y,2)+pow(beta.z,2)));

      // Required for nomalization
      double ralpha = 1.0/((alpha>0)?(alpha+1e-30):(alpha-1e-30));
      // Required for chain rule to eliminate alpha
      double dAlphadBetaX = (-ralpha)*beta.x;
      double dAlphadBetaY = (-ralpha)*beta.y;
      double dAlphadBetaZ = (-ralpha)*beta.z;
      // Components of the Jacobian of the rotation matrix
      tens3d dRdAlpha( vect3d(  2*alpha,  -2*beta.z,  2*beta.y ),
                       vect3d(  2*beta.z,  2*alpha,  -2*beta.z ),
                       vect3d( -2*beta.y,  2*beta.x,  2*alpha  ) );
      tens3d dRdBetaX( vect3d(  2*beta.x,  2*beta.y,  2*beta.z ),
                       vect3d(  2*beta.y, -2*beta.x, -2*alpha  ),
                       vect3d(  2*beta.z,  2*alpha,  -2*beta.x ) );
      tens3d dRdBetaY( vect3d( -2*beta.y,  2*beta.x,  2*alpha  ),
                       vect3d(  2*beta.x,  2*beta.y,  2*beta.z ),
                       vect3d( -2*alpha,   2*beta.z, -2*beta.y ) );
      tens3d dRdBetaZ( vect3d( -2*beta.z, -2*alpha,   2*beta.x ),
                       vect3d(  2*alpha,  -2*beta.z,  2*beta.y ),
                       vect3d(  2*beta.x,  2*beta.y,  2*beta.z ) );

      // Jacobian matrix of Rotation operator
      // Chain in alpha dependence on beta to eliminate alpha from Jacobian
      vector3d< vector3d<double> > dRdQ[3];
      dRdQ[0] = dRdBetaX + dRdAlpha*dAlphadBetaX;
      dRdQ[1] = dRdBetaY + dRdAlpha*dAlphadBetaY;
      dRdQ[2] = dRdBetaZ + dRdAlpha*dAlphadBetaZ;

      // Construct transpose(A)*A and transpose(A)*b

      // LHS (Jacobian) of the system of equations to be solved
      double M[9] = {
        0, 0, 0,
        0, 0, 0,
        0, 0, 0
      };
      Mat<double> ATA(M, 3);
      // RHS (residual) of the system of equations to be solved
      double ATb[3] = {
        0, 0, 0
      };
      // Sum of the squares of the magnitudes of the residual vectors
      double bTb = 0.0;
      vector< pair<vect3d,vect3d> >::const_iterator itr;
      for (itr=edgeData.begin(); itr!=edgeData.end(); ++itr) {
        const vect3d u  = itr->first;
        const vect3d du = itr->second;
        const vect3d U  = vect3d(dot(R.x,u),dot(R.y,u),dot(R.z,u));
        // RHS: (residual for node k)
        //   -L(x[k]) = -( R(x[k]) - x[k] - dx[k] )
        const vect3d b = (u + du) - U;
        // The magnitude of the residual
        bTb += dot(b,b);
        // System matrix
        for (int i=0; i<3; ++i) {
          for (int j=0; j<3; ++j) {
            ATA[i][j] += dot(vect3d(dot(dRdQ[i].x,u),
                                    dot(dRdQ[i].y,u),
                                    dot(dRdQ[i].z,u)),
                             vect3d(dot(dRdQ[j].x,u),
                                    dot(dRdQ[j].y,u),
                                    dot(dRdQ[j].z,u)));
          }
          ATb[i] += dot(vect3d(dot(dRdQ[i].x,u),
                               dot(dRdQ[i].y,u),
                               dot(dRdQ[i].z,u)), b);
        }
      }
#define VERBOSE
#ifdef VERBOSE
      if (n == maxNewtonIter-1) {
      	cerr << "Max Newton iterations exceeded in finding rotation" << endl ;
	int lsz = edgeData.size() ;
	cerr << lsz << endl ;
	cerr.precision(16) ;
	for(int i=0;i<lsz;++i)
	  cerr << edgeData[i].first << ' ' << edgeData[i].second << endl ;
      }
#endif

      double x[3] = {
        0, 0, 0
      };

      pivot_type piv[6];
      ATA.decompose_lu_pivot(piv);
      double mindiag = min(min(fabs(ATA[0][0]),fabs(ATA[1][1])),
                           fabs(ATA[2][2])) ;
      if(mindiag < 1e-10) { // matrix singular, abort with no rotation
        beta = vect3d(0,0,0) ;
        alpha = 1 ;
        break ;
      }
        
      ATA.solve_lu_pivot(ATb, x, piv);

      double del = x[0]*x[0]+x[1]*x[1]+x[2]*x[2] ;
      double err = 0 ;
      double scale = 1.0 ;

      // Don't allow too much change in one step
      if(del > 0.01) 
	scale = 0.1/sqrt(del) ;

      // Now do a line search to make sure that we head towards the minimum
      int nsearch = 16 ;
      for(int i=0;i<nsearch;++i) {
	vect3d betat(beta.x+scale*x[0],beta.y+scale*x[1],beta.z+scale*x[2]) ;
	double beta2t = dot(betat,betat);
	double alphat = 0.0 ;
	if (beta2t > 1.) {
	  double rbeta = 1.0/(sqrt(beta2t));
	  betat *= rbeta;
	} else {
	  alphat = sqrt(1.-beta2t);
	}

	tens3d Rt(vect3d(pow(alphat,2)+pow(betat.x,2)-pow(betat.y,2)-pow(betat.z,2),
			2*(betat.x*betat.y - alphat*betat.z),
			2*(betat.x*betat.z + alphat*betat.y)),
		 vect3d(2*(betat.y*betat.x + alphat*betat.z),
			pow(alphat,2)-pow(betat.x,2)+pow(betat.y,2)-pow(betat.z,2),
			2*(betat.y*betat.z - alphat*betat.x)),
		 vect3d(2*(betat.z*betat.x - alphat*betat.y),
			2*(betat.z*betat.y + alphat*betat.x),
			pow(alphat,2)-pow(betat.x,2)-pow(betat.y,2)+pow(betat.z,2)));
	
	double bTbt= 0 ;
	for (itr=edgeData.begin(); itr!=edgeData.end(); ++itr) {
	  const vect3d u  = itr->first;
	  const vect3d du = itr->second;
	  const vect3d U  = vect3d(dot(Rt.x,u),dot(Rt.y,u),dot(Rt.z,u));
	  // RHS: (residual for node k)
	  //   -L(x[k]) = -( R(x[k]) - x[k] - dx[k] )
	  const vect3d b = (u + du) - U;
	  // The magnitude of the residual
	  bTbt += dot(b,b);
	}
	err = fabs(bTbt-bTb) ;
	if(bTbt < bTb)
	  break ;
	scale *= .5 ;
	if(n+1 == nsearch)
	  scale = 0 ;
      }

      // Update beta
      beta.x  += scale*x[0];
      beta.y  += scale*x[1];
      beta.z  += scale*x[2];


      // compute normalized quaternion
      double beta2 = dot(beta,beta);
      if (beta2 > 1) {
        double rbeta = 1.0/(sqrt(beta2)+1e-30);
        beta *= rbeta;
        alpha = 0.0;
      }
      else {
        alpha = sqrt(1-beta2);
      }
#ifdef TEST
      cout << "beta="<< beta << "err=" << err << ",bTb=" << bTb 
      	   << ",scale=" << scale << ",del=" << del << endl ;
#endif
      // If we have converged to a minimum value, exit iteration
      if(del < 1e-12 || err*scale < 1e-8 || bTb < 1e-7 || err < bTb*1e-6)
	break ;
    }
    beta.x = (fabs(beta.x)<1e-9)?0.0:beta.x ;
    beta.y = (fabs(beta.y)<1e-9)?0.0:beta.y ;
    beta.z = (fabs(beta.z)<1e-9)?0.0:beta.z ;
    double beta2 = dot(beta,beta);
    alpha = sqrt(max(1.-beta2,0.));
  }

  /**
   * Input: position and displacement of the current node under
   * consideration, and a list of positions and displacments of the
   * adjacent nodes relative to the current node and its displacement.
   * This list consists of all edges adjacent to a given node.  The
   * positions are the locations of the adjacent nodes relative to the
   * current node, and the displacements are the nodal displacements
   * minus the displacement of the current node.
   *
   * Algorithm: Use nonlinear least-squares fit to find the four
   * quaternion parameters.  (Actually, since the quaternion is
   * normalized, there are only three free parameters.)
   */
  Rotor::Rotor( const vector< pair<vect3d, vect3d> > & edgeData,
                const vect3d & axis )
    : alpha(1), beta(0,0,0) {
    const int N = edgeData.size();
    if (N<1) {
      //      cerr << "No edges provided to Rotor algorithm.\n";
      return ;
    }

    const int maxNewtonIter = 25;
    const double absTol = 1e-14;
    const double relTol = 1e-7;

    // Normalize the desired axis of rotation
    const vect3d ax = axis/norm(axis);
    // Find the unit Rotor which will map ax onto z-axis
    const vect3d zx(0,0,1);
    const vect3d wx = ax+zx/norm(ax+zx);
    const Rotor R0(ax,wx);
    double alp(1), bet(0);
    // Construct an initial guess for alp and bet
    vector< pair<vect3d,vect3d> >::const_iterator itr;
    for (itr=edgeData.begin(); itr!=edgeData.end(); ++itr) {
      // Search until valid data is found
      vect3d u  = itr->first;
      vect3d du = itr->second;
      // Rotate into the z-axis frame
      u  = R0*u;
      du = R0*du;
      // Project onto the xy-plane
      u.z  = 0.0;
      du.z = 0.0;
      // Determine if u and du are linearly independent
      if ( pow(dot(u,du),2) != dot(u,u)*dot(du,du) ) {
        vect3d v = u + du;
        // Normalize the projected edges
        u = u/norm(u);
        v = v/norm(v);
        // Find a vector oriented halfway between u and v
        vect3d w = (u+v)/norm(u+v);
        // Initialize alpha and beta
        alpha = dot(u, w);
        beta = cross(u, w);
        alp = alpha;
        bet = norm(beta);
        break;
      }
    }
    if ( itr==edgeData.end() ) {
      // No valid displacements found - return identity rotor
      return;
    }

    bool converged = false;
    double prevErr = 1e30;
    for (int n=0; n<maxNewtonIter && !converged; ++n) {
      double R[2][2] = {
        { alp, -bet },
        { bet,  alp }
      };

      // Required for nomalization
      double ralp = 1.0/alp;
      // Required for chain rule to eliminate alpha
      double dAlpdBet = (-ralp)*bet;
      // Components of the Jacobian of the rotation matrix
      double dRdAlp[2][2] = {
        {  1,  0 },
        {  0,  1 }
      };
      double dRdBet[2][2] = {
        {  0, -1 },
        {  1,  0 }
      };

      // Jacobian matrix of Rotation operator
      // Chain in alpha dependence on beta to eliminate alpha from Jacobian
      double dRdQ[2][2];
      for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
          dRdQ[i][j] = dRdBet[i][j] + dRdAlp[i][j]*dAlpdBet;
        }
      }

      // Construct transpose(A)*A and transpose(A)*b

      // LHS (Jacobian) of the system of equations to be solved
      double ATA = 0;

      // RHS (residual) of the system of equations to be solved
      double ATb = 0;

      // Sum of the squares of the magnitudes of the residual vectors
      double bTb = 0.0;
      vector< pair<vect3d,vect3d> >::const_iterator itr;
      for (itr=edgeData.begin(); itr!=edgeData.end(); ++itr) {
        vect3d u  = itr->first;
        vect3d du = itr->second;
        // Rotate into the z-axis frame
        u  = R0*u;
        du = R0*du;
        // Project onto the xy-plane
        u.z  = 0.0;
        du.z = 0.0;
        // Apply current rotor to u:  R(u)
        vect3d Ru = vect3d(R[0][0]*u.x+R[0][1]*u.y,
                           R[1][0]*u.x+R[1][1]*u.y,
                           0.0);
        // RHS: (residual for node k)
        //   -L(x[k]) = -( R(x[k]) - x[k] - dx[k] )
        vect3d b = (u + du) - Ru;
        // The magnitude of the residual
        bTb += dot(b,b);
        // System matrix
        ATA += dot(vect3d(dRdQ[0][0]*u.x + dRdQ[0][1]*u.y,
                          dRdQ[1][0]*u.x + dRdQ[1][1]*u.y,
                          0.0),
                   vect3d(dRdQ[0][0]*u.x + dRdQ[0][1]*u.y,
                          dRdQ[1][0]*u.x + dRdQ[1][1]*u.y,
                          0.0));
        ATb += dot(vect3d(dRdQ[0][0]*u.x + dRdQ[0][1]*u.y,
                          dRdQ[1][0]*u.x + dRdQ[1][1]*u.y,
                          0.0), b);
      }

      double absErr = sqrt(bTb);
      double relErr = fabs(prevErr-absErr)/prevErr;

//         cout << "Iteration: " << n << endl;
//         cerr << "Current Rotation Matrix:" << endl;
//         cerr << setw(15) << R[0][0]
//              << setw(15) << R[0][1] << endl;
//         cerr << setw(15) << R[1][0]
//              << setw(15) << R[1][1] << endl;
//         cerr << "Current Rotor:" << endl;
//         cerr << "alpha: "
//              << setw(15) << alp << endl;
//         cerr << "beta:  "
//              << setw(15) << 0.0
//              << setw(15) << 0.0
//              << setw(15) << bet << endl;
//         cerr << "AbsErr: " << absErr << endl;
//         cerr << "RelErr: " << relErr << endl;

      converged = (absErr < absTol) || (relErr < relTol);
      if (converged) {
        break;
      }
      else if (n == maxNewtonIter-1) {
        if(absErr > 1e-4 && relErr > 1e-2) {
          cerr << "Max Newton iterations exceeded: "
               << "iter(" << setw(2) << n << "), "
               << "absErr(" << setw(10) << absErr << "), "
               << "relErr(" << setw(10) << relErr << ")" << endl;
        }
        break;
      }

      prevErr = absErr;

      double x = ATb/ATA;
      bet += x;

      double bet2 = bet*bet;
      if (bet2 > 1) {
        double rbet = 1.0/sqrt(bet2);
        bet *= rbet;
        alp = 0.0;
      }
      else {
        alp = sqrt(1-bet2);
      }
    }
    alpha = alp;
    beta = (~R0)*vect3d(0,0,bet);
  }

  Rotor Rotor::findOrthogonalRotor( const vect3d & axis,
                                    const vect3d & v1,
                                    const vect3d & v2 ) {
    const vect3d Z = axis/norm(axis);
    const vect3d w1 = v1 - dot(v1,Z)*Z;
    const vect3d w2 = v2 - dot(v2,Z)*Z;
    //    const double w1Sq = dot(w1,w1);
    //    const double w2Sq = dot(w2,w2);
//     if (w1Sq == 0.0 || w2Sq == 0.0) {
//       cerr << "Failed to project vectors onto plane of rotation\n";
//       cerr << "axis: " << axis << endl;
//       cerr << "v1:   " << v1 << endl;
//       cerr << "v2:   " << v2 << endl;
//       cerr << "w1:   " << w1 << endl;
//       cerr << "w2:   " << w2 << endl;
//       cerr << "w1Sq: " << w1Sq << endl;
//       cerr << "w2Sq: " << w2Sq << endl;
//     }
    return Rotor(w1,w2);
  }

}

#ifdef TEST
using namespace std ;
int main() {
  int npnts = 0 ;
  cin >> npnts ;

  vector< pair<gridMotion::vect3d,gridMotion::vect3d> > vp(npnts) ;
  for(int i=0;i<npnts;++i) {
    cin >> vp[i].first >> vp[i].second ;
  }

  gridMotion::Rotor R(vp) ;
}
#endif
