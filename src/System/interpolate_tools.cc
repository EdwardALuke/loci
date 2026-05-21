//#############################################################################
//#
//# Copyright 2015-2025, Mississippi State University
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
#include <interpolate.h>
#include <fact_db.h>

namespace Loci {
  using std::vector ;
  using std::pair ;
  typedef vector3d<double> vect3d ;
 
  // Compute stencil using nearest point in 8 octants
  vector<int> get_stencil(const kdTree::KDTree<float> &kd,
			  kdTree::coord3df ccenter,
			  double deltai) {
    vector<int> neighbors(8) ;
    double delta = 1.4142*deltai ;
    double rmin = delta*delta ; // Note rmin is radius squared.
    double rmin_ref = rmin ;
    
    int id = kd.find_closest(ccenter,rmin) ;

    if(id < 0) {
      // If no points in stencil radius, return closest point
      vector<int> n ;
      id = kd.find_closest(ccenter) ;
      if(id >=0 && id !=std::numeric_limits<int>::max())
        n.push_back(id) ;
      return n ;
    }
    if(rmin <= 1e-30) {
      vector<int> n ;
      if(id != std::numeric_limits<int>::max())
	n.push_back(id) ;
      return n ;
    }

    
    // First gather postive z quadrants
    kdTree::KDTree<float>::bounds box ;

    for(int i=0;i<3;++i) {
      box.minc[i] = ccenter[i] ;
      box.maxc[i] = ccenter[i]+delta ;
    }

    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[0] = id;

    box.maxc[0] = ccenter[0] ;
    box.minc[0] = ccenter[0]-delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[1] = id;

    box.maxc[1] = ccenter[1] ;
    box.minc[1] = ccenter[1]-delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[2] = id;

    box.minc[0] = ccenter[0] ;
    box.maxc[0] = ccenter[0]+delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[3] = id;

    // Now gather negative z quadrants
    for(int i=0;i<3;++i) {
      box.minc[i] = ccenter[i] ;
      box.maxc[i] = ccenter[i]+delta ;
    }
    box.maxc[2] = ccenter[2] ;
    box.minc[2] = ccenter[2]-delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[4] = id;

    box.maxc[0] = ccenter[0] ;
    box.minc[0] = ccenter[0]-delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[5] = id ;

    box.maxc[1] = ccenter[1] ;
    box.minc[1] = ccenter[1]-delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[6] = id ;

    box.minc[0] = ccenter[0] ;
    box.maxc[0] = ccenter[0]+delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[7] = id ;

    return neighbors ;
  }

  // Compute stencil using nearest point in 8 octants 
  vector<int> get_stencil(const kdTree::KDTree<double> &kd,
			    kdTree::coord3d ccenter,
			    double deltai) {
    vector<int> neighbors(8) ;
    double delta = 1.4142*deltai ;
    double rmin = delta*delta ; // Note rmin is radius squared.
    double rmin_ref = rmin ;
    
    int id = kd.find_closest(ccenter,rmin) ;

    if(id < 0 ) {
      // If no points in stencil radius, return closest point
      vector<int> n ;
      id = kd.find_closest(ccenter) ;
      if(id >=0 && id != std::numeric_limits<int>::max())
        n.push_back(id) ;
      return n ;
    }
    if(rmin <= 1e-30) {
      vector<int> n ;
      if(id != std::numeric_limits<int>::max())
	n.push_back(id) ;
      return n ;
    }

    
    // First gather postive z quadrants
    kdTree::KDTree<double>::bounds box ;

    for(int i=0;i<3;++i) {
      box.minc[i] = ccenter[i] ;
      box.maxc[i] = ccenter[i]+delta ;
    }

    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[0] = id;

    box.maxc[0] = ccenter[0] ;
    box.minc[0] = ccenter[0]-delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[1] = id;

    box.maxc[1] = ccenter[1] ;
    box.minc[1] = ccenter[1]-delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[2] = id;

    box.minc[0] = ccenter[0] ;
    box.maxc[0] = ccenter[0]+delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[3] = id;

    // Now gather negative z quadrants
    for(int i=0;i<3;++i) {
      box.minc[i] = ccenter[i] ;
      box.maxc[i] = ccenter[i]+delta ;
    }
    box.maxc[2] = ccenter[2] ;
    box.minc[2] = ccenter[2]-delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[4] = id;

    box.maxc[0] = ccenter[0] ;
    box.minc[0] = ccenter[0]-delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[5] = id ;

    box.maxc[1] = ccenter[1] ;
    box.minc[1] = ccenter[1]-delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[6] = id ;

    box.minc[0] = ccenter[0] ;
    box.maxc[0] = ccenter[0]+delta ;
    rmin = rmin_ref ;
    id = kd.find_closest_box(ccenter,box,rmin) ;
    if(id == std::numeric_limits<int>::max())
      id = -1 ;
    neighbors[7] = id ;

    return neighbors ;
  }

  // Compute stencil weights from a set of candidate neighbors
  // The code outputs an updated neighbors list and stencil weights
  // This routine uses the following algorithm to determine the stencil:
  //  1) Since 8 points come from 8 ocatants, they form a hexahedron.  From
  //     this hexahedron form the two 5 tetrahedra decompositions and find out
  //     which tetrahedra which has the closest four points and contains the
  //     interpolating point.  If a tetrahedra is found, use compute
  //     barycentric coordinates to obtain weights.
  //  2) If no tetrahedra is formed, then some octants must not have contained
  //     points and the case should degenerate to a 2-D case.  For this we
  //     check all possible triangles  that include closest point and select
  //     the triangle with the smallest perimiter that we can project our
  //     interpolant onto.  Use the barycentric corrdinates for this triangle
  //     to compute weights.
  //  3) If no triangle is found from step 2 select edge formed from two
  //     closest points.  Use projected point distances to compute weights.

  // Hexahedron from out of get_stencil will be organized as so:
  //
  //                       2----------3
  //                      /|         /|
  //                     / |        / |
  //                    /  |       /  |
  //                   6---+------7   |
  //                   |   |      |   |
  //                   |   1------+---0
  //                   |  /       |  /
  //                   | /        | /
  //                   |/         |/
  //                   5----------4
  //
  
  void stencil_weights(std::vector<double> &w,
                       std::vector<int> &neighbors,
                       const store<vect3d> &loc,
                       vect3d ipnt) {
    using std::swap ;
    // Find closest point, set it to n0
    int sz = neighbors.size() ;
    if(sz == 0 ) 
      return ;
    if(sz == 1) {
      vector<double> W(1) ;
      W[0] = 1.0 ;
      w.swap(W) ;
      return ;
    }
    FATAL(sz != 8) ;
    // Now search for corner tetrahedra
    const int corner_tets[8][3]
      = { {1, 3, 4}, {0, 5, 2}, {1, 6, 3}, {0, 2, 7},
          {0, 7, 5}, {1, 4, 6}, {2, 5, 7}, {3, 6, 4} } ;

    int corner = -1 ;
    double dist = 1e33 ;
    for(int i=0;i<8;++i) {
      const int p0 = neighbors[i] ;
      const int p1 = neighbors[corner_tets[i][0]] ;
      const int p2 = neighbors[corner_tets[i][1]] ;
      const int p3 = neighbors[corner_tets[i][2]] ;
      // Check if tet is valid
      if(p0 >=0 && p1 >= 0 && p2 >= 0 && p3 >= 0) {
        // Now test to see if point is in this corner tet.
        vect3d v1 = loc[p1]-loc[p2] ;
        vect3d v2 = loc[p3]-loc[p2] ;
        vect3d vp = ipnt-loc[p2] ;
        if(dot(vp,cross(v1,v2)) < 0.0) {
            double sum = dot(vp,vp) ;
            vp = ipnt - loc[i] ;
            sum += dot(vp,vp) ;
            vp = ipnt - loc[p1] ;
            sum += dot(vp,vp) ;
            vp = ipnt - loc[p3] ;
            sum += dot(vp,vp) ;
            if(sum < dist) {
              dist = sum ;
              corner = i ;
            }
        }
      } 
    }

    int t0=-1,t1=-1,t2=-1,t3=-1 ;
    if(corner >= 0) {
      const int i = corner ;
      t0 = neighbors[i] ;
      t1 = neighbors[corner_tets[i][0]] ;
      t2 = neighbors[corner_tets[i][1]] ;
      t3 = neighbors[corner_tets[i][2]] ;
    } else {
      // not a corner, must be in center.  Select best of two center ones.
      int p0 = neighbors[1] ;
      int p1 = neighbors[3] ;
      int p2 = neighbors[4] ;
      int p3 = neighbors[6] ;

      if(p0 >=0 && p1 >= 0 && p2 >= 0 && p3 >= 0) {
        t0 = p0 ;
        t1 = p1 ;
        t2 = p2 ;
        t3 = p3 ;
        dist =
          dot(ipnt-loc[t0],ipnt-loc[t0])+
          dot(ipnt-loc[t1],ipnt-loc[t1])+
          dot(ipnt-loc[t2],ipnt-loc[t2])+
          dot(ipnt-loc[t3],ipnt-loc[t3]) ;
      }
      p0 = neighbors[0] ;
      p1 = neighbors[2] ;
      p2 = neighbors[5] ;
      p3 = neighbors[7] ;
      if(p0 >=0 && p1 >= 0 && p2 >= 0 && p3 >= 0) {
        double dist2 =
          dot(ipnt-loc[p0],ipnt-loc[p0])+
          dot(ipnt-loc[p1],ipnt-loc[p1])+
          dot(ipnt-loc[p2],ipnt-loc[p2])+
          dot(ipnt-loc[p3],ipnt-loc[p3]) ;
        if(dist2 < dist) {
          dist = dist2 ;
          t0 = p0 ;
          t1 = p1 ;
          t2 = p2 ;
          t3 = p3 ;
        }
      }
    }
    if(t0 >=0) { // found tetrahedra
      // Compute tetrahedra barycentric coordinates for weights
      // and modify stencil to only include the four points that define
      // the bounding tetrahedra.
      vector<int> N(4) ;
      vector<double> W(4) ;
      N[0] = t0 ;
      N[1] = t1 ;
      N[2] = t2 ;
      N[3] = t3 ;
      vect3d v0 = loc[N[0]] ;
      vect3d v1 = loc[N[1]]-v0 ;
      vect3d v2 = loc[N[2]]-v0 ;
      vect3d v3 = loc[N[3]]-v0 ;
      vect3d vc = ipnt-v0 ;
      
      double w3 = dot(cross(v2,v1),vc) ;
      double w2 = dot(cross(vc,v1),v3) ;
      double w1 = dot(cross(v2,vc),v3) ;
      double vol = dot(cross(v2,v1),v3) ;
      double w0 = vol-w1-w2-w3 ;
      if(vol >= 1e-30 && w0>=0.0 && w1>= 0.0 && w2 >= 0.0 && w3 >= 0.0 ) {
        double rvol = 1./(vol) ;
        double aspect6 = dist*dist*dist*rvol*rvol ;
        if(aspect6 < 1.6e10) { // aspect ratio less than 50
          W[0] = w0*rvol ;
          W[1] = w1*rvol ;
          W[2] = w2*rvol ;
          W[3] = w3*rvol ;
          w.swap(W) ;
          neighbors.swap(N) ;
          return ;
        } 
      }
      //      cerr << "degenerate tet, corner = "<< corner << ",w0="
      //	   << w0 << ",w1="<<w1 << ",w2=" <<w2 << ",w3=" 
      //	   << ",vol=" << vol << endl ;
    }
    // Got here, so it means we are degenerate.
    int nt[8] ;
    sz = 0 ;
    double mind = 1e30 ;
    for(int i=0;i<8;++i)
      if(neighbors[i] >= 0) {
	nt[sz] = neighbors[i] ;
	sz++ ;
	vect3d dv = ipnt-loc[neighbors[i]] ;
	const double d = dot(dv,dv) ;
	if(d<mind) {
	  mind = d ;
	  swap(nt[sz-1],nt[0]) ;
	}
      }
    if(sz == 0) {
      vector<int> N ;
      vector<double> W ;
      neighbors.swap(N) ;
      w.swap(W) ;
      return ;
    }
    if(sz == 1) {
      vector<int> N(1) ;
      vector<double> W(1) ;
      N[0] = nt[0] ;
      W[0] = 1.0 ;
      neighbors.swap(N) ;
      w.swap(W) ;
      return ;
    }      
    vect3d v0 = loc[nt[0]] ;
    vect3d vc = ipnt-v0 ;
    // Now find best triangle fit
    mind = 1e33 ;
    int n1=-1, n2 = -1 ;
    for(int i=1;i<sz;++i) {
      const vect3d v1 = loc[nt[i]]-v0 ;
      const double d1 = dot(loc[nt[i]]-ipnt,loc[nt[i]]-ipnt) ;
      for(int j=i+1;j<sz;++j) {
        const vect3d dv = loc[nt[j]]-ipnt;
        const vect3d v2 = loc[nt[j]]-v0 ;
        const double d2 = dot(dv,dv) ;
        // find smallest perimiter
        if(mind > d1+d2) {
          const vect3d n = cross(v1,v2) ; // Compute face normal
          const double a2 = dot(n,n) ;
          const double d = dot(vc,n) ; // projected distance to triangle
          const double ra2 = 1./(a2+1e-60) ;
          // Compute projected point on triangle plane
          vect3d pnt = (vc-d*ra2*n) ;
          // Now check to see if pnt is in triangle
          // Compute barycentric coordinates and make sure all are positive
          double c1 = dot(n,cross(v1,pnt)) ;
          double c2 = dot(n,cross(pnt,v2)) ;
          double c3 = dot(n,cross(v1-pnt,v2-pnt)) ;
      
          if(c1 >= 0. && c2 >= 0. && c3 >= 0.) {
            mind = d1+d2 ;
            n2 = j ;
            n1 = i ;
          }
        }
      }
    }

    // found triangle, setup weights    
    if(n2 != -1) {
      swap(nt[1],nt[n1]) ;
      swap(nt[2],nt[n2]) ;
      n1 = 1 ;
      n2 = 2 ;
      vect3d v1 = loc[nt[1]]-v0 ;
      vect3d v2 = loc[nt[2]]-v0 ;
      vect3d n = cross(v1,v2) ;
      double ra2 = 1./(dot(n,n)+1e-30) ;
      double d = dot(vc,n) ; // projected distance to triangle
      vect3d pnt = (vc-d*ra2*n) ;
      double w2 = dot(n,cross(v1,pnt)) ;
      double w1 = dot(n,cross(pnt,v2)) ;
      double w0 = dot(n,cross(v1-pnt,v2-pnt)) ;
      if(w1>=0 && w2>=0 && w0 >=0 && w0+w1+w2 > 1e-30) {
	vector<int> N(3) ;
	vector<double> W(3) ;
	N[0] = nt[0] ;
	N[1] = nt[1] ;
	N[2] = nt[2] ;
        double wsr = 1./(w0+w1+w2) ;
	W[0] = w0*wsr ;
	W[1] = w1*wsr ;
	W[2] = w2*wsr ;
	w.swap(W) ;
	neighbors.swap(N) ;
	return ;
      }
      //      cerr << "degenerate triangle, w0="<<w0 << "w1="<<w1<< "w2="<<w2 << endl ;
    }
    
    // Now we are reduced to a line segment, project onto it to get weights
    // Find second closest point
    mind = dot(ipnt-loc[nt[1]],ipnt-loc[nt[1]]) ;
    for(int i=2;i<sz;++i) {
      vect3d dv = ipnt-loc[nt[i]] ;
      const double d = dot(dv,dv) ;
      if(d<mind) {
        mind = d ;
        swap(nt[i],nt[1]) ;
      }
    }
    // compute weights
    vector<int> N(2) ;
    vector<double> W(2) ;
    N[0] = nt[0] ;
    N[1] = nt[1] ;
    vect3d v1 = loc[nt[1]]-v0 ;
    W[1] = max(double(0.0),min(double(1.0),dot(vc,v1)/(dot(v1,v1)+1e-30))) ;
    W[0] = 1.-W[1] ;
    w.swap(W) ;
    neighbors.swap(N) ;
    return ;
  }
  
  int collectPointsSizes(const kdTree::KDTree<float> &kd,
			 kdTree::KDTree<float>::bounds bnd) {
    MPI_Comm comm = get_exec_comm() ;
    MEMORY_PROFILE(collectPointsBegin) ;
    // Communicate bounds request to other processors
    using namespace kdTree ;    
    int p = get_exec_size() ;
    vector<KDTree<float>::bounds> bnd_req(p) ;
    MPI_Allgather(&bnd,6,MPI_FLOAT,&bnd_req[0],6,MPI_FLOAT,comm) ;


    // Now communicate points that processors need to build stencil
    vector<int> scounts(p,0) ;
    vector<KDTree<float>::coord_info> pntlist ;

    for(int i=0;i<p;++i) { // Find points in i'th processor bounding box
      scounts[i] = kd.count_box(bnd_req[i]) ;
    }
    vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;
    int cnt = rcounts[0] ;
    for(int i=1;i<p;++i)
      cnt += rcounts[i] ;
    return cnt ;
  }


  void getStencilBoundingBox2(kdTree::KDTree<float>::bounds &bnd,
			      double &delta,
			      const kdTree::KDTree<float> &kd,
			      const vector3d<double> pnts[],int start, int end) {
    double deltain = delta ;
    for(int d=0;d<3;++d) {
      bnd.minc[d] = .25*std::numeric_limits<float>::max() ;
      bnd.maxc[d] = -.25*std::numeric_limits<float>::max() ;
    }

    for(int i=start;i<end;++i) {
      bnd.maxc[0] = max(bnd.maxc[0],realToFloat(pnts[i].x)) ;
      bnd.maxc[1] = max(bnd.maxc[1],realToFloat(pnts[i].y)) ;
      bnd.maxc[2] = max(bnd.maxc[2],realToFloat(pnts[i].z)) ;
      bnd.minc[0] = min(bnd.minc[0],realToFloat(pnts[i].x)) ;
      bnd.minc[1] = min(bnd.minc[1],realToFloat(pnts[i].y)) ;
      bnd.minc[2] = min(bnd.minc[2],realToFloat(pnts[i].z)) ;
    }
    int npnts = collectPointsSizes(kd,bnd) ;

    if(npnts == 0) {
      // just estimate based on what we have on this processor
      for(int i=start;i<end;++i) {
	double rmin = 1e30 ;
	kdTree::coord3df pt ;
	pt[0] = realToFloat(pnts[i].x) ;
	pt[1] = realToFloat(pnts[i].y) ;
	pt[2] = realToFloat(pnts[i].z) ;

	kd.find_closest(pt,rmin) ;
	delta = max(delta,sqrt(rmin)) ;
      }
    } else {
      // Compute a delta to expand the bounding box
      double d1 = bnd.maxc[0]-bnd.minc[0] ;
      double d2 = bnd.maxc[1]-bnd.minc[1] ;
      double d3 = bnd.maxc[2]-bnd.minc[2] ;
      if(d1<d2)
	std::swap(d1,d2) ;
      if(d1<d3)
	std::swap(d1,d3) ;
      if(d2<d3)
	std::swap(d2,d3) ;
      
      // Compute mean distance
      double rnpnts = 1./double(npnts) ;
      double dist = max(max(d1*rnpnts,sqrt(d1*d2*rnpnts)),
		      pow(d1*d2*d3*rnpnts,0.333333333)) ;

      // Over estimate to try to get the distance required to find stencil pnts
      //      dist *= 10 ;
      dist *= 5 ;
      dist = max(dist,delta) ;
      delta = dist ;
    }
#ifdef VERBOSE
    debugout << "BoundingBox2: dx=" << bnd.maxc[0]-bnd.minc[0]
                   << ",dy=" << bnd.maxc[1]-bnd.minc[1]
                   << ",dz=" << bnd.maxc[2]-bnd.minc[2]
                   << ",delta="<< delta << endl ;
#endif
    for(int d=0;d<3;++d) {
      bnd.maxc[d] += delta ;
      bnd.minc[d] -= delta ;
    }
    int npnts2 = collectPointsSizes(kd,bnd) ;
    
    if(npnts2 > 1000000) {
      double f = pow(double(npnts)/double(npnts2+1),0.33333) ;
      
#ifdef VERBOSE
      debugout << "npnts2="<<npnts2 <<",f=" << f << endl ;
#endif
      for(int d=0;d<3;++d) {
        bnd.maxc[d] -= (1.-f)*(delta-deltain) ;
        bnd.minc[d] += (1.-f)*(delta-deltain) ;
      }
    }
    double max_delta = 0 ;
    for(int d=0;d<3;++d) {
      max_delta = max(max_delta,double(bnd.maxc[d]-bnd.minc[d])) ;
    }
    max_delta *= 1e-3 ;
    for(int d=0;d<3;++d) {
      if(bnd.maxc[d]-bnd.minc[d] < max_delta) {
        bnd.maxc[d] += max_delta ;
        bnd.minc[d] -= max_delta ;
      }
    }
    
  }

  void getStencilBoundingBox(kdTree::KDTree<float>::bounds &bnd,
                             double &delta,
                             const const_store<vector3d<double> > &pnts,
                             entitySet dom) {
    MPI_Comm comm = get_exec_comm() ;
    //    MEMORY_PROFILE(getStencilBoundingBoxBegin) ;
    for(int d=0;d<3;++d) {
      bnd.minc[d] = .25*std::numeric_limits<float>::max() ;
      bnd.maxc[d] = -.25*std::numeric_limits<float>::max() ;
    }
    FORALL(dom,cc) {
      bnd.maxc[0] = max(bnd.maxc[0],realToFloat(pnts[cc].x)) ;
      bnd.maxc[1] = max(bnd.maxc[1],realToFloat(pnts[cc].y)) ;
      bnd.maxc[2] = max(bnd.maxc[2],realToFloat(pnts[cc].z)) ;
      bnd.minc[0] = min(bnd.minc[0],realToFloat(pnts[cc].x)) ;
      bnd.minc[1] = min(bnd.minc[1],realToFloat(pnts[cc].y)) ;
      bnd.minc[2] = min(bnd.minc[2],realToFloat(pnts[cc].z)) ;
    } ENDFORALL ;

    kdTree::KDTree<float>::bounds bndall ;
    MPI_Allreduce(&bnd.maxc[0],&bndall.maxc[0],3,MPI_FLOAT,MPI_MAX,
                  comm) ;
    MPI_Allreduce(&bnd.minc[0],&bndall.minc[0],3,MPI_FLOAT,MPI_MIN,
                  comm) ;
    int domsize = dom.size() ;
    int npnts = 0 ;
    MPI_Allreduce(&domsize,&npnts,1,MPI_INT,MPI_SUM,comm) ;
    // Compute a delta to expand the bounding box
    double d1 = bndall.maxc[0]-bndall.minc[0] ;
    double d2 = bndall.maxc[1]-bndall.minc[1] ;
    double d3 = bndall.maxc[2]-bndall.minc[2] ;
    if(d1<d2)
      std::swap(d1,d2) ;
    if(d1<d3)
      std::swap(d1,d3) ;
    if(d2<d3)
      std::swap(d2,d3) ;

    // Compute mean distance
    double rnpnts = 1./double(npnts) ;
    double dist = max(max(d1*rnpnts,sqrt(d1*d2*rnpnts)),
                    pow(d1*d2*d3*rnpnts,0.333333333)) ;

    // Over estimate to try to get the distance required to find stencil pnts
    //    dist *= 5 ;
    dist *= 2.5 ;
    dist = max(dist,delta) ;
    delta = dist ;
    for(int d=0;d<3;++d) {
      bnd.maxc[d] += dist ;
      bnd.minc[d] -= dist ;
    }
    //    MEMORY_PROFILE(getStencilBoundingBoxEnd) ;
  }

  void collectPoints(vector<kdTree::KDTree<float>::coord_info> &pout,
                     const kdTree::KDTree<float> &kd,
                     kdTree::KDTree<float>::bounds bnd) {
    MPI_Comm comm = get_exec_comm() ;
    MEMORY_PROFILE(collectPointsBegin) ;
    // Communicate bounds request to other processors
    using namespace kdTree ;
    int p = get_exec_size() ;
    vector<KDTree<float>::bounds> bnd_req(p) ;
    MPI_Allgather(&bnd,6,MPI_FLOAT,&bnd_req[0],6,MPI_FLOAT,comm) ;


    // Now communicate points that processors need to build stencil
    vector<int> scounts(p,0) ;
    vector<KDTree<float>::coord_info> pntlist ;

    for(int i=0;i<p;++i) { // Find points in i'th processor bounding box
      int bsz = pntlist.size() ;
      kd.find_box(pntlist,bnd_req[i]) ;
      scounts[i] = pntlist.size()-bsz ;
    }

    for(size_t i=0;i<scounts.size();++i)
      scounts[i]*=sizeof(KDTree<float>::coord_info) ;

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;

    vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,comm) ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }

    int result_size = (rdispls[p-1]+rcounts[p-1])/sizeof(KDTree<float>::coord_info) ;
#ifdef VERBOSE
    debugout << "results size = " << result_size << endl ;
#endif
    vector<KDTree<float>::coord_info> pcollect(result_size) ;

    MPI_Alltoallv(&pntlist[0],&scounts[0],&sdispls[0],MPI_BYTE,
                  &pcollect[0],&rcounts[0],&rdispls[0],MPI_BYTE,
                  comm) ;
    MEMORY_PROFILE(collectPointsEnd) ;
    pout.swap(pcollect) ;
  }

  void getCommSchedFromStencil(vector<int> &send_info_out,
                               vector<int> &req_sizes_out,
                               vector<int> &snd_sizes_out,
                               vector<int> &access_out,
                               const vector<Array<int,4> > &stencil,
                               const store<int> &ids,
                               const vector<int> &distribution) {
    MPI_Comm comm = get_exec_comm() ;
    MEMORY_PROFILE("pre_get_comm_sched") ;
    vector<int> stmp(stencil.size()*4) ;
    for(size_t i=0;i<stencil.size();++i)
      for(int j=0;j<4;++j)
	stmp[i*4+j] = stencil[i][j] ;

    std::sort(stmp.begin(),stmp.end()) ;
    vector<int>::const_iterator se = std::unique(stmp.begin(),stmp.end()) ;
    vector<int> access(se-stmp.begin()) ;
    
    int cnt = 0 ;

    for(vector<int>::const_iterator ii=stmp.begin();ii!=se;++ii) {
      access[cnt++] = ids[*ii] ;
      WARN(ids[*ii] < 0) ;
    }
    
    std::sort(access.begin(),access.end()) ;

    WARN(access.size()>0 && access[0] < 0) ;

    const int p = get_exec_size() ;
    // Now communicate the accessed info
    vector<int> req_sizes(p,0) ;
    // Count accesses to each processor
    cnt = 0 ;
    for(size_t i=0;i<access.size();++i) {
      while(access[i] >= distribution[cnt+1]&& cnt<p)
        cnt++ ;
      req_sizes[cnt]++ ;
    }

    vector<int> snd_sizes(p,0) ;
    MPI_Alltoall(&req_sizes[0],1,MPI_INT,&snd_sizes[0],1,MPI_INT,
                 comm) ;

    int snd_tot_size = 0 ;
    for(int i=0;i<p;++i)
      snd_tot_size += snd_sizes[i] ;
    vector<int> send_info(snd_tot_size) ;

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+req_sizes[i-1] ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+snd_sizes[i-1] ;
    }
    MPI_Alltoallv(&access[0],&req_sizes[0],&sdispls[0],MPI_INT,
                  &send_info[0],&snd_sizes[0],&rdispls[0],MPI_INT,
                  comm) ;

    MEMORY_PROFILE(recv_access) ;
    send_info_out.swap(send_info) ;
    req_sizes_out.swap(req_sizes) ;
    snd_sizes_out.swap(snd_sizes) ;
    access_out.swap(access) ;
    MEMORY_PROFILE(leaving_comm_stencil) ;
  }

  void remapStencil(vector<Array<int,4> > &stencils,
                    const vector<int> &access,
                    const store<int> &ids) {
    MEMORY_PROFILE(remapStencilStart) ;
    if(stencils.size() == 0)
      return ;
    entitySet locdom = ids.domain() ;

    vector<pair<int,int> > idmap(locdom.size()) ;

    FORALL(locdom,ii) {
      idmap[ii].first = ids[ii] ;
      idmap[ii].second = ii ;
    } ENDFORALL ;

    std::sort(idmap.begin(),idmap.end()) ;
    vector<int> acmap(locdom.size(),-1) ;
    size_t cnt = 0 ;
    for(size_t i=0;i<locdom.size();++i) {
      if(cnt >= access.size())
        break ;
      if(access[cnt] == idmap[i].first)
        acmap[idmap[i].second] = cnt++ ;
    }
    for(size_t i=0;i<stencils.size();++i) {
      for(int j=0;j<4;++j) {
        WARN(stencils[i][j] < 0 || size_t(stencils[i][j]) >=locdom.size()) ;
	WARN(acmap[stencils[i][j]] < 0) ;
	stencils[i][j] = acmap[stencils[i][j]] ;
      }
    }
    //    MEMORY_PROFILE(remapStencilEnd) ;
  }

#define COUNT_SIZE 100
#define SMALLEST_SPLIT 2048
#define STOP_SPLIT 10240
  
  inline double split_histogram(const int counts[COUNT_SIZE]) {
    int tot = counts[0] ;
    int mxcount = counts[0] ;
    int mxcid = 0 ;
    for(int i=1;i<COUNT_SIZE;++i) {
      tot += counts[i] ;
      if(mxcount < counts[i]) {
        mxcid = i ;
        mxcount = counts[i] ;
      }
    }
    int mean = tot/COUNT_SIZE ;

    // If there are zero counts in the histogram, split here.
    int tot2 = counts[0];
    int maxdist = 0 ;
    int sp = -1 ;
    for(int i=1;i<COUNT_SIZE;++i) {
      tot2 += counts[i] ;
      if(counts[i] == 0)
        if(maxdist < min(tot2,tot-tot2)) {
          sp = i ;
          maxdist = max(maxdist,min(tot2,tot-tot2)) ;
        }
    }
    if(maxdist > SMALLEST_SPLIT)
      return (double(sp)+.5)/double(COUNT_SIZE) ;

    if(10*mxcount < 12*mean)  // If less than 20% variation in histogram
      return .5 ;             // split down the middle.
    
    int s1 ;
    for(s1=mxcid-1;s1>=0;--s1)
      if(counts[s1] < mean)
        break ;
    int s2 ;
    for(s2=mxcid+1;s2<COUNT_SIZE;++s2)
      if(counts[s2] < mean)
        break ;
    
    int c1 = 0 ;
    int c2 = 0 ;
    for(int i=s1;i>=0;--i)
      if(counts[i] < mean)
        c1++ ;
    for(int i=s2;i<COUNT_SIZE;++i)
      if(counts[i] < mean)
        c2++ ;

    // find split loc
    int cx = 0 ;
    if(c1 > c2)
      cx = s1+1 ;
    else
      cx = s2 ;
    // Check to see if enough points to split.
    int sum = 0 ;
    for(int i=0;i<cx;++i)
      sum += counts[i] ;
    if(min(sum,tot-sum) < SMALLEST_SPLIT)
      return 0.5 ;
    return double(cx)/double(COUNT_SIZE) ;
  }

  void histogram_part(vect3d vecs[], int ids[],int start,int end,int depth,vector<int> &sizes,int &splits) {
    if(depth == 0 ||splits == 0) {
      sizes.push_back(end-start) ;
      return ;
    }
    if(start == end) {
      return ;
    }
    if(end-start < STOP_SPLIT) {
      sizes.push_back(end-start) ;
      return ;
    }
      
    int counts[COUNT_SIZE] ;
    for(int i=0;i<COUNT_SIZE;++i)
      counts[i] = 0 ;
    vect3d mx=vecs[start],mn=vecs[start] ;
    for(int i=start+1;i<end;++i) {
      mx.x = max(mx.x,vecs[i].x) ;
      mx.y = max(mx.y,vecs[i].y) ;
      mx.z = max(mx.z,vecs[i].z) ;
      mn.x = min(mn.x,vecs[i].x) ;
      mn.y = min(mn.y,vecs[i].y) ;
      mn.z = min(mn.z,vecs[i].z) ;
    }
    double dx = realToDouble(mx.x-mn.x) ;
    double dy = realToDouble(mx.y-mn.y) ;
    double dz = realToDouble(mx.z-mn.z) ;
    int top = start ;
    int bot = end-1 ;

    if(dx > dy && dx > dz) { // x coord split
      for(int i=start;i<end;++i) {
        double t = realToDouble(vecs[i].x-mn.x)/dx ;
        int ind = max(min(int(floor(t*COUNT_SIZE)),COUNT_SIZE-1),0) ;
        counts[ind]++ ;
      }
      double t = split_histogram(counts) ;

      double xs = realToDouble(mn.x) + t*dx ;
      while(top <= bot) {
        if(vecs[top].x > xs) {
          std::swap(vecs[top],vecs[bot]) ;
          std::swap(ids[top],ids[bot]) ;
          bot-- ;
        } else
          top++ ;
      }
    } else if(dy > dz) { // y coord split
      for(int i=start;i<end;++i) {
        double t = realToDouble(vecs[i].y-mn.y)/dy ;
        int ind = max(min(int(floor(t*COUNT_SIZE)),COUNT_SIZE-1),0) ;
        counts[ind]++ ;
      }

      double t = split_histogram(counts) ;

      double ys = realToDouble(mn.y) + t*dy ;
      while(top <= bot) {
        if(vecs[top].y > ys) {
          std::swap(vecs[top],vecs[bot]) ;
          std::swap(ids[top],ids[bot]) ;
          bot-- ;
        } else
          top++ ;
      }
    } else {      // z coord split
      for(int i=start;i<end;++i) {
        double t = realToDouble(vecs[i].z-mn.z)/dz ;
        int ind = max(min(int(floor(t*COUNT_SIZE)),COUNT_SIZE-1),0) ;
        counts[ind]++ ;
      }

      double t = split_histogram(counts) ;

      double zs = realToDouble(mn.z) + t*dz ;
      while(top <= bot) {
        if(vecs[top].z > zs) {
          std::swap(vecs[top],vecs[bot]) ;
          std::swap(ids[top],ids[bot]) ;
          bot-- ;
        } else
          top++ ;
      }
    }
    top = min(top,end) ;
    splits-- ;
    if(min(top-start,end-top) < SMALLEST_SPLIT) {
      sizes.push_back(end-start) ; // If split ineffective, just return chunk
    } else {
      histogram_part(vecs,ids,start,top,depth-1,sizes,splits) ;
      histogram_part(vecs,ids,top,end,depth-1,sizes,splits) ;
    }
  }


  inline bool vpairxcmp(const vpair &x1, const vpair &x2) {
    return x1.first.x < x2.first.x ;
  }
  inline bool vpairycmp(const vpair &x1, const vpair &x2) {
    return x1.first.y < x2.first.y ;
  }
  inline bool vpairzcmp(const vpair &x1, const vpair &x2) {
    return x1.first.z < x2.first.z ;
  }

  // Group points into a collection of bounding boxes (no more than 2^levels)
  // For each bounding box, find the smallest reciprocal point spacing
  // and assign this spacing to the centroid of the bounding box.
  void collectGroups(vector<vpair> &results,
                     vector<vpair> &inputs, int start, int end,
                     int levels, int min_pnts) {
    float max_x=-1e33,min_x=1e33,max_y=-1e33,min_y=1e33,max_z=-1e33,min_z=1e33;
    for(int i=start;i<end;++i) {
      if(inputs[i].second > 0) {
        max_x = max(max_x,inputs[i].first.x) ;
        min_x = min(min_x,inputs[i].first.x) ;
        max_y = max(max_y,inputs[i].first.y) ;
        min_y = min(min_y,inputs[i].first.y) ;
        max_z = max(max_z,inputs[i].first.z) ;
        min_z = min(min_z,inputs[i].first.z) ;
      }
    }
    if(levels == 0 || end-start < min_pnts) {
      vpair v ;
      v.first = .5*vector3d<float>(max_x+min_x,max_y+min_y,max_z+min_z) ;
      int min_space = 2100000000 ;
      for(int i=start;i<end;++i)
        if(inputs[i].second > 0)
          min_space = min(min_space,inputs[i].second) ;
      v.second = min_space ;
      results.push_back(v) ;
      return ;
    }
    levels-- ;
    float dx = max_x-min_x ;
    float dy = max_y-min_y ;
    float dz = max_z-min_z ;
    if(dx > dy && dx > dz) {
      std::sort(&inputs[start],&inputs[end],vpairxcmp) ;
    } else if(dy > dz) {
      std::sort(&inputs[start],&inputs[end],vpairycmp) ;
    } else {
      std::sort(&inputs[start],&inputs[end],vpairzcmp) ;
    }

    // Split space into 4 equal sized parts
    int dinc = (end-start)/4 ;
    for(int i=0;i<3;++i) {
      collectGroups(results,inputs,start,start+dinc,levels,min_pnts) ;
      start += dinc ;
    }
    collectGroups(results,inputs,start,end,levels,min_pnts) ;
  }

  // Input, pnts are the points to be decomposed, this routine sorts these
  // points into bouding box regions
  // boxes are the list of bounding boxes created
  // split is the value to split this level of the tree
  // dim is the dimension of the split (0,1,2 -- x,y,z)
  // start and stop are the indicies into pnts that this call will decompose
  void getBoundingBoxDecomp(vector<kdTree::KDTree<float>::coord_info> &pnts, 
                            vector<bound_info> &boxes, 
                            double split,
                            int dim,
                            int levels,
                            int start,
                            int stop,
                            double delta_lim,
                            int split_lim) {
#ifdef VERBOSE
    debugout << "bbox decomp, dim = " << dim << ",split = " << split << ", levels = " << levels << endl ;
#endif
    kdTree::KDTree<float>::bounds bndl,bndr ;
    for(int d=0;d<3;++d) {
      bndl.minc[d] = .25*std::numeric_limits<float>::max() ;
      bndl.maxc[d] = -.25*std::numeric_limits<float>::max() ;
      bndr.minc[d] = .25*std::numeric_limits<float>::max() ;
      bndr.maxc[d] = -.25*std::numeric_limits<float>::max() ;
    }
    int i=start,j=stop ;
    while (i<=j) {
      if(pnts[i].coords[dim] < split) {
        for(int d=0;d<3;++d) {
          bndl.minc[d] = min(bndl.minc[d],pnts[i].coords[d]) ;
          bndl.maxc[d] = max(bndl.maxc[d],pnts[i].coords[d]) ;
        }
        i++ ;
      } else {
        for(int d=0;d<3;++d) {
          bndr.minc[d] = min(bndr.minc[d],pnts[i].coords[d]) ;
          bndr.maxc[d] = max(bndr.maxc[d],pnts[i].coords[d]) ;
        }
        std::swap(pnts[i],pnts[j]) ;
        j-- ;
      }
    }
    // If number of points to small, why subdivide?
    const int threshold = 16 ;
    // process left branch
    int nleft = i-start;
    double delta = bndl.maxc[0]-bndl.minc[0] ;
    int ndim = 0 ;
    double nsplit = bndl.minc[0] + 0.5*delta ;
    double delta1 = bndl.maxc[1]-bndl.minc[1] ;
    if(delta1 > delta) {
      delta = delta1 ;
      ndim = 1 ;
      nsplit = bndl.minc[1] + 0.5*delta1 ;
    }
    double delta2 = bndl.maxc[2]-bndl.minc[2] ;
    if(delta2 > delta) {
      delta = delta2 ;
      ndim=2 ;
      nsplit = bndl.minc[2] + 0.5*delta2 ;
    }
    if(levels > 1 && nleft > threshold &&
       (delta > delta_lim || nleft > split_lim) ) { // split further
      getBoundingBoxDecomp(pnts,boxes,nsplit,ndim,levels-1,start,i-1,
                           delta_lim,split_lim) ;
    } else {
      if(nleft != 0) {
        bound_info bi ;
        bi.bnd = bndl ;
        bi.start = start ;
        bi.stop = i-1 ;
        boxes.push_back(bi) ;
      }
    }
    // process right branch
    int nright = stop-i+1 ;
    delta = bndr.maxc[0]-bndr.minc[0] ;
    ndim = 0 ;
    nsplit = bndr.minc[0] + 0.5*delta ;
    delta1 = bndr.maxc[1]-bndr.minc[1] ;
    if(delta1 > delta) {
      delta = delta1 ;
      ndim = 1 ;
      nsplit = bndr.minc[1] + 0.5*delta1 ;
    }
    delta2 = bndr.maxc[2]-bndr.minc[2] ;
    if(delta2 > delta) {
      delta = delta2 ;
      ndim=2 ;
      nsplit = bndr.minc[2] + 0.5*delta2 ;
    }
    if(levels > 1 && nright > threshold &&
       (delta > delta_lim || nright > split_lim) ) { // split further
      getBoundingBoxDecomp(pnts,boxes,nsplit,ndim,levels-1,i,stop,
                           delta_lim,split_lim) ;
    } else {
      if(nright != 0) {
        bound_info bi ;
        bi.bnd = bndr ;
        bi.start = i ;
        bi.stop = stop ;
        boxes.push_back(bi) ;
      }
    }
  }

  void getBoundingBoxes(vector<kdTree::KDTree<float>::coord_info> &pnts, 
                        vector<bound_info> &boxes, int levels, double delta_lim,
                        int split_lim) {
    boxes.clear() ;
    int sz = pnts.size() ;
    if(sz == 0)
      return ;
    kdTree::KDTree<float>::bounds bnd ;
    bnd.minc[0] = pnts[0].coords[0] ;
    bnd.maxc[0] = pnts[0].coords[0] ;
    bnd.minc[1] = pnts[0].coords[1] ;
    bnd.maxc[1] = pnts[0].coords[1] ;
    bnd.minc[2] = pnts[0].coords[2] ;
    bnd.maxc[2] = pnts[0].coords[2] ;

    // Sample to get first split.
    int bdelta = (sz/20)+1 ;
    for(int i=1;i<sz;i+=bdelta) {
      bnd.minc[0] = min(bnd.minc[0],pnts[i].coords[0]) ;
      bnd.minc[1] = min(bnd.minc[1],pnts[i].coords[1]) ;
      bnd.minc[2] = min(bnd.minc[2],pnts[i].coords[2]) ;
      bnd.maxc[0] = max(bnd.maxc[0],pnts[i].coords[0]) ;
      bnd.maxc[1] = max(bnd.maxc[1],pnts[i].coords[1]) ;
      bnd.maxc[2] = max(bnd.maxc[2],pnts[i].coords[2]) ;
    }

    int ndim = 0 ;
    double nsplit = 0.5*(bnd.minc[0]+bnd.maxc[0]) ;
    double delta = bnd.maxc[0]-bnd.minc[0] ;
    double delta1 = bnd.maxc[1]-bnd.minc[1] ;
    double delta2 = bnd.maxc[2]-bnd.minc[2] ;
    if(delta < delta1) {
      delta = delta1 ;
      nsplit = 0.5*(bnd.minc[1]+bnd.maxc[1]) ;
      ndim = 1 ;
    }
    if(delta < delta2) {
      delta = delta2 ;
      nsplit = 0.5*(bnd.minc[2]+bnd.maxc[2]) ;
      ndim = 2 ;
    }
    // Check bounds, if under limit, decompose, otherwise return one box.
    if(delta > delta_lim || sz > split_lim)
      getBoundingBoxDecomp(pnts,boxes,nsplit,ndim,levels,0,sz-1,delta_lim, split_lim) ;
    else {
      bound_info bi ;
      bi.bnd = bnd ;
      bi.start = 0 ;
      bi.stop = sz-1 ;
      boxes.push_back(bi) ;
    }
  }


  // Define bounding box distribution where boxes are numbered over the
  // concatenation of all boxes.
  // Output:
  //  box_sp: source processsor of a box
  //  box_tp: target processor of a box
  void boundingBoxDistribution(vector<int> &box_sp,
                               vector<int> &box_tp,
                               vector<int> &b_sizes,
                               vector<int> &send_block_id,
                               vector<int> &recv_block_id,
                               const vector<bound_info> &boxes,
                               const MPI_Comm &comm) {

    int p ;
    MPI_Comm_size(comm,&p) ;
    int r ;
    MPI_Comm_rank(comm,&r) ;
    
    int bsize = boxes.size() ;
    
    tmp_array<int> bsize_list(p)  ;
    MPI_Allgather(&bsize,1,MPI_INT,&bsize_list[0],1,MPI_INT,comm) ;

#ifdef VERBOSE
    debugout << "bsizes = " << endl ;
    for(int i=0;i<p;++i)
      debugout << ' ' << bsize_list[i] ;
    debugout  << endl ;
#endif

    vector<int> lb_sizes(bsize) ;
    for(int i=0;i<bsize;++i) {
      lb_sizes[i] = boxes[i].stop-boxes[i].start+1 ;
    }

    allGatherVec(b_sizes,lb_sizes,comm) ;

    vector<int> box_sp_tmp(b_sizes.size()) ; // source processor of block
    box_sp.swap(box_sp_tmp) ;
    
    int cnt = 0 ;
    for(int i=0;i<p;++i) 
      for(int j=0;j<bsize_list[i];++j)
        box_sp[cnt++] = i ;

    
    // Allocate target processors
    int tot = 0 ;
    int threshold = 0 ;
    for(size_t i=0;i<b_sizes.size();++i) {
      tot += b_sizes[i] ;
      //    threshold = max(threshold,b_sizes[i]+1) ;
    }
    
    vector<int> box_tp_tmp(b_sizes.size()) ; // target processor of block
    box_tp.swap(box_tp_tmp) ;
    
    threshold = max(threshold,(tot / p) + 1) ;
    //  threshold += threshold/16 ;

#ifdef VERBOSE
    debugout << "threshold = " << threshold << endl ;
#endif
    vector<int> psizes(p,0) ;
    cnt = 0 ;
    vector<pair<int,int> > block_list_sort(b_sizes.size()) ;
    for(size_t i=0;i<b_sizes.size();++i) {
      block_list_sort[i].first = b_sizes[i] ;
      block_list_sort[i].second = i ;
      box_tp[i] = -1 ;
    }
   
    std::sort(block_list_sort.begin(),block_list_sort.end()) ;
    // First assign largest blocks to owning processor if less than
    // threshold
    for(size_t i=0;i<b_sizes.size();++i) {
      int bk = block_list_sort[b_sizes.size()-i-1].second ;
      if(psizes[box_sp[bk]] == 0 ||
         psizes[box_sp[bk]] < threshold) {
        box_tp[bk] = box_sp[bk] ;
        psizes[box_sp[bk]] += b_sizes[bk] ;
      }
    }
    
    vector<pair<int,int> > vpsizes(p) ;
    for(int i=0;i<p;++i) {
      vpsizes[i].first = psizes[i] ;
      vpsizes[i].second = i ;
    }
    std::sort(vpsizes.begin(),vpsizes.end()) ;
    int pcnt = 0 ;

    int mode = 1 ;
    for(size_t i=0;i<b_sizes.size();++i) {
      int bk = block_list_sort[b_sizes.size()-i-1].second ;
      if(box_tp[bk] < 0) { // unassigned
        while(psizes[vpsizes[pcnt].second]  + b_sizes[bk]*mode > threshold) {
          pcnt++ ;
          if(pcnt == p) {
            pcnt = 0 ;
            for(int i=0;i<p;++i) {
              vpsizes[i].first = psizes[i] ;
              vpsizes[i].second = i ;
            }
            mode = 0 ;
            std::sort(vpsizes.begin(),vpsizes.end()) ;
          }
        }
        box_tp[bk] = vpsizes[pcnt].second ;
        psizes[vpsizes[pcnt].second] += b_sizes[bk] ;
        //        pcnt++ ;
      }
    }
    
#ifdef VERBOSE
    debugout << "psizes =" ;
    for(int i=0;i<p;++i)
      debugout << ' ' << psizes[i] ;
    debugout << endl ;
#endif

#ifdef VERBOSE        
    debugout << "block sizes = " ;
    for(size_t i=0;i<b_sizes.size();++i) {
      debugout << ' ' << b_sizes[i]
                     << '[' << box_sp[i] << '>' << box_tp[i] <<']' ;
    }
    debugout << endl ;
#endif
    recv_block_id.clear() ;
    send_block_id.clear() ;
    for(size_t i=0;i<b_sizes.size();++i) {
      if(box_sp[i] == r) {
        send_block_id.push_back(i) ;
#ifdef VERBOSE
        debugout << "sending block " << i << " to " << box_tp[i] << endl ;
#endif
      }
      if(box_tp[i] == r) {
        recv_block_id.push_back(i) ;
#ifdef VERBOSE
        debugout << "recving block " << i << " from " << box_sp[i] << endl ;
#endif
      }
    }

  }                             

  using namespace kdTree ;

  void recieveTargetPoints(vector<vector<KDTree<float>::coord_info>  > &recv_targets,
                           const kdTree::KDTree<float> *kd,
                           const vector<bound_info> &boxes_g,
                           const vector<int> &box_sp,
                           const vector<int> &box_tp,
                           const vector<int> &b_sizes,
                           const vector<int> &send_block_id,
                           const vector<int> &recv_block_id,
                           const MPI_Comm &comm) {
    
    MEMORY_PROFILE(RecieveTargetPointsBegin) ;
    int p ;
    MPI_Comm_size(comm,&p) ;
    int r ;
    MPI_Comm_rank(comm,&r) ;

    // loop over boxes and collect target points to send
    vector<KDTree<float>::coord_info> pntlist ;
    vector<int> scounts(boxes_g.size(),0) ;
    vector<int> block_sends(p,0) ;
    vector<vector<int> > proc_send_list(p) ;
    
    for(size_t i=0;i<boxes_g.size();++i) {
      int bsz = pntlist.size() ;
      (kd)->find_box(pntlist,boxes_g[i].bnd) ;
      scounts[i] = pntlist.size()-bsz ;
      
#ifdef VERBOSE
      if(scounts[i] != 0) {

        debugout << "sending " << scounts[i] << " points to " << box_tp[i] << endl ;
        kdTree::KDTree<float>::bounds btmp ;

        for(int d=0;d<3;++d) {
          btmp.minc[d] = pntlist[bsz].coords[d] ;
          btmp.maxc[d] = pntlist[bsz].coords[d] ;
        }
        for(int k=0;k<scounts[i];++k) {
          for(int d=0;d<3;++d) {
            btmp.minc[d] = min(btmp.minc[d],pntlist[bsz+k].coords[d]) ;
            btmp.maxc[d] = max(btmp.minc[d],pntlist[bsz+k].coords[d]) ;
          }
        }
          
        debugout << "boxesg="
                       << boxes_g[i].bnd.minc[0] << ","
                       <<boxes_g[i].bnd.maxc[0] << " "
                       << boxes_g[i].bnd.minc[1] << ","
                       <<boxes_g[i].bnd.maxc[1] << " "
                       << boxes_g[i].bnd.minc[2] << ","
                       <<boxes_g[i].bnd.maxc[2]
                       << endl ;
        debugout << "sent="
                       << btmp.minc[0] << "," <<btmp.maxc[0] << " "
                       << btmp.minc[1] << "," <<btmp.maxc[1] << " "
                       << btmp.minc[2] << "," <<btmp.maxc[2]
                       << endl ;
        
      }
#endif
      if(scounts[i] != 0) {
        FATAL(box_tp[i] < 0 || box_tp[i] >= p) ;
        proc_send_list[box_tp[i]].push_back(i) ;
        block_sends[box_tp[i]]++ ;
      }
    }

    vector<int> soffsets(boxes_g.size(),0) ;
    for(size_t i=1;i<boxes_g.size();++i) {
      soffsets[i] = soffsets[i-1]+scounts[i-1] ;
    }
    // Now we know what to send, we need to send it.

    // find out how many blocks we are recieving and from which processor
    vector<int> block_recvs(p,0) ;
    MPI_Alltoall(&block_sends[0],1,MPI_INT,&block_recvs[0],1,MPI_INT, comm) ;

    // Now recieve the block sizes
    int num_recv = 0 ;
    for(int i=0;i<p;++i) 
      num_recv += block_recvs[i] ;
    
    vector<pair<int,int> > recv_info(num_recv) ;
    vector<MPI_Request> req_queue(num_recv) ;

    int cnt = 0 ;
    for(int i=0;i<p;++i) 
      for(int j=0;j<block_recvs[i];++j) {
        MPI_Irecv(&(recv_info[cnt]),2,MPI_INT,i,j,comm,&req_queue[cnt]) ;
        cnt++ ;
      }

    const int p1 = 3917 ; // two primes to help randomize order of
    const int p2 = 4093 ; // sending to reduce bottlenecks
    for(int i=0;i<p;++i) {
      int cp = (i*p1+r*p2)%p ; // processor to send data this iteration

      FATAL(cp < 0 || cp >=p)  ;
      for(size_t j=0;j<proc_send_list[cp].size();++j) {
        pair<int,int> buf ;
        buf.first = proc_send_list[cp][j] ;
        FATAL(buf.first < 0 || buf.first >= int(scounts.size())) ;
        buf.second = scounts[buf.first] ;

        MPI_Send(&(buf),2,MPI_INT,cp,j,comm) ;
      }
    }
#ifdef VERBOSE
    debugout << "num_recv = " << num_recv << endl ;
#endif
    if(num_recv != 0) {
      vector<MPI_Status> stat_queue(num_recv) ;
      MPI_Waitall(num_recv,&req_queue[0],&stat_queue[0]) ;
    }

    // Add up recieve sizes
    vector<int> recv_tot_size(boxes_g.size(),0) ;
    cnt = 0 ;
    for(int i=0;i<p;++i) {
      for(int j=0;j<block_recvs[i];++j) {
	recv_tot_size[recv_info[cnt].first] += recv_info[cnt].second ;
        cnt++ ;
      }
    }

#ifdef VERBOSE
    debugout << "recv_info:" << endl ;
    cnt = 0 ;
    for(int i=0;i<p;++i) {
      debugout  << "p=" << i << "," ;
      for(int j=0;j<block_recvs[i];++j) {
        debugout << ' ' << recv_info[cnt].first << ',' << recv_info[cnt].second ;

        cnt++ ;
      }
      debugout << endl ;
    }
#endif

    {
      vector<vector<kdTree::KDTree<float>::coord_info>  > tmpx(boxes_g.size()) ;
      recv_targets.swap(tmpx) ;
    }
    MEMORY_PROFILE(RecieveTargetPointsMid) ;

    // allocate points to receive 
    for(size_t i=0;i<recv_tot_size.size();++i) {
      if(recv_tot_size[i] != 0) {
	vector<kdTree::KDTree<float>::coord_info>  tmp(recv_tot_size[i]) ;
	recv_targets[i].swap(tmp) ;
      }
    }
    
    vector<int> recv_offset(recv_tot_size.size(),0) ;
    cnt = 0 ;
    int coord_size = sizeof(kdTree::KDTree<float>::coord_info) ;
    for(int i=0;i<p;++i) 
      for(int j=0;j<block_recvs[i];++j) {
	int bk = recv_info[cnt].first ;
	int offset = recv_offset[bk] ;
	int recv_size = recv_info[cnt].second ;
	recv_offset[bk] += recv_size ;
	MPI_Irecv(&recv_targets[bk][offset],recv_size*coord_size,MPI_BYTE,i,j,
                  comm,&req_queue[cnt]) ;
        cnt++ ;
      }
    
    for(int i=0;i<p;++i) {
      int cp = (i*p1+r*p2)%p ; // processor to send data this iteration
      for(size_t j=0;j<proc_send_list[cp].size();++j) {
        int bk = proc_send_list[cp][j] ;
        int bsz = scounts[bk] ;
        MPI_Send(&pntlist[soffsets[bk]],bsz*coord_size,MPI_BYTE,cp,j,comm) ;
      }
    }

#ifdef VERBOSE
    debugout << "num_recv2 = " << num_recv << endl ;
    debugout << "cnt=" << cnt << endl ;
#endif
    if(num_recv != 0) {
      vector<MPI_Status> stat_queue(num_recv) ;
      MPI_Waitall(num_recv,&req_queue[0],&stat_queue[0]) ;
    }
    MEMORY_PROFILE(RecieveTargetPointsEnd) ;
  }

}
