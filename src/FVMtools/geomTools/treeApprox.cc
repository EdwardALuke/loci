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
#include <utility>
#include "gridMotion/treeApprox.h"
#include <iostream>
using std::cout ;
using std::endl ;
#include <sstream>
#include <string>
using std::string ;
//#include <Loci>

namespace gridMotion {
  
  static void (*reporterFunction)(string s) = 0 ;

  void registerReporterFunction(void (*rfp)(std::string )) {
    reporterFunction = rfp ;
  }
  
  void reportTime(const char *text, double time) {
    if(reporterFunction != 0) {
      std::ostringstream oss ;
      oss << "Time for " << text << ": " << time << endl ;
      string out = oss.str() ;
      (*reporterFunction)(out) ;
    }
  }
  
  using std::min ;
  using std::max ;
//--------------------------------------------------------------------  
inline double frsum_basic(double x[3],double beta[9],double xc[3], int a) {
  double rr1 = sqrt((x[0]-beta[0])*(x[0]-beta[0])+
                    (x[1]-beta[1])*(x[1]-beta[1])+
                    (x[2]-beta[2])*(x[2]-beta[2])) ;
  double rr2 = sqrt((x[0]-beta[3])*(x[0]-beta[3])+
                    (x[1]-beta[4])*(x[1]-beta[4])+
                    (x[2]-beta[5])*(x[2]-beta[5])) ;
  double rr3 = sqrt((x[0]-beta[6])*(x[0]-beta[6])+
                    (x[1]-beta[7])*(x[1]-beta[7])+
                    (x[2]-beta[8])*(x[2]-beta[8])) ;
  double x40 = 4.*xc[0]-beta[0]-beta[3]-beta[6] ;
  double x41 = 4.*xc[1]-beta[1]-beta[4]-beta[7] ;
  double x42 = 4.*xc[2]-beta[2]-beta[5]-beta[8] ;
  
  double rr4 = sqrt((x[0]-x40)*(x[0]-x40)+
                    (x[1]-x41)*(x[1]-x41)+
                    (x[2]-x42)*(x[2]-x42)) ;
  return 0.25*(pow(rr1,-a)+pow(rr2,-a)+ pow(rr3,-a)+pow(rr4,-a)) ;
}

  inline double frsum(double x[3],double beta[9],double xc[3],double J[9],int a) {
    double rr1 = sqrt((x[0]-beta[0])*(x[0]-beta[0])+
		     (x[1]-beta[1])*(x[1]-beta[1])+
		     (x[2]-beta[2])*(x[2]-beta[2])) ;
    double rr2 = sqrt((x[0]-beta[3])*(x[0]-beta[3])+
		     (x[1]-beta[4])*(x[1]-beta[4])+
		     (x[2]-beta[5])*(x[2]-beta[5])) ;
    double rr3 = sqrt((x[0]-beta[6])*(x[0]-beta[6])+
		     (x[1]-beta[7])*(x[1]-beta[7])+
		     (x[2]-beta[8])*(x[2]-beta[8])) ;
    double x40 = 4.*xc[0]-beta[0]-beta[3]-beta[6] ;
    double x41 = 4.*xc[1]-beta[1]-beta[4]-beta[7] ;
    double x42 = 4.*xc[2]-beta[2]-beta[5]-beta[8] ;
  
    double rr4 = sqrt((x[0]-x40)*(x[0]-x40)+
		     (x[1]-x41)*(x[1]-x41)+
		     (x[2]-x42)*(x[2]-x42)) ;
  
    double J4 = .25*pow(rr4,-a-2) ;
    double J1 = .25*pow(rr1,-a-2) ;
    J[0] = -double(a)*((beta[0]-x[0])*J1+(x[0]-x40)*J4) ;
    J[1] = -double(a)*((beta[1]-x[1])*J1+(x[1]-x41)*J4) ;
    J[2] = -double(a)*((beta[2]-x[2])*J1+(x[2]-x42)*J4) ;
    double J2 = .25*pow(rr2,-a-2) ;
    J[3] = -double(a)*((beta[3]-x[0])*J2+(x[0]-x40)*J4) ;
    J[4] = -double(a)*((beta[4]-x[1])*J2+(x[1]-x41)*J4) ;
    J[5] = -double(a)*((beta[5]-x[2])*J2+(x[2]-x42)*J4) ;
    double J3 = .25*pow(rr3,-a-2) ;
    J[6] = -double(a)*((beta[6]-x[0])*J3+(x[0]-x40)*J4) ;
    J[7] = -double(a)*((beta[7]-x[1])*J3+(x[1]-x41)*J4) ;
    J[8] = -double(a)*((beta[8]-x[2])*J3+(x[2]-x42)*J4) ;

    return 0.25*(pow(rr1,-a)+pow(rr2,-a)+ pow(rr3,-a)+pow(rr4,-a)) ;
  }

  //--------------------------------------------------------------------
  // points that sample the surface of a unit sphere
  double testpnts[38][3] = {
    // dodecahedron 
    { 0.57735027, 0.57735027, 0.57735027},
    { 0.57735027, 0.57735027,-0.57735027},
    { 0.57735027,-0.57735027, 0.57735027},
    { 0.57735027,-0.57735027,-0.57735027},
    {-0.57735027, 0.57735027, 0.57735027},
    {-0.57735027, 0.57735027,-0.57735027},
    {-0.57735027,-0.57735027, 0.57735027},
    {-0.57735027,-0.57735027,-0.57735027},
    { 0, 0.343560749722, 0.899453719973},
    { 0, 0.343560749722,-0.899453719973},
    { 0,-0.343560749722, 0.899453719973},
    { 0,-0.343560749722,-0.899453719973},
    { 0.343560749722, 0.899453719973, 0},
    { 0.343560749722,-0.899453719973, 0},
    {-0.343560749722, 0.899453719973, 0},
    {-0.343560749722,-0.899453719973, 0},
    { 0.899453719973, 0, 0.343560749722}, 
    { 0.899453719973, 0,-0.343560749722}, 
    {-0.899453719973, 0, 0.343560749722}, 
    {-0.899453719973, 0,-0.343560749722},
    // isohedron
    { 0,  0.525731112119,  0.850650808352},
    { 0, -0.525731112119,  0.850650808352},
    { 0,  0.525731112119, -0.850650808352},
    { 0, -0.525731112119, -0.850650808352},
    { 0.525731112119,  0.850650808352, 0}, 
    {-0.525731112119,  0.850650808352, 0}, 
    { 0.525731112119, -0.850650808352, 0}, 
    {-0.525731112119, -0.850650808352, 0},
    { 0.850650808352, 0,  0.525731112119},
    {-0.850650808352, 0,  0.525731112119},
    { 0.850650808352, 0, -0.525731112119},
    {-0.850650808352, 0, -0.525731112119},
    // axis
    { 1, 0, 0},
    {-1, 0, 0},
    { 0, 1, 0},
    { 0,-1, 0},
    { 0, 0, 1},
    { 0, 0,-1}
  } ;


  void get_fexact(const vector<nodeDataBase> &nodeDataCollection,
                  int start,
                  int end,
                  int exponent,
                  double radius,
                  vector3d<double> xc,
                  const double testpnts[][3],
                  double fexact[],
                  int npnts) {
  
    for(int j=0;j<npnts;++j) {
      vector3d<double> isopt(testpnts[j][0],testpnts[j][1],testpnts[j][2]) ;
      vector3d<double> pt = vector3d<double>(xc) + double(radius)*isopt ;
      double f = 0 ;
      for(int i=start;i<end;++i) {
	const vector3d<double> npos   = nodeDataCollection[i].pos;
	const double r = norm(pt-npos) ;
	const double wn = nodeDataCollection[i].weight;
	f += wn * pow(r,-exponent) ;
      }
      fexact[j] = f ;
    }
  }
  //--------------------------------------------------------------------  
  template<class T1, class T2, class T3> 
  inline void solve_lu_pivot_approx(const T1 * restrict A,
				    const T2 * restrict b,
				    T3 * restrict x,
				    const pivot_type * restrict pivot,
				    int size) {
    // Perform forward solve Ly = b, note b becomes x after this step
    for(int i=0;i<size;++i) {
      T3 xi = b[pivot[i]] ;
      const T1 * restrict Aj = A ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int j=0;j<i;++j,Aj+=size)
	xi -= Aj[i]*x[j] ;
      x[i] = xi ;
    }
    // Do back solver Ux = y
    const T1 * restrict Ai = A + size*(size-1) ;
    for(int i=size-1;i>=0;--i,Ai-=size) {
      const T1 *restrict Aj = Ai + size ;
      T3 xi = x[i] ;
#ifdef HAVE_IVDEP
#pragma ivdep
#endif
      for(int j=i+1;j<size;++j,Aj+=size)
	xi -= Aj[i]*x[j] ;
      x[i] = xi/(Ai[i]+((Ai[i] < 0)?-1e-30:1e-30)) ;
    }
  }

  void makeQuad(const vector<nodeDataBase> &nodeDataCollection,
		int start,
		int end,
		const double testpnts[][3],
		const double fexact[],
		int nfit,
		int ntest,
		int exponent,
		vector3d<float> xc,
		double radius,
		vector3d<float> &qp1,
		vector3d<float> &qp2,
		vector3d<float> &qp3,
		float &error) {
  
    if((end-start) < 6) { // If too few points abort
      qp1 = xc ;
      qp2 = qp1 ;
      qp3 = qp1 ;
      if(end-start == 2) {
	qp1 = nodeDataCollection[start].pos ;
	qp2 = nodeDataCollection[start].pos ;
	qp3 = nodeDataCollection[start+1].pos ;
      }
      if(end-start == 4) {
	qp1 = nodeDataCollection[start].pos ;
	qp2 = nodeDataCollection[start+1].pos ;
	qp3 = nodeDataCollection[start+2].pos ;
      }
      // compute error
      double maxerr = 0 ;
      for(int jj=0;jj<ntest;++jj) {
	vector3d<float> isopt(testpnts[jj][0],testpnts[jj][1],testpnts[jj][2]) ;
	vector3d<float> pt = xc + radius*isopt ;
	double rr1 = norm(qp1-pt) ;
	double rr2 = norm(qp2-pt) ;
	double rr3 = norm(qp3-pt) ;
	double rr4 = norm(4.*xc-qp1-qp2-qp3-pt) ;
	int a = exponent ;
	double ya=0.25*(pow(rr1,-a)+pow(rr2,-a)+ pow(rr3,-a)+pow(rr4,-a)) ;
	double r = fexact[jj]-ya ;
	maxerr = max(maxerr,fabs(r)/fexact[jj]) ;
      }
      //      cout << "small maxerr = " << maxerr << ", size=" << end-start << endl ;
      error = maxerr ;
      return ;
    }
  
  
    // begin non-linear iteration to find quad points for r^3
    double XC[3] = {xc.x,xc.y,xc.z} ;
    vector3d<double> xcd(xc.x,xc.y,xc.z) ;
    double beta[9] ;
    // initial beta guess
    const int middp = min(start + (end-start)/2,end-1) ;
    beta[0] = nodeDataCollection[start].pos.x ;
    beta[1] = nodeDataCollection[start].pos.y ;
    beta[2] = nodeDataCollection[start].pos.z ;
    beta[3] = nodeDataCollection[end-1].pos.x ;
    beta[4] = nodeDataCollection[end-1].pos.y ;
    beta[5] = nodeDataCollection[end-1].pos.z ;
    beta[6] = nodeDataCollection[middp].pos.x ;
    beta[7] = nodeDataCollection[middp].pos.y ;
    beta[8] = nodeDataCollection[middp].pos.z ;
    
    real fmin = 0 ;
    for(int jj=0;jj<nfit;++jj) 
      fmin = max(fmin,fexact[jj]) ;

    double J[9] ;
    double Data[81] ;
    Mat<double> ATA(Data,9) ;
    double B[9] ;
    // begin non-linear solver iteration
    for(int iter=0;iter<32;++iter) {
      for(int i=0;i<81;++i) // zero out ATA matrix
	Data[i] = 0 ;
      for(int i=0;i<9;++i)
	B[i] = 0 ;
      // Loop over fitting points
      double mxr = 0 ;
      for(int jj=0;jj<nfit;++jj) {
	vector3d<double> isopt(testpnts[jj][0],testpnts[jj][1],testpnts[jj][2]) ;
	vector3d<double> pt = xcd + double(radius)*isopt ;
	double x[3] = {pt.x,pt.y,pt.z} ;
	double ya = frsum(x,beta,XC,J,exponent) ;
	double r = (fexact[jj]-ya) ;
	mxr = max(mxr,fabs(r)) ;
	for(int i=0;i<9;++i) {
	  B[i] += J[i]*r ;
	  for(int j=0;j<9;++j)
	    ATA[i][j] += J[i]*J[j] ;
	}
      }
      for(int i=0;i<9;++i) 
	ATA[i][i] += 1e-20 ;

      pivot_type piv[9] ;
      ATA.decompose_lu_pivot(piv) ;
      double dBeta[9] ;
      //    ATA.solve_lu_pivot(B,dBeta,piv) ;
      solve_lu_pivot_approx(&ATA[0][0],B,dBeta,piv,9) ;
      double maxdbeta = fabs(dBeta[0]) ;
      for(int i=1;i<9;++i)
	maxdbeta = max(maxdbeta,fabs(dBeta[i])) ;
      double alpha = 1.0 ;
    
      if(maxdbeta > 0.1*radius)
	alpha = 0.1*radius/maxdbeta ;
      // Now do line search to find alpha that reduces residual
      int nls = 10 ;
      int ls ;
      for(ls=0;ls<nls;++ls) {
	double btmp[9] ;
	for(int i=0;i<9;++i)
	  btmp[i] = beta[i] + alpha * dBeta[i] ;

	double mxrstep = 0 ;
	for(int jj=0;jj<nfit;++jj) {
	  vector3d<double> isopt(testpnts[jj][0],testpnts[jj][1],testpnts[jj][2]) ;
	  vector3d<double> pt = xcd + double(radius)*isopt ;
	  double x[3] = {pt.x,pt.y,pt.z} ;
	  double ya = frsum_basic(x,btmp,XC,exponent) ;
	  double r = (fexact[jj]-ya) ;
	  mxrstep = max(mxrstep,fabs(r)) ;
	}
	if(mxrstep > mxr)
	  alpha = alpha/double(2+ls) ;
	else {
	  break ;
	}
      }
      if(ls == nls) {
	// line search failed.. probably at minumum so exit
	//	cout << "line search failed, alpha = " << alpha << endl ;
	break ;
      }

      for(int i=0;i<9;++i)
	beta[i] += alpha * dBeta[i] ;
      //      cout << "iter=" << iter << ", maxdbeta = " << maxdbeta
      //	   << ",mxr="<<mxr << ",alpha=" << alpha << ",radius="<< radius << endl ;
      if(mxr < 0.01*fmin || maxdbeta < radius*1e-6)
	break ;
    }
    
    qp1 = vector3d<float>(beta[0],beta[1],beta[2]) ;
    qp2 = vector3d<float>(beta[3],beta[4],beta[5]) ;
    qp3 = vector3d<float>(beta[6],beta[7],beta[8]) ;
  
    // compute error
    double maxerr = 0 ;
    for(int jj=0;jj<ntest;++jj) {
      vector3d<double> isopt(testpnts[jj][0],testpnts[jj][1],testpnts[jj][2]) ;
      vector3d<double> pt = xcd + double(radius)*isopt ;
      double x[3] = {pt.x,pt.y,pt.z} ;
      double ya = frsum(x,beta,XC,J,exponent) ;
      double r = fexact[jj]-ya ;
      maxerr = max(maxerr,fabs(r)/fexact[jj]) ;
    }
    //    cout << "maxerr = " << maxerr <<", size=" << end-start << endl ;
    error = maxerr ;
  }

  int buildTreeRecurse(vector<treeInfoBase> &nodeDataTree,
		       vector<nodeDataBase> &nodeDataSort,
		       int node) {
    treeInfoBase & nodeInfo = nodeDataTree[node] ;
    const int node_min_size = 42 ;
    const int start = nodeInfo.start ;
    const int end = nodeInfo.end ;
    const int node_size = end - start ;
    if(node_size < node_min_size) { // base case
      nodeInfo.left = -1 ;
      nodeInfo.right = -1 ;
      nodeInfo.coord = -1 ;
      nodeInfo.split = 0 ;
      return node ;
    }
  
    // Not base case, so split array.  First find bounding box.
    double xmin = nodeDataSort[start].pos.x ;
    double xmax = xmin ;
    double ymin = nodeDataSort[start].pos.y ;
    double ymax = ymin ;
    double zmin = nodeDataSort[start].pos.z ;
    double zmax = zmin ;

    for(int i=start+1;i<end;++i) {
      vector3d<double> p = nodeDataSort[i].pos ;
      xmin = min(xmin,p.x) ;
      xmax = max(xmax,p.x) ;
      ymin = min(ymin,p.y) ;
      ymax = max(ymax,p.y) ;
      zmin = min(zmin,p.z) ;
      zmax = max(zmax,p.z) ;
    }
    int split = -1 ;
    //  cout << "mxmn="<< xmax << ' ' << xmin << ' ' << ymax << ' ' << ymin
    //       << ' ' << zmax << ' ' << zmin << endl ;


    if((xmax-xmin > ymax-ymin) && (xmax-xmin > zmax-zmin)) {
      // split in x coordinate

      //Set pval to the center of the bounding box x-coordinates
      double pval = 0.5*(xmax+xmin) ;
      //    cout << "xpval = " << pval << ' ' << end-start << endl ;

      nodeInfo.coord = 0 ;
      nodeInfo.split = pval ;
      int p1 = start ;
      int p2 = end-1 ;
      while(p1 <= p2) {
	if(nodeDataSort[p1].pos.x <= pval) { // if less than pivot
	  p1++ ; // leave in first list
	} else { // otherwise move to end of the list
	  while(nodeDataSort[p2].pos.x > pval && p1 < p2)
	    p2-- ;
	  std::swap(nodeDataSort[p1],nodeDataSort[p2]) ;
	  p2-- ; // move end of list pointer
	}
      }
      split = p1 ;
    } else if(ymax-ymin > zmax-zmin) {
      // split in y coordinate

      double pval = 0.5*(ymax+ymin) ;
      //    cout << "ypval = " << pval << ' ' << end-start << endl ;

      nodeInfo.coord = 1 ;
      nodeInfo.split = pval ;
      int p1 = start ;
      int p2 = end-1 ;
      while(p1 <= p2) {
	if(nodeDataSort[p1].pos.y <= pval) { // if less than pivot
	  p1++ ; // leave in first list
	} else { // otherwise move to end of the list
	  while(nodeDataSort[p2].pos.y > pval && p1 < p2)
	    p2-- ;
	  std::swap(nodeDataSort[p1],nodeDataSort[p2]) ;
	  p2-- ; // move end of list pointer
	}
      }
      split = p1 ;
    } else {
      // split in z coordinate

      double pval = 0.5*(zmax+zmin) ;
      //    cout << "zpval = " << pval << ' ' << end-start <<endl ;

      nodeInfo.coord = 2 ;
      nodeInfo.split = pval ;
      int p1 = start ;
      int p2 = end-1 ;
      while(p1 <= p2) {
	if(nodeDataSort[p1].pos.z <= pval) { // if less than pivot
	  p1++ ; // leave in first list
	} else { // otherwise move to end of the list
	  while(nodeDataSort[p2].pos.z > pval && p1 < p2)
	    p2-- ;
	  std::swap(nodeDataSort[p1],nodeDataSort[p2]) ;
	  p2-- ; // move end of list pointer
	}
      }
      split = p1 ;
      if(split == end || split < start) {
	std::cerr << "split failed to split z, zmax=" << zmax<< ",zmin="
	     << zmin << ",pval=" << pval << endl ;
      }
      for(int i=start;i<split;++i)
	if(nodeDataSort[i].pos.z > pval) {
	  std::cerr << "z split fail<, i=" << i 
	       << ' ' <<  nodeDataSort[i].pos.z << ' ' << pval << endl ;
	} 
      for(int i=split;i<end;++i)
	if(nodeDataSort[i].pos.z <= pval) {
	  std::cerr << "zsplit fail>, i=" << i 
	       << ' ' <<  nodeDataSort[i].pos.z << ' ' << pval << endl ;
	}


    }
    //  Now build tree nodes, first left, then right
    int left = nodeDataTree.size() ;
    nodeDataTree[node].left = left ;
    nodeDataTree.push_back(treeInfoBase()) ;
    nodeDataTree[left].start = start ;
    nodeDataTree[left].end = split ;
    buildTreeRecurse(nodeDataTree,nodeDataSort,left) ;

    int right = nodeDataTree.size() ;
    nodeDataTree[node].right = right ;
    nodeDataTree.push_back(treeInfoBase()) ;
    nodeDataTree[right].start = split ;
    nodeDataTree[right].end = end ;
    buildTreeRecurse(nodeDataTree,nodeDataSort,right) ;

    return node ;
  }

  // build the tree assigning branches to processors
  int buildTreeRecurseProc(vector<treeInfoBase> &nodeDataTree,
			   vector<nodeDataBase> &nodeDataSort,
			   vector<pair<int,int> > &proclims,
			   int node) {
    treeInfoBase & nodeInfo = nodeDataTree[node] ;
    const int start = nodeInfo.start ;
    const int end = nodeInfo.end ;
    const int psize = proclims[node].second-proclims[node].first ;
    if(psize==1) { // base case
      nodeInfo.left = -1 ;
      nodeInfo.right = -1 ;
      nodeInfo.coord = -1 ;
      nodeInfo.split = 0 ;
      //      Loci::debugout << "proc=" << proclims[node].first
      //                     << ", alloc=" << start << "," << end << ":"
      //                     << end-start << endl ;
      return node ;
    }
  
    // Not base case, so split array.  First find bounding box.
    double xmin = nodeDataSort[start].pos.x ;
    double xmax = xmin ;
    double ymin = nodeDataSort[start].pos.y ;
    double ymax = ymin ;
    double zmin = nodeDataSort[start].pos.z ;
    double zmax = zmin ;

    for(int i=start+1;i<end;++i) {
      vector3d<double> p = nodeDataSort[i].pos ;
      xmin = min(xmin,p.x) ;
      xmax = max(xmax,p.x) ;
      ymin = min(ymin,p.y) ;
      ymax = max(ymax,p.y) ;
      zmin = min(zmin,p.z) ;
      zmax = max(zmax,p.z) ;
    }
    int split = -1 ;

    //    Loci::debugout << "pmin=" << proclims[node].first 
    //                   << ",pmax=" << proclims[node].second
    //                   << ",xmin=" << xmin <<",xmax="<< xmax
    //                   << ",ymin=" << ymin <<",ymax="<< ymax
    //                   << ",zmin=" << zmin <<",zmax="<< zmax
    //                   << endl ;
    //  cout << "mxmn="<< xmax << ' ' << xmin << ' ' << ymax << ' ' << ymin
    //       << ' ' << zmax << ' ' << zmin << endl ;

    double fac = double(psize/2)/double(psize) ;

    if((xmax-xmin > ymax-ymin) && (xmax-xmin > zmax-zmin)) {
      // split in x coordinate                  
        
      //Set pval to the center of the bounding box x-coordinates
      float pval = xmin + fac*(xmax-xmin) ;
      //    cout << "xpval = " << pval << ' ' << end-start << endl ;

      nodeInfo.coord = 0 ;
      nodeInfo.split = pval ;
      int p1 = start ;
      int p2 = end-1 ;
      while(p1 <= p2) {
	if(nodeDataSort[p1].pos.x <= pval) { // if less than pivot
	  p1++ ; // leave in first list
	} else { // otherwise move to end of the list
	  while(nodeDataSort[p2].pos.x > pval && p1 < p2)
	    p2-- ;
	  std::swap(nodeDataSort[p1],nodeDataSort[p2]) ;
	  p2-- ; // move end of list pointer
	}
      }
      split = p1 ;
    } else if(ymax-ymin > zmax-zmin) {
      // split in y coordinate

      float pval = ymin + fac*(ymax-ymin) ;
      //    cout << "ypval = " << pval << ' ' << end-start << endl ;

      nodeInfo.coord = 1 ;
      nodeInfo.split = pval ;
      int p1 = start ;
      int p2 = end-1 ;
      while(p1 <= p2) {
	if(nodeDataSort[p1].pos.y <= pval) { // if less than pivot
	  p1++ ; // leave in first list
	} else { // otherwise move to end of the list
	  while(nodeDataSort[p2].pos.y > pval && p1 < p2)
	    p2-- ;
	  std::swap(nodeDataSort[p1],nodeDataSort[p2]) ;
	  p2-- ; // move end of list pointer
	}
      }
      split = p1 ;
    } else {
      // split in z coordinate
      float pval = zmin + fac*(zmax-zmin) ;
      //    cout << "zpval = " << pval << ' ' << end-start <<endl ;

      nodeInfo.coord = 2 ;
      nodeInfo.split = pval ;
      int p1 = start ;
      int p2 = end-1 ;
      while(p1 <= p2) {
	if(nodeDataSort[p1].pos.z <= pval) { // if less than pivot
	  p1++ ; // leave in first list
	} else { // otherwise move to end of the list
	  while(nodeDataSort[p2].pos.z > pval && p1 < p2)
	    p2-- ;
	  std::swap(nodeDataSort[p1],nodeDataSort[p2]) ;
	  p2-- ; // move end of list pointer
	}
      }
      split = p1 ;
      if(split == end || split < start) {
	std::cerr << "split failed to split z, zmax=" << zmax<< ",zmin="
                       << zmin << ",pval=" << pval << endl ;
        std::cerr << "fac=" << fac << "," << fac*(zmax+zmin) << endl ;
      }
      for(int i=start;i<split;++i)
	if(nodeDataSort[i].pos.z > pval) {
	  std::cerr << "z split fail<, i=" << i 
	       << ' ' <<  nodeDataSort[i].pos.z << ' ' << pval << endl ;
	} 
      for(int i=split;i<end;++i)
	if(nodeDataSort[i].pos.z <= pval) {
	  std::cerr << "zsplit fail>, i=" << i 
	       << ' ' <<  nodeDataSort[i].pos.z << ' ' << pval << endl ;
	}


    }

    //  Now build tree nodes, first left, then right
    int left = nodeDataTree.size() ;
    nodeDataTree[node].left = left ;
    nodeDataTree.push_back(treeInfoBase()) ;
    proclims.push_back(proclims[node]) ;
    // Split partition among processors according to size
    int psplit = int(double(psize)*double(split-start)/double(end-start)+.5)+proclims[node].first ;
    // make sure a split always assigns a processor
    psplit = max(proclims[node].first+1,psplit) ;
    psplit = min(proclims[node].second-1,psplit) ; 
    //    Loci::debugout << "psplit=" << psplit<<"{"<< proclims[node].first << ","
    //                   << proclims[node].second << "} "
    //                   << "left=" << split-start <<",right=" <<end-split << endl ;
    proclims[left].second = psplit ;
    nodeDataTree[left].start = start ;
    nodeDataTree[left].end = split ;
    buildTreeRecurseProc(nodeDataTree,nodeDataSort,proclims,left) ;

    int right = nodeDataTree.size() ;
    nodeDataTree[node].right = right ;
    nodeDataTree.push_back(treeInfoBase()) ;
    proclims.push_back(proclims[node]) ;
    proclims[right].first = psplit ;
    nodeDataTree[right].start = split ;
    nodeDataTree[right].end = end ;
    buildTreeRecurseProc(nodeDataTree,nodeDataSort,proclims,right) ;

    return node ;
  }


  void computeProcWeightsCentroids(const vector<treeInfoBase> &nodeDataTree,
				   const vector<nodeDataBase> &nodeDataSort,
				   const vector<pair<int,int> > &proclims,
				   vector<double> &weights,
				   vector<vector3d<float> > &centroid,
				   int node) {
    const int left = nodeDataTree[node].left ;
    const int right = nodeDataTree[node].right ;
    if(node > int(proclims.size()) ||
       proclims[node].second-proclims[node].first == 1) {
      return ;
    }
    if(left >= 0) { // follow branches
      computeProcWeightsCentroids(nodeDataTree,nodeDataSort,proclims,
				 weights,centroid,left) ;
      computeProcWeightsCentroids(nodeDataTree,nodeDataSort,proclims,
				 weights,centroid,right) ;
      // Now combine leaves are complete do higher levels of the tree
      weights[node] = weights[left]+weights[right] ;
      centroid[node] = (1./weights[node])*(weights[left]*centroid[left]+
					   weights[right]*centroid[right]) ;
    }
  }

    
  void computeQuadPointsTree(const vector<treeInfoBase> &nodeDataTree,
			     const vector<nodeDataBase> &nodeDataSort,
			     vector<double> &weights,
			     vector<vector3d<float> > &centroid,
			     vector<float> &radius,
			     vector<float> &err,
			     vector<Array<vector3d<float>,3> > &q,
			     int node) {
    const double factor = 3 ;
    const int exponent = 3 ;
    const int ntest = 32 ;
    const int nfit = 20 ;
    const int left = nodeDataTree[node].left ;
    const int right = nodeDataTree[node].right ;
    const int start = nodeDataTree[node].start ;
    const int end = nodeDataTree[node].end ;
    if(left >= 0) { // follow branches

      computeQuadPointsTree(nodeDataTree,nodeDataSort,
			    weights, centroid, radius, err, q,
			    left) ;
      computeQuadPointsTree(nodeDataTree, nodeDataSort,
			    weights, centroid, radius, err, q,
			    right) ;
      // Now combine leaves are complete do higher levels of the tree
      weights[node] = weights[left]+weights[right] ;
      centroid[node] = (1./weights[node])*(weights[left]*centroid[left]+
					   weights[right]*centroid[right]) ;

      float r2 = max(1e-15f,dot(centroid[node],centroid[node])*1e-10f) ;
      for(int i=start;i<end;++i) {
	vector3d<float> dv = nodeDataSort[i].pos-centroid[node] ;
	r2 = max(r2,dot(dv,dv)) ;
      }
      radius[node] = sqrt(r2) ;
      double ftest[ntest] ;
      double rtest = radius[node]*factor ;
      get_fexact(nodeDataSort,start,end,exponent,rtest,centroid[node],
		 testpnts,&ftest[0],ntest) ;
      double rw = 1./weights[node] ;
      for(int i=0;i<ntest;++i)
	ftest[i] *= rw ;
      makeQuad(nodeDataSort,start,end,testpnts,&ftest[0],nfit,ntest,exponent,
	       centroid[node],rtest,
	       q[node][0],q[node][1],q[node][2],err[node]) ;
    } else { // leaf level
      double weight = 0 ;
      vector3d<double> xsum(0,0,0) ;
      for(int i=start;i<end;++i) {
	double w = nodeDataSort[i].weight ;
	weight += w ;
	xsum += w*vector3d<double>(nodeDataSort[i].pos) ;
      }
      weights[node] = weight ;
      xsum *= (1./weight) ;
      centroid[node] = xsum ;
      double r2 = max(1e-15,dot(xsum,xsum)*1e-10) ;
      for(int i=start;i<end;++i) {
	vector3d<float> dv = nodeDataSort[i].pos-vector3d<float>(xsum) ;
	r2 = max(r2,double(dot(dv,dv))) ;
      }
      radius[node] = sqrt(r2) ;
      double ftest[ntest] ;
      double rtest = radius[node]*factor ;
      get_fexact(nodeDataSort,start,end,exponent,rtest,xsum,
		 testpnts,&ftest[0],ntest) ;
      double rw = 1./weights[node] ;
      for(int i=0;i<ntest;++i)
	ftest[i] *= rw ;
      makeQuad(nodeDataSort,start,end,testpnts,&ftest[0],nfit,ntest,exponent,
	       centroid[node],rtest,
	       q[node][0],q[node][1],q[node][2],err[node]) ;
    }
  }


  void deformApproxTree::
  computeDeformationQuadsTree(const vector<double> &weights, int node) {
    const int left = nodeDataTree[node].left ;
    const int right = nodeDataTree[node].right ;
    if(left >= 0) { // follow branches
      computeDeformationQuadsTree(weights,left) ;
      computeDeformationQuadsTree(weights,right) ;
      // right now we compute the higher levels of the tree exactly
      // but it might be more efficient to compute high levels of
      // the tree using approximations, but this is not currently
      // implemented.
    }
    // Compute deformation variables for quad nodes at this level
    const int ntst = 4 ;
    vector3d<float> xc = centroid[node] ;
    vector3d<float> rtestpts[ntst] ;
    rtestpts[0] = qpts[node][0]-xc ;
    rtestpts[1] = qpts[node][1]-xc ;
    rtestpts[2] = qpts[node][2]-xc ;
    rtestpts[3] = 3.f*xc-qpts[node][0]-qpts[node][1]-qpts[node][2]  ;

    for(int i=0;i<ntst;++i) {
      rtestpts[i] *= 1./max(norm(rtestpts[i]),1e-10f) ;
      if(dot(rtestpts[i],rtestpts[i]) < 0.99)
        rtestpts[i] = vector3d<float>(1.,0,0) ;
    }
    
    // extract displacement and rotation matrix at centroid and 
    // compute centroid values for deformation
    vector<vector<vector3d<double> > > defs(4) ;
    for(int i=0;i<4;++i) {
      vector<vector3d<double> > tmp(ntst) ;
      defs[i].swap(tmp) ;
      for(int k=0;k<ntst;++k)
	defs[i][k] = vector3d<double>(0,0,0) ;
    }

    vector3d<double> meanvals[4] ;
    for(int i=0;i<4;++i) 
      meanvals[i] = vector3d<double>(0,0,0) ;

    const int exponent = 3 ;
    const float factor = 3.0 ;
    double wexactsum[ntst] ;
    for(int i=0;i<ntst;++i)
      wexactsum[i] = 0 ;

    if(left >=0) {
      const int branch[2] = {left,right} ;
      for(int i=0;i<2;++i) {
	const int side = branch[i] ;
	vector3d<float> dside = xc-centroid[side] ;
	tensor3d<float> rside = qrot[side][0];
	rside += qrot[side][1] ;
	rside += qrot[side][2] ;
	rside += qrot[side][3] ;
	rside *= 0.25 ;
	vector3d<float> disp = 0.25*(qdisp[side][0]+qdisp[side][1]+
				     qdisp[side][2]+qdisp[side][3]) +
	  dot(rside,dside) ;
	meanvals[0] += vector3d<double>(weights[side]*disp) ;
	meanvals[1] += vector3d<double>(weights[side]*rside.x) ;
	meanvals[2] += vector3d<double>(weights[side]*rside.y) ;
	meanvals[3] += vector3d<double>(weights[side]*rside.z) ;
      }
      for(int i=0;i<4;++i) {
	meanvals[i] *= 1./weights[node] ;
      }
    } else {
      const int start = nodeDataTree[node].start ;
      const int end = nodeDataTree[node].end ;
      for(int i=start;i<end;++i) {
	tens3d rotmat = nodeRotors[i].matrix() ;
	rotmat.x.x -= 1.0 ; // Convert rotation to delta change
	rotmat.y.y -= 1.0 ;
	rotmat.z.z -= 1.0 ;
	vector3d<double> dp(xc-nodeDataSort[i].pos) ;
	vect3d displT = vector3d<double>(nodeDisp[i])+dot(rotmat,dp) ;
	double w = nodeDataSort[i].weight ;
	meanvals[0] += w*displT ;
	meanvals[1] += w*rotmat.x ;
	meanvals[2] += w*rotmat.y ;
	meanvals[3] += w*rotmat.z ;
      }
      for(int i=0;i<4;++i) {
	meanvals[i] *= 1./weights[node] ;
      }
    }

    // Compute exact mapping to test points
    const int start = nodeDataTree[node].start ;
    const int end = nodeDataTree[node].end ;
    for(int i=start;i<end;++i) {
      tens3d rotmat = nodeRotors[i].matrix() ;
      rotmat.x.x -= 1.0 ; // Convert rotation to delta change
      rotmat.y.y -= 1.0 ;
      rotmat.z.z -= 1.0 ;
      vector3d<double> dp(xc-nodeDataSort[i].pos) ;
      vect3d displT = vector3d<double>(nodeDisp[i])+dot(rotmat,dp) ;
      
      for(int j=0;j<ntst;++j) {
	vector3d<float> rtest = rtestpts[j]*radius[node]*factor+xc ;
	const double r = norm(nodeDataSort[i].pos-rtest) ;
	const double wn = nodeDataSort[i].weight ;
	const double wr = wn*pow(r,-exponent) ;
	wexactsum[j] += wr ;
	defs[0][j] += wr*displT ;
	defs[1][j] += wr*rotmat.x ;
	defs[2][j] += wr*rotmat.y ;
	defs[3][j] += wr*rotmat.z ;
      }
    }

      
    for(int i=0;i<4;++i) {
      for(int j=0;j<ntst;++j) {
	defs[i][j] *= 1./wexactsum[j] ;
      }
    }

    const vector3d<double> zero(0.0,0.0,0.0) ;
    float droterr = 0 ; 
    int nrsample = 4 ;
    double W[nrsample][4] ;
    double A00=0,A11=0,A22=0,A01=0,A02=0,A12= 0 ;
    for(int s=0;s<nrsample;++s) {
      // test point
      vector3d<float> pt = rtestpts[s]*radius[node]*factor+xc ;

      // Weights associated with quadripole point for this test point
      for(int i=0;i<4;++i)
	W[s][i] = 0 ;
      double Wsum = 0 ;
      for(int i=0;i<3;++i) {
	const double r = norm(pt-qpts[node][i]) ;
	W[s][i] = pow(r,-exponent) ;
	Wsum += W[s][i] ;
      }
      W[s][3] = pow(norm(pt-4.*xc+qpts[node][0]+qpts[node][1]+qpts[node][2]),
		    -exponent) ;
      Wsum += W[s][3] ;
	  
      double rWsum = 1./Wsum ;
      for(int i=0;i<4;++i)
	W[s][i] *= rWsum ;
      double J[3] = {W[s][0]-W[s][3],W[s][1]-W[s][3],W[s][2]-W[s][3]} ;
	  
      A00 += J[0]*J[0] ;
      A11 += J[1]*J[1] ;
      A22 += J[2]*J[2] ;
      A01 += J[0]*J[1] ;
      A02 += J[0]*J[2] ;
      A12 += J[1]*J[2] ;
    }
    double det = (+A00*A11*A22+A01*A12*A02+A02*A01*A12
		  -A02*A11*A02-A00*A12*A12-A01*A01*A22) ;
    double ATAinv[3][3] ;
    double rdet = 1./(det+((det>0.0)?1e-12:-1e-12)) ;
    ATAinv[0][0] = rdet*(A11*A22 - A12*A12) ;
    ATAinv[0][1] = rdet*(A02*A12 - A01*A22) ;
    ATAinv[0][2] = rdet*(A01*A12 - A02*A11) ;
	  
    ATAinv[1][0] = ATAinv[0][1] ;
    ATAinv[1][1] = rdet*(A00*A22 - A02*A02) ;
    ATAinv[1][2] = rdet*(A02*A01 - A00*A12) ;

    ATAinv[2][0] = ATAinv[0][2] ;
    ATAinv[2][1] = ATAinv[1][2] ;
    ATAinv[2][2] = rdet*(A00*A11 - A01*A01) ;
    
    for(int field=0;field<4;++field) { 
      // Now compute quad point values using least squares
      // Symmetric matrix terms (3x3)
      vector3d<double> Y[3] = {zero,zero,zero} ;
      for(int s=0;s<nrsample;++s) {
	double J[3] = {W[s][0]-W[s][3],W[s][1]-W[s][3],W[s][2]-W[s][3]} ;
	  
	for(int i=0;i<3;++i) {
	  Y[i] += J[i]*(defs[field][s]-4.*W[s][3]*meanvals[field]) ;
	}
      }
      vector3d<double> v1= meanvals[field] ;
      vector3d<double> v2 = v1 ;
      vector3d<double> v3 = v1 ;
      vector3d<double> v4 = v1 ;

      if(fabs(det) > 1e-9) {
	// Compute values at quad points
	v1 = ATAinv[0][0]*Y[0]+ATAinv[0][1]*Y[1]+ATAinv[0][2]*Y[2] ;
	v2 = ATAinv[1][0]*Y[0]+ATAinv[1][1]*Y[1]+ATAinv[1][2]*Y[2] ;
	v3 = ATAinv[2][0]*Y[0]+ATAinv[2][1]*Y[1]+ATAinv[2][2]*Y[2] ;
	v4 = 4.*meanvals[field] -v1-v2-v3;
      }
      double err[4] = {0,0,0,0} ;
      for(int s=0;s<ntst;++s) {
	// test point
	vector3d<float> pt = xc + rtestpts[s]*radius[node]*factor ;
	// Weights associated with quadripole point for this test point
	double W[4] = {0,0,0,0} ;
	double Wsum = 0 ;
	for(int i=0;i<3;++i) {
	  const double r = norm(pt-qpts[node][i]) ;
	  W[i] = pow(r,-exponent) ;
	  Wsum += W[i] ;
	}
	const double r3 = norm(pt-4.*xc+qpts[node][0]+qpts[node][1]+qpts[node][2]) ;
	W[3] = pow(r3,-exponent) ;
	Wsum += W[3] ;
	double rWsum = 1./Wsum ;
	for(int i=0;i<4;++i)
	  W[i] *= rWsum ;
	err[field] = max(err[field],
			 norm(defs[field][s]-
			      (W[0]*v1+W[1]*v2+W[2]*v3+W[3]*v4))) ;
      }
      switch(field) {
      case 0:
	qdisp[node][0] = v1 ;
	qdisp[node][1] = v2 ;
	qdisp[node][2] = v3 ;
	qdisp[node][3] = v4 ;
	droterr = max(droterr,float(err[field]/(max(norm(meanvals[field]),1e-4)))) ;
	break ;
      case 1:
	qrot[node][0].x = v1 ;
	qrot[node][1].x = v2 ;
	qrot[node][2].x = v3 ;
	qrot[node][3].x = v4 ;
	droterr = max(droterr,float(err[field])) ;
	break ;
      case 2:
	qrot[node][0].y = v1 ;
	qrot[node][1].y = v2 ;
	qrot[node][2].y = v3 ;
	qrot[node][3].y = v4 ;
	droterr = max(droterr,float(err[field])) ;
	break ;
      case 3:
	qrot[node][0].z = v1 ;
	qrot[node][1].z = v2 ;
	qrot[node][2].z = v3 ;
	qrot[node][3].z = v4 ;
	droterr = max(droterr,float(err[field])) ;
	break ;
      }
    }
    drot[node]=droterr ;
  }
  
  
  void buildTreeSerialPart(vector<nodeDataBase> &nodeDataSort,
			   vector<treeInfoBase> &nodeDataTree,
			   int &pbuild,
			   vector<pair<int,int> > &proclims,
			   vector<int> &proctreebase,
			   vector<int> &treeptn,
			   MPI_Comm comm) {
    // First gather
    int nsz = nodeDataSort.size() ;
    for(int i=0;i<nsz;++i)  {
      nodeDataSort[i].order = i ;
    }
    int p = 1 ;
    int r = 0 ;
    MPI_Comm_size(comm,&p) ;
    MPI_Comm_rank(comm,&r) ;

    pbuild = max(min(p,nsz/4096),1) ;

    // build tree (serially on each processor)
    nodeDataTree = vector<treeInfoBase>() ;
    proclims = vector<pair<int,int> >() ;
    nodeDataTree.push_back(treeInfoBase()) ;
    nodeDataTree[0].start = 0 ;
    nodeDataTree[0].end = nsz ;
    proclims.push_back(pair<int,int>(0,pbuild)) ;
    buildTreeRecurseProc(nodeDataTree,nodeDataSort,proclims,0) ;
    proctreebase = vector<int>(pbuild,-1) ;
    for(size_t i=0;i<proclims.size();++i)
      if(proclims[i].second-proclims[i].first == 1)
	proctreebase[proclims[i].first] = i ;
    
    treeptn = vector<int>(pbuild+1) ;
    treeptn[0] = nodeDataTree.size() ;
    for(int i=0;i<pbuild;++i) {
      int tp = proctreebase[i] ;
      buildTreeRecurse(nodeDataTree,nodeDataSort,tp) ;
      treeptn[i+1] = nodeDataTree.size() ;
    }
  }

  void buildWeightsApprox(const vector<nodeDataBase> &nodeDataSort,
			  const vector<treeInfoBase> &nodeDataTree,
			  const vector<pair<int,int> > &proclims,
			  const vector<int> &proctreebase,
			  const vector<int> &treeptn,
			  vector<double> &weights,
			  vector<vector3d<float> > &centroid,
			  vector<float> &radius,
			  vector<Array<vector3d<float>,3> > &qpts,
			  vector<float> &qerr,
			  MPI_Comm comm) {
    int p = 1 ;
    int r = 0 ;
    MPI_Comm_size(comm,&p) ;
    MPI_Comm_rank(comm,&r) ;
    int tsz = nodeDataTree.size() ;
    weights = vector<double>(tsz) ;
    centroid = vector<vector3d<float> >(tsz) ;
    radius = vector<float>(tsz) ;
    qpts = vector<Array<vector3d<float>,3> >(tsz) ;
    qerr = vector<float>(tsz) ;

    // Compute quad points assigned to this processor
    computeQuadPointsTree(nodeDataTree, nodeDataSort,
			  weights, centroid, radius, qerr, qpts,
			  proctreebase[r]) ;

    // gather local processor results
    // Collect Centroids
    if(p > 1) {
      vector<vector3d<float> > centbuf(treeptn[r+1]-treeptn[r]) ;
      for(int i=treeptn[r];i<treeptn[r+1];++i)
	centbuf[i-treeptn[r]] = centroid[i] ;
      vector<int> recv_counts(p) ;
      for(int i=0;i<p;++i)
	recv_counts[i] = (treeptn[i+1]-treeptn[i])*3 ;
      vector<int> displs(p) ;
      displs[0] = 0 ;
      for(int i=1;i<p;++i)
	displs[i] = displs[i-1]+recv_counts[i-1] ;
      MPI_Allgatherv(&centbuf[0],centbuf.size()*3,MPI_FLOAT,
		     &centroid[treeptn[0]],
		     &recv_counts[0],&displs[0], MPI_FLOAT,
		     comm) ;
      
      centbuf.resize(p) ;
      vector3d<float> cent = centroid[proctreebase[r]] ;
      MPI_Allgather(&cent,3,MPI_FLOAT,&centbuf[0],3,MPI_FLOAT,comm) ;
      for(int i=0;i<p;++i)
	centroid[proctreebase[i]] = centbuf[i] ;
    }
    if(p > 1) {
      // Collect radius
      vector<float > buf(treeptn[r+1]-treeptn[r]) ;
      for(int i=treeptn[r];i<treeptn[r+1];++i)
	buf[i-treeptn[r]] = radius[i] ;
      vector<int> recv_counts(p) ;
      for(int i=0;i<p;++i)
	recv_counts[i] = (treeptn[i+1]-treeptn[i]) ;
      vector<int> displs(p) ;
      displs[0] = 0 ;
      for(int i=1;i<p;++i)
	displs[i] = displs[i-1]+recv_counts[i-1] ;
      MPI_Allgatherv(&buf[0],buf.size(),MPI_FLOAT,
		     &radius[treeptn[0]],
		     &recv_counts[0],&displs[0], MPI_FLOAT,
		     comm) ;
	// collect qerr
      for(int i=treeptn[r];i<treeptn[r+1];++i)
	buf[i-treeptn[r]] = qerr[i] ;
      MPI_Allgatherv(&buf[0],buf.size(),MPI_FLOAT,
		     &qerr[treeptn[0]],
		     &recv_counts[0],&displs[0], MPI_FLOAT,
		     comm) ;
      buf.resize(p) ;
      float rad = radius[proctreebase[r]] ;
      MPI_Allgather(&rad,1,MPI_FLOAT,&buf[0],1,MPI_FLOAT,comm) ;
      for(int i=0;i<p;++i)
	radius[proctreebase[i]] = buf[i] ;
    }
    if(p > 1) {
      // Collect weights
      vector<double > buf(treeptn[r+1]-treeptn[r]) ;
      for(int i=treeptn[r];i<treeptn[r+1];++i)
	buf[i-treeptn[r]] = weights[i] ;
      vector<int> recv_counts(p) ;
      for(int i=0;i<p;++i)
	recv_counts[i] = (treeptn[i+1]-treeptn[i]) ;
      vector<int> displs(p) ;
      displs[0] = 0 ;
      for(int i=1;i<p;++i)
	displs[i] = displs[i-1]+recv_counts[i-1] ;
      MPI_Allgatherv(&buf[0],buf.size(),MPI_DOUBLE,
		     &weights[treeptn[0]],
		     &recv_counts[0],&displs[0], MPI_DOUBLE,
		     comm) ;
      buf.resize(p) ;
      double w = weights[proctreebase[r]] ;
      
      MPI_Allgather(&w,1,MPI_DOUBLE,&buf[0],1,MPI_DOUBLE,comm) ;
      for(int i=0;i<p;++i)
	weights[proctreebase[i]] = buf[i] ;
    }
    // Collect quad points
    if(p > 1) {
      vector<Array<vector3d<float>,3> > buf(treeptn[r+1]-treeptn[r]) ;
      for(int i=treeptn[r];i<treeptn[r+1];++i)
	buf[i-treeptn[r]] = qpts[i] ;
      vector<int> recv_counts(p) ;
      for(int i=0;i<p;++i)
	recv_counts[i] = (treeptn[i+1]-treeptn[i])*9 ;
      vector<int> displs(p) ;
      displs[0] = 0 ;
      for(int i=1;i<p;++i)
	displs[i] = displs[i-1]+recv_counts[i-1] ;
      MPI_Allgatherv(&buf[0],buf.size()*9,MPI_FLOAT,
		     &qpts[treeptn[0]],
		     &recv_counts[0],&displs[0], MPI_FLOAT,
		     comm) ;
      buf.resize(p) ;
      Array<vector3d<float>,3> qp = qpts[proctreebase[r]] ;
      vector<Array<vector3d<float>,3> > buf2(p) ;
      MPI_Allgather(&qp,9,MPI_FLOAT,&buf2[0],9,MPI_FLOAT,comm) ;
      for(int i=0;i<p;++i)
	qpts[proctreebase[i]] = buf2[i] ;
    }
      
    // Now compute top part of tree distributed among processors
    computeProcWeightsCentroids(nodeDataTree,nodeDataSort,proclims,
				weights,centroid,0) ;
    int start = nodeDataTree[proctreebase[r]].start ;
    int end = nodeDataTree[proctreebase[r]].end ;
    vector<int> treework ;
    for(size_t i=0;i<proclims.size();++i) 
      if(proclims[i].second-proclims[i].first > 1)
	treework.push_back(i) ;
    
    // First compute radius of top tree nodes
    vector<float> radii2(treework.size(),0) ;
    for(size_t i=0;i<treework.size();++i) {
      int node = treework[i] ;
      int sloc = max(start,nodeDataTree[node].start) ;
      int eloc = min(end,nodeDataTree[node].end) ;
      float r2 = max(1e-15f,dot(centroid[node],centroid[node])*1e-10f) ;
      for(int j=sloc;j<eloc;++j) {
	vector3d<float> dv = (nodeDataSort[j].pos-centroid[node]) ;
	r2 = max(r2,dot(dv,dv)) ;
      }
      radii2[i] = r2 ;
    }
    vector<float> radii(treework.size(),0) ;
      
    MPI_Allreduce(&radii2[0],&radii[0],treework.size(),
		  MPI_FLOAT,MPI_MAX,comm) ;
    for(size_t i=0;i<treework.size();++i) {
      int node = treework[i] ;
      radius[node] = sqrt(radii[i]) ;
    }
      
    const int ntest = 32 ;
    const int nfit = 20 ;
    const int exponent = 3 ;
    const float factor = 3.0 ;
    vector<double> ftesti(ntest*treework.size()) ;
    for(size_t i=0;i<treework.size();++i) {
      int node = treework[i] ;
      int sloc = max(start,nodeDataTree[node].start) ;
      double rtest = radius[node]*factor ;
      int eloc = min(end,nodeDataTree[node].end) ;
      if(sloc < eloc)
	get_fexact(nodeDataSort,
		   sloc,eloc,exponent,rtest,
		   centroid[node],
		   testpnts,&ftesti[ntest*i],ntest) ;
    }
    vector<double> ftest(ntest*treework.size()) ;
    MPI_Allreduce(&ftesti[0],&ftest[0],ntest*treework.size(),MPI_DOUBLE,
		  MPI_SUM,comm) ;
    for(size_t i=0;i<treework.size();++i) {
      int node = treework[i] ;
      double rtest = radius[node]*factor ;
      double rw = 1./weights[node] ;
      for(int j=0;j<ntest;++j)
	ftest[ntest*i+j] *= rw ;
      makeQuad(nodeDataSort,
	       nodeDataTree[node].start,
	       nodeDataTree[node].end,
	       testpnts,&ftest[ntest*i],nfit,ntest,exponent,
	       centroid[node],rtest,
	       qpts[node][0],
	       qpts[node][1],
	       qpts[node][2],
	       qerr[node]) ;
    }
  }
  
  

  void buildDefApprox(deformApproxTree &approxTree,
		      const vector<pair<int,int> > &proclims,
		      const vector<int> &proctreebase,
		      const vector<int> &treeptn,
		      const vector<double> &weights,
		      MPI_Comm bcomm) {
    int pbuild = 1 ;
    int r = 0 ;
    MPI_Comm_size(bcomm,&pbuild) ;
    MPI_Comm_rank(bcomm,&r) ;
    int tsz = approxTree.nodeDataTree.size() ;
    approxTree.qdisp = vector<Array<vector3d<float>,4> >(tsz) ;
    approxTree.qrot  = vector<Array<tensor3d<float>,4> >(tsz) ;
    approxTree.drot  = vector<float>(tsz) ;
    approxTree.computeDeformationQuadsTree(weights,
					   proctreebase[r]) ;

    int start = approxTree.nodeDataTree[proctreebase[r]].start ;
    int end = approxTree.nodeDataTree[proctreebase[r]].end ;
    vector<int> treework ;
    for(size_t i=0;i<proclims.size();++i) 
      if(proclims[i].second-proclims[i].first > 1)
	treework.push_back(i) ;
    const int exponent = 3 ;
    const float factor = 3.0 ;
    // now communicate partial tree 
    if(pbuild > 1) {
      vector<float> buf(treeptn[r+1]-treeptn[r]) ;
      for(int i=treeptn[r];i<treeptn[r+1];++i)
	buf[i-treeptn[r]] =approxTree.drot[i] ;
      vector<int> recv_counts(pbuild) ;
      for(int i=0;i<pbuild;++i)
	recv_counts[i] = (treeptn[i+1]-treeptn[i]) ;
      vector<int> displs(pbuild) ;
      displs[0] = 0 ;
      for(int i=1;i<pbuild;++i)
	displs[i] = displs[i-1]+recv_counts[i-1] ;
      MPI_Allgatherv(&buf[0],buf.size(),MPI_FLOAT,
		     &approxTree.drot[treeptn[0]],
		     &recv_counts[0],&displs[0], MPI_FLOAT,
		     bcomm) ;
      float tmp = approxTree.drot[proctreebase[r]] ;
      vector<float> buf2(pbuild) ;
      MPI_Allgather(&tmp,1,MPI_FLOAT,&buf2[0],1,MPI_FLOAT,bcomm) ;
      for(int i=0;i<pbuild;++i)
	approxTree.drot[proctreebase[i]] = buf2[i] ;
    }
    if(pbuild > 1) {
      vector<Array<vector3d<float>, 4> > buf(treeptn[r+1]-treeptn[r]) ;
      for(int i=treeptn[r];i<treeptn[r+1];++i)
	buf[i-treeptn[r]] =approxTree.qdisp[i] ;
      vector<int> recv_counts(pbuild) ;
      for(int i=0;i<pbuild;++i)
	recv_counts[i] = (treeptn[i+1]-treeptn[i])*12 ;
      vector<int> displs(pbuild) ;
      displs[0] = 0 ;
      for(int i=1;i<pbuild;++i)
	displs[i] = displs[i-1]+recv_counts[i-1] ;
      MPI_Allgatherv(&buf[0],buf.size()*12,MPI_FLOAT,
		     &approxTree.qdisp[treeptn[0]],
		     &recv_counts[0],&displs[0], MPI_FLOAT,
		     bcomm) ;
      vector<Array<vector3d<float>, 4> > buf2(pbuild) ;
      Array<vector3d<float>,4> tmp = approxTree.qdisp[proctreebase[r]] ;
      MPI_Allgather(&tmp,12,MPI_FLOAT,&buf2[0],12,MPI_FLOAT,bcomm) ;
      for(int i=0;i<pbuild;++i) 
	approxTree.qdisp[proctreebase[i]] = buf2[i] ;
    }
    if(pbuild > 1) {
      vector<Array<tensor3d<float>, 4> > buf(treeptn[r+1]-treeptn[r]) ;
      for(int i=treeptn[r];i<treeptn[r+1];++i)
	buf[i-treeptn[r]] =approxTree.qrot[i] ;
      vector<int> recv_counts(pbuild) ;
      for(int i=0;i<pbuild;++i)
	recv_counts[i] = (treeptn[i+1]-treeptn[i])*36 ;
      vector<int> displs(pbuild) ;
      displs[0] = 0 ;
      for(int i=1;i<pbuild;++i)
	displs[i] = displs[i-1]+recv_counts[i-1] ;
      MPI_Allgatherv(&buf[0],buf.size()*36,MPI_FLOAT,
		     &approxTree.qrot[treeptn[0]],
		     &recv_counts[0],&displs[0], MPI_FLOAT,
		     bcomm) ;
      vector<Array<tensor3d<float>, 4> > buf2(pbuild) ;
      Array<tensor3d<float>,4> tmp = approxTree.qrot[proctreebase[r]] ;
      MPI_Allgather(&tmp,36,MPI_FLOAT,&buf2[0],36,MPI_FLOAT,bcomm) ;
      for(int i=0;i<pbuild;++i) 
	approxTree.qrot[proctreebase[i]] = buf2[i] ;
    }
    // now compute for parts shared by all processors
    int topsz = treework.size() ;
    // compute mean values for each node
    vector<Array<vector3d<double>, 4> > meanvals(topsz) ;
    for(int j=0;j<topsz;++j)
      for(int i=0;i<4;++i) 
	meanvals[j][i] = vector3d<double>(0,0,0) ;
    for(int j=0;j<topsz;++j) {
      int node = treework[j] ;
      vector3d<float> xc = approxTree.centroid[node] ;
      int sloc = max(start,approxTree.nodeDataTree[node].start) ;
      int eloc = min(end,approxTree.nodeDataTree[node].end) ;
      for(int i=sloc;i<eloc;++i) {
	tens3d rotmat = approxTree.nodeRotors[i].matrix() ;
	rotmat.x.x -= 1.0 ; // Convert rotation to delta change
	rotmat.y.y -= 1.0 ;
	rotmat.z.z -= 1.0 ;
	vector3d<double> dp(xc-approxTree.nodeDataSort[i].pos) ;
	vector3d<double> displT = vector3d<double>(approxTree.nodeDisp[i])+dot(rotmat,dp) ;
	double w = approxTree.nodeDataSort[i].weight ;
	meanvals[j][0] += w*displT ;
	meanvals[j][1] += w*rotmat.x ;
	meanvals[j][2] += w*rotmat.y ;
	meanvals[j][3] += w*rotmat.z ;
      }
    }
    {     
      vector<Array<vector3d<double>, 4> > meanvalsc(topsz) ;

      MPI_Allreduce(&meanvals[0],&meanvalsc[0],topsz*4*3,
		    MPI_DOUBLE,MPI_SUM,bcomm) ;
      meanvals.swap(meanvalsc) ;
    }
    for(int j=0;j<topsz;++j) {
      int node = treework[j] ;
      for(int i=0;i<4;++i) {
	meanvals[j][i] *= 1./weights[node] ;
      }
    }
      
    const int ntst = 4 ;
    vector<Array<vector3d<float>,ntst> > rtestpts(topsz) ;
    for(int j=0;j<topsz;++j) {
      int node = treework[j] ;
      vector3d<float> xc = approxTree.centroid[node] ;
      rtestpts[j][0] = approxTree.qpts[node][0]-xc ;
      rtestpts[j][1] = approxTree.qpts[node][1]-xc ;
      rtestpts[j][2] = approxTree.qpts[node][2]-xc ;
      rtestpts[j][3] = -1.*(rtestpts[j][0]+rtestpts[j][1]+rtestpts[j][2]) ;
      for(int i=0;i<ntst;++i)
	rtestpts[j][i] *= 1./max(norm(rtestpts[j][i]),1e-10f) ;
    }
    vector<Array<Array<vector3d<double>,ntst>,4> > defs(topsz) ;
    vector<Array<double,ntst> > wexactsum(topsz) ;
    for(int j=0;j<topsz;++j) {
      for(int i=0;i<ntst;++i)
	wexactsum[j][i] = 0.0 ;
      for(int i=0;i<ntst;++i)
	for(int k=0;k<4;++k)
	  defs[j][k][i] = vector3d<float>(0.0,0.0,0.0) ;

      int node = treework[j] ;
      vector3d<float> xc = approxTree.centroid[node] ;
      // Compute exact mapping to test points
      int sloc = max(start,approxTree.nodeDataTree[node].start) ;
      int eloc = min(end,approxTree.nodeDataTree[node].end) ;
      for(int i=sloc;i<eloc;++i) {
	tens3d rotmat = approxTree.nodeRotors[i].matrix() ;
	rotmat.x.x -= 1.0 ; // Convert rotation to delta change
	rotmat.y.y -= 1.0 ;
	rotmat.z.z -= 1.0 ;
	vector3d<double> dp(xc-approxTree.nodeDataSort[i].pos) ;
	vect3d displT = vector3d<double>(approxTree.nodeDisp[i])+
	  dot(rotmat,dp) ;
      
	for(int t=0;t<ntst;++t) {
	  vector3d<float> rtest = rtestpts[j][t]*approxTree.radius[node]*factor+xc ;
	  const double r = norm(approxTree.nodeDataSort[i].pos-rtest) ;
	  const double wn = approxTree.nodeDataSort[i].weight ;
	  const double wr = wn*pow(r,-exponent) ;
	  wexactsum[j][t] += wr ;
	  defs[j][0][t] += wr*displT ;
	  defs[j][1][t] += wr*rotmat.x ;
	  defs[j][2][t] += wr*rotmat.y ;
	  defs[j][3][t] += wr*rotmat.z ;
	}
      }
    }
    {
      vector<Array<double,ntst> > wexactsumc(topsz) ;
      MPI_Allreduce(&wexactsum[0],&wexactsumc[0],topsz*ntst,
		    MPI_DOUBLE,MPI_SUM,bcomm) ;
      wexactsum.swap(wexactsumc) ;
    }
    {
      vector<Array<Array<vector3d<double>,ntst>,4> > defsc(topsz) ;
      MPI_Allreduce(&defs[0],&defsc[0],topsz*ntst*4*3,
		    MPI_DOUBLE,MPI_SUM,bcomm) ;
      defs.swap(defsc) ;
    }
    for(int j=0;j<topsz;++j)
      for(int i=0;i<4;++i)
	for(int t=0;t<ntst;++t)
	  defs[j][i][t] *= 1./wexactsum[j][t] ;
      
     
    for(int j=0;j<topsz;++j) {
      int node = treework[j] ;

      vector3d<float> xc = approxTree.centroid[node] ;
      const vector3d<double> zero(0.0,0.0,0.0) ;
      float droterr = 0 ; 
      int nrsample = 4 ;
      double W[nrsample][4] ;
      double A00=0,A11=0,A22=0,A01=0,A02=0,A12= 0 ;
      for(int s=0;s<nrsample;++s) {
	// test point
	vector3d<float> pt = rtestpts[j][s]*approxTree.radius[node]*factor+xc ;

	// Weights associated with quadripole point for this test point
	for(int i=0;i<4;++i)
	  W[s][i] = 0 ;
	double Wsum = 0 ;
	for(int i=0;i<3;++i) {
	  const double r = norm(pt-approxTree.qpts[node][i]) ;
	  W[s][i] = pow(r,-exponent) ;
	  Wsum += W[s][i] ;
	}
	W[s][3] = pow(norm(pt-4.*xc+
			   approxTree.qpts[node][0]+
			   approxTree.qpts[node][1]+
			   approxTree.qpts[node][2]),
		      -exponent) ;
	Wsum += W[s][3] ;
	  
	double rWsum = 1./Wsum ;
	for(int i=0;i<4;++i)
	  W[s][i] *= rWsum ;
	double J[3] = {W[s][0]-W[s][3],W[s][1]-W[s][3],W[s][2]-W[s][3]} ;
	  
	A00 += J[0]*J[0] ;
	A11 += J[1]*J[1] ;
	A22 += J[2]*J[2] ;
	A01 += J[0]*J[1] ;
	A02 += J[0]*J[2] ;
	A12 += J[1]*J[2] ;
      }
      double det = (+A00*A11*A22+A01*A12*A02+A02*A01*A12
		    -A02*A11*A02-A00*A12*A12-A01*A01*A22) ;
      double ATAinv[3][3] ;
      double rdet = 1./(det+((det>0.0)?1e-12:-1e-12)) ;
      ATAinv[0][0] = rdet*(A11*A22 - A12*A12) ;
      ATAinv[0][1] = rdet*(A02*A12 - A01*A22) ;
      ATAinv[0][2] = rdet*(A01*A12 - A02*A11) ;
	
      ATAinv[1][0] = ATAinv[0][1] ;
      ATAinv[1][1] = rdet*(A00*A22 - A02*A02) ;
      ATAinv[1][2] = rdet*(A02*A01 - A00*A12) ;
	
      ATAinv[2][0] = ATAinv[0][2] ;
      ATAinv[2][1] = ATAinv[1][2] ;
      ATAinv[2][2] = rdet*(A00*A11 - A01*A01) ;
    
      for(int field=0;field<4;++field) { 
	// Now compute quad point values using least squares
	// Symmetric matrix terms (3x3)
	vector3d<double> Y[3] = {zero,zero,zero} ;
	for(int s=0;s<nrsample;++s) {
	  double J[3] = {W[s][0]-W[s][3],W[s][1]-W[s][3],W[s][2]-W[s][3]} ;
	    
	  for(int i=0;i<3;++i) {
	    Y[i] += J[i]*(defs[j][field][s]-4.*W[s][3]*meanvals[j][field]) ;
	  }
	}
	vector3d<double> v1= meanvals[j][field] ;
	vector3d<double> v2 = v1 ;
	vector3d<double> v3 = v1 ;
	vector3d<double> v4 = v1 ;

	if(fabs(det) > 1e-9) {
	  // Compute values at quad points
	  v1 = ATAinv[0][0]*Y[0]+ATAinv[0][1]*Y[1]+ATAinv[0][2]*Y[2] ;
	  v2 = ATAinv[1][0]*Y[0]+ATAinv[1][1]*Y[1]+ATAinv[1][2]*Y[2] ;
	  v3 = ATAinv[2][0]*Y[0]+ATAinv[2][1]*Y[1]+ATAinv[2][2]*Y[2] ;
	  v4 = 4.*meanvals[j][field] -v1-v2-v3;
	}
	double err[4] = {0,0,0,0} ;
	for(int s=0;s<ntst;++s) {
	  // test point
	  vector3d<float> pt = xc + rtestpts[j][s]*approxTree.radius[node]*factor ;
	  // Weights associated with quadripole point for this test point
	  double W[4] = {0,0,0,0} ;
	  double Wsum = 0 ;
	  for(int i=0;i<3;++i) {
	    const double r = norm(pt-approxTree.qpts[node][i]) ;
	    W[i] = pow(r,-exponent) ;
	    Wsum += W[i] ;
	  }
	  const double r3 = norm(pt-4.*xc+
				 approxTree.qpts[node][0]+
				 approxTree.qpts[node][1]+
				 approxTree.qpts[node][2]) ;
	  W[3] = pow(r3,-exponent) ;
	  Wsum += W[3] ;
	  double rWsum = 1./Wsum ;
	  for(int i=0;i<4;++i)
	    W[i] *= rWsum ;
	  err[field] = max(err[field],
			   norm(defs[j][field][s]-
				(W[0]*v1+W[1]*v2+W[2]*v3+W[3]*v4))) ;
	}
	switch(field) {
	case 0:
	  approxTree.qdisp[node][0] = v1 ;
	  approxTree.qdisp[node][1] = v2 ;
	  approxTree.qdisp[node][2] = v3 ;
	  approxTree.qdisp[node][3] = v4 ;
	  droterr = max(droterr,float(err[field]/(max(norm(meanvals[j][field]),1e-4)))) ;
	  break ;
	case 1:
	  approxTree.qrot[node][0].x = v1 ;
	  approxTree.qrot[node][1].x = v2 ;
	  approxTree.qrot[node][2].x = v3 ;
	  approxTree.qrot[node][3].x = v4 ;
	  droterr = max(droterr,float(err[field])) ;
	  break ;
	case 2:
	  approxTree.qrot[node][0].y = v1 ;
	  approxTree.qrot[node][1].y = v2 ;
	  approxTree.qrot[node][2].y = v3 ;
	  approxTree.qrot[node][3].y = v4 ;
	  droterr = max(droterr,float(err[field])) ;
	  break ;
	case 3:
	  approxTree.qrot[node][0].z = v1 ;
	  approxTree.qrot[node][1].z = v2 ;
	  approxTree.qrot[node][2].z = v3 ;
	  approxTree.qrot[node][3].z = v4 ;
	  droterr = max(droterr,float(err[field])) ;
	  break ;
	}
      }
      approxTree.drot[node]=droterr ;
    }
  }
  
  void weightApproxTree::build(const vector<vector3d<float> > &pnts,
                                const vector<float> &wts,
                                MPI_Comm comm) {
    // First gather
    int nsz = pnts.size() ;
    // copy in positions and weights
    nodeDataSort = vector<nodeDataBase>(nsz) ;
    for(int i=0;i<nsz;++i)  {
      nodeDataSort[i].pos = pnts[i] ;
      nodeDataSort[i].weight = wts[i] ;
    }
    
    stopWatch st ;
    st.start() ;
    pbuild =  0 ;
    buildTreeSerialPart(nodeDataSort,nodeDataTree,
			pbuild,proclims,proctreebase,treeptn,
			comm) ;
    reportTime("serial part of tree build",st.stop()) ;
    st.start() ;
    int p = 1 ;
    int r = 0 ;
    MPI_Comm_size(comm,&p) ;
    MPI_Comm_rank(comm,&r) ;

    MPI_Comm bcomm ;
    int color = 1 ;
    if(r >=pbuild)
      color=0 ;

    MPI_Comm_split(comm,color,r,&bcomm) ;
    reportTime("data setup before approximations",st.stop()) ;
    st.start() ;

    if(r < pbuild) {    
      vector<double> weights ;
      buildWeightsApprox(nodeDataSort, nodeDataTree,
			 proclims,proctreebase,treeptn,
			 weights,
			 centroid,radius, qpts,qerr,
			 bcomm) ;
      reportTime("parallel weights approx tree",st.stop()) ;
      st.start() ;

      int tsz = nodeDataTree.size() ;

      cweight = vector<float>(tsz) ;
      for(int i=0;i<tsz;++i) {
	cweight[i] = weights[i] ;
      }

    }

    // send tree to processors that did not compute tree
    if(r<pbuild) {
      if(r == 0 && pbuild != p) {
	int tsz = nodeDataTree.size() ;
	MPI_Send(&cweight[0],tsz,MPI_FLOAT,pbuild,21,comm) ;
	MPI_Send(&centroid[0],tsz*3,MPI_FLOAT,pbuild,22,comm) ;
	MPI_Send(&radius[0],tsz,MPI_FLOAT,pbuild,23,comm) ;
	MPI_Send(&qpts[0],tsz*3*3,MPI_FLOAT,pbuild,24,comm) ;
	MPI_Send(&qerr[0],tsz,MPI_FLOAT,pbuild,25,comm) ;
	int flag = 0 ;
	MPI_Status stat ;
	MPI_Recv(&flag,1,MPI_INT,pbuild,29,comm,&stat) ;
      }
    } else {
      int tsz = nodeDataTree.size() ;
      cweight = vector<float>(tsz) ;
      centroid = vector<vector3d<float> >(tsz) ;
      radius = vector<float>(tsz) ;
      qpts = vector<Array<vector3d<float>,3> >(tsz) ;
      qerr = vector<float>(tsz) ;
      if(r == pbuild) {
	MPI_Status stat ;
	MPI_Recv(&cweight[0],tsz,MPI_FLOAT,0,21,comm,&stat) ;
	MPI_Recv(&centroid[0],tsz*3,MPI_FLOAT,0,22,comm,&stat) ;
	MPI_Recv(&radius[0],tsz,MPI_FLOAT,0,23,comm,&stat) ;
	MPI_Recv(&qpts[0],tsz*3*3,MPI_FLOAT,0,24,comm,&stat) ;
	MPI_Recv(&qerr[0],tsz,MPI_FLOAT,0,25,comm,&stat) ;
      }
      MPI_Bcast(&cweight[0],tsz,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&centroid[0],tsz*3,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&radius[0],tsz,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&qpts[0],tsz*3*3,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&qerr[0],tsz,MPI_FLOAT,0,bcomm) ;
      int flag = 1 ;
      MPI_Send(&flag,1,MPI_INT,0,29,comm) ;
    }
    MPI_Comm_free(&bcomm) ;
    reportTime("communication",st.stop()) ;
    st.start() ;
  }
  
  
  void weightApproxTree::
  WeightApprox(double &weight, double &werr,
               const vector3d<double> &pos, int node,
               int a, int b, double L, double alphaL,
               double errpn) const {
    const vector3d<double> dvn = pos-vector3d<double>(centroid[node]) ;
    const double dn = norm(dvn) ;
    //    if(node == 0) 
    //      cout << "errpn=" << errpn << endl ;
    if(dn > 1.5*radius[node]) {
      // evaluate if this node can be approximated

      // Estimate error
      const real rfac = dn/(3.*radius[node]) ;
      double err = (qerr[node])/(rfac*rfac) ;
      if(rfac < .9) // increase error estimate if we are too close
	err *= 10.0 ;
      // estimated weight
      const vector3d<float> posf(pos) ;
      const double r1 = 1./norm(posf-qpts[node][0]) ;
      const double r2 = 1./norm(posf-qpts[node][1]) ;
      const double r3 = 1./norm(posf-qpts[node][2]) ;
      //    qp4 = 4.*xc - qp1 -  qp2 - qp3 ;
      const double r4 = 1./norm(posf-4.*centroid[node]+
				qpts[node][0]+qpts[node][1]+qpts[node][2]) ;
      const double w1 = pow(L*r1,a)+pow(alphaL*r1,b) ;
      const double w2 = pow(L*r2,a)+pow(alphaL*r2,b) ;
      const double w3 = pow(L*r3,a)+pow(alphaL*r3,b) ;
      const double w4 = pow(L*r4,a)+pow(alphaL*r4,b) ;
      const double W = 0.25*cweight[node] ;
      const double wn = W*(w1+w2+w3+w4) ;
      err *= wn ;

      // Check to see if the error is within the acceptable range
      if((err+werr)/(weight+wn) < errpn) { 
	weight += wn ;
	werr += err ;
	return ;
      } 
    }
  
    // approximation has too much error so descend the tree
    const int left = nodeDataTree[node].left ;
    const int right = nodeDataTree[node].right ;
    if(left < 0 || right < 0) { // leaf node so evaluate all points
      int start = nodeDataTree[node].start ;
      int end = nodeDataTree[node].end ;
      for(int i= start;i<end;++i) {
	const vect3d npos   = nodeDataSort[i].pos;
	const double wn     = nodeDataSort[i].weight;
	const vect3d rvec   = (vect3d(pos) - npos);
      
	double rnrvec = 1./max(norm(rvec),1e-20) ;
	const double w1 = pow(L*rnrvec, a) ;
	const double w2 = pow(alphaL*rnrvec, b) ;
      
	const double w = wn*(w1 + w2);
	//	if(fabs(delta.x)+fabs(delta.y)+fabs(delta.z)>1e-12)
	//	  cout << "w=" << w << ",delta =" << delta << endl ;
	weight += w ;
      }
      return ;
    }

    // Otherwise descend tree (first in branch that contains our point)
    double diff = -nodeDataTree[node].split ;
    if(nodeDataTree[node].coord == 0)
      diff += pos.x ;
    if(nodeDataTree[node].coord == 1)
      diff += pos.y ;
    if(nodeDataTree[node].coord == 2)
      diff += pos.z ;

    if(diff<0) {
      WeightApprox(weight,werr, pos,left,  a,b,L,alphaL,errpn) ;
      WeightApprox(weight,werr, pos,right, a,b,L,alphaL,errpn) ;
    } else {
      WeightApprox(weight,werr, pos,right, a,b,L,alphaL,errpn) ;
      WeightApprox(weight,werr, pos,left,  a,b,L,alphaL,errpn) ;
    }
    return ;

  }
  
  void buildDeformApprox(deformApproxTree &approxTree,
			 const vector<NodeData> &nodeDataSort,
			 MPI_Comm comm) {
    stopWatch sttot ;
    sttot.start() ;
    // First gather
    int nsz = nodeDataSort.size() ;
    // copy in positions and weights
    approxTree.nodeDataSort = vector<nodeDataBase>(nsz) ;
    for(int i=0;i<nsz;++i)  {
      approxTree.nodeDataSort[i].pos = nodeDataSort[i].pos ;
      approxTree.nodeDataSort[i].weight = nodeDataSort[i].weight ;
    }
    
    stopWatch st ;
    st.start() ;
    int pbuild =  0 ;
    vector<pair<int,int> > proclims ;
    vector<int> proctreebase ;
    vector<int> treeptn ;
    buildTreeSerialPart(approxTree.nodeDataSort,
			approxTree.nodeDataTree,
			pbuild,proclims,proctreebase,treeptn,
			comm) ;
    reportTime("serial part of tree build",st.stop()) ;
    st.start() ;
    int p = 1 ;
    int r = 0 ;
    MPI_Comm_size(comm,&p) ;
    MPI_Comm_rank(comm,&r) ;

    // Now copy in the deformations into the tree
    approxTree.nodeDisp = vector<vector3d<float> >(nsz) ;
    approxTree.nodeRotors = vector<Rotor>(nsz) ;
    for(int i=0;i<nsz;++i) {
      int access = approxTree.nodeDataSort[i].order ;
      approxTree.nodeDisp[i] = nodeDataSort[access].disp ;
      approxTree.nodeRotors[i] = nodeDataSort[access].rot ;
    }

    MPI_Comm bcomm ;
    int color = 1 ;
    if(r >=pbuild)
      color=0 ;

    MPI_Comm_split(comm,color,r,&bcomm) ;
    reportTime("data setup before approximations",st.stop()) ;
    st.start() ;

    if(r < pbuild) {    
      vector<double> weights ;
      buildWeightsApprox(approxTree.nodeDataSort,
			 approxTree.nodeDataTree,
			 proclims,proctreebase,treeptn,
			 weights,
			 approxTree.centroid,approxTree.radius,
			 approxTree.qpts,approxTree.qerr,
			 bcomm) ;
      reportTime("parallel weights approx tree",st.stop()) ;
      st.start() ;

      buildDefApprox(approxTree,proclims,proctreebase,treeptn,weights,
		     bcomm) ;
      reportTime("deformation approx tree",st.stop()) ;
      st.start() ;
		     

      int tsz = approxTree.nodeDataTree.size() ;

      approxTree.cweight = vector<float>(tsz) ;
      for(int i=0;i<tsz;++i) {
	approxTree.cweight[i] = weights[i] ;
      }

    }

    // send tree to processors that did not compute tree
    if(r<pbuild) {
      if(r == 0 && pbuild != p) {
	int tsz = approxTree.nodeDataTree.size() ;
	MPI_Send(&approxTree.cweight[0],tsz,MPI_FLOAT,pbuild,21,comm) ;
	MPI_Send(&approxTree.centroid[0],tsz*3,MPI_FLOAT,pbuild,22,comm) ;
	MPI_Send(&approxTree.radius[0],tsz,MPI_FLOAT,pbuild,23,comm) ;
	MPI_Send(&approxTree.qpts[0],tsz*3*3,MPI_FLOAT,pbuild,24,comm) ;
	MPI_Send(&approxTree.qerr[0],tsz,MPI_FLOAT,pbuild,25,comm) ;
	MPI_Send(&approxTree.qdisp[0],tsz*3*4,MPI_FLOAT,pbuild,26,comm) ;
	MPI_Send(&approxTree.qrot[0],tsz*9*4,MPI_FLOAT,pbuild,27,comm) ;
	MPI_Send(&approxTree.drot[0],tsz,MPI_FLOAT,pbuild,28,comm) ;
	int flag = 0 ;
	MPI_Status stat ;
	MPI_Recv(&flag,1,MPI_INT,pbuild,29,comm,&stat) ;
      }
    } else {
      int tsz = approxTree.nodeDataTree.size() ;
      approxTree.cweight = vector<float>(tsz) ;
      approxTree.centroid = vector<vector3d<float> >(tsz) ;
      approxTree.radius = vector<float>(tsz) ;
      approxTree.qpts = vector<Array<vector3d<float>,3> >(tsz) ;
      approxTree.qerr = vector<float>(tsz) ;
      approxTree.qdisp = vector<Array<vector3d<float>,4> >(tsz) ;
      approxTree.qrot  = vector<Array<tensor3d<float>,4> >(tsz) ;
      approxTree.drot  = vector<float>(tsz) ;
      if(r == pbuild) {
	MPI_Status stat ;
	MPI_Recv(&approxTree.cweight[0],tsz,MPI_FLOAT,0,21,comm,&stat) ;
	MPI_Recv(&approxTree.centroid[0],tsz*3,MPI_FLOAT,0,22,comm,&stat) ;
	MPI_Recv(&approxTree.radius[0],tsz,MPI_FLOAT,0,23,comm,&stat) ;
	MPI_Recv(&approxTree.qpts[0],tsz*3*3,MPI_FLOAT,0,24,comm,&stat) ;
	MPI_Recv(&approxTree.qerr[0],tsz,MPI_FLOAT,0,25,comm,&stat) ;
	MPI_Recv(&approxTree.qdisp[0],tsz*3*4,MPI_FLOAT,0,26,comm,&stat) ;
	MPI_Recv(&approxTree.qrot[0],tsz*9*4,MPI_FLOAT,0,27,comm,&stat) ;
	MPI_Recv(&approxTree.drot[0],tsz,MPI_FLOAT,0,28,comm,&stat) ;
      }
      MPI_Bcast(&approxTree.cweight[0],tsz,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&approxTree.centroid[0],tsz*3,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&approxTree.radius[0],tsz,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&approxTree.qpts[0],tsz*3*3,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&approxTree.qerr[0],tsz,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&approxTree.qdisp[0],tsz*3*4,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&approxTree.qrot[0],tsz*9*4,MPI_FLOAT,0,bcomm) ;
      MPI_Bcast(&approxTree.drot[0],tsz,MPI_FLOAT,0,bcomm) ;
      int flag = 1 ;
      MPI_Send(&flag,1,MPI_INT,0,29,comm) ;
    }
    MPI_Comm_free(&bcomm) ;
    reportTime("communication",st.stop()) ;
    st.start() ;
    reportTime("Total Tree Build",sttot.stop()) ;
  }    

    
  
  void deformApproxTree::
  TreeApprox(vector3d<double> &disp, double &weight, double &werr,
	     const vector3d<double> &pos, int node,
	     int a, int b, double L, double alphaL,
	     double errpn) const {
    const vector3d<double> dvn = pos-vector3d<double>(centroid[node]) ;
    const double dn = norm(dvn) ;
    //    if(node == 0) 
    //      cout << "errpn=" << errpn << endl ;
    if(dn > 1.5*radius[node]) {
      // evaluate if this node can be approximated

      // Estimate error
      const real rfac = dn/(3.*radius[node]) ;
      //    double err = nodeDataTree[node].err/(rfac*rfac)+nodeDataTree[node].drot/(rfac*rfac) ;
      double err = (qerr[node]+drot[node])/(rfac*rfac) ;
      if(rfac < .9) // increase error estimate if we are too close
	err *= 10.0 ;
      // estimated weight
      const vector3d<float> posf(pos) ;
      const double r1 = 1./norm(posf-qpts[node][0]) ;
      const double r2 = 1./norm(posf-qpts[node][1]) ;
      const double r3 = 1./norm(posf-qpts[node][2]) ;
      //    qp4 = 4.*xc - qp1 -  qp2 - qp3 ;
      const double r4 = 1./norm(posf-4.*centroid[node]+
				qpts[node][0]+qpts[node][1]+qpts[node][2]) ;
      const double w1 = pow(L*r1,a)+pow(alphaL*r1,b) ;
      const double w2 = pow(L*r2,a)+pow(alphaL*r2,b) ;
      const double w3 = pow(L*r3,a)+pow(alphaL*r3,b) ;
      const double w4 = pow(L*r4,a)+pow(alphaL*r4,b) ;
      const double W = 0.25*cweight[node] ;
      const double wn = W*(w1+w2+w3+w4) ;
      err *= wn ;

      // Check to see if the error is within the acceptable range
      if((err+werr)/(weight+wn) < errpn) { 
	vector3d<float> pc = vector3d<float>(pos)-centroid[node] ;
	//	const vector3d<float> def = dot(cRotations[node],pc) + cDisp[node]  ;
	const vector3d<float> disp1 = dot(qrot[node][0],pc)+qdisp[node][0] ;
	const vector3d<float> disp2 = dot(qrot[node][1],pc)+qdisp[node][1] ;
	const vector3d<float> disp3 = dot(qrot[node][2],pc)+qdisp[node][2] ;
	const vector3d<float> disp4 = dot(qrot[node][3],pc)+qdisp[node][3] ;
	const vector3d<float> dispf = W*(w1*disp1+w2*disp2+w3*disp3+w4*disp4) ;
	weight += wn ;
	werr += err ;
	disp += vector3d<double>(dispf) ;
	//	cout << "werrsum = " << werr << " + " << err << endl ;
	return ;
      } 
      //    cout << " err = " << err/wn << ",wn=" << wn
      //	 <<",werr = " << werr << "w=" << weight 
      //	 << ",rfac =" << rfac << ",errpn="<< errpn << endl ;
    }
  
    // approximation has too much error so descend the tree
    const int left = nodeDataTree[node].left ;
    const int right = nodeDataTree[node].right ;
    if(left < 0 || right < 0) { // leaf node so evaluate all points
      int start = nodeDataTree[node].start ;
      int end = nodeDataTree[node].end ;
      //    cout << "cnt=" << cnt << "werr="<< werr << endl ;
      for(int i= start;i<end;++i) {
	const Rotor  R      = nodeRotors[i] ;
	const vect3d npos   = nodeDataSort[i].pos;
	const vect3d ndisp  = nodeDisp[i] ;
	const double wn     = nodeDataSort[i].weight;
	const vect3d rvec   = (vect3d(pos) - npos);
      
	const vect3d delta = ( (R*rvec) + ndisp - rvec );
      
	double rnrvec = 1./max(norm(rvec),1e-20) ;
	const double w1 = pow(L*rnrvec, a) ;
	const double w2 = pow(alphaL*rnrvec, b) ;
      
	const double w = wn*(w1 + w2);
	//	if(fabs(delta.x)+fabs(delta.y)+fabs(delta.z)>1e-12)
	//	  cout << "w=" << w << ",delta =" << delta << endl ;
	disp += w*delta ;
	weight += w ;
      }
      return ;
    }

    // Otherwise descend tree (first in branch that contains our point)
    double diff = -nodeDataTree[node].split ;
    if(nodeDataTree[node].coord == 0)
      diff += pos.x ;
    if(nodeDataTree[node].coord == 1)
      diff += pos.y ;
    if(nodeDataTree[node].coord == 2)
      diff += pos.z ;

    if(diff<0) {
      TreeApprox(disp,weight,werr, pos,left,  a,b,L,alphaL,errpn) ;
      TreeApprox(disp,weight,werr, pos,right, a,b,L,alphaL,errpn) ;
    } else {
      TreeApprox(disp,weight,werr, pos,right, a,b,L,alphaL,errpn) ;
      TreeApprox(disp,weight,werr, pos,left,  a,b,L,alphaL,errpn) ;
    }
    //    if(node == 0) 
    //      cout << "w="<< weight << ",werr="<< werr << "disp=" << disp/weight << endl ;
    return ;

  }

}
