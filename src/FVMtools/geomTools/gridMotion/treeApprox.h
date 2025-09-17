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
//# along with the GridMover software. If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef APPROX_H
#define APPROX_H

#include <istream>
#include <ostream>
#include <iostream>
using std::cerr;
using std::endl;
using std::ostream ;
#include <iomanip>
using std::setw;
#include <vector>
using std::vector;

#include "gridTypes.h"
#include "nodeData.h"
#include "rotor.h"
#include "affine.h"
#include "nodeTree.h"
#include <mpi.h>

namespace gridMotion {
  
  class weightApproxTree {
  public:
    // node Data
    vector<nodeDataBase> nodeDataSort ;
    // tree Data 
    vector<treeInfoBase> nodeDataTree ;
    // parallel info
    int pbuild;
    vector<pair<int,int> > proclims ;
    vector<int> proctreebase,treeptn ;
    
    // weight approximations
    vector<float> cweight ;
    vector<vector3d<float> > centroid ;
    vector<float> radius ;
    vector<Array<vector3d<float>,3> > qpts ;
    vector<float> qerr ;
    void build(const vector<vector3d<float> > &pnts,
               const vector<float> &wts,
               MPI_Comm comm) ; 
    void WeightApprox(double &weight, double &werr,
                      const vector3d<double> &pos, int node, 
                      int a, int b, double L, double alphaL, double errpn) const ;
 } ;    
  class deformApproxTree : public weightApproxTree {
  public:

    // Deformation Info
    vector<vector3d<float> > nodeDisp ;
    vector<Rotor> nodeRotors ;
    // deformation approximations
    vector<Array<vector3d<float>,4> > qdisp ;
    vector<Array<tensor3d<float>,4> > qrot ;
    vector<float> drot ;
      

    void computeDeformationQuadsTree(const vector<double> &weights,
				     int node) ;

    void TreeApprox(vector3d<double> &disp, double &weight, double &werr,
		    const vector3d<double> &pos, int node, 
		    int a, int b, double L, double alphaL, double errpn) const ;
  } ;
    

  void buildDeformApprox(deformApproxTree &approxTree,
			 const vector<NodeData> &nodeDataSort,
			 MPI_Comm comm) ;
      
  // Return the displacement at location pos meeting the absolute error bound
  // (error no greater than the prescribed value).
  inline vector3d<double>
    approxDisplacement(const vector3d<double> &pos,
		       const deformApproxTree &approxTree,
		       int a,
		       int b,
		       double L,
		       double aL,
		       double err_abs) {
    double w = 0 ;
    double werr = 0 ;
    vector3d<double> disp(0,0,0) ;
    vector3d<double> posf(pos.x,pos.y,pos.z) ;
    double errpn = 0.1 ; // ten percent error target to get cheap estimate of
    // displacement
    approxTree.TreeApprox(disp,w,werr,posf,0,a,b,L,aL,errpn) ;
    vector3d<double> d = disp ;
    d *= 1./w ;
    // Now we have the cheap displacement and an estimate of the displacement
    // magnitude, so we can convert the absolute error into a relative
    // error
    double disp_mag = std::max(fabs(d.x)+fabs(d.y)+fabs(d.z),1e-9) ;
    errpn = err_abs/disp_mag ;
    // Now check to see if we already met this relative error
    if(errpn < werr/w) {// no so redo approximation
      w = 0 ;
      werr = 0 ;
      disp = vector3d<double>(0.,0.,0.) ;
      approxTree.TreeApprox(disp,w,werr,posf,0,a,b,L,aL,errpn) ;
      vector3d<double> d = disp ;
      d *= 1./w ;
    }
    return d ;
  }
                           
  void registerReporterFunction(void (*rfp)(std::string )) ;
  void reportTime(const char *text, double time) ;
}
#endif
