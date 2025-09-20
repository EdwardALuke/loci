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
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

namespace geomRep {
  using namespace Loci ;
  
  using std::ios ;
  using std::vector ;
  struct cubicBezierTriangle {
    vector3d<double> b300,b030,b003 ;
    vector3d<double> b210,b120,b021 ;
    vector3d<double> b012,b102,b201 ;
    vector3d<double> b111 ;
    vector3d<double> centroid() { return (1./3.)*(b300+b030+b003) ; }
    double radius() { 
      vector3d<double> c = centroid() ; 
      return max(max(norm(c-b300),norm(c-b030)),norm(c-b003)) ;
    } 
    void projectLinearUV(double &u, double &v, vector3d<double> P) const {
      // triangle normal vector
      vector3d<double> n = cross(b030-b300,b003-b300) ;
      n *= 1./(norm(n)+1e-30) ;
    
      double a0 = dot(cross(b030-P,b003-P),n) ; // Compute projected areas
      double a1 = dot(cross(b003-P,b300-P),n) ;
      double a2 = dot(cross(b300-P,b030-P),n) ;
      if(a0 < 0) {
        // between b030 and b003 so u=0
        v = min(max(dot(P-b003,b030-b003)/
                    (dot(b030-b003,b030-b003)+1e-30),0.0),1.0) ;
        u = 0 ;
      } else if(a1 < 0) {
        // between b003 and b300 so v=0 
        u = min(max(dot(P-b003,b300-b003)/
                    (dot(b300-b003,b300-b003)+1e-30),0.0),1.0) ;
        v = 0 ;
      } else if(a2 < 0) {
        // between b300 and b030 so w=0
        u = min(max(dot(P-b030,b300-b030)/
                    (dot(b300-b030,b300-b030)+1e-30),0.0),1.0) ;
        v = 1.0-u ;
      } else {
        a0 = max(a0,0.0) ;
        a1 = max(a1,0.0) ;
        a2 = max(a2,0.0) ;
        double s = a0+a1+a2 ; // Scale barycentric coordinates to [0,1]
        u = a0/s ;
        v = a1/s ;
      }
    }
    vector3d<double> evaluateLinear(double u, double v) const {
      return u*b300+v*b030+(1.0-u-v)*b003 ;
    }
    vector3d<double> evaluateBezier(double u, double v) const {
      const double w = 1.0-u-v ;
      return (u*u*u*b300 + v*v*v*b030 + w*w*w*b003+ 
              3.*(u*u*v*b210 + u*v*v*b120 + u*u*w*b201 +
                  v*v*w*b021 + u*w*w*b102 + v*w*w*b012) +
              6.*w*u*v*b111);
    }
    vector3d<double> evaluateBezierDu(double u, double v) const {
      const double w = 1.0-u-v ;
      return 3.*(u*u*(b300-b201)  + v*v*(b120-b021) + w*w*(b102-b003))
        + 6.*((b210-b111)*u*v + (b201-b102)*u*w + (b111-b012)*v*w) ;
    }
    vector3d<double> evaluateBezierDv(double u, double v) const {
      const double w = 1.0-u-v ;
      return 3.*(u*u*(b210-b201)+v*v*(b030-b021)+w*w*(b012-b003))
        + 6.*((b120-b111)*u*v + (b111-b102)*u*w + (b021-b012)*v*w) ;
    }

    void projectBezierUV(double &u, double &v, const vector3d<double> &P) const ;

  } ;    

  inline std::istream &operator>>(std::istream &s, cubicBezierTriangle &t) {
    s >> t.b300 >> t.b030 >> t.b003 >> t.b210 >> t.b120 >> t.b021 >> t.b012 >> t.b102 >> t.b201 >> t.b111 ;
    return s ;
  }

  inline std::ostream &operator<<(std::ostream &s, const cubicBezierTriangle &t) {
    s << t.b300 << t.b030 << t.b003 << t.b210 << t.b120 << t.b021 << t.b012 << t.b102 << t.b201 << t.b111 ;
    return s ;
  }

  class surfaceGeometry {
    int npnts ;
    vector<vector3d<double> > pos ;
    int nfaces ;
    vector<Array<int,3> > triNodes ;
    vector<cubicBezierTriangle> triGeom ;
    Loci::kdTree::kd_tree  *kdp ;
    vector<int> searchMap ;
    vector<int> offsetMap ;
  public:
    surfaceGeometry() {
      kdp=0 ;
      npnts= -1 ;
      nfaces= -1 ;
    }
    ~surfaceGeometry() {
      if(kdp)
        delete kdp ;
    }
    bool readASCIIFile(string filename) ;
    void setupSearchMap() ;
    void broadcastGeometry(int root, MPI_Comm comm) ;
    vector3d<double> projectGeometryLinear(vector3d<double> pt, vector3d<double> &normal) const ;
    vector3d<double> projectGeometryBezier(vector3d<double> pt, vector3d<double> &normal) const ;
    int getNumNodes() { return npnts ; }
    int getNumFaces() { return nfaces ; }

  } ;
  
}

