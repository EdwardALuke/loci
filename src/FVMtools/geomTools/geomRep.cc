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
#include <Loci.h>
#include <iostream>
#include <fstream>
#include <string>
#include "geomRep.h"


namespace geomRep {
  using namespace Loci ;
  using namespace Loci::kdTree;
  using std::vector ;
  using std::ifstream ;
  using std::ios ;
  using std::cerr ;
  using std::endl ;
  using std::cout ;
  using std::pair;
  using std::string ;
  typedef  Loci::vector3d<double> vect3d;
  
  inline void normalizeUV(double &u, double &v) {
    double w = max(0.0,1.-u-v) ;
    u = max(0.0,u) ;
    v = max(0.0,v) ;
    double s = 1./(w+u+v) ;
    u *= s ;
    v *= s ;
  }
  
  void cubicBezierTriangle:: projectBezierUV(double &u, double &v,
                                             const vector3d<double> &P) const {
    // Get initial guess
    projectLinearUV(u, v, P) ;
    int MAX_ITER = 32 ;
    int i=0 ;
    for(i=0;i<MAX_ITER;++i) {
      vector3d<double> Ps = evaluateBezier(u,v) ;
      vector3d<double> Psu = evaluateBezierDu(u,v) ;
      double Fu = dot(P-Ps,Psu) ;
      double dFu = -dot(Psu,Psu) ;
      double uold = u ;
      u += -Fu/dFu ;
      vector3d<double> Psv = evaluateBezierDv(u,v) ;
      double Fv = dot(P-Ps,Psv) ;
      double dFv = -dot(Psv,Psv) ;
      double vold = v ;
      v += -Fv/dFv ;
      normalizeUV(u,v) ;
      double err = fabs(uold-u)+fabs(vold-v) ;
      //    cout << "i=" << i << ",u="<< u << ",v=" << v << ",err="<< err << endl ;
      if(err < 1e-5)
        break ;
    }
  }


  bool surfaceGeometry::readASCIIFile(string filename) {
    npnts = -1 ;
    nfaces = -1 ;
    ifstream infile(filename.c_str(),ios::in) ;
    if(infile.fail()) 
      return false ;
    infile >> npnts ;
    {  vector<vector3d<double> > lpos(npnts) ;
      for(int i=0;i<npnts;++i) {
        double x, y, z ;
        infile >> x >> y >> z ;
        if(infile.fail()) 
          return false ;
        lpos[i] = vector3d<double>(x,y,z) ;
      }
      pos.swap(lpos) ;
    }
  
    infile >> nfaces ;

    if(infile.fail()) {
      npnts = -1 ;
      nfaces = -1 ;
      return false ;
    }
    { vector<Array<int,3> > ltriNodes(nfaces) ;
      for(int i=0;i<nfaces;++i) {
        infile >> ltriNodes[i][0] >> ltriNodes[i][1] >> ltriNodes[i][2] ;
        if(infile.fail()) {
          npnts = -1 ;
          nfaces = -1 ;
          return false ;
        }
      }
      triNodes.swap(ltriNodes) ;
    }
    { vector<cubicBezierTriangle> ltriGeom(nfaces) ;
      for(int i=0;i<nfaces;++i) {
        infile >> ltriGeom[i] ;
        if(infile.fail()) {
          npnts = -1 ;
          nfaces = -1 ;
          return false ;
        }
      }
      triGeom.swap(ltriGeom) ;
    } 
    if(infile.fail())
      return false ;

    setupSearchMap() ;

    return true ;
  }

  void surfaceGeometry::setupSearchMap() {
    if(kdp)
      delete kdp ;

    vector<Loci::kdTree::coord3d> fpos(nfaces) ;
    vector<int> fid(nfaces) ;
    for(int i=0;i<nfaces;++i) {
      fpos[i] = triGeom[i].centroid() ;
      fid[i] = i ; ;
    }
    kdp = new Loci::kdTree::kd_tree(fpos,fid) ;

    vector<pair<int,int> > searchPairs ;
    for(int i=0;i<nfaces;++i) {
      Loci::kdTree::coord3d c = triGeom[i].centroid() ;
      double r = triGeom[i].radius() ;
      Loci::kdTree::kd_tree::bounds box ;
      double r2 = 2.05*r ;
      box.maxc[0] = c.x + r2 ;
      box.maxc[1] = c.y + r2 ;
      box.maxc[2] = c.z + r2 ;
      box.minc[0] = c.x - r2 ;
      box.minc[1] = c.y - r2 ;
      box.minc[2] = c.z - r2 ;
      // Find all of the points within a given bounding box
      std::vector< Loci::kdTree::kd_tree::coord_info> found_pts;
      kdp->find_box(found_pts,  box);
    
      for(unsigned int pi = 0; pi < found_pts.size(); pi++){
        vect3d p = vect3d(found_pts[pi].coords[0], found_pts[pi].coords[1], found_pts[pi].coords[2]);
        if(norm(p - c) <= r2){
          pair<int,int> e1(found_pts[pi].id,i) ;
          pair<int,int> e2(i,found_pts[pi].id) ;
          searchPairs.push_back(e1) ;
          searchPairs.push_back(e2) ;
        }
      }
    }
    sort(searchPairs.begin(),searchPairs.end()) ;
    vector<pair<int,int> >::iterator ue = 
      unique(searchPairs.begin(),searchPairs.end()) ;
    int sz = ue-searchPairs.begin() ;
    
    vector<int> lsearchMap(sz) ;
    vector<int> sizes(nfaces,0) ;
    for(int i=0;i<sz;++i) {
      lsearchMap[i] = searchPairs[i].second ;
      sizes[searchPairs[i].first]++ ;
    }
    vector<int> loffsets(nfaces+1,0) ;
    for(int i=0;i<nfaces;++i)
      loffsets[i+1] = loffsets[i]+sizes[i] ;
    searchMap.swap(lsearchMap) ;
    offsetMap.swap(loffsets) ;
  }

  void surfaceGeometry:: broadcastGeometry(int root, MPI_Comm comm) {
    MPI_Bcast(&npnts,1,MPI_INT,root,comm) ;
    MPI_Bcast(&nfaces,1,MPI_INT,root,comm) ;

    int r ;
    MPI_Comm_rank(comm,&r) ;
    if(r != root) {
      vector<vector3d<double> > lpos(npnts) ;
      pos.swap(lpos) ;
      vector<Array<int,3> > ltriNodes(nfaces) ;    
      triNodes.swap(ltriNodes) ;
      vector<cubicBezierTriangle> ltriGeom(nfaces) ;
      triGeom.swap(ltriGeom) ;
    }
    MPI_Bcast(&pos[0],npnts*3,MPI_DOUBLE,root,comm) ;
    MPI_Bcast(&triNodes[0],nfaces*3,MPI_INT,root,comm) ;
    MPI_Bcast(&triGeom[0],nfaces*30,MPI_DOUBLE,root,comm) ;
    if(r != root) {
      setupSearchMap() ;
    }
  }

  vector3d<double> surfaceGeometry::projectGeometryLinear(vector3d<double> pt, vector3d<double> &normal) const {
    Loci::kdTree::coord3d thePoint ;
    thePoint[0] = pt.x ;
    thePoint[1] = pt.y ;
    thePoint[2] = pt.z ;
    double rmin = std::numeric_limits<float>::max() ;
    int id = kdp->find_closest(thePoint,rmin) ;
    double dmin = 1e30 ;
    vect3d surface_pt(0,0,0) ;
    double surface_u=0,surface_v=0 ;
    int surface_t = -1 ;
  
    for(int i=offsetMap[id];i<offsetMap[id+1];++i) {
      int f = searchMap[i] ;
      double u,v ;
      triGeom[f].projectLinearUV(u,v,pt) ;
      vect3d sp = triGeom[f].evaluateLinear(u,v) ;
      double d = norm(sp-pt) ;
      if(d< dmin) {
        dmin = d ;
        surface_pt = triGeom[f].evaluateBezier(u,v) ;
        surface_t = f ;
        surface_u = u ;
        surface_v = v ;
      }
    }
    normal = cross(triGeom[surface_t].b030-triGeom[surface_t].b300,
                   triGeom[surface_t].b003-triGeom[surface_t].b300) ;
    normal *= 1./max(norm(normal),1e-30) ;
    return surface_pt ;
  }
  vector3d<double> surfaceGeometry::projectGeometryBezier(vector3d<double> pt, vector3d<double> &normal) const {
    Loci::kdTree::coord3d thePoint ;
    thePoint[0] = pt.x ;
    thePoint[1] = pt.y ;
    thePoint[2] = pt.z ;
    double rmin = std::numeric_limits<float>::max() ;
    int id = kdp->find_closest(thePoint,rmin) ;
    double dmin = 1e30 ;
    vect3d surface_pt(0,0,0) ;
    double surface_u=0,surface_v=0 ;
    int surface_t = -1 ;
  
    for(int i=offsetMap[id];i<offsetMap[id+1];++i) {
      int f = searchMap[i] ;
      double u,v ;
      triGeom[f].projectBezierUV(u,v,pt) ;
      vect3d sp = triGeom[f].evaluateBezier(u,v) ;
      double d = norm(sp-pt) ;
      if(d< dmin) {
        dmin = d ;
        surface_pt = sp ;
        surface_t = f ;
        surface_u = u ;
        surface_v = v ;
      }
    }
    vector3d<double> dpu = triGeom[surface_t].evaluateBezierDu(surface_u,surface_v) ;
    vector3d<double> dpv = triGeom[surface_t].evaluateBezierDv(surface_u,surface_v) ;
    normal = cross(dpu,dpv) ;
    normal *= 1./max(norm(normal),1e-30) ;
    return surface_pt ;
  }
}

