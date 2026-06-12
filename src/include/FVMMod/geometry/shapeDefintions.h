//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
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
#include <sciTypes.h>

#include <map>
using std::map ;
#include <vector>
using std::vector ;
using std::cout ;
using std::cerr ;
using std::endl ;
using Loci::storeRepP ;
using Loci::rigid_transform ;

#include "geomRep.h"

namespace loci {

  /// Implicit Function Base Class that defines geometry API that computes
  /// distance function from a composition of fundamental geometric elements
  class distFunc : public Loci::CPTR_type {
  public:
    /// Pure virtual function for obtaining distance and gradient of distance
    /// field evaluated at a provided location.  Defined such that it is
    /// evaluated for an array of locations of size sz.
    virtual void getDist(real dist[], vector3d<real_t> n[],
                         const vector3d<real_t> loc[],int sz) const = 0 ;
    /// Transform geometry using rigid body transform
    virtual Loci::CPTR<distFunc> transform(rigid_transform Tv) const = 0 ;
  } ;


  /// The composite function builds a geometry from a set of geometric
  /// components.  You are outside of the geometry ony if you are outside
  /// of all geometric components (e.g. compose as the union of all
  /// geometric elements).  Distance outside of the component is the minimum
  /// distance out of all of the geometric elements.
  class compositeFunc: public distFunc {
    /// List of composed geometries
    std::vector<Loci::CPTR<distFunc> > geomList ;

  public:
    compositeFunc(std::vector<Loci::CPTR<distFunc> > &list) {
      geomList = list ;
      if(geomList.size() == 0)
        cerr << "warning composite function created with no geometry!"
             << endl ;
    }
    virtual void getDist(real dist[], vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
  
    virtual Loci::CPTR<distFunc> transform(rigid_transform Tv) const  ;
  } ;


  /// Compute distance to composite geometry
  void compositeFunc::getDist(real dist[], vector3d<real_t> n[],
                              const vector3d<real_t> loc[], int sz) const  {
    vector<real> dists(sz) ;
    vector<vector3d<real_t>> ns(sz) ;

    /// First get distance of first element in list
    geomList[0]->getDist(dist,n,loc,sz) ;
    /// Loop over remaining geometric elements 
    for(size_t i=1;i<geomList.size();++i) {
      geomList[i]->getDist(&dists[0],&ns[0],loc,sz) ;
      /// Select the closest distance out of the composite list
      for(int j=0;j<sz;++j) {
        /// Define the distance field as the minimum of all distances
        /// in the list of geometric elements
        if(dists[j] < dist[j]) {
          dist[j] = dists[j] ;
          n[j] = ns[j] ;
        }
      }
    }
  }

  /// Apply geometric transform to each element
  Loci::CPTR<distFunc> compositeFunc::transform(rigid_transform Tv) const {
    std::vector<Loci::CPTR<distFunc> > newGeomList(geomList.size()) ;
    for(size_t i=0;i<geomList.size();++i)
      newGeomList[i] = geomList[i]->transform(Tv) ;
    return new compositeFunc(newGeomList) ;
  }


  namespace hexFuncInfo {
    const int triangles[12][3] = {{0, 1, 2}, {3, 2, 1},
                                  {4, 6, 5}, {7, 5, 6},
                                  {1, 0, 5}, {4, 5, 0},
                                  {3, 7, 2}, {6, 2, 7},
                                  {2, 6, 0}, {4, 0, 6},
                                  {3, 1, 7}, {5, 7, 1}} ;
    const int t2e[12][3] = {{ 0, 1, 2},{ 4, 1, 3},
                            { 7, 6, 5},{ 8, 6, 9},
                            { 0,14,13},{ 5,14,10},
                            {12,17, 4},{11,17, 9},
                            {11,15, 2},{10,15, 7},
                            { 3,16,12},{ 8,16,13}} ;
			  
    const int edges[18][2] = {{0,1},{1,2},{2,0},{1,3},{3,2},
                              {4,5},{5,6},{6,4},{5,7},{7,6},
                              {0,4},{2,6},{3,7},{1,5},
                              {0,5},{0,6},{1,7},{2,7}} ;
  }


  /// Define a geometry using 8 points in space.
  class hexFunc: public distFunc {
    vector3d<real_t> corners[8] ;
    vector3d<real_t> fn[12] ;
    vector3d<real_t> ev[18] ;
    real evl[18] ;
    real radius ;
  public:
    hexFunc(vector3d<real_t> pts[8],real r) {
      using namespace hexFuncInfo ;
      for(int i=0;i<8;++i)
        corners[i] = pts[i] ;
      for(int f=0;f<12;++f) {
        vector3d<real_t> dv1 = corners[triangles[f][1]]-corners[triangles[f][0]] ;
        vector3d<real_t> dv2 = corners[triangles[f][2]]-corners[triangles[f][0]] ;
        vector3d<real_t> A = cross(dv1,dv2) ;
        fn[f] = A/max<real>(norm(A),1e-30) ;
      }
      for(int e=0;e<18;++e) {
        vector3d<real_t> dv = corners[edges[e][1]]-corners[edges[e][0]] ;
        evl[e] = norm(dv) ;
        ev[e] = dv/max<real>(evl[e],1e-30) ;
      }
      radius = r ;
    }
    virtual void getDist(real dist[],vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
    virtual Loci::CPTR<distFunc> transform(rigid_transform Tv) const ;
  } ;

  /// Get distance field by projecting onto points, then edges, then triangular
  /// faces of the decomposed hex
  void hexFunc::getDist(real dist[],vector3d<real_t> n[],
                        const vector3d<real_t> loc[], int sz) const {
    using namespace hexFuncInfo ;
    for(int i=0;i<sz;++i) { // loop over points
      vector3d<real_t> pt = loc[i] ;
      // first project onto all of the edges
      real mindist2 = 1e30 ;
      vector3d<real_t> minvec(0,0,0) ;
      for(int e=0;e<18;++e) {
        real w = dot(pt-corners[edges[e][0]],ev[e])/max<real>(evl[e],1e-15) ;
        w = max<real>(0.0,min<real>(1.0,w)) ;
        vector3d<real_t> dv = pt-((1.-w)*corners[edges[e][0]]+w*corners[edges[e][1]]) ;
        real d2 = dot(dv,dv) ;
        if(d2 < mindist2) {
          mindist2 = d2 ;
          minvec = dv ;
        }
      }
      // minimum distance to all edges
      real mindist = sqrt(mindist2) ;
	
      // Now check for face projections
      for(int f=0;f<12;++f) {
        // maps onto all edges, so project onto triangle
        const vector3d<real_t> t0 = corners[triangles[f][0]] ;
        const vector3d<real_t> t1 = corners[triangles[f][1]]-t0 ;
        const vector3d<real_t> t2 = corners[triangles[f][2]]-t0 ;
        real distf = dot(pt-t0,fn[f]) ;
        vector3d<real_t> pp = (pt-distf*fn[f])-t0 ; // point projected on face
        if((dot(fn[f],cross(t1,pp))>= 0 &&
            dot(fn[f],cross(pp,t2))>= 0 &&
            dot(fn[f],cross(t2-t1,pp-t1)) >= 0)
           ) {
          if(fabs(distf) < fabs(mindist)*(1.00001) ) {
            mindist = distf ;
            minvec = fn[f] ;
          }
        }
      }
      dist[i] = mindist-radius ;
      n[i] = (1./max<real>(norm(minvec),1e-30))*minvec ;
    }
  }

  Loci::CPTR<distFunc> hexFunc::transform(rigid_transform Tv) const {
    vector3d<real_t> ncorners[8] ;
    for(int i=0;i<8;++i)
      ncorners[i] = Tv.transform(corners[i]) ;
    return new hexFunc(ncorners,radius) ;
  }


  /// Define geometry of a cylinder given a radius and two points
  class cylinderFunc : public distFunc {
    real radius ;
    vector3d<real_t> p1,p2 ;
  public :
    cylinderFunc(vector3d<real_t> p1i,vector3d<real_t> p2i, real r): radius(r),p1(p1i),p2(p2i) {}
    cylinderFunc(vector3d<real_t> c, real r) {
      radius = r ;
      p1 = vector3d<real_t>(c.x,c.y,-1e30) ;
      p2 = vector3d<real_t>(c.x,c.y,1e30) ;
    }
    virtual void getDist(real dist[],vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
    virtual Loci::CPTR<distFunc> transform(rigid_transform Tv)  const ;
  } ;

  void cylinderFunc::getDist(real dist[], vector3d<real_t> n[],
                             const vector3d<real_t> loc[], int sz) const {
    vector3d<real_t> dv = (p2-p1) ;
    real ndv = norm(dv) ;
    real rdv = 1./max<real_t>(ndv,1e-30) ;
    dv *= rdv ;
    for(int i=0;i<sz;++i) {
      real w = dot(loc[i]-p1,dv)*rdv ;
      vector3d<real_t> paxis = (1.-w)*p1+w*p2 ;
      vector3d<real_t> pr = loc[i]-paxis ;
      real npr = norm(pr) ;
      dist[i] = npr-radius ;
      pr *= 1./max<real_t>(npr,1e-30) ;
      n[i] = pr ;
      if(dist[i]<0.0) { // Inside  so check distance to endcaps
        if(w<= 0) { // outside endcap1
          dist[i] = -w*ndv ;
          n[i] = -1.*dv ;
        } else if(w>=1.0) { // outside endcap2
          dist[i] = (w-1.0)*ndv ;
          n[i] = dv ;
        } else { // inside cylinder
          if(w<0.5) {
            if(dist[i] < -w*ndv) {
              dist[i] = -w*ndv ;
              n[i] = -1.*dv ;
            }
          } else {
            if(dist[i] < (w-1.0)*ndv) {
              dist[i] = (w-1.0)*ndv ;
              n[i] = dv ;
            }
          }
        }
      } else if(w< 0.0) { // end cap
        vector3d<real_t> pedge = p1+n[i]*radius ;
        vector3d<real_t> dvedge = loc[i]-pedge ;
        dist[i] = norm(dvedge) ;
        n[i] = dvedge/max<real_t>(dist[i],1e-30) ;
      } else if(w > 1.0) { // end cap
        vector3d<real_t> pedge = p2+n[i]*radius ;
        vector3d<real_t> dvedge = loc[i]-pedge ;
        dist[i] = norm(dvedge) ;
        n[i] = dvedge/max<real_t>(dist[i],1e-30) ;
      }
    }
  }

  Loci::CPTR<distFunc> cylinderFunc::transform(rigid_transform Tv)  const {
    vector3d<real_t> np1 = Tv.transform(p1) ;
    vector3d<real_t> np2 = Tv.transform(p2) ;

    return new cylinderFunc(np1,np2,radius) ;
  }


  /// Define spherical geometry from center point and radius
  class sphereFunc : public distFunc {
    vector3d<real_t> center;
    real radius ;
  public :
    sphereFunc(vector3d<real_t> c, real r): center(c),radius(r) {}
    virtual void getDist(real dist[], vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
    virtual Loci::CPTR<distFunc> transform(rigid_transform Tv) const ;
  } ;

  void sphereFunc::getDist(real dist[], vector3d<real_t> n[],
                           const vector3d<real_t> loc[], int sz) const {
    for(int i=0;i<sz;++i) {
      vector3d<real_t> nl = loc[i]-center ;
      real nn = norm(nl) ;
      n[i] = (1./max<real>(nn,1e-30))*nl ;
      dist[i]=nn-radius ;
      dist[i]=norm(loc[i]-center)-radius ;
    }
  }

  Loci::CPTR<distFunc> sphereFunc::transform(rigid_transform Tv)  const {
    vector3d<real_t> ct = Tv.transform(center) ;

    return new sphereFunc(ct,radius) ;
  }

  /// Compute coordinate aligned box given two opposite corner points and
  /// an extrusion radii
  class boxFunc : public distFunc {
    vector3d<real_t> pmin,pmax;
    real radius ;
  public :
    boxFunc(vector3d<real_t> c1, vector3d<real_t> c2, real r): radius(r) {
      pmin.x = min(c1.x,c2.x) ;
      pmax.x = max(c1.x,c2.x) ;
      pmin.y = min(c1.y,c2.y) ;
      pmax.y = max(c1.y,c2.y) ;
      pmin.z = min(c1.z,c2.z) ;
      pmax.z = max(c1.z,c2.z) ;
    }
    virtual void getDist(real dist[], vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
    virtual Loci::CPTR<distFunc> transform(rigid_transform Tv) const ;
  } ;

  void boxFunc::getDist(real dist[], vector3d<real_t> n[],
                        const vector3d<real_t> loc[], int sz) const {
    for(int i=0;i<sz;++i) {
      vector3d<real_t> pt = loc[i] ;
      vector3d<real_t> projected = pt ;
      // find projected point
      projected.x = min(pmax.x,max(pmin.x,pt.x)) ;
      projected.y = min(pmax.y,max(pmin.y,pt.y)) ;
      projected.z = min(pmax.z,max(pmin.z,pt.z)) ;
      vector3d<real_t> dv = pt-projected ;
      real dv2 = dot(dv,dv) ;
      if(dv2 == 0) { // inside
        real dst = pt.x-pmin.x ;
        n[i] = vector3d<real_t>(-1.,0,0) ;
        if(dst > pmax.x-pt.x) {
          dst = pmax.x-pt.x ;
          n[i] = vector3d<real_t>(1.,0.,0.) ;
        }
        if(dst > pt.y-pmin.y) {
          dst = pt.y-pmin.y ;
          n[i] = vector3d<real_t>(0.,-1.,0.) ;
        }
        if(dst > pmax.y-pt.y) {
          dst = pmax.y-pt.y ;
          n[i] = vector3d<real_t>(0.,1.,0.) ;
        }
        if(dst > pt.z-pmin.z) {
          dst = pt.z-pmin.z ;
          n[i] = vector3d<real_t>(0.,0.,-1.) ;
        }
        if(dst > pmax.z-pt.z) {
          dst = pmax.z-pt.z ;
          n[i] = vector3d<real_t>(0.,0.,1.) ;
        }
        dist[i] = -dst-radius ;
      } else {
        dist[i] = sqrt(dv2) ;
        n[i] = (1./dist[i])*dv ;
        dist[i] -= radius ;
      }
    }
  }

  Loci::CPTR<distFunc> boxFunc::transform(rigid_transform Tv) const {
    vector3d<real_t> ncorners[8] ;
    ncorners[0] = Tv.transform(vector3d<real_t>(pmin.x,pmin.y,pmax.z)) ;
    ncorners[1] = Tv.transform(vector3d<real_t>(pmax.x,pmin.y,pmax.z)) ;
    ncorners[2] = Tv.transform(vector3d<real_t>(pmin.x,pmax.y,pmax.z)) ;
    ncorners[3] = Tv.transform(vector3d<real_t>(pmax.x,pmax.y,pmax.z)) ;
    ncorners[4] = Tv.transform(vector3d<real_t>(pmin.x,pmin.y,pmin.z)) ;
    ncorners[5] = Tv.transform(vector3d<real_t>(pmax.x,pmin.y,pmin.z)) ;
    ncorners[6] = Tv.transform(vector3d<real_t>(pmin.x,pmax.y,pmin.z)) ;
    ncorners[7] = Tv.transform(vector3d<real_t>(pmax.x,pmax.y,pmin.z)) ;
    return new hexFunc(ncorners,radius) ;
  }

  class conicFunc : public distFunc
{
  real_t   r1;  //!< [m]
  real_t   r2;  //!< [m]
  vector3d<real_t> p1;  //!< [m]
  vector3d<real_t> p2;  //!< [m]
public:

  conicFunc(real_t r1i, real_t r2i, vector3d<real_t> p1i, vector3d<real_t> p2i)
  {
    r1   = r1i;
    r2   = r2i;
    p1   = p1i;
    p2   = p2i;
  } // End of inconeFunc()

    virtual void getDist(real dist[], vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
    virtual Loci::CPTR<distFunc> transform(rigid_transform Tv) const ;
}

  void conicFunc::getDist(real dist[], vector3d<real_t> n[],
                        const vector3d<real_t> loc[], int sz) const {
    vector3d<real_t> dv = (p2-p1) ;
    real ndv = norm(dv) ;
    real rdv = 1./max<real_t>(ndv,1e-30) ;
    dv *= rdv ;
    for(int i=0;i<sz;++i) {
      real w = dot(loc[i]-p1,dv)*rdv ;
      vector3d<real_t> paxis = (1.-w)*p1+w*p2 ;
      vector3d<real_t> pr = loc[i]-paxis ;
      real npr = norm(pr) ;
      dist[i] = npr-r1 ; // 
      pr *= 1./max<real_t>(npr,1e-30) ;
      n[i] = pr ;
      if(dist[i]<0.0) { // Inside  so check distance to endcaps
        if(w<= 0) { // outside endcap1
          dist[i] = -w*ndv ;
          n[i] = -1.*dv ;
        } else if(w>=1.0) { // outside endcap2
          dist[i] = (w-1.0)*ndv ;
          n[i] = dv ;
        } else { // inside cylinder
          if(w<0.5) {
            if(dist[i] < -w*ndv) {
              dist[i] = -w*ndv ;
              n[i] = -1.*dv ;
            }
          } else {
            if(dist[i] < (w-1.0)*ndv) {
              dist[i] = (w-1.0)*ndv ;
              n[i] = dv ;
            }
          }
        }
      } else if(w< 0.0) { // end cap
        vector3d<real_t> pedge = p1+n[i]*(r1-r2) ;
        vector3d<real_t> dvedge = loc[i]-pedge ;
        dist[i] = norm(dvedge) ;
        n[i] = dvedge/max<real_t>(dist[i],1e-30) ;
      } else if(w > 1.0) { // end cap
        vector3d<real_t> pedge = p2+n[i]*(r1-r2) ;
        vector3d<real_t> dvedge = loc[i]-pedge ;
        dist[i] = norm(dvedge) ;
        n[i] = dvedge/max<real_t>(dist[i],1e-30) ;
      }
    }

}