//#############################################################################
//#
//# Copyright 2015-2026, Mississippi State University
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
#ifndef GEOMPRIMITIVES_H
#define GEOMPRIMITIVES_H

#include <Loci.h>
#include <Tools/tools.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
using std::string ;
using std::endl ;
using std::cout ;
using std::cerr ;
using std::ifstream ;
using std::ios ;
#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <set>
using std::set ;

namespace Loci
{

 /// Implicit Function Base Class that defines geometry API that computes
  /// distance function from a composition of fundamental geometric elements
  class geometry_type : public Loci::CPTR_type {
  public:
    /// Pure virtual function for obtaining distance and gradient of distance
    /// field evaluated at a provided location.  Defined such that it is
    /// evaluated for an array of locations of size sz.
    virtual void getDist(real_t dist[], vector3d<real_t> n[],
                         const vector3d<real_t> loc[],int sz) const = 0 ;
    /// Transform geometry using rigid body transform
    virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const = 0 ;
    virtual void inGeometry(real_t dist[],bool inGeom[],int sz) const = 0;
  } ;


  /// The composite function builds a geometry from a set of geometric
  /// components.  You are outside of the geometry ony if you are outside
  /// of all geometric components (e.g. compose as the union of all
  /// geometric elements).  Distance outside of the component is the minimum
  /// distance out of all of the geometric elements.
  class compositeFunc: public geometry_type {
    /// List of composed geometries
    std::vector<Loci::CPTR<geometry_type> > geomList ;
  public:
    compositeFunc(std::vector<Loci::CPTR<geometry_type> > &list) {
      geomList = list ;
      if(geomList.size() == 0)
        cerr << "warning composite function created with no geometry!"
             << endl ;
    }
    virtual void getDist(real_t dist[], vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
  
    virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const  ;
    virtual void inGeometry(real_t dist[],bool inGeom[],int sz) const;
  } ;


  /// Compute distance to composite geometry
  void compositeFunc::getDist(real_t dist[], vector3d<real_t> n[],
                              const vector3d<real_t> loc[], int sz) const  {
    vector<real_t> dists(sz) ;
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

  void compositeFunc::inGeometry(real_t dist[],bool inGeom [], int sz) const
  {
    for (int i = 0; i< sz; i++)
    {
        if (dist[i] < 0) inGeom[i] = true;
        if (dist[i] > 0) inGeom[i] = false;
    }
  }

  /// Apply geometric transform to each element
  Loci::CPTR<geometry_type> compositeFunc::transform(rigid_transform Tv) const {
    std::vector<Loci::CPTR<geometry_type> > newGeomList(geomList.size()) ;
    for(size_t i=0;i<geomList.size();++i)
      newGeomList[i] = geomList[i]->transform(Tv) ;
    return new compositeFunc(newGeomList) ;
  }


class generalCylinderFunc : public geometry_type {
  protected:
    real_t r1, r2;
    vector3d<real_t> p1, p2;
  public: 
    generalCylinderFunc(vector3d<real_t> p1i, vector3d<real_t> p2i, real_t r1i, real_t r2i): 
      r1(r1i),r2(r2i),p1(p1i),p2(p2i) {}
    
    generalCylinderFunc(vector3d<real_t> c, real_t r1i, real_t r2i) {
      r1 = r1i ;
      r2 = r2i ;
      p1 = vector3d<real_t>(c.x,c.y,-1e30) ;
      p2 = vector3d<real_t>(c.x,c.y,1e30) ;
    }
    virtual void getDist(real_t dist[], vector3d<real_t> n[], const vector3d<real_t> loc[], int sz) const;
    virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const;
    virtual void inGeometry(real_t dist[],bool inGeom[],int sz) const;
};

void generalCylinderFunc::getDist(real_t dist[], vector3d<real_t> n[], 
    const vector3d<real_t> loc[], int sz) const
{
    // 1. Calculate the axis geometry
  vector3d<real_t> v = p2-p1;
  real_t H = norm(v);
  real_t oH = 1.0/max<real_t>(1e-30,H);
  v *=oH;
  for (int i=0;i<sz;i++) {
    
    // local vector from point to cap
    vector3d<real_t> w = loc[i] - p1;
    
    // Project p into the 2D local space (r,z)
    real_t axialProj = dot(w, v); 
    vector3d<real_t> radialVec = w - axialProj * v;
    real_t radialDist = norm(radialVec);
    
    // Shift origin to the midpoint of the segment to match the previous local bounds [-h, +h]
    real_t h =  0.5*H;
    real_t r = radialDist;
    real_t z = axialProj - h; 
    
    // 3. Evaluate the 2D profile (Inigo Quilez logic)
    real_t dx = r1 - r2;
    real_t mx = r1 + r2;
    
    real_t k = h * h + dx * dx;
    real_t s = h * r + dx * z - h * mx;
    real_t t = dx * r - h * z + h * (r1 - r2);
    
    real_t invSqrtK = 1.0 / max<real_t>(1e-30,sqrt(k));
    real_t dSide = s * invSqrtK;
    real_t dCap = abs(z) - h;
    
    real_t n2dx = 0.0;
    real_t n2dy = 0.0;
    dist[i] = 0.0;
    
    // Operate in 2D projection space
    if (t < 0.0) 
    {
        dist[i] = max<real_t>(r - r1, dCap);
        if (r - r1 > dCap) { n2dx = 1.0; n2dy = 0.0; } 
        else { n2dx = 0.0; n2dy = (z > 0.0) ? -1.0 : 1.0; }
    } 
    else if (t > 2.0 * k) 
    {
        dist[i] = max<real_t>(r - r2, dCap);
        if (r - r2 > dCap) { n2dx = 1.0; n2dy = 0.0; } 
        else { n2dx = 0.0; n2dy = (z > 0.0) ? -1.0 : 1.0; }
    } 
    else 
    {
        dist[i] = max<real_t>(dSide, dCap);
        if (dSide > dCap) { n2dx = h * invSqrtK; n2dy = dx * invSqrtK; } 
        else { n2dx = 0.0; n2dy = (z > 0.0) ? 1.0 : -1.0; }
    }
    
    // 4. Reconstruct the 3D Normal back to global Space
    // The radial component maps to the normalized radial vector direction
    vector3d<real_t> radialDir = radialVec/max<real_t>(1e-30,radialDist);

    // Combine the transformed 2D normal back into global coordinates
    n[i] = n2dx * radialDir + n2dy * v;
    // normalize 
    real_t oNn = 1.0/max<real_t>(1e-30,norm(n[i]));
    n[i] *=oNn;
  }
}

  Loci::CPTR<geometry_type> generalCylinderFunc::transform(rigid_transform Tv) const {
    vector3d<real_t> np1 = Tv.transform(p1) ;
    vector3d<real_t> np2 = Tv.transform(p2) ;
    
    return new generalCylinderFunc(np1,np2,r1,r2) ;
  }

//   class coneFunc : public generalCylinderFunc {
//     real_t r1, r2;
//     vector3d<real_t> p1, p2;
//   public: 
//     coneFunc(vector3d<real_t> p1i, vector3d<real_t> p2i, real_t ri): 
//       r1(ri),r2(0.0),p1(p1i),p2(p2i) {}
    
//     coneFunc(vector3d<real_t> c, real_t rb) {
//       r1 = rb ;
//       r2 = 0.0;
//       p1 = vector3d<real_t>(c.x,c.y,-1e30) ;
//       p2 = vector3d<real_t>(c.x,c.y,1e30) ;
//     }
//     virtual void getDist(real_t dist[], vector3d<real_t> n[], const vector3d<real_t> loc[], int sz) const;
//     virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const;
// };

// void coneFunc::getDist(real_t dist[], vector3d<real_t> n[], 
//     const vector3d<real_t> loc[], int sz) const
// {
//   // use generalCylinderFunc to compute distance
//   generalCylinderFunc::getDist(dist, n, loc, sz);
// }

//   Loci::CPTR<geometry_type> coneFunc::transform(rigid_transform Tv)  const {
//     vector3d<real_t> np1 = Tv.transform(p1) ;
//     vector3d<real_t> np2 = Tv.transform(p2) ;
    
//     return new coneFunc(np1,np2,r) ;
//   }

class hollowGeneralCylinderFunc : public generalCylinderFunc {
    real_t r1in, r2in;
    real_t thickness;
  public:
    hollowGeneralCylinderFunc(vector3d<real_t> p1i, vector3d<real_t> p2i, real_t r1o, real_t r2o, real_t r1i, real_t r2i, real_t thick): 
      generalCylinderFunc(p1i,p2i,r1o,r2o),r1in(r1i),r2in(r2i), thickness(thick) {}

    hollowGeneralCylinderFunc(vector3d<real_t> c, real_t r1o, real_t r2o, real_t thick)
      : generalCylinderFunc(c,r1o,r2o), r1in(thick*r1o),r2in(thick*r2o),thickness(thick) {}


    virtual void getDist(real_t dist[], vector3d<real_t> n[], const vector3d<real_t> loc[], int sz) const;
    virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const;
    virtual void inGeometry(real_t dist[],bool inGeom[],int sz) const;
};

// For the hollow cylinder, call twice to generalCylinder 
void hollowGeneralCylinderFunc::getDist(real_t dist[], vector3d<real_t> n[], const vector3d<real_t> loc[], int sz) const
{
  // Interior core
  real_t disti[sz]; 
  vector3d<real_t> ni[sz];
  generalCylinderFunc::getDist(disti, ni, loc, sz);

  // Exterior shell
  real_t disto[sz];
  vector3d<real_t> no[sz];
  generalCylinderFunc::getDist(disto,no,loc,sz);

  // Compare -inner to outer
  for (int i=0;i<sz;i++)
  {
    dist[i] = max<real_t>(-disti[i],disto[i]);
    if (disto[i] > disti[i])
    {
      n[i] = no[i];
    }
    else
    {
      n[i] = -1.0*ni[i];
    }
  }
}

  Loci::CPTR<geometry_type> hollowGeneralCylinderFunc::transform(rigid_transform Tv)  const {
    vector3d<real_t> np1 = Tv.transform(p1) ;
    vector3d<real_t> np2 = Tv.transform(p2) ;
    
    return new hollowGeneralCylinderFunc(np1,np2,r1,r2,r1in,r2in,thickness) ;
  }

  // Internal 2D Polygon SDF with analytical normal calculation
real_t sdPolygon2D(vector3d<real_t> q, const std::vector<vector3d<real_t>>& v, vector3d<real_t>& n2d) 
{
    int num = v.size();
    real_t s = 1.0;
    vector3d<real_t> closest = vector3d<real_t>(1.0,0.0,0.0);
    real_t minDistSq = 1e20;

    for (int i = 0, j = num - 1; i < num; j = i++) 
    {
        vector3d<real_t> e = v[j] - v[i];
        vector3d<real_t> w = q - v[i];

        real_t edgeLenSq = dot(e,e);
        real_t oe = 1.0/max<real_t>(1e-30,edgeLenSq);
        real_t h = min<real_t>(1.0,max<real_t>(0.0,dot(w,e)*oe));
        
        vector3d<real_t> d = w - e*h;
        real_t distSq = dot(d,d);
        
        if (distSq < minDistSq) {
            minDistSq = distSq;
            real_t len = sqrt(distSq);
            real_t olen = 1.0/max<real_t>(1e-30,len);
            closest = d*olen;
        }

        // Winding number interior test
        bool cond1 = q.y >= v[i].y;
        bool cond2 = q.y < v[j].y;
        bool cond3 = e.x * w.y > e.y * w.x;
        if ((cond1 && cond2 && cond3) || (!cond1 && !cond2 && !cond3)) {
            s = -s;
        }
    }
    
    n2d = closest;
    return s * sqrt(minDistSq);
}

  class sectorFunc : public geometry_type {
    vector3d<real_t> p1;
    vector3d<real_t> p2;
    vector<vector3d<real_t>> polyVerts;
    real_t angle;
    real_t sectorPosition;
    vector3d<real_t> sectorOrientation;
  public:
    sectorFunc(vector3d<real_t> p1i, vector3d<real_t> p2i, vector<vector3d<real_t>> verts, real_t alpha, real_t per, vector3d<real_t> orient): 
      p1(p1i),p2(p2i),polyVerts(verts),angle(alpha),sectorPosition(per),sectorOrientation(orient) {}
    sectorFunc(vector3d<real_t> c){
        angle = 2.0*M_PI;
        sectorPosition = 0.0;
        sectorOrientation = cross(vector3d<real_t>(c.x,0.0,0.0),vector3d<real_t>(0.0,c.y,0.0));
        polyVerts.push_back(vector3d<real_t>(0.0,0.0,0.0));
        polyVerts.push_back(vector3d<real_t>(c.x,0.0,0.0));
        polyVerts.push_back(vector3d<real_t>(c.x,c.y,0.0));
        polyVerts.push_back(vector3d<real_t>(0.0,c.y,0.0));
    }
    virtual void getDist(real_t dist[], vector3d<real_t> n[], const vector3d<real_t> loc[], int sz) const;
    virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const;
    virtual void inGeometry(real_t dist[],bool inGeom[],int sz) const;
  };

void sectorFunc::getDist(real_t dist[], vector3d<real_t> n[], const vector3d<real_t> loc[], int sz) const
{
    // 1. Calculate the axis geometry
    vector3d<real_t> V = p2 - p1;
    real_t H = norm(V);
    real_t oH = 1.0/max<real_t>(1e-30,H);
    //if (H < 1e-9 || polyVerts.empty()) return 1e10;
    
    V *= oH;
    for (int i = 0; i< sz; i++)
    {
      vector3d<real_t> w = loc[i]-p1;
      
      // Project test point p into 2D cylindrical space
      // CHECK IMPLEMENTATION OF TEST POINT
      vector3d<real_t> q;
      q.y = dot(w, V); // Height along axis
      vector3d<real_t> radialVec = w - q.y*V;
      q.x = norm(radialVec); // Radius away from axis
      
      // 2. Transform the 3D Polygon Vertices into the exact same 2D profile space
      std::vector<vector3d<real_t>> poly2d;
      poly2d.reserve(polyVerts.size());
      
      for (const auto& v3d : polyVerts) {
          vector3d<real_t> vW =  v3d - p1;
          real_t vY = dot(vW, V); // Vertex height along axis
          vector3d<real_t> vRadial = vW - vY * V;
          real_t vX = norm(vRadial); // Vertex radius away from axis
          
          poly2d.push_back({vX,vY,0.0});
      }
      
      // 3. Establish local angular coordinate frame orthogonal to axis vector 'u'
     // vector3d<real_t> referenceDir = (abs(V.y) < 0.9) ? vector3d<real_t>{0.0, 1.0, 0.0} : vector3d<real_t>{1.0, 0.0, 0.0};
      real_t proj = dot(sectorOrientation,V);
      vector3d<real_t> r = sectorOrientation - proj*V;

      // cylindrical coordinate transform
      real_t rdist = norm(r);
      real_t ord = 1.0/max<real_t>(1e-30,rdist);
      r *= ord;
      vector3d<real_t> z = cross(V,r);

      vector3d<real_t> proj2d = vector3d(dot(radialVec,r),dot(radialVec,z),0.0);
      
      // 4. Evaluate the transformed 2D Polygon Profile
      vector3d<real_t> n2d;
      real_t dProfile = sdPolygon2D(q, poly2d, n2d);
      dist[i] = 0.0;
        
      // 5. Handle Angular Sector Bounds
      if (angle >= 2.0 * M_PI - 1e-5) 
      {
          // Full 360-degree revolution
          vector3d<real_t> radialDir = radialVec/max<real_t>(1e-30,q.x);
          n[i] = n2d.x * radialDir + n2d.y * V;
          real_t nn = norm(n[i]);
          real_t oNN = 1.0/max<real_t>(1e-30,nn);
          n[i] *= oNN;

          dist[i] =  dProfile;
          return;
      }
      
      // Partial angular wedge processing
      real_t halfAngle = angle * 0.5;
      real_t sinA = sin(halfAngle);
      real_t cosA = cos(halfAngle);

      // Offset angle relative to orientation and percentage
      real_t offset = (0.5 - sectorPosition) * angle;

      // Angular offset for plane
      real_t cosO = cos(offset);
      real_t sinO = sin(offset);
      vector3d<real_t> r2d; 
      r2d.x = proj2d.x * cosO - proj2d.y * sinO;
      r2d.z = proj2d.x * sinO + proj2d.y * cosO;
      r2d.y = 0.0;
      real_t absZ = abs(r2d.z);
      
      real_t dSector = (r2d.x * sinA - absZ * cosA);
      bool isInsideSector = (r2d.x * cosA + absZ * sinA > 0.0) && (dSector < 0.0);
      
      if (isInsideSector) {
          vector3d<real_t> radialDir = radialVec/max<real_t>(1e-30,q.x);
          n[i] = n2d.x * radialDir + n2d.y * V;
          real_t nn = norm(n[i]);
          real_t oNN = 1.0/max<real_t>(1e-30,nn);
          n[i] *= oNN;
          dist[i] = dProfile;
          return;
      } 
      else 
      {
          // Point is outside the slice; compute distance to the slicing edge planes
          dist[i] = max<real_t>(dProfile, dSector);
          
          if (dProfile > dSector) {
              vector3d<real_t> radialDir = radialVec/max<real_t>(1e-30,q.x);
              n[i] = n2d.x * radialDir + n2d.y * V;
              real_t nn = norm(n[i]);
              real_t oNN = 1.0/max<real_t>(1e-30,nn);
              n[i] *= oNN;
          } else {
              // Normal matches the flat sliced plane faces of the sector wedge
              real_t zSign = (r2d.z >= 0.0) ? 1.0 : -1.0;
              vector3d<real_t> nrot;
              nrot.x = sinA;
              nrot.y = -cosA * zSign;
              nrot.z = 0.0;
              
              // back into the original plane
              n2d.x = nrot.x*cosO + nrot.y*sinO;
              n2d.y = -nrot.x*sinO + nrot.y*cosO;
              
              n[i] = n2d.x * r + n2d.y * z;
              real_t nn = norm(n[i]);
              real_t oNN = 1.0/max<real_t>(1e-30,nn);
              n[i] *= oNN;
          }
      }
    }
  }

  Loci::CPTR<geometry_type> sectorFunc::transform(rigid_transform Tv)  const {
    vector3d<real_t> np1 = Tv.transform(p1) ;
    vector3d<real_t> np2 = Tv.transform(p2) ;
    vector3d<real_t> nori = Tv.transform(sectorOrientation);
    vector<vector3d<real_t>> pv;
    pv.reserve(polyVerts.size());
    for (size_t i = 0; i< polyVerts.size(); i++)
    {
        vector3d<real_t> pl = Tv.transform(polyVerts[i]);
        pv.push_back(pl);
    }
    return new sectorFunc(np1, np2, pv, angle, sectorPosition, nori);
  }

class extrusionFunc : public geometry_type {
    vector3d<real_t> p1; 
    vector3d<real_t> p2;
    vector<vector3d<real_t>> polyVerts; /// verticies for base polygon
    real_t extrusionPercentage;
  public:
    extrusionFunc(vector3d<real_t> p1i, vector3d<real_t> p2i, vector<vector3d<real_t>> verts, real_t per):
        p1(p1i),p2(p2i),polyVerts(verts),extrusionPercentage(per){}
    extrusionFunc(vector3d<real_t> c)
    {
        polyVerts.reserve(4);
        polyVerts.push_back(vector3d<real_t>(0.0,0.0,0.0));
        polyVerts.push_back(vector3d<real_t>(c.x,0.0,0.0));
        polyVerts.push_back(vector3d<real_t>(c.x,c.y,0.0));
        polyVerts.push_back(vector3d<real_t>(0.0,c.y,0.0));
        p1 = vector3d<real_t>(0.0,0.0,0.0);
        p2 = cross(p1,polyVerts[3]);
        extrusionPercentage = 0.0;
    }
    virtual void getDist(real_t dist[], vector3d<real_t> n[], const vector3d<real_t> loc[], int sz) const;
    virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const;
    virtual void inGeometry(real_t dist[],bool inGeom[],int sz) const;    
  };

void extrusionFunc::getDist(real_t dist[], vector3d<real_t> n[], const vector3d<real_t> loc[], int sz) const
{
     // 1. Calculate the main axis geometry (Y-axis of our local frame)
    vector3d<real_t> V = p2 - p1;
    real_t H = norm(V);
    real_t oH = 1.0/max<real_t>(1e-30,H);
    //if (H < 1e-9 || polyVerts.empty()) return 1e10;
    V *= oH;

    for (int i = 0;i<sz;i++)
    {
      vector3d<real_t> w = loc[i] - polyVerts[0];
      vector3d<real_t> wC = loc[i] - p1;
      real_t wClen = dot(wC,V);
      
      // 2. Establish the rest of the 3D coordinate frame
      // We determine the plane normal from the first few vertices of the 3D polygon
      vector3d<real_t> x = polyVerts[1] - polyVerts[0];
      real_t elen = norm(x);
      real_t olene = 1.0/max<real_t>(1e-30,elen);
      x *= olene;


      vector3d<real_t> edge2 = polyVerts[3] - polyVerts[0];
      elen = norm(edge2);
      olene = 1.0/max<real_t>(1e-30,elen);
      edge2 *=olene;
      vector3d<real_t> polyNorm = cross(x,edge2);
      real_t pnlen = norm(polyNorm);
      real_t opnlen = 1.0/max<real_t>(1e-30,pnlen);
      polyNorm *= opnlen;
      
      // // Cross product of edges gives the raw normal of the polygon's plane (Local Z - Extrusion direction)
      vector3d<real_t> y = cross(polyNorm,x);
      real_t leny = norm(y);
      real_t oleny = 1.0/max<real_t>(1e-30,leny);
      y *= oleny;

      // 3. Transform the 3D Polygon Vertices into this flat 2D frame (X=Lateral, Y=Axial)
      std::vector<vector3d<real_t>> poly2d;
      poly2d.reserve(polyVerts.size());
      for (const auto& v3d : polyVerts) {
          vector3d<real_t> vW = v3d - polyVerts[0];
          real_t vX = dot(vW, x);
          real_t vY = dot(vW, y);
          poly2d.push_back({ vX, vY,0.0 });
      }
      
      // 4. Map the test point into the local 3D coordinate values
      vector3d<real_t> p2d = vector3d<real_t>(dot(w,x),dot(w,y),0.0);
      // 5. Evaluate the 2D Polygon profile distance
      vector3d<real_t> n2d;
      real_t dProfile2D = sdPolygon2D(p2d, poly2d, n2d);

      // 6. Evaluate Extrusion Bounds using the user's percentage bias
      // Shift the extrusion center window based on where they anchored the polygon plane
      real_t target = extrusionPercentage * H;
      real_t shift = wClen + target;

      // get distance along the extrusion path
      real_t begin = -shift;
      real_t end = shift - H;
      real_t dExtrusion = max<real_t>(begin,end);
      
      // 7. Combine calculations using the standard bounded 3D Extrusion operator
      // If inside the profile, it finds the closest exit out of the caps or side walls.
      // If both distances are negative, the point is deep inside; check which feature is closer
      if (dProfile2D < 0.0 && dExtrusion < 0.0) {
          dist[i]  = max<real_t>(dProfile2D, dExtrusion);
      } else {
          // Outside regions take the Euclidean combination of the positive distances
          real_t odist = 0.0;
          if (dProfile2D > 0.0) odist += dProfile2D * dProfile2D;
          if (dExtrusion > 0.0) odist += dExtrusion * dExtrusion;
          dist[i] = sqrt(odist);
      }
      
      // 8. Reconstruct the World Normal Vector
      if (dist[i] == dProfile2D || (dProfile2D > dExtrusion && dist[i] <= 0.0)) {
          // Closest feature is the perimeter side walls of the polygon profile
          n[i] = n2d.x*x + n2d.y*y;
      } else if (dist[i] == dExtrusion || dist[i] <= 0.0) {
          // Closest feature is the front or back flat end cap of the extrusion length
          real_t capSign = (shift >= 0.0) ? 1.0 : -1.0;
          n[i] = V*capSign;
      }
      
      // normalization
      real_t nlen = norm(n[i]);
      real_t olen = 1.0/max<real_t>(1e-30,nlen);
      n[i] *= olen;    
  }
}

  Loci::CPTR<geometry_type> extrusionFunc::transform(rigid_transform Tv)  const {
    vector3d<real_t> np1 = Tv.transform(p1) ;
    vector3d<real_t> np2 = Tv.transform(p2) ;

    vector<vector3d<real_t>> pv;
    pv.reserve(polyVerts.size());
    for (size_t i = 0; i< polyVerts.size(); i++)
    {
        vector3d<real_t> pl = Tv.transform(polyVerts[i]);
        pv.push_back(pl);
    }
    return new extrusionFunc(np1, np2, pv, extrusionPercentage);
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
  class hexFunc: public geometry_type {
    vector3d<real_t> corners[8] ;
    vector3d<real_t> fn[12] ;
    vector3d<real_t> ev[18] ;
    real_t evl[18] ;
  public:
    hexFunc(vector3d<real_t> pts[8]) {
      using namespace hexFuncInfo ;
      for(int i=0;i<8;++i)
        corners[i] = pts[i] ;
      for(int f=0;f<12;++f) {
        vector3d<real_t> dv1 = corners[triangles[f][1]]-corners[triangles[f][0]] ;
        vector3d<real_t> dv2 = corners[triangles[f][2]]-corners[triangles[f][0]] ;
        vector3d<real_t> A = cross(dv1,dv2) ;
        fn[f] = A/max<real_t>(norm(A),1e-30) ;
      }
      for(int e=0;e<18;++e) {
        vector3d<real_t> dv = corners[edges[e][1]]-corners[edges[e][0]] ;
        evl[e] = norm(dv) ;
        ev[e] = dv/max<real_t>(evl[e],1e-30) ;
      }
    }
    virtual void getDist(real_t dist[],vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
    virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const ;
    virtual void inGeometry(real_t dist[],bool inGeom[],int sz) const;
  } ;

  /// Get distance field by projecting onto points, then edges, then triangular
  /// faces of the decomposed hex
  void hexFunc::getDist(real_t dist[],vector3d<real_t> n[],
                        const vector3d<real_t> loc[], int sz) const {
    using namespace hexFuncInfo ;
    for(int i=0;i<sz;++i) { // loop over points
      vector3d<real_t> pt = loc[i] ;
      // first project onto all of the edges
      real_t mindist2 = 1e30 ;
      vector3d<real_t> minvec(0,0,0) ;
      for(int e=0;e<18;++e) {
        real_t w = dot(pt-corners[edges[e][0]],ev[e])/max<real_t>(evl[e],1e-15) ;
        w = max<real_t>(0.0,min<real_t>(1.0,w)) ;
        vector3d<real_t> dv = pt-((1.-w)*corners[edges[e][0]]+w*corners[edges[e][1]]) ;
        real_t d2 = dot(dv,dv) ;
        if(d2 < mindist2) {
          mindist2 = d2 ;
          minvec = dv ;
        }
      }
      // minimum distance to all edges
      real_t mindist = sqrt(mindist2) ;
	
      // Now check for face projections
      for(int f=0;f<12;++f) {
        // maps onto all edges, so project onto triangle
        const vector3d<real_t> t0 = corners[triangles[f][0]] ;
        const vector3d<real_t> t1 = corners[triangles[f][1]]-t0 ;
        const vector3d<real_t> t2 = corners[triangles[f][2]]-t0 ;
        real_t distf = dot(pt-t0,fn[f]) ;
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
      dist[i] = mindist;//-radius ;
      n[i] = (1./max<real_t>(norm(minvec),1e-30))*minvec ;
    }
  }

  Loci::CPTR<geometry_type> hexFunc::transform(rigid_transform Tv) const {
    vector3d<real_t> ncorners[8] ;
    for(int i=0;i<8;++i)
      ncorners[i] = Tv.transform(corners[i]) ;
    return new hexFunc(ncorners) ;
  }

 /// Define spherical geometry from center point and radius
  class sphereFunc : public geometry_type {
    vector3d<real_t> center;
    real_t radius ;
  public :
    sphereFunc(vector3d<real_t> c, real_t r): center(c),radius(r) {}
    virtual void getDist(real_t dist[], vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
    virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const ;
    virtual void inGeometry(real_t dist[],bool inGeom[],int sz) const;
  } ;

  void sphereFunc::getDist(real_t dist[], vector3d<real_t> n[],
                           const vector3d<real_t> loc[], int sz) const {
    for(int i=0;i<sz;++i) {
      vector3d<real_t> nl = loc[i]-center ;
      real_t nn = norm(nl) ;
      n[i] = (1./max<real_t>(nn,1e-30))*nl ;
      dist[i]=nn-radius ;
      dist[i]=norm(loc[i]-center)-radius ;
    }
  }

  Loci::CPTR<geometry_type> sphereFunc::transform(rigid_transform Tv)  const {
    vector3d<real_t> ct = Tv.transform(center) ;

    return new sphereFunc(ct,radius) ;
  }

  /// Compute coordinate aligned box given two opposite corner points and
  /// an extrusion radii
  class boxFunc : public geometry_type {
    vector3d<real_t> pmin,pmax;
  public :
    boxFunc(vector3d<real_t> c1, vector3d<real_t> c2) {
      pmin.x = min(c1.x,c2.x) ;
      pmax.x = max(c1.x,c2.x) ;
      pmin.y = min(c1.y,c2.y) ;
      pmax.y = max(c1.y,c2.y) ;
      pmin.z = min(c1.z,c2.z) ;
      pmax.z = max(c1.z,c2.z) ;
    }
    virtual void getDist(real_t dist[], vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
    virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const ;
    virtual void inGeometry(real_t dist[],bool inGeom[],int sz) const;
  } ;

  void boxFunc::getDist(real_t dist[], vector3d<real_t> n[],
                        const vector3d<real_t> loc[], int sz) const {
    for(int i=0;i<sz;++i) {
      vector3d<real_t> pt = loc[i] ;
      vector3d<real_t> projected = pt ;
      // find projected point
      projected.x = min(pmax.x,max(pmin.x,pt.x)) ;
      projected.y = min(pmax.y,max(pmin.y,pt.y)) ;
      projected.z = min(pmax.z,max(pmin.z,pt.z)) ;
      vector3d<real_t> dv = pt-projected ;
      real_t dv2 = dot(dv,dv) ;
      if(dv2 == 0) { // inside
        real_t dst = pt.x-pmin.x ;
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
        dist[i] = -dst; //-radius ;
      } else {
        dist[i] = sqrt(dv2) ;
        n[i] = (1./dist[i])*dv ;
 //       dist[i] -= radius ;
      }
    }
  }

  Loci::CPTR<geometry_type> boxFunc::transform(rigid_transform Tv) const {
    vector3d<real_t> ncorners[8] ;
    ncorners[0] = Tv.transform(vector3d<real_t>(pmin.x,pmin.y,pmax.z)) ;
    ncorners[1] = Tv.transform(vector3d<real_t>(pmax.x,pmin.y,pmax.z)) ;
    ncorners[2] = Tv.transform(vector3d<real_t>(pmin.x,pmax.y,pmax.z)) ;
    ncorners[3] = Tv.transform(vector3d<real_t>(pmax.x,pmax.y,pmax.z)) ;
    ncorners[4] = Tv.transform(vector3d<real_t>(pmin.x,pmin.y,pmin.z)) ;
    ncorners[5] = Tv.transform(vector3d<real_t>(pmax.x,pmin.y,pmin.z)) ;
    ncorners[6] = Tv.transform(vector3d<real_t>(pmin.x,pmax.y,pmin.z)) ;
    ncorners[7] = Tv.transform(vector3d<real_t>(pmax.x,pmax.y,pmin.z)) ;
    return new hexFunc(ncorners) ;
  }


   /// Planar list geometry type specification
  class planelistFunc : public geometry_type {
    std::vector<vector3d<real_t> > pts ;
    std::vector<vector3d<real_t> > normals ;
  public:
    planelistFunc(std::vector<vector3d<real_t> > &ps,
                   std::vector<vector3d<real_t> > &ns) : pts(ps),normals(ns) {}
    virtual void getDist(real_t dist[], vector3d<real_t> n[],
                         const vector3d<real_t> loc[], int sz) const ;
    virtual Loci::CPTR<geometry_type> transform(rigid_transform Tv) const ;
    virtual void inGeometry(real_t dist[],bool inGeom[],int sz) const;
  } ;

  void planelistFunc::getDist(real_t dist[], vector3d<real_t> n[],
                        const vector3d<real_t> loc[], int sz) const {
    real_t maxdist;
    for(int i=0;i<sz;++i) {
      for(size_t j=0;j<pts.size();++j) {
        real_t loc_dist = dot(loc[i]-pts[j],normals[j]);
        maxdist = max<real_t>(maxdist,loc_dist) ;
      }
      dist[i] = maxdist;
      n[i]  = normals[i];
    }
  }

  Loci::CPTR<geometry_type> planelistFunc::transform(rigid_transform Tv) const {
    vector<vector3d<real_t> > ptmp ;
    vector<vector3d<real_t> > ntmp ;
    for(size_t i=0;i<pts.size();++i) {
      vector3d<real_t> pt1 = pts[i] ;
      real_t np1 = norm(pt1)+1 ;
      vector3d<real_t> pt2 = pt1 + np1*normals[i] ;
      vector3d<real_t> p1 = Tv.transform(vector3d<real_t>(pt1.x,pt1.y,pt1.z)) ;
      vector3d<real_t> p2 = Tv.transform(vector3d<real_t>(pt2.x,pt2.y,pt2.z));
      ptmp.push_back(p1) ;
      vector3d<real_t> n1 = (p2-p1)/np1 ;
      ntmp.push_back(n1) ;
    }
    return new planelistFunc(ptmp,ntmp) ;
  }
}

#endif // GEOMPRIMITIVES