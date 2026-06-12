//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
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
#ifndef INCOMPOSITE_H
#define INCOMPOSITE_H
#include <Loci.h>


  /// Implicit Function Base Class that defines geometry API that computes
  /// distance function from a composition of fundamental geometric elements
  class distFunc : public Loci::CPTR_type {
  public:
    /// Pure virtual function for obtaining distance and gradient of distance
    /// field evaluated at a provided location.  Defined such that it is
    /// evaluated for an array of locations of size sz.
    virtual void getDist(real dist[], vect3d n[],
                         const vect3d loc[],int sz) const = 0 ;
    /// Transform geometry using rigid body transform
    virtual Loci::CPTR<distFunc> transform(rigid_transform Tv) const = 0 ;
    virtual real_t levelSet(vector3d<real_t> pt) const = 0; // +inside, -outside
  } ;


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

    /// Compute distance to composite geometry
    void getDist(real dist[], vect3d n[],
                              const vect3d loc[], int sz) const  {
      vector<real> dists(sz) ;
      vector<vect3d> ns(sz) ;

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
   Loci::CPTR<distFunc> transform(rigid_transform Tv) const {
    std::vector<Loci::CPTR<distFunc> > newGeomList(geomList.size()) ;
    for(size_t i=0;i<geomList.size();++i)
      newGeomList[i] = geomList[i]->transform(Tv) ;
    return new compositeFunc(newGeomList) ;
  }
} ;


#endif // INCOMPOSITE_H
 
