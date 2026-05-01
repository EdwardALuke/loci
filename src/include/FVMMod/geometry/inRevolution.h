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
#ifndef INREVOLUTION_H
#define INREVOLUTION_H
#include "geomTest.h"

/// @brief Class that creates a bounding cylinder for two given points and a
///        radius.
 /// Body of revolution geometry type specification
  class inRevolution : public geomTest {
    vector3d<real_t> p1,p2 ; // revolution axis endpoints
    real_t rlen2 ; // reciprocal of len(pt1-pt2)^2
    vector<pair<real_t,real_t> > radius_pairs ;
  public:
    inRevolution() {
      p1=vector3d<real_t>(0,0,0) ;
      p2=vector3d<real_t>(0,0,0) ;
    }

    inRevolution(vector3d<real_t> p1i,
                   vector3d<real_t> p2i,
                   vector<pair<real_t,real_t> > rp ) :
        p1(p1i),p2(p2i),radius_pairs(rp)
    {
      std::sort(radius_pairs.begin(),radius_pairs.end()) ;
      rlen2 = 1./dot(p2-p1,p2-p1) ;
    }

    inRevolution(const options_list &ol){
    // Here get revolution information
      if(ol.optionExists("p1")     ||
      !ol.optionExists("p2") ||
      !ol.optionExists("radius"))
      {
      printf("ERROR::inRevolution: Need two axis points (p1 and p2) and a radii pairs\n");
      Loci::Abort();
      }
      // For a revolution, a series of radius value and offset values are needed
      Loci::option_value_type ovt= ol.getOptionValueType("radius") ;
      if(ovt != Loci::LIST) {
      cerr << "radius for geometry of component '" << tag
              << "' must be a list" << endl ;
      Loci::Abort() ;
      }

      // Get the points set first
      ol.getOptionUnits("p1",     "m", p1);
      ol.getOptionUnits("p2",     "m", p2);

      rlen2 = 1./dot(p2-p1,p2-p1);

      // Get the radii lists
      Loci::options_list::arg_list value_list ;
      ol.getOption("radius",value_list) ;
      int sz = value_list.size() ;
      for(int i=0;i<sz;++i) {
        if(value_list[i].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[i].get_value(vu) ;
          if(!vu.is_compatible("meter")) {
            std::cerr << "wrong type of units for radius in component"
                      << " geometry for component '" << tag << "'"
                      << ": " << vu << std::endl ;
            Loci::Abort() ;
          }
          radius.push_back(vu.get_value_in("meter")) ;
        } else if(value_list[i].type_of() == Loci::REAL) {
          real_t r =0 ;
          value_list[i].get_value(r) ;
          radius.push_back(r) ;
        } else {
          cerr << "wrong type in radius list for geometry of component '"
                << tag << "'" << endl ;
          Loci::Abort() ;
        }
      }

      // Now get the offsets list
      if(ol.optionExists("offsets")) {
        Loci::option_value_type ovt= ol.getOptionValueType("offsets") ;
        if(ovt != Loci::LIST) {
          cerr << "offsets for geometry of component '" << tag
                << "' must be a list" << endl ;
          Loci::Abort() ;
        }
        Loci::options_list::arg_list value_list ;
        ol.getOption("offsets",value_list) ;
        int sz = value_list.size() ;
        for(int i=0;i<sz;++i) {
          if(value_list[i].type_of() == Loci::REAL) {
            real_t o =0 ;
            value_list[i].get_value(o) ;
            offsets.push_back(o) ;
          } else {
            cerr << "wrong type in offset list for geometry of component '"
                  << tag << "'" << endl ;
            Loci::Abort() ;
          }
        }
      }

      // Check data is the same size
      if(radius.size() != offsets.size()) {
        cerr << "radius and offset should have same number of entries in "
            << "geometry for component '" << tag << "'" << endl ;
        Loci::Abort() ;
      }

      // Create the radius-offset pairs
      int sz = radius.size() ;
      radius_pair.resize(sz);
      for(int i=0;i<sz;++i) {
        radius_pair[i].first = offsets[i] ;
        radius_pair[i].second = radius[i] ;
      }
    } // End of inRevolution()

  bool inGeomPt(vector3d<real_t> pt) const {
    vector3d<real_t> p = pt-p1 ;
    real_t t = dot(p,p2-p1)*rlen2 ;
    if(t < 0 || t > 1) {
      return false ;
    }

    vector3d<real_t> v = t*(p2-p1) ;
    real_t r = 0 ;
    for(size_t i=1;i<radius_pairs.size();++i) {
      if(t >= radius_pairs[i-1].first && t < radius_pairs[i].first) {
        real_t t1 = radius_pairs[i-1].first ;
        real_t t2 = radius_pairs[i].first ;
        real_t w = (t1-t)/(t1-t2) ;
        r = radius_pairs[i-1].second*(1.-w)+radius_pairs[i].second*w ;
      }
    }

    real_t r2axis = dot(v-p,v-p) ;
    if(r2axis < r*r) {
      return true ;
    }
    return false ;
  }

  real_t distToSurface(vector3d<real_t> pt) const {
    vector3d<real_t> p = pt-p1 ;
    real_t t = dot(p,p2-p1)*rlen2 ;
    if(t < 0)
      return norm(p) ;
    if(t > 1)
      return norm(pt-p2) ;

    vector3d<real_t> v = t*(p2-p1) ;
    real_t r = 0 ;
    for(size_t i=1;i<radius_pairs.size();++i)
      if(t >= radius_pairs[i-1].first && t < radius_pairs[i].first) {
        real_t t1 = radius_pairs[i-1].first ;
        real_t t2 = radius_pairs[i].first ;
        real_t w = (t1-t)/(t1-t2) ;
        r = radius_pairs[i-1].second*(1.-w)+radius_pairs[i].second*w ;
      }

    real_t r2axis = dot(v-p,v-p) ;
    if(r2axis < r*r)
      return 0 ;
    return sqrt(r2axis)-r ;
  }

}
#endif // INREVOLUTION_H
