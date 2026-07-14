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
#ifndef GMAP_REP_H
#define GMAP_REP_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <gstore_rep.h>

namespace Loci {
  
  class gMapRep : public gStoreRep {
  public:
    virtual ~gMapRep(){} ;

    virtual void set_image_space(gKeySpaceP space) = 0;
    virtual gKeySpaceP get_image_space()const = 0;
    virtual size_t size() const = 0;
    virtual gStoreRepP clone() const = 0;
    virtual gEntitySet image(const gEntitySet &domain) const = 0;
    //remove this function because gMap and gMultiMap have different return value 
    // virtual gEntitySet image(gEntity domain) const = 0;
    virtual gEntitySet image() const = 0;
    virtual std::pair<gEntitySet,gEntitySet>
    preimage(const gEntitySet &codomain) const = 0;
    virtual gStoreRepP get_map()const{
      std::cerr << "gMapRep:get_map() should not be called"
                << std::endl ;
      return gStoreRepP(0);
    }
    //this method remap the SECOND field of this using m
    virtual void inplace_compose(const gMap &m, MPI_Comm comm LOCI_DEFAULT_COMM) = 0 ;
    //this method remap the SECOND field of this using m and return a new map
    // virtual gStoreRepP  recompose(const gMap &m, MPI_Comm comm LOCI_DEFAULT_COMM) = 0 ;
    virtual gstore_type RepType() const {return GMAP;}

    //different from traditional maps, this method is const method
    //dom is the domain after expansion, not out_of_dom
    virtual gStoreRepP expand(gEntitySet &dom, std::vector<gEntitySet> &init_ptn,MPI_Comm comm LOCI_DEFAULT_COMM)const = 0 ;
    virtual gStoreRepP local_inverse() const  = 0;
  } ;

  typedef CPTR<gMapRep> gMapRepP ;
}

#endif
