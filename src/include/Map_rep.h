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
#ifndef MAP_REP_H
#define MAP_REP_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <store_rep.h>

namespace Loci {

  class multiMap ;
  class Map ;
  
  class MapRep : public storeRep {
    int rangeKeySpace ;
  public:
    MapRep() { rangeKeySpace = 0 ; }
    virtual ~MapRep() ;
    virtual entitySet image(const entitySet &domain) const = 0 ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const = 0 ;
    virtual storeRepP get_map() = 0 ;
    virtual void compose(const dMap &m, const entitySet &context) = 0 ;
    virtual store_type RepType() const ;
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn, MPI_Comm comm LOCI_DEFAULT_COMM) = 0 ;
    virtual storeRepP MapRemap(const dMap &dm, const dMap &im) const = 0 ;
    virtual void shift(int_type)
    {std::cerr<<"shift for Map has not been implemented!"<<std::endl ;}
    int getRangeKeySpace() const { return rangeKeySpace ; }
    void setRangeKeySpace(int v) { rangeKeySpace = v ; }
  } ;

  typedef NPTR<MapRep> MapRepP ;
}

#endif
