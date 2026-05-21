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
#ifndef DMAP_H_
#define DMAP_H_

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Tools/debug.h>
#include <Map_rep.h>

#include <Tools/block_hash.h>


namespace Loci {

  class Map ;
  class multiMap ;
  class dMap;
  
  class dMapRepI : public MapRep {
    block_hash<int> attrib_data;
  public:
    dMapRepI() { }
    dMapRepI(const entitySet &p) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual void erase(const entitySet& rm) ;
    virtual void invalidate(const entitySet& valid) ;
    virtual void guarantee_domain(const entitySet& include) ;
    virtual ~dMapRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP MapRemap(const dMap &dm, const dMap &rm) const ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 MPI_Comm comm LOCI_DEFAULT_COMM) ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 const dMap& remap, MPI_Comm comm LOCI_DEFAULT_COMM) ;
    virtual void compose(const dMap &m, const entitySet &context) ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet& e, entitySet& packed) ;
    virtual int pack_size(const entitySet &e) ;
    virtual int estimated_pack_size(const entitySet &e) ;

    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual entitySet domain() const ;

    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
      preimage(const entitySet &codomain) const ;
    virtual storeRepP get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;

#ifdef H5_HAVE_PARALLEL 
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset, hid_t xfer_plist_id) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist_id) const ;
#endif
    
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual block_hash<int> *get_attrib_data() { return &attrib_data; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
  } ;
      
  class dMap : public store_instance {
    friend class const_dMap ;
    typedef dMapRepI MapType ;
    block_hash<int> *attrib_data;
  public:
    // These should be private, but too much code currently depends on it
    // being public.  This code is dangerous because it does a shallow
    // copy which means that results sometimes may be unpredictable when
    // these operations are used.
    dMap(const dMap &var) { setRep(var.Rep()) ; }
    dMap & operator=(const dMap &str) { setRep(str.Rep()) ; return *this ;}

    dMap() { setRep(new MapType) ;}
    dMap(storeRepP &rp) { setRep(rp) ; }

    virtual ~dMap() ;

    virtual void notification() ;
    
    dMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    // this method does erases the domain of the Map that are
    // inside the passed in parameter
    void erase(const entitySet& rm) { Rep()->erase(rm) ; }

    entitySet domain() const { return Rep()->domain() ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }

    int &elem(int indx) { 
      return (*attrib_data)[indx]; 
    }

    const int &const_elem(int indx)  const 
    { 
      return attrib_data->elem(indx); 
    }

    int &operator[](int indx) { return elem(indx); }

    const int &operator[](int indx) const { return const_elem(indx) ; }
    const int &operator()(int indx) const { return const_elem(indx) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

    entitySet image(const entitySet &dom) const {
      return MapRepP(Rep())->image(dom) ;
    }
    std::pair<entitySet,entitySet> preimage(const entitySet &codomain) const {
      return MapRepP(Rep())->preimage(codomain) ;
    }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const dMap &m)
    { return m.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, dMap &m)
    { return m.Input(s) ; }

  class const_dMap : public store_instance {
    typedef dMapRepI MapType ;
    block_hash<int>  *attrib_data;
    const_dMap(const const_dMap &var) {setRep(var.Rep()) ; }
    const_dMap(const dMap &var) {setRep(var.Rep()); }
    const_dMap & operator=(const dMap &str) { setRep(str.Rep()) ; return *this ;}

    const_dMap & operator=(const const_dMap &str)
    { setRep(str.Rep()) ; return *this ;}

  public:
    const_dMap()
    { setRep(new MapType); }
    const_dMap(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_dMap() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_dMap & operator=(storeRepP p) { setRep(p) ; return *this ;}

    entitySet domain() const { return Rep()->domain(); }

    operator MapRepP()
    { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }

    const int &const_elem(int indx)  const {
      return attrib_data->elem(indx); 
    }

    int operator[](int indx) const { return const_elem(indx) ; }
    int operator()(int indx) const { return const_elem(indx) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

        
  inline std::ostream & operator<<(std::ostream &s, const const_dMap &m)
    { return m.Print(s) ; }
 
  void inverseMap(multiMap &result, const dMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage);
}

#endif
