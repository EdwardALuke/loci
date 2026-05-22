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
#ifndef DMULTIMAP_H
#define DMULTIMAP_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Tools/debug.h>
#include <Map_rep.h>
#include <DMap.h>
#include <vector>

#include <Tools/block_hash.h>
#include <Tools/malloc_alloc.h>

#include <Map.h>
#include <store.h>
namespace Loci {

  class dmultiMapRepI : public MapRep {
    block_hash<std::vector<int,malloc_alloc<int> > >  attrib_data;
  public:
    dmultiMapRepI() { }
    dmultiMapRepI(const store<int> &sizes) { allocate(sizes) ; }
    void allocate(const store<int> &sizes) ;
    virtual void allocate(const entitySet &ptn) ;
    virtual void erase(const entitySet& rm) ;
    virtual void invalidate(const entitySet& valid) ;
    virtual void guarantee_domain(const entitySet& include) ;
    virtual ~dmultiMapRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP MapRemap(const dMap &dm, const dMap &rm) const ;
    virtual void compose(const dMap &m, const entitySet &context) ;
    virtual void copy(storeRepP &st, const entitySet &context)  ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context)  ;
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
    
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn, MPI_Comm comm LOCI_DEFAULT_COMM) ;
    virtual storeRepP thaw() ;
    virtual storeRepP freeze() ;
    block_hash<std::vector<int, malloc_alloc<int> > > *get_attrib_data() {return &attrib_data;}
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
  } ;
  
  //***************************************************************************
  
  class dmultiMap : public store_instance {
    friend class const_dmultiMap ;
    typedef dmultiMapRepI MapType ;
    block_hash<std::vector<int,malloc_alloc<int> > > *attrib_data;
    dmultiMap(const dmultiMap &var) { setRep(var.Rep()) ; }
    dmultiMap & operator=(const dmultiMap &str)
    { setRep(str.Rep()) ; return *this ;}

  public:
    dmultiMap() { setRep(new MapType) ; }
        
    dmultiMap(const store<int> &sizes) { setRep( new MapType(sizes) ); }
        

    dmultiMap(storeRepP p) { setRep(p) ; }
    
    virtual ~dmultiMap() ;

    virtual void notification() ;

    dmultiMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    void allocate(const store<int> &sizes) {
      NPTR<MapType> p(Rep()) ;
      p->allocate(sizes) ; 
    }

    entitySet domain() const { return Rep()->domain() ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; 
    }

    std::vector<int,malloc_alloc<int> > &elem(int indx) {
      return (*attrib_data)[indx]; 
    }

    const std::vector<int,malloc_alloc<int> > &const_elem(int indx)  const 
    {
      return attrib_data->elem(indx) ;
    }

    std::vector<int,malloc_alloc<int> > &operator[](int indx) { return elem(indx); }
    
    const std::vector<int,malloc_alloc<int> > &operator[](int indx) const 
    { return const_elem(indx) ; }
    
    int num_elems(int indx) const 
    { return const_elem(indx).size() ; }

    std::ostream &Print(std::ostream &s) const 
    { return Rep()->Print(s) ; }

    std::istream &Input(std::istream &s) 
    { return Rep()->Input(s) ; }
  } ;
  
  //***************************************************************************

  inline std::ostream & operator<<(std::ostream &s, const dmultiMap &m)
  { return m.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, dmultiMap &m)
  { return m.Input(s) ; }

  //***************************************************************************

  class const_dmultiMap : public store_instance {
    typedef dmultiMapRepI      MapType ;
    block_hash<std::vector<int,malloc_alloc<int> > > *attrib_data;
    const_dmultiMap(const const_dmultiMap &var) {  setRep(var.Rep()) ; }
    const_dmultiMap(const dmultiMap &var) { setRep(var.Rep()) ; }
    const_dmultiMap & operator=(const const_dmultiMap &str)
    { setRep(str.Rep()) ; return *this ;}

    const_dmultiMap & operator=(const dmultiMap &str)
    { setRep(str.Rep()) ; return *this ;}

  public:
    const_dmultiMap() { setRep(new MapType) ; }
    
    
    
    const_dmultiMap(storeRepP &rp) { setRep(rp) ; }
    
    virtual ~const_dmultiMap() ;
    virtual void notification() ;
    
    virtual instance_type access() const ;
    
    const_dmultiMap & operator=(storeRepP p) 
    { setRep(p) ; return *this ;}
    
    entitySet domain() const { return Rep()->domain(); }

    operator MapRepP() 
    {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; 
    }

    const std::vector<int,malloc_alloc<int> > &const_elem(int indx)  const {
      return attrib_data->elem(indx) ;
    }

    const std::vector<int,malloc_alloc<int> > operator[](int indx) const 
    { return const_elem(indx) ; }
    
    int num_elems(int indx) const 
    {return const_elem(indx).size() ;}

    std::ostream &Print(std::ostream &s) const 
    { return Rep()->Print(s) ; }
  } ;

  //***************************************************************************

  inline std::ostream & operator<< (std::ostream &s, 
                                    const const_dmultiMap &m)
  { return m.Print(s) ; }

  //***************************************************************************

  void inverseMap(dmultiMap &result, const dMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(dmultiMap &result, const Map &input_map,
                  const entitySet &input_image,
		  const entitySet &input_preimage) ;
  void inverseMap(dmultiMap &result, const const_dMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(dmultiMap &result, const const_Map &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(dmultiMap &result, const dmultiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;
  void inverseMap(dmultiMap &result, const multiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) ;

}

#endif
