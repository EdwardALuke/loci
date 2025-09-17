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
#ifndef DMAPVEC_DEF_H
#define DMAPVEC_DEF_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Tools/hash_map.h>

#include <Tools/debug.h>
#include <store_rep.h>
#include <Map_rep.h>
#include <store.h>
#include <multiMap.h>
#include <Loci_types.h>

namespace Loci {
  template<unsigned int M> class dMapVec ;
  template<unsigned int M> class const_dMapVec ;

  template <unsigned int M> class dMapVecRepI : public MapRep {
  public:
    typedef Array<int,M>   VEC;
  private:
    HASH_MAP(int,VEC)    attrib_data;
  public:
    dMapVecRepI() {}
    dMapVecRepI(const entitySet &p) { allocate(p); }
    virtual void allocate(const entitySet &ptn) ;
    virtual ~dMapVecRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP MapRemap(const dMap &dm, const dMap &rm) const ;
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
    virtual void unpack(void *ptr, int &loc, int &size,  const sequence &seq)  ;
    
    virtual entitySet domain() const ;

    virtual entitySet image(const entitySet &domain) const ;

    virtual std::pair<entitySet,entitySet>
    preimage(const entitySet &codomain) const ;

    virtual storeRepP get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace,hid_t dataset, hsize_t dimension, const char* name, entitySet &en) const ;
#ifdef H5_HAVE_PARALLEL 
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset, hid_t xfer_plist_id) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace,hid_t dataset, hsize_t dimension, const char* name, entitySet &en, hid_t xfer_plist_id) const ;
#endif
    
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
    virtual DatatypeP getType() ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    HASH_MAP(int,VEC) *get_attrib_data() { return &attrib_data; }
    virtual frame_info get_frame_info() ;
  } ;

  template<unsigned int M> class dMapVec : public store_instance 
  {
    friend  class const_dMapVec<M> ;
    typedef dMapVecRepI<M> MapVecType ;
    typedef typename MapVecType::VEC VEC ;
    HASH_MAP(int, VEC)     *attrib_data;
    dMapVec(const dMapVec<M> &var) { setRep(var.Rep()) ; }
    dMapVec & operator=(const dMapVec<M> &str) { 
      setRep(str.Rep()) ; 
      return *this ;
    }
  public:
    dMapVec() { setRep(new MapVecType) ; }
    dMapVec(storeRepP &rp) { setRep(rp) ; }
    virtual ~dMapVec() ;
    virtual void notification() ;
    dMapVec &operator=(storeRepP p) { 
      setRep(p) ; 
      return *this ;
    }
    // -----------------------------------------------------------------
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    // -----------------------------------------------------------------
    entitySet domain() const { 
      return Rep()->domain() ; 
    }
    // -----------------------------------------------------------------
    //    operator storeRepP() { return Rep() ; }
    // -----------------------------------------------------------------
    operator MapRepP() { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }
    // -----------------------------------------------------------------
    VEC &elem(int indx) {
      return (*attrib_data)[indx]; 
    }
    // -----------------------------------------------------------------
    const VEC &const_elem(int indx)  const { 
      typename HASH_MAP(int,VEC)::const_iterator ci;
      ci = attrib_data->find(indx);
      if( ci != attrib_data->end()) return ci->second();
    }
    // -----------------------------------------------------------------
    VEC &operator[](int indx) { 
      return elem(indx); 
    }
    // -----------------------------------------------------------------
    const VEC &operator[](int indx) const { return const_elem(indx) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
    //  entitySet image(const entitySet &dom) const;
  } ;
  
  template<unsigned int M> class const_dMapVec : public store_instance {
    typedef dMapVecRepI<M> MapVecType ;
    typedef typename MapVecType::VEC VEC ;
    HASH_MAP(int,VEC)  *attrib_data;
    const_dMapVec(const const_dMapVec<M> &var) { setRep(var.Rep()) ; } 
    const_dMapVec(const dMapVec<M> &var) { setRep(var.Rep()) ; }
    const_dMapVec & operator=(const const_dMapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
    const_dMapVec & operator=(const dMapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}
  public:
    const_dMapVec() { setRep(new MapVecType) ; }
    const_dMapVec(storeRepP &rp) { setRep(rp) ; }
    virtual ~const_dMapVec() ;
    virtual void notification() ;
    virtual instance_type access() const ;
    const_dMapVec & operator=(storeRepP p) { setRep(p) ; return *this ;}
    entitySet domain() const { return Rep()->domain() ; }
    operator MapRepP() { MapRepP p(Rep()) ; fatal(p==0) ; return p ; }
    const VEC &const_elem(int indx)  const {
      typename HASH_MAP(int,VEC)::const_iterator ci;
      ci = attrib_data->find(indx);
      return( ci->second );
    }
    const VEC &operator[](int indx) const { return const_elem(indx) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;


} // end of namespace Loci

#endif
