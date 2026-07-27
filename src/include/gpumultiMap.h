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
#ifndef LOCI_GPUMULTIMAP_H_
#define LOCI_GPUMULTIMAP_H_

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/debug.h>
#include <Map_rep.h>
#include <gpuMap.h>
#include <store.h>

namespace Loci {

  class gpumultiMapRepI : public gpuMapRep {
    entitySet store_domain ;
    Entity **base_ptr ;
  public:
    gpumultiMapRepI() { base_ptr = 0 ; }
    gpumultiMapRepI(const store<int> &sizes) {
      base_ptr = 0 ;
      allocate(sizes) ; }
    void allocate(const store<int> &sizes) ;
    virtual void allocate(const entitySet &ptn) ;
    virtual ~gpumultiMapRepI() ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP MapRemap(const dMap &dm, const dMap &rm) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
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
    virtual void pack(void *ptr, int &loc,
                      int &size, const entitySet &e, const Map& remap) ;
    virtual void unpack(void *ptr, int &loc,
                        int &size, const sequence &seq, const dMap& remap) ;
    
    virtual entitySet domain() const ;

    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
    preimage(const entitySet &codomain) const ;
    virtual storeRepP get_map() ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
#ifdef H5_HAVE_PARALLEL 
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en, hid_t xfer_plist_id) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist_id) const ;
#endif
    Entity ** get_base_ptr() const {
      Entity ** p = 0 ;
      if(alloc_id>=0)
        p = ((Entity **)GPUstoreAllocateData[alloc_id].alloc_ptr2
             - GPUstoreAllocateData[alloc_id].base_offset) ;
      return p ;
    }
    Entity *begin(int indx) { return base_ptr[indx] ; }
    Entity *end(int indx) { return base_ptr[indx+1] ; }
    const Entity *begin(int indx) const { return base_ptr[indx] ; }
    const Entity *end(int indx) const { return base_ptr[indx+1] ; }
    int vec_size(int indx) const { return end(indx)-begin(indx) ; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
    virtual void copyFrom(const storeRepP &fromMap, entitySet set) ;
    virtual store_type RepType() const ;
  private:
    virtual storeRepP expand(entitySet &out_of_dom, std::vector<entitySet> &init_ptn) ;
  } ;

  template <class T> class multiArrayHelper_const {
    const T *first ;
    const T *last ;
  public:
    GPU_DECL multiArrayHelper_const(const T *f, const T *l) : first(f), last(l){}
    GPU_DECL int size() { return last-first ; }
    GPU_DECL const T &operator[](Entity indx) { return first[indx] ; }
    GPU_DECL const T &operator[](size_t indx) { return first[indx] ; }
    GPU_DECL const T &operator[](unsigned int indx) { return first[indx] ; }
    GPU_DECL  const T &operator[](unsigned char indx) { return first[indx] ; }
    GPU_DECL const T *begin() { return first ; }
    GPU_DECL const T *end() { return last; }
  } ;
  
  template <class T> class constMultiAccessor {
    const T * const * accessor ;
  public:
    void operator=(const T * const * acc) { accessor = acc; } 
    GPU_DECL multiArrayHelper_const<T> operator[](Entity indx) {
      return multiArrayHelper_const<T>(accessor[indx],accessor[indx+1]) ;
    }
    GPU_DECL multiArrayHelper_const<T> operator[](Entity indx) const {
      return multiArrayHelper_const<T>(accessor[indx],accessor[indx+1]) ;
    }
  } ;
    
    
  class gpumultiMap : public store_instance {
    friend class const_gpumultiMap ;
    typedef gpumultiMapRepI MapType ;
    Entity **base_ptr ;
  public:

    class arrayHelper {
      Entity *first ;
      Entity *last ;
    public:
      arrayHelper(Entity *f, Entity *l) : first(f), last(l){}
      int size() { return last-first ; }
      Entity &operator[](int indx) { return first[indx] ; }
      Entity &operator[](size_t indx) { return first[indx] ; }
      Entity &operator[](unsigned int indx) { return first[indx] ; }
      Entity &operator[](unsigned char indx) { return first[indx] ; }
      Entity *begin() { return first ; }
      Entity *end() { return last; }
    } ;
    class arrayHelper_const {
      const Entity *first ;
      const Entity *last ;
    public:
      arrayHelper_const(const Entity *f, const Entity *l) : first(f), last(l){}
      int size() { return last-first ; }
      const Entity &operator[](Entity indx) { return first[indx] ; }
      const Entity &operator[](size_t indx) { return first[indx] ; }
      const Entity &operator[](unsigned int indx) { return first[indx] ; }
      const Entity &operator[](unsigned char indx) { return first[indx] ; }
      const Entity *begin() { return first ; }
      const Entity *end() { return last; }
    } ;

    gpumultiMap(const gpumultiMap &var) { setRep(var.Rep()) ; }
    gpumultiMap & operator=(const gpumultiMap &str)
    { setRep(str.Rep()) ; return *this ;}
    
    gpumultiMap() { setRep(new MapType) ; }

    gpumultiMap(storeRepP p) { setRep(p) ; }
    
    virtual ~gpumultiMap() ;
    virtual void notification() ;

    gpumultiMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    entitySet domain() const { return Rep()->domain() ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }
    arrayHelper elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper(base_ptr[indx],base_ptr[indx+1]) ;
    }
    arrayHelper_const const_elem(int indx)  const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper_const(base_ptr[indx],base_ptr[indx+1]) ;
    }

    arrayHelper operator[](Entity indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper(base_ptr[indx],base_ptr[indx+1]) ;
    }
    arrayHelper_const operator[](Entity indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper_const(base_ptr[indx],base_ptr[indx+1]) ;
    }
    arrayHelper operator[](size_t indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper(base_ptr[indx],base_ptr[indx+1]) ;
    }
    arrayHelper_const operator[](size_t indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper_const(base_ptr[indx],base_ptr[indx+1]) ;
    }
      
    int num_elems(int indx) const {return base_ptr[indx+1]-base_ptr[indx];}
    Entity *begin(int indx) { return base_ptr[indx] ; }
    Entity *end(int indx) { return base_ptr[indx+1] ; }
    const Entity *begin(int indx) const { return base_ptr[indx] ; }
    const Entity *end(int indx) const { return base_ptr[indx+1] ; }
    int vec_size(int indx) const { return end(indx)-begin(indx) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
    int getRangeKeySpace() const { return MapRepP(Rep())->getRangeKeySpace() ; }
    void setRangeKeySpace(int v) { MapRepP(Rep())->setRangeKeySpace(v) ; }
    Entity **ptr() { return base_ptr ; }
    
  } ;
  
  inline std::ostream & operator<<(std::ostream &s, const gpumultiMap &m)
  { return m.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, gpumultiMap &m)
  { return m.Input(s) ; }

  class const_gpumultiMap : public store_instance {
    typedef gpumultiMapRepI MapType ;
    const Entity * const * base_ptr ;
    const_gpumultiMap(const_gpumultiMap &var) {  setRep(var.Rep()) ; }
    const_gpumultiMap & operator=(const const_gpumultiMap &str)
    { setRep(str.Rep()) ; return *this ;}
    const_gpumultiMap & operator=(const gpumultiMap &str)
    { setRep(str.Rep()) ; return *this ;}
  public:
    class arrayHelper_const {
      const Entity *first ;
      const Entity *last ;
    public:
      arrayHelper_const(const int *f, const int *l) : first(f), last(l){}
      int size() { return last-first ; }
      const Entity &operator[](int indx) { return first[indx] ; }
      const Entity &operator[](size_t indx) { return first[indx] ; }
      const Entity &operator[](unsigned int indx) { return first[indx] ; }
      const Entity &operator[](unsigned char indx) { return first[indx] ; }
      const Entity *begin() { return first ; }
      const Entity *end() { return last; }
      
    } ;

    const_gpumultiMap() { setRep(new MapType) ; }
    
    const_gpumultiMap(gpumultiMap &var) { setRep(var.Rep()) ; }
    
    const_gpumultiMap(storeRepP rp) { setRep(rp) ; }
    
    virtual ~const_gpumultiMap() ;
    virtual void notification() ;
    
    virtual instance_type access() const ;
    
    const_gpumultiMap & operator=(storeRepP p) { setRep(p) ; return *this ;}
    
    entitySet domain() const { return Rep()->domain(); }
    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }
    arrayHelper_const const_elem(Entity indx)  const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper_const(base_ptr[indx],base_ptr[indx+1]); }
    arrayHelper_const operator[](Entity indx) const { 
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper_const(base_ptr[indx],base_ptr[indx+1]); }
    arrayHelper_const operator[](size_t indx) const { 
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return arrayHelper_const(base_ptr[indx],base_ptr[indx+1]); }
    int num_elems(int indx) const {return base_ptr[indx+1]-base_ptr[indx];}
    const int *begin(int indx) const { return base_ptr[indx] ; }
    const int *end(int indx) const { return base_ptr[indx+1] ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    int getRangeKeySpace() const { return MapRepP(Rep())->getRangeKeySpace() ; }
    const Entity * const *ptr() const { return base_ptr ; }
  } ;

}


#endif
