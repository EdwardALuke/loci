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
#ifndef GPUMULTISTORE_DEF_H
#define GPUMULTISTORE_DEF_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <Tools/debug.h>
#include <Tools/tools.h>
#include <data_traits.h>
#include <store_rep.h>
#include <gpurep.h>
#include <DMap.h>
#include <multiMap.h>
#include <dist_internal.h>

namespace Loci {

  template<class T> class gpumultiStoreRepI : public gpuRep {
    entitySet store_domain ;
    T **base_ptr ;

    void copyFrom(const storeRepP &p, entitySet set) ;
    void copyTo(storeRepP &p, entitySet set) const ;

    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &en) ;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en) const ;
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en) ;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en) const ;
#ifdef H5_HAVE_PARALLEL
    void hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &en, hid_t xfer_plist_id) ;
    void hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en, hid_t xfer_plist_id) const ;
    void hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en, hid_t xfer_plist_id) ;
    void hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en, hid_t xfer_plist_id) const ;
#endif

    int get_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset) ;
    int get_estimated_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset) ;
    int get_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset) ;
    int get_estimated_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset) ;

    void packdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e) ;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e) ;

    void unpackdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq) ;
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info get_frame_info(IDENTITY_CONVERTER g) ;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g) ;
  public:
    gpumultiStoreRepI()
    { base_ptr = 0 ; }

    gpumultiStoreRepI(const entitySet &p)
    { base_ptr = 0 ; store_domain = p ; }

    gpumultiStoreRepI(const store<int> &sizes)
    { base_ptr = 0 ; allocate(sizes) ; }

    void allocate(const store<int> &sizes) ;
    virtual void shift(int_type offset) ;
    void setSizes(const const_multiMap &mm) ;
    virtual ~gpumultiStoreRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
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

    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
#ifdef H5_HAVE_PARALLEL
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en, hid_t xfer_plist_id) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist_id) const ;
#endif
    T ** get_base_ptr() const {
      T **p = 0 ;
      if(alloc_id>=0)
        p = ((T **)GPUstoreAllocateData[alloc_id].alloc_ptr2) -
          GPUstoreAllocateData[alloc_id].base_offset ;
      return p ;
    }
    T *begin(int indx) { return base_ptr[indx] ; }
    T *end(int indx) { return base_ptr[indx+1] ; }
    const T *begin(int indx) const { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
  } ;


  template<class T> class gpumultiStore : public store_instance {
    typedef gpumultiStoreRepI<T> storeType ;
    T ** base_ptr ;
    gpumultiStore(const gpumultiStore<T> &var) { setRep(var.Rep()) ; }
    gpumultiStore<T> & operator=(const gpumultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef Vect<T> containerType ;
    gpumultiStore() { setRep(new storeType) ; }
    gpumultiStore(storeRepP rp) { setRep(rp) ; }
    virtual ~gpumultiStore() {}
    virtual void notification() ;
    gpumultiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const store<int> &sizes) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->allocate(sizes) ;
    }
    void setSizes(const const_multiMap &m) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->setSizes(m) ;
    }
    
    const entitySet domain() const { return Rep()->domain() ; }
    Vect<T> elem(Entity indx)
    {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ;
    }
    Vect<T> operator[](Entity indx)
    {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ;
    }
    Vect<T> operator[](size_t indx)
    {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ;
    }

    const_Vect<T> operator[](size_t indx) const {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return const_Vect<T>(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ;
    }
    T *begin(int indx) { return base_ptr[indx] ; }
    T *end(int indx) { return base_ptr[indx+1] ; }
    const T *begin(int indx) const { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
    int vec_size(int index) const { return end(index)-begin(index) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

  } ;

  template<class T> class const_gpumultiStore : public store_instance {
    typedef gpumultiStoreRepI<T> storeType ;
    T ** base_ptr ;
    const_gpumultiStore(const const_gpumultiStore<T> &var) { setRep(var.Rep()) ; }
    const_gpumultiStore(const gpumultiStore<T> &var) { setRep(var.Rep()) ; }
    const_gpumultiStore<T> & operator=(const gpumultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_gpumultiStore<T> & operator=(const const_gpumultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef const_Vect<T> containerType ;
    const_gpumultiStore() { setRep(new storeType) ; }
    const_gpumultiStore(storeRepP rp) { setRep(rp) ; }
    virtual ~const_gpumultiStore() {}
    virtual void notification() ;
    virtual instance_type access() const ;
    const_gpumultiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    const entitySet domain() const { return Rep()->domain() ; }
    const T * const * ptr() { return base_ptr; }
    containerType elem(Entity indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return containerType(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }
    containerType operator[](Entity indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return containerType(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }
    containerType operator[](size_t indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL);
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif
      return containerType(base_ptr[indx],base_ptr[indx+1]-base_ptr[indx]) ; }

    const T *begin(int indx) const { return base_ptr[indx] ; }
    const T *end(int indx) const { return base_ptr[indx+1] ; }
    int vec_size(int index) const { return end(index)-begin(index) ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

} // end of namespace Loci

#endif
