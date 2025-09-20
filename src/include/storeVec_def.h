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
#ifndef STOREVEC_DEF_H
#define STOREVEC_DEF_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <hdf5_readwrite.h>
#include <store_rep.h>
#include <vector>
#include <dist_internal.h>
#include <DMap.h>

namespace Loci {
  extern double total_memory_usage ;
  extern std::ofstream debugout ;
  //*******************************************************************/
  template <class T> struct Scalar {
    T val ;
    Scalar(T v) : val(v) { }
  } ;
  template <class T> Scalar<T> mk_Scalar(T v) { return Scalar<T>(v) ;} 
  //*******************************************************************/
  template <class T> class Vect ;
  template <class T> class const_Vect {
  public:
    friend class Vect<T> ;
    const T * restrict ptr ;
#ifdef BOUNDS_CHECK
    int size ;
#endif
  public:
    const_Vect(const T *p ,int sz) {
      ptr = p ;
#ifdef BOUNDS_CHECK
      size = sz ;
#endif
    }
#ifdef BOUNDS_CHECK
    const T &operator[](int idx) const {
      fatal(idx >= size || idx < 0) ;
      return ptr[idx] ;
    }
#else 
    const T &restrict operator[](int idx) const restrict {
      return ptr[idx] ;
    }
#endif
#ifdef BOUNDS_CHECK
    const T &operator[](size_t idx) const {
      fatal(idx >= size_t(size) || idx < 0) ;
      return ptr[idx] ;
    }
#else 
    const T &restrict operator[](size_t idx) const restrict {
      return ptr[idx] ;
    }
#endif
#ifdef BOUNDS_CHECK
    const T &operator[](unsigned int idx) const {
      fatal(idx >= size || idx < 0) ;
      return ptr[idx] ;
    }
#else 
    const T &restrict operator[](unsigned int idx) const restrict {
      return ptr[idx] ;
    }
#endif
#ifdef BOUNDS_CHECK
    const T &operator[](unsigned char idx) const {
      fatal(idx >= size || idx < 0) ;
      return ptr[idx] ;
    }
#else 
    const T &restrict operator[](unsigned char idx) const restrict {
      return ptr[idx] ;
    }
#endif
    operator const T *restrict () const restrict {
      return ptr ;
    }
  } ;
  //******************************************************************/
  template <class T> class Vect {
  public:
    T *ptr ;
    int size ;
  public:
    Vect() {};
    void setSize( int s ) {
      size = s;
    }
    int getSize() const { return size; }
    
    Vect(T *p ,int sz) {
      ptr = p ;
      size = sz ;
    }

    void operator=(const Vect<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ = *p2++ ;
    }


    void operator=(const const_Vect<T> &t) {
      T *p1 = ptr ;
      const T *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ = *p2++ ;
    }


    template <class S> void operator=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ = s.val ;
    }

    template <class S> void operator+=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ += s.val ;
    }
      
    template <class S> void operator*=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ *= s.val ;
    }
      
    template <class S> void operator-=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ -= s.val ;
    }
      
    template <class S> void operator/=(const Scalar<S> &s) {
      T *p1 = ptr ;
      for(int i=0;i<size;++i)
        *p1++ /= s.val ;
    }

    template <class S> void operator+=(const Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ += *p2++ ;
    }

    template <class S> void operator+=(const const_Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ += *p2++ ;
    }

    template <class S> void operator-=(const Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ -= *p2++ ;
    }

    template <class S> void operator-=(const const_Vect<S> &t) {
      T *p1 = ptr ;
      const S *p2 = t.ptr ;
#ifdef BOUNDS_CHECK
      fatal(size != t.size) ;
#endif
      
      for(int i=0;i<size;++i)
        *p1++ -= *p2++ ;
    }

    T &operator[](int idx) {
#ifdef BOUNDS_CHECK
      fatal(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    const T &operator[](int idx) const {
#ifdef BOUNDS_CHECK
      fatal(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }
    T &operator[](size_t idx) {
#ifdef BOUNDS_CHECK
      fatal(idx >= (size_t)size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    const T &operator[](size_t idx) const {
#ifdef BOUNDS_CHECK
      fatal(idx >= (size_t)size || idx < 0) ;
#endif
      return ptr[idx] ;
    }
    T &operator[](unsigned int idx) {
#ifdef BOUNDS_CHECK
      fatal(idx >= (unsigned int)size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    const T &operator[](unsigned int idx) const {
#ifdef BOUNDS_CHECK
      fatal(idx >= (unsigned int)size || idx < 0) ;
#endif
      return ptr[idx] ;
    }
    T &operator[](unsigned char idx) {
#ifdef BOUNDS_CHECK
      fatal(idx >= (unsigned int)size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    const T &operator[](unsigned char idx) const {
#ifdef BOUNDS_CHECK
      fatal(idx >= size || idx < 0) ;
#endif
      return ptr[idx] ;
    }

    operator T*() {
      return ptr ;
    }

    operator const T *() const {
      return ptr ;
    }
  } ;

  //******************************************************************/
  template<class T> class storeVecRepI : public storeRep {
    entitySet    store_domain ;
    int          size ;
    lmutex       mutex ;
    bool         isMat; //if this is a storeMat
    
    int  get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int  get_estimated_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &en) ;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &en) const;
#ifdef H5_HAVE_PARALLEL
    void hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &en, hid_t xfer_plist_id) ;
    void hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &en, hid_t xfer_plist_id) const;
#endif
    
    void packdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e) ;
    void unpackdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int &size,
                    const sequence &seq) ;

    int  get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    int  get_estimated_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en) ;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &en) const;
#ifdef H5_HAVE_PARALLEL
    void hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en, hid_t xfer_plist_id) ;
    void hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &en, hid_t xfer_plist_id) const;
#endif    
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const entitySet &e ) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int &size,
                    const sequence &seq) ;
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info get_frame_info(IDENTITY_CONVERTER g) ;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g) ;
  public:
    storeVecRepI(): size(0),isMat(false) { }
    
    storeVecRepI(const entitySet &p):size(0),isMat(false)  { allocate(p) ; }
    
    virtual ~storeVecRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &eset, const int* p) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual storeRepP thaw(const entitySet& es) const ;
#ifdef DYNAMICSCHEDULING
    virtual storeRepP freeze(const entitySet& es) const {
      std::cerr << "storeRep.freeze(e) is not implemented yet"
                << std::endl ;
      abort() ;
      return storeRepP(0) ;
    }
#endif
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    virtual int pack_size(const entitySet& e, entitySet& packed) ;
    virtual int pack_size(const entitySet &e ) ;
    virtual int estimated_pack_size(const entitySet &e ) ; 
    virtual void pack(void * ptr, int &loc, int &size, const entitySet &e ) ;
    virtual void unpack(void * ptr, int &loc, int &size, const sequence &seq) ;
#ifdef DYNAMICSCHEDULING
    virtual void pack(void* ptr, int& loc,
                      int& size, const entitySet& e, const Map& remap) {
      pack(ptr,loc,size,e) ;
    }
    virtual void unpack(void* ptr, int& loc,
                        int& size, const sequence& seq, const dMap& remap) {
      unpack(ptr,loc,size,seq) ;
    }
#endif
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
    virtual void set_elem_size(int sz) ;
    T *get_alloc_ptr() const { T *p = 0; if(alloc_id>=0) p=(T *)storeAllocateData[alloc_id].base_ptr ; return p ; }
    int get_base_offset() const { int o = 0 ;if(alloc_id>=0) o = storeAllocateData[alloc_id].base_offset ; return o ; }
    int get_size() const { int sz = 0 ; if(alloc_id>=0) sz = storeAllocateData[alloc_id].size ; return sz ; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
    void setIsMat(bool im){isMat=im;}
  } ;

  template<class T> class storeVec : public store_instance {
    typedef storeVecRepI<T> storeType ;
    T* alloc_ptr ;
    int base_offset ;
    int size ;
    storeVec(const storeVec<T> &var) {setRep(var.Rep()) ;}
    storeVec<T> & operator=(const storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef Vect<T> containerType ;
    storeVec() {setRep(new storeType) ;}
    storeVec(storeRepP rp) { setRep(rp) ;}
    virtual ~storeVec() {}
    virtual void notification() ;
    storeVec<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    void setVecSize(int sz) {
      Rep()->set_elem_size(sz) ;
    }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    int vecSize() const { return size ; }
    const entitySet domain() const { return Rep()->domain() ; }
    Vect<T> elem(Entity indx) {
#ifdef BOUNDS_CHECK
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(alloc_ptr+((indx-base_offset)*size),size) ; }
    Vect<T> operator[](Entity indx) {
#ifdef BOUNDS_CHECK
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(alloc_ptr+((indx-base_offset)*size),size) ; }
    Vect<T> operator[](size_t indx) {
#ifdef BOUNDS_CHECK
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Vect<T>(alloc_ptr+((indx-base_offset)*size),size) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
    
  } ;

  template<class T> class const_storeVec : public store_instance {
    typedef storeVecRepI<T> storeType ;
    const T* restrict alloc_ptr ;
    int base_offset ;
    int size ;
    const_storeVec(const const_storeVec<T> &var) {setRep(var.Rep()) ;}
    const_storeVec(const storeVec<T> &var) {setRep(var.Rep()) ;}
    const_storeVec<T> & operator=(const storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_storeVec<T> & operator=(const const_storeVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef const_Vect<T> containerType ;
    const_storeVec() { setRep(new storeType) ; }
    const_storeVec(storeRepP rp) { setRep(rp) ; }
    virtual ~const_storeVec() {}
    virtual void notification() ;
    virtual instance_type access() const ;
    const_storeVec<T> & operator=(storeRepP p) {
      setRep(p) ;
      return *this ;
    }
    int vecSize() const { return size ; }
    const entitySet domain() const { return Rep()->domain() ; }
#ifdef BOUNDS_CHECK
    const_Vect<T> elem(Entity indx) const {
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
      return const_Vect<T>(alloc_ptr+((indx-base_offset)*size),size) ; }
    const_Vect<T> operator[](Entity indx) const {
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
      return const_Vect<T>(alloc_ptr+((indx-base_offset)*size),size) ;
    }
    const_Vect<T> operator[](size_t indx) const {
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
      return const_Vect<T>(alloc_ptr+((indx-base_offset)*size),size) ;
    }
#else
    const_Vect<T> elem(Entity indx) const restrict {
      return const_Vect<T>(alloc_ptr+((indx-base_offset)*size),size) ; }
    const_Vect<T> operator[](Entity indx) const restrict {
      return const_Vect<T>(alloc_ptr+((indx-base_offset)*size),size) ;
    }
    const_Vect<T> operator[](size_t indx) const restrict {
      return const_Vect<T>(alloc_ptr+((indx-base_offset)*size),size) ;
    }
#endif
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }

  } ;


} // end of namespace Loci

#endif
