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
#ifndef DSTORE_DEF_H
#define DSTORE_DEF_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <algorithm>
#include <functional>

#include <mpi.h>

#include <Tools/debug.h>
#include <store_rep.h>
#include <Map.h>
#include <Tools/intervalSet.h>
#include <hdf5_readwrite.h>
#include <Tools/block_hash.h>


namespace Loci {

  template<class T> class dstoreRepI : public storeRep {
    block_hash<T>  attrib_data;

    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &en) ;
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en);

    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en) const;
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en) const;

#ifdef H5_HAVE_PARALLEL 
    void  hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &en, hid_t xfer_plist_id ) ;
    void  hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en, hid_t xfer_plist_id);
    
    void  hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en, hid_t xfer_plist_id) const;
    void  hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en, hid_t xfer_plist_id) const;
#endif

    int   get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int   get_estimated_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    void  packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                   const entitySet &e ) ;
    void  unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size,
                     const sequence &seq) ;
    int  get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    int  get_estimated_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    void  packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                   const entitySet &e ) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                    const sequence &seq) ;
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info get_frame_info(IDENTITY_CONVERTER g) ;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g) ;
  public:
    dstoreRepI(){}
    dstoreRepI(const entitySet &p) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual void erase(const entitySet& rm) ;
    virtual void invalidate(const entitySet& valid) ;
    virtual void guarantee_domain(const entitySet& include) ;
    virtual void shift(int_type offset) ;
    virtual ~dstoreRepI() {}
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 MPI_Comm comm LOCI_DEFAULT_COMM) ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 const dMap& remap, MPI_Comm comm LOCI_DEFAULT_COMM) ;
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
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
#ifdef H5_HAVE_PARALLEL 
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset,hid_t xfer_plist_id ) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist_id) const ;
#endif
    
    virtual entitySet domain() const;
    block_hash<T> *get_attrib_data() { return &attrib_data; }
    const block_hash<T> *get_attrib_data() const { return &attrib_data; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
  } ; 

  template<class T> class dstore : public store_instance {
    typedef dstoreRepI<T>  storeType ;
    block_hash<T>       *attrib_data;
    dstore(const dstore &var) { setRep(var.Rep()) ; }
    dstore<T> & operator=(const dstore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef T containerType ;
    dstore() { setRep(new storeType); }
    dstore(storeRepP &rp) { setRep(rp) ; }
    virtual ~dstore() {}
    virtual void notification() ;
    dstore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    entitySet domain() const { return Rep()->domain(); }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

    T &elem(int indx) {
      return attrib_data->access(indx) ;
    }

    const T &elem(int indx) const { return attrib_data->elem(indx) ; }
  
    T &operator[](int indx) { return elem(indx);}
    const T&operator[](int indx) const { return elem(indx);}
    const T&operator()(int indx) const { return elem(indx) ; }
  } ;

  template<class T> class const_dstore : public store_instance {
    typedef dstoreRepI<T> storeType ;
    block_hash<T>      *attrib_data;
    const_dstore(const dstore<T> &var) { setRep(var.Rep()) ; }
    const_dstore(const const_dstore &var) { setRep(var.Rep()) ; }
    const_dstore<T> & operator=(const const_dstore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_dstore<T> & operator=(const dstore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef T containerType ;
    const_dstore() { setRep(new storeType) ; }
    const_dstore(storeRepP &rp) { setRep(rp) ; }
    virtual ~const_dstore() {}
    virtual void notification() ;
    virtual instance_type access() const  ;
    const_dstore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    entitySet domain() const { return Rep()->domain(); }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }

    const T &elem(int indx) const {
      return attrib_data->elem(indx) ;
    } 
    const T&operator[](int indx) const { return elem(indx);}
    const T&operator()(int indx) const { return elem(indx) ; }
  } ;


} // end of namespace Loci
#endif
