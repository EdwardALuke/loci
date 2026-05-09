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
#ifndef DSTOREVEC_DEF_H
#define DSTOREVEC_DEF_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <Tools/debug.h>
#include <Tools/tools.h>
#include <hdf5_readwrite.h>
#include <store_rep.h>
#include <Tools/lmutex.h>
#include <DMap.h>
#include <vector>
#include <algorithm>
#include <Tools/hash_map.h>

namespace Loci {

  template<class T> class dstoreVecRepI : public storeRep {

    lmutex                    mutex ;
    entitySet                 store_domain ;
    int                       size;
    HASH_MAP(int,std::vector<T>)  attrib_data;
    bool         isMat; //if this is a storeMat
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &en) ;
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &usr) ;

    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en) const;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en) const;

    void hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &en,hid_t xfer_plist_id ) ;
    void hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &usr, hid_t xfer_plist_id) ;

    void hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en, hid_t xfer_plist_id) const;
    void hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en, hid_t xfer_plist_id) const;

    int get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    int get_estimated_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_estimated_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size, const entitySet &e ) ;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size, const entitySet &e ) ;

    void unpackdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int &size, const sequence &seq) ;
    void unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int &size, const sequence &seq) ;
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info get_frame_info(IDENTITY_CONVERTER g) ;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g) ;
  public:
    dstoreVecRepI() {isMat=false;}

    dstoreVecRepI(const entitySet &p) { allocate(p) ; isMat=false; }

    virtual ~dstoreVecRepI() ;
    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
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
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &eset) const ;

    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset, hid_t xfer_plist_id ) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &eset, hid_t xfer_plist_id) const ;

    virtual void set_elem_size(int sz) ;
    
    HASH_MAP(int,std::vector<T>) *get_attrib_data() { return &attrib_data; }
    int get_size() const { return size; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
    void setIsMat(bool im){isMat=im;}
  } ; 

  template<class T> class dstoreVec : public store_instance {
    typedef dstoreVecRepI<T>     storeType ;
    HASH_MAP(int, std::vector<T>) *attrib_data;
    int                          vsize;
    dstoreVec(const dstoreVec<T> &var) {setRep(var.Rep()) ;}
    dstoreVec<T> & operator=(const dstoreVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef std::vector<T> containerType ;
    dstoreVec() {setRep(new storeType) ;}
    dstoreVec(storeRepP &rp) { setRep(rp) ;}
    virtual ~dstoreVec() {}
    virtual void notification() ;
    dstoreVec<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    void setVecSize(int size) {
      Rep()->set_elem_size(size) ;
      vsize = size;
    }
    int getVecSize() {
      return(vsize);
    }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    entitySet domain() const { return Rep()->domain() ; }
    std::vector<T> &elem(int indx) {
      return( (*attrib_data)[indx] );
    }
    std::vector<T> &operator[](int indx) {
      (*attrib_data)[indx].resize(vsize);
      return( elem(indx) );
    }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
  } ;

  template<class T> class const_dstoreVec : public store_instance {
    typedef dstoreVecRepI<T> storeType ;
    HASH_MAP(int, std::vector<T>)  *attrib_data;
    int size ;
    const_dstoreVec(const const_dstoreVec<T> &var) {setRep(var.Rep()) ;}
    const_dstoreVec(const dstoreVec<T> &var) {setRep(var.Rep()) ;}
    const_dstoreVec<T> & operator=(const dstoreVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_dstoreVec<T> & operator=(const const_dstoreVec<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef std::vector<T> containerType ;
    const_dstoreVec() { setRep(new storeType) ; }
    const_dstoreVec(storeRepP &rp) { setRep(rp) ; }
    virtual ~const_dstoreVec() {}
    virtual void notification() ;
    virtual instance_type access() const ;
    const_dstoreVec<T> & operator=(storeRepP p) {
      setRep(p) ;
      return *this ;
    }
    entitySet domain() const { return Rep()->domain() ; }
    std::vector<T> elem(int indx) const {
      std::vector<T>   newVec;
      typename HASH_MAP(int, std::vector<T>)::const_iterator ci;
      ci = attrib_data->find(indx);
      if( ci != attrib_data->end()){
        newVec = ci->second;
        return( newVec );
      } else {
        std::cerr << "Error: DStoreVec, Out range accessing Vector " << std::endl;
      }
      return newVec;
    }
    int vecSize() const { return size ; }
    std::vector<T>  operator[](int indx) const {
      return( elem(indx) );
    }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
  } ;

} // end of namespace Loci

#endif
