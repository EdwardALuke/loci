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
#ifndef DMULTISTORE_DEF_H
#define DMULTISTORE_DEF_H

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
#include <storeVec.h>
#include <DMap.h>
#include <DMultiMap.h>
#include <vector>
#include <algorithm>
#include <Tools/hash_map.h>

namespace Loci {
  template<class T> class dmultiStoreRepI : public storeRep {
    entitySet                 store_domain ;
    HASH_MAP(int,std::vector<T>)  attrib_data;

    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &usr) ;
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &usr);
    
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en) const;
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en) const;

    void  hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &usr, hid_t xfer_plist_id) ;
    void  hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &usr,hid_t xfer_plist_id);
    
    void  hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &en, hid_t xfer_plist_id) const;
    void  hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &en, hid_t xfer_plist_id) const;
    
    int get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_estimated_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset);
    int get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset);
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

    //  Constructor ...
    dmultiStoreRepI() { }

    dmultiStoreRepI(const entitySet &p) { store_domain=p;}

    dmultiStoreRepI(const store<int> &sizes) { allocate(sizes) ; }

    //  Destructors ...
    virtual ~dmultiStoreRepI() ;

    //  Member function ...
    void allocate(const store<int> &sizes) ;
    virtual void shift(int_type offset) ;
    virtual void allocate(const entitySet &ptn) ;
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
    virtual int pack_size(const entitySet &e) ;
    virtual int estimated_pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size,  const sequence &seq ) ;
    		      
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;

    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset, hid_t xfer_plist) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist) const ;

    

    HASH_MAP(int,std::vector<T>) *get_attrib_data(){return &attrib_data; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
  } ;

  template<class T> class dmultiStore : public store_instance {
    typedef dmultiStoreRepI<T>  storeType ;
    HASH_MAP(int, std::vector<T> ) *attrib_data;
    dmultiStore(const dmultiStore<T> &var) {setRep(var.Rep()) ;}
    dmultiStore<T> & operator=(const dmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef std::vector<T> containerType ;
    dmultiStore() {setRep(new storeType) ;}
    dmultiStore(storeRepP &rp) { setRep(rp) ;}
    virtual ~dmultiStore() {}
    virtual void notification() ;
    dmultiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    void allocate(const store<int> &sizes) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->allocate(sizes) ;
    }
    void setSizes(const const_dmultiMap &m) {
      NPTR<storeType> p(Rep()) ;
      fatal(p==0) ;
      p->setSizes(m) ;
    }
    entitySet domain() const { return Rep()->domain() ; }
    std::vector<T> &elem(int indx) {
      return( (*attrib_data)[indx] );
    }
    std::vector<T> &operator[](int indx) {
      return( (*attrib_data)[indx] );
    }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
  } ;

  template<class T> class const_dmultiStore : public store_instance {
    typedef dmultiStoreRepI<T> storeType ;
    HASH_MAP(int, std::vector<T> )   *attrib_data;
    const_dmultiStore(const const_dmultiStore<T> &var) {setRep(var.Rep()) ;}
    const_dmultiStore(const dmultiStore<T> &var) {setRep(var.Rep()) ;}
    const_dmultiStore<T> & operator=(const dmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_dmultiStore<T> & operator=(const const_dmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    //    typedef const_Vect<T> containerType ;
    const_dmultiStore() {setRep(new storeType) ;}
    const_dmultiStore(storeRepP &rp) { setRep(rp) ;}
    virtual ~const_dmultiStore() {}
    virtual void notification() ;
    virtual instance_type access() const ;
    const_dmultiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    const entitySet domain() const { return Rep()->domain() ; }
    std::vector<T> &elem(int indx) {
      return( (*attrib_data)[indx] );
    }
    std::vector<T> &operator[](int indx) {
      return( (*attrib_data)[indx] );
    }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
  } ;

} // end of namespace Loci

#endif
