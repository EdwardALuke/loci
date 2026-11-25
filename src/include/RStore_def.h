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
#ifndef RSTORE_DEF_H
#define RSTORE_DEF_H

#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <algorithm>
#include <functional>
#include <vector>
#include <utility>

#include <mpi.h>

#include <Tools/debug.h>
#include <store_rep.h>
#include <Map.h>
#include <Tools/intervalSet.h>
#include <hdf5_readwrite.h>

namespace Loci {

  // Relational store representation.  Data is stored as a sorted
  // vector of (Entity,T) pairs.  This provides logarithmic-time
  // constant access via binary search and supports efficient bulk
  // insert and erase operations, making it suitable for representing
  // entities with 64-bit keys.
  template<class T> class rstoreRepI : public storeRep {
  public:
    typedef std::pair<Entity,T> value_type ;
    typedef std::vector<value_type> container_type ;
    typedef typename container_type::iterator iterator ;
    typedef typename container_type::const_iterator const_iterator ;
  private:
    container_type attrib_data ;
    bool sorted ;

    void ensure_sorted() {
      if(!sorted) {
        std::sort(attrib_data.begin(), attrib_data.end(),
                  [](const value_type &a, const value_type &b)
                  { return a.first < b.first ; }) ;
        sorted = true ;
      }
    }

    void ensure_sorted() const {
      if(!sorted)
        const_cast<rstoreRepI*>(this)->ensure_sorted() ;
    }

    const_iterator find(Entity e) const {
      ensure_sorted() ;
      value_type key(e, T()) ;
      auto it = std::lower_bound(attrib_data.begin(), attrib_data.end(),
                                 key,
                                 [](const value_type &a, const value_type &b)
                                 { return a.first < b.first ; }) ;
      if(it != attrib_data.end() && it->first == e)
        return it ;
      return attrib_data.end() ;
    }

    iterator find_mut(Entity e) {
      ensure_sorted() ;
      value_type key(e, T()) ;
      auto it = std::lower_bound(attrib_data.begin(), attrib_data.end(),
                                 key,
                                 [](const value_type &a, const value_type &b)
                                 { return a.first < b.first ; }) ;
      if(it != attrib_data.end() && it->first == e)
        return it ;
      return attrib_data.end() ;
    }

    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset,
                   hsize_t dimension, const char* name,
                   IDENTITY_CONVERTER c, frame_info &fi, entitySet &en) ;
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset,
                   hsize_t dimension, const char* name,
                   USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en) ;

    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset,
                    hsize_t dimension, const char* name,
                    IDENTITY_CONVERTER c, const entitySet &en) const ;
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset,
                    hsize_t dimension, const char* name,
                    USER_DEFINED_CONVERTER c, const entitySet &en) const ;

#ifdef H5_HAVE_PARALLEL
    void  hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset,
                    hsize_t dimension, const char* name,
                    IDENTITY_CONVERTER c, frame_info &fi, entitySet &en,
                    hid_t xfer_plist_id) ;
    void  hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset,
                    hsize_t dimension, const char* name,
                    USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en,
                    hid_t xfer_plist_id) ;

    void  hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset,
                     hsize_t dimension, const char* name,
                     IDENTITY_CONVERTER c, const entitySet &en,
                     hid_t xfer_plist_id) const ;
    void  hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset,
                     hsize_t dimension, const char* name,
                     USER_DEFINED_CONVERTER c, const entitySet &en,
                     hid_t xfer_plist_id) const ;
#endif

    int   get_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset) ;
    int   get_estimated_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset) ;
    void  packdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int size,
                   const entitySet &e) ;
    void  unpackdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int size,
                     const sequence &seq) ;
    int   get_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset) ;
    int   get_estimated_mpi_size(USER_DEFINED_CONVERTER c,
                                 const entitySet &eset) ;
    void  packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                   const entitySet &e) ;
    void  unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                     const sequence &seq) ;
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info get_frame_info(IDENTITY_CONVERTER g) ;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g) ;

  public:
    rstoreRepI() : sorted(true) {}
    rstoreRepI(const entitySet &p) : sorted(true) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual ~rstoreRepI() {}
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 MPI_Comm comm=MPI_COMM_WORLD) ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD) ;
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
    virtual void unpack(void *ptr, int &loc, int &size,
                        const sequence &seq) ;

    virtual store_type RepType() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset,
                          hsize_t dimension, const char* name,
                          frame_info &fi, entitySet &user_eset) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset,
                           hsize_t dimension, const char* name,
                           entitySet& en) const ;
#ifdef H5_HAVE_PARALLEL
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset,
                           hsize_t dimension, const char* name,
                           frame_info &fi, entitySet &user_eset,
                           hid_t xfer_plist_id) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset,
                            hsize_t dimension, const char* name,
                            entitySet& en, hid_t xfer_plist_id) const ;
#endif

    virtual entitySet domain() const ;
    container_type *get_attrib_data() { return &attrib_data ; }
    const container_type *get_attrib_data() const { return &attrib_data ; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;

    // Bulk insert: merge a sorted or unsorted vector of pairs into the
    // store.  Existing entries for matching entities are overwritten.
    void bulk_insert(const std::vector<value_type> &entries) ;

    // Bulk erase: remove all entries whose entities are in the set.
    void bulk_erase(const entitySet &eset) ;

    // Const element access via binary search (logarithmic).
    // Undefined behaviour if entity is not present.
    const T &elem(Entity e) const {
      const_iterator it = find(e) ;
      return it->second ;
    }

    // Mutable element access — inserts a default-constructed entry if
    // the entity is not present.
    T &access(Entity e) {
      ensure_sorted() ;
      value_type key(e, T()) ;
      auto it = std::lower_bound(attrib_data.begin(), attrib_data.end(),
                                 key,
                                 [](const value_type &a, const value_type &b)
                                 { return a.first < b.first ; }) ;
      if(it != attrib_data.end() && it->first == e)
        return it->second ;
      it = attrib_data.insert(it, key) ;
      return it->second ;
    }

    // Iterators over the sorted (Entity,T) pairs
    iterator begin() { ensure_sorted() ; return attrib_data.begin() ; }
    iterator end() { return attrib_data.end() ; }
    const_iterator begin() const { ensure_sorted() ; return attrib_data.begin() ; }
    const_iterator end() const { return attrib_data.end() ; }
    size_t size() const { return attrib_data.size() ; }
  } ;

  template<class T> class rstore : public store_instance {
    typedef rstoreRepI<T> storeType ;
    storeType *store_ptr ;
    rstore(const rstore &var) { setRep(var.Rep()) ; }
    rstore<T> & operator=(const rstore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef T containerType ;
    rstore() { setRep(new storeType) ; }
    rstore(storeRepP &rp) { setRep(rp) ; }
    virtual ~rstore() {}
    virtual void notification() ;
    rstore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    entitySet domain() const { return Rep()->domain() ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

    // Bulk operations
    void bulk_insert(const std::vector<std::pair<Entity,T>> &entries) {
      store_ptr->bulk_insert(entries) ;
    }
    void bulk_erase(const entitySet &eset) {
      store_ptr->bulk_erase(eset) ;
    }

    // Const access via binary search (logarithmic)
    const T &operator[](Entity indx) const {
      return store_ptr->elem(indx) ;
    }

    const T &elem(Entity indx) const {
      return store_ptr->elem(indx) ;
    }
  } ;

  template<class T> class const_rstore : public store_instance {
    typedef rstoreRepI<T> storeType ;
    storeType *store_ptr ;
    const_rstore(const rstore<T> &var) { setRep(var.Rep()) ; }
    const_rstore(const const_rstore &var) { setRep(var.Rep()) ; }
    const_rstore<T> & operator=(const const_rstore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_rstore<T> & operator=(const rstore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef T containerType ;
    const_rstore() { setRep(new storeType) ; }
    const_rstore(storeRepP &rp) { setRep(rp) ; }
    virtual ~const_rstore() {}
    virtual void notification() ;
    virtual instance_type access() const ;
    const_rstore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    entitySet domain() const { return Rep()->domain() ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }

    const T &elem(Entity indx) const {
      return store_ptr->elem(indx) ;
    }
    const T &operator[](Entity indx) const { return elem(indx) ; }
  } ;

  // Relational multi-store representation.  Data is stored as a
  // stable-sorted vector of (Entity,T) pairs where an entity may
  // appear zero or more times.  The insertion order of values for
  // a given entity is preserved by using a stable sort.  The domain
  // includes entities that were explicitly allocated even if they
  // have zero associated values.
  template<class T> class rmultiStoreRepI : public storeRep {
  public:
    typedef std::pair<Entity,T> value_type ;
    typedef std::vector<value_type> container_type ;
    typedef typename container_type::iterator iterator ;
    typedef typename container_type::const_iterator const_iterator ;
  private:
    container_type attrib_data ;
    entitySet store_domain ;
    bool sorted ;

    void ensure_sorted() {
      if(!sorted) {
        std::stable_sort(attrib_data.begin(), attrib_data.end(),
                         [](const value_type &a, const value_type &b)
                         { return a.first < b.first ; }) ;
        sorted = true ;
      }
    }

    void ensure_sorted() const {
      if(!sorted)
        const_cast<rmultiStoreRepI*>(this)->ensure_sorted() ;
    }

    // Find the range [first,last) of entries for entity e
    std::pair<const_iterator,const_iterator> equal_range(Entity e) const {
      ensure_sorted() ;
      value_type key(e, T()) ;
      auto cmp = [](const value_type &a, const value_type &b)
                 { return a.first < b.first ; } ;
      return std::equal_range(attrib_data.begin(), attrib_data.end(),
                              key, cmp) ;
    }

    std::pair<iterator,iterator> equal_range_mut(Entity e) {
      ensure_sorted() ;
      value_type key(e, T()) ;
      auto cmp = [](const value_type &a, const value_type &b)
                 { return a.first < b.first ; } ;
      return std::equal_range(attrib_data.begin(), attrib_data.end(),
                              key, cmp) ;
    }

    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset,
                   hsize_t dimension, const char* name,
                   IDENTITY_CONVERTER c, frame_info &fi, entitySet &en) ;
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset,
                   hsize_t dimension, const char* name,
                   USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en) ;

    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset,
                    hsize_t dimension, const char* name,
                    IDENTITY_CONVERTER c, const entitySet &en) const ;
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset,
                    hsize_t dimension, const char* name,
                    USER_DEFINED_CONVERTER c, const entitySet &en) const ;

#ifdef H5_HAVE_PARALLEL
    void  hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset,
                    hsize_t dimension, const char* name,
                    IDENTITY_CONVERTER c, frame_info &fi, entitySet &en,
                    hid_t xfer_plist_id) ;
    void  hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset,
                    hsize_t dimension, const char* name,
                    USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &en,
                    hid_t xfer_plist_id) ;

    void  hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset,
                     hsize_t dimension, const char* name,
                     IDENTITY_CONVERTER c, const entitySet &en,
                     hid_t xfer_plist_id) const ;
    void  hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset,
                     hsize_t dimension, const char* name,
                     USER_DEFINED_CONVERTER c, const entitySet &en,
                     hid_t xfer_plist_id) const ;
#endif

    int   get_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset) ;
    int   get_estimated_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset) ;
    void  packdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int size,
                   const entitySet &e) ;
    void  unpackdata(IDENTITY_CONVERTER c, void *ptr, int &loc, int size,
                     const sequence &seq) ;
    int   get_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset) ;
    int   get_estimated_mpi_size(USER_DEFINED_CONVERTER c,
                                 const entitySet &eset) ;
    void  packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                   const entitySet &e) ;
    void  unpackdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                     const sequence &seq) ;
    DatatypeP getType(IDENTITY_CONVERTER g) ;
    DatatypeP getType(USER_DEFINED_CONVERTER g) ;
    frame_info get_frame_info(IDENTITY_CONVERTER g) ;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g) ;

  public:
    rmultiStoreRepI() : sorted(true) {}
    rmultiStoreRepI(const entitySet &p) : sorted(true) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual ~rmultiStoreRepI() {}
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 MPI_Comm comm=MPI_COMM_WORLD) ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD) ;
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
    virtual void unpack(void *ptr, int &loc, int &size,
                        const sequence &seq) ;

    virtual store_type RepType() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset,
                          hsize_t dimension, const char* name,
                          frame_info &fi, entitySet &user_eset) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset,
                           hsize_t dimension, const char* name,
                           entitySet& en) const ;
#ifdef H5_HAVE_PARALLEL
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset,
                           hsize_t dimension, const char* name,
                           frame_info &fi, entitySet &user_eset,
                           hid_t xfer_plist_id) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset,
                            hsize_t dimension, const char* name,
                            entitySet& en, hid_t xfer_plist_id) const ;
#endif

    virtual entitySet domain() const ;
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;

    // Bulk insert: append entries to the store.  Multiple values per
    // entity are allowed and their relative order is preserved.
    void bulk_insert(const std::vector<value_type> &entries) ;

    // Bulk erase: remove all entries whose entities are in the set
    // and remove those entities from the domain.
    void bulk_erase(const entitySet &eset) ;

    // Number of values associated with entity e.
    int num_elems(Entity e) const {
      auto range = equal_range(e) ;
      return static_cast<int>(range.second - range.first) ;
    }

    // Const access to values for entity e.
    // Returns pointer to first element; count from num_elems().
    const T *elem(Entity e) const {
      auto range = equal_range(e) ;
      if(range.first == range.second)
        return 0 ;
      return &(range.first->second) ;
    }

    // Return the number of elements for a given entity
    // (alias for num_elems for consistency)
    int vec_size(Entity e) const { return num_elems(e) ; }

    // Public const access to equal_range for wrapper classes
    std::pair<const_iterator,const_iterator>
    equal_range_public(Entity e) const { return equal_range(e) ; }

    // Iterators over the sorted (Entity,T) pairs
    iterator begin() { ensure_sorted() ; return attrib_data.begin() ; }
    iterator end() { return attrib_data.end() ; }
    const_iterator begin() const { ensure_sorted() ; return attrib_data.begin() ; }
    const_iterator end() const { return attrib_data.end() ; }
    size_t size() const { return attrib_data.size() ; }
  } ;

  // A lightweight const span over contiguous (Entity,T) pairs for one entity.
  // The values are stored as pair<Entity,T>, so elements are at stride
  // sizeof(pair<Entity,T>).  This provides operator[] over the T values.
  template<class T> class rmultiStore_ref {
    typedef std::pair<Entity,T> value_type ;
    const value_type *ptr ;
    int sz ;
  public:
    rmultiStore_ref(const value_type *p, int n) : ptr(p), sz(n) {}
    int size() const { return sz ; }
    const T &operator[](int i) const { return ptr[i].second ; }
    // Provide begin/end for range-based for loops
    class iterator {
      const value_type *p ;
    public:
      iterator(const value_type *pp) : p(pp) {}
      const T &operator*() const { return p->second ; }
      iterator &operator++() { ++p ; return *this ; }
      bool operator!=(const iterator &o) const { return p != o.p ; }
    } ;
    iterator begin() const { return iterator(ptr) ; }
    iterator end() const { return iterator(ptr + sz) ; }
  } ;

  template<class T> class rmultiStore : public store_instance {
    typedef rmultiStoreRepI<T> storeType ;
    storeType *store_ptr ;
    rmultiStore(const rmultiStore &var) { setRep(var.Rep()) ; }
    rmultiStore<T> & operator=(const rmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef T containerType ;
    rmultiStore() { setRep(new storeType) ; }
    rmultiStore(storeRepP &rp) { setRep(rp) ; }
    virtual ~rmultiStore() {}
    virtual void notification() ;
    rmultiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    entitySet domain() const { return Rep()->domain() ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

    void bulk_insert(const std::vector<std::pair<Entity,T>> &entries) {
      store_ptr->bulk_insert(entries) ;
    }
    void bulk_erase(const entitySet &eset) {
      store_ptr->bulk_erase(eset) ;
    }

    // Const access: returns a reference to the span of values for entity
    rmultiStore_ref<T> operator[](Entity indx) const {
      auto range = store_ptr->equal_range_public(indx) ;
      int n = static_cast<int>(range.second - range.first) ;
      const std::pair<Entity,T> *p = (n > 0) ? &(*range.first) : 0 ;
      return rmultiStore_ref<T>(p, n) ;
    }

    int num_elems(Entity indx) const {
      return store_ptr->num_elems(indx) ;
    }

    int vec_size(Entity indx) const {
      return store_ptr->vec_size(indx) ;
    }
  } ;

  template<class T> class const_rmultiStore : public store_instance {
    typedef rmultiStoreRepI<T> storeType ;
    storeType *store_ptr ;
    const_rmultiStore(const rmultiStore<T> &var) { setRep(var.Rep()) ; }
    const_rmultiStore(const const_rmultiStore &var) { setRep(var.Rep()) ; }
    const_rmultiStore<T> & operator=(const const_rmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_rmultiStore<T> & operator=(const rmultiStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef T containerType ;
    const_rmultiStore() { setRep(new storeType) ; }
    const_rmultiStore(storeRepP &rp) { setRep(rp) ; }
    virtual ~const_rmultiStore() {}
    virtual void notification() ;
    virtual instance_type access() const ;
    const_rmultiStore<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }
    entitySet domain() const { return Rep()->domain() ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }

    rmultiStore_ref<T> operator[](Entity indx) const {
      auto range = store_ptr->equal_range_public(indx) ;
      int n = static_cast<int>(range.second - range.first) ;
      const std::pair<Entity,T> *p = (n > 0) ? &(*range.first) : 0 ;
      return rmultiStore_ref<T>(p, n) ;
    }

    int num_elems(Entity indx) const {
      return store_ptr->num_elems(indx) ;
    }

    int vec_size(Entity indx) const {
      return store_ptr->vec_size(indx) ;
    }
  } ;

} // end of namespace Loci
#endif
