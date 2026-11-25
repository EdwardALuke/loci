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
#ifndef RMULTIMAP_DEF_H
#define RMULTIMAP_DEF_H

#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <algorithm>
#include <vector>
#include <utility>

#include <mpi.h>

#include <Tools/debug.h>
#include <Map_rep.h>
#include <Tools/intervalSet.h>
#include <hdf5_readwrite.h>

namespace Loci {

  // Relational multi-map representation.  Stores a one-to-many mapping
  // from Entity to Entity in a stable-sorted vector of pairs.  An
  // entity may map to zero or more target entities.  Insertion order
  // is preserved via stable sort.  A separate domain set tracks
  // entities with zero mappings.
  class rmultiMapRepI : public MapRep {
  public:
    typedef std::pair<Entity,Entity> value_type ;
    typedef std::vector<value_type> container_type ;
    typedef container_type::iterator iterator ;
    typedef container_type::const_iterator const_iterator ;
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
        const_cast<rmultiMapRepI*>(this)->ensure_sorted() ;
    }

    std::pair<const_iterator,const_iterator> equal_range(Entity e) const {
      ensure_sorted() ;
      value_type key(e, 0) ;
      auto cmp = [](const value_type &a, const value_type &b)
                 { return a.first < b.first ; } ;
      return std::equal_range(attrib_data.begin(), attrib_data.end(),
                              key, cmp) ;
    }

    std::pair<iterator,iterator> equal_range_mut(Entity e) {
      ensure_sorted() ;
      value_type key(e, 0) ;
      auto cmp = [](const value_type &a, const value_type &b)
                 { return a.first < b.first ; } ;
      return std::equal_range(attrib_data.begin(), attrib_data.end(),
                              key, cmp) ;
    }

  public:
    rmultiMapRepI() : sorted(true) {}
    rmultiMapRepI(const entitySet &p) : sorted(true) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual ~rmultiMapRepI() {}
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP MapRemap(const dMap &dm, const dMap &rm) const ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 MPI_Comm comm=MPI_COMM_WORLD) ;
    virtual storeRepP
    redistribute(const std::vector<entitySet>& dom_ptn,
                 const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD) ;
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
    virtual void unpack(void *ptr, int &loc, int &size,
                        const sequence &seq) ;

    virtual entitySet domain() const ;
    virtual entitySet image(const entitySet &domain) const ;
    virtual std::pair<entitySet,entitySet>
    preimage(const entitySet &codomain) const ;
    virtual storeRepP get_map() ;
    virtual storeRepP expand(entitySet &out_of_dom,
                             std::vector<entitySet> &init_ptn) ;
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
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;

    // Bulk insert: append entries (multiple targets per entity allowed).
    void bulk_insert(const std::vector<value_type> &entries) ;

    // Bulk erase: remove all entries for entities in eset and remove
    // those entities from the domain.
    void bulk_erase(const entitySet &eset) ;

    // Number of targets for entity e.
    int num_elems(Entity e) const {
      auto range = equal_range(e) ;
      return static_cast<int>(range.second - range.first) ;
    }

    int vec_size(Entity e) const { return num_elems(e) ; }

    // Public const equal_range for wrapper classes
    std::pair<const_iterator,const_iterator>
    equal_range_public(Entity e) const { return equal_range(e) ; }

    // Iterators
    iterator begin() { ensure_sorted() ; return attrib_data.begin() ; }
    iterator end() { return attrib_data.end() ; }
    const_iterator begin() const { ensure_sorted() ; return attrib_data.begin() ; }
    const_iterator end() const { return attrib_data.end() ; }
    size_t size() const { return attrib_data.size() ; }
  } ;

  // Lightweight const span over the targets of one entity in an rmultiMap
  class rmultiMap_ref {
    typedef std::pair<Entity,Entity> value_type ;
    const value_type *ptr ;
    int sz ;
  public:
    rmultiMap_ref(const value_type *p, int n) : ptr(p), sz(n) {}
    int size() const { return sz ; }
    Entity operator[](int i) const { return ptr[i].second ; }
    class iterator {
      const value_type *p ;
    public:
      iterator(const value_type *pp) : p(pp) {}
      Entity operator*() const { return p->second ; }
      iterator &operator++() { ++p ; return *this ; }
      bool operator!=(const iterator &o) const { return p != o.p ; }
    } ;
    iterator begin() const { return iterator(ptr) ; }
    iterator end() const { return iterator(ptr + sz) ; }
  } ;

  class rmultiMap : public store_instance {
    typedef rmultiMapRepI MapType ;
    MapType *store_ptr ;
    rmultiMap(const rmultiMap &var) { setRep(var.Rep()) ; }
    rmultiMap & operator=(const rmultiMap &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    rmultiMap() { setRep(new MapType) ; }
    rmultiMap(storeRepP &rp) { setRep(rp) ; }
    virtual ~rmultiMap() {}
    virtual void notification() ;
    rmultiMap & operator=(storeRepP p) { setRep(p) ; return *this ; }
    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    entitySet domain() const { return Rep()->domain() ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ;
    }

    void bulk_insert(const std::vector<std::pair<Entity,Entity>> &entries) {
      store_ptr->bulk_insert(entries) ;
    }
    void bulk_erase(const entitySet &eset) {
      store_ptr->bulk_erase(eset) ;
    }

    rmultiMap_ref operator[](Entity indx) const {
      auto range = store_ptr->equal_range_public(indx) ;
      int n = static_cast<int>(range.second - range.first) ;
      const std::pair<Entity,Entity> *p = (n > 0) ? &(*range.first) : 0 ;
      return rmultiMap_ref(p, n) ;
    }

    int num_elems(Entity indx) const {
      return store_ptr->num_elems(indx) ;
    }

    int vec_size(Entity indx) const {
      return store_ptr->vec_size(indx) ;
    }

    entitySet image(const entitySet &dom) const {
      return MapRepP(Rep())->image(dom) ;
    }
    std::pair<entitySet,entitySet> preimage(const entitySet &codomain) const {
      return MapRepP(Rep())->preimage(codomain) ;
    }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const rmultiMap &m)
  { return m.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, rmultiMap &m)
  { return m.Input(s) ; }

  class const_rmultiMap : public store_instance {
    typedef rmultiMapRepI MapType ;
    MapType *store_ptr ;
    const_rmultiMap(const rmultiMap &var) { setRep(var.Rep()) ; }
    const_rmultiMap(const const_rmultiMap &var) { setRep(var.Rep()) ; }
    const_rmultiMap & operator=(const const_rmultiMap &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_rmultiMap & operator=(const rmultiMap &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    const_rmultiMap() { setRep(new MapType) ; }
    const_rmultiMap(storeRepP &rp) { setRep(rp) ; }
    virtual ~const_rmultiMap() {}
    virtual void notification() ;
    virtual instance_type access() const ;
    const_rmultiMap & operator=(storeRepP p) { setRep(p) ; return *this ; }
    entitySet domain() const { return Rep()->domain() ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ;
    }

    rmultiMap_ref operator[](Entity indx) const {
      auto range = store_ptr->equal_range_public(indx) ;
      int n = static_cast<int>(range.second - range.first) ;
      const std::pair<Entity,Entity> *p = (n > 0) ? &(*range.first) : 0 ;
      return rmultiMap_ref(p, n) ;
    }

    int num_elems(Entity indx) const {
      return store_ptr->num_elems(indx) ;
    }

    int vec_size(Entity indx) const {
      return store_ptr->vec_size(indx) ;
    }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const const_rmultiMap &m)
  { return m.Print(s) ; }

} // end of namespace Loci
#endif
