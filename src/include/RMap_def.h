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
#ifndef RMAP_DEF_H
#define RMAP_DEF_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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

  // Relational map representation.  Stores a one-to-one mapping from
  // Entity to Entity in a sorted vector of pairs with logarithmic
  // access time via binary search.  Analogous to rstore<Entity> but
  // extends MapRep to provide image/preimage/compose operations.
  class rMapRepI : public MapRep {
  public:
    typedef std::pair<Entity,Entity> value_type ;
    typedef std::vector<value_type> container_type ;
    typedef container_type::iterator iterator ;
    typedef container_type::const_iterator const_iterator ;
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
        const_cast<rMapRepI*>(this)->ensure_sorted() ;
    }

    const_iterator find(Entity e) const {
      ensure_sorted() ;
      value_type key(e, 0) ;
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
      value_type key(e, 0) ;
      auto it = std::lower_bound(attrib_data.begin(), attrib_data.end(),
                                 key,
                                 [](const value_type &a, const value_type &b)
                                 { return a.first < b.first ; }) ;
      if(it != attrib_data.end() && it->first == e)
        return it ;
      return attrib_data.end() ;
    }

  public:
    rMapRepI() : sorted(true) {}
    rMapRepI(const entitySet &p) : sorted(true) { allocate(p) ; }
    virtual void allocate(const entitySet &ptn) ;
    virtual void shift(int_type offset) ;
    virtual ~rMapRepI() {}
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

    // Bulk insert: merge new entries; existing entities are overwritten.
    void bulk_insert(const std::vector<value_type> &entries) ;

    // Bulk erase: remove entries for entities in eset.
    void bulk_erase(const entitySet &eset) ;

    // Const element access via binary search.
    Entity elem(Entity e) const {
      const_iterator it = find(e) ;
      return it->second ;
    }

    // Mutable element access — inserts 0 if entity not present.
    Entity &access(Entity e) {
      ensure_sorted() ;
      value_type key(e, 0) ;
      auto it = std::lower_bound(attrib_data.begin(), attrib_data.end(),
                                 key,
                                 [](const value_type &a, const value_type &b)
                                 { return a.first < b.first ; }) ;
      if(it != attrib_data.end() && it->first == e)
        return it->second ;
      it = attrib_data.insert(it, key) ;
      return it->second ;
    }

    // Iterators
    iterator begin() { ensure_sorted() ; return attrib_data.begin() ; }
    iterator end() { return attrib_data.end() ; }
    const_iterator begin() const { ensure_sorted() ; return attrib_data.begin() ; }
    const_iterator end() const { return attrib_data.end() ; }
    size_t size() const { return attrib_data.size() ; }
  } ;

  class rMap : public store_instance {
    typedef rMapRepI MapType ;
    MapType *store_ptr ;
    rMap(const rMap &var) { setRep(var.Rep()) ; }
    rMap & operator=(const rMap &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    rMap() { setRep(new MapType) ; }
    rMap(storeRepP &rp) { setRep(rp) ; }
    virtual ~rMap() {}
    virtual void notification() ;
    rMap & operator=(storeRepP p) { setRep(p) ; return *this ; }
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

    // Const access via binary search
    Entity operator[](Entity indx) const {
      return store_ptr->elem(indx) ;
    }

    Entity elem(Entity indx) const {
      return store_ptr->elem(indx) ;
    }

    entitySet image(const entitySet &dom) const {
      return MapRepP(Rep())->image(dom) ;
    }
    std::pair<entitySet,entitySet> preimage(const entitySet &codomain) const {
      return MapRepP(Rep())->preimage(codomain) ;
    }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const rMap &m)
  { return m.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, rMap &m)
  { return m.Input(s) ; }

  class const_rMap : public store_instance {
    typedef rMapRepI MapType ;
    MapType *store_ptr ;
    const_rMap(const rMap &var) { setRep(var.Rep()) ; }
    const_rMap(const const_rMap &var) { setRep(var.Rep()) ; }
    const_rMap & operator=(const const_rMap &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
    const_rMap & operator=(const rMap &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    const_rMap() { setRep(new MapType) ; }
    const_rMap(storeRepP &rp) { setRep(rp) ; }
    virtual ~const_rMap() {}
    virtual void notification() ;
    virtual instance_type access() const ;
    const_rMap & operator=(storeRepP p) { setRep(p) ; return *this ; }
    entitySet domain() const { return Rep()->domain() ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }

    operator MapRepP() {
      MapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ;
    }

    Entity operator[](Entity indx) const {
      return store_ptr->elem(indx) ;
    }
    Entity elem(Entity indx) const {
      return store_ptr->elem(indx) ;
    }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const const_rMap &m)
  { return m.Print(s) ; }

} // end of namespace Loci
#endif
