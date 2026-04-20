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

#ifndef LOCI_GPUMAPVEC_DEF_H_
#define LOCI_GPUMAPVEC_DEF_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Tools/basic_types.h>
#include <Tools/debug.h>
#include <gpurep.h>
#include <gpuMap.h>

namespace Loci {

  template<unsigned int M>
  class gpuMapVecRepI : public gpuMapRep {
  private:
    entitySet store_domain;
    Entity * base_ptr;
    Entity * access_ptr;

  public:
    gpuMapVecRepI() : base_ptr(0), access_ptr(0) { }
    gpuMapVecRepI(entitySet const & p) { allocate(p); }
    ~gpuMapVecRepI();
    virtual void allocate(entitySet const & ptn);
    virtual storeRep * new_store(entitySet const & p) const;
    virtual storeRep * new_store(entitySet const & p, int const * cnt) const;
    virtual storeRepP remap(dMap const & m) const;
    virtual storeRepP MapRemap(dMap const & dm, dMap const & rm) const;
    virtual void compose(dMap const & m, entitySet const & context);
    virtual void copy(storeRepP & st, entitySet const & context);
    virtual void gather(dMap const & m, storeRepP & st, entitySet const & context);
    virtual void scatter(dMap const & m, storeRepP & st, entitySet const & context);
    virtual int pack_size(entitySet const & e, entitySet & packed);
    virtual int pack_size(entitySet const & e);
    virtual int estimated_pack_size(entitySet const & e);
    virtual void pack(void * ptr, int & loc, int & size, entitySet const & e);
    virtual void unpack(void * ptr, int & loc, int & size, sequence const & seq);

#ifdef DYNAMICSCHEDULING
    virtual void pack(void * ptr, int & loc, int & size, entitySet const & e, Map const & remap) {
      pack(ptr, loc, size, e);
    }

    virtual void unpack(void * ptr, int & loc, int & size, sequence const & seq, dMap const & remap) {
      unpack(ptr, loc, size, seq);
    }
#endif

    virtual entitySet domain() const;

    virtual entitySet image(entitySet const & domain) const;
    virtual std::pair<entitySet, entitySet> preimage(entitySet const & codomain) const;
    virtual storeRepP get_map();
    virtual std::ostream & Print(std::ostream & s) const;
    virtual std::istream & Input(std::istream & s);
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, char const * name, frame_info & fi, entitySet & en);
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, char const * name, entitySet & en) const;

#ifdef H5_HAVE_PARALLEL
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, char const * name, frame_info & fi, entitySet & en, hid_t xfer_plist_id);
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, char const * name, entitySet & en, hid_t xfer_plist_id) const;
#endif

    Entity * get_base_ptr() const {
      return base_ptr;
    }

    Entity * get_access_ptr() const {
      return access_ptr;
    }

    virtual storeRepP expand(entitySet & out_of_dom, std::vector<entitySet> & init_ptn);
    virtual storeRepP freeze();
    virtual storeRepP thaw();
    virtual DatatypeP getType();
    virtual frame_info get_frame_info();

#ifdef DYNAMICSCHEDULING
    virtual storeRepP freeze(entitySet const & es) const {
      cerr << "gpuMapVecRepI::freeze is not implemented yet" << endl;
      Loci::Abort();
      return storeRepP(0);
    }

    virtual storeRepP thaw(entitySet const & es) const {
      cerr << "gpuMapVecRepI::thaw is not implemented yet" << endl;
      Loci::Abort();
      return storeRepP(0);
    }
#endif

    void copyFrom(const storeRepP &fromMap, entitySet set) override;

    store_type RepType() const override;
  } ;

  template<unsigned int M>
  class EntityMapVec {
    Entity * ptr ;

  public:
    GPU_DECL
    EntityMapVec(Entity * ptr) : ptr(ptr) { }

    GPU_DECL
    Entity & operator[](int index) {
      return ptr[index];
    }

    GPU_DECL
    Entity const & operator[](int index) const {
      return ptr[index];
    }

    GPU_DECL
    unsigned int size() const {
      return M;
    }
  };

  template<unsigned int M>
  class const_EntityMapVec {
    Entity const * ptr;

  public:
    GPU_DECL
    const_EntityMapVec(Entity const * ptr) : ptr(ptr) { }

    GPU_DECL
    Entity const & operator[](int index) const {
      return ptr[index];
    }

    GPU_DECL
    unsigned int size() const {
      return M;
    }
  };

  template<unsigned int M>
  class MapVecView {
    Entity * access_ptr;

  public:
    GPU_DECL
    MapVecView() : access_ptr(nullptr) { }

    GPU_DECL
    MapVecView(Entity * access_ptr) : access_ptr(access_ptr) { }

    GPU_DECL
    EntityMapVec<M> operator[](Entity e) {
      return EntityMapVec<M>(&access_ptr[e*M]);
    }

    GPU_DECL
    const_EntityMapVec<M> operator[](Entity e) const {
      return EntityMapVec<M>(&access_ptr[e*M]);
    }
  };

  template<unsigned int M>
  class const_MapVecView {
    Entity const * access_ptr;

  public:
    GPU_DECL
    const_MapVecView() : access_ptr(nullptr) { }

    GPU_DECL
    const_MapVecView(int const * access_ptr) : access_ptr(access_ptr) { }

    GPU_DECL
    const_EntityMapVec<M> operator[](Entity e) const {
      return const_EntityMapVec<M>(&access_ptr[e*M]);
    }
  };

  template<unsigned int M> class const_gpuMapVec ;

  template<unsigned int M> class gpuMapVec : public store_instance {
    friend class const_gpuMapVec<M> ;

    Entity * base_ptr ;
    Entity * access_ptr ;

    gpuMapVec(gpuMapVec<M> const & var) {
      setRep(var.Rep()) ;
    }

    gpuMapVec<M> & operator=(gpuMapVec<M> const & src) {
      setRep(src.Rep()) ;
      return *this ;
    }

  public:
    gpuMapVec() {
      setRep(new gpuMapVecRepI<M>()) ;

    }

    gpuMapVec(storeRepP rp) {
      setRep(rp) ;
    }

    virtual ~gpuMapVec() {
    }

    virtual void notification() ;

    gpuMapVec<M> & operator=(storeRepP p) {
      setRep(p) ;
      return *this ;
    }

    void allocate(entitySet const & ptn) {
      Rep()->allocate(ptn);
    }

    entitySet domain() const {
      return Rep()->domain();
    }

    entitySet image(entitySet const & dom) const {
      cerr << "gpuMapVec::image should not be called" << endl;
      Loci::Abort();
      return EMPTY;
    }

    std::pair<entitySet, entitySet> preimage(entitySet const & codomain) const {
      cerr << "gpuMapVec::preimage should not be called" << endl;
      Loci::Abort();
      return std::pair<entitySet, entitySet>(EMPTY, EMPTY);
    }

    operator gpuMapRepP() {
      gpuMapRepP p(Rep());
      fatal(p == 0);
      return p;
    }

    MapVecView<M> view() {
      return MapVecView<M>(access_ptr);
    }

    const_MapVecView<M> view() const {
      return const_MapVecView<M>(access_ptr);
    }

    std::ostream & Print(std::ostream & s) const {
      cerr << "gpuMapVec::Print should not be called" << endl;
      Loci::Abort();
      return s;
    }

    std::istream & Input(std::istream & s) {
      cerr << "gpuMapVec::Input should not be called" << endl;
      Loci::Abort();
      return s;
    }

    int getRangeKeySpace() const {
      return MapRepP(Rep())->getRangeKeySpace();
    }

    void setRangeKeySpace(int v) {
      MapRepP(Rep())->setRangeKeySpace(v);
    }
  } ;

  template<unsigned int M>
  class const_gpuMapVec : public store_instance {

    Entity * base_ptr ;
    Entity * access_ptr ;

    const_gpuMapVec(const_gpuMapVec<M> & src) {
      setRep(src.Rep());
    }

    const_gpuMapVec(gpuMapVec<M> const & src) {
      setRep(src.Rep());
    }

    const_gpuMapVec<M> & operator=(const_gpuMapVec<M> const & src) {
      setRep(src.Rep());
      return *this;
    }

    const_gpuMapVec<M> & operator=(gpuMapVec<M> const & src) {
      setRep(src.Rep());
      return *this;
    }

  public:
    const_gpuMapVec() {
      setRep(new gpuMapVecRepI<M>());
    }

    const_gpuMapVec(storeRepP rp) {
      setRep(rp);
    }

    virtual ~const_gpuMapVec() {
    }

    virtual void notification();

    virtual instance_type access() const;

    const_gpuMapVec<M> operator=(storeRepP p) {
      setRep(p);
      return *this;
    }

    entitySet const domain() const {
      return Rep()->domain();
    }

    entitySet image(entitySet const & dom) const {
      cerr << "const_gpuMapVec::image should not be called" << endl;
      Loci::Abort();
      return EMPTY;
    }

    std::pair<entitySet, entitySet> preimage(entitySet const & codomain) const {
      cerr << "const_gpuMapVec::preimage should not be called" << endl;
      Loci::Abort();
      return std::pair<entitySet, entitySet>(EMPTY, EMPTY);
    }

    operator gpuMapRepP() {
      gpuMapRepP p(Rep());
      fatal(p == 0);
      return p;
    }

    const_MapVecView<M> view() {
      return const_MapVecView<M>(access_ptr);
    }

    std::ostream & Print(std::ostream & s) const {
      return Rep()->Print(s);
    }

    int getRangeKeySpace() const {
      return MapRepP(Rep())->getRangeKeySpace();
    }
  } ;
}

#endif
