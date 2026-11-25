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
#ifndef RMAPVEC_IMPL_H
#define RMAPVEC_IMPL_H

#include <RMapVec_def.h>
#include <MapVec.h>
#include <multiMap.h>
#include <distribute.h>

namespace Loci {

  //------------------------------------------------------------------------
  // bulk_insert
  //------------------------------------------------------------------------
  template<int M>
  void rMapVecRepI<M>::bulk_insert(const std::vector<value_type> &entries) {
    if(entries.empty())
      return ;

    std::vector<value_type> sorted_new(entries) ;
    std::sort(sorted_new.begin(), sorted_new.end(),
              [](const value_type &a, const value_type &b)
              { return a.first < b.first ; }) ;

    ensure_sorted() ;

    container_type merged ;
    merged.reserve(attrib_data.size() + sorted_new.size()) ;

    auto oi = attrib_data.begin() ;
    auto ni = sorted_new.begin() ;

    while(oi != attrib_data.end() && ni != sorted_new.end()) {
      if(oi->first < ni->first) {
        merged.push_back(*oi) ;
        ++oi ;
      } else if(ni->first < oi->first) {
        merged.push_back(*ni) ;
        ++ni ;
      } else {
        merged.push_back(*ni) ;
        ++oi ;
        ++ni ;
      }
    }
    while(oi != attrib_data.end()) {
      merged.push_back(*oi) ;
      ++oi ;
    }
    while(ni != sorted_new.end()) {
      merged.push_back(*ni) ;
      ++ni ;
    }

    attrib_data.swap(merged) ;
    sorted = true ;

    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  // bulk_erase
  //------------------------------------------------------------------------
  template<int M>
  void rMapVecRepI<M>::bulk_erase(const entitySet &eset) {
    if(eset == EMPTY || attrib_data.empty())
      return ;

    ensure_sorted() ;

    container_type kept ;
    kept.reserve(attrib_data.size()) ;

    auto ei = eset.begin() ;
    for(auto &p : attrib_data) {
      while(ei != eset.end() && *ei < p.first)
        ++ei ;
      if(ei != eset.end() && *ei == p.first)
        continue ;
      kept.push_back(p) ;
    }

    attrib_data.swap(kept) ;
    sorted = true ;

    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  // allocate
  //------------------------------------------------------------------------
  template<int M>
  void rMapVecRepI<M>::allocate(const entitySet &ptn) {
    entitySet dom = domain() ;
    entitySet remove = dom - ptn ;
    entitySet add = ptn - dom ;

    bulk_erase(remove) ;

    if(add != EMPTY) {
      std::vector<value_type> new_entries ;
      new_entries.reserve(add.size()) ;
      VEC zero_vec(0) ;
      for(auto ei = add.begin(); ei != add.end(); ++ei) {
        value_type v ;
        v.first = *ei ;
        v.second = zero_vec ;
        new_entries.push_back(v) ;
      }
      bulk_insert(new_entries) ;
    }

    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  template<int M>
  void rMapVecRepI<M>::shift(int_type offset) {
    for(auto &p : attrib_data)
      p.first += offset ;
    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  template<int M>
  entitySet rMapVecRepI<M>::domain() const {
    ensure_sorted() ;
    entitySet dom ;
    for(const auto &p : attrib_data)
      dom += p.first ;
    return dom ;
  }

  //------------------------------------------------------------------------
  template<int M>
  entitySet rMapVecRepI<M>::image(const entitySet &iset) const {
    ensure_sorted() ;
    entitySet dom = domain() ;
    entitySet common = dom & iset ;
    entitySet img ;
    for(auto ei = common.begin(); ei != common.end(); ++ei) {
      auto it = find(*ei) ;
      if(it != attrib_data.end()) {
        for(int j=0;j<M;++j)
          img += it->second[j] ;
      }
    }
    return img ;
  }

  //------------------------------------------------------------------------
  template<int M>
  std::pair<entitySet,entitySet>
  rMapVecRepI<M>::preimage(const entitySet &codomain) const {
    ensure_sorted() ;
    entitySet domaini ;
    entitySet domainu ;
    for(const auto &p : attrib_data) {
      bool vali = true ;
      bool valu = false ;
      for(int j=0;j<M;++j) {
        bool in_set = codomain.inSet(p.second[j]) ;
        vali = vali && in_set ;
        valu = valu || in_set ;
      }
      if(vali)
        domaini += p.first ;
      if(valu)
        domainu += p.first ;
    }
    return std::make_pair(domaini, domainu) ;
  }

  //------------------------------------------------------------------------
  template<int M>
  storeRepP rMapVecRepI<M>::get_map() {
    return getRep() ;
  }

  //------------------------------------------------------------------------
  template<int M>
  void rMapVecRepI<M>::compose(const dMap &m, const entitySet &context) {
    ensure_sorted() ;
    entitySet common = domain() & context ;
    entitySet mdom = m.domain() ;
    FORALL(common, i) {
      auto it = find_mut(i) ;
      if(it != attrib_data.end()) {
        for(int j=0;j<M;++j) {
          if(mdom.inSet(it->second[j]))
            it->second[j] = m[it->second[j]] ;
          else
            it->second[j] = -1 ;
        }
      }
    } ENDFORALL ;
  }

  //------------------------------------------------------------------------
  template<int M>
  storeRepP rMapVecRepI<M>::expand(entitySet &out_of_dom,
                                   std::vector<entitySet> &init_ptn) {
    return getRep() ;
  }

  //------------------------------------------------------------------------
  // Print / Input
  //------------------------------------------------------------------------
  template<int M>
  std::ostream &rMapVecRepI<M>::Print(std::ostream &s) const {
    entitySet dom = domain() ;
    s << '{' << dom << std::endl ;
    ensure_sorted() ;
    for(const auto &p : attrib_data) {
      for(int j=0;j<M;++j)
        s << p.second[j] << ' ' ;
      s << std::endl ;
    }
    s << '}' << std::endl ;
    return s ;
  }

  template<int M>
  std::istream &rMapVecRepI<M>::Input(std::istream &s) {
    entitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '{') {
      std::cerr << "Incorrect Format while reading rMapVec" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    s >> e ;

    std::vector<value_type> entries ;
    entries.reserve(e.size()) ;
    FORALL(e, ii) {
      value_type v ;
      v.first = ii ;
      for(int j=0;j<M;++j)
        s >> v.second[j] ;
      entries.push_back(v) ;
    } ENDFORALL ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '}') {
      std::cerr << "Incorrect Format while reading rMapVec" << std::endl ;
      s.putback(ch) ;
    }

    attrib_data.clear() ;
    bulk_insert(entries) ;

    return s ;
  }

  //------------------------------------------------------------------------
  template<int M>
  storeRep *rMapVecRepI<M>::new_store(const entitySet &p) const {
    return new rMapVecRepI<M>(p) ;
  }

  template<int M>
  storeRep *rMapVecRepI<M>::new_store(const entitySet &p,
                                      const int* cnt) const {
    storeRep* sp = 0 ;
    std::cerr << " This method should not be called for an rMapVec "
              << std::endl ;
    return sp ;
  }

  //------------------------------------------------------------------------
  // remap
  //------------------------------------------------------------------------
  template<int M>
  storeRepP rMapVecRepI<M>::remap(const dMap &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    rMapVec<M> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m, my_store, newdomain) ;
    return s.Rep() ;
  }

  //------------------------------------------------------------------------
  template<int M>
  storeRepP rMapVecRepI<M>::MapRemap(const dMap &dm, const dMap &rm) const {
    entitySet newdomain = dm.domain() & domain() ;
    rMapVec<M> s ;
    s.allocate(dm.image(newdomain)) ;

    ensure_sorted() ;
    std::vector<value_type> entries ;
    FORALL(newdomain, i) {
      auto it = find(i) ;
      if(it != attrib_data.end()) {
        value_type v ;
        v.first = dm[i] ;
        for(int j=0;j<M;++j)
          v.second[j] = rm[it->second[j]] ;
        entries.push_back(v) ;
      }
    } ENDFORALL ;
    s.bulk_insert(entries) ;
    return s.Rep() ;
  }

  //------------------------------------------------------------------------
  // redistribute
  //------------------------------------------------------------------------
  template<int M>
  storeRepP
  rMapVecRepI<M>::redistribute(const std::vector<entitySet>& dom_ptn,
                               MPI_Comm comm) {
    entitySet dom = domain() ;
    int np ;
    MPI_Comm_size(comm, &np) ;
    fatal(size_t(np) != dom_ptn.size()) ;
    entitySet total_dom ;
    std::vector<entitySet> dom_split(np) ;
    for(int i=0;i<np;++i) {
      dom_split[i] = dom & dom_ptn[i] ;
      total_dom += dom_ptn[i] ;
    }

    std::vector<int> send_counts(np,0) ;
    for(int i=0;i<np;++i)
      send_counts[i] = pack_size(dom_split[i]) ;

    std::vector<int> send_displs(np,0) ;
    for(int i=1;i<np;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    std::vector<int> recv_counts(np,0) ;
    MPI_Alltoall(send_counts.data(), 1, MPI_INT,
                 recv_counts.data(), 1, MPI_INT, comm) ;

    std::vector<int> recv_displs(np,0) ;
    for(int i=1;i<np;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    int tot_send_size = 0 ;
    for(int i=0;i<np;++i)
      tot_send_size += send_counts[i] ;

    std::vector<unsigned char> send_buffer(tot_send_size) ;
    std::vector<unsigned char*> send_ptr(np) ;
    send_ptr[0] = send_buffer.data() ;
    for(int i=1;i<np;++i)
      send_ptr[i] = send_ptr[i-1] + send_counts[i-1] ;

    for(int i=0;i<np;++i) {
      int position = 0 ;
      pack(send_ptr[i], position, send_counts[i], dom_split[i]) ;
    }

    std::vector<sequence> pack_seq(np) ;
    for(int i=0;i<np;++i)
      pack_seq[i] = sequence(dom_split[i]) ;

    std::vector<sequence> unpack_seq
      = transpose_sequence(pack_seq, comm) ;
    std::vector<sequence>().swap(pack_seq) ;

    int tot_recv_size = 0 ;
    for(int i=0;i<np;++i)
      tot_recv_size += recv_counts[i] ;

    std::vector<unsigned char> recv_buffer(tot_recv_size) ;
    MPI_Alltoallv(send_buffer.data(), send_counts.data(),
                  send_displs.data(), MPI_PACKED,
                  recv_buffer.data(), recv_counts.data(),
                  recv_displs.data(), MPI_PACKED, comm) ;

    rMapVec<M> ns ;
    entitySet new_domain ;
    for(int i=0;i<np;++i)
      new_domain += entitySet(unpack_seq[i]) ;
    entitySet old_dom = dom - total_dom ;
    ns.allocate(new_domain+old_dom) ;
    storeRepP srp = ns.Rep() ;

    std::vector<unsigned char*> recv_ptr(np) ;
    recv_ptr[0] = recv_buffer.data() ;
    for(int i=1;i<np;++i)
      recv_ptr[i] = recv_ptr[i-1] + recv_counts[i-1] ;

    for(int i=0;i<np;++i) {
      int position = 0 ;
      srp->unpack(recv_ptr[i], position, recv_counts[i], unpack_seq[i]) ;
    }

    ensure_sorted() ;
    for(entitySet::const_iterator ei=old_dom.begin();
        ei!=old_dom.end();++ei) {
      auto it = find(*ei) ;
      if(it != attrib_data.end()) {
        NPTR<rMapVecRepI<M>> nsp(srp) ;
        nsp->access(*ei) = it->second ;
      }
    }

    return srp ;
  }

  template<int M>
  storeRepP
  rMapVecRepI<M>::redistribute(const std::vector<entitySet>& dom_ptn,
                               const dMap& remap, MPI_Comm comm) {
    rMapVec<M> new_store ;
    Loci::Abort() ;
    return new_store.Rep() ;
  }

  //------------------------------------------------------------------------
  // freeze / thaw
  //------------------------------------------------------------------------
  template<int M>
  storeRepP rMapVecRepI<M>::freeze() {
    MapVec<M> static_map ;
    entitySet dom = domain() ;
    static_map.allocate(dom) ;
    ensure_sorted() ;
    for(const auto &p : attrib_data)
      static_map[p.first] = p.second ;
    return static_map.Rep() ;
  }

  template<int M>
  storeRepP rMapVecRepI<M>::thaw() {
    return getRep() ;
  }

  //------------------------------------------------------------------------
  // copy / gather / scatter
  //------------------------------------------------------------------------
  template<int M>
  void rMapVecRepI<M>::copy(storeRepP &st, const entitySet &context) {
    const_rMapVec<M> s(st) ;
    fatal((context-domain()) != EMPTY) ;
    FORALL(context, i) {
      access(i) = s[i] ;
    } ENDFORALL ;
  }

  template<int M>
  void rMapVecRepI<M>::gather(const dMap &m, storeRepP &st,
                              const entitySet &context) {
    const_rMapVec<M> s(st) ;
    fatal(context != EMPTY) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context, i) {
      access(i) = s[m[i]] ;
    } ENDFORALL ;
  }

  template<int M>
  void rMapVecRepI<M>::scatter(const dMap &m, storeRepP &st,
                               const entitySet &context) {
    const_rMapVec<M> s(st) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY) ;
    FORALL(context, i) {
      access(m[i]) = s[i] ;
    } ENDFORALL ;
  }

  //------------------------------------------------------------------------
  // MPI pack/unpack
  //------------------------------------------------------------------------
  template<int M>
  int rMapVecRepI<M>::pack_size(const entitySet &eset) {
    entitySet common = eset & domain() ;
    Entity *p = 0 ;
    return cpypacksize(p, common.size() * M) ;
  }

  template<int M>
  int rMapVecRepI<M>::estimated_pack_size(const entitySet &eset) {
    entitySet common = eset & domain() ;
    Entity *p = 0 ;
    return cpypacksize(p, common.size() * M) ;
  }

  template<int M>
  int rMapVecRepI<M>::pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;
    Entity *p = 0 ;
    return cpypacksize(p, packed.size() * M) ;
  }

  template<int M>
  void rMapVecRepI<M>::pack(void *outbuf, int &position, int &size,
                            const entitySet &usr_eset) {
    entitySet ecommon = domain() & usr_eset ;
    ensure_sorted() ;
    for(auto ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
      auto it = find(*ci) ;
      if(it != attrib_data.end())
        cpypack(outbuf, position, size,
                const_cast<Entity*>(&it->second[0]), M) ;
    }
  }

  template<int M>
  void rMapVecRepI<M>::unpack(void *ptr, int &loc, int &size,
                              const sequence &seq) {
    for(auto ci = seq.begin(); ci != seq.end(); ++ci)
      cpyunpack(ptr, loc, size, &access(*ci)[0], M) ;
  }

  //------------------------------------------------------------------------
  // Type info
  //------------------------------------------------------------------------
  template<int M>
  DatatypeP rMapVecRepI<M>::getType() {
    return DatatypeP(new AtomicType(INT)) ;
  }

  template<int M>
  frame_info rMapVecRepI<M>::get_frame_info() {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }

  //------------------------------------------------------------------------
  // notification / access
  //------------------------------------------------------------------------
  template<int M>
  void rMapVec<M>::notification() {
    NPTR<MapVecType> p(Rep()) ;
    if(p != 0)
      store_ptr = &(*p) ;
    warn(p == 0) ;
  }

  template<int M>
  void const_rMapVec<M>::notification() {
    NPTR<MapVecType> p(Rep()) ;
    if(p != 0)
      store_ptr = &(*p) ;
    warn(p == 0) ;
  }

  template<int M>
  store_instance::instance_type
  const_rMapVec<M>::access() const {
    return READ_ONLY ;
  }

  //------------------------------------------------------------------------
  // HDF5 stubs
  //------------------------------------------------------------------------
  template<int M>
  void rMapVecRepI<M>::readhdf5(hid_t group_id, hid_t dataspace,
                                hid_t dataset, hsize_t dimension,
                                const char* name, frame_info &fi,
                                entitySet &user_eset) {
    warn(true) ;
  }

  template<int M>
  void rMapVecRepI<M>::writehdf5(hid_t group_id, hid_t dataspace,
                                 hid_t dataset, hsize_t dimension,
                                 const char* name,
                                 entitySet& usr_eset) const {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<int M>
  void rMapVecRepI<M>::readhdf5P(hid_t group_id, hid_t dataspace,
                                  hid_t dataset, hsize_t dimension,
                                  const char* name, frame_info &fi,
                                  entitySet &user_eset,
                                  hid_t xfer_plist_id) {
    warn(true) ;
  }

  template<int M>
  void rMapVecRepI<M>::writehdf5P(hid_t group_id, hid_t dataspace,
                                   hid_t dataset, hsize_t dimension,
                                   const char* name, entitySet& usr_eset,
                                   hid_t xfer_plist_id) const {
    warn(true) ;
  }
#endif

} // end of namespace Loci

#endif
