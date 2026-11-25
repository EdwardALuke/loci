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
#ifndef RMAP_IMPL_H
#define RMAP_IMPL_H

#include <RMap_def.h>
#include <Map.h>
#include <multiMap.h>
#include <distribute.h>

namespace Loci {

  //------------------------------------------------------------------------
  // bulk_insert
  //------------------------------------------------------------------------
  inline void rMapRepI::bulk_insert(const std::vector<value_type> &entries) {
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
  inline void rMapRepI::bulk_erase(const entitySet &eset) {
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
  inline void rMapRepI::allocate(const entitySet &ptn) {
    entitySet dom = domain() ;
    entitySet remove = dom - ptn ;
    entitySet add = ptn - dom ;

    bulk_erase(remove) ;

    if(add != EMPTY) {
      std::vector<value_type> new_entries ;
      new_entries.reserve(add.size()) ;
      for(auto ei = add.begin(); ei != add.end(); ++ei)
        new_entries.push_back(value_type(*ei, 0)) ;
      bulk_insert(new_entries) ;
    }

    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  inline void rMapRepI::shift(int_type offset) {
    for(auto &p : attrib_data)
      p.first += offset ;
    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  inline entitySet rMapRepI::domain() const {
    ensure_sorted() ;
    entitySet dom ;
    for(const auto &p : attrib_data)
      dom += p.first ;
    return dom ;
  }

  //------------------------------------------------------------------------
  inline entitySet rMapRepI::image(const entitySet &iset) const {
    ensure_sorted() ;
    entitySet dom = domain() ;
    entitySet common = dom & iset ;
    entitySet img ;
    for(auto ei = common.begin(); ei != common.end(); ++ei) {
      auto it = find(*ei) ;
      if(it != attrib_data.end())
        img += it->second ;
    }
    return img ;
  }

  //------------------------------------------------------------------------
  inline std::pair<entitySet,entitySet>
  rMapRepI::preimage(const entitySet &codomain) const {
    ensure_sorted() ;
    entitySet pre ;
    entitySet img ;
    for(const auto &p : attrib_data) {
      if(codomain.inSet(p.second)) {
        pre += p.first ;
        img += p.second ;
      }
    }
    return std::make_pair(pre, img) ;
  }

  //------------------------------------------------------------------------
  inline storeRepP rMapRepI::get_map() {
    return getRep() ;
  }

  //------------------------------------------------------------------------
  inline void rMapRepI::compose(const dMap &m, const entitySet &context) {
    ensure_sorted() ;
    entitySet common = domain() & context ;
    FORALL(common, i) {
      auto it = find_mut(i) ;
      if(it != attrib_data.end())
        it->second = m[it->second] ;
    } ENDFORALL ;
  }

  //------------------------------------------------------------------------
  inline storeRepP rMapRepI::expand(entitySet &out_of_dom,
                                    std::vector<entitySet> &init_ptn) {
    // Not fully implemented for rMap yet
    return getRep() ;
  }

  //------------------------------------------------------------------------
  // Print / Input
  //------------------------------------------------------------------------
  inline std::ostream &rMapRepI::Print(std::ostream &s) const {
    entitySet dom = domain() ;
    s << '{' << dom << std::endl ;
    ensure_sorted() ;
    for(const auto &p : attrib_data)
      s << p.second << ' ' ;
    s << std::endl << '}' << std::endl ;
    return s ;
  }

  inline std::istream &rMapRepI::Input(std::istream &s) {
    entitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '{') {
      std::cerr << "Incorrect Format while reading rMap" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    s >> e ;

    std::vector<value_type> entries ;
    entries.reserve(e.size()) ;
    FORALL(e, ii) {
      Entity val ;
      s >> val ;
      entries.push_back(value_type(ii, val)) ;
    } ENDFORALL ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '}') {
      std::cerr << "Incorrect Format while reading rMap" << std::endl ;
      s.putback(ch) ;
    }

    attrib_data.clear() ;
    bulk_insert(entries) ;

    return s ;
  }

  //------------------------------------------------------------------------
  inline storeRep *rMapRepI::new_store(const entitySet &p) const {
    return new rMapRepI(p) ;
  }

  inline storeRep *rMapRepI::new_store(const entitySet &p,
                                       const int* cnt) const {
    storeRep* sp = 0 ;
    std::cerr << " This method should not be called for an rMap "
              << std::endl ;
    return sp ;
  }

  //------------------------------------------------------------------------
  // remap
  //------------------------------------------------------------------------
  inline storeRepP rMapRepI::remap(const dMap &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    rMap s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m, my_store, newdomain) ;
    return s.Rep() ;
  }

  //------------------------------------------------------------------------
  inline storeRepP rMapRepI::MapRemap(const dMap &dm, const dMap &rm) const {
    entitySet newdomain = dm.domain() & domain() ;
    rMap s ;
    s.allocate(dm.image(newdomain)) ;

    ensure_sorted() ;
    std::vector<value_type> entries ;
    FORALL(newdomain, i) {
      auto it = find(i) ;
      if(it != attrib_data.end())
        entries.push_back(value_type(dm[i], rm[it->second])) ;
    } ENDFORALL ;
    s.bulk_insert(entries) ;
    return s.Rep() ;
  }

  //------------------------------------------------------------------------
  // redistribute
  //------------------------------------------------------------------------
  inline storeRepP
  rMapRepI::redistribute(const std::vector<entitySet>& dom_ptn,
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

    rMap ns ;
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

    // copy old_dom entries
    ensure_sorted() ;
    for(entitySet::const_iterator ei=old_dom.begin();
        ei!=old_dom.end();++ei) {
      auto it = find(*ei) ;
      if(it != attrib_data.end()) {
        NPTR<rMapRepI> nsp(srp) ;
        nsp->access(*ei) = it->second ;
      }
    }

    return srp ;
  }

  inline storeRepP
  rMapRepI::redistribute(const std::vector<entitySet>& dom_ptn,
                         const dMap& remap, MPI_Comm comm) {
    rMap new_store ;
    Loci::Abort() ;
    return new_store.Rep() ;
  }

  //------------------------------------------------------------------------
  // freeze / thaw
  //------------------------------------------------------------------------
  inline storeRepP rMapRepI::freeze() {
    Map static_map ;
    entitySet dom = domain() ;
    static_map.allocate(dom) ;
    ensure_sorted() ;
    for(const auto &p : attrib_data)
      static_map[p.first] = p.second ;
    return static_map.Rep() ;
  }

  inline storeRepP rMapRepI::thaw() {
    return getRep() ;
  }

  //------------------------------------------------------------------------
  // copy / gather / scatter
  //------------------------------------------------------------------------
  inline void rMapRepI::copy(storeRepP &st, const entitySet &context) {
    const_rMap s(st) ;
    fatal((context-domain()) != EMPTY) ;
    FORALL(context, i) {
      access(i) = s[i] ;
    } ENDFORALL ;
  }

  inline void rMapRepI::gather(const dMap &m, storeRepP &st,
                                const entitySet &context) {
    const_rMap s(st) ;
    fatal(context != EMPTY) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context, i) {
      access(i) = s[m[i]] ;
    } ENDFORALL ;
  }

  inline void rMapRepI::scatter(const dMap &m, storeRepP &st,
                                 const entitySet &context) {
    const_rMap s(st) ;
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
  inline int rMapRepI::pack_size(const entitySet &eset) {
    Entity *p = 0 ;
    return cpypacksize(p, (eset & domain()).size()) ;
  }

  inline int rMapRepI::estimated_pack_size(const entitySet &eset) {
    Entity *p = 0 ;
    return cpypacksize(p, (eset & domain()).size()) ;
  }

  inline int rMapRepI::pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;
    Entity *p = 0 ;
    return cpypacksize(p, packed.size()) ;
  }

  inline void rMapRepI::pack(void *outbuf, int &position, int &size,
                              const entitySet &usr_eset) {
    entitySet ecommon = domain() & usr_eset ;
    ensure_sorted() ;
    for(auto ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
      auto it = find(*ci) ;
      if(it != attrib_data.end())
        cpypack(outbuf, position, size,
                const_cast<Entity*>(&it->second), 1) ;
    }
  }

  inline void rMapRepI::unpack(void *ptr, int &loc, int &size,
                                const sequence &seq) {
    for(auto ci = seq.begin(); ci != seq.end(); ++ci)
      cpyunpack(ptr, loc, size, &access(*ci), 1) ;
  }

  //------------------------------------------------------------------------
  // Type info
  //------------------------------------------------------------------------
  inline DatatypeP rMapRepI::getType() {
    typedef data_schema_traits<Entity> traits_type ;
    return traits_type::get_type() ;
  }

  inline frame_info rMapRepI::get_frame_info() {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }

  //------------------------------------------------------------------------
  // notification / access
  //------------------------------------------------------------------------
  inline void rMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p != 0)
      store_ptr = &(*p) ;
    warn(p == 0) ;
  }

  inline void const_rMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p != 0)
      store_ptr = &(*p) ;
    warn(p == 0) ;
  }

  inline store_instance::instance_type
  const_rMap::access() const {
    return READ_ONLY ;
  }

  //------------------------------------------------------------------------
  // HDF5 stubs
  //------------------------------------------------------------------------
  inline void rMapRepI::readhdf5(hid_t group_id, hid_t dataspace,
                                  hid_t dataset, hsize_t dimension,
                                  const char* name, frame_info &fi,
                                  entitySet &user_eset) {
    warn(true) ;
  }

  inline void rMapRepI::writehdf5(hid_t group_id, hid_t dataspace,
                                   hid_t dataset, hsize_t dimension,
                                   const char* name,
                                   entitySet& usr_eset) const {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  inline void rMapRepI::readhdf5P(hid_t group_id, hid_t dataspace,
                                    hid_t dataset, hsize_t dimension,
                                    const char* name, frame_info &fi,
                                    entitySet &user_eset,
                                    hid_t xfer_plist_id) {
    warn(true) ;
  }

  inline void rMapRepI::writehdf5P(hid_t group_id, hid_t dataspace,
                                    hid_t dataset, hsize_t dimension,
                                    const char* name, entitySet& usr_eset,
                                    hid_t xfer_plist_id) const {
    warn(true) ;
  }
#endif

} // end of namespace Loci

#endif
