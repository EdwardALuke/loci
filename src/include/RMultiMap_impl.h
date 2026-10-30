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
#ifndef RMULTIMAP_IMPL_H
#define RMULTIMAP_IMPL_H

#include <RMultiMap_def.h>
#include <Map.h>
#include <multiMap.h>
#include <store.h>
#include <distribute.h>

namespace Loci {

  //------------------------------------------------------------------------
  // bulk_insert
  //------------------------------------------------------------------------
  inline void rmultiMapRepI::bulk_insert(
    const std::vector<value_type> &entries)
  {
    if(entries.empty())
      return ;

    for(const auto &e : entries)
      store_domain += e.first ;

    attrib_data.insert(attrib_data.end(), entries.begin(), entries.end()) ;
    sorted = false ;
    ensure_sorted() ;

    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  // bulk_erase
  //------------------------------------------------------------------------
  inline void rmultiMapRepI::bulk_erase(const entitySet &eset) {
    if(eset == EMPTY)
      return ;

    store_domain -= eset ;

    if(attrib_data.empty())
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
  inline void rmultiMapRepI::allocate(const entitySet &ptn) {
    entitySet remove = store_domain - ptn ;
    bulk_erase(remove) ;
    store_domain = ptn ;
    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  inline void rmultiMapRepI::shift(int_type offset) {
    for(auto &p : attrib_data)
      p.first += offset ;
    entitySet new_domain ;
    for(auto ei = store_domain.begin(); ei != store_domain.end(); ++ei)
      new_domain += (*ei + offset) ;
    store_domain = new_domain ;
    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  inline entitySet rmultiMapRepI::domain() const {
    return store_domain ;
  }

  //------------------------------------------------------------------------
  inline entitySet rmultiMapRepI::image(const entitySet &iset) const {
    ensure_sorted() ;
    entitySet common = store_domain & iset ;
    entitySet img ;
    for(auto ei = common.begin(); ei != common.end(); ++ei) {
      auto range = equal_range(*ei) ;
      for(auto it = range.first; it != range.second; ++it)
        img += it->second ;
    }
    return img ;
  }

  //------------------------------------------------------------------------
  inline std::pair<entitySet,entitySet>
  rmultiMapRepI::preimage(const entitySet &codomain) const {
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
  inline storeRepP rmultiMapRepI::get_map() {
    return getRep() ;
  }

  //------------------------------------------------------------------------
  inline void rmultiMapRepI::compose(const dMap &m, const entitySet &context) {
    ensure_sorted() ;
    entitySet common = domain() & context ;
    for(auto &p : attrib_data) {
      if(common.inSet(p.first))
        p.second = m[p.second] ;
    }
  }

  //------------------------------------------------------------------------
  inline storeRepP rmultiMapRepI::expand(entitySet &out_of_dom,
                                         std::vector<entitySet> &init_ptn) {
    return getRep() ;
  }

  //------------------------------------------------------------------------
  // Print / Input
  //------------------------------------------------------------------------
  inline std::ostream &rmultiMapRepI::Print(std::ostream &s) const {
    entitySet dom = domain() ;
    s << '{' << dom << std::endl ;
    ensure_sorted() ;
    for(auto ei = dom.begin(); ei != dom.end(); ++ei) {
      auto range = equal_range(*ei) ;
      int cnt = static_cast<int>(range.second - range.first) ;
      s << cnt ;
      for(auto it = range.first; it != range.second; ++it)
        s << ' ' << it->second ;
      s << std::endl ;
    }
    s << '}' << std::endl ;
    return s ;
  }

  inline std::istream &rmultiMapRepI::Input(std::istream &s) {
    entitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '{') {
      std::cerr << "Incorrect Format while reading rmultiMap" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    s >> e ;

    attrib_data.clear() ;
    store_domain = e ;

    FORALL(e, ii) {
      int cnt ;
      s >> cnt ;
      for(int j=0;j<cnt;++j) {
        Entity val ;
        s >> val ;
        attrib_data.push_back(value_type(ii, val)) ;
      }
    } ENDFORALL ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '}') {
      std::cerr << "Incorrect Format while reading rmultiMap" << std::endl ;
      s.putback(ch) ;
    }

    sorted = false ;
    ensure_sorted() ;

    return s ;
  }

  //------------------------------------------------------------------------
  inline storeRep *rmultiMapRepI::new_store(const entitySet &p) const {
    return new rmultiMapRepI(p) ;
  }

  inline storeRep *rmultiMapRepI::new_store(const entitySet &p,
                                            const int* cnt) const {
    storeRep* sp = 0 ;
    std::cerr << " This method should not be called for an rmultiMap "
              << std::endl ;
    return sp ;
  }

  //------------------------------------------------------------------------
  // remap
  //------------------------------------------------------------------------
  inline storeRepP rmultiMapRepI::remap(const dMap &m) const {
    entitySet newdomain = m.domain() & domain() ;
    rmultiMap s ;
    s.allocate(m.image(newdomain)) ;

    ensure_sorted() ;
    std::vector<value_type> new_entries ;
    FORALL(newdomain, i) {
      auto range = equal_range(i) ;
      for(auto it = range.first; it != range.second; ++it)
        new_entries.push_back(value_type(m[i], it->second)) ;
    } ENDFORALL ;
    s.bulk_insert(new_entries) ;
    return s.Rep() ;
  }

  //------------------------------------------------------------------------
  inline storeRepP rmultiMapRepI::MapRemap(const dMap &dm,
                                           const dMap &rm) const {
    entitySet newdomain = dm.domain() & domain() ;
    rmultiMap s ;
    s.allocate(dm.image(newdomain)) ;

    ensure_sorted() ;
    std::vector<value_type> entries ;
    FORALL(newdomain, i) {
      auto range = equal_range(i) ;
      for(auto it = range.first; it != range.second; ++it)
        entries.push_back(value_type(dm[i], rm[it->second])) ;
    } ENDFORALL ;
    s.bulk_insert(entries) ;
    return s.Rep() ;
  }

  //------------------------------------------------------------------------
  // redistribute
  //------------------------------------------------------------------------
  inline storeRepP
  rmultiMapRepI::redistribute(const std::vector<entitySet>& dom_ptn,
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

    rmultiMap ns ;
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

    if(old_dom != EMPTY) {
      ensure_sorted() ;
      NPTR<rmultiMapRepI> nsp(srp) ;
      std::vector<value_type> old_entries ;
      FORALL(old_dom, ei) {
        auto range = equal_range(ei) ;
        for(auto it = range.first; it != range.second; ++it)
          old_entries.push_back(*it) ;
      } ENDFORALL ;
      if(!old_entries.empty())
        nsp->bulk_insert(old_entries) ;
    }

    return srp ;
  }

  inline storeRepP
  rmultiMapRepI::redistribute(const std::vector<entitySet>& dom_ptn,
                              const dMap& remap, MPI_Comm comm) {
    rmultiMap new_store ;
    Loci::Abort() ;
    return new_store.Rep() ;
  }

  //------------------------------------------------------------------------
  // freeze / thaw
  //------------------------------------------------------------------------
  inline storeRepP rmultiMapRepI::freeze() {
    store<int> counts ;
    entitySet dom = domain() ;
    counts.allocate(dom) ;
    ensure_sorted() ;
    FORALL(dom, ei) {
      counts[ei] = num_elems(ei) ;
    } ENDFORALL ;

    multiMap mm ;
    mm.allocate(counts) ;

    FORALL(dom, ei) {
      auto range = equal_range(ei) ;
      int idx = 0 ;
      for(auto it = range.first; it != range.second; ++it)
        mm[ei][idx++] = it->second ;
    } ENDFORALL ;

    return mm.Rep() ;
  }

  inline storeRepP rmultiMapRepI::thaw() {
    return getRep() ;
  }

  //------------------------------------------------------------------------
  // copy / gather / scatter
  //------------------------------------------------------------------------
  inline void rmultiMapRepI::copy(storeRepP &st, const entitySet &context) {
    const_rmultiMap s(st) ;
    bulk_erase(context) ;
    store_domain += context ;
    std::vector<value_type> entries ;
    FORALL(context, i) {
      rmultiMap_ref vals = s[i] ;
      for(int j=0;j<vals.size();++j)
        entries.push_back(value_type(i, vals[j])) ;
    } ENDFORALL ;
    if(!entries.empty())
      bulk_insert(entries) ;
  }

  inline void rmultiMapRepI::gather(const dMap &m, storeRepP &st,
                                    const entitySet &context) {
    const_rmultiMap s(st) ;
    bulk_erase(context) ;
    store_domain += context ;
    std::vector<value_type> entries ;
    FORALL(context, i) {
      rmultiMap_ref vals = s[m[i]] ;
      for(int j=0;j<vals.size();++j)
        entries.push_back(value_type(i, vals[j])) ;
    } ENDFORALL ;
    if(!entries.empty())
      bulk_insert(entries) ;
  }

  inline void rmultiMapRepI::scatter(const dMap &m, storeRepP &st,
                                     const entitySet &context) {
    const_rmultiMap s(st) ;
    entitySet targets = m.image(context) ;
    bulk_erase(targets) ;
    store_domain += targets ;
    std::vector<value_type> entries ;
    FORALL(context, i) {
      rmultiMap_ref vals = s[i] ;
      for(int j=0;j<vals.size();++j)
        entries.push_back(value_type(m[i], vals[j])) ;
    } ENDFORALL ;
    if(!entries.empty())
      bulk_insert(entries) ;
  }

  //------------------------------------------------------------------------
  // MPI pack/unpack — pack count + targets for each entity
  //------------------------------------------------------------------------
  inline int rmultiMapRepI::pack_size(const entitySet &eset) {
    ensure_sorted() ;
    int numBytes = 0 ;
    Entity *p = 0 ;
    entitySet common = eset & store_domain ;
    for(auto ei = common.begin(); ei != common.end(); ++ei) {
      int cnt = num_elems(*ei) ;
      numBytes += cpypacksize(&cnt, 1) ;
      numBytes += cpypacksize(p, cnt) ;
    }
    return numBytes ;
  }

  inline int rmultiMapRepI::estimated_pack_size(const entitySet &eset) {
    entitySet common = eset & store_domain ;
    int numBytes = common.size() * (sizeof(int) + 4*sizeof(Entity)) ;
    return numBytes ;
  }

  inline int rmultiMapRepI::pack_size(const entitySet& e,
                                      entitySet& packed) {
    packed = domain() & e ;
    return pack_size(packed) ;
  }

  inline void rmultiMapRepI::pack(void *outbuf, int &position, int &size,
                                  const entitySet &usr_eset) {
    entitySet ecommon = domain() & usr_eset ;
    ensure_sorted() ;
    for(auto ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
      auto range = equal_range(*ci) ;
      int cnt = static_cast<int>(range.second - range.first) ;
      cpypack(outbuf, position, size, &cnt, 1) ;
      for(auto it = range.first; it != range.second; ++it)
        cpypack(outbuf, position, size,
                const_cast<Entity*>(&it->second), 1) ;
    }
  }

  inline void rmultiMapRepI::unpack(void *ptr, int &loc, int &size,
                                    const sequence &seq) {
    std::vector<value_type> entries ;
    for(auto ci = seq.begin(); ci != seq.end(); ++ci) {
      int cnt = 0 ;
      cpyunpack(ptr, loc, size, &cnt, 1) ;
      store_domain += *ci ;
      for(int j=0;j<cnt;++j) {
        Entity val ;
        cpyunpack(ptr, loc, size, &val, 1) ;
        entries.push_back(value_type(*ci, val)) ;
      }
    }
    if(!entries.empty())
      bulk_insert(entries) ;
  }

  //------------------------------------------------------------------------
  // Type info
  //------------------------------------------------------------------------
  inline DatatypeP rmultiMapRepI::getType() {
    typedef data_schema_traits<Entity> traits_type ;
    return traits_type::get_type() ;
  }

  inline frame_info rmultiMapRepI::get_frame_info() {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }

  //------------------------------------------------------------------------
  // notification / access
  //------------------------------------------------------------------------
  inline void rmultiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p != 0)
      store_ptr = &(*p) ;
    warn(p == 0) ;
  }

  inline void const_rmultiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p != 0)
      store_ptr = &(*p) ;
    warn(p == 0) ;
  }

  inline store_instance::instance_type
  const_rmultiMap::access() const {
    return READ_ONLY ;
  }

  //------------------------------------------------------------------------
  // HDF5 stubs
  //------------------------------------------------------------------------
  inline void rmultiMapRepI::readhdf5(hid_t group_id, hid_t dataspace,
                                      hid_t dataset, hsize_t dimension,
                                      const char* name, frame_info &fi,
                                      entitySet &user_eset) {
    warn(true) ;
  }

  inline void rmultiMapRepI::writehdf5(hid_t group_id, hid_t dataspace,
                                       hid_t dataset, hsize_t dimension,
                                       const char* name,
                                       entitySet& usr_eset) const {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  inline void rmultiMapRepI::readhdf5P(hid_t group_id, hid_t dataspace,
                                       hid_t dataset, hsize_t dimension,
                                       const char* name, frame_info &fi,
                                       entitySet &user_eset,
                                       hid_t xfer_plist_id) {
    warn(true) ;
  }

  inline void rmultiMapRepI::writehdf5P(hid_t group_id, hid_t dataspace,
                                        hid_t dataset, hsize_t dimension,
                                        const char* name,
                                        entitySet& usr_eset,
                                        hid_t xfer_plist_id) const {
    warn(true) ;
  }
#endif

} // end of namespace Loci

#endif
