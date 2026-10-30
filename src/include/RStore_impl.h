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
#ifndef RSTORE_IMPL_H
#define RSTORE_IMPL_H

#include <RStore_def.h>
#include <store_def.h>
#include <multiStore.h>
#include <distribute.h>

namespace Loci {

  //------------------------------------------------------------------------
  // bulk_insert: merge new entries into the sorted vector.
  // Existing entries with matching entities are overwritten.
  //------------------------------------------------------------------------
  template<class T>
  void rstoreRepI<T>::bulk_insert(const std::vector<value_type> &entries) {
    if(entries.empty())
      return ;

    // Sort the incoming entries by entity
    std::vector<value_type> sorted_new(entries) ;
    std::sort(sorted_new.begin(), sorted_new.end(),
              [](const value_type &a, const value_type &b)
              { return a.first < b.first ; }) ;

    ensure_sorted() ;

    // Merge: walk both sorted sequences
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
        // Same entity: new value overwrites old
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
  // bulk_erase: remove all entries whose entities are in eset.
  //------------------------------------------------------------------------
  template<class T>
  void rstoreRepI<T>::bulk_erase(const entitySet &eset) {
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
  // allocate: set the domain to ptn, default-constructing new entries
  // and removing entries not in ptn.
  //------------------------------------------------------------------------
  template<class T>
  void rstoreRepI<T>::allocate(const entitySet &ptn) {
    entitySet dom = domain() ;
    entitySet remove = dom - ptn ;
    entitySet add = ptn - dom ;

    bulk_erase(remove) ;

    if(add != EMPTY) {
      std::vector<value_type> new_entries ;
      new_entries.reserve(add.size()) ;
      for(auto ei = add.begin(); ei != add.end(); ++ei)
        new_entries.push_back(value_type(*ei, T())) ;
      bulk_insert(new_entries) ;
    }

    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  template<class T>
  void rstoreRepI<T>::shift(int_type offset) {
    for(auto &p : attrib_data)
      p.first += offset ;
    // Still sorted if offset is uniform
    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  template<class T>
  std::ostream &rstoreRepI<T>::Print(std::ostream &s) const {
    entitySet dom = domain() ;
    s << '{' << dom << std::endl ;
    ensure_sorted() ;
    for(const auto &p : attrib_data)
      Loci::streamoutput(&p.second, 1, s) ;
    s << '}' << std::endl ;
    return s ;
  }

  //------------------------------------------------------------------------
  template<class T>
  std::istream &rstoreRepI<T>::Input(std::istream &s) {
    entitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    s >> e ;

    std::vector<value_type> entries ;
    entries.reserve(e.size()) ;
    FORALL(e, ii) {
      T val ;
      Loci::streaminput(&val, 1, s) ;
      entries.push_back(value_type(ii, val)) ;
    } ENDFORALL ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }

    attrib_data.clear() ;
    bulk_insert(entries) ;

    return s ;
  }

  //------------------------------------------------------------------------
  template<class T>
  entitySet rstoreRepI<T>::domain() const {
    ensure_sorted() ;
    entitySet dom ;
    for(const auto &p : attrib_data)
      dom += p.first ;
    return dom ;
  }

  //------------------------------------------------------------------------
  template<class T>
  storeRep *rstoreRepI<T>::new_store(const entitySet &p) const {
    return new rstoreRepI<T>(p) ;
  }

  template<class T>
  storeRep *rstoreRepI<T>::new_store(const entitySet &p,
                                     const int* cnt) const {
    storeRep* sp = 0 ;
    std::cerr << " This method should not be called for an rstore "
              << std::endl ;
    return sp ;
  }

  //------------------------------------------------------------------------
  template<class T>
  store_type rstoreRepI<T>::RepType() const {
    return STORE ;
  }

  //------------------------------------------------------------------------
  template<class T>
  void rstore<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      store_ptr = &(*p) ;
    warn(p == 0) ;
  }

  template<class T>
  inline std::ostream & operator<<(std::ostream &s, const rstore<T> &t)
  { return t.Print(s) ; }

  template<class T>
  inline std::istream & operator>>(std::istream &s, rstore<T> &t)
  { return t.Input(s) ; }

  //------------------------------------------------------------------------
  template<class T>
  void const_rstore<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      store_ptr = &(*p) ;
    warn(p == 0) ;
  }

  template<class T>
  store_instance::instance_type
  const_rstore<T>::access() const {
    return READ_ONLY ;
  }

  //------------------------------------------------------------------------
  // remap
  //------------------------------------------------------------------------
  template<class T>
  storeRepP rstoreRepI<T>::remap(const dMap &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    rstore<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m, my_store, newdomain) ;
    return s.Rep() ;
  }

  //------------------------------------------------------------------------
  // redistribute
  //------------------------------------------------------------------------
  template<class T>
  storeRepP
  rstoreRepI<T>::redistribute(const std::vector<entitySet>& dom_ptn,
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

    rstore<T> ns ;
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

    // copy old_dom if any
    ensure_sorted() ;
    for(entitySet::const_iterator ei=old_dom.begin();
        ei!=old_dom.end();++ei) {
      auto it = find(*ei) ;
      if(it != attrib_data.end()) {
        NPTR<rstoreRepI<T>> nsp(srp) ;
        nsp->access(*ei) = it->second ;
      }
    }

    return srp ;
  }

  template<class T>
  storeRepP
  rstoreRepI<T>::redistribute(const std::vector<entitySet>& dom_ptn,
                              const dMap& remap, MPI_Comm comm) {
    std::vector<P2pCommInfo> send, recv ;
    get_p2p_comm(dom_ptn, domain(), 0, 0, comm, recv, send) ;
    rstore<T> new_store ;
    Loci::Abort() ;
    return new_store.Rep() ;
  }

  //------------------------------------------------------------------------
  // freeze / thaw
  //------------------------------------------------------------------------
  template<class T>
  storeRepP rstoreRepI<T>::freeze() {
    store<T> static_store ;
    entitySet tmp_dom = domain() ;
    static_store.allocate(tmp_dom) ;
    ensure_sorted() ;
    for(const auto &p : attrib_data)
      static_store[p.first] = p.second ;
    return static_store.Rep() ;
  }

  template<class T>
  storeRepP rstoreRepI<T>::thaw() {
    return getRep() ;
  }

  //------------------------------------------------------------------------
  // copy / gather / scatter
  //------------------------------------------------------------------------
  template<class T>
  void rstoreRepI<T>::copy(storeRepP &st, const entitySet &context) {
    const_rstore<T> s(st) ;
    fatal((context-domain()) != EMPTY) ;
    FORALL(context, i) {
      access(i) = s[i] ;
    } ENDFORALL ;
  }

  template<class T>
  void rstoreRepI<T>::gather(const dMap &m, storeRepP &st,
                             const entitySet &context) {
    const_rstore<T> s(st) ;
    fatal(context != EMPTY) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context, i) {
      access(i) = s[m[i]] ;
    } ENDFORALL ;
  }

  template<class T>
  void rstoreRepI<T>::scatter(const dMap &m, storeRepP &st,
                              const entitySet &context) {
    const_rstore<T> s(st) ;
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
  template<class T>
  int rstoreRepI<T>::pack_size(const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    schema_converter traits_type ;
    return get_mpi_size(traits_type, eset) ;
  }

  template<class T>
  int rstoreRepI<T>::estimated_pack_size(const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    schema_converter traits_type ;
    return get_estimated_mpi_size(traits_type, eset) ;
  }

  template<class T>
  int rstoreRepI<T>::pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    schema_converter traits_type ;
    return get_mpi_size(traits_type, packed) ;
  }

  template<class T>
  int rstoreRepI<T>::get_mpi_size(IDENTITY_CONVERTER c,
                                  const entitySet &eset) {
    T *p = 0 ;
    return cpypacksize(p, eset.size()) ;
  }

  template<class T>
  int rstoreRepI<T>::get_estimated_mpi_size(IDENTITY_CONVERTER c,
                                            const entitySet &eset) {
    T *p = 0 ;
    return cpypacksize(p, eset.size()) ;
  }

  template<class T>
  int rstoreRepI<T>::get_mpi_size(USER_DEFINED_CONVERTER c,
                                  const entitySet &eset) {
    int numBytes = 0 ;
    for(auto ci = eset.begin(); ci != eset.end(); ++ci) {
      auto it = find(*ci) ;
      if(it != attrib_data.end()) {
        typedef data_schema_traits<T> converter_traits ;
        typename converter_traits::Converter_Type cvtr(
          const_cast<T&>(it->second)) ;
        int size = cvtr.getSize() ;
        numBytes += cpypacksize(&size, 1) ;
        typename converter_traits::Converter_Base_Type *p = 0 ;
        numBytes += cpypacksize(p, size) ;
      }
    }
    return numBytes ;
  }

  template<class T>
  int rstoreRepI<T>::get_estimated_mpi_size(USER_DEFINED_CONVERTER c,
                                            const entitySet &eset) {
    int numBytes = 0 ;
    numBytes = eset.size()*50*sizeof(double) ;
    numBytes += eset.size()*sizeof(int) ;
    return numBytes ;
  }

  template<class T>
  void rstoreRepI<T>::pack(void *outbuf, int &position, int &size,
                           const entitySet &usr_eset) {
    entitySet ecommon ;
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    schema_converter traits_type ;
    ecommon = domain() & usr_eset ;
    packdata(traits_type, outbuf, position, size, ecommon) ;
  }

  template<class T>
  void rstoreRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
                               int &position, int outcount,
                               const entitySet &eset) {
    for(auto ci = eset.begin(); ci != eset.end(); ++ci) {
      auto it = find(*ci) ;
      if(it != attrib_data.end())
        cpypack(outbuf, position, outcount,
                const_cast<T*>(&it->second), 1) ;
    }
  }

  template<class T>
  void rstoreRepI<T>::packdata(USER_DEFINED_CONVERTER c, void *outbuf,
                               int &position, int outcount,
                               const entitySet &eset) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype ;
    std::vector<dtype> inbuf ;

    for(auto ci = eset.begin(); ci != eset.end(); ++ci) {
      auto it = find(*ci) ;
      if(it != attrib_data.end()) {
        typename schema_traits::Converter_Type cvtr(
          const_cast<T&>(it->second)) ;
        int stateSize = cvtr.getSize() ;
        if(stateSize > static_cast<int>(inbuf.size()))
          inbuf.resize(stateSize) ;
        cvtr.getState(inbuf.data(), stateSize) ;
        cpypack(outbuf, position, outcount, &stateSize, 1) ;
        cpypack(outbuf, position, outcount, inbuf.data(), stateSize) ;
      }
    }
  }

  template<class T>
  void rstoreRepI<T>::unpack(void *ptr, int &loc, int &size,
                             const sequence &seq) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    schema_converter traits_type ;
    unpackdata(traits_type, ptr, loc, size, seq) ;
  }

  template<class T>
  void rstoreRepI<T>::unpackdata(IDENTITY_CONVERTER c, void *inbuf,
                                 int &position, int insize,
                                 const sequence &seq) {
    for(auto ci = seq.begin(); ci != seq.end(); ++ci)
      cpyunpack(inbuf, position, insize, &access(*ci), 1) ;
  }

  template<class T>
  void rstoreRepI<T>::unpackdata(USER_DEFINED_CONVERTER c, void *inbuf,
                                 int &position, int insize,
                                 const sequence &seq) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype ;
    std::vector<dtype> outbuf ;

    for(auto ci = seq.begin(); ci != seq.end(); ++ci) {
      int stateSize = 0 ;
      cpyunpack(inbuf, position, insize, &stateSize, 1) ;
      if(stateSize > static_cast<int>(outbuf.size()))
        outbuf.resize(stateSize) ;
      cpyunpack(inbuf, position, insize, outbuf.data(), stateSize) ;
      typename schema_traits::Converter_Type cvtr(access(*ci)) ;
      cvtr.setState(outbuf.data(), stateSize) ;
    }
  }

  //------------------------------------------------------------------------
  // Type info / frame_info
  //------------------------------------------------------------------------
  template<class T>
  frame_info rstoreRepI<T>::get_frame_info() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    return get_frame_info(schema_converter()) ;
  }

  template<class T>
  frame_info rstoreRepI<T>::get_frame_info(IDENTITY_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }

  template<class T>
  frame_info rstoreRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }

  template<class T>
  DatatypeP rstoreRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    return getType(schema_converter()) ;
  }

  template<class T>
  DatatypeP rstoreRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type ;
    return traits_type::get_type() ;
  }

  template<class T>
  DatatypeP rstoreRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype ;
    typedef data_schema_traits<dtype> traits_type ;
    return traits_type::get_type() ;
  }

  //------------------------------------------------------------------------
  // HDF5 (stubs — same as dstore)
  //------------------------------------------------------------------------
  template<class T>
  void rstoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace,
                               hid_t dataset, hsize_t dimension,
                               const char* name, frame_info &fi,
                               entitySet &user_eset) {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rstoreRepI<T>::readhdf5P(hid_t group_id, hid_t dataspace,
                                hid_t dataset, hsize_t dimension,
                                const char* name, frame_info &fi,
                                entitySet &user_eset, hid_t xfer_plist_id) {
    warn(true) ;
  }
#endif

  template<class T>
  void rstoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace,
                               hid_t dataset, hsize_t dimension,
                               const char* name, IDENTITY_CONVERTER c,
                               frame_info &fi, entitySet &eset) {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rstoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace,
                                hid_t dataset, hsize_t dimension,
                                const char* name, IDENTITY_CONVERTER c,
                                frame_info &fi, entitySet &eset,
                                hid_t xfer_plist_id) {
    warn(true) ;
  }
#endif

  template<class T>
  void rstoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace,
                               hid_t dataset, hsize_t dimension,
                               const char* name, USER_DEFINED_CONVERTER c,
                               frame_info &fi, entitySet &eset) {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rstoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace,
                                hid_t dataset, hsize_t dimension,
                                const char* name, USER_DEFINED_CONVERTER c,
                                frame_info &fi, entitySet &eset,
                                hid_t xfer_plist_id) {
    warn(true) ;
  }
#endif

  template<class T>
  void rstoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace,
                                hid_t dataset, hsize_t dimension,
                                const char* name, entitySet& usr_eset) const {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rstoreRepI<T>::writehdf5P(hid_t group_id, hid_t dataspace,
                                 hid_t dataset, hsize_t dimension,
                                 const char* name, entitySet& usr_eset,
                                 hid_t xfer_plist_id) const {
    warn(true) ;
  }
#endif

  template<class T>
  void rstoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace,
                                hid_t dataset, hsize_t dimension,
                                const char* name, IDENTITY_CONVERTER c,
                                const entitySet &eset) const {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rstoreRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace,
                                 hid_t dataset, hsize_t dimension,
                                 const char* name, IDENTITY_CONVERTER c,
                                 const entitySet &en,
                                 hid_t xfer_plist_id) const {
    warn(true) ;
  }
#endif

  template<class T>
  void rstoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace,
                                hid_t dataset, hsize_t dimension,
                                const char* name, USER_DEFINED_CONVERTER c,
                                const entitySet &eset) const {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rstoreRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace,
                                 hid_t dataset, hsize_t dimension,
                                 const char* name, USER_DEFINED_CONVERTER c,
                                 const entitySet &en,
                                 hid_t xfer_plist_id) const {
    warn(true) ;
  }
#endif

  //========================================================================
  // rmultiStoreRepI implementation
  //========================================================================

  //------------------------------------------------------------------------
  // bulk_insert: append entries preserving insertion order within each
  // entity (stable sort).
  //------------------------------------------------------------------------
  template<class T>
  void rmultiStoreRepI<T>::bulk_insert(const std::vector<value_type> &entries) {
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
  // bulk_erase: remove all entries for entities in eset, and remove
  // those entities from the domain.
  //------------------------------------------------------------------------
  template<class T>
  void rmultiStoreRepI<T>::bulk_erase(const entitySet &eset) {
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
  // allocate: set the domain.  Entities in ptn that are not already in
  // the domain are added (with zero values).  Entities not in ptn are
  // removed along with their values.
  //------------------------------------------------------------------------
  template<class T>
  void rmultiStoreRepI<T>::allocate(const entitySet &ptn) {
    entitySet remove = store_domain - ptn ;
    bulk_erase(remove) ;
    store_domain = ptn ;
    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  template<class T>
  void rmultiStoreRepI<T>::shift(int_type offset) {
    for(auto &p : attrib_data)
      p.first += offset ;
    entitySet new_domain ;
    for(auto ei = store_domain.begin(); ei != store_domain.end(); ++ei)
      new_domain += (*ei + offset) ;
    store_domain = new_domain ;
    dispatch_notify() ;
  }

  //------------------------------------------------------------------------
  template<class T>
  std::ostream &rmultiStoreRepI<T>::Print(std::ostream &s) const {
    entitySet dom = domain() ;
    s << '{' << dom << std::endl ;
    ensure_sorted() ;
    for(auto ei = dom.begin(); ei != dom.end(); ++ei) {
      auto range = equal_range(*ei) ;
      int cnt = static_cast<int>(range.second - range.first) ;
      s << cnt << std::endl ;
      for(auto it = range.first; it != range.second; ++it)
        Loci::streamoutput(&it->second, 1, s) ;
    }
    s << '}' << std::endl ;
    return s ;
  }

  //------------------------------------------------------------------------
  template<class T>
  std::istream &rmultiStoreRepI<T>::Input(std::istream &s) {
    entitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '{') {
      std::cerr << "Incorrect Format while reading rmultiStore" << std::endl ;
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
        T val ;
        Loci::streaminput(&val, 1, s) ;
        attrib_data.push_back(value_type(ii, val)) ;
      }
    } ENDFORALL ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '}') {
      std::cerr << "Incorrect Format while reading rmultiStore" << std::endl ;
      s.putback(ch) ;
    }

    sorted = false ;
    ensure_sorted() ;

    return s ;
  }

  //------------------------------------------------------------------------
  template<class T>
  entitySet rmultiStoreRepI<T>::domain() const {
    return store_domain ;
  }

  //------------------------------------------------------------------------
  template<class T>
  storeRep *rmultiStoreRepI<T>::new_store(const entitySet &p) const {
    return new rmultiStoreRepI<T>(p) ;
  }

  template<class T>
  storeRep *rmultiStoreRepI<T>::new_store(const entitySet &p,
                                          const int* cnt) const {
    storeRep* sp = 0 ;
    std::cerr << " This method should not be called for an rmultiStore "
              << std::endl ;
    return sp ;
  }

  //------------------------------------------------------------------------
  template<class T>
  store_type rmultiStoreRepI<T>::RepType() const {
    return STORE ;
  }

  //------------------------------------------------------------------------
  template<class T>
  void rmultiStore<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      store_ptr = &(*p) ;
    warn(p == 0) ;
  }

  template<class T>
  inline std::ostream & operator<<(std::ostream &s, const rmultiStore<T> &t)
  { return t.Print(s) ; }

  template<class T>
  inline std::istream & operator>>(std::istream &s, rmultiStore<T> &t)
  { return t.Input(s) ; }

  //------------------------------------------------------------------------
  template<class T>
  void const_rmultiStore<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      store_ptr = &(*p) ;
    warn(p == 0) ;
  }

  template<class T>
  store_instance::instance_type
  const_rmultiStore<T>::access() const {
    return READ_ONLY ;
  }

  //------------------------------------------------------------------------
  // remap
  //------------------------------------------------------------------------
  template<class T>
  storeRepP rmultiStoreRepI<T>::remap(const dMap &m) const {
    entitySet newdomain = m.domain() & domain() ;
    rmultiStore<T> s ;
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
  // redistribute
  //------------------------------------------------------------------------
  template<class T>
  storeRepP
  rmultiStoreRepI<T>::redistribute(const std::vector<entitySet>& dom_ptn,
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

    rmultiStore<T> ns ;
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
    if(old_dom != EMPTY) {
      ensure_sorted() ;
      NPTR<rmultiStoreRepI<T>> nsp(srp) ;
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

  template<class T>
  storeRepP
  rmultiStoreRepI<T>::redistribute(const std::vector<entitySet>& dom_ptn,
                                   const dMap& remap, MPI_Comm comm) {
    std::vector<P2pCommInfo> send, recv ;
    get_p2p_comm(dom_ptn, domain(), 0, 0, comm, recv, send) ;
    rmultiStore<T> new_store ;
    Loci::Abort() ;
    return new_store.Rep() ;
  }

  //------------------------------------------------------------------------
  // freeze / thaw
  //------------------------------------------------------------------------
  template<class T>
  storeRepP rmultiStoreRepI<T>::freeze() {
    // Convert to a traditional multiStore
    store<int> counts ;
    entitySet dom = domain() ;
    counts.allocate(dom) ;
    ensure_sorted() ;
    FORALL(dom, ei) {
      counts[ei] = num_elems(ei) ;
    } ENDFORALL ;

    multiStore<T> ms ;
    ms.allocate(counts) ;

    FORALL(dom, ei) {
      auto range = equal_range(ei) ;
      int idx = 0 ;
      for(auto it = range.first; it != range.second; ++it)
        ms[ei][idx++] = it->second ;
    } ENDFORALL ;

    return ms.Rep() ;
  }

  template<class T>
  storeRepP rmultiStoreRepI<T>::thaw() {
    return getRep() ;
  }

  //------------------------------------------------------------------------
  // copy / gather / scatter
  //------------------------------------------------------------------------
  template<class T>
  void rmultiStoreRepI<T>::copy(storeRepP &st, const entitySet &context) {
    const_rmultiStore<T> s(st) ;
    // Remove existing entries for context entities, then copy in
    bulk_erase(context) ;
    store_domain += context ;
    std::vector<value_type> entries ;
    FORALL(context, i) {
      rmultiStore_ref<T> vals = s[i] ;
      for(int j=0;j<vals.size();++j)
        entries.push_back(value_type(i, vals[j])) ;
    } ENDFORALL ;
    if(!entries.empty())
      bulk_insert(entries) ;
  }

  template<class T>
  void rmultiStoreRepI<T>::gather(const dMap &m, storeRepP &st,
                                  const entitySet &context) {
    const_rmultiStore<T> s(st) ;
    bulk_erase(context) ;
    store_domain += context ;
    std::vector<value_type> entries ;
    FORALL(context, i) {
      rmultiStore_ref<T> vals = s[m[i]] ;
      for(int j=0;j<vals.size();++j)
        entries.push_back(value_type(i, vals[j])) ;
    } ENDFORALL ;
    if(!entries.empty())
      bulk_insert(entries) ;
  }

  template<class T>
  void rmultiStoreRepI<T>::scatter(const dMap &m, storeRepP &st,
                                   const entitySet &context) {
    const_rmultiStore<T> s(st) ;
    entitySet targets = m.image(context) ;
    bulk_erase(targets) ;
    store_domain += targets ;
    std::vector<value_type> entries ;
    FORALL(context, i) {
      rmultiStore_ref<T> vals = s[i] ;
      for(int j=0;j<vals.size();++j)
        entries.push_back(value_type(m[i], vals[j])) ;
    } ENDFORALL ;
    if(!entries.empty())
      bulk_insert(entries) ;
  }

  //------------------------------------------------------------------------
  // MPI pack/unpack — pack count + values for each entity
  //------------------------------------------------------------------------
  template<class T>
  int rmultiStoreRepI<T>::pack_size(const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    schema_converter traits_type ;
    return get_mpi_size(traits_type, eset) ;
  }

  template<class T>
  int rmultiStoreRepI<T>::estimated_pack_size(const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    schema_converter traits_type ;
    return get_estimated_mpi_size(traits_type, eset) ;
  }

  template<class T>
  int rmultiStoreRepI<T>::pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    schema_converter traits_type ;
    return get_mpi_size(traits_type, packed) ;
  }

  template<class T>
  int rmultiStoreRepI<T>::get_mpi_size(IDENTITY_CONVERTER c,
                                       const entitySet &eset) {
    ensure_sorted() ;
    int numBytes = 0 ;
    T *p = 0 ;
    entitySet common = eset & store_domain ;
    for(auto ei = common.begin(); ei != common.end(); ++ei) {
      int cnt = num_elems(*ei) ;
      numBytes += cpypacksize(&cnt, 1) ;
      numBytes += cpypacksize(p, cnt) ;
    }
    return numBytes ;
  }

  template<class T>
  int rmultiStoreRepI<T>::get_estimated_mpi_size(IDENTITY_CONVERTER c,
                                                 const entitySet &eset) {
    entitySet common = eset & store_domain ;
    int numBytes = common.size() * (sizeof(int) + 4*sizeof(T)) ;
    return numBytes ;
  }

  template<class T>
  int rmultiStoreRepI<T>::get_mpi_size(USER_DEFINED_CONVERTER c,
                                       const entitySet &eset) {
    ensure_sorted() ;
    int numBytes = 0 ;
    entitySet common = eset & store_domain ;
    for(auto ci = common.begin(); ci != common.end(); ++ci) {
      auto range = equal_range(*ci) ;
      int cnt = static_cast<int>(range.second - range.first) ;
      numBytes += cpypacksize(&cnt, 1) ;
      for(auto it = range.first; it != range.second; ++it) {
        typedef data_schema_traits<T> converter_traits ;
        typename converter_traits::Converter_Type cvtr(
          const_cast<T&>(it->second)) ;
        int size = cvtr.getSize() ;
        numBytes += cpypacksize(&size, 1) ;
        typename converter_traits::Converter_Base_Type *bp = 0 ;
        numBytes += cpypacksize(bp, size) ;
      }
    }
    return numBytes ;
  }

  template<class T>
  int rmultiStoreRepI<T>::get_estimated_mpi_size(USER_DEFINED_CONVERTER c,
                                                 const entitySet &eset) {
    entitySet common = eset & store_domain ;
    int numBytes = common.size() * (sizeof(int) + 50*sizeof(double)) ;
    return numBytes ;
  }

  template<class T>
  void rmultiStoreRepI<T>::pack(void *outbuf, int &position, int &size,
                                const entitySet &usr_eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    schema_converter traits_type ;
    entitySet ecommon = domain() & usr_eset ;
    packdata(traits_type, outbuf, position, size, ecommon) ;
  }

  template<class T>
  void rmultiStoreRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
                                    int &position, int outcount,
                                    const entitySet &eset) {
    ensure_sorted() ;
    for(auto ci = eset.begin(); ci != eset.end(); ++ci) {
      auto range = equal_range(*ci) ;
      int cnt = static_cast<int>(range.second - range.first) ;
      cpypack(outbuf, position, outcount, &cnt, 1) ;
      for(auto it = range.first; it != range.second; ++it)
        cpypack(outbuf, position, outcount,
                const_cast<T*>(&it->second), 1) ;
    }
  }

  template<class T>
  void rmultiStoreRepI<T>::packdata(USER_DEFINED_CONVERTER c, void *outbuf,
                                    int &position, int outcount,
                                    const entitySet &eset) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype ;
    std::vector<dtype> inbuf ;

    ensure_sorted() ;
    for(auto ci = eset.begin(); ci != eset.end(); ++ci) {
      auto range = equal_range(*ci) ;
      int cnt = static_cast<int>(range.second - range.first) ;
      cpypack(outbuf, position, outcount, &cnt, 1) ;
      for(auto it = range.first; it != range.second; ++it) {
        typename schema_traits::Converter_Type cvtr(
          const_cast<T&>(it->second)) ;
        int stateSize = cvtr.getSize() ;
        if(stateSize > static_cast<int>(inbuf.size()))
          inbuf.resize(stateSize) ;
        cvtr.getState(inbuf.data(), stateSize) ;
        cpypack(outbuf, position, outcount, &stateSize, 1) ;
        cpypack(outbuf, position, outcount, inbuf.data(), stateSize) ;
      }
    }
  }

  template<class T>
  void rmultiStoreRepI<T>::unpack(void *ptr, int &loc, int &size,
                                  const sequence &seq) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    schema_converter traits_type ;
    unpackdata(traits_type, ptr, loc, size, seq) ;
  }

  template<class T>
  void rmultiStoreRepI<T>::unpackdata(IDENTITY_CONVERTER c, void *inbuf,
                                      int &position, int insize,
                                      const sequence &seq) {
    std::vector<value_type> entries ;
    for(auto ci = seq.begin(); ci != seq.end(); ++ci) {
      int cnt = 0 ;
      cpyunpack(inbuf, position, insize, &cnt, 1) ;
      store_domain += *ci ;
      for(int j=0;j<cnt;++j) {
        T val ;
        cpyunpack(inbuf, position, insize, &val, 1) ;
        entries.push_back(value_type(*ci, val)) ;
      }
    }
    if(!entries.empty())
      bulk_insert(entries) ;
  }

  template<class T>
  void rmultiStoreRepI<T>::unpackdata(USER_DEFINED_CONVERTER c, void *inbuf,
                                      int &position, int insize,
                                      const sequence &seq) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype ;
    std::vector<dtype> outbuf ;
    std::vector<value_type> entries ;

    for(auto ci = seq.begin(); ci != seq.end(); ++ci) {
      int cnt = 0 ;
      cpyunpack(inbuf, position, insize, &cnt, 1) ;
      store_domain += *ci ;
      for(int j=0;j<cnt;++j) {
        int stateSize = 0 ;
        cpyunpack(inbuf, position, insize, &stateSize, 1) ;
        if(stateSize > static_cast<int>(outbuf.size()))
          outbuf.resize(stateSize) ;
        cpyunpack(inbuf, position, insize, outbuf.data(), stateSize) ;
        T val ;
        typename schema_traits::Converter_Type cvtr(val) ;
        cvtr.setState(outbuf.data(), stateSize) ;
        entries.push_back(value_type(*ci, val)) ;
      }
    }
    if(!entries.empty())
      bulk_insert(entries) ;
  }

  //------------------------------------------------------------------------
  // Type info / frame_info
  //------------------------------------------------------------------------
  template<class T>
  frame_info rmultiStoreRepI<T>::get_frame_info() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    return get_frame_info(schema_converter()) ;
  }

  template<class T>
  frame_info rmultiStoreRepI<T>::get_frame_info(IDENTITY_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }

  template<class T>
  frame_info rmultiStoreRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }

  template<class T>
  DatatypeP rmultiStoreRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter ;
    return getType(schema_converter()) ;
  }

  template<class T>
  DatatypeP rmultiStoreRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type ;
    return traits_type::get_type() ;
  }

  template<class T>
  DatatypeP rmultiStoreRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype ;
    typedef data_schema_traits<dtype> traits_type ;
    return traits_type::get_type() ;
  }

  //------------------------------------------------------------------------
  // HDF5 (stubs)
  //------------------------------------------------------------------------
  template<class T>
  void rmultiStoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace,
                                    hid_t dataset, hsize_t dimension,
                                    const char* name, frame_info &fi,
                                    entitySet &user_eset) {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rmultiStoreRepI<T>::readhdf5P(hid_t group_id, hid_t dataspace,
                                     hid_t dataset, hsize_t dimension,
                                     const char* name, frame_info &fi,
                                     entitySet &user_eset,
                                     hid_t xfer_plist_id) {
    warn(true) ;
  }
#endif

  template<class T>
  void rmultiStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace,
                                    hid_t dataset, hsize_t dimension,
                                    const char* name, IDENTITY_CONVERTER c,
                                    frame_info &fi, entitySet &eset) {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rmultiStoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace,
                                     hid_t dataset, hsize_t dimension,
                                     const char* name, IDENTITY_CONVERTER c,
                                     frame_info &fi, entitySet &eset,
                                     hid_t xfer_plist_id) {
    warn(true) ;
  }
#endif

  template<class T>
  void rmultiStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace,
                                    hid_t dataset, hsize_t dimension,
                                    const char* name, USER_DEFINED_CONVERTER c,
                                    frame_info &fi, entitySet &eset) {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rmultiStoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace,
                                     hid_t dataset, hsize_t dimension,
                                     const char* name, USER_DEFINED_CONVERTER c,
                                     frame_info &fi, entitySet &eset,
                                     hid_t xfer_plist_id) {
    warn(true) ;
  }
#endif

  template<class T>
  void rmultiStoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace,
                                     hid_t dataset, hsize_t dimension,
                                     const char* name,
                                     entitySet& usr_eset) const {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rmultiStoreRepI<T>::writehdf5P(hid_t group_id, hid_t dataspace,
                                      hid_t dataset, hsize_t dimension,
                                      const char* name, entitySet& usr_eset,
                                      hid_t xfer_plist_id) const {
    warn(true) ;
  }
#endif

  template<class T>
  void rmultiStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace,
                                     hid_t dataset, hsize_t dimension,
                                     const char* name, IDENTITY_CONVERTER c,
                                     const entitySet &eset) const {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rmultiStoreRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace,
                                      hid_t dataset, hsize_t dimension,
                                      const char* name, IDENTITY_CONVERTER c,
                                      const entitySet &en,
                                      hid_t xfer_plist_id) const {
    warn(true) ;
  }
#endif

  template<class T>
  void rmultiStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace,
                                     hid_t dataset, hsize_t dimension,
                                     const char* name, USER_DEFINED_CONVERTER c,
                                     const entitySet &eset) const {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T>
  void rmultiStoreRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace,
                                      hid_t dataset, hsize_t dimension,
                                      const char* name, USER_DEFINED_CONVERTER c,
                                      const entitySet &en,
                                      hid_t xfer_plist_id) const {
    warn(true) ;
  }
#endif

} // end of namespace Loci

#endif
