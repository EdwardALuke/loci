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

#ifndef LOCI_GPUMAPVEC_IMPL_H_
#define LOCI_GPUMAPVEC_IMPL_H_

#include <gpuMapVec_def.h>

namespace Loci {

  using std::cerr;
  using std::endl;

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::allocate(entitySet const & ptn) {
    if(alloc_id < 0) {
      alloc_id = getGPUStoreAllocateID();
    }

    GPUstoreAllocateData[alloc_id].template allocBasic<Array<Entity, M>>(ptn, 1);
    store_domain = GPUstoreAllocateData[alloc_id].allocset;
    base_ptr = (Array<Entity, M> *)GPUstoreAllocateData[alloc_id].base_ptr;
    base_ptr -= GPUstoreAllocateData[alloc_id].base_offset;
    dispatch_notify();
  }

  template<unsigned int M>
  inline gpuMapVecRepI<M>::~gpuMapVecRepI() {
    if(alloc_id >= 0) {
      GPUstoreAllocateData[alloc_id].template release<Array<Entity, M>>();
      releaseGPUStoreAllocateID(alloc_id);
      alloc_id = -1;
    }
  }

  template<unsigned int M>
  inline storeRep * gpuMapVecRepI<M>::new_store(entitySet const & p) const {
    return new gpuMapVecRepI<M>(p);
  }

  template<unsigned int M>
  inline storeRep * gpuMapVecRepI<M>::new_store(entitySet const & p, int const *) const {
    cerr << "gpuMapVecRepI<M>::new_store is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline storeRepP gpuMapVecRepI<M>::remap(dMap const & m) const {
    cerr << "gpuMapVecRepI<M>::remap is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline storeRepP gpuMapVecRepI<M>::MapRemap(dMap const & dm, dMap const & rm) const {
    cerr << "gpuMapVecRepI<M>::MapRemap is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::compose(dMap const & m, entitySet const & context) {
    cerr << "gpuMapVecRepI<M>::compose is not implemented" << endl;
    Loci::Abort();
  }

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::copy(storeRepP & st, entitySet const & context) {
    cerr << "gpuMapVecRepI<M>::copy is not implemented" << endl;
    Loci::Abort();
  }

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::gather(dMap const & m, storeRepP & st, entitySet const & context) {
    cerr << "gpuMapVecRepI<M>::gather is not implemented" << endl;
    Loci::Abort();
  }

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::scatter(dMap const & m, storeRepP & st, entitySet const & context) {
    cerr << "gpuMapVecRepI<M>::scatter is not implemented" << endl;
    Loci::Abort();
  }

  template<unsigned int M>
  inline int gpuMapVecRepI<M>::pack_size(entitySet const & e, entitySet & packed) {
    cerr << "gpuMapVecRepI<M>::pack_size is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline int gpuMapVecRepI<M>::pack_size(entitySet const & e) {
    cerr << "gpuMapVecRepI<M>::pack_size is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline int gpuMapVecRepI<M>::estimated_pack_size(entitySet const & e) {
    cerr << "gpuMapVecRepI<M>::estimated_pack_size is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::pack(void * ptr, int & loc, int & size, entitySet const & e) {
    cerr << "gpuMapVecRepI<M>::pack is not implemented" << endl;
    Loci::Abort();
  }

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::unpack(void * ptr, int & loc, int & size, sequence const & seq) {
    cerr << "gpuMapVecRepI<M>::unpack is not implemented" << endl;
    Loci::Abort();
  }

  template<unsigned int M>
  inline entitySet gpuMapVecRepI<M>::domain() const {
    return store_domain;
  }

  template<unsigned int M>
  inline entitySet gpuMapVecRepI<M>::image(entitySet const & domain) const {
    return defermap->image(domain);
  }

  template<unsigned int M>
  inline std::pair<entitySet, entitySet> gpuMapVecRepI<M>::preimage(entitySet const & codomain) const {
    return defermap->preimage(codomain);
  }

  template<unsigned int M>
  inline storeRepP gpuMapVecRepI<M>::get_map() {
    cerr << "gpuMapVecRepI<M>::get_map is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline std::ostream & gpuMapVecRepI<M>::Print(std::ostream & s) const {
    cerr << "gpuMapVecRepI<M>::Print should not be called for gpuMapVec" << endl;
    Loci::Abort();
    return s;
  }

  template<unsigned int M>
  inline std::istream & gpuMapVecRepI<M>::Input(std::istream & s) {
    cerr << "gpuMapVecRepI<M>::Input is not implemented" << endl;
    Loci::Abort();
    return s;
  }

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::readhdf5(
    hid_t group_id, hid_t daaspace, hid_t dataset, hsize_t dimension,
    char const * name, frame_info & fi, entitySet & en
  ) {
    cerr << "gpuMapVecRepI<M>::readhdf5 is not implemented" << endl;
    Loci::Abort();
  }

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::writehdf5(
    hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
    char const * name, entitySet & en
  ) const {
    cerr << "gpuMapVecRepI<M>::writehdf5 is not implemented" << endl;
    Loci::Abort();
  }

#ifdef H5_HAVE_PARALLEL
  template<unsigned int M>
  inline void gpuMapVecRepI<M>::readhdf5P(
    hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
    char const * name, frame_info & fi, entitySet & en, hid_t xfer_plist_id
  ) {
    cerr << "gpuMapVecRepI<M>::readhdf5P is not implemented" << endl;
    Loci::Abort();
  }

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::writehdf5P(
    hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
    char const * name, entitySet & en, hid_t xfer_plist_id
  ) const {
    cerr << "gpuMapVecRepI<M>::writehdf5P is not implemented" << endl;
    Loci::Abort();
  }
#endif

  template<unsigned int M>
  inline storeRepP gpuMapVecRepI<M>::expand(
    entitySet & out_of_dom, std::vector<entitySet> & init_ptn
  ) {
    cerr << "gpuMapVecRepI<M>::expand is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline storeRepP gpuMapVecRepI<M>::freeze() {
    cerr << "gpuMapVecRepI<M>::freeze is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline storeRepP gpuMapVecRepI<M>::thaw() {
    cerr << "gpuMapVecRepI<M>::thaw is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline DatatypeP gpuMapVecRepI<M>::getType() {
    cerr << "gpuMapVecRepI<M>::getType is not implemented" << endl;
    Loci::Abort();
    return 0;
  }

  template<unsigned int M>
  inline frame_info gpuMapVecRepI<M>::get_frame_info() {
    cerr << "gpuMapVecRepI<M>::get_frame_info is not implemented" << endl;
    Loci::Abort();
    return frame_info();
  }

  template<unsigned int M>
  inline void gpuMapVecRepI<M>::copyFrom(const storeRepP &fromMap, entitySet set) {
#ifdef USE_CUDA_RT
    int num_intervals = set.num_intervals();
    MapVec<M> m;
    m.setRep(fromMap);
    Array<Entity, M> * gpu_base_ptr = get_base_ptr();
    for(int i = 0; i < num_intervals; ++i) {
      int start = set[i].first;
      int end = set[i].second;
      int sz = end-start+1;

      cudaError_t err = cudaMemcpy(
        &gpu_base_ptr[start], &m[start], sizeof(Array<Entity, M>)*sz,
        cudaMemcpyHostToDevice
      );
      cudaDeviceSynchronize();
      if(err != cudaSuccess) {
        cerr << "cudaMemcpy failed in gpuMapVecRepI<M>::copyFrom" << endl;
        cerr << "error message: " << cudaGetErrorString(err) << endl;
        Loci::Abort();
      }
    }
#endif
  }

  template<unsigned int M>
  inline store_type gpuMapVecRepI<M>::RepType() const  {
    return GPUMAP;
  }

  template<unsigned int M>
  inline void gpuMapVec<M>::notification() {
    NPTR<gpuMapVecRepI<M>> p(Rep());
    if(p != 0) {
      base_ptr = p->get_base_ptr();
    }
    warn(p == 0);
  }

  template<unsigned int M>
  inline store_instance::instance_type const_gpuMapVec<M>::access() const {
    return READ_ONLY;
  }

  template<unsigned int M>
  inline void const_gpuMapVec<M>::notification() {
    NPTR<gpuMapVecRepI<M>> p(Rep());
    if(p != 0) {
      base_ptr = p->get_base_ptr();
    }
    warn(p == 0);
  }
}

#endif
