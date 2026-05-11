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
#include <store_rep.h>

#include <Map.h>

using std::istream ;
using std::ostream ;

namespace Loci {

  std::vector<storeAllocateInfo> storeAllocateData ;
  std::vector<int> storeAllocateFreeList ;

  int getStoreAllocateID() {
    // allocate slot in storeAllocateData
    int id = storeAllocateData.size() ;
    if(!storeAllocateFreeList.empty()) {
      id = storeAllocateFreeList.back() ;
      storeAllocateFreeList.pop_back() ;
    } else {
      storeAllocateData.push_back(storeAllocateInfo()) ;
    }
    storeAllocateData[id].alloc_ptr1 = 0 ;
    storeAllocateData[id].alloc_ptr2 = 0 ;
    storeAllocateData[id].base_ptr = 0 ;
    storeAllocateData[id].base_offset = 0 ;
    storeAllocateData[id].size = 0 ;
    storeAllocateData[id].allocated_size = 0 ;
    storeAllocateData[id].allocated = true ;
    storeAllocateData[id].allocset = EMPTY ;
    return id ;
  }
    
  void releaseStoreAllocateID(int id) {
    storeAllocateData[id].alloc_ptr1 = 0 ;
    storeAllocateData[id].alloc_ptr2 = 0 ;
    storeAllocateData[id].base_ptr = 0 ;
    storeAllocateData[id].base_offset = 0 ;
    storeAllocateData[id].size = 0 ;
    storeAllocateData[id].allocated_size = 0 ;
    storeAllocateData[id].allocated = false ;
    storeAllocateData[id].allocset = EMPTY ;
    storeAllocateFreeList.push_back(id) ;
  }

  
  storeRep::~storeRep() {}
  void storeRep::set_elem_size(int sz) {} //{ warn(true) ; }

  storeRepP storeRep::getRep() { return storeRepP(this) ; }
  storeRepP storeRep::getRep() const { return storeRepP(const_cast<storeRep *>(this)) ; }
  store_instance::instance_type store_instance::access() const {
    return READ_WRITE ;
  }

  store_ref::~store_ref() {}

  int store_ref::get_alloc_id() const { return Rep()->get_alloc_id() ; }
  
  void store_ref::allocate(const entitySet &ptn) {
    Rep()->allocate(ptn) ;
  }

  void store_ref::shift(int_type offset) {
    Rep()->shift(offset) ;
  }

  storeRep *store_ref::new_store(const entitySet &p) const  {
    return Rep()->new_store(p) ;
  }
  storeRep *store_ref::new_store(const entitySet &p, const int* cnt) const  {
    return Rep()->new_store(p, cnt) ;
  }
  store_type store_ref::RepType() const {
    return Rep()->RepType() ;
  }

  storeRepP store_ref::remap(const dMap &m) const {
    return Rep()->remap(m) ;
  }

  storeRepP store_ref::freeze() {
    return Rep()->freeze() ;
  }

  storeRepP store_ref::thaw() {
    return Rep()->thaw() ;
  }

  void store_ref::copy(storeRepP &st, const entitySet &context) {
    Rep()->copy(st,context) ;
  }

  void store_ref::fast_copy(storeRepP &st, const entitySet &context) {
    Rep()->fast_copy(st,context) ;
  }

  void store_ref::gather(const dMap &m,storeRepP &st,
                         const entitySet &context) {
    Rep()->gather(m,st,context) ;
  }

  void store_ref::scatter(const dMap &m, storeRepP &st,
                          const entitySet &context) {
    Rep()->scatter(m,st,context) ;
  }  
  int store_ref::pack_size(const entitySet &e )  {
    return(Rep()->pack_size(e)) ;
  }
  int store_ref::estimated_pack_size(const entitySet &e )  {
    return(Rep()->estimated_pack_size(e)) ;
  }
  int store_ref::pack_size(const entitySet& e, entitySet& packed) {
    return Rep()->pack_size(e,packed) ;
  }
  void store_ref::pack(void * ptr, int &loc, int &size, const entitySet &e )  {
    Rep()->pack(ptr, loc, size, e) ;
  }
  void store_ref::unpack(void * ptr, int &loc, int &size, const sequence &seq )  {
     Rep()->unpack(ptr, loc, size, seq) ;
   }

  ostream &store_ref::Print(ostream &s) const {
    return Rep()->Print(s) ;
  }

  istream &store_ref::Input(istream &s) {
    return Rep()->Input(s) ;
  }

  entitySet store_ref::domain() const {
    return Rep()->domain() ;
  }

  storeRepP store_ref::getRep() {
    return Rep()->getRep() ;
  }

  storeRepP store_ref::getRep() const {
    return Rep()->getRep() ;
  }
    
  void store_ref::notification() {
    dispatch_notify() ;
  }

  void store_ref::
  erase(const entitySet& rm) { Rep()->erase(rm) ;}
  
  void store_ref::
  invalidate(const entitySet& valid) { Rep()->invalidate(valid) ;}
  
  void store_ref::
  guarantee_domain(const entitySet& include)
  { Rep()->guarantee_domain(include) ;}
  
  storeRepP store_ref::
  redistribute(const std::vector<entitySet>& dom_ptn,MPI_Comm comm)
  { return Rep()->redistribute(dom_ptn,comm) ;}
  
  storeRepP store_ref::
  redistribute(const std::vector<entitySet>& dom_ptn,
               const dMap& remap, MPI_Comm comm)
  { return Rep()->redistribute(dom_ptn,remap,comm) ;}
  
}
