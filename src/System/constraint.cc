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
#include <istream>
#include <ostream>
#include <iostream>

#include <constraint.h>
#include <Map.h>
#include <DMap.h>

using std::cerr ;
using std::endl ;
using std::ostream ;
using std::istream ;

namespace Loci {

  constraintRep::constraintRep()
  {
  }

  constraintRep::constraintRep(const entitySet &p)
  {
    constraint_set = p ;
  }

  constraintRep::~constraintRep()
  {
  }

  void constraintRep::allocate(const entitySet &p)
  {
    constraint_set = p ;
    dispatch_notify() ;
  }

  void constraintRep::shift(int_type offset)
  {
    constraint_set >>= offset ;
    dispatch_notify() ;
  }

  storeRep *constraintRep::new_store(const entitySet &p) const
  {
    return new constraintRep(p) ;
  }
  storeRep *constraintRep::new_store(const entitySet &p, const int* cnt) const
  {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a constraint " << endl ;
    return sp ;
  }

  storeRepP constraintRep::remap(const dMap &m) const {
    
    entitySet newconstraint = ~EMPTY ;
    if(constraint_set != ~EMPTY)
      newconstraint = m.image(m.domain()&constraint_set) ;

    constraint r ;
    r = newconstraint ;
    return r.Rep() ;
  }

  storeRepP constraintRep::freeze() {
    return getRep() ;
  }

  storeRepP constraintRep::thaw() {
    return getRep() ;
  }

  void constraintRep::copy(storeRepP &st, const entitySet &context) {
    constraint cs(st) ;
    entitySet sent,tent ;
    sent = *cs ;
    tent = constraint_set ;
    tent -= context ;
    tent += sent & context ;
    constraint_set = tent ;
    dispatch_notify() ;
  }

  void constraintRep::gather(const dMap &m, storeRepP &st,
                             const entitySet &context) {
    constraint cs(st) ;
    entitySet tent = constraint_set ;
    tent -= context ;
    entitySet img = *cs ;
    FORALL(context,i) {
      if(img.inSet(m[i]))
        tent+=i ;
    } ENDFORALL ;
    constraint_set = tent ;
    dispatch_notify() ;
  }

  void constraintRep::scatter(const dMap &m, storeRepP &st,
                              const entitySet &context) {
    constraint cs(st) ;
    entitySet map_image = m.image(context) ;
    entitySet tent = constraint_set ;
    tent -= map_image ;
    entitySet img = *cs ;
    tent += m.image(context&img) ;
    constraint_set = tent ;
    dispatch_notify() ;
  }
  
  int constraintRep::pack_size(const entitySet &e) {
    warn(true) ;
    return 0 ;
  }
 int constraintRep::estimated_pack_size(const entitySet &e) {
    warn(true) ;
    return 0 ;
  }

  

  int constraintRep::
  pack_size(const entitySet& e, entitySet& packed) {
    warn(true) ;
    return 0 ;
  }
  
  void constraintRep::pack(void *ptr, int &loc, int &size, const entitySet&e) {
    warn(true) ;
  }
  
  void constraintRep::unpack(void *ptr, int &loc, int &size, const sequence &seq) {
    warn(true);
  }
  
  store_type constraintRep::RepType() const
  {
    return CONSTRAINT ;
  }

  entitySet constraintRep::domain() const {
    return constraint_set ;
  }

  ostream &constraintRep::Print(ostream &s) const {
    s << constraint_set << endl ;
    return s ;
  }
  DatatypeP constraintRep::getType() {
    return DatatypeP(new AtomicType(INT)) ;
  }

  frame_info constraintRep::get_frame_info() {
    warn(true) ; 
    frame_info fi ;
    return fi ;
  }
  
  void constraintRep::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en){
    warn(true) ;

  }

  void constraintRep::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const{
    warn(true) ;
  }
  
#ifdef H5_HAVE_PARALLEL
  void constraintRep::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en, hid_t xfer_plist_id){
    warn(true) ;
  }
  
  void constraintRep::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist_id) const{
   warn(true) ;
  }
#endif
  
  istream &constraintRep::Input(istream &s) {
    entitySet e ;
    s >> e ;
    allocate(e) ;
    return s ;
  }

  constraint::constraint()
  {
    setRep(new constraintType) ;
  }
    
  constraint::constraint(const constraint &var)
  {
      *data = *(var.data) ;
  }

  constraint::~constraint()
  {
  }

  void constraint::notification()
  {
    NPTR<constraintType> p(Rep());
    if(p!=0)
      data = p->get_constraint() ;
    warn(p==0);
  }

  // a duplicate of the constraint class
  Constraint::Constraint()
  {
    setRep(new constraintType) ;
  }
    
  Constraint::Constraint(const Constraint &var)
  {
    *data = *(var.data) ;
  }

  Constraint::~Constraint()
  {}

  void Constraint::notification()
  {
    NPTR<constraintType> p(Rep());
    if(p!=0)
      data = p->get_constraint() ;
    warn(p==0);
  }

}
