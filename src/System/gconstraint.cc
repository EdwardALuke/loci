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

#include <gconstraint.h>
#include <gmap.h>
#include <constraint.h>
#include <distribute_long.h>
#include <frame_info.h>
using std::cerr ;
using std::endl ;
using std::ostream ;
using std::istream ;

namespace Loci {

  //given send_split, return recv_split;
  //or given recv_split, return send_split
  //allow overlap between processes
  extern std::vector<gEntitySet> transposePtn(const std::vector<gEntitySet> &ptn, MPI_Comm comm);

                                       
  //**************************************************************************/

  gStoreRepP gConstraintRep::remap(const gMap &m) const {
    gEntitySet newConstraint = m.image(m.domain()&constraint_set) ;
    gConstraint r ;
    r = newConstraint ;
    r.set_domain_space(domain_space);
    return r.Rep() ;
  }
 
  //**************************************************************************/
  
  int gConstraintRep::pack_size(const gEntitySet &e)const {
    warn(true) ;
    return 0 ;
  }

  //**************************************************************************/
  
  void gConstraintRep::pack(void *ptr, int &loc, int size, const gEntitySet&e)const {
    warn(true) ;
  }
  
  //**************************************************************************/
 
  void gConstraintRep::unpack(const void *ptr, int &loc, int size) {
    warn(true);
  }
  
  //**************************************************************************/
  
  ostream &gConstraintRep::Print(ostream &s) const {
    s << constraint_set << endl ;
    return s ;
  }

  //**************************************************************************/

  DatatypeP gConstraintRep::getType()const {
    if(sizeof(int_type)==sizeof(gEntity))
      return DatatypeP(new AtomicType(INT)) ;
    else  return DatatypeP(new AtomicType(LONG)) ;
  }

  //**************************************************************************/
  
  istream &gConstraintRep::Input(istream &s) {
    gEntitySet e ;
    s >> e ;
    allocate(e) ;
    return s ;
  }

  //**************************************************************************/

  gStoreRepP gConstraintRep::
  split_redistribute(const std::vector<gEntitySet>& dom_ptn, MPI_Comm comm)const {
   
    gEntitySet dom = constraint_set ;
    gEntitySet old_all = g_all_collect_entitySet<gEntity>(dom, comm); 
    gEntitySet new_all ;
    for(size_t i=0;i<dom_ptn.size();++i)
      new_all += dom_ptn[i] ;
    gEntitySet out = dom - new_all ;

    
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;

    // get the new domain
    gEntitySet new_dom = old_all & dom_ptn[rank] ;
    new_dom += out ;

    gConstraint new_store;
    new_store = new_dom;
    new_store.set_domain_space(domain_space);
    return new_store.Rep() ;
  }

  //**************************************************************************/

  gStoreRepP gConstraintRep::
  redistribute(const std::vector<gEntitySet>& dom_split,
               MPI_Comm comm)const{
   
    std::vector<gEntitySet> recv_split = transposePtn(dom_split, comm);
    gEntitySet new_dom;
    for(size_t i=0;i<recv_split.size();++i)new_dom += recv_split[i];
    gEntitySet old_all = g_all_collect_entitySet<gEntity>(constraint_set, comm);
    new_dom = new_dom&old_all;
    gConstraint new_store;
    new_store = new_dom;
    new_store.set_domain_space(domain_space);
    return new_store.Rep() ;
  }
  
  //**************************************************************************/

  gStoreRepP gConstraintRep::
  redistribute(const std::vector<gEntitySet>& dom_split,
               const gMap& remap, MPI_Comm comm)const{
    gConstraint new_store;
    new_store = redistribute(dom_split);
    new_store.set_domain_space(domain_space);
    return new_store.remap(remap);
  }

  //**************************************************************************/

#ifdef COPY2STORE  
  storeRepP gConstraintRep::copy2store()const {
    constraint cs;
    cs = gconvert(constraint_set); 
    return cs.Rep(); 
  }
#endif
  
  //**************************************************************************/

  frame_info gConstraintRep::get_frame_info() const{
    warn(true) ; 
    frame_info fi ;
    return fi ;
  }

  //**************************************************************************/
 
  void gConstraintRep::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                const char* name, frame_info &fi, const gEntitySet &en){
    warn(true);
  }

  //**************************************************************************/  
  void gConstraintRep::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                 const char* name, const gEntitySet &en) const {
    warn(true);
  }

  //**************************************************************************/
}
