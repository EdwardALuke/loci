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
#include <fstream>
#include <algorithm>

#include <gmap.h>
#include <gmultimap.h>	
#include <hdf5_readwrite_long.h>
#include <gkeyspace.h>
#include <distribute_long.h>
#include <Map.h>
#include <partition.h>
#include <distributed_inverse_long.h>
using std::cerr ;
using std::endl ;
using std::ostream ;
using std::istream ;
using std::ofstream ;

namespace Loci {
  
  extern int GLOBAL_OR(int b, MPI_Comm comm); 
  extern ofstream debugout ;
  
  using std::pair ;
  using std::make_pair ;
  using std::vector ;
  using std::stable_sort ;


  namespace {
    inline  bool gmap_porder1(const std::pair<gEntity, gEntity> &i1, 
                              const std::pair<gEntity, gEntity> &i2) {
      return i1.first < i2.first ;
    }
    inline  bool gmap_porder2(const std::pair<gEntity, gEntity> &i1, 
                              const std::pair<gEntity, gEntity> &i2) {
      return i1.second < i2.second ;
    } 
  }
  
  gMapRepI::gMapRepI(){
    sorted= true;
    domain_space = 0;
    image_space = 0;
  }
  gMapRepI::~gMapRepI(){}
  void  gMapRepI::set_domain_space(gKeySpaceP space){domain_space = &(*space);}
  gKeySpaceP  gMapRepI::get_domain_space()const{return domain_space;}
  void  gMapRepI::set_image_space(gKeySpaceP space){image_space = &(*space);}
  gKeySpaceP  gMapRepI::get_image_space()const{return image_space;}

  void gMap::set_domain_space(gKeySpaceP space){static_cast<gMapRepP>(Rep())->set_domain_space(space);}
  gKeySpaceP gMap::get_domain_space()const{return static_cast<gMapRepP>(Rep())->get_domain_space() ;}
  void gMap::set_image_space(gKeySpaceP space){static_cast<gMapRepP>(Rep())->set_image_space(space);}
  gKeySpaceP gMap::get_image_space()const{return static_cast<gMapRepP>(Rep())->get_image_space();}
  
  void gMapRepI::local_sort(){
    std::sort(attrib_data.begin(), attrib_data.end(), gmap_porder1);
    sorted = true;
  }
  
  void gMapRepI::local_sort2(){
    std::stable_sort(attrib_data.begin(), attrib_data.end(), gmap_porder2);
    sorted = false;
  }
  
  //parSampleSort oot working, because balanceDisttribution will split one entity inside a multiMap
  // or a multiStore 
 
                         

  
  //**************************************************************************/
  gStoreRepP gMapRepI::expand(gEntitySet &dom, std::vector<gEntitySet> &ptn, MPI_Comm comm)const {
    int np;
    MPI_Comm_size(comm, &np);
    std::vector<gEntitySet> recv_request(np);
    for(int i = 0; i < np; ++i) {
      recv_request[i] = dom & ptn[i] ; //recv request
    }
    vector<gEntitySet> send_ptn = transposePtn(recv_request, comm);
    return redistribute(send_ptn, comm);
  }

  gStoreRepP gMapRepI::local_inverse()const {
    gMap result;
    for(const_iterator itr = begin(); itr != end(); itr++){
      result.insert(itr->second, itr->first);
    }
    result.local_sort();
    return result.Rep();
  }
  
 
  gStoreRepP gMapRepI::remap(const gMap &m) const 
  {
    fatal(!(sorted && m.sorted()));
    fatal(m.domain()-domain() != GEMPTY);
    gMap s ;
    
    const_iterator itr2 = m.begin();
    for(const_iterator itr1 = begin(); itr1 != end(); itr1++){
      while(itr2 != m.end() && itr2->first < itr1->first)itr2++ ;
      if(itr2 != m.end() && itr2->first == itr1->first) s.insert(itr2->second, itr1->second);
    }
    s.local_sort();
    s.set_domain_space(domain_space);
    s.set_image_space(image_space);
    return s.Rep();
  }
  
  void gMapRepI::inplace_remap(const gMap &m) 
  {
    fatal(!(sorted && m.sorted()));
    gMap::const_iterator itr2 = m.begin();
    for( iterator itr1 = begin(); itr1 != end(); itr1++){
      while(itr2 != m.end() && itr2->first < itr1->first)itr2++ ;
      if(itr2 != m.end() && itr2->first == itr1->first) itr1->first = itr2->second;
    }
    local_sort();
  }
  
  gStoreRepP gMapRepI::distributed_inverse(gEntitySet global_image,
                                           gEntitySet global_preimage,
                                           const std::vector<gEntitySet> &init_ptn) const{
    gMultiMap result;
    vector<pair<gEntity,gEntity> > input;
 
    for(const_iterator itr = begin(); itr != end(); itr++){
      if(global_preimage.inSet(itr->first)) input.push_back( pair<gEntity,gEntity>(itr->first, itr->second));
    }
 
    vector<gEntity> recv_store;
    distributed_inverse_map<gEntity>(recv_store,
                                     input, 
                                     init_ptn);


    gEntitySet local_input_image = global_image ;
    local_input_image &= init_ptn[MPI_rank] ;
    
    
    for(size_t i=0;i<recv_store.size();++i) {
      gEntity indx = recv_store[i++] ;
      gEntity refer = recv_store[i] ;
      if(local_input_image.inSet(indx))result.insert(indx, refer);
    }
    result.local_sort();
    result.set_vdom(local_input_image - result.domain());
    result.set_domain_space(domain_space);
    result.set_image_space(image_space);
    return result.Rep();
  }
 
  gStoreRepP
  gMapRepI::redistribute(const std::vector<gEntitySet>& dom_split,const gMap& remap,
                         MPI_Comm comm)const {
    gMap result;
    result = redistribute(dom_split, comm);
    result.set_domain_space(domain_space);
    result.set_image_space(image_space);
    return result.remap(remap);
  }
       
  gStoreRepP
  gMapRepI::redistribute(const std::vector<gEntitySet>& dom_split,
                         MPI_Comm comm)const {
    

    // compute the send_counts
    int np;
    MPI_Comm_size(comm, &np);
    std::vector<int> send_counts(np,0) ;
    for(int i=0;i<np;++i){
      send_counts[i] = pack_size(dom_split[i]) ;
    }
    // communicate to get the receive counts
    std::vector<int> recv_counts(np,0) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, comm) ;
    // then pack the buffer
    int tot_send_size = 0 ;
    int tot_recv_size = 0 ;
    for(int i=0;i<np;++i){
      tot_send_size += send_counts[i] ;
      tot_recv_size += recv_counts[i] ;
    }
    
  
    
    int *send_displs = new int[MPI_processes] ;
    int *recv_displs = new int[MPI_processes] ;
    send_displs[0] = 0 ;
    recv_displs[0] = 0 ;
    for(int i=1;i<np;++i){
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    }

  
    unsigned char* send_buffer = new unsigned char[tot_send_size] ;
    int loc_pack = 0 ;
    for(int i = 0; i <  np; ++i)
      pack(send_buffer, loc_pack, tot_send_size, dom_split[i]) ;
  
    // we now proceed to communicate the store contents
    unsigned char* recv_buffer = new unsigned char[tot_recv_size] ;
    MPI_Alltoallv(&send_buffer[0], &send_counts[0],
                  &send_displs[0], MPI_PACKED,
                  &recv_buffer[0], &recv_counts[0],
                  &recv_displs[0], MPI_PACKED, comm) ;

  

    gMap ns ;
    gStoreRepP srp = ns.Rep();
    // unpack
    loc_pack = 0 ;
    for(int i=0;i<np;++i) {
      srp->unpack(recv_buffer, loc_pack, tot_recv_size) ;
    }
    srp->local_sort();
    ns.set_domain_space(domain_space);
    ns.set_image_space(image_space);
    delete[] recv_displs ;
    delete[] send_displs ;
    delete[] send_buffer ;
    delete[] recv_buffer ;
    return srp ;
  }
  //split domain first, then redistribute
  gStoreRepP
  gMapRepI::split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                               MPI_Comm comm)const {
    gEntitySet dom = domain() ;
    // figure out how the domain is split to send to others
    int np ;
    MPI_Comm_size(comm, &np) ;
    fatal(np != dom_ptn.size()) ;
    std::vector<gEntitySet> dom_split(np) ;
    for(int i=0;i<np;++i) {
      dom_split[i] = dom & dom_ptn[i] ;
    }
    return redistribute(dom_split, comm);
  }
  
    
  //**************************************************************************/
  
  void gMapRepI::inplace_compose(const gMap &newmap,  MPI_Comm comm) 
  {
    //check if expand is needed
    gEntitySet old_dom = newmap.domain();
    gEntitySet dom = image();
    gEntitySet out_dom = dom -  old_dom;
    int need_expand = (out_dom == GEMPTY ? 0:1);
    bool global_need_expand = GLOBAL_OR(need_expand, comm);
    gMap m;
    if(global_need_expand){
      std::vector<gEntitySet> old_ptn = g_all_collect_vectors<gEntity>(old_dom, comm) ;
      m.setRep(newmap.expand(dom, old_ptn, comm));
    }else{
      m.setRep(newmap.Rep());
    }

    
    fatal(!(m.sorted()));
    local_sort2();
    gMap::const_iterator itr2 = m.begin();
    for( iterator itr1 = begin(); itr1 != end(); itr1++){
      while(itr2 != m.end() && itr2->first < itr1->second)itr2++ ;
      if(itr2 != m.end() && itr2->first == itr1->second) itr1->second = itr2->second;
    }
    local_sort();
  }
  
  //**************************************************************************/
  gStoreRepP gMapRepI::recompose(const gMap &newmap,  MPI_Comm comm)const 
  {
    //check if expand is needed
    gEntitySet old_dom = newmap.domain();
    gEntitySet dom = image();
    gEntitySet out_dom = dom -  old_dom;
    int need_expand = (out_dom == GEMPTY ? 0:1);
    bool global_need_expand = GLOBAL_OR(need_expand, comm);
    gMap m;
    if(global_need_expand){
      std::vector<gEntitySet> old_ptn = g_all_collect_vectors<gEntity>(old_dom, comm) ;
      m.setRep(newmap.expand(dom, old_ptn, comm));
    }else{
      m.setRep(newmap.Rep());
    }

    
    gMap result;
    fatal(!(m.sorted()));
    const_cast <gMapRepI*> (this)->local_sort2();
    gMap::const_iterator itr2 = m.begin();
    for( const_iterator itr1 = begin(); itr1 != end(); itr1++){
      while(itr2 != m.end() && itr2->first < itr1->second)itr2++ ;
      if(itr2 != m.end() && itr2->first == itr1->second) result.insert(itr1->first, itr2->second);
    }
    const_cast <gMapRepI*> (this)->local_sort();
    result.set_domain_space(domain_space);
    result.set_image_space(image_space); result.local_sort();
    return result.Rep();
  }

  //**************************************************************************/
  gStoreRepP gMapRepI::recompose(const gMultiMap &newmap,  MPI_Comm comm)const 
  {
    //check if expand is needed
    gEntitySet old_dom = newmap.domain();
    gEntitySet dom = image();
    gEntitySet out_dom = dom -  old_dom;
    int need_expand = (out_dom == GEMPTY ? 0:1);
    bool global_need_expand = GLOBAL_OR(need_expand, comm);
    gMultiMap m;
    if(global_need_expand){
      std::vector<gEntitySet> old_ptn = g_all_collect_vectors<gEntity>(old_dom, comm) ;
      m.setRep(newmap.expand(dom, old_ptn, comm));
    }else{
      m.setRep(newmap.Rep());
    }

    
    gMultiMap result;
    fatal(!(m.sorted()));
    const_cast <gMapRepI*> (this)->local_sort2();
    gMultiMap::const_iterator itr2 = m.begin();
    for( const_iterator itr1 = begin(); itr1 != end(); itr1++){
      while(itr2 != m.end() && itr2->first < itr1->second)itr2++ ;
      while(itr2 != m.end() && itr2->first == itr1->second){
        result.insert(itr1->first, itr2->second);
        itr2++;
      }
    }
    const_cast <gMapRepI*> (this)->local_sort();
    result.set_domain_space(domain_space);
    result.set_image_space(image_space);
    result.local_sort();
    return result.Rep();
  }
 
  
  //**************************************************************************/
  gStoreRepP gMapRepI::recompose(gStoreRepP &newmap,  MPI_Comm comm)const {
    if(newmap->RepType()==GMAP){
      gMap m(newmap);
      return recompose(m, comm);
    }else  if(newmap->RepType()==GMULTIMAP){
      gMultiMap m(newmap);
      return recompose(m, comm);
    }else{
      cerr << "ERROR: uncognized RepType() in recompose() method" << endl;
      return gStoreRepP(0);
    }
  }
    
  //**************************************************************************/
  
  int gMapRepI::pack_size(const gEntitySet &eset)const {
    fatal(!sorted);
    int result = 0;
    const_iterator itr1 = attrib_data.begin();
    for( gEntitySet :: const_iterator ci = eset.begin(); ci != eset.end(); ci++){
      while( itr1 != attrib_data.end()  && itr1->first < *ci )itr1++;
      while( itr1 != attrib_data.end() && itr1->first == *ci ){
        result += 2*sizeof(gEntity);
        itr1++;
      }
    }
    return result;
  }
  

 
  
  //**************************************************************************/
  
  void gMapRepI::pack(void *outbuf, int &position, int outsize, const gEntitySet &eset) const
  {
    warn(!sorted);
     
    const_iterator itr1 = attrib_data.begin();
    for(  gEntitySet :: const_iterator ci = eset.begin(); ci != eset.end(); ci++){
      while(itr1 != attrib_data.end() && itr1->first < *ci)itr1++;
      while(itr1 != attrib_data.end() && itr1->first == *ci){
        MPI_Pack( const_cast<gEntity*>(&(itr1->first)), sizeof(gEntity), MPI_BYTE, outbuf,outsize, 
                  &position, MPI_COMM_WORLD) ;
        MPI_Pack( const_cast<gEntity*>(&(itr1->second)), sizeof(gEntity), MPI_BYTE, outbuf,outsize, 
                  &position, MPI_COMM_WORLD) ;
        itr1++;
      }
    }
  }
  
    
  //**************************************************************************/

  void gMapRepI::unpack(const void *inbuf, int &position, int insize) 
  {
    gEntity ei;
    gEntity val;
    while(position != insize){
      MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &ei,
                  sizeof(gEntity), MPI_BYTE, MPI_COMM_WORLD) ;
      MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &val,
                  sizeof(gEntity), MPI_BYTE, MPI_COMM_WORLD) ;
      insert(ei, val);
    }
  }
   
  gEntity gMapRepI::image( gEntity iset) const 
  {
    
    gEntity codomain = 0 ;
    std::pair<gEntity, gEntity> p = std::make_pair(iset, gEntity(0));
    const_iterator itr = std::lower_bound(begin(), end(), p, gmap_porder1);
    if(itr != end()){
      codomain = itr->second;
    }
    return codomain ;
  }
  //**************************************************************************/
  
  gEntitySet gMapRepI::image(const gEntitySet &iset) const 
  {
    
    gEntitySet codomain ;
    GFORALL(iset, ei){
      codomain += image(ei);
    }ENDGFORALL;
    return codomain ;
  }
  
  gEntitySet gMapRepI::image() const 
  {
    gEntitySet codomain ;
    const_iterator itr1 = begin();
    while(itr1 != end()){
      codomain += itr1->second;
      itr1++;
    }
    return codomain ;
  }

  
  
  //**************************************************************************/
  
  pair<gEntitySet,gEntitySet>
  gMapRepI::preimage(const gEntitySet &codomain) const 
  {
   
    const_cast<gMapRepI&>(*this).local_sort2();
    gEntitySet domain ;
    const_iterator itr = begin();
    GFORALL(codomain, ei){
      
      while(itr != end() && itr->second < ei) itr++;
      while(itr != end() && itr->second == ei){
        domain += itr->first;
        itr++;
      }
    }ENDGFORALL;
    const_cast<gMapRepI&>(*this).local_sort();
    return make_pair(domain,domain);
  }
  
  

  //**************************************************************************/

  ostream &gMapRepI::Print(ostream &s) const 
  {
    s << '{' << domain() << std::endl ;
    for(const_iterator itr = begin(); itr!= end(); itr++){
      s << itr->second << std::endl ;
    }
    s << '}' << std::endl ;
    return s ;
  }

  //**************************************************************************/

  istream &gMapRepI::Input(istream &s) 
  {
    gEntitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading gStore" << endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
        
    GFORALL(e,ii) {
      gEntity img;
      s >> img ;
      insert(ii, img);
    } ENDGFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading gStore" << endl ;
      s.putback(ch) ;
    }
    local_sort();
    return s ;
  }

  DatatypeP gMapRepI::getType()const {
    return DatatypeP(new AtomicType(INT)) ;
  }

  //copy gMap to traditional Loci Map 
#ifdef COPY2STORE  
  storeRepP gMapRepI::copy2store() const{
    // if(!sorted)throw(bad_access);
    entitySet dom = gconvert(domain());
    fatal(size() != dom.size()) ;
    Map s;
    s.allocate(dom);

    for(const_iterator itr = begin(); itr!= end(); itr++){
      s[itr->first] = int(itr->second);
    } 
    return s.Rep();
  }
#endif
    
  frame_info gMapRepI::get_frame_info()const {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }

  //**************************************************************************/

  void gMapRepI::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                          const char* name, frame_info &fi, const gEntitySet &usr_eset)  {
    warn(true) ;
  } 

  //**************************************************************************/
  
  void gMapRepI::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                           const char* name, const gEntitySet& usr_eset) const
  {
    warn(true) ;
    
  } 

  
 

}
