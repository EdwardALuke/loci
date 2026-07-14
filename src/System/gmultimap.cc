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

#include <store.h>
#include <multiMap.h>
#include <gmultimap.h>
#include <hdf5_readwrite_long.h>
#include <distribute_long.h>
#include <distributed_inverse_long.h>
#include <field_sort.h>
#include <fact_db.h>
using std::cerr ;
using std::endl ;
using std::ostream ;
using std::istream ;
using std::ofstream ;

namespace Loci {
 
  extern int GLOBAL_OR(int b, MPI_Comm comm);

  //given send_split, return recv_split;
  //or given recv_split, return send_split
  //allow overlap between processes
  extern std::vector<gEntitySet> transposePtn(const std::vector<gEntitySet> &ptn, MPI_Comm comm);
  
  extern ofstream debugout ;
  namespace {
  
   
    
  }
  using std::pair ;
  using std::make_pair ;
  using std::vector ;
  using std::sort ;

  void gMultiMapRepI::local_sort(){
    std::stable_sort(attrib_data.begin(), attrib_data.end(), fieldSort1<gEntity>);
    sorted = true;
  }
  
  
  
  //**************************************************************************/
  gStoreRepP gMultiMapRepI::expand(gEntitySet &dom, std::vector<gEntitySet> &ptn, MPI_Comm comm)const {
    int np;
    MPI_Comm_size(comm, &np);
    std::vector<gEntitySet> recv_request(np);
    for(int i = 0; i < np; ++i) {
      recv_request[i] = dom & ptn[i] ; //recv request
    }
    vector<gEntitySet> send_ptn = transposePtn(recv_request,comm);
    return redistribute(send_ptn, comm);
  }

  gStoreRepP gMultiMapRepI::local_inverse()const {
    gMultiMap result;
    for(const_iterator itr = begin(); itr != end(); itr++){
      result.insert(itr->second, itr->first);
    }
    result.local_sort();
    return result.Rep();
  }
  

  gStoreRepP gMultiMapRepI::remap(const gMap &m) const 
  {
    fatal(!(sorted && m.sorted()));
    fatal(m.domain()-domain() != GEMPTY);
    gMultiMap s ;
    
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
  
  void gMultiMapRepI::inplace_remap(const gMap &m) 
  {
    fatal(!(sorted && m.sorted()));
    gMap::const_iterator itr2 = m.begin();
    for( iterator itr1 = begin(); itr1 != end(); itr1++){
      while(itr2 != m.end() && itr2->first < itr1->first)itr2++ ;
      if(itr2 != m.end() && itr2->first == itr1->first) itr1->first = itr2->second;
    }
    local_sort();
  }

 
  
  gStoreRepP
  gMultiMapRepI::redistribute(const std::vector<gEntitySet>& dom_split,const gMap& remap,
                              MPI_Comm comm)const {
    gMultiMap result;
    result = redistribute(dom_split, comm);
    result.set_domain_space(domain_space);
    result.set_image_space(image_space);
    return result.remap(remap);
  }
       
  gStoreRepP
  gMultiMapRepI::redistribute(const std::vector<gEntitySet>& dom_split,
                              MPI_Comm comm)const {
    
    fatal(!sorted);
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
    
  
    
    int *send_displs = new int[np] ;
    int *recv_displs = new int[np] ;
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

  

    gMultiMap ns ;
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
  gMultiMapRepI::split_redistribute(const std::vector<gEntitySet>& dom_ptn,
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
  
  int gMultiMapRepI::pack_size(const gEntitySet &eset)const {
   
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
  
  void gMultiMapRepI::pack(void *outbuf, int &position, int outsize, const gEntitySet &eset) const
  {
    fatal(!sorted);
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

  void gMultiMapRepI::unpack(const void *inbuf, int &position, int insize) 
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
  
  
  
 
  gEntitySet gMultiMapRepI::image( gEntity iset) const 
  {
    fatal(!sorted);
    gEntitySet codomain ;
    std::pair<gEntity, gEntity> p(iset, gEntity(0));
    std::pair<const_iterator, const_iterator> range = std::equal_range(begin(), end(), p, fieldSort1<gEntity>);
    if(range.first != range.second){
      for(const_iterator itr= range.first; itr != range.second; itr++) codomain += itr->second;
    }
    return codomain ;
  }
  
  //or use count function?
  int gMultiMapRepI::num_elems( gEntity iset) const 
  {
    fatal(!sorted);
    std::pair<gEntity, gEntity> p(iset, gEntity(0));
    std::pair<const_iterator, const_iterator> range = std::equal_range(begin(), end(), p, fieldSort1<gEntity>);
    return std::distance(range.first, range.second);
  }
  
  std::pair<gMultiMapRepI::const_iterator, gMultiMapRepI::const_iterator> gMultiMapRepI::range(gEntity iset)const{
    fatal(!sorted);
    std::pair<gEntity, gEntity> p(iset, gEntity(0));
    std::pair<const_iterator, const_iterator> range = std::equal_range(begin(), end(), p, fieldSort1<gEntity>);
    return range;
  }
  
  //**************************************************************************/
  
  gEntitySet gMultiMapRepI::image(const gEntitySet &iset) const 
  {
    fatal(!sorted);
    gEntitySet codomain ;
    GFORALL(iset, ei){
      codomain += image(ei);
    }ENDGFORALL;
    return codomain ;
  }
  
  gEntitySet gMultiMapRepI::image() const 
  {
    fatal(!sorted);
    gEntitySet codomain ;
    const_iterator itr1 = begin();
    while(itr1 != end()){
      codomain += itr1->second;
      itr1++;
    }
    return codomain ;
  }

  
  

  
  //**************************************************************************/
  
  gStoreRepP gMultiMapRepI::get_map() const 
  {
      
    gMultiMap result ;
    result.set_domain_space(domain_space);
    result.set_image_space(image_space);
    result.set_vdom(vdom);
    for( const_iterator itr = attrib_data.begin(); itr != attrib_data.end(); itr++){
      result.insert(itr->first, itr->second);
    }
    return result.Rep() ;
  }

  //**************************************************************************/

  ostream &gMultiMapRepI::Print(ostream &s) const 
  {
    s << '{' << domain() << endl ;
    s << vdom <<endl;
    
    for( const_iterator itr = attrib_data.begin(); itr != attrib_data.end(); itr++){
      s << itr->first <<' '<< itr->second << endl;
    }
    s << '}' << endl ;
    return s ;
  }

  //**************************************************************************/

  istream &gMultiMapRepI::Input(istream &s) 
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
    s >> vdom;

    gEntity val1 = 0;
    gEntity val2 = 0;
    GFORALL(e,ii) {
      s >> val1 ;
      s>> val2;
      insert(val1, val2);
    } ENDGFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading gStore" << endl ;
      s.putback(ch) ;
    }
    local_sort();
    return s ;
  }

  DatatypeP gMultiMapRepI::getType()const {
    return DatatypeP(new AtomicType(INT)) ;
  }
 
  frame_info gMultiMapRepI::get_frame_info()const {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  //**************************************************************************/

  void gMultiMapRepI::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                               const char* name, frame_info &fi, const gEntitySet &usr_eset)  {
    warn(true) ;
     
  } 

  //**************************************************************************/

  void gMultiMapRepI::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                const char* name, const gEntitySet& usr_eset) const
  {
    warn(true) ;
  } 
  
  //**************************************************************************/

  void gMultiMapRepI::inplace_compose(const gMap &newmap,  MPI_Comm comm) 
  {
    fatal(!sorted);
    //check if expand is needed
    gEntitySet old_dom = newmap.domain();
    gEntitySet dom = image();
    
    gEntitySet out_dom = dom -  old_dom;
    int need_expand = (out_dom == GEMPTY ? 0:1);
    int global_need_expand = GLOBAL_OR(need_expand, comm);
    gMap expanded_map;
    if(global_need_expand){
      std::vector<gEntitySet> old_ptn = g_all_collect_vectors<gEntity>(old_dom, comm) ;
      expanded_map.setRep(newmap.expand(dom, old_ptn, comm));
    }else{
      expanded_map.setRep(newmap.Rep());
    }

    //copy attrib_datat into a vector
    vector<ms_comp<gEntity> > temp_vec;
    {
      short ind = 0;
      const_iterator previous = attrib_data.begin();
      temp_vec.push_back(ms_comp<gEntity>(previous->first, previous->second, ind));
      const_iterator itr = attrib_data.begin();
      itr++;
      for(; itr != attrib_data.end(); itr++){
        previous = itr;
        previous--;
        if(itr->first == previous->first){
          ind++;
          temp_vec.push_back(ms_comp<gEntity>(itr->first, itr->second, ind));
        }else{
          ind = 0;
          temp_vec.push_back(ms_comp<gEntity>(itr->first, itr->second, ind));
        }
      }
    }
    //sort temp_vec according to the image field
    std::sort(temp_vec.begin(), temp_vec.end(), field_sort_img);
    
    //equalJoin and sort according to domain field and index field 
    {
      int i = 0;
      for(gMap::const_iterator itr2 = expanded_map.begin(); itr2!= expanded_map.end(); itr2++){
        while(temp_vec[i].img < itr2->first)i++ ;
        while(itr2->first == temp_vec[i].img){
          temp_vec[i].img =  itr2->second;
          i++;
        }
      }
      std::sort(temp_vec.begin(), temp_vec.end(), field_sort_dom<gEntity>);
    }
    
    //copy back
    {
      iterator itr = begin();
      for(unsigned int i = 0; i < temp_vec.size(); i++){
        itr->second = temp_vec[i].img;
        itr++;
      }
    }
  }


  //**************************************************************************/
  
  gStoreRepP gMultiMapRepI::recompose(const gMap &newmap,  MPI_Comm comm)const 
  {
    fatal(!sorted);
    //check if expand is needed
    gEntitySet old_dom = newmap.domain();
    gEntitySet dom = image();
    
    gEntitySet out_dom = dom -  old_dom;
    int need_expand = (out_dom == GEMPTY ? 0:1);
    bool global_need_expand = GLOBAL_OR(need_expand, comm);
    gMap expanded_map;
    if(global_need_expand){
      std::vector<gEntitySet> old_ptn = g_all_collect_vectors<gEntity>(old_dom, comm) ;
      expanded_map.setRep(newmap.expand(dom, old_ptn, comm));
    }else{
      expanded_map.setRep(newmap.Rep());
    }

    //copy attrib_datat into a vector
    vector<ms_comp<gEntity> > temp_vec;
    {
      short ind = 0;
      const_iterator previous = attrib_data.begin();
      temp_vec.push_back(ms_comp<gEntity>(previous->first, previous->second, ind));
      const_iterator itr = attrib_data.begin();
      itr++;
      for(; itr != attrib_data.end(); itr++){
        previous = itr;
        previous--;
        if(itr->first == previous->first){
          ind++;
          temp_vec.push_back(ms_comp<gEntity>(itr->first, itr->second, ind));
        }else{
          ind = 0;
          temp_vec.push_back(ms_comp<gEntity>(itr->first, itr->second, ind));
        }
      }
    }
    //sort temp_vec according to the image field
    std::sort(temp_vec.begin(), temp_vec.end(), field_sort_img);
    
    //equalJoin and sort according to domain field and index field 
    {
      int i = 0;
      for(gMap::const_iterator itr2 = expanded_map.begin(); itr2!= expanded_map.end(); itr2++){
        while(temp_vec[i].img < itr2->first)i++ ;
        while(itr2->first == temp_vec[i].img){
          temp_vec[i].img =  itr2->second;
          i++;
        }
      }
      std::sort(temp_vec.begin(), temp_vec.end(), field_sort_dom<gEntity>);
    }
    //copy to result
    gMultiMap result;
    for(unsigned int i = 0; i < temp_vec.size(); i++){
      result.insert(temp_vec[i].dom, temp_vec[i].img);
    }
    
    //result.local_sort();
    result.set_domain_space(domain_space);
    result.set_image_space(image_space);
    return result.Rep();
  }


  //**************************************************************************/
  
  gStoreRepP gMultiMapRepI::recompose(const gMultiMap &newmap,  MPI_Comm comm)const 
  {
    fatal(!sorted);
    //check if expand is needed
    gEntitySet old_dom = newmap.domain();
    gEntitySet dom = image();
    
    gEntitySet out_dom = dom -  old_dom;
    int need_expand = (out_dom == GEMPTY ? 0:1);
    bool global_need_expand = GLOBAL_OR(need_expand, comm);
    gMap expanded_map;
    if(global_need_expand){
      std::vector<gEntitySet> old_ptn = g_all_collect_vectors<gEntity>(old_dom, comm) ;
      expanded_map.setRep(newmap.expand(dom, old_ptn, comm));
    }else{
      expanded_map.setRep(newmap.Rep());
    }

    //copy attrib_data into a vector
    vector<pair<gEntity, gEntity> > temp_vec;
    for(const_iterator itr = begin(); itr != end(); itr++){
      temp_vec.push_back(make_pair(itr->first, itr->second));
    }
    //sort temp_vec according to the image field
    std::sort(temp_vec.begin(), temp_vec.end(), fieldSort2);
    
    gMultiMap result;
    gEntity previous = temp_vec[0].second;
    gMultiMap::const_iterator itr_begin = expanded_map.begin();
    gMultiMap::const_iterator itr_end = itr_begin;
    while(itr_end != expanded_map.end() && itr_end->first == previous){
      result.insert(temp_vec[0].first, itr_end->second);
      itr_end++;
    }
    for(unsigned int i = 1; i<temp_vec.size(); i++){
      previous =  temp_vec[i-1].second;
      if(temp_vec[i].second == previous){
        for( gMultiMap::const_iterator itr =itr_begin; itr!=itr_end; itr++){
          result.insert(temp_vec[i].first, itr->second);
        }
      }else{
        itr_begin = itr_end;
        while(itr_end != expanded_map.end() && itr_end->first ==temp_vec[i].second ){
          result.insert(temp_vec[i].first, itr_end->second);
          itr_end++;
        }
      }
    }
    
   
    result.remove_duplication();
    result.set_domain_space(get_domain_space());
    result.set_image_space(gMapRepP(newmap.Rep())->get_image_space());
   
    return result.Rep();
  }

  gStoreRepP gMultiMapRepI::recompose( gStoreRepP &newmap,  MPI_Comm comm)const {
    
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
  std::pair<gEntitySet,gEntitySet>
  gMultiMapRepI::preimage(const gEntitySet &codomain) const {
    fatal(!sorted);
    
    gEntitySet store_domain = domain();
    gEntitySet  domaini,domainu ;
    GFORALL(store_domain,i) {
      std::pair<gEntity, gEntity> p(i, gEntity(0));
      std::pair<const_iterator, const_iterator> range = std::equal_range(begin(), end(), p, fieldSort1<gEntity>);
      if(range.first != range.second){
        bool vali = true ;
        bool valu = false;
        for(const_iterator itr= range.first; itr != range.second; itr++){
          bool in_set = codomain.inSet(itr->second) ;
          vali = vali && in_set ;
          valu = valu || in_set ;
        }
        if(vali)
          domaini += i ;
        if(valu)
          domainu += i ;
      }
    } ENDGFORALL ;
     
    return  make_pair(domaini,domainu) ;
  }
  
 
  
  gStoreRepP gMultiMapRepI::distributed_inverse(gEntitySet global_image,
                                                gEntitySet global_preimage,
                                                const std::vector<gEntitySet> &init_ptn) const{
    fatal(!sorted); 
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
    local_input_image &= init_ptn[get_exec_rank()] ;
    
    
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

  void gMultiMapRepI::remove_duplication(){
    std::sort(attrib_data.begin(), attrib_data.end(), fieldSort1_unique<gEntity>);
    vector<pair<gEntity,gEntity> >::iterator uend ;
    uend = unique(attrib_data.begin(), attrib_data.end());
    attrib_data.erase(uend, attrib_data.end());
  }
  
  void gMultiMap::remove_duplication(){
    CPTR<MapType> p(Rep()) ;
    if(p != 0)p->remove_duplication();
    warn(p==0);
  }



  
  //copy gMultiMap to traditional Loci multiMap
#ifdef COPY2STORE  
  storeRepP gMultiMapRepI::copy2store() const{
    fatal(!sorted);
    gEntitySet dom = domain();
    entitySet pdom = gconvert(dom);
    entitySet v_dom = gconvert(vdom);
    store<int> count;
    count.allocate(pdom+v_dom);
    const_iterator itr = begin();
    FORALL(pdom, ei){
      int cnt = 0;
      while(itr!=end() && itr->first < ei)itr++;
      while(itr!=end() && itr->first == ei){
        cnt++;
        itr++;
      }
      count[ei] = cnt;
    }ENDGFORALL;
    
    FORALL(v_dom, ei){
      count[ei] = 0;
    }ENDGFORALL;
    
    multiMap s;
    s.allocate(count);
    itr = begin();
    GFORALL(dom, ei){
      int cnt = 0;
      while(itr!=end() && itr->first < ei)itr++;
      while(itr!=end() && itr->first == ei){
        s[ei][cnt++] = int(itr->second);
        itr++;
      }
    }ENDGFORALL;
    
    return s.Rep();
  }
#endif

  
}
