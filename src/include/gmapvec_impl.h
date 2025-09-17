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

#ifndef GMAPVEC_IMPL_H_
#define GMAPVEC_IMPL_H_

#include <istream>
#include <ostream>
#include <iostream>
#include <fstream>
#include <MapVec.h>
#include <gmapvec_def.h>
#include <gmap.h>
#include <gmultimap.h>
#include <hdf5_readwrite_long.h>
#include <distribute_long.h>
#include <distributed_inverse_long.h>
#include <field_sort.h>
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
  
  template<unsigned int M> 
  void gMapVecRepI<M>::local_sort(){
    std::stable_sort(attrib_data.begin(), attrib_data.end(), fieldSort1<gEntity>);
    sorted = true;
  }
  
  
  
  //**************************************************************************/
  template<unsigned int M> 
  gStoreRepP gMapVecRepI<M>::expand(gEntitySet &dom, std::vector<gEntitySet> &ptn, MPI_Comm comm)const {
    int np;
    MPI_Comm_size(comm, &np);
    std::vector<gEntitySet> recv_request(np);
    for(int i = 0; i < np; ++i) {
      recv_request[i] = dom & ptn[i] ; //recv request
    }
    vector<gEntitySet> send_ptn = transposePtn(recv_request,comm);
    return redistribute(send_ptn, comm);
  }
  
  //**************************************************************************/
  //the result of local_inverse is a gMultiMap
  template<unsigned int M> 
  gStoreRepP gMapVecRepI<M>::local_inverse()const {
    gMultiMap result;
    for(const_iterator itr = begin(); itr != end(); itr++){
      result.insert(itr->second, itr->first);
    }
    result.local_sort();
    return result.Rep();
  }
  
  //**************************************************************************/
  template<unsigned int M> 
  gStoreRepP gMapVecRepI<M>::remap(const gMap &m) const 
  {
    fatal(!(sorted && m.sorted()));
    fatal(m.domain()-domain() != GEMPTY);
    gMapVec<M> s ;
    
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

  //**************************************************************************/
  template<unsigned int M>
  void gMapVecRepI<M>::inplace_remap(const gMap &m) 
  {
    fatal(!(sorted && m.sorted()));
    gMap::const_iterator itr2 = m.begin();
    for( iterator itr1 = begin(); itr1 != end(); itr1++){
      while(itr2 != m.end() && itr2->first < itr1->first)itr2++ ;
      if(itr2 != m.end() && itr2->first == itr1->first) itr1->first = itr2->second;
    }
    local_sort();
  }

 
  //**************************************************************************/
  template<unsigned int M>
  gStoreRepP gMapVecRepI<M>::redistribute(const std::vector<gEntitySet>& dom_split,const gMap& remap,
                                          MPI_Comm comm)const {
    gMapVec<M> result;
    result = redistribute(dom_split, comm);
    result.set_domain_space(domain_space);
    result.set_image_space(image_space);
    return result.remap(remap);
  }
  
  //**************************************************************************/
  template<unsigned int M>     
  gStoreRepP gMapVecRepI<M>::redistribute(const std::vector<gEntitySet>& dom_split,
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

  

    gMapVec<M> ns ;
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

  
 
  //**************************************************************************/
  //split domain first, then redistribute
  template<unsigned int M>
  gStoreRepP gMapVecRepI<M>::split_redistribute(const std::vector<gEntitySet>& dom_ptn,
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
  template<unsigned int M>
  int gMapVecRepI<M>::pack_size(const gEntitySet &eset)const {
    int result = 0;
    warn(sorted);
   
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
  template<unsigned int M> 
  void gMapVecRepI<M>::pack(void *outbuf, int &position, int outsize, const gEntitySet &eset) const
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
  template<unsigned int M>
  void gMapVecRepI<M>::unpack(const void *inbuf, int &position, int insize) 
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
  
  
  
  //**************************************************************************/
  template<unsigned int M>
  gEntitySet gMapVecRepI<M>::image( gEntity iset) const 
  {
    gEntitySet codomain ;
    std::pair<gEntity, gEntity> p(iset, gEntity(0));
    std::pair<const_iterator, const_iterator> range = std::equal_range(begin(), end(), p, fieldSort1<gEntity>);
    if(range.first != range.second){
      for(const_iterator itr= range.first; itr != range.second; itr++) codomain += itr->second;
    }
    return codomain ;
  }

  //**************************************************************************/
  template<unsigned int M>
  gEntitySet gMapVecRepI<M>::image(const gEntitySet &iset) const 
  {
    
    gEntitySet codomain ;
    GFORALL(iset, ei){
      codomain += image(ei);
    }ENDGFORALL;
    return codomain ;
  }

  //**************************************************************************/
  template<unsigned int M>
  gEntitySet gMapVecRepI<M>::image() const 
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
  template<unsigned int M>
  gStoreRepP gMapVecRepI<M>::get_map() const 
  {
    gMapVec<M> result ;
    result.set_domain_space(domain_space);
    result.set_image_space(image_space);
    for( const_iterator itr = attrib_data.begin(); itr != attrib_data.end(); itr++){
      result.insert(itr->first, itr->second);
    }
    return result.Rep() ;
  }

  //**************************************************************************/
  template<unsigned int M>
  ostream &gMapVecRepI<M>::Print(ostream &s) const 
  {
    s << '{' << domain() << std::endl ;
        
    for( const_iterator itr = attrib_data.begin(); itr != attrib_data.end(); itr++){
      s << itr->second<<endl;
    }
    s << '}' << endl ;
    return s ;
  }

  //**************************************************************************/
  template<unsigned int M>
  istream &gMapVecRepI<M>::Input(istream &s) 
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
      gEntity val;
      for(unsigned int j=0;j<M;++j){
        s >> val;
        insert(ii, val);
      }
    } ENDGFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading gStore" << endl ;
      s.putback(ch) ;
    }
    local_sort();
    return s ;
  }

  //**************************************************************************/
  template<unsigned int M>
  DatatypeP gMapVecRepI<M>::getType()const {
    return DatatypeP(new AtomicType(INT)) ;
  }
  //**************************************************************************/
  template<unsigned int M> 
  frame_info gMapVecRepI<M>::get_frame_info()const {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  //**************************************************************************/
  template<unsigned int M>
  void gMapVecRepI<M>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                const char* name, frame_info &fi, const gEntitySet &usr_eset)  {
    warn(true) ;
     
  } 

  //**************************************************************************/
  template<unsigned int M>
  void gMapVecRepI<M>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                 const char* name, const gEntitySet& usr_eset) const
  {
    warn(true) ;
  } 


  //this is the working version with std::map
  // //**************************************************************************/
  
  // void gMapVecRepI<M>::inplace_compose(const gMap &newmap,  MPI_Comm comm) 
  // {
   
  //   //check if expand is needed
  //   gEntitySet old_dom = newmap.domain();
  //   gEntitySet dom = image();
    
  //   gEntitySet out_dom = dom -  old_dom;
  //   int need_expand = (out_dom == GEMPTY ? 0:1);
  //   bool global_need_expand = GLOBAL_OR(need_expand, comm);
  //   std::map<gEntity, gEntity> amap;
  //   if(global_need_expand){
  //     std::vector<gEntitySet> old_ptn = g_all_collect_vectors<gEntity>(old_dom, comm) ;
  //     gMap temp;
  //     temp.setRep(newmap.expand(dom, old_ptn, comm));
  //     for(gMap::const_iterator itr = temp.begin(); itr!=temp.end(); itr++)
  //       amap[itr->first] = itr->second;
  //   }else{
  //     for(gMap::const_iterator itr = newmap.begin(); itr!=newmap.end(); itr++)
  //       amap[itr->first] = itr->second;
  //   }
 
  //   for( iterator itr1 = begin(); itr1 != end(); itr1++){
  //     itr1->second = amap[itr1->second];
  //   }
  //  Print(debugout);
  // }

  // //**************************************************************************/
  
  // gStoreRepP gMapVecRepI<M>::recompose(const gMap &newmap,  MPI_Comm comm) 
  // {
   
  //   //check if expand is needed
  //   gEntitySet old_dom = newmap.domain();
  //   gEntitySet dom = image();
    
  //   gEntitySet out_dom = dom -  old_dom;
  //   int need_expand = (out_dom == GEMPTY ? 0:1);
  //   bool global_need_expand = GLOBAL_OR(need_expand, comm);
  //   std::map<gEntity, gEntity> amap;
  //   if(global_need_expand){
  //     std::vector<gEntitySet> old_ptn = g_all_collect_vectors<gEntity>(old_dom, comm) ;
  //     gMap temp;
  //     temp.setRep(newmap.expand(dom, old_ptn, comm));
  //     for(gMap::const_iterator itr = temp.begin(); itr!=temp.end(); itr++)
  //       amap[itr->first] = itr->second;
  //   }else{
  //     for(gMap::const_iterator itr = newmap.begin(); itr!=newmap.end(); itr++)
  //       amap[itr->first] = itr->second;
  //   }

  //   gMultiMap result;
  //   for( iterator itr1 = begin(); itr1 != end(); itr1++){
  //     result.insert(itr1->first, amap[itr1->second]);
  //   }
  //   result.local_sort();
  //   result.set_domain_space(domain_space);
  //   result.set_image_space(image_space);
  //   return result.Rep();
  // }








  
  //**************************************************************************/
  template<unsigned int M>
  void gMapVecRepI<M>::inplace_compose(const gMap &newmap,  MPI_Comm comm) 
  { //check if expand is needed
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
  template<unsigned int M>
  gStoreRepP gMapVecRepI<M>::recompose(const gMap &newmap,  MPI_Comm comm)const 
  {
   
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
    gMapVec<M> result;
    for(unsigned int i = 0; i < temp_vec.size(); i++){
      result.insert(temp_vec[i].dom, temp_vec[i].img);
    }
    
    //result.local_sort();
    result.set_domain_space(domain_space);
    result.set_image_space(image_space);
    return result.Rep();
  }


  //**************************************************************************/
  template<unsigned int M> 
  gStoreRepP gMapVecRepI<M>::recompose(const gMultiMap &newmap,  MPI_Comm comm)const 
  {
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

  //**************************************************************************/
  template<unsigned int M> 
  gStoreRepP gMapVecRepI<M>::recompose( gStoreRepP &newmap,  MPI_Comm comm)const {
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
  template<unsigned int M> 
  std::pair<gEntitySet,gEntitySet> gMapVecRepI<M>::preimage(const gEntitySet &codomain) const {
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
  
 
  //**************************************************************************/
  template<unsigned int M>
  gStoreRepP gMapVecRepI<M>::distributed_inverse(gEntitySet global_image,
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
  
  // //**************************************************************************/
  // template<unsigned int M>
  // void gMapVecRepI<M>::remove_duplication(){
  //   std::sort(attrib_data.begin(), attrib_data.end(), fieldSort1_unique<gEntity>);
  //   vector<pair<Entity,Entity> >::iterator uend ;
  //   uend = unique(attrib_data.begin(), attrib_data.end());
  //   attrib_data.erase(uend, attrib_data.end());
  // }
  // //**************************************************************************/
  // template<unsigned int M>
  // void gMapVec<M>::remove_duplication(){
  //   CPTR<MapType> p(Rep()) ;
  //   if(p != 0)p->remove_duplication();
  //   warn(p==0);
  // }



  
  
  //**************************************************************************/
  //copy gMultiMap to traditional Loci multiMap
#ifdef COPY2STORE  
  template<unsigned int M> 
  storeRepP gMapVecRepI<M>::copy2store() const{
    entitySet dom = gconvert(domain());
    MapVec<M> s ;
    s.allocate(dom) ;
    int ind = 0;
    for(const_iterator itr = begin(); itr!= end(); itr++){
      int cnt = ind%M;
      Entity ei = itr->first;
      s[ei][cnt] = int(itr->second);
      ind++;
    }
    return s.Rep();
  }
#endif

  
}
#endif
