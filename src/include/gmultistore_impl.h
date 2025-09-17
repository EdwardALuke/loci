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

#ifndef GMULTISTORE_IMPL_H
#define GMULTISTORE_IMPL_H 1
#include <data_traits.h>
#include <Tools/except.h>
#include <parSampleSort.h>
#include <vector>
#include <sstream>
#include <iterator>
#include <gstore_rep.h>
#include <hdf5_readwrite_long.h>
#include <mpi.h>
#include <string.h>
#include <dist_internal.h>
#include <gmap.h>
#include <gmultimap.h>
#include <distribute_long.h>
#include <field_sort.h>

using std::vector;

namespace Loci {
  
  template<class T> void gMultiStore<T>::make_consistent(){
    debugout<< "WARNING: void gMultiStore<T>::make_consistent() not implemented yet" << endl;
  }
  
  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const gMultiStore<T> &t)
  { return t.Print(s) ; }


  template<class T> void gMultiStoreRepI<T>::local_sort(){
    std::stable_sort(attrib_data.begin(), attrib_data.end(), fieldSort1<T>);
    sorted = true;
  }

  // //parSampleSort not working, because balanceDisttribution will split one entity inside a multiMap
  // // or a multiStore 
  // template<class T> void gMultiStoreRepI<T>::sort(){
  //   cerr<<" WARNING: gMultiStoreRepI<T>::sort() not implemented yet" << endl; 
  // }

  // template<class T> void gMultiStoreRepI<T>::parSplitSort(const std::vector<gEntitySet>& ptn){
  //   cerr<<" WARNING:gMultiStoreRepI<T>::parSplitSort not implemented yet" << endl; 
  // }
 
  template<class T>
  gStoreRepP gMultiStoreRepI<T>::remap(const gMap &m) const{
    //    fatal( !m.sorted());
    fatal(m.domain()-domain() != GEMPTY);
    fatal(!sorted);
    gMultiStore<T> s;
    typename gRep::const_iterator itr1 = begin();
    gMap::const_iterator itr2 = m.begin();
    for(; itr1 != end(); itr1++){
      while(itr2 != m.end() && itr2->first < itr1->first)itr2++ ;
      if(itr2 != m.end() && itr2->first == itr1->first) s.insert(itr2->second, itr1->second);
    }
    s.local_sort();
    s.set_domain_space(domain_space);
    return s.Rep();
  }

  //*********************************************************************/

  template<class T> 
  std::ostream &gMultiStoreRepI<T>::Print(std::ostream &s) const 
  {
    s << '{' << domain() << endl ;
    s << vdom <<endl;
    for( const_iterator itr = attrib_data.begin(); itr != attrib_data.end(); itr++){ 
      s << itr->first <<' ';
      Loci::streamoutput(&(itr->second),1,s) ;
    }
    s << '}' << endl ;
    return s ;
  }

  //************************************************************************/

  template<class T> 
  std::istream &gMultiStoreRepI<T>::Input(std::istream &s) 
  {
   
    gEntitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    s >> e ;
    s >> vdom;
    
    gEntity val1 = 0;
    T val2;
    GFORALL(e,ii) {
      s >> val1; 
      Loci::streaminput(&val2,1,s) ;
      insert(val1, val2);
    } ENDGFORALL ;
        
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    local_sort();
    return s ;
  }

  
  

  //*************************************************************************/
  
  template<class T> 
  gstore_type gMultiStoreRepI<T>::RepType() const 
  {
    return GMULTISTORE ;
  }

  

  //This method does not consider vdom
  //For stores,recompose will compose a store whose domain is the domain of m,
  //whose data is the data of the SECOND field of m.
  //for example, pos.recompose(face2node) will produce the positions  for each face  
  template<class T>
  gStoreRepP  gMultiStoreRepI<T>::recompose(gStoreRepP& m, MPI_Comm comm)const{
    cerr<<"WARNING: gMultiStoreRepI::recompose() not implemented yet" << endl;
    
    // gMultiMap remap(m);
    
    // gEntitySet new_dom = remap.image();
    // gEntitySet old_dom = domain();
    // vector<gEntitySet> ptn = g_all_collect_vectors<gEntity>(new_dom);
    // gMultiStore<T> expanded_store;//the value in expanded store is unique in each process
    // expanded_store = split_redistribute(ptn, comm);
    // expanded_store.local_sort();
    
    // //put expanded store into std::map
    // std::map<gEntity, T> amap;
    // for(const_iterator itr = expanded_store.begin(); itr!=expanded_store.end(); itr++)
    //   amap[itr->first] = itr->second;

    gMultiStore<T> result;
    // for(gMultiMap::const_iterator mi = remap.begin(); mi!= remap.end(); mi++){
    //   result.insert(mi->first, amap[mi->second]);
    // }
    result.set_domain_space(domain_space);
    return result.Rep();
  }


  //This method does not consider vdom
  template<class T>
  gStoreRepP gMultiStoreRepI<T>::redistribute(const std::vector<gEntitySet>& dom_split,
                                              MPI_Comm comm)const {
    fatal(!sorted);
    int np = 1;
    MPI_Comm_size(comm, &np);
    // compute the send_counts
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
    for(int i = 0; i <  MPI_processes; ++i)
      pack(send_buffer, loc_pack, tot_send_size, dom_split[i]) ;
    
    // we now proceed to communicate the store contents
    unsigned char* recv_buffer = new unsigned char[tot_recv_size] ;
    MPI_Alltoallv(&send_buffer[0], &send_counts[0],
                  &send_displs[0], MPI_PACKED,
                  &recv_buffer[0], &recv_counts[0],
                  &recv_displs[0], MPI_PACKED, comm) ;

    
    gMultiStore<T> ns ;
    gStoreRepP srp = ns.Rep() ;
    // unpack
    loc_pack = 0 ;
    srp->unpack(recv_buffer, loc_pack, tot_recv_size) ;
    srp->local_sort();
    srp->set_domain_space(domain_space);
    
    delete[] recv_displs ;
    delete[] send_displs ;
    delete[] send_buffer ;
    delete[] recv_buffer ;
    return srp ;
  }
  
  
  //split domain first, then redistribute
  //This method does not consider vdom
  template<class T>
  gStoreRepP gMultiStoreRepI<T>::split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                                                    MPI_Comm comm)const {
    fatal(!sorted);

    gEntitySet my_dom = domain() ;
    // figure out how the domain is split to send to others
    int np  = 1;
    MPI_Comm_size(comm, &np) ;
    fatal(np != dom_ptn.size()) ;
    std::vector<gEntitySet> dom_split(np) ;
    for(int i=0;i<np;++i) {
      dom_split[i] = my_dom & dom_ptn[i] ;
    }
    return redistribute(dom_split, comm);
  }

  //This method does not consider vdom
  template<class T> gStoreRepP
  gMultiStoreRepI<T>::redistribute(const std::vector<gEntitySet>& dom_split,
                                   const gMap& remap, MPI_Comm comm)const {
    fatal(!sorted);

    gMultiStore<T> result;
    result.setRep(redistribute( dom_split,comm));
    result.local_sort();
    result.set_domain_space(domain_space);
    return result.remap(remap);
  }

  //   template<class T> storeRepP
  //   gMultiStoreRepI<T>::redistribute_omd(const std::vector<gEntitySet>& dom_ptn,
  //                                   const dMap& remap, MPI_Comm comm) {
  //   //   // this is a push operation, thus the send, recv are reversed
  // //     std::vector<P2pCommInfo> send, recv ;
  // //     get_p2p_comm(dom_ptn, domain(), 0, 0, comm, recv, send) ;
  // //     gMultiStore<T> new_store ;
  // //     fill_store_omd(getRep(), 0, new_store.Rep(), &remap, send, recv, comm) ;
  // //     return new_store.Rep() ;
  //   }

 

  //*************************************************************************/
  template <class T> 
  int gMultiStoreRepI<T>::pack_size( const gEntitySet &eset)const {
   
    fatal(eset - domain() != GEMPTY);
   
    fatal(!sorted);
   
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
   
    return get_mpi_size( traits_type, eset );
  }
  
  

  //*******************************************************************/
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T>
  int gMultiStoreRepI<T>::get_mpi_size( IDENTITY_CONVERTER c,
                                        const gEntitySet &eset)const{
    int result = 0;
    const_iterator itr1 = attrib_data.begin();
    GFORALL(eset, ci){
      while(itr1 != attrib_data.end() && itr1->first < ci  )itr1++;
      while(itr1 != attrib_data.end()&& itr1->first == ci ){
        result += (sizeof(gEntity)+sizeof(T));
        itr1++;
      }
    }ENDGFORALL;
    return result;
  }
  

  
  
  //*******************************************************************/
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T>
  int gMultiStoreRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c,
                                        const gEntitySet &eset)const
  {
   
    int result = 0;
    int       size, numBytes=0;
    gEntitySet  ecommon;
    typedef data_schema_traits<T> converter_traits;
    const_iterator itr1 = attrib_data.begin();
    GFORALL(eset, ci){
      while(itr1 != attrib_data.end() && itr1->first < ci  )itr1++;
      while(itr1 != attrib_data.end() && itr1->first == ci){
        typename converter_traits::Converter_Type cvtr(const_cast<T&>(itr1->second));
        size      = cvtr.getSize();
        numBytes = size*sizeof(typename converter_traits::Converter_Base_Type) ;
        result += (sizeof(gEntity)+sizeof(int)+ numBytes);
        itr1++; 
      }
    }ENDGFORALL;
    return result;
  }
  
  //*******************************************************************/
  //This method does not consider vdom
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T> 
  void gMultiStoreRepI<T>::pack( void *outbuf, int &position, int size, 
                                 const gEntitySet &usr_eset ) const 
  {
    fatal(usr_eset - domain() != GEMPTY);
    fatal(!sorted);
    gEntitySet ecommon;
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
    ecommon = domain()&usr_eset;
    packdata( traits_type, outbuf, position, size, ecommon);
  }
 
  //*******************************************************************/
  //This method does not consider vdom 
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T> 
  void gMultiStoreRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
                                    int &position, int outcount,
                                    const gEntitySet &eset )  const
  {
    const_iterator itr1 = attrib_data.begin();
    GFORALL(eset, ci){
      
      while(itr1 != attrib_data.end() && itr1->first < ci) itr1++;
      while(itr1 != attrib_data.end() && itr1->first == ci){
        MPI_Pack( const_cast<gEntity*>(&ci), sizeof(gEntity), MPI_BYTE, outbuf,outcount, 
                  &position, MPI_COMM_WORLD) ;
        MPI_Pack( const_cast<T*>(&(itr1->second)), sizeof(T), MPI_BYTE, outbuf,outcount, 
                  &position, MPI_COMM_WORLD) ;
        itr1++;
      }
    }ENDGFORALL;
  }

  //*******************************************************************/
  /* this function applys to both  gStore and gMultiStore */
  //This method does not consider vdom
  template <class T>
  void gMultiStoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                                     int &position, int outcount, 
                                     const gEntitySet &eset )const
  {
    warn(sorted);
    int  stateSize, incount;
    const_iterator itr1 = attrib_data.begin();
    typedef  data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;
    std::vector<dtype>  inbuf;
    int typesize = sizeof(dtype);
    //-------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-------------------------------------------------------------------
    for(gEntitySet :: const_iterator ci = eset.begin(); ci != eset.end(); ci++) {
      while(itr1 != attrib_data.end() && itr1->first < *ci) itr1++;
      while(itr1 != attrib_data.end()&& itr1->first == *ci ){
        typename schema_traits::Converter_Type cvtr( const_cast<T&>(itr1->second));
        stateSize    = cvtr.getSize();
        if( stateSize > static_cast<int>(inbuf.size()) ) inbuf.resize(stateSize);
        cvtr.getState( &inbuf[0], stateSize);
        
        MPI_Pack( const_cast<gEntity*>(&(itr1->first)), sizeof(gEntity), MPI_BYTE, outbuf,outcount, 
                  &position, MPI_COMM_WORLD) ;
        MPI_Pack(&stateSize, 1, MPI_INT, outbuf, outcount,&position,
                 MPI_COMM_WORLD);
        
        incount =  stateSize*typesize;
        MPI_Pack(&inbuf[0], incount, MPI_BYTE, outbuf, outcount, &position, 
                 MPI_COMM_WORLD) ;
        itr1++;
      }
    }
  }

  //*******************************************************************/
  //This method does not consider vdom
  template <class T> 
  void gMultiStoreRepI<T>::unpack(const void *ptr, int &loc, int size ) 
  {
    
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    unpackdata(traits_type, ptr, loc, size);
    local_sort();
  }

  //*********************************************************************/
  //This method does not consider vdom 
  template <class T> 
  void gMultiStoreRepI<T>::unpackdata(IDENTITY_CONVERTER c, const void *inbuf,
                                      int &position,  int insize
                                      ) 
  {
    gEntity ei;
    T val;
    while(position != insize){
      MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &ei,
                  sizeof(gEntity), MPI_BYTE, MPI_COMM_WORLD) ;
      MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &val,
                  sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
      insert(ei, val);
    }
  }
  //*********************************************************************/
  //This method does not consider vdom
  template <class T> 
  void gMultiStoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, const void *inbuf, 
                                       int &position, int insize) 
  {
    int  stateSize, outcount;
   

    typedef  data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;
    std::vector<dtype>  outbuf;
    gEntity ei;
    int typesize = sizeof(dtype);
    while(position != insize){
      MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &ei,
                  sizeof(gEntity), MPI_BYTE, MPI_COMM_WORLD) ;
      MPI_Unpack(const_cast<void*>(inbuf), insize, &position, &stateSize, 1, 
                 MPI_INT, MPI_COMM_WORLD) ;
      if( stateSize > static_cast<int>(outbuf.size()) ) outbuf.resize(stateSize);
      outcount = stateSize*typesize;
      MPI_Unpack(const_cast<void*>(inbuf),insize, &position, &outbuf[0], outcount, 
                 MPI_BYTE, MPI_COMM_WORLD) ;
      T temp;
      typename schema_traits::Converter_Type cvtr( temp );
      cvtr.setState( &outbuf[0], stateSize);
      insert(ei, temp);
    }

  }

 
  //*********************************************************************/ 
  //this method computes the mean value for each entity in dom
  //and put them in a gStore
  template <class T> 
  gStoreRepP gMultiStore<T>::get_simple_center(const gEntitySet& dom)const{
    gEntitySet mydom = domain();
    gEntitySet iset = dom & mydom;
    gStore<T> result;
    const_iterator itr = begin();
    for(gEntitySet::const_iterator ei = iset.begin(); ei != iset.end(); ei++){
      while(itr !=end() && itr->first < *ei) itr++;
      T pnt = T();
      size_t sz = 0;
      while(itr !=end() && itr->first == *ei){
        if(sz==0) pnt = itr->second;
        else pnt += itr->second;
        sz++;
        itr++;
      }
      if(sz != 0){
        pnt *= 1./double(sz);
        result.insert(*ei, pnt);
      }
    }
    return result.Rep();
  }

  //*********************************************************************/ 
  //this method computes the weighted average value for each entity
  //and put them in a gStore, the weight is the length of neighboring nodes
  template <class T> 
  gStoreRepP gMultiStore<T>::get_wireframe_center(const gEntitySet& dom)const{
    gEntitySet iset = dom & domain();
    gStore<T> result;
    const_iterator itr = begin(), first_itr = begin();
    for(gEntitySet::const_iterator ei = iset.begin(); ei != iset.end(); ei++){
      while(itr !=end() && itr->first < *ei) itr++;
      T nodesum = T(0.0, 0.0, 0.0);
      double lensum = 0;
      
      while(itr !=end() && itr->first == *ei){
        if(lensum == 0) first_itr = itr;
        const_iterator next = itr;
        next++;
        if(next != end() && next->first == *ei){
          T edge_loc = 0.5*(itr->second + next->second) ;
          T edge_vec = itr->second - next->second;
          double len = norm(edge_vec);
          nodesum += len*edge_loc ;
          lensum += len; 
        }else{
          next = first_itr;
          T edge_loc = 0.5*(itr->second + next->second) ;
          T edge_vec = itr->second - next->second;
          double len = norm(edge_vec);
          nodesum += len*edge_loc ;
          lensum += len;
          result.insert(*ei, nodesum/lensum);
          nodesum = T(0.0, 0.0, 0.0);
          lensum = 0; 
        }
        itr++;
      }
    }
    return result.Rep();
  }
  
  //*********************************************************************/ 
  //this method computes the area for each entity in dom
  //and put them in a gStore
 
  template <class T> 
  gStoreRepP gMultiStore<T>::get_area(const gStore<T>& center)const{
    gEntitySet dom = center.domain();
    gEntitySet iset = dom & domain();
    gStore<double> result;
    const_iterator itr = begin();
    const_iterator first_itr = begin();
    typename gStore<T>::const_iterator center_itr = center.begin();
    for(gEntitySet::const_iterator ei = iset.begin(); ei != iset.end(); ei++, center_itr++){
      while(itr !=end() && itr->first < *ei) itr++;
      T sum = T(0.0, 0.0, 0.0);
      int sz = 0;
      while(itr !=end() && itr->first == *ei){
        if(sz == 0) first_itr = itr;
        const_iterator next = itr;
        next++;
        if(next != end() && next->first == *ei){
          sum += cross((itr->second)- (center_itr->second), (next->second) - (center_itr->second)) ;
          sz++;
        }else{
          next = first_itr;
          sum += cross((itr->second)- (center_itr->second), (next->second) - (center_itr->second)) ;
          sz++;
          result.insert(*ei, 0.5*norm(sum));
          sum = T(0.0, 0.0, 0.0);
        }
        itr++;
      }
    }
    return result.Rep();
  }

  //*********************************************************************/ 
  //this method computes the area for each entity in dom
  //and put them in a gStore
  template <class T> 
  gStoreRepP gMultiStore<T>::get_area(const gEntitySet& dom)const{
    gStore<T> center;
    center = get_wireframe_center(dom);
    return get_area(center);
  }  

  //*********************************************************************/ 
  //this method computes the  weighted mean value for each entity in gMultiStore
  //and put them in a gStore, the weight is given in area
  //prerequest: area and this have the same domain, and for each entity in domain
  // the number of elements in area and the number of elements in this are the same 
  //example: this is the recomposed facecenter, and area is the recomposed facearea
  // return value is the centroid of each cell
  template <class T> 
  gStoreRepP gMultiStore<T>::get_weighted_center(const gMultiStore<double>& area)const{
     fatal(area.size() != size());
     fatal(area.domain()- domain() != GEMPTY);
     gStore<T> result;
     const_iterator itr = begin();
     gStore<double>::const_iterator area_itr = area.begin();
     gEntity current = itr->first;
     T nodesum = (itr->second)*(area_itr->second) ;
     
     double  areasum = area_itr->second;
     itr++;
     area_itr++;
     for( ; itr != end(); itr++, area_itr++){
      if(itr->first == current){
        nodesum += (itr->second) * (area_itr->second);
        areasum += (area_itr->second);
      }else{
        nodesum *= 1./areasum ;
        result.insert(current, nodesum) ;
        current = itr->first;
        nodesum = (itr->second)*(area_itr->second);
        areasum = area_itr->second;
      }
     }
     nodesum *= 1./areasum ;
     result.insert(current, nodesum);
     return result.Rep();
  }
    
  

  //*********************************************************************/
  template<class T> 
  frame_info gMultiStoreRepI<T>::get_frame_info()const {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }

  //*********************************************************************/
  template<class T> 
  frame_info gMultiStoreRepI<T>::get_frame_info(IDENTITY_CONVERTER g)const {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  
  //*********************************************************************/
  template<class T> 
  frame_info gMultiStoreRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g)const {
    warn(true) ;
    frame_info fi ;
    return fi ;
  } 
  
  //*********************************************************************/
  template<class T> 
  void gMultiStoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                    const char* name, frame_info &fi,const gEntitySet &user_eset)
  {
    warn(true) ;
     
  }

  //*************************************************************************/

  template<class T>
  void  gMultiStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                     const char* name, IDENTITY_CONVERTER c, frame_info &fi, const gEntitySet &eset)
  {
    warn(true) ;
     
  }
  //*************************************************************************/

  template<class T>
  void  gMultiStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                     const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, const gEntitySet &eset)
  {
    warn(true) ;
      
  }
  //*************************************************************************/

  template<class T> 
  void gMultiStoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                     const char* name, const gEntitySet& usr_eset) const
  {
    warn(true) ;
     
  }

  //*************************************************************************/
  template <class T> 
  void gMultiStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                     const char* name, IDENTITY_CONVERTER c, const gEntitySet &eset)  const
  {
    warn(true) ;
      
  }

  //*************************************************************************/

  template <class T> 
  void gMultiStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                     const char* name, USER_DEFINED_CONVERTER c, const gEntitySet &eset)  const
  {   
    warn(true) ;
     
  }
  //*************************************************************************/
}
#endif
