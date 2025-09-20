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

#ifndef GSTORE_IMPL_H
#define GSTORE_IMPL_H 1
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
#include <gmultistore.h>
#include <distribute_long.h>
#include <field_sort.h>
using std::vector;

namespace Loci {
  

  template<class T> void gStore<T>::make_consistent(){
    debugout<< "WARNING: void gStore<T>::make_consistent() not implemented yet" << endl;
  }
  
  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const gStore<T> &t)
  { return t.Print(s) ; }


  template<class T> void gStoreRepI<T>::local_sort(){
    std::sort(attrib_data.begin(), attrib_data.end(), fieldSort1<T>);
    sorted = true;
  }

  // //parSampleSort not working, because balanceDisttribution will split one entity inside a multiMap
  // // or a multiStore 
  // template<class T> void gStoreRepI<T>::sort(){
  //   cerr<<" WARNING: gStoreRepI<T>::sort() not implemented yet" << endl; 
  // }

  // template<class T> void gStoreRepI<T>::parSplitSort(const std::vector<gEntitySet>& ptn){
  //   cerr<<" WARNING:gStoreRepI<T>::parSplitSort not implemented yet" << endl; 
  // }
 
                         


  template<class T>
  gStoreRepP gStoreRepI<T>::remap(const gMap &m) const{
    fatal(!sorted);
    //    fatal( !m.sorted());
    fatal(m.domain()-domain() != GEMPTY);
    gStore<T> s;
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
  std::ostream &gStoreRepI<T>::Print(std::ostream &s) const 
  {
    s << '{' << domain() << std::endl ;
    for(const_iterator itr = begin(); itr != end(); itr++){
      Loci::streamoutput(&(itr->second),1,s) ;
    }
    s << '}' << std::endl ;
    return s ;
  }

  //************************************************************************/

  template<class T> 
  std::istream &gStoreRepI<T>::Input(std::istream &s) 
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
        
    GFORALL(e,ii) {
      T val;
      Loci::streaminput(&val,1,s) ;
      insert(ii, val);
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
  gstore_type gStoreRepI<T>::RepType() const 
  {
    return GSTORE ;
  }

  template<class T>
  gStoreRepP gStoreRepI<T>::recompose( gStoreRepP &newmap,  MPI_Comm comm)const {
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

  
  //For stores,recompose will compose a store whose domain is the domain of m,
  //whose data is the data of the SECOND field of m.
  //for example, pos.recompose(face2node) will produce the positions  for each face  
  template<class T>
  gStoreRepP  gStoreRepI<T>::recompose(const gMultiMap &remap, MPI_Comm comm)const{
    gEntitySet new_dom = remap.image();
    gEntitySet old_dom = domain();
    vector<gEntitySet> ptn = g_all_collect_vectors<gEntity>(new_dom);
    gStore<T> expanded_store;//the value in expanded store is unique in each process
    expanded_store = split_redistribute(ptn, comm);
    expanded_store.local_sort();

    //copy remap into a vector
    vector<ms_comp<gEntity> > temp_vec;
    {
      short ind = 0;
      gMultiMap::const_iterator previous = remap.begin();
      temp_vec.push_back(ms_comp<gEntity>(previous->first, previous->second, ind));
      gMultiMap::const_iterator itr = remap.begin();
      itr++;
      for(; itr != remap.end(); itr++){
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

    vector<ms_comp<T> > temp_vec2;
    typename gStore<T>::const_iterator itr = expanded_store.begin();
    for(size_t i = 0; i < temp_vec.size(); i++){
      while(itr->first < temp_vec[i].img )itr++ ;
      if(itr->first == temp_vec[i].img) temp_vec2.push_back( ms_comp<T>(temp_vec[i].dom, itr->second, temp_vec[i].ind));
    }
    
    //sort temp_vec2 according to the dom and ind field
    std::sort(temp_vec2.begin(), temp_vec2.end(), field_sort_dom<T>);                                       

    gMultiStore<T> result;
    for(size_t i = 0; i < temp_vec2.size(); i++){
      result.insert(temp_vec2[i].dom, temp_vec2[i].img);
    }
    result.set_domain_space(domain_space);
    return result.Rep();
  }

  template<class T>
  gStoreRepP  gStoreRepI<T>::recompose(const gMap &remap, MPI_Comm comm)const{
        
    gEntitySet new_dom = remap.image();
    gEntitySet old_dom = domain();
    vector<gEntitySet> ptn = g_all_collect_vectors<gEntity>(new_dom);
    gStore<T> expanded_store;//the value in expanded store is unique in each process
    expanded_store = split_redistribute(ptn, comm);
    expanded_store.local_sort();

    const_cast<gMap &>(remap).local_sort2();

    gStore<T> result;
    typename gStore<T>::const_iterator itr2 = expanded_store.begin();
    for(gMap::const_iterator itr1 = remap.begin(); itr1!= remap.end(); itr1++){
      while(itr2 != expanded_store.end() && itr2->first < itr1->second)itr2++ ;
      if(itr2 != end() && itr2->first == itr1->second) result.insert(itr1->first, itr2->second);
    }
    result.local_sort();
    result.set_domain_space(domain_space);
    return result.Rep();
  }
  
  template<class T>
  gStoreRepP gStoreRepI<T>::redistribute(const std::vector<gEntitySet>& dom_split,
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

    
    gStore<T> ns ;
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
  template<class T>
  gStoreRepP gStoreRepI<T>::split_redistribute(const std::vector<gEntitySet>& dom_ptn,
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
  
  template<class T> gStoreRepP
  gStoreRepI<T>::redistribute(const std::vector<gEntitySet>& dom_split,
                              const gMap& remap, MPI_Comm comm)const {
    fatal(!sorted);

    gStore<T> result;
    result.setRep(redistribute( dom_split,comm));
    result.local_sort();
    result.set_domain_space(domain_space);
    return result.remap(remap);
  }

  //   template<class T> storeRepP
  //   gStoreRepI<T>::redistribute_omd(const std::vector<gEntitySet>& dom_ptn,
  //                                   const dMap& remap, MPI_Comm comm) {
  //   //   // this is a push operation, thus the send, recv are reversed
  // //     std::vector<P2pCommInfo> send, recv ;
  // //     get_p2p_comm(dom_ptn, domain(), 0, 0, comm, recv, send) ;
  // //     gStore<T> new_store ;
  // //     fill_store_omd(getRep(), 0, new_store.Rep(), &remap, send, recv, comm) ;
  // //     return new_store.Rep() ;
  //   }

 

  //*************************************************************************/
  template <class T> 
  int gStoreRepI<T>::pack_size( const gEntitySet &eset)const {
   
    fatal(eset - domain() != GEMPTY);
   
    fatal(!sorted);
   
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
   
    return get_mpi_size( traits_type, eset );
  }
  
  

  //*******************************************************************/
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T>
  int gStoreRepI<T>::get_mpi_size( IDENTITY_CONVERTER c,
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
  int gStoreRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c,
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
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T> 
  void gStoreRepI<T>::pack( void *outbuf, int &position, int size, 
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
  /* this function applys to both  gStore and gMultiStore */ 
  template <class T> 
  void gStoreRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
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
  template <class T>
  void gStoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
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

  template <class T> 
  void gStoreRepI<T>::unpack(const void *ptr, int &loc, int size ) 
  {
    
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    unpackdata(traits_type, ptr, loc, size);
    local_sort();
  }

  //*********************************************************************/
  template <class T> 
  void gStoreRepI<T>::unpackdata(IDENTITY_CONVERTER c, const void *inbuf,
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
  template <class T> 
  void gStoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, const void *inbuf, 
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

  //*******************************************************************/
  template<class T>
  frame_info gStoreRepI<T>::get_frame_info()const {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }
  //*******************************************************************/ 
  template<class T>
  frame_info gStoreRepI<T>::get_frame_info(IDENTITY_CONVERTER g)const {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = 1 ;
    return fi ;
  }
  //*******************************************************************/
  template<class T>
  frame_info gStoreRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g)const {

    gEntitySet dom = domain() ;
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = 1 ;
    int stateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    for(const_iterator ci = begin(); ci != end(); ++ci){
      typename schema_traits::Converter_Type cvtr(const_cast<T&>(ci->second));
      stateSize = cvtr.getSize();
      fi.second_level.push_back(stateSize) ;
    }
    return fi ;
  }

  
 
  //*********************************************************************/
  template<class T> 
  void gStoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                               const char* name, frame_info &fi, const gEntitySet &eset)
  {
     
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset) ;
  }

  //*************************************************************************/

  template<class T>
  void  gStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                const char* name, IDENTITY_CONVERTER c, frame_info &fi, const gEntitySet &eset)
  {
    if(dimension != 0) {
      gStoreRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                          H5P_DEFAULT, tmp_array) ;
      if(err < 0) {
        std::string error = "H5Dread: failed" ;
        throw StringError(error) ;
      }
      for(gEntitySet::const_iterator si = eset.begin(); si != eset.end();++si)
        insert(*si, tmp_array[tmp++]) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
      local_sort();
    } 
  }
  //*************************************************************************/

  template<class T>
  void  gStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, const gEntitySet &eset)
  {
    typedef data_schema_traits<T> schema_traits ;
    if(dimension != 0) {
      gStoreRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      std::vector<int> vint = fi.second_level ;
      typedef typename schema_traits::Converter_Base_Type dtype;
      dtype* tmp_array = new dtype[dimension] ;
      hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
			  H5P_DEFAULT, tmp_array) ;
      if(err < 0) {
	std::string error = "H5Dread: failed" ;
	throw StringError(error) ;
      }
      size_t tmp = 0 ;
      int bucsize ;
      size_t indx = 0 ;
      for(gEntitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) {
        T temp;
        typename schema_traits::Converter_Type cvtr( temp );
        bucsize = vint[indx++] ;
	cvtr.setState(tmp_array+tmp, bucsize) ;
	tmp += bucsize ;
        insert(*si, temp);
      }
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
      local_sort();
    } 
  }
  //*************************************************************************/
  
  template<class T> 
  void gStoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                const char* name, const gEntitySet& usr_eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, usr_eset) ;
  }

  //*************************************************************************/
  template <class T> 
  void gStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                const char* name, IDENTITY_CONVERTER c, const gEntitySet &eset)  const
  {
    if(dimension != 0) {
      gStoreRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      const_iterator itr = begin();
      for(gEntitySet::const_iterator si = eset.begin(); si != eset.end();++si){
        while(itr!= end() && itr->first < *si) itr++;
        if(itr!=end() && itr->first==*si)tmp_array[tmp++] = itr->second ;
      }
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
     
  }

  //*************************************************************************/

  template <class T> 
  void gStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                const char* name, USER_DEFINED_CONVERTER c, const gEntitySet &eset)  const
  {   
    typedef data_schema_traits<T> schema_traits ;
    if(dimension != 0) {
      gStoreRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      typedef typename schema_traits::Converter_Base_Type dtype;
      dtype* tmp_array = new dtype[dimension] ;
      size_t tmp = 0 ;
      int stateSize = 0 ;
      const_iterator itr = begin();
      for(gEntitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
        while(itr!= end() && itr->first < *si) itr++;
        if(itr!=end() && itr->first==*si){
          typename schema_traits::Converter_Type cvtr(const_cast<T&>(itr->second));
          cvtr.getState(tmp_array+tmp, stateSize) ;
          tmp +=stateSize ;
        }
      }
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
     
  }
}

#endif
