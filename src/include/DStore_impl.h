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
#ifndef DSTORE_IMPL_H
#define DSTORE_IMPL_H

#include <DStore_def.h>
#include <store_def.h>
#include <distribute.h>

namespace Loci {
  template<class T> 
  void dstoreRepI<T>::allocate(const entitySet &eset)
  {
    entitySet :: const_iterator ci;

    entitySet dom = domain() ;
    entitySet remove = dom - eset ;
    entitySet add = eset - dom ;

    attrib_data.erase_set(remove) ;

    for( ci = add.begin(); ci != add.end(); ++ci)
      attrib_data.access(*ci) ;
  
    dispatch_notify() ;
  }

  template<class T>
  void dstoreRepI<T>::shift(int_type offset) {
    entitySet new_domain = domain() ;
    new_domain >>= offset ;
    allocate(new_domain) ;
  }

  template<class T>
  void dstoreRepI<T>::erase(const entitySet& rm) {
    entitySet valid = domain() & rm ;
    attrib_data.erase_set(valid) ;
    dispatch_notify() ;
  }

  template<class T>
  void dstoreRepI<T>::invalidate(const entitySet& valid) {
    entitySet redundant = domain() - valid ;
    erase(redundant) ;
  }
  template<class T>
  void dstoreRepI<T>::guarantee_domain(const entitySet& include) {
    entitySet new_set = include - domain() ;
    for(entitySet::const_iterator ei=new_set.begin();
        ei!=new_set.end();++ei)
      attrib_data.access(*ei) ;
    dispatch_notify() ;
  }

  //*********************************************************************/

  template<class T> 
  std::ostream &dstoreRepI<T>::Print(std::ostream &s) const 
  {
    entitySet :: const_iterator it;
    entitySet dom = domain() ;

    s << '{' << dom << std::endl ;

    for(it = dom.begin();it!=dom.end();++it) {
      Loci::streamoutput(&(attrib_data.elem(*it)),1,s) ;
    }

    s << '}' << std::endl ;

    return s ;
  }

  //************************************************************************/

  template<class T> 
  std::istream &dstoreRepI<T>::Input(std::istream &s) 
  {
    entitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    s >> e ;
        
    FORALL(e,ii) {
      Loci::streaminput(&attrib_data.access(ii),1,s) ;
    } ENDFORALL ;
        
    do ch = s.get(); while(ch==' ' || ch=='\n') ;

    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }

    return s ;
  }

  //*************************************************************************/
  template<class T>  
  entitySet dstoreRepI<T>::domain() const 
  {
    return attrib_data.domain() ;
  }

  //*************************************************************************/

  template<class T>
  storeRep *dstoreRepI<T>::new_store(const entitySet &p) const 
  {
    return new dstoreRepI<T>(p) ;
  }
  template<class T>
  storeRep *dstoreRepI<T>::new_store(const entitySet &p, const int* cnt) const 
  {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a dstore " << endl ;
    return sp ;
  }
  //*************************************************************************/
  
  template<class T> 
  store_type dstoreRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //************************************************************************/
  template<class T> 
  void dstore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //*************************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const dstore<T> &t)
  { return t.Print(s) ; }

  //*************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, dstore<T> &t)
  { return t.Input(s) ; }

  //************************************************************************/
  template<class T> 
  void const_dstore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //*************************************************************************/

  template<class T> 
  store_instance::instance_type
  const_dstore<T>::access() const 
  { return READ_ONLY; }
        
  //*************************************************************************/

  template<class T> 
  storeRepP dstoreRepI<T>::remap(const dMap &m) const 
  {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    dstore<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    
    return s.Rep() ;
  }

  template<class T> storeRepP
  dstoreRepI<T>::redistribute(const std::vector<entitySet>& dom_ptn,
                              MPI_Comm comm) {
    entitySet dom = domain() ;
    // figure out how the domain is split to send to others
    int np ;
    MPI_Comm_size(comm, &np) ;
    fatal(size_t(np) != dom_ptn.size()) ;
    entitySet total_dom ;
    std::vector<entitySet> dom_split(np) ;
    for(int i=0;i<np;++i) {
      dom_split[i] = dom & dom_ptn[i] ;
      total_dom += dom_ptn[i] ;
    }

    // compute the send_counts
    std::vector<int> send_counts(np,0) ;
    for(int i=0;i<np;++i)
      send_counts[i] = pack_size(dom_split[i]) ;

    std::vector<int> send_displs(np,0) ;
    for(int i=1;i<np;++i)
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;

    // communicate to get the receive counts
    std::vector<int> recv_counts(np,0) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, comm) ;

    std::vector<int> recv_displs(np,0) ;
    for(int i=1;i<np;++i)
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;

    // then pack the buffer
    int tot_send_size = 0 ;
    for(int i=0;i<np;++i)
      tot_send_size += send_counts[i] ;

    unsigned char* send_buffer = new unsigned char[tot_send_size] ;
    unsigned char** send_ptr = new unsigned char*[np] ;
    send_ptr[0] = send_buffer ;
    for(int i=1;i<np;++i)
      send_ptr[i] = send_ptr[i-1] + send_counts[i-1] ;

    for(int i=0;i<np;++i) {
      int position = 0 ;
      pack(send_ptr[i], position, send_counts[i], dom_split[i]) ;
    }
    delete[] send_ptr ;

    ///////////////////////////////////////////////////////////////
    // compute the receive sequence for proper unpacking
    // this involves another MPI_Alltoallv
    std::vector<sequence> pack_seq(np) ;
    for(int i=0;i<np;++i)
      pack_seq[i] = sequence(dom_split[i]) ;

    std::vector<sequence> unpack_seq
      = transpose_sequence(pack_seq, comm) ;
    // done, release useless buffers
    std::vector<sequence>().swap(pack_seq) ;
    ///////////////////////////////////////////////////////////////

    // we now proceed to communicate the store contents
    int tot_recv_size = 0 ;
    for(int i=0;i<np;++i)
      tot_recv_size += recv_counts[i] ;

    unsigned char* recv_buffer = new unsigned char[tot_recv_size] ;
    MPI_Alltoallv(&send_buffer[0], &send_counts[0],
                  &send_displs[0], MPI_PACKED,
                  &recv_buffer[0], &recv_counts[0],
                  &recv_displs[0], MPI_PACKED, comm) ;

    delete[] send_buffer ;

    dstore<T> ns ;
    entitySet new_domain ;
    for(int i=0;i<np;++i)
      new_domain += entitySet(unpack_seq[i]) ;
    // if there's any old domain not being redistributed,
    // we will add it also
    entitySet old_dom = dom - total_dom ;
    ns.allocate(new_domain+old_dom) ;
    storeRepP srp = ns.Rep() ;

    unsigned char** recv_ptr = new unsigned char*[np] ;
    recv_ptr[0] = recv_buffer ;
    for(int i=1;i<np;++i)
      recv_ptr[i] = recv_ptr[i-1] + recv_counts[i-1] ;

    // unpack
    for(int i=0;i<np;++i) {
      int position = 0 ;
      srp->unpack(recv_ptr[i], position, recv_counts[i], unpack_seq[i]) ;
    }
    delete[] recv_ptr ;
    delete[] recv_buffer ;

    // copy old_dom if any
    for(entitySet::const_iterator ei=old_dom.begin();
        ei!=old_dom.end();++ei)
      ns[*ei] = attrib_data[*ei] ;
    
    return srp ;
  }
  
  template<class T> storeRepP
  dstoreRepI<T>::redistribute(const std::vector<entitySet>& dom_ptn,
                              const dMap& remap, MPI_Comm comm) {
    // this is a push operation, thus the send, recv are reversed
    std::vector<P2pCommInfo> send, recv ;
    get_p2p_comm(dom_ptn, domain(), 0, 0, comm, recv, send) ;
    dstore<T> new_store ;
#ifdef DYNAMICSCHEDULING
    fill_store2(getRep(), 0, new_store.Rep(), &remap, send, recv, comm) ;
#else
    Loci::Abort();
#endif
    return new_store.Rep() ;
  }

#ifdef DYNAMICSCHEDULING
  template<class T> storeRepP
  dstoreRepI<T>::redistribute_omd(const std::vector<entitySet>& dom_ptn,
                                  const dMap& remap, MPI_Comm comm) {
    // this is a push operation, thus the send, recv are reversed
    std::vector<P2pCommInfo> send, recv ;
    get_p2p_comm(dom_ptn, domain(), 0, 0, comm, recv, send) ;
    dstore<T> new_store ;
    fill_store_omd(getRep(), 0, new_store.Rep(), &remap, send, recv, comm) ;
    return new_store.Rep() ;
  }
#endif
  
  // ********************************************************************/

  template<class T>
  storeRepP dstoreRepI<T>::freeze() {
    store<T> static_store ;
    entitySet tmp_dom = domain();
    static_store.allocate(tmp_dom) ;
    for(entitySet::const_iterator ei = tmp_dom.begin();
        ei != tmp_dom.end(); ++ei)
      static_store[*ei] = attrib_data[*ei] ;
    return static_store.Rep() ;
  }

  template<class T>
  storeRepP dstoreRepI<T>::thaw() {
    return getRep() ;
  }
  // **********************************************************************/
  
  template<class T> 
  void dstoreRepI<T>::copy(storeRepP &st, const entitySet &context)  
  {
    const_dstore<T> s(st) ;
    fatal((context-domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i] = s[i] ;
    } ENDFORALL ;

  }

  //*************************************************************************/

  template<class T> 
  void dstoreRepI<T>::gather(const dMap &m, storeRepP &st, const entitySet &context) 
  {
    const_dstore<T> s(st) ;

    fatal( context != EMPTY ) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i] = s[m[i]] ;
    } ENDFORALL ;

  }

  //**************************************************************************/

  template<class T> 
  void dstoreRepI<T>::scatter(const dMap &m, storeRepP &st, const entitySet &context) 
  {
    const_dstore<T> s(st) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);

    FORALL(context,i) {
      attrib_data[m[i]] = s[i] ;
    } ENDFORALL ;

  }

  //*************************************************************************/
  template <class T> 
  int dstoreRepI<T>::pack_size( const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size( traits_type, eset );
  }
  
  template <class T> 
  int dstoreRepI<T>::estimated_pack_size( const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
    
    return get_estimated_mpi_size( traits_type, eset );  
  }

  
  template<class T> int dstoreRepI<T>::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;

    typedef typename data_schema_traits<T>::Schema_Converter
      schema_converter;
    schema_converter traits_type;

    return get_mpi_size(traits_type, packed);
  }

  //*******************************************************************/

  template <class T>
  int dstoreRepI<T>::get_mpi_size( IDENTITY_CONVERTER c,
                                   const entitySet &eset)
  {
    T *p = 0 ;
    return cpypacksize(p,eset.size());
  }
  

  //*******************************************************************/

  template <class T>
  int dstoreRepI<T>::get_estimated_mpi_size( IDENTITY_CONVERTER c,
                                             const entitySet &eset)
  {
    T *p = 0 ;
    return cpypacksize(p,eset.size());
  }
  
  //*******************************************************************/
  template <class T>
  int dstoreRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c,
                                   const entitySet &eset)
  {
    entitySet :: const_iterator ci;
    typedef data_schema_traits<T> converter_traits;

    int numBytes = 0 ;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      typename converter_traits::Converter_Type cvtr(attrib_data[*ci]);
      int size      = cvtr.getSize();
      numBytes += cpypacksize(&size,1) ;
      typename converter_traits::Converter_Base_Type *p=0 ;
      numBytes += cpypacksize(p,size) ;
    }
    return numBytes ;
  }
  //*******************************************************************/
  template <class T>
  int dstoreRepI<T>::get_estimated_mpi_size( USER_DEFINED_CONVERTER c,
                                             const entitySet &eset)
  {
    int numBytes = 0 ;
    numBytes = eset.size()*50*sizeof(double);
    numBytes  += eset.size()*sizeof(int);
    return(numBytes) ;
   
  }
  
  //*******************************************************************/
  template <class T> 
  void dstoreRepI<T>::pack( void *outbuf, int &position, int &size, 
                            const entitySet &usr_eset )  
  {
    entitySet ecommon;

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    ecommon = domain()&usr_eset;

    packdata( traits_type, outbuf, position, size, ecommon);

  }
 
  //*******************************************************************/

  template <class T> 
  void dstoreRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
                               int &position, int outcount,
                               const entitySet &eset )  
  {

    entitySet :: const_iterator ci;

    for( ci = eset.begin(); ci != eset.end(); ++ci)
      cpypack(outbuf,position,outcount,&attrib_data[*ci],1) ;
  }

  //*******************************************************************/

  template <class T>
  void dstoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                                int &position, int outcount, 
                                const entitySet &eset )
  {
    entitySet :: const_iterator ci;

    typedef  data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;
    std::vector<dtype>  inbuf;

    //-------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-------------------------------------------------------------------
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( attrib_data[*ci]);

      int stateSize    = cvtr.getSize();
      if( stateSize > static_cast<int>(inbuf.size()) ) inbuf.resize(stateSize);

      cvtr.getState( &inbuf[0], stateSize);
      cpypack(outbuf,position,outcount,&stateSize,1) ;
      cpypack(outbuf,position,outcount,&inbuf[0],stateSize) ;
    }
  }

  //*******************************************************************/

  template <class T> 
  void dstoreRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    unpackdata(traits_type, ptr, loc, size, seq); 
  }

  //*********************************************************************/
  template <class T> 
  void dstoreRepI<T>::unpackdata(IDENTITY_CONVERTER c, void *inbuf,
                                 int &position,  int insize,
                                 const sequence &seq) 
  {
    sequence:: const_iterator ci;

    for( ci = seq.begin(); ci != seq.end(); ++ci)
      cpyunpack(inbuf,position,insize,&attrib_data[*ci],1) ;
  }
  //*********************************************************************/
  template <class T> 
  void dstoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
                                  int &position, int insize, const sequence &seq) 
  {
    sequence :: const_iterator ci;

    typedef  data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;
    std::vector<dtype>  outbuf;

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      int stateSize=0;
      cpyunpack(inbuf,position,insize,&stateSize,1) ;

      if( stateSize > static_cast<int>(outbuf.size()) ) outbuf.resize(stateSize);
 
      cpyunpack(inbuf,position,insize,&outbuf[0],stateSize) ;

      typename schema_traits::Converter_Type cvtr( attrib_data[*ci] );
      cvtr.setState( &outbuf[0], stateSize);
    }

  }

  template<class T> 
  frame_info dstoreRepI<T>::get_frame_info() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }
  template<class T> 
  frame_info dstoreRepI<T>::get_frame_info(IDENTITY_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  template<class T> 
  frame_info dstoreRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  template<class T> 
  DatatypeP dstoreRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T> 
  DatatypeP dstoreRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T> 
  DatatypeP dstoreRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  
  //*********************************************************************/
  template<class T> 
  void dstoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset)
  {
    warn(true) ;
    /*
      typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
      schema_converter traits_type;

      entitySet eset, ecommon;

      HDF5_ReadDomain(group_id, eset);

      ecommon = eset;
      allocate( ecommon );

      hdf5read( group_id, traits_type, eset,  ecommon );
    */
  }
  //*********************************************************************/
#ifdef H5_HAVE_PARALLEL 
  template<class T> 
  void dstoreRepI<T>::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset, hid_t xfer_plist_id) 
  {
    warn(true) ;
  }
#endif
  //*************************************************************************/

  template<class T>
  void  dstoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset)
  {
    warn(true) ;
    /*
      hsize_t dimension;
      size_t indx = 0, arraySize;
      int    rank = 1;

      entitySet::const_iterator ci;
      typedef data_schema_traits<T> traits_type;

      int num_intervals = user_eset.num_intervals();
      interval *it = new interval[num_intervals];

      for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

      //------------------------------------------------------------------------
      // Calculate the offset 
      //------------------------------------------------------------------------

      store<unsigned> offset;
      offset.allocate( eset );

      arraySize = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci)
      offset[*ci] = arraySize++;

      //------------------------------------------------------------------------

      DatatypeP dtype = traits_type::get_type();
      hid_t vDatatype = dtype->get_hdf5_type();

      dimension = arraySize;
      hid_t mDataspace = H5Screate_simple(rank, &dimension, NULL);   // memory  dataspace
      hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
      hid_t vDataset   = H5Dopen( group_id, "VariableData");

      hssize_t  start[]     = {0};  // determines the starting coordinates.
      hsize_t   stride[]    = {1};  // which elements are to be selected.
      hsize_t   block[]     = {1};  // size of element block;
      hssize_t  foffset[]   = {0};  // location (in file) where data is read.
      hsize_t   count[]     = {0};  // how many positions to select from the dataspace

      std::vector<T>  data;
      for( int k = 0; k < num_intervals; k++) {
      count[0] = it[k].second - it[k].first + 1;

      if( count[0] > data.size() ) data.resize(count[0] );

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace, 
      H5P_DEFAULT, &data[0]);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) 
      attrib_data[i] = data[indx++];
      }

      H5Dclose( vDataset   );
      H5Tclose( vDatatype  );
      H5Sclose( mDataspace );
      H5Sclose( vDataspace );
    */
  }

#ifdef H5_HAVE_PARALLEL 
  template<class T>
  void  dstoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset,hid_t xfer_plist_id)
  {
    warn(true) ;
  }
#endif
  //*************************************************************************/

  template<class T>
  void  dstoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset)
  {
    warn(true) ;
    /*
      hsize_t  dimension;
      size_t   indx = 0, arraySize;
      hid_t    vDataset, vDataspace, vDatatype, mDataspace;
      entitySet::const_iterator  ci;

      //---------------------------------------------------------------
      // Size of each Bucket ....
      //---------------------------------------------------------------
      vDatatype  = H5T_NATIVE_INT;
      vDataset   = H5Dopen(group_id,"ContainerSize");
      vDataspace = H5Dget_space(vDataset);
      H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

      std::vector<int> ibuf(dimension);
      H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);
      H5Dclose(vDataset  );
      H5Sclose(vDataspace);

      store<int> container;
      container.allocate( eset );

      indx      = 0;
      arraySize = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      container[*ci] = ibuf[indx++];
      arraySize     += container[*ci];
      }

      //-------------------------------------------------------------------
      // Calculate the offset of each entity in file ....
      //-------------------------------------------------------------------
      store<unsigned>   offset;
      offset.allocate( eset );

      indx     = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = indx;
      indx       += container[*ci];
      }

      //------------------------------------------------------------------
      // Read the data now ....
      //------------------------------------------------------------------
      int num_intervals = usr_eset.num_intervals();
      if( num_intervals == 0) {
      std::cout << "Warning: Number of intervals are zero : " << endl;
      return;
      }

      interval *it = new interval[num_intervals];

      for(int i=0;i< num_intervals;i++) it[i] = usr_eset[i];

      typedef typename data_schema_traits<T>::Converter_Base_Type dtype;
      typedef data_schema_traits<dtype> traits_type;
      DatatypeP atom_type = traits_type::get_type() ;
      vDatatype = atom_type->get_hdf5_type();

      dimension  = arraySize;
      vDataset   = H5Dopen(group_id,"VariableData");
      vDataspace = H5Dget_space(vDataset);
      mDataspace = H5Dget_space(vDataset);
      H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

      hssize_t  start_mem[] = {0};  // determines the starting coordinates.
      hsize_t   stride[]    = {1};  // which elements are to be selected.
      hsize_t   block[]     = {1};  // size of element block;
      hssize_t  foffset[]   = {0};  // location (in file) where data is read.
      hsize_t   count[]     = {0};  // how many positions to select from the dataspace

      std::vector<dtype> data;

      for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
      count[0] +=  container[i];
      if( count[0] > data.size() ) data.resize( count[0] );

      foffset[0] = offset[it[k].first];
      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start_mem, stride,
      count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,   stride,
      count, block);
      H5Dread(vDataset, vDatatype, mDataspace, vDataspace,H5P_DEFAULT, &data[0]);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
      typename data_schema_traits<T>::Converter_Type cvtr( attrib_data[i]);
      cvtr.setState( &data[0]+indx, container[i] );
      indx += container[i];
      }
      }

      H5Tclose(vDatatype );
      H5Dclose(vDataset  );
      H5Sclose(vDataspace);
      H5Sclose(mDataspace);
    */
  }

#ifdef H5_HAVE_PARALLEL 
  template<class T>
  void  dstoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    warn(true) ;
  }
#endif
  //*************************************************************************/

  template<class T> 
  void dstoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& usr_eset) const
  {
    warn(true) ;
    /*
      typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
      schema_converter traits_output_type;

      entitySet eset(usr_eset &domain());

      HDF5_WriteDomain(group_id, eset);

      hdf5write(group_id, traits_output_type, eset);
    */
  }
#ifdef H5_HAVE_PARALLEL 
  template<class T> 
  void dstoreRepI<T>::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& usr_eset, hid_t xfer_plist_id ) const
  {
    warn(true) ;
  }
#endif
  //*************************************************************************/
  template <class T> 
  void dstoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &eset)  const
  {
    warn(true) ;
    /*
      int      rank = 1;
      hsize_t  dimension;

      std::vector<T>   newvec;
      entitySet :: const_iterator  ei;

      int arraySize = eset.size();

      if( arraySize < 1) return;

      //------------------------------------------------------------------------
      // Collect state data from each object and put into 1D array
      //------------------------------------------------------------------------
      std::vector<T>  data(arraySize);

      int indx = 0;
      for( ei = eset.begin(); ei != eset.end(); ++ei) {
      data[indx++]  = attrib_data.elem(*ei) ;
      }

      //-------------------------------------------------------------------------
      // Write (variable) Data into HDF5 format
      //-------------------------------------------------------------------------

      typedef data_schema_traits<T> traits_type;
      DatatypeP dtype = traits_type::get_type();
      hid_t vDatatype = dtype->get_hdf5_type();

      dimension        = arraySize;
      hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
      hid_t vDataset   = H5Dcreate(group_id, "VariableData", vDatatype, vDataspace, H5P_DEFAULT);
      H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

      H5Dclose( vDataset  );
      H5Sclose( vDataspace);
      H5Tclose( vDatatype );
    */
  }
#ifdef H5_HAVE_PARALLEL 
  template <class T> 
  void dstoreRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &eset, hid_t xfer_plist_id )  const
  {
    warn(true) ;
  }
#endif
  //*************************************************************************/

  template <class T> 
  void dstoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &eset)  const
  {   
    warn(true) ;
    /*
      T         Obj;
      hid_t     vDataspace, vDataset, vDatatype;
      int       rank = 1;
      hsize_t   dimension;

      entitySet :: const_iterator ci;

      //-------------------------------------------------------------------------
      // Get the sum of each object size and maximum size of object in the 
      // container for allocation purpose
      //-------------------------------------------------------------------------

      typedef data_schema_traits<T> schema_traits ;

      size_t  arraySize= 0;
      int     stateSize, maxStateSize = 0;

      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      Obj = attrib_data.elem(*ci) ;
      typename schema_traits::Converter_Type cvtr( Obj );
      stateSize    = cvtr.getSize();
      arraySize   += stateSize;
      maxStateSize = max( maxStateSize, stateSize );
      } 

      //-----------------------------------------------------------------------------
      // Collect state data from each object and put into 1D array
      //-----------------------------------------------------------------------------
      typedef typename schema_traits::Converter_Base_Type dtype;

      dtype *data = new dtype[arraySize];

      size_t indx = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      Obj = attrib_data.elem(*ci) ;
      typename schema_traits::Converter_Type cvtr( Obj );
      cvtr.getState( data+indx, stateSize);
      indx +=stateSize ;
      }

      //--------------------------------------------------------------------
      // Write (variable) Data into HDF5 format
      //--------------------------------------------------------------------
      typedef data_schema_traits<dtype> traits_type;
      DatatypeP atom_type = traits_type::get_type() ;
      vDatatype = atom_type->get_hdf5_type();

      dimension  = arraySize;
      vDataspace = H5Screate_simple(rank, &dimension, NULL);
      vDataset   = H5Dcreate(group_id, "VariableData", vDatatype,
      vDataspace, H5P_DEFAULT);
      H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

      H5Dclose( vDataset  );
      H5Sclose( vDataspace);
      H5Tclose( vDatatype );

      delete [] data;

      //--------------------------------------------------------------------
      // Write Container 
      //--------------------------------------------------------------------
      dimension    = eset.size();
      std::vector<int> vbucket(eset.size());

      indx = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      Obj = attrib_data.elem(*ci) ;
      typename schema_traits::Converter_Type cvtr( Obj );
      vbucket[indx++] =  cvtr.getSize();
      }

      vDatatype  = H5T_NATIVE_INT;
      vDataspace = H5Screate_simple(rank, &dimension, NULL);
      vDataset   = H5Dcreate(group_id, "ContainerSize",
      vDatatype, vDataspace, H5P_DEFAULT);
      H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vbucket[0]);

      H5Dclose( vDataset  );
      H5Sclose( vDataspace);
    */
  }

#ifdef H5_HAVE_PARALLEL 
  template <class T> 
  void dstoreRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &eset, hid_t xfer_plist_id)  const
  {   
    warn(true) ;
  }
#endif
  //*************************************************************************/
  
} // end of namespace Loci

#endif
