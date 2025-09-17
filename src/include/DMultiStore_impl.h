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
#ifndef DMULTISTORE_IMPL_H
#define DMULTISTORE_IMPL_H

#include <DMultiStore_def.h>
#include <multiStore_def.h>

namespace Loci {
  
  template <class T> 
  inline std::ostream & operator<<(std::ostream &s, const dmultiStore<T> &m)
  { return m.Print(s) ; }

  //*************************************************************************/
  template<class T> 
  inline std::istream & operator>>(std::istream &s, dmultiStore<T> &m)
  { return m.Input(s) ; }

  //**************************************************************************/
  template<class T> 
  void dmultiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }
  //**************************************************************************/
  template<class T> 
  store_instance::instance_type
  const_dmultiStore<T>::access() const
  { return READ_ONLY ; }

  //**************************************************************************/
  template<class T> 
  void const_dmultiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      attrib_data = p->get_attrib_data() ;
    warn(p == 0) ;
  }

  //**************************************************************************/
  template<class T> 
  void dmultiStoreRepI<T>::allocate(const store<int> &sizes) 
  {
    //------------------------------------------------------------------------/
    // Objective : reserve the memory of multiStore using Store..
    //------------------------------------------------------------------------/

    entitySet eset = sizes.domain() ;
    entitySet :: const_iterator  ci;

    for( ci = eset.begin(); ci != eset.end(); ++ci)
      attrib_data[*ci].resize(sizes[*ci]);

    store_domain = eset ;
    dispatch_notify() ;
  }
  //***************************************************************************/
  
  template<class T> 
  void dmultiStoreRepI<T>::allocate(const entitySet &eset) 
  {
    entitySet redundant, newSet;
    entitySet :: const_iterator  ci ;
    redundant = domain() -  eset;
    newSet    = eset - domain();
    for( ci = redundant.begin(); ci != redundant.end(); ++ci)
      attrib_data.erase(*ci);
    
    std::vector<T>  emptyVec;
    for( ci = newSet.begin(); ci != newSet.end(); ++ci)
      attrib_data[*ci] =  emptyVec;
    
    store_domain = eset;
    dispatch_notify() ;
  }

  template<class T>
  void dmultiStoreRepI<T>::shift(int_type offset) {
    store_domain >>= offset ;
    allocate(store_domain) ;
  }
  //**************************************************************************/


  template<class T> 
  dmultiStoreRepI<T>::~dmultiStoreRepI() 
  {
    attrib_data.clear();
  }

  //**************************************************************************/

  template<class T> 
  storeRep *dmultiStoreRepI<T>::new_store(const entitySet &p) const 
  {
    return new dmultiStoreRepI<T>(p) ;
  }
  template<class T> 
  storeRep *dmultiStoreRepI<T>::new_store(const entitySet &p, const int* cnt) const 
  {
    store<int> count ;
    count.allocate(p) ;
    int t= 0 ;
    FORALL(p, pi) {
      count[pi] = cnt[t++] ; 
    } ENDFORALL ;
    return new dmultiStoreRepI<T>(count) ;
  }
  //**************************************************************************/
  
  template<class T> 
  storeRepP dmultiStoreRepI<T>::remap(const dMap &m) const {
    
    dmultiStore<T> s ;
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  template<class T>
  storeRepP dmultiStoreRepI<T>::freeze() {
    multiStore<T> static_mul ;
    store<int> count ;
    entitySet tmp_dom = domain() ;
    count.allocate(tmp_dom) ;
    for(entitySet::const_iterator ei = tmp_dom.begin();
        ei != tmp_dom.end(); ++ei)
      count[*ei] = attrib_data[*ei].size() ;
    static_mul.allocate(count) ;
    for(entitySet::const_iterator ei = tmp_dom.begin();
        ei != tmp_dom.end(); ++ei) {
      int i = 0 ;
      for(typename std::vector<T>::const_iterator vi=attrib_data[*ei].begin();
          vi != attrib_data[*ei].end(); ++vi,++i) {
	static_mul[*ei][i] = *vi ;
      }
    }
    return static_mul.Rep() ;
  }

  template<class T>
  storeRepP dmultiStoreRepI<T>::thaw() {
    return getRep() ;
  }
  
  //**************************************************************************/
  template<class T> 
  void dmultiStoreRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
    const_dmultiStore<T> s(st) ;
    std::vector<T>    newVec;

    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[i];
      attrib_data[i].reserve( newVec.size() );
      for(typename std::vector<T>::size_type j = 0; j < newVec.size(); j++)
        attrib_data[i].push_back( newVec[j] );
    } ENDFORALL ;

    dispatch_notify() ;

  }

  //**************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::gather(const dMap &m, storeRepP &st, const entitySet &context) 
  {
    const_dmultiStore<T> s(st) ;
    std::vector<T>    newVec;

    FORALL(context,i) {
      attrib_data[i].clear();
      newVec  =   s[m[i]];
      attrib_data[i].reserve( newVec.size() );
      for(typename std::vector<T>::size_type i = 0; i < newVec.size(); i++) 
        attrib_data[i].push_back( newVec[i] );
    } ENDFORALL ;

    dispatch_notify() ;

  }

  //**************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::scatter(const dMap &m, storeRepP &st, const entitySet &context) 
  {
    const_dmultiStore<T> s(st) ;

    fatal((context - s.domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
    
    FORALL(context,i) {
      attrib_data[m[i]] = s[i];
    } ENDFORALL ;

    dispatch_notify() ;
  }

  //**************************************************************************/
 
  template<class T> 
  store_type dmultiStoreRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //***************************************************************************/
  
  template<class T> 
  entitySet dmultiStoreRepI<T>::domain() const 
  {
    typename HASH_MAP(int,std::vector<T> ) :: const_iterator    ci;
    entitySet          storeDomain;
    std::vector<int>        vec;

    for( ci = attrib_data.begin(); ci != attrib_data.end(); ++ci )
      vec.push_back( ci->first ) ;

    std::sort( vec.begin(), vec.end() );

    for(std::vector<int>::size_type i = 0; i < vec.size(); i++)
      storeDomain +=  vec[i];

    return storeDomain ;
  }

  //***************************************************************************/
  
  template<class T> 
  std::ostream &dmultiStoreRepI<T>::Print(std::ostream &s) const 
  {
    s << '{' << domain() << endl ;

    typename HASH_MAP(int,std::vector<T>)::const_iterator ci;

    FORALL(domain(),ii) {
      ci =  attrib_data.find(ii);
      if( ci == attrib_data.end() ) continue;
      s << ci->second.size() << std::endl ;
    } ENDFORALL ;

    FORALL(domain(),ii) {
      ci =  attrib_data.find(ii);
      if( ci == attrib_data.end() ) continue;
      const std::vector<T> &vec =  ci->second;
      Loci::streamoutput(&vec[0],vec.size(),s) ;
    } ENDFORALL ;

    s << '}' << std::endl ;
    return s ;
  }

  //**************************************************************************/

  template<class T> 
  std::istream &dmultiStoreRepI<T>::Input(std::istream &s) 
  {

    entitySet  e ;
    char ch ;

    // Read the opening brackets ...
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    // Read the interval set ....
    s >> e ;

    // Read the size of each entity map.
    store<int> sizes ;
    sizes.allocate(e) ;
    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;
    allocate(sizes) ;

    // Read the attribute data
    std::vector<T>  vec;
    FORALL(e,ii) {
      vec.clear();
      for( int i = 0; i < sizes[ii]; i++) {
        T val ;
        Loci::streaminput(&val,1,s) ;
        vec.push_back(val);
      }
      attrib_data[ii] = vec;
    } ENDFORALL ;
            
    // Close the bracket ..
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }

    return s ;

  }

  //***************************************************************************/
  template <class T> 
  int dmultiStoreRepI<T>::pack_size(const entitySet &e ) 
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size( traits_type, e );
  }
  
  //***************************************************************************/
  template <class T> 
  int dmultiStoreRepI<T>::estimated_pack_size(const entitySet &e ) 
  {
    
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
    
    return get_estimated_mpi_size( traits_type, e );
  }
  

  //***************************************************************************/
  template<class T> int dmultiStoreRepI<T>::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size(traits_type, packed);
  }
  //**************************************************************************/

  template <class T>
  int dmultiStoreRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {
    int packsize = 0 ;
    FORALL(eset,i) {
      int size = attrib_data[i].size() ;
      packsize += cpypacksize(&size,1) ;
      packsize += cpypacksize(&attrib_data[i][0],size);
    } ENDFORALL ;
    
    return packsize ;
  }
  //**************************************************************************/

  template <class T>
  int dmultiStoreRepI<T>::get_estimated_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {
    
    
    return( 4*eset.size()*sizeof(T) + eset.size()*sizeof(int) ) ;
  }
  //**************************************************************************/

  template <class T>
  int dmultiStoreRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {
    std::vector<T>   newVec;
    entitySet :: const_iterator  ei;
    typename HASH_MAP(int,std::vector<T>)::const_iterator ci;

    typedef data_schema_traits<T> schema_traits;

    int packsize = 0 ;
    for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
        newVec      = ci->second;
        int sz = newVec.size() ;
        packsize += cpypacksize(&sz,1) ;
        for( int ivec = 0; ivec < sz; ivec++) {
          typename schema_traits::Converter_Type cvtr(newVec[ivec] );
          int arraySize = cvtr.getSize() ;
          typename schema_traits::Converter_Base_Type *p = 0 ;
          packsize += cpypacksize(&arraySize,1) ;
          packsize += cpypacksize(p,arraySize) ;
        }
      }
    }

    return packsize ;

  }
  //**************************************************************************/

  template <class T>
  int dmultiStoreRepI<T>::get_estimated_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {
  
    int  vsize = 4*eset.size()*50*sizeof(double)+
      4*eset.size()*sizeof(int)  +   // size of each object
      eset.size()*sizeof(int);       // size of each entity

    return(vsize);
  }

  
  //**************************************************************************/
  
  template <class T> 
  void dmultiStoreRepI<T>::pack(void *ptr, int &loc, int &size, const entitySet &eset ) 
  {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
  
    packdata( traits_type, ptr, loc, size, eset);
  }
  
  //**************************************************************************/
  template <class T>
  void dmultiStoreRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf, int &position,
                                     int outcount, const entitySet &eset )
  {
    int vsize ;
    std::vector<T>   inbuf;
    entitySet :: const_iterator   ci;
    typename HASH_MAP(int,std::vector<T> )::const_iterator iter;

    for( ci = eset.begin(); ci != eset.end(); ++ci){
      iter = attrib_data.find(*ci);
      if( iter != attrib_data.end() ) {
        inbuf   = iter->second;
        vsize   = inbuf.size();
        cpypack(outbuf,position,outcount,&vsize,1) ;
        cpypack(outbuf,position,outcount,&inbuf[0],vsize) ;
      } 
    }

  }

  //**************************************************************************/

  template <class T> 
  void dmultiStoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                                     int &position, int outcount, 
                                     const entitySet &eset ) 
  {
    entitySet::const_iterator ci;
    int vsize, stateSize;
    
    typedef data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;

    std::vector<dtype> inbuf;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vsize = attrib_data[*ci].size();
      cpypack(outbuf,position,outcount,&vsize,1) ;
      for( int ivec = 0; ivec < vsize; ivec++){
        typename data_schema_traits<T>::Converter_Type
          cvtr( attrib_data[*ci][ivec] );

        stateSize = cvtr.getSize();
        inbuf.resize(stateSize);

        cvtr.getState( &inbuf[0], stateSize);

        cpypack(outbuf,position,outcount,&stateSize,1) ;
        cpypack(outbuf,position,outcount,&inbuf[0],stateSize) ;
      }
    }
  }
  //**************************************************************************/

  
  template <class T> 
  void dmultiStoreRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq) 
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    entitySet ecommon, ediff,eset(seq);

    ediff = eset - domain();

    if( ediff.size() > 0)
      std::cerr << " Warning: Entities not part of domain, and not unpacked " << ediff <<endl;

    unpackdata( traits_type, ptr, loc, size, seq);
  }
  

  //**************************************************************************/
  template <class T> 
  void dmultiStoreRepI<T>::unpackdata( IDENTITY_CONVERTER c, void *inbuf, int &position, 
                                       int &insize, const sequence &seq)
  {

    unsigned int   vsize;
    std::vector<T>  outbuf;
    sequence :: const_iterator ci;

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      cpyunpack(inbuf,position,insize,&vsize,1) ;

      if( vsize > outbuf.size() ) outbuf.resize(vsize);

      attrib_data[*ci].resize(vsize);
      cpyunpack(inbuf,position,insize,&outbuf[0],vsize) ;

      for(typename std::vector<T>::size_type ivec = 0; ivec < vsize; ivec++) 
        attrib_data[*ci][ivec] = outbuf[ivec];
    }
  }

  //***********************************************************************/
  template <class T> 
  void dmultiStoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
                                       int &position, int &insize, const sequence &seq)
  {

    sequence :: const_iterator ci;
    int  stateSize, vsize;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    std::vector<dtype> outbuf;

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      if( !store_domain.inSet( *ci ) ) {
        std::cerr << "Warning: Entity not present in entityset " << *ci << endl;
        continue;
      }
      cpyunpack(inbuf,position,insize,&vsize,1) ;

      attrib_data[*ci].resize(vsize);

      for( int ivec = 0; ivec < vsize; ivec++) {
        cpyunpack(inbuf,position,insize,&stateSize,1) ;
        if( stateSize > outbuf.size() ) outbuf.resize(stateSize);
        cpyunpack(inbuf,position,insize,&outbuf[0],stateSize) ;

        typename schema_traits::Converter_Type  cvtr( attrib_data[*ci][ivec] );
        cvtr.setState( &outbuf[0], stateSize);
      }
    }

  }
 
  //***************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset) 
  {
    warn(true) ;
    /*
      typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
      schema_converter traits_type;

      entitySet eset, ecommon;

      // Read the entitySet available in file ...
      HDF5_ReadDomain(group_id, eset);

      // Intersection of entityset in file and user defined. Only common entitySet
      // are read from the file ...
      //
      ecommon = eset & user_eset ;   
      allocate( ecommon );

      // Read the common entitities ...
      hdf5read( group_id, traits_type, eset,  ecommon );
    */
  }
  template<class T> 
  void dmultiStoreRepI<T>::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset, hid_t xfer_plist_id) 
  {
    warn(true) ;
  }
  //*************************************************************************/

  template <class T> 
  void dmultiStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &user_eset)
  {
    warn(true) ;
    /*
      hsize_t  dimension;
      hid_t    vDataset, vDataspace, vDatatype, mDataspace;

      int indx = 0, arraySize, vsize;
      int    rank = 1;

      entitySet::const_iterator ci;

      vDatatype  = H5Tcopy(H5T_NATIVE_INT);
      vDataset   = H5Dopen(group_id,"ContainerSize");
      vDataspace = H5Dget_space(vDataset);
      H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

      std::vector<int> ibuf(dimension);
      H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

      //---------------------------------------------------------------------
      // Size of each main container....
      //---------------------------------------------------------------------

      store<int> container;
      container.allocate( eset );

      indx      = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) 
      container[*ci] = ibuf[indx++];

      container.Print(std::cout);

      //-----------------------------------------------------------------------
      // Calculate the offset of each entity in file ....
      //-----------------------------------------------------------------------
      store<unsigned>   offset;
      store<int>        bucket;

      offset.allocate( eset );

      bucket.allocate( user_eset );
      for( ci = user_eset.begin(); ci != user_eset.end(); ++ci) 
      bucket[*ci] = container[*ci];
      allocate( bucket );

      arraySize = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      arraySize  += container[*ci];
      }
      //------------------------------------------------------------------------
      // Read the data now ....
      //------------------------------------------------------------------------
      int num_intervals = user_eset.num_intervals();
      interval *it = new interval[num_intervals];

      for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

      dimension = arraySize;

      typedef data_schema_traits<T> traits_type;

      DatatypeP dtype = traits_type::get_type();
      vDatatype = dtype->get_hdf5_type();

      mDataspace = H5Screate_simple(rank, &dimension, NULL);
      vDataspace = H5Screate_simple(rank, &dimension, NULL);
      vDataset   = H5Dopen( group_id, "VariableData");

      hssize_t  start[]     = {0};  // determines the starting coordinates.
      hsize_t   stride[]    = {1};  // which elements are to be selected.
      hsize_t   block[]     = {1};  // size of element block;
      hssize_t  foffset[]   = {0};  // location (in file) where data is read.
      hsize_t   count[]     = {0};  // how many positions to select from the dataspace

      std::vector<T>  data;
      for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
      count[0] +=  container[i];

      if( count[0] > data.size() ) data.resize( count[0] );

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, &data[0]);

      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
      vsize = container[i];
      attrib_data[i].resize(vsize);
      for( int ivec = 0; ivec < vsize; ivec++)
      attrib_data[i][ivec] = data[indx++];
      }
      }

      H5Dclose( vDataset  );
      H5Sclose( vDataspace);
      H5Sclose( mDataspace);
      H5Tclose( vDatatype );
    */
  }
  template <class T> 
  void dmultiStoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &user_eset, hid_t xfer_plist_id)
  {
    warn(true) ;
  }
  //*************************************************************************/
  template<class T>
  void  dmultiStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &user_eset)
  {
    warn(true) ;
    /*
      hsize_t  dimension;
      hid_t    vDataset, vDataspace, vDatatype, mDataspace;

      int       indx=0, rank=1, vecsize, vsize;
      unsigned  arraySize;

      entitySet::const_iterator ci;

      vDatatype  = H5Tcopy(H5T_NATIVE_INT);
      vDataset   = H5Dopen(group_id,"ContainerSize");
      vDataspace = H5Dget_space(vDataset);
      H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

      std::vector<int> ibuf(dimension);
      H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

      //-----------------------------------------------------------------------
      // Size of each main container....
      //------------------------------------------------------------------------
      store<int> container;
      container.allocate( eset );
      indx      = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci)
      container[*ci] = ibuf[indx++];

      vDatatype  = H5Tcopy(H5T_NATIVE_INT);
      vDataset   = H5Dopen(group_id,"SubContainerSize");
      vDataspace = H5Dget_space(vDataset);
      H5Sget_simple_extent_dims(vDataspace, &dimension, NULL);

      ibuf.resize(dimension);
      H5Dread(vDataset, vDatatype, H5S_ALL,H5S_ALL,H5P_DEFAULT, &ibuf[0]);

      //---------------------------------------------------------------------
      // Calculate the offset of each entity in file ....
      //---------------------------------------------------------------------
      store< unsigned int>   offset;
      dmultiStore<int>  subContainer;
      offset.allocate( eset );

      arraySize = 0;
      indx      = 0;
      for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      vsize       = container[*ci];
      for( int i = 0; i < vsize; i++)  {
      vecsize    =  ibuf[indx++];
      arraySize  += vecsize;
      subContainer[*ci].push_back( vecsize );
      }
      }
      //---------------------------------------------------------------------
      // Read the data now ....
      //---------------------------------------------------------------------
      int num_intervals = user_eset.num_intervals();

      interval *it = new interval[num_intervals];

      for(int i=0;i< num_intervals;i++) it[i] = user_eset[i];

      typedef data_schema_traits<T> converter_traits;
      typename converter_traits::Converter_Base_Type  *data, *buf, dtype;
      DatatypeP base_type =
      data_schema_traits<typename converter_traits::Converter_Base_Type>::get_type() ;
    
      vDatatype = base_type->get_hdf5_type();

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

      for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++){
      vsize = container[i];
      for( int j = 0; j < vsize; j++)
      count[0] +=  subContainer[i][j];
      }

      data = new typename data_schema_traits<T>::Converter_Base_Type[count[0]];

      foffset[0] = offset[it[k].first];
      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start_mem, stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,   stride, count, block);
      H5Dread(vDataset, vDatatype, mDataspace, vDataspace,H5P_DEFAULT, data);

      indx = 0;
      int bucsize;
      for( int i = it[k].first; i <= it[k].second; i++) {
      vsize = container[i];
      attrib_data[i].resize(vsize);
      for( int j = 0; j < vsize; j++) {
      typename data_schema_traits<T>::Converter_Type cvtr( attrib_data[i][j] );
      bucsize = subContainer[i][j];
      cvtr.setState( data+indx, bucsize );
      indx += bucsize ;
      }
      }

      delete[] data;
      }
    */
  }
  template<class T>
  void  dmultiStoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &user_eset, hid_t xfer_plist_id)
  {
    warn(true) ;
  }
  
  template<class T> 
  frame_info dmultiStoreRepI<T>::get_frame_info() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }
  template<class T> 
  frame_info dmultiStoreRepI<T>::get_frame_info(IDENTITY_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  template<class T> 
  frame_info dmultiStoreRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g) {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  template<class T> 
  DatatypeP dmultiStoreRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T> 
  DatatypeP dmultiStoreRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T> 
  DatatypeP dmultiStoreRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }


  //*************************************************************************/

  template<class T> 
  void dmultiStoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &eset) const 
  {
    warn(true) ;
    /*
      typedef typename data_schema_traits<T> ::Schema_Converter schema_converter;
      schema_converter traits_output_type;

      Loci::HDF5_WriteDomain(group_id, eset);

      hdf5write(group_id, traits_output_type, eset);
    */
  }
  template<class T> 
  void dmultiStoreRepI<T>::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &eset,hid_t xfer_plist_id ) const 
  {
    warn(true) ;
  }
  //*************************************************************************/
  template <class T> 
  void dmultiStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &eset) const
  {
    warn(true) ;
    /*
      int      rank = 1;
      hsize_t  dimension;
      hid_t    vDataspace, vDataset, vDatatype;

      std::vector<T>   newvec;
      entitySet :: const_iterator  ei;
      typename HASH_MAP(int,std::vector<T>)::const_iterator ci;
      std::vector<int> container(eset.size());
      size_t  arraySize= 0;

      size_t indx = 0;
      for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
      newvec      = ci->second;
      arraySize  += newvec.size();
      container[indx++] = newvec.size();
      }
      }

      dimension =  eset.size();
      vDataspace = H5Screate_simple(rank, &dimension, NULL);
      vDatatype  = H5Tcopy( H5T_NATIVE_INT);
      vDataset   = H5Dcreate(group_id, "ContainerSize", vDatatype, vDataspace, H5P_DEFAULT);
      H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &container[0]);

      H5Dclose( vDataset  );
      H5Sclose( vDataspace);
      H5Tclose( vDatatype );

      //------------------------------------------------------------------------
      // Write (variable) Data into HDF5 format
      //------------------------------------------------------------------------
      std::vector<T>  data(arraySize);

      indx = 0;
      for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
      newvec  = ci->second;
      for( int j = 0; j < newvec.size(); j++) 
      data[indx++] = newvec[j];
      }
      }

      typedef data_schema_traits<T> traits_type;
      DatatypeP dtype = traits_type::get_type();
      vDatatype = dtype->get_hdf5_type();

      rank      = 1;
      dimension = arraySize;

      vDataspace = H5Screate_simple(rank, &dimension, NULL);
      vDataset   = H5Dcreate(group_id, "VariableData", vDatatype,
      vDataspace, H5P_DEFAULT);
      H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

      H5Dclose( vDataset  );
      H5Sclose( vDataspace);
      H5Tclose( vDatatype );
    */
  }
  template <class T> 
  void dmultiStoreRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &eset,hid_t xfer_plist_id) const
  {
    warn(true) ;
  }
  //***************************************************************************/

  template <class T> 
  void dmultiStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &eset)  const
  {   
    warn(true) ;
    /*
      int      rank = 1, stateSize;
      hsize_t  dimension;
      hid_t    vDataspace, vDataset, vDatatype;
      typedef data_schema_traits<T> schema_traits;

      std::vector<T>   newVec;
      entitySet :: const_iterator  ei;
      typename HASH_MAP(int,std::vector<T>)::const_iterator ci;

      //------------------------------------------------------------------------
      // Get the sum of each object size and maximum size of object in the 
      // container for allocation purpose
      //-------------------------------------------------------------------------
      std::vector<int> container(eset.size());
      std::vector<int> bucketSize;
    
      size_t  arraySize= 0;
      int     count;

      size_t indx = 0;
      for( ei = eset.begin(); ei != eset.end(); ++ei) {
      ci = attrib_data.find(*ei);
      if( ci != attrib_data.end() ) {
      newVec  = ci->second;
      container[indx++] = newVec.size();
      for( int ivec = 0; ivec < newVec.size(); ivec++) {
      typename schema_traits::Converter_Type cvtr(newVec[ivec] );
      bucketSize.push_back(cvtr.getSize());
      arraySize  += cvtr.getSize();
      }
      }
      }

      dimension  = eset.size();
      vDataspace = H5Screate_simple(rank, &dimension, NULL);
      vDatatype  = H5T_NATIVE_INT;
      vDataset   = H5Dcreate(group_id, "ContainerSize", vDatatype, vDataspace, H5P_DEFAULT);
      H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &container[0]);

      H5Dclose( vDataset  );
      H5Sclose( vDataspace);

      dimension  = bucketSize.size();
      vDataspace = H5Screate_simple(rank, &dimension, NULL);
      vDatatype  = H5T_NATIVE_INT;
      vDataset   = H5Dcreate(group_id, "SubContainerSize", vDatatype, vDataspace, H5P_DEFAULT);
      H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bucketSize[0]);

      H5Dclose( vDataset  );
      H5Sclose( vDataspace);

      //-------------------------------------------------------------------------
      // Collect state data from each object and put into 1D array
      //-------------------------------------------------------------------------
      typedef typename schema_traits::Converter_Base_Type dtype;

      std::vector<dtype> data(arraySize);
      indx = 0;
      for( ei = eset.begin(); ei != eset.end(); ++ei) {      ci = attrib_data.find(*ei);
      if( ci == attrib_data.end() ) continue;
      newVec =  ci->second;
      for( int ivec = 0; ivec < newVec.size(); ivec++) {
      typename schema_traits::Converter_Type cvtr( newVec[ivec] );
      cvtr.getState( &data[0]+indx, stateSize);
      indx +=stateSize ;
      }
      }

      //-------------------------------------------------------------------------
      // Write (variable) Data into HDF5 format
      //-------------------------------------------------------------------------
      typedef data_schema_traits<dtype> traits_type;

      DatatypeP atom_type = traits_type::get_type() ;
      vDatatype = atom_type->get_hdf5_type();

      dimension =  arraySize;
      vDataspace = H5Screate_simple(rank, &dimension, NULL);
      vDataset   = H5Dcreate(group_id, "VariableData", vDatatype, vDataspace, 
      H5P_DEFAULT);
      H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

      //-----------------------------------------------------------------------
      // Clean up
      //-----------------------------------------------------------------------
      H5Dclose( vDataset  );
      H5Sclose( vDataspace);
      H5Tclose( vDatatype );
    */
  }
  template <class T> 
  void dmultiStoreRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, const entitySet &eset, hid_t xfer_plist_id)  const
  {   
    warn(true) ;
  }
  //************************************************************************/

} // end of namespace Loci

#endif
