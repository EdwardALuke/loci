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
#ifndef STOREVEC_IMPL_H
#define STOREVEC_IMPL_H

#include <storeVec_def.h>
#include <DStoreVec_def.h>

namespace Loci {
  template<class T> 
  std::ostream &storeVecRepI<T>::Print(std::ostream &s) const
  {
    
    s << '{' << domain() << std::endl ;
    s << size << std::endl ;

    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    FORALL(domain(),ii) {
      T * p = alloc_ptr + (ii-base_offset)*size ;
      Loci::streamoutput(p,size,s) ;
    }ENDFORALL ;
    s << '}' << std::endl ;

    return s ;
  }

  //******************************************************************/

  template<class T> 
  std::istream &storeVecRepI<T>::Input(std::istream &s)
  {

    //------------------------------------------------------------------
    // Objective : Read the storeVec from the input stream.
    //------------------------------------------------------------------
    char ch ;
    
    // Look for the opening brackets ...
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    entitySet e ;
    int sz ;

    s >> e ;               // Read the entitySet intervals.
    s >> sz ;              // Read the size of the vector.
    
    set_elem_size(sz) ;
    allocate(e) ;

    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    FORALL(e,ii) {
      T * p = alloc_ptr + (ii-base_offset)*size ;
      if(!std::is_trivially_default_constructible<T>::value) 
	for(int i=0;i<size;++i)
	  p[i] = T() ;
      
      Loci::streaminput(p,size,s) ;
    } ENDFORALL ;
    
    // Look for the closing brackets ...
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
	  
    return s ;
  }

  //*****************************************************************/
  template<class T> 
  void storeVecRepI<T>::allocate(const entitySet &ptn) {
    if(alloc_id < 0)
      alloc_id = getStoreAllocateID() ;


    if(size != 0) {
      storeAllocateData[alloc_id].template allocBasic<T>(ptn,size) ;
    } else {
      storeAllocateData[alloc_id].template release<T>() ;
    }
    
    store_domain = ptn ;
    
    dispatch_notify() ;
  }
  
  
  template<class T>
  void storeVecRepI<T>::shift(int_type offset) {
    store_domain >>= offset ;
    storeAllocateData[alloc_id].allocset >>=offset ;
    storeAllocateData[alloc_id].base_offset += offset ;
    dispatch_notify() ;
  }

  //*******************************************************************/

  template<class T> 
  storeVecRepI<T>::~storeVecRepI() 
  {
    if(alloc_id >= 0) {
      storeAllocateData[alloc_id].template release<T>() ;
      releaseStoreAllocateID(alloc_id) ;
      alloc_id = -1 ;
    }

  }

  //*******************************************************************/

  template<class T>
  storeRep *storeVecRepI<T>::new_store(const entitySet &p) const 
  {
    storeRep *sp = new storeVecRepI<T>(p)  ;
    sp->set_elem_size(size) ;
    return sp ;
  }
  
  template<class T>
  storeRep *storeVecRepI<T>::new_store(const entitySet &p, const int* cnt) const 
  {
    storeRep* sp = 0;
    cerr << " This method should not be called for a storeVec " << endl ;
    return sp ;
  }

  //*******************************************************************/

  template<class T> 
  store_type storeVecRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //******************************************************************/

  template<class T> 
  entitySet storeVecRepI<T>::domain() const 
  {
    return store_domain ;
  }

  //*****************************************************************/

  template<class T> 
  void storeVecRepI<T>::set_elem_size(int sz) 
  {
    mutex.lock() ;

    //----------------------------------------------------------------
    // Change the size of vector held. It will reclaim the memory used before
    // his call. and allocate new one.
    //----------------------------------------------------------------
    if(alloc_id < 0) 
      alloc_id = getStoreAllocateID() ;
      
    storeAllocateData[alloc_id].size = size ;
    
    if(size != sz) {
      size = sz ;
      storeAllocateData[alloc_id].size = size ;
      allocate(store_domain) ;
    }
    mutex.unlock() ;
  }
  
  //*******************************************************************/
  template<class T> 
  void storeVec<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;

    if(p!=0) {
      alloc_ptr = p->get_alloc_ptr() ;
      base_offset = p->get_base_offset() ;
      size = p->get_size() ;
    }

    warn(p == 0) ;
  }

  //******************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const storeVec<T> &t)
  { return t.Print(s) ; }

  //*******************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, storeVec<T> &t)
  { return t.Input(s) ; }

  //*******************************************************************/
  template<class T> 
  void const_storeVec<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      alloc_ptr = p->get_alloc_ptr() ;
      base_offset = p->get_base_offset() ;
      size = p->get_size() ;
    }
    warn(p == 0) ;
  }

  //******************************************************************/

  template<class T> store_instance::instance_type
  const_storeVec<T>::access() const
  { return READ_ONLY ; }

  //******************************************************************/
  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const const_storeVec<T> &t)
  { return t.Print(s) ; }

  //************************************************************************/
  template <class T> 
  storeRepP storeVecRepI<T>::remap(const dMap &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    storeVec<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  template<class T>
  storeRepP storeVecRepI<T>::freeze() {
    return getRep() ;
  }

  template<class T>
  storeRepP storeVecRepI<T>::thaw() {
    dstoreVec<T> ds ;
    ds.setVecSize(size) ;
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    for(entitySet::const_iterator ei=store_domain.begin();
        ei!=store_domain.end();++ei) {
      T* b = alloc_ptr+( (*ei-base_offset)*size) ;
      T* e = b + size ; 
      ds[*ei] = std::vector<T>(b,e) ;
    }
    return ds.Rep() ;
  }
  
  template<class T>
  storeRepP storeVecRepI<T>::thaw(const entitySet& es) const {
    entitySet shared = domain() & es ;
    entitySet out = es - domain() ;
    
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    dstoreVec<T> ds ;
    ds.setVecSize(size) ;
    for(entitySet::const_iterator ei=shared.begin();
        ei!=shared.end();++ei) {
      T* b = alloc_ptr+( (*ei-base_offset)*size) ;
      T* e = b + size ; 
      ds[*ei] = std::vector<T>(b,e) ;
    }
    
    Entity c = *domain().begin() ;
    T* b = alloc_ptr + ((c-base_offset)*size) ;
    T* e = b + size ;
    for(entitySet::const_iterator ei=out.begin();
        ei!=out.end();++ei)
      ds[*ei] = std::vector<T>(b,e) ;
    
    return ds.Rep() ;
  }
  
  //*************************************************************************/
  template <class T> 
  void storeVecRepI<T>::copy(storeRepP &st, const entitySet &context) {
    const_storeVec<T> s(st) ;
    int sz = s.vecSize() ;
    set_elem_size(sz) ;
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      T *p = alloc_ptr + (i-base_offset)*sz ;
      for(int j=0;j<sz;++j)
        p[j] = s[i][j] ;
    } ENDFORALL ;
  }

  //**************************************************************************/

  template <class T> 
  void storeVecRepI<T>::gather(const dMap &m, storeRepP &st,
                               const entitySet &context) 
  {
    const_storeVec<T> s(st) ;
    int sz = s.vecSize() ;
    set_elem_size(sz) ;
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    fatal(alloc_ptr == 0) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      T *p = alloc_ptr + (i-base_offset)*sz ;
      for(int j=0;j<sz;++j)
        p[j] = s[m[i]][j] ;
    } ENDFORALL ;
  }

  //*************************************************************************/

  template <class T> 
  void storeVecRepI<T>::scatter(const dMap &m, storeRepP &st,
                                const entitySet &context) 
  {
    const_storeVec<T> s(st) ;
    int sz = s.vecSize() ;
    set_elem_size(sz) ;
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    fatal((context != EMPTY) && (alloc_ptr == 0)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
	
    FORALL(context,i) {
      T *p = alloc_ptr + (m[i]-base_offset)*sz ;
      for(int j=0;j<sz;++j)
        p[j] = s[i][j] ;
    } ENDFORALL ;
  }

  //**************************************************************************/

  template <class T>
  inline int storeVecRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {
    int  M = get_size() ;
    return cpypacksize(&M,1)+cpypacksize(get_alloc_ptr(),M*eset.size()) ;
  }
  //**************************************************************************/

  template <class T>
  inline int storeVecRepI<T>::get_estimated_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {
    int mysize;
  
    if( isMat) mysize = sizeof(T) * eset.size() * 25 +sizeof(int);
    else  mysize = sizeof(T) * eset.size() * 5 +sizeof(int) ;
    return (mysize) ;

  }


  
  //**************************************************************************/

  template <class T>
  int storeVecRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {

    int M = get_size() ;
    entitySet::const_iterator ci;

    fatal((eset - domain()) != EMPTY);

    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    entitySet sdom = eset & domain() ;
    int byteCount = cpypacksize(&M,1) ;
    typedef data_schema_traits<T> converter_traits;
    for( ci = sdom.begin(); ci != sdom.end(); ++ci) {
      for( int ivec = 0; ivec < M; ivec++){
        typename converter_traits::Converter_Type cvtr(alloc_ptr[((*ci)-base_offset)*size+ivec]) ;
        int stateSize = cvtr.getSize() ;
        typename data_schema_traits<T>::Converter_Base_Type *p=0 ;
        byteCount += cpypacksize(&stateSize,1)+cpypacksize(p,stateSize) ;
      }
    }

    return byteCount ;
  }
  

  //**************************************************************************/

  template <class T>
  int storeVecRepI<T>::get_estimated_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {
    

    int     arraySize =0, numContainers = 0;
    int estimated_converter_size = 50*sizeof(double);
   
    
   
   
    if(isMat) arraySize = 25*eset.size()*estimated_converter_size;
    else arraySize = 5*eset.size()*estimated_converter_size;
    
    
    
    if(isMat) numContainers =  25*eset.size();
    else  numContainers =  5*eset.size();
    
    
    return(arraySize +
           (numContainers+1)*sizeof(int));
  }




  
  //**************************************************************************/

  template <class T>
  int storeVecRepI<T>::pack_size( const entitySet &eset)
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size( traits_type, eset );
  }

  template <class T>
  int storeVecRepI<T>::estimated_pack_size( const entitySet &eset)
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_estimated_mpi_size( traits_type, eset );
  }

  template<class T> int storeVecRepI<T>::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    return get_mpi_size(traits_type, packed);    
  }
  //*************************************************************************/

  template <class T>
  inline void storeVecRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf, int &position,
                                         int outcount, const entitySet &eset )
  {

    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    const int M = get_size() ;
    for(size_t i = 0; i < eset.num_intervals(); ++i) {
      Loci::int_type indx1 = eset[i].first ;
      Loci::int_type stop  = eset[i].second ;
      T *p = alloc_ptr + M * (indx1 -base_offset);
      const int t = (stop - indx1 + 1) * M ;

      cpypack(outbuf,position,outcount,p,t) ;
    }
  }

  //**************************************************************************/

  template <class T> 
  inline void storeVecRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                                         int &position, int outcount, 
                                         const entitySet &eset ) 
  {
    entitySet::const_iterator ci;

    //-------------------------------------------------------------------------
    // Get the maximum size of container 
    //-------------------------------------------------------------------------
    int stateSize, maxStateSize=0;
    
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++){
        typename data_schema_traits<T>::Converter_Type
          cvtr( alloc_ptr[((*ci)-base_offset)*size+ivec] );
        stateSize = cvtr.getSize();
        maxStateSize = max( maxStateSize, stateSize);
      }
    }

    typedef data_schema_traits<T> schema_traits; 
    typedef typename schema_traits::Converter_Base_Type dtype;

    std::vector<dtype> inbuf(maxStateSize);

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++){
        typename data_schema_traits<T>::Converter_Type
          cvtr( alloc_ptr[((*ci)-base_offset)*size+ivec] );
        cvtr.getState( &inbuf[0], stateSize);

        cpypack(outbuf,position,outcount,&stateSize,1) ;
        cpypack(outbuf,position,outcount,&inbuf[0],stateSize) ;
      }
    }
  }

  //**************************************************************************/

  template <class T>
  void storeVecRepI<T>::pack(void *outbuf, int &position, int &outcount,
                             const entitySet &eset )
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    int M = get_size() ;

    cpypack(outbuf,position,outcount,&M,1) ;
    
    packdata( traits_type, outbuf, position, outcount, eset);
  }

  //**************************************************************************/
  template <class T> 
  inline void storeVecRepI<T>::unpackdata( IDENTITY_CONVERTER c,
                                           void *inbuf,
                                           int &position,
                                           int &insize,
                                           const sequence &seq)
  {
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    const int M = get_size() ;
    for(size_t i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
        const Loci::int_type indx1 = seq[i].first ;
        const Loci::int_type stop =  seq[i].second ;
        for(Loci::int_type indx = indx1; indx != stop-1; --indx) {
          T *p = alloc_ptr + M * (indx-base_offset) ;

          cpyunpack(inbuf,position,insize,p,M) ;
        }
      } else {
        Loci::int_type indx1 = seq[i].first ;
        Loci::int_type stop = seq[i].second ;
        T *p = alloc_ptr + M * (indx1-base_offset) ;
        const int t = (stop - indx1 + 1) * M ;

        cpyunpack(inbuf,position,insize,p,t) ;
      }
    }

  }

  //***********************************************************************/
  template <class T> 
  void storeVecRepI<T>::unpackdata( USER_DEFINED_CONVERTER c,
                                    void *inbuf, 
                                    int &position,
                                    int &insize,
                                    const sequence &seq)
  {

    sequence :: const_iterator ci;
    int  stateSize ;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    std::vector<dtype> outbuf;

    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      for( int ivec = 0; ivec < size; ivec++) {
        cpyunpack(inbuf,position,insize,&stateSize,1) ;

        if( size_t(stateSize) > outbuf.size() ) outbuf.resize(stateSize);

        cpyunpack(inbuf,position,insize,&outbuf[0],stateSize) ;

        typename schema_traits::Converter_Type  cvtr( alloc_ptr[((*ci)-base_offset)*size+ivec] );
        cvtr.setState( &outbuf[0], stateSize);
      }
    }
  }

  //**************************************************************************/
  template <class T> 
  void storeVecRepI<T>::unpack(void *inbuf, int &position, int &insize, const sequence &seq)
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

#ifdef DEBUG
    entitySet ecommon, ediff,eset(seq);

    ediff = eset - domain();
    if( ediff.size() > 0) { 
      std::cout << "Error:Entities not part of domain " << ediff <<endl;
      abort();
    }
#endif
    
    int init_size = get_size() ;
    int M ;
    cpyunpack(inbuf,position,insize,&M,1) ;
    warn(M == 0);
    fatal(M<0) ;
    if(M > init_size) {
      set_elem_size(M) ;
    }
    unpackdata( traits_type, inbuf, position, insize, seq);
  }
  
  template<class T> 
  frame_info storeVecRepI<T>::get_frame_info() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }
  template<class T> 
  frame_info storeVecRepI<T>::get_frame_info(IDENTITY_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = get_size() ;
    return fi ;
  }
  
  template<class T> 
  frame_info storeVecRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g) {
    entitySet dom = domain() ;
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = get_size() ;
    int stateSize = 0;
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    typedef data_schema_traits<T> schema_traits ;
    for(entitySet::const_iterator ci = dom.begin(); ci != dom.end(); ++ci) 
      for(int ivec = 0; ivec < fi.size; ivec++){
        typename schema_traits::Converter_Type cvtr(alloc_ptr[((*ci)-base_offset)*fi.size+ivec] );
        stateSize = cvtr.getSize();
        fi.second_level.push_back(stateSize) ;
      }
    return fi ;
  }
  
  template<class T> 
  DatatypeP storeVecRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T> 
  DatatypeP storeVecRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T> 
  DatatypeP storeVecRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  //*******************************************************************/

  template<class T> 
  void storeVecRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &usr_eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, usr_eset) ;
  }
#ifdef H5_HAVE_PARALLEL
  template<class T> 
  void storeVecRepI<T>::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &usr_eset, hid_t xfer_plist_id) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5writeP(group_id, dataspace, dataset, dimension, name, traits_output_type, usr_eset, xfer_plist_id) ;
  }
#endif  
  //************************************************************************/
  
  template <class T>  
  void storeVecRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &eset) const
  {
    if(dimension != 0) {
      T *alloc_ptr = get_alloc_ptr() ;
      int base_offset = get_base_offset() ;
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      int qs = get_size() ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
        for(int ivec = 0; ivec < qs; ivec++){
          tmp_array[tmp++] = alloc_ptr[((*si)-base_offset)*qs+ivec] ;
        }
      }
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }
#ifdef H5_HAVE_PARALLEL 
  template <class T>  
  void storeVecRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &eset, hid_t xfer_plist_id) const
  {
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    T* tmp_array = new T[dimension] ;
    size_t tmp = 0 ;
    int qs = get_size() ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
      for(int ivec = 0; ivec < qs; ivec++){
        tmp_array[tmp++] = alloc_ptr[((*si)-base_offset)*qs+ivec] ;
      }
    }
    H5Dwrite(dataset, datatype, memspace, dataspace, xfer_plist_id, tmp_array) ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
  }
#endif
  //*************************************************************************/
  
  template <class T>  
  void storeVecRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset) const
  { 
    typedef data_schema_traits<T> schema_traits ;
    if(dimension != 0) {
      T *alloc_ptr = get_alloc_ptr() ;
      int base_offset = get_base_offset() ;
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      typedef typename schema_traits::Converter_Base_Type dtype;
      dtype* tmp_array = new dtype[dimension] ;
      size_t tmp = 0 ;
      int stateSize = 0 ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
        for(int ivec = 0; ivec < size; ivec++){
          typename schema_traits::Converter_Type cvtr(alloc_ptr[((*si)-base_offset)*size+ivec]);
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
#ifdef H5_HAVE_PARALLEL
  template <class T>  
  void storeVecRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset, hid_t xfer_plist_id) const
  {
    typedef data_schema_traits<T> schema_traits ;
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    dtype* tmp_array = new dtype[dimension] ;
    size_t tmp = 0 ;
    int stateSize = 0 ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
      for(int ivec = 0; ivec < size; ivec++){
        typename schema_traits::Converter_Type cvtr(alloc_ptr[((*si)-base_offset)*size+ivec]);
        cvtr.getState(tmp_array+tmp, stateSize) ;
        tmp +=stateSize ;
      }
    }
    H5Dwrite(dataset, datatype, memspace, dataspace, xfer_plist_id, tmp_array) ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
  }
#endif
  //**************************************************************************/

  template<class T> 
  void storeVecRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T> 
  void storeVecRepI<T>::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5readP(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset, xfer_plist_id) ;
  }
#endif
  //**************************************************************************/
  
  template <class T> 
  void storeVecRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER convert, frame_info &fi, entitySet &eset)
  {
    if(dimension != 0) {
      T *alloc_ptr = get_alloc_ptr() ;
      int base_offset = get_base_offset() ;
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      int qs = fi.size ;
      hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                          H5P_DEFAULT, tmp_array) ;
      if(err < 0) {
        std::string error = "H5Dread: failed" ;
        throw StringError(error) ;
      }

      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) 
        for(int ivec = 0; ivec < qs; ivec++) {
          alloc_ptr[((*si)-base_offset)*qs+ivec] = tmp_array[tmp++] ;
        }
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }
#ifdef H5_HAVE_PARALLEL
  template <class T> 
  void storeVecRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER convert, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    T* tmp_array = new T[dimension] ;
    size_t tmp = 0 ;
    int qs = fi.size ;
    hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                        xfer_plist_id, tmp_array) ;
    if(err < 0) {
      std::string error = "H5Dread: failed" ;
      throw StringError(error) ;
    }

    for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) 
      for(int ivec = 0; ivec < qs; ivec++) {
        alloc_ptr[((*si)-base_offset)*qs+ivec] = tmp_array[tmp++] ;
      }
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
  }
#endif
  //************************************************************************/

  template <class T> 
  void storeVecRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset)
  {
    typedef data_schema_traits<T> schema_traits ;
    if(dimension != 0) {
      T *alloc_ptr = get_alloc_ptr() ;
      int base_offset = get_base_offset() ;
      storeRepP qrep = getRep() ;
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
      int qs = fi.size ;
      size_t tmp = 0 ;
      int bucsize ;
      size_t indx = 0 ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) 
        for(int ivec = 0; ivec < qs; ivec++) {
          typename data_schema_traits<T>::Converter_Type cvtr(alloc_ptr[((*si)-base_offset)*qs+ivec]);
          bucsize = vint[indx++] ;
          cvtr.setState(tmp_array+tmp, bucsize) ;
          tmp += bucsize ;
        }
        
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

#ifdef H5_HAVE_PARALLEL
  template <class T> 
  void storeVecRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    typedef data_schema_traits<T> schema_traits ;
    T *alloc_ptr = get_alloc_ptr() ;
    int base_offset = get_base_offset() ;
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    std::vector<int> vint = fi.second_level ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    dtype* tmp_array = new dtype[dimension] ;    
    hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                        xfer_plist_id, tmp_array) ;
    if(err < 0) {
      std::string error = "H5Dread: failed" ;
      throw StringError(error) ;
    }
    int qs = fi.size ;
    size_t tmp = 0 ;
    int bucsize ;
    size_t indx = 0 ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) 
      for(int ivec = 0; ivec < qs; ivec++) {
        typename data_schema_traits<T>::Converter_Type cvtr(alloc_ptr[((*si)-base_offset)*qs+ivec]);
        bucsize = vint[indx++] ;
        cvtr.setState(tmp_array+tmp, bucsize) ;
        tmp += bucsize ;
      }
        
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
  }
#endif
  //******************************************************************/
  
} // end of namespace Loci

#endif
