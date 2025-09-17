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
#ifndef MULTISTORE_IMPL_H
#define MULTISTORE_IMPL_H

#include <multiStore_def.h>
#include <DMultiStore_def.h>

namespace Loci {
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
  template <class T> 
  inline std::ostream & operator<<(std::ostream &s, const multiStore<T> &m)
  { return m.Print(s) ; }

  //************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, multiStore<T> &m)
  { return m.Input(s) ; }
 
  //************************************************************************/
  template<class T> 
  void multiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }
  //************************************************************************/
  template<class T> 
  store_instance::instance_type
  const_multiStore<T>::access() const
  { return READ_ONLY ; }

  //*************************************************************************/
  template<class T> 
  void const_multiStore<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  //*************************************************************************/
  template<class T> 
  void multiStoreRepI<T>::allocate(const store<int> &sizes) 
  {
    if(alloc_id < 0)
      alloc_id = getStoreAllocateID() ;
    //-------------------------------------------------------------------------
    // Objective: Allocate memeory for multiStore data. This call reclaims 
    // all previously held memory
    //-------------------------------------------------------------------------
    // Assign new entitySet ...

    entitySet ptn = sizes.domain() ;
    int cntid = sizes.Rep()->get_alloc_id() ;
    storeAllocateData[alloc_id].template release<T>() ;
    storeAllocateData[alloc_id].template
      allocMulti<T>(storeAllocateData[cntid],ptn) ;

    base_ptr = ((T **)storeAllocateData[alloc_id].alloc_ptr2 -
		storeAllocateData[alloc_id].base_offset) ;
   
    store_domain = ptn ;

    dispatch_notify();
  }


  //*************************************************************************/
   
  template<class T> 
  void multiStoreRepI<T>::setSizes(const const_multiMap &mm)
  {
    //------------------------------------------------------------------------
    // Objective : Set the degree of each entity specified by the map..
    //------------------------------------------------------------------------

    mutex.lock() ;

    store<int> sizes ;
    sizes.allocate(store_domain) ;
    FORALL(store_domain,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    entitySet map_set = mm.domain() & store_domain ;
    FORALL(map_set,i) {
      sizes[i] = (mm.end(i) - mm.begin(i)) ;
    } ENDFORALL ;
    allocate(sizes) ;
    
    mutex.unlock() ;
  }

  //*************************************************************************/
  
  template<class T> 
  void multiStoreRepI<T>::allocate(const entitySet &ptn) 
  {
    
    if(alloc_id < 0)
      alloc_id = getStoreAllocateID() ;
    //------------------------------------------------------------------------
    // Objective : allocate memory specified by the entitySet. Allocation
    // doesn't resize the memory, therefore reclaims previously held memory.
    //------------------------------------------------------------------------
    if(ptn == store_domain)
      return ;

    
    //-------------------------------------------------------------------------
    // Initialize degree of each entity to zero 
    //-------------------------------------------------------------------------

    store<int> count ;
    count.allocate(ptn) ;
    
    FORALL(ptn,i) {
      count[i] = 0 ;
    } ENDFORALL ;
    
    allocate(count) ;

  }

  template<class T>
  void multiStoreRepI<T>::shift(int_type offset) {
    if(alloc_id < 0)
      alloc_id = getStoreAllocateID() ;
    
    store_domain >>= offset ;
    storeAllocateData[alloc_id].allocset >>= offset ;
    storeAllocateData[alloc_id].base_offset += offset ;
    base_ptr = get_base_ptr() ;
    dispatch_notify() ;
  }
  //*************************************************************************/

  template<class T> 
  multiStoreRepI<T>::~multiStoreRepI() 
  {
    if(alloc_id>=0) {
      storeAllocateData[alloc_id].template release<T>() ;
      releaseStoreAllocateID(alloc_id) ;
      alloc_id = -1 ;
    }
  }

  //*************************************************************************/

  template<class T> 
  storeRep *multiStoreRepI<T>::new_store(const entitySet &p) const 
  {
    store<int> count ;
    count.allocate(p) ;
    entitySet ent = p - domain() ;
    
    for(entitySet::const_iterator ei = p.begin(); ei != p.end(); ++ei)
      count[*ei] = base_ptr[*ei+1] - base_ptr[*ei] ;
    
    for(entitySet::const_iterator ei = ent.begin(); ei != ent.end(); ++ei)
      count[*ei] = 0 ;
    
    return new multiStoreRepI<T>(count) ;
  }
  template<class T> 
  storeRep *multiStoreRepI<T>::new_store(const entitySet &p, const int* cnt) const {
    store<int> count ;
    count.allocate(p) ;
    int t= 0 ;
    FORALL(p, pi) {
      count[pi] = cnt[t++] ; 
    } ENDFORALL ;
    return new multiStoreRepI<T>(count) ;
  }
  //*************************************************************************/

  template<class T> 
  storeRepP multiStoreRepI<T>::remap(const dMap &m) const 
  {

    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    multiStore<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  template<class T>
  storeRepP multiStoreRepI<T>::freeze() {
    return getRep() ;
  }

  template<class T>
  storeRepP multiStoreRepI<T>::thaw() {
    dmultiStore<T> ds ;
    for(entitySet::const_iterator ei=store_domain.begin();
        ei!=store_domain.end();++ei) {
      std::vector<T> v ;
      ds[*ei] = v ;
      int t = end(*ei) - begin(*ei) ;
      for(int i=0;i<t;++i)
        ds[*ei].push_back(base_ptr[*ei][i]) ;
    }
    return ds.Rep() ;
  }

  //*************************************************************************/
  
  template<class T> 
  void multiStoreRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
    const_multiStore<T> s(st) ;
    fatal(base_ptr == 0) ;
    fatal((context - domain()) != EMPTY) ;
    fatal((context - s.domain()) != EMPTY) ;
    store<int> count ;
    count.allocate(domain()) ;
    
    FORALL(domain() - context, i) {
      count[i] = base_ptr[i+1] - base_ptr[i] ;
    } ENDFORALL ;
    
    FORALL(context, i) {
      count[i] = s.end(i) - s.begin(i) ;
    } ENDFORALL ;
    
    // Allocate a temporary copy for the copy
    int cntid = count.Rep()->get_alloc_id() ;
    storeAllocateInfo tmp ;
    tmp.template allocMulti<T>(storeAllocateData[cntid],count.domain()) ;

    T **new_base_ptr = ((T **)tmp.alloc_ptr2 -
			     tmp.base_offset) ;

    // Copy data to new allocation
    FORALL(domain()-context,i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;

    FORALL(context,i) {
      for(int j=0;j<count[i];++j)
        new_base_ptr[i][j] = s[i][j] ;
    } ENDFORALL ;
    // release old allocation
    storeAllocateData[alloc_id].template release<T>() ;
    // copy newly copied temporary to main allocatin
    storeAllocateData[alloc_id] = tmp ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  //*************************************************************************/
  
  template<class T> 
  void multiStoreRepI<T>::gather(const dMap &m, storeRepP &st,
                                 const entitySet &context) 
  {
    store<int> count ;
    const_multiStore<T> s(st) ;
    count.allocate(domain()) ;
    
    FORALL(domain()-context,i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    
    FORALL(context,i) {
      count[i] = s.end(m[i])-s.begin(m[i]) ;
    } ENDFORALL ;
    
    // Allocate a temporary copy for the copy
    int cntid = count.Rep()->get_alloc_id() ;
    storeAllocateInfo tmp ;
    tmp.template allocMulti<T>(storeAllocateData[cntid],count.domain()) ;

    T **new_base_ptr = ((T **)tmp.alloc_ptr2 -
			     tmp.base_offset) ;
    FORALL(domain()-context,i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;

    FORALL(context,i) {
      for(int j=0;j<count[i];++j)
        new_base_ptr[i][j] = s[m[i]][j] ;
    } ENDFORALL ;

    // release old allocation
    storeAllocateData[alloc_id].template release<T>() ;
    // copy newly copied temporary to main allocatin
    storeAllocateData[alloc_id] = tmp ;

    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  //*************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::scatter(const dMap &m, storeRepP &st,
                                  const entitySet &context) 
  {
    store<int> count;
    
    const_multiStore<T> s(st) ;
    count.allocate(domain());

    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
    
    FORALL(domain()-m.image(context),i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    
    FORALL(context,i) {
      count[m[i]] = s.end(i)-s.begin(i) ;
    } ENDFORALL ;
    
    // Allocate a temporary copy for the copy
    int cntid = count.Rep()->get_alloc_id() ;
    storeAllocateInfo tmp ;
    tmp.template allocMulti<T>(storeAllocateData[cntid],count.domain()) ;

    T **new_base_ptr = ((T **)tmp.alloc_ptr2 -
			     tmp.base_offset) ;

    FORALL(domain()-m.image(context),i) {
      for(int j=0;j<count[i];++j) 
        new_base_ptr[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;
    FORALL(context,i) {
      for(int j=0;j<count[m[i]];++j) {
        new_base_ptr[m[i]][j] = s[i][j] ;
      }
    } ENDFORALL ;
    // release old allocation
    storeAllocateData[alloc_id].template release<T>() ;
    // copy newly copied temporary to main allocatin
    storeAllocateData[alloc_id] = tmp ;

    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  //*************************************************************************/

  template<class T> 
  store_type multiStoreRepI<T>::RepType() const 
  {
    return STORE ;
  }

  //*************************************************************************/
  
  template<class T> 
  entitySet multiStoreRepI<T>::domain() const 
  {
    return store_domain ;
  }

  //*************************************************************************/
  
  template<class T> 
  std::ostream &multiStoreRepI<T>::Print(std::ostream &s) const 
  {
    //-------------------------------------------------------------------------
    // Objective : Print the multiStore data in the output stream.
    //-------------------------------------------------------------------------

    s << '{' << domain() << endl ;

    //-------------------------------------------------------------------------
    // Write the size of each entity
    //-------------------------------------------------------------------------

    FORALL(domain(),ii) {
      s << end(ii)-begin(ii) << std::endl ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // Write the data of each entity in the domain 
    //-------------------------------------------------------------------------

    FORALL(domain(),ii) {
      Loci::streamoutput(begin(ii),end(ii)-begin(ii),s) ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // close the bracket ...
    //-------------------------------------------------------------------------

    s << '}' << std::endl ;

    return s ;
  }

  //*************************************************************************/

  template<class T> 
  std::istream &multiStoreRepI<T>::Input(std::istream &s) 
  {
    //-------------------------------------------------------------------------
    // Objective : Read the multiStore data from the input stream
    //-------------------------------------------------------------------------

    entitySet e ;
    char ch ;

    // Look for opening bracket ....
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    //-------------------------------------------------------------------------
    // Read the entitySet intervals ....
    //-------------------------------------------------------------------------
    
    s >> e ;

    //-------------------------------------------------------------------------
    // Read the size of each entity in the set ...
    //-------------------------------------------------------------------------
    
    store<int> sizes ;
    sizes.allocate(e) ;

    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // read the data
    //-------------------------------------------------------------------------
    
    allocate(sizes) ;
    FORALL(e,ii) {
      for(T *ip = begin(ii);ip!=end(ii);++ip) 
        *ip = T() ;
      Loci::streaminput(begin(ii),end(ii)-begin(ii),s) ;
    } ENDFORALL ;

    //-------------------------------------------------------------------------
    // Look for closing brackets
    //-------------------------------------------------------------------------
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  //**************************************************************************/

  template <class T> 
  inline int multiStoreRepI<T>::get_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset ) 
  {

    int packsize = 0 ;
    FORALL(eset,i) {
      int sze = end(i)- begin(i) ;
      packsize += cpypacksize(&sze,1) ;
      packsize += cpypacksize(begin(i),sze) ;
    } ENDFORALL ;

    return packsize ;
  }
  //**************************************************************************/

  template <class T> 
  inline int multiStoreRepI<T>::get_estimated_mpi_size(IDENTITY_CONVERTER c, const entitySet &eset ) 
  {
    int sze;
    sze = 4*eset.size()*sizeof(T) ;
    sze += eset.size()*sizeof(int) ;
    return sze ;
  }
  //**************************************************************************/
  template <class T> 
  int multiStoreRepI<T>::get_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset ) 
  {

    entitySet  :: const_iterator ci;
    typedef data_schema_traits<T> schema_traits;

    entitySet sdom = eset & domain() ;

    int packsize = 0 ;
    for( ci = sdom.begin(); ci != sdom.end(); ++ci) {
      int vsize = end(*ci) - begin(*ci) ;
      for( int ivec = 0; ivec < vsize; ivec++){
        typename schema_traits::Converter_Type cvtr(base_ptr[*ci][ivec] );
        int arraySize = cvtr.getSize() ;
        typename schema_traits::Converter_Base_type *p=0 ;
        packsize += cpypacksize(&arraySize,1) ;
        packsize += cpypacksize(p,arraySize) ;
      }
    }

    return packsize ;
  }
  
  //**************************************************************************/
  template <class T> 
  int multiStoreRepI<T>::get_estimated_mpi_size(USER_DEFINED_CONVERTER c, const entitySet &eset ) 
  {
    int  vsize = 4*eset.size()*50*sizeof(double)+
      4*eset.size()*sizeof(int)  +   // size of each object
      eset.size()*sizeof(int);       // size of each entity
    
    return(vsize);
  }
  //**************************************************************************/
  template <class T> 
  int multiStoreRepI<T>::pack_size(const entitySet &eset ) 
  {
    fatal((eset - domain()) != EMPTY);
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    warn(eset-domain() != EMPTY) ;
    return get_mpi_size( traits_type, eset );
  }
  
  //**************************************************************************/
  template <class T> 
  int multiStoreRepI<T>::estimated_pack_size(const entitySet &eset ) 
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
    return get_estimated_mpi_size( traits_type, eset );  
  }
  


  
  template<class T> int multiStoreRepI<T>::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;

    typedef typename data_schema_traits<T>::Schema_Converter
      schema_converter;
    schema_converter traits_type;

    return get_mpi_size(traits_type, packed);    
  }
  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf,
                                    int &position,  int outcount,
                                    const entitySet &eset ) 
  {
    FORALL(eset,i) {
      int vsize   = end(i)-begin(i);
      cpypack(outbuf,position,outcount,&base_ptr[i][0],vsize) ;
    } ENDFORALL ;
  }
  
  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf, 
                                    int &position, int outcount,
                                    const entitySet &eset ) 
  {

    entitySet::const_iterator ci;

    //-------------------------------------------------------------------------
    // Get the maximum size of container
    //-------------------------------------------------------------------------
    int vecsize, stateSize, maxStateSize=0;

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vecsize = end(*ci)-begin(*ci);
      for( int ivec = 0; ivec < vecsize; ivec++){
        typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*ci][ivec]);
        stateSize = cvtr.getSize();
        maxStateSize = max( maxStateSize, stateSize);
      }
    }

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    std::vector<dtype> inbuf(maxStateSize);

    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vecsize = end(*ci)-begin(*ci);
      for( int ivec = 0; ivec < vecsize; ivec++){
        typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*ci][ivec]);
        cvtr.getState( &inbuf[0], stateSize);

        cpypack(outbuf,position,outcount,&stateSize,1) ;
        cpypack(outbuf,position,outcount,&inbuf[0],stateSize) ;
      }
    }
  }

  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::pack( void *outbuf, int &position, int &outcount, 
                                const entitySet &eset ) 
  {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;

    FORALL(eset,ii){
      int size = end(ii) - begin(ii) ;
      cpypack(outbuf,position,outcount,&size,1) ;
    }ENDFORALL ;

    packdata( traits_type, outbuf, position, outcount, eset);
  }

  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::unpack(void *ptr, int &loc, int &size, 
                                 const sequence &seq) 
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_type;
    if(base_ptr == 0) return ;

    entitySet new_dom = domain() | entitySet(seq) ;
    entitySet ent = domain() - entitySet(seq);

    store<int> ecount ;
    ecount.allocate(new_dom) ;

    for(Loci::sequence::const_iterator si = seq.begin(); si != seq.end(); ++si) {
      cpyunpack(ptr,loc,size,&ecount[*si],1) ;
    }

    unpackdata( traits_type, ptr, loc, size, seq); 
  }

  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::unpackdata( IDENTITY_CONVERTER c, void *inbuf, int &position,
                                      int insize, const sequence &seq) 
  {
    sequence :: const_iterator si;

    for(si = seq.begin(); si != seq.end(); ++si) {
      int vsize    = end(*si)-begin(*si);
      cpyunpack(inbuf,position,insize,&base_ptr[*si][0],vsize) ;
    }
  }
  
  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf, 
                                      int &position, int insize,
                                      const sequence &seq) 
  {
    sequence :: const_iterator ci;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;


    std::vector<dtype> outbuf;
    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      int vsize  = end(*ci)-begin(*ci);
      for( int ivec = 0; ivec < vsize; ivec++) {
        int stateSize = 0 ;
        cpyunpack(inbuf,position,insize,&stateSize,1) ;
        if( stateSize > outbuf.size() ) outbuf.resize(stateSize);
        cpyunpack(inbuf,position,insize,&outbuf[0],stateSize) ;

        typename schema_traits::Converter_Type  cvtr( base_ptr[*ci][ivec] );
        cvtr.setState( &outbuf[0], stateSize);
      }
    }
  }
  
  template<class T> 
  frame_info multiStoreRepI<T>::get_frame_info() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }
  template<class T> 
  frame_info multiStoreRepI<T>::get_frame_info(IDENTITY_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = 0 ;
    entitySet dom = domain() ;
    int newsize = 0 ;
    for(entitySet::const_iterator ci = dom.begin(); ci != dom.end(); ++ci) {
      newsize = end(*ci) - begin(*ci);
      fi.first_level.push_back(newsize) ;
    }
    return fi ;
  }
  template<class T> 
  frame_info multiStoreRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g) {
    int vsize ;
    entitySet dom = domain() ;
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = 0 ;
    std::vector<int> fl(dom.size()) ;
    int stateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    for(entitySet::const_iterator ci = dom.begin(); ci != dom.end(); ++ci) {
      vsize = end(*ci) - begin(*ci) ;
      fi.first_level.push_back(vsize) ;
      for(int ivec = 0; ivec < vsize; ivec++){
        typename schema_traits::Converter_Type cvtr(base_ptr[(*ci)][ivec]) ;
        stateSize = cvtr.getSize();
        fi.second_level.push_back(stateSize) ;
      }
    }
    return fi ;
  }
  template<class T> 
  DatatypeP multiStoreRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T> 
  DatatypeP multiStoreRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T> 
  DatatypeP multiStoreRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  //**************************************************************************/
  
  template<class T> 
  void multiStoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset) 
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, user_eset) ;
  }

#ifdef H5_HAVE_PARALLEL
  template<class T> 
  void multiStoreRepI<T>::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &user_eset, hid_t xfer_plist_id) 
  {
    warn(true);
  }
#endif
  //**************************************************************************/

  template <class T> 
  void multiStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset)
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                          H5P_DEFAULT, tmp_array) ;
      if(err < 0) {
        cerr << "H5Dread() failed" << endl ;
      }
      int loc = 0 ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) 
        for(int ivec = 0; ivec < (fi.first_level)[loc]; ivec++) {
          base_ptr[*si][ivec] = tmp_array[tmp++] ;
          loc++ ;
        }
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

#ifdef H5_HAVE_PARALLEL 
  template <class T> 
  void multiStoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset, hid_t cfer_plist_id)
  {
    warn(true);
  }
#endif  
  //**************************************************************************/
  
  template <class T> 
  void multiStoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER, frame_info &fi, entitySet &eset)
  {
    typedef data_schema_traits<T> schema_traits ;
    if(dimension != 0) {
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
        cerr << "H5Dread() failed" << endl ;
      }
      size_t tmp = 0 ;
      int bucsize ;
      size_t indx = 0 ;
      int loc = 0 ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) {
        for(int ivec = 0; ivec < (fi.first_level)[loc]; ivec++) {
          typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*si][ivec]);
          bucsize = vint[indx++] ;
          cvtr.setState(tmp_array+tmp, bucsize) ;
          tmp += bucsize ;
        }
        loc++ ;
      }
      
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }  

#ifdef H5_HAVE_PARALLEL 
  template <class T> 
  void multiStoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    warn(true);
  }
#endif
  //**************************************************************************/

  template<class T> 
  void multiStoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &usr_eset) const 
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, usr_eset) ;
  }
  
#ifdef H5_HAVE_PARALLEL
  template<class T> 
  void multiStoreRepI<T>::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &usr_eset, hid_t xfer_plist_id) const 
  {
    warn(true);
  }
#endif
  
  //**************************************************************************/
  template <class T> 
  void multiStoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &eset) const
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      int newsize = 0 ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
        newsize = end(*si) - begin(*si);
        for(int ivec = 0; ivec < newsize; ivec++){
          tmp_array[tmp++] = base_ptr[*si][ivec] ;
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
  void multiStoreRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &eset, hid_t xfer_plist_id) const
  {
    warn(true);
  }
#endif
  //*************************************************************************/

  template <class T> 
  void multiStoreRepI<T>:: hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset)  const
  {   
    typedef data_schema_traits<T> schema_traits ;
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      typedef typename schema_traits::Converter_Base_Type dtype;
      dtype* tmp_array = new dtype[dimension] ;
      size_t tmp = 0 ;
      int stateSize = 0 ;
      int newsize = 0 ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
        newsize = end(*si) - begin(*si);
        for(int ivec = 0; ivec < newsize; ivec++){
          typename schema_traits::Converter_Type cvtr(base_ptr[*si][ivec]);
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
  template <class T> 
  void multiStoreRepI<T>:: hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset, hid_t xfer_plist_id)  const
  {
#ifndef H5_HAVE_PARALLEL    
    hdf5write(group_id,  dataspace, dataset,  dimension,  name,  g, eset);
#else
    warn(true);
#endif
  }

  //*************************************************************************/
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
} // end of namespace Loci

#endif
