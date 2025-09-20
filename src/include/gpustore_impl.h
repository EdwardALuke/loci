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
#ifndef GPUSTORE_IMPL_H
#define GPUSTORE_IMPL_H

#include <store_def.h>
#include <DStore_def.h>

#ifdef USE_CUDA_RT
#include <cuda_runtime_api.h>
#endif

namespace Loci {
  // Code to copy from cpu container to gpu container
  template<class T> 
  void gpustoreRepI<T>::copyFrom(const storeRepP &p, entitySet set)  {
    int setivals = set.num_intervals() ;
    store<T> v(p) ;
    T *gpu_base_ptr = get_base_ptr() ;
    for(int i=0;i<setivals;++i) {
      int start = set[i].first ;
      int end = set[i].second ;
      int sz = end-start+1 ;
#ifdef USE_CUDA_RT
      cudaError_t err = cudaMemcpy(gpu_base_ptr+start,&v[start],sizeof(T)*sz,
			       cudaMemcpyHostToDevice) ;
      if(err!= cudaSuccess) {
	cerr << "cudaMemcpy failed in gpuMapRepI::copyFrom" << endl ;
	Loci::Abort() ;
      }
      cudaDeviceSynchronize() ;
#endif
    }
    //    FORALL(set,ii) {
    //      base_ptr[ii] = v[ii] ;
    //    } ENDFORALL ;
  }
  
  // code to copy from gpu container to cpu container
  template<class T>
  void gpustoreRepI<T>::copyTo(storeRepP &p, entitySet set) const {
    int setivals = set.num_intervals() ;
    store<T> v(p) ;
    const T *gpu_base_ptr = get_base_ptr() ;
    for(int i=0;i<setivals;++i) {
      int start = set[i].first ;
      int end = set[i].second ;
      int sz = end-start+1 ;
#ifdef USE_CUDA_RT
      cudaDeviceSynchronize() ;
      cudaError_t err = cudaMemcpy(&v[start],gpu_base_ptr+start,sizeof(T)*sz,
			       cudaMemcpyDeviceToHost) ;
      if(err!= cudaSuccess) {
	cerr << "cudaMemcpy failed in gpuMapRepI::copyFrom" << endl ;
	Loci::Abort() ;
      }
#endif
    }
  }

  template<class T> void gpustoreRepI<T>::allocate(const entitySet &eset) {
    if(alloc_id < 0) {
      alloc_id = getGPUStoreAllocateID() ;
      GPUstoreAllocateData[alloc_id].template allocBasic<T>(eset,1) ;
    } else if(eset != GPUstoreAllocateData[alloc_id].allocset) {
      GPUstoreAllocateData[alloc_id].template allocBasic<T>(eset,1) ;
    }      
      
    store_domain = GPUstoreAllocateData[alloc_id].allocset ; ;
    dispatch_notify() ;
    return ;
  }

  template<class T> void gpustoreRepI<T>::shift(int_type offset) {
    store_domain >>= offset ;
    GPUstoreAllocateData[alloc_id].base_offset += offset ;
    dispatch_notify() ;
  }

  template<class T> std::ostream &gpustoreRepI<T>::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    T *base_ptr = get_base_ptr() ;
    FORALL(domain(),ii) {
      Loci::streamoutput(&base_ptr[ii],1,s) ;
    }ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }


  template<class T> std::istream &gpustoreRepI<T>::Input(std::istream &s) {
    entitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    allocate(e) ;

    T *base_ptr = get_base_ptr() ;
    FORALL(e,ii) {
      base_ptr[ii] = T() ;
      Loci::streaminput(&base_ptr[ii],1,s) ;
    } ENDFORALL ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  //*************************************************************************/

  template<class T>  gpustoreRepI<T>::~gpustoreRepI() {
    if(alloc_id>=0) {
      GPUstoreAllocateData[alloc_id].template release<T>() ;
      releaseGPUStoreAllocateID(alloc_id) ;
      alloc_id = -1 ;
    }
  }

  template<class T>  entitySet gpustoreRepI<T>::domain() const {
    return store_domain ;
  }

  template<class T>
  storeRep *gpustoreRepI<T>::new_store(const entitySet &p) const {
    return new gpustoreRepI<T>(p)  ;
  }
  template<class T>
  storeRep *gpustoreRepI<T>::new_store(const entitySet &p, const int* cnt) const {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a store " << endl ;
    return sp ;
  }
  template<class T> store_type gpustoreRepI<T>::RepType() const {
    return GPUSTORE ;
  }

  template<class T> void gpustore<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  template<class T> inline std::ostream & operator<<(std::ostream &s, const gpustore<T> &t)
  { return t.Print(s) ; }

  template<class T> inline std::istream & operator>>(std::istream &s, gpustore<T> &t)
  { return t.Input(s) ; }

  template<class T> void const_gpustore<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p != 0)
      base_ptr = p->get_base_ptr() ;
    warn(p == 0) ;
  }

  template<class T> store_instance::instance_type
  const_gpustore<T>::access() const { return READ_ONLY; }

  template<class T>
  storeRepP gpustoreRepI<T>::remap(const dMap &m) const {
    entitySet newdomain = m.domain() & domain() ;
    entitySet mapimage = m.image(newdomain) ;
    gpustore<T> s ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;

    s.Rep()->scatter(m,my_store,newdomain) ;
    return s.Rep() ;
  }

  template<class T> storeRepP gpustoreRepI<T>::freeze() {
    return getRep() ;
  }

  template<class T> storeRepP gpustoreRepI<T>::thaw() {
    dstore<T> ds ;
    T *base_ptr = get_base_ptr() ;
    for(entitySet::const_iterator ei=store_domain.begin();
        ei!=store_domain.end();++ei)
      ds[*ei] = base_ptr[*ei] ;
    return ds.Rep() ;
  }

  
  template<class T> void gpustoreRepI<T>::copy(storeRepP &st, const entitySet &context)  {
    const_gpustore<T> s(st) ;
    T *base_ptr = get_base_ptr() ;
    fatal((context != EMPTY) && (base_ptr ==0)) ;
    fatal((context-domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[i] ;
    } ENDFORALL ;
  }
  
  template<class T> void gpustoreRepI<T>::gather(const dMap &m, storeRepP &st,
                                              const entitySet &context) {
    const_gpustore<T> s(st) ;
    T *base_ptr = get_base_ptr() ;
    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[m[i]] ;
    } ENDFORALL ;
  }

  template<class T> void gpustoreRepI<T>::scatter(const dMap &m, storeRepP &st,
                                               const entitySet &context) {
    const_gpustore<T> s(st) ;
    T *base_ptr = get_base_ptr() ;
    fatal((context != EMPTY) && (base_ptr == 0)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
    FORALL(context,i) {
      base_ptr[m[i]] = s[i] ;
    } ENDFORALL ;
  }

  //*******************************************************************/

  template <class T>
  inline int gpustoreRepI<T>::get_mpi_size( IDENTITY_CONVERTER c,
                                         const entitySet &eset)
  {
    return ( sizeof(T)*eset.size() );
  }
  //*******************************************************************/
  template <class T>
  inline int gpustoreRepI<T>::get_estimated_mpi_size( IDENTITY_CONVERTER c,
                                                   const entitySet &eset)
  {
    return ( sizeof(T)*eset.size() );
  }
  //*******************************************************************/
  template <class T>
  inline int gpustoreRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c,
                                         const entitySet &eset)
  {
    
    int       size ;
    int numBytes = 0 ;
    entitySet :: const_iterator ci;
    typedef data_schema_traits<T> converter_traits;

    fatal((eset - domain()) != EMPTY);

    T *base_ptr = get_base_ptr() ;
    entitySet sdom = eset & domain() ;
    for( ci = sdom.begin(); ci != sdom.end(); ++ci) {
      typename converter_traits::Converter_Type cvtr(base_ptr[*ci]);
      size      = cvtr.getSize();
      numBytes += size*sizeof(typename converter_traits::Converter_Base_Type) ;
    }

    numBytes  += eset.size()*sizeof(int);
    return(numBytes) ;
  }
  //*******************************************************************/
  template <class T>
  inline int gpustoreRepI<T>::get_estimated_mpi_size( USER_DEFINED_CONVERTER c,
                                                   const entitySet &eset)
  {
    int numBytes = 0 ;
    numBytes = eset.size()*50*sizeof(double);
    numBytes  += eset.size()*sizeof(int);
    return(numBytes) ;
  }

  //*******************************************************************/
  template <class T> int gpustoreRepI<T>::pack_size( const entitySet &eset) {

    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    return get_mpi_size( schema_converter(), eset );
  }

  //*******************************************************************/
  template <class T> int gpustoreRepI<T>::estimated_pack_size( const entitySet &eset) {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    return get_estimated_mpi_size( schema_converter(), eset );
    
  }
  

  template<class T> int gpustoreRepI<T>::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;
    typedef typename data_schema_traits<T>::Schema_Converter
      schema_converter ;

    return get_mpi_size(schema_converter(), packed) ;
  }

  //*******************************************************************/
  template <class T>
  inline void gpustoreRepI<T>::packdata(IDENTITY_CONVERTER c, void *outbuf,
                                     int &position, int outcount,
                                     const entitySet &eset )
  {
    T *base_ptr = get_base_ptr() ;
#ifdef USE_CUDA_RT
    cudaDeviceSynchronize() ;
#endif
    for( size_t i = 0; i < eset.num_intervals(); i++) {
      const Loci::int_type begin = eset[i].first ;
      int t = eset[i].second - eset[i].first + 1 ;

#ifdef USE_CUDA_RT
      cudaError_t err = cudaMemcpy(((char *)outbuf)+position,base_ptr+begin,
				   sizeof(T)*t,cudaMemcpyDeviceToHost) ;
      if(err!= cudaSuccess) {
	cerr << "cudaMemcpy failed in gpustoreRepI::packdata" << endl ;
	Loci::Abort() ;
      }
#endif
      position += t*sizeof(T) ;
      //      MPI_Pack( &base_ptr[begin], t*sizeof(T), MPI_BYTE, outbuf,outcount,
      //                &position, MPI_COMM_WORLD) ;
    }
#ifdef USE_CUDA_RT
    cudaDeviceSynchronize() ;
#endif
    
  }

  //*******************************************************************/

  template <class T>
  inline void gpustoreRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf,
                                      int &position, int outcount,
                                      const entitySet &eset )
  {

    entitySet :: const_iterator ci;
    entitySet  ecommon;

    ecommon = domain()&eset;
    //-------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the
    // container for allocation purpose
    //-------------------------------------------------------------------
    size_t  incount = 0;
    int     stateSize, maxStateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    T *base_ptr = get_base_ptr() ;
    for( ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( base_ptr[*ci] );
      stateSize    = cvtr.getSize();
      incount     += stateSize;
      maxStateSize = max( maxStateSize, stateSize );
    }

    typename schema_traits::Converter_Base_Type *inbuf;

    int typesize = sizeof(typename schema_traits::Converter_Base_Type);

    inbuf = new typename schema_traits::Converter_Base_Type[maxStateSize];

    //-------------------------------------------------------------------
    // Collect state data from each object and put into 1D array
    //-------------------------------------------------------------------

    for( ci = ecommon.begin(); ci != ecommon.end(); ++ci) {
      typename schema_traits::Converter_Type cvtr( base_ptr[*ci]);
      cvtr.getState( inbuf, stateSize);

      incount =  sizeof(int);
      MPI_Pack(&stateSize, 1, MPI_INT, outbuf, outcount,&position,
               MPI_COMM_WORLD);

      incount =  stateSize*typesize;
      MPI_Pack(inbuf, incount, MPI_BYTE, outbuf, outcount, &position,
               MPI_COMM_WORLD) ;
    }
    delete [] inbuf;
  }
  //*******************************************************************/

  template <class T>
  void gpustoreRepI<T>::pack( void *outbuf, int &position, int &size,
                           const entitySet &usr_eset )
  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    warn(usr_eset-domain() != EMPTY) ;

    packdata( schema_converter(), outbuf, position, size, usr_eset);

  }

  //*********************************************************************/
  template <class T>
  inline void gpustoreRepI<T>::unpackdata(IDENTITY_CONVERTER c, void *inbuf,
                                       int &position,  int insize,
                                       const sequence &seq)
  {

    T *base_ptr = get_base_ptr() ;
    for(size_t i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
        const Loci::int_type stop = seq[i].second ;
        for(Loci::int_type indx = seq[i].first; indx != stop-1; --indx) {
#ifdef USE_CUDA_RT
	  cudaError_t err = cudaMemcpy(base_ptr+indx,(char*)inbuf+position, 
				       sizeof(T),cudaMemcpyHostToDevice) ;
	  if(err!= cudaSuccess) {
	    cerr << "cudaMemcpy failed in gpustoreRepI::unpackdata" << endl ;
	    Loci::Abort() ;
	  }
	  position += sizeof(T) ;
	  //          MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
	  //                      sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
#endif
	}
      } else {
	Loci::int_type indx = seq[i].first ;
	int t = seq[i].second - seq[i].first + 1 ;
#ifdef USE_CUDA_RT
	cudaError_t err = cudaMemcpy(base_ptr+indx, (char*)inbuf+position,
				     t*sizeof(T),cudaMemcpyHostToDevice) ;
	if(err!= cudaSuccess) {
	  cerr << "cudaMemcpy failed in gpustoreRepI::unpackdata" << endl ;
	  Loci::Abort() ;
	}
	position += t*sizeof(T) ;
	//        MPI_Unpack( inbuf, insize, &position, &base_ptr[indx],
	//                    t*sizeof(T), MPI_BYTE, MPI_COMM_WORLD) ;
#endif
      }
    }
#ifdef USE_CUDA_RT
    cudaDeviceSynchronize() ;
#endif
    
  }
    
  //*********************************************************************/
  template <class T>
  inline void gpustoreRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf,
                                        int &position, int insize, const sequence &seq)
  {

    sequence :: const_iterator ci;

    //-------------------------------------------------------------------
    // Get the sum of each object size and maximum size of object in the
    // container for allocation purpose
    //-------------------------------------------------------------------
    int  stateSize, outcount;

    typedef  data_schema_traits<T> converter_traits;
    typename converter_traits::Converter_Base_Type *outbuf;

    int typesize = sizeof(typename converter_traits::Converter_Base_Type);
    T *base_ptr = get_base_ptr() ;

    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      if( !store_domain.inSet( *ci ) ) {
        std::cout << "Warning: entity not in entitySet  : "
                  << *ci << std::endl;
        continue;
      }
      outcount = sizeof(int);
      MPI_Unpack(inbuf, insize, &position, &stateSize, 1,
                 MPI_INT, MPI_COMM_WORLD) ;

      outbuf = new typename converter_traits::Converter_Base_Type[stateSize];

      outcount = stateSize*typesize;
      MPI_Unpack(inbuf, insize, &position, outbuf, outcount,
                 MPI_BYTE, MPI_COMM_WORLD) ;

      typename converter_traits::Converter_Type cvtr( base_ptr[*ci] );
      cvtr.setState( outbuf, stateSize);
      delete [] outbuf;
    }

  }


  //*******************************************************************/

  template <class T>
  void gpustoreRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq)
  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    unpackdata(schema_converter(), ptr, loc, size, seq);
  }

  //**********************************************************************/

  template<class T>
  frame_info gpustoreRepI<T>::get_frame_info() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }
  
  template<class T>
  frame_info gpustoreRepI<T>::get_frame_info(IDENTITY_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = 1 ;
    return fi ;
  }

  template<class T>
  frame_info gpustoreRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g) {

    entitySet dom = domain() ;
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = 1 ;
    int stateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    T *base_ptr = get_base_ptr() ;
    for(entitySet::const_iterator ci = dom.begin(); ci != dom.end(); ++ci)
      for(int ivec = 0; ivec < fi.size; ivec++){
        typename schema_traits::Converter_Type cvtr(base_ptr[(*ci)*fi.size+ivec] );
        stateSize = cvtr.getSize();
        fi.second_level.push_back(stateSize) ;
      }
    return fi ;
  }

  
  template<class T>
  DatatypeP gpustoreRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T>
  DatatypeP gpustoreRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T>
  DatatypeP gpustoreRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  //**************************************************************************/

  template<class T>
  void gpustoreRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, eset) ;
  }


#ifdef H5_HAVE_PARALLEL
  template<class T>
  void gpustoreRepI<T>::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& eset, hid_t xfer_plist_id) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5writeP(group_id, dataspace, dataset, dimension, name, traits_output_type, eset, xfer_plist_id) ;
  }
#endif
  
  //**********************************************************************/
  template <class T>
  void gpustoreRepI<T> :: hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &eset)  const
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      size_t tmp = 0 ;
      T *base_ptr = get_base_ptr() ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si)
        tmp_array[tmp++] = base_ptr[*si] ;


      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

#ifdef H5_HAVE_PARALLEL
  template <class T>
  void gpustoreRepI<T> :: hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, const entitySet &eset, hid_t xfer_plist_id)  const
  {
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    T* tmp_array = new T[dimension] ;
    size_t tmp = 0 ;
    T *base_ptr = get_base_ptr() ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si)
      tmp_array[tmp++] = base_ptr[*si] ;


    H5Dwrite(dataset, datatype, memspace, dataspace, xfer_plist_id, tmp_array) ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
  }
#endif  
  //*********************************************************************/

  template <class T>
  void gpustoreRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset) const
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
      T *base_ptr = get_base_ptr() ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
        typename schema_traits::Converter_Type cvtr(base_ptr[*si]);
        cvtr.getState(tmp_array+tmp, stateSize) ;
        tmp +=stateSize ;
      }
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

#ifdef H5_HAVE_PARALLEL
  template <class T>
  void gpustoreRepI<T>::hdf5writeP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset, hid_t xfer_plist_id) const
  {
    typedef data_schema_traits<T> schema_traits ;
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    dtype* tmp_array = new dtype[dimension] ;
    size_t tmp = 0 ;
    int stateSize = 0 ;
    T *base_ptr = get_base_ptr() ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si) {
      typename schema_traits::Converter_Type cvtr(base_ptr[*si]);
      cvtr.getState(tmp_array+tmp, stateSize) ;
      tmp +=stateSize ;
    }
    H5Dwrite(dataset, datatype, memspace, dataspace, xfer_plist_id, tmp_array) ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
  }
#endif
  //*********************************************************************/
  template<class T>
  void gpustoreRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset) ;
  }
#ifdef H5_HAVE_PARALLEL
  template<class T>
  void gpustoreRepI<T>::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5readP(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset, xfer_plist_id) ;
  }
#endif
  
  template<class T>
  void  gpustoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset)
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
        std::string error = "H5Dread: failed" ;
        throw StringError(error) ;
      }
      T *base_ptr = get_base_ptr() ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si)
        base_ptr[*si] = tmp_array[tmp++] ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }
#ifdef H5_HAVE_PARALLEL
  template<class T>
  void  gpustoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER c, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    storeRepP qrep = getRep() ;
    int rank = 1 ;
    DatatypeP dp = qrep->getType() ;
    hid_t datatype = dp->get_hdf5_type() ;
    hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
    T* tmp_array = new T[dimension] ;
    size_t tmp = 0 ;
    hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                        xfer_plist_id, tmp_array) ;
    if(err < 0) {
      std::string error = "H5Dread: failed" ;
      throw StringError(error) ;
    }
    T *base_ptr = get_base_ptr() ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end();++si)
      base_ptr[*si] = tmp_array[tmp++] ;
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
  }
#endif
  //*********************************************************************/

  template<class T>
  void  gpustoreRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset)
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
	std::string error = "H5Dread: failed" ;
	throw StringError(error) ;
      }
      size_t tmp = 0 ;
      int bucsize ;
      size_t indx = 0 ;
      T *base_ptr = get_base_ptr() ;
      for(entitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) {
	typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*si]);
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
  template<class T>
  void  gpustoreRepI<T>::hdf5readP(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension, const char* name, USER_DEFINED_CONVERTER c, frame_info &fi, entitySet &eset,hid_t xfer_plist_id)
  {
    typedef data_schema_traits<T> schema_traits ;
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
    size_t tmp = 0 ;
    int bucsize ;
    size_t indx = 0 ;
    T *base_ptr = get_base_ptr() ;
    for(entitySet::const_iterator si = eset.begin(); si != eset.end(); ++si) {
      typename data_schema_traits<T>::Converter_Type cvtr(base_ptr[*si]);
      bucsize = vint[indx++] ;
      cvtr.setState(tmp_array+tmp, bucsize) ;
      tmp += bucsize ;
    }
    H5Sclose(memspace) ;
    H5Tclose(datatype) ;
    delete [] tmp_array ;
  }
#endif
  //*********************************************************************/

} // end of namespace Loci

#endif
