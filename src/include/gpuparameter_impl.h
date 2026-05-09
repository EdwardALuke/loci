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
#ifndef GPUPARAMETER_IMPL_H
#define GPUPARAMETER_IMPL_H
#include <string>
namespace Loci {
  //**************************************************************************/
 std::vector<entitySet> all_collect_vectors(entitySet &e,MPI_Comm comm) ;
  // a utility function to remap an entityset
   entitySet
  remap_entitySet(const entitySet& es, const dMap& remap) ;



  template<class T> void gpuparamRepI<T>::allocate(const entitySet &p) {
    if(alloc_id < 0) {
      alloc_id = getGPUStoreAllocateID() ;

      entitySet single = interval(0,0) ;
      GPUstoreAllocateData[alloc_id].template allocBasic<T>(single,1) ;
      base_ptr = (T *)GPUstoreAllocateData[alloc_id].base_ptr ;
      //      *base_ptr = defaultData ;
    }
    store_domain = p ;
    dispatch_notify();
    return ;
  }

  template<class T> void gpuparamRepI<T>::erase(const entitySet& rm) {
    store_domain -= rm ;
    dispatch_notify() ;
  }

  template<class T> void gpuparamRepI<T>::shift(int_type offset) {
    store_domain >>= offset ;
    dispatch_notify() ;
  }

  template<class T> void gpuparamRepI<T>::
  guarantee_domain(const entitySet& include) {
    store_domain += include ;
    dispatch_notify() ;
  }

  //**************************************************************************/

  template<class T> gpuparamRepI<T>::~gpuparamRepI() {
    if(alloc_id>=0) {
      GPUstoreAllocateData[alloc_id].template release<T>() ;
      releaseGPUStoreAllocateID(alloc_id) ;
      alloc_id = -1 ;
    }
    return ;
}

  //**************************************************************************/

  template<class T>
  storeRep *gpuparamRepI<T>::new_store(const entitySet &p) const
  {
    return new gpuparamRepI<T>(p) ;
  }

  template<class T>
  storeRep *gpuparamRepI<T>::new_store(const entitySet &p, const int* cnt) const
  {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a parameter " << endl ;
    return sp ;
  }

  //**************************************************************************/

  template<class T>
  store_type gpuparamRepI<T>::RepType() const
  {
    return GPUPARAMETER ;
  }

  //**************************************************************************/

  template<class T> entitySet gpuparamRepI<T>::domain() const {
    return store_domain ;
  }

  //**************************************************************************/

  template<class T>
  std::ostream &gpuparamRepI<T>::Print(std::ostream &s) const
  {
    entitySet dom = domain() ;
    if(dom == ~EMPTY) {
      Loci::streamoutput(base_ptr,1,s) ;
    } else {
      s << '{' << domain() << std::endl ;
      Loci::streamoutput(base_ptr,1,s) ;
      s << '}' << std::endl ;
    }
    return s ;
  }

  //**************************************************************************/

  template<class T>
  std::istream &gpuparamRepI<T>::Input(std::istream &s)
  {
    entitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      s.putback(ch) ;
      e = ~EMPTY ;
      allocate(e) ;
      *base_ptr = T() ;
      Loci::streaminput(base_ptr,1,s) ;
      return s ;
    }

    s >> e ;
    allocate(e) ;

    *base_ptr = T() ;
    Loci::streaminput(base_ptr,1,s) ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading parameter" << std::endl ;
      s.putback(ch) ;
    }

    return s ;
  }

  template<class T>
  frame_info gpuparamRepI<T>::get_frame_info() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }

  template<class T>
  frame_info gpuparamRepI<T>::get_frame_info(IDENTITY_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = 1 ;
    return fi ;
  }
  template<class T>
  frame_info gpuparamRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g) {
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = 1 ;
    int stateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
    typename schema_traits::Converter_Type cvtr(*base_ptr) ; 
    stateSize = cvtr.getSize();
    fi.second_level.push_back(stateSize) ;

    return fi ;
  }


  template<class T>
  DatatypeP gpuparamRepI<T>::getType() {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  template<class T>
  DatatypeP gpuparamRepI<T>::getType(IDENTITY_CONVERTER g) {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  template<class T>
  DatatypeP gpuparamRepI<T>::getType(USER_DEFINED_CONVERTER g) {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }
  //**************************************************************************/
  // Code to copy from cpu container to gpu container
  template<class T> 
  void gpuparamRepI<T>::copyFrom(const storeRepP &p, entitySet set)  {
    const_param<T> v(p) ;
    T *gpu_base_ptr = get_param() ;
#ifdef USE_CUDA_RT
    cudaError_t err = cudaMemcpy(gpu_base_ptr,&(*v),sizeof(T),
				 cudaMemcpyHostToDevice) ;
    if(err!= cudaSuccess) {
      cerr << "cudaMemcpy failed in gpuMapRepI::copyFrom" << endl ;
      Loci::Abort() ;
    }
    cudaDeviceSynchronize() ;
#endif
  }
  
  // code to copy from gpu container to cpu container
  template<class T>
  void gpuparamRepI<T>::copyTo(storeRepP &p, entitySet set) const {
    param<T> v(p) ;
    const T *gpu_base_ptr = get_param() ;
#ifdef USE_CUDA_RT
    cudaDeviceSynchronize() ;
    cudaError_t err = cudaMemcpy(&(*v),gpu_base_ptr,sizeof(T),
				 cudaMemcpyDeviceToHost) ;
    if(err!= cudaSuccess) {
      cerr << "cudaMemcpy failed in gpuMapRepI::copyFrom" << endl ;
      Loci::Abort() ;
    }
#endif
  }
  
  template<class T>
  void gpuparamRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset) ;
  }

  //**************************************************************************/

  template<class T>
  void gpuparamRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, eset) ;
  }
#ifdef H5_HAVE_PARALLEL
  template<class T>
  void gpuparamRepI<T>::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    warn(true);
  }

  template<class T>
  void gpuparamRepI<T>::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &eset,hid_t xfer_plist_id) const
  {
    warn(true);
  }
#endif
  //**************************************************************************/

  template<class T> gpuparam<T>::~gpuparam() {}

  //**************************************************************************/

  template<class T>
  void gpuparam<T>::notification()
  {
    NPTR<paramType> p(Rep());
    if(p!=0) data = p->get_param() ;
    warn(p==0);
  }

  //**************************************************************************/

  template<class T> const_gpuparam<T>::~const_gpuparam() {}

  //**************************************************************************/

  template<class T>
  void const_gpuparam<T>::notification()
  {
    NPTR<paramType> p(Rep());
    if(p!=0) data = p->get_param() ;
    warn(p==0);
  }

  //**************************************************************************/

  template<class T>
  storeRepP gpuparamRepI<T>::remap(const dMap &m) const
  {
    gpuparam<T> r ;
    r.set_entitySet(m.image(m.domain()&domain())) ;
    *r = *base_ptr ;
    return r.Rep() ;
  }

  template<class T>
  storeRepP gpuparamRepI<T>::freeze() {
    return getRep() ;
  }

  template<class T>
  storeRepP gpuparamRepI<T>::thaw() {
    return getRep() ;
  }
  //**************************************************************************/

  template<class T>
  void gpuparamRepI<T>::copy(storeRepP &st, const entitySet &context)
  {
    gpuparam<T> p(st) ;
    *base_ptr = *p ;
    warn((store_domain - context) != EMPTY) ;
    store_domain = context ;
    dispatch_notify() ;
  }

  // note this method can only be used when the gpuparamRepI<T> (i.e., *this)
  // is NOT connected to any of the containers since we don't perform any
  // kind of notification
  template<class T>
  void gpuparamRepI<T>::fast_copy(storeRepP &st, const entitySet &context)
  {
    storeRepP true_rep = st->getRep();
    gpuparamRepI<T>* p = dynamic_cast<gpuparamRepI<T>*>(&(*true_rep));
    fatal(p==0);
    *base_ptr = p->base_ptr[0];
    warn((store_domain - context) != EMPTY) ;
    store_domain = context ;
  }

  //**************************************************************************/

  template<class T>
  void gpuparamRepI<T>::gather(const dMap &m, storeRepP &st,
                            const entitySet &context)
  {
    gpuparam<T> p(st) ;
    fatal((context - store_domain) != EMPTY) ;
    store_domain = context ;
  }

  //**************************************************************************/

  template<class T>
  void gpuparamRepI<T>::scatter(const dMap &m, storeRepP &st,
                             const entitySet &context)
  {

    fatal((context - store_domain) != EMPTY) ;
    store_domain = m.image(context) ;

    gpuparam<T> p(st) ;
    *base_ptr = *p ; 
  }

  //**************************************************************************/

  template <class T>
  int gpuparamRepI<T>::pack_size( const entitySet &eset)
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    return get_mpi_size( schema_converter(), eset );
  }
  //**************************************************************************/
  template <class T>
  int gpuparamRepI<T>::estimated_pack_size( const entitySet &eset)
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    return get_estimated_mpi_size( schema_converter(), eset );
    
  }
  //**************************************************************************/
  template<class T> int gpuparamRepI<T>::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;    

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    return get_mpi_size(schema_converter(), packed);
  }
  //**************************************************************************/

  template <class T>
  int gpuparamRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {

    return( sizeof(T) ) ;
  }
  //**************************************************************************/

  template <class T>
  int gpuparamRepI<T>::get_estimated_mpi_size( IDENTITY_CONVERTER c, const entitySet &eset)
  {

    return( sizeof(T) ) ;
  }
  //**************************************************************************/

  template <class T>
  int gpuparamRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {
    typedef data_schema_traits<T> schema_traits;

    typename schema_traits::Converter_Type cvtr(*base_ptr);
    int arraySize = cvtr.getSize() ;

    return(arraySize*sizeof(typename schema_traits::Converter_Base_Type) + sizeof(int));
  }
  //**************************************************************************/

  template <class T>
  int gpuparamRepI<T>::get_estimated_mpi_size( USER_DEFINED_CONVERTER c, const entitySet &eset)
  {
    return(50*sizeof(double) + sizeof(int));
  } 
  //**************************************************************************/

  template <class T>
  void gpuparamRepI<T>::pack(void *ptr, int &loc, int &size, const entitySet &e )
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    packdata( schema_converter(), ptr, loc, size);
  }

  //**************************************************************************/
  template <class T>
  void gpuparamRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf, int &position,
                               int outcount )
  {
    MPI_Pack( base_ptr, sizeof(T), MPI_BYTE, outbuf, outcount, &position,
              MPI_COMM_WORLD) ;
  }

  //**************************************************************************/

  template <class T>
  void gpuparamRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf,
                               int &position, int outcount )
  {
    int stateSize;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);

    typename data_schema_traits<T>::Converter_Type cvtr( *base_ptr );

    std::vector<dtype> inbuf(cvtr.getSize());
    cvtr.getState( &inbuf[0], stateSize);

    MPI_Pack(&stateSize, 1, MPI_INT, outbuf, outcount,&position,
             MPI_COMM_WORLD);
    int incount =  stateSize*typesize;
    MPI_Pack(&inbuf[0], incount, MPI_BYTE, outbuf, outcount, &position,
             MPI_COMM_WORLD) ;

  }
  //**************************************************************************/

  template <class T>
  void gpuparamRepI<T>::unpack(void *ptr, int &loc, int &size, const sequence &seq)  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    unpackdata( schema_converter(), ptr, loc, size);
  }


  //**************************************************************************/
  template <class T>
  void gpuparamRepI<T>::unpackdata( IDENTITY_CONVERTER c, void *inbuf, int &position,
                                 int &insize)
  {

    /*
      typedef data_schema_traits<T> traits_type;
      DatatypeP    atom_type = traits_type::get_type();
      MPI_Datatype datatype  = atom_type->get_mpi_type();
    */
    MPI_Unpack( inbuf, insize, &position, base_ptr, sizeof(T),
                MPI_BYTE, MPI_COMM_WORLD) ;

  }

  //***********************************************************************/
  template <class T>
  void gpuparamRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, void *inbuf,
                                 int &position, int &insize)
  {

    int  stateSize, outcount;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    MPI_Unpack( inbuf, insize, &position, &stateSize, 1,
                MPI_INT, MPI_COMM_WORLD) ;
    std::vector<dtype> outbuf(stateSize);

    outcount = stateSize*sizeof(dtype);
    MPI_Unpack( inbuf, insize, &position, &outbuf[0], outcount,
                MPI_BYTE, MPI_COMM_WORLD) ;
    typename schema_traits::Converter_Type  cvtr( *base_ptr );
    cvtr.setState( &outbuf[0], stateSize);

  }

  //***********************************************************************/

  template<class T> store_instance::instance_type
  const_gpuparam<T>::access() const
  { return READ_ONLY; }

  //**************************************************************************/

  template<class T>
  inline std::ostream & operator<<(std::ostream &s, const const_gpuparam<T> &t)
  { return t.Print(s) ; }

  //**************************************************************************/
  template <class T>
  void gpuparamRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension, const char* name, IDENTITY_CONVERTER g, const entitySet &eset) const
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;

      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, base_ptr) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
    }
  }

  //*********************************************************************/

  template <class T>
  void gpuparamRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, const entitySet &eset) const
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
      int stateSize = 0 ;
      T tmp = *base_ptr;
      typename schema_traits::Converter_Type cvtr(tmp);
      cvtr.getState(tmp_array, stateSize) ;
      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, tmp_array) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

  //**************************************************************************/
  template <class T>
  void gpuparamRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, IDENTITY_CONVERTER g, frame_info &fi, const entitySet &en)
  {
    if(dimension != 0) {
      storeRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;
      T* tmp_array = new T[dimension] ;
      hid_t err = H5Dread(dataset,  datatype, memspace, dataspace,
                          H5P_DEFAULT, tmp_array) ;
      if(err < 0) {
        cerr << "H5Dread() failed" << endl ;
      }
      *base_ptr = tmp_array[0] ;

      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

  //*************************************************************************/

  template <class T>
  void gpuparamRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, USER_DEFINED_CONVERTER g, frame_info &fi, const entitySet &en)
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
      typename data_schema_traits<T>::Converter_Type cvtr(*base_ptr);
      int bucsize = vint[0] ;
      cvtr.setState(tmp_array, bucsize) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

  template<class T> storeRepP gpuparamRepI<T>::
  redistribute(const std::vector<entitySet>& dom_ptn, MPI_Comm comm) {
    // for a parameter, we just redistribute its domain
    entitySet dom = domain() ;
    entitySet new_all ;
    for(size_t i=0;i<dom_ptn.size();++i)
      new_all += dom_ptn[i] ;
    entitySet out = dom - new_all ;

    std::vector<entitySet> old_dist = all_collect_vectors(dom, comm) ;
    entitySet old_all ;
    for(size_t i=0;i<old_dist.size();++i)
      old_all += old_dist[i] ;

    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;

    // get the new domain
    entitySet new_dom = old_all & dom_ptn[rank] ;
    new_dom += out ;

    gpuparam<T> np ;
    np.set_entitySet(new_dom) ;
    *np = *base_ptr ;

    return np.Rep() ;
  }

  template<class T> storeRepP gpuparamRepI<T>::
  redistribute(const std::vector<entitySet>& dom_ptn,
               const dMap& remap, MPI_Comm comm) {
    // for a parameter, we just redistribute its domain
    entitySet dom = domain() ;
    std::vector<entitySet> old_dist = all_collect_vectors(dom, comm) ;
    entitySet old_all ;
    for(size_t i=0;i<old_dist.size();++i)
      old_all += old_dist[i] ;

    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;

    // get the new domain
    entitySet new_dom = old_all & dom_ptn[rank] ;
    new_dom = remap_entitySet(new_dom, remap) ;

    gpuparam<T> np ;
    np.set_entitySet(new_dom) ;
    *np = *base_ptr ;

    return np.Rep() ;
  }

  //***************************************************************************

}
#endif
