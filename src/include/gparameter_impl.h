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
#ifndef GPARAMETER_IMPL_H
#define GPARAMETER_IMPL_H

namespace Loci {
  //given send_split, return recv_split;
  //or given recv_split, return send_split
  //allow overlap between processes
   std::vector<gEntitySet> transposePtn(const std::vector<gEntitySet> &ptn, MPI_Comm comm);

#ifdef COPY2STORE  
  //**************************************************************************/
  template<class T> storeRepP gParamRepI<T>::copy2store()const{
    param<T> r ;
    r.set_entitySet(gconvert(store_domain));
    *r = attrib_data;
    return r.Rep() ;
  }
#endif
  //**************************************************************************/
  template<class T>
  std::ostream &gParamRepI<T>::Print(std::ostream &s) const
  {
    gEntitySet dom = domain() ;
    if(dom == ~GEMPTY) {
      Loci::streamoutput(&attrib_data,1,s) ;
    } else {
      s << '{' << domain() << std::endl ;
      Loci::streamoutput(&attrib_data,1,s) ;
      s << '}' << std::endl ;
    }
    return s ;
  }

  //**************************************************************************/

  template<class T>
  std::istream &gParamRepI<T>::Input(std::istream &s)
  {
    gEntitySet e ;
    char ch ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      s.putback(ch) ;
      e = ~GEMPTY ;
      allocate(e) ;
      attrib_data = T() ;
      Loci::streaminput(&attrib_data,1,s) ;
      return s ;
    }

    s >> e ;
    allocate(e) ;

    attrib_data = T() ;
    Loci::streaminput(&attrib_data,1,s) ;

    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading gParameter" << std::endl ;
      s.putback(ch) ;
    }

    return s ;
  }
  
  //**************************************************************************/

  template<class T>
  DatatypeP gParamRepI<T>::getType()const {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return getType(schema_converter()) ;
  }
  
  //**************************************************************************/ 

  template<class T>
  DatatypeP gParamRepI<T>::getType(IDENTITY_CONVERTER g)const {
    typedef data_schema_traits<T> traits_type;
    return(traits_type::get_type()) ;
  }
  
  //**************************************************************************/

  template<class T>
  DatatypeP gParamRepI<T>::getType(USER_DEFINED_CONVERTER g)const {
    typedef data_schema_traits<T> schema_traits ;
    typedef typename schema_traits::Converter_Base_Type dtype;
    typedef data_schema_traits<dtype> traits_type;
    return(traits_type::get_type()) ;
  }

  //**************************************************************************/

  template <class T>
  int gParamRepI<T>::pack_size( const gEntitySet &eset)const
  {
    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    return get_mpi_size( schema_converter(), eset );
  }
 
  //**************************************************************************/

  template <class T>
  int gParamRepI<T>::get_mpi_size( IDENTITY_CONVERTER c, const gEntitySet &eset)const
  {

    return( sizeof(T) ) ;
  }
  
  //   //**************************************************************************/

  template <class T>
  int gParamRepI<T>::get_mpi_size( USER_DEFINED_CONVERTER c, const gEntitySet &eset)const
  {
    typedef data_schema_traits<T> schema_traits;

    typename schema_traits::Converter_Type cvtr(const_cast<T&>(attrib_data));
    int arraySize = cvtr.getSize() ;

    return(arraySize*sizeof(typename schema_traits::Converter_Base_Type) + sizeof(int));
  }
  
  //**************************************************************************/

  template <class T>
  void gParamRepI<T>::pack(void *ptr, int &loc, int size, const gEntitySet &e )const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;

    packdata( schema_converter(), ptr, loc, size);
  }

  //**************************************************************************/
  template <class T>
  void gParamRepI<T>::packdata( IDENTITY_CONVERTER c, void *outbuf, int &position,
                                int outcount )const
  {
    MPI_Pack( const_cast<T*>(&attrib_data), sizeof(T), MPI_BYTE, outbuf, outcount, &position,
              MPI_COMM_WORLD) ;
  }

  //**************************************************************************/

  template <class T>
  void gParamRepI<T>::packdata( USER_DEFINED_CONVERTER c, void *outbuf,
                                int &position, int outcount )const
  {
    int stateSize;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    int typesize = sizeof(dtype);

    typename data_schema_traits<T>::Converter_Type cvtr(const_cast<T&>(attrib_data));

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
  void gParamRepI<T>::unpack(const void *ptr, int &loc, int size)  {

    typedef typename
      data_schema_traits<T>::Schema_Converter schema_converter;

    unpackdata( schema_converter(), ptr, loc, size);
  }
  
  //**************************************************************************/
 
  template <class T>
  gStoreRepP gParamRepI<T>::remap(const gMap &m) const{
    gParam<T> r ;
    r.set_entitySet(m.image(m.domain()&domain())) ;
    *r = attrib_data ;
    r.set_domain_space(domain_space);
    return r.Rep() ;
  }
  
  //**************************************************************************/
 
  template <class T>
  gStoreRepP gParamRepI<T>::clone() const{
    gParam<T> result;
    *result = attrib_data;
    result.set_entitySet(store_domain);
    result.set_domain_space(domain_space);
    return result.Rep();
  }

  //**************************************************************************/

  template <class T>
  void gParamRepI<T>::unpackdata( IDENTITY_CONVERTER c, const void *inbuf, int &position,
                                  int insize)
  {

    
    MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &attrib_data, sizeof(T),
                MPI_BYTE, MPI_COMM_WORLD) ;

  }

  //***********************************************************************/

  template <class T>
  void gParamRepI<T>::unpackdata( USER_DEFINED_CONVERTER c, const void *inbuf,
                                  int &position, int insize)
  {

    int  stateSize, outcount;

    typedef data_schema_traits<T> schema_traits;
    typedef typename schema_traits::Converter_Base_Type dtype;

    MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &stateSize, 1,
                MPI_INT, MPI_COMM_WORLD) ;
    std::vector<dtype> outbuf(stateSize);

    outcount = stateSize*sizeof(dtype);
    MPI_Unpack( const_cast<void*>(inbuf), insize, &position, &outbuf[0], outcount,
                MPI_BYTE, MPI_COMM_WORLD) ;
    typename schema_traits::Converter_Type  cvtr( attrib_data );
    cvtr.setState( &outbuf[0], stateSize);

  }

  
  //**************************************************************************/

  // template<class T>
  // inline std::ostream & operator<<(std::ostream &s, const const_gparam<T> &t)
  // { return t.Print(s) ; }

  //**************************************************************************/

  template<class T> gStoreRepP gParamRepI<T>::
  split_redistribute(const std::vector<gEntitySet>& dom_ptn, MPI_Comm comm)const {
    // for a gParameter, we just redistribute its domain
    gEntitySet dom = domain() ;
    gEntitySet new_all ;
    for(size_t i=0;i<dom_ptn.size();++i)
      new_all += dom_ptn[i] ;
    gEntitySet out = dom - new_all ;

    std::vector<gEntitySet> old_dist = g_all_collect_vectors<gEntity>(dom, comm) ;
    gEntitySet old_all ;
    for(size_t i=0;i<old_dist.size();++i)
      old_all += old_dist[i] ;

    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;

    // get the new domain
    gEntitySet new_dom = old_all & dom_ptn[rank] ;
    new_dom += out ;

    gParam<T> np ;
    np.set_entitySet(new_dom) ;
    *np = attrib_data ;
    np.set_domain_space(domain_space);
    return np.Rep() ;
  }
  
  //**************************************************************************/

  template<class T> gStoreRepP gParamRepI<T>::
  redistribute(const std::vector<gEntitySet>& dom_split,
               MPI_Comm comm)const{
    if(store_domain == gEntitySet(gInterval(GUNIVERSE_MIN,GUNIVERSE_MAX))){
      gParam<T> np ;
      np.set_entitySet(store_domain) ;
      *np = attrib_data ;
      np.set_domain_space(domain_space);
      return np.Rep() ;
    }else{
      std::vector<gEntitySet> recv_split = transposePtn(dom_split, comm);
      gEntitySet new_dom;
      for(size_t i=0;i<recv_split.size();++i)new_dom += recv_split[i];
      gEntitySet old_all = g_all_collect_entitySet<gEntity>(store_domain, comm);
      new_dom = new_dom&old_all;
      gParam<T> np ;
      np.set_entitySet(new_dom) ;
      *np = attrib_data ;
      np.set_domain_space(domain_space);
      return  np.Rep() ;
    }
  }
  
  //**************************************************************************/
  
  template<class T> gStoreRepP gParamRepI<T>::
  redistribute(const std::vector<gEntitySet>& dom_split,
               const gMap& remap, MPI_Comm comm)const{
   
    if(store_domain ==  gEntitySet(gInterval(GUNIVERSE_MIN,GUNIVERSE_MAX))){
      gParam<T> np ;
      np.set_entitySet(store_domain) ;
      *np = attrib_data ;
      np.set_domain_space(domain_space);
      return np.Rep() ;
    }else{
      gParam<T> np ;
      np = redistribute(dom_split, comm);
      np.set_domain_space(domain_space);
      return np.remap(remap);
    } 
  }

  //**************************************************************************/
  template<class T>
  frame_info gParamRepI<T>::get_frame_info()const {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    return get_frame_info(schema_converter()) ;
  }
  //**************************************************************************/
  template<class T>
  frame_info gParamRepI<T>::get_frame_info(IDENTITY_CONVERTER g)const {
    frame_info fi ;
    fi.is_stat = 0 ;
    fi.size = 1 ;
    return fi ;
  }
  //**************************************************************************/
  template<class T>
  frame_info gParamRepI<T>::get_frame_info(USER_DEFINED_CONVERTER g)const {
    frame_info fi ;
    fi.is_stat = 1 ;
    fi.size = 1 ;
    int stateSize = 0;
    typedef data_schema_traits<T> schema_traits ;
     T tmp = attrib_data;
    typename schema_traits::Converter_Type cvtr(tmp);
    stateSize = cvtr.getSize();
    fi.second_level.push_back(stateSize) ;
    return fi ;
  }


  //**************************************************************************/
  template <class T>
  void gParamRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension,
                               const char* name, IDENTITY_CONVERTER g, const gEntitySet &eset) const
  {
    if(dimension != 0) {
      gStoreRepP qrep = getRep() ;
      int rank = 1 ;
      DatatypeP dp = qrep->getType() ;
      hid_t datatype = dp->get_hdf5_type() ;
      hid_t memspace = H5Screate_simple(rank, &dimension, NULL) ;

      H5Dwrite(dataset, datatype, memspace, dataspace, H5P_DEFAULT, &attrib_data) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
    }
  }

  //*********************************************************************/

  template <class T>
  void gParamRepI<T>::hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                               const char* name, USER_DEFINED_CONVERTER g, const gEntitySet &eset) const
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
      int stateSize = 0 ;
      T tmp = attrib_data;
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
  void gParamRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                              const char* name, IDENTITY_CONVERTER g, frame_info &fi, const gEntitySet &en)
  {
    if(dimension != 0) {
      gStoreRepP qrep = getRep() ;
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
      attrib_data = tmp_array[0] ;

      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

  //*************************************************************************/

  template <class T>
  void gParamRepI<T>::hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
			       const char* name, USER_DEFINED_CONVERTER g, frame_info &fi, const gEntitySet &gen)
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
        cerr << "H5Dread() failed" << endl ;
      }
      typename data_schema_traits<T>::Converter_Type cvtr(attrib_data);
      int bucsize = vint[0] ;
      cvtr.setState(tmp_array, bucsize) ;
      H5Sclose(memspace) ;
      H5Tclose(datatype) ;
      delete [] tmp_array ;
    }
  }

  //**************************************************************************/
  template<class T>
  void gParamRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, 
			       hsize_t dimension, const char* name, 
			       frame_info &fi, const gEntitySet &eset)
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5read(group_id, dataspace, dataset, dimension, name, traits_output_type, fi, eset) ;
  }
  
  //**************************************************************************/

  template<class T>
  void gParamRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                               const char* name, const gEntitySet &eset) const
  {
    typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
    schema_converter traits_output_type;
    hdf5write(group_id, dataspace, dataset, dimension, name, traits_output_type, eset) ;
  }
 

  
  //**************************************************************************/

  template<class T>
  inline std::ostream & operator<<(std::ostream &s, const gParam<T> &t)
  { return t.Print(s) ; }

  //**************************************************************************/

  template<class T>
  inline std::istream & operator>>(std::istream &s, gParam<T> &t)
  { return t.Input(s) ; }


  
}

#endif
