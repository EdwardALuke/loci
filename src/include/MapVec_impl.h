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
#ifndef MAPVEC_IMPL_H
#define MAPVEC_IMPL_H

#include <MapVec_def.h>
#include <DMapVec_def.h>

namespace Loci {
  
  template<int M> void MapVecRepI<M>::allocate(const entitySet &ptn) {

    if(alloc_id < 0)
      alloc_id = getStoreAllocateID() ;

    storeAllocateData[alloc_id].template allocBasic<VEC>(ptn,1) ;
    store_domain = storeAllocateData[alloc_id].allocset ; ;
    base_ptr = (((VEC *)storeAllocateData[alloc_id].base_ptr) -
		storeAllocateData[alloc_id].base_offset );
    store_domain = ptn ;
    dispatch_notify() ;
    return ;
    
  }

  //*************************************************************************/

  template<int M> MapVecRepI<M>::~MapVecRepI() {
    if(alloc_id>=0) {
      storeAllocateData[alloc_id].template release<VEC>() ;
      releaseStoreAllocateID(alloc_id) ;
      alloc_id = -1 ;
    }
    return ;
  }

  //*************************************************************************/

  template<int M> storeRepP MapVecRepI<M>::get_map()  {
    store<int> sizes ;
    sizes.allocate(store_domain) ;
    FORALL(store_domain,i) {
      sizes[i] = M ;
    } ENDFORALL ;
    multiMap result ;
    result.allocate(sizes) ;
    FORALL(store_domain,i) {
      for(int j=0;j<M;++j) 
        result.begin(i)[j] = base_ptr[i][j] ;
    } ENDFORALL ;
    return result.Rep() ;
  }

  //*************************************************************************/

  template<int M> storeRep *MapVecRepI<M>::new_store(const entitySet &p) const {
    return new MapVecRepI<M>(p) ;
  }
  template<int M> storeRep *MapVecRepI<M>::new_store(const entitySet &p, const int* cnt) const {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a MapVec " << endl ;
    return sp ;
  }
  //*************************************************************************/

  template<int M> entitySet MapVecRepI<M>::domain() const {
    return store_domain ;
  }

  //*************************************************************************/

  template<int M> entitySet MapVecRepI<M>::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    if(d.num_intervals() < IMAGE_THRESHOLD) {
      for(size_t i=0;i<d.num_intervals();++i)
        codomain += Loci::image_section(&base_ptr[d[i].first][0],
                                        &base_ptr[d[i].second+1][0]) ;
    } else {
      std::vector<int> img(d.size()*M) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(size_t i=0;i<d.num_intervals();++i)
        for(int j=d[i].first;j!=d[i].second+1;++j) {
          for(int k=0;k<M;++k) {
            *ins = base_ptr[j][k] ;
            ++ins ;
          }
        }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        codomain += *ins ;
    }
    return codomain ;
  }

  //*************************************************************************/

  template<int M> std::pair<entitySet,entitySet>
  MapVecRepI<M>::preimage(const entitySet &codomain) const {
    entitySet domaini ;
    entitySet domainu ;
    FORALL(store_domain,i) {
      bool vali = true ;
      bool valu = false ;
      for(int j=0;j<M;++j) {
        bool in_set = codomain.inSet(base_ptr[i][j]) ;
        vali = vali && in_set ;
        valu = valu || in_set ;
      }
      if(vali)
        domaini += i ;
      if(valu)
        domainu += i ;
    } ENDFORALL ;
    return std::make_pair(domaini,domainu) ;
  }

  //*************************************************************************/
  
  template<int M> std::ostream &MapVecRepI<M>::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      s << base_ptr[ii][0] ;
      for(int j=1;j<M;++j)
        s << " " << base_ptr[ii][j] ;
      s << std::endl ;
    } ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }

  //*************************************************************************/

  template<int M> std::istream &MapVecRepI<M>::Input(std::istream &s) {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    allocate(e) ;

    FORALL(e,ii) {
      for(int j=0;j<M;++j)
        s >> base_ptr[ii][j] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  //*************************************************************************/
  template<int M> void MapVec<M>::notification() {
    NPTR<MapVecType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  template<int M> inline std::ostream & operator<<(std::ostream &s, const MapVec<M> &m){
    return m.Print(s) ;
  }
    
  template<int M> inline std::istream & operator>>(std::istream &s, MapVec<M> &m) {
    return m.Input(s) ;
  }
  //*************************************************************************/
  template<int M> void const_MapVec<M>::notification() {
    NPTR<MapVecType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  //*************************************************************************/

  template<int M> store_instance::instance_type
  const_MapVec<M>::access() const { return READ_ONLY; }
    
  //*************************************************************************/

  template<int M> inline std::ostream & operator<<(std::ostream &s,
                                                   const const_MapVec<M> &m){
    return m.Print(s) ;
  }

  //*************************************************************************/

  template<int M> storeRepP MapVecRepI<M>::MapRemap(const dMap &dm, const dMap &rm) const {
    entitySet newdomain = dm.domain() & domain() ;
    std::pair<entitySet,entitySet> mappimage = preimage(rm.domain()) ;
    newdomain &= mappimage.first ;
    entitySet mapimage = dm.image(newdomain) ;
    MapVec<M> s ;
    s.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(s.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(dm,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(rm,mapimage) ;
    
    return s.Rep() ;
  }
  //*************************************************************************/
  template<int M> void MapVecRepI<M>::compose(const dMap &m,
                                              const entitySet &context)  {
    fatal((context-store_domain) != EMPTY) ;
    fatal((image(context)-m.domain()) != EMPTY) ;
    entitySet dom = m.domain() ;
    FORALL(context,i) {
      for(int j=0;j<M;++j) {
	if(dom.inSet(base_ptr[i][j]))
          base_ptr[i][j] = m[base_ptr[i][j]] ;
	else
	  base_ptr[i][j] = -1 ;
      }
    } ENDFORALL ;
  }
  //*************************************************************************/
  template<int M> storeRepP MapVecRepI<M>::remap(const dMap &m) const {
    cerr << "Map remap not implemented" << endl ;
    return MapRemap(m,m) ;
  }

  //*************************************************************************/
  template<int M> void MapVecRepI<M>::copy(storeRepP &st,
                                           const entitySet &context)  {
    const_MapVec<M> s(st) ;
    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;
    FORALL(context,i) {
      for(int j=0;j<M;++j)
        base_ptr[i][j] = s[i][j] ;
    } ENDFORALL ;
  }
  //*************************************************************************/
  template<int M> void MapVecRepI<M>::gather(const dMap &m, storeRepP &st,
                                             const entitySet &context)  {
    const_MapVec<M> s(st) ;
    fatal(base_ptr == 0) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ;
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      for(int j=0;j<M;++j)
        base_ptr[i][j] = s[m[i]][j] ;
    } ENDFORALL ;
  }
  //*************************************************************************/
  template<int M> void MapVecRepI<M>::scatter(const dMap &m, storeRepP &st,
                                              const entitySet &context)  {
    const_MapVec<M> s(st) ;
    fatal((base_ptr == 0) &&(context != EMPTY)) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
    FORALL(context,i) {
      for(int j=0;j<M;++j)
        base_ptr[m[i]][j] = s[i][j] ;
    } ENDFORALL ;
  }

  //*************************************************************************/
  template <int M> int MapVecRepI<M>::pack_size( const entitySet &eset) {
    fatal((eset - domain()) != EMPTY);
    int init_size = M;

    int size = cpypacksize(&init_size,1) ;

    entitySet::const_iterator ci;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
      size += cpypacksize(&base_ptr[*ci][0],M) ;

    return size ;
  }

  //*************************************************************************/
  template <int M> int MapVecRepI<M>::estimated_pack_size( const entitySet &eset) {
   
    int size ;
    size = sizeof(int)*eset.size()*M  + sizeof(int);
    return size ; 
  }
   


  
  template<int M> int MapVecRepI<M>::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;

    int init_size = M;

    int size = cpypacksize(&init_size,1) ;

    entitySet::const_iterator ci;
    for( ci = packed.begin(); ci != packed.end(); ++ci)
      size += cpypacksize(&base_ptr[*ci][0],M) ;

    return size ;
  }
  //*************************************************************************/
  template <int M> void MapVecRepI<M>::pack(void *outbuf, int &position, 
                                            int &outcount, const entitySet &eset ) 
  {
    int init_size = M;
    cpypack(outbuf,position,outcount,&init_size,1) ;

    entitySet::const_iterator ci;
    for( ci = eset.begin(); ci != eset.end(); ++ci)
      cpypack(outbuf,position,outcount,&base_ptr[*ci][0],M) ;
  }
  //*************************************************************************/

  template <int M> void MapVecRepI<M>::unpack(void *inbuf, int &position, int &insize, 
                                              const sequence &seq) {

    int init_size;

    cpyunpack(inbuf,position,insize,&init_size,1) ;

    if(init_size != M) {
      std::cerr << " Invalid MapVec container for unpack data " << std::endl;
      Loci::Abort();
    }

    sequence::const_iterator ci;
    for( ci = seq.begin(); ci != seq.end(); ++ci)
      cpyunpack(inbuf,position,insize,&base_ptr[*ci][0],M) ;
  }
  //*************************************************************************/
  
  template<int M> storeRepP
  MapVecRepI<M>::expand(entitySet &out_of_dom,
                        std::vector<entitySet> &ptn) {
    cerr<<"start mapVec expand" << endl;
    int *recv_count = new int[MPI_processes] ;
    int *send_count = new int[MPI_processes] ;
    int *send_displacement = new int[MPI_processes] ;
    int *recv_displacement = new int[MPI_processes] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy(MPI_processes), send_clone(MPI_processes) ;
    for(int i = 0; i < MPI_processes; ++i) {
      entitySet tmp = out_of_dom & ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	copy[i].push_back(*ei) ;
      std::sort(copy[i].begin(), copy[i].end()) ;
      send_count[i] = copy[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	send_clone[i].push_back(recv_buf[j]) ;
      std::sort(send_clone[i].begin(), send_clone[i].end()) ;
    }
    std::vector<entitySet> map_entities(MPI_processes) ;
    for(int i = 0; i < MPI_processes; ++i){
      entitySet dom;
      for(unsigned int j = 0; j <send_clone[i].size(); j++){
        if(store_domain.inSet(send_clone[i][j]))dom += send_clone[i][j];
      }
      map_entities[i] = dom;
    }
    
    
    for(int i = 0; i < MPI_processes; ++i) {
      send_count[i] =  map_entities[i].size()*(1+M) ;
    }
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += send_count[i] ;
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 MPI_COMM_WORLD) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) 
      for(entitySet::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = *miv ;
	++size_send ;
        for(unsigned int k = 0; k < M; k++) { 
	  send_map[size_send] =base_ptr[*miv][k] ;
	  ++size_send ;
	}
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  
    HASH_MAP(int, std::vector<int> ) hm ;
   
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i];) {
        for(int k = 0; k < M; ++k)
          hm[recv_map[j]].push_back(recv_map[j+k+1]);
        j += M + 1 ;
      }
    }
    

    entitySet new_domain = store_domain + out_of_dom ;
    MapVec<M> new_map ;
    new_map.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(new_map.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    new_map.allocate(new_domain) ;
    
    for(ei = store_domain.begin(); ei != store_domain.end(); ++ei)
      for(int p = 0; p < M ; ++p)
	new_map[*ei][p] = base_ptr[*ei][p] ;

    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i];) {
        for(unsigned int k = 0; k < M; ++k)
          new_map[recv_map[j]][k] =recv_map[j+k+1];  
        j += M + 1 ;
      }
    }

    storeRepP sp = new_map.Rep() ;
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_map ;
    delete [] recv_map ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
     cerr<<"end mapVec expand" << endl;
    return sp ;
  }//end of expand
  //*************************************************************************/
  template<int M> storeRepP MapVecRepI<M>::freeze() {
    return getRep() ;
  }
  
  template<int M > storeRepP MapVecRepI<M>::thaw() {
    dMapVec<M> dm ;
    dm.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(dm.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    FORALL(store_domain, i) {
      for(int j=0;j<M;++j)
        dm[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;
    return dm.Rep() ;
  }
  //*************************************************************************/
  
  template<int M> void inverseMap (multiMap &result,
                                   const MapVec<M> &input_map,
                                   const entitySet &input_image,
                                   const entitySet &input_preimage) {
    result.Rep()->setDomainKeySpace(MapRepP(input_map.Rep())->getRangeKeySpace()) ;
    MapRepP(result.Rep())->setRangeKeySpace(input_map.Rep()->getDomainKeySpace()) ;
    store<int> sizes ;
    sizes.allocate(input_image) ;
    FORALL(input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    entitySet preloop = input_preimage & input_map.domain() ;
    FORALL(preloop,i) {
      for(int k=0;k<M;++k)
        if(input_image.inSet(input_map[i][k]))
          sizes[input_map[i][k]] += 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(preloop,i) {
      for(int k=0;k<M;++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) {
          sizes[elem] -= 1 ;
          FATAL(sizes[elem] < 0) ;
          result[elem][sizes[elem]] = i ;
        }
      }
    } ENDFORALL ;
#ifdef DEBUG
    FORALL(input_image,i) {
      FATAL(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
  }      
  template<int M> void inverseMap (dmultiMap &result,
                                   const MapVec<M> &input_map,
                                   const entitySet &input_image,
                                   const entitySet &input_preimage) {
    result.Rep()->setDomainKeySpace(MapRepP(input_map.Rep())->getRangeKeySpace()) ;
    MapRepP(result.Rep())->setRangeKeySpace(input_map.Rep()->getDomainKeySpace()) ;
    store<int> sizes ;
    sizes.allocate(input_image) ;
    FORALL(input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    entitySet preloop = input_preimage & input_map.domain() ;
    FORALL(preloop,i) {
      for(int k=0;k<M;++k)
        if(input_image.inSet(input_map[i][k]))
          sizes[input_map[i][k]] += 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(preloop,i) {
      for(int k=0;k<M;++k) {
        int elem = input_map[i][k] ;
        if(input_image.inSet(elem)) {
          sizes[elem] -= 1 ;
          FATAL(sizes[elem] < 0) ;
          result[elem][sizes[elem]] = i ;
        }
      }
    } ENDFORALL ;
#ifdef DEBUG
    FORALL(input_image,i) {
      FATAL(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
  }      

  
  template<int M> 
    DatatypeP MapVecRepI<M>::getType() {
    return DatatypeP(new AtomicType(INT)) ;
  }
  template<int M> 
    frame_info MapVecRepI<M>::get_frame_info() {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  
  
  //*************************************************************************/
  
  template<int M> 
  void MapVecRepI<M>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset)
  {
    warn(true) ;
    /*
    hsize_t dimension[1];
    int  indx=0, rank=1, size;
 
    Loci::HDF5_ReadDomain(group_id, eset);
    Loci::HDF5_ReadVecSize(group_id, &size);

    allocate( eset );

    entitySet::const_iterator ci;

    //---------------------------------------------------------------------------
    // Calculate the offset of each entity in file ....
    //---------------------------------------------------------------------------
    store<unsigned> offset;
    offset.allocate(eset);

    int arraySize = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      offset[*ci] = arraySize;
      arraySize  += size;
    }

    //---------------------------------------------------------------------------
    // Read the data now ....
    //---------------------------------------------------------------------------
    int num_intervals = eset.num_intervals();
    interval *it = new interval[num_intervals];

    for(int i=0;i< num_intervals;i++) it[i] = eset[i];

    int   *data;
    dimension[0] = size*eset.size();

    hid_t vDatatype  = H5T_NATIVE_INT;
    hid_t mDataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t vDataspace = H5Screate_simple(rank, dimension, NULL);
    hid_t vDataset   = H5Dopen( group_id, "MapVec");

    hssize_t  start[]     = {0};  // determines the starting coordinates.
    hsize_t   stride[]    = {1};  // which elements are to be selected.
    hsize_t   block[]     = {1};  // size of element block;
    hssize_t  foffset[]   = {0};  // location (in file) where data is read.
    hsize_t   count[]     = {0};  // how many positions to select from the dataspace

    for( int k = 0; k < num_intervals; k++) {
      count[0] = 0;
      for( int i = it[k].first; i <= it[k].second; i++)
        count[0] +=  size;

      data = new int[count[0]];

      foffset[0] = offset[it[k].first];

      H5Sselect_hyperslab(mDataspace, H5S_SELECT_SET, start,  stride, count, block);
      H5Sselect_hyperslab(vDataspace, H5S_SELECT_SET, foffset,stride, count, block);
      H5Dread( vDataset, vDatatype, mDataspace, vDataspace, H5P_DEFAULT, data);
 
      indx = 0;
      for( int i = it[k].first; i <= it[k].second; i++) {
        for( int ivec = 0; ivec < size; ivec++)
          base_ptr[i][ivec] = data[indx++];
      }
      delete[] data;
    }

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    H5Sclose( mDataspace);
    */
  }

  //*************************************************************************/
    
  template<int M> 
  void MapVecRepI<M>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& usr_eset) const 
  {
    warn(true) ;
    /*
    int     rank = 1;
    hsize_t dimension;

    entitySet eset(usr_eset &domain());
    int arraySize = M*eset.size();
    if( arraySize < 1) return;

    Loci::HDF5_WriteDomain( group_id, eset);
    Loci::HDF5_WriteVecSize(group_id, M);

    std::vector<int> data(arraySize);
    entitySet::const_iterator ci;

    size_t indx = 0;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      for( int ivec = 0; ivec < M; ivec++)
        data[indx++] =  base_ptr[*ci][ivec];
    }

    dimension        = arraySize;
    hid_t vDataspace = H5Screate_simple(rank, &dimension, NULL);
    hid_t vDatatype  = H5T_NATIVE_INT;
    hid_t vDataset   = H5Dcreate(group_id, "MapVec", vDatatype, vDataspace, 
                                 H5P_DEFAULT);
    H5Dwrite(vDataset, vDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

    H5Dclose( vDataset  );
    H5Sclose( vDataspace);
    */
  }
  //*************************************************************************/
#ifdef H5_HAVE_PARALLEL
  template<int M> 
  void MapVecRepI<M>::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &eset, hid_t xfer_plist_id)
  {
    warn(true) ;
  }
  //*************************************************************************/

  template<int M> 
  void MapVecRepI<M>::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& usr_eset, hid_t xfer_plist_id) const 
  {
    warn(true) ;
  }
#endif

} // end of namespace Loci


#endif
