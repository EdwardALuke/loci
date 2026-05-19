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
#include <istream>
#include <ostream>
#include <iostream>
#include <fstream>

#include <DMap.h>
#include <multiMap.h>
#include <hdf5_readwrite.h>
#include <fact_db.h>
#include <Tools/hash_map.h>
#include <distribute.h>
#include <cstring>

using std::cerr ;
using std::endl ;
using std::ostream ;
using std::istream ;
using std::ofstream ;

namespace Loci {
  extern fact_db *exec_current_fact_db ;

  extern ofstream debugout ;
  
  using std::pair ;
  using std::make_pair ;
  using std::vector ;
  using std::sort ;
  //**************************************************************************/
  storeRepP dMapRepI::expand(entitySet &out_of_dom, std::vector<entitySet> &ptn) {
    MPI_Comm comm = get_exec_comm() ;
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
      sort(copy[i].begin(), copy[i].end()) ;
      send_count[i] = copy[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
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
		  comm) ;  
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	send_clone[i].push_back(recv_buf[j]) ;
      sort(send_clone[i].begin(), send_clone[i].end()) ;
    }

    entitySet dom = domain() ;
    
    std::vector<HASH_MAP(int, int) > map_entities(MPI_processes) ;
    for(int i = 0; i < MPI_processes; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(dom.inSet(*vi))
	  (map_entities[i])[*vi] = attrib_data[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < MPI_processes; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < MPI_processes; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    HASH_MAP(int, int) hm ;
    for(int i = 0; i < MPI_processes; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; j+=2) {
	hm[recv_map[j]] = recv_map[j+1];
      }
    }
    
    for(HASH_MAP(int, int)::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi) 
      attrib_data[hmi->first] = hmi->second ;
    
    dMap dm ;
    dm.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(dm.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    dom = domain() ;
    for(entitySet::const_iterator ei = dom.begin(); ei != dom.end(); ++ei)
      dm[*ei] = attrib_data.elem(*ei) ;
    storeRepP sp = dm.Rep() ;
    delete [] send_buf ;
    delete [] recv_buf ;
    delete [] send_map ;
    delete [] recv_map ;
    delete [] recv_count ;
    delete [] send_count ;
    delete [] send_displacement ;
    delete [] recv_displacement ;
    return sp ;
  }
  
  //****************************************************************************/
  void dMapRepI::allocate(const entitySet &eset)
  {
    entitySet redundant, newSet;
    entitySet :: const_iterator   ci;

    redundant = domain() -  eset;
    newSet    = eset - domain();

    attrib_data.erase_set(redundant) ;
    
    for( ci = newSet.begin(); ci != newSet.end(); ++ci)
      attrib_data[*ci] = 0;
    
    dispatch_notify() ;
  }
  
  //*************************************************************************/
  void dMapRepI::erase(const entitySet &rm) {
    // we need to make sure the erased sets are
    // indeed valid (within the domain)
    entitySet valid = domain() & rm ;
    attrib_data.erase_set(valid) ;
    dispatch_notify() ;
  }

  void dMapRepI::invalidate(const entitySet& valid) {
    entitySet redundant = domain() - valid ;
    erase(redundant) ;
  }
  void dMapRepI::guarantee_domain(const entitySet& include) {
    entitySet new_set = include - domain() ;
    for(entitySet::const_iterator ei=new_set.begin();
        ei!=new_set.end();++ei)
      attrib_data.access(*ei) ;
    dispatch_notify() ;
  }
  
  //************************************************************************/
  
  dMapRepI::~dMapRepI() 
  {}
  
  //**************************************************************************/

  storeRep *dMapRepI::new_store(const entitySet &p) const 
  {
    return new dMapRepI(p)  ;
  }
  storeRep *dMapRepI::new_store(const entitySet &p, const int* count) const 
  {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a dMap " << endl ;
    return sp ;
  }
  //**************************************************************************/
  
  storeRepP dMapRepI::MapRemap(const dMap &dm, const dMap &rm) const 
  {
    dMap s ;
    s.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(s.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    entitySet newdomain = dm.domain() & domain() ;
    pair<entitySet,entitySet> mappimage = preimage(rm.domain()) ;
    newdomain &= mappimage.first ;
    entitySet mapimage = dm.image(newdomain) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(dm,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(rm,mapimage) ;
        
    return s.Rep() ;
  }
  storeRepP dMapRepI::remap(const dMap &newmap) const {
    cerr << "remap should not be called for a DMap!" << endl ;
    return MapRemap(newmap,newmap) ;
  }
  
  storeRepP
  dMapRepI::redistribute(const std::vector<entitySet>& dom_ptn,
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

    ////////////////////////////////////////////////////////////
    // first compute the unpack sequence
    std::vector<sequence> pack_seq(np) ;
    for(int i=0;i<np;++i)
      pack_seq[i] = sequence(dom_split[i]) ;

    std::vector<sequence> unpack_seq =
      transpose_sequence(pack_seq, comm) ;

    std::vector<sequence>().swap(pack_seq) ;
    ////////////////////////////////////////////////////////////

    std::vector<int> send_counts(np,0) ;
    std::vector<int> recv_counts(np,0) ;
    std::vector<int> send_displs(np,0) ;
    std::vector<int> recv_displs(np,0) ;
    int tot_send_size = 0 ;
    int tot_recv_size = 0 ;

    // continue to pack and send the Map
    for(int i=0;i<np;++i)
      send_counts[i] = pack_size(dom_split[i]) ;
    MPI_Alltoall(&send_counts[0], 1, MPI_INT,
                 &recv_counts[0], 1, MPI_INT, comm) ;

    for(int i=1;i<np;++i) {
      send_displs[i] = send_displs[i-1] + send_counts[i-1] ;
      recv_displs[i] = recv_displs[i-1] + recv_counts[i-1] ;
    }
    tot_send_size = send_displs[np-1] + send_counts[np-1] ;
    tot_recv_size = recv_displs[np-1] + recv_counts[np-1] ;

    // pack the buffer
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

    unsigned char* recv_buffer = new unsigned char[tot_recv_size] ;
    MPI_Alltoallv(&send_buffer[0], &send_counts[0],
                  &send_displs[0], MPI_PACKED,
                  &recv_buffer[0], &recv_counts[0],
                  &recv_displs[0], MPI_PACKED, comm) ;
    delete[] send_buffer ;

    dMap nm ;
    nm.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(nm.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    entitySet new_domain ;
    for(int i=0;i<np;++i)
      new_domain += entitySet(unpack_seq[i]) ;
    // if there's any old domain not being redistributed,
    // we will add it also
    entitySet old_dom = dom - total_dom ;
    nm.allocate(new_domain+old_dom) ;
    storeRepP srp = nm.Rep() ;

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

    // copy old_dom contents if any
    for(entitySet::const_iterator ei=old_dom.begin();
        ei!=old_dom.end();++ei)
      nm[*ei] = attrib_data[*ei] ;

    return srp ;
  }
  
  storeRepP dMapRepI::
  redistribute(const std::vector<entitySet>& dom_ptn,
               const dMap& remap, MPI_Comm comm) {
    // this is a push operation, thus the send recv are reversed
    std::vector<P2pCommInfo> send, recv ;
    get_p2p_comm(dom_ptn, domain(), 0, 0, comm, recv, send) ;
    dMap new_map ;
    new_map.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(new_map.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    Loci::Abort();
    return new_map.Rep() ;
  }

  
  // ******************************************************************/
  storeRepP dMapRepI::freeze() {
    Map m ;
    m.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(m.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    m.allocate(domain()) ;
    FORALL(domain(), i) {
      m[i] = attrib_data[i] ;
    } ENDFORALL ;
    return m.Rep() ;
  }
  
  storeRepP dMapRepI::thaw() {
    return getRep() ;
  } 
  
  //**************************************************************************/
  
  void dMapRepI::compose(const dMap &newmap, const entitySet &context) 
  {
    fatal((context-domain()) != EMPTY) ;
    fatal((image(context)-newmap.domain()) != EMPTY) ;

    FORALL(context,i) {
      const int mv = attrib_data.elem(i) ;
      attrib_data[i] = newmap[mv] ;
    } ENDFORALL ;
    
  }
  
  //**************************************************************************/
  
  
  void dMapRepI::copy(storeRepP &st, const entitySet &context) 
  {
    const_dMap s(st) ;

    fatal((context-s.domain()) != EMPTY) ;
    
    FORALL(context,i) {
      attrib_data[i] = s[i] ;
    } ENDFORALL ;
    
  }
  
  //**************************************************************************/

  void dMapRepI::gather(const dMap &m, storeRepP &st, const entitySet &context) 
  {
    const_dMap s(st) ;
    
    fatal((m.image(context) - s.domain()) != EMPTY) ; 
    fatal((context - domain()) != EMPTY) ;

    FORALL(context,i) {
      attrib_data[i] = s[m[i]] ;
    } ENDFORALL ;
  }
  
  //**************************************************************************/
  
  void dMapRepI::scatter(const dMap &m,storeRepP &st, const entitySet &context) 
  {
    const_dMap s(st) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);
    
    FORALL(context,i) {
      attrib_data[m[i]] = s[i] ;
    } ENDFORALL ;
  }
  
  //**************************************************************************/
  
  int dMapRepI::pack_size(const entitySet &e) {
    int *p = 0 ;
    return cpypacksize(p,e.size()) ;
  }

  int dMapRepI::estimated_pack_size(const entitySet &e) {
    
    return e.size()*sizeof(int);
  }

  int dMapRepI::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;
    int *p = 0 ;
    return cpypacksize(p,packed.size()) ;
  }
  
  //**************************************************************************/
  
  void dMapRepI::pack(void *outbuf, int &position, int &outcount, const entitySet &eset) 
  {
    entitySet :: const_iterator ci;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      cpypack(outbuf,position,outcount,&attrib_data[*ci],1) ;
    }
  }
  
  //**************************************************************************/

  void dMapRepI::unpack(void *inbuf, int &position, int &insize, const sequence &seq) 
  {
    sequence:: const_iterator ci;
    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      cpyunpack(inbuf,position,insize,&attrib_data[*ci],1) ;
    }
  }
  
  
  //**************************************************************************/
  
  entitySet dMapRepI::domain() const 
  {
    entitySet dom = attrib_data.domain() ;
    return dom ;
    
  }

  //**************************************************************************/
  
  entitySet dMapRepI::image(const entitySet &iset) const 
  {
    
    entitySet codomain ;
    entitySet :: const_iterator  ei;
    entitySet dom = attrib_data.domain() & iset ;
    int sz = dom.size() ;
    vector<int> mlist(sz) ;
    int i=0 ;
    for( ei = dom.begin(); ei != dom.end(); ++ei){
      mlist[i++] =  attrib_data.elem(*ei) ;
    }
    codomain = create_entitySet(mlist.begin(),mlist.end()) ;
    return codomain ;
  }
  
  //**************************************************************************/
  
  pair<entitySet,entitySet>
  dMapRepI::preimage(const entitySet &codomain) const  
  {
    entitySet domain ;
    entitySet dom = attrib_data.domain() ;
    entitySet::const_iterator  ei;
    for(ei = dom.begin(); ei != dom.end(); ++ei)
      if(codomain.inSet( attrib_data[*ei] ) )
        domain  += *ei ;
    return make_pair(domain,domain);
  }
  
  //**************************************************************************/
  
  storeRepP dMapRepI::get_map() 
  {
    multiMap result ;
    result.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(result.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    store<int> sizes ;
    entitySet storeDomain = attrib_data.domain() ;

    sizes.allocate(storeDomain) ;
    FORALL(storeDomain,i) {
      sizes[i] = 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(storeDomain,i) {
      result.begin(i)[0] = attrib_data[i] ;
    } ENDFORALL ;
    return result.Rep() ;
  }

  //**************************************************************************/

  ostream &dMapRepI::Print(ostream &s) const 
  {
    entitySet dom = domain() ;
    s << '{' << dom << endl ;

    FORALL(dom,ii) {
      s << attrib_data.elem(ii) << endl ;
    } ENDFORALL ;

    s << '}' << endl ;
    return s ;
  }

  //**************************************************************************/

  istream &dMapRepI::Input(istream &s) 
  {

    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    allocate(e) ;

    FORALL(e,ii) {
      s >> attrib_data[ii] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      cerr << "Incorrect Format while reading store" << endl ;
      s.putback(ch) ;
    }

    return s ;
  }

  DatatypeP dMapRepI::getType() {
    return DatatypeP(new AtomicType(INT)) ;
  }
  frame_info dMapRepI::get_frame_info() {
    warn(true) ;
    cerr << "get frame info not implemented for dMapRepI" << endl ;
    frame_info fi ;
    return fi ;
  }

  //**************************************************************************/

  void dMapRepI::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &usr_eset)  {
    warn(true) ;
  }
  
#ifdef H5_HAVE_PARALLEL 
  void dMapRepI::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &usr_eset, hid_t xfer_plist_id)  {
    warn(true) ;
  }
#endif
  
  //**************************************************************************/

  void dMapRepI::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& usr_eset) const
  {
    warn(true) ;
  } 

#ifdef H5_HAVE_PARALLEL 
  void dMapRepI::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& usr_eset, hid_t xfer_plist_id) const
  {
    warn(true) ;
  }
#endif
  
  //**************************************************************************/

  dMap::~dMap() {}

  //**************************************************************************/

  void dMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }    

  //**************************************************************************/

  const_dMap::~const_dMap() {}

  //**************************************************************************/

  void const_dMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      attrib_data = p->get_attrib_data() ;
    warn(p==0) ;
  }

  //**************************************************************************/

  store_instance::instance_type const_dMap::access() const
  { return READ_ONLY ; }

  //****************************************************************************/
  void inverseMap(multiMap &result, const dMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) {
    store<int> sizes ;
    sizes.allocate(input_image) ;
    FORALL(input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    entitySet preloop = input_preimage & input_map.domain() ;
    FORALL(preloop,i) {
      if(input_image.inSet(input_map[i]))
        sizes[input_map[i]] += 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(preloop,i) {
      int elem = input_map[i] ;
      if(input_image.inSet(elem)) {
        sizes[elem] -= 1 ;
        FATAL(sizes[elem] < 0) ;
        result[elem][sizes[elem]] = i ;
      }
    } ENDFORALL ;
#ifdef DEBUG
    FORALL(input_image,i) {
      FATAL(sizes[i] != 0) ;
    } ENDFORALL ;
#endif
  }
}
