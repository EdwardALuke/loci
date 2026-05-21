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
#include <Map.h>
#include <multiMap.h>
#include <DMultiMap.h>
#include <fact_db.h>
#include <Tools/hash_map.h>
#include <iostream>

namespace Loci {
  extern fact_db *exec_current_fact_db ;

  using std::pair ;
  using std::make_pair ;
  
  MapRep::~MapRep() {}
  
  store_type MapRep::RepType() const { return MAP ; }
  
  void MapRepI::allocate(const entitySet &ptn) {
    if(alloc_id < 0)
      alloc_id = getStoreAllocateID() ;

    storeAllocateData[alloc_id].template allocBasic<Entity>(ptn,1) ;
    store_domain = storeAllocateData[alloc_id].allocset ;
    base_ptr = (((Entity *) storeAllocateData[alloc_id].base_ptr) -
		storeAllocateData[alloc_id].base_offset) ;
      
    dispatch_notify() ;
    return ;
  }

  MapRepI::~MapRepI() {
    if(alloc_id>=0) {
      storeAllocateData[alloc_id].template release<Entity>() ;
      releaseStoreAllocateID(alloc_id) ;
      alloc_id = -1 ;
    }
  }

  storeRep *MapRepI::new_store(const entitySet &p) const {
    return new MapRepI(p)  ;
  }
  storeRep *MapRepI::new_store(const entitySet &p, const int* count) const {
    storeRep* sp = 0 ;
    cerr << " This method should not be called for a Map " << endl ;
    return sp ;
  }
  storeRepP MapRepI::MapRemap(const dMap &dm, const dMap &rm) const {
    entitySet newdomain = dm.domain() & domain() ;
    pair<entitySet,entitySet> mappimage = preimage(rm.domain()) ;
    newdomain &= mappimage.first ;
    entitySet mapimage = dm.image(newdomain) ;
    Map s ;
    s.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(s.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(dm,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(rm,mapimage) ;
    return s.Rep() ;
  }

  storeRepP MapRepI::remap(const dMap &m) const {
    cerr << "Map shouldn't use remap!" << endl ;
    return MapRemap(m,m) ;
  }

  void MapRepI::compose(const dMap &m, const entitySet &context) {
    fatal((context-store_domain) != EMPTY) ;
    fatal((image(context)-m.domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = m[base_ptr[i]] ;
    } ENDFORALL ;
  }

  void MapRepI::copy(storeRepP &st, const entitySet &context) {
    const_Map s(st) ;
    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[i] ;
    } ENDFORALL ;
  }

  void MapRepI::gather(const dMap &m, storeRepP &st, const entitySet &context) {
    const_Map s(st) ;
    fatal(base_ptr == 0 && context != EMPTY) ;
    fatal((m.image(context) - s.domain()) != EMPTY) ; 
    fatal((context - domain()) != EMPTY) ;
    FORALL(context,i) {
      base_ptr[i] = s[m[i]] ;
    } ENDFORALL ;
  }

  void MapRepI::scatter(const dMap &m,storeRepP &st, const entitySet &context) {
    const_Map s(st) ;
    fatal(base_ptr == 0 && context != EMPTY) ;
    fatal((context - s.domain()) != EMPTY) ;
    fatal((m.image(context) - domain()) != EMPTY) ;
    fatal((context - m.domain()) != EMPTY);

    FORALL(context,i) {
      base_ptr[m[i]] = s[i] ;
    } ENDFORALL ;
  }
 
  int MapRepI::pack_size(const entitySet &e) {
    fatal((e - domain()) != EMPTY);
    return cpypacksize(base_ptr,e.size());
  }
  int MapRepI::estimated_pack_size(const entitySet &e) {
    return cpypacksize(base_ptr,e.size());
  }
  int MapRepI::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;
    int size = sizeof(Entity) * packed.size() ;
    return size ;
  }
  
  void MapRepI::pack(void *outbuf, int &position, int &outcount, const entitySet &eset) 
  {
    for( size_t i = 0; i < eset.num_intervals(); i++) {
      const Loci::int_type begin = eset[i].first ;
      int t = eset[i].second - eset[i].first + 1 ;

      cpypack(outbuf,position,outcount,&base_ptr[begin],t) ;
    }
  }
  
  void MapRepI::pack(void *outbuf, int &position,
                     int &outcount, const entitySet &eset, const Map& remap) 
  {
    for( size_t i = 0; i < eset.num_intervals(); i++) {
      const Loci::int_type begin = eset[i].first ;
      int t = eset[i].second - eset[i].first + 1 ;
      int* img = new int[t] ;
      for(int k=0;k<t;++k)
        img[k] = remap[base_ptr[begin+k]] ;
      cpypack(outbuf,position,outcount,&img[0],t) ;
      delete[] img ;
    }
  }
  
  void MapRepI::unpack(void *inbuf, int &position, int &insize, const sequence &seq) {

    for(size_t i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
        const Loci::int_type stop = seq[i].second ;
        for(Loci::int_type indx = seq[i].first; indx != stop-1; --indx) {
          cpyunpack(inbuf,position,insize,&base_ptr[indx],1) ;
        }
      } else {
        Loci::int_type indx = seq[i].first ;
        int t = seq[i].second - seq[i].first + 1 ;
        cpyunpack(inbuf,position,insize,&base_ptr[indx],t) ;
      }
    }
  }

  void MapRepI::unpack(void *inbuf, int &position,
                       int &insize, const sequence &seq, const dMap& remap) {

    for(size_t i = 0; i < seq.num_intervals(); ++i) {
      if(seq[i].first > seq[i].second) {
        const Loci::int_type stop = seq[i].second ;
        for(Loci::int_type indx = seq[i].first; indx != stop-1; --indx) {
          cpyunpack(inbuf,position,insize,&base_ptr[indx],1) ;
        }
        // remap
        for(Loci::int_type indx=seq[i].first;indx!=stop-1;--indx)
          base_ptr[indx] = remap[base_ptr[indx]] ;
      } else {
        Loci::int_type indx = seq[i].first ;
        int t = seq[i].second - seq[i].first + 1 ;
        cpyunpack(inbuf,position,insize,&base_ptr[indx],t) ;
        // remap
        for(int k=0;k<t;++k)
          base_ptr[indx+k] = remap[base_ptr[indx+k]] ;
      }
    }
  }

  entitySet MapRepI::domain() const {
    return store_domain ;
  }

  // The integer version of image_section.

  entitySet image_section(const int *start, const int *end) {
    if(start == end)
      return EMPTY ;
    int mx,mn ;
    mx = *start ;
    mn = *start ;
    for(const int *i=start;i!=end;++i) {
      mx = max(mx,*i) ;
      mn = min(mn,*i) ;
    }                                         
    size_t sz = mx-mn+1 ;
    size_t sz2 = end-start ;
    if(sz>2*sz2) {
      // If the image is too sparse then we are better off sorting
      // using standard sort
      entitySet dom ;
      std::vector<int> img(sz2) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(const int *i=start;i!=end;++i) {
        *ins = *i ;
        ++ins ;
      }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        dom += *ins ;
      return dom ;
    }
    std::vector<bool> bits(sz) ;
    for(size_t i=0;i<sz;++i)
      bits[i] = false ;
    for(const int *i=start;i!=end;++i)
      bits[*i-mn] = true ;

    WARN(!bits[0]);

    entitySet result ;
    interval iv(mn,mn) ;
    for(size_t i=0;i<sz;++i)
      if(!bits[i]) {
        iv.second = i+mn-1 ;
        result += iv ;
        for(;i<sz;++i)
          if(bits[i])
            break ;
        iv.first = i+mn ;
      }

    WARN(!bits[sz-1]) ;
    
    iv.second = mx ;
    result += iv ;
    return result ;
  }

#ifdef ENTITY

    // Compare entities using their id numbers.

    Entity max(const Entity e1, const Entity e2){
	int i1 = e1.getid();
	int i2 = e2.getid();
	if (i1 >= i2) {
	    return e1;
	} else {
	    return e2;
	}
    }

    Entity min(const Entity e1, const Entity e2){
	int i1 = e1.getid();
	int i2 = e2.getid();
	if (i1 <= i2) {
	    return e1;
	} else {
	    return e2;
	}
    }

#endif

#ifdef ENTITY

  // The Entity version of image_section.

    entitySet image_section(const Entity *start, const Entity *end) {
    if(start == end)
      return EMPTY ;
    Entity mx,mn ;
    mx = *start ;
    mn = *start ;
    for(const Entity *i=start;i!=end;++i) {
      mx = max(mx,*i) ;
      mn = min(mn,*i) ;
    }        
    int mxi = getEntityIdentity(mx);
    int mni = getEntityIdentity(mn);
    size_t sz = mxi-mni+1 ;
    int si = getEntityIdentity(*start);
    int ei = getEntityIdentity(*end);
    size_t sz2 = ei-si ;
    if(sz>2*sz2) {
      // If the image is too sparse then we are better off sorting
      // using standard sort
      entitySet dom ;
      std::vector<int> img(sz2) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(const Entity *i=start;i!=end;++i) {
        *ins = *i ;
        ++ins ;
      }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        dom += *ins ;
      return dom ;
    }
    std::vector<bool> bits(sz) ;
    for(int i=0;i<sz;++i)
      bits[i] = false ;
    for(const int *i=&si;i!=&ei;++i)
      bits[*i-mn] = true ;

    WARN(!bits[0]);

    entitySet result ;
    interval iv(mn,mn) ;
    for(int i=0;i<sz;++i)
      if(!bits[i]) {
        iv.second = i+mn-1 ;
        result += iv ;
        for(;i<sz;++i)
          if(bits[i])
            break ;
        iv.first = i+mn ;
      }

    WARN(!bits[sz-1]) ;
    
    iv.second = mx ;
    result += iv ;
    return result ;
    }      

#endif

    
  entitySet MapRepI::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    if(d.num_intervals() < IMAGE_THRESHOLD) {
      for(size_t i=0;i<d.num_intervals();++i)
        codomain += image_section(base_ptr+d[i].first,
                                  base_ptr+(d[i].second+1)) ;
    } else {
      std::vector<int> img(d.size()) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(size_t i=0;i<d.num_intervals();++i)
        for(int j=d[i].first;j!=d[i].second+1;++j) {
          *ins = base_ptr[j] ;
          ++ins ;
        }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        codomain += *ins ;
    }
      
    return codomain ;
  }

  pair<entitySet,entitySet>
  MapRepI::preimage(const entitySet &codomain) const  {
    entitySet domain ;
    FORALL(store_domain,i) {
      if(codomain.inSet(base_ptr[i]))
        domain += i ;
    } ENDFORALL ;
    return make_pair(domain,domain) ;
  }
  
  storeRepP MapRepI::expand(entitySet &out_of_dom, std::vector<entitySet> &ptn) {
    storeRepP sp ;
    MPI_Comm comm = get_exec_comm() ;
    const int np = get_exec_size() ;
    int *recv_count = new int[np] ;
    int *send_count = new int[np] ;
    int *send_displacement = new int[np] ;
    int *recv_displacement = new int[np] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy(np), send_clone(np) ;
    for(int i = 0; i < np; ++i) {
      entitySet tmp = out_of_dom & ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	copy[i].push_back(*ei) ;
      std::sort(copy[i].begin(), copy[i].end()) ;
      send_count[i] = copy[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i < np; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < np; ++i)
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < np; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    for(int i = 0; i < np; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	send_clone[i].push_back(recv_buf[j]) ;
      std::sort(send_clone[i].begin(), send_clone[i].end()) ;
    }
    
    std::vector<HASH_MAP(int, int) > map_entities(np) ;
    for(int i = 0; i < np; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) 
	if(store_domain.inSet(*vi))
	  (map_entities[i])[*vi] = base_ptr[*vi] ;
    
    size_send = 0 ;
    for(int i = 0; i < np; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      size_send += send_count[i] ;
    }
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i < np; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < np; ++i) 
      for(HASH_MAP(int, int)::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < np; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    HASH_MAP(int, int) hm ;
    for(int i = 0; i < np; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; j+=2) {
	hm[recv_map[j]] = recv_map[j+1];
      }
    }
    Map dm ;
    dm.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(dm.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    
    entitySet tm = store_domain + out_of_dom ;
    dm.allocate(tm) ;
    for(ei = store_domain.begin(); ei != store_domain.end(); ++ei)
      dm[*ei] = base_ptr[*ei] ;
    for(HASH_MAP(int, int)::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi) 
      dm[hmi->first] = hmi->second ;
    sp = dm.Rep() ;
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

  storeRepP MapRepI::freeze() {
    return getRep() ;
  }
  
  storeRepP MapRepI::thaw() {
    
    dMap dm ;
    dm.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(dm.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    FORALL(store_domain,i) {
      dm[i] = base_ptr[i] ;
    } ENDFORALL ;
    return dm.Rep() ;
  }
  storeRepP MapRepI::get_map() {
    store<int> sizes ;
    sizes.allocate(store_domain) ;
    FORALL(store_domain,i) {
      sizes[i] = 1 ;
    } ENDFORALL ;
    multiMap result ;
    result.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(result.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    result.allocate(sizes) ;
    FORALL(store_domain,i) {
      result.begin(i)[0] = base_ptr[i] ;
    } ENDFORALL ;
    return result.Rep() ;
  }
    
  std::ostream &MapRepI::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      s << base_ptr[ii] << std::endl ;
    }ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }


  std::istream &MapRepI::Input(std::istream &s) {
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
      s >> base_ptr[ii] ;
    } ENDFORALL ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }
  DatatypeP MapRepI::getType() {
    return DatatypeP(new AtomicType(INT)) ;
  }

  frame_info MapRepI::get_frame_info() {
    warn(true) ;
    frame_info fi ;
    return fi ;
   }
  
  void MapRepI::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &usr_eset){
    warn(true) ; 
  } 

#ifdef H5_HAVE_PARALLEL 
  void MapRepI::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &usr_eset, hid_t xfer_plist_id){
    warn(true) ; 
  } 
#endif
  void MapRepI::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &usr_eset) const{
    warn(true) ;
  } 

#ifdef H5_HAVE_PARALLEL 
  void MapRepI::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet &usr_eset, hid_t xfer_plist_id) const{
    warn(true) ;
  } 
#endif  
  Map::~Map() {}

  void Map::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }    

  const_Map::~const_Map() {}

  void const_Map::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  store_instance::instance_type const_Map::access() const
  { return READ_ONLY ; }
    
  void multiMapRepI::allocate(const entitySet &ptn) {
    store<int> count ;
    count.allocate(ptn) ;
    FORALL(ptn,i) {
      count[i] = 0 ;
    } ENDFORALL ;
    allocate(count) ;
  }

  storeRepP multiMapRepI::expand(entitySet &out_of_dom, std::vector<entitySet> &ptn) {
    MPI_Comm comm = get_exec_comm() ;
    const int np = get_exec_size() ;
    int *recv_count = new int[np] ;
    int *send_count = new int[np] ;
    int *send_displacement = new int[np] ;
    int *recv_displacement = new int[np] ;
    entitySet::const_iterator ei ;
    std::vector<int>::const_iterator vi ;
    int size_send = 0 ;
    std::vector<std::vector<int> > copy(np), send_clone(np) ;
    for(int i = 0; i < np; ++i) {
      entitySet tmp = out_of_dom & ptn[i] ;
      for(ei = tmp.begin(); ei != tmp.end(); ++ei)
	copy[i].push_back(*ei) ;
      std::sort(copy[i].begin(), copy[i].end()) ;
      send_count[i] = copy[i].size() ;
      size_send += send_count[i] ; 
    }
    int *send_buf = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i < np; ++i)
      size_send += recv_count[i] ;
    
    int *recv_buf = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < np; ++i)
      for(vi = copy[i].begin(); vi != copy[i].end(); ++vi) {
	send_buf[size_send] = *vi ;
	++size_send ;
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < np; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_buf,send_count, send_displacement , MPI_INT,
		  recv_buf, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    for(int i = 0; i < np; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]; ++j) 
	send_clone[i].push_back(recv_buf[j]) ;
      std::sort(send_clone[i].begin(), send_clone[i].end()) ;
    }
    std::vector<HASH_MAP(int, std::vector<int>) > map_entities(np) ;
    for(int i = 0; i < np; ++i) 
      for(vi = send_clone[i].begin(); vi != send_clone[i].end(); ++vi) {
	if(store_domain.inSet(*vi)) {
	  for(const int* j = begin(*vi); j != end(*vi); ++j)
	    (map_entities[i])[*vi].push_back(*j) ;
	}
      }
    
    for(int i = 0; i < np; ++i) {
      send_count[i] = 2 * map_entities[i].size() ;
      for(HASH_MAP(int, std::vector<int>)::iterator hi = map_entities[i].begin(); hi != map_entities[i].end(); ++hi)
	send_count[i] += hi->second.size() ; 
    }
    size_send = 0 ;
    for(int i = 0; i < np; ++i)
      size_send += send_count[i] ;
    int *send_map = new int[size_send] ;
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT,
		 comm) ; 
    size_send = 0 ;
    for(int i = 0; i < np; ++i)
      size_send += recv_count[i] ;
    int *recv_map = new int[size_send] ;
    size_send = 0 ;
    for(int i = 0; i < np; ++i) 
      for(HASH_MAP(int, std::vector<int> )::const_iterator miv = map_entities[i].begin(); miv != map_entities[i].end(); ++miv) {
	send_map[size_send] = miv->first ;
	++size_send ;
	send_map[size_send] = miv->second.size() ;
	++size_send ;
	for(std::vector<int>::const_iterator vi = miv->second.begin(); vi != miv->second.end(); ++vi) { 
	  send_map[size_send] = *vi ;
	  ++size_send ;
	}
      }
    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i < np; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_count[i-1] ;
    }
    MPI_Alltoallv(send_map,send_count, send_displacement , MPI_INT,
		  recv_map, recv_count, recv_displacement, MPI_INT,
		  comm) ;  
    HASH_MAP(int, std::vector<int> ) hm ;
    std::vector<int> ss ;
    for(int i = 0; i < np; ++i) {
      for(int j = recv_displacement[i]; j <
	    recv_displacement[i]+recv_count[i]-1; ++j) {
	int count = recv_map[j+1] ;
	if(count)
	  for(int k = 0; k < count; ++k)
	    hm[recv_map[j]].push_back(recv_map[j+k+2]);
	else
	  hm[recv_map[j]] = ss ;
	j += count + 1 ;
      }
    }
    multiMap dmul ;
    dmul.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(dmul.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    entitySet tm = store_domain + out_of_dom ;
    store<int> count ;
    count.allocate(tm) ;
    for(ei = store_domain.begin(); ei != store_domain.end(); ++ei)
      count[*ei] = end(*ei) - begin(*ei) ;
    for(HASH_MAP(int, std::vector<int> )::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi)
      count[hmi->first] = hmi->second.size() ;
    dmul.allocate(count) ;
    
    for(ei = store_domain.begin(); ei != store_domain.end(); ++ei)
      for(int p = 0; p < count[*ei] ; ++p)
	dmul[*ei][p] = base_ptr[*ei][p] ;
    
    for(HASH_MAP(int, std::vector<int> )::const_iterator hmi = hm.begin(); hmi != hm.end(); ++hmi) {
      int p = 0 ;
      for(std::vector<int>::const_iterator si = hmi->second.begin(); si != hmi->second.end(); ++si) {
	dmul[hmi->first][p] = *si ;
	p++ ;
      }
    }
    storeRepP sp = dmul.Rep() ;
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

  storeRepP multiMapRepI::freeze() {
    return getRep() ;
  }
  
  storeRepP multiMapRepI::thaw() {
    dmultiMap dm ;
    dm.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(dm.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    FORALL(store_domain, i) {
      size_t sz = end(i)-begin(i) ;
      std::vector<int,malloc_alloc<int> >(sz).swap(dm[i]) ;
      for(size_t j = 0; j < sz; ++j)
        dm[i][j] = base_ptr[i][j] ;
    } ENDFORALL ;
    return(dm.Rep()) ;
  }
  

  void multiMapRepI::allocate(const store<int> &sizes) {
    if(alloc_id < 0)
      alloc_id = getStoreAllocateID() ;

    entitySet ptn = sizes.domain() ;
    int cntid = sizes.Rep()->get_alloc_id() ;
    storeAllocateData[alloc_id].template release<Entity>() ;
    storeAllocateData[alloc_id].template
      allocMulti<Entity>(storeAllocateData[cntid],ptn) ;

    base_ptr = ((Entity **)storeAllocateData[alloc_id].alloc_ptr2 -
		storeAllocateData[alloc_id].base_offset) ;
   
    store_domain = ptn ;
    dispatch_notify() ;
  }
  

  multiMapRepI::~multiMapRepI() {
    if(alloc_id>=0) {
      storeAllocateData[alloc_id].template release<Entity>() ;
      releaseStoreAllocateID(alloc_id) ;
      alloc_id = -1 ;
    }
  }

  storeRep *multiMapRepI::new_store(const entitySet &p) const {
    return new multiMapRepI()  ;
  }
  storeRep *multiMapRepI::new_store(const entitySet &p, const int* cnt) const {
    store<int> count ;
    count.allocate(p) ;
    int t= 0 ;
    FORALL(p, pi) {
      count[pi] = cnt[t++] ; 
    } ENDFORALL ;
    return new multiMapRepI(count)  ;
  }
  storeRepP multiMapRepI::MapRemap(const dMap &dm, const dMap &rm) const {
    entitySet newdomain = dm.domain() & domain() ;
    //    pair<entitySet,entitySet> mappimage = preimage(m.domain()) ;
    //    newdomain &= mappimage.first ;
    entitySet mapimage = dm.image(newdomain) ;
    multiMap s ;
    s.Rep()->setDomainKeySpace(getDomainKeySpace()) ;
    MapRepP(s.Rep())->setRangeKeySpace(getRangeKeySpace()) ;
    
    s.allocate(mapimage) ;
    storeRepP my_store = getRep() ;
    s.Rep()->scatter(dm,my_store,newdomain) ;
    MapRepP(s.Rep())->compose(rm,mapimage) ;
    return s.Rep() ;
  } 
  storeRepP multiMapRepI::remap(const dMap &m) const {
    cerr << "remap should not be used on Map!" << endl ;
    return MapRemap(m,m) ;
  }

  void multiMapRepI::compose(const dMap &m, const entitySet &context) {
    fatal((context-store_domain) != EMPTY) ;
    //fatal((image(context)-m.domain()) != EMPTY) ;
    entitySet dom = m.domain() ;
    FORALL(context,i) {
      for(int *ii = base_ptr[i];ii!=base_ptr[i+1];++ii) {
        if(dom.inSet(*ii))
          *ii = m[*ii] ;
        else
          *ii = -1 ;
      }
    } ENDFORALL ;
  }

    
  void multiMapRepI::copy(storeRepP &st, const entitySet &context) {
    const_multiMap s(st) ;
    fatal((context-domain()) != EMPTY) ;
    fatal((context-s.domain()) != EMPTY) ;

    store<int> count ;
    count.allocate(domain()) ;
    FORALL(domain()-context,i) {
      count[i] = base_ptr[i+1]-base_ptr[i] ;
    } ENDFORALL ;
    FORALL(context,i) {
      count[i] = s.end(i)-s.begin(i) ;
    } ENDFORALL ;

    // Allocate a temporary copy for the copy
    int cntid = count.Rep()->get_alloc_id() ;
    storeAllocateInfo tmp ;
    tmp.template allocMulti<Entity>(storeAllocateData[cntid],count.domain()) ;

    Entity **new_base_ptr = ((Entity **)tmp.alloc_ptr2 -
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
    storeAllocateData[alloc_id].template release<Entity>() ;
    // copy newly copied temporary to main allocatin
    storeAllocateData[alloc_id] = tmp ;
    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }

  void multiMapRepI::gather(const dMap &m, storeRepP &st,
                            const entitySet  &context) {
    store<int> count ;
    const_multiMap s(st) ;
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
    tmp.template allocMulti<Entity>(storeAllocateData[cntid],count.domain()) ;

    Entity **new_base_ptr = ((Entity **)tmp.alloc_ptr2 -
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
    storeAllocateData[alloc_id].template release<Entity>() ;
    // copy newly copied temporary to main allocatin
    storeAllocateData[alloc_id] = tmp ;

    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }
 
  void multiMapRepI::scatter(const dMap &m, storeRepP &st,
                             const entitySet  &context) {
    store<int> count ;
    const_multiMap s(st) ;
    count.allocate(domain()) ;

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
    tmp.template allocMulti<Entity>(storeAllocateData[cntid],count.domain()) ;

    Entity **new_base_ptr = ((Entity **)tmp.alloc_ptr2 -
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
    storeAllocateData[alloc_id].template release<Entity>() ;
    // copy newly copied temporary to main allocatin
    storeAllocateData[alloc_id] = tmp ;

    base_ptr = new_base_ptr ;
    dispatch_notify() ;
  }
  
  int multiMapRepI::pack_size(const  entitySet &eset ) {
    fatal((eset - domain()) != EMPTY);

    int size = 0 ;
    FORALL(eset,i) {
      int cnt = end(i) - begin(i) ;
      size += cpypacksize(&cnt,1) ;
      size += cpypacksize(begin(i),cnt) ;
    } ENDFORALL ;
    
    return size ;
  }
  
  int multiMapRepI::estimated_pack_size(const  entitySet &eset ) {
    return 5*eset.size()*sizeof(Entity);
  }

  int multiMapRepI::
  pack_size(const entitySet& e, entitySet& packed) {
    packed = domain() & e ;
    int size = 0 ;
    FORALL(packed, i) {
      int cnt = end(i) - begin(i) ;
      size += cpypacksize(&cnt,1) ;
      size += cpypacksize(begin(i),cnt) ;
    } ENDFORALL ;

    return size ;
  }

  void multiMapRepI::pack(void *outbuf, int &position, int &outcount, const entitySet &eset) {


    entitySet :: const_iterator ci;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      int vsize    = end(*ci) - begin(*ci);
      cpypack(outbuf,position,outcount,&vsize,1) ;
      cpypack(outbuf,position,outcount,begin(*ci),vsize) ;
    }

  }
  
  void multiMapRepI::pack(void *outbuf, int &position,
                          int &outcount, const entitySet &eset,
                          const Map& remap) {
    int vsize;
    entitySet :: const_iterator ci;
    for( ci = eset.begin(); ci != eset.end(); ++ci) {
      vsize    = end(*ci) - begin(*ci);
      cpypack(outbuf,position,outcount,&vsize,1) ;
      // remap first
      int* img = new int[vsize] ;
      for(int k=0;k<vsize;++k)
        img[k] = remap[base_ptr[*ci][k]] ;
      cpypack(outbuf,position,outcount,&img[0],vsize) ;
      delete[] img ;
    }
  }
  
  void multiMapRepI::unpack(void *inbuf, int &position, int &insize, const sequence &seq) {

    int vsize;
    sequence:: const_iterator ci;
    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      cpyunpack(inbuf,position,insize,&vsize,1) ;
#ifdef DEBUG
      if(vsize != end(*ci)-begin(*ci))
        cerr << "vsize = " << vsize << ",actual = " << end(*ci)-begin(*ci)
             << ",ci = " << *ci << endl ;
#endif
      fatal(vsize != end(*ci)-begin(*ci)) ;
      cpyunpack(inbuf,position,insize,begin(*ci),vsize) ;
    }
  }   
    
  void multiMapRepI::unpack(void *inbuf, int &position,
                            int &insize, const sequence &seq,
                            const dMap& remap) {
    int vsize;
    sequence:: const_iterator ci;
    for( ci = seq.begin(); ci != seq.end(); ++ci) {
      cpyunpack(inbuf,position,insize,&vsize,1) ;
#ifdef DEBUG
      if(vsize != end(*ci)-begin(*ci))
        cerr << "vsize = " << vsize << ",actual = " << end(*ci)-begin(*ci)
             << ",ci = " << *ci << endl ;
#endif
      fatal(vsize != end(*ci)-begin(*ci)) ;
      cpyunpack(inbuf,position,insize,begin(*ci),vsize) ;
      // remap
      for(int k=0;k<vsize;++k)
        base_ptr[*ci][k] = remap[base_ptr[*ci][k]] ;
    }
  }   
    
  entitySet multiMapRepI::domain() const {
    return store_domain ;
  }
    
  entitySet multiMapRepI::image(const entitySet &domain) const {
    entitySet d = domain & store_domain ;
    entitySet codomain ;
    if(d.num_intervals() < IMAGE_THRESHOLD) {
      for(size_t i=0;i<d.num_intervals();++i)
        codomain += image_section(begin(d[i].first),end(d[i].second)) ;
    } else {
      size_t sz = 0 ;
      for(size_t i=0;i<d.num_intervals();++i)
        sz += end(d[i].second)-begin(d[i].first) ;

      std::vector<int> img(sz) ;
      std::vector<int>::iterator ins = img.begin() ;
      for(size_t i=0;i<d.num_intervals();++i)
        for(const int *j=begin(d[i].first);j!=end(d[i].second);++j) {
          *ins = *j ;
          ++ins ;
        }
      std::sort(img.begin(),img.end()) ;
      std::vector<int>::iterator uend = std::unique(img.begin(),img.end());
      for(ins=img.begin();ins!=uend;++ins)
        codomain += *ins ;
    }
    return codomain ;
  }

  pair<entitySet,entitySet>
  multiMapRepI::preimage(const entitySet &codomain) const  {
    entitySet domaini,domainu ;
    FORALL(store_domain,i) {
      bool vali = true ;
      bool valu = false;
      for(const int *ip = begin(i);ip!= end(i);++ip) {
        bool in_set = codomain.inSet(*ip) ;
        vali = vali && in_set ;
        valu = valu || in_set ;
      }
      if(vali)
        domaini += i ;
      if(valu)
        domainu += i ;
    } ENDFORALL ;
    return make_pair(domaini,domainu) ;
  }

  storeRepP multiMapRepI::get_map() {
    return this ;
  }
    
  std::ostream &multiMapRepI::Print(std::ostream &s) const {
    s << '{' << domain() << std::endl ;
    FORALL(domain(),ii) {
      s << end(ii)-begin(ii) << std::endl ;
    } ENDFORALL ;
    FORALL(domain(),ii) {
      for(const int *ip = begin(ii);ip!=end(ii);++ip)
        s << *ip << " " ;
      s << std::endl;
    } ENDFORALL ;
    s << '}' << std::endl ;
    return s ;
  }


  std::istream &multiMapRepI::Input(std::istream &s) {
    entitySet e ;
    char ch ;
    
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '{') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> e ;
    store<int> sizes ;
    sizes.allocate(e) ;
    FORALL(e,ii) {
      s >> sizes[ii] ;
    } ENDFORALL ;

    allocate(sizes) ;
        
    FORALL(e,ii) {
      for(int *ip = begin(ii);ip!=end(ii);++ip)
        s >> *ip  ;
    } ENDFORALL ;
            
    do ch = s.get(); while(ch==' ' || ch=='\n') ;
    if(ch != '}') {
      std::cerr << "Incorrect Format while reading store" << std::endl ;
      s.putback(ch) ;
    }
    return s ;
  }

  DatatypeP multiMapRepI::getType() {
    warn(true) ;
    DatatypeP dp ;
    return dp ;
  }
  frame_info multiMapRepI::get_frame_info() {
    warn(true) ;
    frame_info fi ;
    return fi ;
  }
  
  void multiMapRepI::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &usr_eset) {
    warn(true) ;
  }

#ifdef H5_HAVE_PARALLEL 
  void multiMapRepI::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &usr_eset, hid_t xfer_plist_id) {
    warn(true) ;
  }
#endif
  
  void multiMapRepI::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& eset) const{
    warn(true) ;
  } 
  
#ifdef H5_HAVE_PARALLEL 
  void multiMapRepI::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& eset, hid_t xfer_plist_id) const{
    warn(true) ;
  }  
#endif
  
  multiMap::~multiMap() {}

  void multiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  const_multiMap::~const_multiMap() { }

  void const_multiMap::notification() {
    NPTR<MapType> p(Rep()) ;
    if(p!=0)
      base_ptr = p->get_base_ptr() ;
    warn(p==0) ;
  }

  store_instance::instance_type const_multiMap::access() const
  { return READ_ONLY ; }
    
  void inverseMap(multiMap &result, const Map &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) {
    store<int> sizes ;
    sizes.allocate(input_image) ;
    result.Rep()->setDomainKeySpace(MapRepP(input_map.Rep())->getRangeKeySpace()) ;
    MapRepP(result.Rep())->setRangeKeySpace(input_map.Rep()->getDomainKeySpace()) ;

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
  

  void inverseMap(multiMap &result, const multiMap &input_map,
                  const entitySet &input_image,
                  const entitySet &input_preimage) {
    store<int> sizes ;
    sizes.allocate(input_image) ;
    result.Rep()->setDomainKeySpace(MapRepP(input_map.Rep())->getRangeKeySpace()) ;
    MapRepP(result.Rep())->setRangeKeySpace(input_map.Rep()->getDomainKeySpace()) ;
    
    FORALL(input_image,i) {
      sizes[i] = 0 ;
    } ENDFORALL ;
    entitySet preloop = input_preimage & input_map.domain() ;

    FORALL(preloop,i) {
      for(const int *mi=input_map.begin(i);mi!=input_map.end(i);++mi)
        if(input_image.inSet(*mi))
          sizes[*mi] += 1 ;
    } ENDFORALL ;
    result.allocate(sizes) ;
    FORALL(preloop,i) {
      for(const int *mi=input_map.begin(i);mi!=input_map.end(i);++mi) {
        int elem = *mi ;
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
   
}
