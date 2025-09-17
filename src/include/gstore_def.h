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
#ifndef GSTORE_DEF_H
#define GSTORE_DEF_H

// This header file contains the class definition of
// GStore
// Their corresponding template implementation is in
// gstore_impl.h
// The same design applies to all other container classes.

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>
#include <data_traits.h>
#include <sstream>
#include <iterator>
#include <vector>
#include <gmap.h>
#include <gstore_rep.h>
#include <hdf5_readwrite_long.h>
#include <mpi.h>
#include <string.h>
#include <dist_internal.h>
//for test
#include <store.h>
namespace Loci {
  extern int MPI_processes;
  extern int MPI_rank ;
  
  
  // template<class T>  inline bool equiFF(const std::pair<gEntity,T> &v1,
  //                                         const std::pair<gEntity,T> &v2) {
  //     return v1.first < v2.first ;
  //   }

  //   template<class T>   inline void equiJoinFF(std::vector<std::pair<gEntity,T> > &in1,
  //                                              std::vector<std::pair<gEntity,gEntity> > &in2,
  //                                              std::vector<std::pair<gEntity,T> > &out) {
  //       std::sort(in1.begin(),in1.end(),equiFF) ;
  //       std::sort(in2.begin(),in2.end(),equiFF) ;
    
  //       int p = 0 ;
  //       MPI_Comm_size(MPI_COMM_WORLD,&p) ;

  //       // Sort inputs using same splitters (this will make sure that
  //       // data that needs to be on the same processor ends up on the
  //       // same processor
  //       if(p != 1) {
  //         std::vector<std::pair<gEntity,gEntity> > splitters ;
  //         parGetSplitters(splitters,in1,equiFF,MPI_COMM_WORLD) ;

  //         parSplitSort(in1,splitters,equiFF,MPI_COMM_WORLD) ;
  //         parSplitSort(in2,splitters,equiFF,MPI_COMM_WORLD) ;
  //       }

  //       // Find pairs where first entry are the same and create joined protomap
  //       out.clear() ;
  //       size_t j = 0 ;
  //       for(size_t i=0;i<in1.size();++i) {
  //         while(j<in2.size() && equiFF(in2[j],in1[i]))
  //           ++j ;
  //         size_t k=j ;
  //         while(k<in2.size() && in2[k].first == in1[i].first) {
  //           out.push_back(std::pair<gEntity,T>(in1[i].second,in2[k].second)) ;
  //           k++ ;
  //         }
  //       }

  //       // Remove duplicates from protomap
  //       parSampleSort(out,equiFF,MPI_COMM_WORLD) ;
  //       std::sort(out.begin(),out.end()) ;
  //       out.erase(std::unique(out.begin(),out.end()),out.end()) ;
  //     }

  

  //where storage is allocated
  template<class T> class gStoreRepI : public gStoreRep {
    //storage type
    typedef typename std::vector<std::pair<gEntity,T> > gRep ;
    
  private:
    bool sorted; 
    gRep attrib_data;
    gKeySpaceP domain_space;
  private:
    int  get_mpi_size( IDENTITY_CONVERTER c, const gEntitySet &eset)const;
    void packdata(IDENTITY_CONVERTER, void *ptr, int &loc, int size,
                  const gEntitySet &e )const ;
    void unpackdata(IDENTITY_CONVERTER c, const void *ptr, int &loc, int size) ;
    int  get_mpi_size( USER_DEFINED_CONVERTER c, const gEntitySet &eset)const;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size,
                  const gEntitySet &e )const ;
    void unpackdata(USER_DEFINED_CONVERTER c, const void *ptr, int &loc, int size) ;
    DatatypeP getType(IDENTITY_CONVERTER g)const{
      typedef data_schema_traits<T> traits_type;
      return(traits_type::get_type()) ;
    }
    DatatypeP getType(USER_DEFINED_CONVERTER g)const{
      typedef data_schema_traits<T> schema_traits ;
      typedef typename schema_traits::Converter_Base_Type dtype;
      typedef data_schema_traits<dtype> traits_type;
      return(traits_type::get_type()) ;
    }
    frame_info get_frame_info(IDENTITY_CONVERTER g)const ;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g)const ;
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                   const char* name, IDENTITY_CONVERTER c, frame_info &fi,  const gEntitySet &en) ;
    void  hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                   const char* name, USER_DEFINED_CONVERTER c, frame_info &fi,  const gEntitySet &en);
    
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                    const char* name, IDENTITY_CONVERTER c, const gEntitySet &en) const;
    void  hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                    const char* name, USER_DEFINED_CONVERTER c, const gEntitySet &en) const;
    virtual gStoreRepP recompose(const gMap &m, MPI_Comm comm=MPI_COMM_WORLD )const ;
    virtual gStoreRepP recompose(const gMultiMap &m, MPI_Comm comm=MPI_COMM_WORLD )const ;
  public:
    typedef typename std::vector<std::pair<gEntity,T> >::iterator iterator ;
    typedef typename std::vector<std::pair<gEntity,T> >::const_iterator const_iterator ;
    iterator begin() { return attrib_data.begin(); }
    iterator end() { return attrib_data.end(); }
    const_iterator begin() const{ return attrib_data.begin(); }
    const_iterator end() const { return attrib_data.end(); }
    size_t size() const{return attrib_data.size();}
    void reserve (size_t n){attrib_data.reserve(n);}
    void clear(){attrib_data.clear();}
    gStoreRepI():sorted(true),domain_space(0){}
    void set_domain_space(gKeySpaceP space){domain_space = space;}
    gKeySpaceP get_domain_space()const{return domain_space;}
    
    //For maps, recompose will remap the SECOND field of this using m and create a new map
    //For stores,recompose will compose a store whose domain is the domain of m,
    //whose data is the data of the SECOND field of m.
    //for example, pos.recompose(face2node) will produce the positions  for each face  
    virtual gStoreRepP  recompose( gStoreRepP &m, MPI_Comm comm=MPI_COMM_WORLD)const  ;
    
    
    // this method redistributes the stores according to the split of local domain over a group of process
    //dom_split: the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const;
    
    // this redistribute version takes an additional remap
    // argument, after redistribution, the new store is remapped
    // dom_split:  the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,const gMap &remap,
                 MPI_Comm comm=MPI_COMM_WORLD)const;
    
    // the redistribute takes a vector of gEntitySets as domain
    // distribution over a group of processes and
    // redistributes the container according to the domain partition
    // NOTE: should be pure virtual (= 0), but here we default it
    // to do nothing just to be able to compile the code, all
    // containers will need to implement this method later.
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const;
    
    //binary search function, for random access
    //mixed equality and equivalance
    //should give warning info if search fails
    iterator search(iterator first,
                    iterator last,
                    const gEntity& key){
      // if(!sorted)throw(bad_access);
      if(first->first==key) return first;
      iterator it;
     
      typename  std::iterator_traits<iterator>::difference_type count, step;
      count = distance(first,last);
      while (count>0)
        {
          it = first; step=count/2; advance(it,step);
          if (it->first < key) {                
            first=++it;
            count-=step+1;
          }
          else count=step;
        }
      return first;
    }
    
    void insert(gEntity e, const T &val){
      sorted = false;
      attrib_data.push_back(std::pair<gEntity, T>(e, val));
    }

    void erase(iterator itr1, iterator itr2){
      attrib_data.erase(itr1, itr2);
    }
    gEntitySet domain()const{
      gEntitySet dom = GEMPTY;
      if(!sorted)const_cast<gStoreRepI&>(*this).local_sort();
      for(const_iterator itr = begin(); itr != end(); itr++){
        dom += itr->first;
      }
      return dom;
    }
    void insert(const gEntitySet& seq,  const T* vals){
      sorted = false;
      size_t idx= 0;
      for(typename gEntitySet::const_iterator itr = seq.begin();
          itr!= seq.end(); itr++){
        attrib_data.push_back(std::pair<gEntity, T>(*itr, vals[idx++]));
      }
    }
    
    void local_sort();
    // void parSamplesort();
    // void sort();
    // void parSplitSort(const std::vector<gEntitySet>& ptn);
    
    //binary search function with const access, for random access
    const_iterator search(const_iterator first,
                          const_iterator last,
                          const gEntity& key)const{
      // if(!sorted)throw(bad_access);
      if(first->first==key) return first;
      const_iterator it;
     
      typename  std::iterator_traits<const_iterator>::difference_type count, step;
      count = distance(first,last);
      while (count>0)
        {
          it = first; step=count/2; advance(it,step);
          if (it->first < key) {                
            first=++it;
            count-=step+1;
          }
          else count=step;
        }
      return first;
    }
    //copy gstore to traditional Loci store, storevec or multistore
    //partially implemented
#ifdef COPY2STORE  
    virtual storeRepP copy2store() const{
      fatal(!sorted);
      entitySet dom = gconvert(domain());
      fatal(size() != dom.size()) ;
      store<T> s;
      s.allocate(dom);
      for(const_iterator itr = begin(); itr != end(); itr++){  
        s[itr->first] = itr->second;
      }
      return s.Rep();
    }
#endif
   
    virtual ~gStoreRepI(){}
    virtual gStoreRepP clone() const{return new gStoreRepI(*this);}
    virtual gStoreRepP remap(const gMap &m) const ;
    virtual int pack_size(const gEntitySet &e)const ;
    virtual void pack(void *ptr, int &loc, int size, const gEntitySet &e)const ;
    virtual void unpack(const void *ptr, int &loc, int size) ;
    virtual gstore_type RepType() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
   
    virtual void* get_attrib_data() {return &attrib_data; }
    virtual const void* get_attrib_data() const{return &attrib_data; }
    virtual DatatypeP getType()const{
      typedef typename data_schema_traits<T>::Schema_Converter schema_converter;
      return getType(schema_converter()) ;}
    virtual frame_info get_frame_info()const ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                          const char* name,  frame_info &fi, const gEntitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                           const char* name, const gEntitySet& en) const ;
   
  } ;
  
  
   
  template<class T> class gStore: public gstore_instance {
    typedef typename std::vector<std::pair<gEntity,T> > gRep ; 
    typedef gStoreRepI<T>  storeType ;
    gStore(const gStore &var) { setRep(var.Rep()) ;}
    gStore<T> & operator=(const gStore<T> &str) {
      setRep(str.Rep()) ;
      return *this ;
    }
  public:
    typedef T containerType ;
    gStore() { setRep(new storeType); }
    gStore(gStoreRepP &rp) { setRep(rp) ;}
    
    virtual ~gStore() {}
   
    gStore<T> & operator=(gStoreRepP p) { setRep(p) ; return *this ; }

    void set_domain_space(gKeySpaceP space){Rep()->set_domain_space(space);}
    gKeySpaceP get_domain_space()const{return Rep()->get_domain_space();}
    
    virtual gStoreRepP clone() const{return Rep()->clone();}
    gEntitySet domain() const { return Rep()->domain(); }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

    //iterators, use std iterators so that std functions can be used
    typedef typename std::vector<std::pair<gEntity,T> >::iterator iterator ;
    typedef typename std::vector<std::pair<gEntity,T> >::const_iterator const_iterator ;

    iterator begin() { return static_cast<gRep*>(Rep()->get_attrib_data())->begin(); }
    iterator end() { return static_cast<gRep*>(Rep()->get_attrib_data())->end(); }
    const_iterator begin()const { return static_cast<const gRep*>(Rep()->get_attrib_data())->begin(); }
    const_iterator end()const { return static_cast<const gRep*>(Rep()->get_attrib_data())->end(); }
    size_t size() const{return Rep()->size();}
    void reserve (size_t n){Rep()->reserve(n);}
    void clear(){Rep()->clear();}

    // the remap method merely renumbers the container
    // according to the passed in map
    virtual gStoreRepP remap(const gMap &m) const{return Rep()->remap(m);}
    
    //For maps, recompose will remap the SECOND field of this using m and create a new map
    //For stores,recompose will compose a store whose domain is the domain of m,
    //whose data is the data of the SECOND field of m.
    //for example, pos.recompose(face2node) will produce the positions  for each face  
    virtual gStoreRepP  recompose(gStoreRepP &m, MPI_Comm comm=MPI_COMM_WORLD)const{
      return Rep()->recompose(m, comm);}
    // virtual gStoreRepP recompose(const gMap &m, MPI_Comm comm=MPI_COMM_WORLD )const{
    //   return Rep()->recompose(m, comm);} 
    // virtual gStoreRepP recompose(const gMultiMap &m, MPI_Comm comm=MPI_COMM_WORLD )const{
    //   return Rep()->recompose(m, comm);}


    
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split, comm);}
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 const gMap& remap, MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split,remap,comm);}
    // the redistribute takes a vector of gEntitySets as domain
    // distribution over a group of processes and
    // redistributes the gStores according to the domain partition
    // NOTE: should be pure virtual (= 0), but here we default it
    // to do nothing just to be able to compile the code, all
    // containers will need to implement this method later.
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->split_redistribute(dom_ptn, comm);}
    
    //insert elements into store
    void insert(gEntity e, const T &val){
      CPTR<storeType> p(Rep()) ;
      if(p != 0)
        p->insert(e, val);
      warn(p==0);
    }

    void erase(iterator itr1, iterator itr2){
      CPTR<storeType> p(Rep()) ;
      if(p != 0)
        p->erase(itr1, itr2);
      warn(p==0);
    }
    
    
    void insert(const gEntitySet& seq,  const T* vals){
      CPTR<storeType> p(Rep()) ;
      if(p != 0)
        p->insert(seq, vals);
      warn(p==0);
    }

    
    void local_sort(){
      CPTR<storeType> p(Rep()) ;
      if(p != 0)
        p->local_sort();
      warn(p==0);
    }
    
    void sort(){
      CPTR<storeType> p(Rep()) ;
      if(p != 0)
        p->sort();
      warn(p==0);
    }
    
    
   
    void make_consistent();
    
  } ;

  
   
} // end of namespace Loci

#endif
