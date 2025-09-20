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
#ifndef GMAPVEC_DEF_H_
#define GMAPVEC_DEF_H_

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <Tools/debug.h>
#include <gmap.h>
#include <gkeyspace.h>
namespace Loci {

  class gMap ;
  template<unsigned int M> class gMapVec ;
  template<unsigned int M> class const_gMapVec;
  
  template<unsigned int M> class gMapVecRepI : public gMapRep {
    bool sorted;
    std::vector<std::pair<gEntity,gEntity> > attrib_data;
    
    gKeySpaceP domain_space;
    gKeySpaceP image_space;
  public:
    /**Convenient functions for vertor operations
     */
    typedef  std::vector<std::pair<gEntity,gEntity> >::iterator iterator ;
    typedef  std::vector<std::pair<gEntity,gEntity> >::const_iterator const_iterator ;
    iterator begin() { return attrib_data.begin(); }
    iterator end() { return attrib_data.end(); }
    const_iterator begin()const { return attrib_data.begin(); }
    const_iterator end() const { return attrib_data.end(); }
    size_t size() const{return attrib_data.size();}
    void reserve (size_t n){attrib_data.reserve(n);}
    void clear(){attrib_data.clear();}
    void local_sort();
        
    gMapVecRepI():sorted(true),domain_space(0),image_space(0) {}
    
    void set_domain_space(gKeySpaceP space){domain_space = space;}
    gKeySpaceP get_domain_space()const{return domain_space;}
    void set_image_space(gKeySpaceP space){image_space = space;}
    gKeySpaceP get_image_space()const{return image_space;}
   
    virtual ~gMapVecRepI(){} 
    virtual gStoreRepP clone() const{return new gMapVecRepI(*this); }  
    virtual gStoreRepP remap(const gMap &m) const ;
    virtual void inplace_remap(const gMap &m); 
    
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const ;

    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const ;
    
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_ptn,
                 const gMap& remap, MPI_Comm comm=MPI_COMM_WORLD)const ;
   
    virtual void inplace_compose(const gMap &m, MPI_Comm comm=MPI_COMM_WORLD ) ;
    virtual gStoreRepP recompose(const gMap &m, MPI_Comm comm=MPI_COMM_WORLD )const ;
    virtual gStoreRepP recompose(const gMultiMap &m, MPI_Comm comm=MPI_COMM_WORLD )const ;
    virtual gStoreRepP recompose( gStoreRepP &m, MPI_Comm comm=MPI_COMM_WORLD )const ;
    
    virtual gStoreRepP local_inverse() const;
    virtual gStoreRepP distributed_inverse(gEntitySet global_image,
                                           gEntitySet global_preimage,
                                           const std::vector<gEntitySet> &init_ptn) const;
    virtual int pack_size(const gEntitySet &e)const ;
    virtual void pack(void *ptr, int &loc, int size, const gEntitySet &e)const ;
    virtual void unpack(const void *ptr, int &loc, int size) ;
  
    virtual void insert(gEntity e, gEntity val){
      sorted = false;
      attrib_data.push_back(std::pair<gEntity, gEntity>(e, val));
      
    }

    virtual void insert(const gEntitySet& seq,  const gEntity* vals){
      sorted = false;
      size_t idx= 0;
      for(gEntitySet::const_iterator itr = seq.begin();
          itr!= seq.end(); itr++){
        attrib_data.push_back(std::pair<gEntity, gEntity>(*itr, vals[idx++]));
      }
      
    }
    gEntitySet domain() const {
      if(!sorted) const_cast<gMapVecRepI&>(*this).local_sort();
      gEntitySet dom = GEMPTY;
      for(const_iterator itr = begin(); itr != end(); itr++){
        dom += itr->first;
      }
      return dom;
    }
      
    virtual gstore_type RepType() const {return GMAPVEC;}
    
    virtual bool isSorted() const {return sorted;}
    virtual gEntitySet image(const gEntitySet &domain) const ;
    virtual gEntitySet image( gEntity domain) const ;
    virtual gEntitySet image() const ;
    
    virtual std::pair<gEntitySet,gEntitySet>
    preimage(const gEntitySet &codomain) const ;
    virtual gStoreRepP get_map()const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual frame_info get_frame_info()const ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                          const char* name, frame_info &fi, const gEntitySet &user_eset) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                           const char* name, const gEntitySet& en) const ;

    //different from traditional maps, this method is const method
    //dom is the domain after expansion, not out_of_dom
    virtual gStoreRepP expand(gEntitySet &dom, std::vector<gEntitySet> &init_ptn,MPI_Comm comm=MPI_COMM_WORLD)const ;
  
    virtual void* get_attrib_data() {return &attrib_data; }
    virtual const void* get_attrib_data() const{return &attrib_data; }
#ifdef COPY2STORE  
    virtual storeRepP copy2store() const;
#endif
    virtual DatatypeP getType()const ;
  
  } ;


  
  template<unsigned int M> class gMapVec : public gstore_instance {
    friend class const_gMapVec<M> ;
    typedef gMapVecRepI<M> MapType ;
    typedef std::vector<std::pair<gEntity,gEntity> > gRep ; 
  public:

    typedef  std::vector<std::pair<gEntity,gEntity> >::iterator iterator ;
    typedef  std::vector<std::pair<gEntity,gEntity> >::const_iterator const_iterator ;
    iterator begin()  { return static_cast<gRep*>(Rep()->get_attrib_data())->begin(); }
    iterator end(){ return static_cast<gRep*>(Rep()->get_attrib_data())->end(); }
    const_iterator begin()const  { return static_cast<const gRep*>(Rep()->get_attrib_data())->begin(); }
    const_iterator end() const  { return static_cast<const gRep*>(Rep()->get_attrib_data())->end(); }
    size_t size() const{return Rep()->size();}
    void reserve (size_t n){Rep()->reserve(n);}
    void clear(){Rep()->clear();}
    void local_sort(){Rep()->local_sort();}
    // void remove_duplication();
    
    
    // These should be private, but too much code currently depends on it
    // being public.  This code is dangerous because it does a shallow
    // copy which means that results sometimes may be unpredictable when
    // these operations are used.
    gMapVec(const gMapVec<M> &var) { setRep(var.Rep()) ; }
    gMapVec & operator=(const gMapVec<M> &str) { setRep(str.Rep()) ;  return *this ;}
    gMapVec & operator=(gStoreRepP p) { setRep(p) ;  return *this ;}
    
    gMapVec() {setRep(new MapType) ; }
    gMapVec(gStoreRepP &rp) { setRep(rp) ; }
    virtual ~gMapVec(){}

    void set_vdom(gEntitySet vd){
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        p->set_vdom(vd);
      warn(p==0) ;
    }
    gEntitySet get_vdom()const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->get_vdom();
      warn(p==0) ;
    }
    void set_domain_space(gKeySpaceP space){static_cast<gMapRepP>(Rep())->set_domain_space(space);}
    gKeySpaceP get_domain_space()const{return static_cast<gMapRepP>(Rep())->get_domain_space() ;}
    void set_image_space(gKeySpaceP space){static_cast<gMapRepP>(Rep())->set_image_space(space);}
    gKeySpaceP get_image_space()const{return static_cast<gMapRepP>(Rep())->get_image_space();}
    
    virtual gStoreRepP clone() const{return Rep()->clone() ;}
        
    virtual  void insert(gEntity e, gEntity val){
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        p->insert(e, val);
      warn(p==0) ;}
    virtual void insert(const gEntitySet& seq,  const gEntity* vals){
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        p->insert(seq, vals);
      warn(p==0);
    }
    virtual instance_type access() const{return READ_WRITE ; }
    virtual gStoreRepP remap(const gMap &m) const{return Rep()->remap(m);}
    virtual gStoreRepP local_inverse() const{return static_cast<gMapRepP>(Rep())->local_inverse();}

    virtual gStoreRepP distributed_inverse(gEntitySet global_image,
                                           gEntitySet global_preimage,
                                           const std::vector<gEntitySet> &init_ptn) const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->distributed_inverse(global_image, global_preimage, init_ptn);
      else
        return gStoreRepP(0);
    }
    
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split, comm) ;}
    
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_ptn, comm) ;}

    //different from traditional maps, this method is const method
    //dom is the domain after expansion, not out_of_dom
    virtual gStoreRepP expand(gEntitySet &dom, std::vector<gEntitySet> &init_ptn,MPI_Comm comm=MPI_COMM_WORLD)const{
      return  gMapRepP(Rep())->expand(dom, init_ptn, comm);}

    
    virtual gStoreRepP recompose(const gMap &m, MPI_Comm comm=MPI_COMM_WORLD )const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->recompose(m, comm);
      else return gStoreRepP(0);
    }

    //duplication removal is needed in this method
    virtual gStoreRepP recompose(const gMultiMap &m, MPI_Comm comm=MPI_COMM_WORLD )const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->recompose(m, comm);
      else return gStoreRepP(0);
    }

    virtual gStoreRepP recompose(gStoreRepP &m, MPI_Comm comm=MPI_COMM_WORLD )const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->recompose(m, comm);
      else return gStoreRepP(0);
    }
    
   
    gEntitySet domain() const { return Rep()->domain() ; }
    operator gMapRepP() {
      gMapRepP p(Rep()) ;
      fatal(p==0) ;
      return p ; }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }

    gEntitySet image(const gEntitySet &dom) const {
      return gMapRepP(Rep())->image(dom) ;
    }
    gEntitySet image(gEntity dom) const {
      CPTR<MapType> p(Rep()) ;
      fatal(p==0);
      return p->image(dom);
    }
    gEntitySet image() const {
      return gMapRepP(Rep())->image() ;
    }
   
    std::pair<gEntitySet,gEntitySet> preimage(const gEntitySet &codomain) const {
      return gMapRepP(Rep())->preimage(codomain) ;
    }
  } ;

  template<unsigned int M> 
  inline std::ostream & operator<<(std::ostream &s, const gMapVec<M> &m)
  { return m.Print(s) ; }

  template<unsigned int M> 
  inline std::istream & operator>>(std::istream &s, gMapVec<M> &m)
  { return m.Input(s) ; }

  template<unsigned int M> class const_gMapVec : public gstore_instance {
    typedef gMapVecRepI<M> MapType ;
    typedef std::vector<std::pair<gEntity,gEntity> > gRep ; 
    const_gMapVec(const const_gMapVec<M> &var) {setRep(var.Rep()) ;}
    const_gMapVec(const gMapVec<M> &var) {setRep(var.Rep());}
    const_gMapVec & operator=(const gMapVec<M> &str) { setRep(str.Rep()) ; return *this ;}
    const_gMapVec & operator=(const const_gMapVec<M> &str)
    { setRep(str.Rep()) ; return *this ;}

  public:

    typedef  std::vector<std::pair<gEntity,gEntity> >::const_iterator const_iterator ;
    const_iterator begin()const { return static_cast<const gRep*>(Rep()->get_attrib_data())->begin(); }
    const_iterator end()const { return static_cast<const gRep*>(Rep()->get_attrib_data())->end(); }
    size_t size() const{return Rep()->size();}

    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split, comm) ;}
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 const gMap& remap, MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split,remap, comm) ;}
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_ptn, comm) ;}

    virtual gStoreRepP remap(const gMap &m) const{return Rep()->remap(m);}
    virtual gstore_type RepType() const {return Rep()->RepType();}
    virtual ~const_gMapVec(){}
    
    const_gMapVec(){ setRep(new MapType);  }
    const_gMapVec(gStoreRepP &rp) { setRep(rp) ; }
    
    virtual instance_type access() const{return READ_ONLY ; }
        
    const_gMapVec & operator=(gStoreRepP p) { setRep(p) ; return *this ;}

    gEntitySet domain() const { return Rep()->domain(); }
    
    gEntitySet image(const gEntitySet &dom) const {
      return gMapRepP(Rep())->image(dom) ;
    }
    gEntitySet image(gEntity dom) const {
      CPTR<MapType> p(Rep()) ;
      fatal(p==0);
      return p->image(dom);
    }
    gEntitySet image() const {
      return gMapRepP(Rep())->image() ;
    }
    
    std::pair<gEntitySet,gEntitySet> preimage(const gEntitySet &codomain) const {
      return gMapRepP(Rep())->preimage(codomain) ;
    }
    
    operator gMapRepP()
    { gMapRepP p(Rep()) ; fatal(p==0) ; return p ; }


    virtual gStoreRepP distributed_inverse(gEntitySet global_image,
                                           gEntitySet global_preimage,
                                           const std::vector<gEntitySet> &init_ptn) const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->distributed_inverse(global_image,
                                      global_preimage,
                                      init_ptn);
      else return gStoreRepP(0);
    }

    virtual gStoreRepP recompose(const gMap &m, MPI_Comm comm=MPI_COMM_WORLD )const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->recompose(m, comm);
      else return gStoreRepP(0);
    }

    virtual gStoreRepP recompose(const gMultiMap &m, MPI_Comm comm=MPI_COMM_WORLD )const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->recompose(m, comm);
      else return gStoreRepP(0);
    }

    virtual gStoreRepP recompose( gStoreRepP &m, MPI_Comm comm=MPI_COMM_WORLD )const{
      CPTR<MapType> p(Rep()) ;
      if(p != 0)
        return p->recompose(m, comm);
      else return gStoreRepP(0);
    }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
  } ;

  template<unsigned int M> 
  inline std::ostream & operator<<(std::ostream &s, const_gMapVec<M> &m)
  { return m.Print(s) ; }
    
}

#endif
