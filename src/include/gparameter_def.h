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
#ifndef GPARAMETER_DEF_H
#define GPARAMETER_DEF_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>
#include <mpi.h>
#include <Tools/debug.h>
#include <gstore_rep.h>
#include <data_traits.h>
#include <distribute_long.h>
#include <gkeyspace.h>
#include <string>
namespace Loci {
  
  template<class T> class gParamRepI : public gStoreRep {
    gEntitySet store_domain ;
    T attrib_data ;
    gKeySpaceP domain_space;
  private:
    int get_mpi_size( IDENTITY_CONVERTER c, const gEntitySet &eset)const;
    int get_mpi_size( USER_DEFINED_CONVERTER c, const gEntitySet &eset)const;
    void packdata(IDENTITY_CONVERTER c,     void *ptr, int &loc, int size )const;
    void packdata(USER_DEFINED_CONVERTER c, void *ptr, int &loc, int size )const;
    void unpackdata(IDENTITY_CONVERTER c,     const void *ptr, int &loc, int size);
    void unpackdata(USER_DEFINED_CONVERTER c, const void *ptr, int &loc, int size);
    DatatypeP getType(IDENTITY_CONVERTER g)const ;
    DatatypeP getType(USER_DEFINED_CONVERTER g)const ;
    frame_info get_frame_info(IDENTITY_CONVERTER g)const;
    frame_info get_frame_info(USER_DEFINED_CONVERTER g)const;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset,hsize_t dimension,
                   const char* name, IDENTITY_CONVERTER g, const gEntitySet &eset) const;
    void hdf5write(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                   const char* name, USER_DEFINED_CONVERTER g, const gEntitySet &eset) const;
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                  const char* name, IDENTITY_CONVERTER g, frame_info &fi, const gEntitySet &en);
    void hdf5read(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                  const char* name, USER_DEFINED_CONVERTER g, frame_info &fi, const gEntitySet &en);
  
  public:

    gParamRepI():store_domain(gInterval(GUNIVERSE_MIN,GUNIVERSE_MAX)){
      std::string casename;
      gKeySpaceP space = gKeySpace::get_space("UniverseSpace", casename);
      domain_space = space;
    }
    gParamRepI(const gEntitySet &p):store_domain(p){
      std::string casename;
      gKeySpaceP space = gKeySpace::get_space("UniverseSpace", casename);
      domain_space = space;
    }
    virtual void allocate(const gEntitySet &p){  store_domain = p ; }
    virtual ~gParamRepI() {}
    void set_domain_space(gKeySpaceP space){domain_space = space;}
    gKeySpaceP get_domain_space()const{return domain_space;}
    // this method redistributes the container according to the split of local domain over a group of process
    //dom_split: the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const;

    // the redistribute takes a vector of gEntitySets as domain
    // distribution over a group of processes and
    // redistributes the container according to the domain partition
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const;

    // this redistribute version takes an additional remap
    // argument, after redistribution, the new container is remapped
    // dom_split:  the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 const gMap& remap, MPI_Comm comm=MPI_COMM_WORLD)const;
    //  virtual gStoreRepP
    //         redistribute_omd(const std::vector<gEntitySet>& dom_ptn,
    //                          const dMap& remap, MPI_Comm comm=MPI_COMM_WORLD)  ;

    virtual gstore_type RepType() const { return GPARAMETER ;}
    virtual gEntitySet domain() const { return store_domain ;}
    virtual gStoreRepP remap(const gMap &m) const;
    virtual gStoreRepP clone() const;
    virtual void* get_attrib_data() {
      return static_cast<void*>(&attrib_data);
    }
    virtual const void* get_attrib_data() const{
      return  static_cast<const void*>(&attrib_data);
    }
#ifdef COPY2STORE  
    virtual storeRepP copy2store()const;
#endif
     
    virtual int pack_size(const gEntitySet &e)const ;
    virtual void pack(void *ptr, int &loc, int size, const gEntitySet &e)const ;//only pack attrib_data
    virtual void unpack(const void *ptr, int &loc, int size)  ;//only unpack attribute_data

    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual DatatypeP getType() const ;

    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                          const char* name, frame_info &fi, const gEntitySet &user_eset) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                           const char* name, const gEntitySet& en) const ;
    virtual frame_info get_frame_info()const ;
     T *get_param() { return &attrib_data ; }
  } ;

 
  //**************************************************************************/


  template<class T> class gParam : public gstore_instance {
    typedef gParamRepI<T> gParamType ;
    T * data ;
  private:
    void connect_data(){
      CPTR<gParamType> p(Rep()) ;
      if(p != 0) data = p->get_param();
    }
  public:
    typedef T containerType ;
    gParam() { setRep(new gParamType) ;
      // data = static_cast<T*> (Rep()->get_attrib_data());
      connect_data();
    }
    
    //deep copy
    gParam(const gParam &var) {
      setRep(new gParamType) ;
      connect_data();
      *data = *(var.data) ;
      set_entitySet(var.domain());
    }
    
    gParam(gStoreRepP rp) {
      setRep(rp);
      connect_data();
      //data = static_cast<T*> (Rep()->get_attrib_data());   
    }

    virtual ~gParam() {}
    //deep assignment
    gParam & operator=(const gParam &p) {
      // setRep(p.Rep());
      //data = static_cast<T*> (Rep()->get_attrib_data());
      *data = *(p.data) ;
      set_entitySet(p.domain());
      return *this ;
    }
    
    gParam & operator=(gStoreRepP p) {
      setRep(p) ;
      // data = static_cast<T*> (Rep()->get_attrib_data());
      connect_data();
      return *this ;
    }

    gParam & operator=(const T &v) {
      *data = v ;
      return *this ;
    }
    
    T * operator->() { return data ; }

    const T * operator->() const { return data ; }

    T * operator&() { return data ; }

    const T * operator &() const { return data ; }

    T &operator*() { return *data ; }

    const T &operator*() const { return *data ; }

    void set_domain_space(gKeySpaceP space){Rep()->set_domain_space(space);}
    gKeySpaceP get_domain_space()const{return Rep()->get_domain_space();}
    
    void set_entitySet(const gEntitySet &ptn) {
      CPTR<gParamType> p(Rep()) ;
      if(p != 0) p->allocate(ptn);
      warn(p==0) ; 
    }

    gEntitySet domain() const { return Rep()->domain(); }

    // this method redistributes the container according to the split of local domain over a group of process
    //dom_split: the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split, comm);}

    // this redistribute version takes an additional remap
    // argument, after redistribution, the new container is remapped
    // dom_split:  the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 const gMap& remap, MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->redistribute(dom_split,remap,comm);}


    // the redistribute takes a vector of gEntitySets as domain
    // distribution over a group of processes and
    // redistributes container according to the domain partition
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const{return Rep()->split_redistribute(dom_ptn, comm);}
    

    // the remap method merely renumbers the container
    // according to the passed in map
    virtual gStoreRepP remap(const gMap &m) const{return Rep()->remap(m);} 
    virtual gStoreRepP clone() const{return Rep()->clone();}
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;
   
}

#endif
