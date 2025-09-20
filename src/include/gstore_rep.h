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
#ifndef GSTORE_REP_H
#define GSTORE_REP_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Config/conf.h>
#include <Tools/debug.h>
#include <Tools/cptr.h>
#include <Tools/intervalSet.h>
#include <istream>
#include <ostream>
#include <mpi.h>
#include <data_traits.h>
#include <Tools/nptr.h>
//#include <frame_info.h>
//for test
//also need the declaration of class frame_info
//#include <store_rep.h>
namespace Loci {
  enum gstore_type { GSTORE, GMULTISTORE, GMAP, GMULTIMAP,GMAPVEC, GPARAMETER,GCONSTRAINT, GBLACKBOX } ;
 
 
  
  class gMap ;
  class gMapRep;
  class gKeySpace;
  typedef CPTR<gMapRep> gMapRepP ;
  
  class gStoreRep;
  typedef CPTR<gStoreRep> gStoreRepP ;
  typedef const_CPTR<gStoreRep> const_gStoreRepP ;
  
  class storeRep ;
  typedef NPTR<storeRep> storeRepP ;
  typedef const_NPTR<storeRep> const_storeRepP ;

  class gKeySpace;
  typedef CPTR<gKeySpace> gKeySpaceP ;

  class gStoreRep : public CPTR_type {
  public:
    //the convenient functions for vector operations 

    //for stores and maps , the return value is not the size of domain,
    //but the size of vector storage 
    virtual size_t size()const {
      std::cerr << "gStoreRep:size() should not be called"
                << std::endl ;
      return 0;}
    virtual void reserve(size_t n){}
    virtual void clear(){}
    virtual void local_sort(){}
           
    virtual ~gStoreRep(){}

    virtual void set_domain_space(gKeySpaceP space) = 0;
    virtual gKeySpaceP get_domain_space()const = 0;

    virtual gStoreRepP clone() const = 0 ;
#ifdef COPY2STORE  
    virtual storeRepP copy2store() const = 0;
#endif
    // the remap method merely renumbers the container
    // according to the passed in map
    virtual gStoreRepP remap(const gMap &m) const = 0 ;
    
    
    //For stores, recompose will compose a new store whose domain is the domain of m,
    //whose data is the data of the SECOND field of m.
    //for example, pos.recompose(face2node) will producea  store
    //of the positions  for each face
    //For maps, in necessary, expanded map will be generated from m 
    virtual gStoreRepP  recompose(gStoreRepP &m, MPI_Comm comm=MPI_COMM_WORLD)const{
      std::cerr << "gStoreRep.recompose() is not implemented yet"
                << std::endl ;
      abort() ;
      return gStoreRepP(0) ;
    }


    
    // this method redistributes the stores according to the split of local domain over a group of process
    //dom_split: the send split of local domain
    // NOTE: should be pure virtual (= 0), but here we default it
    // to do nothing just to be able to compile the code, all
    // containers will need to implement this method later.
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,//send split
                 MPI_Comm comm=MPI_COMM_WORLD)const {
      std::cerr << "gStoreRep.redistribute() is not implemented yet"
                << std::endl ;
      abort() ;
      return gStoreRepP(0) ;
    }

    //this method redistributes store according to the partition of the whole domain
    //dom_ptn: the partition of the whole domain
    //first the local domain is splitted according to dom_ptn, then redistribute is performed
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm=MPI_COMM_WORLD)const {
      std::cerr << "gStoreRep.split_redistribute() is not implemented yet"
                << std::endl ;
      abort() ;
      return gStoreRepP(0) ;
    }
   
    
    // this redistribute version takes an additional remap
    // argument, after redistribution, the new store is remapped
    // dom_split:  the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 const gMap& remap, MPI_Comm comm=MPI_COMM_WORLD)const {
      std::cerr << "gStoreRep.redistribute() is not implemented yet"
                << std::endl ;
      abort() ;
      return gStoreRepP(0);
    }
    
    virtual int pack_size(const gEntitySet &e) const= 0;
  
    virtual void pack(void *ptr, int &loc, int size,  const gEntitySet &e)const = 0 ;
    virtual void unpack(const void *ptr, int &loc, int size) = 0 ;

    //return the store type, such as GSTORE, GMAP, GMULTIMAP, etc
    virtual gstore_type RepType() const = 0 ;
    virtual std::ostream &Print(std::ostream &s) const = 0 ;
    virtual std::istream &Input(std::istream &s) = 0 ;
   
    virtual gEntitySet domain() const = 0 ;
    virtual gStoreRepP getRep() { return gStoreRepP(this) ;} 
    virtual gStoreRepP getRep() const{ return gStoreRepP(const_cast<gStoreRep *>(this)) ; } 

    //return the data type of store
    virtual  DatatypeP getType() const = 0 ;
    //return the frame info
    virtual frame_info get_frame_info()const = 0 ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                          const char* name, frame_info &fi, const gEntitySet &en) = 0;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                           const char* name, const gEntitySet &en) const = 0;
   
   
    //return a pointer to storage of stored and maps
    virtual void*  get_attrib_data(){
      std::cerr << "gStoreRep.get_attribute_data() is not implemented yet"
                << std::endl ;
      abort() ;
      return NULL ;
    }

    //return a const pointer to storage of stored and maps
    virtual const void*  get_attrib_data() const{
      std::cerr << "gStoreRep.get_attribute_data() is not implemented yet"
                << std::endl ;
      abort() ;
      return NULL ; 
    }
  } ;


  class gstore_instance {
    gStoreRepP rep ;
  public:
    enum instance_type { READ_WRITE, READ_ONLY } ;
    void setRep(const gStoreRepP &p)
    { rep = p ; }
    gStoreRepP Rep() { return rep->getRep(); }
    gStoreRepP Rep() const { return rep->getRep(); }
    gStoreRepP getRep() {return rep->getRep() ;}
    gStoreRepP getRep() const {return rep->getRep() ;}
    virtual instance_type access() const {return READ_WRITE;}
  } ;
}

#endif
