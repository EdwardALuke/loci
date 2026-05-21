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
#ifndef GCONSTRAINT_H
#define GCONSTRAINT_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/debug.h>
#include <gstore_rep.h>
#include <istream>
#include <ostream>
#include <distribute_long.h>
#include <gkeyspace.h>

namespace Loci {
  
  class gConstraintRep : public gStoreRep {
    gEntitySet constraint_set ;
    gKeySpaceP domain_space;
  public:
    gConstraintRep():domain_space(0){} 
    gConstraintRep(const gEntitySet &p):constraint_set(p),domain_space(0){}
    virtual ~gConstraintRep(){}
    void set_domain_space(gKeySpaceP space){domain_space = space;}
    gKeySpaceP get_domain_space()const{return domain_space;}
    
    virtual void allocate(const gEntitySet &p){
      constraint_set = p;
    }

    // this method redistributes the container according to the split of
    //local domain over a group of process
    //dom_split: the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm LOCI_DEFAULT_COMM)const;

    // the redistribute takes a vector of gEntitySets as domain
    // distribution over a group of processes and
    // redistributes the container according to the domain partition
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm LOCI_DEFAULT_COMM)const;


    // this redistribute version takes an additional remap
    // argument, after redistribution, the new container is remapped
    // dom_split:  the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 const gMap& remap, MPI_Comm comm LOCI_DEFAULT_COMM)const;
    
    
    virtual gStoreRepP remap(const gMap &m) const ;
    virtual int pack_size(const gEntitySet &e)const ;
    virtual void pack(void *ptr, int &loc, int size, const gEntitySet &e)const ;
    virtual void unpack(const void *ptr, int &loc, int size) ;
    virtual size_t size() const{return constraint_set.size();}
    virtual gStoreRepP clone() const{
      return new gConstraintRep(constraint_set); 
    }
#ifdef COPY2STORE  
    virtual storeRepP copy2store()const ;
#endif
    virtual void* get_attrib_data() {
      return static_cast<void*>(&constraint_set);
    }
    virtual const void* get_attrib_data() const{
      return  static_cast<const void*>(&constraint_set);
    }

  
    virtual gstore_type RepType() const{ return GCONSTRAINT  ; }
    virtual gEntitySet domain() const{ return constraint_set ; }
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual DatatypeP getType()const ;
    virtual frame_info get_frame_info() const ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, const gEntitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, const gEntitySet &en) const ;
  } ;

  class gConstraint : public gstore_instance {
    typedef gConstraintRep gConstraintType ;
    gEntitySet *data ;
  public:
    gConstraint() 
    {
      setRep(new gConstraintType) ;
      data = static_cast<gEntitySet*> (Rep()->get_attrib_data());
    }
    
    gConstraint(const gConstraint &var)
    {
      setRep(var.Rep());
      data = static_cast<gEntitySet*> (Rep()->get_attrib_data());
    }
    gConstraint(const gStoreRepP &rp) {
      setRep(rp) ;
      data = static_cast<gEntitySet*> (Rep()->get_attrib_data());
    }
    virtual ~gConstraint(){}

    void set_domain_space(gKeySpaceP space){Rep()->set_domain_space(space);}
    gKeySpaceP get_domain_space()const{return Rep()->get_domain_space();}
    virtual gEntitySet domain() const {return Rep()->domain() ;}

    // this method redistributes the container according to the split of
    //local domain over a group of process
    //dom_split: the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 MPI_Comm comm LOCI_DEFAULT_COMM)const{return Rep()->redistribute(dom_split, comm);}

    // the redistribute takes a vector of gEntitySets as domain
    // distribution over a group of processes and
    // redistributes the container according to the domain partition
    virtual gStoreRepP
    split_redistribute(const std::vector<gEntitySet>& dom_ptn,
                       MPI_Comm comm LOCI_DEFAULT_COMM)const{return Rep()->split_redistribute(dom_ptn, comm);}


    // this redistribute version takes an additional remap
    // argument, after redistribution, the new container is remapped
    // dom_split:  the send split of local domain
    virtual gStoreRepP
    redistribute(const std::vector<gEntitySet>& dom_split,
                 const gMap& remap, MPI_Comm comm LOCI_DEFAULT_COMM)const{return Rep()->redistribute(dom_split,remap,comm);}

    virtual gStoreRepP remap(const gMap &m) const{return Rep()->remap(m) ;}
    gConstraint & operator=(const gConstraint &p)
    { setRep(p.Rep()) ;
      data = static_cast<gEntitySet*> (Rep()->get_attrib_data());
      return *this ;}
    gConstraint & operator=(const gStoreRepP &p)
    { setRep(p) ;
      data = static_cast<gEntitySet*> (Rep()->get_attrib_data());
      return *this ;}
    gConstraint & operator=(const gEntitySet &v)
    { setRep(new gConstraintType(v));
      data = static_cast<gEntitySet*> (Rep()->get_attrib_data());
      return *this ; }

    virtual gStoreRepP clone() const{return Rep()->clone() ;}
    
    gEntitySet &operator*() { return *data ; }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const gConstraint &t)
  { return t.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, gConstraint &t)
  { return t.Input(s) ; }
    
}


#endif
    
