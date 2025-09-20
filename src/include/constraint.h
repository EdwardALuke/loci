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
#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/debug.h>
#include <store_rep.h>
#include <istream>
#include <ostream>

namespace Loci {
    
  class constraintRep : public storeRep {
    entitySet constraint_set ;
  public:
    constraintRep() ;
    constraintRep(const entitySet &p) ;
    virtual ~constraintRep() ;
    virtual void allocate(const entitySet &p) ;
    virtual void shift(int_type offset) ;
    virtual storeRep *new_store(const entitySet &p) const ;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const ;
    virtual storeRepP remap(const dMap &m) const ;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual void copy(storeRepP &st, const entitySet &context) ;
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context)  ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context) ;
    
    virtual int pack_size(const entitySet& e, entitySet& packed) ;
    virtual int pack_size(const entitySet &e) ;
    virtual int estimated_pack_size(const entitySet &e) ;
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e) ;
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual store_type RepType() const ;
    virtual entitySet domain() const ;
    virtual std::ostream &Print(std::ostream &s) const ;
    virtual std::istream &Input(std::istream &s) ;
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
#ifdef H5_HAVE_PARALLEL
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en,hid_t xfer_plist_id ) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist_id) const ;
#endif
    entitySet *get_constraint() { return &constraint_set ; }
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;

#ifdef DYNAMICSCHEDULING
    virtual storeRepP freeze(const entitySet& es) const {
      std::cerr << "storeRep.freeze(e) is not implemented yet"
                << std::endl ;
      abort() ;
      return storeRepP(0) ;
    }
    virtual storeRepP thaw(const entitySet& es) const {
      std::cerr << "storeRep.freeze(e) is not implemented yet"
                << std::endl ;
      abort() ;
      return storeRepP(0) ;
    }
    virtual void pack(void* ptr, int& loc,
                      int& size, const entitySet& e, const Map& remap) {
      pack(ptr,loc,size,e) ;
    }
    virtual void unpack(void* ptr, int& loc,
                        int& size, const sequence& seq, const dMap& remap) {
      unpack(ptr,loc,size,seq) ;
    }
#endif
  } ;

  class constraint : public store_instance {
    typedef constraintRep constraintType ;
      entitySet *data ;
  public:
    constraint() ;
    constraint(const constraint &var) ;
    constraint(const storeRepP &rp) { setRep(rp) ; }
    virtual ~constraint() ;

    constraint & operator=(const constraint &p)
    { *data = *(p.data) ; return *this ;}
    constraint & operator=(const storeRepP &p)
    { setRep(p) ;  return *this ;}
    constraint & operator=(const entitySet &v)
    { *data = v ; return *this ; }

    virtual void notification() ;
    
    entitySet &operator*() { return *data ; }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const constraint &t)
    { return t.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, constraint &t)
    { return t.Input(s) ; }

  // a const version of constraint class
  class const_constraint: public store_instance {
    typedef constraintRep constraintType ;
    entitySet *data ;
  public:
    const_constraint() { setRep(new constraintType) ;}
    const_constraint(const constraint& var) {setRep(var.Rep()) ;}
    const_constraint(const storeRepP& rp) {setRep(rp) ;}
    const_constraint& operator=(const const_constraint& p) {
      setRep(p.Rep()) ;
      return *this ;
    }
    const_constraint& operator=(const storeRepP& p) {
      setRep(p) ;
      return *this ;
    }
    const_constraint& operator=(const entitySet& v) {
      *data = v ;
      return *this ;
    }
    virtual ~const_constraint() {}
    virtual void notification() {
      NPTR<constraintType> p(Rep());
      if(p!=0)
        data = p->get_constraint() ;
      warn(p==0);
    }
    virtual instance_type access() const {return READ_ONLY ;}
    const entitySet& operator*() const {return *data ;}
    std::ostream& Print(std::ostream& s) const {return Rep()->Print(s) ;}
  } ;
  inline std::ostream & operator<<(std::ostream &s, const const_constraint &t)
    { return t.Print(s) ; }
  
  // a duplicate of the constraint class
  class Constraint : public store_instance {
    typedef constraintRep constraintType ;
    entitySet *data ;
  public:
    Constraint() ;
    Constraint(const Constraint &var) ;
    Constraint(const storeRepP &rp) { setRep(rp) ; }
    virtual ~Constraint() ;

    Constraint & operator=(const Constraint &p)
    { setRep(p.Rep()) ; return *this ;}
    Constraint & operator=(const storeRepP &p)
    { setRep(p) ;  return *this ;}
    Constraint & operator=(const entitySet &v)
    { *data = v ; return *this ; }

    virtual void notification() ;
    
    entitySet &operator*() { return *data ; }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s) ; }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const Constraint &t)
    { return t.Print(s) ; }

  inline std::istream & operator>>(std::istream &s, Constraint &t)
    { return t.Input(s) ; }
  
}

#endif
    
