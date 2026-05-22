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
#ifndef GBLACKBOX_H
#define GBLACKBOX_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <mpi.h>
#include <blackbox.h>
#include <Config/conf.h>
#include <Tools/debug.h>
#include <gstore_rep.h>
#include <istream>
#include <ostream>
#include <data_traits.h>

namespace Loci {
  
  template<class T> class gBlackboxRepI : public gStoreRep {
    gEntitySet store_domain;
    T *attrib_data;
    gKeySpaceP domain_space;
  public:
    gBlackboxRepI() { attrib_data=0; store_domain = ~GEMPTY; }
    gBlackboxRepI(const gEntitySet &p) { attrib_data=0; store_domain = p;}
    virtual void allocate(const gEntitySet &p) ;
   
    virtual ~gBlackboxRepI();
    virtual void set_domain_space(gKeySpaceP space){domain_space = space;}
    virtual gKeySpaceP get_domain_space()const{return domain_space;}

    // this method redistributes the container according to the split of local domain over a group of process
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
    virtual gstore_type RepType() const;
    virtual gEntitySet domain() const;
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
    // virtual int pack_size(const gEntitySet& e, gEntitySet& packed) ;
    virtual int pack_size(const gEntitySet &e) const;
    // virtual int estimated_pack_size(const gEntitySet &e);
    
    virtual void pack(void *ptr, int &loc, int &size, const gEntitySet &e)const;
    virtual void unpack(void *ptr, int &loc, int &size) ;
    
    virtual std::ostream &Print(std::ostream &s) const;
    virtual std::istream &Input(std::istream &s);
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset,
                          hsize_t dimension, const char* name, frame_info &fi, const gEntitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset,
                           hsize_t dimension, const char* name, const gEntitySet& en) const ;
    virtual DatatypeP getType()const ;
    virtual frame_info get_frame_info() const ;
    T *get_blackbox() { return attrib_data; }
   
  };

  //**************************************************************************/

  template<class T> void gBlackboxRepI<T>::allocate(const gEntitySet &p) {
    if(p == GEMPTY) {
      if(attrib_data != 0) {
        delete attrib_data ;
        attrib_data = 0 ;
      }
      attrib_data = new T ;
    } else {
      if(attrib_data == 0)
        attrib_data = new T ;
    }
    store_domain = p;
  }

  //**************************************************************************/

  template<class T> gBlackboxRepI<T>::~gBlackboxRepI()
  {
    if(attrib_data != 0)
      delete attrib_data ;
  }

  
  //**************************************************************************/
  
  template<class T> 
  gstore_type gBlackboxRepI<T>::RepType() const 
  {
    return GBLACKBOX;
  }

  //**************************************************************************/

  template<class T> gEntitySet gBlackboxRepI<T>::domain() const {
    return store_domain;
  }

  //**************************************************************************/
        
  template<class T> 
  std::ostream &gBlackboxRepI<T>::Print(std::ostream &s) const 
  {
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    s << '{' << domain() << std::endl;
    s << '}' << std::endl;
    return s;
  }

  //**************************************************************************/

  template<class T> 
  std::istream &gBlackboxRepI<T>::Input(std::istream &s) 
  {
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    return s;
  }

  //**************************************************************************/
  template<class T> 
  void gBlackboxRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                  const char* name, frame_info &fi, const gEntitySet &en) 
  {
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }
    
  //**************************************************************************/
    
  template<class T> 
  void gBlackboxRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension,
                                   const char* name, const gEntitySet& en) const
  {
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }
    
  //**************************************************************************/
    
  template<class T> 
  DatatypeP gBlackboxRepI<T>::getType()const
  {
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    return 0 ;
  }
  //**************************************************************************/
    
    
  template<class T> 
  frame_info gBlackboxRepI<T>::get_frame_info()const
  {
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    frame_info fi ;
    return fi ;
  }

  //**************************************************************************/
  template<class T> gStoreRepP  gBlackboxRepI<T>::
  split_redistribute(const std::vector<gEntitySet>& dom_ptn, MPI_Comm comm)const {
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    return 0 ;
  }
  //**************************************************************************/

  template<class T> gStoreRepP gBlackboxRepI<T>::
  redistribute(const std::vector<gEntitySet>& dom_split,
               MPI_Comm comm)const{
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    return 0 ;
  }

  //**************************************************************************/
  
  template<class T> gStoreRepP gBlackboxRepI<T>::
  redistribute(const std::vector<gEntitySet>& dom_split,
               const gMap& remap, MPI_Comm comm)const{
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    return 0 ;
  }
    
  //**************************************************************************/
    
  template<class T> class gBlackbox : public gstore_instance {
    typedef gBlackboxRepI<T> gBlackboxType;
    T * data;
    gBlackbox(const gBlackbox &var) { setRep(var.Rep()); }
    gBlackbox & operator=(const gBlackbox &p) {setRep(p.Rep()); return *this; }
    void connect_data(){
      CPTR<gBlackboxType> p(Rep()) ;
      if(p != 0) data = p->get_blackbox();
    } 
  public:
    typedef T containerType;
    gBlackbox() { setRep(new gBlackboxType); connect_data(); }

    gBlackbox(gStoreRepP rp) { setRep(rp); connect_data(); }

    virtual ~gBlackbox();


    gBlackbox & operator=(gStoreRepP p) {setRep(p); connect_data(); return *this; }
    gBlackbox & operator=(const T &v) { *data = v;  return *this; }

    
    
    T * operator->() { return data; }
    const T * operator->() const { return data; }
    
    T &operator*() { return *data; }
    const T &operator*() const { return *data; }

    T &operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(data == NULL);
      fatal(!Rep()->domain().inSet(indx));
#endif
      return *data;
    }

    const T &operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(data == NULL);
      fatal(!Rep()->domain().inSet(indx));
#endif
      return *data;
    }

    //    void set_entitySet(const gEntitySet &ptn) {Rep()->allocate(ptn); }

    gEntitySet domain() const { return Rep()->domain(); }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s); }
  };

  //**************************************************************************/

  template<class T> gBlackbox<T>::~gBlackbox() {}
    
 

  //**************************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const gBlackbox<T> &t)
  {
    return t.Print(s);
  }

  //**************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, gBlackbox<T> &t)
  {
    return t.Input(s);
  }

  //  //**************************************************************************/

  //   template<class T> 
  //   class const_gBlackbox : public gstore_instance {
  //     typedef T containerType;
  //     typedef gBlackboxRepI<T> gBlackboxType;
  //     const T * data;
  //     const_gBlackbox(const const_gBlackbox<T> &var) { setRep(var.Rep()); }
  //     const_gBlackbox(const gBlackbox<T> &var) { setRep(var.Rep()); }
  //     const_gBlackbox & operator=(const const_gBlackbox<T> &p)
  //     { setRep(p.Rep); return *this;}
  //     const_gBlackbox & operator=(const gBlackbox<T> &p)
  //     { setRep(p.Rep); return *this;}
  //   public:
  //     const_gBlackbox() { setRep(new gBlackboxType); }
  //     const_gBlackbox(gStoreRepP rp) { setRep(rp); }
    
  //     virtual ~const_gBlackbox();

  //     const_gBlackbox & operator=(gStoreRepP p)
  //     { setRep(p); return *this;}

  //     virtual void notification();
  //     virtual instance_type access() const;
        
  //     const T * operator->() const { return data; }
    
  //     const T &operator*() const { return *data; }

  //     const T &operator[](int indx) const {
  // #ifdef BOUNDS_CHECK
  //       fatal(data == NULL);
  //       fatal(!Rep()->domain().inSet(indx));
  // #endif
  //       return *data;
  //     }

  //     gEntitySet domain() const { return Rep()->domain(); }
    
  //     std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
  //   };

  //   //**************************************************************************/

  //   template<class T> const_gBlackbox<T>::~const_gBlackbox() {}

  
  //**************************************************************************/

  template<class T> 
  gStoreRepP gBlackboxRepI<T>::remap(const gMap &m) const 
  {
    gBlackbox<T> r;
    r.set_entitySet(m.image(m.domain()&domain()));
    if(attrib_data!=0)
      *r = *attrib_data;
    return r.Rep();
  }

  
  //**************************************************************************/

#ifdef COPY2STORE  
  template<class T> 
  storeRepP gBlackboxRepI<T>::copy2store()const 
  {
    blackbox<T> r;
    r.set_entitySet(store_domain) ;
    if(attrib_data != 0)
      *r = *attrib_data ;
    return r.Rep();
  }
#endif
 
  //**************************************************************************/
 
  template <class T> 
  int gBlackboxRepI<T>::pack_size(const gEntitySet &eset)const 
  {
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    return 0;
  }
  
 

  //**************************************************************************/

  template <class T> 
  void gBlackboxRepI<T>::pack(void *ptr, int &loc, int &size,
                              const gEntitySet &e)const 
  {
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }

  //**************************************************************************/

  template <class T> 
  void gBlackboxRepI<T>::unpack(void *ptr, int &loc, int &size)
  {
    cerr << "GBLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }  


  // //**************************************************************************/

  // template<class T>
  // gStore_instance::instance_type const_gBlackbox<T>::access() const
  // {
  //   return READ_ONLY;
  // }

  // //**************************************************************************/
    
  // template<class T> 
  // inline std::ostream & operator<<(std::ostream &s,
  //                                  const const_gBlackbox<T> &t)
  // {
  //   return t.Print(s);
  // }

  //**************************************************************************/
}

#endif
    
