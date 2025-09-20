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
#ifndef BLACKBOX_H
#define BLACKBOX_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <mpi.h>

#include <Config/conf.h>
#include <Tools/debug.h>
#include <store_rep.h>
#include <istream>
#include <ostream>
#include <data_traits.h>

namespace Loci {
  
  template<class T> class blackboxRepI : public storeRep {
    entitySet store_domain;
    T *attrib_data;

  public:
    blackboxRepI() { attrib_data=0; store_domain = interval(UNIVERSE_MIN,UNIVERSE_MAX); }
    blackboxRepI(const entitySet &p) { attrib_data=0; store_domain = p;}
    virtual void allocate(const entitySet &p) ;
    virtual void shift(int_type offset) ;
    virtual ~blackboxRepI();
    virtual store_type RepType() const;
    virtual entitySet domain() const;
    virtual storeRep *new_store(const entitySet &p) const;
    virtual storeRep *new_store(const entitySet &p, const int* cnt) const;
    virtual storeRepP remap(const dMap &m) const;
    virtual storeRepP freeze() ;
    virtual storeRepP thaw() ;
    virtual void copy(storeRepP &st, const entitySet &context);
    virtual void gather(const dMap &m, storeRepP &st,
                        const entitySet &context) ;
    virtual void scatter(const dMap &m, storeRepP &st,
                         const entitySet &context);
    virtual int pack_size(const entitySet& e, entitySet& packed) ;
    virtual int pack_size(const entitySet &e);
    virtual int estimated_pack_size(const entitySet &e);
    
    virtual void pack(void *ptr, int &loc, int &size, const entitySet &e);
    virtual void unpack(void *ptr, int &loc, int &size, const sequence &seq) ;
    
    virtual std::ostream &Print(std::ostream &s) const;
    virtual std::istream &Input(std::istream &s);
    virtual void readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) ;
    virtual void writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const ;
#ifdef H5_HAVE_PARALLEL
    virtual void readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en, hid_t xfer_plist_id) ;
    virtual void writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist_id) const ;
#endif
    virtual DatatypeP getType() ;
    virtual frame_info get_frame_info() ;
    T *get_blackbox() { return attrib_data; }
  };

  //**************************************************************************/

  template<class T> void blackboxRepI<T>::allocate(const entitySet &p) {
    if(p == EMPTY) {
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
    dispatch_notify();
  }

  //**************************************************************************/

  template<class T> blackboxRepI<T>::~blackboxRepI()
  {
    if(attrib_data != 0)
      delete attrib_data ;
  }

  //**************************************************************************/

  template<class T>
  storeRep *blackboxRepI<T>::new_store(const entitySet &p) const
  {
    return new blackboxRepI<T>(p);
  }

  template<class T>
  storeRep *blackboxRepI<T>::new_store(const entitySet &p, const int* cnt) const
  {
    storeRep* sp = 0;
    cerr << " This method should not be called for a dMap " << endl ;
    return sp ; 
  }
  //**************************************************************************/
  
  template<class T> 
  store_type blackboxRepI<T>::RepType() const 
  {
    return BLACKBOX;
  }

  //**************************************************************************/

  template<class T> entitySet blackboxRepI<T>::domain() const {
    return store_domain;
  }

  //**************************************************************************/

  template<class T> void blackboxRepI<T>::shift(int_type offset) {
    cerr << "shift not supported in BLACKBOX!" << endl ;
  }

  //**************************************************************************/
        
  template<class T> 
  std::ostream &blackboxRepI<T>::Print(std::ostream &s) const 
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    s << '{' << domain() << std::endl;
    s << '}' << std::endl;
    return s;
  }

  //**************************************************************************/

  template<class T> 
  std::istream &blackboxRepI<T>::Input(std::istream &s) 
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    return s;
  }

  //**************************************************************************/
  template<class T> 
  void blackboxRepI<T>::readhdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en) 
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }
    
  //**************************************************************************/
    
  template<class T> 
  void blackboxRepI<T>::writehdf5(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en) const
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }
    
 //**************************************************************************/
#ifdef H5_HAVE_PARALLEL
  template<class T> 
  void blackboxRepI<T>::readhdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, frame_info &fi, entitySet &en,hid_t xfer_plist_id ) 
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }
    
  //**************************************************************************/
    
  template<class T> 
  void blackboxRepI<T>::writehdf5P(hid_t group_id, hid_t dataspace, hid_t dataset, hsize_t dimension, const char* name, entitySet& en, hid_t xfer_plist_id) const
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }
#endif  
  //**************************************************************************/
    
  template<class T> 
  DatatypeP blackboxRepI<T>::getType()
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    return 0 ;
  }
  //**************************************************************************/
    
    
  template<class T> 
  frame_info blackboxRepI<T>::get_frame_info()
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    frame_info fi ;
    return fi ;
  }
    
  //**************************************************************************/
    
  template<class T> class blackbox : public store_instance {
    typedef blackboxRepI<T> blackboxType;
    T * data;
    blackbox(const blackbox &var) { setRep(var.Rep()); }
    blackbox & operator=(const blackbox &p) {setRep(p.Rep()); return *this; }
  public:
    typedef T containerType;
    blackbox() { setRep(new blackboxType); }

    blackbox(storeRepP rp) { setRep(rp); }

    virtual ~blackbox();


    blackbox & operator=(storeRepP p) {setRep(p); return *this; }
    blackbox & operator=(const T &v) { *data = v; return *this; }

    virtual void notification();
    
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

    void set_entitySet(const entitySet &ptn) {Rep()->allocate(ptn); }

    entitySet domain() const { return Rep()->domain(); }

    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s); }
  };

  //**************************************************************************/

  template<class T> blackbox<T>::~blackbox() {}
    
  //**************************************************************************/

  template<class T> 
  void blackbox<T>::notification()
  {
    data = 0 ;
    NPTR<blackboxType> p(Rep());
    if(p!=0)
      data = p->get_blackbox();
    warn(p==0);
  }

  //**************************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const blackbox<T> &t)
  {
    return t.Print(s);
  }

  //**************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, blackbox<T> &t)
  {
    return t.Input(s);
  }

  //**************************************************************************/

  template<class T> 
  class const_blackbox : public store_instance {
    typedef T containerType;
    typedef blackboxRepI<T> blackboxType;
    const T * data;
    const_blackbox(const const_blackbox<T> &var) { setRep(var.Rep()); }
    const_blackbox(const blackbox<T> &var) { setRep(var.Rep()); }
    const_blackbox & operator=(const const_blackbox<T> &p)
      { setRep(p.Rep()); return *this;}
    const_blackbox & operator=(const blackbox<T> &p)
      { setRep(p.Rep()); return *this;}
  public:
    const_blackbox() { setRep(new blackboxType); }
    const_blackbox(storeRepP rp) { setRep(rp); }
    
    virtual ~const_blackbox();

    const_blackbox & operator=(storeRepP p)
    { setRep(p); return *this;}

    virtual void notification();
    virtual instance_type access() const;
        
    const T * operator->() const { return data; }
    
    const T &operator*() const { return *data; }

    const T &operator[](int indx) const {
#ifdef BOUNDS_CHECK
      fatal(data == NULL);
      fatal(!Rep()->domain().inSet(indx));
#endif
      return *data;
    }

    entitySet domain() const { return Rep()->domain(); }
    
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
  };

  //**************************************************************************/

  template<class T> const_blackbox<T>::~const_blackbox() {}

  //**************************************************************************/

  template<class T> 
  void const_blackbox<T>::notification() 
  {  
    NPTR<blackboxType> p(Rep());
    if(p!=0) data = p->get_blackbox();
    warn(p==0);
  }
    
  //**************************************************************************/

  template<class T> 
  storeRepP blackboxRepI<T>::remap(const dMap &m) const 
  {
    blackbox<T> r;
    r.set_entitySet(m.image(m.domain()&domain()));
    if(attrib_data!=0)
      *r = *attrib_data;
    return r.Rep();
  }

  template<class T>
  storeRepP blackboxRepI<T>::freeze() {
    return getRep() ;
  }

  template<class T>
  storeRepP blackboxRepI<T>::thaw() {
    return getRep() ;
  }

  //**************************************************************************/

  template<class T> 
  void blackboxRepI<T>::copy(storeRepP &st, const entitySet &context) 
  {
    warn((store_domain - context) != EMPTY);
    blackbox<T> p(st);
    allocate(context) ;
    if(attrib_data != 0)
      *attrib_data = *p ;
    dispatch_notify();
  }

  //**************************************************************************/

  template<class T> 
  void blackboxRepI<T>::gather(const dMap &m, storeRepP &st,
                               const entitySet &context) 
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }

  //**************************************************************************/

  template<class T> 
  void blackboxRepI<T>::scatter(const dMap &m, storeRepP &st,
                                const entitySet &context) 
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }

  //**************************************************************************/
 
  template <class T> 
  int blackboxRepI<T>::pack_size(const entitySet &eset) 
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    return 0;
  }
  
  template <class T> 
  int blackboxRepI<T>::estimated_pack_size(const entitySet &eset) 
  {
    
    return 0;
  }
  
  template <class T> 
  int blackboxRepI<T>::pack_size(const entitySet &e, entitySet& packed) 
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
    return 0;
  }

  //**************************************************************************/

  template <class T> 
  void blackboxRepI<T>::pack(void *ptr, int &loc, int &size,
                             const entitySet &e) 
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }

  //**************************************************************************/

  template <class T> 
  void blackboxRepI<T>::unpack(void *ptr, int &loc, int &size,
                               const sequence &seq)
  {
    cerr << "BLACKBOX " << __FILE__ << "(" << __LINE__ << ")" << endl;
  }  


  //**************************************************************************/

  template<class T>
  store_instance::instance_type const_blackbox<T>::access() const
  {
    return READ_ONLY;
  }

  //**************************************************************************/
    
  template<class T> 
  inline std::ostream & operator<<(std::ostream &s,
                                   const const_blackbox<T> &t)
  {
    return t.Print(s);
  }

  //**************************************************************************/
}

#endif
    
