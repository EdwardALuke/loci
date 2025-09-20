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
#ifndef DSTOREMAT_H
#define DSTOREMAT_H 

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <istream>
#include <ostream>

#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/tools.h>

#include <store_rep.h>

#include <Tools/lmutex.h>

#include <Map.h>
#include <multiMap.h>


namespace Loci {
  template<class T> class dstoreMat : public store_instance {
    typedef storeVecRepI<T> storeType ;
    T* base_ptr ;
    int size_tot ;
    int size_dim ;
    dstoreMat(const storeMat &var) {setRep(var.Rep()) ; }
    dstoreMat<T> & operator=(const dstoreMat<T> &str)
      { setRep(str.Rep()) ; return *this ;}

  public:
    typedef Mat<T> containerType ;
    dstoreMat() {setRep(new storeType) ; Rep()->setIsMat(true);}
    dstoreMat(storeRepP &rp) {setRep(rp) ; Rep()->setIsMat(true); }

    virtual ~dstoreMat() ;
    virtual void notification() ;

    dstoreMat<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void setVecSize(int size) {
      size_dim = size ;
      size_tot = size*size ;
      Rep()->set_elem_size(size_tot) ; 
    }

    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }

    int vecSize() const { return size_dim; }

    const entitySet &domain() const { return Rep()->domain() ; }

    //    operator storeRepP() { return Rep() ; }

    Mat<T> elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Mat<T>(base_ptr+(indx*size_tot),size_dim) ; }

    Mat<T> operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Mat<T>(base_ptr+(indx*size_tot),size_dim) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}

  } ;

  template<class T> dstoreMat<T>::~dstoreMat<T>() { }

  template<class T> void dstoreMat<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      base_ptr = p->get_base_ptr() ;
      size_tot = p->get_size() ;
      size_dim = int(sqrt(float(size_tot))+.5) ;
    }
    warn(p == 0) ;
  }

  template<class T> inline std::ostream & operator<<(std::ostream &s, const dstoreMat<T> &t)
    { return t.Print(s) ; }

  template<class T> inline std::istream & operator>>(std::istream &s, dstoreMat<T> &t)
    { return t.Input(s) ; }

  template<class T> class const_dstoreMat : public store_instance {
    typedef storeVecRepI<T> storeType ;
    const T* base_ptr ;
    int size_tot ;
    int size_dim ;
    const_dstoreMat(const const_dstoreMat<T> &var) { setRep(var.Rep()) ; }
    const_dstoreMat(const dstoreMat<T> &var) { setRep(var.Rep()) ; }
    const_dstoreMat<T> & operator=(const const_dstoreMat<T> &str)
    { setRep(str.Rep()) ; return *this ;}

    const_dstoreMat<T> & operator=(const dstoreMat<T> &str)
    { setRep(str.Rep()) ; return *this ;}

  public:
    typedef const_Mat<T> containerType ;
    const_dstoreMat() { setRep(new storeType) ;  Rep()->setIsMat(true);}

    const_dstoreMat(storeRepP &rp) { setRep(rp) ; Rep()->setIsMat(true); }
    
    virtual ~const_dstoreMat() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_dstoreMat<T> & operator=(storeRepP p)
    { setRep(p) ; return *this ; }

    int vecSize() const { return size_dim; }
    const entitySet &domain() { return Rep()->domain() ; }

    const_Mat<T> elem(int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Mat<T>(base_ptr+(indx*size_tot),size_dim) ; }
    const_Mat<T> operator[](int indx) {
#ifdef BOUNDS_CHECK
      fatal(base_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Mat<T>(base_ptr+(indx*size_tot),size_dim) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
  } ;

  template<class T> const_dstoreMat<T>::~const_dstoreMat<T>() { }

  template<class T> void const_dstoreMat<T>::notification() {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      base_ptr = p->get_base_ptr() ;
      size_tot = p->get_size() ;
      size_dim = int(sqrt(float(size_tot))+.5) ;
    }
    warn(p == 0) ;
  }

  template<class T> store_instance::instance_type
    const_dstoreMat<T>::access() const
    { return READ_ONLY; }
        

  template<class T> inline std::ostream &
    operator<<(std::ostream &s, const const_dstoreMat<T> &t)
    { return t.Print(s) ; }

#endif

