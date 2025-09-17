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
#ifndef STOREMAT_H
#define STOREMAT_H 

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
#include <matrix.h>

namespace Loci {
 
  template<class T> class storeMat : public store_instance {
    typedef storeVecRepI<T> storeType ;
    T* alloc_ptr ;
    int base_offset ;
    int size_tot ;
    int size_dim ;
    storeMat(const storeMat &var) {setRep(var.Rep()) ; }
    storeMat<T> & operator=(const storeMat<T> &str)
    { setRep(str.Rep()) ; return *this ;}

  public:
    typedef Mat<T> containerType ;
    storeMat() {setRep(new storeType) ; Rep()->setIsMat(true);}
    storeMat(storeRepP rp) {setRep(rp); Rep()->setIsMat(true);}

    virtual ~storeMat() ;
    virtual void notification() ;

    storeMat<T> & operator=(storeRepP p) { setRep(p) ; return *this ; }

    void setVecSize(int size) {
      size_dim = size ;
      size_tot = size*size ;
      Rep()->set_elem_size(size_tot) ; 
    }

    void allocate(const entitySet &ptn) { Rep()->allocate(ptn) ; }
    int vecSize() const { return size_dim; }
    const entitySet domain() const { return Rep()->domain() ; }

    Mat<T> elem(Entity indx) {
#ifdef BOUNDS_CHECK
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Mat<T>(alloc_ptr+((indx-base_offset)*size_tot),size_dim) ; }
    Mat<T> operator[](Entity indx) {
#ifdef BOUNDS_CHECK
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Mat<T>(alloc_ptr+((indx-base_offset)*size_tot),size_dim) ; }
    Mat<T> operator[](size_t indx) {
#ifdef BOUNDS_CHECK
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return Mat<T>(alloc_ptr+((indx-base_offset)*size_tot),size_dim) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }
    std::istream &Input(std::istream &s) { return Rep()->Input(s) ;}
  } ;

  //**************************************************************************/

  template<class T> 
  storeMat<T>::~storeMat() { }

  //**************************************************************************/

  template<class T> 
  void storeMat<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      alloc_ptr = p->get_alloc_ptr() ;
      base_offset = p->get_base_offset() ;
      size_tot = p->get_size() ;
      size_dim = int(sqrt(float(size_tot))+.5) ;
    }
    warn(p == 0) ;
  }

  //**************************************************************************/

  template<class T> 
  inline std::ostream & operator<<(std::ostream &s, const storeMat<T> &t)
  { return t.Print(s) ; }

  //**************************************************************************/

  template<class T> 
  inline std::istream & operator>>(std::istream &s, storeMat<T> &t)
  { return t.Input(s) ; }

  //*************************************************************************/

  template<class T> class const_storeMat : public store_instance {
    typedef storeVecRepI<T> storeType ;
    const T* alloc_ptr ;
    int base_offset ;
    int size_tot ;
    int size_dim ;
    const_storeMat(const const_storeMat<T> &var) { setRep(var.Rep()) ; }
    const_storeMat(const storeMat<T> &var) { setRep(var.Rep()) ; }
    const_storeMat<T> & operator=(const const_storeMat<T> &str)
    { setRep(str.Rep()) ; return *this ;}
    const_storeMat<T> & operator=(const storeMat<T> &str)
    { setRep(str.Rep()) ; return *this ;}

  public:
    typedef const_Mat<T> containerType ;
    const_storeMat() { setRep(new storeType) ; Rep()->setIsMat(true); }

    const_storeMat(storeRepP rp) { setRep(rp) ; Rep()->setIsMat(true); }
    
    virtual ~const_storeMat() ;
    virtual void notification() ;

    virtual instance_type access() const ;
        
    const_storeMat<T> & operator=(storeRepP p)
    { setRep(p) ; return *this ; }

    int vecSize() const { return size_dim; }
    const entitySet domain() { return Rep()->domain() ; }
    //    operator storeRepP() { return Rep() ; }

    const_Mat<T> elem(Entity indx) {
#ifdef BOUNDS_CHECK
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Mat<T>(alloc_ptr+((indx-base_offset)*size_tot),size_dim) ; }
    const_Mat<T> operator[](Entity indx) {
#ifdef BOUNDS_CHECK
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Mat<T>(alloc_ptr+((indx-base_offset)*size_tot),size_dim) ; }
    const_Mat<T> operator[](size_t indx) {
#ifdef BOUNDS_CHECK
      fatal(alloc_ptr==NULL); 
      fatal(!((Rep()->domain()).inSet(indx))) ;
#endif 
      return const_Mat<T>(alloc_ptr+((indx-base_offset)*size_tot),size_dim) ; }
    std::ostream &Print(std::ostream &s) const { return Rep()->Print(s); }

  } ;

  //***************************************************************************/

  template<class T> 
  const_storeMat<T>::~const_storeMat() { }

  //***************************************************************************/

  template<class T> 
  void const_storeMat<T>::notification() 
  {
    NPTR<storeType> p(Rep()) ;
    if(p!=0) {
      alloc_ptr = p->get_alloc_ptr() ;
      base_offset = p->get_base_offset() ;
      size_tot = p->get_size() ;
      size_dim = int(sqrt(float(size_tot))+.5) ;
    }
    warn(p == 0) ;
  }

  //***************************************************************************/

  template<class T> 
  store_instance::instance_type
  const_storeMat<T>::access() const
  { return READ_ONLY; }
        
  //***************************************************************************/

  template<class T> inline std::ostream &
  operator<<(std::ostream &s, const const_storeMat<T> &t)
  { return t.Print(s) ; }

  //***************************************************************************/

}

#endif
