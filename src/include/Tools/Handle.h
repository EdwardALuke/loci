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
#ifndef HANDLE_H
#define HANDLE_H (1)

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <typeinfo>

#include <Tools/debug.h>

#include <Tools/lmutex.h>

namespace Loci {
    
  enum Handle_type {NULL_HANDLE} ;

  // Counted pointer template
  template<class T> class HandleGrab {
  private:
    // Reference count, must be > 0 at all times
    int count ;
    // Dont allow grabs to be assigned
    HandleGrab &operator=(const HandleGrab &ptn) 
    { warn(true); GrabItem = ptn.GrabItem ; count = 1 ; return *this ; }
  public:
    HandleGrab() { count = 1 ; }
    HandleGrab(const HandleGrab &g) : GrabItem(g.GrabItem) { count = 1 ; }
#ifndef SUN_CC
    template<class S> explicit HandleGrab(S arg) : GrabItem(arg)
    { count = 1 ; }
#endif
    ~HandleGrab() { fatal(count!=0) ; }
    // add reference to grab
    HandleGrab *LinkGrab()
    {lock.lock(); ++count ; lock.unlock(); return this ;}
    // delete reference to grab
    void UnlinkGrab() // { if(--count==0) delete this; }
    {
      lock.lock(); 
      if(--count == 0) {
        lock.unlock();
        delete this ; 
      } else
        lock.unlock();
    }
    // check if there is only one reference to grab
    bool GrabUnique() const { return count == 1 ; }
    
    T GrabItem ;
    lmutex lock;
  } ;

  template<class T> class ConstHandle ;

  template<class T> class Handle {
  private:
    void ChangeGrab(HandleGrab<T> *gp) {
      if(grab!=gp) {
        HandleGrab<T> *tmp = gp->LinkGrab() ;
        if(grab!=0)
          grab->UnlinkGrab() ;
        grab = tmp ;
      }
    }
  protected:
    HandleGrab<T> *grab ;
  public:
    friend class ConstHandle<T> ;
    // OCF Methods
    Handle() {grab = new HandleGrab<T> ; }
    Handle(Handle_type i) { grab = 0 ; }
#ifndef SUN_CC
    template<class S> explicit Handle(S arg) { grab = new HandleGrab<T>(arg); }
#endif
    Handle(const Handle<T> & hdl)
    { grab = 0 ; if(hdl.grab != 0) grab = hdl.grab->LinkGrab() ;  }
    ~Handle() { if(grab!=0) grab->UnlinkGrab() ; }
    Handle<T> &operator=(const Handle<T> &hr) 
    { ChangeGrab(hr.grab) ; return *this; }
    Handle<T> &operator=(Handle_type i) {
      if(grab!=0) grab->UnlinkGrab() ;
      grab = 0 ;
      return *this;
    }
    bool null() const  { return grab == 0; }
    // Make Handle point to a newly allocated grab
    void New() { HandleGrab<T> *ng = new HandleGrab<T> ;
    ChangeGrab(ng) ; ng->UnlinkGrab() ; }

    // Make Handle point to a uniquely allocated grab
    void MakeUnique() {  if((grab==0)||!grab->GrabUnique()) {
      HandleGrab<T> *ng = new HandleGrab<T>(*grab) ; 
      ChangeGrab(ng) ; ng->UnlinkGrab() ; }}

    // Dereference operator
    T &operator* () { return grab->GrabItem ; }
    const T &operator* () const { return grab->GrabItem ; }

    // Delegation operator
    T * operator-> () { return &grab->GrabItem; }
    const T * operator-> () const { return (const T *)&grab->GrabItem; }

  } ;


  // Handles inherit the relational operators from their template type
  template<class T> bool operator<(const Handle<T> &h1, const Handle<T> &h2)
    { return (*h1 < *h2) ; }
  template<class T> bool operator<=(const Handle<T> &h1, const Handle<T> &h2)
    { return (*h1 <= *h2) ; }
  template<class T> bool operator>(const Handle<T> &h1, const Handle<T> &h2)
    { return (*h1 > *h2) ; }
  template<class T> bool operator>=(const Handle<T> &h1, const Handle<T> &h2)
    { return (*h1 >= *h2) ; }
  template<class T> bool operator==(const Handle<T> &h1, const Handle<T> &h2)
    { return (*h1 == *h2) ; }
  template<class T> bool operator!=(const Handle<T> &h1, const Handle<T> &h2)
    { return (*h1 != *h2) ; }

  template<class T> bool operator==(const Handle<T> &h, Handle_type t)
    { return h.null() ; }
  template<class T> bool operator==(Handle_type t, const Handle<T> &h)
    { return h.null() ; }

  template<class T> bool operator!=(const Handle<T> &h, Handle_type t)
    { return !h.null() ; }
  template<class T> bool operator!=(Handle_type t, const Handle<T> &h)
    { return !h.null() ; }


  template<class T> class ConstHandle {
  private:

    HandleGrab<T> *grab ;

    void ChangeGrab(HandleGrab<T> *gp) {
      if(grab!=gp) {
        HandleGrab<T> *tmp = gp->LinkGrab() ;
        if(grab!=0)
          grab->UnlinkGrab() ;
        grab = tmp ;
      }
    }
    ConstHandle() {grab = new HandleGrab<T> ; }
  public:
    // OCF Methods
    ConstHandle(const Handle<T> & hdl) { grab = hdl.grab->LinkGrab() ; }
    ConstHandle(const ConstHandle<T> & hdl) { grab = hdl.grab->LinkGrab() ; }
    ConstHandle(Handle_type i) { grab = 0 ; }
    ~ConstHandle() {
      if(grab != 0) grab->UnlinkGrab() ; }
    ConstHandle<T> &operator=(const Handle<T> &hr) 
    { ChangeGrab(hr.grab) ; return *this; }
    ConstHandle<T> &operator=(const ConstHandle<T> &hr) 
    { ChangeGrab(hr.grab) ; return *this; }

    bool null() { return grab == 0 ; }

    // Dereference operator
    const T &operator* () const { return grab->GrabItem ; }

    // Delegation operator
    const T * operator-> () const { return (const T *)&grab->GrabItem; }

  } ;

  // ConstHandles inherit the relational operators from their template type
  template<class T> bool operator<(const ConstHandle<T> &h1,
                                   const ConstHandle<T> &h2)
    { return (*h1 < *h2) ; }
  template<class T> bool operator<=(const ConstHandle<T> &h1,
                                    const ConstHandle<T> &h2)
    { return (*h1 <= *h2) ; }
  template<class T> bool operator>(const ConstHandle<T> &h1,
                                   const ConstHandle<T> &h2)
    { return (*h1 > *h2) ; }
  template<class T> bool operator>=(const ConstHandle<T> &h1,
                                    const ConstHandle<T> &h2)
    { return (*h1 >= *h2) ; }
  template<class T> bool operator==(const ConstHandle<T> &h1,
                                    const ConstHandle<T> &h2)
    { return (*h1 == *h2) ; }
  template<class T> bool operator!=(const ConstHandle<T> &h1,
                                    const ConstHandle<T> &h2)
    { return (*h1 != *h2) ; }

  template<class T> bool operator==(const ConstHandle<T> &h, Handle_type t)
    { return h.null() ; }
  template<class T> bool operator==(Handle_type t, const ConstHandle<T> &h)
    { return h.null() ; }

  template<class T> bool operator!=(const ConstHandle<T> &h, Handle_type t)
    { return !h.null() ; }
  template<class T> bool operator!=(Handle_type t, const ConstHandle<T> &h)
    { return !h.null() ; }

}

#endif
