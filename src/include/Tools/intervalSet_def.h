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
#ifndef INTERVALSET_DEF_H
#define INTERVALSET_DEF_H 1

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <cstdint>
#include <Tools/debug.h>
#include <Tools/Handle.h>
#include <Tools/tools.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <iostream>

//this directive should be paired with the typedef of GEntity
#define MPI_GENTITY_TYPE MPI_INT
#define HDF5_GENTITY_TYPE  H5T_NATIVE_INT

namespace Loci {

  typedef int int_type ;
#define LONG_GENTITY

#ifdef LONG_GENTITY
  typedef std::int64_t gEntity;
#else
  typedef int_type gEntity;
#endif

   
  //c++ can not typedef a template class, so interval is replaced by genInterval and moved into genIntervalSet class
  //outside class definition, interval is replaced by std::pair<T, T>. same with pair_vector
  

  //num_intervals return size_t, instead of int
  //use size_t for size of genIntervalSet
  //use size_t for index of genIntervalSet
  //create_intervalSet, create_sequence and FORALL macro are for int_type only
  //GFORALL is added for GEntity
  //add comparison between sequence and intervalSet
  
  template<class T> std::ostream & operator<<(std::ostream &s, const std::pair<T, T> &i) ;
  template<class T> std::istream & operator>>(std::istream &s, std::pair<T, T> &i) ;
 
  template<class T> class genIntervalSet ;
  template<class T> class genSequence ;

  
#ifdef ENTITY
  template<class T> void Union_inplace (Handle<std::vector<std::pair<T, T> > >  &Rep,
                                        const std::pair<T,T> &ivl);
  template<class T> Handle<std::vector<std::pair<T,T> > >
  Union(const Handle<std::vector<std::pair<T,T> > >& Rep1,
        const Handle<std::vector<std::pair<T,T> > >& Rep2);
#endif
        
  template<class T> size_t num_intervals(Handle<std::vector<std::pair<T, T> > > &Rep){
    return Rep->size();
  } 
        
  template<typename T>  class genIntervalSet {
  public:
    
    typedef typename std::pair<T, T> genInterval;
    typedef typename std::vector<std::pair<T,T> > pair_vector ;
    typedef typename std::vector<std::pair<T,T> > genIntervalSetRep ;
    static T  UNIVERSE_MAX; 
    static T  UNIVERSE_MIN;
    static genIntervalSet<T> EMPTY;
  private:
    friend class genSequence<T> ;
    Handle<genIntervalSetRep> Rep ;
    
    struct rep_holder {
      Handle<genIntervalSetRep> Rep ;
    } ;
    
    static rep_holder *rhp ; //the static member of template class potentially dangerous?
    static Handle<genIntervalSetRep> & getEmptyRep() {
      if(rhp==0)
        rhp = new rep_holder ;
      return rhp->Rep ;
    } ;
  public:
    class genIntervalSetIterator {//nesting make iterator aumotically template
      friend class genIntervalSet; //declare friend class
      typename genIntervalSetRep::const_iterator current_interval,end_interval ;//must add class or typename for iterator in template class
      T current_value ;
      void increment() {
        if(++current_value>current_interval->second) {
          ++current_interval ;
          if(current_interval!=end_interval)
            current_value = current_interval->first ;
          else
            current_value = UNIVERSE_MAX ;
        }
      }                
    public:
      genIntervalSetIterator() {}
     
      genIntervalSetIterator(typename genIntervalSetRep::const_iterator r,
                             typename genIntervalSetRep::const_iterator end,T e) :
        current_interval(r),end_interval(end),current_value(e) {}
      T operator*() const { return current_value; }
      genIntervalSetIterator &operator++() {
        increment() ;
        return *this ;
      }
      genIntervalSetIterator operator++(int ) {
        genIntervalSetIterator tmp = *this ;
        increment() ;
        return tmp ;
      }
      bool operator==(const genIntervalSetIterator &i) const {
        return current_value == i.current_value ; }
      bool operator!=(const genIntervalSetIterator &i) const {
        return current_value != i.current_value ; }
    } ;


    genIntervalSet() : Rep(getEmptyRep()) {} ;
    genIntervalSet(const genInterval &ivl) { 
      genInterval i(min(ivl.first,ivl.second),
                    max(ivl.first,ivl.second));
      Rep->push_back(i) ; }
    genIntervalSet(const genIntervalSet<T> &ptn): Rep(ptn.Rep) {}
    explicit genIntervalSet(const genSequence<T> &seq) ;
#ifdef ENTITY
    explicit genIntervalSet(const Handle<pair_vector> &RepIn): Rep(RepIn) {}
#endif
    ~genIntervalSet() {}
    
    typedef genIntervalSetIterator const_iterator ;

    const_iterator begin() const {
      if(Rep->size() != 0)
        return genIntervalSetIterator(Rep->begin(),Rep->end(),(*Rep)[0].first) ;
      else
        return end() ;
    }
    const_iterator end() const {
      return genIntervalSetIterator(Rep->end(),Rep->end(),UNIVERSE_MAX) ;
    }
    

    bool inSet(T indx) const ;

    genIntervalSet & operator=(const genIntervalSet &ptn)  
    { Rep = ptn.Rep ; return *this;}
    

    bool Equal(const genIntervalSet &ptn) const
    { return *(Rep) == *(ptn.Rep) ; }
   

    bool less_than(const genIntervalSet &ptn) const
    { return (*Rep) < *(ptn.Rep) ;}
    bool greater_than(const genIntervalSet &ptn) const
    {return (*Rep) > *(ptn.Rep) ;}
    
    genIntervalSet & operator>>=(T rotval) ;
    genIntervalSet & operator<<=(T rotval) ;

    genIntervalSet operator<<(T rotval) const ;
    genIntervalSet operator>>(T rotval) const ;
   
    //use size_t for size of genIntervalSet?
    size_t size() const {
      size_t size = 0 ;
      typename genIntervalSetRep::const_iterator i ;
      for(i = Rep->begin();i!= Rep->end();++i)
        size += i->second-i->first + 1 ;
      return size ;  }
    //use size_t for num_intervals
    size_t num_intervals() const { return Rep->size() ; }
    //use size_t for index of genIntervalSet?
    const genInterval & operator[](size_t indx) const { 
      fatal(indx>=num_intervals()) ;
      return (*Rep)[indx]; 
    }

#ifdef ENTITY

      
    void Union(const genInterval &ivl){Loci::Union_inplace(Rep,ivl);}

    void Union(const genIntervalSet &ptn){
      size_t psz = ptn.num_intervals() ;
      if(psz == 0)// If ptn == EMPTY, do nothing
        return ;
      if(psz == 1) // only one genInterval, use genInterval update in place
        Loci::Union_inplace(Rep,ptn[0]);
      else
        Rep = Loci::Union(Rep,ptn.Rep) ; // General case
    }
    
    static genIntervalSet Union(const genIntervalSet &set1,
                                const genIntervalSet &set2){
      return genIntervalSet(Loci::Union(set1.Rep,set2.Rep)) ;
    }

#else

    void Union(const genInterval &ivl) ;
    void Union(const genIntervalSet &ptn) ;
    static genIntervalSet Union(const genIntervalSet &set1,
                                const genIntervalSet &set2) ;


#endif


    void Intersection(const genInterval &ivl) ;
    void Intersection(const genIntervalSet &ptn) ;
    static genIntervalSet Intersection(const genIntervalSet &set1,
                                       const genIntervalSet &set2) ;
    void Complement() ;
    static genIntervalSet Complement(const genIntervalSet &set) ;
        
    std::ostream &Print(std::ostream &s) const ;
    void Print() const { Print(std::cout) ; }

    std::istream &Input(std::istream &s) ;

    T Min() const {
      size_t sz = Rep->size();
      return sz==0?UNIVERSE_MAX:(*Rep)[0].first ;}
    T Max() const {
      size_t sz = Rep->size() ;
      return sz==0?UNIVERSE_MIN:(*Rep)[sz-1].second ;}
    void make_unique() { Rep.MakeUnique(); }
  } ;


  template<typename T> T genIntervalSet<T>::UNIVERSE_MAX = std::numeric_limits<T>::max() - 1 ;
  template<typename T> T genIntervalSet<T>::UNIVERSE_MIN = std::numeric_limits<T>::min() + 1 ;
  template<typename T> genIntervalSet<T> genIntervalSet<T>::EMPTY = genIntervalSet<T>() ;

  template<typename T> inline genIntervalSet<T> operator~(const genIntervalSet<T> &e) {
    return genIntervalSet<T>::Complement(e) ;
  }

  template<typename T> inline genIntervalSet<T> & operator+=(genIntervalSet<T> &e, T ival) 
  { e.Union(std::pair<T, T>(ival,ival)) ; return e ; }

  template<typename T> inline genIntervalSet<T> & operator+=(genIntervalSet<T> &e, const std::pair<T, T> &ivl)
  { e.Union(ivl) ; return e; }

  template<typename T> inline genIntervalSet<T> & operator+=(genIntervalSet<T> &e, const genIntervalSet<T> &ptn)
  { e.Union(ptn) ; return e ; }

  template<typename T> inline genIntervalSet<T> & operator|=(genIntervalSet<T> &e, T ival)
  { e.Union(std::pair<T, T>(ival,ival)) ; return e ; }

  template<typename T> inline genIntervalSet<T> & operator|=(genIntervalSet<T> &e, const std::pair<T, T> &ivl)
  { e.Union(ivl) ; return e ; }

  template<typename T> inline genIntervalSet<T> & operator|=(genIntervalSet<T> &e, const genIntervalSet<T> &ptn)
  { e.Union(ptn) ; return e ; }

  template<typename T> inline genIntervalSet<T> & operator&=(genIntervalSet<T> &e, T ival)
  { e.Intersection(std::pair<T, T>(ival,ival)) ;
    return e; }

  template<typename T> inline genIntervalSet<T> & operator&=(genIntervalSet<T> &e, const std::pair<T, T> &ivl)
  { e.Intersection(ivl) ; return e ; }

  template<typename T> inline genIntervalSet<T> & operator&=(genIntervalSet<T> &e, const genIntervalSet<T> &ptn)
  { e.Intersection(ptn) ; return e ; }

  template<typename T> inline genIntervalSet<T> operator+(const genIntervalSet<T> & e, T ival) {
    genIntervalSet<T> retp = e ;
    retp.Union(genInterval(ival,ival)) ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator+(const T &ival, const genIntervalSet<T> &e) {
    genIntervalSet<T> retp = e ;
    retp.Union(genInterval(ival,ival)) ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator+(const genIntervalSet<T> &e, const std::pair<T, T> &ivl) {
    genIntervalSet<T> retp = e ;
    retp.Union(ivl) ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator+(const std::pair<T, T> &ivl, const genIntervalSet<T> &e) {
    genIntervalSet<T> retp = e ;
    retp.Union(ivl) ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator+(const genIntervalSet<T> &e, const genIntervalSet<T> &ptn) {
    genIntervalSet<T> retp = e ;
    retp.Union(ptn) ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator|(const genIntervalSet<T> &e, T ival) {
    genIntervalSet<T> retp = e ;
    retp.Union(genInterval(ival,ival)) ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator|(T ival, const genIntervalSet<T> &e) {
    genIntervalSet<T> retp = e ;
    retp.Union(genInterval(ival,ival)) ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator|(const genIntervalSet<T> &e, const std::pair<T, T> &ivl) {
    genIntervalSet<T> retp = e ;
    retp.Union(ivl) ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator|(const std::pair<T,T> &ivl,const genIntervalSet<T> &e) {
    genIntervalSet<T> retp = e ;
    retp.Union(ivl) ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator|(const genIntervalSet<T> &e, const genIntervalSet<T> &ptn) {
    genIntervalSet<T> retp = e ;
    retp.Union(ptn) ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator&(const genIntervalSet<T> &e, T ival)  {
    return genIntervalSet<T>::Intersection(e,genIntervalSet<T>(std::pair<T, T>(ival,ival))) ;
  }

  template<typename T> inline genIntervalSet<T> operator&(T ival,const genIntervalSet<T> &e)  {
    return genIntervalSet<T>::Intersection(e,genIntervalSet<T>(std::pair<T, T>(ival,ival))) ;
  }

  template<typename T> inline genIntervalSet<T> operator&(const genIntervalSet<T> &e, const std::pair<T, T> &ivl) {
    return genIntervalSet<T>::Intersection(e,genIntervalSet<T>(ivl)) ;
  }

  template<typename T> inline genIntervalSet<T> operator&(const std::pair<T, T> &ivl,const genIntervalSet<T> &e) {
    return genIntervalSet<T>::Intersection(e,genIntervalSet<T>(ivl)) ;
  }

  template<typename T>  inline genIntervalSet<T> operator&(const genIntervalSet<T> &e, const genIntervalSet<T> &ptn) {
    return genIntervalSet<T>::Intersection(e,ptn) ;
  }

  template<typename T> inline genIntervalSet<T> operator-(const genIntervalSet<T> &e, T &ival) {
    return genIntervalSet<T>::Intersection(e,~genIntervalSet<T>(std::pair<T, T>(ival,ival))) ;
  }

  template<typename T> inline genIntervalSet<T> operator-(const genIntervalSet<T> &e, const std::pair<T, T> &ivl) {
    return genIntervalSet<T>::Intersection(e,~genIntervalSet<T>(ivl)) ;
  }

  template<typename T>  inline genIntervalSet<T> operator-(const std::pair<T, T> &ivl, const genIntervalSet<T> &e) {
    return genIntervalSet<T>::Intersection(genIntervalSet<T>(ivl),~e) ;
  }


  template<typename T> inline genIntervalSet<T> operator-(const genIntervalSet<T> &e, const genIntervalSet<T> &ptn) {
    return genIntervalSet<T>::Intersection(e,~ptn) ;
  }

  template<typename T> inline genIntervalSet<T> operator^(const genIntervalSet<T> &e, const genIntervalSet<T> &ptn) {
    genIntervalSet<T> retp = ((e|ptn) & ~(e & ptn)) ;
    return retp ;
  }
  template<typename T> inline genIntervalSet<T> operator^(const genIntervalSet<T> &e, T ival) {
    genIntervalSet<T> ivlp = genInterval(ival,ival) ;
    genIntervalSet<T> retp = e^ivlp ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator^(T ival, const genIntervalSet<T> &e) {
    genIntervalSet<T> ivlp = genInterval(ival,ival) ;
    genIntervalSet<T> retp = e^ivlp ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator^(const genIntervalSet<T> &e, const std::pair<T, T> &ivl) {
    genIntervalSet<T> ivlp = ivl ;
    genIntervalSet<T> retp = e^ivlp ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> operator^(const std::pair<T, T> &ivl, const genIntervalSet<T> &e) {
    genIntervalSet<T> ivlp = ivl ;
    genIntervalSet<T> retp = e^ivlp ;
    return retp ;
  }


  template<typename T> inline bool operator==(const genIntervalSet<T> &es1, const genIntervalSet<T> &es2) {
    return es1.Equal(es2) ;
  }

  template<typename T> inline bool operator!=(const genIntervalSet<T> &es1, const genIntervalSet<T> &es2) {
    return !es1.Equal(es2) ;
  }

  template<typename T> inline bool operator<(const genIntervalSet<T> &es1, const genIntervalSet<T> &es2) {
    return es1.less_than(es2) ;
  }

  template<typename T> inline bool operator>(const genIntervalSet<T> &es1, const genIntervalSet<T> &es2) {
    return es1.greater_than(es2) ;
  }

  template<typename T> inline bool operator<=(const genIntervalSet<T> &es1, const genIntervalSet<T> &es2) {
    return !(es1.greater_than(es2)) ;
  }

  template<typename T> inline bool operator>=(const genIntervalSet<T> &es1, const genIntervalSet<T> &es2) {
    return !(es1.less_than(es2)) ;
  }

  template<typename T> inline genIntervalSet<T> genIntervalSet<T>::operator>>(T rotval) const  {
    genIntervalSet<T> retp = (*this) ;
    retp >>= rotval ;
    return retp ;
  }

  template<typename T> inline genIntervalSet<T> genIntervalSet<T>::operator<<(T rotval) const {
    genIntervalSet<T> retp = (*this) ;
    retp <<= rotval ;
    return retp ;
  }


  template<typename T> inline genIntervalSet<T> & operator^=(genIntervalSet<T> &e, const genIntervalSet<T> &ptn) {
    e =  e^ptn ;
    return e ;
  }
    
  template<typename T> inline genIntervalSet<T> & operator^=(genIntervalSet<T> &e, const std::pair<T, T> &ivl) {
    genIntervalSet<T> ivlp ;
    ivlp = ivl ;
    e ^= ivlp ;
    return e ;
  }

  template<typename T> inline genIntervalSet<T> & operator^=(genIntervalSet<T> &e, T ival) {
    genIntervalSet<T> ivlp ;
    ivlp = genInterval(ival,ival) ;
    e ^= ivlp ;
    return e ;
  }

  template<typename T> inline genIntervalSet<T> & operator-=(genIntervalSet<T> &e, const genIntervalSet<T> &ptn) {
    genIntervalSet<T> ptnc = genIntervalSet<T>::Complement(ptn);
    e &= ptnc ;
    return e ;
  }
    
  template<typename T> inline genIntervalSet<T> & operator-=(genIntervalSet<T> &e, const std::pair<T, T> &ivl) {
    genIntervalSet<T> ptnc ;
    ptnc = ivl ;
    ptnc.Complement() ;
    e &= ptnc ;
    return e ;
  }

  template<typename T> inline genIntervalSet<T> & operator-=(genIntervalSet<T> &e, T ival) {
    genIntervalSet<T> ptnc;
    ptnc = std::pair<T, T>(ival,ival) ;
    ptnc.Complement() ;
    e &= ptnc ;
    return e ;
  }

    
  template<typename T> inline std::ostream & operator<<(std::ostream &s, const genIntervalSet<T> &e)
  { return e.Print(s) ; }
    
  template<typename T> inline std::istream & operator>>(std::istream &s, genIntervalSet<T> &e)
  { return e.Input(s) ;}


    


    
  template<typename T> class genSequence {

  public:
    typedef typename genIntervalSet<T>::genIntervalSetRep genSequenceRep ;//typename is mandatory before a qualified , dependent names 
    typedef typename std::pair<T, T> genInterval;
    typedef typename std::vector<std::pair<T,T> > pair_vector ;
  private:
    
    Handle<genSequenceRep> Rep ;

    struct rep_holder {
      Handle<genSequenceRep> Rep ;
    } ;
    static rep_holder *rhp ;
    static Handle<genSequenceRep> & getEmptyRep() {
      if(rhp==0)
        rhp = new rep_holder ;
      return rhp->Rep ;
    } ;
  public:
    class genSequenceIterator {
      typename genSequenceRep::const_iterator current_genInterval,end_genInterval ;
      T current_value ;
      T dir ;
      void increment() {
        if(current_value == current_genInterval->second) {
          ++current_genInterval ;
          if(current_genInterval!=end_genInterval) {
            current_value = current_genInterval->first ;
            dir = current_genInterval->second>=current_value?1:-1 ;
          } else {
            current_value = genIntervalSet<T>::UNIVERSE_MAX ;
            dir = 0 ;
          }
        } else
          current_value += dir ;
      }                
    public:
      genSequenceIterator() {}
      //      { current_genInterval=0;end_genInterval=0;current_value=0; }
      genSequenceIterator(typename genSequenceRep::const_iterator r,
                          typename genSequenceRep::const_iterator end,T e) :
        current_genInterval(r),end_genInterval(end),current_value(e)
      { dir = ((e!=genIntervalSet<T>::UNIVERSE_MAX)&&(r->second>=r->first))?1:-1 ;}
      T operator*() const { return current_value; }
      genSequenceIterator &operator++() {
        increment() ;
        return *this ;
      }
      genSequenceIterator operator++(int ) {
        genSequenceIterator tmp = *this ;
        increment() ;
        return tmp ;
      }
      bool operator==(const genSequenceIterator &i) const {
        return current_value == i.current_value  &&
          current_genInterval == i.current_genInterval; }
      bool operator!=(const genSequenceIterator &i) const 
      { return !operator==(i) ; }
    } ;

    
    genSequence() : Rep(getEmptyRep()) {} ;
    genSequence(const genSequence<T> &seq): Rep(seq.Rep) {}
    genSequence(const genInterval &ivl) {
      genInterval i(min(ivl.first,ivl.second),
                    max(ivl.first,ivl.second));
      Rep->push_back(i) ;
    }
    genSequence(const genIntervalSet<T> &ptn) : Rep(ptn.Rep) {}
    ~genSequence() {}
    
    typedef genSequenceIterator const_iterator ;

    const_iterator begin() const {
      if(Rep->size() != 0)
        return genSequenceIterator(Rep->begin(),Rep->end(),(*Rep)[0].first) ;
      return end() ;
    }
    const_iterator end() const {
      return genSequenceIterator(Rep->end(),Rep->end(),genIntervalSet<T>::UNIVERSE_MAX) ;
    }
    

    genSequence<T> & operator=(const genSequence<T> &ptn)  
    { Rep = ptn.Rep ; return *this;}
    // genSequence<T> & operator=(const genInterval &ivl) 
    //     { Rep.New() ; Rep->push_back(ivl) ; return *this ; }

    bool Equal(const genSequence<T> &ptn) const
    { return *(Rep) == *(ptn.Rep) ; }
    //add comparison between sequence and genIntervalSet 
    bool Equal(const genIntervalSet<T> &ptn) const
    { return *(Rep) == *(ptn.Rep) ; }


    bool less_than(const genSequence<T> &ptn) const
    { return (*Rep) < *(ptn.Rep) ;}
    bool greater_than(const genSequence<T> &ptn) const
    {return (*Rep) > *(ptn.Rep) ;}
    
    genSequence & operator>>=(T rotval) ;
    genSequence & operator<<=(T rotval) ;

    genSequence operator<<(T rotval) const ;
    genSequence operator>>(T rotval) const ;
   
    size_t size() const {
      size_t size = 0 ;
      typename genSequenceRep::const_iterator i ;
      for(i = Rep->begin();i!= Rep->end();++i)
        size += abs(i->first-i->second) + 1 ;
      return size ;  }
    size_t num_intervals() const { return Rep->size();  }

    const genInterval& operator[](size_t indx) const 
    { fatal(indx>=num_intervals()) ;
      return (*Rep)[indx]; }

    void Append(const std::pair<T, T> &ivl) ;
    void Append(const genSequence &seq) ;
    void Append(const genIntervalSet<T> &ptn) ;

    genSequence &Reverse() ;
    std::ostream &Print(std::ostream &s) const ;
    void Print() const { Print(std::cout) ; }

    std::istream &Input(std::istream &s) ;
    void make_unique() { Rep.MakeUnique(); }
  } ;

  template<typename T> inline genSequence<T> & operator+=(genSequence<T> &e, T ival) 
  { e.Append(std::pair<T, T>(ival,ival)) ; return e ; }

  template<typename T> inline genSequence<T> & operator+=(genSequence<T> &e, const std::pair<T, T> &ivl)
  { e.Append(ivl) ; return e; }

  template<typename T> inline genSequence<T> & operator+=(genSequence<T> &e, const genSequence<T> &seq)
  { e.Append(seq) ; return e ; }

  template<typename T> inline genSequence<T> & operator+=(genSequence<T> &e, const genIntervalSet<T> &ptn)
  { e.Append(ptn) ; return e ; }

  template<typename T> inline genSequence<T> operator+(const genSequence<T> & e, T ival) {
    genSequence<T> retp = e ;
    retp.Append(std::pair<T, T>(ival,ival)) ;
    return retp ;
  }

  template<typename T> inline genSequence<T> operator+(const T &ival, const genSequence<T> &e) {
    genSequence<T> retp = std::pair<T, T>(ival,ival) ;
    retp.Append(e) ;
    return retp ;
  }

  template<typename T> inline genSequence<T> operator+(const genSequence<T> &e, const std::pair<T, T> &ivl) {
    genSequence<T> retp = e ;
    retp.Append(ivl) ;
    return retp ;
  }

  template<typename T> inline genSequence<T> operator+(const std::pair<T, T> &ivl, const genSequence<T> &e) {
    genSequence<T> retp = ivl ;
    retp.Append(e) ;
    return retp ;
  }

  template<typename T> inline genSequence<T> operator+(const genSequence<T> &e, const genSequence<T> &ptn) {
    genSequence<T> retp = e ;
    retp.Append(ptn) ;
    return retp ;
  }

  template<typename T> inline genSequence<T> operator+(const genSequence<T> &e, const genIntervalSet<T> &ptn) {
    genSequence<T> retp = e ;
    retp.Append(ptn) ;
    return retp ;
  }

  template<typename T> inline bool operator==(const genSequence<T> &es1, const genSequence<T> &es2) {
    return es1.Equal(es2) ;
  }

  template<typename T> inline bool operator!=(const genSequence<T> &es1, const genSequence<T> &es2) {
    return !es1.Equal(es2) ;
  }

  //add comparison between sequence and genIntervalSet
  template<typename T> inline bool operator==(const genSequence<T> &es1, const genIntervalSet<T> &es2) {
    return es1.Equal(es2) ;
  }
  
  template<typename T> inline bool operator!=(const genSequence<T> &es1, const genIntervalSet<T> &es2) {
    return !es1.Equal(es2) ;
  }

  
  template<typename T> inline bool operator<(const genSequence<T> &es1, const genSequence<T> &es2) {
    return es1.less_than(es2) ;
  }

  template<typename T> inline bool operator>(const genSequence<T> &es1, const genSequence<T> &es2) {
    return es1.greater_than(es2) ;
  }

  template<typename T> inline bool operator<=(const genSequence<T> &es1, const genSequence<T> &es2) {
    return !(es1.greater_than(es2)) ;
  }

  template<typename T> inline bool operator>=(const genSequence<T> &es1, const genSequence<T> &es2) {
    return !(es1.less_than(es2)) ;
  }

  template<typename T> inline genSequence<T> genSequence<T>::operator>>(T rotval) const  {
    genSequence<T> retp = (*this) ;
    retp >>= rotval ;
    return retp ;
  }

  template<typename T> inline genSequence<T> genSequence<T>::operator<<(T rotval) const {
    genSequence<T> retp = (*this) ;
    retp >>= rotval ;
    return retp ;
  }

  template<typename T> inline std::ostream & operator<<(std::ostream &s, const genSequence<T> &e)
  { return e.Print(s) ; }
    
  template<typename T> inline std::istream & operator>>(std::istream &s, genSequence<T> &e)
  { return e.Input(s) ;}

  template<typename T> inline genIntervalSet<T>::genIntervalSet(const genSequence<T> &seq) {
    for(size_t i=0;i<seq.num_intervals();++i)
      Union(seq[i]) ;
  }

  template<typename T> inline genIntervalSet<int_type> create_intervalSet(T start, T end) {
    if(start==end)
      return genIntervalSet<int_type>::EMPTY ;
    std::sort(start,end) ;
    int_type First = *start ;
    int_type Second = *start ;
    genIntervalSet<int_type> r ;
    for(T p=start;p!=end;++p) {
      if(*p > Second+1) {
        r += std::pair<int_type, int_type>(First,Second) ;
        First = *p ;
      }
      Second = *p ;
    }
    r += std::pair<int_type, int_type>(First,Second) ;
    return r ;
  }

  template<typename T> inline genSequence<int_type> create_sequence(T start, T end) {
    genSequence<int_type> s ;
    for(T p=start;p!=end;++p) {
      s += *p ;
    }
    return s ;
  }

  template<class T> inline genIntervalSet<gEntity> gcreate_intervalSet(T start, T end) {
    if(start==end)
      return genIntervalSet<gEntity>::EMPTY ;
    std::sort(start,end) ;
    gEntity First = *start ;
    gEntity Second = *start ;
    genIntervalSet<gEntity> r ;
    for(T p=start;p!=end;++p) {
      if(*p > Second+1) {
        r += std::pair<gEntity, gEntity>(First,Second) ;
        First = *p ;
      }
      Second = *p ;
    }
    r += std::pair<gEntity, gEntity>(First,Second) ;
    return r ;
  }

  template<class T> inline genSequence<gEntity> gcreate_sequence(T start, T end) {
    genSequence<gEntity> s ;
    for(T p=start;p!=end;++p) {
      s += *p ;
    }
    return s ;
  }
  
}
  
#define FORALL(var,indx)                                                \
  {  const Loci::genIntervalSet<int> &__p = var ;                       \
  for(size_t __ii=0;__ii<__p.num_intervals();++__ii)                       \
    for(Loci::int_type indx=__p[__ii].first;indx<=__p[__ii].second;++indx)
      
#define ENDFORALL }



#define GFORALL(var,indx)                                               \
  {  const Loci::genIntervalSet<gEntity> &__p = var ;                   \
  for(size_t __ii=0;__ii<__p.num_intervals();++__ii)                       \
    for(Loci::gEntity indx=__p[__ii].first;indx<=__p[__ii].second;++indx)
      
#define ENDGFORALL }

#endif 
