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

#ifndef INTERVALSET_IMPL_H
#define INTERVALSET_IMPL_H 1

#include <Tools/debug.h>
#include <Tools/tools.h>
#include <Tools/intervalSet_def.h>

#include <algorithm>
#include <functional>

#include <istream>
#include <ostream>
#include <iostream>
#include <ctype.h>


namespace Loci {
  // provide generic int pair<int,int> output
  std::ostream & operator <<(std::ostream &s, const std::pair<int,int> &p) ;
  std::istream &operator >>(std::istream &s, std::pair<int,int> &p) ;


  template<typename T> typename genIntervalSet<T>::rep_holder* genIntervalSet<T>::rhp = 0 ;
  template<typename T> typename genSequence<T>::rep_holder* genSequence<T>::rhp = 0 ;
 
 
  /* The genIntervalSetRep class does all of the work for the genIntervalSet class.
   * It is referenced through a counted pointer (Handle) template.
   */
  template<typename T> std::ostream &outputInterval(std::ostream &s, const std::pair<T, T> &i) {
    s << '[' ;
    if(i.first==genIntervalSet<T>::UNIVERSE_MIN)
      s << '#' ;
    else 
      s << i.first ;
    s << ',' ;
    if(i.second==genIntervalSet<T>::UNIVERSE_MAX)
      s << '#' ;
    else
      s << i.second ;
    s << ']' ;
    return s;
  }


  template<typename T> std::istream &inputInterval(std::istream &s, std::pair<T, T> &i) {
    char ch ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!='[') {
      std::cerr << "Incorrect format when reading genInterval" << std::endl ;
      std::cerr << "expected a '[' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    if(s.peek() == '#') {
      ch = s.get() ;
      i.first = genIntervalSet<T>::UNIVERSE_MIN ;
    } else
      s >> i.first ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=',') {
      std::cerr << "Incorrect format when reading genInterval" << std::endl ;
      std::cerr << "expected a ',' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    if(s.peek() == '#') {
      ch = s.get() ;
      i.second = genIntervalSet<T>::UNIVERSE_MAX ;
    } else
      s >> i.second ;
    
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=']') {
      std::cerr << "Incorrect format when reading genInterval" << std::endl ;
      std::cerr << "expected a ']' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    return s ;
  }
    
  /*

    GenIntervalSet Class Code starts Here
   
  */

  template<typename T> std::ostream &genIntervalSet<T>::Print(std::ostream &s) const {
    s << "(" ;
    for(size_t i=0;i<num_intervals();++i) {
      outputInterval(s,(*Rep)[i]) ;
    } //      s << (*Rep)[i] ;
    s << ")" ;
    return s ;
  }

  template<typename T> std::istream &genIntervalSet<T>::Input(std::istream &s) {
    Rep.New() ;

    char ch ;
    // Startup ;

    // Kill white space
    do{
      ch = s.get() ;
    } while(isspace(ch)) ;

    // Check for correct startup character
    if(ch != '(') {
      std::cerr << "Incorrect format when reading genIntervalSet" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    // get genIntervals
    for(;;) {
      // Skip whitespace
      while(isspace(s.peek()))
        ch = s.get() ;
      // If next input is not an genInterval, then we are finished.
      if(s.peek()!='[')
        break ;
      // grab genInterval from input
      std::pair<T, T> ivl  ;
      inputInterval(s,ivl) ;
      //      s >> ivl ;
      // Add it up
      Union(ivl) ;
    }
    // Check that everything finished ok, then finish
    ch = s.get() ;
    if(ch!=')') {
      std::cerr << "Incorrect format when reading genIntervalSet" << std::endl ;
      s.putback(ch) ;
    }
      
    return s ;
  }

  // This partial ordering function defines an implicit equivalence
  // relation between genIntervals if they intersect one another
  namespace {
    template<typename T> inline  bool genInterval_porder_intersect(const std::pair<T, T>& i1, 
                                                                const std::pair<T, T>& i2) {
      return (i1.first < i2.first) && (i1.second < i2.first)  ;
    }
  }

  template<typename T> bool genIntervalSet<T>::inSet(T indx) const {
    return std::binary_search(Rep->begin(),Rep->end(),std::pair<T, T>(indx,indx),
                              genInterval_porder_intersect<T>) ;
  }
  
  // This partial ordering function defines an implicit equivalence
  // relation between genIntervals that are ordered lowvalue then high value
  // so that if they intersect or they touch they are in the same equivalence
  // relation.
  namespace {
    template<typename T> inline  bool genInterval_porder_union(const std::pair<T, T> &i1, 
                                                            const std::pair<T, T> &i2) {
      return i1.second+1 < i2.first ;
    }
  }
  
#ifdef ENTITY
  template<typename T> void Union_inplace (Handle<std::vector<std::pair<T, T> > > &Rep, const std::pair<T, T> &ivl) {
    Rep.MakeUnique() ;
    std::pair<T, T> i(min(ivl.first,ivl.second),
                      max(ivl.first,ivl.second)) ;
    if(num_intervals(Rep) == 0) {
      Rep->push_back(i) ;
      return ;
    }
    std::pair< typename std::vector<std::pair<T, T> >::iterator, typename std::vector<std::pair<T, T> >::iterator> range ;
    range = std::equal_range(Rep->begin(),Rep->end(),ivl,
                        genInterval_porder_union<T>) ;
    T range_size = range.second - range.first ;
    FATAL(range_size<0) ;
    switch(range_size) {
    case 0:
      Rep->insert(range.first,i) ;
      break ;
    case 1:
      (*range.first).first = min(i.first,(*range.first).first) ;
      (*range.first).second = max(i.second,(*range.first).second) ;
      break ;
    default:
      T lowval = (*range.first).first ;
      T hival  = (*(range.second-1)).second ;
      (*range.first).first = min(lowval,i.first) ;
      (*range.first).second = max(hival,i.second) ;
      Rep->erase(range.first+1,range.second) ;
      break ;
    }
  }


  template<typename T> Handle<std::vector<std::pair<T, T> > > Union(const Handle<std::vector<std::pair<T, T> > > &Rep1,
                                                                 const Handle<std::vector<std::pair<T, T> > > &Rep2) {
    typename std::vector<std::pair<T, T> >::const_iterator i1 = Rep1->begin() ;
    typename std::vector<std::pair<T, T> >::const_iterator i2 = Rep2->begin() ;
    const typename std::vector<std::pair<T, T> >::const_iterator e1 = Rep1->end() ;
    const typename std::vector<std::pair<T, T> >::const_iterator e2 = Rep2->end() ;
    const size_t size = Rep1->size() + Rep2->size() ;
    Handle<std::vector<std::pair<T, T> > > Rep ;
    Rep->reserve(size) ;
    while(i1!=e1 && i2!=e2) {
      if(genInterval_porder_union(*i1,*i2)) {
        Rep->push_back(*i1) ;
        ++i1 ;
      } else if(genInterval_porder_union(*i2,*i1)) {
        Rep->push_back(*i2) ;
        ++i2 ;
      } else {
        std::pair<T, T> ivl(min((*i1).first,(*i2).first),
                            max((*i1).second,(*i2).second)) ;
        ++i1 ;
        ++i2 ;
        bool flag = true ;
        while(flag) {
          flag = false ;
          if(i1!=e1 && !genInterval_porder_union(ivl,*i1)
             && !genInterval_porder_union(*i1,ivl)) {
            ivl.second = max(ivl.second,(*i1).second) ;
            ++i1 ;
            flag = true ;
          }
          if(i2!=e2 && !genInterval_porder_union(ivl,*i2)
             && !genInterval_porder_union(*i2,ivl)) {
            ivl.second = max(ivl.second,(*i2).second) ;
            ++i2 ;
            flag = true ;
          }
        } 
        Rep->push_back(ivl) ;
      }
    }
    while(i1!=e1) {
      Rep->push_back(*i1) ;
      ++i1 ;
    }
    while(i2!=e2) {
      Rep->push_back(*i2) ;
      ++i2 ;
    }
    return Rep ;
  }
#endif

#ifndef ENTITY
  template<typename T> void genIntervalSet<T>::Union(const std::pair<T, T> &ivl) {
    Rep.MakeUnique() ;
    std::pair<T, T> i(min(ivl.first,ivl.second),
                      max(ivl.first,ivl.second)) ;
    if(num_intervals() == 0) {
      Rep->push_back(i) ;
      return ;
    }
    std::pair<typename genIntervalSetRep::iterator, typename genIntervalSetRep::iterator> range ;
    range = std::equal_range(Rep->begin(),Rep->end(),i,
                        genInterval_porder_union<T>) ;
    T range_size = range.second - range.first ;
    fatal(range_size<0) ;
    switch(range_size) {
    case 0:
      Rep->insert(range.first,i) ;
      break ;
    case 1:
      (*range.first).first = min(i.first,(*range.first).first) ;
      (*range.first).second = max(i.second,(*range.first).second) ;
      break ;
    default:
      T lowval = (*range.first).first ;
      T hival  = (*(range.second-1)).second ;
      (*range.first).first = min(lowval,i.first) ;
      (*range.first).second = max(hival,i.second) ;
      Rep->erase(range.first+1,range.second) ;
      break ;
    }
  }


  template<typename T> void genIntervalSet<T>::Union(const genIntervalSet<T> &ptn) {
    genIntervalSet<T> tmp = genIntervalSet<T>::Union(*this,ptn) ;
    Rep = tmp.Rep ;
  }



  template<typename T> genIntervalSet<T> genIntervalSet<T>::Union(const genIntervalSet<T> &set1,
                                                               const genIntervalSet<T> &set2) {
    typename genIntervalSetRep::const_iterator i1 = set1.Rep->begin() ;
    typename genIntervalSetRep::const_iterator i2 = set2.Rep->begin() ;
    const typename genIntervalSetRep::const_iterator e1 = set1.Rep->end() ;
    const typename genIntervalSetRep::const_iterator e2 = set2.Rep->end() ;
    const size_t size = set1.Rep->size() + set2.Rep->size() ;
    Handle<genIntervalSetRep> Rep ;
    Rep->reserve(size) ;
    while(i1!=e1 && i2!=e2) {
      if(genInterval_porder_union(*i1,*i2)) {
        Rep->push_back(*i1) ;
        ++i1 ;
      } else if(genInterval_porder_union(*i2,*i1)) {
        Rep->push_back(*i2) ;
        ++i2 ;
      } else {
        std::pair<T, T> ivl(min((*i1).first,(*i2).first),
                            max((*i1).second,(*i2).second)) ;
        ++i1 ;
        ++i2 ;
        bool flag = true ;
        while(flag) {
          flag = false ;
          if(i1!=e1 && !genInterval_porder_union(ivl,*i1)
             && !genInterval_porder_union(*i1,ivl)) {
            ivl.second = max(ivl.second,(*i1).second) ;
            ++i1 ;
            flag = true ;
          }
          if(i2!=e2 && !genInterval_porder_union(ivl,*i2)
             && !genInterval_porder_union(*i2,ivl)) {
            ivl.second = max(ivl.second,(*i2).second) ;
            ++i2 ;
            flag = true ;
          }
        } 
        Rep->push_back(ivl) ;
      }
    }
    while(i1!=e1) {
      Rep->push_back(*i1) ;
      ++i1 ;
    }
    while(i2!=e2) {
      Rep->push_back(*i2) ;
      ++i2 ;
    }
    genIntervalSet<T> result ;
    if(Rep->size() != 0)
      result.Rep = Rep ;
    return result ;
  }

#endif

  template<typename T> void genIntervalSet<T>::Complement() {
    Rep.MakeUnique() ;
    if(Rep->size() == 0) {
      Rep->push_back(std::pair<T, T>(genIntervalSet<T>::UNIVERSE_MIN,genIntervalSet<T>::UNIVERSE_MAX)) ;
      return ;
    }

    if((*Rep)[0].first == genIntervalSet<T>::UNIVERSE_MIN) {
      typename genIntervalSetRep::iterator stop = Rep->end()-1 ;
      for(typename genIntervalSetRep::iterator i=Rep->begin();i!=stop;) {
        i->first = i->second+1 ;
        typename genIntervalSetRep::iterator next = i+1 ;
        i->second = next->first-1 ;
        i = next ;
      }
      if((*stop).second == genIntervalSet<T>::UNIVERSE_MAX)
        Rep->erase(stop) ;
      else {
        stop->first = stop->second+1 ;
        stop->second = genIntervalSet<T>::UNIVERSE_MAX ;
      }
      return ;
    }

    typename genIntervalSetRep::iterator i = Rep->begin() ;
    T sec = i->second ;
    T fir = i->first ;
    (*i).first = UNIVERSE_MIN ;
    (*i).second = fir-1 ;
    for(++i;i!=Rep->end();++i) {
      fir = i->first ;
      i->first = sec+1 ;
      sec = i->second ;
      i->second = fir-1 ;
    }
    if(sec != genIntervalSet<T>::UNIVERSE_MAX)
      Rep->push_back(genInterval(sec+1,genIntervalSet<T>::UNIVERSE_MAX)) ;
    return ;
  }

  template<typename T> genIntervalSet<T> genIntervalSet<T>::Complement(const genIntervalSet<T> &set) {
    const size_t size = set.Rep->size() ;
    Handle<genIntervalSetRep> Rep ;
    Rep->reserve(size+1) ;
    if(size == 0)
      Rep->push_back(std::pair<T, T>(genIntervalSet<T>::UNIVERSE_MIN,genIntervalSet<T>::UNIVERSE_MAX)) ;
    else {
      if(set.Rep->front().first != genIntervalSet<T>::UNIVERSE_MIN)
        Rep->push_back(std::pair<T, T>(genIntervalSet<T>::UNIVERSE_MIN,set.Rep->front().first-1)) ;

      T sec = set.Rep->front().second ;
      for(size_t i=1;i<size;++i) {
        const std::pair<T, T> &ci = (*(set.Rep))[i] ;
        Rep->push_back(std::pair<T, T>(sec+1,ci.first-1)) ;
        sec = ci.second ;
      }
      if(sec!=genIntervalSet<T>::UNIVERSE_MAX) 
        Rep->push_back(std::pair<T, T>(sec+1,genIntervalSet<T>::UNIVERSE_MAX)) ;
    }
        
    genIntervalSet<T> result ;
    if(Rep->size() != 0)
      result.Rep = Rep ;
    return result ;
  }

  template<typename T>  genIntervalSet<T> &genIntervalSet<T>::operator>>=(T rotval) {
    Rep.MakeUnique() ;
    for(typename genIntervalSetRep::iterator i = Rep->begin();i!=Rep->end();++i) {
      if(i->first != genIntervalSet<T>::UNIVERSE_MIN)
        i->first += rotval ;
      if(i->second != genIntervalSet<T>::UNIVERSE_MAX)
        i->second += rotval ;
    }
    return *this ;
  }

  template<typename T> genIntervalSet<T> &genIntervalSet<T>::operator<<=(T rotval) {
    Rep.MakeUnique() ;
    for(typename genIntervalSetRep::iterator i = Rep->begin();i!=Rep->end();++i) {
      if(i->first != genIntervalSet<T>::UNIVERSE_MIN)
        i->first -= rotval ;
      if(i->second != genIntervalSet<T>::UNIVERSE_MAX)
        i->second -= rotval ;
    }
    return *this ;
  }


  //  Finding an intersection to an std::pair<T, T> is a simple matter of extracting
  //  the equivalence set for under the partial ordering for overlaping
  //  genIntervals and then truncating the range to the std::pair<T, T> bounds.
  template<typename T> void genIntervalSet<T>::Intersection(const std::pair<T,T> &ivl) {
    std::pair<typename genIntervalSetRep::iterator, typename genIntervalSetRep::iterator> range ;
    range = std::equal_range(Rep->begin(),Rep->end(),ivl,
                        genInterval_porder_intersect<T>) ;
    T range_size = range.second-range.first ;
    fatal(range_size<0) ;
    if(range_size == 0) {
      Rep = EMPTY.Rep ;
      return ;
    }

    // Make new pair vector
    std::vector<std::pair<T, T> > pv(range_size) ;
    size_t cnt = 0 ;
    while(range.first!=range.second) {
      pv[cnt++] = *range.first ;
      ++range.first ;
    }
    const T l = range_size-1 ;
    pv[0].first  = max(ivl.first, pv[0].first) ;
    pv[l].second = min(ivl.second,pv[l].second) ;

    // Make new Rep, and make this set refer to it.
    Handle<std::vector<std::pair<T, T> > > Rpnew ;
    Rpnew->swap(pv) ;
    Rep = Rpnew ;
    return ;
  }


  template<typename T> void genIntervalSet<T>::Intersection(const genIntervalSet<T> &ptn) {
    genIntervalSet<T> tmp = genIntervalSet<T>::Intersection(*this,ptn) ;
    Rep = tmp.Rep ; 
  }
  
 
  template<typename T> genIntervalSet<T> genIntervalSet<T>::Intersection(const genIntervalSet<T> &set1,
                                                                      const genIntervalSet<T> &set2)
  {
    if(set1 == genIntervalSet<T>::EMPTY || set2 == genIntervalSet<T>::EMPTY)
      return genIntervalSet<T>::EMPTY ;
    
    typename genIntervalSetRep::const_iterator i1 = set1.Rep->begin() ;
    typename genIntervalSetRep::const_iterator i2 = set2.Rep->begin() ;
    const typename genIntervalSetRep::const_iterator e1 = set1.Rep->end() ;
    const typename genIntervalSetRep::const_iterator e2 = set2.Rep->end() ;
    const size_t size = set1.Rep->size() + set2.Rep->size() ;
    Handle<genIntervalSetRep> Rep ;
    Rep->reserve(size) ;
    while(i1!=e1 && i2!=e2) {
      if(genInterval_porder_intersect(*i1,*i2))
        ++i1 ;
      else if(genInterval_porder_intersect(*i2,*i1))
        ++i2 ;
      else {
        std::pair<T, T> ivl(max((*i1).first,(*i2).first),
                            min((*i1).second,(*i2).second)) ;
        Rep->push_back(ivl) ;
        if((*i1).second > (*i2).second)
          ++i2 ;
        else
          ++i1 ;
      }
    }
    genIntervalSet<T> result ;
    if(Rep->size() != 0)
      result.Rep = Rep ;
    return result ;
  }
  
  template<typename T>  std::ostream &genSequence<T>::Print(std::ostream &s) const {
    s << "(" ;
    for(size_t i=0;i<num_intervals();++i)
      s << (*Rep)[i] ;
    s << ")" ;
    return s ;
  }
  
  template<typename T> std::istream &genSequence<T>::Input(std::istream &s) {
    Rep.New() ;
    
    char ch ;
    // Startup ;
    
    // Kill white space
    do{
      ch = s.get() ;
    } while(isspace(ch)) ;

    // Check for correct startup character
    if(ch != '(') {
      std::cerr << "Incorrect format when reading genIntervalSet" << std::endl ;
      s.putback(ch) ;
      return s ;
    }

    // get genIntervals
    for(;;) {
      // Skip whitespace
      while(isspace(s.peek()))
        ch = s.get() ;
      // If next input is not an genInterval, then were finished
      if(s.peek()!='[')
        break ;
      // grab std::pair<T, T> from input
      std::pair<T, T> ivl  ;
      s >> ivl ;
      // Add it up
      Append(ivl) ;
    }
    // Check that everything finished ok, then finish
    ch = s.get() ;
    if(ch!=')') {
      std::cerr << "Incorrect format when reading genIntervalSet" << std::endl ;
      s.putback(ch) ;
    }
      
    return s ;
  }
   
  template<typename T> void genSequence<T>::Append(const std::pair<T, T> &ivl) {
    Rep.MakeUnique() ;
    if(Rep->size() == 0) {
      Rep->push_back(ivl) ;
      return ;
    }
            
    std::pair<T, T> &back = Rep->back() ;
    T dir = back.second-back.first ;
    T dirivl = ivl.second-ivl.first ;
    if((ivl.first == (back.second+1)) && dirivl >= 0 && dir >= 0) {
      (Rep->back()).second = ivl.second ;
      return ;
    }
    if((ivl.first == (back.second-1)) && dirivl <= 0 && dir <= 0) {
      (Rep->back()).second = ivl.second ;
      return ;
    }            
    Rep->push_back(ivl) ;
  }

  template<typename T> void genSequence<T>::Append(const genSequence<T> &seq) {
    if(seq.Rep->size() == 0)
      return ;
    Rep.MakeUnique() ;
    typename genSequenceRep::const_iterator i = seq.Rep->begin() ;
    Append(*i) ;
    ++i ;
    Rep->insert(Rep->end(),i,seq.Rep->end()) ;
  }

  template<typename T> void genSequence<T>::Append(const genIntervalSet<T> &ptn) {
    genSequence<T> seq(ptn) ;
    Append(seq) ;
  }

  template<typename T> genSequence<T> & genSequence<T>::Reverse() {
    Rep.MakeUnique() ;
    std::reverse(Rep->begin(),Rep->end()) ;
    for(typename genSequenceRep::iterator i=Rep->begin();i!=Rep->end();++i)
      std::swap(i->first,i->second) ;
    return *this ;
  }
  
  template<typename T> genSequence<T> &genSequence<T>::operator>>=(T rotval) {
    Rep.MakeUnique() ;
    for(typename genSequenceRep::iterator i = Rep->begin();i!=Rep->end();++i) {
      if(i->first != genIntervalSet<T>::UNIVERSE_MIN)
        i->first += rotval ;
      if(i->second != genIntervalSet<T>::UNIVERSE_MAX)
        i->second += rotval ;
    }
    return *this ;
  }

  template<typename T> genSequence<T> &genSequence<T>::operator<<=(T rotval) {
    Rep.MakeUnique() ;
    for(typename genSequenceRep::iterator i = Rep->begin();i!=Rep->end();++i) {
      if(i->first != genIntervalSet<T>::UNIVERSE_MIN)
        i->first -= rotval ;
      if(i->second != genIntervalSet<T>::UNIVERSE_MAX)
        i->second -= rotval ;
    }
    return *this ;
  }
}
    
#endif
