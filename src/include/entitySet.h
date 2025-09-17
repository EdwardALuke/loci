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
#ifndef ENTITYSET_H
#define ENTITYSET_H 1

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/intervalSet.h>


namespace Loci {

  typedef int int_type ;

#ifdef ENTITY

  class Entity {
  private:
    int id;

  public:
    Entity() {}
    Entity(int i) {id = i;}
    Entity(const Entity &e) {id = e.id;} /* copy constructor */
    int getid() const {return id;}
    std::ostream & Print (std::ostream &s) const { 
      s << id ;
      return s;}
    std::istream & Read (std::istream &s) {
      s >> id;
      return s;
    }
    bool operator==(const Entity &e) const {
      return (id == e.getid());
    }
    operator int () const {return getid();}
  };
  
  inline std::ostream &operator<<(std::ostream &s, const Entity &e) {
    return e.Print(s);
  }
  
  inline std::istream &operator>>(std::istream &s, Entity &e) {
    return e.Read(s);
  }
  
  
  inline int getEntityIdentity(const Entity e){return e.getid();}  

  typedef intervalSet entitySet ;


  template<class T> inline entitySet create_entitySet(T start,T end) {
      return create_intervalSet(start,end) ;
  } 

//      class entitySet {
//    public:
//      class entitySetIterator ;
//      typedef pair_vector entitySetRep ;
//    private:
//      friend class sequence ;
//      Handle<entitySetRep> Rep ;
    
//      struct rep_holder {
//        Handle<entitySetRep> Rep ;
//      } ;
    
//      static rep_holder *rhp ;
//      static Handle<entitySetRep> & getEmptyRep() {
//        if(rhp==0)
//          rhp = new rep_holder ;
//        return rhp->Rep ;
//      } ;
//    public:
//      class entitySetIterator {
//        entitySetRep::const_iterator current_interval,end_interval ;
//        int_type current_value ;
//        void increment() {
//          if(++current_value>current_interval->second) {
//            ++current_interval ;
//            if(current_interval!=end_interval)
//              current_value = current_interval->first ;
//            else
//              current_value = UNIVERSE_MAX ;
//          }
//        }                
//      public:
//        entitySetIterator() {}
//          //      {current_interval=0;end_interval=0;current_value=0;}
//        entitySetIterator(entitySetRep::const_iterator r,
//                            entitySetRep::const_iterator end,int_type e) :
//          current_interval(r),end_interval(end),current_value(e) {}
//        int operator*() const { return current_value; }
//        entitySetIterator &operator++() {
//          increment() ;
//          return *this ;
//        }
//        entitySetIterator operator++(int ) {
//          entitySetIterator tmp = *this ;
//          increment() ;
//          return tmp ;
//        }
//        bool operator==(const entitySetIterator &i) const {
//          return current_value == i.current_value ; }
//        bool operator!=(const entitySetIterator &i) const {
//          return current_value != i.current_value ; }
//      } ;


//      entitySet() : Rep(getEmptyRep()) {} ;
//      entitySet(const interval &ivl) { 
//  	interval i(min(ivl.first,ivl.second),
//  		   max(ivl.first,ivl.second));
//  	Rep->push_back(i) ; }
//      entitySet(const entitySet &ptn): Rep(ptn.Rep) {}
//      explicit entitySet(const sequence &seq) ;
//      explicit entitySet(const Handle<pair_vector> &RepIn): Rep(RepIn) {}
//      ~entitySet() {}
    
//      typedef entitySetIterator const_iterator ;

//      const_iterator begin() const {
//        if(Rep->size() != 0)
//          return entitySetIterator(Rep->begin(),Rep->end(),(*Rep)[0].first) ;
//        else
//          return end() ;
//      }
//      const_iterator end() const {
//        return entitySetIterator(Rep->end(),Rep->end(),UNIVERSE_MAX) ;
//      }
    

//      bool inSet(int_type indx) const ;

//      entitySet & operator=(const entitySet &ptn)  
//      { Rep = ptn.Rep ; return *this;}
//      entitySet & operator=(const interval &ivl) 
//      { Rep.New() ; Rep->push_back(ivl) ; return *this ; }

//      bool Equal(const entitySet &ptn) const
//      { return *(Rep) == *(ptn.Rep) ; }

//      bool less_than(const entitySet &ptn) const
//      { return (*Rep) < *(ptn.Rep) ;}
//      bool greater_than(const entitySet &ptn) const
//      {return (*Rep) > *(ptn.Rep) ;}
    
//      entitySet & operator>>=(int rotval) ;
//      entitySet & operator<<=(int rotval) ;

//      entitySet operator<<(int rotval) const ;
//      entitySet operator>>(int rotval) const ;

//      int size() const {
//        int size = 0 ;
//        std::vector<interval>::const_iterator i ;
//        for(i = Rep->begin();i!= Rep->end();++i)
//          size += i->second-i->first + 1 ;
//        return size ;  }

//      int num_intervals() const { return Rep->size() ; }
//      const interval & operator[](int_type indx) const 
//      { fatal(indx<0); fatal(indx>=num_intervals()) ;
//        return (*Rep)[indx]; }

      
//      void Union(const interval &ivl){Loci::Union(Rep,ivl);}
//        void Union(const entitySet &ptn){
//            int psz = ptn.num_intervals() ;
//            if(psz == 0)// If ptn == EMPTY, do nothing
//              return ;
//            if(psz == 1) // only one interval, use interval update in place
//              Loci::Union(Rep,ptn[0]);
//            else
//              Rep = Loci::Union(Rep,ptn.Rep) ; // General case
//        }
      
//      static entitySet Union(const entitySet &set1,
//                                          const entitySet &set2){
//        return entitySet(Loci::Union(set1.Rep,set2.Rep)) ;
//      }

//      void Intersection(const interval &ivl) ;
//      void Intersection(const entitySet &ptn) ;
//      static entitySet Intersection(const entitySet &set1,
//                                      const entitySet &set2) ;
//      void Complement() ;
//      static entitySet Complement(const entitySet &set) ;
        
//      std::ostream &Print(std::ostream &s) const ;
//      void Print() const { Print(std::cout) ; }

//      std::istream &Input(std::istream &s) ;

//      int Min() const {
//        int sz = Rep->size();
//        return sz==0?UNIVERSE_MAX:(*Rep)[0].first ;}
//      int Max() const {
//        int sz = Rep->size() ;
//        return sz==0?UNIVERSE_MIN:(*Rep)[sz-1].second ;}
//    } ;


#else
  typedef Loci::int_type Entity ;
  inline int getEntityIdentity(const int i){return i;}
  typedef intervalSet entitySet ;
  template<class T> inline entitySet create_entitySet(T start,T end) {
      return create_intervalSet(start,end) ;
  }
#endif

  template<class Op> inline void do_loop(const entitySet &iset, Op f) {
    for(int i=0;i<iset.num_intervals();++i) {
      const Loci::int_type stop = iset[i].second ;
      for(Loci::int_type indx=iset[i].first;indx<=stop;++indx)
        f(indx) ;
    }
  }

  template<class Op> inline void do_loop(const sequence &seq, Op f) {
    const int ni = seq.num_intervals() ; 
    for(int i=0;i<ni;++i) {
      const Loci::int_type start = seq[i].first ;
      const Loci::int_type stop = seq[i].second  ;
      if(seq[i].first<=seq[i].second)	
	for(Loci::int_type indx=start;indx!=stop+1;++indx) 
	  f(indx) ;
      else
	for(Loci::int_type indx=start;indx!=stop-1;--indx) 
	  f(indx) ;
    }
  }
  
  template<class T> inline void do_loop(const sequence &seq, T *cp,
                                        void(T::*pmf)(Entity) ) {
    const int ni = seq.num_intervals() ;
    for(int i=0;i<ni;++i) {
      const Loci::int_type i1 = seq[i].first ;
      const Loci::int_type i2 = seq[i].second  ;
      const Loci::int_type start = min(i1,i2) ;
      const Loci::int_type stop = max(i1,i2) ;
      //#ifdef HAVE_IVDEP
      //#pragma ivdep
      //#endif
      for(Loci::int_type indx=start;indx!=stop+1;++indx) 
        (cp->*pmf)(indx) ;
    }
  }

  template<class T> inline void do_loop(const sequence &seq, T *cp) {
    const int ni = seq.num_intervals() ;
    for(int i=0;i<ni;++i) {
      const Loci::int_type i1 = seq[i].first ;
      const Loci::int_type i2 = seq[i].second  ;
      const Loci::int_type start = min(i1,i2) ;
      const Loci::int_type stop = max(i1,i2) ;
      //#ifdef HAVE_IVDEP
      //#pragma ivdep
      //#endif
      for(Loci::int_type indx=start;indx!=stop+1;++indx) 
        cp->calculate(indx) ;
    }
  }

  template<class T> inline void do_segments(const sequence &seq, T *cp,
                                        void(T::*pmf)(Entity) ) {
    const int ni = seq.num_intervals() ;
    for(int i=0;i<ni;++i) {
      const Loci::int_type i1 = seq[i].first ;
      const Loci::int_type i2 = seq[i].second  ;
      const Loci::int_type start = min(i1,i2) ;
      const Loci::int_type stop = max(i1,i2) ;
      (cp->*pmf)(start,stop-start+1) ;
    }
  }

  template<class T> inline void do_segments(const sequence &seq, T *cp) {
    const int ni = seq.num_intervals() ;
    for(int i=0;i<ni;++i) {
      const Loci::int_type i1 = seq[i].first ;
      const Loci::int_type i2 = seq[i].second  ;
      const Loci::int_type start = min(i1,i2) ;
      const Loci::int_type stop = max(i1,i2) ;
      cp->segment(start,stop-start+1) ;
    }
  }

  template<class T> inline void do_loop_seq(const sequence &seq, T *cp,
                                        void(T::*pmf)(Entity) ) {
    const int ni = seq.num_intervals() ;
    for(int i=0;i<ni;++i) {
      const Loci::int_type start = seq[i].first ;
      const Loci::int_type stop = seq[i].second  ;
      if(seq[i].first<=seq[i].second)	
	for(Loci::int_type indx=start;indx!=stop+1;++indx) 
	  (cp->*pmf)(indx) ;
      else
	for(Loci::int_type indx=start;indx!=stop-1;--indx) 
	  (cp->*pmf)(indx) ;
    }
  }

  template<class T> inline void do_loop_seq(const sequence &seq, T *cp) {
    const int ni = seq.num_intervals() ;
    for(int i=0;i<ni;++i) {
      const Loci::int_type start = seq[i].first ;
      const Loci::int_type stop = seq[i].second  ;
      if(seq[i].first<=seq[i].second)	
	for(Loci::int_type indx=start;indx!=stop+1;++indx) 
	  cp->calculate(indx) ;
      else
	for(Loci::int_type indx=start;indx!=stop-1;--indx) 
	  cp->calculate(indx) ;
    }
  }
}

#endif
