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
#include <Tools/parse.h>
#include <Tools/debug.h>
#include <locale>
#ifdef NO_CSTDLIB
#include <stdio.h>
#else
#include <cstdio>
#endif

namespace Loci {
  namespace parse {
    
    using namespace std ;
    void kill_white_space(istream &s) {
      
      bool flushed_comment ;
      do {
	flushed_comment = false ;
	while(!s.eof() && isspace(s.peek()))
	  s.get() ;
	if(s.peek() == '/') { // check for comment
	  s.get() ;
	  if(s.peek() == '/') {
	    while(!s.eof()) {
	      int ch = s.get() ;
	      if(ch=='\n' || ch == '\r')
		break ;
	    }
	    flushed_comment = true ;
	  } else
	    s.putback('/') ;
	}
      } while(!s.eof() && flushed_comment) ;
      
      return ;
    }
    
    bool is_name(istream &s) {
        kill_white_space(s) ;
        int ch = s.peek() ;
        return isalpha(ch) || ch == '_' ;
    }
    
    string get_name(istream &s) {
        if(!is_name(s))
          return "" ;
        string str ;
        while(!s.eof() && (s.peek() != EOF) &&
              (isalnum(s.peek()) || (s.peek() == '_')) )
          str += s.get() ;

        return str ;
    }

    bool is_int(istream &s) {
        kill_white_space(s) ;
        return isdigit(s.peek()) || s.peek()=='-' || s.peek()=='+' ;
    }
    
    long get_int(istream &s) {
        if(!is_int(s))
          return 0 ;
        long l = 0 ;
        s >> l ;
        return l ;
    }

    bool is_real(istream &s) {
        kill_white_space(s) ;
        const char ch = s.peek() ;
        return isdigit(ch) || ch=='-' || ch=='+' || ch =='.' ;
    }
    
    double get_real(istream &s) {
      if(!is_real(s)) {
	return 0.0 ;
      }

      // First grab real into string rval
      string rval ;
      char ch = s.get() ;
      rval += ch ; // since aready passing is_real we know first character
      // is in rval

      bool leading_digit = isdigit(ch) ;

      // any leading digits will go in rval
      while(isdigit(s.peek())) {
	ch = s.get() ;
	leading_digit = true ;
	rval += ch ;
      }
      // If there is a point, then the point and any following digits will
      // go into rval
      if(s.peek() == '.') {
	ch = s.get() ;
	rval += ch ;
	bool trailing_digit = false ;
	while(isdigit(s.peek())) {
	  trailing_digit = true ;
	  ch = s.get() ;
	  rval += ch ; 
	}
	if(!leading_digit && !trailing_digit)  // convert . to .0 
	  rval += '0' ; 
      }
      // If there is an exponent, check to make sure it is followed by a digit
      // if it is then grab the exponent, else put back the character with 
      // unget
      if(s.peek() == 'e' || s.peek() == 'E') {
	ch = s.get() ;
	ch = s.peek() ;
	if(isdigit(ch) || ch=='-' || ch=='+') { // valid exponent
	  rval += 'e' ;
	  ch = s.get() ;
	  rval += ch ;
	  while(isdigit(s.peek())) {
	    ch = s.get() ;
	    rval += ch ; 
	  }
	} else { // invalid exponent, ignore 'e' or 'E'
	  s.unget() ;
	}
      }

      return atof(rval.c_str()) ;
    }

    bool is_string(istream &s) {
        kill_white_space(s) ;
        return s.peek() == '\"' ;
    }
    
    string get_string(istream &s) {
        if(!is_string(s))
          return "" ;
        string str ;
#ifdef DEBUG        
        if(s.eof())
          cerr << "s.eof() true in parse::get_string" << endl  ;
#endif
        s.get() ;
        int ch = s.get() ;
        while(ch != '\"' &&!s.eof()) {
            str += ch ;
            ch = s.get() ;
        }
#ifdef DEBUG
        if(ch!='\"')
          cerr << "no closing \" in parse::get_string" << endl ;
#endif
        return str ;
    }

    bool  is_token(istream &s, const string &token) {
      kill_white_space(s) ;
      const int sz = token.size() ;
      for(int i=0;i<sz;++i) {
        if(s.peek() != token[i]) {
          for(--i;i>=0;--i)
            s.putback(token[i]) ;
          return false ;
        }
        s.get() ;
      }
      for(int i=token.size()-1;i>=0;--i) 
        s.putback(token[i]) ;
      return true ;
    }
    
    bool get_token(istream &s, const string &token) {
      kill_white_space(s) ;
      const int sz = token.size() ;
      for(int i=0;i<sz;++i) {
	if(s.peek() != token[i]) {
	  for(--i;i>=0;--i)
	    s.putback(token[i]) ;
	  return false ;
	}
	s.get() ;
      }
      return true ;
    }
  }    
}
