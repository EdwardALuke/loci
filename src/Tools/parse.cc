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
    namespace {
      // ctype predicates are only defined for EOF and unsigned-char values.
      // This wrapper keeps the stream-facing call sites safe and compact.
      bool is_digit(int ch) {
        return ch != EOF && isdigit(static_cast<unsigned char>(ch)) ;
      }

      // parse identifiers with the same leading-character rule used by
      // is_name(): alphabetic characters and underscores are allowed.
      bool is_name_start(int ch) {
        return ch != EOF &&
          (isalpha(static_cast<unsigned char>(ch)) || ch == '_') ;
      }

      // is_int() accepts signed literals only when the sign is followed by a
      // digit.  Peek ahead without consuming the stream.
      bool signed_int_starts_here(istream &s) {
        const int sign = s.get() ;
        const bool ok = is_digit(s.peek()) ;
        s.putback(static_cast<char>(sign)) ;
        return ok ;
      }

      // is_real() accepts signed literals when the sign is followed by either
      // a digit or a decimal point.  Leave the stream unchanged after probing.
      bool signed_real_starts_here(istream &s) {
        const int sign = s.get() ;
        const int next = s.peek() ;
        s.putback(static_cast<char>(sign)) ;
        return is_digit(next) || next == '.' ;
      }
    }

    void kill_white_space(istream &s) {

      bool flushed_comment ;
      do {
        flushed_comment = false ;
        while(!s.eof() &&
              isspace(static_cast<unsigned char>(s.peek())))
          s.get() ;
        if(s.peek() == '/') { // check for comment
          s.get() ;
          if(s.peek() == '/') {
            while(!s.eof()) {
              int ch = s.get() ;
              if(ch=='\n' || ch == '\r') { break ; }
            }
            flushed_comment = true ;
          } else {
            s.putback('/') ;
          }
        }
      } while(!s.eof() && flushed_comment) ;

      return ;
    }

    bool is_name(istream &s) {
      kill_white_space(s) ;
      int ch = s.peek() ;
      return is_name_start(ch) ;
    }

    string get_name(istream &s) {
      if(!is_name(s)) { return "" ; }
      string str ;
      while(!s.eof() && (s.peek() != EOF) &&
            (isalnum(static_cast<unsigned char>(s.peek())) ||
             (s.peek() == '_')) ) {
        str += s.get() ;
      }

      return str ;
    }

    bool is_int(istream &s) {
      kill_white_space(s) ;
      const int ch = s.peek() ;
      if(is_digit(ch)) {
        return true ;
      }
      return (ch == '-' || ch == '+') && signed_int_starts_here(s) ;
    }

    long get_int(istream &s) {
      if(!is_int(s)) { return 0 ; }
      long l = 0 ;
      s >> l ;
      return l ;
    }

    bool is_real(istream &s) {
      kill_white_space(s) ;
      const int ch = s.peek() ;
      if(is_digit(ch) || ch == '.') {
        return true ;
      }
      return (ch == '-' || ch == '+') && signed_real_starts_here(s) ;
    }

    double get_real(istream &s) {
      if(!is_real(s)) { return 0.0 ; }

      // First grab real into string rval
      string rval ;
      char ch = s.get() ;
      rval += ch ; // since we already passed is_real, the first character
      // is in rval

      bool leading_digit = is_digit(ch) ;

      // any leading digits will go in rval
      while(is_digit(s.peek())) {
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
        while(is_digit(s.peek())) {
          trailing_digit = true ;
          ch = s.get() ;
          rval += ch ;
        }
        if(!leading_digit && !trailing_digit) {  // convert . to .0
          rval += '0' ;
        }
      }
      // If there is an exponent, check to make sure it is followed by a digit
      // after an optional sign.  Otherwise leave the exponent marker untouched
      // so the caller can see the malformed suffix.
      if(s.peek() == 'e' || s.peek() == 'E') {
        const char exp = s.get() ;
        if(is_digit(s.peek())) {
          rval += 'e' ;
          ch = s.get() ;
          rval += ch ;
          while(is_digit(s.peek())) {
            ch = s.get() ;
            rval += ch ;
          }
        } else if(s.peek() == '-' || s.peek() == '+') {
          const char sign = s.get() ;
          if(is_digit(s.peek())) {
            rval += 'e' ;
            rval += sign ;
            while(is_digit(s.peek())) {
              ch = s.get() ;
              rval += ch ;
            }
          } else {
            s.putback(sign) ;
            s.putback(exp) ;
          }
        } else {
          s.putback(exp) ;
        }
      }

      return atof(rval.c_str()) ;
    }

    bool is_string(istream &s) {
      kill_white_space(s) ;
      return s.peek() == '\"' ;
    }

    string get_string(istream &s) {
      if(!is_string(s)) { return "" ; }
      string str ;
#ifdef DEBUG
      if(s.eof()) {
        cerr << "s.eof() true in parse::get_string" << endl  ;
      }
#endif
      s.get() ;
      int ch = s.get() ;
      while(ch != '\"' &&!s.eof()) {
        str += ch ;
        ch = s.get() ;
      }
#ifdef DEBUG
      if(ch!='\"') {
        cerr << "no closing \" in parse::get_string" << endl ;
      }
#endif
      return str ;
    }

    bool  is_token(istream &s, const string &token) {
      kill_white_space(s) ;
      const int sz = token.size() ;
      for(int i=0;i<sz;++i) {
        if(s.peek() != token[i]) {
          for(--i;i>=0;--i) {
            s.putback(token[i]) ;
          }
          return false ;
        }
        s.get() ;
      }
      for(int i=token.size()-1;i>=0;--i) {
        s.putback(token[i]) ;
      }
      return true ;
    }

    bool get_token(istream &s, const string &token) {
      kill_white_space(s) ;
      const int sz = token.size() ;
      for(int i=0;i<sz;++i) {
        if(s.peek() != token[i]) {
          for(--i;i>=0;--i) {
            s.putback(token[i]) ;
          }
          return false ;
        }
        s.get() ;
      }
      return true ;
    }
  }
}
