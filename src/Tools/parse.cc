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
    // -----------------------------------------------------------------------
    /// @brief kill_white_space() advances an input stream past leading white
    /// space and line comments beginning with //
    ///
    /// @param [s] input stream to be positioned at the next non-whitespace,
    /// non-comment character
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

    // -----------------------------------------------------------------------
    /// @brief is_name() checks whether the next token in the stream begins
    /// with a valid identifier lead character
    ///
    /// @param [s] input stream to inspect
    /// @return true if the next non-whitespace character is alphabetic or
    /// underscore
    bool is_name(istream &s) {
      kill_white_space(s) ;
      int ch = s.peek() ;
      return isalpha(ch) || ch == '_' ;
    }

    // -----------------------------------------------------------------------
    /// @brief get_name() extracts an identifier token from the input stream
    ///
    /// @param [s] input stream to consume from
    /// @return identifier made of alphanumeric and underscore characters, or
    /// an empty string if the next token is not a valid name
    string get_name(istream &s) {
      if(!is_name(s)) { return "" ; }
      string str ;
      while(!s.eof() && (s.peek() != EOF) &&
            (isalnum(s.peek()) || (s.peek() == '_')) ) {
        str += s.get() ;
      }

      return str ;
    }

    // -----------------------------------------------------------------------
    /// @brief is_int() checks whether the next token in the stream can begin
    /// an integer literal
    ///
    /// @param [s] input stream to inspect
    /// @return true if the next non-whitespace character is a digit or sign
    bool is_int(istream &s) {
      kill_white_space(s) ;
      return isdigit(s.peek()) || s.peek()=='-' || s.peek()=='+' ;
    }

    // -----------------------------------------------------------------------
    /// @brief get_int() extracts an integer value from the input stream
    ///
    /// @param [s] input stream to consume from
    /// @return parsed integer value, or zero if the next token does not begin
    /// an integer literal
    long get_int(istream &s) {
      if(!is_int(s)) { return 0 ; }
      long l = 0 ;
      s >> l ;
      return l ;
    }

    // -----------------------------------------------------------------------
    /// @brief is_real() checks whether the next token in the stream can begin
    /// a floating-point literal
    ///
    /// @param [s] input stream to inspect
    /// @return true if the next non-whitespace character can start a real
    /// number
    bool is_real(istream &s) {
      kill_white_space(s) ;
      const char ch = s.peek() ;
      return isdigit(ch) || ch=='-' || ch=='+' || ch =='.' ;
    }

    // -----------------------------------------------------------------------
    /// @brief get_real() extracts a floating-point value from the input stream
    ///
    /// @param [s] input stream to consume from
    /// @return parsed floating-point value, or zero if the next token does not
    /// begin a real literal
    double get_real(istream &s) {
      if(!is_real(s)) { return 0.0 ; }

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
        if(!leading_digit && !trailing_digit) {  // convert . to .0
          rval += '0' ;
        }
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

    // -----------------------------------------------------------------------
    /// @brief is_string() checks whether the next token in the stream is a
    /// quoted string literal
    ///
    /// @param [s] input stream to inspect
    /// @return true if the next non-whitespace character is a double quote
    bool is_string(istream &s) {
      kill_white_space(s) ;
      return s.peek() == '\"' ;
    }

    // -----------------------------------------------------------------------
    /// @brief get_string() extracts the contents of a quoted string literal
    ///
    /// @param [s] input stream to consume from
    /// @return characters between the opening and closing double quotes, or an
    /// empty string if the next token is not a string literal
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

    // -----------------------------------------------------------------------
    /// @brief is_token() checks whether the next token in the stream matches a
    /// specific literal without consuming it
    ///
    /// @param [s] input stream to inspect
    /// @param [token] literal token to compare against the input
    /// @return true if the next characters match token exactly
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

    // -----------------------------------------------------------------------
    /// @brief get_token() consumes a specific literal token from the input
    /// stream if it matches exactly
    ///
    /// @param [s] input stream to consume from
    /// @param [token] literal token expected at the current position
    /// @return true if token was matched and consumed, false otherwise
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
