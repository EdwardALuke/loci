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
//#define VERBOSE
#include "lpp.h"
#include "parseAST.h"
#include <ctype.h>
#include <set>
#include <iostream>
#include <sstream>
//#include <sys/timeb.h>
#include <time.h>
#include <vector>

using std::istringstream ;
using std::ostringstream ;

using std::pair ;
using std::list ;
using std::string ;
using std::set ;
using std::map ;

using std::istream ;
using std::ifstream ;
using std::ofstream ;
using std::ostream ;
using std::ios ;
using std::endl ;
using std::cerr ;
using std::cout ;
using std::vector ;
using namespace Loci ;

// This is the Lexical Analysis part of the AST infrastructure

/// Stack of tokens used when parsing to implement return tokens to the input
/// token stream when a parsing rule fails.  This allows alternative rules
/// to be evaluated.
vector<CPTR<AST_Token> >  tokenStack ;

void pushToken(CPTR<AST_Token> &pt) {
#ifdef VERBOSE
  cerr << "pushing token " << pt->text << endl ;
#endif
  tokenStack.push_back(pt) ;
}

CPTR<AST_Token> getNumberToken(std::istream &is, int &linecount) {
  // Process elements that start with 0-9, '.'
  CPTR<AST_Token> AST_data = new AST_Token() ;
  AST_data->lineno = linecount ;
  string numberdata ;
  bool hasPoint = false ;
  if(is.peek() == '.')  { // This might be a floating point number
    numberdata += is.get() ;
    if(is.peek() <'0' || is.peek() > '9') {
      // Not a floating point number
      AST_data->nodeType = AST_type::TK_DOT ;
      AST_data->text = numberdata ;
      return AST_data ;
    }
    hasPoint = true ;
  }
  if((is.peek() >= '0' && is.peek() <='9')) {
    AST_data->nodeType = AST_type::TK_NUMBER ;
    if(is.peek() == '0' && !hasPoint) {
      // This path is either binary, octal, or hexidecimal
      numberdata += is.get() ;
      if(is.peek() == 'x' || is.peek() == 'X') { //hex number
	numberdata += is.get() ;
	while((is.peek() >='0' && is.peek() <= '9') ||
	      (is.peek() >='a' && is.peek() <= 'f') ||
	      (is.peek() >='A' && is.peek() <= 'F')) {
	  numberdata += is.get() ;
	}
	if(numberdata.size() == 2)
	  AST_data->nodeType = AST_type::TK_ERROR ;
	while(is.peek() == 'l' || is.peek() =='L' ||
	      is.peek() == 'u' || is.peek() =='U')
	  numberdata += is.get() ;
	AST_data->text = numberdata ;
	return AST_data ;
      }
      if(is.peek() == 'b' || is.peek() == 'B') { // binary number
	numberdata += is.get() ;
	while(is.peek() >= '0' && is.peek() <= '1') {
	  numberdata += is.get() ;
	}
	if(numberdata.size() == 2 || (is.peek() >= '2' && is.peek() <='9'))
	  AST_data->nodeType = AST_type::TK_ERROR ;
	while(is.peek() == 'l' || is.peek() =='L' ||
	      is.peek() == 'u' || is.peek() =='U')
	  numberdata += is.get() ;
	AST_data->text = numberdata ;
	return AST_data ;
      }
      if(is.peek() >= '0' && is.peek() <= '9') { // octal number
	cout << "is.peek() = " << is.peek() << endl ;
	// octal number
	while(is.peek() >= '0' && is.peek() <= '7') {
	  numberdata += is.get() ;
	}
	if(is.peek() >= '8' && is.peek() <='9')
	  AST_data->nodeType = AST_type::TK_ERROR ;
	while(is.peek() == 'l' || is.peek() =='L' ||
	      is.peek() == 'u' || is.peek() =='U')
	  numberdata += is.get() ;
	AST_data->text = numberdata ;
	return AST_data ;
      }
    }
    while(is.peek() >= '0' && is.peek() <= '9') {
      numberdata += is.get() ;
    }
    if(!hasPoint && is.peek() == '.') {
      numberdata += is.get() ;
      hasPoint = true ;
      while(is.peek() >= '0' && is.peek() <= '9') {
	numberdata += is.get() ;
      }
    }
    bool isFloat = hasPoint ;
    if(is.peek() == 'e' || is.peek() == 'E') {
      isFloat = true ;
      numberdata += is.get() ;
      // get exponent
      if(is.peek()=='+' || is.peek() == '-')
	numberdata += is.get() ;
      while(is.peek() >= '0' && is.peek() <= '9') {
	numberdata += is.get() ;
      }
    }
    if(isFloat) {
      if(is.peek()=='f' || is.peek()=='F' || is.peek()=='l' || is.peek()=='L')
	numberdata += is.get() ;
    } else {
      while(is.peek()=='l' || is.peek() == 'L' || is.peek() == 'u' || is.peek()=='U')
	numberdata += is.get() ;
    }
    AST_data->text=numberdata ;
    return AST_data ;
  }
  AST_data->text=numberdata ;
  AST_data->nodeType = AST_type::TK_ERROR ;
  return AST_data ;
}

struct keywords {
  const char *keyword ;
  AST_type::elementType nodeType;
} ;

keywords keywordDictionary[] = {
  {"alignas", AST_type::TK_ALIGNAS},
  {"alignof", AST_type::TK_ALIGNOF},
  {"and", AST_type::TK_LOGICAL_AND},
  {"and_eq", AST_type::TK_AND_ASSIGN},
  {"asm", AST_type::TK_ASM},
  {"auto", AST_type::TK_AUTO},
  {"bitand", AST_type::TK_AND},
  {"bitor", AST_type::TK_OR},
  {"bool", AST_type::TK_BOOL},
  {"break", AST_type::TK_BREAK},
  {"case", AST_type::TK_CASE},
  {"catch", AST_type::TK_CATCH},
  {"char", AST_type::TK_CHAR},
  {"class", AST_type::TK_CLASS},
  {"const", AST_type::TK_CONST},
  {"continue", AST_type::TK_CONTINUE},
  {"default", AST_type::TK_DEFAULT},
  {"delete", AST_type::TK_DELETE},
  {"do", AST_type::TK_DO},
  {"double", AST_type::TK_DOUBLE},
  {"dynamic_cast", AST_type::TK_DYNAMIC_CAST},
  {"else", AST_type::TK_ELSE},
  {"enum", AST_type::TK_ENUM},
  {"explicit", AST_type::TK_EXPLICIT},
  {"export", AST_type::TK_EXPORT},
  {"extern", AST_type::TK_EXTERN},
  {"false", AST_type::TK_FALSE},
  {"float", AST_type::TK_FLOAT},
  {"for", AST_type::TK_FOR},
  {"friend", AST_type::TK_FRIEND},
  {"goto", AST_type::TK_GOTO},
  {"if", AST_type::TK_IF},
  {"inline", AST_type::TK_INLINE},
  {"int", AST_type::TK_INT},
  {"long", AST_type::TK_LONG},
  {"mutable", AST_type::TK_MUTABLE},
  {"namespace", AST_type::TK_NAMESPACE},
  {"new", AST_type::TK_NEW},
  {"noexcept", AST_type::TK_NOEXCEPT},
  {"not", AST_type::TK_NOT},
  {"not_eq", AST_type::TK_EQUAL},
  {"nullptr", AST_type::TK_NULLPTR},
  {"operator", AST_type::TK_OPERATOR},
  {"or", AST_type::TK_LOGICAL_OR},
  {"or_eq", AST_type::TK_OR_ASSIGN},
  {"private", AST_type::TK_PRIVATE},
  {"protected",AST_type::TK_PROTECTED},
  {"public", AST_type::TK_PUBLIC},
  {"register", AST_type::TK_REGISTER},
  {"reinterpret_cast", AST_type::TK_REINTERPRET_CAST},
  {"return", AST_type::TK_RETURN},
  {"short", AST_type::TK_SHORT},
  {"signed", AST_type::TK_SIGNED},
  {"sizeof", AST_type::TK_SIZEOF},
  {"static", AST_type::TK_STATIC},
  {"static_cast", AST_type::TK_STATIC_CAST},
  {"struct", AST_type::TK_STRUCT},
  {"switch", AST_type::TK_SWITCH},
  {"template", AST_type::TK_TEMPLATE},
  {"this", AST_type::TK_THIS},
  {"throw", AST_type::TK_THROW},
  {"true", AST_type::TK_TRUE},
  {"try", AST_type::TK_TRY},
  {"typedef", AST_type::TK_TYPEDEF},
  {"typeid", AST_type::TK_TYPEID},
  {"typename", AST_type::TK_TYPENAME},
  {"union", AST_type::TK_UNION},
  {"unsigned", AST_type::TK_UNSIGNED},
  {"virtual", AST_type::TK_VIRTUAL},
  {"void", AST_type::TK_VOID},
  {"volatile", AST_type::TK_VOLATILE},
  {"while", AST_type::TK_WHILE},
  {"xor",  AST_type::TK_EXOR},
  {"xor_eq", AST_type::TK_EXOR_ASSIGN}
} ;

inline AST_type::elementType findToken(const string &name) {
  const int ntok = sizeof(keywordDictionary)/sizeof(keywords) ;
  int start = ntok/2 ;
  if(name == keywordDictionary[start].keyword) 
    return keywordDictionary[start].nodeType ;
  int min = 0 ;
  int max = ntok-1 ;
  if(name > keywordDictionary[start].keyword) 
    min = start ;
  else
    max = start ;
  while (max-min > 1) {
    start = (min+max)/2 ;
    if(name == keywordDictionary[start].keyword)
      return keywordDictionary[start].nodeType ;
    if(name > keywordDictionary[start].keyword) 
      min = start ;
    else
      max = start ;
  }
  if(name == keywordDictionary[min].keyword)
    return keywordDictionary[min].nodeType ;
  if(name == keywordDictionary[max].keyword)
    return keywordDictionary[max].nodeType ;
  return AST_type::TK_NAME ;
}

CPTR<AST_Token> getTokenInternal(std::istream &is, int &linecount) {
  killsp(is,linecount) ;
  if(is.fail() || is.eof()) {
    CPTR<AST_Token> AST_data = new AST_Token() ;
    AST_data->lineno = linecount ;
    AST_data->nodeType = AST_type::TK_ERROR ;
    return AST_data ;
  }
  if(is_name(is)) {
    CPTR<AST_Token> AST_data = new AST_Token() ;
    AST_data->text = get_name(is) ;
    AST_data->lineno = linecount ;
    AST_data->nodeType = findToken(AST_data->text) ;
#ifdef VERBOSE
    cerr << "get token NAME("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  }
  if(is.peek() == '.' || (is.peek() >= '0' && is.peek() <='9')) {
    return getNumberToken(is,linecount) ;
  }
  CPTR<AST_Token> AST_data = new AST_Token() ;
  AST_data->lineno = linecount ;
  AST_data->text += is.get() ;
  AST_data->nodeType = AST_type::TK_ERROR ;
  switch(AST_data->text[0]) {
  case '+':
    if(is.peek()=='+') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_INCREMENT ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_PLUS_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_PLUS ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '-':
    if(is.peek()=='-') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_DECREMENT ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_MINUS_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek()=='>') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_ARROW ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
//     if(is.peek() >= '0' && is.peek() <='9') {
//       // this is a number
//       CPTR<AST_Token> num = getNumberToken(is,linecount) ;
//       num->text = string("-")+num->text ;
// #ifdef VERBOSE
//       cerr << "get token NUM(" << num->text << ")"<< endl ;
// #endif
//       return num ;
//     }
    AST_data->nodeType = AST_type::TK_MINUS ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '*':
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_TIMES_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_TIMES ;
#ifdef VERBOSE
    cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '/':
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_DIVIDE_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek()=='/') {
      AST_data->text += is.get() ;
      while(is.peek() != '\n' && !is.eof() && !is.fail()) {
	AST_data->text += is.get() ;
      }
      AST_data->nodeType = AST_type::TK_COMMENT ;
#ifdef VERBOSE
      cerr << "get token COMMENT("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek()=='*') {
      AST_data->text += is.get() ;
      for(;;) {
	if(is.eof() || is.fail())
	  break ;
	if(is.peek() != '*') {
	  AST_data->text += is.get() ;
	  if(is.peek() == '/') {
	    AST_data->text = is.get() ;
	    AST_data->nodeType=AST_type::TK_COMMENT ;
#ifdef VERBOSE
	    cerr << "get token COMMENT("<< AST_data->text<< ")" << endl ;
#endif
	    return AST_data ;
	  }
	} else {
	  if(is.peek() == '\n')
	    linecount++ ;
	  AST_data->text += is.get() ;
	}
      }
      AST_data->nodeType=AST_type::TK_ERROR ;
#ifdef VERBOSE
      cerr << "get token ERROR("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_DIVIDE ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '%':
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_MODULUS_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_MODULUS ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case ',':
    AST_data->nodeType = AST_type::TK_COMMA ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '@':
    AST_data->nodeType = AST_type::TK_AT ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '&':
    if(is.peek() == '&' ) {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_LOGICAL_AND ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek() == '=' ) {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_AND_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_AND ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '|':
    if(is.peek() == '|') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_LOGICAL_OR ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek() == '=' ) {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_OR_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_OR ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '^':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_EXOR_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_EXOR ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '=':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_EQUAL ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_ASSIGN ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '!':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_NOT_EQUAL ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_NOT ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case ':':
    if(is.peek() == ':') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_SCOPE ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }

    AST_data->nodeType = AST_type::TK_COLON ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '<':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_LE ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek() == '<') {
      AST_data->text += is.get() ;
      if(is.peek() == '=') {
	AST_data->text += is.get() ;
	AST_data->nodeType = AST_type::TK_SHIFT_LEFT_ASSIGN ;
#ifdef VERBOSE
	cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
	return AST_data ;
      }
      AST_data->nodeType = AST_type::TK_SHIFT_LEFT ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_LT ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '>':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_GE ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek() == '>') {
      AST_data->text += is.get() ;
      if(is.peek() == '=') {
	AST_data->text += is.get() ;
	AST_data->nodeType = AST_type::TK_SHIFT_RIGHT_ASSIGN ;
#ifdef VERBOSE
	cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
	return AST_data ;
      }
      AST_data->nodeType = AST_type::TK_SHIFT_RIGHT ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_GT ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '~':
    AST_data->nodeType = AST_type::TK_TILDE ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '?':
    AST_data->nodeType = AST_type::TK_QUESTION ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '"':
    while(is.peek() != '"' && !is.eof() && !is.fail()) {
      AST_data->text += is.get() ;
    }
    if(is.peek() == '"') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_STRING ;
#ifdef VERBOSE
      cerr << "get token STRING("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    break ;
  case '\'':
    while(is.peek() != '\'' && !is.eof() && !is.fail()) {
      AST_data->text += is.get() ;
    }
    if(is.peek() == '\'') {
      AST_data->text += is.get() ;
      AST_data->nodeType = AST_type::TK_STRING ;
#ifdef VERBOSE
      cerr << "get token STRING("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    break ;
  case ';':
    AST_data->nodeType = AST_type::TK_SEMICOLON ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '(':
    AST_data->nodeType = AST_type::TK_OPENPAREN ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case ')':
    AST_data->nodeType = AST_type::TK_CLOSEPAREN ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '[':
    AST_data->nodeType = AST_type::TK_OPENBRACKET ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case ']':
    AST_data->nodeType = AST_type::TK_CLOSEBRACKET ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '{':
    AST_data->nodeType = AST_type::TK_OPENBRACE ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '}':
    AST_data->nodeType = AST_type::TK_CLOSEBRACE ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '$':
    // Now we have a Loci variable or a Loci command
    if(is.peek() == '[') { // This is a Loci command
      AST_data->text += is.get() ;
      AST_data->text = "" ;
      while(!is.fail() && !is.eof() && is.peek() != ']') {
	if(is.peek() == '[') {
	  int count = 1 ;
	  while(!is.fail() && !is.eof() && count != 0) {
	    if(is.peek() == '[')
	      count++ ;
	    if(is.peek() == ']')
	      count-- ;
	    AST_data->text += is.get() ;
	  }
	} else 
	  AST_data->text += is.get() ;
      }
      is.get() ;
      AST_data->nodeType = AST_type::TK_LOCI_DIRECTIVE ;
#ifdef VERBOSE
      cerr << "get token DIRECTIVE("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = AST_type::TK_LOCI_VARIABLE ;
    if(is.peek() == '*') {
      AST_data->nodeType = AST_type::TK_LOCI_CONTAINER ;
      AST_data->text += is.peek() ;
    }
    if(!is.fail() && !is.eof() &&
       ((is.peek() >='a' && is.peek() <='z') ||
	(is.peek() >='A' && is.peek() <='Z') ||
	is.peek() == '_' || is.peek() == '@')) {
      // extract variable name
      while(!is.fail() && !is.eof() &&
	    ((is.peek() >='a' && is.peek() <='z') ||
	     (is.peek() >='A' && is.peek() <='Z') ||
	     (is.peek() >='0' && is.peek() <='9') ||
	     is.peek() == '_' || is.peek() == '@' ||
	     is.peek() == '(' || is.peek() == '{')) {
	if(is.peek() == '(') {
	  AST_data->text +=is.get() ;
	  int count = 1 ;
	  while(!is.fail() && !is.eof() && count != 0) {
	    if(is.peek() == '(')
	      count++ ;
	    if(is.peek() == ')')
	      count-- ;
	    AST_data->text += is.get() ;
	  }
	} else if(is.peek() == '{') {
	  AST_data->text +=is.get() ;
	  int count = 1 ;
	  while(!is.fail() && !is.eof() && count != 0) {
	    if(is.peek() == '{')
	      count++ ;
	    if(is.peek() == '}')
	      count-- ;
	    AST_data->text += is.get() ;
	  }
	} else
	  AST_data->text += is.get() ;
      }

    }
#ifdef VERBOSE
    cerr << "get token VAR("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  } 
     
      
  AST_data->nodeType = AST_type::TK_ERROR ;
#ifdef VERBOSE
  cerr << "get token ERR("<< AST_data->text<< ")" << endl ;
#endif
  return AST_data ;

}
  
    
CPTR<AST_Token> getToken(std::istream &is, int &linecount) {
  if(tokenStack.size() != 0) {
    CPTR<AST_Token> tok = tokenStack.back() ;
#ifdef VERBOSE
    cerr << "poping token '" << tok->text << "'" << endl ;
#endif
    tokenStack.pop_back() ;
    return tok ;
  }
  //  return getTokenInternal(is,linecount) ;
  CPTR<AST_Token> tok  = getTokenInternal(is,linecount) ;

  if(tok->nodeType == AST_type::TK_LT) {
    // Disambiguate between '<' operator and template arguments
    // Search forward to find a matching '>'
#ifdef VERBOSE
    cerr << "searching to find if '<' is for template arguments or less than operator" << endl ;
#endif
    vector<int> LT_LIST ;
    vector<int> LT_DEPTH ;
    vector<CPTR<AST_Token> > toklist ;
    toklist.push_back(tok) ;
    LT_LIST.push_back(toklist.size()-1 ) ;
    int parencount = 0 ;
    LT_DEPTH.push_back(parencount) ;
    do {
      tok = getTokenInternal(is,linecount) ;
      toklist.push_back(tok) ;
      switch(tok->nodeType) {
      case AST_type::TK_SEMICOLON:
      case AST_type::TK_LOGICAL_AND:
      case AST_type::TK_LOGICAL_OR:
      case AST_type::TK_OPENBRACKET:
      case AST_type::TK_OPENBRACE:
      case AST_type::TK_CLOSEBRACKET:
      case AST_type::TK_CLOSEBRACE:
      case AST_type::TK_ERROR:
#ifdef VERBOSE
        cerr << "search did not find matching '>', found '" << tok->text << "' instead." << endl ;
#endif
        { // No matching bracket found so  unwind search
          // onto token stack
          while(toklist.size()>1) {
            pushToken(toklist.back()) ;
            toklist.pop_back() ;
          }
          // return TK_LT or TK_OPENTEMPLATE
          tok = toklist[0] ;
          return tok ;
        }
        break ;
      case AST_type::TK_OPENPAREN:
        parencount++ ;
        break ;
      case AST_type::TK_CLOSEPAREN:
        if(parencount == 0) {
#ifdef VERBOSE
        cerr << "search did not find matching '>', found '" << tok->text << "' instead." << endl ;
#endif
          // No matching bracket found so  unwind search
          // onto token stack
          while(toklist.size()>1) {
            pushToken(toklist.back()) ;
            toklist.pop_back() ;
          }
          // return TK_LT or TK_OPENTEMPLATE
          tok = toklist[0] ;
          return tok ;
        }
        parencount-- ;
        break ;
      case AST_type::TK_LT:
        LT_LIST.push_back(toklist.size()-1 ) ;
        LT_DEPTH.push_back(parencount) ;
#ifdef VERBOSE
        cerr << "searching for nested template '<'" << endl ;
#endif
        break;
      case AST_type::TK_GT:
        // found match to '<'
        if(parencount == LT_DEPTH.back()) {
#ifdef VERBOSE
          cerr << "found matching '>'" << endl ;
#endif
          toklist[LT_LIST.back()]->nodeType =
            AST_type::TK_OPENTEMPLATE ;
          toklist[LT_LIST.back()]->text = "<<<" ; ;
          toklist.back()->nodeType =
            AST_type::TK_CLOSETEMPLATE ;
          toklist.back()->text = ">>>" ;
          LT_LIST.pop_back() ;
          LT_DEPTH.pop_back() ;
          if(LT_LIST.size() == 0) {
            while(toklist.size()>1) {
              pushToken(toklist.back()) ;
              toklist.pop_back() ;
            }
            // return TK_OPENTEMPLATE
            tok = toklist[0] ;
            return tok ;
          }
        }
        break ;
      case AST_type::TK_SHIFT_RIGHT:
        // find out if we have a pair of < < in the same paren
        // level to match
        {
#ifdef VERBOSE
          cerr << "found matching '>>'" << endl ;
#endif
          int sz = LT_LIST.size() ;
          if(sz>1 && LT_DEPTH[sz-1] == parencount &&
             LT_DEPTH[sz-2] == parencount) {
            toklist[LT_LIST[sz-1]]->nodeType =
              AST_type::TK_OPENTEMPLATE ;
            toklist[LT_LIST[sz-2]]->nodeType =
              AST_type::TK_OPENTEMPLATE ;
            toklist.back()->nodeType =
              AST_type::TK_CLOSETEMPLATE ;
            toklist.back()->text = ">" ;
            CPTR<AST_Token> AST_data = new AST_Token() ;
            AST_data->text = ">" ; ;
            AST_data->lineno = toklist.back()->lineno ;
            AST_data->nodeType = AST_type::TK_CLOSETEMPLATE ;
            toklist.push_back(AST_data) ;
            LT_LIST.pop_back() ;
            LT_DEPTH.pop_back() ;
            LT_LIST.pop_back() ;
            LT_DEPTH.pop_back() ;
          }
          if(LT_LIST.size() == 0) {
            while(toklist.size()>1) {
              pushToken(toklist.back()) ;
              toklist.pop_back() ;
            }
            // return TK_OPENTEMPLATE
            tok = toklist[0] ;
            return tok ;
          }
        }
      default:
        break ;
      }
    } while(true) ;
  }
  return tok ;


}
