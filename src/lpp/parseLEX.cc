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
  {"using",AST_type::TK_USING},
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
	is.peek() == '_' || is.peek() == '@' ||
        is.peek() == '$')
       ) {
      // extract variable name
      while(!is.fail() && !is.eof() &&
	    ((is.peek() >='a' && is.peek() <='z') ||
	     (is.peek() >='A' && is.peek() <='Z') ||
	     (is.peek() >='0' && is.peek() <='9') ||
	     is.peek() == '_' || is.peek() == '@' ||
	     is.peek() == '(' || is.peek() == '{' ||
             is.peek() == '$')) {
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
/// Convert operator code to corresponding string for printing
string OPtoName(AST_type::elementType val) {
  switch(val) {
  case AST_type::OP_SCOPE:
    return string("OP_SCOPE") ;
  case AST_type::OP_AT:
    return string("OP_AT") ;
  case AST_type::OP_DOT:
    return string("OP_DOT") ;
  case AST_type::OP_ARROW:
    return string("OP_ARROW") ;
  case AST_type::OP_TIMES:
    return string("OP_TIMES") ;
  case AST_type::OP_DIVIDE:
    return string("OP_DIVIDE") ;
  case AST_type::OP_MODULUS:
    return string("OP_MODULUS") ;
  case AST_type::OP_PLUS:
    return string("OP_PLUS") ;
  case AST_type::OP_MINUS:
    return string("OP_MINUS") ;
  case AST_type::OP_SHIFT_RIGHT:
    return string("OP_SHIFT_RIGHT") ;
  case AST_type::OP_SHIFT_LEFT:
    return string("OP_SHIFT_RIGHT") ;
  case AST_type::OP_LT:
    return string("OP_LT") ;
  case AST_type::OP_GT:
    return string("OP_GT") ;
  case AST_type::OP_GE:
    return string("OP_GE") ;
  case AST_type::OP_LE:
    return string("OP_LE") ;
  case AST_type::OP_EQUAL:
    return string("OP_EQUAL") ;
  case AST_type::OP_NOT_EQUAL:
    return string("OP_NOT_EQUAL") ;
  case AST_type::OP_AND:
    return string("OP_AND") ;
  case AST_type::OP_EXOR:
    return string("OP_EXOR") ;
  case AST_type::OP_OR:
    return string("OP_OR") ;
  case AST_type::OP_LOGICAL_AND:
    return string("OP_LOGICAL_AND") ;
  case AST_type::OP_LOGICAL_OR:
    return string("OP_LOGICAL_OR") ;
  case AST_type::OP_TERNARY:
    return string("OP_TERNARY") ;
  case AST_type::OP_ASSIGN:
    return string("OP_ASSIGN") ;
  case AST_type::OP_TIMES_ASSIGN:
    return string("OP_TIMES_ASSIGN") ;
  case AST_type::OP_DIVIDE_ASSIGN:
    return string("OP_DIVIDE_ASSIGN") ;
  case AST_type::OP_MODULUS_ASSIGN:
    return string("OP_MODULUS_ASSIGN") ;
  case AST_type::OP_PLUS_ASSIGN:
    return string("OP_PLUS_ASSIGN") ;
  case AST_type::OP_MINUS_ASSIGN:
    return string("OP_MINUS_ASSIGN") ;
  case AST_type::OP_SHIFT_LEFT_ASSIGN:
    return string("OP_SHIFT_LEFT_ASSIGN") ;
  case AST_type::OP_SHIFT_RIGHT_ASSIGN:
    return string("OP_SHIFT_RIGHT_ASSIGN") ;
  case AST_type::OP_AND_ASSIGN:
    return string("OP_AND_ASSIGN") ;
  case AST_type::OP_OR_ASSIGN:
    return string("OP_OR_ASSIGN") ;
  case AST_type::OP_EXOR_ASSIGN:
    return string("OP_EXOR_ASSIGN") ;
  case AST_type::OP_COMMA:
    return string("OP_COMMA") ;
  case AST_type::OP_COLON:
    return string("OP_COLON") ;
  case AST_type::OP_SEMICOLON:
    return string("OP_SEMICOLON") ;
  case AST_type::OP_NIL:
    return string("OP_NIL") ;
  case AST_type::OP_INCREMENT:
    return string("OP_INCREMENT") ;
  case AST_type::OP_DECREMENT:
    return string("OP_DECREMENT") ;
  case AST_type::OP_POSTINCREMENT:
    return string("OP_POSTINCREMENT") ;
  case AST_type::OP_POSTDECREMENT:
    return string("OP_POSTDECREMENT") ;
  case AST_type::OP_COMMENT:
    return string("OP_COMMENT") ;
  case AST_type::OP_BRACEBLOCK:
    return string("OP_BRACEBLOCK") ;
  case AST_type::OP_NAME:
    return string("OP_NAME") ;
  case AST_type::OP_FUNC:
    return string("OP_FUNC") ;
  case AST_type::OP_ARRAY:
    return string("OP_ARRAY") ;
  case AST_type::OP_NAME_BRACE:
    return string("OP_NAME_BRACE") ;
  case AST_type::OP_FUNC_BRACE:
    return string("OP_FUNC_BRACE") ;
  case AST_type::OP_TEMPLATE:
    return string("OP_TEMPLATE") ;
  case AST_type::OP_TEMPLATE_FUNC:
    return string("OP_TEMPLATE_FUNC") ;
  case AST_type::OP_STRING:
    return string("OP_STRING") ;
  case AST_type::OP_NUMBER:
    return string("OP_NUMBER") ;
  case AST_type::OP_ERROR:
    return string("OP_ERROR") ;
    
  case AST_type::OP_UNARY_PLUS:
    return string("OP_UNARY_PLUS") ;
  case AST_type::OP_UNARY_MINUS:
    return string("OP_UNARY_MINUS") ;
  case AST_type::OP_NOT:
    return string("OP_NOT") ;
  case AST_type::OP_TILDE:
    return string("OP_TILDE") ;
  case AST_type::OP_AMPERSAND:
    return string("OP_AMPERSAND") ;
  case AST_type::OP_DOLLAR:
    return string("OP_DOLLAR") ;
  case AST_type::OP_STAR:
    return string("OP_STAR") ;
  case AST_type::OP_GROUP:
    return string("OP_GROUP") ;
  case AST_type::OP_GROUP_ERROR:
    return string("OP_GROUP_ERROR") ;
  case AST_type::OP_OPENPAREN:
    return string("OP_OPENPAREN") ;
  case AST_type::OP_CLOSEPAREN:
    return string("OP_CLOSEPAREN") ;
  case AST_type::OP_OPENBRACKET:
    return string("OP_OPENBRACKET") ;
  case AST_type::OP_CLOSEBRACKET:
    return string("OP_CLOSEBRACKET") ;
  case AST_type::OP_OPENBRACE:
    return string("OP_OPENBRACE") ;
  case AST_type::OP_CLOSEBRACE:
    return string("OP_CLOSEBRACE") ;
  case AST_type::OP_LOCI_DIRECTIVE:
    return string("OP_LOCI_DIRECTIVE") ;
  case AST_type::OP_LOCI_VARIABLE:
    return string("OP_LOCI_VARIABLE") ;
  case AST_type::OP_LOCI_CONTAINER:
    return string("OP_LOCI_CONTAINER") ;
  case AST_type::OP_TERM:
    return string("OP_TERM") ;
  case AST_type::OP_SPECIAL:
    return string("OP_SPECIAL") ;
  case AST_type::ND_BLOCK:
    return string("ND_BLOCK") ;
  case AST_type::ND_SYNTAXERR:
    return string("ND_SYTNAXERR") ;
  case AST_type::ND_CTRL_IF:
    return string("ND_CTRL_IF") ;
  case AST_type::ND_CTRL_FOR:
    return string("ND_CTRL_FOR") ;
  case AST_type::ND_CTRL_WHILE:
    return string("ND_CTRL_WHILE") ;
  case AST_type::ND_CTRL_DO:
    return string("ND_CTRL_DO") ;
  case AST_type::ND_CTRL_SWITCH:
    return string("ND_CTRL_SWITCH") ;
  case AST_type::ND_SIMPLE_STATEMENT:
    return string("ND_SIMPLE_STATEMENT") ;
  case AST_type::ND_DECL:
    return string("ND_DECL") ;
  case AST_type::ND_TYPE_SPEC:
    return string("ND_TYPE_SPEC") ;
  case AST_type::ND_TERMINAL:
    return string("ND_TERMINAL") ;
  case AST_type::TK_SENTINEL:
    return string("TK_SENTINEL") ;
  case AST_type::TK_BRACEBLOCK:
    return string("TK_BRACEBLOCK") ;
  case AST_type::TK_SCOPE:
    return string("TK_SCOPE") ;
  case AST_type::TK_NAME:
    return string("TK_NAME") ;
  case AST_type::TK_AT:
    return string("TK_AT") ;
  case AST_type::TK_ARROW:
    return string("TK_ARROW") ;
  case AST_type::TK_TIMES:
    return string("TK_TIMES") ;
  case AST_type::TK_DIVIDE:
    return string("TK_DIVIDE") ;
  case AST_type::TK_MODULUS:
    return string("TK_MODULUS") ;
  case AST_type::TK_PLUS:
    return string("TK_PLUS") ;
  case AST_type::TK_MINUS:
    return string("TK_MINUS") ;
  case AST_type::TK_SHIFT_RIGHT:
    return string("TK_SHIFT_RIGHT") ;
  case AST_type::TK_SHIFT_LEFT:
    return string("TK_SHIFT_LEFT") ;
  case AST_type::TK_LT:
    return string("TK_LT") ;
  case AST_type::TK_GT:
    return string("TK_GT") ;
  case AST_type::TK_GE:
    return string("TK_GE") ;
  case AST_type::TK_LE:
    return string("TK_LE") ;
  case AST_type::TK_EQUAL:
    return string("TK_EQUAL") ;
  case AST_type::TK_NOT_EQUAL:
    return string("TK_NOT_EQUAL") ;
  case AST_type::TK_AND:
    return string("TK_AND") ;
  case AST_type::TK_EXOR:
    return string("TK_EXOR") ;
  case AST_type::TK_OR:
    return string("TK_OR") ;
  case AST_type::TK_LOGICAL_AND:
    return string("TK_LOGICAL_AND") ;
  case AST_type::TK_LOGICAL_OR:
    return string("TK_LOGICAL_OR") ;
  case AST_type::TK_ASSIGN:
    return string("TK_ASSIGN") ;
  case AST_type::TK_TIMES_ASSIGN:
    return string("TK_TIMES_ASSIGN") ;
  case AST_type::TK_DIVIDE_ASSIGN:
    return string("TK_DIVIDE_ASSIGN") ;
  case AST_type::TK_MODULUS_ASSIGN:
    return string("TK_MODULUS_ASSIGN") ;
  case AST_type::TK_PLUS_ASSIGN:
    return string("TK_PLUS_ASSIGN") ;
  case AST_type::TK_MINUS_ASSIGN:
    return string("TK_MINUS_ASSIGN") ;
  case AST_type::TK_SHIFT_LEFT_ASSIGN:
    return string("TK_SHIFT_LEFT_ASSIGN") ;
  case AST_type::TK_SHIFT_RIGHT_ASSIGN:
    return string("TK_SHIFT_RIGHT_ASSIGN") ;
  case AST_type::TK_AND_ASSIGN:
    return string("TK_AND_ASSIGN") ;
  case AST_type::TK_OR_ASSIGN:
    return string("TK_OR_ASSIGN") ;
  case AST_type::TK_EXOR_ASSIGN:
    return string("TK_EXOR_ASSIGN") ;
  case AST_type::TK_COMMA:
    return string("TK_COMMA") ;
  case AST_type::TK_DOT:
    return string("TK_DOT") ;
  case AST_type::TK_COLON:
    return string("TK_COLON") ;
  case AST_type::TK_SEMICOLON:
    return string("TK_SEMICOLON") ;
  case AST_type::TK_NIL:
    return string("TK_NIL") ;
  case AST_type::TK_INCREMENT:
    return string("TK_INCREMENT") ;
  case AST_type::TK_DECREMENT:
    return string("TK_DECREMENT") ;
  case AST_type::TK_COMMENT:
    return string("TK_COMMENT") ;
  case AST_type::TK_STRING:
    return string("TK_STRING") ;
  case AST_type::TK_NUMBER:
    return string("TK_NUMBER") ;
  case AST_type::TK_ERROR:
    return string("TK_ERROR") ;
  case AST_type::TK_UNARY_PLUS:
    return string("TK_UNARY_PLUS") ;
  case AST_type::TK_UNARY_MINUS:
    return string("TK_UNARY_MINUS") ;
  case AST_type::TK_NOT:
    return string("TK_NOT") ;
  case AST_type::TK_TILDE:
    return string("TK_TILDE") ;
  case AST_type::TK_QUESTION:
    return string("TK_QUESTION") ;
  case AST_type::TK_AMPERSAND:
    return string("TK_AMPERSAND") ;
  case AST_type::TK_STAR:
    return string("TK_STAR") ;
  case AST_type::TK_OPENPAREN:
    return string("TK_OPENPAREN") ;
  case AST_type::TK_CLOSEPAREN:
    return string("TK_CLOSEPAREN") ;
  case AST_type::TK_OPENBRACKET:
    return string("TK_OPENBRACKET") ;
  case AST_type::TK_CLOSEBRACKET:
    return string("TK_CLOSEBRACKET") ;
  case AST_type::TK_OPENBRACE:
    return string("TK_CLOSEBRACE") ;
  case AST_type::TK_OPENTEMPLATE:
    return string("TK_CLOSETEMPLATE") ;
  case AST_type::TK_LOCI_DIRECTIVE:
    return string("TK_LOCI_DIRECTIVE") ;
  case AST_type::TK_LOCI_VARIABLE:
    return string("TK_LOCI_VARIABLE") ;
  case AST_type::TK_LOCI_CONTAINER:
    return string("TK_LOCI_CONTAINER") ;

  case AST_type::TK_ALIGNAS:
    return string("TK_ALIGNAS") ;

  case AST_type::TK_ALIGNOF:
    return string("TK_ALIGNOF") ;

  case AST_type::TK_ASM: 
    return string("TK_ASM") ; 

  case AST_type::TK_BOOL:
    return string("TK_BOOL") ;

  case AST_type::TK_FALSE:
    return string("TK_FALSE") ;

  case AST_type::TK_TRUE:
    return string("TK_TRUE") ;

  case AST_type::TK_CHAR:
    return string("TK_CHAR") ;

  case AST_type::TK_INT:
    return string("TK_INT") ;

  case AST_type::TK_LONG:
    return string("TK_LONG") ;
      
  case AST_type::TK_SHORT:
    return string("TK_SHORT") ;

  case AST_type::TK_SIGNED:
    return string("TK_SIGNED") ;

  case AST_type::TK_UNSIGNED:
    return string("TK_UNSIGNED") ;

  case AST_type::TK_DOUBLE:
    return string("TK_DOUBLE") ;

  case AST_type::TK_FLOAT:
    return string("TK_FLOAT") ;

  case AST_type::TK_ENUM:
    return string("TK_ENUM") ;

  case AST_type::TK_MUTABLE:
    return string("TK_MUTABLE") ;

  case AST_type::TK_CONST:
    return string("TK_CONST") ;

  case AST_type::TK_STATIC:
    return string("TK_STATIC") ;

  case AST_type::TK_VOLATILE:
    return string("TK_VOLATILE") ;

  case AST_type::TK_AUTO:
    return string("TK_AUTO") ;

  case AST_type::TK_REGISTER:
    return string("TK_REGISTER") ;

  case AST_type::TK_EXPORT:
    return string("TK_EXPORT") ;

  case AST_type::TK_EXTERN:
    return string("TK_EXTERN") ;

  case AST_type::TK_INLINE:
    return string("TK_INLINE") ;

  case AST_type::TK_NAMESPACE:
    return string("TK_NAMESPACE") ;

  case AST_type::TK_EXPLICIT:
    return string("TK_EXPLICIT") ;

  case AST_type::TK_DYNAMIC_CAST:
    return string("TK_DYNAMIC_CAST") ;

  case AST_type::TK_STATIC_CAST:
    return string("TK_STATIC_CAST") ;

  case AST_type::TK_REINTERPRET_CAST:
    return string("TK_REINTERPRET_CAST") ;

  case AST_type::TK_OPERATOR:
    return string("TK_OPERATOR") ;

  case AST_type::TK_PROTECTED:
    return string("TK_PROTECTED") ;

  case AST_type::TK_NOEXCEPT:
    return string("TK_NOEXCEPT") ;

  case AST_type::TK_NULLPTR:
    return string("TK_NULLPTR") ;

  case AST_type::TK_RETURN:
    return string("TK_RETURN") ;

  case AST_type::TK_SIZEOF:
    return string("TK_SIZEOF") ;

  case AST_type::TK_THIS:
    return string("TK_THIS") ;

  case AST_type::TK_TYPEID:
    return string("TK_TYPEID") ;

  case AST_type::TK_SWITCH:
    return string("TK_SWITCH") ;

  case AST_type::TK_CASE:
    return string("TK_CASE") ;

  case AST_type::TK_BREAK:
    return string("TK_BREAK") ;

  case AST_type::TK_DEFAULT:
    return string("TK_DEFAULT") ;

  case AST_type::TK_FOR:
    return string("TK_FOR") ;

  case AST_type::TK_DO:
    return string("TK_DO") ;

  case AST_type::TK_WHILE:
    return string("TK_WHILE") ;

  case AST_type::TK_CONTINUE:
    return string("TK_CONTINUE") ;

  case AST_type::TK_CLASS:
    return string("TK_CLASS") ;

  case AST_type::TK_STRUCT:
    return string("TK_STRUCT") ;

  case AST_type::TK_PUBLIC:
    return string("TK_PUBLIC") ;

  case AST_type::TK_PRIVATE:
    return string("TK_PRIVATE") ;

  case AST_type::TK_FRIEND:
    return string("TK_FRIEND") ;

  case AST_type::TK_UNION:
    return string("TK_UNION") ;

  case AST_type::TK_TYPENAME:
    return string("TK_TYPENAME") ;

  case AST_type::TK_TEMPLATE:
    return string("TK_TEMPLATE") ;

  case AST_type::TK_TYPEDEF:
    return string("TK_TYPEDEF") ;

  case AST_type::TK_VIRTUAL:
    return string("TK_VIRTUAL") ;

  case AST_type::TK_VOID:
    return string("TK_VOID") ;

  case AST_type::TK_TRY:
    return string("TK_TRY") ;

  case AST_type::TK_CATCH:
    return string("TK_CATCH") ;

  case AST_type::TK_THROW:
    return string("TK_THROW") ;

  case AST_type::TK_IF:
    return string("TK_IF") ;

  case AST_type::TK_ELSE:
    return string("TK_ELSE") ;

  case AST_type::TK_GOTO:
    return string("TK_GOTO") ;

  case AST_type::TK_NEW:
    return string("TK_NEW") ;

  case AST_type::TK_DELETE:
    return string("TK_DELETE") ;
    
  default:
    {
      std::ostringstream oss ;
      oss << "ND_UNDEF(0x" << std::hex << int(val) << ")" ;
      return oss.str() ;
    }
  }
  return string("/*error*/") ;
}

