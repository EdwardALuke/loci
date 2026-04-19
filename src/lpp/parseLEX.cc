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
using namespace nodeTypes ;
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
      AST_data->nodeType = TK_DOT ;
      AST_data->text = numberdata ;
      return AST_data ;
    }
    hasPoint = true ;
  }
  if((is.peek() >= '0' && is.peek() <='9')) {
    AST_data->nodeType = TK_NUMBER ;
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
	  AST_data->nodeType = TK_ERROR ;
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
	  AST_data->nodeType = TK_ERROR ;
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
	  AST_data->nodeType = TK_ERROR ;
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
  AST_data->nodeType = TK_ERROR ;
  return AST_data ;
}

struct keywords {
  const char *keyword ;
  AST_type::elementType nodeType;
} ;

keywords keywordDictionary[] = {
  {"alignas", TK_ALIGNAS},
  {"alignof", TK_ALIGNOF},
  {"and", TK_LOGICAL_AND},
  {"and_eq", TK_AND_ASSIGN},
  {"asm", TK_ASM},
  {"auto", TK_AUTO},
  {"bitand", TK_AND},
  {"bitor", TK_OR},
  {"bool", TK_BOOL},
  {"break", TK_BREAK},
  {"case", TK_CASE},
  {"catch", TK_CATCH},
  {"char", TK_CHAR},
  {"class", TK_CLASS},
  {"const", TK_CONST},
  {"continue", TK_CONTINUE},
  {"default", TK_DEFAULT},
  {"delete", TK_DELETE},
  {"do", TK_DO},
  {"double", TK_DOUBLE},
  {"dynamic_cast", TK_DYNAMIC_CAST},
  {"else", TK_ELSE},
  {"enum", TK_ENUM},
  {"explicit", TK_EXPLICIT},
  {"export", TK_EXPORT},
  {"extern", TK_EXTERN},
  {"false", TK_FALSE},
  {"float", TK_FLOAT},
  {"for", TK_FOR},
  {"friend", TK_FRIEND},
  {"goto", TK_GOTO},
  {"if", TK_IF},
  {"inline", TK_INLINE},
  {"int", TK_INT},
  {"long", TK_LONG},
  {"mutable", TK_MUTABLE},
  {"namespace", TK_NAMESPACE},
  {"new", TK_NEW},
  {"noexcept", TK_NOEXCEPT},
  {"not", TK_NOT},
  {"not_eq", TK_EQUAL},
  {"nullptr", TK_NULLPTR},
  {"operator", TK_OPERATOR},
  {"or", TK_LOGICAL_OR},
  {"or_eq", TK_OR_ASSIGN},
  {"private", TK_PRIVATE},
  {"protected",TK_PROTECTED},
  {"public", TK_PUBLIC},
  {"register", TK_REGISTER},
  {"reinterpret_cast", TK_REINTERPRET_CAST},
  {"return", TK_RETURN},
  {"short", TK_SHORT},
  {"signed", TK_SIGNED},
  {"sizeof", TK_SIZEOF},
  {"static", TK_STATIC},
  {"static_cast", TK_STATIC_CAST},
  {"struct", TK_STRUCT},
  {"switch", TK_SWITCH},
  {"template", TK_TEMPLATE},
  {"this", TK_THIS},
  {"throw", TK_THROW},
  {"true", TK_TRUE},
  {"try", TK_TRY},
  {"typedef", TK_TYPEDEF},
  {"typeid", TK_TYPEID},
  {"typename", TK_TYPENAME},
  {"union", TK_UNION},
  {"unsigned", TK_UNSIGNED},
  {"using",TK_USING},
  {"virtual", TK_VIRTUAL},
  {"void", TK_VOID},
  {"volatile", TK_VOLATILE},
  {"while", TK_WHILE},
  {"xor",  TK_EXOR},
  {"xor_eq", TK_EXOR_ASSIGN}
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
  return TK_NAME ;
}

CPTR<AST_Token> getTokenInternal(std::istream &is, int &linecount) {
  killsp(is,linecount) ;
  if(is.fail() || is.eof()) {
    CPTR<AST_Token> AST_data = new AST_Token() ;
    AST_data->lineno = linecount ;
    AST_data->nodeType = TK_ERROR ;
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
  AST_data->nodeType = TK_ERROR ;
  switch(AST_data->text[0]) {
  case '+':
    if(is.peek()=='+') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_INCREMENT ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_PLUS_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_PLUS ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '-':
    if(is.peek()=='-') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_DECREMENT ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_MINUS_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek()=='>') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_ARROW ;
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
    AST_data->nodeType = TK_MINUS ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '*':
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_TIMES_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_TIMES ;
#ifdef VERBOSE
    cerr << "get token OPER("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '/':
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_DIVIDE_ASSIGN ;
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
      AST_data->nodeType = TK_COMMENT ;
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
	    AST_data->nodeType=TK_COMMENT ;
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
      AST_data->nodeType=TK_ERROR ;
#ifdef VERBOSE
      cerr << "get token ERROR("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_DIVIDE ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '%':
    if(is.peek()=='=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_MODULUS_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_MODULUS ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case ',':
    AST_data->nodeType = TK_COMMA ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '@':
    AST_data->nodeType = TK_AT ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '&':
    if(is.peek() == '&' ) {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_LOGICAL_AND ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek() == '=' ) {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_AND_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_AND ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '|':
    if(is.peek() == '|') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_LOGICAL_OR ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek() == '=' ) {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_OR_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_OR ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '^':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_EXOR_ASSIGN ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_EXOR ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '=':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_EQUAL ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_ASSIGN ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '!':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_NOT_EQUAL ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_NOT ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case ':':
    if(is.peek() == ':') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_SCOPE ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }

    AST_data->nodeType = TK_COLON ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '<':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_LE ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek() == '<') {
      AST_data->text += is.get() ;
      if(is.peek() == '=') {
	AST_data->text += is.get() ;
	AST_data->nodeType = TK_SHIFT_LEFT_ASSIGN ;
#ifdef VERBOSE
	cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
	return AST_data ;
      }
      AST_data->nodeType = TK_SHIFT_LEFT ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_LT ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '>':
    if(is.peek() == '=') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_GE ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    if(is.peek() == '>') {
      AST_data->text += is.get() ;
      if(is.peek() == '=') {
	AST_data->text += is.get() ;
	AST_data->nodeType = TK_SHIFT_RIGHT_ASSIGN ;
#ifdef VERBOSE
	cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
	return AST_data ;
      }
      AST_data->nodeType = TK_SHIFT_RIGHT ;
#ifdef VERBOSE
      cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_GT ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '~':
    AST_data->nodeType = TK_TILDE ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '?':
    AST_data->nodeType = TK_QUESTION ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '"':
    {
      while(is.peek() != '"' && !is.eof() && !is.fail()) {
        AST_data->text += is.get() ;
      }
      if(is.peek() != '"') {
        break ;
      }
      AST_data->text += is.get() ;
      // C++: concatenate adjacent string literals ("a" "b" is one literal "ab")
      for(;;) {
        killsp(is, linecount) ;
        if(is.peek() != '"' || is.eof() || is.fail()) {
          break ;
        }
        if(AST_data->text.empty() || AST_data->text.back() != '"') {
          break ;
        }
        AST_data->text.pop_back() ;
        is.get() ; // opening '"' of the next literal
        while(is.peek() != '"' && !is.eof() && !is.fail()) {
          AST_data->text += is.get() ;
        }
        if(is.peek() != '"') {
          AST_data->nodeType = TK_ERROR ;
#ifdef VERBOSE
          cerr << "get token ERR(unclosed string in concatenation)" << endl ;
#endif
          return AST_data ;
        }
        AST_data->text += is.get() ;
      }
      AST_data->nodeType = TK_STRING ;
#ifdef VERBOSE
      cerr << "get token STRING("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
  case '\'':
    while(is.peek() != '\'' && !is.eof() && !is.fail()) {
      AST_data->text += is.get() ;
    }
    if(is.peek() == '\'') {
      AST_data->text += is.get() ;
      AST_data->nodeType = TK_STRING ;
#ifdef VERBOSE
      cerr << "get token STRING("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    break ;
  case ';':
    AST_data->nodeType = TK_SEMICOLON ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '(':
    AST_data->nodeType = TK_OPENPAREN ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case ')':
    AST_data->nodeType = TK_CLOSEPAREN ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '[':
    AST_data->nodeType = TK_OPENBRACKET ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case ']':
    AST_data->nodeType = TK_CLOSEBRACKET ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '{':
    AST_data->nodeType = TK_OPENBRACE ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '}':
    AST_data->nodeType = TK_CLOSEBRACE ;
#ifdef VERBOSE
    cerr << "get token OP("<< AST_data->text<< ")" << endl ;
#endif
    return AST_data ;
  case '#':
    {
      AST_data->text = "" ;
      while(!is.fail() && !is.eof() && is.peek() != '\n') {
        char ch = is.get() ;
        if(ch == '\\' && is.peek() == '\n') {
          AST_data->text += ch ;
          ch = is.get() ;
          AST_data->text += ch ;
          linecount++ ;
          continue ;
        }
        AST_data->text += ch ;
      }
      AST_data->nodeType = TK_MACRO ;
#ifdef VERBOSE
      cerr << "get token MACRO("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
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
      AST_data->nodeType = TK_LOCI_DIRECTIVE ;
#ifdef VERBOSE
      cerr << "get token DIRECTIVE("<< AST_data->text<< ")" << endl ;
#endif
      return AST_data ;
    }
    AST_data->nodeType = TK_LOCI_VARIABLE ;
    if(is.peek() == '*') {
      AST_data->nodeType = TK_LOCI_CONTAINER ;
      AST_data->text += is.get() ;
    }
    AST_data->text = "" ;
    if(!is.fail() && !is.eof() &&
       ((is.peek() >='a' && is.peek() <='z') ||
	(is.peek() >='A' && is.peek() <='Z') ||
	is.peek() == '_' || is.peek() == '@' ||
        is.peek() == '$' )
       ) {
      // extract variable name
      while(!is.fail() && !is.eof() &&
	    ((is.peek() >='a' && is.peek() <='z') ||
	     (is.peek() >='A' && is.peek() <='Z') || 
	     (is.peek() >='0' && is.peek() <='9') ||
	     is.peek() == '_' || is.peek() == '@' ||
             is.peek() == ':' ||
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
	} else if(is.peek() == ':') {
          is.get() ;
          if(is.peek() == ':') {
            AST_data->text += ':' ;
            AST_data->text += is.get() ;
          } else {
            is.putback(':') ;
            break ;
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
     
      
  AST_data->nodeType = TK_ERROR ;
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

  if(tok->nodeType == TK_LT) {
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
      case TK_SEMICOLON:
      case TK_LOGICAL_AND:
      case TK_LOGICAL_OR:
      case TK_OPENBRACKET:
      case TK_OPENBRACE:
      case TK_CLOSEBRACKET:
      case TK_CLOSEBRACE:
      case TK_MACRO:
      case TK_ERROR:
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
      case TK_OPENPAREN:
        parencount++ ;
        break ;
      case TK_CLOSEPAREN:
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
      case TK_LT:
        LT_LIST.push_back(toklist.size()-1 ) ;
        LT_DEPTH.push_back(parencount) ;
#ifdef VERBOSE
        cerr << "searching for nested template '<'" << endl ;
#endif
        break;
      case TK_GT:
        // found match to '<'
        if(parencount == LT_DEPTH.back()) {
#ifdef VERBOSE
          cerr << "found matching '>'" << endl ;
#endif
          toklist[LT_LIST.back()]->nodeType =
            TK_OPENTEMPLATE ;
          toklist[LT_LIST.back()]->text = "<" ; ;
          toklist.back()->nodeType =
            TK_CLOSETEMPLATE ;
          toklist.back()->text = ">" ;
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
      case TK_SHIFT_RIGHT:
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
              TK_OPENTEMPLATE ;
            toklist[LT_LIST[sz-2]]->nodeType =
              TK_OPENTEMPLATE ;
            toklist.back()->nodeType =
              TK_CLOSETEMPLATE ;
            toklist.back()->text = ">" ;
            CPTR<AST_Token> AST_data = new AST_Token() ;
            AST_data->text = ">" ; ;
            AST_data->lineno = toklist.back()->lineno ;
            AST_data->nodeType = TK_CLOSETEMPLATE ;
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
  case OP_SCOPE:
    return string("OP_SCOPE") ;
  case OP_AT:
    return string("OP_AT") ;
  case OP_DOT:
    return string("OP_DOT") ;
  case OP_ARROW:
    return string("OP_ARROW") ;
  case OP_TIMES:
    return string("OP_TIMES") ;
  case OP_DIVIDE:
    return string("OP_DIVIDE") ;
  case OP_MODULUS:
    return string("OP_MODULUS") ;
  case OP_PLUS:
    return string("OP_PLUS") ;
  case OP_MINUS:
    return string("OP_MINUS") ;
  case OP_SHIFT_RIGHT:
    return string("OP_SHIFT_RIGHT") ;
  case OP_SHIFT_LEFT:
    return string("OP_SHIFT_RIGHT") ;
  case OP_LT:
    return string("OP_LT") ;
  case OP_GT:
    return string("OP_GT") ;
  case OP_GE:
    return string("OP_GE") ;
  case OP_LE:
    return string("OP_LE") ;
  case OP_EQUAL:
    return string("OP_EQUAL") ;
  case OP_NOT_EQUAL:
    return string("OP_NOT_EQUAL") ;
  case OP_AND:
    return string("OP_AND") ;
  case OP_EXOR:
    return string("OP_EXOR") ;
  case OP_OR:
    return string("OP_OR") ;
  case OP_LOGICAL_AND:
    return string("OP_LOGICAL_AND") ;
  case OP_LOGICAL_OR:
    return string("OP_LOGICAL_OR") ;
  case OP_TERNARY:
    return string("OP_TERNARY") ;
  case OP_ASSIGN:
    return string("OP_ASSIGN") ;
  case OP_TIMES_ASSIGN:
    return string("OP_TIMES_ASSIGN") ;
  case OP_DIVIDE_ASSIGN:
    return string("OP_DIVIDE_ASSIGN") ;
  case OP_MODULUS_ASSIGN:
    return string("OP_MODULUS_ASSIGN") ;
  case OP_PLUS_ASSIGN:
    return string("OP_PLUS_ASSIGN") ;
  case OP_MINUS_ASSIGN:
    return string("OP_MINUS_ASSIGN") ;
  case OP_SHIFT_LEFT_ASSIGN:
    return string("OP_SHIFT_LEFT_ASSIGN") ;
  case OP_SHIFT_RIGHT_ASSIGN:
    return string("OP_SHIFT_RIGHT_ASSIGN") ;
  case OP_AND_ASSIGN:
    return string("OP_AND_ASSIGN") ;
  case OP_OR_ASSIGN:
    return string("OP_OR_ASSIGN") ;
  case OP_EXOR_ASSIGN:
    return string("OP_EXOR_ASSIGN") ;
  case OP_COMMA:
    return string("OP_COMMA") ;
  case OP_COLON:
    return string("OP_COLON") ;
  case OP_SEMICOLON:
    return string("OP_SEMICOLON") ;
  case OP_NIL:
    return string("OP_NIL") ;
  case OP_INCREMENT:
    return string("OP_INCREMENT") ;
  case OP_DECREMENT:
    return string("OP_DECREMENT") ;
  case OP_POSTINCREMENT:
    return string("OP_POSTINCREMENT") ;
  case OP_POSTDECREMENT:
    return string("OP_POSTDECREMENT") ;
  case OP_COMMENT:
    return string("OP_COMMENT") ;
  case OP_BRACEBLOCK:
    return string("OP_BRACEBLOCK") ;
  case OP_NAME:
    return string("OP_NAME") ;
  case OP_FUNC:
    return string("OP_FUNC") ;
  case OP_ARRAY:
    return string("OP_ARRAY") ;
  case OP_NAME_BRACE:
    return string("OP_NAME_BRACE") ;
  case OP_FUNC_BRACE:
    return string("OP_FUNC_BRACE") ;
  case OP_TEMPLATE:
    return string("OP_TEMPLATE") ;
  case OP_TEMPLATE_FUNC:
    return string("OP_TEMPLATE_FUNC") ;
  case OP_STRING:
    return string("OP_STRING") ;
  case OP_NUMBER:
    return string("OP_NUMBER") ;
  case OP_ERROR:
    return string("OP_ERROR") ;
    
  case OP_UNARY_PLUS:
    return string("OP_UNARY_PLUS") ;
  case OP_UNARY_MINUS:
    return string("OP_UNARY_MINUS") ;
  case OP_NOT:
    return string("OP_NOT") ;
  case OP_TILDE:
    return string("OP_TILDE") ;
  case OP_AMPERSAND:
    return string("OP_AMPERSAND") ;
  case OP_DOLLAR:
    return string("OP_DOLLAR") ;
  case OP_STAR:
    return string("OP_STAR") ;
  case OP_CAST:
    return string("OP_CAST") ;
  case OP_GROUP:
    return string("OP_GROUP") ;
  case OP_GROUP_ERROR:
    return string("OP_GROUP_ERROR") ;
  case OP_OPENPAREN:
    return string("OP_OPENPAREN") ;
  case OP_CLOSEPAREN:
    return string("OP_CLOSEPAREN") ;
  case OP_OPENBRACKET:
    return string("OP_OPENBRACKET") ;
  case OP_CLOSEBRACKET:
    return string("OP_CLOSEBRACKET") ;
  case OP_OPENBRACE:
    return string("OP_OPENBRACE") ;
  case OP_CLOSEBRACE:
    return string("OP_CLOSEBRACE") ;
  case OP_LOCI_DIRECTIVE:
    return string("OP_LOCI_DIRECTIVE") ;
  case OP_LOCI_VARIABLE:
    return string("OP_LOCI_VARIABLE") ;
  case OP_LOCI_CONTAINER:
    return string("OP_LOCI_CONTAINER") ;
  case OP_TERM:
    return string("OP_TERM") ;
  case OP_SPECIAL:
    return string("OP_SPECIAL") ;
  case ND_BLOCK:
    return string("ND_BLOCK") ;
  case ND_SYNTAXERR:
    return string("ND_SYTNAXERR") ;
  case ND_CTRL_IF:
    return string("ND_CTRL_IF") ;
  case ND_CTRL_FOR:
    return string("ND_CTRL_FOR") ;
  case ND_CTRL_WHILE:
    return string("ND_CTRL_WHILE") ;
  case ND_CTRL_DO:
    return string("ND_CTRL_DO") ;
  case ND_CTRL_SWITCH:
    return string("ND_CTRL_SWITCH") ;
  case ND_SIMPLE_STATEMENT:
    return string("ND_SIMPLE_STATEMENT") ;
  case ND_DECL:
    return string("ND_DECL") ;
  case ND_TYPE_SPEC:
    return string("ND_TYPE_SPEC") ;
  case ND_TERMINAL:
    return string("ND_TERMINAL") ;
  case TK_SENTINEL:
    return string("TK_SENTINEL") ;
  case TK_BRACEBLOCK:
    return string("TK_BRACEBLOCK") ;
  case TK_SCOPE:
    return string("TK_SCOPE") ;
  case TK_NAME:
    return string("TK_NAME") ;
  case TK_AT:
    return string("TK_AT") ;
  case TK_ARROW:
    return string("TK_ARROW") ;
  case TK_TIMES:
    return string("TK_TIMES") ;
  case TK_DIVIDE:
    return string("TK_DIVIDE") ;
  case TK_MODULUS:
    return string("TK_MODULUS") ;
  case TK_PLUS:
    return string("TK_PLUS") ;
  case TK_MINUS:
    return string("TK_MINUS") ;
  case TK_SHIFT_RIGHT:
    return string("TK_SHIFT_RIGHT") ;
  case TK_SHIFT_LEFT:
    return string("TK_SHIFT_LEFT") ;
  case TK_LT:
    return string("TK_LT") ;
  case TK_GT:
    return string("TK_GT") ;
  case TK_GE:
    return string("TK_GE") ;
  case TK_LE:
    return string("TK_LE") ;
  case TK_EQUAL:
    return string("TK_EQUAL") ;
  case TK_NOT_EQUAL:
    return string("TK_NOT_EQUAL") ;
  case TK_AND:
    return string("TK_AND") ;
  case TK_EXOR:
    return string("TK_EXOR") ;
  case TK_OR:
    return string("TK_OR") ;
  case TK_LOGICAL_AND:
    return string("TK_LOGICAL_AND") ;
  case TK_LOGICAL_OR:
    return string("TK_LOGICAL_OR") ;
  case TK_ASSIGN:
    return string("TK_ASSIGN") ;
  case TK_TIMES_ASSIGN:
    return string("TK_TIMES_ASSIGN") ;
  case TK_DIVIDE_ASSIGN:
    return string("TK_DIVIDE_ASSIGN") ;
  case TK_MODULUS_ASSIGN:
    return string("TK_MODULUS_ASSIGN") ;
  case TK_PLUS_ASSIGN:
    return string("TK_PLUS_ASSIGN") ;
  case TK_MINUS_ASSIGN:
    return string("TK_MINUS_ASSIGN") ;
  case TK_SHIFT_LEFT_ASSIGN:
    return string("TK_SHIFT_LEFT_ASSIGN") ;
  case TK_SHIFT_RIGHT_ASSIGN:
    return string("TK_SHIFT_RIGHT_ASSIGN") ;
  case TK_AND_ASSIGN:
    return string("TK_AND_ASSIGN") ;
  case TK_OR_ASSIGN:
    return string("TK_OR_ASSIGN") ;
  case TK_EXOR_ASSIGN:
    return string("TK_EXOR_ASSIGN") ;
  case TK_COMMA:
    return string("TK_COMMA") ;
  case TK_DOT:
    return string("TK_DOT") ;
  case TK_COLON:
    return string("TK_COLON") ;
  case TK_SEMICOLON:
    return string("TK_SEMICOLON") ;
  case TK_NIL:
    return string("TK_NIL") ;
  case TK_INCREMENT:
    return string("TK_INCREMENT") ;
  case TK_DECREMENT:
    return string("TK_DECREMENT") ;
  case TK_COMMENT:
    return string("TK_COMMENT") ;
  case TK_MACRO:
    return string("TK_MACRO") ;
  case TK_STRING:
    return string("TK_STRING") ;
  case TK_NUMBER:
    return string("TK_NUMBER") ;
  case TK_ERROR:
    return string("TK_ERROR") ;
  case TK_UNARY_PLUS:
    return string("TK_UNARY_PLUS") ;
  case TK_UNARY_MINUS:
    return string("TK_UNARY_MINUS") ;
  case TK_NOT:
    return string("TK_NOT") ;
  case TK_TILDE:
    return string("TK_TILDE") ;
  case TK_QUESTION:
    return string("TK_QUESTION") ;
  case TK_AMPERSAND:
    return string("TK_AMPERSAND") ;
  case TK_STAR:
    return string("TK_STAR") ;
  case TK_OPENPAREN:
    return string("TK_OPENPAREN") ;
  case TK_CLOSEPAREN:
    return string("TK_CLOSEPAREN") ;
  case TK_OPENBRACKET:
    return string("TK_OPENBRACKET") ;
  case TK_CLOSEBRACKET:
    return string("TK_CLOSEBRACKET") ;
  case TK_OPENBRACE:
    return string("TK_CLOSEBRACE") ;
  case TK_OPENTEMPLATE:
    return string("TK_CLOSETEMPLATE") ;
  case TK_LOCI_DIRECTIVE:
    return string("TK_LOCI_DIRECTIVE") ;
  case TK_LOCI_VARIABLE:
    return string("TK_LOCI_VARIABLE") ;
  case TK_LOCI_CONTAINER:
    return string("TK_LOCI_CONTAINER") ;

  case TK_ALIGNAS:
    return string("TK_ALIGNAS") ;

  case TK_ALIGNOF:
    return string("TK_ALIGNOF") ;

  case TK_ASM: 
    return string("TK_ASM") ; 

  case TK_BOOL:
    return string("TK_BOOL") ;

  case TK_FALSE:
    return string("TK_FALSE") ;

  case TK_TRUE:
    return string("TK_TRUE") ;

  case TK_CHAR:
    return string("TK_CHAR") ;

  case TK_INT:
    return string("TK_INT") ;

  case TK_LONG:
    return string("TK_LONG") ;
      
  case TK_SHORT:
    return string("TK_SHORT") ;

  case TK_SIGNED:
    return string("TK_SIGNED") ;

  case TK_UNSIGNED:
    return string("TK_UNSIGNED") ;

  case TK_DOUBLE:
    return string("TK_DOUBLE") ;

  case TK_FLOAT:
    return string("TK_FLOAT") ;

  case TK_ENUM:
    return string("TK_ENUM") ;

  case TK_MUTABLE:
    return string("TK_MUTABLE") ;

  case TK_CONST:
    return string("TK_CONST") ;

  case TK_STATIC:
    return string("TK_STATIC") ;

  case TK_VOLATILE:
    return string("TK_VOLATILE") ;

  case TK_AUTO:
    return string("TK_AUTO") ;

  case TK_REGISTER:
    return string("TK_REGISTER") ;

  case TK_EXPORT:
    return string("TK_EXPORT") ;

  case TK_EXTERN:
    return string("TK_EXTERN") ;

  case TK_INLINE:
    return string("TK_INLINE") ;

  case TK_NAMESPACE:
    return string("TK_NAMESPACE") ;

  case TK_EXPLICIT:
    return string("TK_EXPLICIT") ;

  case TK_DYNAMIC_CAST:
    return string("TK_DYNAMIC_CAST") ;

  case TK_STATIC_CAST:
    return string("TK_STATIC_CAST") ;

  case TK_REINTERPRET_CAST:
    return string("TK_REINTERPRET_CAST") ;

  case TK_OPERATOR:
    return string("TK_OPERATOR") ;

  case TK_PROTECTED:
    return string("TK_PROTECTED") ;

  case TK_NOEXCEPT:
    return string("TK_NOEXCEPT") ;

  case TK_NULLPTR:
    return string("TK_NULLPTR") ;

  case TK_RETURN:
    return string("TK_RETURN") ;

  case TK_SIZEOF:
    return string("TK_SIZEOF") ;

  case TK_THIS:
    return string("TK_THIS") ;

  case TK_TYPEID:
    return string("TK_TYPEID") ;

  case TK_SWITCH:
    return string("TK_SWITCH") ;

  case TK_CASE:
    return string("TK_CASE") ;

  case TK_BREAK:
    return string("TK_BREAK") ;

  case TK_DEFAULT:
    return string("TK_DEFAULT") ;

  case TK_FOR:
    return string("TK_FOR") ;

  case TK_DO:
    return string("TK_DO") ;

  case TK_WHILE:
    return string("TK_WHILE") ;

  case TK_CONTINUE:
    return string("TK_CONTINUE") ;

  case TK_CLASS:
    return string("TK_CLASS") ;

  case TK_STRUCT:
    return string("TK_STRUCT") ;

  case TK_PUBLIC:
    return string("TK_PUBLIC") ;

  case TK_PRIVATE:
    return string("TK_PRIVATE") ;

  case TK_FRIEND:
    return string("TK_FRIEND") ;

  case TK_UNION:
    return string("TK_UNION") ;

  case TK_TYPENAME:
    return string("TK_TYPENAME") ;

  case TK_TEMPLATE:
    return string("TK_TEMPLATE") ;

  case TK_TYPEDEF:
    return string("TK_TYPEDEF") ;

  case TK_VIRTUAL:
    return string("TK_VIRTUAL") ;

  case TK_VOID:
    return string("TK_VOID") ;

  case TK_TRY:
    return string("TK_TRY") ;

  case TK_CATCH:
    return string("TK_CATCH") ;

  case TK_THROW:
    return string("TK_THROW") ;

  case TK_IF:
    return string("TK_IF") ;

  case TK_ELSE:
    return string("TK_ELSE") ;

  case TK_GOTO:
    return string("TK_GOTO") ;

  case TK_NEW:
    return string("TK_NEW") ;

  case TK_DELETE:
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

