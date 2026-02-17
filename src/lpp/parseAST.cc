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

/// Convert operator code to corresponding string for printing
string OPtoString(AST_type::elementType val) {
  switch(val) {
  case AST_type::OP_SCOPE:
    return string("::") ;
  case AST_type::OP_AT:
    return string("@") ;
  case AST_type::OP_ARROW:
    return string("->") ;
  case AST_type::OP_TIMES:
    return string("*") ;
  case AST_type::OP_DIVIDE:
    return string("/") ;
  case AST_type::OP_MODULUS:
    return string("%") ;
  case AST_type::OP_PLUS:
    return string("+") ;
  case AST_type::OP_MINUS:
    return string("-") ;
  case AST_type::OP_SHIFT_RIGHT:
    return string(">>") ;
  case AST_type::OP_SHIFT_LEFT:
    return string("<<") ;
  case AST_type::OP_LT:
    return string("<") ;
  case AST_type::OP_GT:
    return string(">") ;
  case AST_type::OP_GE:
    return string(">=") ;
  case AST_type::OP_LE:
    return string("<=") ;
  case AST_type::OP_EQUAL:
    return string("==") ;
  case AST_type::OP_NOT_EQUAL:
    return string("!=") ;
  case AST_type::OP_AND:
    return string("&") ;
  case AST_type::OP_EXOR:
    return string("^") ;
  case AST_type::OP_OR:
    return string("|") ;
  case AST_type::OP_LOGICAL_AND:
    return string("&&") ;
  case AST_type::OP_LOGICAL_OR:
    return string("||") ;
  case AST_type::OP_ASSIGN:
    return string("=") ;
  case AST_type::OP_TIMES_ASSIGN:
    return string("*=") ;
  case AST_type::OP_DIVIDE_ASSIGN:
    return string("/=") ;
  case AST_type::OP_MODULUS_ASSIGN:
    return string("%=") ;
  case AST_type::OP_PLUS_ASSIGN:
    return string("+=") ;
  case AST_type::OP_MINUS_ASSIGN:
    return string("-=") ;
  case AST_type::OP_SHIFT_LEFT_ASSIGN:
    return string("<<=") ;
  case AST_type::OP_SHIFT_RIGHT_ASSIGN:
    return string(">>=") ;
  case AST_type::OP_AND_ASSIGN:
    return string("&=") ;
  case AST_type::OP_OR_ASSIGN:
    return string("|=") ;
  case AST_type::OP_EXOR_ASSIGN:
    return string("^=") ;
  case AST_type::OP_COMMA:
    return string(",") ;
  case AST_type::OP_DOT:
    return string(".") ;
  case AST_type::OP_COLON:
    return string(":") ;
  case AST_type::OP_SEMICOLON:
    return string(";") ;
  case AST_type::OP_INCREMENT:
    return string(" ++") ;
  case AST_type::OP_DECREMENT:
    return string(" --") ;
  case AST_type::OP_POSTINCREMENT:
    return string("++ ") ;
  case AST_type::OP_POSTDECREMENT:
    return string("-- ") ;
  case AST_type::OP_UNARY_PLUS:
    return string("+") ;
  case AST_type::OP_UNARY_MINUS:
    return string("-") ;
  case AST_type::OP_NOT:
    return string("!") ;
  case AST_type::OP_TILDE:
    return string("~") ;
  case AST_type::OP_AMPERSAND:
    return string("&") ;
  case AST_type::OP_TERNARY:
    return string("?") ;
  case AST_type::OP_DOLLAR:
    return string("$") ;
  case AST_type::OP_STAR:
    return string("*") ;
  default:
    return string("/*error*/") ;
  }
  return string("/*error*/") ;
}

bool isTypeDecl(CPTR<AST_Token> p, varmap &typemap) {
  switch(p->nodeType) {
  case AST_type::TK_CHAR:
  case AST_type::TK_FLOAT:
  case AST_type::TK_DOUBLE:
  case AST_type::TK_INT:
  case AST_type::TK_BOOL:
  case AST_type::TK_SHORT:
  case AST_type::TK_LONG:
  case AST_type::TK_SIGNED:
  case AST_type::TK_UNSIGNED:
  case AST_type::TK_CONST:
  case AST_type::TK_AUTO:
    return true ;
  default:
    return false ;
  }
}

AST_type::ASTP parseTypeSpecifier(std::istream &is, int &linecount,
                                  const string &fileName,
                                  varmap &typemap) {
  AST_type::ASTList type_spec ;
  bool istype = true ;
  bool builtin_type = false ;
  bool defined_type = false ;

  CPTR<AST_Token> token = getToken(is,linecount) ;
  do {
#ifdef VERBOSE
    cerr << "in parseTypeSpecifier, token = " << token->text << endl ;
#endif
    switch(token->nodeType) {
    case AST_type::TK_CONST:
      type_spec.push_back(AST_type::ASTP(token)) ;
#ifdef VERBOSE
      cerr << "in parseTypeSpecifier const variable" << endl ;
#endif
      break ;
    case AST_type::TK_EXTERN:
      type_spec.push_back(AST_type::ASTP(token)) ;
#ifdef VERBOSE
      cerr << "in parseTypeSpecifier extern variable" << endl ;
#endif
      break ;
    case AST_type::TK_CHAR:
    case AST_type::TK_FLOAT:
    case AST_type::TK_DOUBLE:
    case AST_type::TK_INT:
    case AST_type::TK_BOOL:
    case AST_type::TK_SHORT:
    case AST_type::TK_LONG:
    case AST_type::TK_SIGNED:
    case AST_type::TK_UNSIGNED:
    case AST_type::TK_AUTO:
#ifdef VERBOSE
      cerr << "in parseTypeSpecifier, builtin type" << endl ;
#endif
      builtin_type = true ;
      if(defined_type) { // This is a syntax error
        AST_type::ASTP p =
          AST_type::ASTP(new AST_syntaxError("Mixed concrete and defined types",
                                             token->lineno,fileName)) ;
        type_spec.push_back(p) ;
        istype = false ;
      }
      type_spec.push_back(AST_type::ASTP(token)) ;
      break ;
    case AST_type::TK_SCOPE:
    case AST_type::TK_NAME:
#ifdef VERBOSE
      cerr << "in parseTypeSpecifier, named type" << endl ;
#endif
      if(builtin_type || defined_type ) {
        istype = false ;
      } else {
        pushToken(token) ;
        
        AST_type::ASTP typedec =
          parseIdentifier(is,linecount,fileName,typemap) ;
        if(typedec->nodeType == AST_type::OP_FUNC) {
          // not a proper type name
          AST_type::ASTP p =
            AST_type::ASTP(new AST_syntaxError("Mixed concrete and defined types",
                                               token->lineno,fileName)) ;
          type_spec.push_back(p) ;
        } else {

          // right now only consider unscoped non-templated type names
          defined_type = true ;
          type_spec.push_back(typedec) ;
        }
      }
      break ;
    default:
#ifdef VERBOSE
      cerr << "in parseTypeSpecifier, finished type parsing" << endl ;
#endif
      istype = false ;
      break ;
    }
    if(istype)
      token = getToken(is,linecount) ;
  } while(istype) ;

  pushToken(token) ;
  CPTR<AST_typeSpec> p = new AST_typeSpec ;
  p->nodeType = AST_type::ND_TYPE_SPEC ;
  p->type_spec = type_spec ;
  return AST_type::ASTP(p) ;
}

static int parseIDctr = 0 ;
AST_type::AST_type() {
  nodeType = OP_ERROR ;
  id = parseIDctr++ ;
}

/// Acceptor method AST node, passes node to visitor object
void AST_Token::accept(AST_visitor &v) {  v.visit(*this) ; }

AST_type::ASTP AST_Token::clone() const {
  CPTR<AST_Token> tok = new AST_Token ;
  tok->nodeType = nodeType ;
  tok->text = text ;
  tok->lineno = lineno ;
  return ASTP(tok) ;
}

/// Acceptor method AST node, passes node to visitor object
void AST_syntaxError::accept(AST_visitor &v) {  v.visit(*this) ; }

AST_type::ASTP AST_syntaxError::clone() const {
  CPTR<AST_syntaxError> p = new AST_syntaxError(error,lineno,fileName) ;
  p->nodeType = nodeType ;
  return ASTP(p) ;
}

/// Acceptor method AST node, passes node to visitor object
void AST_SimpleStatement::accept(AST_visitor &v) {  v.visit(*this) ; }

AST_type::ASTP AST_SimpleStatement::clone() const {
  CPTR<AST_SimpleStatement> p = new AST_SimpleStatement(exp->clone(),
                                                        Terminal->clone()) ;
  p->nodeType = nodeType ;
  return ASTP(p) ;
}

/// Acceptor method AST node, passes node to visitor object
void AST_Block::accept(AST_visitor &v) {  v.visit(*this) ; }

inline void cloneList(AST_type::ASTList &l, const AST_type::ASTList &in) {
  AST_type::ASTList tmp(in.size()) ;
  for(size_t i=0;i<in.size();++i)
    tmp[i] = in[i]->clone() ;
  l.swap(tmp) ;
}

AST_type::ASTP AST_Block::clone() const {
  CPTR<AST_Block> p = new AST_Block ;
  p->nodeType = nodeType ;
  cloneList(p->elements,elements) ;
  p->identifiers = identifiers ;
  return ASTP(p) ;
}

/// Acceptor method AST node, passes node to visitor object
void AST_typeSpec::accept(AST_visitor &v) {  v.visit(*this) ; }


AST_type::ASTP AST_typeSpec::clone() const {
  CPTR<AST_typeSpec> p = new AST_typeSpec ;
  p->nodeType = nodeType ;
  cloneList(p->type_spec,type_spec) ;
  return ASTP(p) ;
}

/// Acceptor method AST node, passes node to visitor object
void AST_declaration::accept(AST_visitor &v) {  v.visit(*this) ; }


AST_type::ASTP AST_declaration::clone() const {
  CPTR<AST_declaration> p = new AST_declaration ;
  p->nodeType = nodeType ;
  cloneList(p->type_decl,type_decl) ;
  cloneList(p->decls,decls) ;
  p->decl_varnames = decl_varnames ;
  return ASTP(p) ;
}
/// Acceptor method AST node, passes node to visitor object
void AST_exprOper::accept(AST_visitor &v) {  v.visit(*this) ; }

AST_type::ASTP AST_exprOper::clone() const {
  CPTR<AST_exprOper> p = new AST_exprOper ;
  p->nodeType = nodeType ;
  cloneList(p->terms,terms) ;
  return ASTP(p) ;
}

/// Acceptor method AST node, passes node to visitor object
void AST_controlStatement::accept(AST_visitor &v) {  v.visit(*this) ; }

AST_type::ASTP AST_controlStatement::clone() const {
  CPTR<AST_controlStatement> p = new AST_controlStatement ;
  p->nodeType = nodeType ;
  p->controlType = controlType->clone() ;
  cloneList(p->parts,parts) ;
  p->identifiers = identifiers ;
  return ASTP(p) ;
}

/// Check if token is a terminal symbol
inline bool isTerm(AST_type::elementType e) {
  switch(e) {
  case AST_type::TK_STRING:
  case AST_type::TK_NAME:
  case AST_type::TK_NUMBER:
  case AST_type::TK_TRUE:
  case AST_type::TK_FALSE:
  case AST_type::TK_LOCI_VARIABLE:
  case AST_type::TK_LOCI_CONTAINER:
    return true ;
  default:
    return false ;
  }
  return false ;
}

/// Check if token is a terminal symbol
inline bool isTerm(const CPTR<AST_Token> &p) {
  return isTerm(p->nodeType) ;
}

/// Parse terminal character
AST_type::ASTP parseTerm(std::istream &is, int &linecount) {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseTerm, token = " << token->text << endl ;
#endif
  if(isTerm(token))
    return AST_type::ASTP(token) ;
  else {
    pushToken(token) ;
    AST_type::ASTP(new AST_syntaxError("expecting token",token->lineno,"")) ;    
    return AST_type::ASTP(0) ;
  }
}




/// Parse operator token, convert from TK_* node type to OP_* nodetype
AST_type::ASTP parseOperator(std::istream &is, int &linecount) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseOperator, token = " << openToken->text << endl ;
#endif
  switch(openToken->nodeType) {
  case AST_type::TK_SCOPE:
    openToken->nodeType = AST_type::OP_SCOPE ;
    return AST_type::ASTP(openToken) ;
    //-----------------------------------
    // For using @ to separate namespaces
  case AST_type::TK_AT:
    openToken->nodeType = AST_type::OP_AT ;
    return AST_type::ASTP(openToken) ;
    //-----------------------------------
    // Traditional C operators
  case AST_type::TK_ARROW:
    openToken->nodeType = AST_type::OP_ARROW ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_TIMES:
    openToken->nodeType = AST_type::OP_TIMES ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_DIVIDE:
    openToken->nodeType = AST_type::OP_DIVIDE ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_MODULUS:
    openToken->nodeType = AST_type::OP_MODULUS ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_PLUS:
    openToken->nodeType = AST_type::OP_PLUS ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_MINUS:
    openToken->nodeType = AST_type::OP_MINUS ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_SHIFT_RIGHT:
    openToken->nodeType = AST_type::OP_SHIFT_RIGHT ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_SHIFT_LEFT:
    openToken->nodeType = AST_type::OP_SHIFT_LEFT ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_LT:
    openToken->nodeType = AST_type::OP_LT ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_GT:
    openToken->nodeType = AST_type::OP_GT ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_GE:
    openToken->nodeType = AST_type::OP_GE ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_LE:
    openToken->nodeType = AST_type::OP_LE ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_EQUAL:
    openToken->nodeType = AST_type::OP_EQUAL ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_NOT_EQUAL:
    openToken->nodeType = AST_type::OP_NOT_EQUAL ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_AND:
    openToken->nodeType = AST_type::OP_AND ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_EXOR:
    openToken->nodeType = AST_type::OP_EXOR ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_OR:
    openToken->nodeType = AST_type::OP_OR ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_LOGICAL_AND:
    openToken->nodeType = AST_type::OP_LOGICAL_AND;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_LOGICAL_OR:
    openToken->nodeType = AST_type::OP_LOGICAL_OR ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_ASSIGN:
    openToken->nodeType = AST_type::OP_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_TIMES_ASSIGN:
    openToken->nodeType = AST_type::OP_TIMES_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_DIVIDE_ASSIGN:
    openToken->nodeType = AST_type::OP_DIVIDE_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_MODULUS_ASSIGN:
    openToken->nodeType = AST_type::OP_MODULUS_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_PLUS_ASSIGN:
    openToken->nodeType = AST_type::OP_PLUS_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_MINUS_ASSIGN:
    openToken->nodeType = AST_type::OP_MINUS_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_SHIFT_LEFT_ASSIGN:
    openToken->nodeType = AST_type::OP_SHIFT_LEFT_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_SHIFT_RIGHT_ASSIGN:
    openToken->nodeType = AST_type::OP_SHIFT_RIGHT_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_AND_ASSIGN:
    openToken->nodeType = AST_type::OP_AND_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_OR_ASSIGN:
    openToken->nodeType = AST_type::OP_OR_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_EXOR_ASSIGN:
    openToken->nodeType = AST_type::OP_EXOR_ASSIGN ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_COMMA:
    openToken->nodeType = AST_type::OP_COMMA ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_QUESTION:
    openToken->nodeType = AST_type::OP_TERNARY ;
    return AST_type::ASTP(openToken) ;
  case AST_type::TK_DOT:
    openToken->nodeType = AST_type::OP_DOT ;
    return AST_type::ASTP(openToken) ;
  default:
    pushToken(openToken) ;
    return AST_type::ASTP(0) ;
  }
}

/// Check for input token that is a unary operator
inline bool checkUnaryToken(AST_type::elementType e) {
  switch(e) {
  case AST_type::TK_PLUS:
  case AST_type::TK_MINUS:
  case AST_type::TK_NOT:
  case AST_type::TK_AND:
  case AST_type::TK_TIMES:
  case AST_type::TK_INCREMENT:
  case AST_type::TK_DECREMENT:
    return true ;
  default:
    return false ;
  }
  return false ;
}

inline bool checkUnaryToken(const CPTR<AST_Token> &x) {
  return checkUnaryToken(x->nodeType) ;
}

/// Parse Unary operator, convert to TK_* node type to OP_* node type.
AST_type::elementType unaryOperator(AST_type::elementType e) {
  switch(e) {
  case AST_type::TK_PLUS:
    return AST_type::OP_UNARY_PLUS ;
  case AST_type::TK_MINUS:
    return AST_type::OP_UNARY_MINUS ;
  case AST_type::TK_NOT:
    return AST_type::OP_NOT ;
  case AST_type::TK_AND:
    return AST_type::OP_AMPERSAND ;
  case AST_type::TK_TIMES:
    return AST_type::OP_STAR ;
  case AST_type::TK_INCREMENT:
    return AST_type::OP_INCREMENT ;
  case AST_type::TK_DECREMENT:
    return AST_type::OP_DECREMENT ;
  default:
    return AST_type::OP_ERROR ;
  }
}

/// Check for postfix token '++' or '--' 
bool checkPostFixToken(AST_type::elementType e) {
  return (e == AST_type::TK_INCREMENT ||
	  e == AST_type::TK_DECREMENT) ;
}

/// Pares postfix operator convert node type from TK_* to appropriate OP_* 
AST_type::elementType postFixOperator(AST_type::elementType e) {
  switch(e) {
  case AST_type::TK_INCREMENT:
    return AST_type::OP_POSTINCREMENT ;
  case AST_type::TK_DECREMENT:
    return AST_type::OP_POSTDECREMENT ;
  default:
    return AST_type::OP_ERROR ;
  }
  return AST_type::OP_ERROR ;
}

/// Apply postfix operator to previous expression passed in expr
AST_type::ASTP applyPostFixOperator(AST_type::ASTP expr,
				    std::istream &is, int &linecount,
				    string fileName,
				    varmap &typemap) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in applyPostFixOperator, token = " << openToken->text << endl ;
#endif
  if(checkPostFixToken(openToken->nodeType)) {
    CPTR<AST_exprOper> post = new AST_exprOper ;
    post->nodeType = postFixOperator(openToken->nodeType) ;
    post->terms.push_back(expr) ;
    openToken = getToken(is,linecount) ;
    pushToken(openToken) ;
    if(checkPostFixToken(openToken->nodeType) ||
       ASTEqual(openToken,AST_type::TK_OPENBRACKET))
      return applyPostFixOperator(AST_type::ASTP(post),is,linecount,
				  fileName,typemap) ;
    else
      return AST_type::ASTP(post) ;
  }
  if(ASTEqual(openToken,AST_type::TK_OPENBRACKET)) {
    AST_type::ASTP index = parseExpression(is,linecount,fileName,typemap) ;
    openToken = getToken(is,linecount) ;
    if(openToken->nodeType != AST_type::TK_CLOSEBRACKET) {
      pushToken(openToken) ;
      return AST_type::ASTP(new AST_syntaxError("expecting ']'",openToken->lineno,fileName)) ;
    }
    CPTR<AST_exprOper> array = new AST_exprOper ;
    array->nodeType = AST_type::OP_ARRAY ;
    array->terms.push_back(expr) ;
    array->terms.push_back(index) ;
    
    openToken = getToken(is,linecount) ;
    pushToken(openToken) ;
    if(checkPostFixToken(openToken->nodeType) ||
       ASTEqual(openToken,AST_type::TK_OPENBRACKET))
      return applyPostFixOperator(AST_type::ASTP(array),is,linecount,fileName,
				  typemap) ;
    else
      return AST_type::ASTP(array) ;
  }
  pushToken(openToken) ;
    
  return expr ;
}

AST_type::ASTP parseScopedObject(std::istream &is, int &linecount,
                                 const string &fileName,
                                 varmap &typemap) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "parseScopedObject, parsing named object with scope, "
       << OPtoName(openToken->nodeType) ;
    if(ASTEqual(openToken,AST_type::TK_NAME))
      cerr << " name=" << openToken->text ;
    cerr << endl ;
#endif
  AST_type::ASTList terms ;
  if(ASTEqual(openToken,AST_type::TK_SCOPE)) {
    CPTR<AST_Token> tmp = new AST_Token ;
    tmp->text = "" ;
    tmp->lineno = openToken->lineno ;
    tmp->nodeType = AST_type::TK_NIL ;
    terms.push_back(AST_type::ASTP(tmp)) ;
  } else {
    terms.push_back(AST_type::ASTP(openToken)) ;
    openToken = getToken(is,linecount) ;
#ifdef VERBOSE
    cerr << "searching token = " << openToken->text << endl ;
#endif
    if(!ASTEqual(openToken,AST_type::TK_SCOPE)) {
#ifdef VERBOSE
      cerr << "parseScopeObject, did not find scope operator, pushing token" << endl ;
#endif
        pushToken(openToken) ;
      }
  }
  AST_type::ASTP objname = 0 ;

  // If no scope operator then just process the name
  if(ASTEqual(openToken,AST_type::TK_SCOPE)) {
    // Process scope operators to get real name
#ifdef VERBOSE
    cerr << "parseScopedObject, found scope operator" << endl ;
#endif
    openToken = getToken(is,linecount) ;
    while(ASTEqual(openToken,AST_type::TK_NAME)) {
#ifdef VERBOSE
      cerr << "parseScopedObject, found token = " << openToken->text
           << endl ;
#endif
      terms.push_back(AST_type::ASTP(openToken)) ;
      openToken = getToken(is,linecount) ;
      if(ASTEqual(openToken,AST_type::TK_SCOPE)) {
        openToken = getToken(is,linecount) ;
#ifdef VERBOSE
        cerr << "found scope operator, token after = " << openToken->text
             << endl ;
#endif
        if(!ASTEqual(openToken,AST_type::TK_NAME)) {
          CPTR<AST_exprOper> tmp = new AST_exprOper ;
          tmp->nodeType = AST_type::OP_ERROR ;
          terms.push_back(AST_type::ASTP(tmp)) ;
        }
      }
    }
    pushToken(openToken) ;
    CPTR<AST_exprOper> tmp = new AST_exprOper ;
    tmp->nodeType = AST_type::OP_SCOPE ;
    tmp->terms = terms ;
    objname = AST_type::ASTP(tmp) ;
  } else {
    objname = AST_type::ASTP(terms.back()) ;
  }
  return objname ;
}

AST_type::ASTP parseTemplateArguments(AST_type::ASTP objname,
                                      std::istream &is, int &linecount,
                                      const string &fileName,
                                      varmap &typemap) {
  // Template Arguments
#ifdef VERBOSE
  cerr << "parseTemplateArguments: found term template open '<'" 
           << endl ;
#endif
  CPTR<AST_Token> token = getToken(is,linecount) ;
  AST_type::ASTP arg = AST_type::ASTP(token) ;
  if(!ASTEqual(token,AST_type::TK_NUMBER)) {
    pushToken(token) ;
    AST_type::ASTP arg = parseTypeSpecifier(is,linecount,fileName,typemap) ;
  }
  
#ifdef VERBOSE
  cerr << "parseTemplateArguments: parsed type specification " 
       << endl ;
#endif    
  CPTR<AST_Token> closeToken = getToken(is,linecount) ;
  if(ASTEqual(closeToken,AST_type::TK_CLOSETEMPLATE)) {
#ifdef VERBOSE
  cerr << "parseTemplateArguments: only single argument " 
       << endl ;
#endif    
    // only a single type in template parameters
    CPTR<AST_exprOper> func = new AST_exprOper ;
    func->nodeType = AST_type::OP_TEMPLATE ;
    func->terms.push_back(objname) ;
    func->terms.push_back(arg) ;
    return AST_type::ASTP(func) ;
  }
  if(ASTEqual(closeToken,AST_type::TK_COMMA)) {

#ifdef VERBOSE
  cerr << "parseTemplateArguments: parsing comma separated list of types " 
       << endl ;
#endif    
    
    CPTR<AST_exprOper> tlist = new AST_exprOper ;
    tlist->nodeType = AST_type::OP_COMMA ;    tlist->terms.push_back(arg) ;
    while(ASTEqual(closeToken,AST_type::TK_COMMA)) {
#ifdef VERBOSE
      cerr << "parseTemplateArgument, got " << closeToken->text<< endl ;
#endif
      CPTR<AST_Token> token = getToken(is,linecount) ;
      AST_type::ASTP arg = AST_type::ASTP(token) ;
      if(!ASTEqual(token,AST_type::TK_NUMBER)) {
        pushToken(token) ;
        AST_type::ASTP arg = parseTypeSpecifier(is,linecount,fileName,typemap) ;
      }
      
      tlist->terms.push_back(arg) ;
      closeToken = getToken(is,linecount) ;
    }
    if(!ASTEqual(closeToken,AST_type::TK_CLOSETEMPLATE)) {
      tlist->terms.push_back(AST_type::ASTP(new AST_syntaxError("expecting '>'",closeToken->lineno,fileName))) ;
    }
    CPTR<AST_exprOper> func = new AST_exprOper ;
    func->nodeType = AST_type::OP_TEMPLATE ;
    func->terms.push_back(objname) ;
    func->terms.push_back(AST_type::ASTP(tlist)) ;
    return AST_type::ASTP(func) ;
  }
  
  AST_type::ASTP p =
    AST_type::ASTP(new AST_syntaxError("Syntax error parsing template arguments",
                                       closeToken->lineno,fileName)) ;
  return p ;
}

AST_type::ASTP parseFunctionArguments(AST_type::ASTP objname,
                                      std::istream &is, int &linecount,
                                      const string &fileName,
                                      varmap &typemap) {
#ifdef VERBOSE
  cerr << "parseFunctionArguments: found term template open '('" 
       << endl ;
#endif    
  AST_type::ASTP args = parseExpression(is,linecount,fileName,typemap) ;
  CPTR<AST_Token> closeToken = getToken(is,linecount) ;
  CPTR<AST_exprOper> func = new AST_exprOper ;
  func->nodeType = AST_type::OP_FUNC ;
  func->terms.push_back(AST_type::ASTP(objname)) ;
  if(args != 0)
    func->terms.push_back(args) ;
  if(!ASTEqual(closeToken,AST_type::TK_CLOSEPAREN)) {
    CPTR<AST_Token> err = new AST_Token ;
    *err = *closeToken ;
    err->nodeType = AST_type::OP_ERROR ;
    func->terms.push_back(AST_type::ASTP(err)) ;
    pushToken(closeToken) ;
  }
  return AST_type::ASTP(func) ;
}
/// Parse identifier or function call
AST_type::ASTP parseIdentifier(std::istream &is, int &linecount,
				      const string &fileName,
				      varmap &typemap) {
#ifdef VERBOSE
  cerr << "in parseIdentifier" << endl ;
#endif 

  AST_type::ASTP objname = parseScopedObject(is,linecount,fileName,typemap) ;

  CPTR<AST_Token> openToken = getToken(is,linecount) ;
  if(ASTEqual(openToken,AST_type::TK_OPENTEMPLATE)) {
    objname = parseTemplateArguments(objname,is,linecount,fileName,typemap) ;
    openToken = getToken(is,linecount) ;
    if(ASTEqual(openToken,AST_type::TK_SCOPE)) {
      CPTR<AST_exprOper> tmp = new AST_exprOper ;
      tmp->nodeType = AST_type::OP_SCOPE ;
      tmp->terms.push_back(objname) ;
      openToken = getToken(is,linecount) ;
      tmp->terms.push_back(AST_type::ASTP(openToken)) ;
      if(ASTEqual(openToken,AST_type::TK_NAME)) {
        return AST_type::ASTP(tmp) ;
      } else {
        AST_type::ASTP p =
          AST_type::ASTP(new AST_syntaxError("error parsing scope after template type",openToken->lineno,fileName)) ;
        tmp->terms.push_back(p) ;
        return AST_type::ASTP(tmp) ;
      }
    } else
      pushToken(openToken) ;
  } else
    pushToken(openToken) ;
    
  openToken = getToken(is,linecount) ;
  if(ASTEqual(openToken,AST_type::TK_OPENPAREN)) {
    objname = parseFunctionArguments(objname,is,linecount,fileName,typemap) ;
  } else
    pushToken(openToken) ;
  return objname ;
}

/// Parse the part of an expression that is not a binary or ternary operator
AST_type::ASTP parseExpressionPartial(std::istream &is, int &linecount,
				      const string &fileName,
				      varmap &typemap) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseExpressionPartial, token = " << openToken->text << endl ;
#endif

  // Check for a type cast
  if(ASTEqual(openToken,AST_type::TK_DOUBLE) ||
     ASTEqual(openToken,AST_type::TK_FLOAT) ||
     ASTEqual(openToken,AST_type::TK_INT) ||
     ASTEqual(openToken,AST_type::TK_SHORT) ||
     ASTEqual(openToken,AST_type::TK_CHAR) ||
     ASTEqual(openToken,AST_type::TK_LONG)) {
    CPTR<AST_Token> nextToken = getToken(is,linecount) ;
    pushToken(nextToken) ;
    if(ASTEqual(nextToken,AST_type::TK_OPENPAREN)) {
      CPTR<AST_declaration> AST_data = new AST_declaration ;
      AST_data->type_decl.push_back(AST_type::ASTP(openToken)) ;
      AST_data->type_decl.push_back(parseExpressionPartial(is,linecount,
                                                           fileName,
                                                           typemap)) ;
      return AST_type::ASTP(AST_data) ;
    } 
  }
  // Check for a parenthesized group, if found parse it.
  if(ASTEqual(openToken,AST_type::TK_OPENPAREN)) {
#ifdef VERBOSE
  cerr << "parseExpressionPartial: Found '(' parsing new expression" << endl ;
#endif    
    AST_type::ASTP exp = parseExpression(is,linecount,fileName,typemap) ;
    CPTR<AST_Token> closeToken = getToken(is,linecount) ;
    CPTR<AST_exprOper> group = new AST_exprOper ;
    group->nodeType = AST_type::OP_GROUP ;
    group->terms.push_back(AST_type::ASTP(exp)) ;
    if(closeToken->nodeType != AST_type::TK_CLOSEPAREN) {
      pushToken(closeToken) ;
      return AST_type::ASTP(new AST_syntaxError("expecting ')'",closeToken->lineno,fileName)) ;
      
    }
    return applyPostFixOperator(AST_type::ASTP(group),is,linecount,fileName,
				typemap) ;
  }
  // Check for unary operators
  if(checkUnaryToken(openToken)) {
#ifdef VERBOSE
    cerr << "parseExpressionPartial: found unary operator '"
         << openToken->text << " parsing unary target." 
         << endl ;
#endif    
    //check for unary operators
    AST_type::ASTP expr = parseExpressionPartial(is,linecount,fileName,
						 typemap) ;
    if(expr!= 0) {
      CPTR<AST_exprOper> unary = new AST_exprOper ;
      unary->nodeType = unaryOperator(openToken->nodeType) ;
      unary->terms.push_back(expr) ;
#ifdef VERBOSE
      cerr << "parseExpressionPartial: creating unary operator, "
           << "search for postfix operators ++,--, or []" 
           << endl ;
#endif    
      return applyPostFixOperator(AST_type::ASTP(unary),is,linecount,fileName,
				  typemap) ;
    }
  }
  // Check if we are parsing a scoped name
  if(ASTEqual(openToken,AST_type::TK_NAME) ||
     ASTEqual(openToken,AST_type::TK_SCOPE)) {
    pushToken(openToken) ;
    AST_type::ASTP objname = parseIdentifier(is,linecount,fileName,typemap) ;
    return applyPostFixOperator(objname,is,linecount,fileName, typemap) ;
    
    //    AST_type::ASTP objname = parseScopedObject(is,linecount,fileName,typemap) ;

    openToken = getToken(is,linecount) ;
    if(ASTEqual(openToken,AST_type::TK_OPENTEMPLATE)) {
      objname = parseTemplateArguments(objname,is,linecount,fileName,typemap) ;
    } else
      pushToken(openToken) ;
    
    openToken = getToken(is,linecount) ;
    if(ASTEqual(openToken,AST_type::TK_OPENPAREN)) {
      objname = parseFunctionArguments(objname,is,linecount,fileName,typemap) ;
    } else
      pushToken(openToken) ;

    return applyPostFixOperator(AST_type::ASTP(objname),is,linecount,fileName,
				  typemap) ;
    
  }
      
  if(isTerm(openToken)) {
#ifdef VERBOSE
    cerr << "parseExpressionPartial: found term '" 
         << openToken->text << "', checking for following parenthesis"
         << endl ;
#endif    
    pushToken(openToken) ;
    AST_type::ASTP exp = parseTerm(is,linecount);


    return applyPostFixOperator(exp,is,linecount,fileName,typemap) ;
    //    return exp ;
  }
  pushToken(openToken) ;
  return 0 ;
}

/// get precedence of operator code
inline int getPrecedence(AST_type::elementType op) {
  // The operator precedence is stored in the lowest 8 bits of the
  // operator code.  The mask is used to extract the precidence of the
  // operator code
  const unsigned int mask = ~0x7f ;
  FATAL(op > AST_type::OP_NIL) ;
  return op & mask ;
}

/// Interfce for general ASTP pointer
inline int getPrecedence(AST_type::ASTP op) {
  return getPrecedence(op->nodeType) ;
}
/// Interface for expression AST node
inline int getPrecedence(CPTR<AST_exprOper> op) {
  return getPrecedence(op->nodeType) ;
}

AST_type::ASTP parseExpressionOperator(AST_type::ASTP expr,
                                       std::istream &is, int &linecount,
                                       const string &fileName,
                                       varmap &typemap,
                                       AST_type::elementType
                                       term) {
  // exprStack holds already processed binary operators
  // create stack and push sentinel on the stack

  // Note that if we have repeated operators that have default left
  // associativity (e.g. a+b+c+d = (a+(b+(c+d))) then this will be structured
  // as a single plus node in the tree, e.g. (+:a b c d)
  
  vector<CPTR<AST_exprOper> > exprStack ;
  {
    CPTR<AST_exprOper> tmp = new AST_exprOper ;
    tmp->nodeType = AST_type::OP_NIL ;
    tmp->terms.push_back(expr) ;
    exprStack.push_back(tmp) ;
#ifdef VERBOSE
    cerr << "push op = " << OPtoName(exprStack.back()->nodeType)  << endl ;
#endif
  }
  // After getting the first term we are in a loop of searching for operators
  do {
    // Check if we detect a terminating token
    CPTR<AST_Token> openToken = getToken(is,linecount) ;
    if(openToken->nodeType == term) {
      // We found it, restore token and break out of loop
      pushToken(openToken) ;
      break ;
    }
    // restore token for parseOperator
    pushToken(openToken) ;
    // binary operator check
    AST_type::ASTP op = parseOperator(is, linecount) ;
    if(op == 0)
      break ;

#ifdef VERBOSE
    cerr << "parseExpression operator op = " << OPtoString(op->nodeType) << endl ;
#endif
    // Special case for ?: operator or array
    if(ASTEqual(op,AST_type::OP_TERNARY)) 
      expr = parseExpression(is,linecount,fileName,typemap) ;
    else 
      expr = parseExpressionPartial(is,linecount,fileName,typemap) ;
    
    if(expr == 0) {
      return AST_type::ASTP(new AST_syntaxError("Expecting expression after binary operator",linecount,fileName)) ;
    }
    
    if(ASTEqual(exprStack.back(),AST_type::OP_NIL)) {
      // If no operator is parsed yet, we just get started and initilize
      // the left branch
      exprStack.back()->nodeType = op->nodeType ;
      exprStack.back()->terms.push_back(expr) ;
#ifdef VERBOSE
      cerr << "replacing OP_NIL on exprStack with "
           << OPtoName(exprStack.back()->nodeType)
           << " line " << __LINE__
           << endl ;
#endif
    } else {
      // Now we reorder the tree based on operator precedence
      while(exprStack.size() > 1 &&
            getPrecedence(op) >= getPrecedence(exprStack.back())) {
#ifdef VERBOSE
	if(ASTEqual(op,AST_type::OP_TERNARY)) {
	  cerr << "pop stack when OP_TERNARY" << endl ;
	}
#endif
#ifdef VERBOSE
    cerr << "pop off exprStack:" << OPtoName(exprStack.back()->nodeType)  << endl ;
#endif
	exprStack.pop_back() ;
#ifdef VERBOSE
    cerr << "after pop top of stack is:" << OPtoName(exprStack.back()->nodeType)
         << ", size=" << exprStack.size() << endl ;
#endif
        
      }
      
      if(ASTEqual(op,exprStack.back())) {
	// If operator is the same, just chain the terms
	exprStack.back()->terms.push_back(expr) ;
#ifdef VERBOSE
        cerr << "adding term to operator: "
             << OPtoName(exprStack.back()->nodeType)
             << " stack size=" << exprStack.size() 
             << endl ;
#endif
      } else if(getPrecedence(op) < getPrecedence(exprStack.back())) {
        // if operator is lower precedence than top of stack, then we need
        // to put this operator on the rightmost term
        CPTR<AST_exprOper> np = new AST_exprOper ;
        np->nodeType = op->nodeType ;
#ifdef VERBOSE
        cerr << "top of stack is" << OPtoName(exprStack.back()->nodeType)
             << endl ;
        cerr << "editing rightmost term: " << OPtoName(exprStack.back()->terms.back()->nodeType)
             << endl ;
#endif
        np->terms.push_back(exprStack.back()->terms.back()) ;
        np->terms.push_back(expr) ;
        exprStack.back()->terms.back() = AST_type::ASTP(np) ;
#ifdef VERBOSE
        cerr << "pushing " << OPtoName(exprStack.back()->nodeType)
             << "to stack,  line " << __LINE__
             << endl ;
#endif
        if(ASTEqual(op,AST_type::OP_TERNARY)) {
          CPTR<AST_Token> op2 = getToken(is,linecount) ;
#ifdef VERBOSE
          cerr << " detected OP_TERNARY, op2 =" << op2->text << endl ;
#endif
          if(ASTEqual(op2,AST_type::TK_COLON)) {
            expr = parseExpressionPartial(is,linecount,fileName,typemap) ;
            np->terms.push_back(expr) ;
          } else {
            pushToken(op2) ;
            return AST_type::ASTP(new AST_syntaxError("expecting ':' in tertiary operator",op2->lineno,fileName)) ;
          }
        }
        exprStack.push_back(np) ;
#ifdef VERBOSE
        cerr << "pushing on exprStack:" << OPtoName(exprStack.back()->nodeType)
             << " stack size=" << exprStack.size() << endl ;
#endif
      } else {
	CPTR<AST_exprOper> np = new AST_exprOper ;
	np->nodeType = op->nodeType ;
#ifdef VERBOSE
        cerr << "creating lower precidence op="
             << OPtoName(np->nodeType)  << endl ;
#endif        
	np->terms.push_back(AST_type::ASTP(exprStack.back())) ;
	np->terms.push_back(expr) ;
	if(ASTEqual(op,AST_type::OP_TERNARY)) {
	  CPTR<AST_Token> op2 = getToken(is,linecount) ;
	  if(!ASTEqual(op2,AST_type::TK_COLON)) {
	    pushToken(op2) ;
	    expr = AST_type::ASTP(new AST_syntaxError("unexpected ':' in tertiary operator",linecount,fileName)) ;
	  } else {
	      expr = parseExpressionPartial(is,linecount,fileName,typemap) ;
	  }
	  np->terms.push_back(expr) ;
#ifdef VERBOSE
	  cerr << " processed OP_TERNARY" << endl ;
#endif
	  
	} 
	exprStack.back() = np ;
#ifdef VERBOSE
        cerr << "changing exprStack to" << OPtoName(exprStack.back()->nodeType)
             << " line " << __LINE__
             << endl ;
#endif
      }
    }
  } while(true) ;

  AST_type::ASTP e = CPTR<AST_type>(exprStack.front()) ;
  if(ASTEqual(exprStack.front(),AST_type::OP_NIL)) 
    e = CPTR<AST_type>(exprStack.front()->terms.front()) ;

  return applyPostFixOperator(AST_type::ASTP(e),is,linecount,fileName,typemap) ;
}

/// Parse general expression.  This routine mainly handles operator
/// precedence for binary operators, it passes other work to
/// parseExpressionPartial
AST_type::ASTP parseExpression(std::istream &is, int &linecount,
			       const string &fileName,
			       varmap &typemap) {
#ifdef VERBOSE
  cerr << "in parseExpression"<< endl ;
#endif
  AST_type::ASTP expr = parseExpressionPartial(is,linecount,fileName,typemap) ;
  if(expr == 0) // If no valid expression then return null
    return expr ;
  return parseExpressionOperator(expr,is,linecount,fileName,typemap) ;
  
}


AST_type::ASTP parseCaseStatement(std::istream &is, int &linecount,
				  const string &fileName,
				  varmap &typemap) {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseCaseStatement, token = " << token->text << endl ;
#endif
  if(!ASTEqual(token,AST_type::TK_CASE) &&
     !ASTEqual(token,AST_type::TK_DEFAULT)) {
    cerr << "internal error parsing switch statement on line " << linecount
	 << " of file " << fileName << endl ;
  }
  CPTR<AST_Block> AST_data = new AST_Block ;
  AST_data->nodeType = AST_type::TK_CASE ;
  AST_data->identifiers = typemap ;
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  if(ASTEqual(token,AST_type::TK_CASE)) {
    AST_type::ASTP expr =
      parseExpression(is,linecount,fileName,AST_data->identifiers) ;
    AST_data->elements.push_back(expr) ;
  }
  
  token = getToken(is,linecount) ;
  if(!ASTEqual(token,AST_type::TK_COLON)) {
    pushToken(token) ; 
    return AST_type::ASTP(new AST_syntaxError("expecting ':' in case statement",token->lineno,fileName)) ;
  }
  
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  return AST_type::ASTP(AST_data) ;
}

AST_type::ASTP parseSwitchStatement(std::istream &is, int &linecount,
				    const string &fileName,
				    varmap &typemap)  {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseSwitchStatement, token = " << token->text << endl ;
#endif
  AST_type::ASTP switchtok(token) ;
  CPTR<AST_controlStatement> ctrl = new AST_controlStatement ;
  ctrl->identifiers = typemap  ;
  ctrl->controlType = AST_type::ASTP(token) ;
  ctrl->nodeType = AST_type::ND_CTRL_SWITCH ;

  token = getToken(is,linecount) ;
  if(!ASTEqual(token,AST_type::TK_OPENPAREN)) {
    pushToken(token) ;
    return AST_type::ASTP(new AST_syntaxError("expecting '(' in switch statment",token->lineno,fileName)) ;
  }
  pushToken(token) ;
  AST_type::ASTP conditional =
    parseExpression(is,linecount,fileName,ctrl->identifiers) ;
  // switch statement is just a list of statements in order
  ctrl->parts.push_back(conditional) ;
  token = getToken(is,linecount) ;
  if(!ASTEqual(token,AST_type::TK_OPENBRACE)) {
    pushToken(token) ;
    return AST_type::ASTP(new AST_syntaxError("expecting '{' in switch statment",token->lineno,fileName)) ;
  }
  ctrl->parts.push_back(AST_type::ASTP(token)) ;
  token = getToken(is,linecount) ;
  while(!ASTEqual(token,AST_type::TK_CLOSEBRACE) &&
	!ASTEqual(token,AST_type::TK_ERROR)) {
    if(ASTEqual(token,AST_type::TK_CASE) ||
       ASTEqual(token,AST_type::TK_DEFAULT)) {
      pushToken(token) ;
      ctrl->parts.push_back(parseCaseStatement(is,linecount,fileName,
                                               ctrl->identifiers)) ;
    } else {
      pushToken(token) ;
      ctrl->parts.push_back(parseStatement(is,linecount,fileName,
                                           ctrl->identifiers)) ;
    }
    token = getToken(is,linecount) ;
  }
			
  if(ASTEqual(token,AST_type::TK_ERROR)) 
    return AST_type::ASTP(new AST_syntaxError("failed to find closing brace in swithc statement",linecount,fileName)) ;

  ctrl->parts.push_back(AST_type::ASTP(token)) ;
  return AST_type::ASTP(ctrl) ;

}

AST_type::ASTP parseIfStatement(std::istream &is, int &linecount,
				const string &fileName,
				varmap &typemap)  {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseIfStatement, token = " << token->text << endl ;
#endif
  if(!ASTEqual(token,AST_type::TK_IF)) {
    pushToken(token) ;
    return AST_type::ASTP(new AST_syntaxError("confused in if statement",token->lineno,fileName)) ;
  }
  AST_type::ASTP iftok = AST_type::ASTP(token) ;
  
  token = getToken(is,linecount) ;
  if(!ASTEqual(token,AST_type::TK_OPENPAREN)) {
    pushToken(token) ;
    return AST_type::ASTP(new AST_syntaxError("if expecting '('",token->lineno,fileName)) ;
  }
  pushToken(token) ;
  AST_type::ASTP conditional = parseExpression(is,linecount,fileName,typemap) ;
  if(conditional == 0) 
    return AST_type::ASTP(new AST_syntaxError("malformed if conditional",token->lineno,fileName)) ;

  AST_type::ASTP body = parseStatement(is,linecount,fileName,typemap) ;

  AST_type::ASTP elsetok = 0 ;
  AST_type::ASTP ebody = 0 ;
  token = getToken(is,linecount) ;
  
  if(ASTEqual(token,AST_type::TK_ELSE)) {
    elsetok = AST_type::ASTP(token) ;
    ebody = parseStatement(is,linecount,fileName,typemap) ;
  } else {
    pushToken(token) ;
  }
  CPTR<AST_controlStatement> ctrl = new AST_controlStatement ;
  ctrl->identifiers = typemap ;
  ctrl->constructIf(iftok,conditional,body,elsetok,ebody) ;
  return AST_type::ASTP(ctrl) ;
}			 

AST_type::ASTP parseLoopStatement(std::istream &is, int &linecount,
				  const string &fileName,
				  varmap &typemap) {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseLoopStatement, token = " << token->text << endl ;
#endif
  AST_type::ASTP loop = AST_type::ASTP(token) ;
  CPTR<AST_controlStatement> ctrl = new AST_controlStatement ;
  ctrl->identifiers = typemap ;
  ctrl->controlType = loop ;
  switch(token->nodeType) {
  case AST_type::TK_FOR:
    {

      //--------------------------------------------------------------------
      // parsing for loop
      //
      // Note for the for loop the parts will be
      // TK_OPENPAREN
      // initializer (semicolon included)
      // conditional
      // TK_SEMICOLON
      // advance
      // TK_CLOSEPAREN
      // for body
      
      ctrl->nodeType = AST_type::ND_CTRL_FOR ;
      
      token = getToken(is,linecount) ;
      if(!ASTEqual(token,AST_type::TK_OPENPAREN)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("for expecting '('",token->lineno,fileName)) ;
     } 
      ctrl->parts.push_back(AST_type::ASTP(token)) ;
      token = getToken(is,linecount) ;
#ifdef VERBOSE
      cerr << "for loop parsed TK_OPEN_PAREN, next token = " << OPtoName(token->nodeType) << endl ;
#endif
      
      AST_type::ASTP initializer = 0 ;
      if(ASTEqual(token,AST_type::TK_SEMICOLON)) {
        initializer = AST_type::ASTP(token) ;
      } else {
        CPTR<AST_Token> token2 = getToken(is,linecount) ;

        bool isDeclaration = true ;
        if(ASTEqual(token,AST_type::TK_NAME) &&
           !ASTEqual(token2,AST_type::TK_SCOPE)) {
          auto ii = typemap.find(token->text) ;
          if(ii != typemap.end() && ii->second.isLocalIdentifier())
            isDeclaration = false ;
        }

        pushToken(token2) ;
	pushToken(token) ;

        if(isDeclaration)
          initializer =
            parseDeclaration(is,linecount,fileName,ctrl->identifiers) ;
        else
          initializer =
            parseSimpleStatement(is,linecount,fileName,ctrl->identifiers) ;
      }
#ifdef VERBOSE
      cerr << "initializer = " << OPtoName(initializer->nodeType) << endl ;
#endif      
      ctrl->parts.push_back(initializer) ;
	
#ifdef VERBOSE
      cerr << "for: parsing conditional expression" << endl;
#endif
      AST_type::ASTP conditional = parseExpression(is,linecount,fileName,
						   ctrl->identifiers) ;
      
#ifdef VERBOSE
      cerr << "conditional = " << OPtoName(conditional->nodeType) << endl ;
#endif      
      ctrl->parts.push_back(conditional) ;
      token = getToken(is,linecount) ;
      if(!ASTEqual(token,AST_type::TK_SEMICOLON)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("No ';' after conditional in for loop",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;
      AST_type::ASTP advance = parseExpression(is,linecount,fileName,
                                               ctrl->identifiers) ;
      ctrl->parts.push_back(advance) ;
      token = getToken(is,linecount) ;
      if(!ASTEqual(token,AST_type::TK_CLOSEPAREN)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting ')' in for loop",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;
      AST_type::ASTP body = parseStatement(is,linecount,fileName,
                                           ctrl->identifiers) ;
      ctrl->parts.push_back(body) ;
      return AST_type::ASTP(ctrl) ;
    }
    break ;
  case AST_type::TK_WHILE:
    {
      //--------------------------------------------------------------------
      // parsing while loop
      //
      // Note parts will be
      // conditional
      // loop body
      ctrl->nodeType = AST_type::ND_CTRL_WHILE ;

      token = getToken(is,linecount) ;
      if(!ASTEqual(token,AST_type::TK_OPENPAREN)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting '(' in while loop",token->lineno,fileName)) ;
      }
      pushToken(token) ;
      AST_type::ASTP conditional = parseExpression(is,linecount,fileName,
						   ctrl->identifiers) ;
      ctrl->parts.push_back(conditional) ;

      AST_type::ASTP body = parseStatement(is,linecount,fileName,
                                           ctrl->identifiers) ;
      ctrl->parts.push_back(body) ;
      return AST_type::ASTP(ctrl) ;
    }
    break ;
  case AST_type::TK_DO:
    {
      //--------------------------------------------------------------------
      // parsing do loop
      //
      // Note parts will be
      // loop body
      // TK_WHILE
      // conditional
      // TK_SEMICOLON

      ctrl->nodeType = AST_type::ND_CTRL_DO ;
      AST_type::ASTP body = parseStatement(is,linecount,fileName,
                                           ctrl->identifiers) ;
      ctrl->parts.push_back(body) ;
      token = getToken(is,linecount) ;
      if(!ASTEqual(token,AST_type::TK_WHILE)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting 'while' in do loop",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;

      token = getToken(is,linecount) ;
      if(!ASTEqual(token,AST_type::TK_OPENPAREN)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting '(' in do loop",token->lineno,fileName)) ;
      }
      pushToken(token) ;
      AST_type::ASTP conditional = parseExpression(is,linecount,fileName,
                                                   ctrl->identifiers) ;
      ctrl->parts.push_back(conditional) ;

      token = getToken(is,linecount) ;
      if(!ASTEqual(token,AST_type::TK_SEMICOLON)) {
	pushToken(token) ;
	
	ctrl->parts.push_back(AST_type::ASTP(new AST_syntaxError("Expecting ';' in do loop",token->lineno,fileName))) ;
      } else
	ctrl->parts.push_back(AST_type::ASTP(token)) ;
      return AST_type::ASTP(ctrl) ;
    }
    break ;
  default:
    return 0 ;
  }
  
  return AST_type::ASTP(getToken(is,linecount)) ;
}


string getIdentifierName(std::istream &is,
                         int &linecount,
                         const string &fileName,
                         bool &isFunc) {
#ifdef VERBOSE
  cerr << "in getIdentifierName" << endl ;
#endif
  vector<CPTR<AST_Token>>  token_stack ;
  isFunc = false ;
  string name ;
  CPTR<AST_Token> token = getToken(is,linecount) ;
  // parse the scoped name part
  if(ASTEqual(token,AST_type::TK_SCOPE)) {
    name += token->text ;
    token_stack.push_back(token) ;
    token = getToken(is,linecount) ;
  }
  while(ASTEqual(token,AST_type::TK_NAME)) {
    name += token->text ;
    token_stack.push_back(token) ;
    token = getToken(is,linecount) ;
    if(ASTEqual(token,AST_type::TK_SCOPE)) {
      name += token->text ;
      token_stack.push_back(token) ;
      token = getToken(is,linecount) ;
    } else
      break ;
  }
#ifdef VERBOSE
  cerr << "ident name = " << name << endl ;
#endif
  // if template, parse template arguments
  if(ASTEqual(token,AST_type::TK_OPENTEMPLATE)) {
#ifdef VERBOSE
    cerr << "identifier is template" << endl ;
#endif
    int count = 1 ;
    name += "<" ;
    token_stack.push_back(token) ;
    token = getToken(is,linecount) ;
    while(count > 0) {
      if(ASTEqual(token,AST_type::TK_OPENTEMPLATE)) {
        count++ ;
        name += "<" ;
      } else if(ASTEqual(token,AST_type::TK_CLOSETEMPLATE)) {
        count-- ;
        name += ">" ;
      } else {
        name += token->text ;
      }
      token_stack.push_back(token) ;
      token = getToken(is,linecount) ;
    }
  }
  if(ASTEqual(token,AST_type::TK_OPENPAREN)) {
    isFunc=true ;
  }
  token_stack.push_back(token) ;
  // Push tokens so that they can be parsed
  for(auto ii = token_stack.rbegin();ii!=token_stack.rend();++ii)
    pushToken(*ii) ;
#ifdef VERBOSE
  cerr << "return ident name = " << name << endl ;
#endif
  return name ;
}
    
      
    
/// Determine if statement is a declration or a simple statement by
/// looking up potential typename in varmap.  We assume that the first
/// token is either an TK_NAME or TK_SCOPE as a prerequisite
AST_type::ASTP parseDeclarationOrSimpleStatement(std::istream &is,
                                                 int &linecount,
                                                 const string &fileName,
                                                 varmap &typemap) {
#ifdef VERBOSE
  cerr << "in parseDeclarationOrSimpleStatement"<< endl;
#endif
  bool isFunc = false ;
  string ident = getIdentifierName(is,linecount,fileName,isFunc) ;
#ifdef VERBOASE
  cerr << "ident = " << ident << "isFunc = " << isFunc << endl ;
#endif
  if(isFunc)
    return parseSimpleStatement(is,linecount,fileName,typemap) ;
  auto ii = typemap.find(ident) ;
  if(ii != typemap.end() && ii->second.isLocalIdentifier())
    return parseSimpleStatement(is,linecount,fileName,typemap) ;
  return parseDeclaration(is,linecount,fileName,typemap) ;
}

// Parse type declaration when we know we are expecting a type declaration
AST_type::ASTP parseDeclaration(std::istream &is, int &linecount,
				const string &fileName,
				varmap &typemap) {

  CPTR<AST_declaration> AST_data = new AST_declaration ;

  CPTR<AST_Token> token = getToken(is,linecount) ;

  bool istype = true ;
  bool builtin_type = false ;
  bool defined_type = false ;

  do {
#ifdef VERBOSE
  cerr << "in parseDeclaration, token = " << token->text << endl ;
#endif
    switch(token->nodeType) {
    case AST_type::TK_CONST:
      // TK_REGISTER,TK_VOLATILE,TK_MUTABLE add later
      AST_data->type_decl.push_back(AST_type::ASTP(token)) ;
#ifdef VERBOSE
      cerr << "in parseDeclaration, const variable" << endl ;
#endif
      break ;
    case AST_type::TK_EXTERN:
      AST_data->type_decl.push_back(AST_type::ASTP(token)) ;
#ifdef VERBOSE
      cerr << "in parseDeclaration, extern variable" << endl ;
#endif
      break ;
    case AST_type::TK_CHAR:
    case AST_type::TK_FLOAT:
    case AST_type::TK_DOUBLE:
    case AST_type::TK_INT:
    case AST_type::TK_BOOL:
    case AST_type::TK_SHORT:
    case AST_type::TK_LONG:
    case AST_type::TK_SIGNED:
    case AST_type::TK_UNSIGNED:
    case AST_type::TK_AUTO:
#ifdef VERBOSE
      cerr << "in parseDeclaration, builtin type" << endl ;
#endif
      builtin_type = true ;
      if(defined_type) { // This is a syntax error
        AST_type::ASTP p =
          AST_type::ASTP(new AST_syntaxError("Mixed concrete and defined types",
                                             token->lineno,fileName)) ;
        AST_data->type_decl.push_back(p) ;
        istype = false ;
      }
      AST_data->type_decl.push_back(AST_type::ASTP(token)) ;
      break ;
    case AST_type::TK_SCOPE:
    case AST_type::TK_NAME:
#ifdef VERBOSE
      cerr << "in parseDeclaration, named type" << endl ;
#endif
      if(builtin_type || defined_type ) {
        istype = false ;
      } else {
        pushToken(token) ;
        
        AST_type::ASTP typedec =
          parseIdentifier(is,linecount,fileName,typemap) ;
        if(typedec->nodeType == AST_type::OP_FUNC) {
          // not a proper type name
          AST_type::ASTP p =
            AST_type::ASTP(new AST_syntaxError("Mixed concrete and defined types",
                                               token->lineno,fileName)) ;
          AST_data->type_decl.push_back(p) ;
        } else {

          // right now only consider unscoped non-templated type names
          defined_type = true ;
          AST_data->type_decl.push_back(typedec) ;
        }
      }
      break ;
    default:
#ifdef VERBOSE
      cerr << "in parseDeclaration, finished type parsing" << endl ;
#endif
      istype = false ;
      break ;
    }
    if(istype)
      token = getToken(is,linecount) ;
  } while(istype) ;

#ifdef VERBOSE
  cerr << "in parseDeclaration, processing declarations: " << token->text
       << " " << OPtoName(token->nodeType) << endl;
#endif
  while(ASTEqual(token,AST_type::TK_NAME) ||
        ASTEqual(token,AST_type::TK_TIMES) ||
        ASTEqual(token,AST_type::TK_AND)) {
#ifdef VERBOSE
    cerr << "in parseDeclaration, found name, * or & " << endl ;
#endif

    if(ASTEqual(token,AST_type::TK_TIMES)) {
      AST_data->decls.push_back(AST_type::ASTP(token)) ;
      token = getToken(is,linecount) ;
      if(ASTEqual(token,AST_type::TK_CONST)) {
        AST_data->decls.push_back(AST_type::ASTP(token)) ;
        token = getToken(is,linecount) ;
      }
    }
    if(ASTEqual(token,AST_type::TK_AND)) {
      AST_data->decls.push_back(AST_type::ASTP(token)) ;
      token = getToken(is,linecount) ;
    }
    
#ifdef VERBOSE
  cerr << "in parseDeclaration, processed pointer or reference: " << token->text << endl ;
#endif
    if(ASTEqual(token,AST_type::TK_NAME)) {
      typemap[token->text] = localIdentifier() ;
      AST_data->decl_varnames.push_back(token->text) ;
      pushToken(token) ;
      
      AST_type::ASTP expr = parseExpressionPartial(is,linecount,fileName,typemap) ;
#ifdef VERBOSE
  cerr << "in parseDeclaration, parsing expression " <<  endl ;
#endif
      if(expr != 0) { // If no valid expression then return null
        expr = parseExpressionOperator(expr,is,linecount,fileName,typemap,
                                       AST_type::TK_COMMA) ;
      } else {
#ifdef VERBOSE
  cerr << "in parseDeclaration, parseExpressionPartial returned nil " <<  endl ;
#endif
        
        AST_type::ASTP p =
          AST_type::ASTP(new AST_syntaxError("Syntax error variable declaration",
                                         token->lineno,fileName)) ;
        AST_data->decls.push_back(p) ;
      }
      AST_data->decls.push_back(expr) ;
      token = getToken(is,linecount) ;
    }
    if(ASTEqual(token,AST_type::TK_COMMA)) {
      AST_data->decls.push_back(AST_type::ASTP(token)) ;
      token = getToken(is,linecount) ;
    }
  }
  if(ASTEqual(token,AST_type::TK_SEMICOLON)) {
    AST_data->decls.push_back(AST_type::ASTP(token)) ;
  } else {
    AST_type::ASTP p =
      AST_type::ASTP(new AST_syntaxError("Expected ';' while parsing variable declaration",
                                         token->lineno,fileName)) ;
    AST_data->decls.push_back(p) ;
  }

  return AST_type::ASTP(AST_data) ;
}


/// Parse return, continue, or break control statements
AST_type::ASTP parseSpecialControlStatement(std::istream &is, int &linecount,
					    const string &fileName,
					    varmap &typemap) {
  
  CPTR<AST_Block> AST_data = new AST_Block ;
  AST_data->identifiers = typemap ;
  AST_data->nodeType = AST_type::OP_SPECIAL ;
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseSpecialControlStatement, token = " << token->text << endl ;
#endif
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  token = getToken(is,linecount) ;
  if(!ASTEqual(token,AST_type::TK_SEMICOLON)) {
    if(ASTEqual(AST_data->elements.back(),AST_type::TK_RETURN)) {
      pushToken(token) ;
      AST_type::ASTP exp =
        parseExpression(is,linecount,fileName,AST_data->identifiers) ;
      AST_data->elements.push_back(exp) ;
      token = getToken(is,linecount) ;
      if(!ASTEqual(token,AST_type::TK_SEMICOLON)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting ';' after return",token->lineno,fileName)) ;
      }
    } else {
      pushToken(token) ;
      return AST_type::ASTP(new AST_syntaxError("unexpected ';'" ,token->lineno,fileName)) ;
    }
  }
  AST_data->elements.push_back(AST_type::ASTP(token)) ;

  return AST_type::ASTP(AST_data) ;
}

AST_type::ASTP parseSimpleStatement(std::istream &is, int &linecount,
                                    const string &fileName,
                                    varmap &typemap) {
#ifdef VERBOSE
  cerr << "in parseSimpleStatement" << endl ;
#endif
  AST_type::ASTP exp = parseExpression(is,linecount,fileName,typemap) ;
  
  CPTR<AST_Token> termToken = getToken(is,linecount) ;
  AST_type::ASTP term = AST_type::ASTP(termToken) ;
#ifdef VERBOSE
  cerr << "terminal token for simpleStatement is " << OPtoName(term->nodeType)
       << endl ;
#endif
  if(!ASTEqual(term,AST_type::TK_SEMICOLON)) {
    pushToken(termToken) ;
    return AST_type::ASTP(new AST_syntaxError("Expecting ';' ",
                                              termToken->lineno,fileName)) ;
  }
  AST_type::ASTP stat = new AST_SimpleStatement(exp,term) ;
  return stat ;
}

// Parse an inut statement, this could be a code block denoted by braces,
// A type declaration
// A loop statement
// A if statement
// A switch statement
// A special control statement (e.g. break, control, or return
// A simple statement (expression) followed by a semicolon
// An empty statement
AST_type::ASTP parseStatement(std::istream &is, int &linecount,
			      const string &fileName,
			      varmap &typemap) {
  CPTR<AST_Token> firstToken = getToken(is,linecount) ;

#ifdef VERBOSE
  cerr << "in parseStatement, token = " << firstToken->text << endl ;
#endif
  pushToken(firstToken) ;
  switch(firstToken->nodeType) {
  case AST_type::TK_OPENBRACE:
    return parseBlock(is,linecount,fileName,typemap) ;
  case AST_type::TK_CHAR:
  case AST_type::TK_FLOAT:
  case AST_type::TK_DOUBLE:
  case AST_type::TK_INT:
  case AST_type::TK_BOOL:
  case AST_type::TK_LONG:
  case AST_type::TK_SIGNED:
  case AST_type::TK_UNSIGNED:
  case AST_type::TK_CONST:
  case AST_type::TK_EXTERN:
  case AST_type::TK_AUTO:
    return parseDeclaration(is,linecount,fileName,typemap) ;
  case AST_type::TK_FOR:
  case AST_type::TK_WHILE:
  case AST_type::TK_DO:
    return parseLoopStatement(is,linecount,fileName,typemap) ;
  case AST_type::TK_IF:
    return parseIfStatement(is,linecount,fileName,typemap) ;
  case AST_type::TK_SWITCH:
    return parseSwitchStatement(is,linecount,fileName,typemap) ;
  case AST_type::TK_BREAK:
  case AST_type::TK_CONTINUE:
  case AST_type::TK_RETURN:
    return parseSpecialControlStatement(is,linecount,fileName,typemap) ;
  case AST_type::TK_NAME:
  case AST_type::TK_SCOPE:
    return parseDeclarationOrSimpleStatement(is,linecount,fileName,typemap) ;
    {
      CPTR<AST_Token> tok1 = getToken(is,linecount) ;
      CPTR<AST_Token> tok2 = getToken(is,linecount) ;
      bool isDeclaration = true ;
      if(ASTEqual(tok1,AST_type::TK_NAME) &&
         !ASTEqual(tok2,AST_type::TK_SCOPE)) {
        auto ii = typemap.find(tok1->text) ;
        if(ii != typemap.end() && ii->second.isLocalIdentifier())
          isDeclaration = false ;
      }
      
      pushToken(tok2) ;
      pushToken(tok1) ;
      if(isDeclaration)
	return parseDeclaration(is,linecount,fileName,typemap) ;
      return parseSimpleStatement(is,linecount,fileName,typemap) ;
      AST_type::ASTP exp = parseExpression(is,linecount,fileName,typemap) ;
      
      CPTR<AST_Token> termToken = getToken(is,linecount) ;
      AST_type::ASTP term = AST_type::ASTP(termToken) ;
      if(!ASTEqual(term,AST_type::TK_SEMICOLON)) {
	pushToken(termToken) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting ';' ",
						  termToken->lineno,fileName)) ;
      }
      AST_type::ASTP stat = new AST_SimpleStatement(exp,term) ;
      return stat ;
    }

  case AST_type::TK_USING:
    {
      firstToken = getToken(is,linecount) ;
      CPTR<AST_declaration> AST_data = new AST_declaration ;
      AST_data->type_decl.push_back(AST_type::ASTP(firstToken)) ;
      
      bool isFunc = false ;
      string s = getIdentifierName(is,linecount,fileName,isFunc) ;
      if(!isFunc) {
        typemap[s] = localIdentifier() ;
      }
      AST_type::ASTP typedec =
        parseIdentifier(is,linecount,fileName,typemap) ;
      AST_data->type_decl.push_back(typedec) ;
      firstToken = getToken(is,linecount) ;
      if(!ASTEqual(firstToken,AST_type::TK_SEMICOLON)) {
	pushToken(firstToken) ;
        AST_type::ASTP err =
          AST_type::ASTP(new AST_syntaxError("Expecting ';' ",
                                             firstToken->lineno,fileName)) ;
        AST_data->type_decl.push_back(err) ;
      } else {
        AST_data->type_decl.push_back(AST_type::ASTP(firstToken)) ;
      }
      return AST_type::ASTP(AST_data) ;
    }
    break ;
    
  case AST_type::TK_OPENPAREN:
  case AST_type::TK_PLUS:
  case AST_type::TK_MINUS:
  case AST_type::TK_NOT:
  case AST_type::TK_AND:
  case AST_type::TK_TIMES:
  case AST_type::TK_INCREMENT:
  case AST_type::TK_DECREMENT:
  case AST_type::TK_LOCI_VARIABLE:
    {
      return parseSimpleStatement(is,linecount,fileName,typemap) ;
    } 
    
  case AST_type::TK_SEMICOLON:
    firstToken = getToken(is,linecount) ;
    return AST_type::ASTP(firstToken) ;
  default:
    break ;
  }
  
  firstToken = getToken(is,linecount) ;
  return AST_type::ASTP(new AST_syntaxError("Expecting statement or ';' ",
					    firstToken->lineno,
					    fileName)) ;
}

/// Parse a block of statements enclosed in braces
AST_type::ASTP parseBlock(std::istream &is, int &linecount,
			  const string &fileName,
			  varmap &typemap) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseBlock, token = " << openToken->text << endl ;
#endif

  AST_type::elementType closeType = AST_type::TK_CLOSEBRACE ;
  switch(openToken->nodeType) {
  case AST_type::TK_OPENBRACE:
    closeType = AST_type::TK_CLOSEBRACE ;
    break ;
  default:
    return AST_type::ASTP(openToken) ;
  }

  CPTR<AST_Block> AST_data = new AST_Block ;
  AST_data->identifiers = typemap ;
  AST_data->nodeType = AST_type::TK_BRACEBLOCK ;
  AST_data->elements.push_back(AST_type::ASTP(openToken)) ;
  CPTR<AST_Token> token = getToken(is,linecount) ;
  while(!ASTEqual(token,closeType)) {
    pushToken(token) ;
    CPTR<AST_type> statement =
      parseStatement(is,linecount,fileName,AST_data->identifiers) ;
    AST_data->elements.push_back(statement) ;
    token = getToken(is,linecount) ;
    if(is.fail() || is.eof()) 
      break ;
  }
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  return AST_type::ASTP(AST_data) ;
}

void AST_visitor::visit(AST_SimpleStatement &s) {
  if(s.exp!=0)
    s.exp->accept(*this) ;
  if(s.Terminal!=0) 
    s.Terminal->accept(*this) ;
}
void AST_visitor::visit(AST_Block &s) {
  for(auto ii=s.elements.begin();ii!=s.elements.end();++ii)
    if(*ii!=0)
      (*ii)->accept(*this) ;
}

void AST_visitor::visit(AST_typeSpec &s) {
  for(auto ii=s.type_spec.begin();ii!=s.type_spec.end();++ii)
    if(*ii != 0)
      (*ii)->accept(*this) ;
}

void AST_visitor::visit(AST_declaration &s) {
  for(auto ii=s.type_decl.begin();ii!=s.type_decl.end();++ii)
    if(*ii != 0)
      (*ii)->accept(*this) ;

  for(auto ii=s.decls.begin();ii!=s.decls.end();++ii)
    if(*ii != 0)
      (*ii)->accept(*this) ;
}
void AST_visitor::visit(AST_exprOper &s) {
  for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
    if(*ii != 0)
      (*ii)->accept(*this) ;

}
void AST_visitor::visit(AST_controlStatement &s) {
  s.controlType->accept(*this) ;
  for(auto ii=s.parts.begin();ii!=s.parts.end();++ii) {
    if(*ii != 0)
      (*ii)->accept(*this) ;
  }
}

void AST_errorCheck::visit(AST_syntaxError &s) {
  cerr << s.fileName << ":" << s.lineno << ": syntax, " <<s.error << endl ;
  error_count++ ;
}  

void AST_collectAccessInfo::visit(AST_Token &s) {
  if(ASTEqual(s,AST_type::TK_LOCI_VARIABLE)) {
    Loci::variable v(s.text) ;
    vmap_info vmap ;
    vmap.var += v ;
    accessed.insert(vmap) ;
    id2var[s.id] = v ;
  }
}
AST_type::ASTP getArrayVar(AST_type::ASTP expr) {
  while(ASTEqual(expr,AST_type::OP_ARRAY)) {
    CPTR<AST_exprOper> p(expr) ;
    expr = AST_type::ASTP(p->terms.front()) ;
  }
  return expr ;
}
bool isExprLociVariable(AST_type::ASTP expr) {
  expr = getArrayVar(expr) ;
  
  return ASTEqual(expr,AST_type::TK_LOCI_VARIABLE) ;
}

  
void AST_collectAccessInfo::visit(AST_exprOper &s) {
  switch(s.nodeType) {
  case AST_type::OP_ASSIGN:
  case AST_type::OP_TIMES_ASSIGN:
  case AST_type::OP_DIVIDE_ASSIGN:
  case AST_type::OP_MODULUS_ASSIGN:
  case AST_type::OP_PLUS_ASSIGN:
  case AST_type::OP_MINUS_ASSIGN:
  case AST_type::OP_SHIFT_LEFT_ASSIGN:
  case AST_type::OP_SHIFT_RIGHT_ASSIGN:
  case AST_type::OP_AND_ASSIGN:
  case AST_type::OP_OR_ASSIGN:
  case AST_type::OP_EXOR_ASSIGN:
    {
      auto ii=s.terms.begin();
      AST_collectAccessInfo first ;
      if(ii!=s.terms.end() && *ii != 0)
	(*ii)->accept(first) ;
      for(auto fi=first.accessed.begin();fi!= first.accessed.end();++fi)
	writes.insert(*fi) ;
      for(auto mi=first.id2var.begin();mi!=first.id2var.end();++mi)
	id2var[mi->first] = mi->second ;
      for(auto mi=first.id2vmap.begin();mi!=first.id2vmap.end();++mi)
	id2vmap[mi->first] = mi->second ;
      for(++ii;ii!=s.terms.end();++ii)
	if(*ii != 0) 
          (*ii)->accept(*this) ;
    }
    break ;
  case AST_type::OP_ARROW:
    {
      // Work on this to fill in the maps
      bool allLociVars = true ;
      for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
	if(*ii != 0) 
	  allLociVars = allLociVars && isExprLociVariable(*ii) ;
      if(allLociVars) {
	Loci::vmap_info vm ;
	
	for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
	  if(*ii != 0) {
	    AST_type::ASTP vp = getArrayVar(*ii) ;
	    CPTR<AST_Token> tok(vp) ;

	    Loci::variable v(tok->text) ;

	    variableSet vset ;
	    vset += v ;
	    if(ii+1 == s.terms.end()) {
	      vm.var += vset ;
	    } else {
	      vm.mapping.push_back(vset) ;
	    }
	  }
	id2vmap[s.id] = vm ;
	accessed.insert(vm) ;
	AST_collectAccessInfo base ;
	for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
	  if(*ii != 0)
	    (*ii)->accept(base) ;
	for(auto mi=base.id2var.begin();mi!=base.id2var.end();++mi)
	  id2var[mi->first] = mi->second ;
	for(auto mi=base.id2vmap.begin();mi!=base.id2vmap.end();++mi)
	  id2vmap[mi->first] = mi->second ;

      } else {
	for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
	  if(*ii != 0)
	    (*ii)->accept(*this) ;
      }
    }
    break ;
  default:
    for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
      if(*ii != 0)
	(*ii)->accept(*this) ;
    break ;
  }
}
  
void AST_simplePrint::visit(AST_exprOper &s) {
  auto t = id2rename.find(s.id) ;
  if(t != id2rename.end()) {
    out << t->second ;
    return ;
  }
  switch (s.nodeType) {
  case AST_type::OP_GROUP:
    out << '(' ;
    for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
      if(*ii != 0)
	(*ii)->accept(*this) ;
    out << ')' ;
    break ;
  case AST_type::OP_FUNC:
    {
      AST_type::ASTList::iterator ii=s.terms.begin() ;
      FATAL(ii == s.terms.end()) ;
      (*ii)->accept(*this) ;
      ++ii ;
      out << '(' ;
      if(ii != s.terms.end()) {
	(*ii)->accept(*this) ;
        ++ii ;
      }
      out << ')' ;
      if(ii!=s.terms.end()) {
	cerr << "internal error processing func" ;
        (*ii)->accept(*this) ;
      }
    }
    break ;
  case AST_type::OP_TEMPLATE:
    {
      AST_type::ASTList::iterator ii=s.terms.begin() ;
      FATAL(ii == s.terms.end()) ;
      (*ii)->accept(*this) ;
      ++ii ;
      out << '<' ;
      if(ii != s.terms.end()) {
	(*ii)->accept(*this) ;
        ++ii ;
      }
      out << '>' ;
      if(ii!=s.terms.end()) {
	cerr << "internal error processing func" ;
        (*ii)->accept(*this) ;
      }
    }
    break ;
  case AST_type::OP_ARRAY:
    {
      AST_type::ASTList::iterator ii=s.terms.begin() ;
      if(*ii != 0)
	(*ii)->accept(*this) ;
      ++ii ;
      out << '[' ;
      if(*ii != 0)
	(*ii)->accept(*this) ;
      out << ']' ;
      ++ii ;
      if(ii!=s.terms.end()) {
	cerr << "internal error processing array" ;
	if(*ii != 0)
	  (*ii)->accept(*this);
      }
    }
    break ;
  case AST_type::OP_TERNARY:
    {
      AST_type::ASTList::iterator ii=s.terms.begin() ;
      if(*ii != 0)
	(*ii)->accept(*this) ;
      ++ii ;
      out << '?' ;
      if(*ii != 0)
	(*ii)->accept(*this) ;
      out << ':' ;
      if(ii==s.terms.end()) {
	cerr << "internal error on tertiary operator" << endl ;
      } else
	++ii ;
      if(*ii != 0)
	(*ii)->accept(*this) ;

    }
    break ;
      
  case AST_type::OP_UNARY_PLUS:
  case AST_type::OP_UNARY_MINUS:
  case AST_type::OP_NOT:
  case AST_type::OP_AMPERSAND:
  case AST_type::OP_STAR:
  case AST_type::OP_INCREMENT:
  case AST_type::OP_DECREMENT:
    {
      string op = OPtoString(s) ;
      out << op ;
      for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
	if(*ii != 0)
	  (*ii)->accept(*this) ;
    }
    break ;
  case AST_type::OP_POSTINCREMENT:
  case AST_type::OP_POSTDECREMENT:
    {
      string op = OPtoString(s) ;
      for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
	if(*ii != 0)
	  (*ii)->accept(*this) ;
      out << op ;
    }
    break ;
  default:
    {
      string op = OPtoString(s) ;
      for(auto ii=s.terms.begin();ii!=s.terms.end();) {
	if(*ii != 0)
	  (*ii)->accept(*this) ;
	++ii ;
	if(ii != s.terms.end())
	  out << op ;
      }
    }
    break ;
  }
}

void AST_simplePrint::visit(AST_Token &s) {
  if(s.lineno >0 && lineno != s.lineno) {
    out << endl ;
    if(!prettyPrint)
      if(lineno < 0 || lineno+1 != s.lineno)
	out << "#line " << s.lineno << endl ;

    lineno = s.lineno ;
  }
  auto t = id2rename.find(s.id) ;
  if(t != id2rename.end()) {
    out << t->second ;
  } else {
    if(ASTEqual(s,AST_type::TK_LOCI_DIRECTIVE)) {
      out << "$[" << s.text << "]" ;
    } else if(ASTEqual(s,AST_type::TK_LOCI_CONTAINER)) {
      out << "$*" << s.text ;
    } else if(ASTEqual(s,AST_type::TK_LOCI_VARIABLE)) {
      out << "$" << s.text ;
    } else 
      out <<s.text << ' ' ;
  }
}

bool condenseOp(AST_type::elementType et) {
  switch(et) {
  case AST_type::OP_ARROW:
  case AST_type::OP_SCOPE:
  case AST_type::OP_TIMES:
  case AST_type::OP_PLUS:
  case AST_type::OP_COMMA:
    return true ;
  default:
    return false ;
  }
}
    
void AST_condenseLeftAssociative::visit(AST_exprOper &e) {
  AST_type::elementType op = e.nodeType ;

  if(condenseOp(op) && e.terms[0]->nodeType == op) {
    // push in reverse order e.terms onto terms.
    AST_type::ASTList terms = e.terms ;
    std::reverse(terms.begin(),terms.end()) ;
    // keep expanding last term until we no longer can
    while(terms.back()->nodeType == op) {
      CPTR<AST_exprOper> p = CPTR<AST_exprOper>(terms.back()) ;
      terms.pop_back() ;
      for(auto ii=p->terms.rbegin();ii!=p->terms.rend();++ii)
        terms.push_back(*ii) ;
    }
    // reverse order back to normal
    std::reverse(terms.begin(),terms.end()) ;
    e.terms.swap(terms) ;
  }
  for(auto ii=e.terms.begin();ii!=e.terms.end();++ii)
    if(*ii != 0)
      (*ii)->accept(*this) ;
}
