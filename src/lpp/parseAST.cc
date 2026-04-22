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
using Loci::vmap_info ;
using Loci::variableSet ;
using namespace nodeTypes ;

/// Convert operator code to corresponding string for printing
string OPtoString(AST_type::elementType val) {
  switch(val) {
  case OP_SCOPE:
    return string("::") ;
  case OP_AT:
    return string("@") ;
  case OP_ARROW:
    return string("->") ;
  case OP_TIMES:
    return string("*") ;
  case OP_DIVIDE:
    return string("/") ;
  case OP_MODULUS:
    return string("%") ;
  case OP_PLUS:
    return string("+") ;
  case OP_MINUS:
    return string("-") ;
  case OP_SHIFT_RIGHT:
    return string(">>") ;
  case OP_SHIFT_LEFT:
    return string("<<") ;
  case OP_LT:
    return string("<") ;
  case OP_GT:
    return string(">") ;
  case OP_GE:
    return string(">=") ;
  case OP_LE:
    return string("<=") ;
  case OP_EQUAL:
    return string("==") ;
  case OP_NOT_EQUAL:
    return string("!=") ;
  case OP_AND:
    return string("&") ;
  case OP_EXOR:
    return string("^") ;
  case OP_OR:
    return string("|") ;
  case OP_LOGICAL_AND:
    return string("&&") ;
  case OP_LOGICAL_OR:
    return string("||") ;
  case OP_ASSIGN:
    return string("=") ;
  case OP_TIMES_ASSIGN:
    return string("*=") ;
  case OP_DIVIDE_ASSIGN:
    return string("/=") ;
  case OP_MODULUS_ASSIGN:
    return string("%=") ;
  case OP_PLUS_ASSIGN:
    return string("+=") ;
  case OP_MINUS_ASSIGN:
    return string("-=") ;
  case OP_SHIFT_LEFT_ASSIGN:
    return string("<<=") ;
  case OP_SHIFT_RIGHT_ASSIGN:
    return string(">>=") ;
  case OP_AND_ASSIGN:
    return string("&=") ;
  case OP_OR_ASSIGN:
    return string("|=") ;
  case OP_EXOR_ASSIGN:
    return string("^=") ;
  case OP_COMMA:
    return string(",") ;
  case OP_DOT:
    return string(".") ;
  case OP_COLON:
    return string(":") ;
  case OP_SEMICOLON:
    return string(";") ;
  case OP_INCREMENT:
    return string(" ++") ;
  case OP_DECREMENT:
    return string(" --") ;
  case OP_POSTINCREMENT:
    return string("++ ") ;
  case OP_POSTDECREMENT:
    return string("-- ") ;
  case OP_UNARY_PLUS:
    return string("+") ;
  case OP_UNARY_MINUS:
    return string("-") ;
  case OP_NOT:
    return string("!") ;
  case OP_TILDE:
    return string("~") ;
  case OP_AMPERSAND:
    return string("&") ;
  case OP_TERNARY:
    return string("?") ;
  case OP_DOLLAR:
    return string("$") ;
  case OP_STAR:
    return string("*") ;
  case OP_CAST:
    return string(" ") ;
  case OP_TEMPLATE_CAST:
    return string(" ") ;
  case OP_BRACEBLOCK:
    return string(" ") ;
  default:
    return string("/*error*/") ;
  }
  return string("/*error*/") ;
}

// Parse a type specification used in either a type declaration statement
// or a template argument.
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
    cerr << "in parseTypeSpecifier, token = " << token->text
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    switch(token->nodeType) {
    case TK_CONST:
      type_spec.push_back(AST_type::ASTP(token)) ;
#ifdef VERBOSE
      cerr << "in parseTypeSpecifier const variable"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      break ;
    case TK_EXTERN:
      type_spec.push_back(AST_type::ASTP(token)) ;
#ifdef VERBOSE
      cerr << "in parseTypeSpecifier extern variable"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      break ;
    case TK_CHAR:
    case TK_FLOAT:
    case TK_DOUBLE:
    case TK_INT:
    case TK_VOID:
    case TK_BOOL:
    case TK_SHORT:
    case TK_LONG:
    case TK_SIGNED:
    case TK_UNSIGNED:
    case TK_AUTO:
#ifdef VERBOSE
      cerr << "in parseTypeSpecifier, builtin type"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
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
    case TK_SCOPE:
    case TK_NAME:
#ifdef VERBOSE
      cerr << "in parseTypeSpecifier, named type"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      if(builtin_type || defined_type ) {
        istype = false ;
      } else {
        pushToken(token) ;

        AST_type::ASTP typedec =
          parseIdentifier(is,linecount,fileName,typemap) ;
        if(typedec->nodeType == OP_FUNC) {
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
      cerr << "in parseTypeSpecifier, finished type parsing"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      istype = false ;
      break ;
    }
    if(istype)
      token = getToken(is,linecount) ;
  } while(istype) ;

  pushToken(token) ;
  CPTR<AST_typeSpec> p = new AST_typeSpec ;
  p->nodeType = ND_TYPE_SPEC ;
  p->type_spec = type_spec ;
  return AST_type::ASTP(p) ;
}

// AST_type constructor enumerates all AST_type objects to give them a unique
// identifier (id)
static int parseIDctr = 0 ;
AST_type::AST_type() {
  nodeType = OP_ERROR ;
  id = parseIDctr++ ;
}

/// Acceptor method AST node, passes node to visitor object for AST_token
void AST_Token::accept(AST_visitor &v) {  v.visit(*this) ; }

/// Code to clone an AST_Token
AST_type::ASTP AST_Token::clone() const {
  CPTR<AST_Token> tok = new AST_Token ;
  tok->nodeType = nodeType ;
  tok->text = text ;
  tok->lineno = lineno ;
  return ASTP(tok) ;
}

/// Acceptor method AST node, passes node to visitor object of
/// AST_syntaxError.
void AST_syntaxError::accept(AST_visitor &v) {  v.visit(*this) ; }

/// Code to clone an AST_syntaxError
AST_type::ASTP AST_syntaxError::clone() const {
  CPTR<AST_syntaxError> p = new AST_syntaxError(error,lineno,fileName) ;
  p->nodeType = nodeType ;
  return ASTP(p) ;
}

/// Acceptor method AST node, passes node to visitor object for
/// AST_simpleStatement
void AST_SimpleStatement::accept(AST_visitor &v) {  v.visit(*this) ; }

/// Code to clone an AST_simpleStatement
AST_type::ASTP AST_SimpleStatement::clone() const {
  CPTR<AST_SimpleStatement> p = new AST_SimpleStatement(exp->clone(),
                                                        Terminal->clone()) ;
  p->nodeType = nodeType ;
  return ASTP(p) ;
}

/// Acceptor method AST node, passes node to visitor object for AST_Block
void AST_Block::accept(AST_visitor &v) {  v.visit(*this) ; }

/// General helper code to clone an ASTList
inline void cloneList(AST_type::ASTList &l, const AST_type::ASTList &in) {
  AST_type::ASTList tmp(in.size()) ;
  for(size_t i=0;i<in.size();++i)
    tmp[i] = in[i]->clone() ;
  l.swap(tmp) ;
}

// Code to clone the codeblock AST_Block
AST_type::ASTP AST_Block::clone() const {
  CPTR<AST_Block> p = new AST_Block ;
  p->nodeType = nodeType ;
  cloneList(p->elements,elements) ;
  p->identifiers = identifiers ;
  return ASTP(p) ;
}

/// Acceptor method AST node, passes node to visitor object for
/// AST_typeSpec node.
void AST_typeSpec::accept(AST_visitor &v) {  v.visit(*this) ; }

/// Code to clone AST_typeSpec node.
AST_type::ASTP AST_typeSpec::clone() const {
  CPTR<AST_typeSpec> p = new AST_typeSpec ;
  p->nodeType = nodeType ;
  cloneList(p->type_spec,type_spec) ;
  return ASTP(p) ;
}

/// Acceptor method AST node, passes node to visitor object for
/// AST_declaration node
void AST_declaration::accept(AST_visitor &v) {  v.visit(*this) ; }

/// Code to clone AST_declaration node
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

///Code to clone an expression operator node.  (AST_exprOper)
AST_type::ASTP AST_exprOper::clone() const {
  CPTR<AST_exprOper> p = new AST_exprOper ;
  p->nodeType = nodeType ;
  cloneList(p->terms,terms) ;
  return ASTP(p) ;
}

/// Acceptor method AST node, passes node to visitor object for
/// AST_controlStatement
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
  case TK_STRING:
  case TK_NAME:
  case TK_NUMBER:
  case TK_TRUE:
  case TK_FALSE:
  case TK_LOCI_VARIABLE:
  case TK_LOCI_CONTAINER:
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
  cerr << "in parseTerm, token = " << token->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  if(isTerm(token))
    return AST_type::ASTP(token) ;
  else {
    pushToken(token) ;
    AST_type::ASTP(new AST_syntaxError("expecting token",token->lineno,"")) ;
    return AST_type::ASTP(0) ;
  }
}


AST_type::elementType tokenToOperator(AST_type::elementType nodeType) {
  switch(nodeType) {
  case TK_SCOPE:
    return OP_SCOPE ;
    //-----------------------------------
    // For using @ to separate namespaces
  case TK_AT:
    return OP_AT ;
    //-----------------------------------
    // Traditional C operators
  case TK_ARROW:
    return OP_ARROW ;
  case TK_TIMES:
    return OP_TIMES ;
  case TK_DIVIDE:
    return OP_DIVIDE ;
  case TK_MODULUS:
    return OP_MODULUS ;
  case TK_PLUS:
    return OP_PLUS ;
  case TK_MINUS:
    return OP_MINUS ;
  case TK_SHIFT_RIGHT:
    return OP_SHIFT_RIGHT ;
  case TK_SHIFT_LEFT:
    return OP_SHIFT_LEFT ;
  case TK_LT:
    return OP_LT ;
  case TK_GT:
    return OP_GT ;
  case TK_GE:
    return OP_GE ;
  case TK_LE:
    return OP_LE ;
  case TK_EQUAL:
    return OP_EQUAL ;
  case TK_NOT_EQUAL:
    return OP_NOT_EQUAL ;
  case TK_AND:
    return OP_AND ;
  case TK_EXOR:
    return OP_EXOR ;
  case TK_OR:
    return OP_OR ;
  case TK_LOGICAL_AND:
    return OP_LOGICAL_AND;
  case TK_LOGICAL_OR:
    return OP_LOGICAL_OR ;
  case TK_ASSIGN:
    return OP_ASSIGN ;
  case TK_TIMES_ASSIGN:
    return OP_TIMES_ASSIGN ;
  case TK_DIVIDE_ASSIGN:
    return OP_DIVIDE_ASSIGN ;
  case TK_MODULUS_ASSIGN:
    return OP_MODULUS_ASSIGN ;
  case TK_PLUS_ASSIGN:
    return OP_PLUS_ASSIGN ;
  case TK_MINUS_ASSIGN:
    return OP_MINUS_ASSIGN ;
  case TK_SHIFT_LEFT_ASSIGN:
    return OP_SHIFT_LEFT_ASSIGN ;
  case TK_SHIFT_RIGHT_ASSIGN:
    return OP_SHIFT_RIGHT_ASSIGN ;
  case TK_AND_ASSIGN:
    return OP_AND_ASSIGN ;
  case TK_OR_ASSIGN:
    return OP_OR_ASSIGN ;
  case TK_EXOR_ASSIGN:
    return OP_EXOR_ASSIGN ;
  case TK_COMMA:
    return OP_COMMA ;
  case TK_QUESTION:
    return OP_TERNARY ;
  case TK_COLON:
    return OP_COLON ; // works with OP_TERNARY
  case TK_DOT:
    return OP_DOT ;
  default:
    return OP_ERROR ;
  }
  return OP_ERROR ;
}
/// Parse operator token, convert from TK_* node type to OP_* nodetype
AST_type::ASTP parseOperator(std::istream &is, int &linecount) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseOperator, token = " << openToken->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  AST_type::elementType op = tokenToOperator(openToken->nodeType) ;
  if(op == OP_ERROR) {
    pushToken(openToken) ;
    return AST_type::ASTP(0) ;
  }
  openToken->nodeType = op ;
  return AST_type::ASTP(openToken) ;
}

/// Check for input token that is a unary operator
inline bool checkUnaryToken(AST_type::elementType e) {
  switch(e) {
  case TK_PLUS:
  case TK_MINUS:
  case TK_NOT:
  case TK_AND:
  case TK_TIMES:
  case TK_INCREMENT:
  case TK_DECREMENT:
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
  case TK_PLUS:
    return OP_UNARY_PLUS ;
  case TK_MINUS:
    return OP_UNARY_MINUS ;
  case TK_NOT:
    return OP_NOT ;
  case TK_AND:
    return OP_AMPERSAND ;
  case TK_TIMES:
    return OP_STAR ;
  case TK_INCREMENT:
    return OP_INCREMENT ;
  case TK_DECREMENT:
    return OP_DECREMENT ;
  default:
    return OP_ERROR ;
  }
}

/// Check for postfix token '++' or '--'
bool checkPostFixToken(AST_type::elementType e) {
  return (e == TK_INCREMENT ||
	  e == TK_DECREMENT) ;
}

/// Pares postfix operator convert node type from TK_* to appropriate OP_*
AST_type::elementType postFixOperator(AST_type::elementType e) {
  switch(e) {
  case TK_INCREMENT:
    return OP_POSTINCREMENT ;
  case TK_DECREMENT:
    return OP_POSTDECREMENT ;
  default:
    return OP_ERROR ;
  }
  return OP_ERROR ;
}

/// Apply postfix operator to previous expression passed in expr
AST_type::ASTP applyPostFixOperator(AST_type::ASTP expr,
				    std::istream &is, int &linecount,
				    string fileName,
				    varmap &typemap) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in applyPostFixOperator, token = " << openToken->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  if(checkPostFixToken(openToken->nodeType)) {
    CPTR<AST_exprOper> post = new AST_exprOper ;
    post->nodeType = postFixOperator(openToken->nodeType) ;
    post->terms.push_back(expr) ;
    openToken = getToken(is,linecount) ;
    pushToken(openToken) ;
    if(checkPostFixToken(openToken->nodeType) ||
       ASTEqual(openToken,TK_OPENBRACKET))
      return applyPostFixOperator(AST_type::ASTP(post),is,linecount,
				  fileName,typemap) ;
    else
      return AST_type::ASTP(post) ;
  }
  if(ASTEqual(openToken,TK_OPENPAREN)) {
    return parseFunctionArguments(expr,is,linecount,fileName,typemap) ;
  }
  if(ASTEqual(openToken,TK_OPENBRACKET)) {
    AST_type::ASTP index = parseExpression(is,linecount,fileName,typemap) ;
    openToken = getToken(is,linecount) ;
    if(openToken->nodeType != TK_CLOSEBRACKET) {
      pushToken(openToken) ;
      return AST_type::ASTP(new AST_syntaxError("expecting ']'",openToken->lineno,fileName)) ;
    }
    CPTR<AST_exprOper> array = new AST_exprOper ;
    array->nodeType = OP_ARRAY ;
    array->terms.push_back(expr) ;
    array->terms.push_back(index) ;

    openToken = getToken(is,linecount) ;
    pushToken(openToken) ;
    if(checkPostFixToken(openToken->nodeType) ||
       ASTEqual(openToken,TK_OPENBRACKET))
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
  if(ASTEqual(openToken,TK_NAME))
    cerr << " name=" << openToken->text ;
  cerr << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  AST_type::ASTList terms ;
  // Special case for global scope (eg ::name)
  if(ASTEqual(openToken,TK_SCOPE)) {
    CPTR<AST_Token> tmp = new AST_Token ;
    tmp->text = "" ;
    tmp->lineno = openToken->lineno ;
    tmp->nodeType = TK_NIL ;
    terms.push_back(AST_type::ASTP(tmp)) ;
    openToken = getToken(is,linecount) ;
  }
  // Now openToken should be a name, so push it onto the terms list
  if(!ASTEqual(openToken,TK_NAME)) {
    terms.push_back(AST_type::ASTP(new AST_syntaxError("expecting name following '::'",openToken->lineno,fileName))) ;
  } else
    terms.push_back(AST_type::ASTP(openToken)) ;

  // Now get more terms as long as we have a scope operator
  openToken = getToken(is,linecount) ;

#ifdef VERBOSE
  cerr << "searching token = " << openToken->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  while(ASTEqual(openToken,TK_SCOPE)) {
    openToken = getToken(is,linecount) ;
    if(!ASTEqual(openToken,TK_NAME)) {
      terms.push_back(AST_type::ASTP(new AST_syntaxError("expecting name following '::'",openToken->lineno,fileName))) ;
    } else {
      terms.push_back(AST_type::ASTP(openToken)) ;
    }
    // Get next token, continue chain if is a scope token
    openToken = getToken(is,linecount) ;
  }
  // We are done parsing, so return this to the token input
  pushToken(openToken) ;

  if(terms.size() > 1) {
    CPTR<AST_exprOper> tmp = new AST_exprOper ;
    tmp->nodeType = OP_SCOPE ;
    tmp->terms = terms ;
    return AST_type::ASTP(tmp) ;
  }
  return AST_type::ASTP(terms.back()) ;
}

AST_type::ASTP parseTemplateArguments(AST_type::ASTP objname,
                                      std::istream &is, int &linecount,
                                      const string &fileName,
                                      varmap &typemap) {
  // Template Arguments
#ifdef VERBOSE
  cerr << "parseTemplateArguments: found term template open '<'"
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  CPTR<AST_Token> token = getToken(is,linecount) ;
  AST_type::ASTP arg = AST_type::ASTP(token) ;
  if(!ASTEqual(token,TK_NUMBER)) {
    pushToken(token) ;
    arg = parseTypeSpecifier(is,linecount,fileName,typemap) ;
  }

#ifdef VERBOSE
  cerr << "parseTemplateArguments: parsed type specification "
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  CPTR<AST_Token> closeToken = getToken(is,linecount) ;
  if(ASTEqual(closeToken,TK_CLOSETEMPLATE)) {
#ifdef VERBOSE
    cerr << "parseTemplateArguments: only single argument "
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    // only a single type in template parameters
    CPTR<AST_exprOper> func = new AST_exprOper ;
    func->nodeType = OP_TEMPLATE ;
    func->terms.push_back(objname) ;
    func->terms.push_back(arg) ;
    return AST_type::ASTP(func) ;
  }
  if(ASTEqual(closeToken,TK_COMMA)) {

#ifdef VERBOSE
    cerr << "parseTemplateArguments: parsing comma separated list of types "
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif

    CPTR<AST_exprOper> tlist = new AST_exprOper ;
    tlist->nodeType = OP_COMMA ;    tlist->terms.push_back(arg) ;
    while(ASTEqual(closeToken,TK_COMMA)) {
#ifdef VERBOSE
      cerr << "parseTemplateArgument, got " << closeToken->text
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      CPTR<AST_Token> token = getToken(is,linecount) ;
      AST_type::ASTP arg = AST_type::ASTP(token) ;
      if(!ASTEqual(token,TK_NUMBER)) {
        pushToken(token) ;
        arg = parseTypeSpecifier(is,linecount,fileName,typemap) ;
      }

      tlist->terms.push_back(arg) ;
      closeToken = getToken(is,linecount) ;
    }
    if(!ASTEqual(closeToken,TK_CLOSETEMPLATE)) {
      tlist->terms.push_back(AST_type::ASTP(new AST_syntaxError("expecting '>'",closeToken->lineno,fileName))) ;
    }
    CPTR<AST_exprOper> func = new AST_exprOper ;
    func->nodeType = OP_TEMPLATE ;
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
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  AST_type::ASTP args = parseExpression(is,linecount,fileName,typemap) ;
  CPTR<AST_Token> closeToken = getToken(is,linecount) ;
  CPTR<AST_exprOper> func = new AST_exprOper ;
  func->nodeType = OP_FUNC ;
  func->terms.push_back(AST_type::ASTP(objname)) ;
  if(args != 0)
    func->terms.push_back(args) ;
  if(!ASTEqual(closeToken,TK_CLOSEPAREN)) {
    CPTR<AST_Token> err = new AST_Token ;
    *err = *closeToken ;
    err->nodeType = OP_ERROR ;
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
  cerr << "in parseIdentifier"
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif

  AST_type::ASTP objname = parseScopedObject(is,linecount,fileName,typemap) ;

  CPTR<AST_Token> openToken = getToken(is,linecount) ;
  while(ASTEqual(openToken,TK_OPENTEMPLATE)) {
    objname = parseTemplateArguments(objname,is,linecount,fileName,typemap) ;
    openToken = getToken(is,linecount) ;
    if(!ASTEqual(openToken,TK_SCOPE))
      break ;
    CPTR<AST_exprOper> tmp = new AST_exprOper ;
    tmp->nodeType = OP_SCOPE ;
    tmp->terms.push_back(objname) ;
    openToken = getToken(is,linecount) ;
    if(!ASTEqual(openToken,TK_NAME)) {
      AST_type::ASTP p =
        AST_type::ASTP(new AST_syntaxError("error parsing scope after template type",openToken->lineno,fileName)) ;
      tmp->terms.push_back(AST_type::ASTP(openToken)) ;
      tmp->terms.push_back(p) ;
      return AST_type::ASTP(tmp) ;
    }
    tmp->terms.push_back(AST_type::ASTP(openToken)) ;
    objname = AST_type::ASTP(tmp) ;
    openToken = getToken(is,linecount) ;
  }
  pushToken(openToken) ;

  openToken = getToken(is,linecount) ;
  if(ASTEqual(openToken,TK_OPENPAREN)) {
    objname = parseFunctionArguments(objname,is,linecount,fileName,typemap) ;
  } else
    pushToken(openToken) ;
  return objname ;
}

// Scan forward in tokens to determine if an openparen forms a C style
// type cast operator.  This will scan forward skipping over names,
// builtin types, or pointer/reference identifiers until a close paren
// is found, then scan after that to determine if the format is consistent
// with this being a type operator (e.g. a name or similar following).  This
// will indicate that this is a C style casting operator.
bool scanForCStyleCast(std::istream &is, int &linecount,varmap &typemap) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
  if(!ASTEqual(openToken,TK_OPENPAREN)) {
    pushToken(openToken) ;
    return false ;
  }
  vector<CPTR<AST_Token>> tokenStack ;
  // Now scan forward looking for the close paren
  do {
    tokenStack.push_back(openToken) ;
    openToken = getToken(is,linecount) ;
    // If we find a close paren, then break out of loop
    // to check for
    if(ASTEqual(openToken,TK_CLOSEPAREN) &&
       tokenStack.size() > 1)
      break ;
    bool builtintype = false ;
    bool usertype = false ;
    switch(openToken->nodeType) {
    case TK_CHAR:
    case TK_FLOAT:
    case TK_DOUBLE:
    case TK_INT:
    case TK_BOOL:
    case TK_SHORT:
    case TK_LONG:
    case TK_SIGNED:
    case TK_UNSIGNED:
    case TK_VOID:
      builtintype = true ;
      if(usertype) { // If user type is already defined, then this
        // is not a well formed type definition, so this is not a C style
        // cast.
        pushToken(openToken) ;
        while(tokenStack.size() > 0) {
          pushToken(tokenStack.back()) ;
          tokenStack.pop_back() ;
        }
        return false ;
      }
      continue ;
    case TK_CONST:
      break ;
    case TK_SCOPE:
    case TK_OPENTEMPLATE:
    case TK_CLOSETEMPLATE:
    case TK_NAME:
      usertype = true ;
      if(builtintype) { // Can't mix builtin types and named types
        return false ;
      }
      // Probably should build identifier here, and check if it is a
      // local variable.  Future work.
      continue ;
    case TK_AMPERSAND:
    case TK_TIMES:
      // These need to be at the very end e.g. (real *) or (real * const)
      tokenStack.push_back(openToken) ;
      openToken = getToken(is,linecount) ;
      if(ASTEqual(openToken,TK_CONST)) {
        tokenStack.push_back(openToken) ;
        openToken = getToken(is,linecount) ;
      }
      // We have enough information to return
      pushToken(openToken) ;
      while(tokenStack.size() > 0) {
        pushToken(tokenStack.back()) ;
        tokenStack.pop_back() ;
      }
      return ASTEqual(openToken,TK_CLOSEPAREN) ;

    default:
      pushToken(openToken) ;
      while(tokenStack.size() > 0) {
        pushToken(tokenStack.back()) ;
        tokenStack.pop_back() ;
      }
      return false ;
    }
  } while(true) ;
  // found closing paren, search what follows to determine if this is a
  // C style cast.
  tokenStack.push_back(openToken) ;
  openToken = getToken(is,linecount) ;
  bool castop=false ;
  switch(openToken->nodeType) {
  case TK_OPENPAREN:
  case TK_NAME:
  case TK_LOCI_VARIABLE:
  case TK_LOCI_CONTAINER:
  case TK_NUMBER:
  case TK_STRING:
  case TK_TRUE:
  case TK_FALSE:
  case TK_PLUS:
  case TK_MINUS:
  case TK_NOT:
  case TK_TILDE:
  case TK_STAR:
  case TK_INCREMENT:
  case TK_DECREMENT:
    castop = true ;
  default:
    break ;
  }

  pushToken(openToken) ;
  while(tokenStack.size() > 0) {
    pushToken(tokenStack.back()) ;
    tokenStack.pop_back() ;
  }
  return castop ;
}

/// Parse the part of an expression that is not a binary or ternary operator
AST_type::ASTP parseExpressionPartial(std::istream &is, int &linecount,
				      const string &fileName,
				      varmap &typemap,
                                      AST_type::operatorPrecedence prec) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseExpressionPartial, token = " << openToken->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif

  // Check for a parenthesized group 
  if(ASTEqual(openToken,TK_OPENPAREN)) {
    pushToken(openToken) ;
    // Grab '('
    openToken = getToken(is, linecount) ;
#ifdef VERBOSE
    cerr << "parseExpressionPartial: Found '(' parsing new expression"
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    // parse expression contained within parenthesized group
    AST_type::ASTP exp = parseExpression(is,linecount,fileName,typemap) ;
#ifdef VERBOSE
    cerr << "parseExpressionPartial: creating group"
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
    {
      AST_simplePrint printer(cerr,-2) ;
      exp->accept(printer) ;
      cerr << endl ;
    }
#endif

    CPTR<AST_Token> closeToken = getToken(is,linecount) ;
    CPTR<AST_exprOper> group = new AST_exprOper ;
    group->nodeType = OP_GROUP ;
    group->terms.push_back(AST_type::ASTP(exp)) ;
    if(closeToken->nodeType != TK_CLOSEPAREN) {
      pushToken(closeToken) ;
      return AST_type::ASTP(new AST_syntaxError("expecting ')'",closeToken->lineno,fileName)) ;

    }
    return AST_type::ASTP(group) ;
  }
  if(ASTEqual(openToken,TK_CONST_CAST) ||
     ASTEqual(openToken,TK_DYNAMIC_CAST) ||
     ASTEqual(openToken,TK_REINTERPRET_CAST) ||
     ASTEqual(openToken,TK_STATIC_CAST)) {
#ifdef VERBOSE
    cerr << "in template cast, token = " << openToken->text
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    // parse templated cast
    CPTR<AST_exprOper> cast = new AST_exprOper ;
    cast->nodeType = OP_TEMPLATE_CAST ;
    cast->terms.push_back(AST_type::ASTP(openToken)) ;
    // how scan for template brace
    openToken=getToken(is,linecount) ;
#ifdef VERBOSE
    cerr << "in template cast parser, token = " << openToken->text
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    if(!ASTEqual(openToken,TK_OPENTEMPLATE)) {
      cerr << "expecting '<'" << endl ;
      cast->terms.push_back(AST_type::ASTP(new AST_syntaxError("expecting '<'",openToken->lineno,fileName))) ;
      return AST_type::ASTP(cast) ;
    }
    cast->terms.push_back(AST_type::ASTP(openToken)) ;
    AST_type::ASTP objname =
      parseTypeSpecifier(is,linecount,fileName,typemap) ;
#ifdef VERBOSE
    cerr << "in template cast parser, got template arguments"
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    cast->terms.push_back(objname) ;
    openToken=getToken(is,linecount) ;
    if(ASTEqual(openToken,TK_TIMES) ||
       ASTEqual(openToken,TK_AMPERSAND)) {
      cast->terms.push_back(AST_type::ASTP(openToken)) ;
      openToken = getToken(is,linecount) ;
    }
    if(!ASTEqual(openToken,TK_CLOSETEMPLATE)) {
      cerr << "expecting '>'" << endl ;
      cast->terms.push_back(AST_type::ASTP(new AST_syntaxError("expecting '>'",openToken->lineno,fileName))) ;
      return AST_type::ASTP(cast) ;
    }
    cast->terms.push_back(AST_type::ASTP(openToken)) ;
    openToken=getToken(is,linecount) ;
    if(!ASTEqual(openToken,TK_OPENPAREN)) {
      cerr << "expecting '('" << endl ;
      cast->terms.push_back(AST_type::ASTP(new AST_syntaxError("expecting '('",openToken->lineno,fileName))) ;
      return AST_type::ASTP(cast) ;
    }
    cast->terms.push_back(AST_type::ASTP(openToken)) ;
    cast->terms.push_back(AST_type::ASTP(parseExpression(is,linecount,fileName,typemap))) ;
    openToken=getToken(is,linecount) ;
    if(!ASTEqual(openToken,TK_CLOSEPAREN)) {
      cerr << "expecting ')'" << endl ;
      cast->terms.push_back(AST_type::ASTP(new AST_syntaxError("expecting ')'",openToken->lineno,fileName))) ;
      return AST_type::ASTP(cast) ;
    }
    cast->terms.push_back(AST_type::ASTP(openToken)) ;
    return AST_type::ASTP(cast) ;
  }
  // Brace-enclosed initializer list: { expr, ... }  (commas inside do not end
  // the outer declaration; parseExpression stops before '}').
  if(ASTEqual(openToken,TK_OPENBRACE)) {
    CPTR<AST_Token> afterBrace = getToken(is,linecount) ;
    if(ASTEqual(afterBrace,TK_CLOSEBRACE)) {
      CPTR<AST_exprOper> br = new AST_exprOper ;
      br->nodeType = OP_BRACEBLOCK ;
      return AST_type::ASTP(br) ;
    }
    pushToken(afterBrace) ;
    AST_type::ASTP inner =
      parseExpression(is,linecount,fileName,typemap) ;
    CPTR<AST_Token> closeTok = getToken(is,linecount) ;
    if(closeTok->nodeType != TK_CLOSEBRACE) {
      pushToken(closeTok) ;
      return AST_type::ASTP(new AST_syntaxError(
                                                "expecting '}' in brace initializer", closeTok->lineno,fileName)) ;
    }
    CPTR<AST_exprOper> br = new AST_exprOper ;
    br->nodeType = OP_BRACEBLOCK ;
    br->terms.push_back(inner) ;
    return AST_type::ASTP(br) ;
  }
  // Check for unary operators
  if(checkUnaryToken(openToken)) {
#ifdef VERBOSE
    cerr << "parseExpressionPartial: found unary operator '"
         << openToken->text << " parsing unary target."
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    //check for unary operators
    AST_type::ASTP expr =
      parseExpression(is,linecount,fileName,typemap,
                      AST_type::operatorPrecedence::PREC_UNARY_PREFIX) ;
    if(expr!= 0) {
      CPTR<AST_exprOper> unary = new AST_exprOper ;
      unary->nodeType = unaryOperator(openToken->nodeType) ;
      unary->terms.push_back(expr) ;
#ifdef VERBOSE
      cerr << "parseExpressionPartial: creating unary operator, "
           << "search for postfix operators ++,--, or []"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      return AST_type::ASTP(unary) ;
    }
  }
  // Check if we are parsing a scoped name
  if(ASTEqual(openToken,TK_NAME) ||
     ASTEqual(openToken,TK_SCOPE)) {
    pushToken(openToken) ;
    AST_type::ASTP objname = parseIdentifier(is,linecount,fileName,typemap) ;
    return objname ;
  }

  if(isTerm(openToken)) {
#ifdef VERBOSE
    cerr << "parseExpressionPartial: found term '"
         << openToken->text << "', checking for following parenthesis"
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    pushToken(openToken) ;
    AST_type::ASTP exp = parseTerm(is,linecount);


    return exp ;
  }
  pushToken(openToken) ;
  return 0 ;
}

/// get precedence of operator code (larger = tighter binding).
inline AST_type::operatorPrecedence getPrecedence(AST_type::elementType op) {
  // Note going from highest precedence to lowest
  switch(op) {
  case OP_SCOPE:
  case OP_AT:
    return AST_type::operatorPrecedence::PREC_SCOPE ;
  case OP_POSTINCREMENT:
  case OP_POSTDECREMENT:
  case OP_ARRAY:
  case OP_DOT:
  case OP_ARROW:
    return AST_type::operatorPrecedence::PREC_UNARY_POSTFIX ;
  case OP_INCREMENT:
  case OP_DECREMENT:
  case OP_UNARY_PLUS:
  case OP_UNARY_MINUS:
  case OP_NOT:
  case OP_TILDE:
  case OP_STAR:
  case OP_AMPERSAND:
  case OP_DOLLAR:
    return AST_type::operatorPrecedence::PREC_UNARY_PREFIX ;
    // Not implementing member access right now
    //        return AST_type::operatorPrecedence::PREC_MEMBER_ACCESS ;
  case OP_TIMES:
  case OP_DIVIDE:
  case OP_MODULUS:
    return AST_type::operatorPrecedence::PREC_MULTIPLICATIVE ;
  case OP_PLUS:
  case OP_MINUS:
    return AST_type::operatorPrecedence::PREC_ADDITIVE ;
  case OP_SHIFT_RIGHT:
  case OP_SHIFT_LEFT:
    return AST_type::operatorPrecedence::PREC_SHIFT ;
  case OP_LT:
  case OP_GT:
  case OP_GE:
  case OP_LE:
    return AST_type::operatorPrecedence::PREC_RELATIONAL ;
  case OP_EQUAL:
  case OP_NOT_EQUAL:
    return AST_type::operatorPrecedence::PREC_EQUALITY ;
  case OP_AND:
    return AST_type::operatorPrecedence::PREC_BITWISE_AND ;
  case OP_EXOR:
    return AST_type::operatorPrecedence::PREC_BITWISE_XOR ;
  case OP_OR:
    return AST_type::operatorPrecedence::PREC_BITWISE_OR ;
  case OP_LOGICAL_AND:
    return AST_type::operatorPrecedence::PREC_LOGICAL_AND ;
  case OP_LOGICAL_OR:
    return AST_type::operatorPrecedence::PREC_LOGICAL_OR ;
  case OP_TERNARY:
  case OP_COLON:
  case OP_ASSIGN:
  case OP_TIMES_ASSIGN:
  case OP_DIVIDE_ASSIGN:
  case OP_MODULUS_ASSIGN:
  case OP_PLUS_ASSIGN:
  case OP_MINUS_ASSIGN:
  case OP_SHIFT_LEFT_ASSIGN:
  case OP_SHIFT_RIGHT_ASSIGN:
  case OP_AND_ASSIGN:
  case OP_OR_ASSIGN:
  case OP_EXOR_ASSIGN:
    return AST_type::operatorPrecedence::PREC_ASSIGNMENT ;
  case OP_COMMA:
    return AST_type::operatorPrecedence::PREC_COMMA ;
  case OP_NIL:
    return AST_type::operatorPrecedence::PREC_ERROR ;
  default:
    return AST_type::operatorPrecedence::PREC_ERROR ;
  }
}

/// Interfce for general ASTP pointer
inline  AST_type::operatorPrecedence getPrecedence(AST_type::ASTP op) {
  return getPrecedence(op->nodeType) ;
}
/// Interface for expression AST node
inline AST_type::operatorPrecedence getPrecedence(CPTR<AST_exprOper> op) {
  return getPrecedence(op->nodeType) ;
}

inline bool rightToLeftAssociativePrecedence(AST_type::operatorPrecedence op) {
  return (op == AST_type::operatorPrecedence::PREC_ASSIGNMENT ||
          op == AST_type::operatorPrecedence::PREC_UNARY_PREFIX) ;
}
inline bool precedenceCompare(AST_type::operatorPrecedence prec_op,
                              AST_type::operatorPrecedence prec_opStack) {
  return ((prec_op == prec_opStack &&
           rightToLeftAssociativePrecedence(prec_op)) ||
          prec_op > prec_opStack) ;
}  

inline bool precedenceCompare(AST_type::ASTP &op,
                              CPTR<AST_exprOper> &opStack) {

  AST_type::operatorPrecedence prec_op = getPrecedence(op);
  AST_type::operatorPrecedence prec_opStack = getPrecedence(opStack) ;
  return precedenceCompare(prec_op,prec_opStack) ;
}

void printStack(vector<CPTR<AST_exprOper> > &exprStack,
                ostream &out) {
  out << "Stack:" << endl ;
  for(size_t i=0;i<exprStack.size();++i) {
    out << i << " - " ;
    AST_simplePrint printer(out,-2) ;
    exprStack[i]->accept(printer) ;
    out << endl ;
  }
}

AST_type::ASTP parseExpressionOperator(AST_type::ASTP expr,
                                       std::istream &is, int &linecount,
                                       const string &fileName,
                                       varmap &typemap,
                                       AST_type::operatorPrecedence
                                       prec) {
  // exprStack holds already processed binary operators
  // create stack and push sentinel on the stack

  // Note that if we have repeated operators that have default left
  // associativity (e.g. a+b+c+d = (((a+b)+c)+d) then this will be structured
  // as a single plus node in the tree, e.g. (+:a b c d)

#ifdef VERBOSE
  cerr << "in parseExpressionOperator" ;
  if(expr == 0)
    cerr << " with nil expr" ;
  cerr << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  vector<CPTR<AST_exprOper> > exprStack ;
  {
    CPTR<AST_exprOper> tmp = new AST_exprOper ;
    tmp->nodeType = OP_NIL ;
    if(expr!=0)
      tmp->terms.push_back(expr) ;

    exprStack.push_back(tmp) ;
#ifdef VERBOSE
    cerr << "push op = " << OPtoName(exprStack.back()->nodeType)
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
  }
  // After getting the first term we are in a loop of searching for operators
  do {
#ifdef VERBOSE
    if(exprStack.size() > 0) {
      printStack(exprStack,cerr) ;
    }
#endif
    // Check if we detect a terminating token
    CPTR<AST_Token> openToken = getToken(is,linecount) ;
    AST_type::elementType opcode = tokenToOperator(openToken->nodeType) ;
    AST_type::operatorPrecedence opprec = getPrecedence(opcode) ;
    if(opprec <= prec) {
      // We found it, restore token and break out of loop
      pushToken(openToken) ;
#ifdef VERBOSE
      cerr << "exiting expression parsing due to precedence for token"
           << OPtoName(openToken->nodeType) << endl ;
#endif
      break ;
    }
    // restore token for parseOperator
    pushToken(openToken) ;
    // binary operator check
    AST_type::ASTP op = parseOperator(is, linecount) ;
#ifdef VERBOSE
    if(op == 0) {
      cerr << "exiting expression parsing due to token"
           << OPtoName(openToken->nodeType) << endl ;
    }
#endif
    if(op == 0)
      break ;

#ifdef VERBOSE
    cerr << "parseExpression operator op = " << OPtoString(op->nodeType)
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    expr = parseExpressionPartial(is,linecount,fileName,typemap) ;
    expr = applyPostFixOperator(expr,is,linecount,fileName,typemap) ;
    
    if(expr == 0) {
      return AST_type::ASTP(new AST_syntaxError("Expecting expression after binary operator",linecount,fileName)) ;
    }

    if(ASTEqual(exprStack.back(),OP_NIL)) {
      // If no operator is parsed yet, we just get started and initilize
      // the left branch
      exprStack.back()->nodeType = op->nodeType ;
      exprStack.back()->terms.push_back(expr) ;
#ifdef VERBOSE
      cerr << "replacing OP_NIL on exprStack with "
           << OPtoName(exprStack.back()->nodeType)
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
      if(exprStack.size() > 0) {
        printStack(exprStack,cerr) ;
      }
#endif

    } else {
      // Now we reorder the tree based on operator precedence
      while(exprStack.size() > 1 &&
            !precedenceCompare(op,exprStack.back())) {
#ifdef VERBOSE
        cerr << "pop off exprStack:" << OPtoName(exprStack.back()->nodeType)
             << ", file: " << __FILE__ << ":" << __LINE__
             << endl ;
#endif
	exprStack.pop_back() ;
#ifdef VERBOSE
        cerr << "after pop top of stack is:" << OPtoName(exprStack.back()->nodeType)
             << ", size=" << exprStack.size()
             << ", file: " << __FILE__ << ":" << __LINE__
             << endl ;
#endif
#ifdef VERBOSE
        if(exprStack.size() > 0) {
          printStack(exprStack,cerr) ;
        }
#endif

      }

      if(ASTEqual(op,OP_COLON)) {
#ifdef VERBOSE
        cerr << "FOUND COLON" << endl ;
        printStack(exprStack,cerr) ;
#endif
        if(exprStack.back()->nodeType != OP_TERNARY) {
          cerr << "expected ':' to be working with '?' in TERNARY operator" << endl ;
          CPTR<AST_Token> optoken = CPTR<AST_Token>(op) ;
          AST_type::ASTP err =
            AST_type::ASTP(new AST_syntaxError("expecting ':' to be paired with '?'",
                                               optoken->lineno,fileName)) ;
          exprStack.back()->terms.push_back(err) ;
        }
      }
      if(!ASTEqual(op,OP_TERNARY) && ASTEqual(op,exprStack.back())) {
	// If operator is the same, just chain the terms
	exprStack.back()->terms.push_back(expr) ;
#ifdef VERBOSE
        cerr << "adding term to operator: "
             << OPtoName(exprStack.back()->nodeType)
             << " stack size=" << exprStack.size()
             << ", file: " << __FILE__ << ":" << __LINE__
             << endl ;
#endif
#ifdef VERBOSE
        if(exprStack.size() > 0) {
          printStack(exprStack,cerr) ;
        }
#endif
      } else if(precedenceCompare(op,exprStack.back())) {
        // if operator is lower precedence than top of stack, then we need
        // to put this operator on the rightmost term
        CPTR<AST_exprOper> np = new AST_exprOper ;
        np->nodeType = op->nodeType ;
#ifdef VERBOSE
        cerr << "top of stack is" << OPtoName(exprStack.back()->nodeType)
             << endl ;
        cerr << "editing rightmost term: " << OPtoName(exprStack.back()->terms.back()->nodeType)
             << ", file: " << __FILE__ << ":" << __LINE__
             << endl ;
#endif
        np->terms.push_back(exprStack.back()->terms.back()) ;
        np->terms.push_back(expr) ;
        exprStack.back()->terms.back() = AST_type::ASTP(np) ;
#ifdef VERBOSE
        cerr << "pushing " << OPtoName(exprStack.back()->nodeType)
             << "to stack,  line " << __LINE__
             << ", file: " << __FILE__ << ":" << __LINE__
             << endl ;
#endif
#ifdef VERBOSE
        if(exprStack.size() > 0) {
          printStack(exprStack,cerr) ;
        }
#endif
        exprStack.push_back(np) ;
#ifdef VERBOSE
        cerr << "pushing on exprStack:" << OPtoName(exprStack.back()->nodeType)
             << " stack size=" << exprStack.size()
             << ", file: " << __FILE__ << ":" << __LINE__
             << endl ;
#endif
#ifdef VERBOSE
        if(exprStack.size() > 0) {
          printStack(exprStack,cerr) ;
        }
#endif
      } else {
	CPTR<AST_exprOper> np = new AST_exprOper ;
	np->nodeType = op->nodeType ;
#ifdef VERBOSE
        cerr << "creating lower precidence op="
             << OPtoName(np->nodeType)
             << ", file: " << __FILE__ << ":" << __LINE__
             << endl ;
#endif
	np->terms.push_back(AST_type::ASTP(exprStack.back())) ;
	np->terms.push_back(expr) ;
	exprStack.back() = np ;
#ifdef VERBOSE
        cerr << "changing exprStack to" << OPtoName(exprStack.back()->nodeType)
             << ", file: " << __FILE__ << ":" << __LINE__
             << endl ;
#endif
#ifdef VERBOSE
        if(exprStack.size() > 0) {
          printStack(exprStack,cerr) ;
        }
#endif
      }
    }
  } while(true) ;


#ifdef VERBOSE
  cerr << "parseExpressionOperator exit loop stack" << endl ;
  if(exprStack.size() > 0) {
    printStack(exprStack,cerr) ;
  }
#endif

  AST_type::ASTP e = CPTR<AST_type>(exprStack.front()) ;
  if(ASTEqual(exprStack.front(),OP_NIL))
    e = CPTR<AST_type>(exprStack.front()->terms.front()) ;

  return applyPostFixOperator(AST_type::ASTP(e),is,linecount,fileName,typemap) ;
}

/// Parse general expression.  This routine mainly handles operator
/// precedence for binary operators, it passes other work to
/// parseExpressionPartial
AST_type::ASTP parseExpression(std::istream &is, int &linecount,
			       const string &fileName,
			       varmap &typemap,
                               AST_type::operatorPrecedence prec) {
#ifdef VERBOSE
  cerr << "in parseExpression"
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif

  AST_type::ASTP expr = parseExpressionPartial(is,linecount,fileName,typemap,prec) ;
  expr = applyPostFixOperator(expr,is,linecount,fileName,typemap) ;
  if(expr == 0) // If no valid expression then return null
    return expr ;

  return parseExpressionOperator(expr,is,linecount,fileName,typemap,prec) ;

}


AST_type::ASTP parseCaseStatement(std::istream &is, int &linecount,
				  const string &fileName,
				  varmap &typemap) {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseCaseStatement, token = " << token->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  if(!ASTEqual(token,TK_CASE) &&
     !ASTEqual(token,TK_DEFAULT)) {
    cerr << "internal error parsing switch statement on line " << linecount
	 << " of file " << fileName << endl ;
  }
  CPTR<AST_Block> AST_data = new AST_Block ;
  AST_data->nodeType = TK_CASE ;
  AST_data->identifiers = typemap ;
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  if(ASTEqual(token,TK_CASE)) {
    AST_type::ASTP expr =
      parseExpression(is,linecount,fileName,AST_data->identifiers,
                      AST_type::operatorPrecedence::PREC_ASSIGNMENT) ;
    AST_data->elements.push_back(expr) ;
  }

  token = getToken(is,linecount) ;
  if(!ASTEqual(token,TK_COLON)) {
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
  cerr << "in parseSwitchStatement, token = " << token->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  AST_type::ASTP switchtok(token) ;
  CPTR<AST_controlStatement> ctrl = new AST_controlStatement ;
  ctrl->identifiers = typemap  ;
  ctrl->controlType = AST_type::ASTP(token) ;
  ctrl->nodeType = ND_CTRL_SWITCH ;

  token = getToken(is,linecount) ;
  if(!ASTEqual(token,TK_OPENPAREN)) {
    pushToken(token) ;
    return AST_type::ASTP(new AST_syntaxError("expecting '(' in switch statment",token->lineno,fileName)) ;
  }
  pushToken(token) ;
  AST_type::ASTP conditional =
    parseExpression(is,linecount,fileName,ctrl->identifiers) ;
  // switch statement is just a list of statements in order
  ctrl->parts.push_back(conditional) ;
  token = getToken(is,linecount) ;
  if(!ASTEqual(token,TK_OPENBRACE)) {
    pushToken(token) ;
    return AST_type::ASTP(new AST_syntaxError("expecting '{' in switch statment",token->lineno,fileName)) ;
  }
  ctrl->parts.push_back(AST_type::ASTP(token)) ;
  token = getToken(is,linecount) ;
  while(!ASTEqual(token,TK_CLOSEBRACE) &&
	!ASTEqual(token,TK_ERROR)) {
    if(ASTEqual(token,TK_CASE) ||
       ASTEqual(token,TK_DEFAULT)) {
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

  if(ASTEqual(token,TK_ERROR))
    return AST_type::ASTP(new AST_syntaxError("failed to find closing brace in swithc statement",linecount,fileName)) ;

  ctrl->parts.push_back(AST_type::ASTP(token)) ;
  return AST_type::ASTP(ctrl) ;

}

AST_type::ASTP parseIfStatement(std::istream &is, int &linecount,
				const string &fileName,
				varmap &typemap)  {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseIfStatement, token = " << token->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  if(!ASTEqual(token,TK_IF)) {
    pushToken(token) ;
    return AST_type::ASTP(new AST_syntaxError("confused in if statement",token->lineno,fileName)) ;
  }
  AST_type::ASTP iftok = AST_type::ASTP(token) ;

  token = getToken(is,linecount) ;
  if(!ASTEqual(token,TK_OPENPAREN)) {
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

  if(ASTEqual(token,TK_ELSE)) {
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
  cerr << "in parseLoopStatement, token = " << token->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  AST_type::ASTP loop = AST_type::ASTP(token) ;
  CPTR<AST_controlStatement> ctrl = new AST_controlStatement ;
  ctrl->identifiers = typemap ;
  ctrl->controlType = loop ;
  switch(token->nodeType) {
  case TK_FOR:
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

      ctrl->nodeType = ND_CTRL_FOR ;

      token = getToken(is,linecount) ;
      if(!ASTEqual(token,TK_OPENPAREN)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("for expecting '('",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;
      token = getToken(is,linecount) ;
#ifdef VERBOSE
      cerr << "for loop parsed TK_OPEN_PAREN, next token = " << OPtoName(token->nodeType)
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif

      AST_type::ASTP initializer = 0 ;
      if(ASTEqual(token,TK_SEMICOLON)) {
        initializer = AST_type::ASTP(token) ;
      } else {
        CPTR<AST_Token> token2 = getToken(is,linecount) ;

        bool isDeclaration = true ;
        if(ASTEqual(token,TK_NAME) &&
           !ASTEqual(token2,TK_SCOPE)) {
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
      cerr << "initializer = " << OPtoName(initializer->nodeType)
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      ctrl->parts.push_back(initializer) ;

#ifdef VERBOSE
      cerr << "for: parsing conditional expression"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl;
#endif
      AST_type::ASTP conditional = parseExpression(is,linecount,fileName,
						   ctrl->identifiers) ;

#ifdef VERBOSE
      cerr << "conditional = " << OPtoName(conditional->nodeType)
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      ctrl->parts.push_back(conditional) ;
      token = getToken(is,linecount) ;
      if(!ASTEqual(token,TK_SEMICOLON)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("No ';' after conditional in for loop",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;
      AST_type::ASTP advance = parseExpression(is,linecount,fileName,
                                               ctrl->identifiers) ;
      ctrl->parts.push_back(advance) ;
      token = getToken(is,linecount) ;
      if(!ASTEqual(token,TK_CLOSEPAREN)) {
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
  case TK_WHILE:
    {
      //--------------------------------------------------------------------
      // parsing while loop
      //
      // Note parts will be
      // conditional
      // loop body
      ctrl->nodeType = ND_CTRL_WHILE ;

      token = getToken(is,linecount) ;
      if(!ASTEqual(token,TK_OPENPAREN)) {
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
  case TK_DO:
    {
      //--------------------------------------------------------------------
      // parsing do loop
      //
      // Note parts will be
      // loop body
      // TK_WHILE
      // conditional
      // TK_SEMICOLON

      ctrl->nodeType = ND_CTRL_DO ;
      AST_type::ASTP body = parseStatement(is,linecount,fileName,
                                           ctrl->identifiers) ;
      ctrl->parts.push_back(body) ;
      token = getToken(is,linecount) ;
      if(!ASTEqual(token,TK_WHILE)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting 'while' in do loop",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;

      token = getToken(is,linecount) ;
      if(!ASTEqual(token,TK_OPENPAREN)) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting '(' in do loop",token->lineno,fileName)) ;
      }
      pushToken(token) ;
      AST_type::ASTP conditional = parseExpression(is,linecount,fileName,
                                                   ctrl->identifiers) ;
      ctrl->parts.push_back(conditional) ;

      token = getToken(is,linecount) ;
      if(!ASTEqual(token,TK_SEMICOLON)) {
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
  cerr << "in getIdentifierName"
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  vector<CPTR<AST_Token>>  token_stack ;
  isFunc = false ;
  string name ;
  CPTR<AST_Token> token = getToken(is,linecount) ;
  // parse the scoped name part
  if(ASTEqual(token,TK_SCOPE)) {
    name += token->text ;
    token_stack.push_back(token) ;
    token = getToken(is,linecount) ;
  }
  while(ASTEqual(token,TK_NAME)) {
    name += token->text ;
    token_stack.push_back(token) ;
    token = getToken(is,linecount) ;
    if(ASTEqual(token,TK_SCOPE)) {
      name += token->text ;
      token_stack.push_back(token) ;
      token = getToken(is,linecount) ;
    } else
      break ;
  }
#ifdef VERBOSE
  cerr << "ident name = " << name
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  // if template, parse template arguments
  if(ASTEqual(token,TK_OPENTEMPLATE)) {
#ifdef VERBOSE
    cerr << "identifier is template"
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    int count = 1 ;
    name += "<" ;
    token_stack.push_back(token) ;
    token = getToken(is,linecount) ;
    while(count > 0) {
      if(ASTEqual(token,TK_OPENTEMPLATE)) {
        count++ ;
        name += "<" ;
      } else if(ASTEqual(token,TK_CLOSETEMPLATE)) {
        count-- ;
        name += ">" ;
      } else {
        name += token->text ;
      }
      token_stack.push_back(token) ;
      token = getToken(is,linecount) ;
    }
  }
  if(ASTEqual(token,TK_OPENPAREN)) {
    isFunc=true ;
  }
  token_stack.push_back(token) ;
  // Push tokens so that they can be parsed
  for(auto ii = token_stack.rbegin();ii!=token_stack.rend();++ii)
    pushToken(*ii) ;
#ifdef VERBOSE
  cerr << "return ident name = " << name
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
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
  cerr << "in parseDeclarationOrSimpleStatement"
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl;
#endif
  // If starts with a unary operator or '(' then it is a simple statement
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
  pushToken(openToken) ;
  if(checkUnaryToken(openToken) || ASTEqual(openToken,TK_OPENPAREN)) {
#ifdef VERBOSE
    cerr << "Starts with operator or grouping so call parseSimpleStatement"
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    return parseSimpleStatement(is,linecount,fileName,typemap) ;
  }
  if(ASTEqual(openToken,TK_STRUCT) || ASTEqual(openToken,TK_CLASS))
    return parseDeclaration(is,linecount,fileName,typemap) ;

  bool isFunc = false ;
  string ident = getIdentifierName(is,linecount,fileName,isFunc) ;
#ifdef VERBOASE
  cerr << "ident = " << ident << "isFunc = " << isFunc << endl ;
#endif
  if(isFunc)
    return parseSimpleStatement(is,linecount,fileName,typemap) ;

  auto ii = typemap.find(ident) ;
#ifdef VERBOSE
  if( ii != typemap.end()) {
    cerr << "typemap defines " << ident
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
  }
#endif
  if(ii != typemap.end() && ii->second.isLocalIdentifier()) {
    return parseSimpleStatement(is,linecount,fileName,typemap) ;
  }
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
    cerr << "in parseDeclaration, token = " << token->text
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    switch(token->nodeType) {
    case TK_CONST:
      // TK_REGISTER,TK_VOLATILE,TK_MUTABLE add later
      AST_data->type_decl.push_back(AST_type::ASTP(token)) ;
#ifdef VERBOSE
      cerr << "in parseDeclaration, const variable"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      break ;
    case TK_EXTERN:
      AST_data->type_decl.push_back(AST_type::ASTP(token)) ;
#ifdef VERBOSE
      cerr << "in parseDeclaration, extern variable"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      break ;
    case TK_CHAR:
    case TK_FLOAT:
    case TK_DOUBLE:
    case TK_INT:
    case TK_VOID:
    case TK_BOOL:
    case TK_SHORT:
    case TK_LONG:
    case TK_SIGNED:
    case TK_UNSIGNED:
    case TK_AUTO:
    case TK_STRUCT:
    case TK_CLASS:
#ifdef VERBOSE
      cerr << "in parseDeclaration, builtin type"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
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
    case TK_SCOPE:
    case TK_NAME:
#ifdef VERBOSE
      cerr << "in parseDeclaration, named type"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      if(builtin_type || defined_type ) {
        istype = false ;
      } else {
        pushToken(token) ;

        AST_type::ASTP typedec =
          parseIdentifier(is,linecount,fileName,typemap) ;
        if(typedec->nodeType == OP_FUNC) {
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
      cerr << "parsed token " << OPtoName(token->nodeType)
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      istype = false ;
      break ;
    }
    if(istype)
      token = getToken(is,linecount) ;
  } while(istype) ;

#ifdef VERBOSE
  cerr << "in parseDeclaration, processing declarations: " << token->text
       << " " << OPtoName(token->nodeType)
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl;
#endif
  while(ASTEqual(token,TK_NAME) ||
        ASTEqual(token,TK_TIMES) ||
        ASTEqual(token,TK_AND)) {
#ifdef VERBOSE
    cerr << "in parseDeclaration, found name, * or & "
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif

    if(ASTEqual(token,TK_TIMES)) {
      AST_data->decls.push_back(AST_type::ASTP(token)) ;
      token = getToken(is,linecount) ;
      if(ASTEqual(token,TK_CONST)) {
        AST_data->decls.push_back(AST_type::ASTP(token)) ;
        token = getToken(is,linecount) ;
      }
    }
    if(ASTEqual(token,TK_AND)) {
      AST_data->decls.push_back(AST_type::ASTP(token)) ;
      token = getToken(is,linecount) ;
    }

#ifdef VERBOSE
    cerr << "in parseDeclaration, processed pointer or reference: " << token->text
         << ", file: " << __FILE__ << ":" << __LINE__
         << endl ;
#endif
    if(ASTEqual(token,TK_NAME)) {
      typemap[token->text] = localIdentifier() ;
      AST_data->decl_varnames.push_back(token->text) ;
      pushToken(token) ;

      AST_type::ASTP expr = parseExpressionPartial(is,linecount,fileName,typemap) ;
      expr = applyPostFixOperator(expr,is,linecount,fileName,typemap) ;
#ifdef VERBOSE
      cerr << "in parseDeclaration, parsing expression "
           << ", file: " << __FILE__ << ":" << __LINE__
           <<  endl ;
#endif
      if(expr != 0) { // If no valid expression then return null
        expr = parseExpressionOperator(expr,is,linecount,fileName,typemap,
                                       AST_type::operatorPrecedence::PREC_COMMA) ;
      } else {
#ifdef VERBOSE
        cerr << "in parseDeclaration, parseExpressionPartial returned nil "
             << ", file: " << __FILE__ << ":" << __LINE__
             <<  endl ;
#endif

        AST_type::ASTP p =
          AST_type::ASTP(new AST_syntaxError("Syntax error variable declaration",
                                             token->lineno,fileName)) ;
        AST_data->decls.push_back(p) ;
      }
      AST_data->decls.push_back(expr) ;
      token = getToken(is,linecount) ;
    }
    if(ASTEqual(token,TK_COMMA)) {
#ifdef VERBOSE
      cerr << "pushing TK_COMMA onto decls"
           << ", file: " << __FILE__ << ":" << __LINE__
           << endl ;
#endif
      AST_data->decls.push_back(AST_type::ASTP(token)) ;
      token = getToken(is,linecount) ;
    }
  }
  if(ASTEqual(token,TK_SEMICOLON)) {
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
  AST_data->nodeType = OP_SPECIAL ;
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseSpecialControlStatement, token = " << token->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  token = getToken(is,linecount) ;
  if(!ASTEqual(token,TK_SEMICOLON)) {
    if(ASTEqual(AST_data->elements.back(),TK_RETURN)) {
      pushToken(token) ;
      AST_type::ASTP exp =
        parseExpression(is,linecount,fileName,AST_data->identifiers) ;
      AST_data->elements.push_back(exp) ;
      token = getToken(is,linecount) ;
      if(!ASTEqual(token,TK_SEMICOLON)) {
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
  cerr << "in parseSimpleStatement"
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  AST_type::ASTP exp = parseExpression(is,linecount,fileName,typemap) ;

  CPTR<AST_Token> termToken = getToken(is,linecount) ;
  AST_type::ASTP term = AST_type::ASTP(termToken) ;
#ifdef VERBOSE
  cerr << "terminal token for simpleStatement is " << OPtoName(term->nodeType)
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  if(!ASTEqual(term,TK_SEMICOLON)) {
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
  cerr << "in parseStatement, token = " << firstToken->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif
  pushToken(firstToken) ;
  switch(firstToken->nodeType) {
  case TK_OPENBRACE:
    return parseBlock(is,linecount,fileName,typemap) ;
  case TK_CHAR:
  case TK_FLOAT:
  case TK_DOUBLE:
  case TK_INT:
  case TK_VOID:
  case TK_BOOL:
  case TK_SHORT:
  case TK_LONG:
  case TK_SIGNED:
  case TK_UNSIGNED:
  case TK_CONST:
  case TK_EXTERN:
  case TK_AUTO:
  case TK_STRUCT:
  case TK_CLASS:
    return parseDeclaration(is,linecount,fileName,typemap) ;
  case TK_FOR:
  case TK_WHILE:
  case TK_DO:
    return parseLoopStatement(is,linecount,fileName,typemap) ;
  case TK_IF:
    return parseIfStatement(is,linecount,fileName,typemap) ;
  case TK_SWITCH:
    return parseSwitchStatement(is,linecount,fileName,typemap) ;
  case TK_BREAK:
  case TK_CONTINUE:
  case TK_RETURN:
    return parseSpecialControlStatement(is,linecount,fileName,typemap) ;
  case TK_NAME:
  case TK_SCOPE:
    return parseDeclarationOrSimpleStatement(is,linecount,fileName,typemap) ;

  case TK_USING:
    {
      firstToken = getToken(is,linecount) ;
      CPTR<AST_declaration> AST_data = new AST_declaration ;
      AST_data->type_decl.push_back(AST_type::ASTP(firstToken)) ;

      firstToken = getToken(is,linecount) ;
      if(ASTEqual(firstToken,TK_NAMESPACE)) {
        AST_data->type_decl.push_back(AST_type::ASTP(firstToken)) ;
      } else {
        pushToken(firstToken) ;
      }
      bool isFunc = false ;
      string s = getIdentifierName(is,linecount,fileName,isFunc) ;
      if(!isFunc) {
        typemap[s] = localIdentifier() ;
      }
      AST_type::ASTP typedec =
        parseIdentifier(is,linecount,fileName,typemap) ;
      AST_data->type_decl.push_back(typedec) ;
      firstToken = getToken(is,linecount) ;
      if(!ASTEqual(firstToken,TK_SEMICOLON)) {
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

  case TK_OPENPAREN:
  case TK_PLUS:
  case TK_MINUS:
  case TK_NOT:
  case TK_AND:
  case TK_TIMES:
  case TK_INCREMENT:
  case TK_DECREMENT:
  case TK_LOCI_VARIABLE:
  case TK_LOCI_CONTAINER:
    {
      return parseSimpleStatement(is,linecount,fileName,typemap) ;
    }

  case TK_MACRO:
    firstToken = getToken(is,linecount) ;
    return AST_type::ASTP(firstToken) ;
  case TK_SEMICOLON:
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
  cerr << "in parseBlock, token = " << openToken->text
       << ", file: " << __FILE__ << ":" << __LINE__
       << endl ;
#endif

  AST_type::elementType closeType = TK_CLOSEBRACE ;
  switch(openToken->nodeType) {
  case TK_OPENBRACE:
    closeType = TK_CLOSEBRACE ;
    break ;
  default:
    return AST_type::ASTP(openToken) ;
  }

  CPTR<AST_Block> AST_data = new AST_Block ;
  AST_data->identifiers = typemap ;
  AST_data->nodeType = TK_BRACEBLOCK ;
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
  if(ASTEqual(s,TK_LOCI_VARIABLE)) {
    Loci::variable v(s.text) ;
    vmap_info vmap ;
    vmap.var += v ;
    accessed.insert(vmap) ;
    id2var[s.id] = v ;
  }
}
AST_type::ASTP getArrayVar(AST_type::ASTP expr) {
  while(ASTEqual(expr,OP_ARRAY)) {
    CPTR<AST_exprOper> p(expr) ;
    expr = AST_type::ASTP(p->terms.front()) ;
  }
  return expr ;
}
bool isExprLociVariable(AST_type::ASTP expr) {
  expr = getArrayVar(expr) ;

  return ASTEqual(expr,TK_LOCI_VARIABLE) ;
}


void AST_collectAccessInfo::visit(AST_exprOper &s) {
  switch(s.nodeType) {
  case OP_ASSIGN:
  case OP_TIMES_ASSIGN:
  case OP_DIVIDE_ASSIGN:
  case OP_MODULUS_ASSIGN:
  case OP_PLUS_ASSIGN:
  case OP_MINUS_ASSIGN:
  case OP_SHIFT_LEFT_ASSIGN:
  case OP_SHIFT_RIGHT_ASSIGN:
  case OP_AND_ASSIGN:
  case OP_OR_ASSIGN:
  case OP_EXOR_ASSIGN:
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
  case OP_ARROW:
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
  case OP_GROUP:
    out << '(' ;
    for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
      if(*ii != 0)
	(*ii)->accept(*this) ;
    out << ')' ;
    break ;
  case OP_CAST:
    out << '(' ;
    if(s.terms.size() >= 1 && s.terms[0] != 0)
      s.terms[0]->accept(*this) ;
    out << ')' ;
    if(s.terms.size() >= 2 && s.terms[1] != 0)
      s.terms[1]->accept(*this) ;
    break ;
  case OP_TEMPLATE_CAST:
    for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
      if(*ii != 0)
	(*ii)->accept(*this) ;
    break ;
  case OP_BRACEBLOCK:
    out << '{' ;
    if(s.terms.size() >= 1 && s.terms[0] != 0)
      s.terms[0]->accept(*this) ;
    out << '}' ;
    break ;
  case OP_FUNC:
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
  case OP_TEMPLATE:
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
  case OP_ARRAY:
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
  case OP_TERNARY:
    {
      AST_type::ASTList::iterator ii=s.terms.begin() ;
      
      if(ii!=s.terms.end() && *ii != 0) {
	(*ii)->accept(*this) ;
      }
      if(ii!=s.terms.end()) 
        ++ii ;
      out << '?' ;
      if(ii!=s.terms.end() && *ii != 0)
	(*ii)->accept(*this) ;
    }
    break ;

  case OP_UNARY_PLUS:
  case OP_UNARY_MINUS:
  case OP_NOT:
  case OP_AMPERSAND:
  case OP_STAR:
  case OP_INCREMENT:
  case OP_DECREMENT:
    {
      string op = OPtoString(s) ;
      out << op ;
      for(auto ii=s.terms.begin();ii!=s.terms.end();++ii)
	if(*ii != 0)
	  (*ii)->accept(*this) ;
    }
    break ;
  case OP_POSTINCREMENT:
  case OP_POSTDECREMENT:
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
  if(lineno != -2 && s.lineno >0 && lineno != s.lineno) {
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
    if(ASTEqual(s,TK_LOCI_DIRECTIVE)) {
      out << "$[" << s.text << "]" ;
    } else if(ASTEqual(s,TK_LOCI_CONTAINER)) {
      out << "$*" << s.text ;
    } else if(ASTEqual(s,TK_LOCI_VARIABLE)) {
      out << "$" << s.text ;
    } else if(ASTEqual(s,TK_MACRO)) {
      out << endl << '#' << s.text << endl ;
    } else
      out <<s.text << ' ' ;
  }
}

bool condenseOp(AST_type::elementType et) {
  switch(et) {
  case OP_ARROW:
  case OP_SCOPE:
  case OP_TIMES:
  case OP_PLUS:
  case OP_COMMA:
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
