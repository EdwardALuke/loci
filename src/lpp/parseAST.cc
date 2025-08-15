//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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

static int parseIDctr = 0 ;
AST_type::AST_type() {
  nodeType = OP_ERROR ;
  id = parseIDctr++ ;
}

void AST_Token::accept(AST_visitor &v) {  v.visit(*this) ; }

void AST_syntaxError::accept(AST_visitor &v) {  v.visit(*this) ; }

void AST_SimpleStatement::accept(AST_visitor &v) {  v.visit(*this) ; }

void AST_Block::accept(AST_visitor &v) {  v.visit(*this) ; }

void AST_declaration::accept(AST_visitor &v) {  v.visit(*this) ; }

void AST_exprOper::accept(AST_visitor &v) {  v.visit(*this) ; }

void AST_controlStatement::accept(AST_visitor &v) {  v.visit(*this) ; }

AST_type::ASTP parseTerm(std::istream &is, int &linecount) {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseTerm, token = " << token->text << endl ;
#endif
  switch(token->nodeType) {
  case AST_type::TK_NAME:
  case AST_type::TK_TRUE:
  case AST_type::TK_FALSE:
  case AST_type::TK_NUMBER:
  case AST_type::TK_STRING:
  case AST_type::TK_LOCI_VARIABLE:
  case AST_type::TK_LOCI_CONTAINER:
    return AST_type::ASTP(token) ;
  default:
    return AST_type::ASTP(0) ;
  }
  return AST_type::ASTP(0) ;
}
  
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
#ifdef VERBOSE
    cerr << "returning OP_TERNARY" << endl ;
#endif
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

bool checkUnaryToken(AST_type::elementType e) {
  return (e == AST_type::TK_PLUS ||
	  e == AST_type::TK_MINUS ||
	  e == AST_type::TK_NOT ||
	  e == AST_type::TK_AND ||
	  e == AST_type::TK_TIMES ||
	  e == AST_type::TK_INCREMENT ||
	  e == AST_type::TK_DECREMENT) ;
}

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

bool checkPostFixToken(AST_type::elementType e) {
  return (e == AST_type::TK_INCREMENT ||
	  e == AST_type::TK_DECREMENT) ;
}
  
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

AST_type::ASTP applyPostFixOperator(AST_type::ASTP expr,
				    std::istream &is, int &linecount,
				    string fileName,
				    const varmap &typemap) {
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
       openToken->nodeType == AST_type::TK_OPENBRACKET)
      return applyPostFixOperator(AST_type::ASTP(post),is,linecount,
				  fileName,typemap) ;
    else
      return AST_type::ASTP(post) ;
  }
  if(openToken->nodeType == AST_type::TK_OPENBRACKET) {
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
       openToken->nodeType == AST_type::TK_OPENBRACKET)
      return applyPostFixOperator(AST_type::ASTP(array),is,linecount,fileName,
				  typemap) ;
    else
      return AST_type::ASTP(array) ;
  }
  pushToken(openToken) ;
    
  return expr ;
}


AST_type::ASTP parseExpressionPartial(std::istream &is, int &linecount,
				      const string &fileName,
				      const varmap &typemap) {
  CPTR<AST_Token> openToken = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseExpressionPartial, token = " << openToken->text << endl ;
#endif
  if(openToken->nodeType == AST_type::TK_OPENPAREN) {
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
  if(checkUnaryToken(openToken->nodeType)) {
    //check for unary operators
    AST_type::ASTP expr = parseExpressionPartial(is,linecount,fileName,
						 typemap) ;
    if(expr!= 0) {
      CPTR<AST_exprOper> unary = new AST_exprOper ;
      unary->nodeType = unaryOperator(openToken->nodeType) ;
      unary->terms.push_back(expr) ;
      return applyPostFixOperator(AST_type::ASTP(unary),is,linecount,fileName,
				  typemap) ;
    }
  }
  if(isTerm(openToken->nodeType)) {
    pushToken(openToken) ;
    AST_type::ASTP exp = parseTerm(is,linecount);

    openToken = getToken(is,linecount) ;

    if(openToken->nodeType == AST_type::TK_OPENPAREN) {
      // Function
      AST_type::ASTP args = parseExpression(is,linecount,fileName,typemap) ;
      CPTR<AST_Token> closeToken = getToken(is,linecount) ;
      CPTR<AST_exprOper> func = new AST_exprOper ;
      func->nodeType = AST_type::OP_FUNC ;
      func->terms.push_back(AST_type::ASTP(exp)) ;
      func->terms.push_back(args) ;
      if(closeToken->nodeType != AST_type::TK_CLOSEPAREN) {
	CPTR<AST_Token> err = new AST_Token ;
	*err = *closeToken ;
	err->nodeType = AST_type::OP_ERROR ;
	func->terms.push_back(AST_type::ASTP(err)) ;
	pushToken(closeToken) ;
      }
      return applyPostFixOperator(AST_type::ASTP(func),is,linecount,fileName,
				  typemap) ;
    }
    pushToken(openToken) ;
    return applyPostFixOperator(exp,is,linecount,fileName,typemap) ;
    //    return exp ;
  }
  pushToken(openToken) ;
  return 0 ;
}
  


AST_type::ASTP parseExpression(std::istream &is, int &linecount,
			       const string &fileName,
			       const varmap &typemap) {
#ifdef VERBOSE
  cerr << "in parseExpression"<< endl ;
#endif
  AST_type::ASTP expr = parseExpressionPartial(is,linecount,fileName,typemap) ;
  if(expr == 0) // If no valid expression then return null
    return expr ;
  vector<CPTR<AST_exprOper> > exprStack ;
  {
    CPTR<AST_exprOper> tmp = new AST_exprOper ;
    tmp->nodeType = AST_type::OP_NIL ;
    tmp->terms.push_back(expr) ;
    exprStack.push_back(tmp) ;
  }
  const unsigned int mask = ~0x7f ;
  // After getting the first term we are in a loop of searching for operators
  do {
    // binary operator check
    AST_type::ASTP op = parseOperator(is, linecount) ;
    if(op == 0)
      break ;

    if(op->nodeType == AST_type::OP_TERNARY ) 
      expr = parseExpression(is,linecount,fileName,typemap) ;
    else
      expr = parseExpressionPartial(is,linecount,fileName,typemap) ;
    
    if(expr == 0) {
      return AST_type::ASTP(new AST_syntaxError("Expecting expression after binary operator",linecount,fileName)) ;
    }
    
    if(exprStack.back()->nodeType == AST_type::OP_NIL) {
      // If no operator is parsed yet, we just get started and initilize
      // the left branch
      exprStack.back()->nodeType = op->nodeType ;
      exprStack.back()->terms.push_back(expr) ;
    } else {
      // Now we reorder the tree based on operator precedence
      while(exprStack.size() >1 &&
	    ((op->nodeType&mask) >= (mask&exprStack.back()->nodeType))) {
#ifdef VERBOSE
	if(op->nodeType == AST_type::OP_TERNARY) {
	  cerr << "popping stack when OP_TERNARY" << endl ;
	}
#endif
	exprStack.pop_back() ;
      }
      if(op->nodeType == exprStack.back()->nodeType) {
	// If operator is the same, just chain the terms
	exprStack.back()->terms.push_back(expr) ;
      } else if(((op->nodeType)&mask) < ((exprStack.back()->nodeType)&mask)) {
        // if operator is lower precedence
	CPTR<AST_exprOper> np = new AST_exprOper ;
	np->nodeType = op->nodeType ;
	np->terms.push_back(exprStack.back()->terms.back()) ;
	np->terms.push_back(expr) ;
	exprStack.back()->terms.back() = AST_type::ASTP(np) ;
	if(op->nodeType == AST_type::OP_TERNARY ) {
	  CPTR<AST_Token> op2 = getToken(is,linecount) ;
#ifdef VERBOSE
	  cerr << " detected OP_TERNARY, op2 =" << op2->text << endl ;
#endif
	  if(op2->nodeType == AST_type::TK_COLON) {
	    expr = parseExpressionPartial(is,linecount,fileName,typemap) ;
	    np->terms.push_back(expr) ;
	  } else {
	    pushToken(op2) ;
	    return AST_type::ASTP(new AST_syntaxError("expecting ':' in tertiary operator",op2->lineno,fileName)) ;
	  }
	}
	exprStack.push_back(np) ;
      } else {
	CPTR<AST_exprOper> np = new AST_exprOper ;
	np->nodeType = op->nodeType ;
	np->terms.push_back(AST_type::ASTP(exprStack.back())) ;
	np->terms.push_back(expr) ;
	if(op->nodeType == AST_type::OP_TERNARY) {
	  CPTR<AST_Token> op2 = getToken(is,linecount) ;
	  if(op2->nodeType != AST_type::TK_COLON) {
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
      }
    }
  } while(true) ;

  return applyPostFixOperator(AST_type::ASTP(exprStack.front()),is,linecount,fileName,typemap) ;
}


AST_type::ASTP parseCaseStatement(std::istream &is, int &linecount,
				  const string &fileName,
				  const varmap &typemap) {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseCaseStatement, token = " << token->text << endl ;
#endif
  if(token->nodeType != AST_type::TK_CASE &&
     token->nodeType != AST_type::TK_DEFAULT) {
    cerr << "internal error parsing switch statement on line " << linecount
	 << " of file " << fileName << endl ;
  }
  CPTR<AST_Block> AST_data = new AST_Block ;
  AST_data->nodeType = AST_type::TK_CASE ;
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  if(token->nodeType == AST_type::TK_CASE ) {
    AST_type::ASTP expr = parseExpression(is,linecount,fileName,typemap) ;
    AST_data->elements.push_back(expr) ;
  }
  
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_COLON) {
    pushToken(token) ; 
    return AST_type::ASTP(new AST_syntaxError("expecting ':' in case statement",token->lineno,fileName)) ;
  }
  
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  return AST_type::ASTP(AST_data) ;
}

AST_type::ASTP parseSwitchStatement(std::istream &is, int &linecount,
				    const string &fileName,
				    const varmap &typemap)  {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseSwitchStatement, token = " << token->text << endl ;
#endif
  AST_type::ASTP switchtok(token) ;
  CPTR<AST_controlStatement> ctrl = new AST_controlStatement ;
  ctrl->controlType = AST_type::ASTP(token) ;
  ctrl->nodeType = AST_type::ND_CTRL_SWITCH ;

  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_OPENPAREN) {
    pushToken(token) ;
    return AST_type::ASTP(new AST_syntaxError("expecting '(' in switch statment",token->lineno,fileName)) ;
  }
  pushToken(token) ;
  AST_type::ASTP conditional = parseExpression(is,linecount,fileName,typemap) ;
  // switch statement is just a list of statements in order
  ctrl->parts.push_back(conditional) ;
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_OPENBRACE) {
    pushToken(token) ;
    return AST_type::ASTP(new AST_syntaxError("expecting '{' in switch statment",token->lineno,fileName)) ;
  }
  ctrl->parts.push_back(AST_type::ASTP(token)) ;
  token = getToken(is,linecount) ;
  while(token->nodeType != AST_type::TK_CLOSEBRACE &&
	token->nodeType != AST_type::TK_ERROR) {
    if(token->nodeType == AST_type::TK_CASE ||
       token->nodeType == AST_type::TK_DEFAULT) {
      pushToken(token) ;
      ctrl->parts.push_back(parseCaseStatement(is,linecount,fileName,typemap)) ;
    } else {
      pushToken(token) ;
      ctrl->parts.push_back(parseStatement(is,linecount,fileName,typemap)) ;
    }
    token = getToken(is,linecount) ;
  }
			
  if(token->nodeType == AST_type::TK_ERROR) 
    return AST_type::ASTP(new AST_syntaxError("failed to find closing brace in swithc statement",linecount,fileName)) ;

  ctrl->parts.push_back(AST_type::ASTP(token)) ;
  return AST_type::ASTP(ctrl) ;

}

AST_type::ASTP parseIfStatement(std::istream &is, int &linecount,
				const string &fileName,
				const varmap &typemap)  {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseIfStatement, token = " << token->text << endl ;
#endif
  if(token->nodeType != AST_type::TK_IF) {
    pushToken(token) ;
    return AST_type::ASTP(new AST_syntaxError("confused in if statement",token->lineno,fileName)) ;
  }
  AST_type::ASTP iftok = AST_type::ASTP(token) ;
  
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_OPENPAREN) {
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
  
  if(token->nodeType == AST_type::TK_ELSE) {
    elsetok = AST_type::ASTP(token) ;
    ebody = parseStatement(is,linecount,fileName,typemap) ;
  } else {
    pushToken(token) ;
  }
  CPTR<AST_controlStatement> ctrl = new AST_controlStatement ;
  ctrl->constructIf(iftok,conditional,body,elsetok,ebody) ;
  return AST_type::ASTP(ctrl) ;
}			 

bool isTypeDecl(CPTR<AST_Token> p, const varmap &typemap) {
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
  case AST_type::TK_NAME:
    {
      auto t = typemap.find(p->text) ;
      if(t == typemap.end())
	return false ;
      return t->second.isType ;
    }
  default:
    return false ;
  }
}


AST_type::ASTP parseLoopStatement(std::istream &is, int &linecount,
				  const string &fileName,
				  const varmap &typemap) {
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseLoopStatement, token = " << token->text << endl ;
#endif
  AST_type::ASTP loop = AST_type::ASTP(token) ;
  CPTR<AST_controlStatement> ctrl = new AST_controlStatement ;
  ctrl->controlType = loop ;
  switch(token->nodeType) {
  case AST_type::TK_FOR:
    {

      //--------------------------------------------------------------------
      // parsing for loop
      //
      // Note for the for loop the parts will be
      // TK_OPENPAREN
      // initializer
      // TK_SEMICOLON
      // conditional
      // TK_SEMICOLON
      // advance
      // TK_CLOSEPAREN
      // for body
      
      ctrl->nodeType = AST_type::ND_CTRL_FOR ;
      
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_OPENPAREN) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("for expecting '('",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;
      
      token = getToken(is,linecount) ;
      AST_type::ASTP initializer = 0 ;
      if(isTypeDecl(token,typemap)) {
	pushToken(token) ;
	initializer = parseDeclaration(is,linecount,fileName,typemap) ;
      } else {
	pushToken(token) ;
	initializer = parseExpression(is,linecount,fileName,typemap) ;
      }
      ctrl->parts.push_back(initializer) ;
	
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_SEMICOLON) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("No ';' after initializer in for loop",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;
      
      AST_type::ASTP conditional = parseExpression(is,linecount,fileName,
						   typemap) ;
      ctrl->parts.push_back(conditional) ;
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_SEMICOLON) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("No ';' after conditional in for loop",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;
      AST_type::ASTP advance = parseExpression(is,linecount,fileName,typemap) ;
      ctrl->parts.push_back(advance) ;
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_CLOSEPAREN) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting ')' in for loop",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;
      AST_type::ASTP body = parseStatement(is,linecount,fileName,typemap) ;
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
      if(token->nodeType != AST_type::TK_OPENPAREN) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting '(' in while loop",token->lineno,fileName)) ;
      }
      pushToken(token) ;
      AST_type::ASTP conditional = parseExpression(is,linecount,fileName,
						   typemap) ;
      ctrl->parts.push_back(conditional) ;

      AST_type::ASTP body = parseStatement(is,linecount,fileName,typemap) ;
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
      AST_type::ASTP body = parseStatement(is,linecount,fileName,typemap) ;
      ctrl->parts.push_back(body) ;
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_WHILE) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting 'while' in do loop",token->lineno,fileName)) ;
      }
      ctrl->parts.push_back(AST_type::ASTP(token)) ;

      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_OPENPAREN) {
	pushToken(token) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting '(' in do loop",token->lineno,fileName)) ;
      }
      pushToken(token) ;
      AST_type::ASTP conditional = parseExpression(is,linecount,fileName,typemap) ;
      ctrl->parts.push_back(conditional) ;

      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_SEMICOLON) {
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

AST_type::ASTP parseType(std::istream &is, int &linecount,
			 const varmap &typemap) {
  //  cout << "parsing type" << endl ;
  
  CPTR<AST_declaration> AST_data = new AST_declaration ;

  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseType, token = " << token->text << endl ;
#endif
  while(isTypeDecl(token,typemap)) {
    AST_data->type_decl.push_back(AST_type::ASTP(token)) ;
    token = getToken(is,linecount) ;
  }
  pushToken(token) ;
  
  return AST_type::ASTP(AST_data) ;
}

AST_type::ASTP parseDeclaration(std::istream &is, int &linecount,
				const string &fileName,
				const varmap &typemap) {

  //  cout << "parsing declaration" << endl ;
  
  CPTR<AST_declaration> AST_data = new AST_declaration ;

  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseDeclaration, token = " << token->text << endl ;
#endif
  while(isTypeDecl(token,typemap)) {
    AST_data->type_decl.push_back(AST_type::ASTP(token)) ;
    token = getToken(is,linecount) ;
  }
  pushToken(token) ;
  AST_data->decls = parseExpression(is,linecount,fileName,typemap) ;
  
  return AST_type::ASTP(AST_data) ;
}

AST_type::ASTP parseSpecialControlStatement(std::istream &is, int &linecount,
					    const string &fileName,
					    const varmap &typemap) {
  
  CPTR<AST_Block> AST_data = new AST_Block ;
  AST_data->nodeType = AST_type::OP_SPECIAL ;
  CPTR<AST_Token> token = getToken(is,linecount) ;
#ifdef VERBOSE
  cerr << "in parseSpecialControlStatement, token = " << token->text << endl ;
#endif
  AST_data->elements.push_back(AST_type::ASTP(token)) ;
  token = getToken(is,linecount) ;
  if(token->nodeType != AST_type::TK_SEMICOLON) {
    if(AST_data->elements.back()->nodeType == AST_type::TK_RETURN) {
      pushToken(token) ;
      AST_type::ASTP exp = parseExpression(is,linecount,fileName,typemap) ;
      AST_data->elements.push_back(exp) ;
      token = getToken(is,linecount) ;
      if(token->nodeType != AST_type::TK_SEMICOLON) {
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

AST_type::ASTP parseStatement(std::istream &is, int &linecount,
			      const string &fileName,
			      const varmap &typemap) {
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
    {
      CPTR<AST_Token> tok1 = getToken(is,linecount) ;
      CPTR<AST_Token> tok2 = getToken(is,linecount) ;
      switch(tok2->nodeType) {
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
      case AST_type::TK_NAME:
	pushToken(tok2) ;
	pushToken(tok1) ;
	return parseDeclaration(is,linecount,fileName,typemap) ;
      default:
	pushToken(tok2) ;
	pushToken(tok1) ;
	break ;
      }
      if(isTypeDecl(firstToken,typemap))
	return parseDeclaration(is,linecount,fileName,typemap) ;
      
      AST_type::ASTP exp = parseExpression(is,linecount,fileName,typemap) ;

      CPTR<AST_Token> termToken = getToken(is,linecount) ;
      AST_type::ASTP term = AST_type::ASTP(termToken) ;
      if(term->nodeType != AST_type::TK_SEMICOLON) {
	pushToken(termToken) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting ';' ",
						  termToken->lineno,fileName)) ;
      }
      AST_type::ASTP stat = new AST_SimpleStatement(exp,term) ;
      return stat ;
    }

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
      AST_type::ASTP exp = parseExpression(is,linecount,fileName,typemap) ;

      CPTR<AST_Token> termToken = getToken(is,linecount) ;
      AST_type::ASTP term = AST_type::ASTP(termToken) ;
      
      if(term->nodeType != AST_type::TK_SEMICOLON) {
	pushToken(termToken) ;
	return AST_type::ASTP(new AST_syntaxError("Expecting ';' ",
						  termToken->lineno,fileName)) ;
      }
      AST_type::ASTP stat = new AST_SimpleStatement(exp,term) ;
      return stat ;
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
  //  cerr << "got to end with token " << firstToken->text << " on line " << firstToken->lineno << endl ;
  //  return AST_type::ASTP(firstToken) ;
}

AST_type::ASTP parseBlock(std::istream &is, int &linecount,
			  const string &fileName,
			  const varmap &typemap) {
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
  AST_data->nodeType = AST_type::TK_BRACEBLOCK ;
  AST_data->elements.push_back(AST_type::ASTP(openToken)) ;
  CPTR<AST_Token> token = getToken(is,linecount) ;
  while(token->nodeType != closeType) {
    pushToken(token) ;
    CPTR<AST_type> statement = parseStatement(is,linecount,fileName,typemap) ;
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
  for(AST_type::ASTList::iterator ii=s.elements.begin();ii!=s.elements.end();++ii)
    if(*ii!=0)
      (*ii)->accept(*this) ;
}

void AST_visitor::visit(AST_declaration &s) {
  for(AST_type::ASTList::iterator ii=s.type_decl.begin();ii!=s.type_decl.end();++ii)
    if(*ii != 0)
      (*ii)->accept(*this) ;

  if(s.decls != 0)
    s.decls->accept(*this) ;
}
void AST_visitor::visit(AST_exprOper &s) {
  for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();++ii)
    if(*ii != 0)
      (*ii)->accept(*this) ;

}
void AST_visitor::visit(AST_controlStatement &s) {
  s.controlType->accept(*this) ;
  for(AST_type::ASTList::iterator ii=s.parts.begin();ii!=s.parts.end();++ii) {
    if(*ii != 0)
      (*ii)->accept(*this) ;
  }
}

void AST_errorCheck::visit(AST_syntaxError &s) {
  cerr << s.fileName << ":" << s.lineno << ": syntax, " <<s.error << endl ;
  error_count++ ;
}  

void AST_collectAccessInfo::visit(AST_Token &s) {
  if(s.nodeType == AST_type::TK_LOCI_VARIABLE) {
    Loci::variable v(s.text.substr(1,s.text.size()-1)) ;
    vmap_info vmap ;
    vmap.var += v ;
    accessed.insert(vmap) ;
    id2var[s.id] = v ;
  }
}
AST_type::ASTP getArrayVar(AST_type::ASTP expr) {
  while(expr->nodeType == AST_type::OP_ARRAY) {
    CPTR<AST_exprOper> p(expr) ;
    expr = AST_type::ASTP(p->terms.front()) ;
  }
  return expr ;
}
bool isExprLociVariable(AST_type::ASTP expr) {
  expr = getArrayVar(expr) ;
  
  return(expr->nodeType == AST_type::TK_LOCI_VARIABLE) ;
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
      AST_type::ASTList::iterator ii=s.terms.begin();
      AST_collectAccessInfo first ;
      if(ii!=s.terms.end() && *ii != 0)
	(*ii)->accept(first) ;
      for(std::set<Loci::vmap_info>::iterator fi=first.accessed.begin();
	  fi!= first.accessed.end();++fi)
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
	    
	    Loci::variable v(tok->text.substr(1,tok->
text.size()-1)) ;
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
    for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();++ii)
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
    for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();++ii)
      if(*ii != 0)
	(*ii)->accept(*this) ;
    out << ')' ;
    break ;
  case AST_type::OP_FUNC:
    {
      AST_type::ASTList::iterator ii=s.terms.begin() ;
      if(*ii != 0)
	(*ii)->accept(*this) ;
      ++ii ;
      out << '(' ;
      if(*ii != 0)
	(*ii)->accept(*this) ;
      out << ')' ;
      ++ii ;
      if(ii!=s.terms.end()) {
	cerr << "internal error processing func" ;
	if(*ii != 0)
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
      string op = OPtoString(s.nodeType) ;
      out << op ;
      for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();++ii)
	if(*ii != 0)
	  (*ii)->accept(*this) ;
    }
    break ;
  case AST_type::OP_POSTINCREMENT:
  case AST_type::OP_POSTDECREMENT:
    {
      string op = OPtoString(s.nodeType) ;
      for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();++ii)
	if(*ii != 0)
	  (*ii)->accept(*this) ;
      out << op ;
    }
    break ;
  default:
    {
      string op = OPtoString(s.nodeType) ;
      for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();) {
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
  if(lineno != s.lineno) {
    out << endl ;
    if(!prettyPrint)
      if(lineno < 0 || lineno+1 != s.lineno)
	out << "#line " << s.lineno << endl ;

    lineno = s.lineno ;
  }
  auto t = id2rename.find(s.id) ;
  if(t != id2rename.end()) {
    out << t->second ;
  } else
    out <<s.text << ' ' ;
}
