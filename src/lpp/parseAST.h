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
#ifndef PARSEAST_H
#define PARSEAST_H

#include <Tools/cptr.h>
using Loci::CPTR ;
using Loci::CPTR_type ;
#include <list>
using std::list ;
#include <vector>
using std::vector ;
#include <string>
using std::string ;

#include <map>
using std::map ;

#include <iostream>
using std::istream ;
using std::ifstream ;
using std::ofstream ;
using std::ostream ;
using std::ios ;
using std::endl ;
using std::cerr ;
using std::cout ;

#include "Tools/variable.h"

// Setup for facilties for parsing and creating Abstract Syntax Trees (AST)

class AST_visitor ;

struct varinfo {
  /// Identifier is a type definition
  bool statusType ;
  /// Identifier is a local variable
  bool statusLocal ;
  /// Type is a templated type
  bool statusTemplate ;
  /// identifier is a function call
  bool statusFunction ;
  /// identifier is defined external to the rule
  bool statusExternal ;
  /// identifier is undefined, assumed to be external
  bool statusUndef ;
  varinfo(): statusType(true), statusLocal(false), statusTemplate(false),
             statusFunction(false),statusExternal(true), statusUndef(true)  {}
  varinfo(bool v, bool t) : statusType(v),statusLocal(false),statusTemplate(t),
                            statusFunction(false),statusExternal(true),
                            statusUndef(true) {}

  bool isTypeName() { return statusType ; }
  bool isLocalIdentifier() { return statusLocal ; }
} ;

/// create varinfo type for local variable 
inline varinfo localIdentifier() {
  varinfo v ;
  v.statusType=false ;
  v.statusLocal = true ;
  v.statusTemplate = false ;
  v.statusFunction = false ;
  v.statusExternal = false ;
  v.statusUndef = false ;
  return v ;
}
/// create varinfo for identifier that is assumed to be a type name as
/// a default
inline varinfo defaultTypeName() {
  varinfo v ;
  v.statusType     = true ;
  v.statusLocal    = false ;
  v.statusTemplate = false ;
  v.statusFunction = false ;
  v.statusExternal = true ;
  v.statusUndef    = true ;
  return v ;
}
/// create varinfo for templated identifier that is assumed to be a type name
/// as a default
inline varinfo defaultTemplateTypeName() {
  varinfo v ;
  v.statusType     = true ;
  v.statusLocal    = false ;
  v.statusTemplate = true ;
  v.statusFunction = false ;
  v.statusExternal = true ;
  v.statusUndef    = true ;
  return v ;
}

typedef map<std::string,varinfo> varmap ;

// Abstract Syntax Tree
class AST_type : public CPTR_type {
public:
  virtual ~AST_type() {}
  AST_type() ;
  typedef CPTR<AST_type> ASTP ;
  typedef std::vector<ASTP> ASTList ;
  enum elementType {
		    OP_SCOPE=0x000,
		    OP_AT=0x080, // For using @ to separate namespaces
		    // Traditional C operators
                    OP_DOT=0x100,
		    OP_ARROW,
		    OP_TIMES = 0x300, OP_DIVIDE, OP_MODULUS,
		    OP_PLUS  = 0x400, OP_MINUS, 
		    OP_SHIFT_RIGHT = 0x500, OP_SHIFT_LEFT,
		    OP_LT = 0x600, OP_GT, OP_GE, OP_LE,
		    OP_EQUAL = 0x700, OP_NOT_EQUAL, 
		    OP_AND=0x800, OP_EXOR=0x900, OP_OR=0xa00,
		    OP_LOGICAL_AND=0xb00, OP_LOGICAL_OR=0xc00,
		    OP_TERNARY,
		    OP_ASSIGN=0xd00,
		    OP_TIMES_ASSIGN,
		    OP_DIVIDE_ASSIGN,
		    OP_MODULUS_ASSIGN,
		    OP_PLUS_ASSIGN,
		    OP_MINUS_ASSIGN,
		    OP_SHIFT_LEFT_ASSIGN,
		    OP_SHIFT_RIGHT_ASSIGN,
		    OP_AND_ASSIGN,
		    OP_OR_ASSIGN,
		    OP_EXOR_ASSIGN,
		    OP_COMMA=0xe00, 
		    OP_COLON=0xf00,
		    OP_SEMICOLON,
		    // terminal for empty statement
		    OP_NIL=0x1000,
		    // terminals for variable name, function, array or name{args}
		    OP_INCREMENT,OP_DECREMENT,
		    OP_POSTINCREMENT, OP_POSTDECREMENT,
		    OP_COMMENT,
		    OP_BRACEBLOCK,
		    OP_NAME, OP_FUNC, OP_ARRAY, OP_NAME_BRACE, OP_FUNC_BRACE,
                    OP_TEMPLATE,OP_TEMPLATE_FUNC,
		    // terminal for string, integer, or unspecified error condition
		    OP_STRING, OP_NUMBER, OP_ERROR,
		    // Unary operations
		    OP_UNARY_PLUS, OP_UNARY_MINUS, OP_NOT, OP_TILDE,
		    OP_AMPERSAND, OP_DOLLAR, OP_STAR,
		    OP_GROUP,OP_GROUP_ERROR,
		    OP_OPENPAREN,OP_CLOSEPAREN,OP_OPENBRACKET,OP_CLOSEBRACKET,
		    OP_OPENBRACE,OP_CLOSEBRACE,
		    OP_LOCI_DIRECTIVE,OP_LOCI_VARIABLE,OP_LOCI_CONTAINER,
		    OP_TERM, OP_SPECIAL,
		    TK_BRACEBLOCK=0x2000,
		    TK_SCOPE,
		    TK_AT, // For using @ to separate namespaces
		    // Traditional C operators
		    TK_ARROW, 
		    TK_TIMES, TK_DIVIDE, TK_MODULUS,
		    TK_PLUS, TK_MINUS, 
		    TK_SHIFT_RIGHT, TK_SHIFT_LEFT,
		    TK_LT, TK_GT, TK_GE, TK_LE,
		    TK_EQUAL, TK_NOT_EQUAL, 
		    TK_AND, TK_EXOR, TK_OR,
		    TK_LOGICAL_AND, TK_LOGICAL_OR, 
		    TK_ASSIGN,
		    TK_TIMES_ASSIGN,
		    TK_DIVIDE_ASSIGN,
		    TK_MODULUS_ASSIGN,
		    TK_PLUS_ASSIGN,
		    TK_MINUS_ASSIGN,
		    TK_SHIFT_LEFT_ASSIGN,
		    TK_SHIFT_RIGHT_ASSIGN,
		    TK_AND_ASSIGN,
		    TK_OR_ASSIGN,
		    TK_EXOR_ASSIGN,
		    TK_COMMA, TK_DOT,
		    TK_COLON,
		    TK_SEMICOLON,
		    // terminal for empty statement
		    TK_NIL,
		    // terminals for variable name, function, array or name{args}
		    TK_INCREMENT,TK_DECREMENT,TK_COMMENT,
		    TK_NAME, 
		    // terminal for string, integer, or unspecified error condition
		    TK_STRING, TK_NUMBER, TK_ERROR,
		    // Unary operations
		    TK_UNARY_PLUS, TK_UNARY_MINUS, TK_NOT, TK_TILDE,
		    TK_QUESTION, TK_AMPERSAND, TK_STAR,
		    TK_OPENPAREN,TK_CLOSEPAREN,TK_OPENBRACKET,TK_CLOSEBRACKET,
		    TK_OPENBRACE,TK_CLOSEBRACE,
                    TK_OPENTEMPLATE,TK_CLOSETEMPLATE,
		    TK_LOCI_DIRECTIVE,TK_LOCI_VARIABLE,TK_LOCI_CONTAINER,
		    // Now the keywords
		    TK_ALIGNAS, TK_ALIGNOF, TK_ASM, 
		    TK_BOOL,TK_FALSE,TK_TRUE,TK_CHAR,TK_INT,TK_LONG,
		    TK_SHORT,TK_SIGNED,TK_UNSIGNED,TK_DOUBLE,TK_FLOAT,TK_ENUM,
		    TK_MUTABLE,TK_CONST,TK_STATIC,TK_VOLATILE,TK_AUTO,
		    TK_REGISTER,TK_EXPORT,TK_EXTERN,TK_INLINE,TK_NAMESPACE,
		    TK_EXPLICIT,TK_DYNAMIC_CAST,TK_STATIC_CAST,
		    TK_REINTERPRET_CAST,
		    TK_OPERATOR,TK_PROTECTED,TK_NOEXCEPT,TK_NULLPTR,
		    TK_RETURN,TK_SIZEOF,TK_THIS,TK_TYPEID,
		    TK_SWITCH,TK_CASE,TK_BREAK,TK_DEFAULT,
		    TK_FOR,TK_DO,TK_WHILE,TK_CONTINUE,
		    TK_CLASS,TK_STRUCT,TK_PUBLIC,TK_PRIVATE,TK_FRIEND,
		    TK_UNION,TK_TYPENAME,TK_TEMPLATE,TK_TYPEDEF,
		    TK_VIRTUAL,TK_VOID,TK_TRY,TK_CATCH,TK_THROW,
		    TK_IF,TK_ELSE,TK_GOTO,TK_NEW,TK_DELETE,
		    ND_SYNTAXERR,
		    ND_CTRL_IF,ND_CTRL_FOR,ND_CTRL_WHILE, ND_CTRL_DO,
		    ND_CTRL_SWITCH, ND_SIMPLE_STATEMENT,ND_BLOCK,
		    ND_DECL,ND_TERMINAL,
		    TK_SENTINEL 
		    
  } ;
  int id ;
  elementType nodeType ;
  virtual void accept(AST_visitor &v) = 0 ;
  virtual ASTP clone() const = 0 ;
  
} ;

extern std::string OPtoString(AST_type::elementType val) ;

class AST_syntaxError: public AST_type {
 public:
  string error ;
  int lineno ;
  string fileName ;
  AST_syntaxError(string e,int l, string f) :error(e),lineno(l),fileName(f)
    {nodeType=AST_type::ND_SYNTAXERR;}
  void accept(AST_visitor &v) ;
  ASTP clone() const ;
} ;

/// Token in AST structure typically will be a terminal symbol
class AST_Token : public AST_type {
public:
  /// string that comprises this token
  string text ;
  /// line number where this token was read from file
  int lineno ;
  /// visitor interface
  void accept(AST_visitor &v) ;
  ASTP clone() const ;
  /// Default constructor is an error, will be overridden during parsing
  AST_Token() {nodeType = AST_type::OP_ERROR; ; lineno=-1; }
} ;

/// A statement that is typically terminated by a semicolon
class AST_SimpleStatement: public AST_type {
public:
  AST_SimpleStatement(ASTP e, ASTP t) : exp(e),Terminal(t)
    {nodeType = AST_type::ND_SIMPLE_STATEMENT ;}
  /// Expression of statement.
  ASTP exp ;
  /// Statement Termination
  ASTP Terminal ;
  /// visitor interface
  void accept(AST_visitor &v) ;
  ASTP clone() const ;
} ;

/// A block of statements (enclosed by braces)
class AST_Block : public AST_type {
public:
  ASTList elements ;
  varmap identifiers ;
  void accept(AST_visitor &v) ;
  ASTP clone() const ;
  AST_Block() {nodeType = AST_type::ND_BLOCK; }
} ;

/// Variable declaration statement
class AST_declaration : public AST_type {
public:
  ASTList type_decl ;
  ASTList decls ;
  std::vector<std::string> decl_varnames ;
  void accept(AST_visitor &v) ;
  ASTP clone() const ;
  AST_declaration() { nodeType = AST_type::ND_DECL; }
} ;

/// Operator in expression
class AST_exprOper : public AST_type {
public:
  ASTList terms ;
  void accept(AST_visitor &v) ;
  ASTP clone() const ;
  AST_exprOper() {nodeType = AST_type::OP_ERROR; }
} ;

/// Control statement (if else 
class AST_controlStatement: public AST_type {
 public:
  ASTP controlType ;
  ASTList parts ;
  varmap identifiers ;
  void accept(AST_visitor &v) ;
  ASTP clone() const ;
  AST_controlStatement() {nodeType = AST_type::OP_ERROR; }
  void constructIf(ASTP tok, ASTP C, ASTP IFBLOCK, ASTP ELSE, ASTP ELSEBLOCK) {
    nodeType = AST_type::ND_CTRL_IF ;
    controlType = tok ;
    parts.push_back(C) ;
    parts.push_back(IFBLOCK) ;
    if(ELSE != 0) {
      parts.push_back(ELSE) ;
      parts.push_back(ELSEBLOCK) ;
    }
  }
} ;


/// Visitor abstract base class
class AST_visitor {
 public :
  virtual ~AST_visitor() {} ;
  virtual void visit(AST_SimpleStatement &)  ;
  virtual void visit(AST_Block &)  ;
  virtual void visit(AST_declaration &)  ;
  virtual void visit(AST_exprOper &)  ;
  virtual void visit(AST_controlStatement &) ;
  virtual void visit(AST_Token &) {}
  virtual void visit(AST_syntaxError &) {}
} ;

/// Visitor class that scans an AST to determine if any systax errors were
/// created during parsing
class AST_errorCheck : public AST_visitor {
 public:
  int error_count ;
  AST_errorCheck() {error_count = 0  ; }
  bool hasErrors() const {return error_count != 0; }
  virtual ~AST_errorCheck() {}
  virtual void visit(AST_syntaxError &) ;
} ;



/// Visitor that prints an AST using a simple substitution map
class AST_simplePrint : public AST_visitor {
 public:
  ostream &out ;
  int lineno ;
  map<int,std::string> id2rename ;
  bool prettyPrint ;
  AST_simplePrint(ostream &s, int line=-1,bool pp=true): out(s),lineno(line),prettyPrint(pp) {} 
  virtual void visit(AST_exprOper &)  ;
  virtual void visit(AST_Token &) ;
} ;

class AST_collectAccessInfo: public AST_visitor {
public:
  std::set<Loci::vmap_info> accessed ;
  std::set<Loci::vmap_info> writes ;
  std::map<int,Loci::variable> id2var ;
  std::map<int,Loci::vmap_info> id2vmap ;
  AST_collectAccessInfo() {} ;
  virtual void visit(AST_exprOper &) ;
  virtual void visit(AST_Token &) ;
} ;
  

/// Parse a general expression
extern AST_type::ASTP parseExpression(std::istream &is, int &linecount,
				      const string &fileName,
				      varmap &typemap) ;
/// Parse terms excluding operator
extern AST_type::ASTP parseExpressionPartial(std::istream &is, int &linecount,
					     const string &fileName,
					     varmap &typemap) ;
/// Parse a statement after the first identifier is parsed
extern AST_type::ASTP parseExpressionOperator(AST_type::ASTP expr,
                                              std::istream &is, int &linecount,
                                              const string &fileName,
                                              varmap &typemap,
                                              AST_type::elementType
                                              term=AST_type::TK_SENTINEL) ;
/// Parse operator token, convert from TK_* node type to OP_* nodetype
extern AST_type::ASTP parseOperator(std::istream &is, int &linecount) ;
/// Apply postfix operator to previous expression passed in expr
extern AST_type::ASTP applyPostFixOperator(AST_type::ASTP expr,
                                           std::istream &is, int &linecount,
                                           string fileName,
                                           varmap &typemap) ;
/// Parse a named object using the scope operator
extern AST_type::ASTP parseScopedObject(std::istream &is, int &linecount,
                                        const string &fileName,
                                        varmap &typemap) ;
/// Parse template instantiation arguments
extern AST_type::ASTP parseTemplateArguments(AST_type::ASTP objname,
                                             std::istream &is, int &linecount,
                                             const string &fileName,
                                             varmap &typemap) ;
/// Parse function evocation arguments
extern AST_type::ASTP parseFunctionArguments(AST_type::ASTP objname,
                                             std::istream &is, int &linecount,
                                             const string &fileName,
                                             varmap &typemap) ;

/// Parse a declaration statement
extern AST_type::ASTP parseDeclaration(std::istream &is, int &linecount,
				       const string &fileName,
				       varmap &typemap) ;

/// Parse an inut statement, this could be a code block denoted by braces,
/// type declaration, loop control statement, an if statement,
/// a switch statement, a special control statement (e.g. break, control,
/// or return, a simple statement (expression) followed by a semicolon,
/// or an empty statement.
extern AST_type::ASTP parseStatement(std::istream &is, int &linecount,
				     const string &fileName,
				     varmap &typemap) ;
/// Parse a looping control statement
extern AST_type::ASTP parseLoopStatement(std::istream &is, int &linecount,
					 const string &fileName,
					 varmap &typemap) ;
/// parse an if (else) control statement
extern AST_type::ASTP parseIfStatement(std::istream &is, int &linecount,
				       const string &fileName,
				       varmap &typemap) ;
/// parse a switch control statement
extern AST_type::ASTP parseSwitchStatement(std::istream &is, int &linecount,
					   const string &fileName,
					   varmap &typemap)  ;
/// parse a case statement within a switch statement
extern AST_type::ASTP parseCaseStatement(std::istream &is, int &linecount,
                                         const string &fileName,
                                         varmap &typemap) ;
/// Parse return, continue, or break control statements
extern AST_type::ASTP parseSpecialControlStatement(std::istream &is,
                                                   int &linecount,
                                                   const string &fileName,
                                                   varmap &typemap) ;

/// Parse a brace enclosed block of statements
extern AST_type::ASTP parseBlock(std::istream &is, int &linecount,
				 const string &fileName,
				 varmap &typemap) ;
/// Parse a terminal symbol (name, number, or string
extern AST_type::ASTP parseTerm(std::istream &is, int &linecount) ;

/// Get a token from the lexical analyzer from input stream and update
/// line count while parsing input stream.
extern CPTR<AST_Token> getToken(std::istream &is, int &linecount) ;

/// Push a token (effectively undo reading a token from the input stream
extern void pushToken(CPTR<AST_Token> &pt) ;


extern istream &killsp(istream &s, int &lines) ;
extern bool is_name(istream &s) ;
extern string get_name(istream &s) ;


/// Convert operator code to corresponding string for printing
inline string OPtoString(AST_type::elementType val) {
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

inline string OPtoString(const AST_type &val)
{ return OPtoString(val.nodeType) ; }

/// Convert operator code to corresponding string for printing
inline string OPtoName(AST_type::elementType val) {
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
    return string("OP_UNARY_MINUD") ;
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
  default:
    {
      std::ostringstream oss ;
      oss << "ND_UNDEF(0x" << std::hex << int(val) << ")" ;
      return oss.str() ;
    }
  }
  return string("/*error*/") ;
}

inline string OPtoName(const AST_type &val)
{ return OPtoName(val.nodeType) ; }

/// Check if two AST nodes are equal
template<class S, class T> inline bool
ASTEqual(const CPTR<S> &t1, const CPTR<T> &t2) {
  return t1->nodeType == t2->nodeType ;
}
/// Check if two AST nodes are equal
template<class S> inline bool
ASTEqual(const CPTR<S> &t1, AST_type::elementType t2) {
  return t1->nodeType == t2 ;
}
//. Check if type of AST node
inline bool ASTEqual(const AST_type &s, AST_type::elementType t) {
  return s.nodeType == t ;
}

#endif
