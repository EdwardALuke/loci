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

// Setup for facilties for parsing and creating abstract syntax trees (AST)

class AST_visitor ;

struct varinfo {
  bool isType ;
  bool isTemplate ;
varinfo(): isType(true), isTemplate(false) {}
varinfo(bool v, bool t) : isType(v),isTemplate(t) {}
} ;

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
		    OP_ARROW=0x100, 
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
		    OP_COMMA=0xe00, OP_DOT,
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
} ;
  
class AST_Token : public AST_type {
public:
  string text ;
  int lineno ;
  void accept(AST_visitor &v) ;
  AST_Token() {nodeType = AST_type::OP_ERROR; }
} ;

class AST_SimpleStatement: public AST_type {
public:
  AST_SimpleStatement(ASTP e, ASTP t) : exp(e),Terminal(t)
    {nodeType = AST_type::ND_SIMPLE_STATEMENT ;}
  ASTP exp ;
  ASTP Terminal ;
  void accept(AST_visitor &v) ;
} ;

class AST_Block : public AST_type {
public:
  ASTList elements ;
  void accept(AST_visitor &v) ;
  AST_Block() {nodeType = AST_type::ND_BLOCK; }
} ;

class AST_declaration : public AST_type {
public:
  ASTList type_decl ;
  ASTP decls ;
  void accept(AST_visitor &v) ;
  AST_declaration() { nodeType = AST_type::ND_DECL; }
} ;

class AST_exprOper : public AST_type {
public:
  ASTList terms ;
  void accept(AST_visitor &v) ;
  AST_exprOper() {nodeType = AST_type::OP_ERROR; }
} ;

class AST_controlStatement: public AST_type {
 public:
  ASTP controlType ;
  ASTList parts ;
  void accept(AST_visitor &v) ;
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

class AST_errorCheck : public AST_visitor {
 public:
  int error_count ;
  AST_errorCheck() {error_count = 0  ; }
  bool hasErrors() const {return error_count != 0; }
  virtual ~AST_errorCheck() {}
  virtual void visit(AST_syntaxError &) ;
} ;

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
  

extern AST_type::ASTP parseExpression(std::istream &is, int &linecount,
				      const string &fileName,
				      const varmap &typemap) ;
extern AST_type::ASTP parseExpressionPartial(std::istream &is, int &linecount,
					     const string &fileName,
					     const varmap &typemap) ;
extern AST_type::ASTP parseDeclaration(std::istream &is, int &linecount,
				       const string &fileName,
				       const varmap &typemap) ;
extern AST_type::ASTP parseStatement(std::istream &is, int &linecount,
				     const string &fileName,
				     const varmap &typemap) ;
extern AST_type::ASTP parseLoopStatement(std::istream &is, int &linecount,
					 const string &fileName,
					 const varmap &typemap) ;
extern AST_type::ASTP parseIfStatement(std::istream &is, int &linecount,
				       const string &fileName,
				       const varmap &typemap) ;
extern AST_type::ASTP parseSwitchStatement(std::istream &is, int &linecount,
					   const string &fileName,
					   const varmap &typemap)  ;
extern AST_type::ASTP parseBlock(std::istream &is, int &linecount,
				 const string &fileName,
				 const varmap &typemap) ;
extern AST_type::ASTP parseTerm(std::istream &is, int &linecount) ;

extern CPTR<AST_Token> getToken(std::istream &is, int &linecount) ;

extern void pushToken(CPTR<AST_Token> &pt) ;
extern bool isTerm(AST_type::elementType e) ;

extern istream &killsp(istream &s, int &lines) ;
extern bool is_name(istream &s) ;
extern string get_name(istream &s) ;




#endif
