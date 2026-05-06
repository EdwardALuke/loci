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
#include "lpp.h"
#include "parseAST.h"
#include <ctype.h>
#include <set>
#include <iostream>
#include <sstream>
//#include <sys/timeb.h>
#include <time.h>
#include <vector>
#include <unistd.h>

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

using Loci::variable ;
using Loci::variableSet ;
using Loci::vmap_info ;
using Loci::exprList ;
using Loci::exprP ;
using Loci::expression ;
using Loci::exprError ;

bool is_name(istream &s) {
  int ch = s.peek() ;
  return isalpha(ch) || ch == '_' ;
}
    
string get_name(istream &s) {
  if(!is_name(s))
    throw parseError("expected name") ;
  string str ;
  while(!s.eof() && (s.peek() != EOF) &&
        (isalnum(s.peek()) || (s.peek() == '_')) )
    str += s.get() ;
  
  return str ;
}

bool is_string(istream &s) {
  return s.peek() == '\"' ;
}
    
string get_string(istream &s) {
  if(!is_string(s))
    throw parseError("expected string") ;
  string str ;
  if(s.eof())
    throw parseError("unexpected EOF") ;
  s.get() ;
  int ch = s.get() ;
  while(ch != '\"' &&!s.eof()) {
    str += ch ;
    ch = s.get() ;
  }
  if(ch!='\"')
    throw parseError("no closing \" for string") ;
  return str ;
}

bool is_comment(istream &s) {
  if(s.peek() != '/')
    return false ;

  s.get() ;
  char c = s.peek() ;
  s.unget() ;
  if(c == '/' || c == '*')
    return true ;
  return false ;
}

istream &killComment(istream &s, int & lines) {
  s.get() ;
  char c = s.get()  ;
  if(c == '/') { // read to end of line
    while(s.peek() != EOF && s.peek() !='\n') {
      s.get() ;
    }
    if(s.peek() == '\n') {
      lines++ ;
      s.get() ;
    }
    return s ;
  }
  for(;;) {
    if(s.peek() == EOF)
      break ;
    char c = s.get() ;
    if(c == '\n')
      lines++ ;
    if(c == '*') {
      if(s.peek() == '/') {
        s.get() ;
        break ;
      }
    }
  }
  return s ;
}
    
istream &killsp(istream &s, int &lines) {

  bool foundstuff = false ;
  do {
    foundstuff = false ;
    while(s.peek() == ' ' || s.peek() == '\t' || s.peek() == '\n'
          || s.peek() == '\r') {
      if(s.peek() == '\n') lines++ ;
      s.get();
      foundstuff = true ;
    }
    if(is_comment(s)) {
      killComment(s,lines) ;
      foundstuff = true ;
    }
  } while(foundstuff) ;
  return s ;
}


string killCommentOut(istream &s, int & lines,ostream &out) {
  string current_comment ;
  s.get() ;
  out << '/' ;
  char c = s.get()  ;
  out << c ;
  if(c == '/') { // read to end of line
    if(s.peek() == '/' || s.peek() == '!') {
      // if JavaDoc style comment record string, otherwise ignore
      c = s.get() ;
      out << c ;
      while(s.peek() != EOF && s.peek() !='\n') {
	char c = s.get() ;
	current_comment += c ;
	out << c ;
      }
      current_comment += '\n' ;
    } else {
      while(s.peek() != EOF && s.peek() !='\n') {
	char c = s.get() ;
	out << c ;
      }
    }
    if(s.peek() == '\n') {
      lines++ ;
      s.get() ;
      out << '\n' ;
    }
    return current_comment ;
  }
  // Now check for javadoc style comment
  bool javadoc = false ;
  if(s.peek() == '*') {
    char c = s.get() ;
    out << c ;
    if(s.peek() == ' ' || s.peek() == '\t') {
      char c = s.get() ;
      out << c ;
      javadoc = true ;
    }
  }
  for(;;) {
    if(s.peek() == EOF)
      break ;
    char c = s.get() ;
    out << c ;
    if(c == '\n')
      lines++ ;
    if(c == '*') {
      if(s.peek() == '/') {
        out << '/' ;
        s.get() ;
        break ;
      }
    }
    if(javadoc)
      current_comment += c ;
  }
  return current_comment ;
}
    
string killspOut(istream &s, int &lines, ostream &out) {

  string current_comment ;
  bool foundstuff = false ;
  do {
    foundstuff = false ;
    while(s.peek() == ' ' || s.peek() == '\t' || s.peek() == '\n'
          || s.peek() == '\r') {
      if(s.peek() == '\n') lines++ ;
      char c = s.get();
      out << c ;
      foundstuff = true ;
    }
    if(is_comment(s)) {
      current_comment += killCommentOut(s,lines,out) ;
      foundstuff = true ;
    }
  } while(foundstuff) ;
  return current_comment ;
}

inline bool spaceChar(char c) {
  return (c == ' ' || c == '\t' || c == '\n' || c == '*') ;
}

string cleanupComment(const string &s) {
  string cleancomment ;
  size_t i = 0 ;
  // skip initial spaces
  while(i<s.size() && spaceChar(s[i]))
    i++ ;

  for(;i<s.size();++i) {
    if(spaceChar(s[i])) {
      // skip over spaces, replace with single space
      while(i+1<s.size() && spaceChar(s[i+1])) 
	i++ ;
      cleancomment += ' ' ;
    } else if(s[i] == '\\') {
      cleancomment += "\\\\" ;
    } else if(s[i] == '"') {
      cleancomment += "\\\"" ;
    } else if(s[i] >= ' ' &&  s[i] <='~') // valid ascii character
      cleancomment += s[i] ;
  }
  return cleancomment ;
}


int parseFile::killsp() {
  int l = line_no ;
  ::killsp(is,line_no) ;
  return l-line_no ;
}

string parseFile::killspout(std::ostream &outputFile) {
  return ::killspOut(is,line_no,outputFile) ;
}

class parsebase {
public:
  int lines ;
  parsebase(): lines(0) { }
  istream &killsp(istream &s) {
    ::killsp(s,lines) ;
    return s ;
  }
} ;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
template<class T> class funclist : public parsebase {
public:
  list<T> flist ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(s.peek() != '(')
      return s ;
    char c = s.get();
    parsebase::killsp(s) ;
    for(;;) {
      T tmp ;
      tmp.get(s) ;
      flist.push_back(tmp) ;
      parsebase::killsp(s) ;
      if(s.peek() == ')') {
        c = s.get() ;
        return s ;
      }
      if(s.peek() != ',') {
        throw parseError("syntax error") ;
      }
      s.get(); // get comma
    }
  }
  string str() const {
    string s ;
    if(flist.begin() != flist.end()) {
      s += "(" ;
      auto ii = flist.begin() ;
      s+= ii->str() ;
      ++ii ;
      for(;ii!=flist.end();++ii) {
        s+="," ;
        s+= ii->str() ;
      }
      s += ")" ;
    }
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    for(auto ii=flist.begin();ii!=flist.end();++ii) {
      i+= ii->num_lines() ;
    }
    return i ;
  }
} ;
  
#pragma GCC diagnostic pop

template<class T> class templlist : public parsebase {
public:
  list<T> flist ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(s.peek() != '<')
      return s ;
    char c ;
    s.get(c);
    parsebase::killsp(s) ;
    for(;;) {
      T tmp ;
      tmp.get(s) ;
      flist.push_back(tmp) ;
      parsebase::killsp(s) ;
      if(s.peek() == '>') {
        s.get() ;
        return s ;
      }
      if(s.peek() != ',') {
        throw parseError("syntax error, expected comma!") ;
      }
      s.get(); // get comma
    }
  }
  string str() const {
    string s ;
    if(flist.begin() != flist.end()) {
      s += "<" ;
      auto ii = flist.begin() ;
      s+= ii->str() ;
      ++ii ;
      for(;ii!=flist.end();++ii) {
        s+="," ;
        s+= ii->str() ;
      }
      s += "> " ;
    }
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    for(auto ii=flist.begin();ii!=flist.end();++ii) {
      i+= ii->num_lines() ;
    }
    return i ;
  }
} ;

class typestuff : public parsebase {
public:
  string name ;
  templlist<typestuff> templ_args ;
  string scopedPostfix ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(isalpha(s.peek()) || s.peek() == '_') {
      char c = s.peek() ;
      while(isalpha(c = s.peek()) || isdigit(c) || c == '_' ||  c == ':')
        name += s.get() ;
    } else if(isdigit(s.peek())) {
      while(isdigit(s.peek())) {
        name += s.get() ;
      }
    } else
      throw parseError("syntax error") ;
    templ_args.get(s) ;
    if(s.peek() == ':') {
      char c = s.peek() ;
      while(isalpha(c = s.peek()) || isdigit(c) || c == '_' ||  c == ':')
        scopedPostfix += s.get() ;
    }
    return s ;
  }
  string str() const {
    string s ;
    s+= name ;
    s+= templ_args.str() ;
    s+= scopedPostfix ;
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    i+= templ_args.num_lines() ;
    return i ;
  }
} ;
class bracestuff : public parsebase {
public:
  string stuff ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(s.peek() == '{') {
      char c = s.get() ;
      while(s.peek() != EOF && s.peek() != '}') {
        c = s.get() ;
        if(c == '{')
          throw parseError("syntax error") ;
        stuff += c ;
        parsebase::killsp(s) ;
      }
      if(s.peek() == EOF)
        throw parseError("unexpected EOF") ;
      c = s.get() ;
      parsebase::killsp(s) ;
    }
    return s ;
  }
    
  string str() const {
    string s ;
    if(stuff == "")
      return s ;
    s += "{" ;
    s += stuff ;
    s += "}" ;
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    return i ;
  }
} ;
  

class var : public parsebase {
public:
  bool isdollar ;
  string name ;
  list<string> prio_list ;
  list<string> nspace_list ;
  funclist<var> param_args ;
  bracestuff bs ;
  var() : isdollar(false) {}
  
  istream &get(istream &s) {
    isdollar = false ;
    parsebase::killsp(s) ;
    if(s.peek() == '$') {
      s.get() ;
      isdollar=true ;
    }
    if(!is_name(s))
      throw parseError("syntax error: expecting name after '$'") ;
    name = get_name(s) ;
    parsebase::killsp(s) ;
    if(s.peek() == ':') {
      while(s.peek() == ':') {
        s.get() ;
        if(s.peek() != ':') {
	  string err = "syntax error, improper trailing colon, use parenthesis around variable '"+ name+"' to fix." ;
	  throw parseError(err.c_str()) ; 
	}
        s.get() ;
        parsebase::killsp(s) ;
        prio_list.push_back(name);
        if(!is_name(s)) 
          throw parseError("syntax error near ':'") ;
        name = get_name(s) ;
        parsebase::killsp(s) ;
      }
    }
    if(s.peek() == '@') {
      while(s.peek() == '@') {
        s.get() ;
        parsebase::killsp(s) ;
        nspace_list.push_back(name);
        if(!is_name(s)) 
          throw parseError("syntax error near '@'") ;
        name = get_name(s) ;
        parsebase::killsp(s) ;
      }
    }
    
    param_args.get(s) ;
    bs.get(s) ;

    return s ;
  }
  string str() const {
    string s ;
    for(auto li=prio_list.begin();li!=prio_list.end();++li)
      s+= *li + "::" ;
    if(isdollar)
      s+="$" ;
    for(auto li=nspace_list.begin();li!=nspace_list.end();++li)
      s += *li + "@" ;
    s+=name ;
    s+= param_args.str() ;
    s+= bs.str() ;
    return s ;
  }
  int num_lines() const {
    int i = lines ;
    i += param_args.num_lines() ;
    i += bs.num_lines() ;
    return i ;
  }
} ;

class nestedparenstuff : public parsebase {
public:
  string paren_contents ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(s.peek() != '(')
      throw parseError("syntax error, expecting '('") ;
    s.get() ;
    int open_parens = 0 ;
    parsebase::killsp(s) ;
    while(s.peek() != ')' || open_parens != 0) {
      if(s.peek() == '"') { // grab string
	paren_contents += s.get() ;
	while(s.peek() != '"') {
	  if(s.peek() == EOF) {
	    throw parseError("unexpected EOF parsing string") ;
	  }
	  if(s.peek() == '\n' ) {
	    lines++ ;
	  }
	  paren_contents += s.get() ;
	}
	paren_contents += s.get() ;
	continue ;
      }
      if(s.peek() == EOF)
        throw parseError("unexpected EOF") ;
      if(s.peek() == '(')
        open_parens++ ;
      if(s.peek() == ')')
        open_parens-- ;
      if(s.peek() == '\n') {
        s.get() ;
        lines++ ;
        continue ;
      }
      paren_contents += s.get() ;
      parsebase::killsp(s) ;
    }
    s.get() ;
    parsebase::killsp(s) ;
    return s ;
  }
  string str() {
    return paren_contents ;
  }
  int num_lines() {
    return lines ;
  }
} ;

class nestedbracketstuff : public parsebase {
public:
  string bracket_contents ;
  istream &get(istream &s) {
    parsebase::killsp(s) ;
    if(s.peek() != '[')
      throw parseError("syntax error, expecting '['") ;
    s.get() ;
    parsebase::killsp(s) ;
    int open_brackets = 0 ;
    while(s.peek() != ']' || open_brackets != 0) {
      if(s.peek() == EOF)
        throw parseError("unexpected EOF") ;
      if(s.peek() == '[')
        open_brackets++ ;
      if(s.peek() == ']')
        open_brackets-- ;
      if(s.peek() == '\n' ) {
        s.get() ;
        lines++ ;
        continue ;
      }
      bracket_contents += s.get() ;
      parsebase::killsp(s) ;
    }
    s.get() ;
    return s ;
  }
  string str() {
    return bracket_contents ;
  }
  int num_lines() {
    return lines ;
  }
} ;

variable convertVariable(variable v) {
  variable::info vinfo = v.get_info() ;
  vinfo.priority = std::vector<std::string>() ;
  for(size_t i=0;i<vinfo.v_ids.size();++i) {
    std::ostringstream ss ;
    ss << 'X' << i << endl ;
    variable xi = variable(ss.str()) ;
    vinfo.v_ids[i] = xi.ident() ;
  }
  return variable(vinfo) ;
}

void parseFile::setup_Type(std::ostream &outputFile,const string &comment) {
  var vin ;
  vin.get(is) ;
  typedoc tdoc ;
  tdoc.filename = filename ;
  tdoc.lineno = line_no ;
  tdoc.comment = comment ;
  typestuff tin ;
  tin.get(is) ;
  while(is.peek() == ' ' || is.peek() == '\t') 
    is.get() ;
  if(is.peek() != ';')
    throw parseError("syntax error, missing ';'") ;
  is.get() ;
  variable v(vin.str()) ;
  outputFile << "// $type " << v << ' ' << tin.str() ;
  int nl = vin.num_lines()+tin.num_lines() ;
  line_no += nl ;
  for(int i=0;i<nl;++i)
    outputFile << endl ;

  if(type_map.find(v) != type_map.end()) {
    // check to see if the type is changing
    auto mi = type_map.find(v) ;
    if(mi->second.container != tin.name ||
       mi->second.container_args !=tin.templ_args.str()) {
      cerr << filename << ":" << line_no << ":1: warning: variable " << v << " retyped!" <<endl << "Did you intend to change the type of this variable?  If so, use $untype " << v << "; to silence warning" << endl ;
    }
  }
  tdoc.container = tin.name ;
  tdoc.container_args = tin.templ_args.str() ;
  tdoc.v = v ;
  v = convertVariable(v) ;
  type_map[v] = tdoc ;
}

void parseFile::setup_Untype(std::ostream &outputFile) {
  var vin ;
  vin.get(is) ;
  while(is.peek() == ' ' || is.peek() == '\t') 
    is.get() ;
  if(is.peek() != ';')
    throw parseError("syntax error, missing ';'") ;
  is.get() ;
  variable v(vin.str()) ;
  v = convertVariable(v) ;
  outputFile << "// $untype " << v << ";";
  int nl = vin.num_lines() ;
  line_no += nl ;
  for(int i=0;i<nl;++i)
    outputFile << endl ;
  auto mi = type_map.find(v) ;
  if(mi == type_map.end()) {
    cerr << filename << ":" << line_no << ":1: warning: variable " << v << " not defined for untype directive!" <<endl ;
  } else
    type_map.erase(mi) ;
}


namespace {
  inline void fill_descriptors(set<vmap_info> &v, const exprList &in) {

    using namespace Loci ;
    for(auto i = in.begin();i!=in.end();++i) {
      // This needs to be improved to use an actual variable syntax
      // certification.  This test will just get the blindingly obvious
      if((*i)->op != OP_ARROW &&
         (*i)->op != OP_NAME &&
         (*i)->op != OP_FUNC &&
         (*i)->op != OP_NAME_BRACE &&
         (*i)->op != OP_FUNC_BRACE &&
         (*i)->op != OP_SCOPE &&
         (*i)->op != OP_AT &&
         (*i)->op != OP_DOLLAR) {
        cerr << "malformed descriptor: " ;
        (*i)->Print(cerr) ;
        cerr << endl ;
        throw parseError("rule signature error") ;
      }
      vmap_info di(*i) ;
      if(v.find(di) != v.end())
        cerr << "Warning, duplicate variable in var set." << endl ;
      else
        v.insert(di) ;
    }
  }
}

void parseFile::process_SpecialCommand(std::ostream &outputFile,
                                       const map<variable,string> &vnames,
                                       int &openbrace) {
  is.get() ; // get leading [
  string name = get_name(is) ;
  if(is.peek() != ']') {
    cerr << "expecting ']' to close special command '" << name << "'" << endl ;
    throw parseError("syntax error") ;
  }
  is.get() ;

  int nsz = name.size() ;
  for(int i=0;i<nsz;++i)
    if(name[i] >= 'A' || name[i] <= 'Z')
      name[i] = std::tolower(name[i]) ;
  
  if(name == "once") {
    killsp() ;
    if(is.peek() != '{') {
      cerr << "expecting '{' after $[Once] command" << endl ;
      cerr << "found " << char(is.peek()) << " instead." <<endl ;
      throw parseError("syntax error") ;
    }
    outputFile << "if(Loci::is_leading_execution()) " ;

  } else if(name == "atomic") {
    killsp() ;
    if(is.peek() != '{') {
      cerr << "expecting '{' after $[Atomic] command" << endl ;
      cerr << "found " << char(is.peek()) << " instead." <<endl ;
      throw parseError("syntax error") ;
    }
    is.get() ;
    openbrace++ ;
    outputFile << "{ Loci::atomic_region_helper L__ATOMIC_REGION ; " << endl ;
  } else {
    cerr << "unknown special command '[" << name << "]' !" << endl ;
    throw parseError("syntax error") ;
  }
}

void parseFile::process_Prelude(std::ostream &outputFile,
                                const map<variable,string> &vnames) {
  if(is.peek() == '{')  // eat open brace
     is.get() ;
  else
    throw parseError("expected open brace") ;
  
  outputFile << "    virtual void prelude(const Loci::sequence &seq) { "
             << endl ;
  syncFile(outputFile) ;
    
  
  int openbrace = 1 ;
  for(;;) {
    killspout(outputFile) ;
    if(is.peek() == EOF)
      throw parseError("unexpected EOF") ;
      
    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      
      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      openbrace++ ;
      continue ;
    }
    if(is.peek() == '$') {
      string name ;
      variable v ;
      is.get() ;
      if(is.peek() == '[') {
        process_SpecialCommand(outputFile,vnames,openbrace) ;
        continue ;
      } 
      var vin ;
      vin.get(is) ;
      line_no += vin.num_lines() ;
      v = variable(vin.str()) ;
        
      auto vmi = vnames.find(v) ;
      if(vmi == vnames.end()) {
        cerr << "variable " << v << " is unknown to this rule!" << endl ;
        throw parseError("type error") ;
      }
      outputFile << vmi->second  ;
    }
  
    char c = is.get() ;
    if(c == '\n')
      line_no++ ;
    outputFile << c ;
  } ;
}

void parseFile::process_Compute(std::ostream &outputFile,
                                const map<variable,string> &vnames) {
  outputFile << "    void compute(const Loci::sequence &seq) { " ;
  is.get() ;
  
  int openbrace = 1 ;
  for(;;) {
    killspout(outputFile) ;
    if(is.peek() == EOF)
      throw parseError("unexpected EOF") ;
      
    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      
      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      openbrace++ ;
      continue ;
    }
    if(is.peek() == '"') {
      is.get() ;
      outputFile << '"' ;
      while(is.peek() != '"' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '"' ;
      continue ;
    }
    if(is.peek() == '\'') {
      is.get() ;
      outputFile << '\'' ;
      while(is.peek() != '\'' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '\'' ;
      continue ;
    }      
    if(is.peek() == '$') {
      variable v ;
      is.get() ;
      if(is.peek() == '[') {
        process_SpecialCommand(outputFile,vnames,openbrace) ;
        continue ;
      }
      bool deref = true ;
      if(is.peek() == '*') {
        is.get() ;
        deref = false ;
      }
      

      var vin ;
      vin.get(is) ;
      v = variable(vin.str()) ;
      line_no += vin.num_lines() ;

      auto vmi = vnames.find(v) ;
      if(vmi == vnames.end()) {
        cerr << "variable " << v << " is unknown to this rule!" << endl ;
        throw parseError("type error") ;
      }
      auto mi = lookupVarType(v) ;
      if(checkTypeValid(mi) &&
               (mi->second.container == "Constraint" || !deref)) {
        outputFile << vmi->second ;
      } else {
        outputFile << "(*" << vmi->second << ')' ;
      }
      
    }
    char c = is.get() ;
    if(c == '\n')
      line_no++ ;
    outputFile << c ;
  } 
}

string getNumber(std::istream &is) {
  string num ;
  while(isdigit(is.peek()))
    num+= is.get();
  if(is.peek() == '.') {
    num += is.get() ;
    while(isdigit(is.peek()))
      num+= is.get();
  }
  if(is.peek() == 'e' || is.peek() == 'E') {
    num += is.get() ;
    if(is.peek() == '-' || is.peek() == '+')
      num += is.get() ;
    while(isdigit(is.peek()))
      num += is.get() ;
  }
  return num ;
}

string parseFile::process_String(string in,
                                 const map<variable,string> &vnames,
				 const set<list<variable> > &validate_set) {
  ostringstream outputFile ;
  istringstream is(in) ;

  int line_no = 0 ;

  for(;;) {
    ::killspOut(is,line_no,outputFile) ;

    if(is.peek() == EOF)
      break ;
      
    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      continue ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      continue ;
    }
    if(is.peek() == '"') {
      is.get() ;
      outputFile << '"' ;
      while(is.peek() != '"' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '"' ;
      continue ;
    }
    if(is.peek() == '\'') {
      is.get() ;
      outputFile << '\'' ;
      while(is.peek() != '\'' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '\'' ;
      continue ;
    }      
    if(is.peek() == '/') {
      is.get() ;
      outputFile << '/' ;
      if(is.peek() == '/') { // comment line
        is.get() ;
        outputFile << '/' ;
        while(is.peek() != '\n') {
          char c = is.get() ;
          outputFile << c ;
        }
        ::killspOut(is,line_no,outputFile) ;
      }
      continue ;
    }
          
    if(is.peek() == '#') {
      is.get() ;
      outputFile << '#' ;
      while(is.peek() != '\n') {
        char c = is.get() ;
        outputFile << c ;
      }
      ::killspOut(is,line_no,outputFile) ;
      continue ;
    }

    if(isdigit(is.peek())) {
      outputFile << getNumber(is) ;
      continue ;
    }

    if(is_name(is) || is.peek() == '$') {
      bool first_name = is_name(is) ;
      string name ;
      variable v ;
      string brackets ;
      if(first_name) 
        name = get_name(is) ;
      else {
        is.get() ;
        if(is.peek() == '*') {
          is.get() ;
          var vin ;
          vin.get(is) ;
          line_no += vin.num_lines() ;
          variable v(vin.str()) ;
          auto  vmi = vnames.find(v) ;
          if(vmi == vnames.end()) {
            cerr << "variable " << v << " is unknown to this rule!" << endl ;
            throw parseError("type error") ;
          }
          
          outputFile << vmi->second ;
          continue ;
        }
        
        var vin ;
        vin.get(is) ;
        line_no += vin.num_lines() ;
        v = variable(vin.str()) ;
        ::killspOut(is,line_no,outputFile) ;
        if(is.peek() == '[') {
          nestedbracketstuff nb ;
          nb.get(is) ;
          string binfo = process_String(nb.str(),vnames,validate_set) ;
          brackets = "[" + binfo + "]" ;
        }
      }
      list<variable> vlist ;
      list<string> blist ;
      bool dangling_arrow = false ;

      for(;;) { // scan for ->$ chain
        ::killspOut(is,line_no,outputFile) ;
        if(is.peek() != '-')
          break ;
        char c=is.get() ;
        if(c== '-' && is.peek() == '>') {
          c=is.get() ;
          ::killspOut(is,line_no,outputFile) ;
          if(is.peek() == '$') {
            is.get() ;
            var vin ;
            vin.get(is) ;
            vlist.push_back(variable(vin.str())) ;
            line_no += vin.num_lines() ;
            string brk ;
            ::killspOut(is,line_no,outputFile) ;
            if(is.peek() == '[') {
              nestedbracketstuff nb ;
              nb.get(is) ;
              string binfo = process_String(nb.str(),vnames,validate_set) ;
              brk = "[" + binfo +"]";
            }
            blist.push_back(brk) ;
          } else {
            dangling_arrow = true ;
            break ;
          }
        } else {
          is.unget() ;
          break ;
        }
      }
      if(dangling_arrow && first_name) {
        outputFile << name << " ->" ;
        continue ;
      }
      if(dangling_arrow)
        throw parseError("syntax error, near '->' operator") ;

      validate_VariableAccess(v,vlist,first_name,vnames,validate_set) ;
      
      if(first_name && (vlist.size() == 0)) {
        outputFile << name << ' ' ;
        continue ;
      }
      for(auto ri=vlist.rbegin();ri!=vlist.rend();++ri) {
        auto vmi = vnames.find(*ri) ;
        if(vmi == vnames.end()) {
          cerr << "variable " << *ri << " is unknown to this rule!" << endl ;
          throw parseError("type error") ;
        }
        outputFile << vmi->second << '[' ;
      }
      if(first_name) {
        outputFile << '*' << name ;
      } else {
        auto vmi = vnames.find(v) ;
        if(vmi == vnames.end()) {
          cerr << "variable " << v << " is unknown to this rule!" << endl ;
          throw parseError("type error: is this variable in the rule signature?") ;
        }
        if(prettyOutput)
          outputFile << vmi->second << "[e]" ;
        else
          outputFile << vmi->second << "[_e_]" ;
      }

      outputFile << brackets ;
      for(auto rbi=blist.begin();rbi!=blist.end();++rbi) {
        outputFile << ']' << *rbi ;
      }

    }
    if(is.peek() != EOF) {
      char c = is.get() ;
      outputFile << c ;
    }
  } 

  
  return outputFile.str() ;
}


void parseFile::validate_VariableAccess(variable v, const list<variable> &vlist,
					bool first_name,
					const map<variable,string> &vnames,
					const set<list<variable> > &validate_set) {

  list<variable> vlistall ;
  for(auto vitmp=vlist.begin();vitmp!=vlist.end();++vitmp) {
    variable vt = *vitmp ;
    while(vt.get_info().priority.size() != 0)
      vt = vt.drop_priority() ;
    vlistall.push_back(vt) ;
  }
  variable vt = v ;
  while(vt.get_info().priority.size() != 0)
    vt = vt.drop_priority() ;
  vlistall.push_front(vt) ;
  
  if(!first_name && !vlistall.empty()
     && validate_set.find(vlistall) == validate_set.end()) {
    ostringstream msg ;
    msg << "variable access " ;
    for(auto lvi=vlistall.begin();lvi!=vlistall.end();) {
      msg << *lvi ;
      ++lvi ;
      if(lvi!=vlistall.end())
	msg << "->" ;
    }
    msg << " not consistent with rule signature!" ;
    throw parseError(msg.str()) ;
  }
  
  for(auto ri=vlist.rbegin();ri!=vlist.rend();++ri) {
    auto vmi = vnames.find(*ri) ;
    if(vmi == vnames.end()) {
      cerr << "variable " << *ri << " is unknown to this rule!" << endl ;
      throw parseError("type error") ;
    }
  }
  if(!first_name) {
    auto vmi = vnames.find(v) ;
    if(vmi == vnames.end()) {
      cerr << "variable " << v << " is unknown to this rule!" << endl ;
      throw parseError("type error: is this variable in the rule signature?") ;
    }
  }
}




void parseFile::process_Calculate(std::ostream &outputFile,
                                  const map<variable,string> &vnames,
                                  const set<list<variable> > &validate_set) {
  if(prettyOutput)
    outputFile << "    void calculate(Loci::Entity e) { " << endl ;
  else
    outputFile << "    void calculate(Loci::Entity _e_) { " << endl ;
  is.get() ;
  while(is.peek() == ' ' || is.peek() == '\t')
    is.get() ;
  if(is.peek() == '\n') {
    is.get() ;
    line_no++ ;
  }
  syncFile(outputFile) ;
  killspout(outputFile) ;
  int openbrace = 1 ;
  for(;;) {
    killspout(outputFile) ;
    if(is.peek() == EOF)
      throw parseError("unexpected EOF in process_Calculate") ;
      
    if(is.peek() == '}') {
      is.get() ;
      outputFile << '}' ;
      
      openbrace-- ;
      if(openbrace == 0)
        break ;
    }
    if(is.peek() == '{') {
      is.get() ;
      outputFile << '{' ;
      openbrace++ ;
      continue ;
    }
    if(is.peek() == '"') {
      is.get() ;
      outputFile << '"' ;
      while(is.peek() != '"' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '"' ;
      continue ;
    }
    if(is.peek() == '\'') {
      is.get() ;
      outputFile << '\'' ;
      while(is.peek() != '\'' && is.peek() != EOF) {
        char c = is.get() ;
        outputFile << c ;
      }
      is.get() ;
      outputFile << '\'' ;
      continue ;
    }      
    if(is.peek() == '/') {
      is.get() ;
      outputFile << '/' ;
      if(is.peek() == '/') { // comment line
        is.get() ;
        outputFile << '/' ;
        while(is.peek() != '\n') {
          char c = is.get() ;
          outputFile << c ;
        }
        killspout(outputFile) ;
      }
      continue ;
    }
          
    if(is.peek() == '#') {
      is.get() ;
      outputFile << '#' ;
      while(is.peek() != '\n') {
        char c = is.get() ;
        outputFile << c ;
      }
      killspout(outputFile) ;
      continue ;
    }

    if(isdigit(is.peek())) {
      outputFile << getNumber(is) ;
      continue ;
    }
    
    if(is_name(is) || is.peek() == '$') {
      int lcount = 0 ;
      bool first_name = is_name(is) ;
      if(!first_name) {
        is.get() ;
        if(is.peek() == '[') { // special command
          process_SpecialCommand(outputFile,vnames,openbrace) ;
          continue ;
        }
      }
      string name ;
      variable v ;
      string brackets ;
      if(first_name) 
        name = get_name(is) ;
      else {
        if(is.peek() == '*') {
          is.get() ;
          var vin ;
          vin.get(is) ;
          line_no += vin.num_lines() ;
          lcount += vin.num_lines() ;
          
          variable v(vin.str()) ;
          auto vmi = vnames.find(v) ;
          if(vmi == vnames.end()) {
            cerr << "variable " << v << " is unknown to this rule!" << endl ;
            throw parseError("type error") ;
          }
          
          outputFile << vmi->second ;
          continue ;
        }
        
        var vin ;
        vin.get(is) ;
        line_no += vin.num_lines() ;
        lcount += vin.num_lines();
        v = variable(vin.str()) ;
        killsp() ;
        if(is.peek() == '[') {
          nestedbracketstuff nb ;
          nb.get(is) ;
          string binfo = process_String(nb.str(),vnames,validate_set) ;
          brackets = "[" + binfo + "]" ;
          line_no += nb.num_lines() ;
          lcount += nb.num_lines() ;
        }
      }
      list<variable> vlist ;
      list<string> blist ;
      bool dangling_arrow = false ;

      for(;;) { // scan for ->$ chain
        lcount += killsp() ;
        if(is.peek() != '-')
          break ;
        char c=is.get() ;
        if(c== '-' && is.peek() == '>') {
          c=is.get() ;
          lcount += killsp() ;
          if(is.peek() == '$') {
            is.get() ;
            var vin ;
            vin.get(is) ;
            line_no += vin.num_lines() ;
            vlist.push_back(variable(vin.str())) ;
            string brk ;
            lcount += killsp() ;
            if(is.peek() == '[') {
              nestedbracketstuff nb ;
              nb.get(is) ;
              string binfo = process_String(nb.str(),vnames,validate_set) ;
              brk = "[" + binfo +"]";
              line_no += nb.num_lines() ;
              lcount += nb.num_lines() ;
            }
            blist.push_back(brk) ;
          } else {
            dangling_arrow = true ;
            break ;
          }
        } else {
          is.unget() ;
          break ;
        }
      }
      if(dangling_arrow && first_name) {
        outputFile << name << " ->" ;
        continue ;
      }
      if(dangling_arrow)
        throw parseError("syntax error, near '->' operator") ;

      if(first_name && (vlist.empty())) {
        outputFile << name << ' ' ;
        continue ;
      }

      validate_VariableAccess(v,vlist,first_name,vnames,validate_set) ;

      for(auto ri=vlist.rbegin();ri!=vlist.rend();++ri) {
        auto vmi = vnames.find(*ri) ;
        if(vmi == vnames.end()) {
          cerr << "variable " << *ri << " is unknown to this rule!" << endl ;
          throw parseError("type error") ;
        }
        outputFile << vmi->second << '[' ;
      }
      if(first_name) {
        outputFile << '*' << name ;
      } else {
        auto vmi = vnames.find(v) ;
        if(vmi == vnames.end()) {
          cerr << "variable " << v << " is unknown to this rule!" << endl ;
          throw parseError("type error: is this variable in the rule signature?") ;
	}
        if(prettyOutput)
          outputFile << vmi->second << "[e]" ;
        else
          outputFile << vmi->second << "[_e_]" ;
      }

      outputFile << brackets ;
      for(auto rbi=blist.begin();rbi!=blist.end();++rbi) {
        outputFile << ']' << *rbi ;
      }

      for(int i=0;i<lcount;++i)
        outputFile << endl ;
      continue ;
      
    }
    
    char c = is.get() ;
    if(c == '\n')
      line_no++ ;
    outputFile << c ;

  }
}

/// Visitor that prints an AST using a simple substitution map
class AST_printTree : public AST_visitor {
 public:
  ostream &out ;
  int indent_level ;
  void indent() {
    for(int i=0;i<indent_level;++i)
      out << "  " ;
  }
  void pushindent(AST_type &s) { out << endl ;
    indent() ;
    out << "[[" << OPtoName(s) << " ";
    indent_level++ ; }
  void popindent() { indent_level-- ; out <<"]]"<< endl ; indent() ; }
  AST_printTree(ostream &s): out(s),indent_level(0) {} 
  virtual void visit(AST_exprOper &)  ;
  virtual void visit(AST_Token &) ;
  virtual void visit(AST_Block &) ;
  virtual void visit(AST_typeSpec &) ;
  virtual void visit(AST_declaration &) ;
  virtual void visit(AST_SimpleStatement &) ;
  virtual void visit(AST_controlStatement &) ;
} ;

void AST_printTree::visit(AST_exprOper &s) {
  using namespace nodeTypes ;
  switch (s.nodeType) {
  case OP_GROUP:

    pushindent(s) ;
    out << '(' ;
    for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();++ii) {
      if(*ii != 0)
	(*ii)->accept(*this) ;
    }
    out << ')' ;
    popindent() ;
    
    break ;
  case OP_CAST:
    {
      pushindent(s) ;
      out << '(' ;
      AST_type::ASTList::iterator ii = s.terms.begin() ;
      if(ii != s.terms.end() && *ii != 0)
        (*ii)->accept(*this) ;
      out << ')' ;
      ++ii ;
      if(ii != s.terms.end() && *ii != 0)
        (*ii)->accept(*this) ;
      popindent() ;
    }
    break ;
  case OP_TEMPLATE_CAST:
    pushindent(s) ;
    for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();++ii) {
      if(*ii != 0)
	(*ii)->accept(*this) ;
    }
    popindent() ;
    break ;
  case OP_BRACEBLOCK:
    {
      pushindent(s) ;
      out << '{' ;
      for(AST_type::ASTList::iterator ii = s.terms.begin();
          ii != s.terms.end(); ++ii)
        if(*ii != 0)
          (*ii)->accept(*this) ;
      out << '}' ;
      popindent() ;
    }
    break ;
  case OP_FUNC:
    {
      AST_type::ASTList::iterator ii=s.terms.begin() ;
      pushindent(s) ;
      FATAL(ii == s.terms.end()) ;
      out << "[[" ;
      (*ii)->accept(*this) ;
      out << "]]" ;
      ++ii ;
      out << '(' ;
      out << "[[" ;
      if(ii != s.terms.end()) {
	(*ii)->accept(*this) ;
        ++ii ;
      }
      out << "]]" ;
      out << ')' ;
      popindent() ;
      
      if(ii!=s.terms.end()) {
	cerr << "internal error processing func" ;
        (*ii)->accept(*this) ;
      }
    }
    break ;
  case OP_TEMPLATE:
    {
      AST_type::ASTList::iterator ii=s.terms.begin() ;
      pushindent(s) ;
      FATAL(ii == s.terms.end()) ;
      out << "[[" ;
      (*ii)->accept(*this) ;
      out << "]]" ;
      ++ii ;
      out << '<' ;
      out << "[[" ;
      if(ii != s.terms.end()) {
	(*ii)->accept(*this) ;
        ++ii ;
      }
      out << "]]" ;
      out << '>' ;
      popindent() ;
      if(ii!=s.terms.end()) {
	cerr << "internal error processing func" ;
        (*ii)->accept(*this) ;
      }
    }
    break ;
  case OP_ARRAY:
    {
      AST_type::ASTList::iterator ii=s.terms.begin() ;
      pushindent(s) ;
      if(*ii != 0)
	(*ii)->accept(*this) ;
      ++ii ;
      out << '[' ;
      if(*ii != 0)
	(*ii)->accept(*this) ;
      out << ']' ;
      popindent() ;
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
      pushindent(s) ;
      if(*ii != 0)
	(*ii)->accept(*this) ;
      ++ii ;
      out << '?' ;
      if(*ii != 0)
	(*ii)->accept(*this) ;
      popindent() ;
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
      pushindent(s) ;
      for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();++ii)
	if(*ii != 0)
	  (*ii)->accept(*this) ;
      popindent() ;
    }
    break ;
  case OP_POSTINCREMENT:
  case OP_POSTDECREMENT:
    {
      pushindent(s) ; 
      for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();++ii)
	if(*ii != 0)
	  (*ii)->accept(*this) ;
      popindent() ;
    }
    break ;
  default:
    {
      pushindent(s) ;
      for(AST_type::ASTList::iterator ii=s.terms.begin();ii!=s.terms.end();) {
	if(*ii != 0)
	  (*ii)->accept(*this) ;
	++ii ;
      }
      popindent() ;
    }

    break ;
  }
}

void AST_printTree::visit(AST_Token &s) {
  using namespace nodeTypes ;
  if(ASTEqual(s,TK_LOCI_DIRECTIVE)) {
    out << "$[" << s.text << "] " ;
  } else if(ASTEqual(s,TK_LOCI_CONTAINER)) {
    out << "$*" << s.text << " " ;
  } else if(ASTEqual(s,TK_LOCI_VARIABLE)) {
    out << "$" << s.text  << " " ;
  } else if(ASTEqual(s,TK_MACRO)) {
    out << "#" << s.text << endl ;
  } else 
    out <<s.text << ' ' ;
}

void AST_printTree::visit(AST_Block &s) {
  pushindent(s) ;
  for(auto ii=s.elements.begin();ii!=s.elements.end();++ii)
    if(*ii!=0)
      (*ii)->accept(*this) ;
  popindent() ;
}

void AST_printTree::visit(AST_typeSpec &s) {
  pushindent(s) ;
  out << "[[" ;
  for(auto ii=s.type_spec.begin();ii!=s.type_spec.end();++ii)
    if(*ii != 0)
      (*ii)->accept(*this) ;
  out << "]]" ;
  popindent() ;
}

void AST_printTree::visit(AST_declaration &s) {
  pushindent(s) ;
  out << "[[" ;
  for(auto ii=s.type_decl.begin();ii!=s.type_decl.end();++ii)
    if(*ii != 0)
      (*ii)->accept(*this) ;
  out << "]][[" ;
  for(auto ii=s.decls.begin();ii!=s.decls.end();++ii)
    if(*ii != 0)
      (*ii)->accept(*this) ;
  out << "]]" ;
  popindent() ;
}

void AST_printTree::visit(AST_SimpleStatement &s) {
  pushindent(s) ;
  if(s.exp!=0)
    s.exp->accept(*this) ;
  if(s.Terminal!=0) 
    s.Terminal->accept(*this) ;
  popindent() ;
}
void AST_printTree::visit(AST_controlStatement &s) {
  pushindent(s) ;
  s.controlType->accept(*this) ;
  for(auto ii=s.parts.begin();ii!=s.parts.end();++ii) {
    if(*ii != 0)
      (*ii)->accept(*this) ;
  }
  popindent() ;
}

class AST_editLociMapArrayAccess : public AST_visitor {
public:
  virtual void visit(AST_exprOper &) ;
} ;  

void AST_editLociMapArrayAccess::visit(AST_exprOper &op) {
  using namespace nodeTypes ;
  
  const int sz = op.terms.size() ;
  // Check to see if this is a Loci mapping operator that
  // ends in an array. In this case, the arrow needs to bind
  // most tightly, so rearrange expression tree
  if(sz>0 &&
     ASTEqual(op,OP_ARROW) && 
     ASTEqual(op.terms[sz-1],OP_ARRAY)) {

    CPTR<AST_exprOper> last = CPTR<AST_exprOper>(op.terms[sz-1]) ;
    // rearrange tree so that array operator is moved to the top
    // and the mapping operator applies to the last variable
    std::swap(op.nodeType,last->nodeType) ;
    std::swap(op.terms,last->terms) ;
    std::swap(last->terms[sz-1],op.terms[0]) ;
  }
  // Now visit children
  for(size_t i=0;i<op.terms.size();++i) {
    if(op.terms[i] !=0)
      op.terms[i]->accept(*this) ;
  }
}

class AST_editLociVariableAccess : public AST_visitor {
public:
  const std::map<variable,std::string> &vnames ;
  AST_type::ASTP entityIndex ;
  AST_type::ASTP convertLociVar(AST_type::ASTP var) {
    CPTR<AST_Token> p = CPTR<AST_Token>(var) ;
    variable v(p->text) ;
    //    while(v.get_info().priority.size() != 0)
    //      v = v.drop_priority() ;
    
    auto vmi = vnames.find(v) ;
    if(vmi == vnames.end()) {
      cerr << "variable " << v << " is unknown to this rule!" << endl ;
      throw parseError("type error: is this variable in the rule signature?") ;
    }
    CPTR<AST_Token> np = new AST_Token ;
    np->lineno = p->lineno ;
    np->text = vmi->second ;
    np->nodeType = nodeTypes::TK_NAME ;
    return AST_type::ASTP(np) ;
  }
  AST_type::ASTP arrayAccess(AST_type::ASTP var, AST_type::ASTP index) {
    CPTR<AST_exprOper> e = new AST_exprOper ;
    e->nodeType = nodeTypes::OP_ARRAY ;
    e->terms.push_back(var) ;
    e->terms.push_back(index) ;
    return AST_type::ASTP(e) ;
  }
  AST_type::ASTP addEntityIndex(AST_type::ASTP var) {
    return arrayAccess(var,entityIndex) ;
  }
  
  AST_editLociVariableAccess(const std::map<variable,std::string> &vnames_in):
    vnames(vnames_in) {
    CPTR<AST_Token> e = new AST_Token ;
    e->lineno = -1 ;
    e->text = "_e_" ;
    e->nodeType = nodeTypes::TK_NAME ;
    entityIndex = AST_type::ASTP(e) ;
  }
  virtual void visit(AST_exprOper &) ;
} ;



void AST_editLociVariableAccess::visit(AST_exprOper &op) {
  using namespace nodeTypes ;
  
  const int sz = op.terms.size() ;
  if(op.nodeType == OP_ARROW) {
    // Check to see if this is a Loci mapping operator
    if(ASTEqual(op.terms[sz-1],TK_LOCI_VARIABLE)) {
      // It is so we need to edit create a tree of array accessor operations
      // First create the root of the tree which starts at the beginning
      CPTR<AST_exprOper> rootptr = new AST_exprOper ;
      if(ASTEqual(op.terms[0], TK_NAME)) {
        // This is the special case of a pointer type (sometimes used to
        // iterate over multiMaps (may need to be deprecated in the future
        // as this exposes the memory layout of the multiMap data structure
        // which may need to change on GPGPUs
        rootptr->nodeType = OP_STAR ;
        rootptr->terms.push_back(op.terms[0]) ;
      } else if(ASTEqual(op.terms[0],TK_LOCI_VARIABLE)) {
        // base map just add entity index operator
        rootptr->nodeType = OP_ARRAY ;
        rootptr->terms.push_back(convertLociVar(op.terms[0])) ;
        rootptr->terms.push_back(entityIndex) ;
      } else if(ASTEqual(op.terms[0],OP_ARRAY)) {
        // base map is a multiMap, still need to insert the entity index
        // operator
        CPTR<AST_exprOper> mapaccess= CPTR<AST_exprOper>(op.terms[0]) ;
        if(mapaccess->terms.size() != 2 ||
           mapaccess->terms[0]->nodeType != TK_LOCI_VARIABLE) {
          cerr << "invalid map at base of Loci mapping operator" << endl;
          throw parseError("invalid map at base of Loci mapping operator") ;
        }
        rootptr->nodeType = OP_ARRAY ;
        AST_type::ASTP p = addEntityIndex(convertLociVar(mapaccess->terms[0])) ;
        rootptr->terms.push_back(p) ;
        rootptr->terms.push_back(mapaccess->terms[1]) ;
      } else {
        cerr << "syntax error in Loci mapping operator" << endl ;
        throw parseError("Invalid Loci mapping operator") ;
      }

      // Now we have the root pointer start building the access tree
      for(int i=1;i<sz;++i) {
        CPTR<AST_exprOper> newroot = 0 ;
        if(ASTEqual(op.terms[i],TK_LOCI_VARIABLE)) {
          newroot = CPTR<AST_exprOper>(arrayAccess(convertLociVar(op.terms[i]),
                                                   AST_type::ASTP(rootptr))) ;
        } else if(ASTEqual(op.terms[i],OP_ARRAY)) {
          // base map is a multiMap, still need to insert the entity index
          // operator
          CPTR<AST_exprOper> mapaccess= CPTR<AST_exprOper>(op.terms[i]) ;

          if(mapaccess->terms.size() != 2 ||
             mapaccess->terms[0]->nodeType != TK_LOCI_VARIABLE) {
            cerr << "invalid map at base of Loci mapping operator" << endl;
            throw parseError("invalid map at base of Loci mapping operator") ;
          }
          AST_type::ASTP var = convertLociVar(mapaccess->terms[0]) ;
          newroot =
            CPTR<AST_exprOper>(arrayAccess(arrayAccess(var,
                                                       AST_type::ASTP(rootptr)),
                                           mapaccess->terms[1])) ;
        } else {
          cerr << "invalid Loci mapping operator" << endl ;
          throw parseError("Invalid Loci mapping operator") ;
        }
        if(newroot != 0)
          rootptr = newroot ;
      }
      op.nodeType = rootptr->nodeType ;
      op.terms = rootptr->terms ;
    }
  }
  for(size_t i=0;i<op.terms.size();++i) {
    if(ASTEqual(op.terms[i],TK_LOCI_VARIABLE)) {
      op.terms[i] = addEntityIndex(convertLociVar(op.terms[i])) ;
    } else if(ASTEqual(op.terms[i],TK_LOCI_CONTAINER)) {
      op.terms[i] = convertLociVar(op.terms[i]); 
    } else {
      op.terms[i]->accept(*this) ;
    }
  }
  
}

class AST_editLociVariableAccess2 : public AST_visitor {
public:
  const std::map<variable,std::string> &vnames ;
  const std::map<variable,std::string> &vtypes ;
  AST_type::ASTP entityIndex ;

  AST_type::ASTP convertLociVar(AST_type::ASTP var) {
    CPTR<AST_Token> p = CPTR<AST_Token>(var) ;
    variable v(p->text) ;
    //    while(v.get_info().priority.size() != 0)
    //      v = v.drop_priority() ;
    
    auto vmi = vnames.find(v) ;
    if(vmi == vnames.end()) {
      cerr << "variable " << v << " is unknown to this rule!" << endl ;
      throw parseError("type error: is this variable in the rule signature?") ;
    }
    CPTR<AST_Token> np = new AST_Token ;
    np->lineno = p->lineno ;
    np->text = vmi->second ;
    np->nodeType = nodeTypes::TK_NAME ;
    return AST_type::ASTP(np) ;
  }

  AST_type::ASTP arrayAccess(AST_type::ASTP var, AST_type::ASTP index) {
    CPTR<AST_exprOper> e = new AST_exprOper ;
    e->nodeType = nodeTypes::OP_ARRAY ;
    e->terms.push_back(var) ;
    e->terms.push_back(index) ;
    return AST_type::ASTP(e) ;
  }

  AST_type::ASTP addEntityIndex(AST_type::ASTP var) {
    return arrayAccess(var,entityIndex) ;
  }
  
  AST_editLociVariableAccess2(
                              const std::map<variable,std::string> &vnames_in,
                              const std::map<variable,std::string> &vtypes_in):
    vnames(vnames_in), vtypes(vtypes_in) {
    CPTR<AST_Token> e = new AST_Token ;
    e->lineno = -1 ;
    e->text = "_e_" ;
    e->nodeType = nodeTypes::TK_NAME ;
    entityIndex = AST_type::ASTP(e) ;
  }

  virtual void visit(AST_exprOper &) ;
} ;



void AST_editLociVariableAccess2::visit(AST_exprOper &op) {
  using namespace nodeTypes ;
  
  const int sz = op.terms.size() ;
  if(op.nodeType == OP_ARROW) {
    // Check to see if this is a Loci mapping operator
    if(ASTEqual(op.terms[sz-1],TK_LOCI_VARIABLE)) {
      // It is so we need to edit create a tree of array accessor operations
      // First create the root of the tree which starts at the beginning
      CPTR<AST_exprOper> rootptr = new AST_exprOper ;
      if(ASTEqual(op.terms[0], TK_NAME)) {
        // This is the special case of a pointer type (sometimes used to
        // iterate over multiMaps (may need to be deprecated in the future
        // as this exposes the memory layout of the multiMap data structure
        // which may need to change on GPGPUs
        rootptr->nodeType = OP_STAR ;
        rootptr->terms.push_back(op.terms[0]) ;
      } else if(ASTEqual(op.terms[0],TK_LOCI_VARIABLE)) {
        // base map just add entity index operator
        rootptr->nodeType = OP_ARRAY ;
        rootptr->terms.push_back(convertLociVar(op.terms[0])) ;
        rootptr->terms.push_back(entityIndex) ;
      } else if(ASTEqual(op.terms[0],OP_ARRAY)) {
        // base map is a multiMap, still need to insert the entity index
        // operator
        CPTR<AST_exprOper> mapaccess= CPTR<AST_exprOper>(op.terms[0]) ;
        if(mapaccess->terms.size() != 2 ||
           mapaccess->terms[0]->nodeType != TK_LOCI_VARIABLE) {
          cerr << "invalid map at base of Loci mapping operator" << endl;
          throw parseError("invalid map at base of Loci mapping operator") ;
        }
        rootptr->nodeType = OP_ARRAY ;
        AST_type::ASTP p = addEntityIndex(convertLociVar(mapaccess->terms[0])) ;
        rootptr->terms.push_back(p) ;
        rootptr->terms.push_back(mapaccess->terms[1]) ;
      } else {
        cerr << "syntax error in Loci mapping operator" << endl ;
        throw parseError("Invalid Loci mapping operator") ;
      }

      // Now we have the root pointer start building the access tree
      for(int i=1;i<sz;++i) {
        CPTR<AST_exprOper> newroot = 0 ;
        if(ASTEqual(op.terms[i],TK_LOCI_VARIABLE)) {
          newroot = CPTR<AST_exprOper>(arrayAccess(convertLociVar(op.terms[i]),
                                                   AST_type::ASTP(rootptr))) ;
        } else if(ASTEqual(op.terms[i],OP_ARRAY)) {
          // base map is a multiMap, still need to insert the entity index
          // operator
          CPTR<AST_exprOper> mapaccess= CPTR<AST_exprOper>(op.terms[i]) ;

          if(mapaccess->terms.size() != 2 ||
             mapaccess->terms[0]->nodeType != TK_LOCI_VARIABLE) {
            cerr << "invalid map at base of Loci mapping operator" << endl;
            throw parseError("invalid map at base of Loci mapping operator") ;
          }
          AST_type::ASTP var = convertLociVar(mapaccess->terms[0]) ;
          newroot =
            CPTR<AST_exprOper>(arrayAccess(arrayAccess(var,
                                                       AST_type::ASTP(rootptr)),
                                           mapaccess->terms[1])) ;
        } else {
          cerr << "invalid Loci mapping operator" << endl ;
          throw parseError("Invalid Loci mapping operator") ;
        }
        if(newroot != 0)
          rootptr = newroot ;
      }
      op.nodeType = rootptr->nodeType ;
      op.terms = rootptr->terms ;
    }
  }

  for(size_t i=0;i<op.terms.size();++i) {
    if(ASTEqual(op.terms[i],TK_LOCI_VARIABLE)) {
      bool is_param = false ;
      CPTR<AST_Token> tok(op.terms[i]) ;
      variable v(tok->text) ;
      auto t = vtypes.find(v) ;
      if(t != vtypes.end()) {
        if(t->second == "param") {
          is_param = true ;
        }
      }

      if(is_param) {
        CPTR<AST_exprOper> param_access = new AST_exprOper ;
        param_access->nodeType = nodeTypes::OP_STAR ;
        param_access->terms.push_back(convertLociVar(op.terms[i])) ;

        CPTR<AST_exprOper> param_group = new AST_exprOper ;
        param_group->nodeType = nodeTypes::OP_GROUP ;
        param_group->terms.push_back(AST_type::ASTP(param_access)) ;

        op.terms[i] = AST_type::ASTP(param_group) ;
      } else {
        op.terms[i] = addEntityIndex(convertLociVar(op.terms[i])) ;
      }
    } else if(ASTEqual(op.terms[i],TK_LOCI_CONTAINER)) {
      op.terms[i] = convertLociVar(op.terms[i]);
    } else {
      op.terms[i]->accept(*this) ;
    }
  }
}

void parseFile::process_Calculate2(std::ostream &outputFile,
                                   const map<variable,string> &vnames,
                                   const set<list<variable> > &validate_set,
                                   const parseSharedInfo &parseInfo) {
  varmap typemap ;
  typemap["cerr"] = localIdentifier() ;
  typemap["std::cerr"] = localIdentifier() ;
  typemap["cout"] = localIdentifier() ;
  typemap["std::cout"] = localIdentifier() ;
  typemap["debugout"] = localIdentifier() ;
  typemap["Loci::debugout"] = localIdentifier() ;
      
  if(is.peek() != '{')
    throw parseError("syntax error, expecting '{'") ;
      
      
  CPTR<AST_type> ap = parseBlock(is,line_no,filename,typemap) ;

  AST_condenseLeftAssociative condenseOps ;
  ap->accept(condenseOps) ;

  // This is sort of a hack because the precedence of the mapping
  // operator (->) is context senstive
  AST_editLociMapArrayAccess mapEditOps ;
  ap->accept(mapEditOps) ;
  
  if(parseInfo.diag_level > 0) {
    AST_printTree diagout(cerr) ;
    ap->accept(diagout) ;
  }

  AST_errorCheck syntaxChecker ;
  ap->accept(syntaxChecker) ;
  if(syntaxChecker.hasErrors()) {
#ifdef VERBOSE
    AST_simplePrint printer(cerr,-1,false) ;
    ap->accept(printer) ;
#endif
    throw parseError("syntax error") ;
  }

  //  cerr << "vnames = " << endl ;
  //  for(auto ii=vnames.begin();ii!=vnames.end();++ii) {
  //    cerr << ii->first << " " << ii->second << endl ;
  //  }
  AST_editLociVariableAccess AST_editor(vnames) ;
  ap->accept(AST_editor) ;
  
  outputFile << "    void calculate(Loci::Entity _e_) " << endl ;
  AST_simplePrint printer(outputFile,-1,prettyOutput) ;
  ap->accept(printer) ;
  
  // AST_collectAccessInfo varaccess ;
  // ap->accept(varaccess) ;
  // //  cout << "variables = " << varaccess.accessed << endl ;
  // //  cout << "write variables = " << varaccess.writes << endl ;

  // variableSet readvars ;
  // variableSet writevars ;
      
  // for(auto i=varaccess.accessed.begin();i!=varaccess.accessed.end();++i) {
  //   readvars += i->var ;
  //   for(size_t j=0;j<i->mapping.size();++j)
  //     readvars += i->mapping[j] ;
  // }
  // for(auto i=varaccess.writes.begin();i!=varaccess.writes.end();++i) {
  //   writevars += i->var ;
  //   for(size_t j=0;j<i->mapping.size();++j)
  //     readvars += i->mapping[j] ;
  // }

  // cerr << "writevars=" << writevars << endl;
  // cerr << "readvars=" << readvars << endl ;
  // readvars -= writevars ;
  
  // // Now remove and save the open and close braces in the parseBlock
  // CPTR<AST_Block> bigblock = CPTR<AST_Block>(ap) ;
  // CPTR<AST_type> open = bigblock->elements[0] ;
  // int bsz = bigblock->elements.size() ;
  // CPTR<AST_type> close = bigblock->elements[bsz-1] ;
  // for(int i=0;i<bsz-1;++i)
  //   bigblock->elements[i] = bigblock->elements[i+1] ;
  // bigblock->elements.pop_back() ;
  // bigblock->elements.pop_back() ;
  
  // AST_simplePrint printer(outputFile,-1,prettyOutput) ;
  // map<string,string> maplist ;
  // for(auto i = varaccess.id2vmap.begin();i!=varaccess.id2vmap.end();++i) {
  //   auto p = vnames.find(*(i->second.var.begin())) ;
  //   string mapaccess = p->second ;
  //   string mapvar ;
  //   string mapsurrogate = "M_";
  //   mapaccess += "[" ;
  //   for(auto j = i->second.mapping.rbegin(); j!=i->second.mapping.rend();++j) {
  //     p = vnames.find(*(j->begin())) ;
  //     string mv = p->second ;
  //     mapvar += mv+"[" ;
  //     if(prettyOutput)
  //       mapsurrogate += mv ;
  //     else
  //       mapsurrogate += mv.substr(2,mv.size()-2) ;
  //   }
  //   mapvar += "_e_" ;
  //   for(auto j = i->second.mapping.rbegin(); j!=i->second.mapping.rend();++j) 
  //     mapvar +="]" ;
  //   maplist[mapsurrogate] = mapvar ;
  //   mapaccess += mapsurrogate + "]" ;

  //   printer.id2rename[i->first] = mapaccess ;
  // }
      
  // if(prettyOutput)
  //   outputFile << "    void calculate(Loci::Entity e) { " << endl ;
  // else
  //   outputFile << "    void calculate(Loci::Entity _e_) { " << endl ;

  // ap->accept(printer) ;
  // close->accept(printer) ;
}
string var2name(variable v) {
  string vn = v.str() ;
  string name ;
  if(!prettyOutput)
    name += "L_" ;
  for(size_t si=0;si!=vn.size();++si) {
    if(isalpha(vn[si]) || isdigit(vn[si]) || vn[si] == '_')
      name += vn[si] ;
    if(vn[si]=='{' || vn[si] == '}')
      name += '_' ;
    if(vn[si]=='=')
      name += "_EQ_" ;
    if(vn[si]=='+')
      name += "_P_" ;
    if(vn[si]=='-')
      name += "_M_" ;
  }
  if(!prettyOutput)
    name += "_" ;
  return name ;
}

// expand mapping list into all possible map strings
std::vector<list<variable> > expand_mapping(std::vector<variableSet> vset) {
  // if we have sliced off all of the variable sets, then the list is empty
  if(vset.size() == 0) {
    return std::vector<list<variable> >() ;
  }
  // get map set for the last item in the list
  variableSet mlast = vset.back() ;
  vset.pop_back() ;

  // expand remainder of list
  std::vector<list<variable> > tmp  = expand_mapping(vset) ;

  // Now build list by enumerating all maps from this level
  std::vector<list<variable> > tmp2 ;
  int tsz = tmp.size() ;
  if(tmp.size() == 0) {
    for(auto vi=mlast.begin();vi!=mlast.end();++vi) {
      list<variable> l1 ;
      l1.push_back(*vi) ;
      tmp2.push_back(l1) ;
    }
  } else {
    for(int i=0;i<tsz;++i) {
      for(auto vi=mlast.begin();vi!=mlast.end();++vi) {
        list<variable> l1= tmp[i] ;
        l1.push_back(*vi) ;
        tmp2.push_back(l1) ;
      }
    }
  }
  return tmp2 ;
}


void parseFile::setup_cudaRule(std::ostream &outputFile, const string &comment,
                               const parseSharedInfo &parseInfo) {
  killsp() ;
  string rule_type ;
  if(is_name(is)) {
    rule_type = get_name(is) ;
  } else 
    throw parseError("syntax error") ;
  if(rule_type != "pointwise")
    throw parseError("cudaRule only valid for pointwise") ;
  nestedparenstuff signature ;

  signature.get(is) ;
  line_no += signature.num_lines() ;
  nestedbracketstuff apply_op ;
  killsp() ;

  string constraint, conditional ;
  string parametric_var ;
  list<string> options ;
  list<string> comments ;
  list<pair<variable,variable> > inplace ;
  
  bool use_prelude = false ;
  bool is_specialized = false ;
  using namespace Loci ;
  while(is.peek() == ',') {
    is.get() ;
    killsp() ;
    if(!is_name(is))
      throw parseError("syntax error") ;

    string s = get_name(is) ;
    if(s == "constraint") {
      nestedparenstuff con ;
      con.get(is) ;
      if(constraint == "")
        constraint = con.str() ;
      else 
        constraint += "," + con.str() ;
      line_no += con.num_lines() ;
    } else if(s == "parametric") {
      nestedparenstuff con ;
      con.get(is) ;
      if(parametric_var != "") {
        throw parseError("syntax error: canot specify more than one parametric variable") ;
      }
        
      parametric_var = con.str() ;
      line_no += con.num_lines() ;
    } else if(s == "conditional") {
      nestedparenstuff con ;
      con.get(is) ;
      if(conditional != "") {
        throw parseError("syntax error: canot specify more than one conditional variable") ;
      }
      conditional = con.str() ;
      line_no += con.num_lines() ;
      // Check variable
      variable cond(conditional) ;
      auto mi  = lookupVarType(cond) ;

      if(!checkTypeValid(mi)) {
        cerr << filename << ':' << line_no << ":0: warning: type of conditional variable '" << cond << "' not found!"  << endl  ;
      } else {
        // clean up type string
        string val = mi->second.container + mi->second.container_args ;
        string val2 ;
        int valsz = val.size() ;
        for(int i=0;i<valsz;++i)
          if(val[i] != ' ' && val[i] != '\t' && val[i] != '\r' && val[i] != '\n')
            val2 += val[i] ;
        
        if(val2 != "param<bool>") {
          throw(parseError("conditional variable must be typed as a param<bool>")) ;
        }
      }
    } else if(s == "inplace") {
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      exprP p = expression::create(ip.str()) ;
      exprList l = collect_associative_op(p,OP_OR) ;
      if(l.size() != 2) 
        throw parseError("inplace needs two variables with a '|' separator") ;
        
      auto i = l.begin() ;
      variable v1(*i) ;
      ++i ;
      variable v2(*i) ;
      inplace.push_back(pair<variable,variable>(v1,v2)) ;
    
    } else if(s == "prelude") {
      use_prelude=true ;
      killsp() ;
      continue ;
    } else if(s == "specialized") {
      is_specialized = true ;
    } else if(s == "option") {
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      options.push_back(ip.str()) ;
    } else if(s == "comments") {
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      comments.push_back(ip.str()) ;
    } else {
      throw parseError("unknown rule modifier") ;
    }
    killsp() ;
  }

  string sig = signature.str() ;
  string heads,bodys ;
  exprP head=0,body=0 ;
  for(size_t i=0;i<sig.size()-1;++i) {
    if(sig[i]=='<' && sig[i+1]=='-') {
      heads = sig.substr(0,i) ;
      bodys = sig.substr(i+2,sig.size()) ;
      head = expression::create(heads) ;
      body = expression::create(bodys) ;
    }
  }
  if(head == 0) {
    heads = sig ;
    head = expression::create(heads) ;
    if(constraint == "") {
      throw parseError("rules without bodies should have a defined constraint as input!") ;
    }
  }
  
  string class_name = "file_" ;
  for(size_t i=0;i<filename.size();++i) {
    char c = filename[i] ;
    if(isalpha(c) || isdigit(c) || c=='_')
      class_name += c ;
    if(c == '.')
      break ;
  }
  class_name += '0' + (cnt/100)%10 ;
  class_name += '0' + (cnt/10)%10 ;
  class_name += '0' + (cnt)%10 ;

  //  timeb tdata ;
  //  ftime(&tdata) ;
  timespec tdata ;
  clock_gettime(CLOCK_MONOTONIC,&tdata) ;
  
  
  ostringstream tss ;
  tss << '_' << tdata.tv_sec << 'm' << tdata.tv_nsec/1000000 ;
  
  class_name += tss.str() ;
  cnt++ ;
  
  set<vmap_info> sources ;
  set<vmap_info> targets ;
  if(body != 0)
    fill_descriptors(sources,collect_associative_op(body,OP_COMMA)) ;
  fill_descriptors(targets,collect_associative_op(head,OP_COMMA)) ;

  variableSet input,output ;
  for(auto i=sources.begin();i!=sources.end();++i) {
    for(size_t j=0;j<i->mapping.size();++j)
      input += i->mapping[j] ;
    input += i->var ;
  }

  for(auto i=targets.begin();i!=targets.end();++i) {
    for(size_t j=0;j<i->mapping.size();++j)
      input += i->mapping[j] ;
    output += i->var ;
  }

  set<std::list<variable> > validate_set ;
  for(auto i=sources.begin();i!=sources.end();++i) {
    if(i->mapping.size() == 0) {
      for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
        std::list<variable> vbasic ;
      
        vbasic.push_back(*vi) ;
        validate_set.insert(vbasic) ;
      }
    } else {
      std::vector<std::list<variable> > maplist = expand_mapping(i->mapping) ;
      int msz = maplist.size() ;
      for(int j=0;j<msz;++j) {
        std::list<variable> mapping_list = maplist[j] ;
        validate_set.insert(mapping_list) ;
        for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
          std::list<variable> mapping_list2 = maplist[j] ;
          mapping_list2.push_back(*vi) ;
          validate_set.insert(mapping_list2) ;
        }
        mapping_list.pop_back() ;
        while(!mapping_list.empty()) {
          validate_set.insert(mapping_list) ;
          mapping_list.pop_back() ;
        }
      }
    }
  }

  for(auto i=targets.begin();i!=targets.end();++i) {
    if(i->mapping.size() == 0) {
      for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
        std::list<variable> vbasic ;
        variable vt = *vi ;
        while(vt.get_info().priority.size() != 0)
          vt = vt.drop_priority() ;
        vbasic.push_back(vt) ;
        validate_set.insert(vbasic) ;
      }
    } else {
      std::vector<std::list<variable> > maplist = expand_mapping(i->mapping) ;
      int msz = maplist.size() ;
      for(int j=0;j<msz;++j) {
        std::list<variable> mapping_list = maplist[j] ;
        validate_set.insert(mapping_list) ;
        for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
          std::list<variable> mapping_list2 = maplist[j] ;
          variable vt = *vi ;
          while(vt.get_info().priority.size() != 0)
            vt = vt.drop_priority() ;
          mapping_list2.push_back(vt) ;
          validate_set.insert(mapping_list2) ;
        }
        mapping_list.pop_back() ;
        while(!mapping_list.empty()) {
          validate_set.insert(mapping_list) ;
          mapping_list.pop_back() ;
        }
      }
    }
  }

  

  map<variable,string> vnames ;
  variableSet all_vars = input;
  all_vars += output ;

  for(auto vi=input.begin();vi!=input.end();++vi) {
    if(vi->get_info().priority.size() != 0) {
      ostringstream oss ;
      oss<< "improper use of priority annotation on rule input, var=" << *vi << endl ;
      throw parseError(oss.str()) ;
    }
  }

  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    vnames[*vi] = var2name(*vi) ;
    if(vi->get_info().priority.size() != 0) {
      variable v = *vi ;
      while(v.get_info().priority.size() != 0)
        v = v.drop_priority() ;
      vnames[v] = vnames[*vi] ;
    }
  }


  variableSet checkset ;
  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    variable v = *vi ;
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;
    checkset += v ;
  }
  for(auto ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    vnames[ipi->first] = vnames[ipi->second] ;
    variable v = ipi->first ;
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;
    if(!checkset.inSet(v)) {
      ostringstream oss ;
      oss << "inplace variable '"<< ipi->first << "' not input or output variable!" ;
      throw parseError(oss.str()) ;
    }
    v = ipi->second ;
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;
    if(!checkset.inSet(v)) {
      ostringstream oss ;
      oss << "inplace variable '"<< ipi->second << "' not input or output variable!" ;
      throw parseError(oss.str()) ;
    }
    
  }
  
  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      string s ;
      s = "unable to determine type of variable " ;
      s += (*vi).str() ;
      throw parseError(s) ;
    }
  }


  //  if(!prettyOutput)
  //    outputFile << "namespace {" ;
  outputFile << "class " << class_name << " : public Loci::" << rule_type << "_rule" ;
  if(rule_type == "pointwise") {
    for(auto vi=output.begin();vi!=output.end();++vi) {
      auto mi = lookupVarType(*vi) ;
      if(mi->second.container == "param" && vi->get_info().name != "OUTPUT") {
        throw(parseError("pointwise rule cannot compute param, use singleton")) ;
      }
    }
  }

  outputFile << " {" << endl ;
  syncFile(outputFile) ;

  variableSet outs = output ;
  for(auto ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    outs -= ipi->first ;
    outs += ipi->second ;
  }
  variableSet ins = input ;
  ins -= outs ;
  map<variable,string> typetable ;
  map<variable,string> ctypetable ;
  for(auto vi=ins.begin();vi!=ins.end();++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    ctypetable[*vi] = mi->second.container ;
    if(mi->second.container == "Map") {
      typetable[*vi] = "int" ;
    } else if(mi->second.container == "MapVec") {
      string scratch = mi->second.container_args ;
      if(scratch.size() > 2) {
        typetable[*vi] = string("Array<Entity,") + scratch.substr(1,scratch.size()-3) + ">" ;
      } else {
        cerr << "unexpected loci variable type!" << endl ;
      }
    } else {
      string scratch = mi->second.container_args ;
      if(scratch.size() > 2) {
        typetable[*vi] = scratch.substr(1,scratch.size()-3) ;
      } else {
        cerr << "unexpected loci variable type!" << endl ;
      }
    }

    if(!prettyOutput) {
      outputFile << "    Loci::const_gpu" << mi->second.container
                 <<  mi->second.container_args ;
    } else {
      outputFile << "    const_gpu" << mi->second.container
                 <<  mi->second.container_args ;
    }

    outputFile << " " << vnames[*vi] << " ; " << endl ;
    syncFile(outputFile) ;
  }

  for(auto vi=outs.begin();vi!=outs.end();++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    ctypetable[*vi] = mi->second.container ;
    if(mi->second.container == "Map") {
      typetable[*vi] = "int" ;
    } else if(mi->second.container == "MapVec") {
      string scratch = mi->second.container_args ;
      if(scratch.size() > 2) {
        typetable[*vi] = string("Array<Entity,") + scratch.substr(1,scratch.size()-3) + ">" ;
      } else {
        cerr << "unexpected loci variable type!" << endl ;
      }
    } else {
      string scratch = mi->second.container_args ;
      if(scratch.size() > 2) {
        typetable[*vi] = scratch.substr(1,scratch.size()-3) ;
      } else {
        cerr << "unexpected loci variable type!" << endl ;
      }
    }
    if(!prettyOutput)
      outputFile << "    Loci::gpu" << mi->second.container
		 <<  mi->second.container_args ;
    else
      outputFile << "    gpu" << mi->second.container <<
	mi->second.container_args ;
    outputFile << " " << vnames[*vi] << " ; " << endl ;
    syncFile(outputFile) ;
  }
  outputFile << "public:" << endl ;
  syncFile(outputFile) ;
  outputFile <<   "    " << class_name << "() {" << endl ;
  syncFile(outputFile) ;
  for(auto ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    all_vars -= ipi->first ;
  }

  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    outputFile << "       name_store(\"" << *vi << "\","
               << vnames[*vi] << ") ;" << endl ;
    syncFile(outputFile) ;
  }
  if(bodys != "") {
    outputFile <<   "       input(\"" << bodys << "\") ;" << endl ;
    syncFile(outputFile) ;
  }

  for(auto i=targets.begin();i!=targets.end();++i) {
    outputFile <<   "       output(\"" ;
    for(size_t j=0;j<i->mapping.size();++j)
      outputFile << i->mapping[j] << "->" ;

    // Output target variables, adding inplace notation if needed
    if(i->var.size() > 1)
      outputFile << '(' ;
    for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
      if(vi != i->var.begin())
        outputFile << ',' ;
      auto ipi = inplace.begin() ;
      for(;ipi!=inplace.end();++ipi) {
        if((ipi->first) == *vi)
          break ;
      }
      if(ipi!=inplace.end()) {
        if(i->mapping.size() == 0 || i->var.size() > 1)
          outputFile << ipi->first << "=" << ipi->second ;
        else
          outputFile << '('<<ipi->first << "=" << ipi->second <<')';
      } else
        outputFile << *vi ;
    }
    if(i->var.size() > 1)
      outputFile << ')' ;

    outputFile <<  "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  outputFile << "disable_threading() ;" << endl ;
  syncFile(outputFile) ;
  //  outputFile <<   "       output(\"" << heads << "\") ;" << endl ;
  //  syncFile(outputFile) ;

  if(constraint!="") {
    // Check to see that the constraint is typed
    exprP C = expression::create(constraint) ;
    set<vmap_info> Cdigest ;
    fill_descriptors(Cdigest,collect_associative_op(C,OP_COMMA)) ;
    variableSet constraint_vars ;
    for(auto i=Cdigest.begin();i!=Cdigest.end();++i) {
      for(size_t j=0;j<i->mapping.size();++j)
	constraint_vars += i->mapping[j] ;
      constraint_vars += i->var ;
    }

    for(auto vi=constraint_vars.begin();vi!=constraint_vars.end();++vi) {
      auto mi = lookupVarType(*vi) ;

      if(!checkTypeValid(mi)) {

        cerr << filename << ':' << line_no << ":0: warning: type of constraint variable '" << *vi << "' not found!"  << endl  ;

      }
    }

    outputFile <<   "       constraint(\"" << constraint << "\") ;" << endl ;
    syncFile(outputFile) ;
  }

  if(parametric_var != "") {
    outputFile <<   "       set_parametric_variable(\""
               << parametric_var << "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  if(is_specialized) {
    outputFile <<   "       set_specialized() ; " << endl ;
    syncFile(outputFile) ;
  }
  if(conditional!="") {
    outputFile <<   "       conditional(\"" << conditional << "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  for(auto lsi=options.begin();lsi!=options.end();++lsi) {
    string s = *lsi ;
    bool has_paren = false ;
    for(size_t i = 0;i<s.size();++i)
      if(s[i] == '(')
        has_paren = true ;
    outputFile <<   "       " << s ;
    if(!has_paren)
      outputFile << "()" ;
    outputFile << " ;" << endl;
    syncFile(outputFile) ;
  }
  if(comments.size() == 0) {
    // check to see if there is a javadoc compatible comment before
    if(comment.size() > 0)
      comments.push_back(comment) ;
    else if(rule_type=="optional" || rule_type=="default") {
      auto mi = lookupVarType(*output.begin()) ;
      comments.push_back(mi->second.comment) ;
    }
  }
  for(auto lsi=comments.begin();lsi!=comments.end();++lsi) {
    outputFile <<   "       comments(" << *lsi << ") ;" << endl ;
    syncFile(outputFile) ;
  }
  // document file location of rule
  outputFile << "       set_file(\"" << filename << ":" << line_no << "\") ;" << endl ;
  syncFile(outputFile) ;
  outputFile <<   "    }" << endl ;
  syncFile(outputFile) ;

  //  bool use_compute = true ;

  if(use_prelude) {
    throw parseError("prelude not compatible with cuda rules") ;
  }

  //  if(use_compute && is.peek() != '{')
  //    throw parseError("syntax error, expecting '{'") ;

  bool sized_outputs = false;
  variableSet outsmi = outs ;
  outsmi -= input ;
  for(auto vi=outsmi.begin();vi!=outsmi.end();++vi) {
    auto mi = lookupVarType(*vi) ;
    const string &ot  = mi->second.container ;
    if(ot == "storeVec" || ot == "storeMat" || ot == "multiStore")
      sized_outputs = true ;
  }
  if(sized_outputs)
    throw parseError("cuda rules currently incompatible with storeVec, storemat or multiStore types") ;

  //  if(use_compute)
  //    process_Calculate(outputFile,vnames,validate_set) ;

  outputFile <<   "    void compute(const Loci::sequence &seq) ;" << endl ;
  syncFile(outputFile) ;
  outputFile <<   "} ;" << endl ;
  syncFile(outputFile) ;

  //  process_Calculate(outputFile,vnames,validate_set) ;
  varmap typemap ;

  if(is.peek() != '{')
    throw parseError("syntax error, expecting '{'") ;

  int startline = line_no ;

  CPTR<AST_type> ap = parseBlock(is,line_no,filename,typemap) ;
  //    outputFile << "Parsed TEST:" << endl ;

  AST_condenseLeftAssociative condenseOps ;
  ap->accept(condenseOps) ;

  AST_editLociMapArrayAccess mapEditOps ;
  ap->accept(mapEditOps) ;

  if(parseInfo.diag_level > 0) {
    AST_printTree diagout(cerr) ;
    ap->accept(diagout) ;
  }

  AST_errorCheck syntaxChecker ;
  ap->accept(syntaxChecker) ;
  if(syntaxChecker.hasErrors()) {
#ifdef VERBOSE
    AST_simplePrint printer(cerr,-1,false) ;
    ap->accept(printer) ;
#endif
    throw parseError("syntax error") ;
  }
  
  AST_collectAccessInfo varaccess ;
  ap->accept(varaccess) ;
  //cerr << "variables = " << varaccess.accessed << endl ;
  //cerr << "write variables = " << varaccess.writes << endl ;
  //for(auto i = varaccess.id2var.begin();i!=varaccess.id2var.end();++i) {
  //  cerr << "id2var[" << i->first << "] = " << i->second << endl ;
  //}
  //for(auto i = varaccess.id2vmap.begin();i!=varaccess.id2vmap.end();++i) {
  //  cerr << "id2vmap[" << i->first << "] = " << i->second << endl ;
  //}

  variableSet readvars ;
  variableSet writevars ;
  
  for(auto i=varaccess.accessed.begin();i!=varaccess.accessed.end();++i) {
    readvars += i->var ;
    for(size_t j=0;j<i->mapping.size();++j)
      readvars += i->mapping[j] ;
  }
  for(auto i=varaccess.writes.begin();i!=varaccess.writes.end();++i) {
    writevars += i->var ;
    for(size_t j=0;j<i->mapping.size();++j)
      readvars += i->mapping[j] ;
  }

  readvars -= writevars ;
  
  // Now remove and save the open and close braces in the parseBlock
  CPTR<AST_Block> bigblock = CPTR<AST_Block>(ap) ;
  CPTR<AST_type> open = bigblock->elements[0] ;
  int bsz = bigblock->elements.size() ;
  CPTR<AST_type> close = bigblock->elements[bsz-1] ;
  for(int i=0;i<bsz-1;++i)
    bigblock->elements[i] = bigblock->elements[i+1] ;
  bigblock->elements.pop_back() ;
  bigblock->elements.pop_back() ;
  
  AST_simplePrint printer(outputFile,-1,prettyOutput) ;

  if(!prettyOutput)
    outputFile << "#line " << startline << endl ;
  outputFile << "__global__" << endl ;
  if(!prettyOutput)
    outputFile << "#line " << startline << endl ;
  outputFile << "void "<< class_name << "_kernel(" ;
  for(auto i=writevars.begin();i!=writevars.end();) {
    outputFile << typetable[*i]<< " *" << vnames[*i] ;
    ++i ;
    if(i!=writevars.end())
      outputFile << "," ;
  }
  if(readvars!=EMPTY)
    outputFile << "," ;
  for(auto i=readvars.begin();i!=readvars.end();) {
    outputFile << "const " << typetable[*i] << " *" << vnames[*i] ;
    ++i ;
    if(i!=readvars.end())
      outputFile << "," ;
  }
  outputFile << ", int _start_, int _end_, int _blksz_)" << endl ;

  open->accept(printer) ;

  if(!prettyOutput)
    outputFile << endl << "#line " << printer.lineno  << endl ;
  outputFile << "   int _e_ = _start_+blockIdx.x*_blksz_+threadIdx.x ;" <<endl;
  if(!prettyOutput)
    outputFile << "#line " << printer.lineno  << endl ;
  outputFile << "   if(_e_ <= _end_) {" << endl ;
  if(!prettyOutput)
    outputFile <<  "#line " << printer.lineno << endl  ;

  AST_editLociVariableAccess2 AST_editor(vnames, ctypetable) ;
  ap->accept(AST_editor) ;

  ap->accept(printer) ;

  close->accept(printer) ;
  close->accept(printer) ;
  outputFile << endl ;
  syncFile(outputFile) ;
  outputFile << "void " << class_name << "::compute(const Loci::sequence &seq) {" << endl ;
  syncFile(outputFile) ;
  outputFile << "  const int NTHREADS=64 ;" << endl ;
  syncFile(outputFile) ;
  outputFile << "  const int ni = seq.num_intervals() ;"<< endl ;
  syncFile(outputFile) ;
  outputFile << "  for(int i=0;i<ni;++i) {" << endl ;
  syncFile(outputFile) ;
  outputFile << "    const Loci::int_type i1 = seq[i].first ;" << endl ;
  syncFile(outputFile) ;
  outputFile << "    const Loci::int_type i2 = seq[i].second  ;" << endl ; 
  syncFile(outputFile) ;
  outputFile << "    const Loci::int_type start = min(i1,i2) ;" << endl ;
  syncFile(outputFile) ;
  outputFile << "    const Loci::int_type stop = max(i1,i2) ;" << endl ;
  syncFile(outputFile) ;
  outputFile << "    const int nblks = (stop-start+NTHREADS)/NTHREADS ;" << endl ;
  syncFile(outputFile) ;

  if(parseInfo.debug_info > 0) {
    ostringstream oss ;
    for(auto i=targets.begin();i!=targets.end();) {
      for(size_t j=0;j<i->mapping.size();++j)
        oss << i->mapping[j] << "->" ;
      // Output target variables, adding inplace notation if needed
      if(i->var.size() > 1)
        oss << '(' ;
      for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
        if(vi != i->var.begin())
          oss << ',' ;
        oss << *vi ;
      }
      if(i->var.size() > 1)
        oss << ')' ;
      ++i;
      if(i != targets.end())
        oss << "," ;
    }
    oss << "<-" ;
    oss << bodys ;
    if(constraint!="")
      oss << ",constraint(" << constraint<<")" ;
    outputFile << "nvtxRangePush(\"" << oss.str() << "\");" << endl ;
    syncFile(outputFile) ;
  }
  outputFile <<"    " <<class_name << "_kernel<<<nblks,NTHREADS,0,Loci::getGPUStream()>>>(";
  for(auto i=writevars.begin();i!=writevars.end();) {
    outputFile << vnames[*i] << ".ptr()";
    ++i ;
    if(i!=writevars.end())
      outputFile << "," ;
  }
  if(readvars!=EMPTY)
    outputFile << "," ;
  for(auto i=readvars.begin();i!=readvars.end();) {
    outputFile <<vnames[*i]<< ".ptr()" ;
    ++i ;
    if(i!=readvars.end())
      outputFile << "," ;
  }
  outputFile << ",start, stop, NTHREADS) ;"  ;
  if(parseInfo.debug_info > 0) {
    outputFile << "nvtxRangePop();" ;
  }
  outputFile<<endl ;
  syncFile(outputFile) ;
  syncFile(outputFile) ;
  outputFile << "}" << endl; // end of for loop
  syncFile(outputFile) ;
  outputFile << "}" << endl; // end of compute method
  syncFile(outputFile) ;
  
  if(!prettyOutput)
    outputFile << "Loci::register_rule<"<<class_name<<"> register_"<<class_name
               << " ;" << endl ;
  else
    outputFile << "register_rule<"<<class_name<<"> register_"<<class_name
               << " ;" << endl ;
  syncFile(outputFile) ;

  //  if(!prettyOutput) {
  //    outputFile << "}" << endl ;
  //    syncFile(outputFile) ;
  //  }

  
  if(!use_prelude && sized_outputs && (rule_type != "apply")) 
    throw parseError("need prelude to size output type!") ;
}
// rule_type
void parseFile::setup_Rule(std::ostream &outputFile, const string &comment,
                           const parseSharedInfo &parseInfo) {
  killsp() ;
  string rule_type ;
  if(is_name(is)) {
    rule_type = get_name(is) ;
  } else 
    throw parseError("syntax error") ;
  nestedparenstuff signature ;

  signature.get(is) ;
  line_no += signature.num_lines() ;
  nestedbracketstuff apply_op ;
  killsp() ;
  if(rule_type == "apply") {
    if(is.peek() != '[') 
      throw parseError("apply rule missing '[operator]'") ;
    apply_op.get(is) ;
    line_no += apply_op.num_lines() ;
    killsp() ;
  }
  

  string constraint, conditional ;
  string parametric_var ;
  list<string> options ;
  list<string> comments ;
  list<pair<variable,variable> > inplace ;
  
  bool use_prelude = false ;
  bool is_specialized = false ;
  while(is.peek() == ',') {
    is.get() ;
    killsp() ;
    if(!is_name(is))
      throw parseError("syntax error") ;

    string s = get_name(is) ;
    if(s == "constraint") {
      nestedparenstuff con ;
      con.get(is) ;
      if(constraint == "")
        constraint = con.str() ;
      else 
        constraint += "," + con.str() ;
      line_no += con.num_lines() ;
    } else if(s == "parametric") {
      nestedparenstuff con ;
      con.get(is) ;
      if(parametric_var != "") {
        throw parseError("syntax error: canot specify more than one parametric variable") ;
      }
        
      parametric_var = con.str() ;
      line_no += con.num_lines() ;
    } else if(s == "conditional") {
      nestedparenstuff con ;
      con.get(is) ;
      if(conditional != "") {
        throw parseError("syntax error: canot specify more than one conditional variable") ;
      }
      conditional = con.str() ;
      line_no += con.num_lines() ;
      // Check variable
      variable cond(conditional) ;
      
      auto mi = lookupVarType(cond) ;
      if(!checkTypeValid(mi)) {
        cerr << filename << ':' << line_no << ":0: warning: type of conditional variable '" << cond << "' not found!"  << endl  ;
      } else {
        // clean up type string
        string val = mi->second.container + mi->second.container_args ;
        string val2 ;
        int valsz = val.size() ;
        for(int i=0;i<valsz;++i)
          if(val[i] != ' ' && val[i] != '\t' && val[i] != '\r' && val[i] != '\n')
            val2 += val[i] ;
        
        if(val2 != "param<bool>") {
          throw(parseError("conditional variable must be typed as a param<bool>")) ;
        }
      }
    } else if(s == "inplace") {
      using namespace Loci ;
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      exprP p = expression::create(ip.str()) ;
      exprList l = collect_associative_op(p,OP_OR) ;
      if(l.size() != 2) 
        throw parseError("inplace needs two variables with a '|' separator") ;
        
      auto i = l.begin() ;
      variable v1(*i) ;
      ++i ;
      variable v2(*i) ;
      inplace.push_back(pair<variable,variable>(v1,v2)) ;
    
    } else if(s == "prelude") {
      use_prelude=true ;
      killsp() ;
      continue ;
    } else if(s == "specialized") {
      is_specialized = true ;
    } else if(s == "option") {
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      options.push_back(ip.str()) ;
    } else if(s == "comments") {
      nestedparenstuff ip ;
      ip.get(is) ;
      line_no += ip.num_lines() ;
      comments.push_back(ip.str()) ;
    } else {
      throw parseError("unknown rule modifier") ;
    }
    killsp() ;
  }

  string sig = signature.str() ;
  string heads,bodys ;
  exprP head=0,body=0 ;
  for(size_t i=0;i<sig.size()-1;++i) {
    if(sig[i]=='<' && sig[i+1]=='-') {
      heads = sig.substr(0,i) ;
      bodys = sig.substr(i+2,sig.size()) ;
      head = expression::create(heads) ;
      body = expression::create(bodys) ;
      if(rule_type == "optional" || rule_type == "default") {
	throw parseError("'optional' or 'default' rules should not have a body (defined by '<-' operator)!") ;
      }
    }
  }
  if(head == 0) {
    heads = sig ;
    head = expression::create(heads) ;
    if(rule_type == "optional" || rule_type == "default") {
      if(constraint != "") 
	throw parseError("'optional' or 'default' rules should not have a constraint!") ;      
    } else {
      if(constraint == "") {
	throw parseError("rules without bodies should have a defined constraint as input!") ;
      }
    }
  }
  
  string class_name = "file_" ;
  for(size_t i=0;i<filename.size();++i) {
    char c = filename[i] ;
    if(isalpha(c) || isdigit(c) || c=='_')
      class_name += c ;
    if(c == '.')
      break ;
  }
  class_name += '0' + (cnt/100)%10 ;
  class_name += '0' + (cnt/10)%10 ;
  class_name += '0' + (cnt)%10 ;

  //  timeb tdata ;
  //  ftime(&tdata) ;
  timespec tdata ;
  clock_gettime(CLOCK_MONOTONIC,&tdata) ;
  
  
  ostringstream tss ;
  //  tss <<  '_' << tdata.time << 'm'<< tdata.millitm;
  tss << '_' << tdata.tv_sec << 'm' << tdata.tv_nsec/1000000 ;
  
  class_name += tss.str() ;
  cnt++ ;
#ifdef OLD
  class_name += "_rule_" ;
  if(conditional != "")
    sig += "_" + conditional ;
  if(constraint != "")
    sig += "_" + constraint ;
  for(size_t i=0;i<sig.size();++i) {
    if(isalpha(sig[i]) || isdigit(sig[i]))
      class_nam e+= sig[i] ;
    if(sig[i] == ',' || sig[i] == '-' || sig[i] == '>' || sig[i] == '('||
       sig[i] == ')' || sig[i] == '{' || sig[i] == '}' || sig[i] == '='||
       sig[i] == '+' || sig[i] == '_')
      class_name += '_' ;
  }
#endif
  using namespace Loci ;  
  set<vmap_info> sources ;
  set<vmap_info> targets ;
  if(body != 0)
    fill_descriptors(sources,collect_associative_op(body,OP_COMMA)) ;
  fill_descriptors(targets,collect_associative_op(head,OP_COMMA)) ;

  variableSet input,output ;
  for(auto i=sources.begin();i!=sources.end();++i) {
    for(size_t j=0;j<i->mapping.size();++j)
      input += i->mapping[j] ;
    input += i->var ;
  }

  for(auto i=targets.begin();i!=targets.end();++i) {
    for(size_t j=0;j<i->mapping.size();++j)
      input += i->mapping[j] ;
    output += i->var ;
  }

  set<std::list<variable> > validate_set ;
  for(auto i=sources.begin();i!=sources.end();++i) {
    if(i->mapping.size() == 0) {
      for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
        std::list<variable> vbasic ;
      
        vbasic.push_back(*vi) ;
        validate_set.insert(vbasic) ;
      }
    } else {
      std::vector<std::list<variable> > maplist = expand_mapping(i->mapping) ;
      int msz = maplist.size() ;
      for(int j=0;j<msz;++j) {
        std::list<variable> mapping_list = maplist[j] ;
        validate_set.insert(mapping_list) ;
        for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
          std::list<variable> mapping_list2 = maplist[j] ;
          mapping_list2.push_back(*vi) ;
          validate_set.insert(mapping_list2) ;
        }
        mapping_list.pop_back() ;
        while(!mapping_list.empty()) {
          validate_set.insert(mapping_list) ;
          mapping_list.pop_back() ;
        }
      }
    }
  }

  for(auto i=targets.begin();i!=targets.end();++i) {
    if(i->mapping.size() == 0) {
      for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
        std::list<variable> vbasic ;
        variable vt = *vi ;
        while(vt.get_info().priority.size() != 0)
          vt = vt.drop_priority() ;
        vbasic.push_back(vt) ;
        validate_set.insert(vbasic) ;
      }
    } else {
      std::vector<std::list<variable> > maplist = expand_mapping(i->mapping) ;
      int msz = maplist.size() ;
      for(int j=0;j<msz;++j) {
        std::list<variable> mapping_list = maplist[j] ;
        validate_set.insert(mapping_list) ;
        for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
          std::list<variable> mapping_list2 = maplist[j] ;
          variable vt = *vi ;
          while(vt.get_info().priority.size() != 0)
            vt = vt.drop_priority() ;
          mapping_list2.push_back(vt) ;
          validate_set.insert(mapping_list2) ;
        }
        mapping_list.pop_back() ;
        while(!mapping_list.empty()) {
          validate_set.insert(mapping_list) ;
          mapping_list.pop_back() ;
        }
      }
    }
  }

  

  map<variable,string> vnames ;
  variableSet all_vars = input;
  all_vars += output ;

  for(auto vi=input.begin();vi!=input.end();++vi) {
    if(vi->get_info().priority.size() != 0) {
      ostringstream oss ;
      oss<< "improper use of priority annotation on rule input, var=" << *vi << endl ;
      throw parseError(oss.str()) ;
    }
  }

  if(rule_type != "pointwise" && rule_type != "default") {
    for(auto vi=output.begin();vi!=output.end();++vi) {
      if(vi->get_info().priority.size() != 0) {
        ostringstream oss ;
        oss << "only pointwise rules can use priority annotation, var="<< *vi << endl ;
        throw parseError(oss.str()) ;
      }
    }
  }    
  
  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    vnames[*vi] = var2name(*vi) ;
    if(vi->get_info().priority.size() != 0) {
      variable v = *vi ;
      while(v.get_info().priority.size() != 0)
        v = v.drop_priority() ;
      vnames[v] = vnames[*vi] ;
    }
  }


  variableSet checkset ;
  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    variable v = *vi ;
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;
    checkset += v ;
  }
  for(auto ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    vnames[ipi->first] = vnames[ipi->second] ;
    variable v = ipi->first ;
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;
    if(!checkset.inSet(v)) {
      ostringstream oss ;
      oss << "inplace variable '"<< ipi->first << "' not input or output variable!" ;
      throw parseError(oss.str()) ;
    }
    v = ipi->second ;
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;
    if(!checkset.inSet(v)) {
      ostringstream oss ;
      oss << "inplace variable '"<< ipi->second << "' not input or output variable!" ;
      throw parseError(oss.str()) ;
    }
    
  }
  
  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      string s ;
      s = "unable to determine type of variable " ;
      s += (*vi).str() ;
      throw parseError(s) ;
    }
  }

  if(!prettyOutput)
    outputFile << "namespace {" ;
  outputFile << "class " << class_name << " : public Loci::" << rule_type << "_rule" ;
  if(rule_type == "pointwise") {
    for(auto vi=output.begin();vi!=output.end();++vi) {
      auto mi = lookupVarType(*vi) ;
      if(mi->second.container == "param" && vi->get_info().name != "OUTPUT") {
        throw(parseError("pointwise rule cannot compute param, use singleton")) ;
      }
    }
  }
  if(rule_type == "singleton") {
    for(auto vi=output.begin();vi!=output.end();++vi) {
      auto mi = lookupVarType(*vi) ;
      const string &t = mi->second.container ;
      if(t == "store" || t == "storeVec" || t == "multiStore") {
        throw(parseError("singleton rule cannot compute store's, use pointwise")) ;
      }
    }
  }
    
              
  bool singletonApply = false ;
  if(rule_type == "apply") {
    if(output.size() != 1) 
      throw parseError("apply rule should have only one output variable") ;
    variable av = *(output.begin()) ;
    typedoc tinfo = lookupVarType(av)->second ;
    outputFile << "< " << tinfo.container << tinfo.container_args <<","
               << apply_op.str() ;
    if(tinfo.container == "storeVec") {
      outputFile << "<Vect" << tinfo.container_args <<" > " ;
    } else if(tinfo.container == "storeMat") {
      outputFile << "<Mat" << tinfo.container_args <<" > " ;
    } else {
      outputFile << tinfo.container_args ;
    }
    if(tinfo.container == "param") {
      bool allparam = true ;
      for(auto vi=input.begin();vi!=input.end();++vi) {
        typedoc tinfo2 = lookupVarType(*vi)->second ;
        if(tinfo2.container != "param") {
          allparam = false ;
        }
      }
      if(allparam)
        singletonApply = true ;
    }
    outputFile << "> " ;
  }
  outputFile << " {" << endl ;
  syncFile(outputFile) ;

  variableSet outs = output ;
  for(auto ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    outs -= ipi->first ;
    outs += ipi->second ;
  }
  variableSet ins = input ;
  ins -= outs ;
  for(auto vi=ins.begin();vi!=ins.end();++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    if(!prettyOutput) 
      outputFile << "    Loci::const_" << mi->second.container
		 <<  mi->second.container_args ;
    else 
      outputFile << "    const_" << mi->second.container
		 <<  mi->second.container_args ;
    outputFile << " " << vnames[*vi] << " ; " << endl ;
    syncFile(outputFile) ;
  }
  bool output_param = false ;
  for(auto vi=outs.begin();vi!=outs.end();++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }
    if(vi->get_info().name != "OUTPUT" && mi->second.container == "param") {
      output_param= true ;
    }
    if(!prettyOutput)
      outputFile << "    Loci::" << mi->second.container
		 <<  mi->second.container_args ;
    else
      outputFile << "    " << mi->second.container
		 <<  mi->second.container_args ;
    outputFile << " " << vnames[*vi] << " ; " << endl ;
    syncFile(outputFile) ;
  }
  outputFile << "public:" << endl ;
  syncFile(outputFile) ;
  outputFile <<   "    " << class_name << "() {" << endl ;
  syncFile(outputFile) ;
  for(auto ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    all_vars -= ipi->first ;
  }

  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    outputFile << "       name_store(\"" << *vi << "\","
               << vnames[*vi] << ") ;" << endl ;
    syncFile(outputFile) ;
    auto mi = access_map.find(lookupVarType(*vi)->second.getFileLoc()) ;
    if(mi != access_map.end()) {
      outputFile << "       store_info_id(\"" << *vi << "\","
               << mi->second << ") ;" << endl ;
      syncFile(outputFile) ;
    }
  }
  if(bodys != "") {
    outputFile <<   "       input(\"" << bodys << "\") ;" << endl ;
    syncFile(outputFile) ;
  }

  for(auto i=targets.begin();i!=targets.end();++i) {
    outputFile <<   "       output(\"" ;
    for(size_t j=0;j<i->mapping.size();++j)
      outputFile << i->mapping[j] << "->" ;

    // Output target variables, adding inplace notation if needed
    if(i->var.size() > 1)
      outputFile << '(' ;
    for(auto vi=i->var.begin();vi!=i->var.end();++vi) {
      if(vi != i->var.begin())
        outputFile << ',' ;
      auto ipi = inplace.begin() ;
      for(;ipi!=inplace.end();++ipi) {
        if((ipi->first) == *vi)
          break ;
      }
      if(ipi!=inplace.end()) {
        if(i->mapping.size() == 0 || i->var.size() > 1)
          outputFile << ipi->first << "=" << ipi->second ;
        else
          outputFile << '('<<ipi->first << "=" << ipi->second <<')';
      } else
        outputFile << *vi ;
    }
    if(i->var.size() > 1)
      outputFile << ')' ;

    outputFile <<  "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  //  outputFile <<   "       output(\"" << heads << "\") ;" << endl ;
  //  syncFile(outputFile) ;

  if(constraint!="") {
    // Check to see that the constraint is typed
    exprP C = expression::create(constraint) ;
    set<vmap_info> Cdigest ;
    fill_descriptors(Cdigest,collect_associative_op(C,OP_COMMA)) ;
    variableSet constraint_vars ;
    for(auto i=Cdigest.begin();i!=Cdigest.end();++i) {
      for(size_t j=0;j<i->mapping.size();++j)
	constraint_vars += i->mapping[j] ;
      constraint_vars += i->var ;
    }

    variableSet undoc = constraint_vars ;
    undoc -= all_vars ;
    for(auto vi=undoc.begin();vi!=undoc.end();++vi) {
      if(checkTypeValid(lookupVarType(*vi))) {
        auto mi = access_map.find(lookupVarType(*vi)->second.getFileLoc()) ;
        if(mi != access_map.end()) {
          outputFile << "       store_info_id(\"" << *vi << "\","
                     << mi->second << ") ;" << endl ;
          syncFile(outputFile) ;
        }
      }
    }
    for(auto vi=constraint_vars.begin();vi!=constraint_vars.end();++vi) {
      auto mi = lookupVarType(*vi) ;
      
      if(!checkTypeValid(mi)) {
        cerr << filename << ':' << line_no << ":0: warning: type of constraint variable '" << *vi << "' not found!"  << endl  ;

      } 
    }
    
    

    outputFile <<   "       constraint(\"" << constraint << "\") ;" << endl ;
    syncFile(outputFile) ;
  }

  if(parametric_var != "") {
    outputFile <<   "       set_parametric_variable(\""
               << parametric_var << "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  if(is_specialized) {
    outputFile <<   "       set_specialized() ; " << endl ;
    syncFile(outputFile) ;
  }
  if(conditional!="") {
    outputFile <<   "       conditional(\"" << conditional << "\") ;" << endl ;
    syncFile(outputFile) ;
  }
  for(auto lsi=options.begin();lsi!=options.end();++lsi) {
    string s = *lsi ;
    bool has_paren = false ;
    for(size_t i = 0;i<s.size();++i)
      if(s[i] == '(')
        has_paren = true ;
    outputFile <<   "       " << s ;
    if(!has_paren)
      outputFile << "()" ;
    outputFile << " ;" << endl;
    syncFile(outputFile) ;
  }
  
  if(comments.size() == 0) {
    // check to see if there is a javadoc compatible comment before
    if(comment.size() > 0)
      comments.push_back("\""+cleanupComment(comment)+"\"") ;
    else if(rule_type=="optional" || rule_type=="default") {
      auto mi = lookupVarType(*output.begin()) ;
      comments.push_back("\""+cleanupComment(mi->second.comment)+"\"") ;
    }
  }
  for(auto lsi=comments.begin();lsi!=comments.end();++lsi) {
    outputFile <<   "       comments(" << *lsi << ") ;" << endl ;
    syncFile(outputFile) ;
  }
  // document file location of rule
  outputFile << "       setvardoc(" << docvarname << ") ;" << endl ;
  
  outputFile << "       set_file(\"" << filename << ":" << line_no << "\") ;" << endl ;

  syncFile(outputFile) ;
  outputFile <<   "    }" << endl ;
  syncFile(outputFile) ;

  bool use_compute = true ;

  if(use_prelude) {
    process_Prelude(outputFile,vnames) ;
    killsp() ;
    if(is.peek() == ';') {
      is.get() ;
      use_compute = false ;
    }
    if(is_name(is)) {
      string s = get_name(is) ;
      if(s != "compute") {
        throw parseError("syntax error, expecting 'compute'") ;
      }
    }
    killsp() ;
  }

  
  if(use_compute && is.peek() != '{')
    throw parseError("syntax error, expecting '{'") ;

  bool sized_outputs = false;
  variableSet outsmi = outs ;
  outsmi -= input ;
  for(auto vi=outsmi.begin();vi!=outsmi.end();++vi) {
    const string &ot = lookupVarType(*vi)->second.container ;
    if(ot == "storeVec" || ot == "storeMat" || ot == "multiStore")
      sized_outputs = true ;
  }

    
    

  if(rule_type == "singleton" ||
     rule_type == "optional"  ||
     rule_type == "default" ||
     rule_type == "constraint" ||
     (output_param && rule_type != "apply" ) ) {
    if(use_prelude) {
      string error = "inappropriate prelude on " + rule_type + " rule." ;
      throw parseError(error) ;
    }
    process_Compute(outputFile,vnames) ;
  } else {
    if(use_compute) {
      if(parseInfo.test_parse) 
        process_Calculate2(outputFile,vnames,validate_set,parseInfo) ;
      else
        process_Calculate(outputFile,vnames,validate_set) ;
    }
    outputFile <<   "    void compute(const Loci::sequence &seq) { " << endl ;
    syncFile(outputFile) ;
//     if(use_prelude) {
//       outputFile <<   "      prelude(seq) ;" << endl ;
//       syncFile(outputFile) ;
//     }
    if(use_compute) {
      if(singletonApply) {
        cerr << "NOTE: parameter only apply rule on '" << output << "' now executes single instance." << endl ;
        // Note, this is better than before, but if rank 0 owns no entity
        // we still may get an out of bounds error with bounds checking turned
        // on.  Technically this isn't wrong except for the edge case that
        // the rule is applied over an empty set.  This probably will require
        // some work on the scheduling to fix, but doesn't impact any current
        // use cases.
        outputFile <<   "      if(Loci::MPI_rank == 0) calculate(seq.num_intervals()>0?seq[0].first:0)  ; " << endl ;
        syncFile(outputFile) ;
      } else {
        outputFile <<   "      do_loop(seq,this) ;" << endl ;
        syncFile(outputFile) ;
      }
    }
    outputFile <<   "    }" << endl ;
    syncFile(outputFile) ;
  }
  outputFile <<   "} ;" << endl ;
  syncFile(outputFile) ;

  if(!prettyOutput)
    outputFile << "Loci::register_rule<"<<class_name<<"> register_"<<class_name
               << " ;" << endl ;
  else
    outputFile << "register_rule<"<<class_name<<"> register_"<<class_name
               << " ;" << endl ;
  syncFile(outputFile) ;

  if(!prettyOutput) {
    outputFile << "}" << endl ;
    syncFile(outputFile) ;
  }

  if(!use_prelude && sized_outputs && (rule_type != "apply")) 
    throw parseError("need prelude to size output type!") ;
}


void parseFile::skip_lpp_conditional(std::ostream &outputFile) {
  int nesting = 1 ;
  char c ;

  bool elseblock = false ;
  do {
    killsp() ;
    while(is.peek() != '$') {
      is.get(c) ;
      killsp() ;
      if(is.eof())
        throw parseError("Unexpected EOF processing $if directives") ;
    }
    is.get(c) ;
    if(is_name(is)) {
      string key = get_name(is) ;
      if((nesting == 1 && key == "else") || key == "endif") {
        nesting-- ;
      }
      if(key == "if" || key == "ifdef" || key == "ifndef") {
        nesting++ ;
      }
      elseblock = (key == "else") ;
    }
    killsp() ;
  } while(nesting != 0) ;
  if(elseblock) {
    outputFile << "// $else block output" << endl;
  }
}


void parseFile::processFile(string file, ostream &outputFile,
			    parseSharedInfo &parseInfo,int level) {
  bool error = false ;
  filename = file ;
  line_no = 1 ;
  ostringstream ss ;
  ss << "docvar_" ;
  for(size_t i=0;i<file.size();++i) {
    char c = file[i] ;
    if(isalpha(c) || isdigit(c) || c=='_')
      ss << c ;
    if(c == '.')
      break ;
  }
  pid_t pid = getpid() ;
  ss << "_" << pid ;
  
  docvarname = ss.str() ;

  parseInfo.fileNameStack.push_back(file) ;
  is.open(file.c_str(),ios::in) ;
  if(is.fail()) {
    for(auto li=include_dirs.begin();li!=include_dirs.end();++li) {
      string s = *li + "/" + file ;
      is.clear() ;
      is.open(s.c_str(),ios::in) ;
      if(!is.fail()) {
	parseInfo.dependFileList.push_back(s) ;
        break ;
      }
    }
    
    if(is.fail()) {
      string s = "can't open include file '" ;
      s += file ;
      s += "'" ;
      parseInfo.fileNameStack.pop_back() ;
      throw parseError(s) ;
    }
  } else {
    parseInfo.dependFileList.push_back(file) ;
  }
  char c ;
  
  if(level==0) {
    if(!parseInfo.no_cuda && parseInfo.debug_info>0) {
      outputFile << "#include <nvtx3/nvtx3.hpp>" << endl ;
      syncFile(outputFile) ;
    }
    outputFile << "extern const char *" << docvarname << "[] ;" << endl ;
  }
  syncFile(outputFile) ;
  bool mapset = false ;
  std::map<std::string,int> systemvarmap ;
    
  do {
    string comment = killspout(outputFile) ;
    //    cout << "comment:"<< comment << endl ;
    try {
      if(is.peek() == '$') { // Loci specific code!
        is.get(c) ; // get the $
        if(is.peek() == '[') {
          map<variable,string> vnames ;
          int openbrace = 0 ;
          process_SpecialCommand(outputFile,vnames,openbrace) ;
        } else  if(is_name(is)) {
          std::string key = get_name(is) ;
          if(key == "type") {
            setup_Type(outputFile,comment) ;
          } else if(key == "untype") {
            setup_Untype(outputFile) ;
          } else if(key == "rule") {
            if(level != 0) {
              throw parseError("$rule is not allowed in include file!") ;
            }
            setup_Rule(outputFile,comment,parseInfo) ;
	  } else if(key == "cudarule") {
            if(level != 0) {
              throw parseError("$rule is not allowed in include file!") ;
            }
	    if(parseInfo.no_cuda)
	      setup_Rule(outputFile,comment,parseInfo) ;
	    else
	      setup_cudaRule(outputFile,comment,parseInfo) ;
          } else if(key == "include") {
            killsp() ;
            if(!is_string(is)) {
              throw parseError("syntax error") ;
	    }
            string newfile = get_string(is) ;
	    // check to see if the file was already included, if so then
	    // don't perform the include repeatedly
            if(parseInfo.includedFiles.find(newfile) ==
	       parseInfo.includedFiles.end()) {
	      parseInfo.includedFiles.insert(newfile) ;
	      parseFile parser ;

	      parser.processFile(newfile,outputFile,parseInfo,level+1) ;
	      syncFile(outputFile) ;

              for(auto mi=parser.type_map.begin();mi!=parser.type_map.end();++mi)
                type_map[mi->first] = mi->second ;
              access_types = parser.access_types ;
              for(auto mi = parser.access_map.begin();mi!=parser.access_map.end();++mi)
              access_map[mi->first] = mi->second ;



	      for(auto mi=parser.type_map.begin();mi!=parser.type_map.end();++mi)
		type_map[mi->first] = mi->second ;
            } else {
	      outputFile << "//$include \"" << newfile << '"' << endl ;
	    }
          } else if(key == "define") {
            killsp() ;
            if(!is_name(is))
              throw parseError("syntax error parsing $define") ;
            string var = get_name(is) ;
            while(is.peek() == ' ' || is.peek() == '\t')
              is.get(c) ;
            string val ;
            if(is_string(is)) {
              val = get_string(is) ;
            }
            
            while((is.peek() == ' ') || (is.peek() == '\t')) {
              is.get(c) ;
            }
            if(is.peek() != '\r' && is.peek() != '\n')
              throw parseError("extraneous text after $define statement") ;
            
            setSystemVar(var,val) ;
            mapset = false ;
            outputFile << "// $define " << var << " \"" << val <<"\"" << endl ;
          } else if(key == "if" || key == "ifdef" || key == "ifndef") {
            killsp() ;
            if(key == "ifdef" && !is_name(is))
              throw parseError("syntax error parsing $ifdef") ;
            if(key == "ifndef" && !is_name(is))
              throw parseError("syntax error parsing $ifndef") ;
            bool check = false ;
            string var ;
            if(key == "ifdef") {
              var = get_name(is) ;
              check = checkSystemVar(var) ;
            } else if(key == "ifndef") {
              var = get_name(is) ;
              check = !checkSystemVar(var) ;
            } else {
              while(is.peek() != '\n' && is.peek() != '\r') {
                is.get(c) ;
                if(is.eof())
                  break ;
                // check if it might be a continuation line
                if(c == '\\') {
                  // it is continuation line, so eat \n \r characters
                  if(is.peek() == '\n' || is.peek() == '\r') {
                    is.get(c) ;
                    if(c == '\n') {
                      line_no++ ;
                    } else {
                      // Handle case of \r\n
                      if(is.peek() == '\n') {
                        is.get(c) ;
                        line_no++ ;
                      }
                    }
                  } else // otherwise add \ to string
                    var += '\\' ;
                } else { // not part of continuation line, add to expression
                  var += c ;
                }

              }
              try {
                exprP C = expression::create(var) ;
                if(!mapset)
                  evalSystemVars(systemvarmap) ;
                mapset = true ;
                check = C->evaluate(systemvarmap);
              } catch(exprError &err) {
                err.Print(cerr) ;
                throw parseError("syntax error parsing $if expression") ;
              }
            }

            if(!check) {
              // if check not true, parse forward until else or endif
              // keeping track of nesting
              skip_lpp_conditional(outputFile) ;
            } else {
              outputFile << "// $" << key << " " << var << " true" <<endl ;
            }
              syncFile(outputFile) ;
              
          } else if(key == "else") {
            killsp() ;
            // We got to the else which means the first part was true, so
            // skip nested part
            skip_lpp_conditional(outputFile) ;
            
            syncFile(outputFile) ;

          } else if(key == "endif") {
            killsp() ;
            outputFile << "// $endif"<< endl  ;
            syncFile(outputFile) ;
          } else {
            throw parseError("syntax error: unknown key") ;
          }
        } else {
          throw parseError("unable to process '$' command") ;
        }
      } else {
	bool foundComment = false ;
        while(is.peek() != '\n' && is.peek() != EOF) {
	  if(is_comment(is)) {
	    killCommentOut(is,line_no,outputFile) ;
	    foundComment = true ;
	    break ;
	  }

          is.get(c) ;
          outputFile << c ;
        }
	if(!foundComment) {
	  is.get(c) ;
	  outputFile << endl ;
	  line_no++ ;
	}
      }
      if(is.peek() == EOF)
        is.get(c) ;
    }
    catch(parseError pe) {
      cerr << filename << ':' << line_no << ": " << pe.error_type << endl ;
      char buf[512] ;
      is.getline(buf,512) ;
      line_no++ ;
      //      cerr << "remaining line = " << buf << endl ;
      error = true ;
    }

  } while(!is.eof()) ;
  parseInfo.fileNameStack.pop_back() ;
  if(error) 
    throw parseError("syntax error") ;


  if(access_types.size() > 0) {
    outputFile << endl << "const char *" << docvarname << "[] = {" << endl ;
    for(auto mi=access_types.begin();mi!=access_types.end();++mi) {
      outputFile << "\"" << mi->getFileLoc() << "\\000"
                 << mi->v << "\\000"
                 << mi->container << mi->container_args << "\\000"
                 << cleanupComment(mi->comment)<< "\","<< endl ;
    }
    outputFile << "\"\\000\\000\\000\" } ;" << endl ; 
  }
  
  if(parseInfo.fileNameStack.empty()) {
    outputFile << endl << "//DEPEND:" ;
    
    for(size_t i=1;i<parseInfo.dependFileList.size();++i)
      outputFile << ' ' << parseInfo.dependFileList[i] ;
    outputFile << endl ;
  }
}
