
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
      current_comment += ' ' ;
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
    if(c == '\n')
      c = ' ' ;
      
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

string cleanupCommentQuoted(const std::string& s) {
  auto begin = find_if_not(s.begin(), s.end(), spaceChar) ;
  auto end = find_if_not(s.rbegin(), s.rend(), spaceChar).base() ;

  while(begin < end) {
     if(*begin == '"') {
       ++begin ;
     } else if(*begin == '\\' && begin+1 < end && begin[1] == '"') {
       begin+=2 ;
     } else {
       break ;
     }
  }

  while(begin < end) {
    if(*(end-1) == '"') {
      --end ;
    } else if(end-2 >= begin && *(end-2) == '\\' && *(end-1) == '"') {
      end -= 2 ;
    } else {
      break ;
    }
  }

  string result ;
  result.reserve(distance(begin, end)) ;

  bool pendingSpace = false ;

  for(; begin < end; ++begin) {
    if(spaceChar(*begin)) {
      pendingSpace = true ;
      continue ;
    }

    if(pendingSpace && !result.empty())
      result += ' ';

    pendingSpace = false ;

    switch (*begin) {
    case '\\':
      result += R"(\\)" ;
      break ;
    case '"':
      result += R"(\")" ;
      break ;
    default:
      if (*begin >= ' ' && *begin <= '~')
        result += *begin ;
      break ;
    }
  }

  return result ;
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

bool get_block(
  istream & is, bool include_guard,
  char const open_char, char const close_char, string & contents,
  int & lines, ostream & err
) {
  enum class State {
    Code, String, Char, LineComment, BlockComment
  } ;

  if(is.peek() != open_char) {
    err << "expecting block start '" << open_char << "'" ;
    return false ;
  }

  int open_count = 1 ;
  contents.clear() ;
  lines = 0 ;

  char ch = is.get() ;
  if(include_guard) {
    contents += ch ;
  }

  State state = State::Code ;

  while(is) {
    ch = is.get() ;
    switch(state) {
    case State::Code:
      if(ch == EOF) {
        err << "unexpected end of block" ;
        return false ;
      } else if(ch == '"') {
        contents += ch ;
        state = State::String ;
      } else if(ch == '\'') {
        contents += ch ;
        state = State::Char ;
      } else if(ch == '/') {
        contents += ch ;
        if(is.peek() == '/') {
          ch = is.get() ;
          contents += ch ;
          state = State::LineComment ;
        } else if(is.peek() == '*') {
          ch = is.get() ;
          contents += ch ;
          state = State::BlockComment ;
        }
      } else if(ch == open_char) {
        contents += ch ;
        ++open_count ;
      } else if(ch == close_char) {
        --open_count ;
        if(open_count == 0) {
          if(include_guard) {
            contents += close_char ;
          }
          return true ;
        }
        contents += close_char ;
      } else if(ch == '\n') {
        contents += ch ;
        ++lines ;
      } else {
        contents += ch ;
      }
      break ;
    case State::String:
      contents += ch ;
      if(contents.back() == '\\') {
        if(is.peek() == EOF) {
          err << "unexpected end of block" ;
          return false ;
        }
        ch = is.get() ;
        contents += ch ;
      } else if(contents.back() == '"') {
        state = State::Code ;
      }
      break ;
    case State::Char:
      contents += ch ;
      if(contents.back() == '\\') {
        if(is.peek() == EOF) {
          err << "unexpected end of block" ;
          return false ;
        }
        contents += ch ;
      } else if(contents.back() == '\'') {
        state = State::Code ;
      }
      break ;
    case State::LineComment:
      contents += ch ;
      if(contents.back() == '\n') {
        ++lines ;
        state = State::Code ;
      }
      break ;
    case State::BlockComment:
      contents += ch ;
      if(contents.back() == '*' && is.peek() == '/') {
        ch = is.get() ;
        contents += ch ;
        state = State::Code ;
      }
      break ;
    }
  }

  err << "unexpected end of block" ;
  return false ;
}

class nestedbracestuff : public parsebase {
  string brace_contents ;

public:
  istream & get(istream & s) {
    parsebase::killsp(s) ;
    ostringstream err ;
    err << "error parsing { ... } block: " ;
    if(!get_block(s, true, '{', '}', brace_contents, lines, err)) {
      throw parseError(err.str()) ;
    }
    return s ;
  }

  string str() {
    return brace_contents ;
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
  virtual void visit(AST_BlockRaw &) ;
  virtual void visit(AST_typeSpec &) ;
  virtual void visit(AST_declaration &) ;
  virtual void visit(AST_SimpleStatement &) ;
  virtual void visit(AST_controlStatement &) ;
  virtual void visit(AST_LociDirective &) ;
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

void AST_printTree::visit(AST_BlockRaw &s) {
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

void AST_printTree::visit(AST_LociDirective &s) {
  pushindent(s) ;
  out << "[[" ;
  s.type->accept(*this) ;
  out << "]][[" ;
  s.body->accept(*this) ;
  out << "]]" ;
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
  
  AST_editLociVariableAccess(
    const std::map<variable,std::string> &vnames_in,
    const std::map<variable,std::string> &vtypes_in
  ) : vnames(vnames_in), vtypes(vtypes_in) {
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
      bool is_param_like = false ;
      CPTR<AST_Token> tok(op.terms[i]) ;
      variable v(tok->text) ;
      auto t = vtypes.find(v) ;
      if(t != vtypes.end()) {
        if(t->second == "param" ||
           t->second == "blackbox" ||
           t->second == "Constraint" ||
           t->second == "constraint") {
          is_param_like = true ;
        }
      }

      if(is_param_like) {
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

class AST_editGPULociVariableAccess : public AST_visitor {
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
  
  AST_editGPULociVariableAccess(
    const std::map<variable,std::string> &vnames_in,
    const std::map<variable,std::string> &vtypes_in
  ) : vnames(vnames_in), vtypes(vtypes_in) {
    CPTR<AST_Token> e = new AST_Token ;
    e->lineno = -1 ;
    e->text = "_e_" ;
    e->nodeType = nodeTypes::TK_NAME ;
    entityIndex = AST_type::ASTP(e) ;
  }

  virtual void visit(AST_exprOper &) ;
} ;



void AST_editGPULociVariableAccess::visit(AST_exprOper &op) {
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

class AST_editLociDirective : public AST_visitor {
  string filename_ ;

public:
  AST_editLociDirective(string const & filename) : filename_(filename) { }
  virtual void visit(AST_Block & s) ;
  virtual void visit(AST_controlStatement & s) ;
  AST_type::ASTP createLociDirectiveReplacement(AST_type::ASTP e) ;
} ;

AST_type::ASTP AST_editLociDirective::createLociDirectiveReplacement(AST_type::ASTP e) {
  CPTR<AST_LociDirective> dir(e) ;
  CPTR<AST_Token> dir_type_tok(dir->type) ;
  CPTR<AST_Block> dir_body(dir->body) ;

  string type ;
  for(auto c: dir_type_tok->text) {
    type += tolower(c) ;
  }

  if(type == "once") {
    CPTR<AST_Token> ns = new AST_Token ;
    ns->nodeType = nodeTypes::TK_NAME ;
    ns->text = "Loci" ;
    ns->lineno = dir_type_tok->lineno ;

    CPTR<AST_Token> funname = new AST_Token ;
    funname->nodeType = nodeTypes::TK_NAME ;
    funname->text = "is_leading_execution" ;
    funname->lineno = dir_type_tok->lineno ;

    CPTR<AST_exprOper> scoped = new AST_exprOper ;
    scoped->nodeType = nodeTypes::OP_SCOPE ;
    scoped->terms.push_back(AST_type::ASTP(ns)) ;
    scoped->terms.push_back(AST_type::ASTP(funname)) ;

    CPTR<AST_exprOper> funcall = new AST_exprOper ;
    funcall->nodeType = nodeTypes::OP_FUNC ;
    funcall->terms.push_back(AST_type::ASTP(scoped)) ;

    CPTR<AST_exprOper> conditional = new AST_exprOper ;
    conditional->nodeType = nodeTypes::OP_GROUP ;
    conditional->terms.push_back(AST_type::ASTP(funcall)) ;

    CPTR<AST_Token> if_tok = new AST_Token ;
    if_tok->nodeType = nodeTypes::TK_IF ;
    if_tok->text = "if" ;
    if_tok->lineno = dir_type_tok->lineno ;

    AST_type::ASTP else_tok = 0 ;
    AST_type::ASTP else_body = 0 ;

    CPTR<AST_controlStatement> stmt = new AST_controlStatement ;
    stmt->constructIf(
      AST_type::ASTP(if_tok),AST_type::ASTP(conditional),
      dir->body,else_tok,else_body
    ) ;

    return AST_type::ASTP(stmt) ;
  } else if(type == "atomic") {
    CPTR<AST_Token> type_scope1 = new AST_Token ;
    type_scope1->nodeType = nodeTypes::TK_NAME ;
    type_scope1->text = "Loci" ;
    type_scope1->lineno = dir_type_tok->lineno ;

    CPTR<AST_Token> type_scope2 = new AST_Token ;
    type_scope2->nodeType = nodeTypes::TK_NAME ;
    type_scope2->text = "atomic_region_helper" ;
    type_scope2->lineno = dir_type_tok->lineno ;

    CPTR<AST_exprOper> type = new AST_exprOper ;
    type->nodeType = nodeTypes::OP_SCOPE ;
    type->terms.push_back(AST_type::ASTP(type_scope1)) ;
    type->terms.push_back(AST_type::ASTP(type_scope2)) ;

    CPTR<AST_Token> varname = new AST_Token ;
    varname->nodeType = nodeTypes::TK_NAME ;
    varname->text = "L__ATOMIC_REGION" ;
    varname->lineno = dir_type_tok->lineno ;

    CPTR<AST_Token> term = new AST_Token ;
    term->nodeType = nodeTypes::TK_SEMICOLON ;
    term->text = ";" ;
    term->lineno = dir_type_tok->lineno ;

    CPTR<AST_declaration> atomic_region = new AST_declaration ;
    atomic_region->type_decl.push_back(AST_type::ASTP(type)) ;
    atomic_region->decls.push_back(AST_type::ASTP(varname)) ;
    atomic_region->decls.push_back(AST_type::ASTP(term)) ;

    dir_body->identifiers[varname->text] = localIdentifier() ;

    dir_body->elements.insert(
      dir_body->elements.begin()+1, AST_type::ASTP(atomic_region)
    ) ;
    return AST_type::ASTP(dir_body) ;
  }

  return e ;
}

void AST_editLociDirective::visit(AST_Block & s) {
  size_t esz = s.elements.size() ;
  for(size_t i = 0; i < esz; ++i) {
    if(ASTEqual(s.elements[i], nodeTypes::OP_LOCI_DIRECTIVE)) {
      s.elements[i] = createLociDirectiveReplacement(s.elements[i]) ;
    } else {
      s.elements[i]->accept(*this) ;
    }
  }
}

void AST_editLociDirective::visit(AST_controlStatement & s) {
  size_t psz = s.parts.size() ;
  for(size_t i = 0; i < psz; ++i) {
    if(ASTEqual(s.parts[i], nodeTypes::OP_LOCI_DIRECTIVE)) {
      s.parts[i] = createLociDirectiveReplacement(s.parts[i]) ;
    } else {
      s.parts[i]->accept(*this) ;
    }
  }
}

class AST_editPrelude : public AST_visitor {
public:
  std::string filename ;
  const std::map<variable,std::string> &vnames ;
  AST_editPrelude(
    std::string const & filename_in,
    const std::map<variable,std::string> &vnames_in
  ) : filename(filename_in), vnames(vnames_in) { }
  virtual void visit(AST_Token &) ;
  virtual void visit(AST_BlockRaw &) ;
} ;

void AST_editPrelude::visit(AST_Token &s) {
  if(ASTEqual(s, nodeTypes::TK_LOCI_VARIABLE)) {
    variable v(s.text) ;
    auto vmi = vnames.find(v) ;
    if(vmi == vnames.end()) {
      ostringstream ss ;
      ss << "type error: variable " << v << " is unknown to this rule!" << endl ;
      throw parseError(ss.str()) ;
    }

    s.text = vmi->second ;
    s.nodeType = nodeTypes::TK_NAME ;
  }
}

void AST_editPrelude::visit(AST_BlockRaw &s) {
  using namespace nodeTypes ;

  AST_type::ASTList new_elements ;
  for(auto ii = s.elements.begin(); ii != s.elements.end(); ++ii) {
    if((*ii)->nodeType == OP_LOCI_DIRECTIVE) {
      CPTR<AST_LociDirective> dir(*ii) ;

      if(!ASTEqual(dir->body, ND_BLOCK_RAW)) {
        ostringstream ss ;
        ss << "error processing prelude at " << filename << ":"
           << dir->type->lineno
           << ": expected raw block for body of Loci directive" ;
        throw parseError(ss.str()) ;
      }

      CPTR<AST_BlockRaw> dir_body(dir->body) ;

      string type ;
      for(auto c : dir->type->text) {
        type += tolower(c) ;
      }

      if(type == "once") {
        CPTR<AST_Token> if_tok = new AST_Token ;
        if_tok->nodeType = TK_IF ;
        if_tok->text = "if" ;
        if_tok->lineno = dir->type->lineno ;

        CPTR<AST_Token> open_tok = new AST_Token ;
        open_tok->nodeType = TK_OPENPAREN ;
        open_tok->text = "(" ;
        open_tok->lineno = dir->type->lineno ;

        CPTR<AST_Token> loci_tok = new AST_Token ;
        loci_tok->nodeType = TK_NAME ;
        loci_tok->text = "Loci" ;
        loci_tok->lineno = dir->type->lineno ;

        CPTR<AST_Token> scope_tok = new AST_Token ;
        scope_tok->nodeType = TK_SCOPE ;
        scope_tok->text = "::" ;
        scope_tok->lineno = dir->type->lineno ;

        CPTR<AST_Token> le_tok = new AST_Token ;
        le_tok->nodeType = TK_NAME ;
        le_tok->text = "is_leading_execution" ;
        le_tok->lineno = dir->type->lineno ;

        CPTR<AST_Token> close_tok = new AST_Token ;
        close_tok->nodeType = TK_CLOSEPAREN ;
        close_tok->text = ")" ;
        close_tok->lineno = dir->type->lineno ;

        dir->body->accept(*this) ;

        new_elements.push_back(AST_type::ASTP(if_tok)) ;
        new_elements.push_back(AST_type::ASTP(open_tok)) ;
        new_elements.push_back(AST_type::ASTP(loci_tok)) ;
        new_elements.push_back(AST_type::ASTP(scope_tok)) ;
        new_elements.push_back(AST_type::ASTP(le_tok)) ;
        new_elements.push_back(open_tok->clone()) ;
        new_elements.push_back(close_tok->clone()) ;
        new_elements.push_back(AST_type::ASTP(close_tok)) ;
        new_elements.push_back(dir->body) ;
      } else if(type == "atomic") {
        CPTR<AST_Token> loci_tok = new AST_Token ;
        loci_tok->nodeType = TK_NAME ;
        loci_tok->text = "Loci" ;
        loci_tok->lineno = dir->type->lineno ;

        CPTR<AST_Token> scope_tok = new AST_Token ;
        scope_tok->nodeType = TK_SCOPE ;
        scope_tok->text = "::" ;
        scope_tok->lineno = dir->type->lineno ;

        CPTR<AST_Token> arh_tok = new AST_Token ;
        arh_tok->nodeType = TK_NAME ;
        arh_tok->text = "atomic_region_helper" ;
        arh_tok->lineno = dir->type->lineno ;

        CPTR<AST_Token> var_tok = new AST_Token ;
        var_tok->nodeType = TK_NAME ;
        var_tok->text = "L__ATOMIC_REGION" ;
        var_tok->lineno = dir->type->lineno ;

        CPTR<AST_Token> term_tok = new AST_Token ;
        term_tok->nodeType = TK_SEMICOLON ;
        term_tok->text = ";" ;
        term_tok->lineno = dir->type->lineno ;

        CPTR<AST_BlockRaw> new_block = new AST_BlockRaw ;
        size_t sz = dir_body->elements.size() ;
        new_block->elements.push_back(dir_body->elements[0]) ;
        new_block->elements.push_back(AST_type::ASTP(loci_tok)) ;
        new_block->elements.push_back(AST_type::ASTP(scope_tok)) ;
        new_block->elements.push_back(AST_type::ASTP(arh_tok)) ;
        new_block->elements.push_back(AST_type::ASTP(var_tok)) ;
        new_block->elements.push_back(AST_type::ASTP(term_tok)) ;
        for(size_t i = 1; i < sz; ++i) {
          new_block->elements.push_back(dir_body->elements[i]) ;
        }
        (*ii) = AST_type::ASTP(new_block) ;

        (*ii)->accept(*this) ;
        new_elements.push_back(*ii) ;
      } else {
        ostringstream ss ;
        ss << "error processing prelude at " << filename << ":"
           << dir->type->lineno
           << ": unsupported Loci directive '" << dir->type->text << "'" ;
        throw parseError(ss.str()) ;
      }
    } else {
      (*ii)->accept(*this) ;
      new_elements.push_back(*ii) ;
    }
  }

  s.elements.swap(new_elements) ;
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
  AST_editLociVariableAccess AST_editor(vnames, vnames) ;
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
  int rule_type_line_no = 0 ;
  int signature_line_no = 0 ;
  int apply_op_line_no = 0 ;
  int constraint_line_no = 0 ;
  int parametric_line_no = 0 ;
  int conditional_line_no = 0 ;
  int specialized_line_no = 0 ;
  vector<int> comments_line_no ;
  int compute_line_no = 0 ;

  killsp() ;
  string rule_type ;
  if(is_name(is)) {
    rule_type_line_no = line_no ;
    rule_type = get_name(is) ;
  } else 
    throw parseError("syntax error") ;

  nestedparenstuff signature ;
  signature.get(is) ;
  signature_line_no = line_no ;
  line_no += signature.num_lines() ;

  nestedbracketstuff apply_op ;
  killsp() ;
  if(rule_type == "apply") {
    if(is.peek() != '[')
      throw parseError("apply rule missing '[operator]'") ;
    apply_op.get(is) ;
    apply_op_line_no = line_no ;
    line_no += apply_op.num_lines() ;
    killsp() ;
  }

  string constraint, conditional ;
  string parametric_var ;
  list<string> options ;
  vector<string> comments ;
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
      constraint_line_no = line_no ;
      line_no += con.num_lines() ;
    } else if(s == "parametric") {
      nestedparenstuff con ;
      con.get(is) ;
      if(parametric_var != "") {
        throw parseError("syntax error: cannot specify more than one parametric variable") ;
      }

      parametric_var = con.str() ;
      parametric_line_no = line_no ;
      line_no += con.num_lines() ;
    } else if(s == "conditional") {
      nestedparenstuff con ;
      con.get(is) ;
      if(conditional != "") {
        throw parseError("syntax error: cannot specify more than one conditional variable") ;
      }
      conditional = con.str() ;
      conditional_line_no = line_no ;
      line_no += con.num_lines() ;
      // Check variable
      variable cond(conditional) ;
      auto mi  = lookupVarType(cond) ;

      if(!checkTypeValid(mi)) {
        cerr << filename << ':' << line_no
             << ":0: warning: type of conditional variable '"
             << cond << "' not found!"  << endl  ;
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
      specialized_line_no = line_no ;
      line_no += ip.num_lines() ;
      options.push_back(ip.str()) ;
    } else if(s == "comments") {
      nestedparenstuff ip ;
      ip.get(is) ;
      comments_line_no.push_back(line_no) ;
      line_no += ip.num_lines() ;
      comments.push_back(ip.str()) ;
    } else {
      throw parseError("unknown rule modifier") ;
    }
    killsp() ;
  }

  if(use_prelude) {
    throw parseError("prelude not compatible with cuda rules") ;
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

  // Inputs cannot have priority.
  for(auto vi=input.begin();vi!=input.end();++vi) {
    if(vi->get_info().priority.size() != 0) {
      ostringstream oss ;
      oss<< "improper use of priority annotation on rule input, var=" << *vi << endl ;
      throw parseError(oss.str()) ;
    }
  }

  // Only pointwise rules can have priority over output variables.
  if(rule_type != "pointwise" && rule_type != "default") {
    for(auto vi = output.begin(); vi != output.end(); ++vi) {
      if(vi->get_info().priority.size() != 0) {
        ostringstream oss ;
        oss << "only pointwise rules can use priority annotation, var=" << *vi << endl ;
        throw parseError(oss.str()) ;
      }
    }
  }

  // Set of inputs + outputs.
  variableSet all_vars = input;
  all_vars += output ;

  // Catch undeclared Loci variable.
  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      string s ;
      s = "unable to determine type of variable " ;
      s += (*vi).str() ;
      throw parseError(s) ;
    }
  }

  // Create C++ variable names for Loci variables.
  map<variable,string> vnames ;
  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    vnames[*vi] = var2name(*vi) ;
    if(vi->get_info().priority.size() != 0) {
      variable v = *vi ;
      while(v.get_info().priority.size() != 0)
        v = v.drop_priority() ;
      vnames[v] = vnames[*vi] ;
    }
  }
  for(auto ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    vnames[ipi->first] = vnames[ipi->second] ;
  }

  // Check if variables paired by inplace specification are specified as either
  // input or output to the rule.
  variableSet checkset ;
  for(auto vi=all_vars.begin();vi!=all_vars.end();++vi) {
    variable v = *vi ;
    while(v.get_info().priority.size() != 0)
      v = v.drop_priority() ;
    checkset += v ;
  }
  for(auto ipi=inplace.begin();ipi!=inplace.end();++ipi) {
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

  if(rule_type == "pointwise") {
    for(auto vi=output.begin();vi!=output.end();++vi) {
      auto mi = lookupVarType(*vi) ;
      if(mi->second.container == "param" && vi->get_info().name != "OUTPUT") {
        throw(parseError("pointwise rule cannot compute param, use singleton")) ;
      }
    }
  }

  bool paramUnit = 0 ;
  if(rule_type == "unit") {
    if(output.size() != 1) {
      throw parseError("unit rule should have only one output variable") ;
    }

    variable v = *(output.begin()) ;
    typedoc tinfo = lookupVarType(v)->second ;
    if(tinfo.container != "param") {
      throw parseError("unit rule should have only param output variable") ;
    }

    if(tinfo.container == "param") {
      paramUnit = 1 ;
    }
  }

  int paramApply = 0 ;
  int singletonApply = 0 ;
  if(rule_type == "apply") {
    if(output.size() != 1) {
      throw parseError("apply rule should have only one output variable") ;
    }

    variable v = *(output.begin()) ;
    typedoc tinfo = lookupVarType(v)->second ;
    if(tinfo.container != "param") {
      throw parseError("apply rule should have param as output variable") ;
    }

    bool allparam = true ;
    for(auto vi = input.begin(); vi != input.end(); ++vi) {
      typedoc tinfo2 = lookupVarType(*vi)->second ;
      if(tinfo2.container != "param") {
        allparam = false ;
      }
    }

    if(tinfo.container == "param") {
      paramApply = 1 ;
    }

    if(allparam) {
      singletonApply = 1 ;
    }
  }

  if(singletonApply) {
    cerr << "NOTE: parameter only apply rule on '"
         << *(output.begin()) << "' now executes single instance." << endl ;
  }

  // Extract constraint variables and check if they are typed.
  variableSet constraint_vars ;
  set<vmap_info> constraints ;
  if(constraint != "") {
    exprP C = expression::create(constraint) ;
    fill_descriptors(constraints, collect_associative_op(C,OP_COMMA)) ;

    for(auto i = constraints.begin(); i != constraints.end(); ++i) {
      for(size_t j = 0; j < i->mapping.size(); ++j)
        constraint_vars += i->mapping[j] ;
      constraint_vars += i->var ;
    }

    for(auto vi = constraint_vars.begin(); vi!=constraint_vars.end(); ++vi) {
      auto mi = lookupVarType(*vi) ;

      if(!checkTypeValid(mi)) {
        cerr << filename << ':' << constraint_line_no
             << ":0: warning: type of constraint variable '"
             << *vi << "' not found!"  << endl  ;
      }
    }
  }

  varmap typemap ;

  if(is.peek() != '{') {
    throw parseError("syntax error, expecting '{'") ;
  }

  compute_line_no = line_no ;

  CPTR<AST_type> ap = parseBlock(is, line_no, filename, typemap) ;

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
    AST_simplePrint printer(cerr, -1, false) ;
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

  for(auto i = varaccess.accessed.begin(); i != varaccess.accessed.end(); ++i) {
    readvars += i->var ;
    for(size_t j = 0; j < i->mapping.size(); ++j) {
      readvars += i->mapping[j] ;
    }
  }
  for(auto i = varaccess.writes.begin(); i != varaccess.writes.end(); ++i) {
    writevars += i->var ;
    for(size_t j = 0; j < i->mapping.size(); ++j) {
      readvars += i->mapping[j] ;
    }
  }

  readvars -= writevars ;

  // Now remove and save the open and close braces in the parseBlock
  CPTR<AST_Block> bigblock = CPTR<AST_Block>(ap) ;
  CPTR<AST_type> open = bigblock->elements[0] ;
  int bsz = bigblock->elements.size() ;
  CPTR<AST_type> close = bigblock->elements[bsz-1] ;
  for(int i = 0; i < bsz-1; ++i) {
    bigblock->elements[i] = bigblock->elements[i+1] ;
  }
  bigblock->elements.pop_back() ;
  bigblock->elements.pop_back() ;

  if(rule_type == "apply") {
    std::stringstream ss ;
    int linecount = compute_line_no ;

    ss << "loci_reduction_t _reduce_op_ ;" ;
    linecount = compute_line_no ;
    AST_type::ASTP reduce_op_decl = parseDeclaration(
      ss, linecount, filename, typemap
    ) ;

    ss.str("") ;
    ss.clear() ;
    ss << "value_t _local_ = _reduce_op_.identity() ;" ;
    linecount = compute_line_no ;
    AST_type::ASTP local_value_decl = parseDeclaration(
      ss, linecount, filename, typemap
    ) ;

    ss.str("") ;
    ss.clear() ;
    ss << "return _local_ ;" ;
    linecount = compute_line_no ;
    AST_type::ASTP return_local_value = parseSpecialControlStatement(
      ss, linecount, filename, typemap
    ) ;

    bigblock->elements.insert(bigblock->elements.begin(), local_value_decl) ;
    bigblock->elements.insert(bigblock->elements.begin(), reduce_op_decl) ;
    bigblock->elements.insert(bigblock->elements.end(), return_local_value) ;

    AST_editJoin edit_join ;
    ap->accept(edit_join) ;
  }

  //AST_printObjectTree treeout(cerr) ;
  //ap->accept(treeout) ;

  variableSet outs = output ;
  for(auto ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    outs -= ipi->first ;
    outs += ipi->second ;
  }
  variableSet ins = input ;
  ins -= outs ;

  map<variable,string> typetable ;
  map<variable,string> ctypetable ;
  map<variable,string> cargtable ;
  for(auto vi=ins.begin();vi!=ins.end();++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }

    if(mi->second.container != "param" &&
       mi->second.container != "Map" &&
       mi->second.container != "MapVec" &&
       mi->second.container != "multiMap" &&
       mi->second.container != "store" &&
       mi->second.container != "storeVec" &&
       mi->second.container != "storeMat" &&
       mi->second.container != "multiStore") {
       cerr << "unknown container type '" << mi->second.container
            << "' for variable " << *vi << endl ;
       throw parseError("unsupported Loci container") ;
    }

    ctypetable[*vi] = mi->second.container ;

    if(mi->second.container == "Map") {
      typetable[*vi] = "int const *" ;
      cargtable[*vi] = "" ;
    } else if(mi->second.container == "multiMap") {
      typetable[*vi] = "Loci::constMultiAccessor<Loci::Entity> " ;
      cargtable[*vi] = "" ;
    } else if(mi->second.container == "MapVec") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = string("Array<Entity,") + arg + "> const *" ;
        cargtable[*vi] = arg ;
      } else {
        cerr << "unexpected loci variable type!" << endl ;
      }
    } else if(mi->second.container == "multiStore") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = "Loci::constMultiAccessor<" + arg + "> ";
        cargtable[*vi] = arg ;
      } else {
        cerr << "unexpected loci variable type!" << endl ;
      }
    } else {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = arg + " const *";
        cargtable[*vi] = arg ;
      } else {
        cerr << "unexpected loci variable type!" << endl ;
      }
    }
  }

  for(auto vi=outs.begin();vi!=outs.end();++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      cerr << "unknown type for variable " << *vi << endl ;
      throw parseError("untyped Loci variable") ;
    }

    if(mi->second.container != "param" &&
       mi->second.container != "Map" &&
       mi->second.container != "MapVec" &&
       mi->second.container != "store" &&
       mi->second.container != "storeVec" &&
       mi->second.container != "storeMat" &&
       mi->second.container != "multiStore") {
       cerr << "unknown container type '" << mi->second.container
            << "' for variable " << *vi << endl ;
       throw parseError("unsupported Loci container") ;
    }

    ctypetable[*vi] = mi->second.container ;

    if(mi->second.container == "Map") {
      typetable[*vi] = "int *" ;
      cargtable[*vi] = "" ;
    } else if(mi->second.container == "MapVec") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = string("Array<Entity,") + arg + ">  *" ;
        cargtable[*vi] = arg ;
      } else {
        cerr << "unexpected loci variable type!" << endl ;
      }
    } else {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = arg + " *";
        cargtable[*vi] = arg ;
      } else {
        cerr << "unexpected loci variable type!" << endl ;
      }
    }
  }

  for(auto ipi=inplace.begin();ipi!=inplace.end();++ipi) {
    all_vars -= ipi->first ;
  }

  if(constraint != "") {
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
        cerr << filename << ':' << line_no
             << ":0: warning: type of constraint variable '"
             << *vi << "' not found!"  << endl  ;
      }
    }
  }

  if(comments.size() == 0) {
    // check to see if there is a javadoc compatible comment before
    if(comment.size() > 0) {
      comments.push_back(comment) ;
      comments_line_no.push_back(rule_type_line_no) ;
    } else if(rule_type=="optional" || rule_type=="default") {
      auto mi = lookupVarType(*output.begin()) ;
      comments.push_back(mi->second.comment) ;
      comments_line_no.push_back(mi->second.lineno) ;
    }
  }

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
    throw parseError("cuda rules currently incompatible with storeVec, storeMat or multiStore types") ;

  AST_editGPULociVariableAccess AST_editor(vnames, ctypetable) ;
  ap->accept(AST_editor) ;

  string rule_debug_name ;
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

    rule_debug_name = oss.str() ;
  }

  DictionaryTemplateValue rule_ctx ;

  rule_ctx["type"] = rule_type ;
  rule_ctx["class"] = class_name ;
  rule_ctx["file"] = filename ;
  rule_ctx["line_number"] = rule_type_line_no ;
  rule_ctx["debug_name"] = rule_debug_name ;

  {
    DictionaryTemplateValue signature_ctx ;
    signature_ctx["line_number"] = signature_line_no ;
    rule_ctx["signature"] = signature_ctx ;
  }

  {
    ArrayTemplateValue input_stores_ctx ;
    for(auto vi = ins.begin(); vi != ins.end(); ++vi) {
      DictionaryTemplateValue ctx ;
      ctx["name"] = (*vi).str() ;
      ctx["vname"] = vnames[*vi] ;
      ctx["ctype"] = ctypetable[*vi] ;
      ctx["vtype"] = typetable[*vi] ;
      ctx["carg"] = cargtable[*vi] ;
      input_stores_ctx.append(ctx) ;
    }
    rule_ctx["input_stores"] = input_stores_ctx ;
  }

  {
    ArrayTemplateValue output_stores_ctx ;
    for(auto vi = outs.begin(); vi != outs.end(); ++vi) {
      DictionaryTemplateValue ctx ;
      ctx["name"] = (*vi).str() ;
      ctx["vname"] = vnames[*vi] ;
      ctx["ctype"] = ctypetable[*vi] ;
      ctx["vtype"] = typetable[*vi] ;
      ctx["carg"] = cargtable[*vi] ;
      output_stores_ctx.append(ctx) ;
    }
    rule_ctx["output_stores"] = output_stores_ctx ;
  }

  {
    ArrayTemplateValue name_stores_ctx ;
    for(auto vi = all_vars.begin(); vi != all_vars.end(); ++vi) {
      DictionaryTemplateValue ctx ;
      ctx["name"] = (*vi).str() ;
      ctx["vname"] = vnames[*vi] ;

      auto mi = access_map.find(lookupVarType(*vi)->second.getFileLoc()) ;
      if(mi != access_map.end()) {
        ctx["has_info_id"] = 1 ;
        ctx["info_id"] = mi->second ;
      } else {
        ctx["has_info_id"] = 0 ;
      }

      name_stores_ctx.append(ctx) ;
    }
    rule_ctx["name_stores"] = name_stores_ctx ;
  }

  {
    ArrayTemplateValue inputs_ctx ;
    for(auto i = sources.begin(); i != sources.end(); ++i) {
      ostringstream ss ;
      if(i->mapping.size() > 1) {
        ss << '(' ;
      }
      for(auto vi = i->mapping.begin(); vi != i->mapping.end(); ++vi) {
        if(vi != i->mapping.begin()) {
          ss << ',' ;
        }
        ss << *vi ;
      }
      if(i->mapping.size() > 1) {
        ss << ')' ;
      }
      if(i->mapping.size() > 0) {
        ss << "->" ;
      }

      if(i->var.size() > 1) {
        ss << '(' ;
      }
      for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
        if(vi != i->var.begin()) {
          ss << ',' ;
        }
        ss << *vi ;
      }
      if(i->var.size() > 1) {
        ss << ')' ;
      }
      DictionaryTemplateValue ctx ;
      ctx["str"] = ss.str() ;
      inputs_ctx.append(ctx) ;
    }
    rule_ctx["inputs"] = inputs_ctx ;
  }

  {
    ArrayTemplateValue outputs_ctx ;
    for(auto i = targets.begin(); i != targets.end(); ++i) {
      ostringstream ss ;
      if(i->mapping.size() > 1) {
        ss << '(' ;
      }
      for(auto vi = i->mapping.begin(); vi != i->mapping.end(); ++vi) {
        if(vi != i->mapping.begin()) {
          ss << ',' ;
        }
        ss << *vi ;
      }
      if(i->mapping.size() > 1) {
        ss << ')' ;
      }
      if(i->mapping.size() > 0) {
        ss << "->" ;
      }

      if(i->var.size() > 1) {
        ss << '(' ;
      }
      for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
        if(vi != i->var.begin()) {
          ss << ',' ;
        }

        auto ipi = inplace.begin() ;
        while(ipi != inplace.end()) {
          if(ipi->first == *vi) break ;
          ++ipi ;
        }
        if(ipi != inplace.end()) {
          if(i->mapping.size() == 0 || i->var.size() > 1) {
            ss << ipi->first << '=' << ipi->second ;
          } else {
            ss << '(' << ipi->first << '=' << ipi->second << ')' ;
          }
        } else {
          ss << *vi ;
        }
      }
      if(i->var.size() > 1) {
        ss << ')' ;
      }

      DictionaryTemplateValue ctx ;
      ctx["str"] = ss.str() ;
      outputs_ctx.append(ctx) ;
    }

    rule_ctx["outputs"] = outputs_ctx ;
  }

  {
    ArrayTemplateValue constraint_spec_ctx ;
    for(auto i = constraints.begin(); i != constraints.end(); ++i) {
      ostringstream ss ;
      if(i->mapping.size() > 1) {
        ss << '(' ;
      }
      for(auto vi = i->mapping.begin(); vi != i->mapping.end(); ++vi) {
        if(vi != i->mapping.begin()) {
          ss << ',' ;
        }
        ss << *vi ;
      }
      if(i->mapping.size() > 1) {
        ss << ')' ;
      }
      if(i->mapping.size() > 0) {
        ss << "->" ;
      }

      if(i->var.size() > 1) {
        ss << '(' ;
      }
      for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
        if(vi != i->var.begin()) {
          ss << "," ;
        }
        ss << *vi ;
      }
      if(i->var.size() > 1) {
        ss << ')' ;
      }

      DictionaryTemplateValue ctx ;
      ctx["str"] = ss.str() ;
      constraint_spec_ctx.append(ctx) ;
    }

    DictionaryTemplateValue constraints_ctx ;
    constraints_ctx["spec"] = constraint_spec_ctx ;
    constraints_ctx["line_number"] = constraint_line_no ;
    rule_ctx["constraints"] = constraints_ctx ;
  }

  rule_ctx["option_disable_threading"] = 1 ;

  {
    DictionaryTemplateValue parametric_ctx ;
    if(parametric_var != "") {
      rule_ctx["is_parametric"] = 1 ;
      parametric_ctx["spec"] = parametric_var ;
      parametric_ctx["line_number"] = parametric_line_no ;
    } else {
      rule_ctx["is_parametric"] = 0 ;
    }
    rule_ctx["parametric"] = parametric_ctx ;
  }

  {
    DictionaryTemplateValue specialized_ctx ;
    if(is_specialized) {
      specialized_ctx["line_number"] = specialized_line_no ;
    }
    rule_ctx["is_specialized"] = is_specialized ;
    rule_ctx["specialized"] = specialized_ctx ;
  }

  {
    DictionaryTemplateValue conditional_ctx ;

    if(conditional != "") {
      rule_ctx["is_conditional"] = 1 ;
      conditional_ctx["spec"] = conditional ;
      conditional_ctx["line_number"] = conditional_line_no ;
    } else {
      rule_ctx["is_conditional"] = 0 ;
    }

    rule_ctx["conditional"] = conditional_ctx ;
  }

  {
    ArrayTemplateValue comments_ctx ;

    size_t size = comments.size() ;
    for(size_t i = 0; i < size; ++i) {
      DictionaryTemplateValue ctx ;
      ctx["str"] = comments[i] ;
      ctx["line_number"] = comments_line_no[i] ;
      comments_ctx.append(ctx) ;
    }

    rule_ctx["comments"] = comments_ctx ;
  }

  if(rule_type == "pointwise") {
    rule_ctx["is_pointwise"] = 1 ;
    rule_ctx["is_unit"] = 0 ;
    rule_ctx["is_apply"] = 0 ;
  } else if(rule_type == "unit") {
    variable unit_var = *(output.begin()) ;

    DictionaryTemplateValue unit_ctx ;
    unit_ctx["target_name"] = unit_var.str() ;
    unit_ctx["target_vname"] = vnames[unit_var] ;
    unit_ctx["container"] = ctypetable[unit_var] ;
    unit_ctx["container_args"] = cargtable[unit_var] ;
    unit_ctx["is_param"] = paramUnit ;

    rule_ctx["is_pointwise"] = 0 ;
    rule_ctx["is_unit"] = 1 ;
    rule_ctx["is_apply"] = 0 ;
    rule_ctx["unit"] = unit_ctx ;
  } else if(rule_type == "apply") {
    variable apply_var = *(output.begin()) ;

    DictionaryTemplateValue apply_ctx ;
    apply_ctx["target_name"] = apply_var.str() ;
    apply_ctx["target_vname"] = vnames[apply_var] ;
    apply_ctx["container"] = ctypetable[apply_var] ;
    apply_ctx["container_args"] = cargtable[apply_var] ;
    apply_ctx["operator"] = apply_op.str() ;
    apply_ctx["line_number"] = apply_op_line_no ;
    apply_ctx["is_param"] = paramApply ;
    apply_ctx["is_singleton"] = singletonApply ;

    rule_ctx["is_pointwise"] = 0 ;
    rule_ctx["is_unit"] = 0 ;
    rule_ctx["is_apply"] = 1 ;
    rule_ctx["apply"] = apply_ctx ;
  }

  {
    ostringstream compute_ss ;
    AST_simplePrint printer(compute_ss, -1, prettyOutput) ;
    ap->accept(printer) ;

    DictionaryTemplateValue compute_ctx ;
    compute_ctx["line_number"] = compute_line_no ;
    compute_ctx["spec"] = compute_ss.str() ;
    rule_ctx["compute"] = compute_ctx ;
  }

  DictionaryTemplateValue root_ctx ;
  root_ctx["pln"] = !prettyOutput ;
  root_ctx["debug_info"] = parseInfo.debug_info ;
  root_ctx["rule"] = rule_ctx ;

  string rule_text ;
  try {
    char const * rule_template_name = nullptr ;
    if(rule_type == "pointwise") {
      rule_template_name = "pointwise_rule" ;
    } else if(rule_type == "unit") {
      if(paramUnit) {
        rule_template_name = "param_unit_rule" ;
      }
    } else if(rule_type == "apply") {
      if(paramApply) {
        if(singletonApply) {
          rule_template_name = "singleton_param_apply_rule" ;
        } else {
          rule_template_name = "param_apply_rule" ;
        }
      }
    }
    
    if(rule_template_name == nullptr) {
      ostringstream ss ;
      ss << "unsupported cudarule type " << rule_type ;
      throw std::runtime_error(ss.str()) ;
    }

    rule_text = cuda_templates.render(rule_template_name, root_ctx,
      [this, rule_type_line_no](
        std::ostream & s, char const * partial
      ) {
        if(!prettyOutput) {
          s << std::endl << "#line " << rule_type_line_no << " \""
            << filename << "\"" << std::endl ;
        } else {
          s << std::endl ;
        }
      }) ;
  } catch(std::runtime_error const & error) {
    ostringstream ss ;
    string const message = error.what() ;
    ss << "error rendering rule: " << error.what() ;
    throw parseError(ss.str()) ;
  }
  outputFile << rule_text << endl ;

  if(!use_prelude && sized_outputs && (rule_type != "apply")) 
    throw parseError("need prelude to size output type!") ;
}

void parseRuleLineInfo::clear() {
  rule_type = 0 ;
  signature = 0 ;
  applyop = 0 ;
  constraint = 0 ;
  parametric = 0 ;
  conditional = 0 ;
  specialized = 0 ;
  options.clear() ;
  inplace.clear() ;
  comments.clear() ;
  prelude = 0 ;
  compute = 0 ;
}

void parseRuleInfo::clear() {
  is_gpu = 0 ;
  rule_type.clear() ;
  signature.clear() ;
  applyop.clear() ;
  constraint.clear() ;
  parametric.clear() ;
  conditional.clear() ;
  is_specialized = 0 ;
  options.clear() ;
  inplace.clear() ;
  comments.clear() ;
  use_prelude = 0 ;
  prelude.clear() ;
  use_compute = 1 ;
  template_name.clear() ;
  compute.clear() ;
  lines.clear() ;
}

void parseFile::render_rule(
  ostream & outputFile,
  parseSharedInfo const & parseInfo,
  parseRuleInfo const & ruleInfo
) {
  using namespace Loci ;
  string rule_text ;
  try {    
    if(ruleInfo.template_name.empty()) {
      ostringstream ss ;
      ss << "unsupported rule type " << ruleInfo.rule_type ;
      throw std::runtime_error(ss.str()) ;
    }

    if(ruleInfo.is_gpu) {
      parseRuleLineInfo const & lines = ruleInfo.lines ;
      rule_text = cuda_templates.render(
        ruleInfo.template_name, ruleInfo.ctx,
        [this, lines](
          std::ostream & s, char const * partial
        ) {
          if(!prettyOutput) {
            s << std::endl << "#line " << lines.rule_type << " \""
              << filename << "\"" << std::endl ;
          } else {
             s << std::endl ;
          }
        }
      ) ;
    } else {
      parseRuleLineInfo const & lines = ruleInfo.lines ;
      rule_text = cpu_templates.render(
        ruleInfo.template_name, ruleInfo.ctx,
        [this, lines](
          std::ostream & s, char const * partial
        ) {
          if(!prettyOutput) {
            s << std::endl << "#line " << lines.rule_type << " \""
              << filename << "\"" << std::endl ;
          } else {
            s << std::endl ;
          }
        }
      ) ;
    }
  } catch(std::runtime_error const & error) {
    ostringstream ss ;
    string const message = error.what() ;
    ss << "error rendering rule: " << error.what() ;
    throw parseError(ss.str()) ;
  }
  outputFile << rule_text << endl ;
}

void parseFile::process_and_validate_rule_info(
  string const & comment, string const & docvarname,
  parseSharedInfo const & parseInfo,
  parseRuleInfo & ruleInfo
) {
  using namespace Loci ;

  // if conditional variable is defined, it must be of type param<bool>
  if(!ruleInfo.conditional.empty()) {
    variable cond(ruleInfo.conditional) ;
    auto mi = lookupVarType(cond) ;
    if(!checkTypeValid(mi)) {
      cerr << filename << ':' << ruleInfo.lines.conditional
           << ":0: warning: type of conditional variable '"
           << cond << "' not found!" << endl ;
    } else {
      string val = mi->second.container + mi->second.container_args ;
      val.erase(std::remove_if(val.begin(), val.end(),
        [](unsigned char c) { return isspace(c) ; }), val.end()) ;
      if(val != "param<bool>") {
        ostringstream ss ;
        ss << "conditional variable must be typed as a param<bool> at line "
           << ruleInfo.lines.conditional ;
        throw parseError(ss.str()) ;
      }
    }
  }

  // check for correctness of inplace specification and create inplace pairs
  vector<pair<variable, variable>> inplace_pairs ;
  if(!ruleInfo.inplace.empty()) {
    for(size_t i = 0; i < ruleInfo.inplace.size(); ++i) {
      exprP p = expression::create(ruleInfo.inplace[i]) ;
      exprList l = collect_associative_op(p, OP_OR) ;
      if(l.size() != 2) {
        ostringstream ss ;
        ss << "inplace needs two variables with a '|' separator, at line "
           << ruleInfo.lines.inplace[i] ;
        throw parseError(ss.str()) ;
      }

      auto iter = l.begin() ;
      variable v1(*iter++) ;
      variable v2(*iter++) ;
      inplace_pairs.push_back(std::make_pair(v1, v2)) ;
    }
  }

  // separate signature into head and body expressions
  exprP head = 0, body = 0 ;
  {
    string::size_type pos = ruleInfo.signature.find("<-") ;
    if(pos == string::npos) {
      head = expression::create(ruleInfo.signature) ;
      body = 0 ;
    } else {
      head = expression::create(ruleInfo.signature.substr(0, pos)) ;
      body = expression::create(ruleInfo.signature.substr(pos+2)) ;
      if(ruleInfo.signature.find("<-", pos+2) != string::npos) {
        ostringstream ss ;
        ss << "rule signature must have a single head and a single body separated by '<-', at line "
           << ruleInfo.lines.signature;
        throw parseError(ss.str()) ;
      }
    }
  }

  // default and optional rules cannot have inputs
  if(ruleInfo.rule_type == "optional" || ruleInfo.rule_type == "default") {
    if(body != 0) {
      ostringstream ss ;
      ss << "'optional' or 'default' rules should not have a body "
         << "(defined by '<-' operator)!, at line " << ruleInfo.lines.signature ;
      throw parseError(ss.str()) ;
    }
  }

  // default and optional rules cannot have constraints
  // for other rule types, lack of inputs means the constraints are mandatory
  if(ruleInfo.rule_type == "optional" || ruleInfo.rule_type == "default") {
    if(!ruleInfo.constraint.empty()) {
      ostringstream ss ;
      ss << "'optional' and 'default' rules should not have a constraint, at line "
         << ruleInfo.lines.constraint ;
      throw parseError(ss.str()) ;
    }
  } else {
    if(body == 0 && ruleInfo.constraint.empty()) {
      ostringstream ss ;
      ss << "rules without bodies should have a defined constraint as input, at line "
         << ruleInfo.lines.constraint ;
      throw parseError(ss.str()) ;
    }
  }

  // collect target, source, and constraint variable mappings
  set<vmap_info> sources ;
  set<vmap_info> targets ;
  set<vmap_info> constraints ;
  if(body != 0) {
    fill_descriptors(sources,collect_associative_op(body,OP_COMMA)) ;
  }
  fill_descriptors(targets,collect_associative_op(head,OP_COMMA)) ;
  if(!ruleInfo.constraint.empty()) {
    exprP C = expression::create(ruleInfo.constraint) ;
    fill_descriptors(constraints,collect_associative_op(C, OP_COMMA)) ;
  }

  // extract input, output, and constraint variables
  variableSet input_vars ;
  variableSet output_vars ;
  variableSet constraint_vars ;
  for(auto i = sources.begin(); i != sources.end(); ++i) {
    for(size_t j = 0; j < i->mapping.size(); ++j) {
      input_vars += i->mapping[j] ;
    }
    input_vars += i->var ;
  }
  for(auto i = targets.begin(); i != targets.end(); ++i) {
    for(size_t j = 0; j < i->mapping.size(); ++j) {
      input_vars += i->mapping[j] ;
    }
    output_vars += i->var ;
  }
  for(auto i = constraints.begin(); i != constraints.end(); ++i) {
    for(size_t j = 0; j < i->mapping.size(); ++j) {
      constraint_vars += i->mapping[j] ;
    }
    constraint_vars += i->var ;
  }

  // TODO: is validate_set used anywhere?
  set<list<variable>> validate_set ;
  for(auto i = sources.begin(); i != sources.end(); ++i) {
    if(i->mapping.size() == 0) {
      for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
        list<variable> vbasic ;
        vbasic.push_back(*vi) ;
        validate_set.insert(vbasic) ;
      }
    } else {
      vector<list<variable>> maplist = expand_mapping(i->mapping) ;
      size_t msz = maplist.size() ;
      for(size_t j = 0; j < msz; ++j) {
        list<variable> mapping_list = maplist[j] ;
        validate_set.insert(mapping_list) ;
        for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
          list<variable> mapping_list2 = maplist[j] ;
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
  for(auto i = targets.begin(); i != targets.end(); ++i) {
    if(i->mapping.size() == 0) {
      for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
        list<variable> vbasic ;
        variable vt = *vi ;
        while(vt.get_info().priority.size() != 0) {
          vt = vt.drop_priority() ;
        }
        vbasic.push_back(vt) ;
        validate_set.insert(vbasic) ;
      }
    } else {
      vector<list<variable>> maplist = expand_mapping(i->mapping) ;
      size_t msz = maplist.size() ;
      for(size_t j = 0; j < msz; ++j) {
        list<variable> mapping_list = maplist[j] ;
        validate_set.insert(mapping_list) ;
        for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
          list<variable> mapping_list2 = maplist[j] ;
          variable vt = *vi ;
          while(vt.get_info().priority.size() != 0) {
            vt = vt.drop_priority() ;
          }
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

  // input variables cannot have priority
  for(auto vi = input_vars.begin(); vi != input_vars.end(); ++vi) {
    if(vi->get_info().priority.size() != 0) {
      ostringstream ss ;
      ss << "improper use of priority annotation on rule input, var=" << *vi ;
      throw parseError(ss.str()) ;
    }
  }

  // only pointwise and default rules can have priority over output variables
  if(ruleInfo.rule_type != "pointwise" && ruleInfo.rule_type != "default") {
    for(auto vi = output_vars.begin(); vi != output_vars.end(); ++vi) {
      if(vi->get_info().priority.size() != 0) {
        ostringstream ss ;
        ss << "only pointwise and default rules can use priority annotation, var=" << *vi ;
        throw parseError(ss.str()) ;
      }
    }
  }

  // catch undelcared input, output, and constraint variables
  for(auto vi = input_vars.begin(); vi != input_vars.end(); ++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      ostringstream ss ;
      ss << "unable to determine type of Loci variable '" << *vi << "'" ;
      throw parseError(ss.str()) ;
    }
  }
  for(auto vi = output_vars.begin(); vi != output_vars.end(); ++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      ostringstream ss ;
      ss << "unable to determine type of Loci variable '" << *vi << "'" ;
      throw parseError(ss.str()) ;
    }
  }
  for(auto vi = constraint_vars.begin(); vi != constraint_vars.end(); ++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      ostringstream ss ;
      ss << "unable to determine type of Loci variable '" << *vi << "'" ;
      throw parseError(ss.str()) ;
    }
  }

  // check if variables paired by inplace specification are specified as either
  // input or output to the rule
  variableSet checkset ;
  for(auto vi = input_vars.begin(); vi != input_vars.end(); ++vi) {
    variable v = *vi ;
    while(v.get_info().priority.size() != 0) {
      v = v.drop_priority() ;
    }
    checkset += v ;
  }
  for(auto vi = output_vars.begin(); vi != output_vars.end(); ++vi) {
    variable v = *vi ;
    while(v.get_info().priority.size() != 0) {
      v = v.drop_priority() ;
    }
    checkset += v ;
  }
  for(auto ipi = inplace_pairs.begin(); ipi != inplace_pairs.end(); ++ipi) {
    variable v = ipi->first ;
    while(v.get_info().priority.size() != 0) {
      v = v.drop_priority() ;
    }
    if(!checkset.inSet(v)) {
      ostringstream ss ;
      ss << "inplace variable '" << ipi->first << "' not input or output variable" ;
      throw parseError(ss.str()) ;
    }
    v = ipi->second ;
    while(v.get_info().priority.size() != 0) {
      v = v.drop_priority() ;
    }
    if(!checkset.inSet(v)) {
      ostringstream ss ;
      ss << "inplace variable '" << ipi->first << "' not input or output variable" ;
      throw parseError(ss.str()) ;
    }
  }

  // pointwise rules cannot compute param
  if(ruleInfo.rule_type == "pointwise") {
    for(auto vi = output_vars.begin(); vi != output_vars.end(); ++vi) {
      auto mi = lookupVarType(*vi) ;
      if(mi->second.container == "param" && vi->get_info().name != "OUTPUT") {
        throw parseError("pointwise rule cannot compute param, use singleton") ;
      }
    }
  }

  // singleton rules cannot compute stores
  if(ruleInfo.rule_type == "singleton") {
    for(auto vi = output_vars.begin(); vi != output_vars.end(); ++vi) {
      auto mi = lookupVarType(*vi) ;
      if(mi->second.container == "store" ||
         mi->second.container == "storeVec" ||
         mi->second.container == "multiStore") {
        throw parseError("singleton rule cannot compute stores, use pointwise") ;
      }
    }
  }

  bool singletonApply = false ;
  if(ruleInfo.rule_type == "apply") {
    if(output_vars.size() != 1) {
      throw parseError("apply rule should have only one output variable") ;
    }
    variable av = *(output_vars.begin()) ;
    typedoc tinfo = lookupVarType(av)->second ;
    if(tinfo.container == "param") {
      bool allparam = true ;
      for(auto vi = input_vars.begin(); vi != input_vars.end(); ++vi) {
        typedoc tinfo2 = lookupVarType(*vi)->second ;
        if(tinfo2.container != "param") {
          allparam = false ;
        }
      }
      if(allparam) {
        singletonApply = true ;
        cerr << "NOTE: parameter-only apply rule on '"
             << av << "' now executes single instance" << endl ;
      }
    }
  }

  variableSet output_stores = output_vars ;
  for(auto ipi = inplace_pairs.begin(); ipi != inplace_pairs.end(); ++ipi) {
    output_stores -= ipi->first ;
    output_stores += ipi->second ;
  }

  variableSet input_stores = input_vars ;
  input_stores -= output_stores ;

  variableSet named_stores = input_stores ;
  named_stores += output_stores ;
  for(auto ipi = inplace_pairs.begin(); ipi != inplace_pairs.end(); ++ipi) {
    named_stores -= ipi->first ;
  }

  bool paramOutput = false ;
  for(auto vi = output_stores.begin(); vi != output_stores.end(); ++vi) {
    auto mi = lookupVarType(*vi) ;
    if(vi->get_info().name != "OUTPUT" && mi->second.container == "param") {
      paramOutput = true ;
    }
  }

  // C++ names of variables
  map<variable, string> vnames ;

  // type returned by store_instance's operator[]
  map<variable, string> typetable ;

  // container type
  map<variable, string> ctypetable ;

  // container arguments
  map<variable, string> cargtable ;

  // reduction operator template argument type
  map<variable, string> rargtable ;

  for(auto vi = input_stores.begin(); vi != input_stores.end(); ++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      ostringstream ss ;
      ss << "unknown type for Loci variable '" << *vi << "'" ;
      throw parseError(ss.str()) ;
    }

    if(mi->second.container != "param" &&
       mi->second.container != "Map" &&
       mi->second.container != "MapVec" &&
       mi->second.container != "multiMap" &&
       mi->second.container != "store" &&
       mi->second.container != "storeVec" &&
       mi->second.container != "storeMat" &&
       mi->second.container != "multiStore" &&
       mi->second.container != "Constraint" &&
       mi->second.container != "constraint" &&
       mi->second.container != "blackbox") {
       ostringstream ss ;
       ss << "unsupported container '" << mi->second.container
          << "' for Loci variable '" << *vi << "'" ;
       throw parseError(ss.str()) ;
    }

    vnames[*vi] = var2name(*vi) ;

    ctypetable[*vi] = mi->second.container ;

    if(mi->second.container == "Map") {
      typetable[*vi] = "int const *" ;
      cargtable[*vi] = "" ;
    } else if(mi->second.container == "multiMap") {
      if(ruleInfo.is_gpu) {
        typetable[*vi] = "Loci::constMultiAccessor<Loci::Entity> " ;
      } else {
        typetable[*vi] = "Loci::multiMap::arrayHelper_const " ;
      }
      cargtable[*vi] = "" ;
      rargtable[*vi] = "" ;
    } else if(mi->second.container == "MapVec") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = string("Array<Entity,") + arg + "> const *" ;
        cargtable[*vi] = arg ;
        rargtable[*vi] = "" ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    } else if(mi->second.container == "multiStore") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        if(ruleInfo.is_gpu) {
          typetable[*vi] = "Loci::constMultiAccessor<" + arg + "> ";
        } else {
          typetable[*vi] = "Loci::const_Vect<" + arg + "> " ;
        }
        cargtable[*vi] = arg ;
        rargtable[*vi] = "Loci::Vect<" + arg + "> " ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    } else if(mi->second.container == "storeVec") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        if(ruleInfo.is_gpu) {
          ostringstream ss ;
          ss << "storeVec is not supported on gpu rules" ;
          throw parseError(ss.str()) ;
        } else {
          typetable[*vi] = "Loci::const_Vect<" + arg + "> " ;
        }
        cargtable[*vi] = arg ;
        rargtable[*vi] = "Loci::Vect<" + arg + "> " ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    } else if(mi->second.container == "storeMat") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        if(ruleInfo.is_gpu) {
          ostringstream ss ;
          ss << "storeMat is not supported on gpu rules" ;
          throw parseError(ss.str()) ;
        } else {
          typetable[*vi] = "Loci::const_Mat<" + arg + "> " ;
        }
        cargtable[*vi] = arg ;
        rargtable[*vi] = "Loci::Mat<" + arg + "> " ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    } else if(mi->second.container == "blackbox") {
      if(ruleInfo.is_gpu) {
        ostringstream ss ;
        ss << "unsupported container '" << mi->second.container
           << "' for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }

      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = arg + " const *" ;
        cargtable[*vi] = arg ;
        rargtable[*vi] = arg ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    } else if(mi->second.container == "Constraint" || mi->second.container == "constraint") {
      if(ruleInfo.is_gpu) {
        ostringstream ss ;
        ss << "unsupported container '" << mi->second.container
           << "' for Loci variable '" << *vi << "' in gpu rule" ;
        throw parseError(ss.str()) ;
      }

      typetable[*vi] = "entitySet const &" ;
      cargtable[*vi] = "" ;
      rargtable[*vi] = "" ;
    } else {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = arg + " const *" ;
        cargtable[*vi] = arg ;
        rargtable[*vi] = arg ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    }

    if(vi->get_info().priority.size() != 0) {
      variable vnp = *vi ;
      while(vnp.get_info().priority.size() != 0) {
        vnp = vnp.drop_priority() ;
      }

      vnames[vnp] = vnames[*vi] ;
      ctypetable[vnp] = ctypetable[*vi] ;
      typetable[vnp] = typetable[*vi] ;
      cargtable[vnp] = cargtable[*vi] ;
      rargtable[vnp] = rargtable[*vi] ;
    }
  }

  for(auto vi = output_stores.begin(); vi != output_stores.end(); ++vi) {
    auto mi = lookupVarType(*vi) ;
    if(!checkTypeValid(mi)) {
      ostringstream ss ;
      ss << "unknown type for Loci variable '" << *vi << "'" ;
      throw parseError(ss.str()) ;
    }

    if(mi->second.container != "param" &&
       mi->second.container != "Map" &&
       mi->second.container != "MapVec" &&
       mi->second.container != "multiMap" &&
       mi->second.container != "store" &&
       mi->second.container != "storeVec" &&
       mi->second.container != "storeMat" &&
       mi->second.container != "multiStore" &&
       mi->second.container != "Constraint" &&
       mi->second.container != "constraint" &&
       mi->second.container != "blackbox") {
       ostringstream ss ;
       ss << "unsupported container '" << mi->second.container
          << "' for Loci variable '" << *vi << "'" ;
       throw parseError(ss.str()) ;
    }

    vnames[*vi] = var2name(*vi) ;

    ctypetable[*vi] = mi->second.container ;

    if(mi->second.container == "Map") {
      typetable[*vi] = "int * " ;
      cargtable[*vi] = "" ;
    } else if(mi->second.container == "multiMap") {
      if(ruleInfo.is_gpu) {
        typetable[*vi] = "Loci::MultiAccessor<Loci::Entity> " ;
      } else {
        typetable[*vi] = "Loci::multiMap::arrayHelper " ;
      }
      cargtable[*vi] = "" ;
      rargtable[*vi] = "" ;
    } else if(mi->second.container == "MapVec") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = string("Array<Entity,") + arg + "> * " ;
        cargtable[*vi] = arg ;
        rargtable[*vi] = "" ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    } else if(mi->second.container == "multiStore") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        if(ruleInfo.is_gpu) {
          typetable[*vi] = "Loci::MultiAccessor<" + arg + "> ";
        } else {
          typetable[*vi] = "Loci::Vect<" + arg + "> " ;
        }
        cargtable[*vi] = arg ;
        rargtable[*vi] = "Loci::Vect<" + arg + "> " ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    } else if(mi->second.container == "storeVec") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        if(ruleInfo.is_gpu) {
          ostringstream ss ;
          ss << "Loci variable of type storeVec is not supported on gpu rules" ;
          throw parseError(ss.str()) ;
        } else {
          typetable[*vi] = "Loci::Vect<" + arg + "> " ;
        }
        cargtable[*vi] = arg ;
        rargtable[*vi] = "Loci::Vect<" + arg + "> " ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    } else if(mi->second.container == "storeMat") {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        if(ruleInfo.is_gpu) {
          ostringstream ss ;
          ss << "Loci variable of type storeMat is not supported on gpu rules" ;
          throw parseError(ss.str()) ;
        } else {
          typetable[*vi] = "Loci::Mat<" + arg + "> " ;
        }
        cargtable[*vi] = arg ;
        rargtable[*vi] = "Loci::Mat<" + arg + "> " ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    } else if(mi->second.container == "blackbox") {
      if(ruleInfo.is_gpu) {
        ostringstream ss ;
        ss << "unsupported container '" << mi->second.container
           << "' for Loci variable '" << *vi << "' in gpu rule" ;
        throw parseError(ss.str()) ;
      }

      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = arg + " *" ;
        cargtable[*vi] = arg ;
        rargtable[*vi] = arg ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    } else if(mi->second.container == "Constraint" || mi->second.container == "constraint") {
      if(ruleInfo.is_gpu) {
        ostringstream ss ;
        ss << "unsupported container '" << mi->second.container
           << "' for Loci variable '" << *vi << "' in gpu rule" ;
        throw parseError(ss.str()) ;
      }

      typetable[*vi] = "entitySet &" ;
      cargtable[*vi] = "" ;
      rargtable[*vi] = "" ;
    } else {
      string scratch = mi->second.container_args ;
      string::size_type start = scratch.find('<') ;
      string::size_type end = scratch.rfind('>') ;
      if(start != string::npos && end != string::npos && start+1 < end) {
        string arg = scratch.substr(start+1, end-start-1) ;
        typetable[*vi] = arg + " *";
        cargtable[*vi] = arg ;
        rargtable[*vi] = arg ;
      } else {
        ostringstream ss ;
        ss << "unexpected type for Loci variable '" << *vi << "'" ;
        throw parseError(ss.str()) ;
      }
    }

    if(vi->get_info().priority.size() != 0) {
      variable vnp = *vi ;
      while(vnp.get_info().priority.size() != 0) {
        vnp = vnp.drop_priority() ;
      }

      vnames[vnp] = vnames[*vi] ;
      ctypetable[vnp] = ctypetable[*vi] ;
      typetable[vnp] = typetable[*vi] ;
      cargtable[vnp] = cargtable[*vi] ;
      rargtable[vnp] = rargtable[*vi] ;
    }
  }

  // for inplace variables, the target variable type information is same as
  // the source variable ;
  for(auto ipi = inplace_pairs.begin(); ipi != inplace_pairs.end(); ++ipi) {
    vnames[ipi->first] = vnames[ipi->second] ;
    typetable[ipi->first] = typetable[ipi->second] ;
    ctypetable[ipi->first] = ctypetable[ipi->second] ;
    cargtable[ipi->first] = cargtable[ipi->second] ;
    rargtable[ipi->first] = rargtable[ipi->second] ;
  }

  // calculate name of the rule class
  string class_name = "file_" ;
  for(auto c : filename) {
    if(std::isalpha(c) || std::isdigit(c) || c == '_')
      class_name += c ;
    if(c == '.') break ;
  }
  class_name += '0' + (cnt/100)%10 ;
  class_name += '0' + (cnt/10)%10 ;
  class_name += '0' + (cnt)%10 ;
  timespec tdata ;
  clock_gettime(CLOCK_MONOTONIC, &tdata) ;
  ostringstream tss ;
  tss << '_' << tdata.tv_sec << 'm' << tdata.tv_nsec/1000000 ;
  class_name += tss.str() ;

  // calculate debug name of the rule
  string rule_debug_name ;
  if(parseInfo.debug_info > 0) {
    ostringstream oss ;
    for(auto i = targets.begin(); i != targets.end();) {
      for(size_t j = 0; j < i->mapping.size(); ++j)
        oss << i->mapping[j] << "->" ;
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
    for(auto i = sources.begin(); i != sources.end();) {
      for(size_t j = 0; j < i->mapping.size(); ++j)
        oss << i->mapping[j] << "->" ;
      if(i->var.size() > 1)
        oss << "(" ;
      for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
        if(vi != i->var.begin())
          oss << "," ;
        oss << *vi ;
      }
      if(i->var.size() > 1)
        oss << ")" ;
      ++i ;
      if(i != sources.end())
        oss << "," ;
    }
    if(!ruleInfo.constraint.empty())
      oss << ",constraint(" << ruleInfo.constraint << ")" ;

    rule_debug_name = oss.str() ;
  }

  // add comment if rule does not have comments specified
  if(ruleInfo.comments.empty()) {
    if(comment.empty()) {
      auto mi = lookupVarType(*output_vars.begin()) ;
      ruleInfo.comments.push_back(cleanupCommentQuoted(mi->second.comment)) ;
      ruleInfo.lines.comments.push_back(ruleInfo.lines.rule_type) ;
    } else {
      ruleInfo.comments.push_back(cleanupCommentQuoted(comment)) ;
      ruleInfo.lines.comments.push_back(ruleInfo.lines.rule_type) ;
    }
  }

  // determine if outputs of the rule are variable-sized stores
  bool sizedOutputs = false ;
  variableSet sized_output_stores = output_stores ;
  sized_output_stores -= input_vars ;
  for(auto vi = sized_output_stores.begin(); vi != sized_output_stores.end(); ++vi) {
    auto mi = lookupVarType(*vi) ;
    string const & ct = mi->second.container ;
    if(ct == "storeVec" || ct == "storeMat" || ct == "multiStore") {
      sizedOutputs = true ;
    }
  }

  // gpu rules cannot have sized outputs
  if(ruleInfo.is_gpu && sizedOutputs) {
    ostringstream ss ;
    ss << "gpu rules cannot have variable sized outputs" ;
    throw parseError(ss.str()) ;
  }

  if(!ruleInfo.use_prelude && sizedOutputs && ruleInfo.rule_type != "apply") {
    ostringstream ss ;
    ss << "need prelude to size outputs of rule at line "
       << ruleInfo.lines.rule_type ;
    throw parseError(ss.str()) ;
  }

  if(ruleInfo.rule_type == "singleton" ||
     ruleInfo.rule_type == "optional" ||
     ruleInfo.rule_type == "default" ||
     ruleInfo.rule_type == "constraint" ||
     (paramOutput && ruleInfo.rule_type != "apply")) {
    if(ruleInfo.use_prelude) {
      ostringstream ss ;
      ss << "inappropriate prelude on " << ruleInfo.rule_type << " rule" ;
      throw parseError(ss.str()) ;
    }
  }

  // process prelude block
  string prelude_body ;
  if(ruleInfo.use_prelude) {
    varmap typemap ;
    if(!ruleInfo.is_gpu) {
      typemap["cerr"] = localIdentifier() ;
      typemap["std::cerr"] = localIdentifier() ;
      typemap["cout"] = localIdentifier() ;
      typemap["std::cout"] = localIdentifier() ;
      typemap["debugout"] = localIdentifier() ;
      typemap["Loci::debugout"] = localIdentifier() ;
      typemap["EMPTY"] = localIdentifier() ;
    }

    istringstream is(ruleInfo.prelude) ;
    int prelude_line_no = ruleInfo.lines.prelude ;

    AST_type::ASTP prelude_ast = parseBlockRaw(
      is, prelude_line_no, filename, typemap
    ) ;
    //cerr << endl << "prelude AST" << endl ;
    //AST_printObjectTree treeout(cerr) ;
    //prelude_ast->accept(treeout) ;

    AST_editPrelude edit_prelude(filename, vnames) ;
    prelude_ast->accept(edit_prelude) ;

    // Now remove and save the open and close braces in the parseBlockRaw
    CPTR<AST_BlockRaw> bigblock = CPTR<AST_BlockRaw>(prelude_ast) ;
    CPTR<AST_type> open = bigblock->elements[0] ;
    int bsz = bigblock->elements.size() ;
    CPTR<AST_type> close = bigblock->elements[bsz-1] ;
    for(int i = 0; i < bsz-1; ++i) {
      bigblock->elements[i] = bigblock->elements[i+1] ;
    }
    bigblock->elements.pop_back() ;
    bigblock->elements.pop_back() ;

    ostringstream prelude_oss ;
    AST_simplePrint printer(prelude_oss, -1, prettyOutput) ;
    prelude_ast->accept(printer) ;

    prelude_body = prelude_oss.str() ;
  }

  // process compute block
  string compute_body ;
  if(ruleInfo.use_compute) {
    varmap typemap ;
    if(!ruleInfo.is_gpu) {
      typemap["cerr"] = localIdentifier() ;
      typemap["std::cerr"] = localIdentifier() ;
      typemap["cout"] = localIdentifier() ;
      typemap["std::cout"] = localIdentifier() ;
      typemap["debugout"] = localIdentifier() ;
      typemap["Loci::debugout"] = localIdentifier() ;
      typemap["EMPTY"] = localIdentifier() ;
    }

    int compute_line_no = ruleInfo.lines.compute ;
    istringstream compute_ss(ruleInfo.compute) ;

    CPTR<AST_type> compute_ast = parseBlock(
      compute_ss, compute_line_no, filename, typemap
    ) ;
    //cerr << endl << "AST" << endl ;
    //AST_printObjectTree treeout(cerr) ;
    //compute_ast->accept(treeout) ;

    AST_condenseLeftAssociative condenseOps ;
    compute_ast->accept(condenseOps) ;

    // This is sort of a hack because the precedence of the mapping
    // operator (->) is context senstive
    AST_editLociMapArrayAccess mapEditOps ;
    compute_ast->accept(mapEditOps) ;

    if(parseInfo.diag_level > 0) {
      AST_printTree diagout(cerr) ;
      compute_ast->accept(diagout) ;
    }

    AST_errorCheck syntaxChecker ;
    compute_ast->accept(syntaxChecker) ;
    if(syntaxChecker.hasErrors()) {
#ifdef VERBOSE
      AST_simplePrint printer(cerr, -1, false) ;
      compute_ast->accept(printer) ;
#endif
      throw parseError("syntax error in compute body") ;
    }

    if(ruleInfo.is_gpu) {
      AST_collectAccessInfo varaccess ;
      compute_ast->accept(varaccess) ;
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

      for(auto i = varaccess.accessed.begin(); i != varaccess.accessed.end(); ++i) {
        readvars += i->var ;
        for(size_t j = 0; j < i->mapping.size(); ++j) {
          readvars += i->mapping[j] ;
        }
      }
      for(auto i = varaccess.writes.begin(); i != varaccess.writes.end(); ++i) {
        writevars += i->var ;
        for(size_t j = 0; j < i->mapping.size(); ++j) {
          readvars += i->mapping[j] ;
        }
      }

      readvars -= writevars ;

      // Now remove and save the open and close braces in the parseBlock
      CPTR<AST_Block> bigblock = CPTR<AST_Block>(compute_ast) ;
      CPTR<AST_type> open = bigblock->elements[0] ;
      int bsz = bigblock->elements.size() ;
      CPTR<AST_type> close = bigblock->elements[bsz-1] ;
      for(int i = 0; i < bsz-1; ++i) {
        bigblock->elements[i] = bigblock->elements[i+1] ;
      }
      bigblock->elements.pop_back() ;
      bigblock->elements.pop_back() ;

      if(ruleInfo.rule_type == "apply") {
        std::stringstream ss ;

        ss << "loci_reduction_t _reduce_op_ ;" ;
        compute_line_no = ruleInfo.lines.compute ;
        AST_type::ASTP reduce_op_decl = parseDeclaration(
          ss, compute_line_no, filename, typemap
        ) ;

        ss.str("") ;
        ss.clear() ;
        ss << "value_t _local_ = _reduce_op_.identity() ;" ;
        compute_line_no = ruleInfo.lines.compute ;
        AST_type::ASTP local_value_decl = parseDeclaration(
          ss, compute_line_no, filename, typemap
        ) ;

        ss.str("") ;
        ss.clear() ;
        ss << "return _local_ ;" ;
        compute_line_no = ruleInfo.lines.compute ;
        AST_type::ASTP return_local_value = parseSpecialControlStatement(
          ss, compute_line_no, filename, typemap
        ) ;

        bigblock->elements.insert(bigblock->elements.begin(), local_value_decl) ;
        bigblock->elements.insert(bigblock->elements.begin(), reduce_op_decl) ;
        bigblock->elements.insert(bigblock->elements.end(), return_local_value) ;

        AST_editJoin edit_join ;
        compute_ast->accept(edit_join) ;
      }

      AST_editGPULociVariableAccess AST_vareditor(vnames, ctypetable) ;
      compute_ast->accept(AST_vareditor) ;

      //AST_printObjectTree treeout(cerr) ;
      //compute_ast->accept(treeout) ;
    } else {
      AST_editLociDirective AST_directive(filename) ;
      compute_ast->accept(AST_directive) ;

      AST_editLociVariableAccess AST_editor(vnames, ctypetable) ;
      compute_ast->accept(AST_editor) ;

      CPTR<AST_Block> bigblock = CPTR<AST_Block>(compute_ast) ;
      CPTR<AST_type> open = bigblock->elements[0] ;
      int bsz = bigblock->elements.size() ;
      CPTR<AST_type> close = bigblock->elements[bsz-1] ;
      for(int i = 0; i < bsz-1; ++i) {
        bigblock->elements[i] = bigblock->elements[i+1] ;
      }
      bigblock->elements.pop_back() ;
      bigblock->elements.pop_back() ;
    }

    ostringstream compute_oss ;
    AST_simplePrint printer(compute_oss, -1, prettyOutput) ;
    compute_ast->accept(printer) ;

    compute_body = compute_oss.str() ;
  }

  // set name of the rule template
  if(ruleInfo.rule_type == "pointwise") {
    ruleInfo.template_name = "pointwise_rule" ;
  } else if(ruleInfo.rule_type == "unit") {
    if(paramOutput) {
      ruleInfo.template_name = "param_unit_rule" ;
    } else {
      ruleInfo.template_name = "unit_rule" ;
    }
  } else if(ruleInfo.rule_type == "apply") {
    if(paramOutput) {
      if(singletonApply) {
        ruleInfo.template_name = "singleton_param_apply_rule" ;
      } else {
        ruleInfo.template_name = "param_apply_rule" ;
      }
    } else {
      ruleInfo.template_name = "apply_rule" ;
    }
  } else if(ruleInfo.rule_type == "singleton") {
    ruleInfo.template_name = "singleton_rule" ;
  } else if(ruleInfo.rule_type == "optional") {
    ruleInfo.template_name = "optional_rule" ;
  } else if(ruleInfo.rule_type == "default") {
    ruleInfo.template_name = "default_rule" ;
  } else if(ruleInfo.rule_type == "constraint") {
    ruleInfo.template_name = "constraint_rule" ;
  } else if(ruleInfo.rule_type == "blackbox") {
    ruleInfo.template_name = "blackbox_rule" ;
  }

  if(ruleInfo.template_name.empty()) {
    ostringstream ss ;
    ss << "no template set for rule: " << rule_debug_name ;
    throw parseError(ss.str()) ;
  }

  DictionaryTemplateValue rule_ctx ;

  rule_ctx["type"] = ruleInfo.rule_type ;
  rule_ctx["parent_class"] = ruleInfo.rule_type + "_rule" ;
  rule_ctx["class"] = class_name ;
  rule_ctx["file"] = filename ;
  rule_ctx["docvar"] = docvarname ;
  rule_ctx["line_number"] = ruleInfo.lines.rule_type ;
  rule_ctx["debug_name"] = rule_debug_name ;

  {
    DictionaryTemplateValue signature_ctx ;
    signature_ctx["line_number"] = ruleInfo.lines.signature ;
    rule_ctx["signature"] = signature_ctx ;
  }

  {
    ArrayTemplateValue input_stores_ctx ;
    for(auto vi = input_stores.begin(); vi != input_stores.end(); ++vi) {
      DictionaryTemplateValue ctx ;
      ctx["name"] = (*vi).str() ;
      ctx["vname"] = vnames[*vi] ;
      ctx["ctype"] = ctypetable[*vi] ;
      ctx["vtype"] = typetable[*vi] ;
      ctx["carg"] = cargtable[*vi] ;
      ctx["rarg"] = rargtable[*vi] ;
      input_stores_ctx.append(ctx) ;
    }
    rule_ctx["input_stores"] = input_stores_ctx ;
  }

  {
    ArrayTemplateValue output_stores_ctx ;
    for(auto vi = output_stores.begin(); vi != output_stores.end(); ++vi) {
      DictionaryTemplateValue ctx ;
      ctx["name"] = (*vi).str() ;
      ctx["vname"] = vnames[*vi] ;
      ctx["ctype"] = ctypetable[*vi] ;
      ctx["vtype"] = typetable[*vi] ;
      ctx["carg"] = cargtable[*vi] ;
      ctx["rarg"] = rargtable[*vi] ;
      output_stores_ctx.append(ctx) ;
    }
    rule_ctx["output_stores"] = output_stores_ctx ;

  }

  {
    ArrayTemplateValue name_stores_ctx ;
    for(auto vi = named_stores.begin(); vi != named_stores.end(); ++vi) {
      DictionaryTemplateValue ctx ;
      ctx["name"] = (*vi).str() ;
      ctx["vname"] = vnames[*vi] ;

      auto mi = access_map.find(lookupVarType(*vi)->second.getFileLoc()) ;
      if(mi != access_map.end()) {
        ctx["has_info_id"] = 1 ;
        ctx["info_id"] = mi->second ;
      } else {
        ctx["has_info_id"] = 0 ;
      }

      name_stores_ctx.append(ctx) ;
    }
    rule_ctx["name_stores"] = name_stores_ctx ;
  }

  {
    ArrayTemplateValue inputs_ctx ;
    for(auto i = sources.begin(); i != sources.end(); ++i) {
      ostringstream ss ;
      for(size_t j = 0; j < i->mapping.size(); ++j) {
        ss << i->mapping[j] << "->" ;
      }

      if(i->var.size() > 1) {
        ss << '(' ;
      }
      for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
        if(vi != i->var.begin()) {
          ss << ',' ;
        }
        ss << *vi ;
      }
      if(i->var.size() > 1) {
        ss << ')' ;
      }
      DictionaryTemplateValue ctx ;
      ctx["str"] = ss.str() ;
      inputs_ctx.append(ctx) ;
    }
    rule_ctx["inputs"] = inputs_ctx ;
  }

  {
    ArrayTemplateValue outputs_ctx ;
    for(auto i = targets.begin(); i != targets.end(); ++i) {
      ostringstream ss ;
      for(size_t j = 0; j < i->mapping.size(); ++j) {
        ss << i->mapping[j] << "->" ;
      }

      if(i->var.size() > 1) {
        ss << '(' ;
      }
      for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
        if(vi != i->var.begin()) {
          ss << ',' ;
        }

        auto ipi = inplace_pairs.begin() ;
        while(ipi != inplace_pairs.end()) {
          if(ipi->first == *vi) break ;
          ++ipi ;
        }
        if(ipi != inplace_pairs.end()) {
          if(i->mapping.size() == 0 || i->var.size() > 1) {
            ss << ipi->first << '=' << ipi->second ;
          } else {
            ss << '(' << ipi->first << '=' << ipi->second << ')' ;
          }
        } else {
          ss << *vi ;
        }
      }
      if(i->var.size() > 1) {
        ss << ')' ;
      }

      DictionaryTemplateValue ctx ;
      ctx["str"] = ss.str() ;
      outputs_ctx.append(ctx) ;
    }

    rule_ctx["outputs"] = outputs_ctx ;
  }

  {
    ArrayTemplateValue constraint_spec_ctx ;
    for(auto i = constraints.begin(); i != constraints.end(); ++i) {
      ostringstream ss ;
      for(size_t j = 0; j < i->mapping.size(); ++j) {
        ss << i->mapping[j] << "->" ;
      }

      if(i->var.size() > 1) {
        ss << '(' ;
      }
      for(auto vi = i->var.begin(); vi != i->var.end(); ++vi) {
        if(vi != i->var.begin()) {
          ss << "," ;
        }
        ss << *vi ;
      }
      if(i->var.size() > 1) {
        ss << ')' ;
      }

      DictionaryTemplateValue ctx ;
      ctx["str"] = ss.str() ;
      constraint_spec_ctx.append(ctx) ;
    }

    DictionaryTemplateValue constraints_ctx ;
    constraints_ctx["spec"] = constraint_spec_ctx ;
    constraints_ctx["line_number"] = ruleInfo.lines.constraint ;
    rule_ctx["constraints"] = constraints_ctx ;
  }

  if(ruleInfo.is_gpu) {
    rule_ctx["option_disable_threading"] = 1 ;
  } else {
    auto opt = find(ruleInfo.options.begin(),
      ruleInfo.options.end(), "disable_threading") ;
    if(opt == ruleInfo.options.end()) {
      rule_ctx["option_disable_threading"] = 0 ;
    } else {
      rule_ctx["option_disable_threading"] = 1 ;
    }
  }

  {
    DictionaryTemplateValue parametric_ctx ;
    if(!ruleInfo.parametric.empty()) {
      rule_ctx["is_parametric"] = 1 ;
      parametric_ctx["spec"] = ruleInfo.parametric ;
      parametric_ctx["line_number"] = ruleInfo.lines.parametric ;
    } else {
      rule_ctx["is_parametric"] = 0 ;
    }
    rule_ctx["parametric"] = parametric_ctx ;
  }

  {
    DictionaryTemplateValue specialized_ctx ;
    if(ruleInfo.is_specialized) {
      specialized_ctx["line_number"] = ruleInfo.lines.specialized ;
    }
    rule_ctx["is_specialized"] = ruleInfo.is_specialized ;
    rule_ctx["specialized"] = specialized_ctx ;
  }

  {
    DictionaryTemplateValue conditional_ctx ;

    if(!ruleInfo.conditional.empty()) {
      rule_ctx["is_conditional"] = 1 ;
      conditional_ctx["spec"] = ruleInfo.conditional ;
      conditional_ctx["line_number"] = ruleInfo.lines.conditional ;
    } else {
      rule_ctx["is_conditional"] = 0 ;
    }

    rule_ctx["conditional"] = conditional_ctx ;
  }

  {
    ArrayTemplateValue comments_ctx ;

    size_t size = ruleInfo.comments.size() ;
    for(size_t i = 0; i < size; ++i) {
      DictionaryTemplateValue ctx ;
      ctx["str"] = ruleInfo.comments[i] ;
      ctx["line_number"] = ruleInfo.lines.comments[i] ;
      comments_ctx.append(ctx) ;
    }

    rule_ctx["comments"] = comments_ctx ;
  }

  DictionaryTemplateValue prelude_ctx ;
  prelude_ctx["line_number"] = ruleInfo.lines.prelude ;
  prelude_ctx["is_specialized"] = ruleInfo.use_prelude ;
  prelude_ctx["spec"] = prelude_body ;
  rule_ctx["prelude"] = prelude_ctx ;

  DictionaryTemplateValue compute_ctx ;
  compute_ctx["line_number"] = ruleInfo.lines.compute ;
  compute_ctx["is_specialized"] = ruleInfo.use_compute ;
  compute_ctx["spec"] = compute_body ;
  rule_ctx["compute"] = compute_ctx ;

  if(ruleInfo.rule_type == "pointwise") {
    rule_ctx["is_pointwise"] = 1 ;
    rule_ctx["is_unit"] = 0 ;
    rule_ctx["is_apply"] = 0 ;
  } else if(ruleInfo.rule_type == "unit") {
    variable unit_var = *(output_stores.begin()) ;

    DictionaryTemplateValue unit_ctx ;
    unit_ctx["target_name"] = unit_var.str() ;
    unit_ctx["target_vname"] = vnames[unit_var] ;
    unit_ctx["container"] = ctypetable[unit_var] ;
    unit_ctx["container_args"] = cargtable[unit_var] ;
    unit_ctx["reduction_args"] = rargtable[unit_var] ;
    unit_ctx["is_param"] = paramOutput ;

    rule_ctx["is_pointwise"] = 0 ;
    rule_ctx["is_unit"] = 1 ;
    rule_ctx["is_apply"] = 0 ;
    rule_ctx["unit"] = unit_ctx ;
  } else if(ruleInfo.rule_type == "apply") {
    variable apply_var = *(output_stores.begin()) ;

    DictionaryTemplateValue apply_ctx ;
    apply_ctx["target_name"] = apply_var.str() ;
    apply_ctx["target_vname"] = vnames[apply_var] ;
    apply_ctx["container"] = ctypetable[apply_var] ;
    apply_ctx["container_args"] = cargtable[apply_var] ;
    apply_ctx["reduction_args"] = rargtable[apply_var] ;
    apply_ctx["operator"] = ruleInfo.applyop ;
    apply_ctx["line_number"] = ruleInfo.lines.applyop ;
    apply_ctx["is_param"] = paramOutput ;
    apply_ctx["is_singleton"] = singletonApply ;

    rule_ctx["is_pointwise"] = 0 ;
    rule_ctx["is_unit"] = 0 ;
    rule_ctx["is_apply"] = 1 ;
    rule_ctx["apply"] = apply_ctx ;
  }

  ruleInfo.ctx["pln"] = !prettyOutput ;
  ruleInfo.ctx["debug_info"] = parseInfo.debug_info ;
  ruleInfo.ctx["rule"] = rule_ctx ;
}

void parseFile::parse_rule_info(
  bool is_gpu, parseSharedInfo const & parseInfo, parseRuleInfo & ruleInfo
) {
  ++cnt ;

  ruleInfo.clear() ;

  ruleInfo.is_gpu = is_gpu ;

  killsp() ;

  if(is_name(is)) {
    ruleInfo.lines.rule_type = line_no ;
    ruleInfo.rule_type = get_name(is) ;
  } else {
    throw parseError("syntax error") ;
  }

  nestedparenstuff signature ;
  signature.get(is) ;
  ruleInfo.signature = signature.str() ;
  ruleInfo.lines.signature = line_no ;
  line_no += signature.num_lines() ;

  killsp() ;

  if(ruleInfo.rule_type == "apply") {
    if(is.peek() != '[') {
      throw parseError("apply rule missing '[operator]'") ;
    }
    nestedbracketstuff applyop ;
    applyop.get(is) ;
    ruleInfo.applyop = applyop.str() ;
    ruleInfo.lines.applyop = line_no ;
    line_no += applyop.num_lines() ;
    killsp() ;
  }

  while(is.peek() == ',') {
    is.get() ;

    killsp() ;

    if(!is_name(is)) {
      throw parseError("syntax error") ;
    }

    string s = get_name(is) ;
    if(s == "constraint") {
      nestedparenstuff stuff ;
      stuff.get(is) ;
      if(ruleInfo.constraint.empty()) {
        ruleInfo.constraint = stuff.str() ;
        ruleInfo.lines.constraint = line_no ;
      } else {
        ruleInfo.constraint += "," + stuff.str() ;
      }
      line_no += stuff.num_lines() ;
    } else if(s == "parametric") {
      nestedparenstuff stuff ;
      stuff.get(is) ;
      if(!ruleInfo.parametric.empty()) {
        throw parseError("syntax error: cannot specify more than one parametric variable") ;
      }
      ruleInfo.parametric = stuff.str() ;
      ruleInfo.lines.parametric = line_no ;
      line_no += stuff.num_lines() ;
    } else if(s == "conditional") {
      nestedparenstuff stuff ;
      stuff.get(is) ;
      if(!ruleInfo.conditional.empty()) {
        throw parseError("syntax error: cannot specify more than one conditional variable") ;
      }
      ruleInfo.conditional = stuff.str() ;
      ruleInfo.lines.conditional = line_no ;
      line_no += stuff.num_lines() ;
    } else if(s == "specialized") {
      ruleInfo.is_specialized = 1 ;
      ruleInfo.lines.specialized = line_no ;
    } else if(s == "option") {
      nestedparenstuff stuff ;
      stuff.get(is) ;
      ruleInfo.options.push_back(stuff.str()) ;
      ruleInfo.lines.options.push_back(line_no) ;
      line_no += stuff.num_lines() ;
    } else if(s == "inplace") {
      nestedparenstuff stuff ;
      stuff.get(is) ;
      ruleInfo.inplace.push_back(stuff.str()) ;
      ruleInfo.lines.inplace.push_back(line_no) ;
      line_no += stuff.num_lines() ;
    } else if(s == "comments") {
      nestedparenstuff stuff ;
      stuff.get(is) ;
      ruleInfo.comments.push_back(cleanupCommentQuoted(stuff.str())) ;
      ruleInfo.lines.comments.push_back(line_no) ;
      line_no += stuff.num_lines() ;
    } else if(s == "prelude") {
      ruleInfo.use_prelude = 1 ;
    } else {
      throw parseError("unknown rule modifier") ;
    }

    killsp() ;
  }

  if(ruleInfo.use_prelude) {
    nestedbracestuff stuff ;
    stuff.get(is) ;
    ruleInfo.prelude = stuff.str() ;
    ruleInfo.lines.prelude = line_no ;
    line_no += stuff.num_lines() ;

    killsp() ;
    if(is.peek() == ';') {
      is.get() ;
      ruleInfo.use_compute = 0 ;
    }
    if(is_name(is)) {
      string s = get_name(is) ;
      if(s != "compute") {
        throw parseError("syntax error, expecting 'compute'") ;
      }
    }
    killsp() ;
  }

  if(ruleInfo.use_compute) {
    nestedbracestuff stuff ;
    stuff.get(is) ;
    ruleInfo.compute = stuff.str() ;
    ruleInfo.lines.compute = line_no ;
    line_no += stuff.num_lines() ;
  }
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
        throw parseError("syntax error: cannot specify more than one parametric variable") ;
      }
        
      parametric_var = con.str() ;
      line_no += con.num_lines() ;
    } else if(s == "conditional") {
      nestedparenstuff con ;
      con.get(is) ;
      if(conditional != "") {
        throw parseError("syntax error: cannot specify more than one conditional variable") ;
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

void parseFile::initialize() {
  // Template for rule constructor.
  char const * rule_ctor = R"(
  $[rule.class]$() {
${each rule.name_stores}$
    name_store("$[name]$", $[vname]$) ;
  ${if has_info_id}$
    store_info_id("$[name]$", $[info_id]$) ;
  ${endif}$
${endeach}$
${each rule.inputs}$
    input("$[str]$") ;
${endeach}$
${each rule.outputs}$
    output("$[str]$") ;
${endeach}$
${each rule.constraints.spec}$
    constraint("$[str]$") ;
${endeach}$
${if rule.option_disable_threading}$
    disable_threading() ;
${endif}$
${if rule.is_parametric}$
    set_parametric_variable("$[rule.parametric.spec]$") ;
${endif}$
${if rule.is_specialized}$
    set_specialized() ;
${endif}$
${if rule.is_conditional}$
    conditional("$[rule.conditional.spec]$") ;
${endif}$
${each rule.comments}$
    comments("$[str]$") ;
${endeach}$
    setvardoc($[rule.docvar]$) ;
    set_file("$[rule.file]$:$[rule.line_number]$") ;
  })" ;
  cuda_templates.define("rule_ctor", rule_ctor) ;
  cpu_templates.define("rule_ctor", rule_ctor) ;

  // Template for rule store_instance declarations.
  cuda_templates.define("rule_store_decl", R"(
${each rule.input_stores}$
  Loci::const_gpu$[ctype]$$[if carg]$<$[carg]$>$[endif]$ $[vname]$ ;
${endeach}$
${each rule.output_stores}$
  Loci::gpu$[ctype]$$[if carg]$<$[carg]$>$[endif]$ $[vname]$ ;
${endeach}$)") ;

  cpu_templates.define("rule_store_decl", R"(
${each rule.input_stores}$
  Loci::const_$[ctype]$$[if carg]$<$[carg]$>$[endif]$ $[vname]$ ;
${endeach}$
${each rule.output_stores}$
  Loci::$[ctype]$$[if carg]$<$[carg]$>$[endif]$ $[vname]$ ;
${endeach}$)") ;

  // Template for pointwise rule.
  cuda_templates.define("pointwise_rule", R"(
class $[rule.class]$ : public Loci::$[rule.parent_class]$ {
  ${> rule_store_decl}$

public:
  ${> rule_ctor}$

  typedef struct {
${each rule.input_stores}$
    $[vtype]$ $[vname]$ ;
${endeach}$
${each rule.output_stores}$
    $[vtype]$ $[vname]$ ;
${endeach}$
    GPU_DECL void operator()(Entity _e_) {
$[rule.compute.spec]$
    }
  } compute_t ;

  void compute(Loci::sequence const & seq) override ;
} ;

__global__ void $[rule.class]$_kernel(
  int start, int stop, $[rule.class]$::compute_t cop
) {
  int _e_ = blockIdx.x*blockDim.x + threadIdx.x + start ;
  if(_e_ < stop) {
    cop(_e_) ;
  }
}

void $[rule.class]$::compute(Loci::sequence const & seq) {
  size_t const ni = seq.num_intervals() ;

  if(ni == 0) return ;

  compute_t cop ;
${each rule.input_stores}$
  cop.$[vname]$ = $[vname]$.ptr() ;
${endeach}$
${each rule.output_stores}$
  cop.$[vname]$ = $[vname]$.ptr() ;
${endeach}$
  
  for(size_t i = 0; i < ni; ++i) {
    Entity start = seq[i].first, stop = seq[i].second+1 ;
    int minGridSize, blockSize ;
    if(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, $[rule.class]$_kernel) != cudaSuccess) {
      std::cerr << "could not determine CUDA block size" << std::endl ;
      Loci::Abort() ;
    }
    int gridSize = ((stop-start) + blockSize - 1) / blockSize ;
${if debug_info}$
    nvtxRangePush("$[rule.debug_name]$") ;
${endif}$
    $[rule.class]$_kernel<<<gridSize, blockSize>>>(start, stop, cop) ;
${if debug_info}$
    nvtxRangePop() ;
${endif}$
  }
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;
)") ;

    // Template for parameter unit rule.
    cuda_templates.define("param_unit_rule", R"(
class $[rule.class]$ : public Loci::$[rule.parent_class]$ {
  ${> rule_store_decl}$

public:
  ${> rule_ctor}$

  typedef $[rule.unit.container_args]$ value_t ;

  typedef struct {
${each rule.input_stores}$
    $[vtype]$ $[vname]$ ;
${endeach}$
${each rule.output_stores}$
    $[vtype]$ $[vname]$ ;
${endeach}$

    GPU_DECL void operator()() {
$[rule.compute.spec]$
    }
  } compute_t ;

  void compute(Loci::sequence const & seq) override ;
} ;

__global__ void $[rule.class]$_kernel($[rule.class]$::compute_t cop) {
  if(threadIdx.x == 0)
    cop() ;
}

void $[rule.class]$::compute(Loci::sequence const & seq) {
  size_t const ni = seq.num_intervals() ;

  if(ni == 0) return ;
  
  compute_t cop ;
${each rule.input_stores}$
  cop.$[vname]$ = $[vname]$.ptr() ;
${endeach}$
${each rule.output_stores}$
  cop.$[vname]$ = $[vname]$.ptr() ;
${endeach}$

${if debug_info}$
  nvtxRangePush("$[rule.debug_name]$") ;
${endif}$
  $[rule.class]$_kernel<<<1, 1>>>(cop) ;
${if debug_info}$
  nvtxRangePop() ;
${endif}$
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;
)") ;

  cuda_templates.define("param_apply_rule", R"(
class $[rule.class]$ : public Loci::$[rule.parent_class]$<
  Loci::gpu$[rule.apply.container]$<$[rule.apply.container_args]$>,
  $[rule.apply.operator]$<$[rule.apply.reduction_args]$>
> {
  ${> rule_store_decl}$
public:
  ${> rule_ctor}$

  typedef $[rule.apply.container_args]$ value_t ;

  typedef $[rule.apply.operator]$<value_t> loci_reduction_t ;

  typedef struct {
    GPU_DECL
    value_t operator()(value_t const & lhs, value_t const & rhs) {
      value_t tmp = lhs ;
      loci_reduction_t op ;
      op(tmp, rhs) ;
      return tmp ;
    }

    GPU_DECL
    value_t identity() const {
      loci_reduction_t tmp ;
      return tmp.identity() ;
    }
  } reduction_t ;

  typedef struct {
${each rule.input_stores}$
    $[vtype]$ $[vname]$ ;
${endeach}$
${each rule.output_stores}$
    $[vtype]$ $[vname]$ ;
${endeach}$

    GPU_DECL
    value_t operator()(Entity _e_) {
$[rule.compute.spec]$
    }
  } compute_t ;

  void compute(Loci::sequence const & seq) override ;
} ;

__global__ void $[rule.class]$_reducevar_kernel(
  $[rule.class]$::value_t * res,
  $[rule.class]$::value_t const * part
) {
  if(threadIdx.x == 0) {
    $[rule.class]$::loci_reduction_t op ;
    op(*res, *part) ;
  }
}

void $[rule.class]$::compute(Loci::sequence const & seq) {
${if debug_info}$
  nvtxRangePush("$[rule.debug_name]$") ;
${endif}$
  size_t const ni = seq.num_intervals() ;
  if(ni == 0) return ;

  thrust::device_vector<value_t> rdata(ni+1) ;
  value_t * rptr = thrust::raw_pointer_cast(rdata.data()) ;

  thrust::host_vector<Entity> hint(ni*2) ;
  for(size_t i = 0; i < ni; ++i) {
    hint[i] = seq[i].first ;
    hint[i+ni] = seq[i].second+1 ;
  }
  thrust::device_vector<Entity> dint(hint) ;

  compute_t cop ;
${each rule.input_stores}$
  cop.$[vname]$ = $[vname]$.ptr() ;
${endeach}$

  reduction_t rop ;
  value_t const unit_value = rop.identity() ;

  thrust::counting_iterator eiter = thrust::make_counting_iterator(0) ;
  thrust::transform_iterator citer = thrust::make_transform_iterator(eiter, cop) ;

  void * tmpptr = nullptr ;
  size_t tmpsize = 0 ;

  cub::DeviceSegmentedReduce::Reduce(
    tmpptr, tmpsize, citer, rptr, ni,
    dint.begin(), dint.begin()+ni, rop, unit_value
  ) ;

  thrust::device_vector<std::uint8_t> dtmp(tmpsize) ;
  tmpptr = thrust::raw_pointer_cast(dtmp.data()) ;

  cub::DeviceSegmentedReduce::Reduce(
    tmpptr, tmpsize, citer, rptr, ni,
    dint.begin(), dint.begin()+ni, rop, unit_value
  ) ;

  tmpptr = nullptr ;
  tmpsize = 0 ;

  cub::DeviceReduce::Reduce(tmpptr, tmpsize, rptr, rptr+ni, ni, rop, unit_value) ;

  dtmp.resize(tmpsize) ;
  tmpptr = thrust::raw_pointer_cast(dtmp.data()) ;

  cub::DeviceReduce::Reduce(tmpptr, tmpsize, rptr, rptr+ni, ni, rop, unit_value) ;

  $[rule.class]$_reducevar_kernel<<<1, 1>>>($[rule.apply.target_vname]$.ptr(), rptr+ni) ;

${if debug_info}$
  nvtxRangePop() ;
${endif}$
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;
)") ;

  cuda_templates.define("singleton_param_apply_rule", R"(
class $[rule.class]$ : public Loci::$[rule.parent_class]$<
  Loci::gpu$[rule.apply.container]$<$[rule.apply.container_args]$>,
  $[rule.apply.operator]$<$[rule.apply.reduction_args]$>
> {
  ${> rule_store_decl}$
public:
  ${> rule_ctor}$

  typedef $[rule.apply.container_args]$ value_t ;

  typedef $[rule.apply.operator]$<value_t> loci_reduction_t ;

  typedef struct {
    GPU_DECL
    value_t operator()(value_t const & lhs, value_t const & rhs) {
      value_t tmp = lhs ;
      loci_reduction_t op ;
      op(tmp, rhs) ;
      return tmp ;
    }

    GPU_DECL
    value_t identity() const {
      loci_reduction_t tmp ;
      return tmp.identity() ;
    }
  } reduction_t ;

  typedef struct {
${each rule.input_stores}$
    $[vtype]$ $[vname]$ ;
${endeach}$
${each rule.output_stores}$
    $[vtype]$ $[vname]$ ;
${endeach}$
  } compute_t ;

  void compute(Loci::sequence const & seq) override ;
} ;

__global__ void $[rule.class]$_computevar_kernel(
  $[rule.class]$::value_t * target,
  $[rule.class]$::compute_t cop
) {
  if(threadId.x == 0) {
    $[rule.class]$::loci_reduction_t rop ;
    $[rule.class]$::value_t const part = cop() ;
    rop(*target, part) ;
  }
}

void $[rule.class]$::compute(Loci::sequence const & seq) {
  if(Loci::MPI_rank == 0) {
    compute_t cop ;
${each rule.input_stores}$
    cop.$[vname]$ = $[vname]$.ptr() ;
${endeach}$
    $[rule.class]$_computevar_kernel<<<1, 1>>>(
      $[rule.apply.target_vname]$.ptr(), cop
    ) ;
  }
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;
)") ;

  // Template for CPU pointwise rule.
  cpu_templates.define("pointwise_rule", R"(
namespace {

class $[rule.class]$ : public Loci::$[rule.parent_class]$ {
  ${> rule_store_decl}$

public:
  ${> rule_ctor}$

${if rule.prelude.is_specialized}$
  void prelude(Loci::sequence const & seq) override ;
${endif}$

${if rule.compute.is_specialized}$
  void calculate(Loci::Entity _e_) ;
${endif}$

  void compute(Loci::sequence const & seq) override ;
} ;

${if rule.prelude.is_specialized}$
void $[rule.class]$::prelude(Loci::sequence const & seq) {
$[rule.prelude.spec]$
}
${endif}$

${if rule.compute.is_specialized}$
void $[rule.class]$::calculate(Loci::Entity _e_) {
$[rule.compute.spec]$
}
${endif}$

void $[rule.class]$::compute(Loci::sequence const & seq) {
${if rule.compute.is_specialized}$
  do_loop(seq, this) ;
${endif}$
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;

} // end: anonymous namespace of rule
)") ;

  cpu_templates.define("unit_rule", R"(
namespace {

class $[rule.class]$ : public Loci::$[rule.parent_class]$ {
  ${> rule_store_decl}$

public:
  ${> rule_ctor}$

${if rule.prelude.is_specialized}$
  void prelude(Loci::sequence const & seq) override ;
${endif}$

${if rule.compute.is_specialized}$
  void calculate(Loci::Entity _e_) ;
${endif}$

  void compute(Loci::sequence const & seq) override ;
} ;

${if rule.prelude.is_specialized}$
void $[rule.class]$::prelude(Loci::sequence const & seq) {
$[rule.prelude.spec]$
}
${endif}$

${if rule.compute.is_specialized}$
void $[rule.class]$::calculate(Loci::Entity _e_) {
$[rule.compute.spec]$
}
${endif}$

void $[rule.class]$::compute(Loci::sequence const & seq) {
${if rule.compute.is_specialized}$
  do_loop(seq, this) ;
${endif}$
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;

} // end: anonymous namespace of rule
)") ;

  cpu_templates.define("apply_rule", R"(
namespace {

class $[rule.class]$ : public Loci::$[rule.parent_class]$<
  Loci::$[rule.apply.container]$<$[rule.apply.container_args]$>,
  $[rule.apply.operator]$<$[rule.apply.reduction_args]$>
> {
  ${> rule_store_decl}$

public:
  ${> rule_ctor}$

${if rule.prelude.is_specialized}$
  void prelude(Loci::sequence const & seq) override ;
${endif}$

${if rule.compute.is_specialized}$
  void calculate(Loci::Entity _e_) ;
${endif}$

  void compute(Loci::sequence const & seq) override ;
} ;

${if rule.prelude.is_specialized}$
void $[rule.class]$::prelude(Loci::sequence const & seq) {
$[rule.prelude.spec]$
}
${endif}$

${if rule.compute.is_specialized}$
void $[rule.class]$::calculate(Loci::Entity _e_) {
$[rule.compute.spec]$
}
${endif}$

void $[rule.class]$::compute(Loci::sequence const & seq) {
${if rule.compute.is_specialized}$
  do_loop(seq, this) ;
${endif}$
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;

} // end: anonymous namespace of rule
)") ;

  cpu_templates.define("param_unit_rule", R"(
namespace {

class $[rule.class]$ : public Loci::$[rule.parent_class]$ {
  ${> rule_store_decl}$

public:
  ${> rule_ctor}$

${if rule.prelude.is_specialized}$
  void prelude(Loci::sequence const & seq) override ;
${endif}$

  void compute(Loci::sequence const & seq) override ;
} ;

${if rule.prelude.is_specialized}$
void $[rule.class]$::prelude(Loci::sequence const & seq) {
$[rule.prelude.spec]$
}
${endif}$

void $[rule.class]$::compute(Loci::sequence const & seq) {
$[rule.compute.spec]$
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;

} // end: anonymous namespace of rule
)") ;

  cpu_templates.define("param_apply_rule", R"(
namespace {
class $[rule.class]$ : public Loci::$[rule.parent_class]$<
  Loci::$[rule.apply.container]$<$[rule.apply.container_args]$>,
  $[rule.apply.operator]$<$[rule.apply.reduction_args]$>
> {
  ${> rule_store_decl}$

public:
  ${> rule_ctor}$

${if rule.prelude.is_specialized}$
  void prelude(Loci::sequence const & seq) override ;
${endif}$

${if rule.compute.is_specialized}$
  void calculate(Loci::Entity _e_) ;
${endif}$

  void compute(Loci::sequence const & seq) override ;
} ;

${if rule.prelude.is_specialized}$
void $[rule.class]$::prelude(Loci::sequence const & seq) {
$[rule.prelude.spec]$
}
${endif}$

${if rule.compute.is_specialized}$
void $[rule.class]$::calculate(Loci::Entity _e_) {
$[rule.compute.spec]$
}
${endif}$

void $[rule.class]$::compute(Loci::sequence const & seq) {
${if rule.compute.is_specialized}$
  do_loop(seq, this) ;
${endif}$
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;

} // end: anonymous namespace of rule
)") ;

  cpu_templates.define("singleton_param_apply_rule", R"(
namespace {
class $[rule.class]$ : public Loci::$[rule.parent_class]$<
  Loci::$[rule.apply.container]$<$[rule.apply.container_args]$>,
  $[rule.apply.operator]$<$[rule.apply.reduction_args]$>
> {
  ${> rule_store_decl}$

public:
  ${> rule_ctor}$

${if rule.prelude.is_specialized}$
  void prelude(Loci::sequence const & seq) override ;
${endif}$

  void compute(Loci::sequence const & seq) override ;
} ;

${if rule.prelude.is_specialized}$
void $[rule.class]$::prelude(Loci::sequence const & seq) {
$[rule.prelude.spec]$
}
${endif}$

void $[rule.class]$::compute(Loci::sequence const & seq) {
${if rule.compute.is_specialized}$
  if(Loci::MPI_rank == 0) {
$[rule.compute.spec]$
  }
${endif}$
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;

} // end: anonymous namespace of rule
)") ;

  cpu_templates.define("blackbox_rule", R"(
namespace {
class $[rule.class]$ : public Loci::$[rule.parent_class]$ {
  ${> rule_store_decl}$

public:
  ${> rule_ctor}$

${if rule.prelude.is_specialized}$
  void prelude(Loci::sequence const & seq) override ;
${endif}$

  void compute(Loci::sequence const & seq) override ;
} ;

${if rule.prelude.is_specialized}$
void $[rule.class]$::prelude(Loci::sequence const & seq) {
$[rule.prelude.spec]$
}
${endif}$

void $[rule.class]$::compute(Loci::sequence const & seq) {
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;
} // end: anonymous namespace of rule
)") ;


  char const * sodc_rule = R"(
namespace {

class $[rule.class]$ : public Loci::$[rule.parent_class]$ {
  ${> rule_store_decl}$

public:
  ${> rule_ctor}$

  void compute(Loci::sequence const & seq) override ;
} ;

void $[rule.class]$::compute(Loci::sequence const & seq) {
$[rule.compute.spec]$
}

Loci::register_rule<$[rule.class]$> register_$[rule.class]$ ;

} // end: anonyous namespace of rule
)" ;

  cpu_templates.define("singleton_rule", sodc_rule) ;
  cpu_templates.define("optional_rule", sodc_rule) ;
  cpu_templates.define("default_rule", sodc_rule) ;
  cpu_templates.define("constraint_rule", sodc_rule) ;

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
            if(parseInfo.test_parse) {
              parseRuleInfo ruleInfo ;
              parse_rule_info(false, parseInfo, ruleInfo) ;
              process_and_validate_rule_info(
                comment, docvarname, parseInfo, ruleInfo
              ) ;
              render_rule(outputFile, parseInfo, ruleInfo) ;
            } else {
              setup_Rule(outputFile,comment,parseInfo) ;
            }
	  } else if(key == "cudarule") {
            if(level != 0) {
              throw parseError("$rule is not allowed in include file!") ;
            }
	    if(parseInfo.no_cuda) {
              if(parseInfo.test_parse) {
                parseRuleInfo ruleInfo ;
                parse_rule_info(false, parseInfo, ruleInfo) ;
                process_and_validate_rule_info(
                  comment, docvarname, parseInfo, ruleInfo
                ) ;
                render_rule(outputFile, parseInfo, ruleInfo) ;
              } else {
                setup_Rule(outputFile,comment,parseInfo) ;
              }
	    } else {
              parseRuleInfo ruleInfo ;
              parse_rule_info(true, parseInfo, ruleInfo) ;
              process_and_validate_rule_info(
                comment, docvarname, parseInfo, ruleInfo
              ) ;
              render_rule(outputFile, parseInfo, ruleInfo) ;
            }
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
