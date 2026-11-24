#ifndef TOOLS_NAMELIST_H
#define TOOLS_NAMELIST_H

#include <Config/conf.h>
#include <Tools/expr.h>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace Loci {
  class nameList {
    std::vector<std::string> names ;
  public:
    nameList() {}

    void clear() { names.clear() ; }
    std::vector<std::string>::const_iterator begin() const
    { return names.begin() ; }
    std::vector<std::string>::const_iterator end() const
    { return names.end() ; }
    bool empty() const { return names.size() == 0 ; }
    
    std::istream &Input(std::istream &s) {
      names.clear() ;
      parse::kill_white_space(s) ;
      bool openbracket = false ;
      if(s.peek() == '[') {
        openbracket = true ;
        s.get() ;
      }
      parse::kill_white_space(s) ;
      if(openbracket && s.peek() == ']') {
        s.get() ;
        return s ;
      }

      if(parse::is_name(s)) {
        std::string name = parse::get_name(s) ;
        names.push_back(name) ;
        parse::kill_white_space(s) ;
        while(s.peek() == ',') {
          s.get() ;
          parse::kill_white_space(s) ;
          if(parse::is_name(s)) {
            name = parse::get_name(s) ;
            names.push_back(name) ;
          } else {
            throw exprError(std::string("Syntax Error"),
                            std::string("expected name in nameList"),
                            Loci::exprErrorType::ERR_SYNTAX) ;
          }
        }
        if(openbracket) {
          parse::kill_white_space(s) ;
          if(s.peek() != ']')
            throw exprError(std::string("Syntax Error"),
                            std::string("expected ']' in nameList"),
                            Loci::exprErrorType::ERR_SYNTAX) ;
          else
            s.get() ;
        }
      } else {
        throw exprError(std::string("Syntax Error"),
                        std::string("expected name in nameList"),
                        Loci::exprErrorType::ERR_SYNTAX) ;
      }
      return s ;
    }
    std::ostream &Print(std::ostream &s) const {
      if(names.size() == 0) {
        s << "[]" << std::endl ;
      } else {
        s << "[" << names[0] ;
        for(size_t i=1;i<names.size();++i)
          s << ", " << names[i] ;
        s <<"]" << std::endl ;
      }
      return s ;
    }
  } ;

  inline std::istream & operator >> (std::istream &s, nameList &n) {
    return n.Input(s) ;
  }
  inline std::ostream & operator << (std::ostream &s, const nameList &n) {
    return n.Print(s) ;
  }

}

#endif
