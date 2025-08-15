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
#ifndef LPP_H
#define LPP_H

#include "Tools/variable.h"
#include <list>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <map>

extern bool prettyOutput ;

extern Loci::variable convertVariable(Loci::variable v) ;

struct parseError {
  std::string error_type ;
  parseError(std::string errs) : error_type(errs) {}
} ;

struct parseSharedInfo {
  std::set<std::string> includedFiles ;
  std::vector<std::string> fileNameStack ;
  std::vector<std::string> dependFileList ;
  bool no_cuda ;
  parseSharedInfo() { no_cuda = true ; }
  
} ;
  

class parseFile {
  int cnt ;
  std::string filename ;
  int line_no ;
  std::ifstream is ;
  std::string docvarname ;

  struct typedoc {
    std::string container ;
    std::string container_args ;
    std::string filename ;
    std::string comment ;
    int lineno ;
    Loci::variable v ;
    std::string getFileLoc() const {
      std::ostringstream buf ;
      buf << filename << ':' << lineno ;
      return buf.str() ;
    }

    typedoc() { lineno = -1 ; }
  } ;
  std::map<Loci::variable,typedoc> type_map ;
  std::vector<typedoc> access_types ;
  std::map<std::string,int> access_map ;

    
  void addAccess(const typedoc &doc) {
    std::string key = doc.getFileLoc() ;
    if(access_map.find(key) == access_map.end()) {
      int id = access_types.size() ;
      access_types.push_back(doc) ;
      access_map[key] = id ;
    }
  }
  
  bool checkTypeValid(std::map<Loci::variable,typedoc>::const_iterator mi) {
    return(mi != type_map.end()) ;
  }
  std::map<Loci::variable,typedoc>::const_iterator
  lookupVarType(Loci::variable v) {
    v = convertVariable(v) ;
    auto mi  = type_map.find(v) ;
    if(!checkTypeValid(mi)) {
      v = v.new_offset(0) ;
      v = v.drop_assign() ;
      while(v.time() != Loci::time_ident())
        v = v.parent() ;

      mi = type_map.find(v) ;
      if(!checkTypeValid(mi)) {
        while(v.get_info().namespac.size() != 0)
          v = v.drop_namespace() ;
        mi = type_map.find(v) ;
      }
    }
    if(checkTypeValid(mi))
      addAccess(mi->second) ;
    else if(v.is_time_variable()) {
      typedoc data ;
      data.container = "param" ;
      data.container_args = "<int>" ;
      data.v = v ;
      type_map[v] = data ;
      return lookupVarType(v) ;
    }
    return mi ;
  }
  int killsp() ;
  std::string killspout(std::ostream &outputFile) ;

  void syncFile(std::ostream &outputFile) {
    if(!prettyOutput)
      outputFile << "#line " << line_no << " \"" << filename << "\"" << std::endl ;
  }

  void validate_VariableAccess(Loci::variable v, 
			       const std::list<Loci::variable> &vlist,
			       bool first_name,
			       const std::map<Loci::variable,std::string> &vnames,
			       const std::set<std::list<Loci::variable> > &validate_set) ;

  std::string process_String(std::string instring,
			     const std::map<Loci::variable,std::string> &vnames,
			     const std::set<std::list<Loci::variable> > &validate_set) ;
  
  void process_SpecialCommand(std::ostream &outputFile,
                              const std::map<Loci::variable,std::string> &vnames,
                              int &openbrace) ;
  void process_Prelude(std::ostream &outputFile,
                       const std::map<Loci::variable,std::string> &vnames) ;
  void process_Compute(std::ostream &outputFile,
                       const std::map<Loci::variable,std::string> &vnames) ;
  void process_Calculate(std::ostream &outputFile,
                         const std::map<Loci::variable,std::string> &vnames,
                         const std::set<std::list<Loci::variable> > & validate_set) ;

  void setup_Type(std::ostream &outputFile, const std::string &comment) ;
  void setup_Untype(std::ostream &outputFile) ;
  void setup_Rule(std::ostream &outputFile,const std::string &comment) ;
  void setup_cudaRule(std::ostream &outputFile,const std::string &comment) ;

  void setup_Test(std::ostream &outputFile) ;
public:
  parseFile() {
    line_no = 0 ;
    cnt = 0 ;
    Loci::variable OUTPUT("OUTPUT") ;
    typedoc data ;
    data.container = "param" ;
    data.container_args = "<bool>" ;
    data.v = OUTPUT ;
    data.filename = "SYSTEM" ;
    data.lineno = 1 ;
    
    type_map[OUTPUT] = data ;
    Loci::variable EMPTY("EMPTY") ;
    data.container = "Constraint" ;
    data.container_args = "" ;
    data.v = EMPTY ;
    data.lineno = 2 ;
    type_map[EMPTY] = data ;
    Loci::variable UNIVERSE("UNIVERSE") ;
    data.v = UNIVERSE ;
    data.lineno = 3 ;
    type_map[UNIVERSE] = data ;
  }
  void processFile(std::string file, std::ostream &outputFile,
		   parseSharedInfo &parseInfo,int level = 0) ;
} ;

extern std::list<std::string> include_dirs ;

#endif
