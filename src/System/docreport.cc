//#############################################################################
//#
//# Copyright 2025, Loci Consulting Services, LLC
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
#include <Loci>
#include <Tools/debug.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include <stdlib.h>
#include <string.h>

using std::cout ;
using std::cerr ;
using std::endl ;
using std::ios ;
using std::ifstream ;
using std::ostream ;
using std::ofstream ;
using std::string ;
using std::ostringstream ;
using std::vector ;
using std::map ;

namespace Loci {
  void prettyPrintString(string i1, string i2,string s, ostream &o) {
    o << i1<<i2  << ": " ;

    size_t initial_space = i2.size() + 2 ;
    size_t count = initial_space ;
    string current_word ;
    size_t loc = 0 ;
    do {
      current_word = "" ;
      while(loc < s.size() && s[loc] != ' ')
        current_word += s[loc++] ;
      if(loc < s.size()) 
        loc++ ;
    
      if(count+current_word.size() >= 77) {
        o << endl ;
        o << i1 ;
        for(size_t t=0;t<initial_space; ++t)
          o << ' ' ;
        count = initial_space +i1.size() ;
      }
      o << current_word << ' ' ;
      count += current_word.size()+1 ;
    } while (loc < s.size()) ;
    o << endl ;
  }
  void prettyPrintString(string i, string s, ostream &o) {
    prettyPrintString("",i,s,o) ;
  }

  void describeInputs(const rule_db &rdb,ostream &o) {
    fact_db local ;
    // first of all, we need to process the default and optional rules
    ruleSet special_rules = rdb.get_default_rules() ;
    // first we process the default rules
    o << "------------------------------------------------------------------------------" << endl ;
    for(ruleSet::const_iterator ri=special_rules.begin();
        ri!=special_rules.end();++ri) {
      // first we need to create the facts in the fact_db
      variableSet targets = ri->targets() ;
      rule_implP rp = ri->get_rule_implP() ;
      for(variableSet::const_iterator vi=targets.begin();
          vi!=targets.end();++vi) {
        // we need to get the storeRep for this variable
        storeRepP srp = rp->get_store(*vi) ;
        if(srp == 0) {
          cerr << "rule " << *ri << " unable to provide type for " << *vi
               << endl ;
          exit(-1) ;
        }
        local.create_fact(*vi,srp) ;
      }
      // then we need to call the compute method to set
      // the default value for this variable
      rp->initialize(local) ;
      rp->compute(sequence(EMPTY)) ;
      for(variableSet::const_iterator vi=targets.begin();
          vi!=targets.end();++vi) {
        storeRepP srp = local.get_variable(*vi) ;
        o << *vi << ": " ;
        srp->Print(o) ;
        string comment = rp->get_comments() ;
        if(comment.size() != 0)
          prettyPrintString("comment",comment,o) ;
        o << "------------------------------------------------------------------------------" << endl ;
      }
    }
    // then we process the optional rules
    special_rules = rdb.get_optional_rules() ;
    for(ruleSet::const_iterator ri=special_rules.begin();
        ri!=special_rules.end();++ri) {
      // first we need to create the facts in the fact_db
      variableSet targets = ri->targets() ;
      rule_implP rp = ri->get_rule_implP() ;
      for(variableSet::const_iterator vi=targets.begin();
          vi!=targets.end();++vi) {
        o << *vi << ": NO DEFAULT VALUE" << endl ;
        string comment = rp->get_comments() ;
        if(comment.size() != 0)
          prettyPrintString("comment",comment,o) ;
        o << "------------------------------------------------------------------------------" << endl ;
      }
    }
  }

  enum ruleCategories {
    RULE_POINTWISE=0x1, RULE_UNIT=0x2, RULE_APPLY=0x4,
    RULE_SINGLETON=0x8, RULE_DEFAULT=0x10, RULE_OPTIONAL=0x20,
    RULE_CONSTRAINT=0x40, RULE_BLACKBOX=0x80, RULE_SYSTEM=0x100
  } ;
  ruleCategories ruleType(rule_impl::rule_impl_type t) {
    switch(t) {
    case rule_impl::POINTWISE:
      return RULE_POINTWISE ;
    case rule_impl::SINGLETON:
      return RULE_SINGLETON ;
    case rule_impl::UNIT:
      return RULE_UNIT ;
    case rule_impl::APPLY:
      return RULE_APPLY ;
    case rule_impl::DEFAULT:
      return RULE_DEFAULT ;
    case rule_impl::OPTIONAL:
      return RULE_OPTIONAL ;
    case rule_impl::CONSTRAINT_RULE:
      return RULE_CONSTRAINT ;
    case rule_impl::BLACKBOX_RULE:
      return RULE_BLACKBOX ;
    default:
      return RULE_SYSTEM ;
    }
  }
  

    
  variable basicVariable(variable v) {
    variable::info vinfo = v.get_info() ;
    // remove priority strings
    vinfo.priority = std::vector<std::string>() ;
    // If this is a parametric variable, then rename them to generic
    // terms
    for(size_t i=0;i<vinfo.v_ids.size();++i) {
      std::ostringstream ss ;
      ss << 'X' << i << endl ;
      variable xi = variable(ss.str()) ;
      vinfo.v_ids[i] = xi.ident() ;
    }
    // ignore time spefier
    vinfo.time_id = time_ident() ;
    vinfo.offset = 0 ;
    vinfo.assign = false ;
    return variable(vinfo) ;
  }

  void prettyPrintVarInfo(variable v, rule_implP rp, ostream &o) {
    ostringstream oss ;
    oss << v ;
    const char * vardoc = rp->getvardoc(v) ;
    const char *varfile = vardoc ;
    vardoc += strlen(vardoc)+1 ;
    //const char *varname = vardoc ;
    vardoc += strlen(vardoc)+1 ;
    const char *vartype = vardoc ;
    vardoc += strlen(vardoc)+1 ;
    string comment = string(vardoc) ;
    if(comment.size()>0)
      prettyPrintString("== ",oss.str(),comment,o) ;
    if(strlen(varfile)>0) {
      if(strlen(vartype) > 0)
        comment = "$type " + string(vartype) ;
      comment += " defined in " + string(varfile)  ;
      prettyPrintString("== ",oss.str(),comment,o) ;
    }
  }

  bool checkParametric(rule r, variable &parametric_var) {
    rule_implP rp = r.get_rule_implP() ;
    bool is_parametric = false ;
    variableSet targets = r.targets() ;
    if(targets.size()==1 && targets.begin()->get_info().v_ids.size()>0) {
      if(!rp->is_specialized()) {
        is_parametric = true ;
        parametric_var = *targets.begin() ;
      }
    }

    if(rp->is_parametric_provided()) {
      is_parametric = true ;
      parametric_var = rp->get_parametric_variable() ;
    }
    return is_parametric ;
  }
  
  void prettyPrintRuleDoc(rule r,ostream &o) {
    o << "== ____________________________________________________________________________" << endl ;
    rule_implP rp = r.get_rule_implP() ;
    switch(rp->get_rule_class()) {
    case rule_impl::POINTWISE:
      o<< "== pointwise ";
      break ;
    case rule_impl::SINGLETON:
      o<< "== singleton " ;
      break ;
    case rule_impl::UNIT:
      o<< "== unit " ;
      break ;
    case rule_impl::APPLY:
      o<< "== apply " ;
      break ;
    case rule_impl::DEFAULT:
      o<< "== default " ;
      break ;
    case rule_impl::OPTIONAL:
      o<< "== optional " ;
      break ;
    case rule_impl::CONSTRAINT_RULE:
      o<< "== constraint " ;
      break ;
    case rule_impl::BLACKBOX_RULE:
      o<< "== blackbox " ;
      break ;
    default:
      o<< "== system " ;
      break ;
    }
    o << "rule signature:" << endl ;
    o << "== " << r << endl ;
    bool is_specialized = rp->is_specialized() ;
    variable parametric_var ;
    variableSet targets = r.targets() ;
    bool is_parametric = checkParametric(r,parametric_var) ;
                            
    if(is_parametric) 
      o << "== Rule is parametric on " << parametric_var  << endl ;

    if(is_specialized) 
      o << "== Rule is a specialized instantiation" << endl ;

    string fileloc = rp->get_fileloc() ;
    if(fileloc.size()>0) 
      o << "== Defined in " << fileloc << endl ;
    string comments = rp->get_comments() ;
    if(comments.size()>0)
      prettyPrintString("=="," Comment",comments,o) ;
    
    if(strlen(rp->getvardoc(*targets.begin())) > 0) {
      // We have variable documentation
      if(targets.size()>1)
        o << "== >>>>>>>>>>>>>>>>>>>>>>>>>>>> outputs >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl ;
      else
        o << "== >>>>>>>>>>>>>>>>>>>>>>>>>>>> output >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl ;
      for(auto vi=targets.begin();vi!=targets.end();++vi) 
        prettyPrintVarInfo(*vi,rp,o) ;

      variableSet sources = r.sources() ;
      if(sources.size()>0) {
        if(sources.size()>1)
          o << "== <<<<<<<<<<<<<<<<<<<<<<<<<<<< inputs <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl ;
        else
          o << "== <<<<<<<<<<<<<<<<<<<<<<<<<<<< input <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl ;
        for(auto vi=sources.begin();vi!=sources.end();++vi)
          prettyPrintVarInfo(*vi,rp,o) ;
      }
    }
  }
  
  void ruleDocumentation(const rule_db &rdb, variableSet vars, ostream &o) {
    ruleSet rules = rdb.all_rules() ;
    if(vars == EMPTY) {
      // The empty set implies all of the variables
      for(auto ri = rules.begin();ri!=rules.end();++ri) {
        vars += ri->sources() ;
        vars += ri->targets() ;
        vars += ri->constraints() ;
      }
    }

    variableSet basicset ;
    for(auto vi=vars.begin();vi!=vars.end();++vi)
      basicset += basicVariable(*vi) ;

    // sort out the parametric variables to process separately
    map<variable,ruleSet> parametric_rules ;
    variableSet parametric_targets ;
    for(auto ri = rules.begin();ri!=rules.end();++ri) {
      variable  parametric_var ;
      bool is_parametric = checkParametric(*ri,parametric_var) ;
      parametric_var = basicVariable(parametric_var) ;
      if(ri->targets().size()==1 &&
         ri->targets().begin()->get_info().v_ids.size()>0) {
        is_parametric=true ;
        parametric_var = basicVariable(*ri->targets().begin()) ;
      }
      if(is_parametric && basicset.inSet(parametric_var)) {
        parametric_rules[parametric_var] += *ri ;
        parametric_targets += ri->targets() ;
      }
    }
    for(auto vi=parametric_targets.begin();vi!=parametric_targets.end();++vi) {
      variable v = basicVariable(*vi) ;
      basicset -= v ;
    }
    
    map<string,ruleSet> printSet ;
    for(auto ri = rules.begin();ri!=rules.end();++ri) {
      variableSet targets = ri->targets() ;
      for(auto vi = targets.begin();vi != targets.end();++vi) {
        variable v = basicVariable(*vi) ;
        if(basicset.inSet(v)) {
          // Ok collate the rules for this variable
          ostringstream oss ;
          oss << v ; //<< ":" << vardoc ;
          printSet[oss.str()] += *ri ;
        }
      }
    }
    for(auto mi= parametric_rules.begin();mi!=parametric_rules.end();++mi) {
      variable v = mi->first ;
      ostringstream oss ;
      oss << "0"<< v ; 
      printSet[oss.str()] += mi->second ;
    }
    
    for(auto mi=printSet.begin();mi!=printSet.end();++mi) {
      ruleSet rs = mi->second ;
      variable vt ;
      if(mi->first[0] != '0')
        vt = variable(mi->first) ;
      else
        vt = variable(mi->first.substr(1)) ;
      variableSet vars ;
      for(auto ri=rs.begin();ri!=rs.end();++ri) {
        variableSet targets = ri->targets() ;
        
        for(auto vi=targets.begin();vi!=targets.end();++vi)
          if(vt == basicVariable(*vi))
            vars += *vi ;
        variable parametric_var ;
        if(checkParametric(*ri,parametric_var))
          vars += parametric_var ;
      }
      FATAL(vars.size() == 0) ;
      variable vselect = *vars.begin() ;
      for(auto vi=vars.begin();vi!=vars.end();++vi) {
        if(vi->get_info().priority.size() < vselect.get_info().priority.size())
          vselect = *vi ;
      }
      int ruleTypes = 0 ;
      // fall back if it doesn't find
      const char *dummy = "\000\000\000\000" ;
      const char *vardoc = dummy;
      rule unit ; //= *rs.begin() ;
      for(auto ri=rs.begin();ri!=rs.end();++ri) {
        rule_implP rp = ri->get_rule_implP() ;
        int type = ruleType(rp->get_rule_class()) ;
        ruleTypes |= type ;
        if(type ==ruleCategories::RULE_UNIT)
          unit = *ri ;
        if(ri->targets().inSet(vselect)) {
          vardoc=rp->getvardoc(vselect) ;
          if(*vardoc != '\0')
            break ;
        }
      }

      const char *varfile = vardoc ;
      vardoc += strlen(vardoc)+1 ;
      const char *varname = vardoc ;
      vardoc += strlen(vardoc)+1 ;
      const char *vartype = vardoc ;
      vardoc += strlen(vardoc)+1 ;
      o << "===============================================================================" << endl ;
      if(varfile==dummy) {
        o << "== Parametric Group based on unbound parametric variable: " ;
      } else if(ruleTypes & (ruleCategories::RULE_UNIT|
                             ruleCategories::RULE_APPLY)) 
        o << "== Reduction Variable: " ;
      else if(ruleTypes &(ruleCategories::RULE_DEFAULT|
                          ruleCategories::RULE_OPTIONAL))
        o << "== Vars File Input Variable: " ;
      else
        o << "== Variable: " ;
      if(*varname!='\0')
        o << varname <<  endl ;
      else
        o << vselect << endl ;
      if(*vartype != '\0')
        o << "== Type: " << vartype << endl ;
      if(*varfile != '\0')
         o << "== Defined in " << varfile << endl ;
      if(*vardoc != '\0') {
        prettyPrintString("=="," Comment",vardoc,o) ;
      }
      // Now document rules computing variable
      o << "==-----------------------------------------------------------------------------" << endl ;
      if(vselect.get_info().v_ids.size()>0)
        o << "== Rules parameterized by variable " << vselect << endl ;
      else
        o << "== Rules computing variable:"  << endl ;
      if(ruleTypes & ruleCategories::RULE_UNIT) {
        rs -= unit ;
        prettyPrintRuleDoc(unit,o) ;
      }
      map<string,rule> sortrules ;
      for(auto ri=rs.begin();ri!=rs.end();++ri) {
        rule_implP rp = ri->get_rule_implP() ;
        sortrules[rp->get_name()] = *ri ;
      }
      for(auto ri=sortrules.begin();ri!=sortrules.end();++ri) {
        prettyPrintRuleDoc(ri->second,o) ;
      }
    }
    o << "===============================================================================" << endl ;
    
  }

}
