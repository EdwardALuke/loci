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
#include <param_rule.h>

using std::string ;
using std::ostringstream ;
#include <distribute.h>

//#define VERBOSE

namespace Loci {
  using std::pair ;
  using std::string ;
  using std::map ;
  using std::vector ;
  
  typedef pair<string,int> param_rule_key ;
  typedef map<param_rule_key, ruleSet> param_rule_db ;

  string var2key(variable v) {
    if(v.get_info().namespac.size() == 0) {
      return v.get_info().name ;
    }
    vector<string> ns = v.get_info().namespac ;
    string n ;
    for(size_t i = 0;i<ns.size();++i)
      n+= ns[i] + '@' ;
    n += v.get_info().name ;
    return n ;
  }
  
  // Return set of parametric variables that were accessed by the given ruleSet
  variableSet scanRulesForParametrics(ruleSet rs) {
    ruleSet::const_iterator ri ;
    variableSet pvars ;
    for(ri=rs.begin();ri!=rs.end();++ri) {
      variableSet sources = ri->sources() ;
      variableSet::const_iterator vi ;
      for(vi=sources.begin();vi!=sources.end();++vi)
        if(variable(*vi).get_arg_list().size() > 0) {
          variable v = *vi ;
          while(v.time() != time_ident())
            v=v.parent() ;
          while(v.get_info().priority.size() != 0)
            v = v.drop_priority() ;

          pvars += v ;
        }
    }
    return pvars ;
  }

  // -----------------------------------------------------------------------
  ///
  /// @brief substituteString() is used to do substitution of underscore
  ///   delimited substrings (or the name itself) for parametric rule
  ///   instantiation
  ///
  /// @param [name] input string to be substituted
  /// @param [transform_map] input transformation map that describes
  ///   the substitutions needed to instatiate the parametric rule
  /// @param [subst] output of substitution variables applied that may be 
  ///   used in later stages if substituted variable has parametric
  ///   arguments
  /// @return string that has been transformed according to transform_map
  ///
  // -----------------------------------------------------------------------
  string substituteString(string name,
                          map<string,variable> &transform_map,
                          vector<variable> &subst) {
    vector<string> vname ;
    vname.push_back(string()) ;
    for(size_t i=0;i<name.size();++i) {
      if(name[i] == '_')
        vname.push_back(string()) ;
      else
        vname.back() += name[i] ;
    }

    for(size_t i=0;i<vname.size();++i) {
      map<string,variable>::const_iterator ii = transform_map.find(vname[i]) ;
      if(ii != transform_map.end()) {
        subst.push_back(ii->second) ;
        vname[i] = ii->second.get_info().name ;
      }
    }
    string sname =vname[0] ;
    for(size_t i=1;i<vname.size();++i)
      sname += '_' + vname[i] ;
    return sname ;
  }
    
  // -----------------------------------------------------------------------
  ///
  /// @brief applySubstitution() is used to modify a variable according to
  ///   the parametric instation rules
  ///
  /// @param [transform_map] input transformation map that describes
  ///   the substitutions needed to instatiate the parametric rule
  /// @param [v] input variable to be transformed
  /// @return variable that has been transformed according to transform_map
  ///
  // -----------------------------------------------------------------------
  variable applySubstitution(map<string,variable> &transform_map, variable v) {
#ifdef VERBOSE2
    {
      debugout << "applySubstitution([ " ;
      map<string,variable>::const_iterator ii ;
      for(ii=transform_map.begin();ii!=transform_map.end();++ii)
	debugout << ii->first << ">"<<ii->second<<" ";
      debugout << "], " << v << ")" << endl ;
    }
#endif
    
      
    variable::info vinfo = v.get_info() ;
    // Substitute variable name with underscore separation
    string name = vinfo.name ;
    vector<variable> subst ;
    vinfo.name = substituteString(name,transform_map,subst) ;

    // Now substitute variable arguments, e.g. A substituted for grad(T) should
    // become grad(T).  The first step changes A to grad, this next step
    // adds the (T) argument
    vector<int> args = vinfo.v_ids ;
    for(size_t i=0;i<subst.size();++i) {
      if(subst[i].get_info().v_ids.size() != 0) {
        // If it was A(T) for example, then the substitution doesn't make sense
        if(args.size() != 0) {
          cerr << "problem with substitution argument conflict when processing variable " << v << endl ;
          cerr << "parametric rule may have irregular substitution" << endl;
        }
        args = subst[i].get_info().v_ids ;
      }
    }

    // If a parametric variable, then also apply substitution to arguments
    // of parametric variable, e.g. if X->temperature, then grad(X) should
    // become grad(temperature)
    if(v.get_info().v_ids.size()  != 0) {
      for(size_t i=0;i<args.size();++i) {
	variable tmp(args[i]) ;
	tmp = applySubstitution(transform_map,tmp) ;
	args[i] = tmp.ident() ;
      }
    } 
    vinfo.v_ids = args ;

    // Now update time variables
    // If we have a time annotation on a variable with underscores, substitute
    // here also, e.g. A{it_X} with X->temperature, should transform to
    // A{X_temperature}
    time_ident time_id = vinfo.time_id ;
    if(time_id != time_ident()) {
      vector<string> tlevel ;
      /// pull out time levels, check if any contain an underscore character
      bool has_underscore = false ;
      while(time_id != time_ident()) {
        tlevel.push_back(time_id.level_name()) ;
        if(tlevel.back().find('_') != std::string::npos)
          has_underscore = true ;
        time_id = time_id.parent() ;
      }
      if(has_underscore) {
        /// If there are underscores, then apply substitution to time identifier
        time_ident time_id_new ;
        for(auto i=tlevel.rbegin();i!=tlevel.rend();++i) {
          vector<variable> subst ;
          time_id_new = time_ident(time_id_new,
                                   substituteString(*i,transform_map,subst)) ;
        }
        vinfo.time_id = time_id_new ;
      }
    }
    // create variable based on edited vinfo
    return variable(vinfo) ;
  }
  
  rule instantiateRule(rule r, variable v) {
    variableSet target = r.targets() ;
    variable vt ;
    rule_implP rp = r.get_rule_implP() ;
    FATAL(rp == 0) ;
    if(rp->is_parametric_provided()) {
      vt = rp->get_parametric_variable() ;
    } else {
      if(target.size() != 1) {
        cerr << "can't handle parametric rules with multiple outputs!" << endl ;
        cerr << "rule = " << r << "variable " << v << endl ;
      }
      vt = *target.begin() ;
    }
    // Find substitution rules
    map<string,variable> transform_map ;
    std::vector<int> alt = vt.get_arg_list() ;
    std::vector<int> aln = v.get_arg_list() ;

    if(alt.size() != aln.size()) {
      cerr << "arg lists for parametric rule " << r << " don't match variable " << v << endl ;
      cerr << "this should not happen!" << endl ;
      Loci::Abort() ;
    }
#ifdef VERBOSE
    debugout << "instantiateRule("<<r<< ":[" << v<< "])" << endl ;
    debugout << "transform_map =" ;
#endif
    for(size_t i=0;i<alt.size();++i) {
      variable tv(alt[i]) ;
      variable sv(aln[i]) ;
      transform_map[tv.get_info().name] = sv ;
#ifdef VERBOSE
      debugout << " ("<< tv.get_info().name << "," << sv << ")" << endl ;
#endif
    }

#ifdef VERBOSE
    debugout << "generating variable renaming map" << endl ;
#endif
    // variable renaming map
    map<variable,variable> vm ;
    // Now loop over variables and apply substitution rules
    variableSet all_vars = target ;
    all_vars += r.sources() ;
    variableSet:: const_iterator vi ;

    for(vi=all_vars.begin();vi!=all_vars.end();++vi) {
#ifdef VERBOSE
      debugout << "result from applySubstitution("<< *vi << ") = " ;
      debugout.flush() ;
#endif
      
      variable sv = applySubstitution(transform_map,*vi) ;
#ifdef VERBOSE
      debugout << sv << endl ;
#endif
      if(sv != *vi) {
        vm[*vi] = sv ;
      }
    }

    rp->rename_vars(vm) ;
    return(rule(rp)) ;
  }

  ruleSet instantiateParametrics(param_rule_db &prule_db, variableSet pvars) {
    variableSet::const_iterator vi ;
    ruleSet newrules ;
    for(vi=pvars.begin();vi!=pvars.end();++vi) {
      variable v = *vi ;
      int psize = v.get_arg_list().size() ;
      param_rule_key rk(var2key(v),psize) ;
      param_rule_db::const_iterator pri = prule_db.find(rk) ;
      if(pri==prule_db.end()) {
        //        if(MPI_rank == 0)
        //          cerr << "WARNING: Unable to instantiate parametric rules for variable "
        //               << v << endl ;
        debugout << "Unable to instantiate parametric rules for variable "
                 << v << endl ;
        continue ;
      }
      ruleSet::const_iterator ri ;
      for(ri=pri->second.begin();ri!=pri->second.end();++ri) {
        newrules += instantiateRule(*ri,*vi) ;
      }
    }
    return newrules ;
  }
  
  rule_db parametric_rdb(const rule_db &rdb,variableSet query_vars) {
#ifdef VERBOSE
   debugout << "in parametric_rdb call with query_vars = "
	    << query_vars << endl ;
#endif
    Loci::rule_db par_rdb ;
    param_rule_db prule_db ; // Parametric rule db ;

    // put parametric rules in prule_db, all else in par_rdb
    ruleSet rset = rdb.all_rules() ;
    for(ruleSet::const_iterator rsi = rset.begin(); rsi != rset.end();++rsi) {
      variableSet target = rsi->targets() ;
      bool param=false ;
      rule_implP rp = rsi->get_rule_implP() ;
      if(rp == 0)
        continue ;
      if(!rp->is_specialized()) {
        if(rp->is_parametric_provided()) {
          variable pv = rp->get_parametric_variable() ;
#ifdef VERBOSE
	  debugout << "parametric variable provided = " << pv << endl ;
#endif
	  vector<int> vlist = pv.get_arg_list() ;
          int psize = vlist.size() ;
	  param_rule_key rk(var2key(pv),psize) ;
	  prule_db[rk] += *rsi ;
          param = true ;
        } else {
          for(variableSet::const_iterator vsi = target.begin(); vsi !=
                target.end(); ++vsi) {
	    variable pv = *vsi ;
	    vector<int> vlist = pv.get_arg_list() ;
	    int psize = vlist.size() ;
	    if(psize > 0) {
	      param_rule_key rk(var2key(*vsi),psize) ;
	      prule_db[rk] += *rsi ;
              param = true ;
            }
          }
        }
      }
      if(!param)
        par_rdb.add_rule(*rsi) ;
    }

    variableSet parvars = scanRulesForParametrics(par_rdb.all_rules()) ;
    for(variableSet::const_iterator vi=query_vars.begin();vi!=query_vars.end();++vi) {
      if(variable(*vi).get_arg_list().size() > 0) {
        parvars += *vi ;
      }
    }
    
#ifdef VERBOSE
    debugout << " entering instantiating loop: parvars = " << parvars << endl ;
#endif
    
    variableSet processed ;
    int cnt = 0 ;
    while(parvars != EMPTY) {
      processed += parvars ;
      ruleSet newrules = instantiateParametrics(prule_db,parvars) ;

#ifdef VERBOSE
      debugout << "generating rules = " << endl << newrules << endl ;
#endif
      parvars = scanRulesForParametrics(newrules) ;
      ruleSet::const_iterator ri ;
      for(ri=newrules.begin();ri!=newrules.end();++ri)
        par_rdb.add_rule(*ri) ;
      parvars -= processed ;
#ifdef VERBOSE
      debugout << "next parvars = " << parvars << endl ;
#endif
      cnt++ ;
      if(cnt == 100) {
        if(MPI_rank == 0) 
          cerr << "Warning, parametric rule instatiation 100 levels deep!" << endl
               << "Probably this is caused by a recursive parametric rule" << endl
               << "currently working on the parametric variables:" << endl
               << parvars << endl ;
      }
    }
    
#ifdef VERBOSE
    debugout << "finished parametric_rdb call" << endl ;
#endif
    return par_rdb ;
  }
}
