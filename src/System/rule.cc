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
#include <istream>
#include <ostream>
#include <iostream>
#include <string>
#include <sstream>

#include "rule.h"
#include <fact_db.h>
#include <constraint.h>

#include <vector>
using std::vector ;
#include <set>
using std::set ;

using std::pair ;

using std::cerr ;
using std::endl ;
using std::istream ;
using std::ostream ;

using std::ostringstream ;
using std::string ;

//#define VERBOSE

namespace Loci {

  namespace {
    variableSet convert_set(const variableSet &vset, time_ident tl) {
      variableSet res ;
      for(variableSet::const_iterator i=vset.begin();i!=vset.end();++i)
        res += variable(*i,tl) ;
      return res ;
    }

    // prepend the time_ident to each of the variable in the Set
    variableSet prepend_set(const variableSet& vset,
                            const time_ident& tl) {
      variableSet res ;
      for(variableSet::const_iterator i=vset.begin();i!=vset.end();++i)
        res += variable(tl,*i) ;

      return res ;
    }

    vmap_info convert_vmap_info(const vmap_info &in, time_ident tl) {
      vmap_info res ;
      for(vector<variableSet>::const_iterator i=in.mapping.begin();
          i!=in.mapping.end(); ++i)
        res.mapping.push_back(convert_set(*i,tl)) ;
      res.var = convert_set(in.var,tl) ;
      for(vector<pair<variable, variable> >::const_iterator i=in.assign.begin();
          i!=in.assign.end();++i)
        res.assign.push_back(std::make_pair(variable(i->first,tl),
                                            variable(i->second,tl))) ;
      return res ;
    }

    // prepend time_ident to vmap_info
    vmap_info prepend_vmap_info(const vmap_info& in,
                                const time_ident& tl) {
      vmap_info res ;
      for(vector<variableSet>::const_iterator i=in.mapping.begin();
          i!=in.mapping.end();++i)
        res.mapping.push_back(prepend_set(in.var,tl)) ;

      res.var = prepend_set(in.var,tl) ;

      for(vector<pair<variable,variable> >::const_iterator
            i=in.assign.begin();i!=in.assign.end();++i)
        res.assign.push_back(std::make_pair(variable(tl,i->first),
                                            variable(tl,i->second))) ;
      return res ;
    }

    variableSet rename_set(const variableSet &vset,
                           const std::map<variable, variable> &rvm ) {
      typedef std::map<variable, variable>::const_iterator map_iter ;
      variableSet res ;
      for(variableSet::const_iterator i=vset.begin();i!=vset.end();++i) {
	map_iter mi ;
        mi = rvm.find(*i) ;
	if(mi != rvm.end())
	  res += mi->second ;
	else
	  res += *i ;
      }
      return res ;
    }
    vmap_info rename_vmap_info(const vmap_info &in,
                               const std::map<variable, variable> &rvm) {
      vmap_info res ;
      for(vector<variableSet>::const_iterator i=in.mapping.begin();
          i!=in.mapping.end(); ++i)
        res.mapping.push_back(rename_set(*i,rvm)) ;
      res.var = rename_set(in.var,rvm) ;
      for(vector<pair<variable, variable> >::const_iterator
            i=in.assign.begin();i!=in.assign.end();++i) {
	variable v1, v2 ;
	std::map<variable, variable>::const_iterator mi = rvm.find(i->first) ;
	if(mi != rvm.end())
	  v1 = mi->second ;
	else
	  v1 = i->first ;

	mi = rvm.find(i->second) ;
	if(mi != rvm.end())
	  v2 = mi->second ;
	else
	  v2 = i->second ;

	res.assign.push_back(std::make_pair(v1, v2)) ;
      }

      return res ;
    }

  }

  string rule_impl::info::rule_identifier() const {
    ostringstream ss ;
    ss << targets << "<-" << sources ;
    if(constraints.begin() != constraints.end()) {
      if(sources.begin() != sources.end())
        ss << "," ;
      ss << "CONSTRAINT(" << constraints << ")" ;
    }

    if(!dynamic_constraints.empty()) {
      if(!constraints.empty())
        ss << "," ;
      ss << "DYNAMIC_CONSTRAINT(" << dynamic_constraints << ")" ;
    }

    if(conditionals != EMPTY)
      ss << ",CONDITIONAL(" << conditionals << ")" ;
    return ss.str() ;
  }

  variableSet rule_impl::info::input_vars() const {
    variableSet input = conditionals ;
    set<vmap_info>::const_iterator i ;
    for(i=sources.begin();i!=sources.end();++i) {
      for(size_t j=0;j<i->mapping.size();++j)
        input += i->mapping[j] ;
      input += i->var ;
    }
    for(i=constraints.begin();i!=constraints.end();++i) {
      for(size_t j=0;j<i->mapping.size();++j)
        input += i->mapping[j] ;
      input += i->var ;
    }
    for(i=targets.begin();i!=targets.end();++i)
      for(size_t j=0;j<i->mapping.size();++j)
        input += i->mapping[j] ;

    return input ;
  }

  variableSet rule_impl::info::output_vars() const {
    variableSet output ;
    set<vmap_info>::const_iterator i ;

    for(i=targets.begin();i!=targets.end();++i)
      output += i->var ;

    return output ;
  }

  rule_impl::rule_impl() {
    name = "UNNAMED" ;
    rule_impl_class = UNKNOWN ;
    rule_threading = true ;
    gpgpu_kernel = false ;
    use_dynamic_schedule = false ;
    relaxed_recursion = false ;
    specialized_parametric = false ;
    use_parametric_variable = false ;
    use_prelude = false;
    use_postlude = false;
    vardoc = 0 ;
    // setting the keyspace tag to "main"
    // this is a temporary solution. in the
    // future, we would need a more systematic
    // way to categorize rule in keyspaces.
    // a more consistent and systematic implementation
    // also requires that the fact database to be
    // reorganized to support "true" keyspace partition
    space_tag = "main" ;
    space_dist = false ;
  }

  void rule_impl::rule_name(const string &fname) {
    name = fname ;
  }

  string rule_impl::get_name() const {
    if(name == "UNNAMED")
      name = typeid(*this).name() ;
    return name ;
  }

  rule_implP rule_impl::add_namespace(const string& n) const {
    rule_implP with_namespace = new_rule_impl();

    variableSet vars = with_namespace->get_var_list() ;
    std::map<variable,variable> new_vars;
    for(variableSet::variableSetIterator i=vars.begin();i!=vars.end();++i) {
      new_vars[*i] = i->add_namespace(n);
    }

    with_namespace->rename_vars(new_vars);
    return with_namespace ;
  }

  namespace {
    inline void fill_descriptors(set<vmap_info> &v, const exprList &in) {

      for(exprList::const_iterator i = in.begin();i!=in.end();++i) {
        vmap_info di(*i) ;
        if(v.find(di) != v.end()) {
          cerr << "Warning, duplicate variable in var set." << endl ;
	  cerr << "i=" << *i << endl ;
          cerr << "expr = " << di << endl ;
	  throw int(-1) ;
        } else
          v.insert(di) ;
      }
    }

    class NULL_RULE_IMPL: public rule_impl {
    public:
      NULL_RULE_IMPL() {}
      rule_implP new_rule_impl() const { return new NULL_RULE_IMPL ; }
      void compute(const sequence &seq) {}
      virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
    } ;
  }

  void rule_impl::source(const string &invar) {
    exprP p = expression::create(invar) ;
    try {
      fill_descriptors(rule_info.sources,
		       collect_associative_op(p,OP_COMMA)) ;
    } catch(int i) {
      cerr << "parsing source string " << invar ;
    }
  }

  void rule_impl::target(const string &outvar) {
    exprP p = expression::create(outvar) ;
    try {
      fill_descriptors(rule_info.targets,
		       collect_associative_op(p,OP_COMMA)) ;
    } catch(int i) {
      cerr << "parsing target string " << outvar << endl ;
    }
  }

  void rule_impl::constraint(const string &constrain) {
    exprP p = expression::create(constrain) ;
    try {
      fill_descriptors(rule_info.constraints,
		       collect_associative_op(p,OP_COMMA)) ;
    } catch(int i) {
      cerr << "parsing constraint string = " << constrain << endl ;
    }
  }

  void rule_impl::conditional(const string &cond) {
    exprP p = expression::create(cond) ;
    rule_info.conditionals += variableSet(p) ;
  }

  void rule_impl::set_store(variable v, const storeRepP &p) {
    typedef storeIMap::iterator SI ;
    std::pair<SI, SI> sip = var_table.equal_range(v) ;
    for(SI si = sip.first; si != sip.second; ++si) {
      fatal(si == var_table.end()) ;
      si->second->setRep(p) ;
    }
  }

  storeRepP rule_impl::get_store(variable v) const {
    //Print(cout) ;
    typedef storeIMap::const_iterator SI ;

    std::pair<SI, SI> sip = var_table.equal_range(v) ;
    SI sp = sip.first ;
    if(sip.first == sip.second) {
      set<vmap_info>::const_iterator vmsi ;
      for(vmsi=rule_info.targets.begin();
          vmsi!=rule_info.targets.end();
          ++vmsi) {
        for(size_t i=0;i<vmsi->assign.size();++i) {
	  if(vmsi->assign[i].first == v) {
	    sip = var_table.equal_range(vmsi->assign[i].second) ;
	    sp = sip.first ;
	    if(sip.first != sip.second)
	      return sp->second->Rep() ;
	  }
	}
      }
      return storeRepP(0) ;
    }
    return sp->second->Rep() ;
  }

  void rule_impl::prot_rename_vars(std::map<variable,variable>  &rvm){
#ifdef VERBOSE
    debugout << " ***********************************************" << endl ;
    debugout << "Before Prot_rename_vars " << endl ;
    Print(debugout) ;
    debugout << "rvm = " ;
    std::map<variable,variable>::const_iterator ip ;
    for(ip=rvm.begin();ip!=rvm.end();++ip) {
      debugout << '(' << ip->first << ',' << ip->second << ") "
               << '[' << ip->first.ident() << ' '
               << ip->first.time().ident() << ','
               << ip->second.ident() << ' '
               << ip->second.time().ident() << "] " ;
    }
    debugout << endl ;
#endif
    typedef storeIMap::iterator smap_iter ;
    typedef std::map<variable, variable>::const_iterator map_iter ;
    storeIMap tmp_var_table ;
    for(smap_iter si = var_table.begin(); si != var_table.end(); ++si) {
      map_iter mi = rvm.find(si->first) ;
      if(mi != rvm.end()) {
	std::pair<smap_iter, smap_iter> sip = var_table.equal_range(mi->first) ;
	for(smap_iter i = sip.first; i != sip.second; ++i) {
          //	  std::pair<variable, store_instance*> vsp = make_pair(mi->second, i->second) ;
          //	  tmp_var_table.insert(vsp) ;
	  storeRepP sp = (i->second)->Rep() ;
	  tmp_var_table.insert(std::pair<const variable, store_instance *>(mi->second,i->second)) ;
	}
      } else {
#ifdef VERBOSE
        debugout << "not renaming " << si->first
                 << " [" << si->first.ident() << ','
                 << si->first.time().ident() << "]" << endl ;
#endif
	tmp_var_table.insert(*si) ;
      }
    }
    var_table.swap(tmp_var_table) ;
    std::set<vmap_info>::const_iterator i ;
    std::set<vmap_info> tmp ;
    for(i = rule_info.sources.begin(); i != rule_info.sources.end(); ++i)
      tmp.insert(rename_vmap_info(*i, rvm)) ;
    rule_info.sources.swap(tmp) ;
    tmp.clear() ;
    for(i = rule_info.targets.begin(); i != rule_info.targets.end(); ++i)
      tmp.insert(rename_vmap_info(*i, rvm)) ;
    rule_info.targets.swap(tmp) ;
    tmp.clear() ;
    for(i=rule_info.constraints.begin();
        i!=rule_info.constraints.end();++i)
      tmp.insert(rename_vmap_info(*i, rvm)) ;
    rule_info.constraints.swap(tmp) ;
    tmp.clear() ;
    rule_info.conditionals = rename_set(rule_info.conditionals,rvm) ;
    if(use_parametric_variable) {
      std::map<variable,variable>::const_iterator mi ;
      mi = rvm.find(ParametricVariable) ;

      if(mi != rvm.end()) {
	ParametricVariable = mi->second ;
      }
    }

#ifdef VERBOSE
    debugout << " After prot rename vars " << endl ;
    Print(debugout) ;
    debugout << " *************************************************" << endl ;
#endif
  }

  void rule_impl::name_store(const string &nm, store_instance &si) {
    variable v(expression::create(nm)) ;
    var_table.insert(std::pair<const variable, store_instance *>(v,&si)) ;
  }


  bool rule_impl::check_perm_bits() const {

    variableSet read_set,write_set ;
    set<vmap_info>::const_iterator i ;

    for(i=rule_info.sources.begin();i!=rule_info.sources.end();++i){
      for(size_t j=0;j<i->mapping.size();++j)
        read_set += i->mapping[j] ;
      read_set += i->var ;
    }
    for(i=rule_info.targets.begin();i!=rule_info.targets.end();++i){
      for(size_t j=0;j<i->mapping.size();++j)
        read_set += i->mapping[j] ;
      write_set += i->var ;
      for(size_t j=0;j<i->assign.size();++j) {
        write_set -= i->assign[j].first ;
        write_set += i->assign[j].second ;
      }
    }

    bool retval = true ;
    for(storeIMap::const_iterator i=var_table.begin();
        i!=var_table.end();++i) {
      if(write_set.inSet(i->first)) {
        if(i->second->access() != store_instance::READ_WRITE) {
          cerr << "WARNING! read-only var '" << i->first
               << "' in target list for rule "
               << get_name() << endl ;
          retval = false ;
        }
      } else if(read_set.inSet(i->first)) {
        if(i->second->access() != store_instance::READ_ONLY) {
          cerr << "WARNING! read-write var '" << i->first
               << "' only in source list for rule "
               << get_name() << endl ;
          retval = false ;
        }
      } else {
        cerr << "WARNING! var '" << i->first
             << "' not in source or target lists for rule "
             << get_name() << endl ;
        retval = false ;
      }
    }

    variableSet::const_iterator si,sri ;
    for(si=read_set.begin();si!=read_set.end();++si) {
      if(var_table.find(*si) == var_table.end()) {
        cerr << "WARNING! var '" << *si << "' has not been named in"
             << " rule " << get_name() << endl ;
        retval = false ;
      }
    }

    for(si=write_set.begin();si!=write_set.end();++si) {
      storeIMap::const_iterator mi ;
      if((mi=var_table.find(*si)) != var_table.end()) {
        switch(rule_impl_class) {
        case POINTWISE:
          if(!isSTORE(mi->second->Rep()) &&
             ((si->get_info()).name != "OUTPUT")) {
            cerr << "-------------------------------------------------"<<endl;
            cerr << "Pointwise rule should have targets of store type."<<endl;
            cerr << "perhaps this rule should be a singleton_rule, or"<<endl;
            cerr << "apply_rule."<< endl ;
            cerr << "error occured for rule " << get_name()
                 << " and variable " << *si << endl ;
            cerr << "-------------------------------------------------"<<endl;
            retval = false ;
          }
          break;
        case DEFAULT:
        case OPTIONAL:
          if(!isPARAMETER(mi->second->Rep())) {
            cerr << "-------------------------------------------------"<<endl;
            cerr << "Default and optional rule should have targets" << endl;
            cerr << " of param. Perhaps this rule should be a" << endl;
	    cerr << "pointwise_rule, or apply_rule."<< endl ;
            cerr << "Error occured for rule " << get_name()
		 << " and variable " << *si << endl ;
            cerr << "-------------------------------------------------"<<endl;
            retval = false ;
          }
          break ;
        case DELETION:
	  if(!isPARAMETER(mi->second->Rep())) {
            cerr << "-------------------------------------------------"<<endl ;
            cerr << "Deletion rule should have only one target of " << endl ;
            cerr << "parameter. Error occured for rule "<<get_name()<<endl ;
            cerr << " and variable " << *si << endl ;
            cerr << "-------------------------------------------------"<<endl;
            retval = false ;
          }
          break ;
        case CONSTRAINT_RULE:
          if(mi->second->Rep()->RepType() != CONSTRAINT) {
            cerr << "-------------------------------------------------"<<endl;
            cerr << "Constraint rule should have targets" << endl;
            cerr << " of constraint. Perhaps this rule should be a" << endl;
	    cerr << "pointwise_rule, or apply_rule."<< endl ;
            cerr << "Error occured for rule " << get_name()
		 << " and variable " << *si << endl ;
            cerr << "-------------------------------------------------"<<endl;
            retval = false ;
          }
          break ;
        case MAP_RULE:
          if(mi->second->Rep()->RepType() != MAP) {
            cerr << "-------------------------------------------------"<<endl;
            cerr << "Map rule should have targets" << endl;
            cerr << " of map. Perhaps this rule should be a" << endl;
	    cerr << "pointwise_rule, or apply_rule."<< endl ;
            cerr << "Error occured for rule " << get_name()
		 << " and variable " << *si << endl ;
            cerr << "-------------------------------------------------"<<endl;
            retval = false ;
          }
          break ;
        case BLACKBOX_RULE:
          if(mi->second->Rep()->RepType() != BLACKBOX) {
            cerr << "-------------------------------------------------"<<endl;
            cerr << "Blackbox rule should have targets" << endl;
            cerr << " of blackbox. Perhaps this rule should be a" << endl;
	    cerr << "pointwise_rule, or apply_rule."<< endl ;
            cerr << "Error occured for rule " << get_name()
		 << " and variable " << *si << endl ;
            cerr << "-------------------------------------------------"<<endl;
            retval = false ;
          }
          break ;
        case SINGLETON:
          if(!isPARAMETER(mi->second->Rep()) &&
	     !isBLACKBOX(mi->second->Rep())) {
            cerr << "-------------------------------------------------"<<endl;
            cerr << "Singleton rule should have targets of param or" << endl;
            cerr << "blackbox type.  Perhaps this rule should be a" << endl;
	    cerr << "pointwise_rule, or apply_rule."<< endl ;
            cerr << "Error occured for rule " << get_name()
		 << " and variable " << *si << endl ;
            cerr << "-------------------------------------------------"<<endl;
            retval = false ;
          }
          for(sri=read_set.begin();sri!=read_set.end();++sri) {
            if((mi=var_table.find(*sri)) != var_table.end() &&
               !isPARAMETER(mi->second->Rep()) &&
               !isMAP(mi->second->Rep()) &&
	       !isBLACKBOX(mi->second->Rep())) {
	      cerr << "-------------------------------------------------"<<endl;
	      cerr << "Singleton rule should have sources of param or" << endl;
	      cerr << "blackbox type.  Perhaps this rule should be a" << endl;
	      cerr << "pointwise_rule, or apply_rule."<< endl ;
	      cerr << "Error occured for rule " << get_name()
		   << " and variable " << *sri << endl ;
	      cerr << "-------------------------------------------------"<<endl;
	      retval = false ;
            }
          }
          if(si->get_info().priority.size() != 0) {
            cerr << "-------------------------------------------------"<<endl;
            cerr << "Singleton rule cannot use priority override" << endl;
            cerr << "Error occured for rule " << get_name()
                 << " and variable " << *si << endl ;
            cerr << "-------------------------------------------------"<<endl;

          }
          break;
        default:
          break ;
        }
      } else {
        cerr << "WARNING! var '" << *si << "' has not been named in"
             << " rule " << get_name() << endl ;
        retval = false ;
      }
    }

    return retval ;
  }

  void rule_impl::initialize(fact_db &facts) {
    storeIMap::iterator sp ;

    for(sp=var_table.begin();sp!=var_table.end();++sp) {
      storeRepP srp = facts.get_variable(sp->first) ;
      if(srp == 0) {
        cerr << "ERROR!: rule_impl::initialize unable to extract '"
             << sp->first << "' from store data base."
             << endl ;
        cerr << "Error occured in rule '"
             << typeid(*this).name() << "'" << endl ;
        Loci::Abort() ;
      }
      sp->second->setRep(srp) ;
    }
  }

  void
  rule_impl::replace_map_constraints(fact_db& facts) {
    // NOTE: because the rule_info utilizes a set<vmap_info> as
    // the basic components, it is therefore more complicated
    // to replace any contents inside, since the iterator
    // to an std::set data-structure is always const (regardless
    // of "const_iterator" or "iterator"), we'll actually need to
    // insert a new vmap_info and erase the original one instead.

    // check to see if the constraints are Maps
    vector<vmap_info> new_vmap_info ;

    set<vmap_info>::iterator si=rule_info.constraints.begin() ;
    set<vmap_info>::iterator si_bak ;

    while(si != rule_info.constraints.end()) {
      // for the constraint data, the "mapping," "assign" fields
      // within the vmap_info structure should all be empty
      // If not, there is likely a problem there, and in our case
      // we wouldn't want to continue.
      if(!si->mapping.empty() || !si->assign.empty()) {
        ++si ;
        continue ;
      }
      // get all the constraint variable names
      variableSet maps ;
      variableSet new_constraints ;
      for(variableSet::const_iterator vi=si->var.begin();
          vi!=si->var.end();++vi) {
        // get the storeRepP in the fact_db for the var
        storeRepP srp = facts.get_variable(*vi) ;
        if(srp == 0)
          continue ;
        // see if it is a Map
        if(!isMAP(srp))
          continue ;
        // it is a map, we need to create a constraint
        // and replace it in the rule
        // get the name of *vi
        string map_name = (vi->get_info()).name ;
        // generate a new name
        string new_name = "__" + map_name + "_MAP_constraint__" ;
        // create a real constraint
        // NOTE: other rules might have already created it
        // then in this case, we don't create it
        storeRepP crp = facts.get_variable(new_name) ;
        if(crp == 0) {
          Loci::constraint mapc ;
          *mapc = srp->domain() ;
	  mapc.Rep()->setDomainKeySpace(srp->getDomainKeySpace()) ;
          // install it in fact_db
          facts.create_fact(new_name, mapc) ;
        }
        // record the changes
        maps += *vi ;
        new_constraints += variable(new_name) ;
      }
      if(maps != EMPTY) {
        // we know we need to replace the vmap_info
        vmap_info nv(*si) ;
        // modify the vmap_info var field
        for(variableSet::const_iterator vi=maps.begin();
            vi!=maps.end();++vi)
          nv.var -= *vi ;

        for(variableSet::const_iterator vi=new_constraints.begin();
            vi!=new_constraints.end();++vi)
          nv.var += *vi ;

        new_vmap_info.push_back(nv) ;

        // erase the original copy
        si_bak = si ;
        ++si_bak ;
        rule_info.constraints.erase(si) ;
        si = si_bak ;
      } else
        ++si ;
    }
    // insert new vmap_info objects
    for(vector<vmap_info>::const_iterator vi=new_vmap_info.begin();
        vi!=new_vmap_info.end();++vi)
      rule_info.constraints.insert(*vi) ;

  }

  void
  rule_impl::split_constraints(const variableSet& dc) {

    vector<vmap_info> new_constraints ;

    set<vmap_info>::iterator si=rule_info.constraints.begin() ;
    set<vmap_info>::iterator si_bak ;

    // split the constraints field into static & dynamic ones
    while(si != rule_info.constraints.end()) {
      variableSet local_dc ;
      // get the constraints variables
      for(variableSet::const_iterator vi=si->var.begin();
          vi!=si->var.end();++vi) {
        // since a constraints variables may appear in other forms
        // e.g., having an offset, or an assign operator, e.g.,
        // if we have a dynamic constraints A{n}, it may appear in
        // a rule as "A{n=0}" or "A{n+1}". therefore, we need to
        // drop possible assigns and offsets.
        variable nv = (vi->new_offset(0)).drop_assign() ;
        if(dc.inSet(nv))
          local_dc += *vi ;
      }
      if(local_dc != EMPTY) {
        // make a copy of vmap_info
        vmap_info dc_copy(*si) ;
        // modify its var field
        dc_copy.var = local_dc ;
        // insert it into the dynamic_constraints set
        rule_info.dynamic_constraints.insert(dc_copy) ;

        // construct a new copy of static constraints
        variableSet local_sc = variableSet(si->var - local_dc) ;
        if(local_sc != EMPTY) {
          vmap_info sc_copy(*si) ;
          sc_copy.var = local_sc ;
          new_constraints.push_back(sc_copy) ;
        }
        // then erase original copy in constraints set
        si_bak = si ;
        ++si_bak ;
        rule_info.constraints.erase(si) ;
        si = si_bak ;
      } else
        ++si ;
    }
    // finally push the new static constraints ones to the set
    rule_info.constraints.insert(new_constraints.begin(),
                                 new_constraints.end()) ;
  }


  variableSet rule_impl::get_var_list() {
    storeIMap::iterator sp ;
    set<vmap_info>::const_iterator i ;
    variableSet vset ;
    for(sp = var_table.begin(); sp != var_table.end(); ++sp)
      vset += sp->first ;

    for(i=rule_info.sources.begin();i!=rule_info.sources.end();++i) {
      for(vector<variableSet>::const_iterator vi = i->mapping.begin();
	  vi != i->mapping.end(); ++vi)
	vset += *vi ;
      vset += i->var ;
      for(vector<pair<variable, variable> >::const_iterator pi=i->assign.begin();
	  pi != i->assign.end();++pi) {
	vset += pi->first ;
	vset += pi->second ;
      }
    }
    for(i=rule_info.targets.begin();i!=rule_info.targets.end();++i) {
      for(vector<variableSet>::const_iterator vi = i->mapping.begin();
	  vi != i->mapping.end(); ++vi)
	vset += *vi ;
      vset += i->var ;
      for(vector<pair<variable, variable> >::const_iterator pi=i->assign.begin();
	  pi != i->assign.end();++pi) {
	vset += pi->first ;
	vset += pi->second ;
      }
    }
    for(i=rule_info.constraints.begin();
	i!=rule_info.constraints.end();++i) {
      for(vector<variableSet>::const_iterator vi = i->mapping.begin();
	  vi != i->mapping.end(); ++vi)
	vset += *vi ;
      vset += i->var ;
      for(vector<pair<variable, variable> >::const_iterator pi=i->assign.begin();
	  pi != i->assign.end();++pi) {
	vset += pi->first ;
	vset += pi->second ;
      }
    }
    vset += rule_info.conditionals ;

    return vset ;
  }

  void rule_impl::set_variable_times(time_ident tl) {
    set<vmap_info>::const_iterator i ;
    set<vmap_info> tmp ;
    for(i=rule_info.sources.begin();i!=rule_info.sources.end();++i)
      tmp.insert(convert_vmap_info(*i,tl)) ;
    rule_info.sources.swap(tmp) ;
    tmp.clear() ;
    for(i=rule_info.targets.begin();i!=rule_info.targets.end();++i)
      tmp.insert(convert_vmap_info(*i,tl)) ;
    rule_info.targets.swap(tmp) ;
    tmp.clear() ;
    for(i=rule_info.constraints.begin();
        i!=rule_info.constraints.end();++i)
      tmp.insert(convert_vmap_info(*i,tl)) ;
    rule_info.constraints.swap(tmp) ;
    tmp.clear() ;
    rule_info.conditionals = convert_set(rule_info.conditionals,tl) ;
    storeIMap tmp2 ;
    storeIMap::iterator j ;
    for(j=var_table.begin();j!=var_table.end();++j) {
      variable v(j->first,tl) ;
      //      std::pair<variable, store_instance*>  vsp = make_pair(v, j->second) ;
      //      tmp2.insert(vsp) ;
      tmp2.insert(std::pair<const variable, store_instance*>(v,j->second)) ;
    }
    var_table.swap(tmp2) ;
    name = rule_info.rule_identifier() ;
  }

  void rule_impl::copy_store_from(rule_impl &f) {
    storeIMap::iterator sp ;

    for(sp=var_table.begin();sp!=var_table.end();++sp)
      sp->second->setRep(f.get_store(sp->first)) ;
  }

  void rule_impl::Print(ostream &s) const {
    s << "------------------------------------------------" << endl;
    s << "--- rule " << get_name() << ", class = " ;
    switch(rule_impl_class) {
    case POINTWISE:
      s << "POINTWISE" ;
      break ;
    case SINGLETON:
      s << "SINGLETON" ;
      break ;
    case UNIT:
      s << "UNIT" ;
      break ;
    case APPLY:
      s << "APPLY" ;
      break ;
    case DEFAULT:
      s << "DEFAULT" ;
      break ;
    case OPTIONAL:
      s << "OPTIONAL" ;
      break ;
    case CONSTRAINT_RULE:
      s << "CONSTRAINT_RULE" ;
      break ;
    case BLACKBOX_RULE:
      s << "BLACKBOX_RULE" ;
      break ;
    case INSERTION:
      s << "INSERTION" ;
      break ;
    case DELETION:
      s << "DELETION" ;
      break ;
    default:
      s << "ERROR" ;
      break ;
    }
    s << endl ;

    string rule_impl_key = rule_info.rule_identifier() ;
    s << rule_impl_key << endl ;
    s << "rule_info.sources = " << rule_info.sources << endl ;
    s << "rule_info.targets = " << rule_info.targets << endl ;
    s << "rule_info.constraints = " << rule_info.constraints << endl ;
    s << "rule_info.conditionals = " << rule_info.conditionals << endl ;

    typedef storeIMap::const_iterator smap_iter ;
    for(smap_iter si = var_table.begin(); si != var_table.end(); ++si) {
      s << "var_table[" << si->first << "]  = " << endl ;
      std::pair<smap_iter, smap_iter> sip = var_table.equal_range(si->first) ;
      for(smap_iter i = sip.first; i != sip.second; ++i) {
	storeRepP sp = (i->second)->Rep() ;
      }
    }

    s << "------------------------------------------------" << endl;

  }

  vector<string>
  rule_impl::gather_keyspace_dist() const {
    vector<string> ks ;
    // right now, it is now the space the rule is in
    ks.push_back(space_tag) ;
    return ks ;
  }

  rule::rule_db *rule::rdb = 0 ;


  rule::info::info(const rule_implP &fp) {
    rule_impl = fp ;
    desc = fp->get_info() ;
    auto &fpdr = *fp ;
    impl_name = typeid(fpdr).name() ;
    rule_ident.append(impl_name) ;
    rule_ident.append("#") ;
    rule_ident.append(desc.rule_identifier()) ;

    set<vmap_info>::const_iterator i ;
    variableSet svars,tvars,tvar_types ;
    for(i=desc.sources.begin();i!=desc.sources.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      source_vars += (*i).var ;
      svars += (*i).var ;
    }
    for(i=desc.constraints.begin();i!=desc.constraints.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      constraint_vars += (*i).var ;
      source_vars += (*i).var ;
    }
    source_vars += desc.conditionals ;
    for(i=desc.targets.begin();i!=desc.targets.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      if((*i).assign.size() != 0) {
        variableSet v = (*i).var ;
        for(variableSet::const_iterator vi=v.begin();vi!=v.end();++vi) {
          bool t = true ;
          for(size_t k=0;k<(*i).assign.size();++k)
            if(*vi == (*i).assign[k].first) {
              tvar_types += (*i).assign[k].second ;
              t = false ;
            }
          if(t)
            tvar_types += *vi ;
        }
      } else
        tvar_types += (*i).var ;

      target_vars += (*i).var ;
      tvars += (*i).var ;
    }

    time_ident source_time,target_time ;
    //Changes svars to source_vars here.
    for(variableSet::const_iterator i=source_vars.begin();
        i!=source_vars.end();
        ++i) {
      source_time =  source_time.before((*i).get_info().time_id)
        ?(*i).get_info().time_id:source_time ;
    }
    int target_offset = 0 ;
    bool target_asgn = 0;

    for(variableSet::const_iterator i=tvars.begin();i!=tvars.end();++i) {
      if(i==tvars.begin()) {
        target_time = (*i).get_info().time_id ;
        target_offset = (*i).get_info().offset ;
        target_asgn = (*i).get_info().assign ;
      } else
        if(target_time != (*i).get_info().time_id ||
           target_asgn != (*i).get_info().assign ||
           (!target_asgn &&
            (target_offset != (*i).get_info().offset))) {
          cerr << "targets not all at identical time level in rule : "
               << endl ;
          rule_impl->Print(cerr) ;
        }
    }

    output_is_parameter = false ;
    for(variableSet::const_iterator
          i=tvar_types.begin();i!=tvar_types.end();++i) {
      //      cerr << "i=" << *i << endl ;
      if(rule_impl->get_store(*i)==0) {
        cerr << "unable to access variable " << *i << endl ;
        cerr << " error occured in rule " ;
        rule_impl->Print(cerr) ;
      } else if(isPARAMETER(rule_impl->get_store(*i))) {
        if(i!=tvars.begin() && !output_is_parameter) {
          cerr << "can't mix parameters and stores in target" << endl
               << "error occured in rule ";
          rule_impl->Print(cerr) ;
        }
        output_is_parameter = true ;
      } else {
        if(i!=tvars.begin() && output_is_parameter) {
          cerr << "can't mix parameters and stores in target" << endl
               << "error occured in rule ";
          rule_impl->Print(cerr) ;
        }
      }
    }

    source_level = source_time ;
    target_level = target_time ;

    rule_class = TIME_SPECIFIC ;

    if(source_time == target_time) {
      if(source_time == time_ident())
        rule_class = GENERIC ;
    } else if(target_time.before(source_time)) {
      rule_class = COLLAPSE ;
    } else if(source_time.before(target_time)) {
      rule_class = BUILD ;
    } else {
      cerr << "unable to infer time hierarchy from rule :" << endl ;
      rule_impl->Print(cerr) ;
    }

    time_advance = false ;
    if(1 == target_offset)
      time_advance = true ;
    if(time_advance && rule_class == BUILD)
      cerr << "can not advance time in build rule " << rule_impl->get_name()
           << endl;
    if(target_offset > 2)
      cerr << "invalid target offset in rule "
           << rule_impl->get_name() << endl;
  }

  rule::info::info(const info &fi, time_ident tl) {
    if(fi.rule_class == INTERNAL) {
      *this = fi ;
      set<vmap_info>::const_iterator i ;
      set<vmap_info> tmp ;
      for(i=desc.sources.begin();i!=desc.sources.end();++i)
        tmp.insert(convert_vmap_info(*i,tl)) ;
      desc.sources.swap(tmp) ;
      tmp.clear() ;
      for(i=desc.targets.begin();i!=desc.targets.end();++i)
        tmp.insert(convert_vmap_info(*i,tl)) ;
      desc.targets.swap(tmp) ;
      tmp.clear() ;
      for(i=desc.constraints.begin();
          i!=desc.constraints.end();++i)
        tmp.insert(convert_vmap_info(*i,tl)) ;
      desc.constraints.swap(tmp) ;
      tmp.clear() ;
      desc.conditionals = convert_set(desc.conditionals,tl) ;
      rule_impl = new NULL_RULE_IMPL ;
      rule_ident = internal_qualifier + ":" + desc.rule_identifier() ;

      source_vars = variableSet() ;
      target_vars = variableSet() ;
      map_vars = variableSet() ;
      constraint_vars = variableSet() ;

      variableSet svars,tvars ;
      for(i=desc.sources.begin();i!=desc.sources.end();++i) {
        for(size_t j=0;j<(*i).mapping.size();++j) {
          source_vars += (*i).mapping[j] ;
          map_vars += (*i).mapping[j] ;
        }
        source_vars += (*i).var ;
        svars += (*i).var ;
      }
      for(i=desc.constraints.begin();i!=desc.constraints.end();++i) {
        for(size_t j=0;j<(*i).mapping.size();++j) {
          source_vars += (*i).mapping[j] ;
          map_vars += (*i).mapping[j] ;
        }
        constraint_vars += (*i).var ;
        source_vars += (*i).var ;
      }
      source_vars += desc.conditionals ;
      for(i=desc.targets.begin();i!=desc.targets.end();++i) {
        for(size_t j=0;j<(*i).mapping.size();++j) {
          source_vars += (*i).mapping[j] ;
          map_vars += (*i).mapping[j] ;
        }
        target_vars += (*i).var ;
        tvars += (*i).var ;
      }

      time_ident source_time,target_time ;

      for(variableSet::const_iterator i=svars.begin();i!=svars.end();++i) {
        source_time =  source_time.before((*i).get_info().time_id)
          ?(*i).get_info().time_id:source_time ;
      }
      int target_offset = 0 ;
      bool target_asgn = 0 ;

      for(variableSet::const_iterator i=tvars.begin();i!=tvars.end();++i) {
        if(i==tvars.begin()) {
          target_time = (*i).get_info().time_id ;
          target_offset = (*i).get_info().offset ;
          target_asgn = (*i).get_info().assign ;
        } else
          if(rule_class != rule::INTERNAL &&
             (target_time != (*i).get_info().time_id ||
              target_asgn != (*i).get_info().assign ||
              (!target_asgn &&
               (target_offset != (*i).get_info().offset)))) {
            cerr << "targets not all at identical time level in rule : "
                 << endl ;
            rule_impl->Print(cerr) ;
          }
      }

      source_level = source_time ;
      target_level = target_time ;

      time_advance = false ;
      if(1 == target_offset)
        time_advance = true ;

      return ;
    }
    rule_impl = fi.rule_impl->new_rule_impl() ;
    //rule_impl->set_variable_times(tl) ;
    variableSet vset = rule_impl->get_var_list() ;
    std::map<variable, variable> rm ;
    for(variableSet::const_iterator vsi = vset.begin(); vsi != vset.end(); ++vsi) {
      rm[variable(*vsi)] = variable(variable(*vsi), tl) ;
    }
    rule_impl->rename_vars(rm) ;

    warn(fi.rule_class != GENERIC) ;
    source_level = tl ;
    target_level = tl ;
    rule_class = GENERIC ;
    output_is_parameter = fi.output_is_parameter ;
    time_advance = false ;
    desc = rule_impl->get_info() ;
    rule_ident = desc.rule_identifier() ;

    set<vmap_info>::const_iterator i ;
    variableSet svars,tvars ;
    for(i=desc.sources.begin();i!=desc.sources.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      source_vars += (*i).var ;
      svars += (*i).var ;
    }
    for(i=desc.constraints.begin();i!=desc.constraints.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      source_vars += (*i).var ;
      constraint_vars += (*i).var ;
    }
    source_vars += desc.conditionals ;
    for(i=desc.targets.begin();i!=desc.targets.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      target_vars += (*i).var ;
      tvars += (*i).var ;
    }

  }

  // prepend time_ident to the info
  rule::info::info(time_ident tl, const info& fi) {
    if(fi.rule_class == INTERNAL) {
      *this = fi ;
      set<vmap_info>::const_iterator i ;
      set<vmap_info> tmp ;
      for(i=desc.sources.begin();i!=desc.sources.end();++i)
        tmp.insert(prepend_vmap_info(*i,tl)) ;
      desc.sources.swap(tmp) ;
      tmp.clear() ;
      for(i=desc.targets.begin();i!=desc.targets.end();++i)
        tmp.insert(prepend_vmap_info(*i,tl)) ;
      desc.targets.swap(tmp) ;
      tmp.clear() ;
      for(i=desc.constraints.begin();
          i!=desc.constraints.end();++i)
        tmp.insert(prepend_vmap_info(*i,tl)) ;
      desc.constraints.swap(tmp) ;
      tmp.clear() ;
      desc.conditionals = prepend_set(desc.conditionals,tl) ;
      rule_impl = new NULL_RULE_IMPL ;
      rule_ident = internal_qualifier + ":" + desc.rule_identifier() ;

      source_vars = variableSet() ;
      target_vars = variableSet() ;
      map_vars = variableSet() ;
      constraint_vars = variableSet() ;

      variableSet svars,tvars ;
      for(i=desc.sources.begin();i!=desc.sources.end();++i) {
        for(size_t j=0;j<(*i).mapping.size();++j) {
          source_vars += (*i).mapping[j] ;
          map_vars += (*i).mapping[j] ;
        }
        source_vars += (*i).var ;
        svars += (*i).var ;
      }
      for(i=desc.constraints.begin();i!=desc.constraints.end();++i) {
        for(size_t j=0;j<(*i).mapping.size();++j) {
          source_vars += (*i).mapping[j] ;
          map_vars += (*i).mapping[j] ;
        }
        constraint_vars += (*i).var ;
        source_vars += (*i).var ;
      }
      source_vars += desc.conditionals ;
      for(i=desc.targets.begin();i!=desc.targets.end();++i) {
        for(size_t j=0;j<(*i).mapping.size();++j) {
          source_vars += (*i).mapping[j] ;
          map_vars += (*i).mapping[j] ;
        }
        target_vars += (*i).var ;
        tvars += (*i).var ;
      }

      time_ident source_time,target_time ;

      for(variableSet::const_iterator i=svars.begin();i!=svars.end();++i) {
        source_time =  source_time.before((*i).get_info().time_id)
          ?(*i).get_info().time_id:source_time ;
      }
      int target_offset = 0 ;
      bool target_asgn = 0 ;

      for(variableSet::const_iterator i=tvars.begin();i!=tvars.end();++i) {
        if(i==tvars.begin()) {
          target_time = (*i).get_info().time_id ;
          target_offset = (*i).get_info().offset ;
          target_asgn = (*i).get_info().assign ;
        } else
          if(rule_class != rule::INTERNAL &&
             (target_time != (*i).get_info().time_id ||
              target_asgn != (*i).get_info().assign ||
              (!target_asgn &&
               (target_offset != (*i).get_info().offset)))) {
            cerr << "targets not all at identical time level in rule : "
                 << endl ;
            rule_impl->Print(cerr) ;
          }
      }

      source_level = source_time ;
      target_level = target_time ;

      time_advance = false ;
      if(1 == target_offset)
        time_advance = true ;

      return ;
    }
    rule_impl = fi.rule_impl->new_rule_impl() ;
    //rule_impl->set_variable_times(tl) ;
    variableSet vset = rule_impl->get_var_list() ;
    std::map<variable, variable> rm ;
    for(variableSet::const_iterator vsi = vset.begin();
        vsi != vset.end(); ++vsi) {
      rm[variable(*vsi)] = variable(tl,variable(*vsi)) ;
    }

    rule_impl->rename_vars(rm) ;

    //warn(fi.rule_class != GENERIC) ;
    //source_level = prepend_time(tl,source_level) ;
    //target_level = prepend_time(tl,target_level) ;
    //rule_class = GENERIC ;
    //rule_class = TIME_SPECIFIC ;
    rule_class = fi.type() ;
    output_is_parameter = fi.output_is_parameter ;
    //time_advance = false ;
    desc = rule_impl->get_info() ;
    rule_ident = desc.rule_identifier() ;

    set<vmap_info>::const_iterator i ;
    variableSet svars,tvars ;
    for(i=desc.sources.begin();i!=desc.sources.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      source_vars += (*i).var ;
      svars += (*i).var ;
    }
    for(i=desc.constraints.begin();i!=desc.constraints.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      source_vars += (*i).var ;
      constraint_vars += (*i).var ;
    }
    source_vars += desc.conditionals ;
    for(i=desc.targets.begin();i!=desc.targets.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      target_vars += (*i).var ;
      tvars += (*i).var ;
    }

    /////// experiment code
    time_ident source_time, target_time ;
    for(variableSet::const_iterator i=svars.begin();i!=svars.end();++i) {
      source_time =  source_time.before((*i).get_info().time_id)
        ?(*i).get_info().time_id:source_time ;
    }
    int target_offset = 0 ;
    bool target_asgn = 0 ;

    for(variableSet::const_iterator i=tvars.begin();i!=tvars.end();++i) {
      if(i==tvars.begin()) {
        target_time = (*i).get_info().time_id ;
        target_offset = (*i).get_info().offset ;
        target_asgn = (*i).get_info().assign ;
      } else
        if(rule_class != rule::INTERNAL &&
           (target_time != (*i).get_info().time_id ||
            target_asgn != (*i).get_info().assign ||
            (!target_asgn &&
             (target_offset != (*i).get_info().offset)))) {
          cerr << "targets not all at identical time level in rule : "
               << endl ;
          rule_impl->Print(cerr) ;
        }
    }
    source_level = source_time ;
    target_level = target_time ;

    time_advance = false ;
    if(1 == target_offset)
      time_advance = true ;

  }

  // rename function that renames variables in the rule
  // according to the rename map passed in. It is the
  // general interface for rule promotion.
  rule rule::rename_vars(std::map<variable,variable>& rvm) const {
    if(type() != rule::INTERNAL) {
      rule_implP rp = get_rule_implP() ;
      rp->rename_vars(rvm) ;
      return rule(rp) ;
    }else {
      rule::info newinfo = get_info() ;

      std::set<vmap_info>::const_iterator i ;
      std::set<vmap_info> tmp ;
      for(i = newinfo.desc.sources.begin();
          i != newinfo.desc.sources.end(); ++i)
        tmp.insert(rename_vmap_info(*i, rvm)) ;
      newinfo.desc.sources.swap(tmp) ;
      tmp.clear() ;
      for(i = newinfo.desc.targets.begin();
          i != newinfo.desc.targets.end(); ++i)
        tmp.insert(rename_vmap_info(*i, rvm)) ;
      newinfo.desc.targets.swap(tmp) ;
      tmp.clear() ;
      for(i=newinfo.desc.constraints.begin();
          i!=newinfo.desc.constraints.end();++i)
        tmp.insert(rename_vmap_info(*i, rvm)) ;
      newinfo.desc.constraints.swap(tmp) ;
      tmp.clear() ;
      newinfo.desc.conditionals = rename_set(newinfo.desc.conditionals,rvm) ;

      newinfo.rule_impl = new NULL_RULE_IMPL ;
      newinfo.rule_ident = newinfo.internal_qualifier + ":"
        + newinfo.desc.rule_identifier() ;

      newinfo.source_vars = variableSet() ;
      newinfo.target_vars = variableSet() ;
      newinfo.map_vars = variableSet() ;
      newinfo.constraint_vars = variableSet() ;

      variableSet svars,tvars ;
      for(i=newinfo.desc.sources.begin();
          i!=newinfo.desc.sources.end();++i) {
        for(size_t j=0;j<(*i).mapping.size();++j) {
          newinfo.source_vars += (*i).mapping[j] ;
          newinfo.map_vars += (*i).mapping[j] ;
        }
        newinfo.source_vars += (*i).var ;
        svars += (*i).var ;
      }
      for(i=newinfo.desc.constraints.begin();
          i!=newinfo.desc.constraints.end();++i) {
        for(size_t j=0;j<(*i).mapping.size();++j) {
          newinfo.source_vars += (*i).mapping[j] ;
          newinfo.map_vars += (*i).mapping[j] ;
        }
        newinfo.constraint_vars += (*i).var ;
        newinfo.source_vars += (*i).var ;
      }
      newinfo.source_vars += newinfo.desc.conditionals ;
      for(i=newinfo.desc.targets.begin();
          i!=newinfo.desc.targets.end();++i) {
        for(size_t j=0;j<(*i).mapping.size();++j) {
          newinfo.source_vars += (*i).mapping[j] ;
          newinfo.map_vars += (*i).mapping[j] ;
        }
        newinfo.target_vars += (*i).var ;
        tvars += (*i).var ;
      }

      time_ident source_time,target_time ;

      for(variableSet::const_iterator i=svars.begin();i!=svars.end();++i) {
        source_time =  source_time.before((*i).get_info().time_id)
          ?(*i).get_info().time_id:source_time ;
      }
      int target_offset = 0 ;
      bool target_asgn = 0 ;

      for(variableSet::const_iterator i=tvars.begin();i!=tvars.end();++i) {
        if(i==tvars.begin()) {
          target_time = (*i).get_info().time_id ;
          target_offset = (*i).get_info().offset ;
          target_asgn = (*i).get_info().assign ;
        } else
          if(newinfo.rule_class != rule::INTERNAL &&
             (target_time != (*i).get_info().time_id ||
              target_asgn != (*i).get_info().assign ||
              (!target_asgn &&
               (target_offset != (*i).get_info().offset)))) {
            cerr << "targets not all at identical time level in rule : "
                 << endl ;
            newinfo.rule_impl->Print(cerr) ;
          }
      }

      newinfo.source_level = source_time ;
      newinfo.target_level = target_time ;

      newinfo.time_advance = false ;
      if(1 == target_offset)
        newinfo.time_advance = true ;

      return rule(newinfo) ;
    }
  }

  rule_implP rule::info::get_rule_implP() const {
    rule_implP fp = rule_impl->new_rule_impl() ;
    return fp ;
  }

  void
  rule::rename(const std::string& s) {
    // get the rule info.
    rule::info info = rdb->fiv[-(id+1)] ;
    // save the old name
    std::string old_name = info.name() ;
    // replace the "rule_ident" inside ;
    info.rule_ident = s ;
    // Get new id
    id = rdb->get_id(info) ;
  }

  rule rule::get_rule_by_name(std::string &name) {
    int myId = rdb->query_name(name);
    rule myRule(myId);
    return myRule;

  }

  namespace {
    // utility function
    inline variableSet get_rule_var_list(const rule& r) {
      variableSet vset ;
      if(r.type() == rule::INTERNAL) {
        vset += r.sources() ;
        vset += r.targets() ;
      }else{
        vset = r.get_rule_implP()->get_var_list() ;
      }

      return vset ;
    }
  } // end of unnamed namespace

  // rule promotion
  rule promote_rule(const rule& r, const time_ident& t) {
    std::map<variable,variable> rvm ;
    variableSet vset = get_rule_var_list(r) ;

    for(variableSet::const_iterator vi=vset.begin();
        vi!=vset.end(); ++vi) {
      rvm[*vi] = variable(*vi,t) ;
    }
    return r.rename_vars(rvm) ;
  }

  // prepend a time to the rule
  rule prepend_rule(const rule& r, const time_ident& t) {
    std::map<variable,variable> rvm ;
    variableSet vset = get_rule_var_list(r) ;

    for(variableSet::const_iterator vi=vset.begin();
        vi!=vset.end(); ++vi) {
      rvm[*vi] = variable(t,*vi) ;
    }
    return r.rename_vars(rvm) ;
  }

  rule_implP rule_impl::new_rule_impl() const {
    cerr << "rule_impl::new_rule_impl() was called:" << endl
         << "this rule should never be called, use copy_rule_impl<> to"<< endl
         << "generate a copyable rule" << endl ;
    abort() ;
    return new NULL_RULE_IMPL ;
  }

  void rule_impl::rename_vars(std::map<variable,variable> &rvm) {
    cerr << "rule_impl::rename_vars(std::map<variable,variable> &rv) was called:" << endl << "this rule should never be called, use copy_rule_impl<> to"<< endl
         << "rename the variables" << endl ;
    abort() ;
  }



  rule::info::info(const string &s) {
    rule_class = rule::INTERNAL ;
    exprP p = expression::create(s) ;
    exprList l = collect_associative_op(p,OP_COMMA) ;
    exprList::const_iterator li ;
    for(li=l.begin();li!=l.end();++li) {
      const exprP &f = *li ;
      if(f->op == OP_FUNC) {
        if(f->name == "source") {
          fill_descriptors(desc.sources,f->expr_list) ;
        } else if(f->name == "target") {
          fill_descriptors(desc.targets,f->expr_list) ;
        } else if(f->name == "constraint") {
          fill_descriptors(desc.constraints,f->expr_list) ;
        } else if(f->name == "conditional") {
          desc.conditionals = variableSet(f->expr_list.front()) ;
        } else if(f->name == "qualifier") {
          internal_qualifier = f->expr_list.front()->name ;
        } else {
          cerr << "unable to interpret internal rule representation"
               << endl ;
          cerr << "rule was given " << s << endl ;
          Loci::Abort() ;
        }
      } else {
        cerr << "syntax error parsing internal rule representation"
             << endl ;
        cerr << "rule was given " << s << endl ;
        Loci::Abort() ;
      }


    }

    rule_impl = new NULL_RULE_IMPL ;

    rule_ident = internal_qualifier + ":" + desc.rule_identifier() ;

    set<vmap_info>::const_iterator i ;
    variableSet svars,tvars ;
    for(i=desc.sources.begin();i!=desc.sources.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      source_vars += (*i).var ;
      svars += (*i).var ;
    }
    for(i=desc.constraints.begin();i!=desc.constraints.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      source_vars += (*i).var ;
      constraint_vars += (*i).var ;
    }
    source_vars += desc.conditionals ;
    for(i=desc.targets.begin();i!=desc.targets.end();++i) {
      for(size_t j=0;j<(*i).mapping.size();++j) {
        source_vars += (*i).mapping[j] ;
        map_vars += (*i).mapping[j] ;
      }
      target_vars += (*i).var ;
      tvars += (*i).var ;
    }

    time_ident source_time,target_time ;

    for(variableSet::const_iterator i=svars.begin();i!=svars.end();++i) {
      source_time =  source_time.before((*i).get_info().time_id)
        ?(*i).get_info().time_id:source_time ;
    }
    int target_offset = 0 ;
    bool target_asgn = 0 ;

    for(variableSet::const_iterator i=tvars.begin();i!=tvars.end();++i) {
      if(i==tvars.begin()) {
        target_time = (*i).get_info().time_id ;
        target_offset = (*i).get_info().offset ;
        target_asgn = (*i).get_info().assign ;
      } else
        if(rule_class != rule::INTERNAL &&
           (target_time != (*i).get_info().time_id ||
            target_asgn != (*i).get_info().assign ||
            (!target_asgn &&
             (target_offset != (*i).get_info().offset)))) {
          cerr << "targets not all at identical time level in rule : "
               << endl ;
          rule_impl->Print(cerr) ;
        }
    }

    source_level = source_time ;
    target_level = target_time ;

    time_advance = false ;
    if(1 == target_offset)
      time_advance = true ;

    if(target_offset > 2)
      cerr << "invalid target offset in rule "
           << rule_impl->get_name() << endl;

  }

  ostream &ruleSet::Print(ostream &s) const {
    for(ruleSet::const_iterator i=begin();i!=end();++i)
      s << *i << endl ;
    return s;
  }
  //Definition of global rule lists
  register_rule_impl_list register_rule_list ;
  rule_impl_list global_rule_list ;

  std::string register_module::load_module() const { return string("") ; }
  std::string register_module::using_nspace() const { return string("") ; }
  std::string register_module::load_nspace() const { return string("") ; }
  std::string register_module::disable_compute_vars() const { return string("") ; }

  rule_impl_list::~rule_impl_list() {
    rule_list_ent *p,*v ;
    for(p=list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
  }
  void rule_impl_list::clear() {
    rule_list_ent *p,*v ;
    for(p=list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
    list = 0 ;
  }
  void rule_impl_list::push_rule(register_rule_type *p) {
    rule_list_ent *flp = new rule_list_ent(p,list) ;
    list = flp ;
  }

  void rule_impl_list::copy_rule_list(const rule_impl_list& rl) {
    rule_list_ent *p, *v ;
    for(p = rl.list; p != 0; p=v) {
      push_rule(p->rr) ;
      v = p->next ;
    }
  }
  void rule_impl_list::copy_rule_list(const register_rule_impl_list& rl) {
    rule_list_ent *p, *v ;
    for(p = rl.global_list; p != 0; p=v) {
      push_rule(p->rr) ;
      v = p->next ;
    }
  }

  //Declaration of static variable global_list
  rule_impl_list::rule_list_ent *register_rule_impl_list::global_list = 0 ;

  register_rule_impl_list::~register_rule_impl_list() {
    rule_list_ent *p,*v ;
    for(p=global_list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
  }
  void register_rule_impl_list::clear() {
    rule_list_ent *p,*v ;
    for(p=global_list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
    global_list = 0 ;
  }
  bool register_rule_impl_list::empty() {
    return (global_list == 0) ;
  }
  void register_rule_impl_list::push_rule(register_rule_type *p) {
    rule_list_ent *flp = new rule_list_ent(p,global_list) ;
    global_list = flp ;
  }

  const ruleSet rule_db::EMPTY_RULE ;

  void rule_db::add_rule(const rule_implP &fp) {
    // Check for rule consistency
    fp->check_perm_bits() ;
    add_rule(rule(fp)) ;
  }

  void rule_db::add_rule(rule f) {
    string fname = f.get_info().rule_impl->get_name() ;
    rule_map_type::const_iterator fmti = name2rule.find(fname) ;
    if(fmti != name2rule.end()) {
      fname = f.get_info().name() ;
    }
    name2rule[fname] = f ;
    // default and optional rules need to be in a different ruleSet
    if(f.get_info().rule_impl->get_rule_class() == rule_impl::DEFAULT) {
      if(default_rules.inSet(f)) {
        cerr << "Warning, adding duplicate rule to rule database"
             << endl
             << " Rule = " << f << endl ;
      }else
        default_rules += f ;
      return ;
    }
    if(f.get_info().rule_impl->get_rule_class() == rule_impl::OPTIONAL) {
      if(optional_rules.inSet(f)) {
        cerr << "Warning, adding duplicate rule to rule database"
             << endl
             << " Rule = " << f << endl ;
      }else
        optional_rules += f ;
      return ;
    }

    if(known_rules.inSet(f)) {
      cerr << "Warning, adding duplicate rule to rule database"
           << endl
           << " Rule = " << f << endl ;
    } else {
      // Now link all rule sources
      variableSet svars = f.sources() ;
      variableSet tvars = f.targets() ;
      for(variableSet::const_iterator i=svars.begin();i!=svars.end();++i)
        srcs2rule[*i] += f ;
      for(variableSet::const_iterator i=tvars.begin();i!=tvars.end();++i) {
        variable v = *i ;
        while(v.get_info().priority.size() != 0) {
          trgt2rule[v] += f ;
          v = v.drop_priority() ;
        }
        trgt2rule[v] += f ;
      }
      if(svars == EMPTY)
        cerr << "WARNING, rule " << f << " has no sources" << endl ;
      if(tvars == EMPTY)
        cerr << "WARNING, rule " << f << " has no targets" << endl ;
      // now add rules to corresponding keyspace partition
      keyspace2rule[f.get_info().rule_impl->get_keyspace_tag()] += f ;
      known_rules += f ;
    }
  }

  void rule_db::add_rules(rule_impl_list &gfl) {
    for(rule_impl_list::iterator i=gfl.begin();i!=gfl.end();++i)
      if(!(i.get_p())->rr->is_module_rule())
	add_rule(*i) ;
  }
  void rule_db::add_rules(register_rule_impl_list &gfl) {
    for(rule_impl_list::iterator i=gfl.begin();i!=gfl.end();++i)
      if(!(i.get_p())->rr->is_module_rule())
	add_rule(*i) ;
  }

  void rule_db::remove_rule(rule f) {
    string fname = f.get_info().rule_impl->get_name() ;
    rule_map_type::const_iterator fmti = name2rule.find(fname) ;
    if(fmti != name2rule.end()) {
      fname = f.get_info().name() ;
    }
    name2rule.erase(fname) ;
    // remove rules in the appropriate categories
    if(f.get_info().rule_impl->get_rule_class() == rule_impl::DEFAULT) {
      default_rules -= f ;
      return ;
    }
    if(f.get_info().rule_impl->get_rule_class() == rule_impl::OPTIONAL) {
      optional_rules -= f ;
      return ;
    }

    known_rules -= f ;
    variableSet svars = f.sources() ;
    variableSet tvars = f.targets() ;
    for(variableSet::const_iterator i=svars.begin();i!=svars.end();++i)
      srcs2rule[*i] -= f ;
    for(variableSet::const_iterator i=tvars.begin();i!=tvars.end();++i) {
      variable v = *i ;
      while(v.get_info().priority.size() != 0) {
        trgt2rule[v] += f ;
        v = v.drop_priority() ;
      }
      trgt2rule[v] -= f ;
    }
    keyspace2rule[f.get_info().rule_impl->get_keyspace_tag()] -= f ;
  }

  void rule_db::remove_rules(const ruleSet& rs) {
    for(ruleSet::const_iterator ri=rs.begin();ri!=rs.end();++ri)
      remove_rule(*ri) ;
  }

}


