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
#ifndef RULE_H
#define RULE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/debug.h>
#include <entitySet.h>
#include <store_rep.h>
#include <key_manager.h>

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <iterator>

#include <Tools/cptr.h>
#include <Tools/expr.h>
#include <Tools/variable.h>
#include <Map.h>
#include <parameter.h>
#include <blackbox.h>
#include <storeVec.h>
#include <storeMat.h>
#include <multiStore.h>
#include <DMultiStore.h>

namespace Loci {
  
  class fact_db ;
  class sched_db ;
  
  class joiner : public CPTR_type {
  public:
    virtual CPTR<joiner> clone() = 0 ;
    virtual storeRepP getTargetRep() = 0 ;
    virtual void SetArgs(storeRepP &target, storeRepP &source) = 0 ;
    virtual void Join(const sequence &seq) = 0 ;
    virtual void Join(Map &t2s, const sequence &seq) = 0 ;
  } ;  
  
  class rule;
  class rule_impl : public CPTR_type {
  public:
    struct info {
      std::set<vmap_info> sources, targets, constraints ;
      std::set<vmap_info> dynamic_constraints ;
      variableSet conditionals ;
      std::string rule_identifier() const ;
      variableSet input_vars() const ;
      variableSet output_vars() const ;
    } ;
    typedef CPTR<rule_impl> rule_implP ;
    enum rule_impl_type {POINTWISE,SINGLETON,UNIT,APPLY,DEFAULT,OPTIONAL,CONSTRAINT_RULE,MAP_RULE,BLACKBOX_RULE,INSERTION,DELETION,ERASE,SUPER_RULE,UNKNOWN} ;
  private:
    rule_impl_type rule_impl_class ;
    bool rule_threading ;
    bool gpgpu_kernel ;
    bool use_dynamic_schedule ;
    bool relaxed_recursion ;
    bool specialized_parametric ;
    bool use_parametric_variable ;
    variable ParametricVariable ;
    mutable std::string name ;
    info rule_info ;
    typedef std::multimap<variable, store_instance *> storeIMap ;
    storeIMap var_table ;
    std::map<variable,variable> rvmap ;
    std::map<variable,int> varInfoId ;
    const char **vardoc ;

    void source(const std::string &invar) ;
    void target(const std::string &outvar) ;
    std::string rule_comments ; // the comments for a rule_impl
    std::string fileloc ;
    // keyspace tag that dictates which keyspace the rule
    // is currently assigned to, defaults to empty string,
    // which currently means that the rule should be scheduled
    // under the default "metis" key space.
    std::string space_tag ;
    // bit indicates whether keyspace (the one this rule is in)
    // distribution should be considered after this rule's execution
    bool space_dist ;
    // these flags record whether the rule has pre- and postlude methods
    // defaults would be false for both
    bool use_prelude;
    bool use_postlude;
  protected:
    rule_impl(rule_impl &f) { fatal(true) ; }
    void rule_class(rule_impl_type ft) { rule_impl_class = ft ; }
    void disable_threading() { rule_threading = false ; }
    void gpgpu() { gpgpu_kernel = true ; }
    void enable_dynamic_scheduling() { use_dynamic_schedule = true ; }
    void load_balance() { use_dynamic_schedule = true ; }
    void set_relaxed_recursion() { relaxed_recursion = true ; }
    void set_specialized() { specialized_parametric = true ; }
    void set_parametric_variable(const std::string name) {
      use_parametric_variable = true ;
      ParametricVariable = variable(name) ;
    }
    // these should be called if the pre- and postlude methods are present
    void enable_prelude() { use_prelude = true; }
    void enable_postlude() { use_postlude = true; }
    void rule_name(const std::string &name) ;
    void name_store(const std::string &name,store_instance &si) ;
    void input(const std::string &invar) { source(invar) ; }
    void output(const std::string &outvar) { target(outvar) ; }
    void constraint(const std::string &constrain) ;
    void conditional(const std::string &cond) ;
    void comments(const char* c) {rule_comments += c ;}
    void set_file(const char* c) {fileloc = c ; }
    void store_info_id(const std::string &var, int id) {
      variable v(var) ;
      varInfoId[v] = id ;
    }
    void setvardoc(const char **docarray) { vardoc = docarray ; }
    // set the keyspace tag
    void keyspace_tag(const std::string& t) {space_tag = t ;}
    // set the space_dist bit
    void keyspace_dist_hint() {space_dist = true ;}
  public:
    const char *getvardoc(variable v) const {
      if(vardoc) {
        auto mi = varInfoId.find(v) ;
        if(mi==varInfoId.end())
          return "\0\0\0\0" ;
        else
          return vardoc[mi->second] ;
      }  else
        return "\0\0\0\0" ;
    }
    rule_impl() ;
    bool check_perm_bits() const ;
    bool thread_rule() const { return rule_threading; }
    bool is_gpu_kernel() const { return gpgpu_kernel ; }
    bool dynamic_schedule_rule() const { return use_dynamic_schedule; }
    bool is_relaxed() const { return relaxed_recursion ; }
    bool is_specialized() { return specialized_parametric; }
    bool is_parametric_provided() { return use_parametric_variable ; }
    // check if pre- and postlude methods are present
    bool has_prelude() const { return use_prelude; }
    bool has_postlude() const { return use_postlude; }
    variable get_parametric_variable() { return ParametricVariable ; }
    void initialize(fact_db &facts) ;
    // this method returns all the keyspaces involved in a rule
    // that need to be considered for distribution
    std::vector<std::string>
    gather_keyspace_dist() const ;
    // this method returns whether this rule affects keyspace
    // distribution
    bool
    affect_keyspace_dist() const {return space_dist ;}

    std::string
    get_keyspace_tag() const {return space_tag ;}
    
    std::string get_name() const ;
    rule_impl_type get_rule_class() const { return rule_impl_class ; }
    const info &get_info() const { return rule_info ; }
    void set_store(variable v, const storeRepP &p) ;
    void set_store(const std::string &nm, const storeRepP &p) 
      { set_store(variable(expression::create(nm)),p) ; }
    
    storeRepP get_store(variable v) const ;
    storeRepP get_store(const std::string &nm) const
      { return get_store(variable(expression::create(nm))) ; }
    
    // This function checks to
    // see if the rule contains constraints that are of type Map. In case
    // of yes, it creates identical constraints
    // (whose value equal the Map domain) in the fact_db and substitutes
    // the Map constraints in the rule with the real constraints.
    //
    // The motivation of this function is that in the parallel code,
    // a Map constraint will often not be expanded enough to include
    // the clone region, thus causing problems in the rule execution
    // schedule. Since the context of the rule including the clone region
    // will often exceed the domain of the constraint, this will either
    // cause the clone region not being computed properly, or in the case
    // of unit/apply rule, cause conflicts in the existential analysis.
    //
    // We could expand the Maps used in rule contraints to include the
    // clone region. However doing so would often require duplicating
    // the Map on all processes, thus incurring a memory cost, or else,
    // we would be allocating the Map on domains that do not have
    // meaningful values for the Map.
    //
    // Thus, here we are doing a substitution to replace all Maps
    // in rule constraint as real constraint variable. If substitution
    // have happened, "facts" may include newly created constraints;
    // the rules may have its "vmap_info" structure modified
    // to reflect the substitution of constraints for maps.
    void replace_map_constraints(fact_db& facts) ;

    void split_constraints(const variableSet& dc) ;
    
    void set_variable_times(time_ident tl) ;
    void copy_store_from(rule_impl &f) ;
    void Print(std::ostream &s) const ;
   
    void prot_rename_vars(std::map<variable, variable> &rvm) ;
    virtual void rename_vars(std::map<variable, variable>  &rvm) ;
    variableSet get_var_list() ;

    virtual rule_implP new_rule_impl() const ;
    // the compute method handles the kernel computation
    virtual void compute(const sequence &) = 0 ;
    // the prelude method is intended to be used to setup
    // global parameters (such as the vec size of a store)
    // it is not required to define such method when creating
    // a rule and the preprocessor might optimize away this
    // method if a user does not define it in a rule
    virtual void prelude(const sequence&) {}
    // the postlude is similarly defined as the prelude method
    virtual void postlude(const sequence&) {}
    virtual CPTR<joiner> get_joiner() = 0 ;
    virtual rule_implP add_namespace(const std::string& n) const ;
    std::string get_comments() const {return rule_comments ;}
    std::string get_fileloc() const { return fileloc; }
  } ;
  
  typedef rule_impl::rule_implP rule_implP ;
  
  
  template <class TCopyRuleImpl> class copy_rule_impl : public TCopyRuleImpl {
    typedef std::list<std::map<variable, variable> > rename_varList ;
    typedef std::list<std::map<variable, variable> >::const_iterator list_iter;
    rename_varList rvlist ;
  public:
    virtual rule_implP new_rule_impl() const ;
    virtual void rename_vars(std::map<variable, variable> &rvm) ;
  } ; 
  
  template <class TCopyRuleImpl> rule_implP copy_rule_impl<TCopyRuleImpl>::new_rule_impl() const {
    rule_implP realrule_impl = new copy_rule_impl<TCopyRuleImpl> ;
    for(list_iter li = rvlist.begin(); li != rvlist.end(); ++li) { 
      std::map<variable, variable> rvm = *li;
      realrule_impl->rename_vars(rvm) ;
    }
    return realrule_impl ;
  }

  template <class TCopyRuleImpl> 
    void copy_rule_impl<TCopyRuleImpl>::rename_vars(std::map<variable,variable> &rvm) {
    rvlist.push_back(rvm) ;
    rule_impl::prot_rename_vars(rvm) ;
  }


  // this is the new rule type for setting default
  // values in the facts database. It should not have
  // any inputs
  class default_rule: public rule_impl {
  protected:
    default_rule() { rule_class(DEFAULT) ; }
    void name_store(const std::string &nm, store_instance &si)
      { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar) {
      std::cerr << "Warning: a DEFAULT rule should not have any inputs!"
                << endl ;
    }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
      { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  } ;

  // this is the new rule type for setting optional
  // values in the facts database. It should not have
  // any inputs
  class optional_rule: public rule_impl {
  protected:
    optional_rule() { rule_class(OPTIONAL) ; }
    void name_store(const std::string &nm, store_instance &si)
      { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar) {
      std::cerr << "Warning: an OPTIONAL rule should not have any inputs!"
                << endl ;
    }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
      { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  } ;

  class constraint_rule: public rule_impl {
  protected:
    constraint_rule() { rule_class(CONSTRAINT_RULE) ; }
    void name_store(const std::string &nm, store_instance &si)
      { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    // do we allow constraint in a constraint rule???
    // I don't think so currently --- so we disable it for now.
//     void constraint(const std::string &constrain)
//     { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
      { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  } ;

  class map_rule: public rule_impl {
  protected:
    map_rule() { rule_class(MAP_RULE) ; }
    void name_store(const std::string& nm, store_instance& si)
      {rule_impl::name_store(nm,si) ; }
    void input(const std::string& invar)
      {rule_impl::input(invar) ;}
    void output(const std::string& outvar)
      {rule_impl::output(outvar) ;}
    void constraint(const std::string& constrain)
      {rule_impl::constraint(constrain) ;}
    void conditional(const std::string& cond)
      {rule_impl::conditional(cond) ;}
    virtual CPTR<joiner> get_joiner() {return CPTR<joiner>(0) ;}
  } ;
  
  class blackbox_rule : public rule_impl {
  protected:
    blackbox_rule() { rule_class(BLACKBOX_RULE) ; disable_threading(); }
    void name_store(const std::string &nm, store_instance &si)
    { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  } ;

#ifdef DYNAMICSCHEDULING
  // this one is purely for interface purpose
  class insertion_rule_interface: public rule_impl {
  public:
    virtual void set_key_manager(KeyManagerP kp) = 0 ;
    virtual KeyManagerP get_key_manager() const = 0 ;
    virtual const KeySet& get_keys_inserted() const = 0 ;
  } ;
  typedef CPTR<insertion_rule_interface> insertion_rule_interfaceP ;

  // we require that "SequentialContainer" provides an iterator
  // interface, and a "value_type" definition. Good examples of
  // "SequentialContainer" are "std::vector", "std::list" etc.
  template<class SequentialContainer>
  class insertion_rule: public insertion_rule_interface {
    KeyManagerP key_manager ;
    KeySet keys_inserted ;

    int input_num ;
  protected:
    insertion_rule():key_manager(0),keys_inserted(EMPTY),input_num(0)
    {rule_class(INSERTION) ;}

    void name_store(const std::string &nm, store_instance &si)
    { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar) {
      ++input_num ;
      if(input_num >= 2) {
        std::cerr << "Warning: an INSERTION rule should only have ONE input!"
                  << endl ;
      } else
        rule_impl::input(invar) ;
    }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain) {
      std::cerr << "Warning: an INSERTION rule should not have "
                << "any constraints!" << endl ;
    }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }

    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }

    // method to be defined by users
    // supplied by two arguments: the Entity e is the newly
    // created entity for the value passed in, users would
    // just need to decide how to fill the value to their
    // target facts
    virtual void
    insert(const typename SequentialContainer::value_type& value,
           Entity e) = 0 ;
  public:
    virtual void set_key_manager(KeyManagerP kp) {
      key_manager = kp ;
    }
    virtual KeyManagerP get_key_manager() const {
      return key_manager ;
    }
    virtual const KeySet& get_keys_inserted() const {
      return keys_inserted ;
    }
    
    virtual void compute(const sequence& seq) {
      // at the beginning, we'll need to clear the keys inserted
      keys_inserted = EMPTY ;
      // we assume that there are only
      // one input to this rule and its type
      // is blackbox<SequentialContainer>
      
      // lets get the input's rep first
      variable rule_input = *( (get_info().sources.begin())->var.begin()) ;
      storeRepP rule_input_rep = get_store(rule_input) ;
      const_blackbox<SequentialContainer> input(rule_input_rep) ;
      const SequentialContainer& input_sequence = *input ;
      // now we can know how many keys to be created
      size_t num_of_keys = std::distance(input_sequence.begin(),
                                         input_sequence.end()) ;
      keys_inserted = key_manager->generate_key(num_of_keys) ;
      KeySet::const_iterator ki = keys_inserted.begin() ;
      typename SequentialContainer::const_iterator b, e ;
      for(b=input_sequence.begin(),e=input_sequence.end();b!=e;++b,++ki)
        insert(*b, *ki) ;
    }
  } ;

  class deletion_rule: public rule_impl {
    KeyManagerP key_manager ;
    KeySet keys_deleted ;
    bool if_destroy_keys ;
  protected:
    deletion_rule():key_manager(0),
                    keys_deleted(EMPTY),if_destroy_keys(false) {
      rule_class(DELETION) ;
    }

    void name_store(const std::string &nm, store_instance &si)
    { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }

    // this one tells if the deleted keys will also be destroyed
    void destroy_keys() {if_destroy_keys = true ;}

    void delete_key(Key k) {keys_deleted += k ;}

    // defined by the users
    virtual void evaluate_key(Key k) = 0 ;
  public:
    virtual void set_key_manager(KeyManagerP kp) {
      key_manager = kp ;
    }
    virtual KeyManagerP get_key_manager() const {
      return key_manager ;
    }
    const KeySet& get_keys_deleted() const {return keys_deleted ;}
    bool destroy_deleted_keys() const {return if_destroy_keys ;}

    void compute(const sequence& seq) {
      // first clear the keys to be destroyed
      keys_deleted = EMPTY ;
      do_loop(seq,this,&deletion_rule::evaluate_key) ;
    }
  } ;
  typedef CPTR<deletion_rule> deletion_ruleP ;

  class erase_rule: public rule_impl {
    KeySet record_erased ;
  protected:
    erase_rule():record_erased(EMPTY) {
      rule_class(ERASE) ;
    }

    void name_store(const std::string &nm, store_instance &si)
    { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
    void erase_record(Key k) {record_erased += k ;}
    // defined by the users
    virtual void evaluate_record(Key k) = 0 ;
  public:
    const KeySet& get_erased_record() const {return record_erased ;}
    void compute(const sequence& seq) {
      // first clear the keys to be destroyed
      record_erased = EMPTY ;
      do_loop(seq,this,&erase_rule::evaluate_record) ;
    }
  } ;
  typedef CPTR<erase_rule> erase_ruleP ;
#endif
  
  // This rule is 
  class super_rule : public rule_impl {
  protected:
    super_rule() { rule_class(SUPER_RULE) ; disable_threading(); }
    void name_store(const std::string &nm, store_instance &si)
    { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  public:
    virtual void process_existential(rule r, fact_db &facts, sched_db &scheds) = 0 ;
    virtual void process_requests(rule r, fact_db &facts, sched_db &scheds) = 0 ;
  } ;

  class pointwise_rule : public rule_impl {
  protected:
    pointwise_rule() { rule_class(POINTWISE) ; }
    void name_store(const std::string &nm, store_instance &si)
      { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
      { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  } ;

  class singleton_rule : public rule_impl {
  protected:
    singleton_rule() { rule_class(SINGLETON) ; }
    void name_store(const std::string &nm, store_instance &si)
    { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  } ;
  
  class unit_rule : public rule_impl {
   protected:
    unit_rule() { rule_class(UNIT) ; }
    void name_store(const std::string &nm, store_instance &si)
    { rule_impl::name_store(nm,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  } ;

  template <class T, class Op> class joinOp : public joiner {
    Op join ;
    T t,s ;
  public:
    virtual CPTR<joiner> clone() ;
    virtual storeRepP getTargetRep() ;
    virtual void SetArgs(storeRepP &target, storeRepP &source) ;
    virtual void Join(const sequence &seq) ;
    virtual void Join(Map &t2s, const sequence &seq)  ;
  } ;

  template <class T, class Op> CPTR<joiner> joinOp<T,Op>::clone() {
    return CPTR<joiner>(new joinOp<T,Op> );
  }

  template <class T, class Op>  storeRepP joinOp<T,Op>::getTargetRep() {
    T st ;
    return st.Rep();
  }
  
  template <class T, class Op>
    void joinOp<T,Op>::SetArgs(storeRepP &target, storeRepP &source)
  { s.setRep(source) ; t.setRep(target) ; }
  
  template <class T, class Op> void joinOp<T,Op>::Join(const sequence &seq) {
    for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) 
      join(t[*i],s[*i]) ;
  }
  
  
  
  template <class T, class Op>
    void joinOp<T,Op>::Join(Map &t2s, const sequence &seq){ 
    for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
      join(t[*i],s[t2s[*i]]) ;
    }
  }
  
  template<class Type,class Op> class joinOp<blackbox<Type>,Op> : public joiner {
    Op join ;
  public:
    CPTR<joiner> clone() 
      { return CPTR<joiner>(new joinOp<blackbox<Type>,Op> ); }
    
    virtual void SetArgs(storeRepP &target, storeRepP &source)
    { std::cerr << "Blackbox Joiner should not be called" << std::endl; }
    
    virtual storeRepP getTargetRep()
      { return 0; }
    
    virtual void Join(const sequence &seq)
      {std::cerr << "Blackbox Joiner should not be called" << std::endl; }
    
    virtual void Join(Map &t2s, const sequence &seq)
      {std::cerr << "Blackbox Joiner should not be called" << std::endl; }
  } ;

  template<class Type,class Op> class joinOp<param<Type>,Op> : public joiner {
    Op join ;
    param<Type> s,t ;
  public:
    CPTR<joiner> clone() 
      { return CPTR<joiner>(new joinOp<param<Type>,Op> ); }
    
    virtual void SetArgs(storeRepP &target, storeRepP &source)
      { s.setRep(source) ; t.setRep(target) ; }
    
    virtual storeRepP getTargetRep()
      { param<Type> st ; storeRepP rep = st.Rep(); return rep; }
    
    virtual void Join(const sequence &seq)
      {join(*t,*s) ;}
    
    virtual void Join(Map &t2s, const sequence &seq)
      {join(*t,*s) ;}
  } ;

  template<class Type,class Op> class joinOp<storeVec<Type>,Op> : public joiner {
    Op join ;
    storeVec<Type> s,t ;
  public:
    virtual CPTR<joiner> clone()
      { return CPTR<joiner>(new joinOp<storeVec<Type>,Op> ); }
    
    virtual void SetArgs(storeRepP &target, storeRepP &source)
      { s.setRep(source) ; t.setRep(target) ; }
    
    virtual storeRepP getTargetRep()
      { storeVec<Type> st ; storeRepP rep = st.Rep(); return rep; }
    
    virtual void Join(const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        Vect<Type> m = t[*i] ;
        join(m,s[*i]) ;
      }
    }

    virtual void Join(Map &t2s, const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        Vect<Type> m = t[*i] ;
        join(m,s[t2s[*i]]) ;
      }
    }
  } ;

  template<class Type,class Op> class joinOp<storeMat<Type>,Op> : public joiner {
    Op join ;
    storeMat<Type> s,t ;
  public:
    virtual CPTR<joiner> clone()
    { return CPTR<joiner>(new joinOp<storeMat<Type>,Op> ); }

    virtual void SetArgs(storeRepP &target, storeRepP &source)
    { s.setRep(source) ; t.setRep(target) ; }

    virtual storeRepP getTargetRep()
    { storeMat<Type> st ; storeRepP rep = st.Rep(); return rep; }

    virtual void Join(const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        Mat<Type> m = t[*i] ;
        join(m,s[*i]) ;
      }
    }

    virtual void Join(Map &t2s, const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        Mat<Type> m = t[*i] ;
        join(m,s[t2s[*i]]) ;
      }
    }
  } ;

  template<class Type, class Op> class joinOp<multiStore<Type>,Op> :
  public joiner {
    Op join ;
    multiStore<Type> s,t ;
  public:
    virtual CPTR<joiner> clone()
    { return CPTR<joiner>(new joinOp<multiStore<Type>,Op> ); }

    virtual void SetArgs(storeRepP &target, storeRepP &source)
    { s.setRep(source) ; t.setRep(target) ; }

    virtual storeRepP getTargetRep()
    { multiStore<Type> st ; storeRepP rep = st.Rep(); return rep; }

    virtual void Join(const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        Vect<Type> m = t[*i] ;
        join(m,s[*i]) ;
      }
    }

    virtual void Join(Map &t2s, const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        Vect<Type> m = t[*i] ;
        join(m,s[t2s[*i]]) ;
      }
    }
  } ;
    
  template <class T, class Op > class apply_rule : public rule_impl {
  protected:
    apply_rule() { rule_class(APPLY) ; }
    void name_store(const std::string &n, store_instance &si)
    { rule_impl::name_store(n,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    void join(typename T::containerType &t1,
              const typename T::containerType &t2) {
      Op f ;
      f(t1,t2) ;
    }
    
    template<class U> void join(Vect<U> t1,
                                const_Vect<U> t2) {
      Op f ;
      f(t1,t2) ;
    }

    template<class U> void join(Mat<U> t1,
                                const_Mat<U> t2) {
      Op f ;
      f(t1,t2) ;
    }
  public:
    virtual CPTR<joiner> get_joiner() {
      CPTR<joinOp<T,Op> > jp = new joinOp<T,Op> ;
      return CPTR<joiner>(jp) ;
    }

  } ;
  
  template <class T> struct NullOp {
    void operator()(T &res, const T &arg)
    { std::cerr << "join should not be called for NullOp" << std::endl; }
    void operator()(Vect<T> &res, const Vect<T> &arg)
    { std::cerr << "join should not be called for NullOp" << std::endl; }
  } ;

  template <class T> struct Summation {
    void operator()(T &res, const T &arg)
    { res += arg ; }
    template <class U> void operator()(T &res, const U &arg)
    { res += arg ; }
  } ;

  template <class T> struct Product {
    void operator()(T &res, const T &arg)
    { res *= arg ; }
    template <class U> void operator()(T &res, const U &arg)
    { res *= arg ; }

  } ;

  template <class T> struct Maximum {
    void operator()(T &res ,const T &arg)
    { res = max(res,arg) ; }
    template <class U> void operator()(T &res, const U &arg)
    { res = max(res,arg) ; }
  } ;
  template <class T> struct Maximum<Vect<T> > {
    template <class U> void operator()(Vect<T> &res ,const U &arg)
    {
      int vs = res.getSize() ;
      for(int i=0;i<vs;++i)
        res[i] = max(res[i],arg[i]) ;
    }
  } ;  

  template <class T> struct Minimum {
    void operator()(T &res, const T &arg)
    { res = min(res,arg) ; }
    template <class U> void operator()(T &res, const U &arg)
    { res = min(res,arg) ; }
  } ;
  template <class T> struct Minimum<Vect<T> > {
    template <class U> void operator()(Vect<T> &res ,const U &arg)
    {
      int vs = res.getSize() ;
      for(int i=0;i<vs;++i)
        res[i] = min(res[i],arg[i]) ;
    }
  } ;  

  
  class rule {
  public:
    enum rule_type {BUILD=0,COLLAPSE=1,GENERIC=2,TIME_SPECIFIC=3,INTERNAL=4} ;
    struct info {
      rule_implP rule_impl ;
      
      rule_impl::info desc ;
      
      std::string rule_ident ;
      time_ident source_level, target_level ;
      variableSet source_vars, target_vars, map_vars, constraint_vars ;
      
      rule_type rule_class ;
      bool output_is_parameter ;
      bool time_advance ;
      std::string internal_qualifier ;
      std::string impl_name ;
      const std::string &name() const { return rule_ident ; }
      info() { rule_ident = "NO_RULE" ;}
      info(const rule_implP &fp) ;
      info(const info &fi, time_ident tl) ;
      // prepend time_ident to info
      info(time_ident tl, const info& fi) ;
      info(const std::string &s) ;

      const info &get_info() const { return *this ; }
      const variableSet &sources() const { return source_vars ; }
      const variableSet &targets() const { return target_vars ; }
      const variableSet &maps() const { return map_vars ; }
      const variableSet &constraints() const {return constraint_vars; }
      
      const std::string &qualifier() const { return internal_qualifier ; }
      rule_type type() const { return rule_class ; }
      time_ident target_time() const { return target_level ; }
      time_ident source_time() const { return source_level ; }
      time_ident time() const { return source_level ; }
      int ident() const { return rule::rdb->get_id(*this) ; }
      
      rule_implP get_rule_implP() const ;
      // SH - namespace support
      rule_implP add_namespace(const std::string& n) const { 
	return rule_impl->add_namespace(n) ;
      }
    } ;
  private:
    friend struct rule::info ;
    friend rule promote_rule(const rule&, const time_ident&) ;
    friend rule prepend_rule(const rule&, const time_ident&) ;
    struct rule_db {
      std::vector<info> fiv ;
      std::map<std::string,int> fmap ;
      rule_db() {
        info no_rule ;
        fiv.push_back(no_rule) ;
        fmap[no_rule.name()] = -1 ;
      }
      int get_id(const info &fi) {
        std::map<std::string,int>::iterator fmi ;
        std::string fname = fi.name() ;
        if((fmi = fmap.find(fname)) == fmap.end()) {
          fiv.push_back(fi) ;
          int id = - int(fiv.size()) ;
          fmap.insert(std::pair<std::string,int>(fname,id)) ;
          return id ;
        }
        return fmi->second ;
      }
      int query_name(std::string &name) {
	std::map<std::string,int>::iterator fmi ;
	if((fmi = fmap.find(name)) == fmap.end()) {
	  fmi = fmap.find("NO_RULE");
	}

	return fmi->second;
      }
      const info &get_info(int id) const { return fiv[-(id+1)] ; }
    } ;
    
    static rule_db *rdb ;
    int id ;
    void create_rdb() {
      if(0==rdb) { 
	rdb = new rule::rule_db ;
	rule();
      }
    }
  protected:
    rule(const rule::info& ri)
      { create_rdb(); id = rdb->get_id(ri) ; }
  public:
    static void rdb_cleanup() { if(rdb) { delete rdb ; rdb = 0; } }
    static int rdb_size()
    { WARN(rdb->fiv.size() != rdb->fmap.size()) ; return rdb->fiv.size() ; }
    rule() { create_rdb() ; id = rdb->get_id(info()) ;}
    explicit rule(int i)
      { create_rdb(); id = i ; }
    rule(const rule_implP &fp)
      { create_rdb(); id = rdb->get_id(info(fp)) ;}
    rule(rule f, time_ident tl)
      { create_rdb(); id = rdb->get_id(info(rdb->get_info(f.id),tl)) ; }
    // prepend time_ident to rule f
    rule(time_ident tl, rule f)
      { create_rdb(); id = rdb->get_id(info(tl,rdb->get_info(f.id))) ; }
    rule(const std::string &s)
      { create_rdb(); id = rdb->get_id(info(s)) ; }
      
    // SH - Namespace support
    // We use a new rule_implP identical to the original but with namespace'd variables to construct a new rule (rule(rule_implP&))
    // get_rule_implP() gives us our rule_implP for this rule, the rule_implP adds the namespace to a copy of itself
    rule add_namespace(const std::string& n) const { 
      return rule(get_rule_implP()->add_namespace(n));
    } 
    rule parent() const { return rule(*this,time().parent()) ; }

    /* 
    std::ostream &Print(std::ostream &s) const
      { s << rdb->get_info(id).name() ; return s ; }
    */
    // remove the stuffs before "#" if any.
    std::ostream &Print(std::ostream &s) const {
      std::string name = rdb->get_info(id).name() ;
      std::string::iterator pos ;
      pos = std::find(name.begin(),name.end(),'#') ;
      s << std::string( (pos==name.end()?name.begin():pos+1),name.end()) ;
      return s ;
    }
    // this function is used to rename a rule
    // i.e., modify the corresponding string inside
    void rename(const std::string&) ;

    static rule get_rule_by_name(std::string &name);
        
    bool operator<(const rule &f) const { return id < f.id ; }
    bool operator==(const rule &f) const { return id == f.id ; }
    bool operator!=(const rule &f) const { return id != f.id ; }

    int ident() const { return id ; }
    const rule::info &get_info() const { return rdb->get_info(id) ; }

    const variableSet &sources() const { return rdb->get_info(id).sources(); }
    const variableSet &targets() const { return rdb->get_info(id).targets(); }

    // rename function that renames variables in the rule
    // according to the rename map passed in. It is the
    // general interface for rule promotion.
    rule rename_vars(std::map<variable,variable>& rvm) const ;
    rule_type type() const { return rdb->get_info(id).type() ; }
    time_ident target_time() const { return rdb->get_info(id).target_time() ;}
    time_ident source_time() const { return rdb->get_info(id).source_time() ;}
    time_ident time() const { return rdb->get_info(id).time() ; }
    rule_implP get_rule_implP() const { return rdb->get_info(id).get_rule_implP() ;}
  } ;

  inline std::ostream &operator<<(std::ostream &s, const rule &f)
    { return f.Print(s) ; }

  // global rule promotion functions
  // rule promotion
  rule promote_rule(const rule& r, const time_ident& t) ;
  // time prepend to a rule
  rule prepend_rule(const rule& r, const time_ident& t) ;

  class ruleSet : public intervalSet {
  public:
    ruleSet() {}
    explicit ruleSet(const exprP &e) ;
    explicit ruleSet(const intervalSet &v)
      {*(static_cast<intervalSet *>(this)) = v ;}
    ruleSet &operator=(const intervalSet &v)
      {*(static_cast<intervalSet *>(this)) = v ; return *this ;}
    ruleSet &operator+=(const rule &v)
      { *this += v.ident() ; return *this ; }
    ruleSet &operator-=(const rule &v)
      { *this -= v.ident() ; return *this ; }
    bool inSet(const rule &v) const
      { return intervalSet::inSet(v.ident()) ; }
    class ruleSetIterator {
      intervalSet::const_iterator ii ;
    public:
      ruleSetIterator() {}
      ruleSetIterator(const intervalSet::const_iterator &i)
        { ii = i ; }
      rule operator*() const { return rule(*ii) ; }
      const rule::info *operator->() const
        { return &(rule(*ii).get_info()) ; }
      ruleSetIterator &operator++() { ++ii ; return *this ;}
      ruleSetIterator operator++(int )
        { return ruleSetIterator(ii++); }
      bool operator==(const ruleSetIterator &i) const 
        { return ii == i.ii ; }
      bool operator!=(const ruleSetIterator &i) const
        { return ii != i.ii ; } ;
    } ;
    typedef ruleSetIterator const_iterator ;
    const_iterator begin() const {
      return const_iterator(intervalSet::begin()) ; }
    const_iterator end() const {
      return const_iterator(intervalSet::end()) ; }
    std::ostream &Print(std::ostream &s) const ;

  } ;

  inline std::ostream &operator<<(std::ostream &s, const ruleSet& v)
    { return v.Print(s) ; }

  class register_rule_type {
  public:
      virtual ~register_rule_type() {}
      virtual rule_implP get_func() const = 0 ;
      virtual bool is_module_rule() const = 0 ;
      
  } ; 

  class register_rule_impl_list ;
  
  class rule_impl_list {
  public:
    class rule_list_iterator ;
    friend class rule_list_iterator ;
    friend class register_rule_impl_list ;
    class rule_list_ent {
    public:
      rule_list_ent(register_rule_type *p, rule_list_ent *nxt) :
        rr(p), next(nxt) {}
      register_rule_type *rr ;
      rule_list_ent *next ;
    } ;
    rule_list_ent *list ;
  public:
    class rule_list_iterator {
      rule_list_ent *p ;
    public:
      rule_list_iterator(rule_list_ent *ptr) : p(ptr) {}
      rule_implP operator*() { return p->rr->get_func() ; }
      rule_list_iterator &operator++() {
        p = p->next ;
        return *this ;
      }
      rule_list_iterator operator++(int ) {
        rule_list_iterator tmp(p) ;
        p = p->next ;
        return tmp ;
      }
      rule_list_ent* get_p() { return p; } 
      bool operator==(const rule_list_iterator &i) { return i.p == p ; }
      bool operator!=(const rule_list_iterator &i) { return i.p != p ; }
    } ;
    typedef rule_list_iterator iterator ;
    
    rule_impl_list() {list = 0 ; }
    ~rule_impl_list() ;
    
    void push_rule(register_rule_type *rr) ;
    iterator begin() { return iterator(list) ; }
    iterator end() { return iterator(0) ; }
    void clear() ;
    void copy_rule_list(const rule_impl_list &rl) ;
    void copy_rule_list(const register_rule_impl_list &rl) ;
    rule_impl_list(const rule_impl_list &x) {
      list = 0 ;
      copy_rule_list(x) ;
    }
    rule_impl_list &operator=(const rule_impl_list &x) {
      list = 0 ;
      copy_rule_list(x) ;
      return *this ;
    }
  } ;
  class register_rule_impl_list : public rule_impl_list {
  public:
    static rule_list_ent *global_list ;
    register_rule_impl_list() {}
    ~register_rule_impl_list() ;
    void clear() ;
    bool empty() ;
    void push_rule(register_rule_type *p) ;
    iterator begin() { return iterator(global_list) ; }
    iterator end() { return iterator(0) ; }
  } ;
  extern register_rule_impl_list register_rule_list ;    
  extern rule_impl_list global_rule_list ;    
  template<class T> class register_rule : public register_rule_type {
  public:
    register_rule() { register_rule_list.push_rule(this) ; }
    virtual bool is_module_rule()  const{ return false; }
    virtual rule_implP get_func() const { return new copy_rule_impl<T> ; }
  } ;
  
  class register_module : public register_rule_type {
  public:
    register_module() { register_rule_list.push_rule(this) ; }
    virtual ~register_module() {}
    virtual bool is_module_rule() const{ return true ; }
    virtual rule_implP get_func() const { return 0 ; }
    virtual std::string load_module() const ;
    virtual std::string using_nspace() const ;
    virtual std::string input_vars() const = 0 ;
    virtual std::string output_vars() const = 0 ;
    virtual std::string load_nspace() const ;
    virtual std::string disable_compute_vars() const ;
  } ;
  
  class rule_db {
    typedef std::map<variable,ruleSet> varmap ;
    typedef varmap::const_iterator vc_iterator ;
    typedef std::map<std::string,rule> rule_map_type ;
      
    static const ruleSet EMPTY_RULE ;
    ruleSet known_rules ;
    // rules set for default and optional rules
    ruleSet default_rules ;
    ruleSet optional_rules ;

    varmap srcs2rule,trgt2rule ;
    rule_map_type name2rule ;
    // partition rules according their keyspace
    std::map<std::string,ruleSet> keyspace2rule ;

  public:
    void add_rule(const rule_implP &fp) ;
    void add_rule(rule f) ;
    void add_rules(rule_impl_list &gfl) ;
    void add_rules(register_rule_impl_list &gfl) ;
    void remove_rule(rule f) ;
    void remove_rules(const ruleSet& rs) ;
    rule_implP get_rule(const std::string &name) {
      return name2rule[name].get_info().rule_impl->new_rule_impl() ;
    }
    
    const ruleSet &all_rules() const { return known_rules ; }
    // return all the rules in "keyspace_tag"
    const ruleSet& all_rules(const std::string& keyspace_tag) const {
      std::map<std::string,ruleSet>::const_iterator
        mi = keyspace2rule.find(keyspace_tag) ;
      if(mi == keyspace2rule.end())
        return EMPTY_RULE ;
      else
        return mi->second ;
    }
    const ruleSet& get_default_rules() const {return default_rules ;}
    const ruleSet& get_optional_rules() const {return optional_rules ;}
    const ruleSet &rules_by_source(variable v) const {
      vc_iterator vmi = srcs2rule.find(v) ;
      if(vmi == srcs2rule.end())
        return EMPTY_RULE ;
      else
        return vmi->second;
    }
    const ruleSet &rules_by_target(variable v) const {
      vc_iterator vmi = trgt2rule.find(v) ;
      if(vmi == trgt2rule.end())
        return EMPTY_RULE ;
      else
        return vmi->second;
    }
      
  } ;

}

#endif
