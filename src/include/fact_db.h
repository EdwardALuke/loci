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
#ifndef FACT_DB_H
#define FACT_DB_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/except.h>

#include <string>
#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <store_rep.h>
#include <store.h>
#include <Tools/variable.h>
#include <Map_rep.h>
#include <Map.h>
#include <DMap.h>
#include <key_manager.h>
#include <keyspace.h>

#define LOCI_COMPAT_MODE1
#define LOCI_HAS_STATIC_KEYSPACE 1

namespace Loci {
  class rule_db ;
  
  class fact_db {
  public:
    
    struct distribute_info : public CPTR_type {
      struct dist_data {
        int proc ;
        entitySet entities ;
        int size ;
        dist_data(int p, entitySet &e) :proc(p),entities(e)
        { size = e.size() ; }
      } ;
      
      int myid ;
      int isDistributed ;
      Map l2g ; // local numbering to global numbering
      store<unsigned char> key_domain ; // local numbering to key domain
      Map l2f ; // local numbering to file numbering

#ifdef LOCI_COMPAT_MODE1
      dMap g2l,g2f ;
#endif
     
      std::vector<dMap> g2lv ; // global numbering to local numbering
      // indexed by key domain
      entitySet my_entities ;
      entitySet comp_entities;//local numbering of global_comp_entities
      
      std::vector<dist_data> copy ; // Data describing which processors own
      // the entities in the clone region
      std::vector<dist_data> xmit ; // Data describing which processors to
      // send entities that this processor owns
      int copy_total_size ;
      int xmit_total_size ;
      //      dMap remap ;
      std::vector<dMap> g2fv ; // Global to file numbering
      distribute_info() {} ;
    }  ;
    std::vector<std::vector<entitySet> > init_ptn ;
    /// Global numbering partition indexed by key space
    std::vector<int> gmax_alloc ;
    std::vector<std::string> keyDomainName ;

    //Related to Duplication: global_comp_entities set defines
    //entities that can be computed on a processor
    entitySet global_comp_entities; 

    typedef CPTR<distribute_info> distribute_infoP ;
    distribute_infoP distributed_info ;

    std::set<std::vector<variableSet> > intensive_output_maps;

#ifdef DYNAMICSCHEDULING
    // experimental inclusion of keyspace partitions
    std::map<std::string,KeySpaceP> keyspace ;
    // experimental inclusion of a key manager
    KeyManagerP key_manager ;
#endif
    
    variable remove_synonym(variable v) const {
      std::map<variable,variable>::const_iterator mi ;
      while((mi=synonyms.find(v)) != synonyms.end())
        v = mi->second ;
      return v ;
    }
  private:
    struct fact_info {
      store_refP data_rep ;
    } ;
    
    /*! all variables that point to the same storeRepP, the second of the pair at the end of chain
      is the variable suppose to appear in fmap*/
    std::map<variable,variable> synonyms ;
    
    // a copy function
    void copy_all_from(const fact_db& f) ;

    std::pair<entitySet, entitySet> get_dist_alloc(int size, size_t kd) ;
    
    int maximum_allocated ;
    int dist_from_start ;
    
    std::map<variable,fact_info> fmap ;
    std::map<variable,storeRepP> tmap ;
    std::vector<std::string> nspace_vec ;
    // support for multiple queries and experimental
    // extensions to the fact_db to distinguish
    // extensional facts and intensional facts
    variableSet extensional_facts ;

    variable add_namespace(variable v) const ;
    // this is the basic method that creats a fact
    // in the fact_db. It is served as the basis
    // for create_fact & the create_intensional_fact methods
    void create_pure_fact(const variable& v, storeRepP st) ;
  public:
    fact_db() ;
    ~fact_db() ;
    // copy constructor
    fact_db(const fact_db &f) {
      copy_all_from(f) ;
    }
    // the assignment operator
    fact_db &operator=(const fact_db &f) {
      if(&f != this) {
        copy_all_from(f) ;
      }
      return *this ;
    }
    void set_maximum_allocated(int i) ;
    int get_max_alloc(int kd ) ;

    int getKeyDomain(std::string name) ;

    void set_variable_type(variable v, storeRepP st) ;
    void set_variable_type(std::string vname, storeRepP st)
    { set_variable_type(variable(vname),st) ;}
    void set_variable_type(variable v, store_instance &si)
    { set_variable_type(v,si.Rep()) ; }
    void set_variable_type(std::string vname, store_instance &si)
    { set_variable_type(variable(vname),si) ; }

    storeRepP get_variable_type(variable v) const ;
    storeRepP get_variable_type(std::string vname) const
    { return get_variable_type(variable(vname)) ;}

      
    
    // create_fact now defaults to create an extensional fact
    // as this is the primary interface for users of Loci
    void create_fact(const variable& v, storeRepP st) {
      variable tmp_v =add_namespace(v) ;
      create_pure_fact(tmp_v,st) ;
      extensional_facts += tmp_v ;
    }
    void create_fact(const std::string& vname, storeRepP st) {
      variable v = variable(vname) ;
      create_fact(v,st) ;
    }
    void create_fact(const variable& v, store_instance &si) {
      variable tmp_v = add_namespace(v) ;
      create_pure_fact(tmp_v,si.Rep()) ;
      si.setRep(get_variable(v)) ;
      extensional_facts += tmp_v ;
    }
    void create_fact(const std::string& vname, store_instance &si) {
      variable v = add_namespace(variable(vname)) ;
      create_pure_fact(v,si.Rep()) ;
      si.setRep(get_variable(vname)) ;
      extensional_facts += v ;
    }
    
    void update_fact(variable v, storeRepP st) ;
    void update_fact(std::string vname, storeRepP st)
    { update_fact(variable(vname),st) ;}
    void update_fact(variable v, store_instance &si)
    { update_fact(v,si.Rep()) ; si.setRep(get_variable(v)) ; }
    void update_fact(std::string vname, store_instance &si)
    { update_fact(variable(vname),si) ; }

    void replace_fact(variable v, storeRepP st)
    { remove_variable(v) ; create_fact(v,st) ; }
    void replace_fact(std::string vname, storeRepP st)
    { replace_fact(variable(vname),st) ;}
    void replace_fact(variable v, store_instance &si)
    { replace_fact(v,si.Rep()) ; si.setRep(get_variable(v)) ; }
    void replace_fact(std::string vname, store_instance &si)
    { replace_fact(variable(vname),si) ; }

    storeRepP get_fact(variable &v) { return get_variable(v); }
    storeRepP get_fact(std::string vname)
    { return get_variable(variable(vname)) ; }

    storeRepP get_variable(variable v) ;
    storeRepP get_variable(std::string vname)
    { return get_variable(variable(vname)) ; }
    
    void remove_variable(variable v) ;
    
    void synonym_variable(variable v, variable synonym) ;

    void rotate_vars(const std::list<variable> &lvars) ;
#ifdef DYNAMICSCHEDULING
    // this method will adjust the rotation vars
    // so that the all the history part will be
    // having the same domain as the current part
    // for example, we have A{n+1}, A{n}, A{n-1}
    // that is rotating. this method will make
    // A{n} and A{n-1} have the same domain as
    // A{n+1}, for those out of A{n+1}'s domain,
    // they are erased, for those missing, they
    // are copied from A{n+1}. this method is
    // mainly for those dynamic rule's rotation
    // list because during the execution, there
    // might have been erase and insertion occurred
    // and we therefore must adjust the history
    // variables accordingly.
    void adjust_rotation_vars(const std::list<variable>& lvars) ;
#endif
    void set_namespace(std::string name_space){
      nspace_vec.insert(nspace_vec.begin(), 1, name_space) ; 
    }
    void remove_namespace() {
      nspace_vec.pop_back() ;
    }  
    void unset_namespace() {
      nspace_vec.clear() ;
    }
    //    entitySet get_allocation(int size) {
    //      entitySet alloc = interval(maximum_allocated,maximum_allocated+size-1) ;
    //      maximum_allocated += size ;
    //      return alloc ;
    //    }

    void update_remap(const std::vector<std::pair<int, int> > &remap_update,
		      size_t kd);

    std::pair<entitySet, entitySet> get_distributed_alloc(int size,size_t kd) ;
    std::pair<entitySet, entitySet> get_distributed_alloc(const std::vector<int> &remap_entities, size_t kd);
    std::pair<entitySet, entitySet> get_distributed_alloc(int size, int offset, size_t kd);
    int is_distributed_start() {return dist_from_start ;}
#ifdef LOCI_COMPAT_MODE1
    std::vector<entitySet>& get_init_ptn() {return init_ptn[0] ;}
    std::pair<entitySet, entitySet> get_distributed_alloc(int size) {
      return get_distributed_alloc(size,0) ;
    }
    std::pair<entitySet, entitySet> get_distributed_alloc(const std::vector<int> &remap_entities) {
      return get_distributed_alloc(remap_entities,0) ;
    }
#endif
    std::vector<entitySet>& get_init_ptn(size_t kd) {return init_ptn[kd] ;}
    void  put_init_ptn(std::vector<entitySet> &t_init, size_t kd ) {init_ptn[kd] = t_init ;}
    int getNumKeyDomains() const { return init_ptn.size() ; }
    
    fact_db::distribute_infoP get_distribute_info() ;
    void put_distribute_info(fact_db::distribute_infoP dp) ;
    bool isDistributed() ;
    
    
    variableSet get_typed_variables() const ;
    std::ostream &write(std::ostream &s) const ;
    std::istream &read(std::istream &s) ;
    void setupDefaults(const rule_db &rdb) ;
    std::istream& read_vars(std::istream& s, const rule_db& rdb) ;

    void read_vars(std::string filename, const rule_db &rdb) {
      std::ifstream ifile(filename.c_str(),std::ios::in) ;
      if(ifile.fail()) {
        std::string error = std::string("Can't open '") +filename + "'" ;
        throw(StringError(error)) ;
      }
      read_vars(ifile,rdb) ;
      ifile.close() ;
    }

    /////////////////////////////////////////////////////////
    // support methods for extensional & intensional facts //
    /////////////////////////////////////////////////////////
    variableSet get_extensional_facts() const {
      return extensional_facts ;
    }
    variableSet get_intensional_facts() const {
      return variableSet(get_typed_variables()-extensional_facts) ;
    }

    // we still provide these methods with explicit name to
    // create extentional facts, they are just as the same as
    // the default create_fact methods
    void create_extensional_fact(const variable& v, storeRepP st) {
      create_fact(v,st) ;
    }
    void create_extensional_fact(const std::string& vname, storeRepP st) {
      create_fact(vname,st) ;
    }
    void create_extensional_fact(const variable& v, store_instance &si) {
      create_fact(v,si) ;
    }
    void create_extensional_fact(const std::string& vname,
                                 store_instance &si) {
      create_fact(vname,si) ;
    }
    // this method will convert all intensional facts (if any) in
    // the fact database to extensional facts
    void make_all_extensional() {
      variableSet intensional_facts = get_intensional_facts() ;
      extensional_facts += intensional_facts ;
    }
    // and then we have the corresponding intensional facts creation
    void create_intensional_fact(const variable& v, storeRepP st) {
      variable v_tmp = add_namespace(v) ;
      create_pure_fact(v_tmp,st) ;
    }
    void create_intensional_fact(const std::string& vname, storeRepP st) {
      create_pure_fact(add_namespace(variable(vname)),st) ;
    }
    void create_intensional_fact(const variable& v, store_instance &si) {
      variable v_tmp = add_namespace(v) ;
      create_pure_fact(v_tmp,si.Rep()) ;
      si.setRep(get_variable(v_tmp)) ;
    }
    void create_intensional_fact(const std::string& vname,
                                 store_instance &si) {
      variable v = add_namespace(variable(vname)) ;
      create_pure_fact(v,si.Rep()) ;
      si.setRep(get_variable(v)) ;
    }
    // this method erases all intensional facts
    void erase_intensional_facts() {
      variableSet intensional_facts = get_intensional_facts() ;
      for(variableSet::const_iterator vi=intensional_facts.begin();
          vi!=intensional_facts.end();++vi)
        remove_variable(*vi) ;
    }
    // this method will convert an intensional fact to
    // a extensional fact
    void make_extensional_fact(const variable& v) ;
    void make_extensional_fact(const std::string& vname) {
      make_extensional_fact(variable(vname)) ;
    }
    // this method will convert a extensional fact
    // to an intensional one
    void make_intensional_fact(const variable& v) ;
    void make_intensional_fact(const std::string& vname) {
      make_intensional_fact(variable(vname)) ;
    }

   
    void write_all_hdf5(const char *filename) ;
    void read_all_hdf5(const char *filename) ;
    void write_hdf5(const char *filename, variableSet &vars) ;
    void read_hdf5(const char *filename, variableSet &vars) ;

   
    void Print_diagnostics() ;

#ifdef DYNAMICSCHEDULING
    // experimental code to create keyspace from the
    // global registered keyspace list
    // returns "true" to indicate the methods succeeded,
    // "false" to indicate an error.
    bool create_keyspace(KeySpaceList& global_list) ;
    KeySpaceP get_keyspace(const std::string& kname) ;
    KeyManagerP get_key_manager() const {return key_manager ;}
    void init_key_manager() ;
#endif
  } ;
  
  void reorder_facts(fact_db &facts, dMap &remap) ;
  void serial_freeze(fact_db &facts) ;
  

  entitySet collect_entitySet(entitySet e, fact_db &facts) ;

  extern fact_db *exec_current_fact_db ;

  inline entitySet collect_entitySet(entitySet e)
  { return collect_entitySet(e,*exec_current_fact_db) ; }

  entitySet all_collect_entitySet(const entitySet &e) ;

  inline entitySet all_collect_entitySet(entitySet localset,fact_db &facts) {
    if(facts.is_distributed_start())
      return Loci::all_collect_entitySet(localset) ;
    return localset ;
  }
}

#endif
