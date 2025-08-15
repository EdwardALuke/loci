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
#ifndef LOCI_MOD_DB_H
#define LOCI_MOD_DB_H
#define LOCI_WITH_MODULES
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <dlfcn.h>

#include <rule.h>
#ifdef DYNAMICSCHEDULING
#include <keyspace.h>
#endif
#include <fact_db.h>
#include <string>
namespace Loci {
  std::string remove_space(const std::string &str) ; 
  void parse_str(const std::string& str, std::vector<std::string> &str_vec) ;
  class mod {
  public: 
    struct mod_info {
      rule_impl_list loaded_rule_list ;
#ifdef DYNAMICSCHEDULING
      KeySpaceList loaded_keyspace_list ;
#endif
      std::string mod_name ;
      void *m_library ;
      void (*m_init_model)(fact_db &facts, rule_db &rdb, const char *problem_name) ;
      std::string name() { return mod_name ; } 
      mod_info(rule_impl_list& rl, std::string str) {
	loaded_rule_list.copy_rule_list(rl) ;
	mod_name = str ;
      }
      mod_info(rule_impl_list& rl) { 
	loaded_rule_list.copy_rule_list(rl) ;
	mod_name = std::string("") ;
      }
      mod_info() {
	mod_name = "" ;
	m_library = 0 ;
	m_init_model = 0 ;
      }
      mod_info(const mod_info& mi) {
	loaded_rule_list.copy_rule_list(mi.loaded_rule_list) ;
#ifdef DYNAMICSCHEDULING
        loaded_keyspace_list.copy_space_list(mi.loaded_keyspace_list) ;
#endif
	mod_name = mi.mod_name ;
        m_library = mi.m_library ;
        m_init_model = mi.m_init_model ;
      }
      ~mod_info() { } 
    } ;
  private: 
    friend struct mod_info ;
    struct mod_db {
      std::map<std::string,mod_info> mod_map ;
      typedef std::map<std::string,mod_info>::iterator MI ;
      void put_info(mod_info& md) {
	MI msi ;
	if((msi = mod_map.find(md.name())) == mod_map.end())
	  mod_map.insert(std::make_pair(md.name(),md ));
      }
      
      mod_info &get_info(const std::string &str) ;
      mod_info &get_info(const std::string &str, const std::string &to_str, const char* problem_name, fact_db &facts) ;
    } ;
  
    static mod_db *mdb ;
    static void cleanup_db() { if(mdb) delete mdb ; }
    void create_mod_db() {if(0==mdb) mdb = new mod::mod_db ; }
  public:
    mod(const std::string& name) {
      create_mod_db() ;
      get_info(name) ;
    }
    mod(const std::string& name, const std::string& to_str, const char* problem_name, fact_db &facts) {
      create_mod_db() ;
      get_info(name, to_str, problem_name, facts) ;
    }
    mod() {
     create_mod_db() ;
     mod_info md = mod_info(global_rule_list) ;
    }
    mod_info &get_info(const std::string name) { 
      return mdb->get_info(name) ; }
    mod_info &get_info(const std::string name, const std::string &to_str, const char* problem_name, fact_db &facts) {
      return mdb->get_info(name, to_str, problem_name, facts); }
    void put_info(mod_info &md) {mdb->put_info(md) ; }
  } ;
  
  void load_module(const std::string from_str, const std::string to_str, rule_db& rdb, std::set<std::string> &str_set) ;
  inline void load_module(std::string module_name,rule_db &rdb) {
    std::set<std::string> str_set ;
    load_module(module_name,"",rdb,str_set) ;
  }
  
  void load_module(const std::string from_str, const std::string to_str, const char* problem_name, fact_db &facts, rule_db& rdb, std::set<std::string> &str_set) ;	 
  void AddModuleSearchDir(std::string dirname) ;
}

#endif
