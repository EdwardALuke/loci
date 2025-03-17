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
#ifndef KEY_SPACE_H
#define KEY_SPACE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#ifdef DYNAMICSCHEDULING
#include <vector>
#include <map>
#include <set>
#include <string>
#include <mpi.h>
#include <Loci_types.h>
#include <Tools/intervalSet.h>
#include <Tools/cptr.h>
#include <Tools/expr.h>
#include <Tools/variable.h>
#include <Map.h>
#include <DMap.h>
#include <store_rep.h>

namespace Loci {

  typedef Entity Key ;
  typedef entitySet KeySet ;

  // this type is used to categorize the dynamism of a keyspace.
  // it is mainly used to help to generate schedules for
  // the rules inside a particular keyspace

  // rightnow, the first twos are well established.
  // STATIC -- means a keyspace can be planned completely
  //   ahead of the actual run, i.e., the old Loci planning
  //   engine can be used to generate an optimal schedule
  // DYNAMIC -- means a keyspace is completely dynamic, that
  //   in general, the best strategy is to NOT preplan
  //   any schedule, rather it is best to rely on a
  //   per-rule evaluation strategy at runtime.
  // SEMI_DYNAMIC -- means that a keyspace is a mixture
  //   of the above two situations, and as of this time,
  //   the most effective strategy is still to be researched
  //   and developed. (we would use the dynamic planning
  //   for now for this type.) also note, that in the future,
  //   we may as well have more types to add to this list.
  //   or, we may only limit the type to be STATIC and DYNAMIC
  //   and the in-between types are represented by the
  //   composition of the keyspaces.

  // Also, these properties are currently supplied as
  // user annotation when creating a keyspace. It is however
  // intended in the future that the scheduler is
  // smart enough to figure out such property on its own.
  // Then at that time, the scheduler is free to choose
  // to adhere to the user annotation or completely
  // ignore them.
  enum KeySpaceDynamism {STATIC, DYNAMIC, SEMI_DYNAMIC} ;

  class rule ;
  class execute_modules ;
  
  class KeySpace: public CPTR_type {
  public:
    typedef CPTR<KeySpace> KeySpaceP ;
  protected:
    // the communicator that determines the group of
    // processess that share this set of keys
    MPI_Comm comm ;
    int rank ;                  // processes rank in comm
    int comm_size ;             // communicator size
    // this is the key distribution across the processess
    // it is in global numbering
    std::vector<KeySet> key_ptn ;
    // the entire key sets in global number
    KeySet keys ;
    // the entire key sets in local numbering,
    // if there is no local numbering, then it is
    // the same as key_ptn[rank] and keys
    KeySet my_keys ;
    dMap g2l ;
    Map l2g ;
  protected:
    // the name of the keyspace
    std::string space_name ;
    // name -> storeRep mapping for the critical structures
    typedef std::map<variable,store_instance*> CSNMap ;
    CSNMap csn_table ;
    // used to indicate whether the csn_table has
    // been hooked with the fact_db's storeRep records
    typedef std::map<variable,bool> CSBMap ;
    CSBMap csn_ready ;

    KeySpaceDynamism dyn_type ;
    std::map<variable, std::string> tunnel_map ;
    std::map<std::string, variableSet> tunnel_map_inverse ;
    variableSet tunnels ;

    // comments may be put for a keyspace
    std::string user_comments ;
    // shadow storeReps are effectively the clone portion
    // of the facts that this keyspace is referring to
    // through its "tunnel" to its target keyspace. it is
    // managed separately by this keyspace because it may
    // choose to communicate whatever subsets of the facts
    // deemed useful and in whatever numbering scheme deemed
    // convenient and helpful.
    std::map<variable, store_refP> shadow ;
    // this is a reference to fact_db's synonyms map
    std::map<variable, variable>* synonyms ;
    variable
    remove_synonym(variable v) const {
      std::map<variable,variable>::const_iterator mi ;
      while( (mi=synonyms->find(v)) != synonyms->end())
        v = mi->second ;
      return v ;
    }
    // all critical structures
    variableSet critical_vars ;
    // all variables generated (managed) by rules in this space
    // they are the ones to redistribute if keyspace needs
    // redistribution
    variableSet control_vars ;
    // flag to indicate whether a local numbering is in use
    bool use_local_number ;
    // // all the rules that belong to this keyspace, filled later
    // ruleSet rules ;

    // these records are used for setting up the run-time
    // variable -> rule execution module mapping, used
    // for invalidate dynamic rule's context and communication.
    // A better implementation is probably to have a visitor
    // class hierarchy defined for all the execution modules,
    // this implementation is a quick method that avoids the
    // visitor classes. First, in the dynamic rule identification
    // stage, the var->rule mapping is constructed. then in
    // the execution module initialization stage, the execution
    // module would register themselves here, and finally in
    // the execution stage, the keyspace init execution module
    // will finally hood up everything properly.
    // We will probably want to revise the design later,
    // perhaps to use the visitor design
    std::map<variable, std::set<std::string> > dcontrol_map ;
    std::map<variable, execute_modules*> dcontrol_vars ;
    std::map<std::string, execute_modules*> dcontrol_rules ;
  protected:    
    void // set space name
    keyspace_name(const std::string& n) {space_name = n ;}

    void // build the critical structure table
    csn_store(const std::string& name,store_instance& si) {
      variable v(expression::create(name)) ;
      csn_table.insert(std::pair<variable,store_instance*>(v,&si)) ;
      csn_ready.insert(std::pair<variable,bool>(v,false)) ;
      critical_vars += v ;
    }

    void
    keyspace_dynamism(KeySpaceDynamism kd) {
      dyn_type = kd ;
    }

    void // set the keyspace's MPI_COMM_GROUP
    MPI_comm_group(MPI_Comm c) {
      comm = c ;
      MPI_Comm_rank(comm, &rank) ;
      MPI_Comm_size(comm, &comm_size) ;
    }

    // define tunnels to other keyspace
    // currently if the target keyspace is not specified,
    // it is assumed to be "main" keyspace
    void
    create_tunnel(const std::string& tunnel_name,
                  const std::string& target_keyspace_name="main") {
      variable v(expression::create(tunnel_name)) ;
      if(tunnels.inSet(v)) {
        std::cerr << "Error: tunnel: " << tunnel_name
                  << " already created in keyspace: "
                  << space_name << std::endl ;
        return ;
      }
      tunnels += v ;
      tunnel_map[v] = target_keyspace_name ;
      tunnel_map_inverse[target_keyspace_name] += v ;
    }

    void
    local_number_scheme(bool b) {
      use_local_number = b ;
    }

    void // add comments
    comments(const std::string& c) { user_comments += c ;}
  private:
    // a copy function used in the copy constructor and
    // assignment operator
    void
    copy_from(const KeySpace& ks) {
      comm = ks.comm ;
      rank = ks.rank ;
      comm_size = ks.comm_size ;
      key_ptn = ks.key_ptn ;
      keys = ks.keys ;
      // deep copy g2l and l2g maps
      entitySet dom = ks.g2l.domain() ;
      g2l.allocate(dom) ;
      for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei)
        g2l[*ei] = ks.g2l[*ei] ;
      dom = ks.l2g.domain() ;
      l2g.allocate(dom) ;
      for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei)
        l2g[*ei] = ks.l2g[*ei] ;

      space_name = ks.space_name ;
      csn_table = ks.csn_table ;
      csn_ready = ks.csn_ready ;
      dyn_type = ks.dyn_type ;
      tunnel_map = ks.tunnel_map ;
      tunnel_map_inverse = ks.tunnel_map_inverse ;
      tunnels = ks.tunnels ;
      user_comments = ks.user_comments ;
      shadow = ks.shadow ;
      synonyms = ks.synonyms ;
      critical_vars = ks.critical_vars ;
      control_vars = ks.control_vars ;
      use_local_number = ks.use_local_number ;

      dcontrol_map = ks.dcontrol_map ;
      dcontrol_vars = ks.dcontrol_vars ;
      dcontrol_rules = ks.dcontrol_rules ;
    }
  public:
    KeySpace():space_name("__@unnamed@keyspace@__"),
               dyn_type(DYNAMIC),use_local_number(false) {
      MPI_comm_group(MPI_COMM_WORLD) ;
      key_ptn.resize(comm_size) ;
    }
    KeySpace(const KeySpace& ks) {copy_from(ks) ;}
    KeySpace& operator=(const KeySpace& ks) {
      if(&ks != this)
        copy_from(ks) ;
      return *this ;
    }

    virtual ~KeySpace() {}

    std::string
    get_comments() const {return user_comments ;}

    std::string
    get_name() const {return space_name ;}

    KeySpaceDynamism
    get_dynamism() const {return dyn_type ;}

    // this one is used to setup storeReps for
    // the critical structures
    void
    setup_csn_rep(variable v, Loci::storeRepP srp) {
      CSNMap::iterator csnii = csn_table.find(v) ;
      if(csnii == csn_table.end()) {
        if(rank == 0) {
          std::cerr << "Error: KeySpace " << space_name
                    << " has no critical structure: " << v << endl ;
        }
        return ;
      }
      csnii->second->setRep(srp) ;
      CSBMap::iterator csbii = csn_ready.find(v) ;
      fatal(csbii == csn_ready.end()) ;
      csbii->second = true ;
    }

    // function for keyspace integrity check
    // returns "true" if the test passes,
    // returns "false" if something is not properly setup,
    // and in this case, the string "err" records the
    // error message.
    bool
    integrity_check(std::string& err) const {
      bool check_pass = true ;
      if(space_name == "__@unnamed@keyspace@__") {
        err += "KeySpace has no name. " ;
        check_pass = false ;
      }
      
      if(csn_ready.empty())
        return check_pass ;
      
      CSBMap::const_iterator csbii = csn_ready.begin() ;

      std::string vars_not_set ;

      if(!csbii->second) {
        vars_not_set = "(" + (csbii->first).get_info().name ;
        check_pass = false ;
      }
      ++csbii ;
      
      for(;csbii!=csn_ready.end();++csbii) {
        if(!csbii->second) {
          std::string s = ", " + (csbii->first).get_info().name ;
          vars_not_set += s ;
          check_pass = false ;
        }
      }

      if(!check_pass) {
        vars_not_set += ")" ;
        err += "KeySpace " + space_name + " critical structures "
          + vars_not_set + " have not been set up yet!" ;
      }

      return check_pass ;
    }

    // return if the space has been named or not
    bool
    named_space() const {
      return space_name != "__@unnamed@keyspace@__" ;
    }

    void
    set_synonyms(std::map<variable,variable>* s) {
      synonyms = s ;
    }

    // setup a shadow variable
    void
    create_shadow(const variable& v, storeRepP sp) ;

    // get the rep information of a shadow variable
    storeRepP
    get_shadow(const variable& v) const ;

    const variableSet&
    get_tunnels() const {
      return tunnels ;
    }

    void
    reset_tunnels(const variableSet& t) {
      tunnels = t ;
    }

    // given a tunnel name, return the space name it extends to
    std::string
    get_tunnel_space(variable t) const {
      t = remove_synonym(t) ;
      std::map<variable,std::string>::const_iterator
        mi=tunnel_map.find(t) ;
      if(mi != tunnel_map.end())
        return mi->second ;
      else                      // return a empty name if not found
        return std::string() ;
    }

    void
    add_control_vars(const variableSet& vs) {
      control_vars += vs ;
    }
    void
    add_control_vars(const variable& v) {
      control_vars += v ;
    }
    const variableSet&
    get_control_vars() const {return control_vars ;}

    const variableSet&
    get_critical_vars() const {return critical_vars ;}

    bool
    has_local_number() const {return use_local_number ;}

    MPI_Comm
    get_mpi_comm() const {return comm ;}

    int
    get_comm_rank() const {return rank ;}

    int
    get_comm_size() const {return comm_size ;}

    bool
    is_distributed() const {return comm_size > 1 ;}

    const std::vector<KeySet>&
    get_key_ptn() const {return key_ptn ;}

    void
    set_key_ptn(const std::vector<KeySet>& ptn) {
      key_ptn = ptn ;
      keys = key_ptn[rank] ;
      my_keys = key_ptn[rank] ;
    }
    
    const KeySet&
    get_keys() const {return keys ;}

    const KeySet&
    get_keys_local() const {return my_keys ;}

    // this method adds new keys into the space
    void
    add_keys(const KeySet& ks) ;

    // while this one removes keys from the space
    void
    remove_keys(const KeySet& ks) ;

    const dMap&
    get_g2l_map() const {return g2l ;}

    const Map&
    get_l2g_map() const {return l2g ;}

    void
    set_dcontrol_map(const variable& v, const rule& r) ;

    void
    register_dcontrol_var(const variable& v, execute_modules* ep) ;

    void
    register_dcontrol_rule(const rule& r, execute_modules* ep) ;

    void
    set_dcontrol() const ;


    // // add, remove, and retrieve rules associated with this keyspace
    // void
    // add_rule(const rule& r) {rules += r ;}
    // void
    // remove_rule(const rule& r) {rules -= r ;}
    // const ruleSet&
    // get_rules() const {return rules ;}

    // generates a new key partition based on the internal rules
    // returns "true" if key_ptn has been updated,
    //         "false" if not
    virtual bool
    update_key_ptn() = 0 ;
    // this should not be called, instead the CopyKeySpace<> should
    // be called, this is modeled after the rule_impl implementation
    virtual KeySpaceP new_keyspace() const ;
    // this one returns a clone of the current keyspace
    //virtual KeySpaceP clone_keyspace() const ;
  } ; // end of class KeySpace

  typedef KeySpace::KeySpaceP KeySpaceP ;

  // ORB keyspace type definition
  // mainly, the ORB keyspace adds several routines that
  // perform the ORB partition and also it defines the
  // critical data-structure
  class OrbKeySpace: public KeySpace {
  protected:
    // first the critical data-structure:
    // a vector of <coord> that describes the
    // spatial location of keys corresponding to
    // the vector <key id>
    std::vector<Loci::vector3d<float> > coords ;
    std::vector<Key> coord_keys ;
  protected:
    // the function that partitions the space according
    // to the critical structure, it returns a keyset partition

    // NOTE: after this function, the critical data-structures
    // "coords" and "keys" may have been changed and may not
    // validly correspond to the current key-data distribution.
    // therefore, every time before using this function, it is
    // necessary to refill the "coords" and "keys" data-structures
    std::vector<KeySet>
    orb_partition() ;
  } ;

  // the following codes provide facilities to set up user defined
  // KeySpaces. the codes mirror the currently design and implementation
  // of rule list and the register rule mechanism.
  template<class TCopyKeySpace>
  class CopyKeySpace: public TCopyKeySpace {
  public:
    virtual KeySpaceP new_keyspace() const {
      return new CopyKeySpace<TCopyKeySpace> ;
    }
    // virtual KeySpaceP clone_keyspace() const {
    //   return new CopyKeySpace<TCopyKeySpace>(*this) ;
    // }
  } ;
  
  class RegisterKeySpaceType {
  public:
    virtual ~RegisterKeySpaceType() {}
    virtual KeySpaceP get_space() const = 0 ;
  } ;
  
  class RegisterKeySpaceList ;

  class KeySpaceList {
  public:
    class KeySpaceListIterator ;
    friend class KeySpaceListIterator ;
    friend class RegisterKeySpaceList ;
    class KeySpaceListEnt {
    public:
      KeySpaceListEnt(RegisterKeySpaceType* p, KeySpaceListEnt* nxt)
        :rr(p), next(nxt) {}
      RegisterKeySpaceType* rr ;
      KeySpaceListEnt* next ;
    } ;
    KeySpaceListEnt* list ;
  public:
    class KeySpaceListIterator {
      KeySpaceListEnt* p ;
    public:
      KeySpaceListIterator(KeySpaceListEnt* ptr):p(ptr) {}
      KeySpaceP operator*() {return p->rr->get_space() ;}
      KeySpaceListIterator& operator++() {
        p = p->next ;
        return *this ;
      }
      KeySpaceListIterator operator++(int) {
        KeySpaceListIterator tmp(p) ;
        p = p->next ;
        return tmp ;
      }
      KeySpaceListEnt* get_p() {return p ;}
      bool operator==(const KeySpaceListIterator& i) {return i.p == p ;}
      bool operator!=(const KeySpaceListIterator& i) {return i.p != p ;}
    } ;
    typedef KeySpaceListIterator Iterator ;

    KeySpaceList() {list = 0 ;}
    ~KeySpaceList() ;

    void push_space(RegisterKeySpaceType* rr) ;
    Iterator begin() {return Iterator(list) ;}
    Iterator end() {return Iterator(0) ;}
    void clear() ;
    void copy_space_list(const KeySpaceList& sl) ;
    void copy_space_list(const RegisterKeySpaceList& sl) ;
    KeySpaceList(const KeySpaceList& x) {
      list = 0 ;
      copy_space_list(x) ;
    }
    KeySpaceList& operator=(const KeySpaceList& x) {
      if(&x != this) {
        clear() ;
        copy_space_list(x) ;
      }
      return *this ;
    }
  } ;

  class RegisterKeySpaceList: public KeySpaceList {
  public:
    static KeySpaceListEnt* global_list ;
    RegisterKeySpaceList() {}
    ~RegisterKeySpaceList() ;
    void clear() ;
    bool empty() ;
    void push_space(RegisterKeySpaceType* p) ;
    Iterator begin() {return Iterator(global_list) ;}
    Iterator end() {return Iterator(0) ;}
  } ;

  extern RegisterKeySpaceList register_key_space_list ;
  extern KeySpaceList global_key_space_list ;

  template<class T>
  class register_key_space: public RegisterKeySpaceType {
  public:
    register_key_space() {register_key_space_list.push_space(this) ;}
    virtual KeySpaceP get_space() const {return new CopyKeySpace<T> ;}
  } ;

} // end of namespace Loci
#endif

#endif
