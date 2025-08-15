//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef READGRID_H
#define READGRID_H

#include <Loci.h>
#include <flowTypes.h>
#include <vector>
#include <string>

namespace flowPsi {
  void read_grid(fact_db &facts, const rule_db &rdb, string filename,
                 std::vector<pair<std::string,std::string> > &command_vars,
                 string gridfile, bool dryrun) ;

  void find_mind_noslip(fact_db &facts) ;
  void find_mind_surf(fact_db &facts) ;
  void find_mind_surf_node(fact_db &facts) ;

  void create_cell_stencil(fact_db &facts) ;

  bool check_boundary_conditions(fact_db &facts) ;
  bool check_scalar_units(const options_list &o, std::string option,
                          std::string unit) ;
  bool check_vector_units(const options_list &ol,std::string vname,
                          std::string units) ;
}

namespace flowPsi {

  class BC_Check : public Loci::CPTR_type {
  public:
    // Returns string containing comma separated list of boundary
    // condition tags that are checked by this checker
    virtual std::string boundaryConditions() = 0 ;
    // Returns string containing variables that are tested by this
    // checker.
    virtual std::string variablesChecked() = 0 ;
    // Check the boundary condition
    virtual bool checkOptions(const options_list &bc_options) = 0 ;
    // Emit Error Message
    virtual std::ostream  &ErrorMessage(std::ostream &s) = 0 ;
  } ;

  typedef Loci::CPTR<BC_Check> BC_implP ;
  
  class register_BC_type {
  public:
    virtual ~register_BC_type() {}
    virtual BC_implP get_BC() const = 0 ;
  } ; 
  
  class register_BC_impl_list ;
  
  class BC_impl_list {
  public:
    class BC_list_iterator ;
    friend class BC_list_iterator ;
    friend class register_BC_impl_list ;
    class BC_list_ent {
    public:
      BC_list_ent(register_BC_type *p, BC_list_ent *nxt) :
        rr(p), next(nxt) {}
      register_BC_type *rr ;
      BC_list_ent *next ;
    } ;
    BC_list_ent *list ;
  public:
    class BC_list_iterator {
      BC_list_ent *p ;
    public:
      BC_list_iterator(BC_list_ent *ptr) : p(ptr) {}
      BC_implP operator*() { return p->rr->get_BC() ; }
      BC_list_iterator &operator++() {
        p = p->next ;
        return *this ;
      }
      BC_list_iterator operator++(int ) {
        BC_list_iterator tmp(p) ;
        p = p->next ;
        return tmp ;
      }
      BC_list_ent* get_p() { return p; } 
      bool operator==(const BC_list_iterator &i) { return i.p == p ; }
      bool operator!=(const BC_list_iterator &i) { return i.p != p ; }
    } ;
    typedef BC_list_iterator iterator ;
    
    BC_impl_list() {list = 0 ; }
    ~BC_impl_list() ;
    
    void push_BC(register_BC_type *rr) ;
    iterator begin() { return iterator(list) ; }
    iterator end() { return iterator(0) ; }
    void clear() ;
    void copy_BC_list(const BC_impl_list &rl) ;
    void copy_BC_list(const register_BC_impl_list &rl) ;
    BC_impl_list(const BC_impl_list &x) {
      list = 0 ;
      copy_BC_list(x) ;
    }
    BC_impl_list &operator=(const BC_impl_list &x) {
      list = 0 ;
      copy_BC_list(x) ;
      return *this ;
    }
  } ;
  class register_BC_impl_list : public BC_impl_list {
  public:
    static BC_list_ent *global_list ;
    register_BC_impl_list() {}
    ~register_BC_impl_list() ;
    void clear() ;
    bool empty() ;
    void push_BC(register_BC_type *p) ;
    iterator begin() { return iterator(global_list) ; }
    iterator end() { return iterator(0) ; }
  } ;
  extern register_BC_impl_list register_BC_list ;    
  //  extern BC_impl_list global_BC_list ;    

  template<class T> class register_BC : public register_BC_type {
  public:
    register_BC() { register_BC_list.push_BC(this) ; }
    virtual BC_implP get_BC() const { return BC_implP(new T) ; }
  } ;
  

  class gridPostStep : public Loci::CPTR_type {
  public:
    virtual void processData(fact_db &facts, rule_db &rdb) const = 0 ;
  } ;

  class gridPostList {
  public:
    struct gridPostListEnt {
      Loci::CPTR<gridPostStep> entry ;
      gridPostListEnt *next ;
    } ;
    static gridPostListEnt *postlist  ;
    static void insertStep(Loci::CPTR<gridPostStep> cp) {
      gridPostListEnt *p = new gridPostListEnt ;
      p->next = gridPostList::postlist ;
      p->entry = cp ;
      gridPostList::postlist = p ;
    }
  } ;
    
  template<class T> class registerGridPostProcess {
  public:
    registerGridPostProcess() {
      Loci::CPTR<gridPostStep> p = new T ;
      gridPostList::insertStep(p) ;
    }
  } ;
      

    
  bool check_scalar_units(const options_list &o, std::string option,
                          std::string unit) ;
  bool check_vector_units(const options_list &ol,std::string vname,
                          std::string units) ;

}

#endif
