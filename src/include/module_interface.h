// Copyright 2026, Streamline Numerics, Inc.
// Distributed under LGPL-3.0-or-later; see COPYING.LESSER.
#ifndef LOCI_MODULE_INTERFACE_H
#define LOCI_MODULE_INTERFACE_H

#include <fact_db.h>
#include <rule.h>
#include <iosfwd>

namespace Loci {
  /// Returned by the optional C symbol loci_module_interface_v1.
  /// Names are local concrete identifiers. All inputs require host bindings.
  struct module_interface_v1 {
    std::vector<std::string> inputs, outputs ;
  } ;

  /// Qualify local names and scope local iteration labels/counters per instance.
  /// Explicit external references retain their names and time annotations.
  variable qualify_module_variable(variable name,
                                   const std::vector<std::string> &scope) ;

  /// Instance identity of an explicit local iteration (before substitution).
  time_ident module_time(const std::vector<std::string> &scope,
                         const time_ident &local_time) ;

  /// Setup access for one instance. The context is valid during its callback.
  /// Read imported data through const Loci containers; this is not a sandbox
  /// against C++ code deliberately mutating a returned representation.
  class module_context {
  public:
    variable resolve(variable name) const ;
    variable resolve(const std::string &name) const {
      return resolve(variable(name)) ;
    }
    /// Look up existing data without scheduling; return null when absent.
    storeRepP get_variable(variable name) const ;
    storeRepP get_variable(const std::string &name) const {
      return get_variable(variable(name)) ;
    }
    /// Required setup data; absence is an interface error.
    storeRepP read(variable name) const ;
    storeRepP read(const std::string &name) const { return read(variable(name)) ; }
    void create(variable name, storeRepP data) ;
    void create(const std::string &name, storeRepP data) { create(variable(name), data) ; }
    void update(variable name, storeRepP data) ;
    void update(const std::string &name, storeRepP data) { update(variable(name), data) ; }
    void erase(variable name) ;
    void erase(const std::string &name) { erase(variable(name)) ; }
    const std::string &instance() const { return origin.instance ; }

  private:
    friend class module_composition ;
    module_context(fact_db &facts, const module_bindings &bindings,
                   const std::vector<std::string> &scope,
                   const module_origin &origin) ;
    variable write_name(variable name) const ;
    fact_db &facts ;
    const module_bindings &bindings ;
    const std::vector<std::string> &scope ;
    const module_origin &origin ;
  } ;

  /// Explicit host composition: declare, prepare, initialize, then query.
  /// No implicit dependencies, nested setup, or application scheduling.
  class module_composition {
  public:
    module_composition() = default ;
    module_composition(const module_composition &) = delete ;
    module_composition &operator=(const module_composition &) = delete ;

    /// Repeating an identical declaration before preparation is harmless.
    /// Instance roots must not overlap; siblings and grouping prefixes are valid.
    /// Optional output bindings precede the optional deck filename.
    void add(const std::string &library, const std::vector<std::string> &scope,
             const module_bindings::binding_map &inputs,
             const module_bindings::binding_map &outputs = {},
             const std::string &vars_file = "") ;

    /// Validate the entire composition before any deck or callback runs.
    void prepare(const rule_db &host_rules = rule_db()) ;

    /// Defaults and decks for all instances precede callbacks in add order.
    /// Repeating initialization on the same database does not rerun callbacks.
    /// An initialization failure is terminal; arbitrary side effects roll back
    /// only if the application itself provides that behavior.
    void initialize(fact_db &facts) ;
    const rule_db &rules() const ;

    /// Inspect final identities and relocated input requests without running setup.
    void print_connections(std::ostream &out) const ;

  private:
    typedef void (*init_function)(module_context &) ;
    struct instance_data {
      module_origin origin ;
      std::string vars_file ;
      std::vector<std::string> scope ;
      module_bindings::binding_map inputs, outputs ;
      module_interface_v1 interface ;
      rule_db declarations ;
      init_function initialize = nullptr ;
    } ;
    struct type_declaration {
      storeRepP type ;
      variable local ;
      rule provider ;
    } ;
    void collect_binding_types() ;
    void check_fact_types(fact_db &facts) const ;
    std::string declaration_location(const type_declaration &type, variable name) const ;
    enum phase { DECLARING, PREPARED, INITIALIZING, READY, FAILED } ;
    phase state = DECLARING ;
    std::vector<instance_data> instances ;
    rule_db prepared_rules ;
    std::map<variable, type_declaration> connected_types ;
    fact_db *initialized_facts = nullptr ;
  } ;
}

#endif
