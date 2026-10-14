// Copyright 2026, Streamline Numerics, Inc.
// Distributed under LGPL-3.0-or-later; see COPYING.LESSER.
#ifndef LOCI_MODULE_INTERFACE_H
#define LOCI_MODULE_INTERFACE_H

#include <fact_db.h>
#include <rule.h>

namespace Loci {
  /// Returned by the optional C symbol loci_module_interface_v1.
  /// Names are local concrete identifiers. All inputs require host bindings.
  struct module_interface_v1 {
    std::vector<std::string> inputs, outputs ;
  } ;

  /// Qualify local names per variable; retain explicit names and built-ins.
  variable qualify_module_variable(variable name,
                                   const std::vector<std::string> &scope) ;

  /// Setup access for one instance. The context is valid during its callback.
  /// Read imported data through const Loci containers; this is not a sandbox
  /// against C++ code deliberately mutating a returned representation.
  class module_context {
  public:
    variable resolve(variable name) const ;
    variable resolve(const std::string &name) const {
      return resolve(variable(name)) ;
    }
    storeRepP read(const std::string &name) const ;
    void create(const std::string &name, storeRepP data) ;
    void update(const std::string &name, storeRepP data) ;
    void erase(const std::string &name) ;
    const std::string &instance() const { return origin ; }

  private:
    friend class module_composition ;
    module_context(fact_db &facts, const module_bindings &bindings,
                   const std::vector<std::string> &scope,
                   const std::string &origin) ;
    variable write_name(variable name) const ;
    fact_db &facts ;
    const module_bindings &bindings ;
    const std::vector<std::string> &scope ;
    const std::string &origin ;
  } ;

  /// Explicit host composition: declare, prepare, initialize, then query.
  /// No implicit dependencies, nested setup, or application scheduling.
  class module_composition {
  public:
    module_composition() = default ;
    module_composition(const module_composition &) = delete ;
    module_composition &operator=(const module_composition &) = delete ;

    /// Repeating an identical declaration before preparation is harmless.
    void add(const std::string &library, const std::vector<std::string> &scope,
             const module_bindings::binding_map &inputs,
             const std::string &vars_file = "") ;

    /// Validate the entire composition before any deck or callback runs.
    void prepare(const rule_db &host_rules = rule_db()) ;

    /// Defaults and decks for all instances precede callbacks in add order.
    /// Repeating initialization on the same database does not rerun callbacks.
    /// An initialization failure is terminal; arbitrary side effects roll back
    /// only if the application itself provides that behavior.
    void initialize(fact_db &facts) ;
    const rule_db &rules() const ;

  private:
    typedef void (*init_function)(module_context &) ;
    struct instance_data {
      std::string library, origin, vars_file ;
      std::vector<std::string> scope ;
      module_bindings::binding_map inputs ;
      module_interface_v1 interface ;
      rule_db rules ;
      init_function initialize = nullptr ;
    } ;
    enum phase { DECLARING, PREPARED, INITIALIZING, READY, FAILED } ;
    phase state = DECLARING ;
    std::vector<instance_data> instances ;
    rule_db prepared_rules ;
    fact_db *initialized_facts = nullptr ;
  } ;
}

#endif
