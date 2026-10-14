// Copyright 2026, Streamline Numerics, Inc.
// Distributed under LGPL-3.0-or-later; see COPYING.LESSER.
#include <module_interface.h>
#include <mod_db.h>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <typeinfo>

namespace Loci {
  namespace {
    bool builtin(variable name) {
      const std::string &s = name.get_info().name ;
      return name.is_time_variable() || s == "UNIVERSE" || s == "EMPTY"
        || s == "OUTPUT" ;
    }

    bool identifier(const std::string &name) {
      if(name.empty() || !(std::isalpha(static_cast<unsigned char>(name[0]))
                           || name[0] == '_'))
        return false ;
      for(unsigned char c : name)
        if(!std::isalnum(c) && c != '_')
          return false ;
      return true ;
    }

    variable base_name(variable value) {
      variable::info info = value.get_info() ;
      info.time_id = time_ident() ;
      info.offset = 0 ;
      info.assign = false ;
      info.priority.clear() ;
      return variable(info) ;
    }

    ruleSet all_rules(const rule_db &rules) {
      ruleSet result = rules.all_rules() ;
      result += rules.get_default_rules() ;
      result += rules.get_optional_rules() ;
      return result ;
    }

    variableSet declarations(const std::vector<std::string> &names,
                             const std::string &origin) {
      variableSet result ;
      for(const auto &name : names) {
        if(!identifier(name) || builtin(variable(name)))
          throw module_binding_error(origin + ": invalid interface name "
                                     + name) ;
        const variable value(name) ;
        if(result.inSet(value))
          throw module_binding_error(origin + ": repeated declaration "
                                     + name) ;
        result += value ;
      }
      return result ;
    }

    void check_type(variable name, storeRepP expected, storeRepP actual) {
      if(expected == 0 || actual == 0)
        return ;
      // Compare empty frozen representations: a dynamic map/store can satisfy
      // the corresponding static declaration without copying numerical data.
      storeRepP declared = expected->new_store(EMPTY)->freeze() ;
      storeRepP supplied = actual->new_store(EMPTY)->freeze() ;
      if(typeid(*declared->getRep()) != typeid(*supplied->getRep()))
        throw module_binding_error("module interface: incompatible type for "
                                   + name.str()) ;
    }

    // Compare known concrete declarations; rule execution retains its normal
    // type checks for types that become known only through instantiation.
    void check_binding_types(const rule_db &rules, fact_db *facts) {
      std::map<variable, storeRepP> types ;
      const auto &bindings = rules.get_module_bindings() ;
      variableSet destinations ;
      for(const auto &input : bindings.entries())
        destinations += input.second ;
      for(rule r : all_rules(rules)) {
        rule_implP impl = r.get_rule_implP() ;
        if(impl == 0)
          continue ;
        for(variable value : impl->get_var_list()) {
          if(!value.get_arg_list().empty() || builtin(value))
            continue ;
          const variable name = bindings.resolve(base_name(value)) ;
          if(!destinations.inSet(name))
            continue ;
          storeRepP data = impl->get_store(value) ;
          if(data == 0)
            continue ;
          auto found = types.find(name) ;
          if(found != types.end())
            check_type(name, found->second, data) ;
          else
            types[name] = data ;
        }
      }
      if(facts != nullptr)
        for(const auto &type : types) {
          storeRepP data = facts->get_variable(type.first) ;
          if(data == 0)
            data = facts->get_variable_type(type.first) ;
          check_type(type.first, type.second, data) ;
        }
    }
  }

  variable qualify_module_variable(variable name,
                                   const std::vector<std::string> &scope) {
    if(builtin(name))
      return name ;
    variable::info info = name.get_info() ;
    if(info.namespac.empty())
      info.namespac = scope ;
    for(int &arg : info.v_ids)
      arg = qualify_module_variable(variable(arg), scope).ident() ;
    return variable(info) ;
  }

  module_context::module_context(fact_db &data, const module_bindings &table,
      const std::vector<std::string> &space, const std::string &owner)
    : facts(data), bindings(table), scope(space), origin(owner) {
    if(!facts.get_namespace().empty())
      throw module_binding_error(origin + ": setup requires unscoped fact_db") ;
  }

  variable module_context::resolve(variable name) const {
    return bindings.resolve(qualify_module_variable(name, scope)) ;
  }

  variable module_context::write_name(variable name) const {
    const variable local = qualify_module_variable(name, scope) ;
    const variable canonical = bindings.resolve(local) ;
    if(bindings.is_input(local) || bindings.is_import(canonical, origin))
      throw module_binding_error(origin + ": cannot write imported input "
                                 + name.str() + " -> " + canonical.str()) ;
    const auto &space = canonical.get_info().namespac ;
    if(space.size() < scope.size()
       || !std::equal(scope.begin(), scope.end(), space.begin()))
      throw module_binding_error(origin + ": cannot write outside instance: "
                                 + canonical.str()) ;
    return canonical ;
  }

  storeRepP module_context::read(const std::string &name) const {
    const variable canonical = resolve(name) ;
    storeRepP data = facts.get_variable(canonical) ;
    if(data == 0)
      throw module_binding_error(origin + ": setup input unavailable: " + name
                                 + " -> " + canonical.str()) ;
    return data ;
  }

  void module_context::create(const std::string &name, storeRepP data) {
    const variable canonical = write_name(variable(name)) ;
    if(data == 0 || facts.get_variable(canonical) != 0)
      throw module_binding_error(origin + ": invalid create: "
                                 + canonical.str()) ;
    check_type(canonical, facts.get_variable_type(canonical), data) ;
    facts.create_fact(canonical, data) ;
  }

  void module_context::update(const std::string &name, storeRepP data) {
    const variable canonical = write_name(variable(name)) ;
    storeRepP current = facts.get_variable(canonical) ;
    if(data == 0 || current == 0)
      throw module_binding_error(origin + ": invalid update: "
                                 + canonical.str()) ;
    check_type(canonical, current, data) ;
    facts.update_fact(canonical, data) ;
  }

  void module_context::erase(const std::string &name) {
    facts.remove_variable(write_name(variable(name))) ;
  }

  void module_composition::add(const std::string &library,
      const std::vector<std::string> &scope,
      const module_bindings::binding_map &inputs,
      const std::string &vars_file) {
    if(state != DECLARING)
      throw module_binding_error("module composition: declarations frozen") ;
    if(library.empty() || scope.empty())
      throw module_binding_error("module composition: library/scope required") ;
    for(const auto &part : scope)
      if(!identifier(part))
        throw module_binding_error("invalid module namespace component: "
                                   + part) ;
    for(const auto &instance : instances)
      if(instance.scope == scope) {
        if(instance.library == library && instance.inputs == inputs
           && instance.vars_file == vars_file)
          return ;
        throw module_binding_error("conflicting declaration: "
                                   + instance.origin) ;
      }
    instance_data instance ;
    instance.library = library ;
    instance.scope = scope ;
    instance.inputs = inputs ;
    instance.vars_file = vars_file ;
    instance.origin = library + "[" ;
    for(const auto &part : scope)
      instance.origin += (instance.origin.back() == '[' ? "" : "@") + part ;
    instance.origin += "]" ;
    instances.push_back(instance) ;
  }

  void module_composition::prepare(const rule_db &host_rules) {
    if(state != DECLARING)
      throw module_binding_error("module composition: already prepared") ;
    try {
      if(!host_rules.get_module_bindings().empty())
        throw module_binding_error("host rules already have module bindings") ;
      rule_db combined(host_rules) ;
      module_bindings::binding_map table ;
      module_bindings::import_map imports ;
      for(auto &instance : instances) {
        mod library ;
        mod::mod_info definition = library.get_info(instance.library) ;
        typedef const module_interface_v1 *(*describe_function)() ;
        auto describe = reinterpret_cast<describe_function>(
          dlsym(definition.m_library, "loci_module_interface_v1")) ;
        if(describe == nullptr)
          throw module_binding_error(instance.origin
                                     + ": missing interface v1") ;
        const module_interface_v1 *interface = describe() ;
        if(interface == nullptr)
          throw module_binding_error(instance.origin + ": null interface") ;
        instance.interface = *interface ;
        instance.initialize = reinterpret_cast<init_function>(
          dlsym(definition.m_library, "loci_init_module_v1")) ;
        const variableSet inputs = declarations(interface->inputs,
                                                instance.origin) ;
        const variableSet outputs = declarations(interface->outputs,
                                                 instance.origin) ;
        if((inputs & outputs) != EMPTY)
          throw module_binding_error(instance.origin
                                     + ": input/output overlap") ;
        for(const auto &input : instance.inputs) {
          if(!inputs.inSet(input.first))
            throw module_binding_error(instance.origin + ": unknown input "
                                       + input.first.str()) ;
          const variable local = qualify_module_variable(input.first,
                                                         instance.scope) ;
          table.insert(std::make_pair(local, input.second)) ;
          imports[instance.origin] += input.second ;
        }
        for(variable input : inputs)
          if(instance.inputs.find(input) == instance.inputs.end())
            throw module_binding_error(instance.origin + ": missing binding "
                                       + input.str()) ;

        for(auto entry = definition.loaded_rule_list.begin();
            entry != definition.loaded_rule_list.end(); ++entry) {
          if(entry.get_p()->rr->is_module_rule())
            throw module_binding_error(instance.origin
              + ": legacy module declarations require explicit host adapters") ;
          rule_implP impl = *entry ;
          if(!impl->get_info().dynamic_constraints.empty())
            throw module_binding_error(instance.origin
              + ": dynamic constraints are not supported") ;
          std::map<variable, variable> names ;
          for(variable name : impl->get_var_list())
            names[name] = qualify_module_variable(name, instance.scope) ;
          impl->rename_vars(names) ;
          rule r(impl) ;
          instance.rules.add_rule(r) ;
          combined.add_rule(r) ;
          combined.add_module_origin(r, instance.origin) ;
        }
      }
      combined.set_module_bindings(module_bindings(table, imports)) ;
      validate_module_bindings(combined, combined.get_module_bindings()) ;
      check_binding_types(combined, nullptr) ;
      prepared_rules = combined ;
      state = PREPARED ;
    } catch(...) {
      state = FAILED ;
      throw ;
    }
  }

  void module_composition::initialize(fact_db &facts) {
    if(state == READY && initialized_facts == &facts)
      return ;
    if(state != PREPARED)
      throw module_binding_error("module composition: cannot initialize") ;
    state = INITIALIZING ;
    try {
      if(!facts.get_namespace().empty())
        throw module_binding_error("module composition needs unscoped facts") ;
      const auto &bindings = prepared_rules.get_module_bindings() ;
      for(const auto &input : bindings.entries())
        if(facts.get_extensional_facts().inSet(input.first))
          throw module_binding_error("fact already initializes imported input "
                                     + input.first.str()) ;
      check_binding_types(prepared_rules, &facts) ;
      facts.setupDefaults(prepared_rules) ;
      for(auto &instance : instances) {
        if(instance.vars_file.empty())
          continue ;
        std::ifstream deck(instance.vars_file.c_str()) ;
        if(!deck)
          throw module_binding_error(instance.origin + ": cannot open "
                                     + instance.vars_file) ;
        module_context context(facts, bindings,
                               instance.scope, instance.origin) ;
        try {
          facts.read_vars(deck, instance.rules, [&](variable name) {
            return context.write_name(name) ;
          }) ;
        } catch(const module_binding_error &error) {
          throw module_binding_error(instance.vars_file + ": " + error.what()) ;
        } catch(const BasicException &error) {
          std::ostringstream message ;
          message << instance.origin << ": " << instance.vars_file << ": " ;
          error.Print(message) ;
          throw module_binding_error(message.str()) ;
        }
      }
      for(auto &instance : instances) {
        module_context context(facts, bindings,
                               instance.scope, instance.origin) ;
        if(instance.initialize != nullptr)
          instance.initialize(context) ;
      }
      check_binding_types(prepared_rules, &facts) ;
      initialized_facts = &facts ;
      state = READY ;
    } catch(...) {
      state = FAILED ;
      throw ;
    }
  }

  const rule_db &module_composition::rules() const {
    if(state != READY)
      throw module_binding_error("module composition: not initialized") ;
    return prepared_rules ;
  }
}
