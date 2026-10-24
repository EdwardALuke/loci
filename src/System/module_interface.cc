// Copyright 2026, Streamline Numerics, Inc.
// Distributed under LGPL-3.0-or-later; see COPYING.LESSER.
#include <module_interface.h>
#include <mod_db.h>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <typeinfo>
#ifdef __GNUG__
#include <cxxabi.h>
#include <cstdlib>
#endif

namespace Loci {
  namespace {
    std::string time_scope(const std::vector<std::string> &scope) {
      // Encode separators too: a@b and a_b remain distinct. Hex contains no
      // underscore-delimited formal label that substitution could replace.
      std::string result = "ns" ;
      const char *hex = "0123456789abcdef" ;
      for(size_t i=0;i<scope.size();i++) {
        if(i != 0) result += "40" ;
        for(unsigned char c : scope[i]) {
          result += hex[c >> 4] ;
          result += hex[c & 15] ;
        }
      }
      return result ;
    }

    time_ident scope_time(const time_ident &time, const std::string &prefix) {
      if(time == time_ident()) return time ;
      return time_ident(scope_time(time.parent(), prefix),
                        prefix + "_" + time.level_name()) ;
    }

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

    std::string type_name(const std::type_info &type) {
#ifdef __GNUG__
      int status = 0 ;
      char *name = abi::__cxa_demangle(type.name(), nullptr, nullptr, &status) ;
      const std::string result = status == 0 && name != nullptr ? name : type.name() ;
      std::free(name) ;
      return result ;
#else
      return type.name() ;
#endif
    }

    // A dynamic map/store can satisfy the corresponding static declaration.
    // These empty representations carry type information, not numerical data.
    storeRepP frozen_type(storeRepP data) {
      return data == 0 ? storeRepP() : data->new_store(EMPTY)->freeze() ;
    }

    bool incompatible_types(storeRepP expected, storeRepP actual) {
      return expected != 0 && actual != 0
        && typeid(*expected->getRep()) != typeid(*actual->getRep()) ;
    }

    module_binding_error type_error(variable name, storeRepP expected, storeRepP actual,
                                    const std::string &expected_at, const std::string &actual_at) {
      std::ostringstream message ;
      message << "incompatible type for " << name
              << "\n  " << expected_at << ": " << type_name(typeid(*expected->getRep()))
              << "\n  " << actual_at << ": " << type_name(typeid(*actual->getRep())) ;
      return module_binding_error(message.str()) ;
    }
  }

  variable qualify_module_variable(variable name,
                                   const std::vector<std::string> &scope) {
    variable::info info = name.get_info() ;
    if(info.tvar) {
      info.name = time_scope(scope) + "_" + info.name ;
      info.time_id = module_time(scope, info.time_id) ;
      return variable(info) ;
    }
    if(info.namespac.empty()) {
      if(!builtin(name)) info.namespac = scope ;
      info.time_id = module_time(scope, info.time_id) ;
    }
    for(int &arg : info.v_ids)
      arg = qualify_module_variable(variable(arg), scope).ident() ;
    return variable(info) ;
  }

  time_ident module_time(const std::vector<std::string> &scope,
                         const time_ident &local_time) {
    return scope_time(local_time, time_scope(scope)) ;
  }

  module_context::module_context(fact_db &data, const module_bindings &table,
      const std::vector<std::string> &space, const module_origin &owner)
    : facts(data), bindings(table), scope(space), origin(owner) {
    if(!facts.get_namespace().empty())
      throw module_binding_error(origin.description() + ": setup requires unscoped fact_db") ;
  }

  variable module_context::resolve(variable name) const {
    return bindings.resolve(qualify_module_variable(name, scope)) ;
  }

  variable module_context::write_name(variable name) const {
    const variable local = qualify_module_variable(name, scope) ;
    const variable canonical = bindings.resolve(local) ;
    if(bindings.is_input(local) || bindings.is_import(canonical, origin.instance))
      throw module_binding_error(origin.description() + ": cannot write imported input "
                                 + name.str() + " -> " + canonical.str()) ;
    const auto &space = canonical.get_info().namespac ;
    if(!bindings.is_export(canonical, origin.instance)
       && (space.size() < scope.size()
           || !std::equal(scope.begin(), scope.end(), space.begin())))
      throw module_binding_error(origin.description() + ": cannot write outside instance: "
                                 + canonical.str()) ;
    return canonical ;
  }

  storeRepP module_context::get_variable(variable name) const {
    return facts.get_variable(resolve(name)) ;
  }

  storeRepP module_context::read(variable name) const {
    const variable canonical = resolve(name) ;
    storeRepP data = facts.get_variable(canonical) ;
    if(data == 0)
      throw module_binding_error(origin.description() + ": setup input unavailable: " + name.str()
                                 + " -> " + canonical.str()) ;
    return data ;
  }

  void module_context::create(variable name, storeRepP data) {
    const variable canonical = write_name(name) ;
    if(data == 0 || facts.get_variable(canonical) != 0)
      throw module_binding_error(origin.description() + ": invalid create: "
                                 + canonical.str()) ;
    const storeRepP expected = frozen_type(facts.get_variable_type(canonical)) ;
    const storeRepP actual = frozen_type(data) ;
    if(incompatible_types(expected, actual))
      throw type_error(canonical, expected, actual,
                       origin.description() + ": declared " + canonical.str(), "create " + name.str()) ;
    facts.create_fact(canonical, data) ;
  }

  void module_context::update(variable name, storeRepP data) {
    const variable canonical = write_name(name) ;
    storeRepP current = facts.get_variable(canonical) ;
    if(data == 0 || current == 0)
      throw module_binding_error(origin.description() + ": invalid update: "
                                 + canonical.str()) ;
    const storeRepP expected = frozen_type(current), actual = frozen_type(data) ;
    if(incompatible_types(expected, actual))
      throw type_error(canonical, expected, actual,
                       origin.description() + ": current " + canonical.str(), "update " + name.str()) ;
    facts.update_fact(canonical, data) ;
  }

  void module_context::erase(variable name) {
    facts.remove_variable(write_name(name)) ;
  }

  void module_composition::add(const std::string &library,
      const std::vector<std::string> &scope,
      const module_bindings::binding_map &inputs,
      const module_bindings::binding_map &outputs,
      const std::string &vars_file) {
    if(state != DECLARING)
      throw module_binding_error("module composition: declarations frozen") ;
    if(library.empty() || scope.empty())
      throw module_binding_error("module composition: library/scope required") ;
    for(const auto &part : scope)
      if(!identifier(part))
        throw module_binding_error("invalid module namespace component: "
                                   + part) ;
    module_origin origin = {"", library} ;
    for(const auto &part : scope)
      origin.instance += (origin.instance.empty() ? "" : "@") + part ;
    for(const auto &instance : instances) {
      if(instance.scope == scope) {
        if(instance.origin.library == library && instance.inputs == inputs
           && instance.outputs == outputs
           && instance.vars_file == vars_file)
          return ;
        throw module_binding_error("module composition declare: conflicting declaration: "
                                   + instance.origin.description() + " and " + origin.description()) ;
      }
      // Separate instances cannot own ancestor/descendant paths. Grouping
      // prefixes and sibling paths remain available without parent authority.
      const size_t common = std::min(scope.size(), instance.scope.size()) ;
      if(std::equal(scope.begin(), scope.begin() + common, instance.scope.begin()))
        throw module_binding_error("module composition declare: overlapping instance roots: "
                                   + instance.origin.description() + " and " + origin.description()) ;
    }
    instance_data instance ;
    instance.scope = scope ;
    instance.inputs = inputs ;
    instance.outputs = outputs ;
    instance.vars_file = vars_file ;
    instance.origin = origin ;
    instances.push_back(instance) ;
  }

  std::string module_composition::declaration_location(const type_declaration &type, variable name) const {
    std::ostringstream location ;
    for(const auto &origin : prepared_rules.get_module_origins(type.provider))
      location << origin.description() << ": " ;
    location << type.local << " -> " << name << ", rule " << type.provider ;
    const std::string file = type.provider.get_rule_implP()->get_fileloc() ;
    if(!file.empty())
      location << " at " << file ;
    return location.str() ;
  }

  // Collect concrete declarations once. Types revealed by instantiation still
  // use normal Loci checks; diagnostics are formatted only on failure.
  void module_composition::collect_binding_types() {
    const auto &bindings = prepared_rules.get_module_bindings() ;
    variableSet destinations ;
    for(const auto &input : bindings.input_entries())
      destinations += input.second ;
    for(const auto &output : bindings.output_entries())
      destinations += output.second ;
    for(const auto &instance : instances)
      for(const auto &output : instance.interface.outputs)
        destinations += bindings.resolve(qualify_module_variable(variable(output), instance.scope)) ;
    for(rule r : all_rules(prepared_rules)) {
      rule_implP impl = r.get_rule_implP() ;
      if(impl == 0)
        continue ;
      for(variable value : impl->get_var_list()) {
        if(!value.get_arg_list().empty() || builtin(value))
          continue ;
        const variable name = bindings.resolve(base_name(value)) ;
        if(!destinations.inSet(name))
          continue ;
        const storeRepP type = frozen_type(impl->get_store(value)) ;
        if(type == 0)
          continue ;
        const type_declaration declaration = {type, value, r} ;
        const auto inserted = connected_types.emplace(name, declaration) ;
        const auto &known = inserted.first->second ;
        if(!inserted.second && incompatible_types(known.type, type))
          throw type_error(name, known.type, type,
                           declaration_location(known, name), declaration_location(declaration, name)) ;
      }
    }
  }

  void module_composition::check_fact_types(fact_db &facts) const {
    for(const auto &type : connected_types) {
      storeRepP data = facts.get_variable(type.first) ;
      if(data == 0)
        data = facts.get_variable_type(type.first) ;
      const storeRepP actual = frozen_type(data) ;
      if(incompatible_types(type.second.type, actual))
        throw type_error(type.first, type.second.type, actual,
                         declaration_location(type.second, type.first), "host fact " + type.first.str()) ;
    }
  }

  void module_composition::prepare(const rule_db &host_rules) {
    if(state != DECLARING)
      throw module_binding_error("module composition: already prepared") ;
    try {
      if(!host_rules.get_module_bindings().empty())
        throw module_binding_error("host rules already have module bindings") ;
      rule_db combined(host_rules) ;
      module_bindings::binding_map table, published ;
      module_bindings::ownership_map imports, exports ;
      std::map<variable, std::string> output_owners ;
      mod library ;
      for(auto &instance : instances) {
        const module_origin &origin = instance.origin ;
        const std::string description = origin.description() ;
        mod::mod_info &definition = library.get_info(origin.library) ;
        typedef const module_interface_v1 *(*describe_function)() ;
        auto describe = reinterpret_cast<describe_function>(
          dlsym(definition.m_library, "loci_module_interface_v1")) ;
        if(describe == nullptr)
          throw module_binding_error(description
                                     + ": missing interface v1") ;
        const module_interface_v1 *interface = describe() ;
        if(interface == nullptr)
          throw module_binding_error(description + ": null interface") ;
        instance.interface = *interface ;
        instance.initialize = reinterpret_cast<init_function>(
          dlsym(definition.m_library, "loci_init_module_v1")) ;
        const variableSet inputs = declarations(interface->inputs,
                                                description) ;
        const variableSet outputs = declarations(interface->outputs,
                                                 description) ;
        if((inputs & outputs) != EMPTY)
          throw module_binding_error(description
                                     + ": input/output overlap") ;
        for(const auto &input : instance.inputs) {
          if(!inputs.inSet(input.first))
            throw module_binding_error(description + ": unknown input "
                                       + input.first.str() + " -> " + input.second.str()) ;
          const variable local = qualify_module_variable(input.first,
                                                         instance.scope) ;
          table.insert(std::make_pair(local, input.second)) ;
          imports[origin.instance] += input.second ;
        }
        for(variable input : inputs)
          if(instance.inputs.find(input) == instance.inputs.end())
            throw module_binding_error(description + ": missing binding "
                                       + input.str()) ;
        for(const auto &output : instance.outputs) {
          if(!outputs.inSet(output.first))
            throw module_binding_error(description + ": unknown output "
                                       + output.first.str() + " -> " + output.second.str()) ;
          published.emplace(qualify_module_variable(output.first, instance.scope),
                            output.second) ;
        }
        for(variable output : outputs) {
          const variable local = qualify_module_variable(output, instance.scope) ;
          const auto bound = instance.outputs.find(output) ;
          const variable destination = bound == instance.outputs.end() ? local : bound->second ;
          const std::string owner = description + " output " + output.str() ;
          const auto inserted = output_owners.emplace(destination, owner) ;
          if(!inserted.second)
            throw module_binding_error("module composition: duplicate output destination "
              + destination.str() + " from " + inserted.first->second + " and " + owner) ;
          exports[origin.instance] += destination ;
        }

        for(auto entry = definition.loaded_rule_list.begin();
            entry != definition.loaded_rule_list.end(); ++entry) {
          if(entry.get_p()->rr->is_module_rule())
            throw module_binding_error(description
              + ": legacy module declarations require explicit host adapters") ;
          rule_implP impl = *entry ;
          if(!impl->get_info().dynamic_constraints.empty()) {
            std::ostringstream message ;
            message << description << ": dynamic constraints are not supported: "
                    << rule(impl) << " at " << impl->get_fileloc() ;
            throw module_binding_error(message.str()) ;
          }
          std::map<variable, variable> names ;
          variableSet variables = impl->get_var_list() ;
          if(impl->is_parametric_provided())
            variables += impl->get_parametric_variable() ;
          for(variable name : variables)
            names[name] = qualify_module_variable(name, instance.scope) ;
          impl->rename_vars(names) ;
          rule r(impl) ;
          if(impl->get_rule_class() == rule_impl::DEFAULT
             || impl->get_rule_class() == rule_impl::OPTIONAL) {
            instance.declarations.add_rule(r) ;
            instance.declarations.add_module_origin(r, origin) ;
          }
          combined.add_rule(r) ;
          combined.add_module_origin(r, origin) ;
        }
      }
      combined.set_module_bindings(module_bindings(table, imports, published, exports)) ;
      validate_module_bindings(combined, combined.get_module_bindings()) ;
      prepared_rules = combined ;
      collect_binding_types() ;
      state = PREPARED ;
    } catch(const module_binding_error &error) {
      state = FAILED ;
      throw module_binding_error("module composition prepare: " + std::string(error.what())) ;
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
      for(const auto &input : bindings.input_entries())
        if(facts.get_extensional_facts().inSet(input.first))
          throw module_binding_error("fact already initializes imported input "
                                     + input.first.str()) ;
      for(const auto &output : bindings.output_entries())
        if(output.first != output.second
           && facts.get_extensional_facts().inSet(output.first))
          throw module_binding_error("fact already initializes relocated output "
                                     + output.first.str()) ;
      check_fact_types(facts) ;
      const rule_db declarations = resolve_module_declarations(prepared_rules, bindings) ;
      facts.setupDefaults(declarations) ;
      facts.setupOptionalTypes(declarations) ;
      for(auto &instance : instances) {
        if(instance.vars_file.empty())
          continue ;
        std::ifstream deck(instance.vars_file.c_str()) ;
        if(!deck)
          throw module_binding_error(instance.origin.description() + ": cannot open "
                                     + instance.vars_file) ;
        module_context context(facts, bindings,
                               instance.scope, instance.origin) ;
        try {
          const rule_db declarations = resolve_module_declarations(instance.declarations, bindings) ;
          facts.read_vars(deck, declarations, [&](variable name) {
            return context.write_name(name) ;
          }) ;
        } catch(const module_binding_error &error) {
          throw module_binding_error(instance.vars_file + ": " + error.what()) ;
        } catch(const BasicException &error) {
          std::ostringstream message ;
          message << instance.origin.description() << ": " << instance.vars_file << ": " ;
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
      check_fact_types(facts) ;
      initialized_facts = &facts ;
      state = READY ;
    } catch(const module_binding_error &error) {
      state = FAILED ;
      throw module_binding_error("module composition initialize: " + std::string(error.what())) ;
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

  void module_composition::print_connections(std::ostream &out) const {
    if(state != PREPARED && state != READY)
      throw module_binding_error("module composition: not prepared") ;
    const auto &bindings = prepared_rules.get_module_bindings() ;
    for(const auto &instance : instances) {
      out << instance.origin.description() << '\n' ;
      for(const auto &name : instance.interface.inputs) {
        const variable local = qualify_module_variable(variable(name), instance.scope) ;
        out << "  input " << local << " -> " << bindings.resolve(local) << '\n' ;
        const variable requested = instance.inputs.at(variable(name)) ;
        if(requested != bindings.resolve(local))
          out << "    requested: " << requested << '\n' ;
      }
      for(const auto &name : instance.interface.outputs) {
        const variable local = qualify_module_variable(variable(name), instance.scope) ;
        out << "  output " << local << " -> " << bindings.resolve(local) << '\n' ;
      }
    }
  }
}
