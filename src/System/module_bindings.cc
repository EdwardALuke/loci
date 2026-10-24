//#############################################################################
//# Copyright 2026, Streamline Numerics, Inc.
//#
//# This file is part of the Loci Framework, distributed under the GNU Lesser
//# General Public License, version 3 or (at your option) any later version.
//# See <http://www.gnu.org/licenses/>.
//#############################################################################
#include <module_bindings.h>
#include <rule.h>
#include <sstream>

namespace Loci {
  namespace {
    variable base_name(variable value) {
      variable::info base ;
      base.name = value.get_info().name ;
      base.namespac = value.get_info().namespac ;
      return variable(base) ;
    }

    bool concrete_name(variable value) {
      const variable::info &info = value.get_info() ;
      return !info.name.empty() && !info.tvar && !info.assign
        && info.offset == 0 && info.time_id == time_ident()
        && info.priority.empty() && info.v_ids.empty()
        && info.name != "UNIVERSE" && info.name != "EMPTY"
        && info.name != "OUTPUT" && info.name != "*NONAME*" ;
    }

    void invalid_binding(const char *message, variable source, variable destination) {
      std::ostringstream out ;
      out << "module bindings: " << message << ": " << source << " -> " << destination ;
      throw module_binding_error(out.str()) ;
    }

    ruleSet all_rules(const rule_db &rules) {
      ruleSet result = rules.all_rules() ;
      result += rules.get_default_rules() ;
      result += rules.get_optional_rules() ;
      return result ;
    }

    std::string rule_location(rule provider, const rule_db &rules) {
      std::ostringstream out ;
      for(const auto &origin : rules.get_module_origins(provider))
        out << origin.description() << ": " ;
      out << provider ;
      const std::string file = provider.get_rule_implP()->get_fileloc() ;
      if(!file.empty()) out << " at " << file ;
      return out.str() ;
    }

    void check_write(variable value, rule provider, const rule_db &rules,
                     const module_bindings &bindings) {
      if(bindings.is_input(value)) {
        std::ostringstream out ;
        out << "module bindings: rule writes imported input " << value
            << " bound to " << bindings.resolve(value)
            << ": " << rule_location(provider, rules) ;
        throw module_binding_error(out.str()) ;
      }
      for(const auto &origin : rules.get_module_origins(provider))
        if(bindings.is_import(bindings.resolve(value), origin.instance)) {
          std::ostringstream out ;
          out << "module bindings: rule writes imported destination " << value
              << " resolved as " << bindings.resolve(value)
              << ": " << rule_location(provider, rules) ;
          throw module_binding_error(out.str()) ;
        }
    }
  }

  module_bindings::module_bindings(const binding_map &values,
                                   const ownership_map &owners,
                                   const binding_map &published,
                                   const ownership_map &writers)
    : inputs(values), outputs(published), imports(owners), exports(writers) {
    for(const auto &entry : outputs) {
      if(!concrete_name(entry.first)
         || entry.first.get_info().namespac.empty())
        invalid_binding("expected qualified concrete output", entry.first, entry.second) ;
      if(!concrete_name(entry.second))
        invalid_binding("expected concrete destination", entry.first, entry.second) ;
      if(inputs.count(entry.first) || inputs.count(entry.second))
        invalid_binding("output aliases an imported input", entry.first, entry.second) ;
      if(entry.first != entry.second && outputs.count(entry.second))
        invalid_binding("output destination is another output alias", entry.first, entry.second) ;
    }
    for(auto &entry : inputs) {
      if(!concrete_name(entry.first)
         || entry.first.get_info().namespac.empty())
        invalid_binding("expected qualified concrete input", entry.first, entry.second) ;
      if(!concrete_name(entry.second))
        invalid_binding("expected concrete destination", entry.first, entry.second) ;
      if(inputs.find(entry.second) != inputs.end())
        invalid_binding("destination is another input", entry.first, entry.second) ;
      const auto output = outputs.find(entry.second) ;
      if(output != outputs.end())
        entry.second = output->second ;
    }
    // Ownership uses final identities, just like rules and setup access.
    for(auto &owner : imports)
      owner.second = resolve(owner.second) ;
    for(auto &owner : exports)
      owner.second = resolve(owner.second) ;
  }

  bool module_bindings::is_import(variable value,
                                 const std::string &origin) const {
    if(value.is_time_variable() || !value.get_arg_list().empty())
      return false ;
    const auto owner = imports.find(origin) ;
    return owner != imports.end() && owner->second.inSet(base_name(value)) ;
  }

  bool module_bindings::is_input(variable value) const {
    if(value.is_time_variable() || !value.get_arg_list().empty())
      return false ;
    return inputs.find(base_name(value)) != inputs.end() ;
  }

  bool module_bindings::is_export(variable value,
                                 const std::string &origin) const {
    if(value.is_time_variable() || !value.get_arg_list().empty())
      return false ;
    const auto owner = exports.find(origin) ;
    return owner != exports.end() && owner->second.inSet(base_name(value)) ;
  }

  variable module_bindings::resolve(variable value) const {
    if(empty() || value.is_time_variable() || !value.get_arg_list().empty())
      return value ;
    const auto input = inputs.find(base_name(value)) ;
    const auto output = outputs.find(base_name(value)) ;
    if(input == inputs.end() && output == outputs.end())
      return value ;
    const variable destination = input != inputs.end() ? input->second : output->second ;
    variable::info info = value.get_info() ;
    info.name = destination.get_info().name ;
    info.namespac = destination.get_info().namespac ;
    return variable(info) ;
  }

  variableSet module_bindings::resolve(const variableSet &values) const {
    variableSet result ;
    for(variable value : values)
      result += resolve(value) ;
    return result ;
  }

  void validate_module_bindings(const rule_db &rules,
                               const module_bindings &bindings) {
    if(bindings.empty())
      return ;
    for(rule provider : all_rules(rules)) {
      rule_implP impl = provider.get_rule_implP() ;
      if(impl == 0)
        continue ;
      const auto &info = impl->get_info() ;
      // Existing rule renaming does not cover dynamic constraints.
      if(!info.dynamic_constraints.empty())
        throw module_binding_error(
          "module bindings: dynamic constraints are not supported: "
          + rule_location(provider, rules)) ;
      if(impl->get_rule_class() == rule_impl::OPTIONAL)
        continue ;
      for(variable output : info.output_vars())
        check_write(output, provider, rules, bindings) ;
      for(const auto &target : info.targets)
        for(const auto &alias : target.assign) {
          check_write(alias.first, provider, rules, bindings) ;
          check_write(alias.second, provider, rules, bindings) ;
        }
    }
  }

  rule_db resolve_module_bindings(const rule_db &rules,
                                 const module_bindings &bindings) {
    if(bindings.empty())
      return rules ;
    validate_module_bindings(rules, bindings) ;
    rule_db result ;
    for(rule provider : all_rules(rules)) {
      rule_implP impl = provider.get_rule_implP() ;
      if(impl == 0) {
        result.add_rule(provider) ;
        continue ;
      }
      std::map<variable, variable> names ;
      for(variable value : impl->get_var_list()) {
        const variable resolved = bindings.resolve(value) ;
        if(resolved != value)
          names[value] = resolved ;
      }
      const rule resolved = names.empty() ? provider
        : provider.rename_vars(names) ;
      result.add_rule(resolved) ;
      result.copy_module_origins(resolved, rules, provider) ;
    }
    return result ;
  }

  rule_db resolve_module_declarations(const rule_db &rules,
                                     const module_bindings &bindings) {
    rule_db declarations ;
    ruleSet selected = rules.get_default_rules() ;
    selected += rules.get_optional_rules() ;
    for(rule r : selected) {
      declarations.add_rule(r) ;
      declarations.copy_module_origins(r, rules, r) ;
    }
    return resolve_module_bindings(declarations, bindings) ;
  }
}
