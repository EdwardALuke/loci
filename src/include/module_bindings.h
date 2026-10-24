//#############################################################################
//# Copyright 2026, Streamline Numerics, Inc.
//#
//# This file is part of the Loci Framework, distributed under the GNU Lesser
//# General Public License, version 3 or (at your option) any later version.
//# See <http://www.gnu.org/licenses/>.
//#############################################################################
#ifndef LOCI_MODULE_BINDINGS_H
#define LOCI_MODULE_BINDINGS_H

#include <Tools/variable.h>
#include <map>
#include <stdexcept>

namespace Loci {
  class rule_db ;

  /// Invalid composition; no execution schedule may use these bindings.
  class module_binding_error : public std::invalid_argument {
  public:
    explicit module_binding_error(const std::string &message)
      : std::invalid_argument(message) {} ;
  } ;

  /// Namespace identity is separate from the library label used in diagnostics.
  struct module_origin {
    std::string instance, library ;
    std::string description() const { return library + "[" + instance + "]" ; }
    bool operator<(const module_origin &other) const {
      return instance < other.instance
        || (instance == other.instance && library < other.library) ;
    }
  } ;

  /// Concrete interface identities resolved after parametric rule expansion.
  /// This is a low-level composition primitive, not a module loader API.
  class module_bindings {
  public:
    typedef std::map<variable, variable> binding_map ;
    typedef std::map<std::string, variableSet> ownership_map ;

    module_bindings() {} ;

    /// Keys must be qualified concrete names; destinations may be at root.
    /// Both must be stationary, unprioritized names. Inputs may name an output
    /// alias; preparation resolves that single hop. Other chains are rejected.
    /// Imports and exports record destinations accessible by each rule origin.
    explicit module_bindings(const binding_map &inputs,
                             const ownership_map &imports = ownership_map(),
                             const binding_map &outputs = binding_map(),
                             const ownership_map &exports = ownership_map()) ;

    bool empty() const { return inputs.empty() && outputs.empty() ; }
    const binding_map &input_entries() const { return inputs ; }
    const binding_map &output_entries() const { return outputs ; }
    bool operator==(const module_bindings &other) const {
      return inputs == other.inputs && imports == other.imports
        && outputs == other.outputs && exports == other.exports ;
    }

    /// Resolve a whole fact occurrence, retaining its time and priority.
    /// Parameter argument labels are unchanged, including qualified labels.
    variable resolve(variable value) const ;
    variableSet resolve(const variableSet &values) const ;

    /// True for an imported local identity, including its time variants.
    bool is_input(variable value) const ;

    /// Check a canonical destination against its originating module's imports.
    bool is_import(variable value, const std::string &origin) const ;

    /// True for a declared output destination writable during instance setup.
    bool is_export(variable value, const std::string &origin) const ;

  private:
    binding_map inputs, outputs ;
    ownership_map imports, exports ;
  } ;

  /// Reject writes to imported local identities before setup or execution.
  /// Call again after expansion to check substitution-generated targets.
  void validate_module_bindings(const rule_db &rules,
                               const module_bindings &bindings) ;

  /// Resolve expanded rules without modifying their original implementations.
  /// Existing rule renaming handles maps, constraints and store references.
  rule_db resolve_module_bindings(const rule_db &rules,
                                 const module_bindings &bindings) ;

  /// Resolve defaults and optional types for setup, without expanding families.
  rule_db resolve_module_declarations(const rule_db &rules,
                                     const module_bindings &bindings) ;
}

#endif
