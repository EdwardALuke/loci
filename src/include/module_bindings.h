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

  /// Concrete input identities resolved after parametric rule expansion.
  /// This is a low-level composition primitive, not a module loader API.
  class module_bindings {
  public:
    typedef std::map<variable, variable> binding_map ;

    module_bindings() {} ;

    /// Keys must be qualified concrete names; destinations may be at root.
    /// Both must be stationary, unprioritized names. Chains are rejected.
    explicit module_bindings(const binding_map &inputs) ;

    bool empty() const { return inputs.empty() ; }
    const binding_map &entries() const { return inputs ; }

    /// Resolve a whole fact occurrence, retaining its time and priority.
    /// Parameter argument labels are unchanged, including qualified labels.
    variable resolve(variable value) const ;
    variableSet resolve(const variableSet &values) const ;

    /// True for an imported local identity, including its time variants.
    bool is_input(variable value) const ;

  private:
    binding_map inputs ;
  } ;

  /// Reject writes to imported local identities before setup or execution.
  /// Call again after expansion to check substitution-generated targets.
  void validate_module_bindings(const rule_db &rules,
                               const module_bindings &bindings) ;

  /// Resolve expanded rules without modifying their original implementations.
  /// Existing rule renaming handles maps, constraints and store references.
  rule_db resolve_module_bindings(const rule_db &rules,
                                 const module_bindings &bindings) ;
}

#endif
