// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_PARAMETRIC_PLAN_H
#define LOCI_FVMADAPT_PARAMETRIC_PLAN_H
#include "gridInterface.h"

namespace Loci {
namespace FVMAdapt {
  /// Balance file-based requests and retain history without constructing the output mesh.
  /// tagFile uses the current mesh's file order; an empty name selects geometric input only.
  bool markParametricPlan(rule_db &rules,RefinementState &state,const std::string &mesh,
      const RefinementOptions &options,const std::string &tagFile,std::string &error) ;
  /// Mark, balance and write a mesh in one grid load. Retain the resulting history in state.
  /// Optional donor pairs refer to the input state, or the original cells when starting fresh.
  bool refineParametricFile(rule_db &rules,RefinementState &state,const std::string &mesh,const RefinementOptions &options,
      const std::string &tagFile,const std::string &output,const std::string &cellParents,std::string &error) ;
  /// Reconstruct a saved plan. Optional donor pairs refer to parent, or the original grid if null.
  bool writeParametricMesh(rule_db &rules,const RefinementState &state,const RefinementState *parent,
      const RefinementOptions &options,const std::string &output,const std::string &cellParents,std::string &error) ;
  /// Collective plan I/O. Files retain geometry controls and distributed histories, not request files.
  /// Parallel rank counts may change; serial/parallel transitions are unsupported.
  bool readParametricPlan(const std::string &file,const std::string &mesh,RefinementState &state,
      RefinementOptions &options,std::string &error) ;
  bool writeParametricPlan(const std::string &file,const RefinementState &state,const RefinementOptions &options,std::string &error) ;
}
}
#endif
