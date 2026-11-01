// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMTOOLS_ADAPT_INPUT_H
#define LOCI_FVMTOOLS_ADAPT_INPUT_H
#include <FVMAdapt/parametric_plan.h>

namespace FVMTools {
  /// Files, controls and optional saved history for marker and refine.
  struct AdaptInput {
    std::string mesh,tagFile,outputMesh,outputPlan,cellParents ;
    Loci::RefinementOptions options ;
    Loci::RefinementState state ;
  } ;
  /// Parse requests and restore controls from -r before applying command-line overrides.
  AdaptInput readAdaptInput(int argc,char **argv,bool meshOutput) ;
  void showAdaptHelp(bool meshOutput) ;
}
#endif
