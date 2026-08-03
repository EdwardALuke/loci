//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# This program is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//#############################################################################
#ifndef FVMADAPT_REFINEMENT_STATE_H
#define FVMADAPT_REFINEMENT_STATE_H

#include <store.h>

namespace Loci {

  /// Values stored in the public adaptResult fact.
  namespace adapt_result {
    enum value {
      derefined = -1,
      retained = 0,
      refined = 1
    } ;
  }

  /// Per-cell state carried with a newly generated grid.
  struct refinedCellState {
    store<int> refinementDepth ;
    store<int> rootCellFileNumber ;
    store<int> adaptResult ;
    bool hasAdaptResult ;

    refinedCellState() : hasAdaptResult(false) {}
  } ;
}

#endif
