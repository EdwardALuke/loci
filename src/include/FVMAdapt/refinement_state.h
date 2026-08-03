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

#include "defines.h"

#include <utility>
#include <vector>

class Cell ;
class HexCell ;
class Prism ;

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

  /// Gather accepted leaves by their plan-assigned index and report tree depth.
  bool getLeafRefinementDepths(const Cell* root,
                               std::vector<int>& depths) ;
  bool getLeafRefinementDepths(const HexCell* root,
                               std::vector<int>& depths) ;
  bool getLeafRefinementDepths(const Prism* root,
                               std::vector<int>& depths) ;

  /// Classify each new leaf from the cardinality of the old/new relation.
  bool classifyAdaptResult(
    const std::vector<std::pair<int32, int32> >& cell2parent,
    int firstNewCell,
    int numberOfNewCells,
    std::vector<int>& result) ;

  /// Flatten rule-derived state and place it on the generated-cell partition.
  bool collectRefinedCellState(refinedCellState& state,
                               fact_db& facts,
                               const std::vector<entitySet>& localCells,
                               int cellBase,
                               bool includeAdaptResult) ;
}

#endif
