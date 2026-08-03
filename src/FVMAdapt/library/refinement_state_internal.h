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
#ifndef FVMADAPT_REFINEMENT_STATE_INTERNAL_H
#define FVMADAPT_REFINEMENT_STATE_INTERNAL_H

#include <FVMAdapt/defines.h>
#include <FVMAdapt/refinement_state.h>

#include <utility>
#include <vector>

class Cell ;
class HexCell ;
class Prism ;

namespace Loci {
  namespace detail {

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
    bool collectRefinedCellState(
      refinedCellState& state,
      const const_store<std::vector<int> >& fineDepth,
      const const_store<std::vector<int> >* fineResult,
      const const_store<int>& cellOffset,
      const const_store<int>& rootFileNumber,
      const entitySet& sourceCells,
      const std::vector<entitySet>& localCells,
      Entity cellBase) ;

  }
}

#endif
