// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef LOCI_FVMADAPT_CELL_STATE_INTERNAL_H
#define LOCI_FVMADAPT_CELL_STATE_INTERNAL_H

namespace Loci {
  class refinedGridData ;
  class fact_db ;
  /// Distribute optional root-owned leaf fields onto the generated cells.
  void collectCellFields(refinedGridData &, fact_db &, bool cellState, bool edgeLengths) ;
}
#endif
