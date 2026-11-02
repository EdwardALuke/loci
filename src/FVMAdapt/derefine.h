// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.
#ifndef FVMADAPT_DEREFINE_H
#define FVMADAPT_DEREFINE_H

class Cell ;
class HexCell ;
class Prism ;

namespace Loci {
  enum class DerefineTags { Nodes, Cells } ;

  /// Collapse eligible parents selected from the current leaves.
  void derefine_marked_cells(Cell& cell, DerefineTags tags) ;
  void derefine_marked_cells(HexCell& cell, DerefineTags tags) ;
  void derefine_marked_cells(Prism& cell, DerefineTags tags) ;
}
#endif
