// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.

#include <list>
#include <set>
#include "../derefine.h"
#include "hexcell.h"
#include "prism.h"

namespace Loci {
  namespace {
    template<class CellType>
    void derefine_leaves(CellType& cell, DerefineTags tags) {
      std::list<CellType*> leaves ;
      cell.sort_leaves(leaves) ;
      std::set<CellType*> parents ;
      for(CellType* leaf : leaves) {
        const char tag = tags == DerefineTags::Nodes ? leaf->get_tagged() : leaf->getTag() ;
        if(tag != 2) continue ;
        CellType* parent = leaf->getParentCell() ;
        if(parent != 0 && (tags == DerefineTags::Nodes ?
                          parent->needDerefine() : parent->needDerefine_ctag()))
          parents.insert(parent) ;
      }
      // Select all parents before deleting any children.
      for(CellType* parent : parents) parent->derefine() ;
    }
  }

  void derefine_marked_cells(HexCell& cell, DerefineTags tags) {
    derefine_leaves(cell,tags) ;
  }

  void derefine_marked_cells(Prism& cell, DerefineTags tags) {
    derefine_leaves(cell,tags) ;
  }
}
