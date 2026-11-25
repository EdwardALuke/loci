// Copyright 2026, Mississippi State University. LGPL-3.0-or-later.

#include <list>
#include <set>
#include "../derefine.h"
#include "diamondcell.h"
#include "hexcell.h"
#include "prism.h"

namespace Loci {
  namespace {
    template<class Leaf, class Root>
    bool derefine_leaves(Root& cell, DerefineTags tags) {
      std::list<Leaf*> leaves ;
      cell.sort_leaves(leaves) ;
      std::set<Leaf*> parents ;
      bool check_root = false ;
      for(Leaf* leaf : leaves) {
        const char tag = tags == DerefineTags::Nodes ? leaf->get_tagged() : leaf->getTag() ;
        if(tag != 2) continue ;
        Leaf* parent = leaf->getParentCell() ;
        if(parent == 0) check_root = true ;
        if(parent != 0 && (tags == DerefineTags::Nodes ?
                          parent->needDerefine() : parent->needDerefine_ctag()))
          parents.insert(parent) ;
      }
      // Select all parents before deleting any children.
      for(Leaf* parent : parents) parent->derefine() ;
      return check_root ;
    }
  }

  void derefine_marked_cells(Cell& cell, DerefineTags tags) {
    // Direct DiamondCell children have no DiamondCell parent. Check the Cell last.
    if(derefine_leaves<DiamondCell>(cell,tags) &&
       (tags == DerefineTags::Nodes ? cell.needDerefine() : cell.needDerefine_ctag()))
      cell.derefine() ;
  }

  void derefine_marked_cells(HexCell& cell, DerefineTags tags) {
    derefine_leaves<HexCell>(cell,tags) ;
  }

  void derefine_marked_cells(Prism& cell, DerefineTags tags) {
    derefine_leaves<Prism>(cell,tags) ;
  }
}
