//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#include <FVMAdapt/parametric/history_types.h>

#include <algorithm>
#include <tuple>

namespace Loci {
namespace FVMAdapt {
  namespace {
    bool childLess(const Child &a, const Child &b) {
      return std::tie(a.shape,a.placement.origin,a.placement.columns) <
             std::tie(b.shape,b.placement.origin,b.placement.columns) ;
    }
    bool stepLess(const CellStep &a, const CellStep &b) {
      if(*a.construction < *b.construction) return true ;
      if(*b.construction < *a.construction) return false ;
      return a.child < b.child ;
    }
  }

  std::size_t SplitConstruction::childCount() const {
    if(planar) return planar->children.size() ;
    return vertexVolume ? vertexVolume->children.size() : children.size() ;
  }

  bool operator==(const SplitConstruction &a, const SplitConstruction &b) {
    return !(a < b) && !(b < a) ;
  }

  bool operator<(const SplitConstruction &a, const SplitConstruction &b) {
    const auto first = std::tie(a.parentShape,a.family,a.version) ;
    const auto second = std::tie(b.parentShape,b.family,b.version) ;
    if(first != second) return first < second ;
    if(bool(a.vertexVolume) != bool(b.vertexVolume)) return bool(a.vertexVolume) < bool(b.vertexVolume) ;
    if(a.vertexVolume) return a.vertexVolume->points < b.vertexVolume->points ;
    if(bool(a.planar) != bool(b.planar)) return bool(a.planar) < bool(b.planar) ;
    if(a.planar) return a.planar->points < b.planar->points ;
    return std::lexicographical_compare(a.children.begin(),a.children.end(),
                                        b.children.begin(),b.children.end(),childLess) ;
  }

  bool operator==(const CellKey &a, const CellKey &b) {
    return a.root == b.root && a.lineage.size() == b.lineage.size() &&
      std::equal(a.lineage.begin(),a.lineage.end(),b.lineage.begin(),
        [](const CellStep &first, const CellStep &second) {
          return first.child == second.child && *first.construction == *second.construction ;
        }) ;
  }

  bool operator!=(const CellKey &a, const CellKey &b) { return !(a == b) ; }

  bool operator<(const CellKey &a, const CellKey &b) {
    if(a.root != b.root) return a.root < b.root ;
    return std::lexicographical_compare(a.lineage.begin(),a.lineage.end(),
                                        b.lineage.begin(),b.lineage.end(),stepLess) ;
  }

}
}
