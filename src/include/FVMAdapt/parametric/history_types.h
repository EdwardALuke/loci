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

#ifndef LOCI_FVMADAPT_HISTORY_TYPES_H
#define LOCI_FVMADAPT_HISTORY_TYPES_H

#include "subdivision.h"
#include "vertex_volume_split.h"
#include "planar_split.h"

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <optional>
#include <vector>

namespace Loci {
namespace FVMAdapt {

  /// Index in one history snapshot; the root is zero. Compaction invalidates indices.
  typedef std::size_t TreeNodeId ;
  constexpr TreeNodeId noTreeNode = static_cast<TreeNodeId>(-1) ;

  /// Child ordinals from the root; a path alone does not identify a replaced construction.
  typedef std::vector<unsigned> CellPath ;

  /// Contiguous child indices in one history snapshot.
  class NodeRange {
    TreeNodeId first_ ;
    std::size_t size_ ;
  public:
    class iterator {
      TreeNodeId id_ ;
    public:
      typedef std::input_iterator_tag iterator_category ;
      typedef TreeNodeId value_type ;
      typedef std::ptrdiff_t difference_type ;
      typedef void pointer ;
      typedef TreeNodeId reference ;
      explicit iterator(TreeNodeId id = 0) : id_(id) {}
      TreeNodeId operator*() const { return id_ ; }
      iterator &operator++() { ++id_ ; return *this ; }
      iterator operator++(int) { iterator previous = *this ; ++*this ; return previous ; }
      bool operator==(iterator other) const { return id_ == other.id_ ; }
      bool operator!=(iterator other) const { return id_ != other.id_ ; }
    } ;
    NodeRange(TreeNodeId first = 0, std::size_t size = 0) : first_(first), size_(size) {}
    std::size_t size() const { return size_ ; }
    bool empty() const { return size_ == 0 ; }
    /// Indexing requires index < size(); front() requires a nonempty range.
    TreeNodeId operator[](std::size_t index) const { return first_+index ; }
    TreeNodeId front() const { return first_ ; }
    iterator begin() const { return iterator(first_) ; }
    iterator end() const { return iterator(first_+size_) ; }
  } ;

  /// Explicit construction tags; values are stable in saved records.
  enum class SplitFamily : unsigned { tensor = 1, triangle = 2, prism = 3, tet = 4, pyramid = 5, vertexVolume = 6, planar = 7 } ;

  /// Accepted local placements, shared by child identities independently of table indices.
  struct SplitConstruction {
    /// Reference shape of the parent.
    Shape parentShape ;
    SplitFamily family ;
    unsigned version = 1 ;
    /// Affine placements for native splits; empty for mapped constructions.
    std::vector<Child> children ;
    /// Accepted nonlinear reference charts for a vertex-volume split.
    std::shared_ptr<const VertexVolumeSplit> vertexVolume ;
    /// Accepted physical polygons for a planar split.
    std::shared_ptr<const PlanarSplit> planar ;
    std::size_t childCount() const ;
  } ;

  bool operator==(const SplitConstruction &a, const SplitConstruction &b) ;
  bool operator<(const SplitConstruction &a, const SplitConstruction &b) ;

  /// Nonnull accepted construction and the child's ordinal within it.
  struct CellStep {
    std::shared_ptr<const SplitConstruction> construction ;
    unsigned child ;
  } ;

  /// Persistent root and accepted construction lineage; geometry revisions are separate.
  struct CellKey {
    std::uint64_t root ;
    std::vector<CellStep> lineage ;
  } ;

  bool operator==(const CellKey &a, const CellKey &b) ;
  bool operator!=(const CellKey &a, const CellKey &b) ;
  bool operator<(const CellKey &a, const CellKey &b) ;

}
}

#endif
