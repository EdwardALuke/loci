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

#ifndef LOCI_FVMADAPT_SUBDIVISION_INTERNAL_H
#define LOCI_FVMADAPT_SUBDIVISION_INTERNAL_H

#include <FVMAdapt/parametric/subdivision.h>

namespace Loci {
namespace FVMAdapt {
namespace detail {

  /// Construct validated native children without face or edge correspondence.
  std::vector<Child> splitChildren(Shape shape, const TensorSplit &request) ;
  std::vector<Child> splitChildren(const TriangleSplit &request) ;
  std::vector<Child> splitChildren(const PrismSplit &request) ;
  std::vector<Child> splitChildren(const TetSplit &request) ;
  std::vector<Child> splitChildren(const PyramidSplit &request) ;

  /// Place a positively oriented tet, exchanging the last two vertices when necessary.
  Child tetChild(std::array<ReferencePoint,4> vertices) ;

  /// Derive complete boundary and edge correspondence for a conforming split template.
  void completeSubdivision(SplitResult &result) ;

  /// Validate accepted child regions against the requested native partition in linear work.
  /// Positive chart orientations may differ; accepted placements are never regenerated.
  void validateSubdivision(Shape parent, const std::vector<Child> &children,
                           const SplitOperation &operation, ExactBudget &budget) ;

}
}
}

#endif
