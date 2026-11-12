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

#include <FVMAdapt/parametric/mapped_volume_subdivision.h>
#include <doctest.h>

using namespace Loci::FVMAdapt;

namespace {
  Point center(Shape shape) {
    const auto &element=referenceElement(shape); Point p{{0,0,0}};
    for(unsigned v=0;v<element.vertexCount;++v)
      for(unsigned d=0;d<3;++d) p[d]+=element.vertices[v][d]/element.vertexCount;
    return p;
  }
}

/// The vertex-centered construction gives eight hex or six prism children;
/// splitting every child again gives 64 or 48 leaves.
TEST_CASE("two isotropic levels reach every hex and prism child") {
  for(auto shape:{Shape::hex,Shape::prism}) {
    CAPTURE(int(shape));
    const auto &element=referenceElement(shape);
    auto map=std::make_shared<CornerMap>(shape,
      std::vector<Point>(element.vertices.begin(),element.vertices.begin()+element.vertexCount));
    MappedVolumeSubdivision cell(map);
    for(unsigned level=0;level<2;++level) {
      for(auto leaf:cell.leafIds()) {
        const auto &child=referenceElement(cell.geometry(leaf)->shape());
        std::vector<Point> centers;
        for(unsigned f=0;f<child.boundaryCount;++f) centers.push_back(center(child.boundaries[f].shape));
        cell.splitLeaf(leaf,centers,center(child.shape));
      }
      CHECK(cell.leafIds().size()==(level ? 8u:1u)*element.vertexCount);
    }
  }
}
