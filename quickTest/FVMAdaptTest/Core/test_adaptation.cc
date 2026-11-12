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

#include <Loci.h>
#include <FVMAdapt/parametric/affine_subdivision.h>
#include <FVMAdapt/parametric/refinement_record.h>
#include <FVMAdapt/parametric/face_plan.h>
#include <FVMAdapt/parametric/corner_boundary.h>

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

using namespace Loci::FVMAdapt;

/// Directional hex cuts give two children per selected reference direction.
TEST_CASE("hex directional splits survive history reload") {
  for(unsigned mask=0; mask<8; ++mask) {
    CAPTURE(mask);
    AffineSubdivision cell(Shape::hex);
    TensorSplit split; unsigned count=1;
    for(unsigned direction=0; direction<3; ++direction) if(mask & (1u<<direction)) {
      split.cuts[direction]={Fraction(1,2)}; count*=2;
    }
    cell.splitLeaf(0,split);
    CHECK(cell.leafIds().size()==count);
    AffineSubdivision restored(decodeRefinement(encodeRefinement(cell.history())));
    CHECK(restored.leafIds().size()==count);
  }
}

/// Native prism cuts retain triangular ends: axial, transverse, and combined.
TEST_CASE("prism directional splits survive history reload") {
  for(unsigned mask=0; mask<4; ++mask) {
    CAPTURE(mask);
    AffineSubdivision cell(Shape::prism); PrismSplit split; unsigned count=1;
    if(mask & 1) { split.axialCuts={Fraction(1,2)}; count*=2; }
    if(mask & 2) { split.base=TriangleSplit::uniform4(); count*=4; }
    cell.splitLeaf(0,split);
    CHECK(cell.leafIds().size()==count);
    AffineSubdivision restored(decodeRefinement(encodeRefinement(cell.history())));
    CHECK(restored.leafIds().size()==count);
  }
}

/// Nested refinement can be restored and coarsened one family at a time.
TEST_CASE("nested cell refinement preserves complete families") {
  for(auto shape:{Shape::hex,Shape::prism}) {
    CAPTURE(int(shape));
    AffineSubdivision cell(shape);
    cell.refine(0,isotropicSplit(shape));
    const auto first=cell.children(0)[0];
    cell.refine(first,isotropicSplit(shape));
    const auto saved=encodeRefinement(cell.history());
    AffineSubdivision restored(decodeRefinement(saved));
    CHECK(restored.leafIds().size()==15);
    restored.coarsen(restored.children(0)[0]);
    CHECK(restored.leafIds().size()==8);
    restored.coarsen(0);
    CHECK(restored.leafIds().size()==1);
  }
}

namespace {
  CornerBoundary square() {
    const std::vector<Point> points={{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}}};
    CornerBoundary face; std::vector<BoundaryId> vertices,edges;
    for(const auto &point:points) vertices.push_back(face.addPoint(point));
    for(unsigned i=0;i<4;++i) edges.push_back(face.addEdge(i,(i+1)%4,
      std::make_shared<CornerMap>(Shape::segment,std::vector<Point>{points[i],points[(i+1)%4]})));
    face.addSurface(vertices,edges,std::make_shared<CornerMap>(Shape::quad,points));
    return face;
  }
  void divide(CornerBoundary &face,BoundaryId id) {
    CornerBoundary::Transaction edit(face);
    face.splitSurface(id,CenterMethod::reference,{},edit); edit.committed=true;
  }
}

/// Both neighbors retain their requested face divisions.
TEST_CASE("shared-face plans merge both neighbors refinement") {
  auto left=square(),right=square();
  divide(left,0); divide(right,0);
  divide(left,left.boundaries.faces[0].splits[0].children[0]);
  divide(right,right.boundaries.faces[0].splits[0].children[2]);
  const auto a=encodeFacePlan(left.boundaries,0),b=encodeFacePlan(right.boundaries,0);
  const auto merged=mergeFacePlans(a,b);
  CHECK(merged==mergeFacePlans(b,a));
  CHECK(applyFacePlan(left,0,merged)); CHECK(applyFacePlan(right,0,merged));
  CHECK(encodeFacePlan(left.boundaries,0)==encodeFacePlan(right.boundaries,0));
  for(unsigned child:{0u,2u}) {
    const auto id=left.boundaries.faces[0].splits[0].children[child];
    REQUIRE(left.boundaries.faces[id].splits.size()==1);
    CHECK(left.boundaries.faces[id].splits[0].children.size()==4);
  }
}

int main(int argc, char** argv) {
  Loci::Init(&argc, &argv);
  doctest::Context context; context.applyCommandLine(argc, argv);
  const int result=context.run();
  Loci::Finalize();
  return result;
}
