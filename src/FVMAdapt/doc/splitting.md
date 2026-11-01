# Splitting Cells {#fvmadapt_geometric_splitting}

Refinement divides a cell in local coordinates. Its geometry map places the
new points in physical space. Both neighbors use the same shared-face divisions
and accepted centers.

## Available Splits

| Construction | Result | Center choice |
| --- | --- | --- |
| Hex cuts (`TensorSplit`) | Cuts in selected local directions; three midpoint cuts give eight hexes | Cut positions define the split |
| Prism cuts (`PrismSplit`) | Triangle subdivision and independent axial cuts | Cut positions define the split |
| Tet cuts (`TetSplit`) | Edge bisection or eight tets; the latter has three central diagonal choices | Fixed pattern |
| Pyramid cuts (`PyramidSplit`) | Six pyramids and four tets | Fixed pattern; no directional split |
| Split around each vertex, hex/prism/tet/pyramid | Eight/six/four/five children; pyramid children need not be hexes | Reference, weighted or integrated |
| Planar triangle/quad | Three/four quads, extended through a fixed thickness for 2D meshes | Reference, weighted or integrated |
| General source cell with polygon faces | One child per original vertex | Reference, weighted or integrated |

`TensorSplit` also handles segments and quads. `TriangleSplit` supplies edge
bisection and four-triangle splitting, including the prism's base. Cuts may be
non-midpoint. A child's shape can differ from its parent's.

Choose `splitFamily="mixed"` to use native cuts in hexes and prisms and
one-child-per-vertex splits in tets, pyramids and general cells. Set
`centerMethod="reference"`; `automaticSplit=true` selects long hex/prism
directions. Mixed mode is 3D and supports midpoint cuts at mixed interfaces.
Weighted/integrated mixed interfaces and 2D extrusion are not supported.

Mixed neighbors intersect their face pieces in the original face's local
coordinates. Both sides use the resulting points and the same saved face map.
This matching is temporary module work; it adds no persistent store or MPI
operation to the kernel.

The one-child-per-vertex construction connects edge midpoints, face centers
and a cell center. Standard hex, prism and tet parents produce hex children.
A pyramid produces four hexes and a ten-vertex child at the apex. General
cells use the same construction on their bounding faces.

<img src="figures/faces/quadrilaterial/quadrilateral_face_isotropic_split.svg"
     alt="A quad face split from its edge midpoints to its center" width="520">

<img src="figures/volume_cells/hexahedron/isotropic/hexahedron_midpoint_based_isotropic.svg"
     alt="A hex split using edge midpoints, face centers and a cell center" width="420">

The cell chooses a split; its faces and edges supply the matching divisions.
The library records the accepted points. Repeating a split after coarsening
reuses those points, even if the preferred center method has changed.

For a standard directional cut:

```cpp
AffineSubdivision cells(Shape::hex);
TensorSplit cut;
cut.cuts[0] = {Fraction(1,2)};
cells.splitLeaf(0, cut);
```

For straight 2D meshes, the planar triangle or quad splits first. The same
points are used on both end faces, preserving the extrusion thickness.
Local directions need not align with physical axes.

A valid reference split does not guarantee good physical cells. Check mapped
geometry and output quality before accepting new constructions. See
[cell geometry](@ref fvmadapt_reference_geometry).
