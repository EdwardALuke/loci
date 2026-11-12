# Cell Geometry {#fvmadapt_reference_geometry}

A geometry map converts a point in a simple reference shape into its physical
position. The reference shape and map form a parametric representation.
For example, local coordinates in a unit cube describe positions inside a
distorted hex. Local directions need not align with physical X, Y and Z.

## Geometry from Corners

VOG supplies node positions and face connectivity. These define the following
geometry:

| Shape | Geometry from its corners |
| --- | --- |
| Triangle, tet | Linear |
| Quad | Bilinear; the corners need not be coplanar |
| Polygon with more than four vertices | Rational surface over a regular reference polygon |
| Hex | Trilinear |
| Triangular prism | Linear over the triangle and along its axis |
| Pyramid | Bilinear base joined to the apex |
| General polyhedron | Shared face maps bounding one cell |

A polygon surface interpolates its corners and straight edges; noncoplanar
corners can produce a curved surface. Some faces use triangular patches
sharing one center. Both neighboring cells use the same face geometry.
Corners alone do not define an intended CAD surface. Cubic edge controls can
describe curved boundaries in addition to the original mesh.

## Straight and Cubic Edges

An edge's local coordinate `t` runs from zero at its first endpoint to one at
its second. A straight edge interpolates between these endpoints. A cubic
Bezier edge uses four ordered points: endpoints `P0`, `P3` and controls `P1`,
`P2`. The controls steer the curve and generally lie off it; they are not
additional mesh vertices.

<img src="figures/cubic_edge.svg"
     alt="Cubic edge from P0 to P3, with controls P1 and P2 off the curve and a split vertex at t=1/2"
     width="560">

The dashed lines join the controls. The red point splits the curve at `t = 1/2`.
Each physical coordinate is a polynomial of degree at most three in `t`:

```text
X(t) = (1-t)^3 P0 + 3(1-t)^2 t P1 + 3(1-t)t^2 P2 + t^3 P3
```

Placing the controls one-third and two-thirds along the straight segment gives
ordinary linear interpolation. Reversing an edge reverses all four points:
`P3, P2, P1, P0`, with local coordinate `1-t`.

For a child retaining the interval `[a,b]`, its coordinate `s` evaluates the
original curve at `t = a + (b-a)*s`. A split at `t = 1/2` need not divide the
curve's physical length equally.

## Faces and Cells from Shared Edges

Each face is shaped by the curves along its edges. Faces sharing an edge use
the same curve, even when they follow it in opposite directions.

A cubic triangular face needs one additional control point to define its interior
shape. FVMAdapt computes this point from the vertices and edge controls. It
shapes the surface without changing its edges; it is not a mesh vertex or a
splitting point.

For a quadrilateral face, the four edge curves determine the surface between
them.

These triangular and quadrilateral surfaces form the boundaries of tetrahedra,
hexes, prisms, and pyramids. Each cell type has a geometry map that determines
positions inside the cell while matching its boundary surfaces. Both cells
sharing a face use the same surface.

Straight edges reduce these constructions to the ordinary corner maps.
A curved pyramid's apex has a defined position but no single derivative.

## Keeping Geometry through Refinement

For a cell described by a volume map, each child retains its placement in
the original reference shape:

```text
point in child -> corresponding point in original shape -> physical position
```

Further refinement and coarsening keep this relationship, including for children
that have no standard element shape. Refitting a surface or volume through a
child's corners or edges would generally change its interior. Piecewise faces
retain their original patches as well.

General cells retain their shared face maps. In straight 2D refinement, each
accepted triangle or quad has its own corner map.

Neighbors share the complete face and edge geometry and accepted split points.
Different local orientations must give the same physical positions. Surviving
points stay in place.

## Choosing a Center

A centroid is an area or volume average. A split point is where the new edges
meet. On a curved face, the centroid may be off the surface, while the split
point must lie on it.

| Method | Split point |
| --- | --- |
| Reference | Center of the reference shape, mapped to physical space |
| Weighted | Edge-length weighted face center or face-area weighted cell center, located on the face or inside the cell |
| Integrated | Area or volume centroid, located on the face or inside the cell |

For a general cell, the reference center is its retained interior point.
The weighted method uses mapped edge midpoints, endpoint chord lengths and
corner-polygon area vectors; these approximate curved lengths and areas.
Projection places a computed face center on the retained surface without
reshaping it. Available methods depend on the cell and split; see
[splitting cells](@ref fvmadapt_geometric_splitting).

## General Cells

A general input polyhedron is bounded by shared face maps and uses straight
edges, although its faces may be nonplanar. It requires an interior point from
which the whole boundary is visible without leaving the cell.

Splitting creates one child around each vertex. Internal quadrilateral faces
follow their four boundary curves; children need not have a standard element shape.
Nonstandard children of curved standard cells retain their original parent's
map; the straight-edge restriction applies to general input polyhedra.

## Geometry and Output Quality

The output mesh uses nodes and polygon faces whose straight connections
approximate the retained geometry. Children dividing a fixed mapped parent
reproduce its volume within integration accuracy, but coarse and fine polygonal
output meshes can have different volumes.

Map validity and output quality are separate requirements. The output cells
need positive volumes and acceptable face angles and skewness. A center method
determines split placement and does not guarantee good cell quality.
