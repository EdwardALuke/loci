# Cell Geometry {#fvmadapt_reference_geometry}

A geometry map converts a point in a simple reference shape into its physical
position. Refinement divides that shape and keeps the map.

## Geometry from VOG

VOG supplies node positions and face connectivity. The importer uses them as
follows:

| Shape | Geometry from its corners |
| --- | --- |
| Triangle, tet | Linear |
| Quad | Bilinear; the corners need not be coplanar |
| Polygon with more than four vertices | Rational surface over a regular reference polygon |
| Hex | Trilinear |
| Triangular prism | Linear over the triangle and along its axis |
| Pyramid | Bilinear base joined to the apex |
| General polyhedron | Shared face maps bounding one cell |

`CornerMap` implements the standard shapes. `BezierMap` supports curved segments,
quads and hexes with additional control points. VOG corners alone do not define
an intended CAD surface.

`PolygonMap` interpolates the polygon's corners and straight edges using
[Wachspress coordinates](https://pmc.ncbi.nlm.nih.gov/articles/PMC3549276/).
Its children retain this surface, including its curvature when the corners
are not coplanar.

`referenceElement(shape)` gives the local vertices, edges and faces.
`GeometryMap::evaluate(q)` gives the physical position and, where available,
derivatives. Local directions need not align with physical X, Y and Z.

## Keeping a Surface through Refinement

A child is evaluated through its parent (`RestrictedMap`):

```text
point in child -> corresponding point in parent -> physical position
```

For a warped quad `X(u,v) = (u,v,u*v/4)`, the center of its left-half child has
parent coordinates `(1/4,1/2)` and physical position `(1/4,1/2,1/32)`.
Further refinement uses this same surface, rather than fitting another surface
through the child's corners.

Neighbors share the complete face and edge geometry and the accepted split
points. Different local orientations must give the same physical positions.
Ordinary refinement and coarsening leave surviving points in place.

The planar 2D path instead retains each accepted triangle or quad with straight
edges. Its extrusion keeps the front and back points aligned and the thickness
unchanged.

## Choosing a Center

A centroid is an area or volume average. A split point is where the new edges
meet. On a curved face, these can be different: the centroid may be off the
surface and must be located on it before splitting.

| `CenterMethod` | Split point |
| --- | --- |
| `reference` | Center of the reference shape |
| `weighted` | Edge-length weighted face center or face-area weighted cell center, located on or inside the map |
| `integrated` | Area or volume centroid, located on or inside the map |

The weighted face calculation uses mapped edge midpoints and endpoint chord
lengths. Cell weights use the area vectors of the corner polygons. These are
approximations, not curved-edge lengths or exact curved-face areas.

```cpp
auto q = faceCenterCoordinates(surface, CenterMethod::integrated);
auto splitPoint = surface.evaluate(q).position;
```

Projection and numerical integration have accuracy and work limits. Failure to
find an acceptable point rejects the request. Accepted coordinates are saved;
changing the preferred method does not move a shared split already in use.
Supported methods depend on the split, as listed in the interface reference.

## General Cells

`CornerBoundary::splitCell()` uses a cell's shared face maps and creates one
child around each vertex. Internal faces follow their four edge curves.
`CornerSubdivision` retains the parent/child history for pyramids and general
cells. Children need not have a standard element shape.

The construction requires an interior point from which the whole boundary is
visible without leaving the cell. The importer uses the vertex average and
rejects cells that fail its geometry checks. Arbitrary non-star-shaped cells
are unsupported. Sampled derivatives do not prove validity everywhere.

## Geometry, Output and Quality

`geometryMeasure()` computes length, area or volume and its centroid;
`MeasureOptions` sets accuracy and work limits. Check its convergence result.
`geometryQuality()` and `splitQuality()` assess the mapped shapes.
When a cell is described by its faces, volume integration treats the linear
direction toward its interior point exactly. Numerical integration remains on
the faces.

The output mesh uses nodes and polygon faces. Straight connections can
approximate curved connections differently before and after refinement.
Children dividing a fixed mapped parent must reproduce its volume within
integration accuracy; the coarse and fine output meshes can have different
computed volumes. Field transfer must use the solver's measures consistently.

Check the output cells as well as the maps: positive volumes, face angles,
skewness and solver behavior. Small pieces used only for integration are not
necessarily small solver cells.

Keep the original geometry and accepted history for subsequent adaptation.
Re-importing an output VOG starts new geometry and loses that history.
