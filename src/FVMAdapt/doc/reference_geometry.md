# Cell Geometry {#fvmadapt_reference_geometry}

A geometry map converts a point in a simple reference shape into its physical
position. The reference shape and map form a parametric representation.
For example, local coordinates in a unit cube describe positions inside a
distorted hex. Refinement divides that shape and keeps the map.
Local directions need not align with physical X, Y and Z.

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

The library also supports curved segments, quads and hexes defined by additional
control points. VOG corners alone do not define an intended CAD surface.

A polygon's surface interpolates its corners and straight edges. Its children
retain this surface, including its curvature when the corners are not coplanar.

An imported face may have extra vertices along an edge. When neighboring edges
are nearly collinear, a smooth map can fold at that corner. Such a face instead
retains triangles joining its edges to one center. Both neighbors use the same
construction. Its cell uses the shared-boundary representation for 3D refinement;
original node positions are unchanged. A planar quad keeps its bilinear map
when a positive-interior check permits its straight corner.

## Keeping a Surface through Refinement

A child is evaluated through its parent:

```text
point in child -> corresponding point in parent -> physical position
```

For a warped quad `X(u,v) = (u,v,u*v/4)`, the center of its left-half child has
parent coordinates `(1/4,1/2)` and physical position `(1/4,1/2,1/32)`.
Further refinement uses this same surface, rather than fitting another surface
through the child's corners.

For a piecewise surface, each child keeps its boundary in the source face's
coordinates. Area and volume integration follow the retained pieces.

Neighbors share the complete face and edge geometry and the accepted split
points. Different local orientations must give the same physical positions.
Ordinary refinement and coarsening leave surviving points in place.

The planar 2D path instead retains each accepted triangle or quad with straight
edges. Its extrusion keeps the front and back points aligned and the thickness
unchanged.

## Choosing a Center

A centroid is an area or volume average. A split point is where the new edges
meet. On a curved face, these can be different: the centroid may be off the
surface. The split point must lie on the surface.

| Method | Split point |
| --- | --- |
| Reference | Center of the reference shape, mapped to physical space |
| Weighted | Edge-length weighted face center or face-area weighted cell center, located on the face or inside the cell |
| Integrated | Area or volume centroid, located on the face or inside the cell |

For a general cell, the reference center is its stored interior point.

The weighted face calculation uses mapped edge midpoints and endpoint chord
lengths. Cell weights use the area vectors of the corner polygons. These are
approximations, not curved-edge lengths or exact curved-face areas.

Locating a computed face center on the saved surface is a projection. It chooses
a point on that surface; it does not reshape the face. Evaluating the geometry
map then gives the physical position of the accepted local point.

Projection and numerical integration have accuracy and work limits. Failure to
find an acceptable point rejects the request. Accepted coordinates are saved;
changing the preferred method does not move a shared split already in use.
Supported methods depend on the split, as listed in the
[interface reference](@ref fvmadapt_reference_subdivision).

## General Cells

A general cell is described by its shared face maps and an interior point,
rather than by one standard reference shape. Splitting creates one child around
each vertex. Internal faces follow their four edge curves. Children need not
have a standard element shape.

When the boundary is piecewise or has a straight corner, internal faces join
those edge curves to one center. This preserves the curves without requiring a
smooth tangent at that corner.

The construction requires an interior point from which the whole boundary is
visible without leaving the cell. The importer uses the vertex average and
rejects cells that fail its geometry checks. Arbitrary non-star-shaped cells
are unsupported. Sampled derivatives do not prove validity everywhere.
If a child's vertex average fails, the library tries a point along its original
corner-to-center diagonal. This point must pass the same geometry checks.

## Geometry, Output and Quality

The output mesh uses nodes and polygon faces. Straight connections can
approximate curved connections differently before and after refinement.
Children dividing a fixed mapped parent must reproduce its volume within
integration accuracy; the coarse and fine output meshes can have different
computed volumes. Field transfer must use the solver's measures consistently.

Check the output cells as well as the maps: positive volumes, face angles,
skewness and solver behavior. A center method defines how a split is placed;
it does not by itself guarantee good cell quality.

A map can have zero derivatives at a boundary while describing a valid cell.
These cases require a positive-interior proof; source edges and faces must
remain nonzero. Quality checks still report the boundary's poor conditioning.

Keep the original geometry and accepted history for subsequent adaptation.
Re-importing an output VOG starts new geometry and loses that history.
