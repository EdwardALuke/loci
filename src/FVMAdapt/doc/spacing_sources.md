# Spacing sources

Set `RefinementOptions.sourceFile` or the `parametricSourceFile` fact to a text
file. Use this or an XML region, alongside optional cell/node tags. Sources are
applied after tags and before balancing. They request refinement, not coarsening.
Clear the filename when replaying saved history without new source requests.
Direct rule callers also supply `parametricRequests`, empty for source-only refinement.

The file contains a positive source count followed by eleven numbers per source:

```text
count
x1 y1 z1   x2 y2 z2   r0 s0   r1 s1   a
```

The two points define a segment. At distance `r` from that segment, spacing is:

- `s0` inside radius `r0`, strictly between the segment's end planes;
- `s1` elsewhere within radius `r1`, including the end caps;
- `s1 * (r/r1)^a` farther away.

Coincident endpoints define a point source and use the outer spacing `s1`.
Radii must satisfy `0 <= r0 <= r1` with `r1 > 0`; both spacings must be positive
and `a >= 0`. All numbers must be finite. Lengths use the grid's physical units.

This segment requests spacing 0.05 in its inner cylinder and 0.2 in its outer
capsule, with quadratic growth outside:

```text
1
0 0 -1   0 0 1   0.1 0.05   0.5 0.2   2
```

A cell samples all sources at its own vertices and their average position. It
refines when half its shortest edge is at least the smallest sampled spacing.
The existing small-feature check also refines while the shortest edge exceeds
`r1/2` and either an endpoint lies in the cell's bounding box or the source
midpoint is closer than `shortest edge + r1` to the average vertex position.
This is a proximity estimate, not an exact segment/cell intersection.

Children are checked recursively in the same proposal. The selected split
family and center method apply; automatic directions are optional for native
3D cells. Extrusion ignores thickness in the edge-size check but samples sources
in physical 3D coordinates. `minEdgeLength`, if positive, limits these requests;
balancing may refine further. Invalid files or exhausted depth/node limits
reject the proposal and retain the accepted mesh.
