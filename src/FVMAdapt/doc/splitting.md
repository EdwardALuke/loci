# FVMAdapt Geometric Splitting Reference {#fvmadapt_geometric_splitting}

This page expands the splitting terminology used by the
[FVMAdapt Mesh Adaptation Overview](@ref fvmadapt_overview). It collects the
detailed isotropic and directional split operations so the overview can stay
focused on the main adaptation flow.

The figures use idealized reference shapes to expose topology. They do not
imply that physical mesh cells are regular or orthogonal. The `xi`, `eta`, and
`zeta` directions below are the local directions defined in the
[numbering conventions](@ref fvmadapt_numbering), not global Cartesian axes.

This page uses *face center* and *cell center* for the nodes returned by
FVMAdapt's current `centroid()` helpers. With `CENTROID = 1`, a face center is
the edge-length-weighted average of its edge midpoints, and a cell center is
the face-area-weighted average of its face centers. These points are not
exact geometric centroids.


## Isotropic Face Splitting

The general polygonal `Face` path splits an n-edge face into n quadrilateral
children:

1. Split each boundary edge at its midpoint.
2. Create one face-center node.
3. Connect the face center to every edge midpoint.
4. Each original vertex forms one quadrilateral child face with the two
   adjacent edge midpoints and the face center.

The specialized `QuadFace` class also supports directional splits used by
`HexCell` and `Prism`. The examples in this section show only the isotropic
operation on a general `Face`.

| Parent face | Library handling | Isotropic split result |
| --- | --- | --- |
| Triangle | General polygonal `Face` | 3 quadrilateral child faces |
| Quadrilateral | General polygonal `Face` | 4 quadrilateral child faces |
| Pentagon | General polygonal `Face` | 5 quadrilateral child faces |
| n-sided polygon | General polygonal `Face` | n quadrilateral child faces |

### Triangle Face

<img src="figures/faces/triangle/triangle_face_isotropic_split.svg"
     alt="Triangle face split into three quadrilateral child faces"
     width="300">

The parent has three edges, so the split creates three quadrilateral children.
One child is highlighted below.

<img src="figures/faces/triangle/triangle_face_isotropic_split_quad_face.svg"
     alt="One quadrilateral child face highlighted inside a split triangle"
     width="200">

### Quadrilateral Face

<img src="figures/faces/quadrilaterial/quadrilateral_face_isotropic_split.svg"
     alt="Quadrilateral face split into four quadrilateral child faces"
     width="300">

The parent has four edges, so the split creates four quadrilateral children.
One child is highlighted below.

<img src="figures/faces/quadrilaterial/quadrilaterial_face_isotropic_split_quad_face.svg"
     alt="One quadrilateral child face highlighted inside a split quadrilateral"
     width="200">

### Pentagon Face

<img src="figures/faces/pentagon/pentagon_face_isotropic_split.svg"
     alt="Pentagon face split into five quadrilateral child faces"
     width="300">

The parent has five edges, so the split creates five quadrilateral children.
One child is highlighted below.

<img src="figures/faces/pentagon/pentagon_face_isotropic_split_quad_face.svg"
     alt="One quadrilateral child face highlighted inside a split pentagon"
     width="200">

### General Polygonal Face

<img src="figures/faces/polyhedra/polyhedra_face_isotropic_split.svg"
     alt="General polygonal face split into quadrilateral child faces"
     width="300">

An n-edge parent produces n quadrilateral children. One child is highlighted
below.

<img src="figures/faces/polyhedra/polyhedra_face_isotropic_split_quad_face.svg"
     alt="One quadrilateral child face highlighted inside a split general polygon"
     width="200">

## Isotropic Cell Splitting

FVMAdapt does not send every cell shape through one generic splitter.
Cells classified as hexahedra use `HexCell`, cells classified as prisms use
`Prism`, and other polyhedra use the general `Cell` path. The resulting child
objects and topologies therefore depend on the parent path.

### Hexahedron

An isotropic `HexCell` split uses code `7`, which splits all three local
directions and creates eight `HexCell` children.

<img src="figures/volume_cells/hexahedron/isotropic/hexahedron_cell.svg"
     alt="Hexahedron cell"
     width="300">

Each boundary `QuadFace` is split in both face-local directions. The
construction uses the edge midpoints, six face centers, and one cell center to
form the child connectivity. In the figure below, red points mark the inserted
edge-midpoint and face-center nodes, while the green point marks the cell
center. Red segments connect edge midpoints to face centers, and green segments
connect face centers to the cell center.

<img src="figures/volume_cells/hexahedron/isotropic/hexahedron_midpoint_based_isotropic.svg"
     alt="Hexahedron midpoint-based isotropic split"
     width="300">

The animation separates the eight children.

<img src="figures/volume_cells/hexahedron/isotropic/hexahedron_midpoint_based_isotropic_exploded_motion_yaw_wire.gif"
     alt="Hexahedron midpoint-based isotropic split"
     width="600">

### Tetrahedron

A tetrahedron follows the general `Cell` path. `Cell::split()` creates one
`DiamondCell` at each original vertex. Every tetrahedron vertex has three
incident edges, so the result is four 3-fold `DiamondCell` children. A 3-fold
`DiamondCell` has hexahedral connectivity, but it remains a `DiamondCell`
object rather than a `HexCell`.

<img src="figures/volume_cells/tetrahedron/tetrahedron_cell.svg"
     alt="Tetrahedron cell"
     width="300">

<img src="figures/volume_cells/tetrahedron/tetrahedron_midpoint_based_isotropic.svg"
     alt="Tetrahedron midpoint-based isotropic split"
     width="300">

The animation separates the four children.

<img src="figures/volume_cells/tetrahedron/tetrahedron_midpoint_based_isotropic_exploded_motion_yaw_wire.gif"
     alt="Tetrahedron midpoint-based isotropic split"
     width="600">

### Triangular Prism

A classified triangular prism uses the specialized `Prism` path with
`nfold = 3`. Isotropic split code `3` splits both prism directions and creates
`2*nfold`, or six, children. The implementation represents each child as a
`Prism` with `nfold = 4`; each therefore has hexahedral connectivity.

<img src="figures/volume_cells/prism/isotropic/prism_cell.svg"
     alt="Prism cell"
     width="300">

<img src="figures/volume_cells/prism/isotropic/prism_midpoint_based_isotropic.svg"
     alt="Prism midpoint-based isotropic split"
     width="300">

The animation separates the six children.

<img src="figures/volume_cells/prism/isotropic/prism_midpoint_based_isotropic_exploded_motion_yaw_wire.gif"
     alt="Prism midpoint-based isotropic split"
     width="600">

### Pyramid

A pyramid also follows the general `Cell` path. Its four base vertices each
produce a 3-fold `DiamondCell`, while its four-edge apex produces one 4-fold
`DiamondCell`. The result is five children: four with hexahedral connectivity
and one 4-fold diamond.

<img src="figures/volume_cells/pyramid/pyramid_cell.svg"
     alt="Pyramid cell"
     width="300">

<img src="figures/volume_cells/pyramid/pyramid_midpoint_based_isotropic.svg"
     alt="Pyramid midpoint-based isotropic split"
     width="300">

The animation separates the five children.

<img src="figures/volume_cells/pyramid/pyramid_midpoint_based_isotropic_exploded_motion_yaw_wire.gif"
     alt="Pyramid midpoint-based isotropic split"
     width="600">

### General Cell and DiamondCell Topology

The general `Cell` split first splits every boundary face. It then creates one
cell-center node, connects that node to every face center, and creates one
interior quadrilateral face for every original edge. Finally, it creates one
`DiamondCell` child for every original vertex.

A general parent `Cell` does not have one fold value. Each original vertex
produces a `DiamondCell` whose fold equals that vertex's incident-edge count,
so one split can produce children with different folds.

The fold of a child is the number of original edges incident on its
corresponding vertex. A fold `n` `DiamondCell` has:

- `2*n + 2` nodes
- `4*n` edges
- `2*n` quadrilateral faces

Only the 3-fold case has hexahedral connectivity: eight nodes and six
quadrilateral faces. That topological description does not imply right angles
or identify the object as a `HexCell`.

| Diamond fold | Nodes | Edges | Faces | Example |
| --- | --- | --- | --- | --- |
| 3-fold | 8 | 12 | 6 | Tetrahedron child or pyramid base-vertex child |
| 4-fold | 10 | 16 | 8 | Pyramid apex child |
| n-fold | `2*n + 2` | `4*n` | `2*n` | General-cell child at an n-valent vertex |

### Isotropic Split Summary

| Parent shape | Refinement path | Local split rule | Isotropic child topology |
| --- | --- | --- | --- |
| Hexahedron | `HexCell` | Code `7`: split all three local directions | 8 `HexCell` children |
| Tetrahedron | General `Cell` | One fold-3 child per vertex | 4 `DiamondCell` children with hexahedral connectivity |
| Triangular prism | `Prism` | Code `3`: split both prism directions | 6 `Prism` children with `nfold = 4` |
| Pyramid | General `Cell` | One child per vertex | 4 fold-3 and 1 fold-4 `DiamondCell` children |
| Other polyhedron | General `Cell` | One child per vertex | `DiamondCell` fold follows the vertex valence |


## Directional Cell Splitting

`HexCell` and `Prism` support split codes that refine selected local
directions. This section defines what each code constructs. The policy that
selects a code from tags, geometry, or the configured split mode is separate
from these code semantics.

### HexCell Split Codes

A `HexCell` code is a three-bit mask: bit `0` selects `zeta`, bit `1` selects
`eta`, and bit `2` selects `xi`. For a nonzero code, the number of children is
`2^k`, where `k` is the number of selected directions.

| Code | Binary | Split directions | Children |
| --- | --- | --- | --- |
| 0 | `000` | none | 0; the cell remains a leaf |
| 1 | `001` | `zeta` | 2 `HexCell` children |
| 2 | `010` | `eta` | 2 `HexCell` children |
| 3 | `011` | `eta`, `zeta` | 4 `HexCell` children |
| 4 | `100` | `xi` | 2 `HexCell` children |
| 5 | `101` | `xi`, `zeta` | 4 `HexCell` children |
| 6 | `110` | `xi`, `eta` | 4 `HexCell` children |
| 7 | `111` | `xi`, `eta`, `zeta` | 8 `HexCell` children |

Code `7` is the isotropic `HexCell` split shown earlier. The remaining
directional codes are illustrated below.

#### Code 1: Zeta

Code `1` splits only the `zeta` direction and produces two children.

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_1_zeta.svg"
     alt="Hexahedron directional split in zeta, code 1"
     width="300">

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_1_zeta_exploded_motion_yaw_wire.gif"
     alt="Children of the HexCell zeta split, code 1"
     width="600">

#### Code 2: Eta

Code `2` splits only the `eta` direction and produces two children.

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_2_eta.svg"
     alt="Hexahedron directional split in eta, code 2"
     width="300">

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_2_eta_exploded_motion_yaw_wire.gif"
     alt="Children of the HexCell eta split, code 2"
     width="600">

#### Code 3: Eta and Zeta

Code `3` combines the `eta` and `zeta` directions and produces four children.

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_3_eta_zeta.svg"
     alt="Hexahedron directional split in eta and zeta, code 3"
     width="300">

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_3_eta_zeta_exploded_motion_yaw_wire.gif"
     alt="Children of the HexCell eta-zeta split, code 3"
     width="600">

#### Code 4: Xi

Code `4` splits only the `xi` direction and produces two children.

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_4_xi.svg"
     alt="Hexahedron directional split in xi, code 4"
     width="300">

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_4_xi_exploded_motion_yaw_wire.gif"
     alt="Children of the HexCell xi split, code 4"
     width="600">

#### Code 5: Xi and Zeta

Code `5` combines the `xi` and `zeta` directions and produces four children.

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_5_xi_zeta.svg"
     alt="Hexahedron directional split in xi and zeta, code 5"
     width="300">

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_5_xi_zeta_exploded_motion_yaw_wire.gif"
     alt="Children of the HexCell xi-zeta split, code 5"
     width="600">

#### Code 6: Xi and Eta

Code `6` combines the `xi` and `eta` directions and produces four children.

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_6_xi_eta.svg"
     alt="Hexahedron directional split in xi and eta, code 6"
     width="300">

<img src="figures/volume_cells/hexahedron/anisotropic/hexahedron_split_code_6_xi_eta_exploded_motion_yaw_wire.gif"
     alt="Children of the HexCell xi-eta split, code 6"
     width="600">

### Prism Split Codes

A `Prism` has two corresponding n-sided end faces and `nfold`
quadrilateral side faces. A classified triangular prism has `nfold = 3`.
`Prism::nfold` counts the sides of the prism's end faces; it is distinct from
the vertex-based fold of a `DiamondCell`.

| Code | Split operation | Children |
| --- | --- | --- |
| 0 | none | 0; the prism remains a leaf |
| 1 | axial | 2 `Prism` children with the parent's `nfold` |
| 2 | transverse | `nfold` `Prism` children, each with `nfold = 4` |
| 3 | axial and transverse | `2*nfold` `Prism` children, each with `nfold = 4` |

#### Code 1: Axial

The axial split bisects the edges that join corresponding end-face vertices.
Their midpoint nodes form one new interior n-sided face. A triangular prism
therefore produces two triangular-prism children.

<img src="figures/volume_cells/prism/anisotropic/prism_split_code_1_axial.svg"
     alt="Prism axial directional split, code 1"
     width="300">

<img src="figures/volume_cells/prism/anisotropic/prism_split_code_1_axial_exploded_motion_yaw_wire.gif"
     alt="Children of the prism axial split, code 1"
     width="600">

#### Code 2: Transverse

The transverse split subdivides both end faces, splits each side face in its
transverse direction, and connects the two end-face centers with a new
interior edge. A triangular prism produces three `Prism` children with
`nfold = 4`, so each child has hexahedral connectivity.

<img src="figures/volume_cells/prism/anisotropic/prism_split_code_2_transverse.svg"
     alt="Prism transverse directional split, code 2"
     width="300">

<img src="figures/volume_cells/prism/anisotropic/prism_split_code_2_transverse_exploded_motion_yaw_wire.gif"
     alt="Children of the prism transverse split, code 2"
     width="600">

#### Code 3: Axial and Transverse

Code `3` combines both operations. A triangular prism produces six `Prism`
children with `nfold = 4`. This is the isotropic prism split shown earlier.


## Further Reference

See [FVMAdapt Numbering Conventions](@ref fvmadapt_numbering) for the local
directions and orientation mappings used to interpret these split codes.

See [FVMAdapt Refinement Plans and Balancing](@ref fvmadapt_plans_and_balancing)
for how split codes are stored in plans and coordinated across shared faces
and edges.
