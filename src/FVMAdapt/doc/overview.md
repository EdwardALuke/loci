# FVMAdapt Mesh Adaptation {#fvmadapt_overview}

FVMAdapt refines and coarsens an existing mesh. The solver chooses which cells
to change. FVMAdapt splits or merges them, makes neighboring faces agree,
balances refinement levels, and builds the new mesh.

Refinement selection can be supplied through cell tags, node tags, geometric
source parameters, or XML regions. An upstream application may derive tags
from solution data, but FVMAdapt starts from supplied refinement inputs.

## Mesh Topology and Refinement Trees

Two different relationships appear throughout FVMAdapt. The first describes
how mesh entities are connected:

```text
cell → face → edge → node
```

A cell is bounded by faces, a face by edges, and an edge by two nodes.
A face can be shared by two cells, and an edge can belong to several faces.

The second is the parent-and-child relationship created by refinement:

```text
original cell
├─ child 0
├─ child 1
└─ ...
```

The original cell is the *root*. A current cell that is not split further is
a *leaf*. The leaves form the adapted mesh. The refinement history records
these relationships and the accepted split locations, so later operations can
reconstruct the same cells.

## How Adaptation Works

1. **Keep the original geometry.** Describe each cell and its shared faces and
   edges. A geometry map converts local coordinates into physical positions.
   For example, a unit cube can describe a distorted hex.
2. **Split the requested cells.** Divide them in local coordinates. Use the
   saved geometry to place new points. Each cell type has suitable split patterns.
3. **Match neighbors and balance levels.** Neighbors use the same points and
   pieces of a shared face. Split additional cells where the permitted level
   difference requires it.
4. **Build and check the mesh.** Assign nodes, faces and cells their output
   numbers and preserve boundary labels. The solver then transfers its fields.
5. **Keep the history.** Save the accepted splits and geometry for the next
   adaptation and for restart.

Physical geometry is also evaluated when choosing centers and checking splits.
The final construction stage assembles the output mesh from the accepted result.

Coarsening restores a parent when all its immediate children can be removed.
Face divisions still needed by a finer neighbor remain:

```text
  refined neighbor       coarse cell
  +---------+----------------------+
  | child 0 |                      |
  +---------+      one cell        |
  | child 1 |                      |
  +---------+----------------------+
            two shared face pieces
```

A coarse cell can therefore remain one cell with several smaller faces on its
boundary. It need not split solely to preserve a standard element shape.
Balancing can still require a split to limit the difference in refinement levels.

Existing points stay in place. Both neighbors use the same saved face geometry
and split points, even when their cell types differ. Reloading history uses
those accepted points; it does not choose new centers.

## Splits and Centers

The familiar isotropic construction connects edge midpoints, face centers and
a cell center to form one child around each vertex. Other patterns cut along
local directions or preserve a cell type, such as dividing a tet into eight tets.
These patterns need not produce the same mesh.

For splits that use centers, the caller can choose a reference center, a weighted
center, or an area/volume centroid. A curved face's centroid may lie off its
surface; the split point must lie on it. The geometry map supplies that position.
See [cell geometry](@ref fvmadapt_reference_geometry) and
[splitting](@ref fvmadapt_geometric_splitting).

The solver receives a mesh of nodes and polygon faces, which approximates the
saved curved geometry. Refinement can change that approximation and its computed
volumes. Geometry accuracy and field conservation therefore need separate checks.

## Code Responsibilities

- `library/parametric/`: local C++ geometry, splitting, coarsening and history.
  It contains no MPI communication or Loci stores.
- FVMAdapt rules and helpers: mesh connectivity, neighbor agreement, balancing,
  communication, output and old/new cell and face correspondence.
- Solver: adaptation criteria, field reconstruction and time histories.

The public tools and online calls use the same geometry and history. The
interface reference gives the supported combinations of cells and split patterns.

## Topics

- @subpage fvmadapt_reference_geometry
- @subpage fvmadapt_geometric_splitting
- @subpage fvmadapt_plans_and_balancing
- @subpage fvmadapt_derefinement
- @subpage fvmadapt_numbering
- @subpage fvmadapt_reference_subdivision
