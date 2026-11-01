# FVMAdapt Mesh Adaptation {#fvmadapt_overview}

FVMAdapt refines and coarsens an existing mesh. The solver chooses which cells
to change. FVMAdapt splits or merges them, makes neighboring faces agree,
balances refinement levels, and builds the new mesh.

The public tools and online calls use one parametric implementation.

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
[splitting and the module interface](@ref fvmadapt_reference_subdivision).

The solver receives a mesh of nodes and polygon faces, which approximates the
saved curved geometry. Refinement can change that approximation and its computed
volumes. Geometry accuracy and field conservation therefore need separate checks.

## Code Responsibilities

- `library/parametric/`: local C++ geometry, splitting, coarsening and history.
  It contains no MPI communication or Loci stores.
- FVMAdapt rules and helpers: mesh connectivity, neighbor agreement, balancing,
  communication, output and old/new cell and face correspondence.
- Solver: adaptation criteria, field reconstruction and time histories.

Parallel refinement and coarsening work for standard cells, general cells with
polygon faces and straight 2D extrusions. `marker`, `refmesh` and `refine` use
parametric histories, as do `onlineRefineMesh()` and `initializeGridFromPlan()`.
The interface reference gives the supported combinations.

## Topics

- @subpage fvmadapt_geometric_splitting
- @subpage fvmadapt_plans_and_balancing
- @subpage fvmadapt_derefinement
- @subpage fvmadapt_numbering

- @subpage fvmadapt_reference_geometry
- @subpage fvmadapt_reference_subdivision
