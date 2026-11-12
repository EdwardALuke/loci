# FVMAdapt Mesh Adaptation Overview {#fvmadapt_overview}

FVMAdapt takes an existing mesh and a description of requested refinement and
constructs the topology of an adapted mesh. Refinement selection can be supplied
through cell tags, fine-cell tags, geometric source parameters, or XML regions.
An upstream application may derive tags from solution data, but FVMAdapt starts
from supplied refinement inputs.

This page follows the refinement path. FVMAdapt also supports
[derefinement](@ref fvmadapt_derefinement), which replaces an
eligible family of fine cells with its parent in the refinement tree.

The module has two cooperating layers:

- The lower-level C++ library implements cell, face, and edge splitting,
  geometry evaluation, and refinement histories.
- The Loci rule layer owns the mesh-wide history stores, coordinates shared faces
  and edges, schedules balancing, and invokes mesh construction and data
  transfer.

VOG supplies node positions and connectivity. Optional controls prescribe
cubic edges; the output consists of nodes and polygon faces.


## From Refinement Request to Adapted Mesh

FVMAdapt determines the adapted topology before constructing the final
mesh. This final mesh is referred to as the **fine mesh**. The refinement workflow is:

1. **Describe the original geometry.** A geometry map converts local
   coordinates into physical positions. A unit cube, for example, can
   describe a distorted hex.
2. **Apply the requested splits or merges.** Split cells in local coordinates
   or restore a parent in place of its complete family of leaf children.
3. **Match neighbors and balance levels.** Combine shared face and edge
   divisions, adding cell splits where the grading rules require them.
4. **Construct the output mesh.** Number its nodes, faces and cells, and use
   the saved geometry to place the points. Keep the accepted history for
   subsequent adaptations.

Physical geometry is also evaluated when choosing centers and checking splits.
The final construction stage assembles the output mesh from the accepted result.


## Why General Polyhedral Cells Matter

Consider two cells, `A` and `B`, that share a face. If the subdivision of `A`
divides that face, both sides must use the same fine-face subdivision.
Cell `B` can remain one leaf while its boundary is reconstructed with several
fine faces:

```text
  refined neighbor       coarse cell
  +---------+----------------------+
  | child 0 |                      |
  +---------+      one cell        |
  | child 1 |                      |
  +---------+----------------------+
            two shared face pieces
```

It need not split solely to preserve a standard element shape. Balancing can
still require a split to limit the difference in refinement levels.


## Mesh Topology and Refinement Trees

Two different relationships appear throughout the FVMAdapt implementation.

The first is **mesh incidence**:

```text
cell → face → edge → node
```

This describes how unlike mesh entities are connected. A face can be shared by
two cells, and an edge can belong to several faces.

The second is the parent-and-child relationship created by refinement:

```text
original cell
├─ child 0
├─ child 1
└─ ...
```

The original cell is the *root*. A current cell that is not split further is
a *leaf*. The leaves form the adapted mesh.


## FVMAdapt Documentation

- @subpage fvmadapt_reference_geometry
- @subpage fvmadapt_geometric_splitting
- @subpage fvmadapt_plans_and_balancing
- @subpage fvmadapt_derefinement
- @subpage fvmadapt_numbering
