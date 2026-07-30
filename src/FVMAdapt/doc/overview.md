# FVMAdapt Mesh Adaptation Overview {#fvmadapt_overview}

FVMAdapt takes an existing mesh and a description of requested refinement and
constructs the topology of an adapted mesh. Refinement selection can be supplied
through cell tags, fine-cell tags, geometric source parameters, or XML regions.
An upstream application may derive tags from solution data, but FVMAdapt starts
from supplied refinement inputs.

This page follows the refinement path. FVMAdapt also supports
[plan-based derefinement](@ref fvmadapt_derefinement), which replaces an
eligible family of fine cells with its parent in the refinement tree.

The module has two cooperating layers:

- The lower-level C++ library implements cell, face, and edge splitting,
  temporary refinement trees, plan serialization, and plan replay.
- The Loci rule layer owns the mesh-wide plan stores, coordinates shared faces
  and edges, schedules balancing, and invokes mesh construction and data
  transfer.


## From Refinement Request to Adapted Mesh

FVMAdapt determines the adapted topology before constructing the final
mesh. This final mesh is referred to as the **fine mesh**. The refinement workflow is:

1. Tags or geometric inputs select refinement in the original mesh or in
   previously refined cells.
2. The rule layer creates or updates one cell-plan vector for each original
   cell. An empty vector means that the cell is not split.
3. Each cell plan implies subdivisions of its boundary faces. Contributions
   from adjacent cells are mapped into a common face orientation and merged.
4. The merged face plans imply subdivisions of their boundary edges. Those
   contributions are mapped into a common edge direction and merged.
5. The shared face and edge requirements are fed back into the cell plans.
   Steps 3 through 5 repeat until another pass makes no cell plan change.
6. The stable cell plans and their compatible face and edge plans are replayed
   to construct and number the fine nodes, faces, and cells used by transfer and
   mesh-output operations.

Steps 3 through 5 are *balancing*. Here, balancing means making the shared
topology compatible and enforcing the applicable refinement-grading rules.
Balancing can add cell splits required by those
topological rules; it does not construct the final mesh.

A plan is a sequence of split codes associated with one original cell, face, or
edge. A split code says whether the current tree object is left unsplit or
which split operation it uses. Replaying a plan means reading those codes to
recreate the corresponding parent-and-child objects in memory. An object that
is not split further is a *leaf*. The leaves identify the fine entities that
will be constructed from that original entity.


### A Concrete Cell Plan

Consider the HexCell plan `[7, 1]`. HexCell code `7` splits all three local
directions and creates eight children. Code `1` splits only the local `zeta`
direction and creates two children.

Replay begins with the original cell, called the root:

```text
plan position 0: root consumes code 7
                 → create root children 0 through 7

plan position 1: root child 0 consumes code 1
                 → create children 0.0 and 0.1

end of plan:     no more explicit codes
                 → children 0.0, 0.1, and root children 1 through 7 are leaves
```

The plan is read in breadth-first order: the root is processed first, followed
by its children in child order, followed by the next level(of splitting of the children).
When the stored vector ends, the remaining queued objects are treated as unsplit leaves. The
plan can therefore omit trailing no-split codes.

In this example, replay produces nine leaves: the two children beneath root
child 0 and the seven unsplit root children. The values in `[7, 1]` do not
contain coordinates or final mesh cell numbers. They only say which tree
objects split and how. The geometric split operations create the required
midpoint and center nodes during replay, and later construction rules assign
the final mesh numbering.


## Why General Polyhedral Cells Matter

Consider two cells, `A` and `B`, that share a face. If the plan for `A`
subdivides that face, both sides must use the same fine-face subdivision:

```text
cell A requests a split
        ↓
A's boundary view subdivides the shared face
        ↓
the shared-face plan carries that subdivision to cell B
```

For a cell handled by the general `Cell` path, `B` can sometimes remain one
leaf while its boundary is reconstructed with several fine faces. In that case,
`B` becomes a general polyhedral fine cell with a more detailed boundary; it
does not need a cell split solely to preserve a standard element shape.

This does not guarantee that refinement stops at `B`. If the boundary trees
violate a refinement-grading rule, balancing expands `B`'s cell plan. That
change can subdivide another face and affect another neighbor on the next
balancing pass. General polyhedral cells can therefore limit refinement
propagation caused only by element-shape restrictions, while compatibility and
grading requirements can still propagate refinement.


## Mesh Topology and Refinement Trees

Two different relationships appear throughout the FVMAdapt implementation.

The first is **mesh incidence**:

```text
cell → face → edge → node
```

This describes how unlike mesh entities are connected. A face can be shared by
two cells, and an edge can belong to several faces.

The second is **refinement lineage**:

```text
original entity
├─ child 0
├─ child 1
└─ child 2
```

This describes parent and child objects of the same kind. A cell has child
cells, a face has child faces, and an edge has child edges. This relationship between
a mesh entity and the children/grandchildren forms a tree. These temporary trees are
reconstructed when FVMAdapt interprets or modifies a plan. The plans preserve the
topology between rule operations, allowing FVMAdapt to communicate and adjust it
without keeping the full refinement trees in memory.

Recognized hexahedra and prisms use specialized refinement paths. Tetrahedra,
pyramids, and other polyhedra use the general `Cell` path. The exact child
topology and split-code meanings depend on that path.


## FVMAdapt Documentation

This overview is the entry point for the FVMAdapt documentation:

1. @subpage fvmadapt_geometric_splitting
   shows the children created by the face and cell split operations.
2. @subpage fvmadapt_plans_and_balancing
   develops the plan representation, orientation mappings, balancing loop, and
   reconstruction process in more detail.
3. @subpage fvmadapt_derefinement
   explains how coarsening requests are evaluated, how eligible sibling
   families are collapsed, and how data is transferred to coarsened cells.
4. @subpage fvmadapt_numbering
   defines the local node, edge, face, and direction conventions used by the
   specialized split codes.
