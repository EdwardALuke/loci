# FVMAdapt Mesh Adaptation Overview {#fvmadapt_overview}

FVMAdapt takes refinement requests for an existing mesh and turns them into a
topologically compatible refined mesh. Requests may come from tags, regions,
or solution-based indicators, while mesh-focused rules account for the effects
of each requested refinement on neighboring entities.

The module first describes the requested topology with refinement plans. These
plans record which cells, faces, and edges are split and how their children are
ordered; they do not contain the physical coordinates of the refined mesh. The
plans are balanced across the mesh and then replayed to construct the final
nodes, faces, and cells.

This overview focuses on refinement. FVMAdapt also supports plan-based
derefinement.


## Why General Polyhedral Cells Matter

A central principle of FVMAdapt is that general polyhedral cells can help keep
the effects of refinement local. When one cell subdivides a shared face, its
neighbor must agree with that new face topology. Requiring every affected cell
to remain a standard element could force refinement to cascade through several
layers of neighboring cells.

FVMAdapt can instead represent affected cells as general polyhedra when needed.
This lets a neighboring cell incorporate a subdivided shared face without being
split only to preserve a standard cell shape. Balancing may still propagate
refinement where compatibility or refinement-depth rules require it, but
general cells reduce propagation caused solely by element-shape restrictions.


## Geometric Model

The basic mesh-incidence chain is:

```text
cell → face → edge → node
```

Cells drive the refinement process, while their faces and edges determine how
the resulting topology must match neighboring cells. FVMAdapt handles standard
cell shapes such as hexahedra, prisms, tetrahedra, and pyramids, as well as
general polyhedral cells.

Many refinement paths use midpoint-based construction: edges receive midpoint
nodes, faces receive center nodes, and cells receive center nodes as required by
the split. The exact nodes, child faces, and child cells depend on the entity
type and split code.

<div class="term-box">
<div class="term-title">Splitting</div>
A split subdivides an entity into child entities. Depending on the entity and
split code, this may create new nodes, edges, faces, or cells.
</div>


## Refinement Trees and Plans

The lower-level refinement library can be understood as a collection of
replayable trees for cells, faces, and edges. A tree begins with a root object
representing an entity in the original mesh. Splitting that object creates its
children, and splitting a child extends the tree to another refinement level.
An object with no children is a leaf. After replay, the leaves represent the
fine entities produced by that refinement tree.

These trees are temporary working objects rather than the permanent mesh-wide
representation. Compact plans describe how to reconstruct them:

- A cell plan describes a cell refinement tree.
- A face plan describes a face refinement tree.
- An edge plan describes an edge refinement tree.

Plans are compact sequences of split codes stored in breadth-first order.
Breadth-first means that the plan describes one level of an entity's tree before
continuing to the next level. The plan is therefore a recipe for replaying the
tree, not the tree itself and not the final refined mesh.


### A Cell Plan Example

Consider an illustrative cell whose first split creates four children:

```text
cell
├─ child 0
├─ child 1
├─ child 2
└─ child 3
```

If child 1 is split again, the tree becomes:

```text
cell
├─ child 0
├─ child 1
│  ├─ child 1.0
│  └─ child 1.1
├─ child 2
└─ child 3
```

The child counts in this example are illustrative; the actual count depends on
the cell type and split code. A cell plan stores the split codes needed to
recreate this structure in breadth-first order. Creating a plan serializes an
existing temporary tree, while replaying a plan reconstructs that tree. Each
cell type defines the split codes that its plans can contain.


## From Refinement Request to Adapted Mesh

The module has two cooperating layers. The lower-level C++ library creates,
splits, serializes, and replays individual entity trees. The Loci rule layer
coordinates plans across the mesh, balances neighboring entities, and invokes
the construction and data-transfer operations that produce the refined mesh.

The overall flow is:

1. Tags, region inputs, or solution-based indicators identify where refinement
   is requested.
2. Cell plans are created or updated to represent those requests.
3. Cell plans are balanced to satisfy mesh-wide compatibility requirements.
4. Face plans are derived from cell plans and merged on shared faces.
5. Edge plans are derived from face plans and merged on shared edges.
6. The compatible plans are replayed to generate fine nodes, faces, cells,
   numbering, transfer data, and data used by later mesh-output operations.

The plan representation itself is topological and does not store physical
coordinates. Geometric measurements can still influence which splits are
selected, and coordinates and mesh numbering are used when the plans are
replayed. The collection of plans across all cells carries the topological
information needed for balancing and final mesh generation.


## Why Balancing Is Necessary

An initial cell plan represents requested refinement for one cell, but that
request cannot be applied independently of its neighbors. Cells that share a
face must agree on how that face is subdivided, and the faces that share an edge
must agree on the edge subdivision. Each cell may view a shared face with a
different local orientation, so the orientation mapping must identify which
child entities correspond.

Balancing adjusts requested cell plans to satisfy these shared-topology and
refinement-depth requirements, such as limits on refinement-depth differences
between neighboring entities. During this process, temporary face and edge
information is derived from the current cell plans and used to determine
whether those cell plans must expand. The process continues until the cell
plans no longer need to change. Compatible face and edge plans are then derived
from the balanced cell plans.

The result of balancing is a compatible set of cell, face, and edge plans, not
the fine mesh itself. Those plans are replayed later during construction,
numbering, transfer, and output.


## Further Reference

See [FVMAdapt Geometric Splitting Reference](@ref fvmadapt_geometric_splitting)
for the isotropic and anisotropic split behavior of faces and cells, including
child topology and split codes.

See [FVMAdapt Refinement Plans and Balancing](@ref fvmadapt_plans_and_balancing)
for the relationship between temporary refinement trees, cell, face, and edge
plans, and the iterative balancing process.

See [FVMAdapt Numbering Conventions](@ref fvmadapt_numbering) for the local node,
edge, and face numbering used to interpret split codes, child ordering, and
face orientation.
