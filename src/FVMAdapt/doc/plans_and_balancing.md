# FVMAdapt Refinement Plans and Balancing {#fvmadapt_plans_and_balancing}

FVMAdapt records the topology it intends to construct before it creates the
fine mesh. The recorded plans let the module rebuild and adjust temporary
refinement trees without keeping all of those trees in memory between Loci rule
operations.

This page follows the refinement path.
[Derefinement](@ref fvmadapt_derefinement) uses related plan machinery, but it
has different request, sibling-family, and edge-depth conditions. Here,
*balancing* means making the subdivisions of shared faces and edges compatible
and enforcing refinement-grading rules. It does not mean processor load
balancing.

The [mesh-adaptation overview](@ref fvmadapt_overview) introduces the complete
workflow. The
[geometric splitting reference](@ref fvmadapt_geometric_splitting) shows what
the split operations create, and the
[numbering conventions](@ref fvmadapt_numbering) define the local positions
used by those operations.


## What Is Stored

The Loci rule layer owns the mesh-wide plan stores. A plan is a variable-length
vector of split codes associated with one original mesh entity. The C++ library
does not keep a persistent `Cell`, `HexCell`, or `Prism` object for every plan.
Instead, it receives one plan vector at a time and reconstructs the temporary
objects needed for that operation.

| Plan | Associated with | What replay reconstructs |
| --- | --- | --- |
| Cell plan | One original cell | Fine cells inside the original cell |
| Face plan | One original face | Fine faces covering the original face |
| Edge plan | One original edge | Fine edge segments covering the original edge |

Cell plans describe volume subdivision. Face and edge plans describe the
shared boundaries on which independently reconstructed cells must agree. An
original face or edge provides the common reference even when neighboring
cells describe it with different local numbering or orientation.

Faces and edges created inside a split cell are consequences of the cell plan.
They do not each receive another independently stored mesh-wide plan.

```text
plan vector in a Loci store
        ↓ replay for one original entity
temporary parent-and-child objects
        ↓ inspect, balance, or construct
updated plan or fine mesh entities
        ↓
discard the temporary objects
```

The vectors therefore preserve the planned topology between rule operations
without retaining the more memory-intensive object trees.


## How One Plan Reconstructs a Tree

A plan records split codes in breadth-first order. Replay begins with the
original entity as the root of a queue. Each stored code applies to the first
object in the queue. If that code splits the object, its children are appended
to the queue.

Code `0` leaves the current object unsplit. When the stored vector ends, every
object still in the queue is also treated as unsplit. A compact plan can
therefore omit trailing zeroes.

The HexCell plan `[7, 1]`, introduced in the overview, produces the following
queue operations:

| Queue before the entry | Code | Operation | Queue after the entry |
| --- | --- | --- | --- |
| root | `7` | Split the root into children 0 through 7 | 0, 1, 2, 3, 4, 5, 6, 7 |
| 0, 1, 2, 3, 4, 5, 6, 7 | `1` | Split child 0 into children 0.0 and 0.1 | 1, 2, 3, 4, 5, 6, 7, 0.0, 0.1 |

The plan then ends, so the nine objects remaining in the queue are leaves.
Those leaves represent the fine cells that can later be constructed inside the
original cell.

The values `7` and `1` select HexCell split operations. The values used for
child positions, such as `0` and `0.1`, are not mesh cell identifiers. Other
cell, face, and edge types have their own split codes and child-numbering
conventions.

### Refinement Lineage and Mesh Incidence

The parent-and-child relationship reconstructed by a plan connects objects of
the same kind:

```text
cell → child cells
face → child faces
edge → child edges
```

That refinement lineage forms a tree. It is different from mesh incidence,
which connects unlike entities:

```text
cell → boundary faces → boundary edges → endpoint nodes
```

A split cell is an internal object in its cell tree. A cell that is not split
further is a cell-tree leaf. The same terminology applies separately to face
and edge trees.

### Local Positions Are Not Mesh IDs

Local positions identify roles inside one cell, face, edge, or child array.
They do not provide persistent global identity. For example, two adjacent
cells can refer to their shared face using different cell-local face
positions. Each cell can also number and orient the children of that face
differently.

Before FVMAdapt combines the two descriptions, it maps both into the original
face's orientation. The mapped child positions then refer to the same parts of
the shared face. Final mesh numbers are assigned later, when the fine
connectivity is constructed; they are not encoded in the plan.


## Deriving Shared Face and Edge Plans

A cell plan determines how each boundary face appears from that cell's local
view. FVMAdapt extracts that subdivision as a face plan. For an interior face,
it extracts one contribution from each adjacent cell.

The contributions are first mapped into the shared face's orientation. They
are then merged so that the shared-face plan contains every subdivision
required by either cell:

```text
left cell plan  --extract and orient--┐
                                      ├──> one shared-face plan
right cell plan --extract and orient--┘
```

### Combining Different Face Subdivisions

The merge operates on the subdivisions after orientation mapping, not on the
raw cell or face split codes. For the supported plan forms, it produces a
common refinement: the shared plan retains every distinct cut required by
either side.

| Left requirement | Right requirement after orientation | Shared result |
| --- | --- | --- |
| No split | Direction A | Direction A |
| Direction A | Direction A | Direction A |
| Direction A | Direction B | Directions A and B |
| Coarse subdivision | Same subdivision with a deeper child split | Coarse and deeper subdivisions |

For example, suppose a HexCell and a prism share a quadrilateral face. If the
hex contribution cuts that face in one face-local direction and the prism
contribution cuts it in the perpendicular direction, the merged root code is
QuadFace code `3`. Replaying that code creates four child faces. If the two raw
directional codes instead map to the same cut in the shared face orientation,
the duplicate cut is retained only once.

A quadrilateral face contributed by the general `Cell` path does not introduce
an unrelated subdivision pattern. Its isotropic general-face split is
converted to QuadFace code `3` before it is merged with a HexCell or prism
contribution. In the face-plan topology, that full split already contains both
supported directional subdivisions. FVMAdapt later reconstructs the shared
face once from the merged plan, so the adjacent cells do not create competing
face-center nodes or incompatible cuts independently.

For a general polygonal face, the merge is the corresponding tree operation:
if either input splits a tree position, the shared plan splits that position
and retains the descendant subdivisions from both inputs.

A boundary face has only its cell-side contribution. An interior face has
contributions from both sides, even when the adjacent cells use different cell
types.

The same process continues from faces to edges. A face plan implies
subdivisions of its boundary edges. Contributions from all incident faces are
mapped into the shared edge direction and merged:

```text
incident face plans
        ↓ extract edge subdivisions
map each contribution into the shared edge direction
        ↓ merge
one shared-edge plan
```

The result is one shared topological subdivision for each original face and
edge. Neighboring local numbering does not need to match. The orientation
mappings only need to make the local child positions refer to the same parts of
the shared entity.


## What Balancing Enforces

Balancing has two related jobs: shared-topology compatibility and refinement
grading.

### Shared-Topology Compatibility

Every cell touching a shared face must be compatible with the same face plan.
Every face touching a shared edge must be compatible with the same edge plan.
Without this agreement, replaying the plans independently could produce
different connectivity on opposite sides of the same mesh entity.

Compatibility does not require adjacent cells to have identical cell trees.
For a cell handled by the general `Cell` path, one cell-tree leaf can sometimes
be reconstructed with several fine faces on its boundary. The resulting fine
cell has a more detailed polyhedral boundary, but the cell does not need to
split solely to retain a standard element shape.

### Refinement Grading

Compatibility alone would allow an unsplit cell to border a boundary tree of
arbitrary depth. The grading checks limit that difference.

For example, suppose a boundary edge supplied to a cell-tree leaf is split,
and one of those child edges is split again. Relative to the leaf, the boundary
edge tree is two levels deep. The default edge-depth check requires that leaf
cell to split. The cell implementation selects an appropriate split operation
for its topology and, for specialized cells, the permitted split directions.

The optional `balance_option` policies can impose additional constraints.
They include splitting when more than half of a cell's faces are already
subdivided and, in cell paths that implement the test, splitting for certain
subdivided opposing-face combinations. These are separate topology-regularity
policies. Their precise split decisions depend on the cell type and split mode;
they are not one uniform rule applied identically to every cell representation.

A cell classified as indivisible is not expanded by the ordinary local
balancing operation. Its original boundary is still reconstructed using the
shared face and edge plans so that the resulting connectivity is compatible.


## The Iterative Balancing Loop

Shared face and edge plans are not merely outputs of cell balancing. They also
supply the boundary trees used to reconsider the cell plans:

```text
cell plans at iteration n
        ↓
extract, orient, and merge shared-face plans
        ↓
extract, orient, and merge shared-edge plans
        ↓
reconstruct each original cell with those boundary trees
replay and locally balance its current cell plan
        ↓
cell plans at iteration n + 1
        ↓
are all cell plans unchanged?
        ├── no: repeat
        └── yes: preserve the balanced plans
```

The implementation reaches stability at two scales. Within one reconstructed
cell, the library repeatedly checks its cell-tree leaves against the supplied
face and edge trees. Local balancing stops when another check adds no cell
split.

Across the mesh, the Loci rule layer regenerates the shared-face and
shared-edge requirements, reconstructs and balances the cells, and compares
each resulting cell plan with its plan from the preceding iteration. Mesh-wide
balancing stops only when every cell plan is unchanged.

The unchanged plans form a stable plan set: another balancing iteration would
produce the same cell plans. Final compatible face and edge plans can then be
derived from those balanced cell plans.


## How Refinement Propagates

Consider the cells `A`, `B`, and `C`, where `A` shares a face with `B` and `B`
shares another face with `C`.

As described in the overview, a requested split in `A` can subdivide the
`A`-`B` face. Cell `B` is then reconstructed with that shared-face subdivision.
If `B` is handled by the general `Cell` path, it may remain one cell with a
more detailed polyhedral boundary. If a grading or other active balancing
policy requires a cell split, however, the plan for `B` is expanded.

That change can affect the other side of `B`:

```text
A's plan subdivides the A-B face
        ↓
B is reconstructed and reconsidered
        ↓
B's plan expands
        ↓
the B-C face gains a subdivision
        ↓
C is reconsidered during a later mesh-wide iteration
```

General polyhedral cells can limit propagation caused only by the need to
preserve a standard element shape. They do not stop propagation required by
shared-topology compatibility, refinement grading, or another enabled
balancing policy.


## Balanced Plans and Fine-Mesh Construction

A balanced result provides:

- a stable plan for each original cell
- one compatible subdivision for every original shared face and edge
- cell plans that satisfy the applicable grading checks
- enough topological information to reconstruct the fine mesh.

It does not require adjacent cells to have identical cell trees, local
numbering, or orientation. It also does not imply uniform refinement,
coordinate smoothing, or general mesh-quality optimization.

After balancing, FVMAdapt replays the compatible plans to construct the fine
mesh:

```text
balanced cell, face, and edge plans
        ↓ replay each plan
temporary refinement trees and their leaves
        ↓ create and number the fine nodes, faces, and cells
numbered fine nodes, faces, and cells
        ↓ transfer field data and write the adapted mesh
adapted mesh with transferred field data
```

Coordinates and final mesh numbering enter during this construction stage.
Later rules handle geometric placement, connectivity numbering, field
transfer, and mesh output.
