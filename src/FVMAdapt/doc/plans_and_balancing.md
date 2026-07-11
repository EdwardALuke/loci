# FVMAdapt Refinement Plans and Balancing {#fvmadapt_plans_and_balancing}

FVMAdapt decides the topology of an adapted mesh before it constructs that
mesh. Temporary refinement trees describe possible subdivisions of individual
mesh entities. Compact plans record how to reconstruct those trees, and
balancing coordinates the plans so independently reconstructed parts of the
mesh agree.

This page develops those ideas in more detail. The
[mesh-adaptation overview](@ref fvmadapt_overview) gives the shorter workflow,
the [geometric splitting reference](@ref fvmadapt_geometric_splitting) shows
what each split creates, and the
[numbering conventions](@ref fvmadapt_numbering) describe the local positions
used to interpret the plans.


## Mesh Incidence and Refinement Lineage

Two different relationships are important in FVMAdapt.

The first is **mesh incidence**:

```text
cell --bounded by--> faces --bounded by--> edges --end at--> nodes
```

This is not a refinement tree. A face can be shared by two cells, an edge can
belong to several faces, and a node can belong to several edges. These are
connections between different kinds of mesh entities.

The second relationship is **refinement lineage**:

```text
original entity
├── child 0
├── child 1
│   ├── child 1.0
│   └── child 1.1
└── child 2
```

This is a tree. The original entity is the root. A split entity is an internal
tree object, and an entity that is not split further is a leaf. Its children
are entities of the same kind: a cell has child cells, a face has child faces,
and an edge has child edges.

In this discussion, a *tree object* means an entity occupying a position in a
refinement tree. It should not be confused with a mesh `Node`, which represents
a geometric vertex.


## A Forest of Related Trees

FVMAdapt does not construct one global refinement tree for the entire mesh. A
useful conceptual model is a forest of related trees rooted at original mesh
entities. The implementation reconstructs the trees temporarily when it needs
to interpret or modify their plans.

| Plan | Root entity | What its leaves represent |
| --- | --- | --- |
| Cell plan | One original cell | Fine cells inside that original cell |
| Face plan | One original face | Fine faces covering that original face |
| Edge plan | One original edge | Fine edge segments covering that original edge |

The cell trees describe volume subdivision. Face and edge trees describe the
shared boundaries through which independently reconstructed cell trees must
agree. The original face or edge provides the common point of reference even
when neighboring cells use different local orientations.

These mesh-wide face and edge plans govern boundaries rooted at original mesh
faces and edges. New interior faces and edges created by a cell split are
consequences of that cell tree; they do not each receive another independently
stored mesh-wide plan.

The trees used during replay are temporary working objects. They are created
while a rule examines or constructs part of the mesh and are discarded when
that work is complete. The plans, rather than those temporary objects, are the
compact topological state carried between operations.


## Plans Are Replay Recipes

A plan is a sequence of split codes stored in breadth-first order. It records
which tree objects split and which split operation each one uses. It does not
store coordinates, temporary object addresses, or global node, edge, face, and
cell identifiers.

Breadth-first order processes all queued objects at one tree level before
objects at the next level:

```text
start with the root in a queue
        ↓
read the next split code for the first queued object
        ↓
if it splits, append its children to the queue
        ↓
remove the processed object and continue
```

Code `0` means that the current object remains a leaf. When a compact plan ends,
the remaining queued objects are also treated as leaves. Trailing no-split
codes therefore do not need to be stored.

For example, consider the HexCell plan `[7, 1]`:

```text
plan position 0: the root consumes code 7
                 → create root children 0 through 7

plan position 1: root child 0 consumes code 1
                 → create two children beneath root child 0

end of plan:     root children 1 through 7 and both new grandchildren
                 remain leaves
```

The numbers `0` through `7` in this example are ordered child positions, not
mesh cell identifiers. Code `7` and code `1` have the HexCell meanings given in
the splitting reference. Other entity types have their own split codes and
child counts.

The same plan can be viewed in either direction:

```text
temporary tree --serialize--> compact plan
compact plan    --replay-----> temporary tree
```

After replay, the tree leaves can be walked to construct the fine entities.
The plan is therefore neither the temporary tree nor the final fine mesh.


## Local Positions and Identity

Every child uses the local numbering convention of its entity type. A HexCell
child, for example, has its own local nodes, edges, and faces. Those numbers
describe positions within that child; they are not inherited global identities
from the parent.

Existing corner objects can be reused by a child, and newly created midpoint
or center objects can be shared by siblings. The same geometric object may
therefore occupy different local positions in different cells or faces.
Orientation mappings identify the corresponding physical positions.

The important identity layers are:

```text
original cell, face, or edge
        → is the root addressed by a plan
local position
        → identifies a role inside one cell, face, edge, or child array
temporary object
        → connects one reconstructed topology during replay
final mesh number
        → identifies an entity in the generated mesh connectivity
```

Final mesh numbering is assigned later. It is not encoded in the refinement
plan and is not the identity used to balance the trees.


## How Cell, Face, and Edge Plans Meet

A cell plan determines how each boundary face appears from that cell's local
view. FVMAdapt extracts that boundary subdivision as a face plan. For an
interior face, it performs this operation for both adjacent cells.

Before the two views can be compared, each one is mapped into the shared face's
local orientation. The mapped plans are then merged so the shared face contains
every subdivision required by either cell:

```text
left cell plan  --extract and orient--┐
                                      ├──> one shared-face plan
right cell plan --extract and orient--┘
```

A boundary face has only the cell-side contribution. An interior face has both
cell-side contributions, even if the adjacent cells have different element
types.

The same idea continues from faces to edges. Each face plan implies split
positions on its boundary edges. Contributions from every incident face are
mapped into the edge direction and merged:

```text
incident face plans
        ↓ extract edge subdivisions
map every contribution into the shared edge direction
        ↓ merge
one shared-edge plan
```

The result is one physical subdivision of each original shared face and edge.
Neighboring local numbering does not need to be identical; the orientation
mappings only need to make the local views refer to the same physical
subentities.


## What Balancing Enforces

Balancing has two related jobs.

### Shared-Topology Compatibility

Every cell touching a shared face must be compatible with the same face
subdivision. Every face touching a shared edge must be compatible with the same
edge subdivision. Without this agreement, independently replayed plans could
describe different connectivity on the two sides of the same physical entity.

Compatibility does not require adjacent cells to have identical cell trees. A
cell may remain one leaf while its boundary is represented by several fine
faces. The resulting fine cell simply has a more detailed polyhedral boundary.
This is how general polyhedral cells can absorb some neighboring refinement
without forcing the cell itself to split.

### Refinement Grading

A splittable cell cannot remain arbitrarily coarse while the trees on its
boundary become arbitrarily deep. In the default grading check, a leaf cell is
split when one of its incident edge trees extends more than one additional
refinement level beneath that leaf's boundary edge.

This condition is more precise than saying that every pair of neighboring
cells must have refinement levels that differ by at most one. The
implementation examines the boundary edge trees of each leaf and handles the
required split according to the cell type and available split directions.

The library also contains stronger grading policies that can split a cell when
too many of its faces, or certain opposing faces, are already subdivided. These
are additional topology-regularity policies rather than the basic requirement
that both sides of a shared entity describe compatible connectivity.

A cell classified as indivisible is not expanded by the ordinary local
balancing step. It can still be reconstructed with the compatible subdivision
of its original boundary.


## The Balancing Feedback Loop

Balancing is not a single pass from cell plans to face plans to edge plans.
Temporary face and edge plans feed their requirements back into the cell plans:

```text
current cell plans at iteration n
        ↓
extract each cell's boundary-face views
map neighboring views into shared-face orientation and merge them
        ↓
temporary shared-face plans
        ↓
extract, orient, and merge their edge subdivisions
        ↓
temporary shared-edge plans
        ↓
reconstruct each original cell with those boundary trees
replay its current cell plan and locally balance its leaves
        ↓
cell plans at iteration n + 1
        ↓
are all cell plans unchanged?
        ├── no: repeat with the expanded plans
        └── yes: preserve the balanced cell plans
```

On the refinement path, the shared face and edge plans accumulate the
subdivisions required for compatibility. Cell plans preserve requested
refinement and expand when the grading checks require a cell split. A change in
one cell can therefore affect a neighbor on the next pass, and that neighbor
can affect another cell through a different face.

Once the cell plans stop changing, the stable plans form a fixed point: another
balancing pass would produce the same cell plans. Final compatible face and
edge plans are then derived from those balanced cell plans.


## Local and Mesh-Wide Convergence

The implementation reaches the fixed point at two scales.

Within one reconstructed splittable cell, local balancing repeatedly examines
its leaves against the currently supplied face and edge trees. It adds required
child splits until another local pass makes no change.

Across the mesh, the rule layer repeats the larger exchange. It regenerates
temporary shared-face and shared-edge requirements, reconstructs and balances
the affected cell trees, and compares the resulting cell plans with the plans
from the previous iteration. Mesh-wide balancing finishes only when every cell
plan is unchanged.

```text
local convergence:     leaves within one reconstructed cell stop changing
mesh-wide convergence: plans for all original cells stop changing
```

The two scales explain why a refinement request can propagate through several
neighbors over several balancing iterations.


## A Two-Cell Propagation Example

Consider two cells, `A` and `B`, sharing a face:

```text
requested split in cell A
        ↓
A's cell plan implies a finer subdivision of the shared face
        ↓
the shared-face and shared-edge plans carry that subdivision to cell B
        ↓
B is reconstructed with the finer boundary
        ↓
does B satisfy the grading rules without a cell split?
        ├── yes: B remains one cell with a subdivided polyhedral boundary
        └── no:  B's cell plan is expanded
```

If `B` is expanded, its other boundary faces can acquire new subdivisions.
Those changes may reach another neighbor `C` during the next mesh-wide
iteration:

```text
A requests refinement
        → shared boundary changes B
        → B may change another shared boundary
        → C is reconsidered on a later pass
```

General polyhedral cells limit this propagation because a cell can often
accept a subdivided boundary without being split solely to retain a standard
element shape. The shared boundary plans still carry compatibility
requirements, but cell-split propagation continues only where grading or
another active topology policy requires it.


## What a Balanced Result Is and Is Not

A balanced result is:

- a stable set of cell plans
- one compatible subdivision for every shared face and edge
- cell plans stable under the applicable grading checks
- sufficient topological information to reconstruct the adapted mesh

A balanced result does not require:

- identical cell trees on opposite sides of a face
- identical local numbering or orientation in neighboring cells
- uniform refinement throughout the mesh.

Here, *balancing* means topological compatibility and refinement grading. It is
not processor load balancing, coordinate smoothing, or general mesh-quality
optimization. It produces compatible plans, not the final mesh coordinates and
connectivity.


## From Balanced Plans to the Fine Mesh

After balancing, FVMAdapt replays the compatible plans to construct the fine
mesh:

```text
balanced cell, face, and edge plans
        ↓ replay
temporary trees and their leaves
        ↓ construct and number
fine nodes, faces, and cells
        ↓ transfer and output
adapted mesh and associated field data
```

Coordinates and final mesh numbering enter during this construction stage.
This separation allows the balancing process to reason primarily about
topology while later rules handle geometric placement, connectivity numbering,
field transfer, and mesh output.
