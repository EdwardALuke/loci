# History and Balancing {#fvmadapt_plans_and_balancing}

FVMAdapt records the topology it intends to construct before it creates the
fine mesh. Here, *balancing* means making the subdivisions of shared faces and
edges compatible and enforcing refinement-grading rules. It does not mean
processor load balancing.

## What Is Stored

Each original cell has a refinement history. It records the parent-and-child
relationships, the accepted splits and their locations, and which descendants
form the current mesh. The saved geometry places those descendants in physical
space. The history describes the complete refinement below the original cell,
not merely the changes requested in the latest cycle.

A plan is a stored description from which a subdivision can be reconstructed;
it need not contain the complete output mesh.

Cell histories describe volume subdivision. Face and edge divisions describe
the shared boundaries on which independently reconstructed cells must agree.
An original face or edge provides the common reference even when neighboring
cells describe it with different local numbering or orientation.

Faces and edges created inside a split cell are consequences of the cell's
subdivision. They do not each receive another independently stored mesh-wide plan.

## A Cell History

Consider a hex split into eight children, with child 0 then split into two:

```text
root
├─ 0
│  ├─ 0.0  leaf
│  └─ 0.1  leaf
├─ 1       leaf
├─ 2       leaf
├─ ...
└─ 7       leaf
```

There are nine leaves: the two children beneath child 0 and the seven unsplit
root children. These are the cells used to construct the current mesh.
The labels describe positions in the tree, not output mesh numbers.

Coarsening children 0.0 and 0.1 makes their parent, child 0, a leaf again.
Their accepted construction remains saved. Repeating the same split reuses
its points rather than choosing new centers. Loading saved history likewise
restores accepted choices.

## History and Mesh Correspondence {#fvmadapt_mesh_correspondence}

Refinement history describes the descendants of each original cell and persists
between adaptations.

Mesh correspondence describes one adaptation: how cells and faces in the
preceding mesh relate to those in the resulting mesh. A refined child comes
from one preceding cell; a coarsened parent receives contributions from several.
Face correspondence also records orientation and area fractions.

## Deriving Shared Face and Edge Divisions

Each cell's subdivision determines how its boundary faces appear from that
cell's local view. For an interior face, there is one contribution from each
adjacent cell. The contributions are first mapped into the shared face's
orientation, then combined:

```text
left cell  --extract and orient-┐
                                ├──> one shared-face subdivision
right cell --extract and orient-┘
```

Edges shared by several faces also need one consistent set of divisions.
Both sides must use the same accepted split points. Incompatible choices are
rejected, rather than moving an already accepted point to make them agree.
See [numbering conventions](@ref fvmadapt_numbering) for orientation.

Compatibility does not require adjacent cells to have identical refinement
trees. A coarse cell can retain several fine faces on its boundary without
being split. Grading rules determine when further cell refinement is needed.

## What Balancing Enforces

Refinement grading limits the difference in cell split depth across shared
faces, edges or vertices. The default allows a difference of one level across
vertices. Depth counts cell splits, not physical size ratios or separate
directional depths.

Optional policies also refine a cell when more than half its faces are divided,
or when two divided faces share no vertex. The planar path counts edges.
Only divisions used by current neighbors count; unused saved splits do not
prevent coarsening.

The usual response is an isotropic split, preserving a requested extrusion
axis. Directional balancing can instead split hexes and prisms along the
required local directions. Missing direction information uses a full split.

## The Iterative Balancing Loop

```text
apply requested refinements and coarsenings to a candidate history
        ↓
combine shared face and edge requirements
        ↓
check neighbor depths and enabled balancing policies
        ↓
add required cell splits or restore a refused coarsening
        ↓
has any cell changed?
        ├── yes: repeat the shared-boundary and grading checks
        └── no: accept if enabled physical-mesh quality limits also pass
```

A split in cell A can require its neighbor B to split. That split can then
affect B's other neighbor C, which is why more than one pass may be needed.
The kernel reconstructs and edits individual cell histories; the Loci rules
combine shared-boundary requirements and coordinate the mesh-wide iteration.
A rejected candidate leaves the accepted histories intact.

A balanced result provides compatible shared boundaries and the requested
refinement grading. It does not imply uniform refinement, coordinate smoothing,
or general mesh-quality optimization.
