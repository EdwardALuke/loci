# FVMAdapt Refinement Plans {#fvmadapt_refinement_plans}

This page documents the compact plan encodings used by the existing
`FVMAdapt` implementation. It is deliberately source-backed: each numeric
claim names the current header, function, or table that defines the behavior.

## What Is A Refinement Plan? {#fvmadapt_what_is_a_refinement_plan}

A refinement plan is a compact recipe for reconstructing the refinement tree
for one mesh entity. It is not the refined mesh itself: it does not store child
coordinates, fine-cell connectivity, or object pointers. Instead, it stores the
sequence of split decisions that the C++ tree builders replay when they need to
recreate the same edge, face, or cell hierarchy.

The tree starts at the original entity:

- The root is the original cell, face, or edge.
- A non-zero split code says how that tree node is refined.
- A `0` split code says that tree node is a leaf.
- The children created by a split are queued and read in breadth-first order.

This is why the implementation can store most plans as `std::vector<char>`.
The `char` values are compact numeric byte codes, not printable letters. For
example, byte value `1` means "split" for an edge or general face; it is not
the printable ASCII digit 1.

For an edge, the plan can be read as a tiny binary-tree recipe:

| Plan | Meaning |
| ---- | ------- |
| `[]` | Leave the original edge unchanged. |
| `[1]` | Split the original edge once; both child edges are leaves. |
| `[1, 1]` | Split the original edge, then split the first child edge. The missing trailing entry for the second child is treated like `0`. |
| `[1, 0, 1]` | Split the original edge, leave the first child edge alone, and split the second child edge. |

Faces and cells use the same broad idea, but their non-zero codes carry
entity-specific meaning. For example, `Face::make_faceplan()` writes only `0`
and `1`, `QuadFace` uses directional face split codes `1`, `2`, and `3`, and
`HexCell` uses a three-bit local-direction mask from `1` through `7`.

The word "plan" is therefore literal: the vector records what split should
happen at each tree node if that entity is replayed later. Loci rules can pass,
balance, merge, extract, write, and reload those recipes without keeping the
full in-memory tree alive between stages.

## Plan Storage And Traversal

Most adaptation plans are `std::vector<char>` values. The vector is consumed as
a breadth-first traversal of a refinement tree: each entry describes the split
at the next tree node. Several tree walkers treat missing trailing entries as
`0`, and many plan writers trim trailing `0` entries before returning.

Important examples:

- `Edge::resplit()` in `node_edge.h` walks an edge tree and interprets the
  edge plan.
- `Face::resplit()` and `Face::make_faceplan()` in `face.h` walk polygonal
  face trees.
- `QuadFace::resplit()` in `quadface.h` walks directional quadrilateral face
  trees.
- `HexCell::empty_resplit()`, `HexCell::resplit()`, and
  `HexCell::make_cellplan()` in `hexcell.h` walk hexahedral cell trees.
- `Prism::empty_resplit()`, `Prism::resplit()`, and
  `Prism::make_cellplan()` in `prism.h` walk prism cell trees.
- `DiamondCell` and `Cell` in `diamondcell.h` handle general-cell
  decomposition and traversal.

An empty vector usually means that the entity is unchanged by the current plan.
Do not assume that it is identical to every internal no-split case: callers
also accept explicit `0` entries inside non-empty plans, and extraction tables
use propagation code `8` for a different purpose.

One serialization convention uses a printable marker: `cellplan_io.loci` writes
an empty `balancedCellPlan` as a one-entry vector containing `C`. The readers
`get_cellPlan`, `get_DBcellPlan`, `get_balanced_cellPlan`,
`get_balanced_cellPlanDB`, `get_parentPlan`, and `get_parentPlanDB` normalize
that sentinel back to an empty vector. Treat this as an I/O marker for "no cell
plan", not as a normal split code in the refinement tree.

## Common Fact Names

`cellPlan`
: Cell-tree plan. The code meaning depends on the cell type.

`facePlan`
: Face-tree plan. Face plans may be extracted from cell plans, merged from
  neighboring cells, or converted between general-face and quad-face
  conventions.

`edgePlan`
: Edge-tree plan. Edge plans are commonly extracted from face plans.

`balancedCellPlan`, `balancedFacePlan`, `balancedEdgePlan`
: Post-balancing facts used by fine-grid generation and numbering rules after
  neighboring cell, face, and edge plans have been made compatible.

`balancedCellPlan1`, `balancedFacePlan1`, `balancedEdgePlan1`
: Related derefinement/restart-side facts used while building
  `priority::restart::balancedCellPlan`. The suffix `1` is historical; check
  the rule path before treating it as a time level.

`balancedPointSet1`
: Intermediate edge split-coordinate set used inside the general derefinement
  edge path. It is built from `balancedFacePlan1` and serialized into
  `balancedEdgePlan1`; multiple faces can contribute to the same edge, so the
  apply rule merges contributions with `Loci::SetUnion`.

`nodeTag`, `cellNodeTag`, `fineCellTag`
: Tag vectors used by resplit/derefinement paths. These tags are not
  interchangeable with plan vectors; check the specific builder before changing
  semantics.

## Edge Plans

Defined by `Edge` in `node_edge.h` and implemented in `library/node_edge.cc`.

| Code | Meaning |
| ---- | ------- |
| `0`  | Leaf edge, no split at this tree node. |
| `1`  | Split the edge at its midpoint and enqueue the two children. |

`Edge::resplit(edgePlan, needReverse, node_list)` uses `needReverse` to control
child traversal order when the edge is being replayed for an owning face.

## General Face Plans

Defined by `Face` in `face.h` and implemented in `library/face.cc`.

| Code | Meaning |
| ---- | ------- |
| `0`  | Leaf face, no split at this tree node. |
| `1`  | Split the polygonal face. |

`Face::split()` first ensures the parent edges are split, then creates one
quadrilateral child face per parent edge by connecting edge midpoints through a
face-center node. This is isotropic polygon-face refinement and is distinct
from `QuadFace` directional split codes.

## Quad Face Plans

Defined by `QuadFace` in `quadface.h` and implemented in
`library/quadface.cc`.

| Code | Meaning |
| ---- | ------- |
| `0`  | Leaf quad face, no split at this tree node. |
| `1`  | Split in the face-local y direction; `childy` is populated. |
| `2`  | Split in the face-local x direction; `childx` is populated. |
| `3`  | Split in both face-local directions; `child`, `childx`, and `childy` participate in the child layout. |

The local quad layout is documented in `quadface.h`: edge `0` runs along the
bottom, edge `1` along the right side, edge `2` along the top, and edge `3`
along the left side. The child id diagram in that header is the source of the
current child ordering.

Orientation helpers in `quadface.h` convert between face-local and cell-local
frames:

- `orient_splitCode_f2c()`
- `orient_childID_f2c()`
- `orient_edgeID_f2c()`
- `orient_edgeID_c2f()`

The numeric meanings of `orientCode` values are intentionally not documented
here yet; they should be recorded only after each conversion table has been
audited.

## Hex Cell Plans

Defined by `HexCell` in `hexcell.h` and implemented in `library/hexcell.cc`.

`HexCell::mySplitCode` is a three-bit local-direction mask. The comment in
`hexcell.h` states that the bit order is `xyz`; code `0` means no split.
`HexCell::numChildren()` gives the current child-count interpretation:

| Code | Split directions | Children |
| ---- | ---------------- | -------- |
| `0`  | none             | `0` |
| `1`  | z                | `2` |
| `2`  | y                | `2` |
| `4`  | x                | `2` |
| `3`  | y and z          | `4` |
| `5`  | x and z          | `4` |
| `6`  | x and y          | `4` |
| `7`  | x, y, and z      | `8` |

This is the hex-specific interpretation of the `xyz` bit mask: bit value `4`
is local x, bit value `2` is local y, and bit value `1` is local z. Do not
apply this table to prism, quad-face, general-face, or edge plans.

Hex cell plans imply quad face plans through `extract_hex_face()` and the
`faceCodeTable`/`faceIDTable` entries in `tables.h`.

## Prism Cell Plans

Defined by `Prism` in `prism.h` and implemented in `library/prism.cc`.

`Prism::numChildren()` gives the current child-count interpretation:

| Code | Split mode | Children |
| ---- | ---------- | -------- |
| `0`  | none | `0` |
| `1`  | two-child split; `split()` calls `QuadFace::split(1)` on each side face and creates one new end face | `2` |
| `2`  | ring split; `split()` calls `QuadFace::split(2)` on each side face and splits the two end faces | `nfold` |
| `3`  | combined split; `split()` calls `QuadFace::split(3)`, splits the two end faces, and creates center topology | `2 * nfold` |

In this implementation, `nfold` is used as the count for the side-face array
and for loops over the vertices/edges of the two end faces. The default/root
construction uses `nfold == 3`. Split paths that create children around the
side-face ring construct those child `Prism` objects with `nfold == 4`.

This is not the same formalism as `HexCell::mySplitCode`. Do not describe
prism codes as the hex `xyz` bit mask. The geometric direction names and the
best terminology for `nfold == 4` children remain documentation questions until
the prism split, extraction, and merge paths are audited together.

Prism face handling is split between triangular end faces (`Face`) and
quadrilateral side faces (`QuadFace`). `get_c1_prism()` uses different
strategies for these face families, and `extract_prism_face()`/merge helpers
in the prism face rules should be treated as part of the prism plan contract.

## General Cell Plans

General cells are represented through `Cell` and `DiamondCell` in
`diamondcell.h`, with helper implementations in `library/diamondcell.cc` and
`library/build_general_cell.cc`.

General-cell refinement decomposes a cell into diamond cells. A `DiamondCell`
with fold `n` has `2n + 2` nodes, `2n` faces, and `4n` edges according to the
header comment in `diamondcell.h`. This path uses isotropic polygon-face
splitting through `Face`.

## Propagation Code 8

`tables.h` defines a special code `8` in both `faceCodeTable` and
`edgeCodeTable`.

For cell-to-face extraction, `8` means that the parent cell is split but the
selected lower-dimensional face is not split at that step. It is not the same
as a plain no-split `0`, because the table comment says the face can still be
split in a later step and that tree construction treats `8` as a one-child
propagation case.

For face-to-edge extraction, the same idea applies: a parent face is split but
the selected edge is not split at that step, and the edge may still split later.

Do not replace `8` with `0` or trim it as if it were a trailing no-op without
auditing every extraction and merge path that consumes these tables.

## Documentation Questions

The following details should stay marked as open until they are verified
directly from the implementation:

- Complete numeric meaning of each `orientCode`.
- Exact child ordering for all mixed prism/general-cell cases.
- Best terminology for prism `nfold == 4` children: the code still represents
  them with `Prism` objects, but the intended mesh-element description should
  be verified with the split/extraction/merge paths before documenting it.
- Which intermediate plan facts should be considered stable restart-facing
  facts rather than temporary Loci facts. Current restart builders publish
  `priority::restart::balancedCellPlan`.
- Which derefinement tag values are local implementation details versus stable
  file or fact semantics.
