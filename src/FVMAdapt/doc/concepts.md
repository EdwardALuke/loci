# FVMAdapt Concepts and Terminology {#fvmadapt_concepts}

This page is a developer note for terms and concepts used inside the existing
`FVMAdapt` implementation. It is intentionally implementation-facing: prefer
small, source-backed statements here over broad claims about mesh adaptation
design.

## Scope

`FVMAdapt` is split between Loci rules and C++ helper code:

- `src/FVMAdapt/*.loci` contains Loci rules that create, balance, merge, and
  consume adaptation facts.
- `src/FVMAdapt/library/*.cc` contains C++ helper routines for tree operations,
  plan conversion, face and edge extraction, and VOG output.
- `src/include/FVMAdapt/*.h` declares the helper classes and low-level tables
  used by the module.
- `src/FVMAdapt/doc/README_xml` documents the XML region input accepted by one
  marker mode.

The notes below describe the current code vocabulary. They are not a proposed
new API.

Related pages:

- \ref fvmadapt_developer_reference
- \ref fvmadapt_refinement_plans
- \ref fvmadapt_diagram_requests
- \ref fvmadapt_testing
- \ref fvmadapt_refactor_ledger

## Loci Rule Styles In FVMAdapt

FVMAdapt uses two rule-writing styles in `.loci` files:

- hand-written C++ classes derived from rule base classes such as
  `pointwise_rule`, `unit_rule`, and `apply_rule`
- newer `$rule ...` syntax that the Loci preprocessor lowers into generated
  C++ rule classes

These two styles are closer than they first look. For a typical store-valued
pointwise rule, lpp generates the same broad shape as the hand-written class:

- `$type name store<T>;` supplies the store type that a hand-written class would
  otherwise declare with a member such as `store<T>` or `const_store<T>`.
- The `$rule` signature supplies the constructor-level `input(...)`,
  `output(...)`, and `constraint(...)` declarations.
- The ordinary `$rule` body becomes a generated `calculate(Entity _e_)` method.
- The generated `compute(const sequence &seq)` wrapper calls
  `do_loop(seq, this)`.
- lpp also emits the `register_rule<...>` call.

For example, the pointwise rule

```cpp
$type newCellPlan store<std::vector<char> >;
$type cellUnchanged store<bool>;
$type tmpCellPlan store<std::vector<char> >;

$rule pointwise(cellUnchanged{n=0}, tmpCellPlan{n=0} <- newCellPlan),
  constraint(gnrlcells) {
  $cellUnchanged{n=0} = ($newCellPlan.size() == 0);
  $tmpCellPlan{n=0} = $newCellPlan;
  reduce_vector($tmpCellPlan{n=0});
}
```

is lowered into a generated `pointwise_rule` class with `name_store(...)` calls
for `newCellPlan`, `cellUnchanged{n=0}`, and `tmpCellPlan{n=0}`; constructor
calls for `input("newCellPlan")`, `output(...)`, and `constraint("gnrlcells")`;
a generated `calculate(Entity)` containing the assignments; and a generated
`compute(seq)` that loops over the scheduled sequence.

The word `compute` is therefore overloaded:

- In hand-written C++ rules, `compute(const sequence &seq)` is the sequence-level
  virtual method called by the runtime. If it calls `do_loop(seq, this)`, the
  per-entity work is in `calculate(Entity)`.
- In `$rule pointwise`, `$rule unit`, and `$rule apply` syntax, the body after
  `compute { ... }` is normally transformed into `calculate(Entity)`, while lpp
  still emits the sequence-level `compute(seq)` wrapper.
- A `$rule prelude { ... }` block becomes `prelude(seq)`. It runs once for the
  scheduled sequence before the kernel loop, so it is the right place for
  output sizing and setup that should not be repeated per entity.
- Singleton, constraint, optional, default, and some parameter-output rules are
  different: their rule body is emitted as sequence-level `compute(seq)` rather
  than per-entity `calculate(Entity)`.

When translating a hand-written FVMAdapt class rule to `$rule` syntax, first
move any C++ member type information into `$type` declarations. Without those
types, lpp cannot resolve variables that were previously known only through the
class members.

### C++ Class To `$rule` Conversion Checklist

For small conversions, start locally and keep the translated rule beside the
old rule until lpp accepts it:

1. Check for existing `$type` declarations in shared `.lh` files before adding
   local declarations. Standard mesh and topology facts such as `lower`,
   `upper`, `boundary_map`, `face2node`, `face2edge`, `edge2node`, `pos`,
   `geom_cells`, `faces`, `interior_faces`, and `fileNumber(X)` come from
   `FVM.lh`; FVMAdapt-specific facts should live in `FVMAdapt/fvmadapt.lh`
   once the translation is stable.
2. If a conversion is still experimental, it is fine to start with nearby
   `$type` declarations. Move them into an `.lh` file before finishing if
   another rule file or external module needs the same fact types.
3. Use the exact fact names from `name_store(...)`, not necessarily the C++
   member names. Preserve time and priority labels such as `{n=0}`, `{n+1}`,
   and `priority::...`. For example,
   `name_store("tmpCellPlan{n+1}", tmpCellPlann1)` means the `$rule`
   signature and body should use `tmpCellPlan{n+1}`; occurrences of the C++
   member `tmpCellPlann1` in `compute(...)` or `calculate(...)` become
   `$tmpCellPlan{n+1}`.
4. Convert constructor declarations into the `$rule` signature:
   `input(...)` becomes the source side, `output(...)` becomes the target side,
   and `constraint(...)` stays a rule modifier.
5. For a pointwise class whose `compute(seq)` only calls `do_loop(seq, this)`,
   drop the wrapper and translate `calculate(Entity e)` into the `$rule` body.
   A guard such as `if (seq.size() != 0)` around only `do_loop(seq, this)` is
   still just a no-op wrapper around the loop.
6. Replace array-style member access with `$` variables. For example,
   `cellUnchanged[cc] = ...` becomes `$cellUnchanged = ...`, and
   `tmpCellPlan[cc]` becomes `$tmpCellPlan`.
   When the old body passed the current `Entity` itself to a helper, use `_e_`
   in the `$rule` body; lpp emits the body inside `calculate(Entity _e_)`.
   For pointwise parameter inputs, a hand-written dereference such as
   `*split_mode_par` becomes `$split_mode_par`, not `$*split_mode_par`.
7. Preserve comments from the C++ `compute(...)` or `calculate(...)` body during
   the first transcription when they still describe the algorithm. Rewrite or
   remove only comments that become stale under the `$rule` form.
8. Use `$*fact` when the old C++ body passed the whole store or map container to
   a helper function. Plain `$fact` means the value at the current rule entity.
9. For a multimap row such as `lower[cc]`, `$lower` is already the current row;
   use `$lower.begin()` and `$lower.size()` instead of the old
   `lower[cc].begin()` and `lower.num_elems(cc)` pair.
10. Translate chained map access one hop at a time. The `$rule` signature
   declares the traversal path, and the body uses `->$fact` to read the next
   fact at the entity reached by the previous hop. For a rule over faces,
   `face2node[f][i]` becomes `$face2node[i]`, and
   `pos[face2node[f][i]]` becomes `$face2node[i]->$pos`. A deeper C++ access
   such as `edge2node[face2edge[f][i]][0]` becomes
   `$face2edge[i]->$edge2node[0]`; adding the endpoint position gives
   `$face2edge[i]->$edge2node[0]->$pos`. Likewise, when a cell rule loops over
   face rows, `pos[face2node[lower[cc][fi]][j]]` becomes
   `$lower[fi]->$face2node[j]->$pos`.
11. If the class `compute(seq)` does real work before or after `do_loop`, that
   work may belong in `prelude { ... }`, `postlude { ... }`, or may mean the
   rule is not a simple pointwise conversion.
12. For `apply_rule<Container, Op>` classes, put the reduction operator in the
   `$rule apply(...)[OpTemplate]` brackets and keep the body as the per-entity
   `join(...)` call. lpp instantiates bracketed apply operators with the output
   value type. Prefer an existing templated Loci reducer when one matches, such
   as `Loci::LogicalAnd`, `Loci::LogicalOr`, `Loci::SetUnion`,
   `Loci::Maximum`, or `Loci::Summation`.
13. Be cautious with `disable_threading()`, direct `fact_db` access, file IO,
   MPI/distribution work, and module databases. Those rules often rely on
   sequence-level side effects and should not be mechanically rewritten.
14. Run lpp on the file before trusting the translation. A warning about a
   runtime-created constraint can be acceptable in existing FVMAdapt code, but
   an unknown fact type usually means a missing `$type`.

## Core Formalism

The existing implementation is easiest to read if refinement is treated as a
set of typed trees:

- The original cell, face, or edge is the root of the tree.
- A split creates children under the current tree node.
- A leaf is a tree node with no children in the current plan.
- A parent is the tree node that owns a child pointer array such as
  `Edge::child`, `Face::child`, `QuadFace::child`, `HexCell::childCell`, or
  `Prism::childCell`.
- A plan is the compact replay recipe for that tree.

Most tree objects store both topology and implementation state. For example,
`Edge` stores `head`, `tail`, child edges, and a refinement `level`; `Face`
stores boundary edges plus a child-face tree; `HexCell` and `Prism` store
faces, child cells, and cell-local split codes.

The word "sibling" needs care. In an ordinary tree discussion, siblings share
the same parent. Some FVMAdapt neighbor helpers use "sibling neighbor" more
loosely for same-size face neighbors, and comments in `HexCell` and `Prism`
explicitly say those neighbors do not necessarily have the same parent.

## Split, Resplit, And Empty Replay

The code uses a few related verbs:

`split`
: Mutate the current tree node by creating geometric children. A split usually
  creates new nodes, edges, faces, or cells and appends those allocations to
  cleanup lists supplied by the caller.

`empty_split`
: Create the child tree shape without creating the geometric node/edge/face
  data. This is used when code needs to replay plan structure, count leaves, or
  compare tree shape without rebuilding full geometry.

`resplit`
: Replay an existing plan vector into a tree. Resplitting is how a stored
  `edgePlan`, `facePlan`, or `cellPlan` becomes an in-memory helper tree again.

`make_*plan`
: Traverse an already-built tree and serialize its split decisions back into a
  compact plan vector.

`derefine`
: Low-level tree removal used by derefinement/coarsening paths. This word is
  used inside the existing FVMAdapt implementation; higher-level solver-facing
  documents may prefer "coarsen" when describing requested adaptation behavior.

## Directional Frames

Plan codes are local to the entity type that consumes them. They are not a
single global convention:

- `Edge` has only one local direction, from `head` to `tail`.
- `Face` handles general polygon faces and uses isotropic split code `1`.
- `QuadFace` uses a face-local two-dimensional frame. Code `1` is local y,
  code `2` is local x, and code `3` is both.
- `HexCell` uses a three-bit local `xyz` mask. Code `4` splits x, code `2`
  splits y, code `1` splits z, and combinations split multiple directions.
- `Prism` uses a two-bit prism mask: code `1` splits the axial/z direction,
  code `2` splits the polygonal xy cross-section, and code `3` splits both.
- `DiamondCell` and general-cell helpers do not use the hex `xyz` bit-mask
  formalism. They use isotropic diamond-cell decomposition and general polygon
  face splitting.

Orientation helpers bridge these local frames. For example, quad-face helpers
convert split codes, child ids, and edge ids between face-local ordering and
cell-local ordering. `transfer_f2c()` and `transfer_c2f()` do similar work for
integer face ranges.

## Centers And Centroids

The implementation uses names such as `simple_center()`, `wireframe()`, and
`centroid()` for the points used to split faces and cells. These names should
not be read as a general guarantee of an exact mathematical centroid.

`defines.h` currently sets `CENTROID = 1`, which selects the wireframe-style
center in several helper classes. In that mode, face and cell centers are often
computed with `weighted_center()` over edge centers or face centers. The
unweighted helper is `point_center()`, whose source comment calls it the
"unweighted mass center" of a set of points.

For documentation and diagrams, use "split center" or "computed center node"
when the exact formula does not matter. Use "center of mass" only when the
text is referring to `point_center()` or `weighted_center()` behavior directly,
and reserve "centroid" for references to the actual `centroid()` helper
functions.

## What Is A Plan?

A plan is a compact recipe for reconstructing the refinement tree of one mesh
entity. It is usually stored as a `std::vector<char>`, where each `char` is a
small numeric split code for the next tree node visited in breadth-first order.
The entries are byte values, not readable text: code `1` is the numeric value
one, not the printable ASCII digit 1.

The root of the tree is the original edge, face, or cell. A non-zero split code
creates children and queues them for later entries in the vector. Code `0`
means the current tree node is a leaf. The replayed tree is then used to build
fine nodes, faces, cells, numbering, and transfer information.

Several helpers treat missing trailing entries as no-split entries. For
example, `Edge::resplit()`, `Face::make_faceplan()`,
`QuadFace::empty_resplit()`, and the `get_c1_*()` helpers all use breadth-first
tree walks where running out of plan entries means code `0` for the current
node. Many plan-writing routines remove trailing `0` values before returning.

An empty plan vector normally means the entity is unchanged at that level. This
convention is common, but callers should still check the specific helper they
are using because some routines also accept explicit `0` entries inside a
non-empty plan.

The detailed split-code tables and examples are in
\ref fvmadapt_refinement_plans.

## Common Plan Names

`cellPlan`
: Plan for a cell tree. Hex cells, prisms, and general cells each interpret
  their cell split codes through their own class methods.

`facePlan`
: Plan for a face tree. A face plan may be extracted from a cell plan, merged
  from neighboring cells, converted between general-face and quad-face
  conventions, or replayed to create fine faces.

`edgePlan`
: Plan for an edge tree. Edge plans are often extracted from face plans so that
  face boundaries remain consistent.

`balancedCellPlan`, `balancedFacePlan`, `balancedEdgePlan`
: Loci fact names used after balancing rules have made neighboring cell, face,
  and edge plans compatible.

`balancedCellPlan1`, `balancedFacePlan1`, `balancedEdgePlan1`
: Related fact names used in derefinement and restart reprocessing paths. These
  are not a separate split-code formalism from the unsuffixed balanced plans.

`nodeTag`, `cellNodeTag`, `fineCellTag`
: Tag vectors used by selected resplit and derefinement paths. These are not
  general substitutes for plans; inspect the specific builder or Loci rule
  before changing their meaning.

## Split Codes

The current implementation stores split decisions as small `char` codes. Avoid
renaming or reinterpreting them without checking every table and tree traversal
that consumes the code.

### Edge Plans

An edge is a binary tree:

- `0`: no split at this tree node.
- `1`: split the edge at its midpoint and enqueue the two children.

`Edge::resplit(edgePlan, needReverse, node_list)` uses `needReverse` to queue
children in the order needed by the owning face.

### General Face Plans

`Face` represents a polygonal face. `Face::make_faceplan()` writes:

- `0`: leaf face.
- `1`: split the face.

`Face::split()` creates one quadrilateral child face for each parent edge by
connecting edge midpoints to a face-center node. This is different from the
directional split codes used by `QuadFace`.

### Quad Face Plans

`QuadFace` supports directional split codes:

- `0`: no split.
- `1`: split in the face-local y direction, stored through `childy`.
- `2`: split in the face-local x direction, stored through `childx`.
- `3`: split in both directions, stored through `child`.

The comments in `quadface.h` define the local edge and child layout. Helper
functions such as `orient_splitCode_f2c()`, `orient_childID_f2c()`,
`orient_edgeID_c2f()`, and `orient_edgeID_f2c()` convert between face-local and
cell-local orderings.

### Hex Cell Plans

`HexCell::mySplitCode` is a three-bit local-direction mask as documented in
`hexcell.h`: `0` means no split, and a value such as `3` (`011`) means the y
and z coordinates are split. The number of children is derived from that code:
one split direction creates 2 children, two split directions create 4 children,
and three split directions create 8 children.

Face plans for hex cells can be extracted from cell plans with
`extract_hex_face()`. The extraction uses `faceCodeTable` and `faceIDTable` in
`tables.h`.

### Prism Cell Plans

`Prism::mySplitCode` is interpreted by `Prism` methods rather than by the hex
tables. The current `Prism::numChildren()` implementation uses:

- `0`: no split.
- `1`: 2 children.
- `2`: `nfold` children.
- `3`: `2 * nfold` children.

Prism faces need separate handling for triangular end faces and quadrilateral
side faces. In current code, `get_c1_prism()` handles side faces with 2D
integer ranges when `faceID >= 2`, and triangular faces with `Face` tree
traversal when `faceID < 2`.

### Special Propagation Code `8`

`tables.h` uses code `8` when a parent cell or face is split but the selected
lower-dimensional entity is not split at that step. It acts as a propagation
case during extraction and merging. Do not treat it as the same thing as a
plain `0`: comments in `tables.h` note that `8` can allow the face or edge to
be split in a later step.

## Orientation

Many helpers distinguish between a face as ordered by global `face2node` data
and the same face as seen from a particular cell. The `orientCode` values encode
that relationship.

Known, source-backed facts:

- `transfer_f2c()` maps integer face-local coordinates into a cell-local face
  frame.
- `transfer_c2f()` maps in the opposite direction.
- `orient_*_f2c()` helpers convert split, child, or edge ids from face-local to
  cell-local order.
- `orient_edgeID_c2f()` converts quad edge ids from cell-local to face-local
  order.

The numeric meanings of all orientation codes should be verified in the
specific helper before documenting them as a public contract.

## Integer Face Ranges

Several routines use `Range2d` over `Point2d` integer coordinates to compare
fine faces without building the full geometric tree. The coordinate box often
uses:

```text
maxX = int64(1) << MAXLEVEL
maxY = int64(1) << MAXLEVEL
```

`MAXLEVEL` is defined in `defines.h`. Comments there warn that deeper
refinement would exceed the integer coordinate limit.

`contain_2d()` is a common helper for this style: it maps each fine face range
to the containing cell-side range and returns the associated cell id.

## Cell, Face, And Edge Trees

The C++ helpers model refinement with tree-like objects:

- `Node` stores a point, an output/input index, and a tag.
- `Edge` is a binary tree with `head`, `tail`, optional children, and a
  refinement level.
- `Face` is a polygonal face tree. Splitting creates one quadrilateral child per
  parent edge.
- `QuadFace` is a quadrilateral face tree with directional child storage.
- `HexCell`, `Prism`, and `DiamondCell` model cell refinement trees.
- `Cell` in `diamondcell.h` handles general cells through diamond-cell
  decomposition.

Ownership is manual. Destructors often delete child arrays, while many builder
paths also collect temporary `Node`, `Edge`, `Face`, or `QuadFace` objects in
lists that are later passed to `cleanup_list()`. When changing code, identify
which object owns each pointer before moving allocation or deletion.

## Balancing And Consistency

In FVMAdapt, balancing means expanding a requested plan until adjacent cells,
faces, and edges can be replayed consistently. Balancing does not directly
build the fine mesh; it produces compatible plan facts that later construction,
numbering, transfer, and restart rules consume.

Plan consistency is maintained across dimensions:

- Cell plans imply face plans on each cell face.
- Face plans from neighboring cells are merged into a shared face plan.
- Edge plans are extracted from face plans.
- Later rules generate fine nodes, fine faces, cell numbering, and transfer
  data from the balanced plans.

Plan derivation is the lower-dimensional part of that process. Cell plans are
used to derive face plans, and face plans are used to derive edge plans or
intermediate edge split-coordinate sets. Existing source names usually say
`extract`, `merge`, or `make_*plan`; use "derive" or "project" only as
explanatory prose, not as a separate API term.

The difficult part is not only whether an entity is split, but whether the
ordering of generated children, faces, and nodes remains consistent between
neighboring cells and across restart or derefinement paths.

`balancedCellPlan`, `balancedFacePlan`, and `balancedEdgePlan` are the ordinary
post-balancing plan facts used by fine-grid generation and numbering rules.
`balancedCellPlan1`, `balancedFacePlan1`, and `balancedEdgePlan1` are the
corresponding derefinement/restart-side plan facts used while building
`priority::restart::balancedCellPlan`. The suffix `1` is historical; do not
infer a time level from the name without checking the rule path.

General derefinement edge balancing may use an intermediate set of edge split
coordinates rather than an edge plan vector directly. For example,
`balancedPointSet1` accumulates split coordinates contributed by faces and is
then serialized into `balancedEdgePlan1`.

## Documentation Practice

When adding to this page:

- Prefer source-backed descriptions and include the function, header, or table
  name that proves the statement.
- Mark uncertain behavior as an open question instead of guessing.
- Keep numeric code tables tied to the current implementation.
- Avoid documenting a desired future AMR API here; this page describes the
  existing `FVMAdapt` module.
