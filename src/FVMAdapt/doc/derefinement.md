# FVMAdapt Derefinement {#fvmadapt_derefinement}

Derefinement, also called *coarsening*, removes existing refinement. In
FVMAdapt, this does not mean merging an arbitrary collection of neighboring
fine cells. It means replacing one complete family of leaf cells with their
immediate parent in the refinement tree.

```text
before                                      after

parent                                      parent
├─ child 0  leaf, requests derefinement     └─ leaf
├─ child 1  leaf, requests derefinement
└─ ... all remaining children
```

The parent in this diagram is not another live mesh cell stored alongside its
children. FVMAdapt reconstructs the parent-and-child objects by replaying the
cell plan. If the family can be collapsed, it removes the children from that
temporary tree and serializes the shorter tree as a new plan. The adapted mesh
is then constructed from the leaves of the new plan.


## Derefinement Requests

The tag data used by the refinement and derefinement paths has three states:

| Value | Requested action |
| --- | --- |
| `0` | Leave the current entity unchanged |
| `1` | Refine |
| `2` | Derefine |

A value of `2` is a request, not permission to remove a cell. Derefinement is
accepted only for an eligible sibling family.

FVMAdapt can associate these states with the current leaves in two ways. A
fine-cell tag applies directly to one current fine cell. With node-based tags,
a leaf cell requests derefinement only when all of its nodes have value `2`.
If any node has value `1`, the cell requests refinement; other mixtures leave
the cell unchanged. Requiring all nodes for a node-based derefinement request
keeps a partially tagged cell from being interpreted as a request to remove
the whole cell.


## The Sibling-Family Test

Suppose a split tree object `P` has direct children `C0`, `C1`, and so on.
FVMAdapt can collapse that split only when all of the following are true:

1. `P` is currently split.
2. Every direct child of `P` is a leaf. None of the children can still have
   children of its own.
3. Every direct child requests derefinement.
4. None of the edges bounding `P` has refinement extending more than one level
   below `P`.

The last condition preserves the edge-depth grading assumed by the
derefinement code. The check is made against reconstructed boundary trees, so
it accounts for the subdivisions required on shared faces and edges. If any
condition fails, the family remains refined during that adaptation cycle.

This is why the useful unit of derefinement is a *sibling family*, not an
individual leaf. Removing only one child would leave a hole: the other
children and the parent overlap rather than partitioning the same volume.
Replacing all of the children with their parent preserves the volume covered
by that branch of the tree.


## A Concrete Plan Example

The overview uses the HexCell plan `[7, 1]`. Code `7` splits the original cell
into eight children, and code `1` splits child 0 into children 0.0 and 0.1:

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

If 0.0 and 0.1 both request derefinement and the edge-depth test passes, their
parent, child 0, becomes a leaf. The resulting compact plan is `[7]`.

Returning all the way to the original cell requires another eligible
sibling-family collapse. On a later adaptation cycle, if children 0 through 7
all request derefinement and the boundary checks pass, the root can become a
leaf and the plan becomes empty.

One derefinement pass can collapse many families in different branches, but it
collapses at most one level along a particular lineage. A parent is tested
while its children still reflect the candidate tree for that pass; it is not
tested again after one of those children has been collapsed. Successive
adaptation cycles can therefore walk back through several refinement levels.
The original mesh cell is the root of the plan and is the coarsest state to
which that plan can return.


## Where Derefinement Fits in Plan Processing

Derefinement does not simply erase split codes from the plan supplied at the
start of a cycle. FVMAdapt first constructs a candidate topology and makes its
shared boundaries consistent:

```text
plan from the previous cycle
        ↓ replay and associate request tags with current leaves
add requested refinements
        ↓
balance the candidate cell plans
        ↓
derive compatible candidate face and edge plans
        ↓
replay the candidate tree and test requested sibling families
        ↓
collapse eligible parents and serialize the final cell plan
        ↓
derive final shared plans and construct the adapted mesh
```

The intermediate balanced plan matters because a new refinement request or a
neighbor's shared-boundary requirements may change the candidate tree before
derefinement is considered. The derefinement check then sees the applicable
cell, face, and edge structure rather than acting only on the raw tags.

After eligible cell branches are collapsed, the final cell plans again drive
the shared face and edge plans. A neighboring cell can therefore retain a
subdivision on a shared boundary even when a cell on the other side has been
coarsened. Derefinement removes an eligible cell split; it does not
independently discard subdivisions still required by the final shared plans.

For the refinement-side plan merge and balancing loop, see
[Refinement Plans and Balancing](@ref fvmadapt_plans_and_balancing).


## Plan History and Data Transfer

The cell plan remains rooted at an original mesh cell across adaptation
cycles. It is cumulative in the sense that it describes the complete current
refinement tree below that original cell, not merely the splits added in the
latest cycle. Refinement grows branches of this plan, while derefinement prunes
branches from it.

FVMAdapt also compares the previous and final plans to construct a
previous-mesh-to-current-mesh cell mapping. This mapping has a different
purpose from refinement lineage:

- The **cell plan** describes the current tree below an original cell and is
  carried across cycles.
- The **transfer mapping**, exposed internally through `cell2parent`/`c2p`
  data, relates cells from the immediately previous mesh to cells in the new
  mesh. It is rebuilt for each adaptation transition.

For the `[7, 1]` to `[7]` example, previous cells 0.0 and 0.1 both map to the
new cell 0. The grid-interface interpolation path recognizes that several
previous fine cells feed one coarsened cell and forms a volume-weighted average
of their cell data. The transfer mapping should therefore not be read as the
persistent genealogy of a cell; the plan carries the longer-lived topological
history.


## Summary

The main derefinement ideas to keep separate are:

- **Derefinement request:** tag value `2`; it asks FVMAdapt to remove
  refinement but does not guarantee that the request can be accepted.
- **Sibling family:** all direct children of one split parent; this is the unit
  that can be collapsed.
- **Eligible parent:** a split parent whose children are all requested leaf
  cells and whose boundary edges pass the depth check.
- **Collapse:** deletion of that parent's children from the temporary
  refinement tree, making the parent a leaf in the new plan.
- **Original-cell root:** the coarsest state represented by the plan; FVMAdapt
  cannot derefine that lineage past the original mesh cell.
