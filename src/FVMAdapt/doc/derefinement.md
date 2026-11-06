# Coarsening {#fvmadapt_derefinement}

Derefinement, also called *coarsening*, removes existing refinement. In
FVMAdapt, this does not mean merging an arbitrary collection of neighboring
fine cells. It means replacing one complete family of leaf cells with their
immediate parent in the refinement tree.

```text
before                                      after

parent                                      parent (leaf)
├─ child 0  leaf, requests coarsening
├─ child 1  leaf, requests coarsening
└─ ... all remaining children
```

The parent in this diagram is not another live mesh cell stored alongside its
children. It is retained in the history. Coarsening selects the parent instead
of its children for the output mesh, keeping the accepted split locations for
later reuse.

## The Sibling-Family Test

FVMAdapt can accept a coarsening only when all of the following are true:

1. The parent is currently split.
2. Every direct child is a leaf. None is still split into further children.
3. Every direct child requests coarsening.
4. The resulting cell satisfies the enabled balancing policies.

Replacing just one child with the parent would overlap the remaining children.
Simply removing that child would leave a hole. The complete sibling family is
therefore the unit of coarsening.

Requested merges are applied to a candidate history before balancing. If a
neighbor's refinement requires the family to remain split, balancing restores
its previous split. A coarsening request does not guarantee a coarser output cell.

## Coarsening through Several Levels

In the [cell-history example](@ref fvmadapt_plans_and_balancing), children 0.0
and 0.1 can be replaced by their parent, child 0. On a later adaptation cycle,
children 0 through 7 can be replaced by the root if they all request coarsening
and the balancing checks pass.

One tag-based adaptation can collapse many families in different branches,
but removes at most one generation along a lineage. The original mesh cell is
the coarsest state to which that lineage can return.

The parent recovers its saved geometry. Shared face divisions needed by a finer
neighbor remain in the output. Surviving points stay in place. Repeating an
accepted split after coarsening reuses its children and centers.

## Data Transfer

The refinement history remains rooted at the original mesh cell across cycles.
The transfer mapping instead relates cells in the immediately previous mesh to
cells in the new mesh. It is rebuilt for each adaptation transition.

New-cell donor lists contain the old leaves covered by each restored parent.
The solver uses these lists and consistent cell measures to transfer fields.
Geometry history does not store solver values or interpolation weights.
