# FVMAdapt Derefinement {#fvmadapt_derefinement}

Derefinement, also called *coarsening*, removes existing refinement. In
FVMAdapt, this does not mean merging an arbitrary collection of neighboring
fine cells. It means replacing one complete family of leaf cells with their
immediate parent in the refinement tree.

```text
before                                      after

parent                                      parent (leaf)
├─ child 0  leaf, requests derefinement
├─ child 1  leaf, requests derefinement
└─ ... all remaining children
```

The parent in this diagram is not another live mesh cell stored alongside its
children. It is retained in the history. Coarsening selects the parent instead
of its children for the output mesh, keeping the accepted split locations for
later reuse.

## The Sibling-Family Test

Suppose a split tree object `P` has direct children `C0`, `C1`, and so on.
For tag-based adaptation, FVMAdapt can collapse that split only when all of
the following are true:

1. `P` is currently split.
2. Every direct child of `P` is a leaf. None is still split into further children.
3. Every direct child requests derefinement.
4. The resulting cell satisfies the enabled balancing policies.

An explicit coarsening request names the parent directly; it still replaces
the complete family of leaf children.

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

One tag-based derefinement pass can collapse many families in different branches, but it
collapses at most one level along a particular lineage. Successive
adaptation cycles can therefore walk back through several refinement levels.
The original mesh cell is the root of the history and is the coarsest state to
which that history can return.

The parent recovers its saved geometry. Shared face divisions needed by a finer
neighbor remain in the output. Surviving points stay in place. Repeating an
accepted split after coarsening reuses its children and centers.
