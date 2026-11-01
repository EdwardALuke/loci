# Coarsening {#fvmadapt_derefinement}

Coarsening restores a parent whose immediate children are all active leaves.
Solver tags must request the complete sibling family. Deeper trees coarsen from
the bottom up; unrelated cells cannot be merged.

The parent recovers its saved geometry. Shared face divisions needed by a finer
neighbor remain in the output. Surviving points stay in place.

```cpp
auto children = cells.refine(0, CenterMethod::weighted);
cells.coarsen(0);
cells.refine(0); // reuse the accepted children and centers
```

Balancing may refuse a merge and restore the previous split. Any rank's proposal
failure preserves the preceding accepted histories. Output uses the framework's
normal failure handling after acceptance.

New-cell donor lists contain the old leaves covered by each restored parent.
The solver uses these lists and consistent cell measures to transfer fields.
Geometry history does not store solver values or interpolation weights.

`RefinementHistory` can discard unused constructions through `compactHistory()`.
Compaction is separate from ordinary coarsening. The boundary-defined
`CornerSubdivision` path retains its accepted constructions.
