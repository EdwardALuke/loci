# FVMAdapt quick tests

The suite separates direct library behavior from complete adaptation workflows:

- `Core` checks refinement-tree and topology behavior through the C++ library.
- `Module` runs the Loci tools on representative HexCell, Prism, mixed-cell,
  XML, parameter-file, and tag-file inputs.
- `Illustrations` creates optional human-inspection artifacts and is not part of
  the pass/fail suite.

## Behavioral reference data

Each Module scenario writes a small `*.actual` manifest and compares it with a
checked-in file under that scenario's `dats/` directory. A manifest records the
observable mesh behavior that should survive an implementation refactor:

```text
cells=512
faces=1728
nodes=729
convexity=pass
volume.Main=64
```

The references intentionally do not contain raw plan bytes, entity ids, HDF5
layout, schedule text, or output ordering. Those are implementation details and
may change without changing the refinement behavior.

Run the same checked-in suite in each worktree, using that worktree's own build:

```sh
make -C quickTest/FVMAdaptTests \
  LOCI_BASE="$PWD/OBJ" \
  TEST_BASE="$PWD/quickTest"
```

For a refactoring branch, base the branch on the test commit (or create a local
integration branch with that ancestry). Both worktrees then compare their own
outputs with exactly the same reference files; no build products are shared.

## Regression intent

Two tests document defects found while developing this suite:

- the Core level-refinement test requires every HexCell and Prism branch to
  reach the requested depth;
- the three-rank UniformHex case requires an uneven MPI tag distribution to
  preserve the final refinement request.

These tests may fail on an unfixed baseline. They should pass on a branch that
contains the corresponding production fixes.

`Module/Prism` also provides an `empty-partition` diagnostic target. It sends a
one-cell tag file through three ranks. Tag reading reaches the later plan-output
stage, but the current implementation cannot write a plan when ranks own no
cells. Keep that diagnostic separate from the default suite until the
empty-domain workflow is fixed.

`empty-partition-known-failure` is an opt-in wrapper for the current defect. It
passes only when the marker reaches cell-plan output and reproduces the known
crash. An unexpected success or a different error fails the wrapper. The
positive `empty-partition` target remains the behavior required from a future
fix; neither target is part of the default suite.
