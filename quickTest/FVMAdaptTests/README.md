# FVMAdapt quick tests

The suite separates direct library behavior from complete adaptation workflows:

- `Core` checks refinement-tree and topology behavior through the C++ library.
- `Module` runs the Loci tools on representative HexCell, Prism, mixed-cell,
  XML, parameter-file, tag-file, and thin extruded-mesh inputs.
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

`Module/ExtrudedModes` applies localized cell tags to thin hexahedral and
prismatic meshes. Its references compare modes 0, 1, and 2 using cell-type
counts and the number of coordinate planes in each direction. This documents
the observable difference between edge-length-driven, no-z, and isotropic
refinement through the complete `marker`/`refmesh` workflow.

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

Three tests document defects found while developing this suite:

- the Core level-refinement test requires every HexCell and Prism branch to
  reach the requested depth;
- the three-rank UniformHex case requires an uneven MPI tag distribution to
  preserve the final refinement request;
- the three-rank one-prism case requires `marker` and `refmesh` to operate when
  some ranks have empty cell ranges.

These tests may fail on an unfixed baseline. They should pass on a branch that
contains the corresponding production fixes.
