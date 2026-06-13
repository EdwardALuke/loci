# FVMAdapt Existing Tests {#fvmadapt_testing}

The current checked-in smoke test for FVMAdapt is
`quickTest/FVMAdaptTest/Makefile`. It exercises the command-line adaptation
path through `marker`, `refmesh`, `vogcheck`, and `extract`, then compares
sorted coordinate output against checked-in reference data.

## Test Flow

The default target builds and checks three refinement outputs:

1. `testGrid.plan`
   - `vogcheck testGrid` validates the input grid.
   - `extract -ascii testGrid 0 x` creates a tag file.
   - `marker -g testGrid.vog -tag testGrid.tag -o testGrid.plan` creates the
     first adaptation plan.

2. `testGridRef.vog`
   - `refmesh -g testGrid.vog -r testGrid.plan -o testGridRef.vog` applies the
     first plan.
   - `vogcheck testGridRef` validates the refined grid.
   - `extract -ascii testGridRef 0 x y z` writes coordinates for comparison
     against `quickTest/FVMAdaptTest/dats/testGridRef.dat`.

3. `testGridRef2.vog` and `testGridRef3.vog`
   - The test creates additional tags and plans from the original and refined
     grids.
   - The resulting coordinate outputs are compared against
     `testGridRef2.dat` and `testGridRef3.dat`.

The comparisons are intentionally coordinate-order tolerant: both generated
and reference coordinate files are sorted before `diff`.

## Running The Test

After configuring and building Loci, run the quick-test target through the
repo's normal test machinery, or run the makefile directly from a configured
build tree where `LOCI_BASE` and `TEST_BASE` are defined.

For direct inspection, the relevant file is:

```bash
quickTest/FVMAdaptTest/Makefile
```

## Coverage Notes

This smoke test is useful because it exercises the public command-line route
through plan generation, remeshing, output validation, and coordinate
comparison. It is not a complete behavioral specification for the adaptation
library.

Known limitations:

- It does not individually assert every plan byte code.
- It does not isolate every cell type or orientation case.
- It does not document all derefinement paths.
- It does not prove that a refactor preserved subtle child ordering invariants
  unless those changes affect the tested coordinate outputs.

For documentation work, treat this test as the current smoke gate, not as proof
that every internal helper is covered.
