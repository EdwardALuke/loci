# FVMAdapt tests

Small tool workflows and C++ `TEST_CASE` tests for the original FVMAdapt.
The original `RefMesh` grid and references are unchanged; additional cases
come from branch `369-fvmadapt_tests`. No FVMAdapt2 APIs or VTK output are used.

| Folder | Behavior checked |
| --- | --- |
| `RefMesh` | Original node-tag regression, continuing from a saved plan and from a generated mesh. |
| `PlanInputs` | XML-region and parameter-file refinement through `marker` and `refmesh`. |
| `UniformHex` | Refine/derefine 64 hexes; separately check the last cell's tag across three MPI ranks. |
| `Prism` | Refine/derefine one prism; separately check refinement with more ranks than cells. |
| `Mixed` | Refine either side of a shared hex/prism face. |
| `ExtrudedModes` | Thin hex/prism meshes in modes 0, 1, and 2, including split directions. |
| `Core` | Split plans, nested refinement/derefinement, shared faces, edge ordering, and two-level geometric refinement. |
| `FaceOutput` | Query the face writer with two empty ranks; preserve face and cell counts. |
| `CellParents` | Map old tetrahedral-cell descendants to new cells after full or partial collapse, in serial and on three ranks. |

Tool cases check mesh counts, convexity, and volume. `ExtrudedModes` also
checks cell types and coordinate planes. The original `RefMesh` test compares
cell-centroid coordinates. `CellParents` checks local old/new cell indices,
not global entity IDs or the order of mapping entries. Tool comparisons do not
require identical entity IDs or plan bytes. These tests are not exhaustive
coverage of balancing or online adaptation.

## Run and clean

From the repository root, select a built or installed Loci tree:

```sh
make -C quickTest/FVMAdaptTest LOCI_BASE=/path/to/loci/OBJ
make -C quickTest/FVMAdaptTest PlanInputs LOCI_BASE=/path/to/loci/OBJ
make -C quickTest/FVMAdaptTest/UniformHex final_cell_tag LOCI_BASE=/path/to/loci/OBJ
make -C quickTest/FVMAdaptTest/CellParents serial LOCI_BASE=/path/to/loci/OBJ
make -C quickTest/FVMAdaptTest clean
```

The full run attempts every case and prints a final `TestResults` summary,
even after failures; it returns nonzero if any case fails. No `-k` is needed.
The suite also participates in the normal `quickTest` target.

Each run starts in a fresh `<folder>/work/<case>/`. Commands and diagnostics
are saved in `run.log` there; meshes, plans, and schedules stay beside it for
inspection. `clean` removes all these generated files and results, not the
inputs or references, and does not require `LOCI_BASE`.

Each folder's Makefile contains its test commands. `common.mk` only handles
setup, isolated working directories, reporting, and cleanup. MPI commands
come from `quickTest/test.conf`; the empty-rank and final-tag cases explicitly
use three ranks. The selected Loci/module paths are printed before running.

Tools must be built; XML tests require libxml2 support. Put `h5dump`, the MPI
launcher, and `timeout` on `PATH`. Overrides include `H5DUMP=/path/to/h5dump`
and `TIMEOUT=gtimeout` on macOS. Each executable has a `TIME_LIMIT=120s`
default to bound hangs; increase it for slower machines.

## Known failures

On `dev` at `b87ecb23`, these remain ordinary failing tests:

- `UniformHex/final_cell_tag`: 64 cells instead of 71; the last tag is lost.
- `Prism/empty_ranks`: three-rank `marker` crashes during grid reading.
- `Core/plans`: edge ordering and two-level hex/prism refinement fail;
  the other six C++ tests pass.
- `FaceOutput/empty_ranks`: empty face ownership creates an invalid face
  range in the writer. This is separate from the prism reader failure.
- `CellParents/serial`: full collapse reports one new cell, but its mapping
  retains child indices. Refinement and partial-collapse controls pass.
  The three-rank case also encounters the grid-reader empty-rank defect.

Source fixes belong in separate changes. Do not regenerate references or
turn these failures into expected successes to hide defects.
