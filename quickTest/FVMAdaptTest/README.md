# FVMAdapt tests

Tests for mesh refinement, derefinement, and cell mappings.

| Folder | Behavior checked |
| --- | --- |
| `RefMesh` | Node-tag refinement from saved plans and generated meshes; unit-cube geometry and retained points. |
| `PlanInputs` | XML-region and parameter-file refinement. |
| `UniformHex` | Hex refinement/derefinement and tag distribution across MPI ranks. |
| `Prism` | Prism refinement/derefinement, including empty MPI ranks. |
| `Mixed` | Refinement across a shared hex/prism face. |
| `ExtrudedModes` | Refinement modes and split directions on thin meshes. |
| `Core` | Directional and isotropic splits, history reload, coarsening, and shared faces. |
| `FaceOutput` | Face output with empty MPI ranks. |
| `CellParents` | Old-to-new cell mappings after refinement and derefinement. |

## Run and clean

Requires Loci tools built with libxml2 support, an MPI compiler and launcher,
`h5dump`, and `timeout`.

From the repository root, set `LOCI_BASE` to your Loci build or installation:

```sh
# Run all tests
make -C quickTest/FVMAdaptTest LOCI_BASE=/path/to/loci/OBJ
# Run one group
make -C quickTest/FVMAdaptTest PlanInputs LOCI_BASE=/path/to/loci/OBJ
# Remove generated files
make -C quickTest/FVMAdaptTest clean
```

To keep results outside the source checkout, create a work directory and run:

```sh
make -C /path/to/work -f /path/to/loci/quickTest/FVMAdaptTest/Makefile \
  LOCI_BASE=/path/to/loci/OBJ
```

The work directory receives the same group folders and result summary. Set
`H5DUMP=/path/to/h5dump` if the HDF5 tool is not on `PATH`.

Results from the source-directory command are summarized in `quickTest/FVMAdaptTest/TestResults`; each case's
log is in `<folder>/work/<case>/run.log`. Any failure makes `make` return
nonzero. Cleanup preserves test inputs and references.

