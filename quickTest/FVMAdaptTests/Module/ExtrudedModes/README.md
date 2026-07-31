# Extruded-mesh split modes

This test sends localized cell tags through `marker` and `refmesh` on two
meshes that are one cell thick in z:

- `thin_hex.ugrid` is a 3-by-3 hexahedral slab. Its center column is tagged,
  and each cell is longer in x than y and much thinner in z.
- `thin_prism.ugrid` is a square formed by two extruded triangles. One prism
  is tagged.

The checked manifests document these outcomes:

| Mesh | Mode | Observed split directions | Cells | Cell types after refinement |
| --- | --- | --- | ---: | --- |
| Hex | 0, edge-length driven | x and y | 18 | 12 hex, 6 general |
| Hex | 1, no z refinement | x only | 12 | 12 hex |
| Hex | 2, isotropic | x, y, and z | 30 | 24 hex, 6 general |
| Prism | 0, edge-length driven | x and y | 4 | 4 hex |
| Prism | 1, no z refinement | x and y | 4 | 4 hex |
| Prism | 2, isotropic | x, y, and z | 7 | 6 hex, 1 general |

The element types above come from the topology written by `vogcheck`. The
coordinate-plane counts in `dats/` independently verify which directions
acquired new nodes. Volume and convexity are checked for every result.

From the repository root, run only this workflow with:

```sh
make -C quickTest/FVMAdaptTests/Module/ExtrudedModes \
  LOCI_BASE="$PWD/OBJ" \
  TEST_BASE="$PWD/quickTest"
```

Remove its generated results with the corresponding `clean` target.
