# History and Balancing {#fvmadapt_plans_and_balancing}

A history records parents, children and accepted split locations. Each parent
selects either no split or one saved split. A selected cell with no children
is a leaf. Coarsening changes this selection without deleting the saved points.

## Local Histories

| Type | Retained state |
| --- | --- |
| `RefinementHistory` | Standard cuts and mapped child placements |
| `HexPlan` | Compact eight-child hex splits and accepted centers |
| `CornerSubdivision` | Cells defined by their shared faces, including pyramids and general cells |

`AffineSubdivision`, `MappedVolumeSubdivision` and `PlanarSubdivision` provide
cell operations over `RefinementHistory`. The module's `ParametricCellData`
presents these representations through one rule-facing interface.
`applyRequests()` and `balance()` edit independent candidates.
`workingCell()` supplies read-only working geometry; `workingCopy()` supplies
an editable copy. Expanding compact state can be expensive, so ordinary queries
use the direct cell-data methods.

`CellKey` identifies an accepted construction lineage. A `CellPath` is only the
sequence of child ordinals. Local array indices are not output mesh numbers and
may change when records are reordered or compacted.

Saving retains the geometry and accepted choices. Loading validates them;
it does not recompute centers. Unknown or invalid records reject before replacing
accepted state. Keep saved histories with their original grid.

## Balancing in Loci

1. Cells contribute their accepted face and edge divisions.
2. Loci reductions combine contributions on shared entities. Conflicting centers reject.
3. Cells read neighbor depths and request any additional splits.
4. Repeat until no cell changes, or reject if a limit is reached.

The kernel performs local edits. The module owns communication, reductions,
iteration and collective acceptance.

`parametricLevelDifference` bounds the difference in split depth across faces,
edges or vertices, selected by `parametricGrading`. Depth counts cell splits,
not physical size ratios or separate directional depths.

`parametricFaceBalance="majority"` also refines a leaf when more than half its
faces are divided. `"opposite"` also handles two divided faces sharing no vertex.
The planar path counts edges. Only divisions used by current neighbors count;
unused saved splits do not prevent coarsening.

The usual response is isotropic, preserving a requested extrusion axis.
`parametricDirectionalBalance=true` lets native hexes and prisms respond in
the divided directions. Missing direction information uses a full split.
A refused merge restores its previous accepted split.
