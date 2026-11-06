# Module Interface {#fvmadapt_reference_subdivision}

See [splitting](@ref fvmadapt_geometric_splitting),
[history and balancing](@ref fvmadapt_plans_and_balancing), and
[cell geometry](@ref fvmadapt_reference_geometry) for the local operations.

## Distributed Loci Rules

Include `FVMAdapt/parametric_adapt.lh`.
The rules apply requests, match shared faces, balance levels and build output.
Facts live on the original cells:

| Fact | Meaning |
| --- | --- |
| `parametricCellHistory` | Accepted saved history; empty initializes from VOG corners |
| `parametricRequests` | Encoded edits; empty means none |
| `parametricAccepted` | All ranks accepted balancing and any enabled quality limits |
| `parametricError` | Local rejection reason, otherwise empty |
| `balancedParametricHistory` | Accepted new history, or the preceding history after rejection |
| `parametricDonorLeaves` | Old leaf indices contributing to each new leaf |
| `parametricParentLeaf` | Previous leaf index for each output leaf; -1 for a merge |

The donor lists use `ParametricCell::leafIds()` order. Initialize once with
empty histories and requests. Query acceptance, histories and errors before
requesting mesh output. Carry `balancedParametricHistory` into the next
`parametricCellHistory`, keeping the original connectivity and file IDs.

Inside the module, `ParametricCellData` owns accepted per-cell state. Its query
and edit methods choose the appropriate compact plan or working geometry.
`applyRequests()` and `balance()` edit independent candidates. `workingCell()`
supplies read-only working geometry; `workingCopy()` supplies an editable copy.
Expanding compact state can be expensive, so ordinary queries use the direct
cell-data methods.
The rules retain responsibility for neighbor reductions and iteration.
Shared face plans contain accepted local centers and edge divisions, without
physical geometry maps. Both sides use the original grid's shared face geometry.

To encode one isotropic request:

```cpp
auto tree = makeParametricCell(topology, corners, saved,
                              parametricHistoryLimits(10, 1000000), "corner", -1);
std::vector<ParametricLeafRequest> edits = {
  {tree.leafIds().front(), TensorSplit(), ParametricAction::isotropic}
};
auto request = encodeParametricRequests(tree, edits);
```

A coarsen request names the parent of a complete leaf family. An explicit cut
uses `ParametricAction::split` and a `SplitOperation` with the native family.
Stale or overlapping requests reject.

## Controls

| Control | Default and choices |
| --- | --- |
| `parametricSplitFamily` | `"native"`: standard cuts; `"corner"`: split around each vertex, including general cells; `"extruded"`: planar splitting of hex/prism; `"mixed"`: native hex/prism and corner splits elsewhere |
| `parametricAutomaticSplit` | `false`; choose native directions from current leaf edges |
| `parametricSplitAspectRatio` | `2`; finite and greater than one |
| `parametricMinEdgeLength` | `0`; automatic/geometric half-edge cutoff, in grid units; XML requires a positive value |
| `parametricRegionFile` | Empty; optional [XML region](README_xml) |
| `parametricSourceFile` | Empty; optional [spacing sources](spacing_sources.md); choose this or an XML region |
| `parametricCenterMethod` | `"reference"`; corner/extruded also accept `"weighted"` and `"integrated"` |
| `parametricExtrusionAxis` | `-1`; 0/1/2 preserves straight X/Y/Z thickness for hex/prism |
| `parametricMaxDepth` | `10`, at most `30` |
| `parametricMaxNodes` | `1000000` retained nodes per original cell |
| `parametricLevelDifference` | `1` |
| `parametricGrading` | `"vertices"`; also `"edges"` or `"faces"` |
| `parametricFaceBalance` | `"none"`; also `"majority"` or `"opposite"` |
| `parametricDirectionalBalance` | `false`; native hex/prism balance splits follow divided-face directions |
| `parametricMaxRounds` | `32` balancing rounds |
| `parametricRefinementLevels` | `1` subdivision per refine tag; coarsening always merges one generation |

Mixed mode uses native cuts in hexes and prisms and one-child-per-vertex splits
in tets, pyramids and general cells. It requires reference centers and full 3D
refinement, with midpoint cuts at mixed interfaces. Weighted/integrated mixed
interfaces and 2D extrusion are not supported.

Automatic selection compares mean edge lengths in each local direction. It
splits directions longer than the smallest mean by the aspect ratio; otherwise
it splits all directions. A direction is skipped when half its shortest edge
is below the cutoff. Tets and pyramids keep their isotropic patterns.
Explicit cuts and balancing do not use this selector.

XML regions select a leaf when any own vertex is inside. Regions and spacing
sources also visit newly created children until their size requirement is met.
Hanging points from neighbors do not affect selection. Isotropic requests stop
when half the shortest own edge is below the cutoff; extrusion ignores thickness.
These controls guide refinement, not a guaranteed minimum output edge length.
See [XML regions](README_xml) and [spacing sources](spacing_sources.md).

Histories have a 64 MiB limit per original cell. Parallel point records require
signed-int reference numerators and denominators. Unsupported coordinates or
exhausted limits reject the proposal on every rank.

## Mesh Quality

Quality checking is optional and defaults to `off`. It measures the physical
polygons after balancing, including divisions imposed by neighbors.

| Control | Default and meaning |
| --- | --- |
| `parametricQualityMode` | `"off"`; `"report"` measures quality; `"reject"` also rejects a proposal exceeding any enabled limit |
| `parametricQualityCentroid` | `"wireframe"`; also `"exact"`, matching Cartesian FVM metrics |
| `parametricMaxNonOrthogonality` | `-1`; maximum angle in degrees between a face normal and the line joining cell centers; boundaries use the cell-to-face line |
| `parametricMaxFaceWarpage` | `-1`; maximum angle in degrees between a polygon's area vector and its center-to-edge triangle normals |
| `parametricMaxChildVolumeRatio` | `-1`; largest/smallest immediate child volume at any split, summing active descendants |

`-1` disables a limit. Angle limits lie in `[0,180]`; a volume-ratio limit is
at least one. These are Cartesian mesh measures, not axisymmetric swept volumes.
The metric centroid choice is separate from the method used to place split points.

`parametricQualityReport` gives each maximum and its original cell ID, current
root-local leaf index, and original face ID when applicable. It also identifies
invalid output geometry, which rejects in either enabled mode. The usual
topology and subdivision checks still apply with quality off.

Online callers set `RefinementOptions::quality` and may pass a
`MeshQualityReport*` as the final `onlineRefineMesh()` argument. Rejection leaves
the accepted grid, histories and transfer database unchanged. No alternative
split or local retry is selected automatically. Base-face recording and accepted
checkpoint restoration report quality without rejecting on optional limits.

The command-line tools accept `-quality off|report|reject`,
`-quality-centroid wireframe|exact`, `-max-nonorthogonality`,
`-max-face-warpage` and `-max-child-volume-ratio`. Saved plans retain these
controls; command-line values override them.

## Solver Connection

Include `FVMAdapt/gridInterface.h` and load the `fvmadapt` rules. The collective
`Loci::onlineRefineMesh()` accepts the original VOG case name, `RefinementState`, and
file-order cell tags: 0 keeps, 1 refines, 2 requests coarsening. Only complete
tagged sibling families merge. Null tags replay saved history unless a geometric
input file requests refinement.
`RefinementOptions` supplies the controls above. Set `tagLocation="nodes"` to
mark the current mesh's vertices in output file order. Any own vertex tagged 1
refines a leaf; all own vertices must be tagged 2 to request coarsening.
Hanging vertices from finer neighbors do not mark the coarse leaf. The rule
interface also accepts `tagfile_par` through the existing node-tag reader
(text values: zero keeps, positive refines, negative coarsens). Use one tag input
per adaptation.

Set `options.regionFile="region.xml"` and a positive `options.minEdgeLength`
for XML refinement. It may accompany cell/node tags; region requirements
are enforced after those requests.

Set `options.sourceFile="sources.par"` for [segment-based spacing](spacing_sources.md).
It uses the same recursive refinement and acceptance rules. The default zero
`minEdgeLength` leaves resolution to the source spacing. Choose one geometric
input file per proposal; either may accompany cell/node tags.

The result is `refinedGridData`, with cell correspondence in `DataXFER_DB` as `c2p`
for the existing transfer rules. A rejected proposal retains the accepted mesh
and state. The caller owns field reconstruction and checkpoint storage.

For restart, load histories into a fresh `RefinementState`, then call
`Loci::initializeGridFromPlan(grid,rules,state,casename,options,error)`.
It rebuilds the mesh without applying tags or geometric request files.
`readParametricPlan()` reads files produced by the command-line tools; solvers
may store `state.histories` with their own checkpoint metadata. The preceding
face state belongs in `state.faces`, not in the options. Original `cellPlan`
checkpoints and the old mode/level signatures are no longer supported.

Optional cell fields use the existing interface:

```cpp
RefinementOptions options;
options.splitFamily = "corner";
options.cellState = true;   // refinementDepth, rootCellFileNumber, cellChange
options.edgeLengths = true; // maxEdgeLength, maxEdgeLengthXY
```

Pass these options to `onlineRefineMesh()`, then install the returned
`refinedGridData` with `setupFVMGridFromContainer()`. The fields are stores on
the installed `geom_cells`. Both options default to false and are independent.
Use the same options on every rank.

Depth starts at zero. Root IDs use the original mesh's file numbering, not
the installed entity IDs. Change flags describe the final accepted split or
merge, including balancing. Replaying saved history reports `CELL_CHANGE_UNKNOWN`.
Edge lengths are straight distances between each leaf's own edge endpoints;
neighboring subdivisions do not shorten them. XY lengths project those edges
onto physical XY. They are not curved-edge arc lengths or cell diameters.

For any split family, set `RefinementOptions::faceHistory` and call
`recordRefinementBaseFaces()` before the first adaptation. Install the returned
grid with the `refinedGridData` overload of `setupFVMGridFromContainer()`, then
retain `getFaceState(facts)` in `state.faces`. `getFaceHistory(facts)` gives
oriented old-face donors and identifies new and removed internal faces.

Extruded faces use flat overlap areas. Native and corner faces use areas
integrated on the saved surfaces. Native internal faces are identified by their
exact plane in the root's reference coordinates; refinement does not rename them.
Physical fractions distribute each old face's value completely.

Mapped face transfer requires nested pieces: each overlap must be a whole old
or new piece. Crossing old/new cuts are rejected during face-history installation.
Area integration failure rejects the adaptation before output.

MPI-to-MPI restart with a changed processor count is tested for all three paths. Reload saved histories with `state.weights` unset, as Stream
does; replay rebuilds weights for subsequent adaptations. Restoring arbitrarily
partitioned weight stores is not supported by the current grid reader.
One-rank/MPI transitions remain unsupported because the
current Loci modes assign different original entity IDs.

## Marking and Mesh Files

`marker` balances requests and saves the accepted history. It does not construct
the output mesh. `refmesh` evaluates that history and writes a VOG:

```sh
mpirun -np 2 marker -g base.vog -ctag cells.txt -o refined.plan
mpirun -np 4 refmesh -g base.vog -r refined.plan -o refined.vog -c2p parents.txt
```

`refine` combines both operations in one grid load, with the same controls:

```sh
mpirun -np 2 refine -g base.vog -ctag cells.txt -o refined.vog -oplan refined.plan -c2p parents.txt
```

`-oplan` and `-c2p` are optional. With `refine -r previous.plan`, donor pairs
refer to that input mesh; without `-r`, they refer to original cells. All three
tools read the same plan format.

Tags use the current mesh's file order: positive refines, negative coarsens,
zero keeps. Use `-tag` for node tags, `-xml` for a region, or `-par` for spacing
sources. XML requires a positive `-tol`. `-levels N` requests N subdivisions
below each tagged leaf; coarsening removes one complete generation.

The default is isotropic splitting around each vertex (`-mode 2`) with weighted
centers. `-mode 1` preserves a Z extrusion; `-mode 0` selects automatic native
cuts. `-family`, `-center` and `-axis` select these controls individually.
Use `-family mixed -automatic 1` for automatic hex/prism cuts and isotropic
splits of other cell types. `-automatic 0` disables automatic selection.
`-balance 0|1|2` selects level grading, added majority-face refinement, or added
opposite-face refinement. `-help` lists the controls. The old `-fold` option is
unsupported; `-tol` limits geometric and automatic requests.

For another cycle, give `marker -r refined.plan` and tags for `refined.vog`.
Give `refmesh -pr refined.plan` to map new cells to that previous mesh; otherwise
parents refer to the original grid. `-c2p` writes one-based child/parent pairs,
including every donor of a merged cell.

Plans retain geometry controls and histories, so reconstruction needs no request
files. Parallel rank counts may change; serial/parallel transitions are rejected.
Original plan files must be regenerated. Keep each plan with its original VOG.

## Local C++ Representations

The local geometry and history classes are in `library/parametric/`, without
MPI communication or Loci stores.

`referenceElement(shape)` gives the local vertices, edges and faces.
`GeometryMap::evaluate(q)` gives a physical position and, where available,
derivatives. `CornerMap` implements the standard shapes. `PolygonMap` uses
[Wachspress coordinates](https://pmc.ncbi.nlm.nih.gov/articles/PMC3549276/)
over a regular reference polygon. `BezierMap` supports curved segments,
quads and hexes with additional control points. `RestrictedMap` evaluates
a child through its parent without refitting the parent geometry.

| Type | Retained state |
| --- | --- |
| `RefinementHistory` | Standard cuts and mapped child placements |
| `HexPlan` | Compact eight-child hex splits and accepted centers |
| `CornerSubdivision` | Cells defined by their shared faces, including pyramids and general cells |

`AffineSubdivision`, `MappedVolumeSubdivision` and `PlanarSubdivision` provide
cell operations over `RefinementHistory`. For a standard directional cut:

```cpp
AffineSubdivision cells(Shape::hex);
TensorSplit cut;
cut.cuts[0] = {Fraction(1,2)};
cells.splitLeaf(0, cut);
```

`CellKey` identifies an accepted construction lineage. A `CellPath` is only the
sequence of child ordinals. Local array indices are not output mesh numbers and
may change when records are reordered or compacted. `compactHistory()` discards
unused constructions separately from ordinary coarsening. `CornerSubdivision`
retains its accepted constructions and does not provide this compaction.

`geometryMeasure()` computes length, area or volume and its centroid;
`MeasureOptions` sets accuracy and work limits. Check its convergence result.
`geometryQuality()` and `splitQuality()` assess the mapped shapes. Small pieces
used only for integration are not necessarily small solver cells.
