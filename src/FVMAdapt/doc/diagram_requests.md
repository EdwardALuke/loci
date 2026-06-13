# FVMAdapt Diagram Requests {#fvmadapt_diagram_requests}

This page lists figures that would make the FVMAdapt developer documentation
easier to read. The intent is to draw implementation reference figures, not
polished marketing diagrams. Use simple lines, explicit labels, and enough
spacing that the labels remain readable in generated Doxygen HTML.

Suggested source format: Inkscape SVG. Keep text as text when possible so the
figure remains editable. A good destination for committed source images is
`src/FVMAdapt/doc/figures/`; the generated Doxygen page can later reference
those SVG files directly.

## Drawing Conventions

- Use zero-based labels that match the code: `node 0`, `edge[0]`, `child[0]`,
  `face[0]`, and so on.
- Show local axes whenever a split code depends on direction.
- Label whether a frame is edge-local, face-local, or cell-local.
- Use arrows for stored edge direction, not just geometric adjacency.
- Keep global mesh ids out of these figures unless the figure is explicitly
  about input/output numbering.
- Prefer one concept per figure. Crowded diagrams are harder to review against
  the implementation.

## Figure 1: Plan Tree Replay

Suggested file: `src/FVMAdapt/doc/figures/fvmadapt_plan_tree.svg`

Draw a small breadth-first refinement tree and the matching plan vector.

Required labels:

- root entity
- leaf entity
- split entity
- parent and child
- queue/traversal order
- example plan such as `[1, 0, 1]`

Purpose:

This figure should explain that a plan is a replay recipe for a tree, not the
refined mesh itself. It should be generic enough to apply to edges, faces, and
cells.

## Figure 2: Edge Split

Suggested file: `src/FVMAdapt/doc/figures/fvmadapt_edge_split.svg`

Draw one edge before and after a midpoint split.

Required labels:

- `head`
- `tail`
- midpoint node
- `child[0]`
- `child[1]`
- stored direction from `head` to `tail`
- optional note that `needReverse` can change traversal order for an owning
  face without changing the stored edge geometry

Purpose:

This is the simplest plan example: edge code `0` is a leaf and code `1`
creates two child edges.

## Figure 3: General Face Split

Suggested file: `src/FVMAdapt/doc/figures/fvmadapt_general_face_split.svg`

Draw a polygonal face, preferably a pentagon, split through edge midpoints and
a computed center node.

Required labels:

- original face
- edge midpoints
- computed center node
- one child face per parent edge
- `Face::child[i]`

Purpose:

This should distinguish `Face` splitting from directional `QuadFace`
splitting. General-face code `1` creates a ring of quadrilateral child faces.

## Figure 4: QuadFace Numbering And Split Modes

Suggested file: `src/FVMAdapt/doc/figures/fvmadapt_quadface_numbering.svg`

Draw the local quadrilateral face used by `QuadFace`.

Required labels:

- `node 0`, `node 1`, `node 2`, `node 3`
- `edge[0]` through `edge[3]` placed along the edge middles
- arrows showing edge directions
- local x and local y axes
- split code `1`: local y split, `childy[0]`, `childy[1]`
- split code `2`: local x split, `childx[0]`, `childx[1]`
- split code `3`: four-quadrant layout with `child[0]`, `child[1]`,
  `child[2]`, `child[3]`

Purpose:

This figure should replace the need to mentally parse the ASCII diagram in
`quadface.h`.

## Figure 5: HexCell Local Directions And Split Codes

Suggested file: `src/FVMAdapt/doc/figures/fvmadapt_hex_split_codes.svg`

Draw a hexahedral cell with local x, y, and z directions and the split-code bit
weights.

Required labels:

- local x, y, z axes
- code `4` = split x
- code `2` = split y
- code `1` = split z
- combined examples: `3` = y+z, `6` = x+y, `7` = x+y+z
- optional face labels from `HexCell`: RIGHT, LEFT, FRONT, BACK, UP, DOWN

Purpose:

This figure should make clear that hex cells use a true three-bit local
direction mask.

## Figure 6: Prism Split Modes

Suggested file: `src/FVMAdapt/doc/figures/fvmadapt_prism_split_modes.svg`

Draw a triangular prism and show the two split families.

Required labels:

- polygonal xy cross-section
- axial/z direction
- code `1` = axial/z split into two children
- code `2` = xy cross-section split into `nfold` children
- code `3` = both, producing `2*nfold` children
- `nfold = 3` for the original triangular prism
- optional note that `nfold` can become `4` for child prisms after side-face
  splitting

Purpose:

This figure should make clear that prism codes are not the same as the hex
`xyz` bit mask.

## Figure 7: General Cell To Diamond Cells

Suggested file: `src/FVMAdapt/doc/figures/fvmadapt_diamond_cell_decomposition.svg`

Draw a simple general cell decomposed into diamond cells.

Required labels:

- original general cell
- one diamond cell associated with a source cell node
- diamond apex nodes `0` and `1`
- fold `n`
- `2n + 2` nodes
- `2n` faces
- `4n` edges
- `n` faces sharing node `0` and `n` faces sharing node `1`

Purpose:

This figure should explain why general cells do not follow the hex/prism
directional split-code formalism.

## Figure 8: Plan Flow Across Dimensions

Suggested file: `src/FVMAdapt/doc/figures/fvmadapt_plan_flow.svg`

Draw the high-level consistency pipeline.

Required labels:

- tags or region input
- `cellPlan`
- extracted or merged `facePlan`
- extracted `edgePlan`
- `balancedCellPlan`
- `balancedFacePlan`
- `balancedEdgePlan`
- fine nodes, fine faces, fine cells, and transfer/restart data

Purpose:

This figure should show why changing one local split-code convention can affect
neighbor consistency, lower-dimensional extraction, and output topology.

## Figure 9: Face-Local Versus Cell-Local Orientation

Suggested file: `src/FVMAdapt/doc/figures/fvmadapt_orientation_frames.svg`

Draw the same quadrilateral face once in face-local order and once as seen from
a neighboring cell.

Required labels:

- face-local node/edge order
- cell-local node/edge order
- `orientCode`
- `orient_splitCode_f2c()`
- `orient_childID_f2c()`
- `orient_edgeID_f2c()`
- `orient_edgeID_c2f()`

Purpose:

This is the orientation warning figure. It should communicate that face plans
and child ids are meaningful only after choosing the frame in which they are
being interpreted.
