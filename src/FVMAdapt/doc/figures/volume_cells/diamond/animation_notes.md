# Diamond Animation Notes

This directory keeps the editable TikZ source and generated web assets for the
diamond refinement figures used in the FVMAdapt overview notes.

## Source Model

`diamond_cell.tex` draws an octahedron, also describable as a square
bipyramid. In this context it represents the diamond-cell shape produced by
pyramid-like decompositions: six original vertices, twelve original edges, and
eight triangular faces.

`diamond_face_based_isotropic.tex` draws the same cell with the
face-centroid, edge-midpoint, and volume-centroid split lines.

`diamond_face_based_isotropic_exploded.tex` defines the wireframe exploded
vertex-cell geometry for the combined motion/yaw animation. The view is
controlled by `\fvmadaptYawAngle`, which changes the projected `x` and `y`
basis vectors while leaving `z` vertical.

The explosion amount is controlled by `\fvmadaptExplodeAmount`, where `0`
places the duplicated child cells in their assembled positions and `1` moves
them to the fully exploded positions. Explosion distance can be scaled with
`\fvmadaptExplodeScale`; the combined renderer uses `1.35`.

The child cells are drawn in their natural, un-exploded coordinates and moved
with TikZ `scope` shifts. The current explosion heuristic keeps the lower-apex
child fixed, moves the four equator children outward in the `xi`/`eta` plane,
and lifts the upper-apex child in `zeta`.

The source also sets a fixed canvas with `\path[use as bounding box]`. This is
important for animations: the yaw sweep changes the projected image extents,
and a fixed canvas keeps every rendered frame the same size.

## Regeneration

Run this command from anywhere in the repository:

```sh
bash src/FVMAdapt/doc/figures/volume_cells/diamond/render_diamond_face_based_isotropic_exploded_motion_yaw_wire.sh
```

The combined wire motion/yaw script regenerates:

- `diamond_face_based_isotropic_exploded_motion_yaw_wire.gif`

The GIF is a generated build artifact and is intentionally not tracked in git.

The combined wire motion/yaw GIF starts assembled at the normal yaw angle,
explodes, sweeps left and right while fully exploded, returns to the normal yaw
angle, and then collapses before restarting. It uses 41 frames with
`-delay 20`.

The script intentionally avoids ImageMagick's `-layers Optimize` pass. That
optimization rewrites frames as cropped subimages with offsets, which can look
truncated in some GIF viewers.

## Toolchain

Required local commands:

- `pdflatex` with TikZ/PGF
- `mutool`
- ImageMagick `magick` or `convert`
- `awk`
- `seq`
