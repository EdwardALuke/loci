# Prism Animation Notes

This directory keeps the editable TikZ source and generated web assets for the
prism refinement figures used in the FVMAdapt overview notes.

## Source Model

`prism_cell.tex` draws the prism cell, and `prism_midpoint_based_isotropic.tex`
draws the same cell with the face-centroid, edge-midpoint, and volume-centroid
split lines.

`prism_midpoint_based_isotropic_exploded.tex` defines the wireframe exploded prism
vertex-cell geometry for the combined motion/yaw animation. The view is
controlled by `\fvmadaptYawAngle`, which changes the projected `x` and `y`
basis vectors while leaving `z` vertical.
The explosion amount is controlled by `\fvmadaptExplodeAmount`, where `0`
places the duplicated child cells in their assembled positions and `1` moves
them to the fully exploded positions. Explosion distance can be scaled with
`\fvmadaptExplodeScale`; the combined renderer uses `1.56`.

The child cells are drawn in their natural, un-exploded coordinates and moved
with TikZ `scope` shifts. The current explosion heuristic is origin-anchored:
the child at the `\xi`, `\eta`, and `\zeta` minimum corner stays fixed, while
the other children move mostly along the dominant coordinate direction of their
parent vertices.

The source also sets a fixed canvas with `\path[use as bounding box]`. This is
important for animations: the yaw sweep changes the projected image extents,
and a fixed canvas keeps every rendered frame the same size.

## Regeneration

Run this command from anywhere in the repository:

```sh
bash src/FVMAdapt/doc/figures/volume_cells/prism/isotropic/render_prism_midpoint_based_isotropic_exploded_motion_yaw_wire.sh
```

The combined wire motion/yaw script regenerates:

- `prism_midpoint_based_isotropic_exploded_motion_yaw_wire.gif`

The combined wire motion/yaw GIF starts assembled at the normal yaw angle,
explodes, sweeps left and right while fully exploded, returns to the normal yaw
angle, and then collapses before restarting. It uses 158 frames with
`-delay 10`.

The script intentionally avoids ImageMagick's `-layers Optimize` pass. That
optimization rewrites frames as cropped subimages with offsets, which can look
truncated in some GIF viewers.

## Toolchain

Required local commands:

- `pdflatex`
- `dvisvgm`
- `mutool`
- `convert` from ImageMagick
- `awk`
- `seq`
