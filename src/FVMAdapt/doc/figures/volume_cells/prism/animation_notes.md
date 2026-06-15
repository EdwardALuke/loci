# Prism Animation Notes

This directory keeps the editable TikZ source and generated web assets for the
prism face-refinement animation used in the FVMAdapt overview notes.

## Source Model

`prism_cell_face_to_volume_centroid.tex` defines the prism geometry once. The
view is controlled by `\fvmadaptYawAngle`, which changes the projected `x` and
`y` basis vectors while leaving `z` vertical.

The source also sets a fixed canvas with `\path[use as bounding box]`. This is
important for animations: the yaw sweep changes the projected image extents,
and a fixed canvas keeps every rendered frame the same size.

## Regeneration

Run this command from anywhere in the repository:

```sh
bash src/FVMAdapt/doc/figures/volume_cells/prism/render_prism_face_to_volume_centroid_yaw.sh
```

The script regenerates:

- `prism_cell_face_to_volume_centroid.svg`
- `prism_cell_face_to_volume_centroid_yaw.gif`

The GIF is built from 20 PNG frames rendered at yaw angles from `-10` to `10`
degrees and back. ImageMagick uses `-delay 22` for a slower playback rate and
`-loop 0` so browser playback loops forever.

The script intentionally avoids ImageMagick's `-layers Optimize` pass. That
optimization rewrites frames as cropped subimages with offsets, which can look
truncated in some GIF viewers.

## Toolchain

Required local commands:

- `pdflatex`
- `dvisvgm`
- `mutool`
- `convert` from ImageMagick
