# Hexahedron Anisotropic Split Figures

This directory contains hexahedron split-code figures for the primitive
anisotropic cases. The fully combined `7` case is the all-direction split shown
in `../isotropic/`.

The implementation stores `HexCell::mySplitCode` as a three-bit local-direction
mask in `xi`, `eta`, `zeta` order:

- `0`: no split.
- `1` / `001`: split in the `zeta` direction. This creates two child
  hexahedra.
- `2` / `010`: split in the `eta` direction. This creates two child hexahedra.
- `3` / `011`: split in the `eta` and `zeta` directions. This creates four
  child hexahedra.
- `4` / `100`: split in the `xi` direction. This creates two child hexahedra.
- `5` / `101`: split in the `xi` and `zeta` directions. This creates four
  child hexahedra.
- `6` / `110`: split in the `xi` and `eta` directions. This creates four child
  hexahedra.
- `7` / `111`: split in all three directions. This is treated as the isotropic
  hexahedron case.

## Files

The short `hexahedron_split_code_*.tex` files are wrappers around
`hexahedron_split_code_template.tex`, which draws the static split scaffold for
one code. The matching `*_exploded.tex` wrappers use
`hexahedron_split_code_exploded_template.tex`, which draws the separated
wireframe children.

Run `render_hexahedron_anisotropic_split_codes.sh` from this directory to
regenerate the static SVG files and the assembled, explode, yaw, return, and
collapse GIFs for codes `1` through `6`.
