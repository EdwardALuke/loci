#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/../../../../../../.." && pwd)

tex_file="$script_dir/prism_cell_face_to_volume_centroid.tex"
base_name=prism_cell_face_to_volume_centroid
out_svg="$script_dir/${base_name}.svg"
out_gif="$script_dir/${base_name}_yaw.gif"

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/fvmadapt-prism-yaw.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT

static_dir="$work_dir/static"
frames_dir="$work_dir/frames"
mkdir -p "$static_dir" "$frames_dir"

pdflatex -halt-on-error -interaction=nonstopmode \
  -output-directory="$static_dir" "$tex_file" >/dev/null
dvisvgm --pdf --no-fonts \
  --output="$static_dir/${base_name}.svg" "$static_dir/${base_name}.pdf" >/dev/null
cp "$static_dir/${base_name}.svg" "$out_svg"

idx=0
for angle in -10 -8 -6 -4 -2 0 2 4 6 8 10 8 6 4 2 0 -2 -4 -6 -8; do
  frame_base=$(printf 'frame_%03d' "$idx")
  frame_dir="$work_dir/tex_${frame_base}"
  mkdir -p "$frame_dir"
  pdflatex -halt-on-error -interaction=nonstopmode \
    -jobname="$frame_base" \
    -output-directory="$frame_dir" \
    "\\def\\fvmadaptYawAngle{$angle}\\input{$tex_file}" >/dev/null
  mutool draw -r 180 -o "$frames_dir/${frame_base}.png" \
    "$frame_dir/${frame_base}.pdf" >/dev/null
  idx=$((idx + 1))
done

# Use full-size frames rather than ImageMagick's cropped sub-frame optimization.
# Some GIF viewers display optimized sub-frames as if they were clipped.
convert -dispose Background -delay 22 -loop 0 \
  "$frames_dir"/frame_*.png +repage "$out_gif"

echo "wrote $out_svg"
echo "wrote $out_gif"
