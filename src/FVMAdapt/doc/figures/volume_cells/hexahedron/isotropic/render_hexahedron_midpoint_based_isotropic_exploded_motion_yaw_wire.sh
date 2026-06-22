#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

tex_file="$script_dir/hexahedron_midpoint_based_isotropic_exploded.tex"
base_name=hexahedron_midpoint_based_isotropic_exploded_motion_yaw_wire
out_gif="$script_dir/${base_name}.gif"

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/fvmadapt-hexahedron-midpoint-motion-yaw-wire.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT

frames_dir="$work_dir/frames"
mkdir -p "$frames_dir"

idx=0
render_frame() {
  local amount=$1
  local yaw=$2
  frame_base=$(printf 'frame_%03d' "$idx")
  frame_dir="$work_dir/tex_${frame_base}"
  mkdir -p "$frame_dir"
  pdflatex -halt-on-error -interaction=nonstopmode \
    -jobname="$frame_base" \
    -output-directory="$frame_dir" \
    "\\def\\fvmadaptExplodeAmount{$amount}\\def\\fvmadaptExplodeScale{1.69}\\def\\fvmadaptYawAngle{$yaw}\\input{$tex_file}" >/dev/null 2>&1
  mutool draw -c rgba -r 180 -o "$frames_dir/${frame_base}.png" \
    "$frame_dir/${frame_base}.pdf" >/dev/null 2>&1
  idx=$((idx + 1))
}

ease_amount() {
  awk -v i="$1" -v n="$2" 'BEGIN {
    pi = 3.141592653589793
    t = i / n
    printf "%.5f", 0.5 - 0.5 * cos(pi * t)
  }'
}

interpolate() {
  awk -v a="$1" -v b="$2" -v e="$3" 'BEGIN {
    printf "%.5f", a + (b - a) * e
  }'
}

base_yaw=8
left_yaw=-10
right_yaw=10

# Briefly show the assembled cell before the pieces separate.
for _ in $(seq 1 10); do
  render_frame 0 "$base_yaw"
done

# Explode from the assembled state at the normal view.
for i in $(seq 1 30); do
  render_frame "$(ease_amount "$i" 30)" "$base_yaw"
done

# While fully exploded, yaw left, sweep right, then return to the normal view.
for i in $(seq 1 24); do
  ease=$(ease_amount "$i" 24)
  render_frame 1 "$(interpolate "$base_yaw" "$left_yaw" "$ease")"
done

for i in $(seq 1 40); do
  ease=$(ease_amount "$i" 40)
  render_frame 1 "$(interpolate "$left_yaw" "$right_yaw" "$ease")"
done

for i in $(seq 1 24); do
  ease=$(ease_amount "$i" 24)
  render_frame 1 "$(interpolate "$right_yaw" "$base_yaw" "$ease")"
done

# Collapse back to the assembled state so the loop restarts cleanly.
for i in $(seq 29 -1 0); do
  render_frame "$(ease_amount "$i" 30)" "$base_yaw"
done

# Use full-size frames rather than ImageMagick's cropped sub-frame optimization.
# Some GIF viewers display optimized sub-frames as if they were clipped.
convert -dispose Background -delay 10 -loop 0 \
  "$frames_dir"/frame_*.png +repage "$out_gif"

echo "wrote $out_gif"
