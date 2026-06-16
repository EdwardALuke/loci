#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$script_dir"
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/fvmadapt-hexahedron-anisotropic.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT

specs=(
  "1:zeta"
  "2:eta"
  "3:eta_zeta"
  "4:xi"
  "5:xi_zeta"
  "6:xi_eta"
)

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

render_static_svg() {
  local code=$1
  local name=$2
  local tex_file="$script_dir/hexahedron_split_code_${code}_${name}.tex"
  local out_svg="$script_dir/hexahedron_split_code_${code}_${name}.svg"
  local job="hexahedron_split_code_${code}_${name}"
  local build_dir="$work_dir/static_${code}_${name}"

  mkdir -p "$build_dir"
  pdflatex -halt-on-error -interaction=nonstopmode \
    -jobname="$job" \
    -output-directory="$build_dir" \
    "$tex_file" >/dev/null
  dvisvgm --pdf --no-fonts --exact \
    -o "$out_svg" "$build_dir/${job}.pdf" >/dev/null
}

render_motion_gif() {
  local code=$1
  local name=$2
  local tex_file="$script_dir/hexahedron_split_code_${code}_${name}_exploded.tex"
  local base_name="hexahedron_split_code_${code}_${name}_exploded_motion_yaw_wire"
  local out_gif="$script_dir/${base_name}.gif"
  local frames_dir="$work_dir/${base_name}_frames"
  local idx=0

  mkdir -p "$frames_dir"

  render_frame() {
    local amount=$1
    local yaw=$2
    local frame_base
    local frame_dir
    frame_base=$(printf 'frame_%03d' "$idx")
    frame_dir="$work_dir/${base_name}_${frame_base}"
    mkdir -p "$frame_dir"
    pdflatex -halt-on-error -interaction=nonstopmode \
      -jobname="$frame_base" \
      -output-directory="$frame_dir" \
      "\\def\\fvmadaptExplodeAmount{$amount}\\def\\fvmadaptExplodeScale{1.35}\\def\\fvmadaptYawAngle{$yaw}\\input{$tex_file}" >/dev/null
    mutool draw -r 150 -o "$frames_dir/${frame_base}.png" \
      "$frame_dir/${frame_base}.pdf" >/dev/null 2>&1
    idx=$((idx + 1))
  }

  local base_yaw=8
  local left_yaw=-10
  local right_yaw=10

  for _ in $(seq 1 10); do
    render_frame 0 "$base_yaw"
  done

  for i in $(seq 1 30); do
    render_frame "$(ease_amount "$i" 30)" "$base_yaw"
  done

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

  for i in $(seq 29 -1 0); do
    render_frame "$(ease_amount "$i" 30)" "$base_yaw"
  done

  convert -dispose Background -delay 10 -loop 0 \
    "$frames_dir"/frame_*.png +repage "$out_gif"
}

for spec in "${specs[@]}"; do
  code=${spec%%:*}
  name=${spec#*:}
  render_static_svg "$code" "$name"
  echo "wrote hexahedron_split_code_${code}_${name}.svg"
  render_motion_gif "$code" "$name"
  echo "wrote hexahedron_split_code_${code}_${name}_exploded_motion_yaw_wire.gif"
done
