#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$script_dir"
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/fvmadapt-hexahedron-anisotropic.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT
render_static=1

if [ "${1:-}" = "--gifs-only" ]; then
  render_static=0
  shift
fi

if [ "$#" -ne 0 ]; then
  echo "usage: $0 [--gifs-only]" >&2
  exit 2
fi

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
    "$tex_file" >/dev/null 2>&1
  dvisvgm --pdf --no-fonts --exact \
    -o "$out_svg" "$build_dir/${job}.pdf" >/dev/null 2>&1
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
      "\\def\\fvmadaptExplodeAmount{$amount}\\def\\fvmadaptExplodeScale{1.35}\\def\\fvmadaptYawAngle{$yaw}\\input{$tex_file}" >/dev/null 2>&1
    mutool draw -c rgba -r 120 -o "$frames_dir/${frame_base}.png" \
      "$frame_dir/${frame_base}.pdf" >/dev/null 2>&1
    idx=$((idx + 1))
  }

  local base_yaw=8
  local left_yaw=-10
  local right_yaw=10

  for _ in $(seq 1 3); do
    render_frame 0 "$base_yaw"
  done

  for i in $(seq 1 8); do
    render_frame "$(ease_amount "$i" 8)" "$base_yaw"
  done

  for i in $(seq 1 6); do
    ease=$(ease_amount "$i" 6)
    render_frame 1 "$(interpolate "$base_yaw" "$left_yaw" "$ease")"
  done

  for i in $(seq 1 10); do
    ease=$(ease_amount "$i" 10)
    render_frame 1 "$(interpolate "$left_yaw" "$right_yaw" "$ease")"
  done

  for i in $(seq 1 6); do
    ease=$(ease_amount "$i" 6)
    render_frame 1 "$(interpolate "$right_yaw" "$base_yaw" "$ease")"
  done

  for i in $(seq 7 -1 0); do
    render_frame "$(ease_amount "$i" 7)" "$base_yaw"
  done

  "${FVMADAPT_IMAGEMAGICK:-convert}" -dispose Background -delay 20 -loop 0 \
    "$frames_dir"/frame_*.png +repage "$out_gif"
}

for spec in "${specs[@]}"; do
  code=${spec%%:*}
  name=${spec#*:}
  if [ "$render_static" -eq 1 ]; then
    render_static_svg "$code" "$name"
    echo "wrote hexahedron_split_code_${code}_${name}.svg"
  fi
  render_motion_gif "$code" "$name"
  echo "wrote hexahedron_split_code_${code}_${name}_exploded_motion_yaw_wire.gif"
done
