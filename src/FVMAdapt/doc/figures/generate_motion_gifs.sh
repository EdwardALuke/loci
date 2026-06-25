#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
figures_root=$script_dir
script_path="$script_dir/generate_motion_gifs.sh"
toolchain_checked=0
static_svg_toolchain_checked=0

usage() {
  echo "usage: $0 [--figures-root DIR]" >&2
  exit 2
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --figures-root)
      [ "$#" -ge 2 ] || usage
      figures_root=$(cd "$2" && pwd -P)
      shift 2
      ;;
    *)
      usage
      ;;
  esac
done

render_specs=(
  "volume_cells/diamond/render_diamond_face_based_isotropic_exploded_motion_yaw_wire.sh|volume_cells/diamond/diamond_face_based_isotropic_exploded_motion_yaw_wire.gif"
  "volume_cells/general_polyhedron/render_general_polyhedron_face_based_isotropic_exploded_motion_yaw_wire.sh|volume_cells/general_polyhedron/general_polyhedron_face_based_isotropic_exploded_motion_yaw_wire.gif"
  "volume_cells/hexahedron/isotropic/render_hexahedron_midpoint_based_isotropic_exploded_motion_yaw_wire.sh|volume_cells/hexahedron/isotropic/hexahedron_midpoint_based_isotropic_exploded_motion_yaw_wire.gif"
  "volume_cells/hexahedron/anisotropic/render_hexahedron_anisotropic_split_codes.sh|volume_cells/hexahedron/anisotropic/hexahedron_split_code_1_zeta_exploded_motion_yaw_wire.gif volume_cells/hexahedron/anisotropic/hexahedron_split_code_2_eta_exploded_motion_yaw_wire.gif volume_cells/hexahedron/anisotropic/hexahedron_split_code_3_eta_zeta_exploded_motion_yaw_wire.gif volume_cells/hexahedron/anisotropic/hexahedron_split_code_4_xi_exploded_motion_yaw_wire.gif volume_cells/hexahedron/anisotropic/hexahedron_split_code_5_xi_zeta_exploded_motion_yaw_wire.gif volume_cells/hexahedron/anisotropic/hexahedron_split_code_6_xi_eta_exploded_motion_yaw_wire.gif"
  "volume_cells/prism/anisotropic/render_prism_split_code_1_axial_exploded_motion_yaw_wire.sh|volume_cells/prism/anisotropic/prism_split_code_1_axial_exploded_motion_yaw_wire.gif"
  "volume_cells/prism/anisotropic/render_prism_split_code_2_transverse_exploded_motion_yaw_wire.sh|volume_cells/prism/anisotropic/prism_split_code_2_transverse_exploded_motion_yaw_wire.gif"
  "volume_cells/prism/isotropic/render_prism_midpoint_based_isotropic_exploded_motion_yaw_wire.sh|volume_cells/prism/isotropic/prism_midpoint_based_isotropic_exploded_motion_yaw_wire.gif"
  "volume_cells/pyramid/render_pyramid_midpoint_based_isotropic_exploded_motion_yaw_wire.sh|volume_cells/pyramid/pyramid_midpoint_based_isotropic_exploded_motion_yaw_wire.gif"
  "volume_cells/tetrahedron/render_tetrahedron_midpoint_based_isotropic_exploded_motion_yaw_wire.sh|volume_cells/tetrahedron/tetrahedron_midpoint_based_isotropic_exploded_motion_yaw_wire.gif"
)

require_toolchain() {
  local missing=()
  local tool

  for tool in pdflatex mutool awk seq; do
    if ! command -v "$tool" >/dev/null 2>&1; then
      missing+=("$tool")
    fi
  done

  if ! imagemagick_cmd=$(find_imagemagick); then
    missing+=("ImageMagick (magick or convert)")
  fi

  if command -v pdflatex >/dev/null 2>&1 && ! check_tikz; then
    missing+=("TikZ/PGF for pdflatex")
  fi

  if [ "${#missing[@]}" -ne 0 ]; then
    echo "ERROR: FVMAdapt motion GIF generation needs additional utilities." >&2
    echo "Missing:" >&2
    for tool in "${missing[@]}"; do
      echo "  - $tool" >&2
    done
    echo >&2
    echo "Install the missing utilities and rerun the documentation build." >&2
    echo "Suggested packages:" >&2
    echo "  Debian/Ubuntu: sudo apt install texlive-latex-base texlive-pictures mupdf-tools imagemagick" >&2
    echo "  macOS/Homebrew: install MacTeX or BasicTeX with TikZ/PGF, then brew install mupdf-tools imagemagick" >&2
    echo "Alternatively, install from a tarball-with-docs archive that already contains generated documentation." >&2
    exit 1
  fi

  export FVMADAPT_IMAGEMAGICK="$imagemagick_cmd"
}

require_static_svg_toolchain() {
  local missing=()
  local tool

  for tool in pdflatex dvisvgm; do
    if ! command -v "$tool" >/dev/null 2>&1; then
      missing+=("$tool")
    fi
  done

  if command -v pdflatex >/dev/null 2>&1 && ! check_tikz; then
    missing+=("TikZ/PGF for pdflatex")
  fi

  if [ "${#missing[@]}" -ne 0 ]; then
    echo "ERROR: FVMAdapt static SVG generation needs additional utilities." >&2
    echo "Missing:" >&2
    for tool in "${missing[@]}"; do
      echo "  - $tool" >&2
    done
    echo >&2
    echo "Install the missing utilities and rerun the documentation build." >&2
    echo "Suggested packages:" >&2
    echo "  Debian/Ubuntu: sudo apt install texlive-latex-base texlive-pictures dvisvgm" >&2
    echo "  macOS/Homebrew: install MacTeX or BasicTeX with TikZ/PGF and dvisvgm" >&2
    echo "Alternatively, install from a tarball-with-docs archive that already contains generated documentation." >&2
    exit 1
  fi
}

find_imagemagick() {
  if [ -n "${FVMADAPT_IMAGEMAGICK:-}" ]; then
    if command -v "$FVMADAPT_IMAGEMAGICK" >/dev/null 2>&1; then
      printf '%s\n' "$FVMADAPT_IMAGEMAGICK"
      return 0
    fi
    return 1
  fi

  if command -v magick >/dev/null 2>&1; then
    printf '%s\n' magick
  elif command -v convert >/dev/null 2>&1; then
    printf '%s\n' convert
  else
    return 1
  fi
}

check_tikz() {
  local probe_dir
  local status
  probe_dir=$(mktemp -d "${TMPDIR:-/tmp}/fvmadapt-tikz-probe.XXXXXX")

  cat > "$probe_dir/tikz_probe.tex" <<'EOF'
\documentclass{article}
\usepackage{tikz}
\pagestyle{empty}
\begin{document}
\begin{tikzpicture}
  \draw (0,0) -- (1,0);
\end{tikzpicture}
\end{document}
EOF

  set +e
  pdflatex -halt-on-error -interaction=nonstopmode \
    -output-directory="$probe_dir" \
    "$probe_dir/tikz_probe.tex" >/dev/null 2>&1
  status=$?
  set -e
  rm -rf "$probe_dir"
  return "$status"
}

ensure_toolchain() {
  if [ "$toolchain_checked" -eq 0 ]; then
    require_toolchain
    toolchain_checked=1
  fi
}

ensure_static_svg_toolchain() {
  if [ "$static_svg_toolchain_checked" -eq 0 ]; then
    require_static_svg_toolchain
    static_svg_toolchain_checked=1
  fi
}

newer_input_exists() {
  local output=$1
  local input_dir
  input_dir=$(dirname "$output")

  find -L "$input_dir" -maxdepth 1 \
    \( -name '*.tex' -o -name 'render*.sh' \) \
    -newer "$output" -print -quit | grep -q .
}

is_static_svg_source() {
  local tex=$1
  local base
  base=$(basename "$tex")

  case "$base" in
    *_common.tex|*_exploded.tex|*_template.tex)
      return 1
      ;;
  esac

  return 0
}

needs_static_svg_render() {
  local tex=$1
  local output=$2

  if [ ! -f "$output" ]; then
    return 0
  fi

  if [ "$tex" -nt "$output" ] || [ "$script_path" -nt "$output" ]; then
    return 0
  fi

  if newer_input_exists "$output"; then
    return 0
  fi

  return 1
}

render_static_svg() {
  local tex=$1
  local output=$2
  local tex_dir
  local tex_file
  local job
  local build_dir
  local status

  tex_dir=$(cd "$(dirname "$tex")" && pwd -P)
  tex_file=$(basename "$tex")
  job=${tex_file%.tex}
  build_dir=$(mktemp -d "${TMPDIR:-/tmp}/fvmadapt-static-svg.XXXXXX")

  set +e
  (
    cd "$tex_dir" &&
    pdflatex -halt-on-error -interaction=nonstopmode \
      -output-directory="$build_dir" \
      "$tex_file" >/dev/null 2>&1 &&
    dvisvgm --pdf --no-fonts --exact \
      -o "$output" "$build_dir/$job.pdf" >/dev/null 2>&1
  )
  status=$?
  set -e

  rm -rf "$build_dir"

  if [ "$status" -ne 0 ]; then
    echo "ERROR: failed to render static SVG for $tex" >&2
    exit "$status"
  fi
}

render_static_svgs() {
  local tex
  local output
  local checked_count=0
  local generated_count=0

  while IFS= read -r -d '' tex; do
    if ! is_static_svg_source "$tex"; then
      continue
    fi

    output=${tex%.tex}.svg
    checked_count=$((checked_count + 1))

    if [ -L "$output" ]; then
      rm -f "$output"
    fi

    if needs_static_svg_render "$tex" "$output"; then
      ensure_static_svg_toolchain
      render_static_svg "$tex" "$output"
      generated_count=$((generated_count + 1))
    fi

    if [ ! -f "$output" ]; then
      echo "ERROR: renderer did not create $output" >&2
      exit 1
    fi
  done < <(find -L "$figures_root" -type f -name '*.tex' -print0)

  echo "Checked $checked_count FVMAdapt static SVG(s); generated $generated_count."
}

needs_render() {
  local render_script=$1
  local output=$2

  if [ ! -f "$output" ]; then
    return 0
  fi

  if [ "$render_script" -nt "$output" ] || [ "$script_path" -nt "$output" ]; then
    return 0
  fi

  if newer_input_exists "$output"; then
    return 0
  fi

  return 1
}

generated_count=0
checked_count=0

render_static_svgs

for spec in "${render_specs[@]}"; do
  render_rel=${spec%%|*}
  outputs_rel=${spec#*|}
  render_script="$figures_root/$render_rel"
  render_needed=0

  if [ ! -f "$render_script" ]; then
    echo "ERROR: missing renderer $render_script" >&2
    exit 1
  fi

  for output_rel in $outputs_rel; do
    output="$figures_root/$output_rel"
    if [ -L "$output" ]; then
      rm -f "$output"
    fi
    checked_count=$((checked_count + 1))
    if needs_render "$render_script" "$output"; then
      render_needed=1
    fi
  done

  if [ "$render_needed" -eq 1 ]; then
    ensure_toolchain
    echo "Generating ${render_rel#$figures_root/}"
    case "$render_rel" in
      volume_cells/hexahedron/anisotropic/render_hexahedron_anisotropic_split_codes.sh)
        bash "$render_script" --gifs-only
        ;;
      *)
        bash "$render_script"
        ;;
    esac
  fi

  for output_rel in $outputs_rel; do
    output="$figures_root/$output_rel"
    if [ ! -f "$output" ]; then
      echo "ERROR: renderer did not create $output" >&2
      exit 1
    fi
    if [ "$render_needed" -eq 1 ]; then
      generated_count=$((generated_count + 1))
    fi
  done
done

echo "Checked $checked_count FVMAdapt motion GIF(s); generated $generated_count."
