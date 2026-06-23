#!/usr/bin/env bash
set -euo pipefail

asset_root=${1:-../..}
html_dir=${2:-html}
generator_source_root=${3:-}

asset_root=$(cd "$asset_root" && pwd -P)
if [ -n "$generator_source_root" ]; then
    generator_source_root=$(cd "$generator_source_root" && pwd -P)
fi
mkdir -p "$html_dir"
html_dir=$(cd "$html_dir" && pwd -P)

manifest=$(mktemp "${TMPDIR:-/tmp}/loci-doxygen-assets.XXXXXX")
trap 'rm -f "$manifest"' EXIT

find_figure_dirs() {
    find "$asset_root" \
        \( -path "$asset_root/OBJ" \
        -o -path "$asset_root/ext" \
        -o -path "$asset_root/contrib" \
        -o -path "$asset_root/tmpcopy" \
        -o -path "$asset_root/loci_install" \
        -o -path "$asset_root/docs/doxygen/html" \
        -o -path "$asset_root/docs/doxygen/latex" \) -prune \
        -o -type d -path '*/doc/figures' -print0
}

find_generator() {
    figures_dir=$1
    local_generator="$figures_dir/generate_motion_gifs.sh"
    if [ -f "$local_generator" ]; then
        printf '%s\t%s\n' "$local_generator" ""
        return 0
    fi

    if [ -z "$generator_source_root" ]; then
        return 1
    fi

    rel_path=${figures_dir#"$asset_root"/}
    for source_figures_dir in \
        "$generator_source_root/$rel_path" \
        "$generator_source_root/src/$rel_path"; do
        source_generator="$source_figures_dir/generate_motion_gifs.sh"
        if [ -f "$source_generator" ]; then
            printf '%s\t%s\n' "$source_generator" "$figures_dir"
            return 0
        fi
    done

    return 1
}

generators=0
while IFS= read -r -d '' figures_dir; do
    if generator_record=$(find_generator "$figures_dir"); then
        generator=${generator_record%%$'\t'*}
        output_root=${generator_record#*$'\t'}
        if [ -n "$output_root" ]; then
            bash "$generator" --figures-root "$output_root"
        else
            bash "$generator"
        fi
        generators=$((generators + 1))
    fi
done < <(find_figure_dirs)

if [ "$generators" -gt 0 ]; then
    echo "Ran $generators Doxygen figure generator(s)."
fi

while IFS= read -r -d '' figures_dir; do
    figures_dir=$(cd "$figures_dir" && pwd -P)
    while IFS= read -r -d '' asset; do
        rel_path=${asset#"$figures_dir"/}
        target_path="$html_dir/figures/$rel_path"
        printf '%s\t%s\n' "$target_path" "$asset" >> "$manifest"
    done < <(find -L "$figures_dir" -type f \( -iname '*.svg' -o -iname '*.gif' \) -print0)
done < <(find_figure_dirs)

if [ ! -s "$manifest" ]; then
    echo "No Doxygen figure assets found under $asset_root."
    exit 0
fi

duplicates=$(cut -f1 "$manifest" | sort | uniq -d)
if [ -n "$duplicates" ]; then
    echo "Doxygen figure asset path collision:" >&2
    while IFS= read -r target_path; do
        echo "  ${target_path#$html_dir/}" >&2
        awk -F '\t' -v target="$target_path" '$1 == target { print "    " $2 }' "$manifest" >&2
    done <<< "$duplicates"
    exit 1
fi

count=0
while IFS=$'\t' read -r target_path asset; do
    mkdir -p "$(dirname "$target_path")"
    cp -pL "$asset" "$target_path"
    count=$((count + 1))
done < "$manifest"

echo "Staged $count Doxygen figure asset(s) in ${html_dir#$PWD/}."
