#!/usr/bin/env bash
set -euo pipefail

source_root=${1:-../..}
html_dir=${2:-html}

source_root=$(cd "$source_root" && pwd -P)
mkdir -p "$html_dir"
html_dir=$(cd "$html_dir" && pwd -P)

manifest=$(mktemp "${TMPDIR:-/tmp}/loci-doxygen-assets.XXXXXX")
trap 'rm -f "$manifest"' EXIT

while IFS= read -r -d '' figures_dir; do
    figures_dir=$(cd "$figures_dir" && pwd -P)
    while IFS= read -r -d '' asset; do
        rel_path=${asset#"$figures_dir"/}
        target_path="$html_dir/figures/$rel_path"
        printf '%s\t%s\n' "$target_path" "$asset" >> "$manifest"
    done < <(find -L "$figures_dir" -type f \( -iname '*.svg' -o -iname '*.gif' \) -print0)
done < <(
    find "$source_root" \
        \( -path "$source_root/OBJ" \
        -o -path "$source_root/ext" \
        -o -path "$source_root/contrib" \
        -o -path "$source_root/tmpcopy" \
        -o -path "$source_root/loci_install" \
        -o -path "$source_root/docs/doxygen/html" \
        -o -path "$source_root/docs/doxygen/latex" \) -prune \
        -o -type d -path '*/doc/figures' -print0
)

if [ ! -s "$manifest" ]; then
    echo "No Doxygen figure assets found under $source_root."
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
