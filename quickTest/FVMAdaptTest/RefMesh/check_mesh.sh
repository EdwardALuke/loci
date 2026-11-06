#!/bin/sh
# This fixture fills the unit cube; refinement must preserve its existing points.
set -eu
before=$1
after=$2

awk -F= '
  FNR==NR { if ($1=="cells") old_cells=$2; if ($1=="nodes") old_nodes=$2; next }
  $1=="cells" { cells=$2 }
  $1=="nodes" { nodes=$2 }
  $1=="convexity" { convex=($2=="pass") }
  $1=="volume.Main" { volume=$2; measured=1 }
  END {
    if (!(old_cells>0 && old_nodes>0 && cells>old_cells && nodes>old_nodes && convex && measured &&
          volume>1-1e-10 && volume<1+1e-10)) {
      print "Refinement failed cell count, convexity or unit-cube volume checks" > "/dev/stderr"; exit 1
    }
  }' "$before.summary" "$after.summary"

awk '
  /======Minimum cell volume/ {
    found=1
    if ($NF !~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ ||
        !($NF+0>0 && $NF+0<=1+1e-10)) bad=1
  }
  END { if (!found || bad) { print "Missing or nonpositive cell volume" > "/dev/stderr"; exit 1 } }
  ' "$after.check"

before_nodes=$(awk -F= '$1=="nodes" { print $2 }' "$before.summary")
after_nodes=$(awk -F= '$1=="nodes" { print $2 }' "$after.summary")
awk -v before_nodes="$before_nodes" -v after_nodes="$after_nodes" '
  NF!=3 { bad=1; next }
  { for(i=1;i<=3;++i) if ($i !~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/) bad=1 }
  FNR==NR { original[sprintf("%.12g %.12g %.12g",$1,$2,$3)]=1; count++; next }
  {
    after_count++
    for(i=1;i<=3;++i) if ($i < -1e-12 || $i > 1+1e-12) bad=1
    delete original[sprintf("%.12g %.12g %.12g",$1,$2,$3)]
  }
  END {
    for (point in original) bad=1
    if (count!=before_nodes || after_count!=after_nodes || bad) { print "Incomplete coordinates, moved point or point outside the unit cube" > "/dev/stderr"; exit 1 }
  }' "$before.dat" "$after.dat"
printf '%s\n' "$before -> $after: PASSED"
