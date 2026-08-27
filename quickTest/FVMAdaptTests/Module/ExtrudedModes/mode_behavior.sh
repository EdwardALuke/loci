#!/bin/sh

set -eu

if [ "$#" -ne 6 ]; then
  echo "usage: $0 fixture mode mesh.vog vogcheck.log topology.h5 coordinates.dat" >&2
  exit 2
fi

fixture=$1
mode=$2
mesh=$3
check_log=$4
topology=$5
coordinates=$6

element_count()
{
  dataset=$1
  count=$(
    h5dump -H -d "/elements/$dataset" "$topology" 2>/dev/null |
      awk '/DATASPACE  SIMPLE/ { print $5; exit }'
  )
  if [ -z "$count" ]; then
    count=0
  fi
  printf '%s\n' "$count"
}

coordinate_planes()
{
  column=$1
  awk -v column="$column" '{ print $column }' "$coordinates" |
    sort -gu |
    awk 'END { print NR }'
}

hexes=$(element_count hexahedra)
prisms=$(element_count prism)
tetrahedra=$(element_count tetrahedra)
pyramids=$(element_count pyramid)
general=$(element_count GeneralCellNfaces)
classified=$((hexes + prisms + tetrahedra + pyramids + general))
cells=$(
  h5dump -a "/file_info/numCells" "$mesh" |
    awk '/\(0\):/ { print $2; found = 1; exit }
         END { if (!found) exit 1 }'
)

if [ "$classified" -ne "$cells" ]; then
  echo "classified cell count $classified does not match mesh cell count $cells" >&2
  exit 1
fi

printf 'fixture=%s\n' "$fixture"
printf 'mode=%s\n' "$mode"
sh ../mesh_behavior.sh "$mesh" "$check_log"
printf 'cell_type.hex=%s\n' "$hexes"
printf 'cell_type.prism=%s\n' "$prisms"
printf 'cell_type.tetrahedron=%s\n' "$tetrahedra"
printf 'cell_type.pyramid=%s\n' "$pyramids"
printf 'cell_type.general=%s\n' "$general"
printf 'coordinate_planes.x=%s\n' "$(coordinate_planes 1)"
printf 'coordinate_planes.y=%s\n' "$(coordinate_planes 2)"
printf 'coordinate_planes.z=%s\n' "$(coordinate_planes 3)"
