#!/bin/sh

# Summarize mesh counts and volume after checking convexity.
set -eu

if [ "$#" -ne 2 ]; then
  echo "usage: $0 mesh.vog vogcheck.log" >&2
  exit 2
fi

mesh=$1
check_log=$2

attribute_value()
{
  attribute=$("${H5DUMP:-h5dump}" -a "/file_info/$1" "$mesh")
  printf '%s\n' "$attribute" |
    awk '/\(0\):/ { print $2; found = 1; exit }
         END { if (!found) exit 1 }'
}

grep -Fq "All cells in this mesh pass the convexity check!" "$check_log"

cells=$(attribute_value numCells)
faces=$(attribute_value numFaces)
nodes=$(attribute_value numNodes)
# These fixtures have one volume component, named Main.
volume=$(awk '/^Volume of component / {
                printf "volume.%s=%s\n", $4, $6
                found = 1
              }
              END { if (!found) exit 1 }' "$check_log")

printf 'cells=%s\n' "$cells"
printf 'faces=%s\n' "$faces"
printf 'nodes=%s\n' "$nodes"
printf 'convexity=pass\n'
printf '%s\n' "$volume"
