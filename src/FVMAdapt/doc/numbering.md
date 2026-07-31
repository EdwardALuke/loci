# FVMAdapt Numbering Conventions {#fvmadapt_numbering}

By the time FVMAdapt operates, the mesh topology is represented by nodes,
edges, faces, and cells, together with the connections between them. A cell is
bounded by faces, a face is bounded by edges, and an edge has two endpoint
nodes. All four are mesh entities. In the original mesh topology, an interior
face is one shared face incident to two cells, not one face belonging to each
cell.

An entity's identity is different from its position in a face or cell's local
topology. Suppose an original mesh node `n` lies on a quadrilateral face and a
HexCell. The face may use `n` as its local node 0 while the HexCell uses the
same `n` as its local node 4. The cell on the other side of the face may assign
it another cell-local position. These local numbers describe the role that `n`
plays in each view; they do not identify different mesh nodes.

| Kind of number | Example | What it identifies |
| --- | --- | --- |
| Mesh entity identity | original mesh node `n` | Which mesh entity is being referenced |
| Local position | QuadFace node 0 or HexCell node 4 | The entity's role in one local topology |
| Fine-mesh output number | a numbered fine node | Its number after the adapted mesh is constructed |

FVMAdapt plans use local positions. Split directions, boundary edges, and child
positions are defined relative to the local topology of the planned cell, face,
or edge. Local positions are useful for describing those operations, but the
positions used by two neighboring cells do not have to match.

Shared entities provide the common references needed to compare these local
descriptions. Each original face has one ordered `face2node` list, which
defines its face-local node positions and normal direction. Each original edge
has one ordered pair of `edge2node` endpoints, which defines its stored
direction. A cell also has local positions for its faces, edges, and nodes, but
those positions are derived separately for that cell.

Before FVMAdapt compares or merges the face-plan contributions from adjacent
cells, it maps both contributions into the shared face's local ordering. It
likewise maps subdivisions from incident faces into a shared edge's stored
direction before combining them. This mapping is called *orientation*. It does
not require the neighboring local numberings to become identical. It makes
their node, edge, child, and split-direction references name the same parts of
the shared entity.

The coordinate names below are also local. They identify positions and
directions in a reference shape, not physical Cartesian axes. The mesh entities
may be skewed or stretched, and their edges do not have to be perpendicular or
equal in length.

In the figures, a `q` subscript marks a QuadFace-local number and an `H`
subscript marks a HexCell-local number. For example,
n<sub>q</sub><sup>(0)</sup> means QuadFace node 0, while
e<sub>H</sub><sup>(6)</sup> means HexCell edge 6.

An `Edge` has a stored direction from its head node to its tail node. A face or
cell may traverse the same mesh edge in the opposite direction, in which case
FVMAdapt records that reversal.


## QuadFace

A `QuadFace` has four corner-node positions and four boundary-edge positions.
Both are numbered from zero. The figure shows the face-local view used by
`QuadFace` split codes.

<img src="figures/numbering/quadface_numbering.svg"
     alt="QuadFace node and edge numbering"
     width="620">

QuadFace node 0 is the local origin. Face-local `x` runs from node 0 toward
node 1, and face-local `y` runs from node 0 toward node 3. These names
distinguish the two local edge directions. They do not mean that the physical
edges are perpendicular.

The boundary-edge numbers are:

| QuadFace edge | Directed endpoints |
| --- | --- |
| 0 | node 0 to node 1 |
| 1 | node 1 to node 2 |
| 2 | node 3 to node 2 |
| 3 | node 0 to node 3 |


## HexCell

A `HexCell` has six quad faces. Its local view is a reference hexahedron with
HexCell node 0 at `(xi, eta, zeta) = (0, 0, 0)`. The reference coordinates name
positions in that local shape; they do not require the physical cell to be a
cube.

The diagram shows node numbers without boxes and face numbers inside boxes.
Each face is identified by the local coordinate that remains constant on it.

<img src="figures/numbering/hexcell_numbering.svg"
     alt="HexCell local node and face numbering"
     width="620">

| HexCell face | Local coordinate |
| --- | --- |
| 0 | `xi = 1` |
| 1 | `xi = 0` |
| 2 | `eta = 1` |
| 3 | `eta = 0` |
| 4 | `zeta = 1` |
| 5 | `zeta = 0` |

The local node numbering is:

| HexCell node | `xi` | `eta` | `zeta` |
| --- | --- | --- | --- |
| 0 | 0 | 0 | 0 |
| 1 | 0 | 0 | 1 |
| 2 | 0 | 1 | 0 |
| 3 | 0 | 1 | 1 |
| 4 | 1 | 0 | 0 |
| 5 | 1 | 0 | 1 |
| 6 | 1 | 1 | 0 |
| 7 | 1 | 1 | 1 |

Equivalently, the local node number is `4*xi + 2*eta + zeta` for `xi`, `eta`,
and `zeta` values of zero or one.

The `HexCell` also has 12 local edges. The implementation stores each edge with
a head and tail node:

| HexCell edge | Endpoints (`head` to `tail`) |
| --- | --- |
| 0 | node 0 to node 4 |
| 1 | node 1 to node 5 |
| 2 | node 2 to node 6 |
| 3 | node 3 to node 7 |
| 4 | node 0 to node 2 |
| 5 | node 1 to node 3 |
| 6 | node 4 to node 6 |
| 7 | node 5 to node 7 |
| 8 | node 0 to node 1 |
| 9 | node 2 to node 3 |
| 10 | node 4 to node 5 |
| 11 | node 6 to node 7 |


## Prism

A `Prism` has two n-sided end faces and `nfold` quadrilateral side faces. A
triangular prism enters this path with `nfold = 3`. Refinement can create
`Prism` children with `nfold = 4`.

The local numbering is arranged in groups:

| Item | Local numbers | Meaning |
| --- | --- | --- |
| End faces | 0 and 1 | The two n-sided faces |
| Side faces | 2 through `nfold + 1` | The quadrilateral faces around the prism |
| First end-face nodes | 0 through `nfold - 1` | Nodes on end face 0 |
| Second end-face nodes | `nfold` through `2*nfold - 1` | Corresponding nodes on end face 1 |
| First end-face edges | 0 through `nfold - 1` | Edges around end face 0 |
| Second end-face edges | `nfold` through `2*nfold - 1` | Edges around end face 1 |
| Axial edges | `2*nfold` through `3*nfold - 1` | Edges joining corresponding end-face nodes |

For a local end-face node `i`, node `i + nfold` is the corresponding node on
the other end face. The axial edge with local number `2*nfold + i` joins those
two nodes. This arrangement gives the `Prism` its transverse and axial split
directions.


## Mapping a Face into a Cell-Local View

Each mesh face has an ordered `face2node` list. When FVMAdapt builds a
`HexCell`, it compares that list with the cell's expected order for the same
face and stores the result in `hexOrientCode`. The two views can differ in two
ways:

- They can choose different corners as local node 0.
- They can walk around the four corners in opposite directions.

Changing the starting corner shifts the local numbers around the face.
Changing the direction reverses their order. Neither change moves or rotates
the physical face.

The node order also sets the direction of the mesh face normal. The
cell-to-face relation tells FVMAdapt whether that normal points into or out of
the cell. FVMAdapt uses both pieces of information when it decides whether the
face and cell orders match or are reversed.

To *orient* a face plan means to translate it from one of these local views
into another. This is necessary on a shared face because, for example, a split
in one cell's face-local `x` direction may be a split in the other cell's
face-local `y` direction. FVMAdapt maps both cell contributions into the
shared face's `face2node` orientation before it compares or merges them.

An orientation code records how one cell-local face view maps to the shared face’s
`face2node` ordering. Each adjacent cell uses its own mapping; neither cell’s local
view is chosen as the master. Orienting does not make the two cells' outward normals
point in the same direction. Their outward normals on a shared face point in opposite
directions. It makes the node, edge, child, and split-direction references in the two
plans refer to the same physical parts of the face.

`hexOrientCode` supplies this mapping for a `HexCell`. The translation can
change a directional split code and can remap child and edge numbers. Prism
faces use the same kind of mapping, stored in `prismOrientCode`.

The examples below use HexCell face 0, the `xi = 1` face. The left side of each
figure shows its HexCell node and edge numbers. The projected copy shows one
possible QuadFace-local view of the same face. Its blue labels and arrows are
QuadFace-local node numbers, edge numbers, and edge directions.

<img src="figures/numbering/quadface_orientation_q0_h4.svg"
     alt="HexCell xi equals 1 face projected with QuadFace n_q zero at n_H four"
     width="820">

In the first example, QuadFace node 0 is HexCell node 4. The QuadFace and
HexCell views start at the same corner and proceed around the face in the same
direction.

In the next example, the physical face and the HexCell view are unchanged, but
the `face2node` list begins at another corner. QuadFace node 0 is now HexCell
node 6. The remaining QuadFace node and edge numbers shift with it.

<img src="figures/numbering/quadface_orientation_q0_h6.svg"
     alt="HexCell xi equals 1 face projected with QuadFace n_q zero at n_H six"
     width="820">

In the final example, QuadFace node 0 is again HexCell node 4, but the
`face2node` order runs in the opposite direction. The reversed order also
reverses the face normal shown in the figure.

<img src="figures/numbering/quadface_orientation_reversed.svg"
     alt="HexCell xi equals 1 face projected with reversed QuadFace local order"
     width="820">


## General Cells

The general `Cell` path does not have one fixed node, edge, and face table.
Its local ordering follows the cell's actual connectivity, which can vary with
the number and shape of its faces and with the number of edges at each vertex.
The fixed QuadFace, HexCell, and Prism tables above therefore should not be
applied to an arbitrary polyhedral cell.
