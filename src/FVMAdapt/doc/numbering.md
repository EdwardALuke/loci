# Identities and Output Numbers {#fvmadapt_numbering}

The original grid supplies stable cell, face, edge and vertex IDs. Each shared
edge has one direction; each shared face has one canonical vertex order.
A cell's local orientation maps into that shared order.

An entity's identity is different from its position in a face or cell's local
topology. The same mesh node can be local node 0 on a face and local node 4 in
a cell. The cell on the other side may assign it another local position.
These numbers describe its role in each view; they do not identify different
mesh nodes.

The numbers have different jobs:

| Number | Meaning |
| --- | --- |
| Original entity ID | Ownership and correspondence across adaptation and MPI migration |
| Local position | A vertex, edge or face's role within one cell or face |
| History node ID | A cell within one original cell's retained history |
| Output entity number | A node, face or cell in the current mesh |

## Mapping a Face into a Cell-Local View

Two views of the same face can differ in two ways:

- They can choose different corners as local node 0.
- They can walk around the corners in opposite directions.

Changing the starting corner shifts the local numbers around the face.
Changing the direction reverses their order. Neither change moves or rotates
the physical face. For example, these orders describe the same four vertices:

```text
shared face:       A B C D
one cell's view:   C B A D
```

Before comparing subdivisions, FVMAdapt maps both cell contributions into the
shared face's order. A cut along one cell's first local direction may be a cut
along the other's second direction. Their outward normals remain opposite;
their descriptions must identify the same points and pieces of the face.

## Constructing the Output Mesh

A new point belongs to the smallest original entity containing it: vertex,
edge, face or cell interior. Shared reference coordinates or accepted point
ordinals identify it within that entity. Both neighbors therefore name the
same point, despite different local cell orientations.

The output rules count and number owned points and leaves, then build face
connectivity. Hanging edge points needed by finer neighbors are included.
Physical coordinates are evaluated from the saved geometry.

Old/new cell correspondence uses accepted ancestry. A refined child has one
old ancestor; a merged parent has several old donors. Face correspondence also
retains orientation and area fractions. These relationships are separate from
the output numbering and from solver field reconstruction.

Restart must retain the original grid and its history. Re-importing an adapted
VOG creates a new base grid and does not recover the original ancestry.
