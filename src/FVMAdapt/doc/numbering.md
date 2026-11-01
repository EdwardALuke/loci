# Identities and Output Numbers {#fvmadapt_numbering}

The original grid supplies stable cell, face, edge and vertex IDs. Each shared
edge has one direction; each shared face has one canonical vertex order.
A cell's local orientation maps into that shared order.

Three kinds of numbers have different jobs:

| Number | Meaning |
| --- | --- |
| Original entity ID | Ownership and correspondence across adaptation and MPI migration |
| History node ID | A cell within one original cell's retained history |
| Output entity number | A node, face or cell in the current mesh |

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
