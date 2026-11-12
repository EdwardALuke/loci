# FVMAdapt Numbering Conventions {#fvmadapt_numbering}

By the time FVMAdapt operates, the mesh topology is represented by nodes,
edges, faces, and cells, together with the connections between them. A cell is
bounded by faces, a face is bounded by edges, and an edge has two endpoint
nodes. All four are mesh entities. In the original mesh topology, an interior
face is one shared face incident to two cells, not one face belonging to each
cell.

An entity's identity is different from its position in a face or cell's local
topology. Suppose an original mesh node `n` lies on a quadrilateral face and a
hex. The face may use `n` as its local node 0 while the hex uses the
same `n` as its local node 4. The cell on the other side of the face may assign
it another cell-local position. These local numbers describe the role that `n`
plays in each view; they do not identify different mesh nodes.

| Kind of number | Example | What it identifies |
| --- | --- | --- |
| Mesh entity identity | original mesh node `n` | Which mesh entity is being referenced |
| Local position | quad node 0 or hex node 4 | The entity's role in one local topology |
| History node ID | node ID 5 in one original cell's history | A cell within that retained history |
| Fine-mesh output number | a numbered fine node | Its number after the adapted mesh is constructed |

Shared entities provide the common references needed to compare these local
descriptions. The original vertex IDs determine one common order for each face
and one direction for each edge. A cell also has local positions for its faces,
edges, and nodes, but those positions are derived separately for that cell.

## Mapping a Face into a Cell-Local View

Two views of the same face can differ in two ways:

- They can choose different corners as local node 0.
- They can walk around the corners in opposite directions.

Changing the starting corner shifts the local numbers around the face.
Changing the direction reverses their order. Neither change moves or rotates
the physical face.

For example, these orders describe the same four vertices:

```text
shared face:       A B C D
one cell's view:   C B A D
```

To *orient* a face subdivision means to translate it from one of these local views
into another. This is necessary on a shared face because, for example, a split
in one cell's face-local `x` direction may be a split in the other cell's
face-local `y` direction. FVMAdapt maps both cell contributions into the
shared face's common orientation before it compares or merges them.
Their outward normals remain opposite; their descriptions must identify the
same points and pieces of the face.

## Constructing the Output Mesh

A new point belongs to the smallest original entity containing it: vertex,
edge, face or cell interior. Shared reference coordinates or accepted point
ordinals identify it within that entity. Both neighbors therefore name the
same point, despite different local cell orientations.

Output numbers identify points, faces and leaf cells in the current mesh.
Hanging edge points needed by finer neighbors are included. Physical
coordinates come from the saved geometry.

Old/new correspondence is separate from output numbering; see
[mesh correspondence](@ref fvmadapt_mesh_correspondence).

The retained history belongs to the original grid. An adapted VOG contains the
current output mesh; importing it as a new base grid does not recover that
ancestry.
