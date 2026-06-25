# FVMAdapt Mesh Adaptation Overview {#fvmadapt_overview}

This page is intended to be an entry point for understanding the existing details of the governing the principles and implementation around the mesh adaptation library.

A primary principle here is that general polyhedral elements can localize the effects of mesh refinement. The general element as the basis for the library matters because it gives freedom for the the adaptation to convert various element types to general cells in order to keep the effects of refinment localized. If every face mismatch must be removed immediately, refinement in one cell would trigger a cascade of neighboring refinements.

The idea of this module is to take an existing mesh, and create a refinement plan for how to refine the mesh based on a combination of externally applied markers and mesh-focused heuristic logic(like not allowing cells next to each other to be refined more than a set number of times). These drivers of the refinment must be coordinated and the ripple effects of the refinment need to be accounted for. A **Plan** is the thing that results from this process of taking these input and heuristics and combining them to come up with what the final mesh should look like. This plan isn't concerned with the details of the mesh itself in terms of coordinates. It is a topological plan that is concerned with the relationships and connections between entities such as nodes, edges, faces, and cells. The plan is firstly a cell-plan in the sense that it is concerned with what cells are to be split and what they are refined into. From a cell plan, other plans for lower level entities can be constructed, such as how faces are to be split, or edges split.

There is a general lower level C++ library that is concerned with cells and how they are split, and then there is a higher level Loci portion of the FVMAdapt module that concerned with orchestrating these per-cell splitting behaviors over an entire mesh and mapping the actual mesh details such as numbering and coordinates in order to generate a refined mesh.


### General Flow

Baseline mesh:
  elements are standard or general polyhedra

Sensor step:
  mark elements based on solution features or error

Refinement decision:
  choose which elements to refine

Subdivision step:
  split selected elements using isotropic or anisotropic rules

Hanging-node handling:
  store parent/child edge and face relationships

Final topology update:
  replace parent faces/edges with child faces/edges
  unrefined neighbors become general elements if necessary


## Geometric Concepts

The mesh adaptation library concerns itself with various types of geometric entities, such as
hexahedra, prisms, pyramids, diamonds, and general polyhedral cells. These concerns center around
how to break apart and recombine these shapes. These shapes are of course related to the mesh cell shapes
that occur in meshes that you would find in finite volume grids. The general polyhedral cell is the cornerstone of this library.

The breakdown of the heirarchy that is central is:
  - cell contains faces → face contains edges → edge contains nodes

A principle of the library is that cells are the entities that drive the refinement behaviors. The refinement algorithm utilizes a midpoint-based refinement strategy where edges are split by putting a node at the center of the edge, faces are split by putting a node at the middle of the face, and cells are split by putting a node at the center of the cell.

In this library, general polyhedral cells are supported, and cells can become polyhedral cells during the course of a refinement. This allows adjacent cells to transform from a standard shape into a general polyhedral shape in order to satisfy a refinement that occured on a shared face between the two cells. Once again, this general cell allows for refinements to be locally limited and not propagate widely throughout a mesh.

<div class="term-box">
<div class="term-title">Splitting</div>
Mutate the current tree node by creating geometric children. A split usually creates new nodes, edges, faces, or cells.
</div>

The detailed enumeration of isotropic midpoint split behavior and hexahedron anisotropic split codes is kept on a separate reference page: [FVMAdapt Geometric Splitting Reference](@ref fvmadapt_geometric_splitting). The local node, edge, and face numbering used by the refinement library is kept in [FVMAdapt Numbering Conventions](@ref fvmadapt_numbering).


## The Model For the Library

The lower level library can be thought of as a set of replayable refinement trees. These refinement trees describe how cells, faces, and edges in the mesh topology are split and how many times entities are split.

An edge tree has a root edge. A split inserts a midpoint and creates two child edges. A face tree has a root face. For a general polygonal Face, and isotropic split first splits the parent edges, then inserts a face-center node and creates one quadrilateral child face per parent edge. A general cell has a root cell. The isotropic split first splits its faces, inserts a cell-center node, creates interior faces from edge/face/cell centers, and creates one child `DiamondCell` associated with each original cell node.

The tree is not kept as a permanent object. Instead, compact plans that describe the refinement are passed around

- A cell plan describes the cell refinement tree.
- A face plan describes the face refinement tree.
- An edge plan describes the edge refinement tree.

Most plans are `std::vector<char>` values consumed in breadth-first order. Those vectors are recipes for replaying a tree, not the refined mesh itself. Face and edge plans are constructed and used to ensure that neighboring cells agree on the topology of a shared face or edge.

A useful way to think about the general flow of this library is:

1. Tags or region inputs mark where refinement is requested.
2. Cell plans are created or updated.
3. Cell plans are balanced.
4. Face plans are extracted from neighboring cell plans and merged on shared faces.
5. Edge plans are extracted from face plans and merged on shared edges.
6. Balanced plans are replayed to generate fine nodes, fine faces, fine cells, transfer data, and VOG output.

If we just consider a single cell, the cell plan for that cell is a succint summary of the process that would be involved in that cell for example splitting once, and then one of its children splitting and the other one doesn't split. This tree of the actual objects representing the cell and its children is the "real" thing, and the plan is just a compact summary of the operations that happened to create this particular tree for that particular cell. The collection of all of these plans for all cells is what represents the data that is used to perform the balancing and final construction of the refined mesh.

## What is Balancing?

As part of the refinement process, cells are initially marked to be refined, and once they are refined in the plan, there still remains an accounting that must be done to take into account several factors such as face and edge consistency between neighbor cells, and also other heuristic rules such as a neighbor not being permitted to be more than 1 level refined than its neighbors. This is an iterative process which converges to a set of cell, face, and edge plans that do not change with further iterations. This convergence represents the state when all of the requirements imposed on the mesh refinement have been satisfied.

Balancing is the set of operations that must be performed such that there is agreement between a requested refinement and the mesh topology that must result from it. During this stage, an existing cell plan can be adjusted

A refinement request first becomes a cell plan: a breadth-first vector of split codes describing how to replay one cell’s refinement tree. Balancing adjusts those requested plans so neighboring cells remain topologically compatible. The balanced cell plans then imply face plans on cell faces; neighboring face plans are merged into shared face plans; edge plans are extracted and merged from those faces. The result is not the fine mesh yet, but a compatible set of balanced cell, face, and edge plans that can then be replayed to generate a final refined mesh. The breadth-first nomenclature means that all cells are considered at a certain refinement level before an algorithm moves down to the next refinement level.

Balancing and lower-dimensional plan merging are central to preserving shared topology in the refined grid. For example, the general-cell balancing path checks whether edge refinement depth differs by more than one level and can split the cell to restore compatibility. Face and edge plans are also extracted and merged so neighboring cells agree on the topology of a shared face or edge. At the end of the process all entities that are shared within a grid should agree on what they look like. One cell can not have a 3 sided face and the neighbor think that the face has 4 sides for example.

In balancing, a requested plan is expanded until adjacent cells, faces, and edges can be replayed consistently. Balancing does not directly
build the fine mesh; it produces compatible plan topology that later construction, numbering, transfer, and restart rules consume.

Plan consistency is maintained across dimensions:

- Cell plans imply face plans on each cell face.
- Face plans from neighboring cells are merged into a shared face plan.
- Edge plans are extracted from face plans.
- Later rules generate fine nodes, fine faces, cell numbering, and transfer data from the balanced plans.

Plan derivation is the lower-dimensional part of that process. Cell plans are used to derive face plans, and face plans are used to derive edge plans or intermediate edge split-coordinate sets.  One of the difficult parts of balancing is not only whether an entity is split, but whether the ordering of generated children, faces, and nodes remains consistent between neighboring cells and across restart or derefinement paths. This consistency is what the final result of the balancing provides.


## The Cell Plan

The library mostly operates on one entity tree at a time. A root object represents the original cell, a split creates child cells under the root. Child cells are linked back to
their parent. A leaf is a cell that has not been split further. Fine cells are generated by walking the leaves and child
relationships. Plan vectors are only replay inputs or outputs for these temporary trees.

For a cell, a cell plan would be a log of the splitting information for that cell and the children and those children's children, etc.

Consider one cell object: cell

After one split:

```text
cell
├─ child 0
├─ child 1
├─ child 2
└─ child 3
```

If child 1 is split again:

```text
cell
├─ child 0
├─ child 1
│  ├─ child 1.0
│  └─ child 1.1
├─ child 2
└─ child 3
```

A plan vector is not the tree itself. It is a compact input that can be used used to recreate the tree. Creating a plan walks an existing tree and writes the split codes into a vector. Each entity type defines their own split codes.
