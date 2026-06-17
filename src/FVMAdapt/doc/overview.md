# FVMAdapt Mesh Adaptation Overview {#fvmadapt_overview}

This page is intended to be an entry point for understanding the existing details of that governing the principles and implementation around the mesh adaptation library.

The principle here is that general polyhedral elements can localize the effects of mesh refinement.

## Geometric Concepts

The mesh adaptation library concerns itself with various types of geometric entities, such as
hexahedra, prisms, pyramids, diamonds, and general polyhedral cells. At the highest level, the concerns center around
how to break apart and recombine these shapes. These shapes are of course related to the mesh cell shapes
that occur in a general mesh.

The breakdown of the heirarchy that is at play in this module is:
  - element contains faces → face contains edges → edge contains nodes

A principle here is that cells are the entities that drive the refinement behaviors. The refinement algorithm utilizes a midpoint-based refinement strategy where edges are split by putting a node at the center of the edge, faces are split by putting a node at the middle of the face, and cells are split by putting a node at the center of the cell.

But in the formalism here, general polyhedral cells are supported, and thus adjacent cells can transform from a standard shape into a general polyhedral shape in order to satisfy a refinement that occured on a shared face between it and its neighbor cell. This general cell allows for refinements to be locally limited and not propagate widely throughout a mesh.


For standard canonical cell shapes, if they are split using the face-based isotropic strategy, they break into:

 - Hexahedra -> 8 smaller hexahedral cells

TODO: Place image here

 - Tetrahedra -> 4 hexahedral cells

TODO: Place image here

 - Prism cell -> 6 hexahedral cells

 - Pyriamid cell -> 4 hexahedral cells + 1 diamond cell

With this splitting behavior, there is a tendency towards hexahedral cells as isotropic refinement is applied.

## Face-Based Isotropic Refinement Strategy

For an isotropic face-based refinement strategy, the refinement of the faces determines how the volume elements are refined.
A n-sided face is split into n faces by insering a node a the centroid location of the face, and then connecting that node to the midpoints of the n edges of the face.

A neighborhood heuristic for if an element should be refined is if a neighboring element is two or more levels of refinement. Also, generic
marking of an element for refinment can also be done by a solver for any reason. For an element marked for refinement:

1. Each edge of the element is split by insering a node at the midpoint.
2. For each N-sided face of the element:
   - A node is inerted at the centroid of the face.
   - A node at the midpoint of each edge is connected to the node at the face centroid. This forms N new edges.
   - Each node on the face lies exactly on two edges and forms a four-sided face with the midpoints of the two edges and the centroid of the face so that N four-sided faces are formed.

Below are some diagrams showing this process on a set of simple faces.


** Pictures of the different splitting cases on the 2d faces. **


3. For the N-faced element, that has E edges and V nodes:
  - Inert a new node at the centroid of the element
  - Connect that node to the centroid nodes of all the faces of the element.
  - This will create V new elements from the original element


## Edge-Based Anisotropic Refinement Strategy

Edge-based refinement strategies are called this because the element type *and* the division of the edges determines how an element is split. This is still technically a face-based splitting strategy like the ones discussed above, with the exception that the way that the faces are split isn't a simple centroid approach.

The cell types that can be anisotropically refined are the prism and hexahedral type.


## What is Balancing?

As part of the refinement process, cells are initially marked to be refined, and once they are refined in the plan, there still remains an accounting that must be done to take into account several factors such as face and edge consistency between neighbor cells, and also other heuristic rules such as a neighbor not being permitted to be more than 1 level refined than its neighbors.


## Various Refinement Criteria

## Mesh Quality Concepts

The idea of the quality of a mesh cell is somewhat subjective as the quality of the cell depends on
where that cell is in a domain and the physics that is being solved on that cell. You can a pyramidal element that
can be broken down into tetrahedra and then a standard metric for quality tetrahedra can be applied to the resulting tets.

Another approach is one where there is an assumption that a volume element is a good quality if the faces are of good quality.
This isn't always true though.


## General Flow

Baseline mesh:
  elements are standard or general polyhedra

Sensor step:
  mark nodes/faces/elements based on solution features or error

Refinement decision:
  choose which elements to refine

Subdivision step:
  split selected elements using isotropic, centroidal, or anisotropic rules

Hanging-node handling:
  store parent/child edge and face relationships

Final topology update:
  replace parent faces/edges with child faces/edges
  unrefined neighbors become general elements if necessary

Solver:
  finite-volume fluxes are computed over the resulting set of faces


## Glossary of Terms

* Convex
* Concave
* Divergence Angle
* Isotropic
* Hanging Node

* Folded Face
  - A face is called a folded face if it is non-planar in 3D. Folded can cause issues with solvers, and as such care needs to be taken to not generate any folded faces during a refinement. A metric for marking if a face is a folded face is if the angle between the normals of two triangles on it is larger than a specified threshold. Folded faces often have small edge lengths, so a threshold on the minimum length edge that can be refined also serves to remedy the issue of folded faces.

* Refinement Level
  - A history of an element's refinement state. If an element has been refined many times, that information can be used by the refinement algorithm to refine its neighbors to prevent large cell size variations in a region where refinement is occuring.

* Element Neighbor
  - For an element, a neighbor is defined as two elements that share one or more nodes.

* Plan
  - A plan is a compact replay of data. It tells how to rebuild a refinement tree. It is not a mesh. It is a way to replay recipes for refinement trees.

* Cell Plan

* Edge Plan

* Face Plan

* Balancing

* Split Codes
  - A specific vocabulary used to describe the types of splitting to be used. This is not universal.

* Template
  - A specific approach to splitting cells that is determined by the cell type that is being split.






## Paper To Implementation Crosswalk

The papers are useful background, but the implementation is not a line-by-line
copy of either paper. This crosswalk lists the safest correspondences to use
when reading the code.

| Paper idea | Closest FVMAdapt implementation | Caution |
| --- | --- | --- |
| General elements allow faces and cells with non-template topology. | `Face`, `Cell`, and `DiamondCell` in `src/include/FVMAdapt` and `src/FVMAdapt/library`. | The code still has specialized `QuadFace`, `HexCell`, and `Prism` paths. Do not treat every path as a generic polyhedron path. |
| Face-based isotropic refinement splits an n-sided face through edge midpoints and a face centroid. | `Face::split()` splits parent edges, creates a face-center node, and builds one quadrilateral child face per parent edge. | Orientation and child ordering are implementation contracts and should be checked against extraction and output helpers before changing comments. |
| Cell subdivision inserts an element center and connects it with face and edge centers. | `Cell::split()` splits all faces, creates a cell-center node, builds interior faces, and creates `DiamondCell` children. | `DiamondCell` is the implementation vehicle for general-cell leaves; its child ordering is not explained by the paper alone. |
| Refinement state must survive across cycles. | `cellPlan`, `facePlan`, and `edgePlan` store breadth-first `std::vector<char>` replay recipes. | These vectors are Loci facts and library inputs. They are not the same thing as a full mesh-state file. |
| Local refinement still needs consistency checks. | `balance_general_cell.loci`, `balance_general_face.loci`, and `balance_general_edge.loci` update and merge plans across cells, faces, and edges. | "General elements localize refinement" should not be read as "no balancing is needed." |
| Anisotropic refinement can be driven by edge or direction choices. | Directional concepts appear in `QuadFace`, `HexCell`, `Prism`, and transfer/extraction helpers. | The mapping to the thesis' anisotropic formalism is only partially audited. Keep those code paths documented separately. |
| Solvers can integrate fluxes on general elements. | Outside the current FVMAdapt documentation scope. | This page documents remeshing/adaptation mechanics, not a solver-facing AMR contract. |



## The Central Model

The current implementation is easiest to read as a set of replayable refinement
trees.

An edge tree has a root edge. A split inserts a midpoint and creates two child
edges. A face tree has a root face. For a general polygonal `Face`,
`Face::split()` first splits the parent edges, then inserts a face-center node
and creates one quadrilateral child face per parent edge. A general cell uses
the same pattern one dimension higher: `Cell::split()` splits its faces,
inserts a cell-center node, creates interior faces from edge/face/cell centers,
and creates one child `DiamondCell` associated with each original cell node.

This is the closest match between the papers and the code: the paper-level
face-based isotropic refinement description appears in the library-level
midpoint/centroid construction used by `Face`, `Cell`, and `DiamondCell`.

The tree is not kept as a permanent object between every rule. Instead, rules
pass compact plans:

- `cellPlan` describes the cell refinement tree.
- `facePlan` describes the face refinement tree.
- `edgePlan` describes the edge refinement tree.

Most plans are `std::vector<char>` values consumed in breadth-first order.
Those vectors are recipes for replaying a tree, not the refined mesh itself.
See \ref fvmadapt_refinement_plans for the current source-backed code tables.

Suggested figure: `src/FVMAdapt/doc/figures/fvmadapt_plan_tree.svg` should show
one small tree, its breadth-first queue order, and the matching plan vector.

## General Elements And Locality

The general-element idea matters because it changes what adaptation has to
force. If every face mismatch must be removed immediately, refinement in one
cell can trigger a cascade of neighboring refinements. General elements allow
the code to represent a more local topology after refinement.

That does not mean the implementation ignores consistency. FVMAdapt still has
balancing rules and helper checks. For example, the general-cell balancing path
checks whether edge refinement depth differs by more than one level and can
split the cell to restore compatibility. Face and edge plans are also extracted
and merged so neighboring cells agree on the topology of a shared face or edge.

A useful way to read the rule database is:

1. Tags or region inputs mark where refinement is requested.
2. Cell plans are created or updated.
3. Cell plans are balanced.
4. Face plans are extracted from neighboring cell plans and merged on shared
   faces.
5. Edge plans are extracted from face plans and merged on shared edges.
6. Balanced plans are replayed to generate fine nodes, fine faces, fine cells,
   transfer data, and VOG output.

Suggested figure: `src/FVMAdapt/doc/figures/fvmadapt_plan_flow.svg` should show
this path from refinement request to balanced plans to generated VOG topology.

## Isotropic And Directional Paths

The papers discuss both isotropic and anisotropic refinement. The current code
contains several related, but not identical, formalisms:

- General cells use a face-based, centroid-style isotropic construction through
  `Face`, `Cell`, and `DiamondCell`.
- `QuadFace` has face-local directional split codes.
- `HexCell` has a three-bit local-direction split mask.
- `Prism` has prism-specific split codes and an `nfold` convention that still
  needs a focused documentation audit before the terminology is broadened.

Do not assume these code paths share one universal split-code meaning. A code
value that is directional for a `HexCell` is not automatically directional for
a general face or prism. The current documentation should keep those meanings
separate until each path has been checked against its split, extraction, merge,
and output helpers.

## What Is Stable Enough To Rely On

These concepts are stable enough to use as reading guides:

- Plans are replay recipes for refinement trees.
- General faces and general cells are first-class in the helper library.
- Balancing and lower-dimensional plan merging are central to preserving shared
  topology.
- The rule layer carries facts through the Loci scheduler; the library layer
  implements the geometry/tree mechanics.
- The smoke-test path exercises marker, remeshing, VOG validation, and sorted
  coordinate comparison.

These concepts should still be treated as open or partially documented:

- Exact orientation-code meanings across every helper.
- Exact child ordering for every prism and mixed general/prism/hex case.
- Which restart and transfer facts are stable external semantics rather than
  current internal plumbing.
- Which broad header comments are validated source documentation and which are
  still draft signposts.

## Where To Go Next

Read these pages in this order:

1. \ref fvmadapt_refinement_plans for the current plan-vector vocabulary.
2. \ref fvmadapt_concepts for rule translation notes and implementation terms.
3. \ref fvmadapt_testing for the current smoke-test path.
4. \ref fvmadapt_diagram_requests for figures that would make the tree and
   orientation conventions reviewable.
5. \ref fvmadapt_doxygen_status for the current generated-documentation status
   and known `.loci` parser limitations.
6. \ref fvmadapt_developer_reference for generated Doxygen groups and API
   browsing.

The development notes are useful, but they are not all user-facing
documentation. When a statement matters to users of this module, prefer to
promote it into this overview or \ref fvmadapt_refinement_plans after checking
it against the implementation.
