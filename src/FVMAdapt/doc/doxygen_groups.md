# FVMAdapt Developer Reference {#fvmadapt_developer_reference}

This page is the generated-reference landing page for the existing
`FVMAdapt` implementation. It links the higher-level developer notes to the
Doxygen topic groups used for classes, helper functions, and `.loci` rules.

\defgroup fvmadapt FVMAdapt mesh adaptation
\brief Existing Loci module and helper library for mesh refinement and
coarsening.

`FVMAdapt` is an implementation-facing adaptation subsystem. It combines Loci
rules, compiled helper classes, compact plan vectors, and VOG output utilities.
These pages document the current implementation; they are not a proposed new
AMR API.

Start with:

- \ref fvmadapt_overview for the user-facing entry point and guiding
  principles.
- \ref fvmadapt_concepts for terminology and the high-level data flow.
- \ref fvmadapt_refinement_plans for byte-code plan conventions.
- \ref fvmadapt_diagram_requests for figures that would clarify numbering,
  splitting, and orientation conventions.
- \ref fvmadapt_testing for the existing smoke-test path.
- \ref fvmadapt_doxygen_status for the current generated-documentation wiring
  and known `.loci` parser limitations.
- `src/FVMAdapt/doc/README_xml` for XML region input accepted by selected
  marker workflows.

The generated topic groups are intentionally broad. They are meant to make the
legacy implementation easier to browse while more precise comments are added
to individual rules, classes, and helper functions.

\defgroup fvmadapt_plans Refinement plan representation
\ingroup fvmadapt
\brief `std::vector<char>` plan encodings for cells, faces, and edges.

Plan vectors are interpreted by cell-type-specific and face-type-specific tree
walks. See \ref fvmadapt_refinement_plans before changing numeric plan codes or
removing apparently redundant entries.

\defgroup fvmadapt_elements Mesh element helper trees
\ingroup fvmadapt
\brief C++ helper classes that model refinement trees for nodes, edges, faces,
and cells.

This group covers `Node`, `Edge`, `Face`, `QuadFace`, `HexCell`, `Prism`,
`DiamondCell`, and the general-cell helper `Cell`.

\defgroup fvmadapt_loci_rules Loci adaptation rules
\ingroup fvmadapt
\brief `.loci` rules that create, balance, merge, and consume adaptation facts.

The rules move between tags, temporary plans, balanced plans, node/face/cell
numbering, and output maps.

\defgroup fvmadapt_balancing Balancing and consistency
\ingroup fvmadapt
\brief Helpers and rules that keep neighboring cell, face, and edge plans
compatible.

Balancing is dimension-coupled: cell plans imply face plans, shared face plans
imply edge plans, and child ordering must remain consistent across neighboring
cells.

\defgroup fvmadapt_topology Topology generation
\ingroup fvmadapt
\brief Builders for fine-grid nodes, faces, cells, and VOG output topology.

\defgroup fvmadapt_transfer Transfer and restart support
\ingroup fvmadapt
\brief Helpers that transfer plan and field data across refinement or restart
paths.

\defgroup fvmadapt_input Adaptation input
\ingroup fvmadapt
\brief Tag, source, and XML-region input paths used to seed adaptation plans.

\defgroup fvmadapt_tests Existing tests
\ingroup fvmadapt
\brief Current quick-test coverage and smoke-test commands.
