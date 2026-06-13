# FVMAdapt Documentation Improvement Plan

This plan scopes a documentation-focused pull request for the existing
`FVMAdapt` module and its supporting `fvmadaptfunc` library. The intent is to
make the current implementation easier for developers to navigate, review, and
maintain without changing behavior or proposing a new AMR API.

## Motivation

`FVMAdapt` contains robust but old mesh-adaptation code whose behavior is
encoded across Loci rules, C++ mesh element classes, and compact refinement
plans stored as `std::vector<char>`. Some internal values are effectively byte
codes: cell, face, edge, orientation, and tag codes are interpreted by several
tables and switch statements. The immediate need is better in-code
documentation and higher-level guidance for the framework that already exists,
not a redesign of adaptation semantics.

The PR should therefore be documentation-only except for comments. It should
not rename APIs, replace byte codes, rewrite plan storage, change refinement
behavior, or introduce design requirements for a separate AMR module.

## First PR Scope

1. Add a developer-facing FVMAdapt overview page.
   - Describe the module/library split from `src/FVMAdapt/Makefile`:
     `fvmadapt_m.so` for Loci rules and `libfvmadaptfunc` for C++ helper
     logic.
   - Map the main source areas: `.loci` rules, `library/*.cc`
     implementations, `include/FVMAdapt/*.h` public/internal headers, XML
     region docs, and `quickTest/FVMAdaptTest`.
   - Explain the main data flow: tags or region requests become refinement
     plans; plans are balanced across cells, faces, and edges; refined topology
     and transfer maps are generated.

2. Add Doxygen structure for the existing framework.
   The repository Doxygen setup already scans `*.cc`, `*.h`, `*.loci`, `*.lh`,
   and `*.md`, maps `.loci` as C++, and uses Javadoc autobrief comments. Use
   that existing machinery instead of inventing a parallel docs system.
   - Add `@defgroup fvmadapt FVMAdapt mesh adaptation`.
   - Add subgroups for plan representation, mesh element trees, Loci rules,
     balancing, topology generation, transfer/restart support, and XML region
     input.
   - Add `@file` comments to the FVMAdapt headers and key `.loci` rule files.
   - Add `@ingroup` annotations where they make the generated API reference
     easier to browse.
   - Add a short Doxygen page that links the Markdown guide, code groups, and
     existing XML region documentation.

3. Improve in-code comments and docstrings.
   - Document class responsibilities for `Node`, `Edge`, `Face`, `QuadFace`,
     `HexCell`, `Prism`, `DiamondCell`, and general-cell helpers.
   - Document ownership expectations for tree nodes, faces, edges, child
     arrays, and temporary lists because these are easy to misunderstand.
   - Document ordering assumptions for local cell faces, face edges, generated
     children, and fine-grid output.
   - Use `@warning` for invariants that must not be changed casually, such as
     shared-face/edge consistency and orientation-dependent ordering.
   - Keep comments factual and local. If behavior is unclear, mark it as an
     open documentation question rather than guessing.

4. Document refinement-plan representations.
   - `cellPlan`, `facePlan`, `edgePlan`, and node/fine-cell tag vectors.
   - Empty plan vectors versus explicit no-split codes.
   - Hex split codes `0..7` as local-direction bit codes, as currently
     described in `HexCell::mySplitCode`.
   - Prism split codes `0..3`, including the role of `nfold` for generated
     prism children.
   - Quad-face split codes `0..3`, including x-only, y-only, and x/y split
     behavior.
   - Special propagation code `8` in `tables.h`, currently used for face/edge
     extraction when a parent entity is split but a particular lower-dimensional
     entity is not split at that step.
   - Orientation-code conversion responsibilities, without claiming the numeric
     orientation code meanings until they are verified against the
     implementation.

5. Document existing execution and testing entry points.
   - Reference `quickTest/FVMAdaptTest/Makefile` as the current smoke test:
     `marker`, `refmesh`, `vogcheck`, and `extract` produce sorted coordinate
     comparisons against checked-in reference data.
   - Explain how the existing marker/refmesh path exercises the module at a
     high level.
   - Record test coverage limitations in neutral terms, without turning this
     PR into a test expansion or new AMR design effort.

## Suggested Documentation Files

- `src/FVMAdapt/doc/developer_guide.md`
  Main developer overview, data-flow map, invariants, and glossary.

- `src/FVMAdapt/doc/refinement_plans.md`
  Tables for currently implemented plan byte codes and how cell, face, and edge
  plans interact.

- `src/FVMAdapt/doc/doxygen_groups.md`
  Doxygen landing page for the FVMAdapt groups and higher-level guide links.

- `src/FVMAdapt/doc/testing.md`
  Existing quick-test coverage and how to run the current smoke test.

The files should be Markdown so Doxygen includes them in generated HTML.

## Review Strategy

Keep the first PR easy to review:

- No behavior changes.
- No generated Doxygen output committed unless the project explicitly wants
  generated API-reference artifacts in version control.
- No new AMR-module requirements or API-design use cases in this PR.
- Every numeric-code table should cite the source declaration, table, or switch
  statement used to derive it.
- Ambiguous behavior should be marked as an open question, not guessed.
- Use small Doxygen comments near dense code as signposts, but put longer
  explanations in Markdown pages.

## Follow-Up PRs

After the documentation PR is reviewed, later PRs can safely be scoped around
incremental maintainability work:

1. Add regression tests for documented current behavior and invariants.
2. Replace magic byte codes with named constants or enums where that can be done
   without changing serialized plan compatibility.
3. Add read-only query helpers for refinement state only after the current
   representation is documented well enough to avoid accidental semantic drift.
4. Refactor ownership and tree-walking code only behind passing tests that cover
   the affected cell types and shared-face/edge interactions.
