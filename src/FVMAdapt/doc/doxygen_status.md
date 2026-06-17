# FVMAdapt Doxygen Status {#fvmadapt_doxygen_status}

This page records the current documentation-generation status for the existing
`FVMAdapt` implementation. It is meant to keep review conversations concrete:
what Doxygen can show today, what it cannot parse cleanly, and how to validate
future documentation edits.

## Current Wiring

The repository has two Doxygen configuration files:

- `Doxyfile` writes generated output under `doc`.
- `doxygen/Doxyfile` writes generated output under `docs`.

Both configurations currently:

- enable Markdown support;
- map `.loci` and `.lh` files through the C++ parser with
  `EXTENSION_MAPPING = loci=C++ lh=C++`;
- include C++ sources, headers, Loci rule files, Loci headers, and Markdown
  pages in the scanned file patterns;
- use `README.md` as the project main page.

For FVMAdapt specifically, `src/FVMAdapt/doc/doxygen_groups.md` defines the
current topic groups and links the higher-level Markdown pages. The intended
human entry point is \ref fvmadapt_overview.

## What Works Today

Doxygen is useful for:

- rendering the FVMAdapt Markdown pages into the generated documentation;
- grouping helper classes and functions under broad FVMAdapt topics;
- browsing source for the C++ helper library and `.loci` rule files;
- checking that cross-page references such as \ref fvmadapt_refinement_plans
  and \ref fvmadapt_testing resolve.

## Known Limitations

The `.loci` and `.lh` files are not C++. Mapping them through the C++ parser is
still useful for source browsing, but it creates parser limitations.

The current known warning pattern is for raw Loci `$type` declarations with
template-valued stores, for example `store<std::vector<char> >`. Doxygen reports
these as "documented symbol" warnings even when the nearby comments are ordinary
source comments rather than Doxygen comments. Treat those warnings as a parser
limitation, not as evidence that a `$type` declaration has a reliable generated
API entry.

Until there is a Loci-aware Doxygen filter or parser strategy:

- keep `$type` explanations as ordinary source comments or Markdown prose;
- do not rely on Doxygen to document raw `$type` declarations as API symbols;
- prefer Markdown pages for rule-database concepts, plan vocabulary, and
  source-backed crosswalks;
- reserve Doxygen comments in `.loci` files for rule-level signposts that have
  been checked against the implementation.

## Validation Command

For documentation edits, build into a temporary directory so generated HTML does
not dirty the worktree:

```bash
tmpdir=$(mktemp -d /tmp/loci-fvmadapt-doxy.XXXXXX)
awk '
  /^OUTPUT_DIRECTORY[[:space:]]*=/ {print "OUTPUT_DIRECTORY = " out; next}
  /^INPUT[[:space:]]*=/ {print "INPUT = src/FVMAdapt src/include/FVMAdapt"; next}
  /^USE_MDFILE_AS_MAINPAGE[[:space:]]*=/ {print "USE_MDFILE_AS_MAINPAGE = src/FVMAdapt/doc/overview.md"; next}
  /^QUIET[[:space:]]*=/ {print "QUIET = YES"; next}
  {print}
' out="$tmpdir" Doxyfile > "$tmpdir/Doxyfile"
doxygen "$tmpdir/Doxyfile"
```

A clean-enough run for the current branch is one that completes successfully and
does not introduce new warnings beyond the known raw `$type` parser warnings.
