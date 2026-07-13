# Doxygen-based documentation

This directory contains the Doxygen configuration and HTML theme files for the
Loci API reference documentation. The files here are source files; generated
HTML is written into the configured `OBJ` build tree.

## Building the documentation

### Prerequisites

- Doxygen 1.9.8 or later
- Graphviz (`dot` command), optional. Enable `HAVE_DOT` in `Doxyfile` if you
  want Doxygen's Graphviz-backed include, inheritance, or collaboration graphs.

### Steps

Run the repository docs target after configuring Loci:

```bash
make docs
```

The generated API reference lands in `OBJ/docs/doxygen/html/`. Open
`OBJ/docs/doxygen/html/index.html` in a browser to browse the output.

For a quick source-tree preview without using `OBJ`, run Doxygen from this
directory:

```bash
cd docs/doxygen
doxygen Doxyfile
```

That writes ignored preview output to `docs/doxygen/html/`.

### Cleaning the output

```bash
make distclean
```

`make clean` preserves Doxygen HTML as a final documentation product.
From the repository top level, `make distclean` removes generated Doxygen output
from both the configured `OBJ` tree and from a source tree unpacked from a
`tarball-with-docs` archive.

# Build Configuration

## Files

| File             | Purpose                                               |
|------------------|-------------------------------------------------------|
| `Doxyfile`       | Main Doxygen configuration                            |
| `theme/header.html`    | Custom HTML header injected into every generated page |
| `theme/footer.html`    | Custom HTML footer injected into every generated page |
| `theme/stylesheet.css` | Extra CSS applied on top of the Doxygen default theme |

## Build configuration summary

| Setting                  | Value                                                                    |
|--------------------------|--------------------------------------------------------------------------|
| Project                  | Loci                                                                     |
| Input root               | `../..` from `OBJ/docs/doxygen`, i.e. the configured build tree           |
| File patterns            | `*.cc`, `*.h`, `*.loci`, `*.lh`, `*.md`                                  |
| Recursive                | Yes                                                                      |
| Excluded directories     | `OBJ`, `ext`, `contrib`, `tmpcopy`, and generated Doxygen output          |
| Output directory         | `.` from `OBJ/docs/doxygen`                                               |
| HTML output subdirectory | `OBJ/docs/doxygen/html`                                                   |
| Formats generated        | HTML only (LaTeX, Man, XML disabled)                                     |
| Search engine            | Enabled (client-side)                                                    |
| Extract private members  | Yes                                                                      |
| Extract static members   | Yes                                                                      |
| JavaDoc autobrief        | Yes — first sentence of a `/** */` comment becomes the brief description |
| Preprocessing            | Enabled; macro expansion disabled                                        |

## Writing documentation

### Comment style

The configuration has `JAVADOC_AUTOBRIEF = YES`, so the first sentence of any
`/** */` block is automatically used as the brief description — no explicit
`@brief` tag is needed:

```text
/**
 * This sentence will go in the brief description. Additional text will
 * appear in the detailed description section.
 *
 * @param path  Filesystem path to test.
 * @return true when @p path is a directory.
 */
bool isDirectory(const string &path);
```

Use `///` or `//!` for single-line member documentation:

```text
int rule_count; ///< total number of active rules in the schedule
```

or as an alternative to the multi-line block `/** */` comment shown above:

```text
/// @brief Returns true when @p path exists and is a directory.
///
/// @param path  Filesystem path to test.
/// @return True when @p path is a directory.
bool isDirectory(const string &path);
```
### Documenting `.loci` files

`.loci` source files are included in the input scan. Doxygen treats them as
C-family source, so standard `/** */` and `///` comments work. 

No `.loci` files currently carry a file-level comment. Adding one
provides a visible entry in the file index.

```text
/**
 * @file inviscid_fluxes.loci
 * @brief Rules required for computing invisicid fluxes.
 *
 * This file include rules for computing ...
 */
```

### Grouping related symbols

No `@defgroup` or `@ingroup` annotations currently exist in the
codebase. Adding them would let Doxygen cluster related rules, types,
and functions into named modules regardless of which file they live
in. This is could be potentially useful as the coverage of the
documentation grows.

```text
/**
 * @defgroup fvm_tools FVM Mesh Utilities
 * @brief Tools for finite-volume mesh construction and interrogation.
 */

/** @ingroup fvm_tools */
$rule pointwise(face_area <- pos, face2node) { ... }
```

### Marking work in progress

| Tag           | Meaning                                     |
|---------------|---------------------------------------------|
| `@todo`       | Known improvement needed                    |
| `@bug`        | Known defect                                |
| `@warning`    | Non-obvious constraint callers must respect |
| `@note`       | Important but non-critical remark           |
| `@deprecated` | Symbol scheduled for removal                |

Doxygen collects all `@todo` and `@bug` entries into dedicated index pages,
making them easy to find and triage.
