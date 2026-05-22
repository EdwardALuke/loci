# Doxygen-based documentation

This directory contains the Doxygen configuration and HTML theme files for the
Loci API reference documentation.

## Building the documentation

### Prerequisites

- Doxygen 1.9.8 or later
- Graphviz (`dot` command) — required for collaboration and include graphs

### Steps

Run Doxygen from **this directory** so that relative paths in the `Doxyfile`
resolve correctly:

```bash
cd doxygen
doxygen Doxyfile
```

The generated HTML lands in `docs/html/`. Open `docs/html/index.html`
in a browser to browse the output.

### Cleaning the output

```bash
rm -rf ../docs/html
```

# Build Configuration

## Files

| File             | Purpose                                               |
|------------------|-------------------------------------------------------|
| `Doxyfile`       | Main Doxygen configuration                            |
| `header.html`    | Custom HTML header injected into every generated page |
| `footer.html`    | Custom HTML footer injected into every generated page |
| `stylesheet.css` | Extra CSS applied on top of the Doxygen default theme |

## Build configuration summary

| Setting                  | Value                                                                    |
|--------------------------|--------------------------------------------------------------------------|
| Project                  | Loci                                                                     |
| Input root               | `..` (repository root)                                                   |
| File patterns            | `*.cc`, `*.h`, `*.loci`, `*.lh`, `*.md`                                  |
| Recursive                | Yes                                                                      |
| Excluded directories     | `ext`, `contrib`, `docs`, `tmpcopy`                                       |
| Output directory         | `../docs`                                                                |
| HTML output subdirectory | `../docs/html`                                                           |
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
