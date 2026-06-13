# Loci Documentation Sources

This directory contains editable documentation sources for the Loci repository.
Generated documentation should be written into the configured `OBJ` build tree
or an install prefix, not committed here.

## Layout

- `tutorial/`: Loci tutorial source files and companion example programs.
- `doxygen/`: Doxygen configuration and theme files for the API reference.
- `developer/`: Developer-facing notes that are not part of the tutorial.

After running `./configure`, use the normal build entry point:

```bash
make docs
```

The tutorial PDF is built under `OBJ/docs/tutorial/docs/`, and the Doxygen API
reference is built under `OBJ/docs/doxygen/html/` when Doxygen is available.

Release builders can include prebuilt documentation in a source archive with:

```bash
make tarball-with-docs
```

That target builds the docs in the configured `OBJ` tree, stages a temporary
copy of the source archive, overlays `tutorial.pdf` and the Doxygen HTML output
under `docs/`, and writes `docs/generated-docs-manifest.txt` into the archive.
The ordinary `make tarball` target remains source-only.
