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
