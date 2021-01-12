# SeidarT Documentation

This document describes how to modify and build SeidarT documentation.

## Overview

There are two documentation folders in the root directory:

```
docs/       # <- this is the folder that github pages pays attention to (production)
docsrc/     # <- this is the editing and staging directory
```

`docs/` 

# Requirements

In order to build this documentation, you will need two pieces of Python software:

- Sphinx
- sphinx-rtd-theme

To install these, simply activate your SeidarT environment (or another environment with Python3) and use pip to install the requirements:

```
conda activate SeidarT
pip install sphinx sphinx-rtd-theme
```

# Editing the source files

Sphinx will build all `.rst` files in the `docsrc/` directory. Editing those files will cause sphinx to render changes in the resulting production html files.

# Build

To build the documentation after editing, first change to the `docsrc/` directory (`cd docsrc`).

There are three makefile commands to know:

- `make clean` readies the build area for a new documentation build. This is important because certain javascript elements do not get written on each build. If a change (especially one pertaining to linking between documents) isn't showing up correctly, then this command will clear out the build directory so that Sphinx knows to rebuild everything.

- `make html` will build the documentation under `docsrc/_build/html/`. Use this to preview changes before sending to production.

- `make github` will both build and copy the documentation to `docs/` to stage the changes for production. Use this command, then `git add --all` and `git commit -m 'documentation update'` and finally `git push` to send the updated documentation to Github Pages.
