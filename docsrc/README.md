# SeidarT Documentation

This document describes how to modify and build SeidarT documentation.

## Overview

There are two documentation folders in the root directory:

```
docs/       # <- this is the folder that github pages pays attention to (production)
docsrc/     # <- this is the editing and staging directory
```

All edits to the documentation are made in `docsrc/`. 

# Requirements

In order to build this documentation, you will need two pieces of Python software:

- Sphinx
- sphinx-rtd-theme

To install these, simply activate your SeidarT environment (or another environment with Python3) and use pip or conda to install the requirements:

```
conda activate SeidarT
pip install sphinx sphinx-rtd-theme
```

# Editing the source files

Sphinx will use all `.rst` files in the `docsrc/` directory when building. Editing those files will cause Sphinx to render changes in the resulting production html files.

It may help to have the [reStructuredText Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html) handy while editing.

# Build

To build the documentation after editing, first change to the `docsrc/` directory (`cd docsrc`).

There are three makefile commands to know:

- `make clean` readies the build area for a new documentation build. This is important because certain javascript elements do not get written on each build. If a change (especially one pertaining to linking between documents) isn't showing up correctly, then this command will clear out the build directory so that Sphinx knows to rebuild everything.

- `make html` will build the documentation under `docsrc/_build/html/`. Use this command to preview changes before sending them to production. *note: the `docsrc/_build/` directory is ignored by git, so no changes inside there will be committed*

- `make github` same as `make html` but also copy the documentation to `docs/` to stage the changes for production. Use this command, then `git add --all` and `git commit -m 'documentation update'` and finally `git push` to send the updated documentation to Github Pages.

Updated documentation should be available on Github Pages shortly after pushing.


# Troubleshooting

Usually, an error running `make html` or `make github` means that you are not in an environment that has Sphinx in it. See [Requirements](#requirements).
