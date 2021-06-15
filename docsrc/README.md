# SeidarT Documentation

This document describes how to modify and build SeidarT documentation.

Documentation is currently hosted [here](https://sbernsen.github.io/SeidarT).

## Overview

There are two documentation folders in the root directory:

```
docs/       # <- this is the folder that github pages pays attention to (production)
docsrc/     # <- this is the editing and build directory
```

All edits to the documentation are made in `docsrc/`.


## Requirements

In order to build this documentation, you will need two pieces of Python software:

- Sphinx
- sphinx-rtd-theme

To install these, simply activate your SeidarT environment (or another environment with Python3) and use pip or conda to install the requirements:

```
conda activate SeidarT
pip install sphinx sphinx-rtd-theme rst2pdf
```


## Editing the source files

Sphinx will use all `.rst` files in the `docsrc/` directory when building. Editing those files will cause Sphinx to render changes in the resulting production html files.

It may help to have the [reStructuredText Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html) handy while editing.


## Build

To build the documentation after editing, first change to the `docsrc/` directory (`cd docsrc`).

There are three makefile commands to know:

- `make clean` readies the build area for a new documentation build. This is important because certain javascript elements do not get written on each build. If a change (especially one pertaining to linking between documents) isn't showing up correctly, then this command will clear out the build directory so that Sphinx knows to rebuild everything.

- `make html` will build the documentation under `docsrc/_build/html/`. Use this command to preview changes before sending them to production. *note: the `docsrc/_build/` directory is ignored by git, so no changes inside there will be committed*


## Review

Navigate to the `SeidarT/docsrc/_build/html/` folder in your file browser, then right click the html files and open them with your choice of web browser.

Internal links should be fully navigable.


## Stage for Production

- `make github` is the same as `make html` but will also copy the documentation to `docs/` in order to stage the changes for production.

Use this command, then `git add --all` and `git commit -m 'documentation update'` and finally `git push` to send the updated documentation to Github Pages.

Updated documentation should be available on Github Pages shortly after pushing.


## Troubleshooting

Usually, an error running `make html` or `make github` means that you are not in an environment that has Sphinx in it. See [Requirements](#requirements).

Sphinx 4.0 and above depreciates a function that rst2pdf relies on to build the PDF file. At this time it is necessary to downgrade to Sphinx 3.5.4 in order to build documentation in PDF format.
If you encounter an error stemming from rst2pdf, downgrading usually solves this.
