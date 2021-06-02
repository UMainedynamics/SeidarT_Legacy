# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
CPDIR		  = ../docs

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

pdf:
	@mkdir -p "$(BUILDDIR)/pdf"
	@$(SPHINXBUILD) -b pdf "$(SOURCEDIR)" "$(BUILDDIR)/pdf"
	@echo "Build finished. The PDF files are in $(BUILDDIR)/pdf."

github:
	@make html
	@cp -a "$(BUILDDIR)/html/." "$(CPDIR)"
	@make pdf
	@cp -a "$(BUILDDIR)/pdf/." "$(CPDIR)"


# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)