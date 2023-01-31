# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = source
BUILDDIR      = _build
FILESDIR      = source/_files

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@echo "  pdfdocs     to delete $(FILESDIR) folder, generate pdf and save them into $(FILESDIR)"

.PHONY: help Makefile

# Get and adapt toulbar2/README.md :
# - delete lines 'master:...' and 'cpd:...' that are not ok for LaTeX
# - change some titles
readme:
	cat ../README.md > $(FILESDIR)/README.md
	sed -i '/^master:/d' $(FILESDIR)/README.md
	sed -i '/^cpd:/d' $(FILESDIR)/README.md
	sed -i '1,10 s/# toulbar2//' $(FILESDIR)/README.md
	sed -i '1,10 s/## Exact optimization/# Exact optimization/' $(FILESDIR)/README.md
	sed -i 's/<!-- \(.*\) -->/\1/g' $(FILESDIR)/README.md

# the two first lines (with svg's) are removed (build status and deployment)
readme_latex:
	make readme
	sed -i '/Build Status/d' $(FILESDIR)/README.md
	sed -i '/PyPi version/d' $(FILESDIR)/README.md

readme_html:
	rm -f $(FILESDIR)/README.md
	make readme

# build pdf files, get and adapt toulbar2/README.md
# and install them into FILESDIR
# Note : has to be done before html pages
files:
	rm -f $(FILESDIR)/*.*
	make readme_latex
	make latexpdf
	cp $(BUILDDIR)/latex/*.pdf $(FILESDIR)/.

# builds pdf, html files; pdf into html/pdf
docs:
	make clean
	make files
	make readme_html
	make html
	make files
	mkdir $(BUILDDIR)/html/pdf
	cp $(FILESDIR)/*.pdf $(BUILDDIR)/html/pdf/.

# epub not ok => not done :
# builds epub file; epub into html/epub
#   make epub
#   mkdir $(BUILDDIR)/html/epub
#   cp $(BUILDDIR)/epub/toulbar2.epub $(BUILDDIR)/html/epub/.

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
