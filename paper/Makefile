# Makefile for the structural-completeness paper.
#
# Targets:
#   make             — build PDF
#   make clean       — remove build artefacts
#   make arxiv       — produce arxiv.tar.gz submission bundle
#   make eprint      — produce eprint.zip submission bundle (alternate format)
#   make deanonymize — restore author block from AUTHORS file (camera-ready)

TEX     := structural_completeness.tex
PDF     := $(TEX:.tex=.pdf)
AUX     := $(TEX:.tex=.aux) $(TEX:.tex=.log) $(TEX:.tex=.out) $(TEX:.tex=.toc) $(TEX:.tex=.bbl) $(TEX:.tex=.blg)
ARXIV_TARBALL := arxiv-submission.tar.gz
EPRINT_ZIP    := eprint-submission.zip

.PHONY: all clean arxiv eprint deanonymize check

all: $(PDF)

$(PDF): $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)  # second pass for refs

clean:
	rm -f $(PDF) $(AUX)
	rm -f $(ARXIV_TARBALL) $(EPRINT_ZIP)
	rm -rf submission-staging/

# arXiv expects a single .tar.gz with .tex source(s) and any included
# figures.  No PDF (arXiv compiles it server-side).
arxiv: $(TEX)
	mkdir -p submission-staging/arxiv
	cp $(TEX) submission-staging/arxiv/
	cp abstract.txt submission-staging/arxiv/ 2>/dev/null || true
	cd submission-staging/arxiv && tar czf ../../$(ARXIV_TARBALL) *
	rm -rf submission-staging/arxiv
	@echo "arXiv submission bundle ready: $(ARXIV_TARBALL)"

# IACR ePrint accepts PDF + optional source.  Produce a zip with both.
eprint: $(PDF)
	mkdir -p submission-staging/eprint
	cp $(PDF) submission-staging/eprint/
	cp $(TEX) submission-staging/eprint/
	cp abstract.txt submission-staging/eprint/ 2>/dev/null || true
	cp README.md submission-staging/eprint/ 2>/dev/null || true
	cd submission-staging/eprint && zip -q ../../$(EPRINT_ZIP) *
	rm -rf submission-staging/eprint
	@echo "IACR ePrint submission bundle ready: $(EPRINT_ZIP)"

# Restore author block from AUTHORS file (for non-anonymous releases).
# AUTHORS file should contain LaTeX-formatted \author{...} block.
deanonymize: $(TEX)
	@if [ ! -f AUTHORS ]; then \
	  echo "ERROR: AUTHORS file missing — create it with the author block"; \
	  echo "Example AUTHORS file contents:"; \
	  echo '  \author{Alice \and Bob \\\\ Affiliation}'; \
	  echo '  \date{\today}'; \
	  exit 1; \
	fi
	@cp $(TEX) $(TEX).anon-backup
	@python3 -c 'import re, sys; \
	tex = open("$(TEX)").read(); \
	authors = open("AUTHORS").read(); \
	tex = re.sub(r"\\author\{Anonymous\}\s*\\date\{\}", authors, tex); \
	open("$(TEX)", "w").write(tex); \
	print("De-anonymized $(TEX); backup saved as $(TEX).anon-backup")'

# Sanity check: ensure the paper compiles and key cross-references resolve
check: $(PDF)
	@echo "Checking for unresolved references in $(PDF)..."
	@if grep -q '\?\?' $(PDF) 2>/dev/null; then \
	  echo "WARNING: unresolved references found"; \
	else \
	  echo "All references resolved."; \
	fi
	@echo "Page count:"
	@pdfinfo $(PDF) 2>/dev/null | grep -i pages || echo "  (install pdfinfo to see page count)"
