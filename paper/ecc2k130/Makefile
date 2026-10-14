# Makefile for the ECC2K-130 campaign paper.
#
# Targets:
#   make             — build PDF (pdflatex + bibtex + two more passes)
#   make check       — unresolved refs, overfull boxes, page count
#   make arxiv       — produce arxiv-submission.tar.gz
#   make eprint      — produce eprint-submission.zip
#   make deanonymize — restore author block from AUTHORS
#   make clean       — remove build artefacts (keeps the committed PDF)

TEX     := ecc2k130.tex
PDF     := $(TEX:.tex=.pdf)
SECTIONS := sec-intro.tex sec-history.tex sec-curve.tex sec-walk.tex sec-arith.tex \
            sec-isa.tex sec-blackwell.tex sec-gpu.tex sec-cpu.tex sec-nist.tex sec-othergpu.tex sec-fpga.tex sec-campaign.tex \
            sec-bounds.tex sec-solvetime.tex sec-open.tex sec-related.tex sec-conclusion.tex app-reproduction.tex app-evidence.tex
BIB     := ecc2k130.bib
AUX     := $(TEX:.tex=.aux) $(TEX:.tex=.log) $(TEX:.tex=.out) $(TEX:.tex=.toc) \
           $(TEX:.tex=.bbl) $(TEX:.tex=.blg) $(TEX:.tex=.lof) $(TEX:.tex=.lot)
ARXIV_TARBALL := arxiv-submission.tar.gz
EPRINT_ZIP    := eprint-submission.zip

.PHONY: all clean arxiv eprint deanonymize check figures-export

all: $(PDF)

FIGURES := $(wildcard figures/*.tex)

$(PDF): $(TEX) $(SECTIONS) $(FIGURES) $(BIB) abstract.txt
	pdflatex -interaction=nonstopmode $(TEX) || true
	bibtex $(TEX:.tex=)
	pdflatex -interaction=nonstopmode $(TEX) || true
	pdflatex -interaction=nonstopmode $(TEX)

clean:
	rm -f $(AUX)
	rm -f $(ARXIV_TARBALL) $(EPRINT_ZIP)
	rm -rf submission-staging/

# arXiv expects a single .tar.gz with .tex sources and .bib.
# No PDF (arXiv compiles it server-side).
arxiv: $(TEX) $(SECTIONS) $(BIB)
	mkdir -p submission-staging/arxiv
	mkdir -p submission-staging/arxiv/figures && cp $(FIGURES) submission-staging/arxiv/figures/
	cp $(TEX) $(SECTIONS) $(BIB) abstract.txt AUTHORS.example Makefile README.md \
	  submission-staging/arxiv/
	cd submission-staging/arxiv && tar czf ../../$(ARXIV_TARBALL) *
	rm -rf submission-staging/arxiv
	@echo "arXiv submission bundle ready: $(ARXIV_TARBALL)"

# IACR ePrint accepts PDF + optional source.
eprint: $(PDF)
	mkdir -p submission-staging/eprint/figures
	cp $(FIGURES) submission-staging/eprint/figures/
	cp $(PDF) $(TEX) $(SECTIONS) $(BIB) abstract.txt README.md Makefile \
	  submission-staging/eprint/
	cd submission-staging/eprint && zip -q -r ../../$(EPRINT_ZIP) .
	rm -rf submission-staging/eprint
	@echo "IACR ePrint submission bundle ready: $(EPRINT_ZIP)"

# Restore author block from AUTHORS (for non-anonymous releases).
deanonymize: $(TEX)
	@if [ ! -f AUTHORS ]; then \
	  echo "ERROR: AUTHORS file missing — create it with the author block"; \
	  echo "Example AUTHORS file contents:"; \
	  echo '  \\author{Alice \\and Bob \\\\ Affiliation}'; \
	  echo '  \\date{\\today}'; \
	  exit 1; \
	fi
	@cp $(TEX) $(TEX).anon-backup
	@python3 -c 'import re; \
tex = open("$(TEX)").read(); \
authors = open("AUTHORS").read().strip(); \
pa = re.search(r"^%\s*pdfauthor\s*=\s*(.+)$$", authors, re.M); \
authors = re.sub(r"^%.*\n?", "", authors, flags=re.M).strip(); \
tex = re.sub(r"\\author\{Anonymous\}\s*\\date\{\}", lambda _: authors + "\n", tex, count=1); \
tex = tex.replace("pdfauthor={Anonymous}", "pdfauthor={" + pa.group(1).strip() + "}") if pa else tex; \
open("$(TEX)", "w").write(tex); \
print("De-anonymized $(TEX); backup saved as $(TEX).anon-backup" + ("" if pa else "; pdfauthor left as Anonymous (add a %pdfauthor= line to AUTHORS)"))'

check: $(PDF)
	@echo "Checking $(PDF)..."
	@if grep -q 'LaTeX Warning: There were undefined references' $(TEX:.tex=.log); then \
	  echo "FAIL: unresolved references"; \
	  grep 'LaTeX Warning: Reference' $(TEX:.tex=.log) || true; \
	  exit 1; \
	else \
	  echo "All references resolved."; \
	fi
	@if grep -qE '(LaTeX|Package natbib) Warning: There were undefined citations' $(TEX:.tex=.log); then \
	  echo "FAIL: unresolved citations"; \
	  exit 1; \
	fi
	@if grep -E 'Overfull \\hbox' $(TEX:.tex=.log); then \
	  echo "FAIL: overfull boxes"; \
	  exit 1; \
	else \
	  echo "No overfull boxes."; \
	fi
	@echo "Page count:"
	@pdfinfo $(PDF) | grep -i pages
	@ls -lh $(PDF)

# Every figure as SVG and 150-dpi PNG (repository policy for performance
# figures). FIGDIR defaults to this round's benchmarks entry.
FIGDIR ?= ../../benchmarks/20261014-paper-figures
figures-export: figures-standalone.tex $(FIGURES)
	pdflatex -interaction=nonstopmode figures-standalone.tex >/dev/null
	mkdir -p $(FIGDIR)
	@i=1; for f in $$(sed -n 's/^\\input{figures\/\(.*\)}$$/\1/p' figures-standalone.tex); do \
	  pdftocairo -svg -f $$i -l $$i figures-standalone.pdf $(FIGDIR)/$$f.svg; \
	  pdftoppm -r 150 -png -singlefile -f $$i -l $$i figures-standalone.pdf $(FIGDIR)/$$f; \
	  i=$$((i+1)); done
	cp figures-standalone.pdf $(FIGDIR)/figures.pdf
	rm -f figures-standalone.aux figures-standalone.log figures-standalone.pdf
	@echo "exported $$(ls $(FIGDIR)/*.svg | wc -l) figures to $(FIGDIR)"
