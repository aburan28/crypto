#!/bin/sh
# Document rendering only; no experiment or statistical analysis.
set -eu
cd "$(dirname "$0")"
dot -Tsvg decision.dot -o decision.svg
dot -Tpng -Gdpi=150 decision.dot -o decision.png
pandoc README.md --standalone --pdf-engine=xelatex \
  -V mainfont='DejaVu Sans' -V monofont='DejaVu Sans Mono' \
  -V geometry:margin=22mm -V fontsize=10pt -o report.pdf
