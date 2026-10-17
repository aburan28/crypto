#!/bin/sh
# Document rendering only: no experiment execution.
set -eu
cd "$(dirname "$0")"
dot -Tsvg flow.dot -o flow.svg
dot -Tpng -Gdpi=150 flow.dot -o flow.png
pandoc README.md --standalone --pdf-engine=xelatex \
  -V mainfont='DejaVu Sans' -V monofont='DejaVu Sans Mono' \
  -V geometry:margin=18mm -V fontsize=10pt -o report.pdf
