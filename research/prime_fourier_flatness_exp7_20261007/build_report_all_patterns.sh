#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/exp7-fontconfig-cache}"
mkdir -p "$XDG_CACHE_HOME"

dot -Tsvg pipeline_all_patterns.dot -o pipeline_all_patterns.svg
dot -Tpdf pipeline_all_patterns.dot -o pipeline_all_patterns.pdf
inkscape peak_ratios_all_patterns.svg \
  --export-filename=peak_ratios_all_patterns.pdf >/dev/null 2>&1

report_for_pdf="$(mktemp /tmp/exp7-all-report-pdf-XXXXXX.md)"
trap 'rm -f "$report_for_pdf"' EXIT
sed -e 's/](pipeline_all_patterns.svg)/](pipeline_all_patterns.pdf)/g' \
    -e 's/](peak_ratios_all_patterns.svg)/](peak_ratios_all_patterns.pdf)/g' \
    RESULT_ALL_PATTERNS.md > "$report_for_pdf"
pandoc "$report_for_pdf" --from=gfm --pdf-engine=xelatex \
  --resource-path=. \
  -V geometry:margin=18mm -V fontsize=10pt \
  -V mainfont='DejaVu Serif' -V monofont='DejaVu Sans Mono' \
  -V colorlinks=true \
  -o RESULT_ALL_PATTERNS.pdf
pdfinfo RESULT_ALL_PATTERNS.pdf >/dev/null
