#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/exp7-fontconfig-cache}"
mkdir -p "$XDG_CACHE_HOME"

dot -Tsvg pipeline.dot -o pipeline.svg
dot -Tpdf pipeline.dot -o pipeline.pdf
inkscape peak_ratios.svg --export-filename=peak_ratios.pdf >/dev/null 2>&1

# GitHub displays the SVGs; the PDF embeds the corresponding vector PDFs.
report_for_pdf="$(mktemp /tmp/exp7-report-pdf-XXXXXX.md)"
trap 'rm -f "$report_for_pdf"' EXIT
sed -e 's/](pipeline.svg)/](pipeline.pdf)/g' \
    -e 's/](peak_ratios.svg)/](peak_ratios.pdf)/g' \
    RESULT.md > "$report_for_pdf"
pandoc "$report_for_pdf" --from=gfm --pdf-engine=xelatex \
  --resource-path=. --metadata title='EXP7 Fourier flatness sweep' \
  -V geometry:margin=18mm -V fontsize=10pt \
  -V monofont='DejaVu Sans Mono' -V colorlinks=true \
  -o RESULT.pdf
pdfinfo RESULT.pdf >/dev/null
