#!/bin/sh
# Builds REPORT.pdf from the note and the two SVGs. Thin orchestration only:
#   sh build.sh
# Needs pandoc and a Chromium binary (CHROME, default: the Playwright one).
set -e
here=$(cd "$(dirname "$0")" && pwd)
note="$here/../notes/index-calculus/RESEARCH_PDP_SAT_SMT_NAGAO_SURVEY.md"
chrome="${CHROME:-/opt/pw-browsers/chromium-1194/chrome-linux/chrome}"
tmp="$here/.build"; mkdir -p "$tmp"
# The note links the study directory as ../../pdp_sat_smt_survey_20261010/;
# inside the study directory those links are local, and the figures are
# inlined after §8 so the PDF carries them.
sed -e 's#\.\./\.\./pdp_sat_smt_survey_20261010/##g' "$note" > "$tmp/report.md"
cat >> "$tmp/report.md" <<'MD'

---

## Figures

![Per-target oracle work against the enumeration floor](oracle_ladder.svg)

![Where a solver constant sits in the product law](product_law.svg)
MD
cp "$here"/oracle_ladder.svg "$here"/product_law.svg "$tmp/"
cat > "$tmp/style.css" <<'CSS'
body { font-family: Helvetica, Arial, sans-serif; font-size: 10.5pt; line-height: 1.35; max-width: 100%; margin: 0; }
table { border-collapse: collapse; font-size: 8pt; width: 100%; }
th, td { border: 1px solid #bbb; padding: 2px 4px; vertical-align: top; text-align: left; }
code { font-size: 8.5pt; }
pre { font-size: 8.5pt; background: #f4f4f4; padding: 6px; }
img { max-width: 100%; }
h1 { font-size: 17pt; } h2 { font-size: 13.5pt; margin-top: 1.2em; } h3 { font-size: 11.5pt; }
@page { size: A4 landscape; margin: 12mm; }
CSS
pandoc "$tmp/report.md" -s --metadata title="SAT and SMT solvers, WDSat and Nagao decomposition for the PDP" \
  -c style.css -o "$tmp/report.html"
"$chrome" --headless=new --no-sandbox --disable-gpu --no-pdf-header-footer \
  --print-to-pdf="$here/REPORT.pdf" "file://$tmp/report.html" 2>/dev/null
rm -rf "$tmp"
ls -l "$here/REPORT.pdf"
