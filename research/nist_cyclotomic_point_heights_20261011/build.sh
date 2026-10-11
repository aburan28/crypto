#!/bin/sh
# Builds REPORT.pdf from the note and the two SVGs. Thin orchestration only:
#   sh build.sh
# Needs pandoc and a Chromium binary (CHROME, default: the Playwright one).
set -e
here=$(cd "$(dirname "$0")" && pwd)
note="$here/../notes/cm-isogeny/RESEARCH_NIST_CYCLOTOMIC_POINT_HEIGHTS.md"
chrome="${CHROME:-/opt/pw-browsers/chromium-1194/chrome-linux/chrome}"
tmp="$here/.build"; mkdir -p "$tmp"
# The note links the study directory as ../../nist_cyclotomic_point_heights_20261011/;
# inside the study directory those links are local.
sed -e 's#\.\./\.\./nist_cyclotomic_point_heights_20261011/##g' "$note" > "$tmp/report.md"
cat >> "$tmp/report.md" <<'MD'

---

## Figures

![Where the cyclotomic structure acts, and where it does not](structure.svg)

![Cheon / Brown–Gallant cost against published powers](cheon_cost.svg)
MD
cp "$here"/structure.svg "$here"/cheon_cost.svg "$tmp/"
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
pandoc "$tmp/report.md" -s --metadata title="Cyclotomic factorisation of P-192, P-224, P-256: what reaches the points, and point heights" \
  -c style.css -o "$tmp/report.html"
"$chrome" --headless=new --no-sandbox --disable-gpu --no-pdf-header-footer \
  --print-to-pdf="$here/REPORT.pdf" "file://$tmp/report.html" 2>/dev/null
rm -rf "$tmp"
ls -l "$here/REPORT.pdf"
