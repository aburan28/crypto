#!/bin/sh
# Builds REPORT.pdf from REPORT.md, the protocol, the amendment, both result
# tables and both figures. Thin orchestration only:
#   sh build.sh
# The figures come from the experiment program, from recorded rows only:
#   cargo run --release --example nist_cyclotomic_controls -- --protocol v1 --plot
#   cargo run --release --example nist_cyclotomic_controls -- --protocol v2 --plot
# Needs pandoc and a Chromium binary (CHROME, default: the Playwright one).
set -e
here=$(cd "$(dirname "$0")" && pwd)
chrome="${CHROME:-/opt/pw-browsers/chromium-1194/chrome-linux/chrome}"
tmp="$here/.build"; mkdir -p "$tmp"
{
  cat "$here/REPORT.md"
  printf '\n\n---\n\n'
  sed 's/^#/##/' "$here/RESULTS_v2.md"
  printf '\n\n---\n\n'
  sed 's/^#/##/' "$here/RESULTS_v1.md"
  printf '\n\n---\n\n## Appendix A: PROTOCOL.md (frozen)\n\n'
  sed 's/^#/###/' "$here/PROTOCOL.md"
  printf '\n\n## Appendix B: AMENDMENT_1.md (frozen)\n\n'
  sed 's/^#/###/' "$here/AMENDMENT_1.md"
} > "$tmp/report.md"
cp "$here"/controls_v1.svg "$here"/controls_v2.svg "$tmp/"
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
pandoc "$tmp/report.md" -s --metadata title="Controlled re-test of the NIST cyclotomic / point-height note" \
  -c style.css -o "$tmp/report.html"
"$chrome" --headless=new --no-sandbox --disable-gpu --no-pdf-header-footer \
  --print-to-pdf="$here/REPORT.pdf" "file://$tmp/report.html" 2>/dev/null
rm -rf "$tmp"
ls -l "$here/REPORT.pdf"
