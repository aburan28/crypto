#!/usr/bin/env bash
# README.md -> REPORT.html (pandoc, figures inlined) -> REPORT.pdf (headless Chromium).
cd "$(dirname "$0")/.."
cat > /tmp/report_style_$$.css <<'CSS'
body{font-family:Inter,Helvetica,Arial,sans-serif;max-width:60em;margin:1.5em auto;font-size:10.5pt;line-height:1.4;color:#1a1a1a}
table{border-collapse:collapse;font-size:9pt;margin:.6em 0}th,td{border:1px solid #ccc;padding:2px 6px}
img{max-width:100%}code,pre{font-size:8.5pt;background:#f4f4f4}pre{padding:.5em;overflow-x:auto;white-space:pre-wrap}
h1{font-size:17pt}h2{font-size:14pt;border-bottom:1px solid #ddd}h3{font-size:12pt}
CSS
pandoc README.md -s --embed-resources --metadata title="Item 2: solution-count census" -c /tmp/report_style_$$.css -o REPORT.html
/opt/pw-browsers/chromium-1194/chrome-linux/chrome --headless --no-sandbox --disable-gpu --no-pdf-header-footer \
  --print-to-pdf=REPORT.pdf "file://$PWD/REPORT.html" 2>/dev/null
rm -f /tmp/report_style_$$.css REPORT.html
