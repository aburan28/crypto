#!/usr/bin/env bash
# README.md -> README.html -> README.pdf (pandoc, then headless Chromium).
set -eu
cd "$(dirname "$0")"
CSS=$(mktemp --suffix=.css)
cat > "$CSS" <<'CSS'
body{font-family:system-ui,sans-serif;max-width:1000px;margin:24px auto;font-size:12.5px;line-height:1.45;color:#1a1a1a;background:#fff}
table{border-collapse:collapse;font-size:10.5px;margin:8px 0}th,td{border:1px solid #bbb;padding:3px 5px}
th{background:#f0f0f0}code{font-size:11px}pre{background:#f6f6f6;padding:8px;overflow-x:auto;font-size:9.5px;white-space:pre-wrap}
img{max-width:100%}
CSS
pandoc README.md -s --metadata pagetitle="Item 05: in-engine linear harvest" -c "$CSS" --embed-resources -o README.html
/opt/pw-browsers/chromium-1194/chrome-linux/chrome --headless --no-sandbox --disable-gpu \
  --no-pdf-header-footer --print-to-pdf=README.pdf "file://$PWD/README.html" 2>/dev/null
rm -f "$CSS"
ls -la README.pdf
