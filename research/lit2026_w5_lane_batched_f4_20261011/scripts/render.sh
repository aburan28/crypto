#!/bin/sh
# Render a Typst source to PDF beside it, and its pages to PNG under
# /Volumes/SSD990-2/crypto/targets/lit2026-w5/tmp/pages/ for inspection.
PY=${PY:-/private/tmp/claude-501/-Volumes-SSD990-crypto-elliptic-curve-crypto-research-88c7cc/ab9482d3-a14b-4bea-8cfa-130f6184453a/scratchpad/venv/bin/python}
src=$1
name=$2
pages=/Volumes/SSD990-2/crypto/targets/lit2026-w5/tmp/pages
mkdir -p "$pages"
"$PY" -I -c "import typst,sys; typst.compile(sys.argv[1], output=sys.argv[1][:-4]+'.pdf'); typst.compile(sys.argv[1], output=sys.argv[2]+'/'+sys.argv[3]+'-{p}.png', format='png', ppi=110)" "$src" "$pages" "$name"
ls "$pages" | grep "^$name-"
