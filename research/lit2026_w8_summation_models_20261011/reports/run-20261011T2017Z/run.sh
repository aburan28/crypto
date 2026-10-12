#!/bin/sh
# Thin orchestration only: runs the Rust harness once per part.
B=${B:-/private/tmp/claude-501/-Volumes-SSD990-crypto-elliptic-curve-crypto-research-88c7cc/ab9482d3-a14b-4bea-8cfa-130f6184453a/scratchpad/w8/w8bin}
R=research/lit2026_w8_summation_models_20261011/reports/run-20261011T2017Z/raw
mkdir -p "$R"
for part in "$@"; do
  echo "== $part start $(date -u +%FT%TZ) $(uptime)"
  /usr/bin/time -l "$B" "$part" --rounds 5 --targets 8 $EXTRA --out "$R/$part$SUF.json" 2> "$R/$part$SUF.stderr"
  echo "== $part exit $? $(date -u +%FT%TZ)"
done
