#!/bin/bash
# Reproduce: bash run.sh   (one core; per-cell time caps in seconds)
set -u
cd "$(dirname "$0")"
RAW=../../factor_base_dor_campaign_20261013/raw
B=dorv/target/release/dorv
PIN="taskset -c 0"; command -v taskset >/dev/null || PIN=""
(cd dorv && cargo build --release -q)
mkdir -p out/replay out/delta
stage() { # stage-name cap-per-cell files...
  local name=$1 cap=$2; shift 2
  for f in "$@"; do
    c=$(basename "$f" .json); s=${c##*_s}
    [ -s out/replay/$c.json ] && continue
    echo "[$name] $c"; $PIN $B replay-grid "$f" "$s" out/replay/$c.json "$cap" 2>>out/replay.log
  done
}
date -u +%FT%TZ > out/started
stage R1 120 $RAW/m2_n9_*.json $RAW/m2_n11_*.json $RAW/m2_n13_*.json $RAW/m2_n15_*.json
for n in 11 13; do
  [ -s out/delta/delta_n${n}_m2.json ] || $PIN $B replay-delta $RAW/delta_n${n}_m2.json out/delta/delta_n${n}_m2.json 300 2>>out/replay.log
done
stage R3 120 $RAW/m2_n17_*.json
stage R4 120 $RAW/m2_n19_d[2-7]_*.json
[ -s out/e2.json ] || $PIN $B e2 out/e2.json
[ -s out/null.json ] || $PIN $B null out/replay out/null.json 200 900
$PIN $B analyze $RAW out out/results.json
date -u +%FT%TZ > out/finished
