#!/bin/sh
# lit2026 W11 (b)-(e): the remaining measurements, one after another
# (thin orchestration only; every number is printed by a Rust example).
#   run_rest.sh <bin_dir> <raw_dir>
set -eu
bin=$1; raw=$2
log="$raw/run_rest.log"
: > "$log"
step() { echo "== $1 $(date -u +%FT%TZ) load: $(uptime | sed -E 's/.*load averages?: //')" >> "$log"; }
step binary_sect163
"$bin/lit2026_w11_binary_xonly" sect163 8 2000 7 > "$raw/binary_sect163.jsonl"
step binary_torsion19
"$bin/lit2026_w11_binary_xonly" torsion19 20 9 > "$raw/binary_torsion19.jsonl"
step biext_k2
"$bin/lit2026_w11_biext_k2" 4 400 7 > "$raw/biext_k2.jsonl"
step dp_stats_K0_53
"$bin/lit2026_w11_dp_stats" 0 53 4,8,12 64 > "$raw/dp_stats.jsonl"
step dp_stats_K0_61
"$bin/lit2026_w11_dp_stats" 0 61 4,8,12 32 >> "$raw/dp_stats.jsonl"
step ladder_30_48
"$bin/lit2026_w11_challenge_ladder" solve 30 48 2 8 > "$raw/ladder_solve.jsonl"
step ladder_50_56
"$bin/lit2026_w11_challenge_ladder" solve 50 56 2 4 >> "$raw/ladder_solve.jsonl"
step ladder_58_60
"$bin/lit2026_w11_challenge_ladder" solve 58 60 2 2 >> "$raw/ladder_solve.jsonl"
step done
