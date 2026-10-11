#!/bin/sh
# lit2026 W11 (e): the challenge-ladder rho campaign (thin orchestration).
#   run_ladder.sh <bin_dir> <out.jsonl>
# 64 seeds per level at 30-48 bits, 16 at 50-56, 8 at 58-60.
set -eu
bin=$1; out=$2
: > "$out"
"$bin/lit2026_w11_challenge_ladder" solve 30 48 2 64 >> "$out"
"$bin/lit2026_w11_challenge_ladder" solve 50 56 2 16 >> "$out"
"$bin/lit2026_w11_challenge_ladder" solve 58 60 2 8 >> "$out"
