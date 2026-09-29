#!/bin/bash
# WALK-CONSTANT.md section 11.6: the device harness's seed-to-seed scatter,
# under rule v2 and, as an independent replicate, rule v1.  The command
# lines are the declaration's; about 46 minutes on four cores.
#
# v2 is this tree's device harness; v1 is the same source built from
# 825f3a84, the last tree with rule v1.  Both are compiled exactly as the
# Makefile's walk-constant target compiles walk-constant-host-h8 (that target
# also runs four harness jobs, so it is not called here).
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$(git rev-parse --show-toplevel)
B=$(mktemp -d)
trap 'rm -rf "$B"' EXIT
FLAGS="-O2 -std=c++17 -DECC_WALK_TABLE=1 -DECC_TABLE_PIVOT_BYTES=0 -Wno-unknown-pragmas -Wno-unused-function -pthread -DECC_TABLE_BRANCHES=8"
mkdir -p "$B/v1"
git -C "$ROOT" archive 825f3a84 ecc2k130/src ecc2k130/include ecc2k130/generated | tar -x -C "$B/v1"
# shellcheck disable=SC2086
g++ $FLAGS ../../src/walkconstant.cpp -o "$B/v2-h8"
# shellcheck disable=SC2086
g++ $FLAGS "$B/v1/ecc2k130/src/walkconstant.cpp" -o "$B/v1-h8"
{
  echo "v2: $(git -C "$ROOT" rev-parse HEAD)$(git -C "$ROOT" diff --quiet HEAD -- ecc2k130/src ecc2k130/include ecc2k130/generated || echo ' (uncommitted changes)')"
  echo "v1: 825f3a84"
  sha256sum "$B/v2-h8" "$B/v1-h8" | sed "s|$B/||"
} > scatter-build.txt

# Disjoint seeds: with one seed the two rules would draw the same tables,
# salts and start points trial for trial.  --threads 4 is part of the
# design, since the seed-to-stream map depends on the thread count.
for s in $(seq 240 255); do
  "$B/v2-h8" --n 23 --walk table --dist ecc2k130 --walks 8 --trials 20000 --seed "$s" --threads 4
done > scatter-v2.jsonl 2> scatter-v2.log
for s in $(seq 260 275); do
  "$B/v1-h8" --n 23 --walk table --dist ecc2k130 --walks 8 --trials 20000 --seed "$s" --threads 4
done > scatter-v1.jsonl 2> scatter-v1.log
python3 scatter.py | tee scatter.txt
