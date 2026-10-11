#!/usr/bin/env bash
# Runs every registered cell of PREREGISTRATION.md §4 once, in order, into
# runs/registered/. Thin orchestration only; the harness is
# examples/tree_split_structured.rs.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
out="$here/runs/registered"
mkdir -p "$out"
cd "$root"
cargo build --release --example tree_split_structured
bin="$root/target/release/examples/tree_split_structured"
git rev-parse HEAD > "$out/commit"
sha256sum "$bin" | cut -d' ' -f1 > "$out/sha256"
"$bin" selftest > "$out/selftest.txt"
: > "$out/cells.jsonl"
: > "$out/cells.log"
cell() {  # vtype n a seed ell
  local t0=$(date -u +%s)
  timeout 1800 "$bin" cell --vtype "$1" --n "$2" --a "$3" --seed "$4" --ell "$5" >> "$out/cells.jsonl" \
    && echo "$(date -u +%FT%TZ) ok   $* $(( $(date -u +%s) - t0 ))s" >> "$out/cells.log" \
    || echo "$(date -u +%FT%TZ) FAIL $* exit=$?" >> "$out/cells.log"
}
for a in 0 1; do
  cell subfield 21 $a 20261011 7
  for n in 17 19 21; do for l in 4 5 6 7 8; do cell poly $n $a 20261011 $l; done; done
  for s in 20261011 20261012; do for l in 4 5 6 7 8; do cell frob 21 $a $s $l; done; done
  cell frob 17 $a 20261011 8
  for l in 4 5 6 7 8; do cell rand 21 $a 20261011 $l; done
done
cd "$here"
"$bin" readout runs/registered/cells.jsonl > runs/registered/readout.txt
"$bin" lean runs/registered/cells.jsonl lean/TreeSplit/Data.lean
