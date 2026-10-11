#!/usr/bin/env bash
# Runs every registered cell of PREREGISTRATION.md §4 once, in order, into
# runs/registered/. Thin orchestration; the harness is
# examples/tree_split_structured_v2.rs.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
out="$here/runs/registered"
mkdir -p "$out"
cd "$root"
cargo build --release --example tree_split_structured_v2
bin="$root/target/release/examples/tree_split_structured_v2"
git rev-parse HEAD > "$out/commit"
sha256sum "$bin" | cut -d' ' -f1 > "$out/sha256"
"$bin" selftest > "$out/selftest.txt"
: > "$out/cells.jsonl"
: > "$out/cells.log"
S=20261011
cell() {  # object vtype n a seed ell
  local t0=$(date -u +%s)
  timeout 3600 "$bin" cell --object "$1" --vtype "$2" --n "$3" --a "$4" --seed "$5" --ell "$6" >> "$out/cells.jsonl" \
    && echo "$(date -u +%FT%TZ) ok   $* $(( $(date -u +%s) - t0 ))s" >> "$out/cells.log" \
    || echo "$(date -u +%FT%TZ) FAIL $* exit=$?" >> "$out/cells.log"
}
for a in 0 1; do
  cell T2 subfield 21 $a $S 7
  for s in $S 20261012; do for l in 5 6 7 8; do cell T2 frob2 21 $a $s $l; done; done
  for l in 7 8; do cell T2 frob2 35 $a $S $l; done
  for n in 17 19 21 23; do for l in 4 5 6 7 8; do cell T2 nbw $n $a $S $l; done; done
  for n in 23 25; do for l in 4 5 6 7 8; do cell T2 poly $n $a $S $l; done; done
  for n in 13 15; do for l in 4 5; do cell Z poly $n $a $S $l; cell Z nbw $n $a $S $l; done; done
  for l in 4 5; do cell Z frob2 15 $a $S $l; done
done
cd "$here"
"$bin" readout runs/registered/cells.jsonl > runs/registered/readout.txt
"$bin" lean runs/registered/cells.jsonl lean/TreeSplit2/Data.lean
