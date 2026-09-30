#!/usr/bin/env bash
# Cut-out degree test (PREREGISTRATION.md): every registered cell, one at a time.
# Nothing here is timed (the metric is a degree), so no CPU pinning.  Outputs are
# never overwritten; a cell killed by its CPU or memory limit keeps its exit status
# and is censored, never negative evidence.
#
#   cargo build --release --example gf2_span
#   research/ic_tree_split_cutout_20260930/run.sh RUNS_DIR
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
runs="${1:?RUNS_DIR}"
mkdir -p "$(dirname "$runs")" && mkdir "$runs"
git -C "$here" rev-parse HEAD > "$runs/commit"
sha256sum "$here/cutout.py" "$here/../../target/release/examples/gf2_span" > "$runs/sha256"

cell() {  # object n a ell seed
  local name="$1-n$2a$3-l$4-s$5"
  set +e
  ( ulimit -t 3600; ulimit -v 12000000; exec python3 "$here/cutout.py" \
      --object "$1" --n "$2" --a "$3" --ell "$4" --seed "$5" ) > "$runs/$name.json" 2> "$runs/$name.stderr"
  echo $? > "$runs/$name.exit"
  set -e
}

for n in 17 19; do for a in 0 1; do for seed in 20260930 20261001; do for ell in 4 5 6 7 8; do
  cell T2 "$n" "$a" "$ell" "$seed"
done; done; done; done
for n in 13 15; do for a in 0 1; do for ell in 4 5; do
  cell Z "$n" "$a" "$ell" 20260930
done; done; done
