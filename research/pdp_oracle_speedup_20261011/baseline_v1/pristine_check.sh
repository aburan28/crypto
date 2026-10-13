#!/usr/bin/env bash
# Build the harness against the UNEDITED lib sources of main (4b6e09b78) in a separate target dir,
# then restore the edited sources.  The edits are cfg-gated counters and one visibility change.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
S=${PS0_SCRATCH:-/tmp/ps0_scratch}; mkdir -p "$S"
cp src/cryptanalysis/semaev_decomp.rs $S/semaev_decomp.edited.rs
cp Cargo.toml $S/Cargo.edited.toml
trap 'cp $S/semaev_decomp.edited.rs src/cryptanalysis/semaev_decomp.rs; cp $S/Cargo.edited.toml Cargo.toml' EXIT
git show 4b6e09b78f010d72de866b7eed2236782c47f656:src/cryptanalysis/semaev_decomp.rs > src/cryptanalysis/semaev_decomp.rs
git show 4b6e09b78f010d72de866b7eed2236782c47f656:Cargo.toml > Cargo.toml
# the example's `folds()` use is only under the opcount feature, which main does not declare
CARGO_TARGET_DIR=target/pristine python3 tools/isolated_bench.py busy -- cargo build --release --example pdp_speedup
# then, on the quiet machine:
#   for l in 8 9; do valgrind --tool=callgrind --callgrind-out-file=/tmp/cg_pristine_l$l.out --toggle-collect='*run_one*' \
#     target/pristine/release/examples/pdp_speedup one --l $l --corpus research/pdp_oracle_speedup_20261011/baseline_v1/corpus/corpus_l$l.txt; done
# and compare "Collected" with callgrind/cg_l$l.log (equal: 77504893 and 337453547).
