#!/usr/bin/env bash
# Registered runs of research/notes/ecc2k130/RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md §2:
# the reference (occurring-variable multipliers, the default merged in #690)
# and the candidate (support-local multipliers) on five suites, three
# repetitions each, one binary; every control explicit and recorded.
#
#   cargo build --release --example groebner_stage_bench
#   research/support_local_multipliers_20260924/run.sh [OUT_DIR]
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here}"
bench="${BENCH:-target/release/examples/groebner_stage_bench}"
arms=("reference occurring" "candidate support")
for suite in frozen chain chain-holdout chain-holdout-2 r2-holdout; do
  ladder=()
  [ "$suite" != frozen ] && ladder=(--ladder "$suite")
  for spec in "${arms[@]}"; do
    read -r arm multipliers <<<"$spec"
    for rep in 1 2 3; do
      dir="$out/$suite/$arm/rep$rep"
      [ -e "$dir/stage.json" ] && { echo "exists: $dir (never overwritten)"; continue; }
      KIC_F4_MULTIPLIERS=$multipliers KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete \
        "$bench" --label "$arm" "${ladder[@]}" --out "$dir" > /dev/null
      echo "$suite $arm rep$rep done"
    done
  done
done
mkdir -p "$out/comparisons"
for suite in frozen chain chain-holdout chain-holdout-2 r2-holdout; do
  tag="$suite.cross_tree.candidate-vs-reference"
  python3 research/inherited_f4_20260922/compare_cross_tree.py "$out/$suite/reference" "$out/$suite/candidate" \
    --output "$out/comparisons/$tag.json" > "$out/comparisons/$tag.log" 2>&1
  echo "$tag: exit $?"
done
