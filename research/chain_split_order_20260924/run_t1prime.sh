#!/usr/bin/env bash
# §5 of RESEARCH_CHAIN_SPLIT_ORDER.md: the supplementary holdout T1′, the
# reference and the registered candidate, three repetitions each, every
# control explicit.  Same conventions as run.sh; never overwrites a run.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here}"
bench="${BENCH:-target/release/examples/groebner_stage_bench}"
arms=(
  "reference layout 0 complete"
  "candidate interleaved 1 complete"
)
suite=chain-holdout-2
for spec in "${arms[@]}"; do
  read -r arm order linear drop <<<"$spec"
  for rep in 1 2 3; do
    dir="$out/$suite/$arm/rep$rep"
    [ -e "$dir/stage.json" ] && { echo "exists: $dir (never overwritten)"; continue; }
    KIC_CHAIN_ORDER=$order KIC_LINEAR_ELIM=$linear KIC_F4_DROP=$drop \
      "$bench" --label "$arm" --ladder "$suite" --out "$dir" > /dev/null
    echo "$suite $arm rep$rep done"
  done
done
python3 research/inherited_f4_20260922/compare_cross_tree.py "$out/$suite/reference" "$out/$suite/candidate" \
  --output "$out/comparisons/$suite.cross_tree.candidate-vs-reference.json" \
  > "$out/comparisons/$suite.cross_tree.candidate-vs-reference.log" 2>&1 || true
