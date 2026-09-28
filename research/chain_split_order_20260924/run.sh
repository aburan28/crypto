#!/usr/bin/env bash
# Registered runs of research/notes/ecc2k130/RESEARCH_CHAIN_SPLIT_ORDER.md §2:
# every arm of the O × L × D factorial on the frozen, chain and chain-holdout
# ladders, three repetitions each, on one binary.  Every control is set
# explicitly on every arm and recorded in each stage.json's `policy`.
#
#   cargo build --release --example groebner_stage_bench
#   research/chain_split_order_20260924/run.sh [OUT_DIR]
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here}"
bench="${BENCH:-target/release/examples/groebner_stage_bench}"

# arm  KIC_CHAIN_ORDER  KIC_LINEAR_ELIM  KIC_F4_DROP
arms=(
  "reference layout 0 complete"
  "D layout 0 rebuild"
  "L layout 1 complete"
  "LD layout 1 rebuild"
  "O interleaved 0 complete"
  "OD interleaved 0 rebuild"
  "candidate interleaved 1 complete"
  "OLD interleaved 1 rebuild"
)
for suite in frozen chain chain-holdout; do
  ladder=()
  [ "$suite" != frozen ] && ladder=(--ladder "$suite")
  for spec in "${arms[@]}"; do
    read -r arm order linear drop <<<"$spec"
    for rep in 1 2 3; do
      dir="$out/$suite/$arm/rep$rep"
      [ -e "$dir/stage.json" ] && { echo "exists: $dir (never overwritten)"; continue; }
      KIC_CHAIN_ORDER=$order KIC_LINEAR_ELIM=$linear KIC_F4_DROP=$drop \
        "$bench" --label "$arm" "${ladder[@]}" --out "$dir" > /dev/null
      echo "$suite $arm rep$rep done"
    done
  done
done
