#!/usr/bin/env bash
# Registered whole-logarithm runs of RESEARCH_CHAIN_SPLIT_ORDER.md §4: the
# reference and the candidate on one `ic` binary, the same curve, factor base
# and seed, so both draw the same trial points; `ic run` checks [k]G = Q
# against the planted k.  Every run is kept, whatever its status.
#
#   cargo build --release --bin ic
#   research/chain_split_order_20260924/e2e.sh [OUT_DIR]
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here}/e2e"
ic="${IC:-target/release/ic}"
cells=(
  "0 13 1 2 3 4 5 6 7 8 9 10 101 102 103 104 105"
  "0 9 1 2 3 4 5"
)
for arm in reference candidate; do
  if [ "$arm" = reference ]; then env=(KIC_CHAIN_ORDER=layout KIC_LINEAR_ELIM=0 KIC_F4_DROP=complete)
  else env=(KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete); fi
  mkdir -p "$out/$arm"
  for spec in "${cells[@]}"; do
    read -r a n seeds <<<"$spec"
    for seed in $seeds; do
      f="$out/$arm/K${a}_2^${n}_seed$seed.json"
      [ -e "$f" ] && { echo "exists: $f"; continue; }
      env "${env[@]}" timeout 1800 "$ic" run --degree "$n" --curve-a "$a" --summands 3 \
        --solver groebner --random-target --seed "$seed" --batch 1 --json > "$f" 2> "$f.stderr"
      echo "$arm K_$a/2^$n seed $seed: exit $?"
    done
  done
done
