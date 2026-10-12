#!/usr/bin/env bash
# §4 of RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md: whole logarithms, both arms on
# one `ic` binary, the same curve, base and seed; `ic run` verifies [k]G = Q.
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here}/e2e"
ic="${IC:-target/release/ic}"
cells=(
  "0 13 201 202 203 204 205 206 207 208 209 210 301 302 303 304 305"
  "0 9 201 202 203 204 205"
)
for spec in "reference occurring" "candidate support"; do
  read -r arm multipliers <<<"$spec"
  mkdir -p "$out/$arm"
  for cell in "${cells[@]}"; do
    read -r a n seeds <<<"$cell"
    for seed in $seeds; do
      f="$out/$arm/K${a}_2^${n}_seed$seed.json"
      [ -e "$f" ] && { echo "exists: $f"; continue; }
      KIC_F4_MULTIPLIERS=$multipliers KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete \
        timeout 1800 "$ic" run --degree "$n" --curve-a "$a" --summands 3 --solver groebner \
        --random-target --seed "$seed" --batch 1 --json > "$f" 2> "$f.stderr"
      echo "$arm K_$a/2^$n seed $seed: exit $?"
    done
  done
done
