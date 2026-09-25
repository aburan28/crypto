#!/usr/bin/env bash
# Hoisting the Macaulay row cap out of the per-row loop and caching two
# per-node environment switches: main's binary (232fceec) against the
# change, the default policy on both, interleaved (arm order rotates with the
# repetition), three repetitions, five stage ladders and the twenty
# whole-logarithm seeds of RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md §4.  The
# change touches only environment lookups, so every counter must be equal;
# check.py checks that and tabulates wall time (a practicality note).
#
#   OLD=<dir holding main's ic and examples/groebner_stage_bench> research/env_reads_20260924/run.sh
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
old="${OLD:?set OLD to the directory of the binaries built from main}"
new="${NEW:-target/release}"
arms=("main $old" "change $new")
for rep in 1 2 3; do
  for k in 0 1; do
    read -r arm dir <<<"${arms[$(((k + rep) % 2))]}"
    for suite in frozen chain chain-holdout chain-holdout-2 r2-holdout; do
      ladder=()
      [ "$suite" != frozen ] && ladder=(--ladder "$suite")
      d="$here/stage/$suite/$arm/rep$rep"
      [ -e "$d/stage.json" ] && { echo "exists: $d"; continue; }
      "$dir/examples/groebner_stage_bench" --label "$arm" "${ladder[@]}" --out "$d" > /dev/null
      echo "$suite $arm rep$rep: exit $?"
    done
  done
done
cells=(
  "0 13 201 202 203 204 205 206 207 208 209 210 301 302 303 304 305"
  "0 9 201 202 203 204 205"
)
for rep in 1 2 3; do
  for cell in "${cells[@]}"; do
    read -r a n seeds <<<"$cell"
    for seed in $seeds; do
      for k in 0 1; do
        read -r arm dir <<<"${arms[$(((k + rep + seed) % 2))]}"
        mkdir -p "$here/e2e/$arm/rep$rep"
        f="$here/e2e/$arm/rep$rep/K${a}_2^${n}_seed$seed.json"
        [ -e "$f" ] && { echo "exists: $f"; continue; }
        timeout 1800 "$dir/ic" run --degree "$n" --curve-a "$a" --summands 3 --solver groebner \
          --random-target --seed "$seed" --batch 1 --json > "$f" 2> "$f.stderr"
        echo "e2e $arm rep$rep K_$a/2^$n seed $seed: exit $?"
      done
    done
  done
done
