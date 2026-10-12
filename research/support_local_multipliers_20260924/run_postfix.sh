#!/usr/bin/env bash
# Amendment runs of RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md §6: the registered
# binary (built from 1f0d9751) against the binary with the monomial-mask hash
# fix (d69936e1), both arms on each, interleaved (the arm order rotates with
# the repetition) so that wall ratios are not back-to-back drift.  The hash
# fix touches only uncharged row packing, so every counter must equal the
# registered runs'; check_identity.py checks that.  Nothing registered is
# overwritten: output goes to postfix/ and e2e_postfix/.
#
#   (cd <worktree at 1f0d9751> && CARGO_TARGET_DIR=<old> cargo build --release --example groebner_stage_bench --bin ic)
#   cargo build --release --example groebner_stage_bench --bin ic
#   OLD=<old>/release research/support_local_multipliers_20260924/run_postfix.sh
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
old="${OLD:?set OLD to the target/release directory built from 1f0d9751}"
new="${NEW:-target/release}"
arms=(
  "registered-reference $old occurring"
  "registered-candidate $old support"
  "fixed-reference $new occurring"
  "fixed-candidate $new support"
)
run_arm() { # name dir multipliers -- command...
  local multipliers=$1; shift
  KIC_F4_MULTIPLIERS=$multipliers KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete "$@"
}
for rep in 1 2 3; do
  for k in 0 1 2 3; do
    read -r arm dir multipliers <<<"${arms[$(((k + rep - 1) % 4))]}"
    for suite in frozen chain chain-holdout chain-holdout-2 r2-holdout; do
      ladder=()
      [ "$suite" != frozen ] && ladder=(--ladder "$suite")
      d="$here/postfix/$suite/$arm/rep$rep"
      [ -e "$d/stage.json" ] && { echo "exists: $d (never overwritten)"; continue; }
      run_arm "$multipliers" "$dir/examples/groebner_stage_bench" --label "$arm" "${ladder[@]}" --out "$d" > /dev/null
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
      for k in 0 1 2 3; do
        read -r arm dir multipliers <<<"${arms[$(((k + rep + seed - 1) % 4))]}"
        mkdir -p "$here/e2e_postfix/$arm/rep$rep"
        f="$here/e2e_postfix/$arm/rep$rep/K${a}_2^${n}_seed$seed.json"
        [ -e "$f" ] && { echo "exists: $f"; continue; }
        run_arm "$multipliers" timeout 1800 "$dir/ic" run --degree "$n" --curve-a "$a" --summands 3 \
          --solver groebner --random-target --seed "$seed" --batch 1 --json > "$f" 2> "$f.stderr"
        echo "e2e $arm rep$rep K_$a/2^$n seed $seed: exit $?"
      done
    done
  done
done
