#!/usr/bin/env bash
# Baseline reproduction control (PREREGISTRATION.md §5.1): re-run the frozen chain and
# chain-holdout ladders of research/chain_split_order_20260924 with the reference and
# candidate policies and compare every row with the frozen rep1 stage.json exactly
# (word XORs, F4 calls, verdict digest, decomposed count, matrix rows/cols, reductions,
# splits).  The two cells the survey names are K_0/2^9 m=4 (chain) and K_1/2^15 m=4
# (chain-holdout); the other rows of the two ladders are compared as well.
#
#   BENCH=$WORK/bin/groebner_stage_bench-2809b498 \
#     research/ic_m4_exponent_audit_20260928/baseline.sh OUT_DIR
#
# Existing OUT_DIR/*/stage.json files are never overwritten; they are compared as they are.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
frozen="$(cd "$here/../chain_split_order_20260924" && pwd)"
out="${1:?OUT_DIR}"
bench="${BENCH:?path to the groebner_stage_bench binary}"
for spec in "reference layout 0 complete" "candidate interleaved 1 complete"; do
  read -r arm order linear drop <<<"$spec"
  for suite in chain chain-holdout; do
    short=${suite/chain-holdout/holdout}
    dir="$out/$short-$arm"
    [ -e "$dir/stage.json" ] && { echo "exists, compared as is: $dir"; continue; }
    env -i PATH="$PATH" HOME="$HOME" KIC_CHAIN_ORDER=$order KIC_LINEAR_ELIM=$linear \
      KIC_F4_DROP=$drop "$bench" --label "$arm" --ladder "$suite" --out "$dir" > /dev/null
  done
done
python3 - "$frozen" "$out" <<'EOF'
import json, sys
frozen, out = sys.argv[1], sys.argv[2]
keys = ["word_ops", "f4_calls", "verdict_digest", "decomposed", "matrix_rows", "matrix_cols",
        "reductions", "splits", "infeasible_branches", "exhausted"]
named = {("chain", "K_0/2^9", 4), ("chain-holdout", "K_1/2^15", 4)}
ok_all = ok_named = True
for suite, short in (("chain", "chain"), ("chain-holdout", "holdout")):
    for arm in ("reference", "candidate"):
        f = {(r["curve"], r["m"]): r for r in json.load(open(f"{frozen}/{suite}/{arm}/rep1/stage.json"))["rows"]}
        g = {(r["curve"], r["m"]): r for r in json.load(open(f"{out}/{short}-{arm}/stage.json"))["rows"]}
        for k, a in f.items():
            b = g.get(k)
            if b is None:
                print(f"      {suite:<13} {arm:<9} {k[0]:<9} m={k[1]}  not in the rerun")
                if (suite, *k) in named:
                    ok_named = False
                continue
            diff = [x for x in keys if a.get(x) != b.get(x)]
            tag = "NAMED " if (suite, *k) in named else "      "
            print(f"{tag}{suite:<13} {arm:<9} {k[0]:<9} m={k[1]}  frozen {a['word_ops']:>11,}  "
                  f"rerun {b['word_ops']:>11,}  {'MATCH' if not diff else 'DIFF ' + ','.join(diff)}")
            ok_all = ok_all and not diff
            if (suite, *k) in named:
                ok_named = ok_named and not diff
print(f"\nnamed cells reproduce exactly: {ok_named}; every row of both ladders: {ok_all}")
sys.exit(0 if ok_named else 1)
EOF
