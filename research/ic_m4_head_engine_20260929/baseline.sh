#!/usr/bin/env bash
# Reproduction control for the head-engine m = 4 rerun (PREREGISTRATION.md §3).
# Re-runs the chain and chain-holdout ladders with the reference and candidate
# policies on the head engine and compares every row EXACTLY with the head totals
# the first audit recorded at 1722bad1 (../ic_m4_exponent_audit_20260928/baseline/
# head-1722bad1).  The engine sources changed since 1722bad1 only by test code and
# one unused helper, so an exact match is expected; a mismatch is disclosed, never
# absorbed.  Each ladder runs pinned to one CPU under the benchmark lock.
#
#   BENCH=$WORK/bin/groebner_stage_bench-4ff512f2 CPU=3 \
#     research/ic_m4_head_engine_20260929/baseline.sh OUT_DIR
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(git -C "$here" rev-parse --show-toplevel)"
recorded="$here/../ic_m4_exponent_audit_20260928/baseline/head-1722bad1"
out="${1:?OUT_DIR}"
bench="${BENCH:?path to the groebner_stage_bench binary}"
cpu="${CPU:-3}"
for spec in "reference layout 0 complete" "candidate interleaved 1 complete"; do
  read -r arm order linear drop <<<"$spec"
  for suite in chain chain-holdout; do
    short=${suite/chain-holdout/holdout}
    dir="$out/$short-$arm"
    [ -e "$dir/stage.json" ] && { echo "exists, compared as is: $dir"; continue; }
    python3 "$root/tools/isolated_bench.py" busy -- taskset -c "$cpu" env -i PATH="$PATH" HOME="$HOME" \
      KIC_CHAIN_ORDER=$order KIC_LINEAR_ELIM=$linear KIC_F4_DROP=$drop \
      "$bench" --label "$arm" --ladder "$suite" --out "$dir" > /dev/null
  done
done
python3 - "$recorded" "$out" <<'PY'
import json, sys
recorded, out = sys.argv[1], sys.argv[2]
keys = ["word_ops", "f4_calls", "verdict_digest", "decomposed", "matrix_rows", "matrix_cols",
        "reductions", "splits", "infeasible_branches", "exhausted"]
named = {("chain", "K_0/2^9", 4), ("holdout", "K_1/2^15", 4)}
ok_all = ok_named = True
for short in ("chain", "holdout"):
    for arm in ("reference", "candidate"):
        f = {(r["curve"], r["m"]): r for r in json.load(open(f"{recorded}/{short}-{arm}/stage.json"))["rows"]}
        g = {(r["curve"], r["m"]): r for r in json.load(open(f"{out}/{short}-{arm}/stage.json"))["rows"]}
        for k, a in f.items():
            b = g.get(k)
            if b is None:
                print(f"      {short:<8} {arm:<9} {k[0]:<9} m={k[1]}  not in the rerun")
                ok_all = False
                ok_named = ok_named and (short, *k) not in named
                continue
            diff = [x for x in keys if a.get(x) != b.get(x)]
            tag = "NAMED " if (short, *k) in named else "      "
            print(f"{tag}{short:<8} {arm:<9} {k[0]:<9} m={k[1]}  1722bad1 {a['word_ops']:>11,}  "
                  f"head {b['word_ops']:>11,}  {'MATCH' if not diff else 'DIFF ' + ','.join(diff)}")
            ok_all = ok_all and not diff
            if (short, *k) in named:
                ok_named = ok_named and not diff
print(f"\nnamed cells reproduce exactly: {ok_named}; every row of both ladders: {ok_all}")
sys.exit(0 if ok_named else 1)
PY
