#!/usr/bin/env bash
# The closure-at-scale study (PREREGISTRATION.md).  Thin orchestration only:
# every computation is examples/dreg_ladder (draws) and examples/dreg_closure.
#
#   research/dreg_closure_scale_20261016/run.sh LADDER_BIN CLOSURE_BIN [gen|measure|control]
#
# gen:     write draws/<cell>.jsonl from the frozen seed (exact solution count,
#          XL to degree 3 only).
# measure: one process per unsatisfiable draw, cheapest cell first; a draw whose
#          runs/<cell>.d<K>.jsonl is non-empty is skipped, so a restart loses only
#          the draw in flight.
# control: the first two satisfiable draws of each fresh cell with N <= 36.
set -u
LADDER=$1
CLOSURE=$2
STEP=${3:-measure}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
DRAWS=$HERE/draws
RUNS=$HERE/runs
mkdir -p "$DRAWS" "$RUNS"
export F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-4}
SEED=20261016
# m:n:ell, cheapest first
FRESH="3:16:6 3:18:6 4:12:3 3:19:7 4:14:3 3:21:7 4:13:4"
R=research
COMMITTED="
$R/dreg_surplus_control_20260925/runs/cell-11-5-6.u0.jsonl
$R/dreg_surplus_control_20260925/runs/cell-11-5-6.u1.jsonl
$R/dreg_surplus_control_20260925/runs/cell-11-5-6.u2.jsonl
$R/dreg_surplus_control_20260925/runs/cell-11-5-6.u3.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-13-5-6.u0.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-13-5-6.u1.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-13-5-6.u2.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-13-5-6.u3.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-15-5-6.u0.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-15-5-6.u1.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-15-5-6.u2.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-15-5-6.u3.jsonl
"
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$RUNS/queue.log"; }

one() { # name file draw [extra flags]
  local name=$1 file=$2 k=$3; shift 3
  local out=$RUNS/$name.d$k.jsonl
  if [ -s "$out" ]; then log "$name d$k: done already, skipped"; return; fi
  log "$name d$k: start"
  local start; start=$(date +%s)
  (cd "$REPO" && prlimit --as=15000000000 timeout 7200 \
     "$CLOSURE" --to 4 --max-cols 130000 --draw "$k" "$@" "$file") > "$out" 2> "$RUNS/$name.d$k.log"
  log "$name d$k: exit $? after $(( $(date +%s) - start )) s"
}

case $STEP in
gen)
  for c in $FRESH; do
    IFS=: read -r m n ell <<< "$c"
    (cd "$REPO" && "$LADDER" --cells "$n:$ell:3" --m "$m" --unsat 4 --ffd-max 3 --seed $SEED) \
      > "$DRAWS/m$m-$n-$ell.jsonl" 2> "$DRAWS/m$m-$n-$ell.log"
  done ;;
measure)
  for f in $COMMITTED; do
    name="committed-$(basename "$f" .jsonl)"
    k=$(grep -v satisfiable "$REPO/$f" | grep -o '"draw":[0-9]*' | head -1 | cut -d: -f2)
    one "$name" "$f" "$k"
  done
  for c in $FRESH; do
    IFS=: read -r m n ell <<< "$c"
    f=research/dreg_closure_scale_20261016/draws/m$m-$n-$ell.jsonl
    for k in $(grep -v '"satisfiable"' "$REPO/$f" | grep -o '"draw":[0-9]*' | cut -d: -f2); do
      one "fresh-m$m-$n-$ell" "$f" "$k"
    done
  done ;;
control)
  for c in 3:16:6 3:18:6 4:12:3; do
    IFS=: read -r m n ell <<< "$c"
    f=research/dreg_closure_scale_20261016/draws/m$m-$n-$ell.jsonl
    for k in $(grep '"satisfiable"' "$REPO/$f" | grep -o '"draw":[0-9]*' | cut -d: -f2 | head -2); do
      one "control-m$m-$n-$ell" "$f" "$k" --satisfiable
    done
  done ;;
esac
