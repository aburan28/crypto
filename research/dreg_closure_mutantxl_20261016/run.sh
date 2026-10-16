#!/usr/bin/env bash
# Run the closure-vs-XL study (PREREGISTRATION.md).  Thin orchestration only:
# every computation is examples/dreg_closure.
#
#   research/dreg_closure_mutantxl_20261016/run.sh BINARY
#
# One output file per input file, in the registered order (cheapest first).
# A file whose runs/<name>.jsonl is non-empty is skipped, so a restart loses
# only the file in flight.  Each file gets a 3-hour wall limit.
set -u
BIN=$1
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
RUNS=$HERE/runs
mkdir -p "$RUNS"
export F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-4}

R=research
INPUTS="
$R/dreg_fixed_surplus_20260923/runs/cell-5-2-7.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-7-2-7.jsonl
$R/dreg_ell_grid_20260925/runs/cell-8-2-7.jsonl
$R/dreg_ell_grid_20260925/runs/cell-10-2-7.jsonl
$R/dreg_ell_grid_20260925/runs/cell-12-2-7.jsonl
$R/dreg_ell_grid_20260925/runs/cell-4-3-7.jsonl
$R/dreg_m4_chain_20261006/runs/cell-4-1-7.jsonl
$R/dreg_m4_chain_20261006/runs/cell-5-1-7.jsonl
$R/dreg_m4_chain_20261006/runs/cell-6-1-7.jsonl
$R/dreg_ell_grid_20260925/runs/cell-5-3-7.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-7-3-7.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-9-3-7.jsonl
$R/dreg_m4_chain_20261006/runs/cell-5-2-8.jsonl
$R/dreg_surplus_control_20260925/runs/cell-8-4-7.jsonl
$R/dreg_surplus_control_20260925/runs/cell-7-4-7.jsonl
$R/dreg_m4_chain_20261006/runs/cell-6-2-8.u0.jsonl
$R/dreg_m4_chain_20261006/runs/cell-6-2-8.u1.jsonl
$R/dreg_m4_chain_20261006/runs/cell-6-2-8.u2.jsonl
$R/dreg_m4_chain_20261006/runs/cell-6-2-8.u3.jsonl
$R/dreg_m4_chain_20261006/runs/cell-9-3-6.u0.jsonl
$R/dreg_m4_chain_20261006/runs/cell-9-3-6.u1.jsonl
$R/dreg_m4_chain_20261006/runs/cell-9-3-6.u2.jsonl
$R/dreg_m4_chain_20261006/runs/cell-9-3-6.u3.jsonl
$R/dreg_m4_chain_20261006/runs/cell-7-2-7.u0.jsonl
$R/dreg_m4_chain_20261006/runs/cell-7-2-7.u1.jsonl
$R/dreg_m4_chain_20261006/runs/cell-7-2-7.u2.jsonl
$R/dreg_m4_chain_20261006/runs/cell-7-2-7.u3.jsonl
$R/dreg_m4_chain_20261006/runs/cell-10-3-6.u0.jsonl
$R/dreg_m4_chain_20261006/runs/cell-10-3-6.u1.jsonl
$R/dreg_m4_chain_20261006/runs/cell-10-3-6.u2.jsonl
$R/dreg_m4_chain_20261006/runs/cell-10-3-6.u3.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-11-4-6.jsonl
$R/dreg_m4_chain_20261006/runs/cell-8-2-7.u0.jsonl
$R/dreg_m4_chain_20261006/runs/cell-8-2-7.u1.jsonl
$R/dreg_m4_chain_20261006/runs/cell-8-2-7.u2.jsonl
$R/dreg_m4_chain_20261006/runs/cell-8-2-7.u3.jsonl
$R/dreg_fixed_surplus_20260923/runs/cell-13-4-6.jsonl
$R/dreg_degree7_ell5_20260930/runs/cell-10-5-7.u0.jsonl
$R/dreg_degree7_ell5_20260930/runs/cell-10-5-7.u1.jsonl
$R/dreg_degree7_ell5_20260930/runs/cell-10-5-7.u2.jsonl
$R/dreg_degree7_ell5_20260930/runs/cell-10-5-7.u3.jsonl
"

log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$RUNS/queue.log"; }

for input in $INPUTS; do
  study=$(echo "$input" | cut -d/ -f2)
  name="$(echo "$study" | sed 's/^dreg_//; s/_2026.*//')-$(basename "$input" .jsonl)"
  out=$RUNS/$name.jsonl
  if [ -s "$out" ]; then
    log "$name: done already, skipped"
    continue
  fi
  log "$name: start max_cols=70000 input=$input"
  start=$(date +%s)
  (cd "$REPO" && timeout 10800 "$BIN" --max-cols 70000 "$input") > "$out" 2> "$RUNS/$name.log"
  code=$?
  log "$name: exit $code after $(( $(date +%s) - start )) s"
done
