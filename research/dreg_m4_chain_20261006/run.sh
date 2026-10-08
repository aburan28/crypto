#!/usr/bin/env bash
# Run the m = 4 refutation-degree study (PREREGISTRATION.md), one job at a
# time.  Thin orchestration only: every computation is dreg_ladder's.
#
#   research/dreg_m4_chain_20261006/run.sh BINARY [JOB ...]
#
# A JOB is N-L-D for a whole cell (four unsatisfiable draws, every draw
# recorded) or N-L-D.uK for one unsatisfiable draw in its own process.  With
# no JOB the registered order runs.  A job whose runs/cell-JOB.jsonl is
# non-empty is skipped, so a restart loses only the job in flight.
set -u
BIN=$1
shift
HERE=$(cd "$(dirname "$0")" && pwd)
RUNS=$HERE/runs
mkdir -p "$RUNS"
SEED=20261006
export F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000
export KIC_SPARSE_DENSE_FINISH=1 KIC_SPARSE_F5=1
export KIC_SPARSE_DENSE_BUDGET_MB=${KIC_SPARSE_DENSE_BUDGET_MB:-11000}
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-4}

ORDER="4-1-7 5-1-7 6-1-7 5-2-8"
for cell in 7-2-7 6-2-8 8-2-7; do
  for k in 0 1 2 3; do ORDER="$ORDER $cell.u$k"; done
done

log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$RUNS/queue.log"; }

for job in ${@:-$ORDER}; do
  out=$RUNS/cell-$job.jsonl
  if [ -s "$out" ]; then
    log "$job: done already, skipped"
    continue
  fi
  cell=${job%%.u*}
  spec=${cell//-/:}
  if [ "$cell" = "$job" ]; then
    mode=(--unsat 4 --controls 0)
  else
    mode=(--unsat-index "${job##*.u}")
  fi
  cmd=("$BIN" --m 4 --cells "$spec" --ffd-max 5 --seed "$SEED" "${mode[@]}")
  log "$job: start budget_mb=$KIC_SPARSE_DENSE_BUDGET_MB f5=1 threads=$RAYON_NUM_THREADS cmd=${cmd[*]}"
  start=$(date +%s)
  "${cmd[@]}" > "$out" 2> "$RUNS/cell-$job.log"
  code=$?
  log "$job: exit $code after $(( $(date +%s) - start )) s"
done
