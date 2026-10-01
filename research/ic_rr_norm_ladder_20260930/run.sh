#!/usr/bin/env bash
# RR norm-form refutation-degree ladder (PREREGISTRATION.md): the registered cells in three
# lanes, one lane per CPU (1, 2, 3), the whole run under the benchmark lock.  A degree is
# not a timed quantity, so lanes may share the machine; nothing else runs beside them.
# Outputs are never overwritten; a cell killed by its CPU or memory limit keeps its lines
# (censored, never negative evidence) and its exit status is recorded.
#
#   BIN=$WORK/bin/rr_degree_ladder research/ic_rr_norm_ladder_20260930/run.sh RUNS_DIR [--resume]
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(git -C "$here" rev-parse --show-toplevel)"
runs="${1:?RUNS_DIR}"
bin="${BIN:?set BIN to the rr_degree_ladder binary}"

cell() {  # cpu a n ell  -- up to three attempts of 6000 CPU-s, resumed at draw granularity
  local cpu="$1" a="$2" n="$3" ell="$4" dmax=7 attempt rc
  [ "$ell" -ge 6 ] && dmax=6
  local name="K${a}n${n}l${ell}"
  if [ -e "$runs/$name.exit" ] && [ "$(cat "$runs/$name.exit")" = 0 ]; then return; fi
  for attempt in 1 2 3; do
    [ "$(wc -l < "$runs/$name.attempts" 2>/dev/null || echo 0)" -ge 3 ] && return
    set +e
    ( ulimit -t 6000; ulimit -v 4500000; exec taskset -c "$cpu" env -i PATH="$PATH" F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000 "$bin" \
        --a "$a" --n "$n" --ell "$ell" --unsat 4 --max-draws 256 --d-max "$dmax" --d-max-x4 9 --seed 20260930 \
        --resume --out "$runs/$name.jsonl" ) > /dev/null 2>> "$runs/$name.stderr"
    rc=$?
    set -e
    echo "$(date -u +%FT%TZ) rc=$rc" >> "$runs/$name.attempts"
    echo "$rc" > "$runs/$name.exit"
    [ "$rc" = 0 ] && return
  done
}

if [ "${2:-}" = "--cells" ]; then
  # lanes balanced by the smoke's costs: l = 5 and 6 cells are the heavy ones
  ( cell 1 1 19 6; cell 1 0 13 2; cell 1 0 13 3; cell 1 0 13 4; cell 1 0 13 5 ) &
  ( cell 2 1 17 6; cell 2 1 17 2; cell 2 1 17 3; cell 2 1 17 4; cell 2 1 17 5 ) &
  ( cell 3 1 19 5; cell 3 1 19 2; cell 3 1 19 3; cell 3 1 19 4 ) &
  wait
  exit 0
fi

if [ "${2:-}" = "--resume" ]; then
  # the container may reclaim a run; resume into the same directory, recording the fact
  echo "$(date -u +%FT%TZ) resume at $(git -C "$here" rev-parse HEAD) binary $(sha256sum "$bin" | cut -c1-16)" >> "$runs/resumes.txt"
else
  [ -e "$runs" ] && { echo "$runs exists; never overwritten (use --resume)" >&2; exit 1; }
  mkdir -p "$runs"
  sha256sum "$bin" > "$runs/binary.sha256"
  sha256sum "$root/examples/rr_degree_ladder.rs" > "$runs/source.sha256"
  git -C "$here" rev-parse HEAD > "$runs/branch_head.txt"
  date -u +%FT%TZ > "$runs/started_utc.txt"
fi
BIN="$bin" python3 "$root/tools/isolated_bench.py" busy -- bash "$here/run.sh" "$runs" --cells
date -u +%FT%TZ > "$runs/finished_utc.txt"
python3 "$here/analyze.py" "$runs" | tee "$runs/readout.txt"
