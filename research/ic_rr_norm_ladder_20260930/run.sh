#!/usr/bin/env bash
# RR norm-form refutation-degree ladder (PREREGISTRATION.md): every registered cell, one
# at a time, pinned to CPU 3, the whole run under the benchmark lock.  Outputs are never
# overwritten; a cell killed by its CPU or memory limit keeps its lines (censored, never
# negative evidence) and its exit status is recorded.
#
#   BIN=$WORK/bin/rr_degree_ladder research/ic_rr_norm_ladder_20260930/run.sh RUNS_DIR
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(git -C "$here" rev-parse --show-toplevel)"
runs="${1:?RUNS_DIR}"
bin="${BIN:?set BIN to the rr_degree_ladder binary}"

if [ "${2:-}" = "--cells" ]; then
  for curve in "0 13 5" "1 17 6" "1 19 6"; do
    read -r a n top <<<"$curve"
    for ell in $(seq 2 "$top"); do
      name="K${a}n${n}l${ell}"
      set +e
      ( ulimit -t 3600; ulimit -v 10000000; exec taskset -c 3 env -i PATH="$PATH" F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000 "$bin" \
          --a "$a" --n "$n" --ell "$ell" --unsat 4 --max-draws 256 --d-max 7 --d-max-x4 9 --seed 20260930 \
          --out "$runs/$name.jsonl" ) > /dev/null 2> "$runs/$name.stderr"
      echo $? > "$runs/$name.exit"
      set -e
    done
  done
  exit 0
fi

[ -e "$runs" ] && { echo "$runs exists; never overwritten" >&2; exit 1; }
mkdir -p "$runs"
sha256sum "$bin" > "$runs/binary.sha256"
git -C "$here" rev-parse HEAD > "$runs/branch_head.txt"
date -u +%FT%TZ > "$runs/started_utc.txt"
BIN="$bin" python3 "$root/tools/isolated_bench.py" busy -- bash "$here/run.sh" "$runs" --cells
date -u +%FT%TZ > "$runs/finished_utc.txt"
python3 "$here/analyze.py" "$runs" | tee "$runs/readout.txt"
