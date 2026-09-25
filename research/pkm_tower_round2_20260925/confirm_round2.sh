#!/bin/bash
# The confirmation runs of note §11.4, after the cells of run_round2.sh:
# every m = 2 cell of round 2 with a finished system at D ≥ 6 is re-run
# once, whole, with the degree bound at 5 (--cap 5) and otherwise the same
# flags. §11.6: each such system must then end with pairs above the bound,
# without refuting and without reaching a staircase stop.
#
#     research/pkm_tower_round2_20260925/confirm_round2.sh
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:-$here/../../target/release/examples/pkm_tower_pilot}
cd "$here/runs"
ulimit -v 14000000
common=(--engine tower --budget 7200 --max-nnz 2000000000 --ladder-t none --planted 0 --random 2)
run() {
  name=$1
  shift
  echo "$(date -u +%FT%TZ) start $name: $*" >> ../progress.txt
  "$bin" "${common[@]}" "$@" --out "$name.jsonl" 2> "$name.log"
  # Read the status before anything else runs: `$(date)` would reset it.
  rc=$?
  echo "$(date -u +%FT%TZ) end $name exit $rc" >> ../progress.txt
}
# The cells to confirm, one line each: file stem, prime, kind, control, t,
# and the staircase stop the cell ran with.
python3 - > ../confirm_cells.txt <<'PY'
import glob, json
cells = set()
for path in sorted(glob.glob("K*.jsonl") + glob.glob("I0*.jsonl")):
    for line in open(path):
        r = json.loads(line)
        if "N" not in r or r["m"] != 2 or r["timed_out"]:
            continue
        if r["solving_degree_max"] >= 6:
            cells.add((path[: -len(".jsonl")], r["p"], r["kind"], r["control"], r["t"]))
for stem, p, kind, control, t in sorted(cells):
    stop = "8" if stem.startswith(("K0", "I0")) else "none"
    print(stem, p, kind, control, t, stop)
PY
while read -r stem p kind control t stop; do
  extra=()
  if [ "$stop" != none ]; then
    extra=(--stop-below "$stop")
  fi
  run "C-$stem-t$t" --p "$p" --kinds "$kind" --m 2 --controls "$control" \
    --t-min "$t" --max-t "$t" --cap 5 "${extra[@]}"
done < ../confirm_cells.txt
echo "$(date -u +%FT%TZ) confirmations done" >> ../progress.txt
