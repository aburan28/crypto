#!/bin/bash
# Round 3 of the PKM tower measurement (note §12.4): the cell M4b, its two
# systems one at a time, then the confirmation run if a system finishes at
# D = 7. The example is built from the commit that pre-registered the cell,
# after the identity check of §12.2 has passed. From the repository root:
#
#     cargo build --release --example pkm_tower_pilot
#     research/pkm_tower_round3_20260925/run_round3.sh
#
# The cell writes runs/<name>.jsonl (one row per system) and runs/<name>.log
# (the example's stderr, with each step's trace as it ends), and progress.txt
# records start and end times and exit statuses.
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:-$here/../../target/release/examples/pkm_tower_pilot}
mkdir -p "$here/runs"
cd "$here/runs"
ulimit -v 14000000
P1=2013265921
cell=(--engine tower --p $P1 --kinds kummer --m 4 --controls tower --t-min 4 --max-t-m4 4
      --planted 0 --random 2 --ladder-t none --budget 14400 --max-nnz 2000000000 --trace)
run() {
  name=$1
  shift
  echo "$(date -u +%FT%TZ) start $name: $*" >> ../progress.txt
  "$bin" "$@" --out "$name.jsonl" 2> "$name.log"
  # Read the status before anything else runs: `$(date)` would reset it.
  rc=$?
  echo "$(date -u +%FT%TZ) end $name exit $rc" >> ../progress.txt
}
run M4b-kummer-m4-p1-N16 "${cell[@]}"
# §12.4: a system that finished at D = 7 has its cell re-run whole with the
# degree bound at 6.
if python3 - <<'PY'
import json, sys
rows = [json.loads(l) for l in open("M4b-kummer-m4-p1-N16.jsonl") if l.strip()]
done = [r for r in rows if "N" in r and not r["timed_out"] and r["solving_degree_max"] == 7]
sys.exit(0 if done else 1)
PY
then
  run C-M4b-kummer-m4-p1-N16 "${cell[@]}" --cap 6
fi
echo "$(date -u +%FT%TZ) round 3 done" >> ../progress.txt
