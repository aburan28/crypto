#!/usr/bin/env bash
# Specialisation study (PREREGISTRATION.md).  Thin orchestration only: every
# computation is examples/dreg_closure --fix.
#
#   research/dreg_closure_specialise_20261017/run.sh CLOSURE_BIN
#
# One process per (draw, arm); a run whose output is non-empty is skipped, so a
# restart loses only the run in flight.
set -u
CLOSURE=$1
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
RUNS=$HERE/runs
mkdir -p "$RUNS"
export F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-4}
D=research/dreg_closure_scale_20261016/draws
CELLS="m3-19-7:0,3,5,7 m3-21-7:0,2,3,6"
ARMS="A0:6=0 A1:6=1 B:6=0,13=0,20=0 C:21=0"
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$RUNS/queue.log"; }
for cell in $CELLS; do
  name=${cell%%:*}; draws=${cell#*:}
  for k in ${draws//,/ }; do
    for arm in $ARMS; do
      a=${arm%%:*}; fix=${arm#*:}
      out=$RUNS/$name.d$k.$a.jsonl
      if [ -s "$out" ]; then log "$name d$k $a: done already, skipped"; continue; fi
      log "$name d$k $a: start fix=$fix"
      start=$(date +%s)
      (cd "$REPO" && prlimit --as=15000000000 timeout 7200 "$CLOSURE" --to 4 --max-cols 130000 \
         --draw "$k" --fix "$fix" "$D/$name.jsonl") > "$out" 2> "$RUNS/$name.d$k.$a.log"
      log "$name d$k $a: exit $? after $(( $(date +%s) - start )) s"
    done
  done
done
