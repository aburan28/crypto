#!/bin/bash
# R11 on v6, as declared (PROTOCOL.md): the manifest, the
# candidate's pin, the A/A, the suite rows, the holdouts, the extension and
# the analysis.  The tests ran on the candidate's tree before this script
# (tests.log).  Untimed CPU work runs under the benchmark lock (`busy`);
# timed steps take it per process.  Every step resumes from disk.
# Usage: run.sh <scratch dir holding bin/> <candidate commit>.
set -u
SP=${1:?scratch dir}
CC=${2:?candidate commit}
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
R=$ROOT/research/ic_tool_program/rounds/R11-build-streams
LOG=$SP/r11-run.log
P=$SP/bin/icprog-r11
ISO=$SP/bin/isolated_bench
BUSY="$ISO busy --wait --"
BASE=$SP/bin/ic-r08-v5-3d4c5c0e
CAND=$SP/bin/ic-r11-v6-${CC:0:8}
BC=3d4c5c0ef88057e8ece09dcc91e8fad4ff7e781a
log() { echo "== $(date -u +%T) $*" >> $LOG; }
grep -q 'tests exit 0' $R/tests.log || { log "stopped: no passing tests.log"; exit 1; }
arms=(--root $ROOT --isolate $ISO --base $BASE --cand $CAND --base-commit $BC --cand-commit $CC)
for step in manifest pin; do
  log $step
  $BUSY $P run r11 $step "${arms[@]}" > $SP/r11-$step.out 2>> $LOG
  rc=$?; echo "   exit $rc $(tr '\n' ' ' < $SP/r11-$step.out | head -c 200)" >> $LOG
  [ $rc = 0 ] || { log "stopped at $step"; exit 1; }
done
for step in aa compare holdout extend; do
  log $step
  $P run r11 $step "${arms[@]}" > $SP/r11-$step.out 2>> $LOG
  rc=$?; echo "   exit $rc $(tail -c 200 $SP/r11-$step.out | tr '\n' ' ')" >> $LOG
  [ $rc = 0 ] || { log "stopped at $step"; exit 1; }
done
log analyse
$P analyse r11 --root $ROOT > $R/analysis.json 2>> $LOG; echo "   exit $?" >> $LOG
log "all done"
