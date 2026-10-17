#!/bin/bash
# R08 on v5, as amendment 1 declares (PROTOCOL.md): the manifest, the
# candidate's pin, the A/A, the suite rows, the holdouts, the extension and
# the analysis.  The tests ran on the candidate's tree before this script
# (tests.log).  Untimed CPU work runs under the benchmark lock (`busy`);
# timed steps take it per process.  Every step resumes from disk.
# Usage: run.sh <scratch dir holding bin/> <candidate commit>.
set -u
SP=${1:?scratch dir}
CC=${2:?candidate commit}
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
R=$ROOT/research/ic_tool_program/rounds/R08-scan-fold
LOG=$SP/r08-run.log
P=$SP/bin/icprog-r10
ISO=$SP/bin/isolated_bench
BUSY="$ISO busy --wait --"
BASE=$SP/bin/ic-r06-v4-ab211e39
CAND=$SP/bin/ic-r08-v5-${CC:0:8}
BC=ab211e39733fd279f2bfab2fe4bfc5d1e7bfc564
log() { echo "== $(date -u +%T) $*" >> $LOG; }
grep -q 'tests exit 0' $R/tests.log || { log "stopped: no passing tests.log"; exit 1; }
arms=(--root $ROOT --isolate $ISO --base $BASE --cand $CAND --base-commit $BC --cand-commit $CC)
for step in manifest pin; do
  log $step
  $BUSY $P run r08 $step "${arms[@]}" > $SP/r08-$step.out 2>> $LOG
  rc=$?; echo "   exit $rc $(tr '\n' ' ' < $SP/r08-$step.out | head -c 200)" >> $LOG
  [ $rc = 0 ] || { log "stopped at $step"; exit 1; }
done
for step in aa compare holdout extend; do
  log $step
  $P run r08 $step "${arms[@]}" > $SP/r08-$step.out 2>> $LOG
  rc=$?; echo "   exit $rc $(tail -c 200 $SP/r08-$step.out | tr '\n' ' ')" >> $LOG
  [ $rc = 0 ] || { log "stopped at $step"; exit 1; }
done
log analyse
$P analyse r08 --root $ROOT > $R/analysis.json 2>> $LOG; echo "   exit $?" >> $LOG
log "all done"
