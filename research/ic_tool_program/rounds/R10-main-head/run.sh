#!/bin/bash
# R10, as declared (PROTOCOL.md): the tests, the manifest, both arms' pins,
# the A/A, the suite rows, the holdouts, the extension, the callgrind
# cross-check and the analysis.  Untimed CPU work runs under the benchmark
# lock (`busy`); timed steps take it per process.  Every step resumes from
# disk.  Usage: run.sh <scratch dir holding bin/>.
set -u
SP=${1:?scratch dir}
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
R=$ROOT/research/ic_tool_program/rounds/R10-main-head
LOG=$SP/r10-run.log
P=$SP/bin/icprog-r10
ISO=$SP/bin/isolated_bench
BUSY="$ISO busy --wait --"
BASE=$SP/bin/ic-v3-995ea207-rebuilt
CAND=$SP/bin/ic-main-92246530
BC=995ea2071cc7453877a503d30eb561ae82cddab9
CC=92246530c4f333724374e70897330bc38fa159a1
log() { echo "== $(date -u +%T) $*" >> $LOG; }
arms=(--root $ROOT --isolate $ISO --base $BASE --cand $CAND --base-commit $BC --cand-commit $CC)
if ! grep -q 'tests exit 0' $LOG 2>/dev/null; then
  log tests
  ( cd $ROOT && CARGO_TARGET_DIR=$SP/target-main $BUSY cargo test --release --locked --lib -- koblitz_fast koblitz_index_calculus semaev_decomp ) > $R/tests.log 2>&1
  rc=$?; echo "   tests exit $rc $(grep -E '^test result:' $R/tests.log | tr '\n' ' ')" >> $LOG
  [ $rc = 0 ] || { log "stopped: a test failed"; exit 1; }
fi
for step in manifest pin; do
  log $step
  $BUSY $P run r10 $step "${arms[@]}" > $SP/r10-$step.out 2>> $LOG
  rc=$?; echo "   exit $rc $(tr '\n' ' ' < $SP/r10-$step.out | head -c 200)" >> $LOG
  [ $rc = 0 ] || { log "stopped at $step"; exit 1; }
done
log "base pin"
mkdir -p $R/base-pin
$BUSY $P run r10 pin --root $ROOT --runs $R/base-pin --isolate $ISO --base $BASE --cand $BASE > $SP/r10-base-pin.out 2>> $LOG
rc=$?; echo "   exit $rc" >> $LOG; [ $rc = 0 ] || { log "stopped: the base does not pin"; exit 1; }
for step in aa compare holdout extend callgrind; do
  log $step
  $P run r10 $step "${arms[@]}" > $SP/r10-$step.out 2>> $LOG
  rc=$?; echo "   exit $rc $(tail -c 200 $SP/r10-$step.out | tr '\n' ' ')" >> $LOG
  [ $rc = 0 ] || { log "stopped at $step"; exit 1; }
done
log analyse
$P analyse r10 --root $ROOT > $R/analysis.json 2>> $LOG; echo "   exit $?" >> $LOG
log "all done"
