#!/bin/bash
# The rule's comparison at v6, as PROTOCOL.md declares: the manifest, every
# row (with the pin), the reference check, the claims and the analysis.
# Every step resumes what its run tree holds.  Usage: run.sh <scratch dir>.
set -u
SP=${1:?scratch dir}
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
D=$ROOT/research/ic_tool_program/rule/v6
LOG=$SP/rule6-run.log
P=$SP/bin/icprog-rule6
ISO=$SP/bin/isolated_bench
IC=$SP/bin/ic-r08-v5-3d4c5c0e
FIX=$SP/bin/koblitz_rho_fixture-v6
C=3d4c5c0ef88057e8ece09dcc91e8fad4ff7e781a
log() { echo "== $(date -u +%T) $*" >> $LOG; }
cd $ROOT
step() {
  log "$1"
  shift
  "$@" >> $LOG 2>&1
  rc=$?; echo "   exit $rc" >> $LOG
  [ $rc = 0 ] || { log "stopped"; exit 1; }
}
step manifest $P rule manifest --comparison v6 --ic $IC --ic-commit $C --isolate $ISO
step all $P rule all --comparison v6 --ic $IC --ic-commit $C --isolate $ISO
step reference $P rule reference --comparison v6 --fixture $FIX --fixture-commit $C --isolate $ISO
step claims $P rule claims --comparison v6
log analyse
$P rule analyse --comparison v6 > $D/analysis.json 2>> $LOG; echo "   exit $?" >> $LOG
log "all done"
