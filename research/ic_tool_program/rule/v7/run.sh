#!/bin/bash
# The rule's comparison at v7, as PROTOCOL.md declares: the manifest, every
# row (with the pin), the reference check, the claims and the analysis.
# Every step resumes what its run tree holds.  Usage: run.sh <scratch dir>.
set -u
SP=${1:?scratch dir}
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
D=$ROOT/research/ic_tool_program/rule/v7
LOG=$SP/rule7-run.log
P=$SP/bin/icprog-rule7
ISO=$SP/bin/isolated_bench
IC=$SP/bin/ic-r11-v6-c2af9b7c
FIX=$SP/bin/koblitz_rho_fixture-v7
C=c2af9b7ccd92469830048dc2deaf97e1d194b27a
log() { echo "== $(date -u +%T) $*" >> $LOG; }
cd $ROOT
step() {
  log "$1"
  shift
  "$@" >> $LOG 2>&1
  rc=$?; echo "   exit $rc" >> $LOG
  [ $rc = 0 ] || { log "stopped"; exit 1; }
}
step manifest $P rule manifest --comparison v7 --ic $IC --ic-commit $C --isolate $ISO
step all $P rule all --comparison v7 --ic $IC --ic-commit $C --isolate $ISO
step reference $P rule reference --comparison v7 --fixture $FIX --fixture-commit a475cde33904207fea38d850d7279cc87e0d27b8 --isolate $ISO
step claims $P rule claims --comparison v7
log analyse
$P rule analyse --comparison v7 > $D/analysis.json 2>> $LOG; echo "   exit $?" >> $LOG
log "all done"
