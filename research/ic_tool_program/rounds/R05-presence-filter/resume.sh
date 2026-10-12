#!/bin/bash
# R05's holdouts resumed on the native runner after the container rebuild
# (IC_TOOL_PROGRAM.md §10a); then the declared extension check.
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
B=$SP/bin
ARGS="--root $SP/r05run-wt --runs $SP/r05-runs --base $B/ic-r03-cand-on-c1a2e5f8-30f6c153 --cand $B/ic-r05-cand-on-30f6c153-edcb0bec --isolate $B/isolated_bench-n2"
echo "== resume $(date -u) icprog $(sha256sum $B/icprog-n2 | cut -c1-16) isolated_bench $(sha256sum $B/isolated_bench-n2 | cut -c1-16)"
$B/icprog-n2 run r05 holdout $ARGS
echo "holdout exit $? $(date -u)"
$B/icprog-n2 run r05 extend $ARGS
echo "extend exit $? $(date -u)"
