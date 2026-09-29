#!/bin/bash
# Whole-process retired instructions (callgrind Ir) for one rung at the frozen cell.
set -u
cd "$(dirname "$0")"
R=${1:?rung 0..3}
BIN=${BIN:?set BIN}
uptime > callgrind_rung${R}_load_before.txt
KIC_RHO_RUNG=$R KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=n53-ks-growing-1024-v1 \
  valgrind --tool=callgrind --callgrind-out-file=callgrind_rung${R}.out \
  "$BIN" 53 0 signed_frobenius 1024 531310 \
  > callgrind_rung${R}.stdout.jsonl 2> callgrind_rung${R}.stderr.log
echo "rc=$?" >> callgrind_rung${R}.stderr.log
uptime > callgrind_rung${R}_load_after.txt
echo "CALLGRIND_RUNG${R}_DONE" >> callgrind_rung${R}.stderr.log
