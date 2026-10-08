#!/bin/bash
# Native (uninstrumented) runs of every rung at the frozen n=53 L=1024 K=440 cell.
set -u
cd "$(dirname "$0")"
BIN=${BIN:?set BIN to the strong-rho example binary}
for R in 0 1 2 3; do
  uptime > rung${R}_native_load_before.txt
  START=$(date +%s.%N)
  KIC_RHO_RUNG=$R KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=n53-ks-growing-1024-v1 \
    "$BIN" 53 0 signed_frobenius 1024 531310 > rung${R}_native.jsonl 2> rung${R}_native.stderr.log
  RC=$?
  END=$(date +%s.%N)
  uptime > rung${R}_native_load_after.txt
  echo "rung=$R rc=$RC wall_s=$(echo "$END - $START" | bc)" | tee -a rung_native_walls.txt
done
echo NATIVE_RUNGS_DONE >> rung_native_walls.txt
