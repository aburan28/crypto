#!/bin/sh
# One recorded job: job.sh <out_dir> <args...>; skipped if the --out file exists and is non-empty.
# Round 4 rules: single-threaded, nice, 4.5 GB address-space cap, 60-minute wall limit.
OUT=$1; shift
BIN=${BIN:-../../../target/release/examples/ic_improve4_01_genus3_rel}
eval "last=\${$#}"
[ -s "$last" ] && exit 0
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
S=$(date +%s)
nice timeout 3600 prlimit --as=4500000000 -- "$BIN" "$@" 2>>"$OUT/progress.log" || echo "FAILED rc=$? $*" >>"$OUT/progress.log"
echo "JOB $(( $(date +%s) - S ))s $*" >>"$OUT/progress.log"
