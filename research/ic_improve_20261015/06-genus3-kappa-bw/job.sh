#!/bin/sh
# One recorded job: job.sh <out_dir> <args...>; skipped if the --out file exists and is non-empty.
OUT=$1; shift
BIN=${BIN:-../../../target/release/examples/ic_improve2_06_genus3_kappa_bw}
eval "last=\${$#}"
[ -s "$last" ] && exit 0
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
S=$(date +%s)
nice timeout 3600 prlimit --as=8000000000 -- "$BIN" "$@" 2>>"$OUT/progress.log" || echo "FAILED rc=$? $*" >>"$OUT/progress.log"
echo "JOB $(( $(date +%s) - S ))s $*" >>"$OUT/progress.log"
