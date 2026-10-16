#!/bin/sh
# One recorded job: job.sh <out_dir> <skip_if_nonempty_file|-> <args...>
OUT=$1; SKIP=$2; shift 2
BIN=${BIN:-../../../target/release/examples/ic_improve2_01_genus3_linalg}
[ "$SKIP" != "-" ] && [ -s "$SKIP" ] && exit 0
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
nice timeout 3600 prlimit --as=4000000000 -- "$BIN" "$@" 2>>"$OUT/progress.log" || echo "FAILED rc=$? $*" >>"$OUT/progress.log"
