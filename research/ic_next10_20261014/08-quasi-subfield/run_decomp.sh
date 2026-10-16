#!/usr/bin/env bash
# Arm D: QSP / structured subspace bases vs 4 random subspaces, same instrument.
set -u
cd "$(dirname "$0")/../../.."
B=${QSF_DECOMP_BIN:-target/release/examples/ic_next10_08_qsf_decomp}
OUT=research/ic_next10_20261014/08-quasi-subfield
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 F4_F2_MAX_ROWS=60000 F4_F2_MAX_COLS=120000
run() { # file n k base yield refute dmax seed
  timeout 2700 nice prlimit --as=700000000 -- $B "$2" "$3" "$4" "$5" "$6" "$7" "$8" >> "$OUT/$1" 2>> "$OUT/decomp_stderr.txt" \
    || echo "{\"n\":$2,\"k\":$3,\"base\":\"$4\",\"status\":\"failed_or_timeout\"}" >> "$OUT/decomp_failures.txt"
}
# n = 31 runs at d_max = 3: the degree-4 closure exceeded the 700 MB cap
# (allocation failure during closure readback) on the first attempt.
CFGS=${QSF_CFGS:-"decomp_n7_k4.jsonl 7 4 g=0x17 2000 6 5
decomp_n18_k9.jsonl 18 9 g=0x201 2000 6 5
decomp_n17_k9.jsonl 17 9 g=0x279 2000 6 5
decomp_n31_k15.jsonl 31 15 g=0x808b 1000 6 3
decomp_n31_k16.jsonl 31 16 g=0x10117 1000 6 3"}
while read -r cfg; do
  set -- $cfg
  run "$1" "$2" "$3" "$4" "$5" "$6" "$7" 11
  for s in 1 2 3 4; do run "$1" "$2" "$3" "random=$s" "$5" "$6" "$7" 11; done
done <<< "$CFGS"
echo ALLDONE
