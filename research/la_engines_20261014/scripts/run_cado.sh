#!/usr/bin/env bash
# Thin orchestration: solve one exported core with CADO-NFS bwc (GF(p),
# inhomogeneous, right nullspace) and verify the answer natively.
# usage: run_cado.sh <cols> <weight> <seed> <thr e.g. 2x2> <m> <n> <outdir>
set -euo pipefail
cols=$1 weight=$2 seed=$3 thr=$4 bm=$5 bn=$6 out=$7
CADO=${CADO_BUILD:-/home/user/ext/cado-nfs/build}
BENCH=${BENCH_BIN:-/home/user/crypto/target/release/examples/la_engine_bench}
P=4611686018427387847
dir=$(mktemp -d "${TMPDIR:-/tmp}/cado-core.XXXXXX")
trap 'rm -rf "$dir"' EXIT
"$BENCH" export-core --cols "$cols" --weight "$weight" --seed "$seed" --dir "$dir" > "$dir/export.json"
"$CADO/linalg/bwc/mf_scan2" --withcoeffs "$dir/matrix.bin" > "$dir/scan.log" 2>&1
start=$(date +%s.%N)
"$CADO/linalg/bwc/bwc.pl" :complete prime=$P matrix="$dir/matrix.bin" \
  rhs="$dir/rhs.txt" nullspace=right mpi=1x1 thr="$thr" m="$bm" n="$bn" \
  wdir="$dir/wdir" > "$dir/bwc.log" 2>&1 || { tail -30 "$dir/bwc.log" >&2; exit 1; }
end=$(date +%s.%N)
sol=$(ls "$dir"/wdir/K.sols* | head -1)
"$BENCH" check-cado --cols "$cols" --weight "$weight" --seed "$seed" --dir "$dir" --sol "$sol" > "$dir/check.json"
mkdir -p "$out"
tag="cado-c${cols}-w${weight}-s${seed}-t${thr}-m${bm}n${bn}"
cp "$dir/bwc.log" "$out/$tag.bwc.log"
jq -c --argjson e "$(cat "$dir/export.json")" --argjson c "$(cat "$dir/check.json")" \
  --argjson steps "$(grep -oE 'Timings for [a-z_]+: \(wct\) [0-9.]+' "$dir/bwc.log" | awk '{gsub(":","",$3); t[$3]+=$5} END{printf "{"; s=""; for(k in t){printf "%s\"%s\":%s", s, k, t[k]; s=","} printf "}"}')" \
  --arg thr "$thr" --arg wall "$(awk -v a="$start" -v b="$end" 'BEGIN{print b-a}')" --arg bm "$bm" --arg bn "$bn" \
  -n '{arm:"cado-bwc", cols:$e.cols, weight:$e.weight, seed:$e.seed, core_dim:$e.core_dim,
       core_nonzeros:$e.core_nonzeros, thr:$thr, block_m:($bm|tonumber), block_n:($bn|tonumber),
       solve_s:($wall|tonumber), steps_wct:$steps, verified:$c.verified, rows_ok:$c.rows_ok}'
