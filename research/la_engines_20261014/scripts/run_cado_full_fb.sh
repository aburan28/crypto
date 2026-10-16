#!/usr/bin/env bash
# Thin orchestration: CADO-NFS's own DL filter and solver on the raw
# relation system (purge -> merge-dl -> replay-dl -> bwc), then native
# right-hand-side assembly and log reconstruction, each step timed.
# Amendment 2 variant: if bwc's internal balancing aborts ("expects a column
# permutation replicated on rows ... pass --reorder columns to mf_bal"),
# rebalance with `mf_bal --reorder columns` into a fresh wdir and rerun bwc
# at the same thread split; record reorder_retry. The timed bwc step is
# mf_bal plus the rerun (the aborted first attempt is excluded and logged).
# usage: run_cado_full_fb.sh <cols> <weight> <seed> <threads> <outdir>
set -euo pipefail
cols=$1 weight=$2 seed=$3 t=$4 out=$5
CADO=${CADO_BUILD:-/home/user/ext/cado-nfs/build}
BENCH=${BENCH_BIN:-/home/user/crypto/target/release/examples/la_engine_bench}
P=4611686018427387847
thr=2x2; [[ $t == 1 ]] && thr=1x1
dir=$(mktemp -d "${TMPDIR:-/tmp}/cado-full.XXXXXX")
trap 'rm -rf "$dir"' EXIT
common=(--cols "$cols" --weight "$weight" --seed "$seed" --coeffs positive --dir "$dir")
"$BENCH" cado-export-rels "${common[@]}" > "$dir/export.json"
now() { date +%s.%N; }
lap() { awk -v a="$1" -v b="$2" 'BEGIN{print b-a}'; }
t0=$(now)
"$CADO/filter/purge" -out "$dir/purged" -outdel "$dir/deleted" -keep 160 -t "$t" "$dir/rels" > "$dir/purge.log" 2>&1
t1=$(now)
"$CADO/filter/merge-dl" -mat "$dir/purged" -out "$dir/history" -t "$t" -target_excess 0 > "$dir/merge.log" 2>&1
t2=$(now)
"$CADO/filter/replay-dl" -purged "$dir/purged" -his "$dir/history" -out "$dir/small.bin" \
  -index "$dir/small.index" -ideals "$dir/small.ideal" > "$dir/replay.log" 2>&1
t3=$(now)
"$BENCH" cado-rhs "${common[@]}" > "$dir/rhs.json"
mat=$(ls "$dir"/small*.bin | grep -vE '\.(rw|cw)\.bin$' | head -1)
[[ -f ${mat%.bin}.rw.bin ]] || "$CADO/linalg/bwc/mf_scan2" --withcoeffs "$mat" > "$dir/scan.log" 2>&1
t4=$(now)
"$CADO/linalg/bwc/bwc.pl" :complete prime=$P matrix="$mat" rhs="$dir/rhs.txt" \
  nullspace=right mpi=1x1 thr="$thr" m=4 n=4 wdir="$dir/wdir" > "$dir/bwc.log" 2>&1 || bwc_rc=$?
fallback=false
if [[ ${bwc_rc:-0} -ne 0 ]]; then
  grep -q "column permutation replicated on rows" "$dir/bwc.log" || { tail -30 "$dir/bwc.log" >&2; exit 1; }
  fallback=true; rm -rf "$dir/wdir"; mv "$dir/bwc.log" "$dir/bwc.first.log"
  t4=$(now)
  base=$(basename "${mat%.bin}"); nhv=${thr}
  mkdir -p "$dir/wdir/$base.$nhv"
  "$CADO/linalg/bwc/mf_bal" ${thr%x*} ${thr#*x} "$mat" --withcoeffs --reorder columns \
    --out "$dir/wdir/$base.$nhv/$base.$nhv.bin" > "$dir/mfbal.log" 2>&1
  "$CADO/linalg/bwc/bwc.pl" :complete prime=$P matrix="$mat" rhs="$dir/rhs.txt" \
    nullspace=right mpi=1x1 thr="$thr" m=4 n=4 wdir="$dir/wdir" > "$dir/bwc.log" 2>&1 \
    || { tail -30 "$dir/bwc.log" >&2; exit 1; }
fi
t5=$(now)
sol=$(ls "$dir"/wdir/K.sols* | grep -v txt | head -1)
"$BENCH" cado-reconstruct "${common[@]}" --sol "$sol" > "$dir/recon.json"
mkdir -p "$out"
tag="cadofullfb-c${cols}-w${weight}-s${seed}-t${t}"
for f in purge merge replay bwc; do cp "$dir/$f.log" "$out/$tag.$f.log"; done

jq -nc --argjson r "$(cat "$dir/rhs.json")" --argjson c "$(cat "$dir/recon.json")" \
  --argjson cols "$cols" --argjson w "$weight" --argjson s "$seed" --argjson t "$t" \
  --arg purge "$(lap $t0 $t1)" --arg merge "$(lap $t1 $t2)" --arg replay "$(lap $t2 $t3)" \
  --argjson fb "$fallback" --arg thr "$thr" --arg bwc "$(lap $t4 $t5)" --arg nnz "$(grep -oE 'Weight of the sparse submatrix: [0-9]+' "$dir/replay.log" | awk '{print $NF}')" \
  '{arm:"cado-full", cols:$cols, weight:$w, seed:$s, coeffs:"positive", threads:$t, reorder_retry:$fb, bwc_thr:$thr,
    small_dim:$r.small_rows, small_nonzeros:($nnz|tonumber? // null),
    purge_s:($purge|tonumber), merge_s:($merge|tonumber), replay_s:($replay|tonumber),
    rhs_s:$r.rhs_s, bwc_s:($bwc|tonumber), reconstruct_s:$c.reconstruct_s,
    solve_s:(($purge|tonumber)+($merge|tonumber)+($replay|tonumber)+$r.rhs_s+($bwc|tonumber)+$c.reconstruct_s),
    determined:$c.determined, verified:$c.verified}'
