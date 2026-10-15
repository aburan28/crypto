#!/usr/bin/env bash
# Amendment 3: re-run each contended cell once on the idle host.
set -uo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
MAIN=/home/user/crypto/target/release/examples/la_engine_bench
A1=/tmp/cbt/release/examples/la_engine_bench
out=$here/raw/rerun_contended.jsonl
export TMPDIR=/tmp/cb
while read -r r; do
  pass=$(jq -r .pass <<<"$r"); arm=$(jq -r .arm <<<"$r"); d=$(jq -r .cols <<<"$r")
  w=$(jq -r .weight <<<"$r"); s=$(jq -r .seed <<<"$r"); t=4; [[ $pass == onethread ]] && t=1
  start=$(date +%s.%N)
  case $arm in
    cado-bwc) thr=2x2; [[ $t == 1 ]] && thr=1x1
      line=$(BENCH_BIN=$MAIN "$here/scripts/run_cado.sh" "$d" "$w" "$s" "$thr" 4 4 "$here/raw/cado_logs" 2>/dev/null) ;;
    cado-full) line=$(BENCH_BIN=$A1 "$here/scripts/run_cado_full.sh" "$d" 3 "$s" 4 "$here/raw/cado_full_logs" 2>/dev/null) ;;
    full-bw-pos) line=$("$A1" ours --arm full-bw --coeffs positive --cols "$d" --weight 3 --seed "$s") ;;
    full-bl-pos) line=$("$A1" ours --arm full-bl --coeffs positive --cols "$d" --weight 3 --seed "$s") ;;
    *) line=$(RAYON_NUM_THREADS=$t "$MAIN" ours --arm "$arm" --cols "$d" --weight "$w" --seed "$s") ;;
  esac
  rc=$?
  el=$(awk -v a="$start" -v b="$(date +%s.%N)" 'BEGIN{print b-a}')
  [[ -z $line ]] && line=$(jq -nc '{failed:true}')
  echo "$line" | jq -c --arg pass "$pass" --arg arm "$arm" --argjson t "$t" --argjson el "$el" --argjson rc "$rc" \
    '. + {pass:$pass, arm:$arm, threads_requested:$t, process_elapsed_s:$el, exit:$rc, rerun:true, utc:(now|todate)}' >> "$out"
  echo "rerun $pass $arm D=$d m=$w s=$s rc=$rc ${el}s" >&2
done < "$here/raw/contended_runs.jsonl"
echo done >&2
