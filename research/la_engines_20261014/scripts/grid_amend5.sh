#!/usr/bin/env bash
# Amendment 5 arms (PROTOCOL_AMENDMENT_5.md). Appends to raw/grid_amend5.jsonl.
set -uo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
BENCH=${BENCH_BIN:-/home/user/crypto/target/release/examples/la_engine_bench}
out=$here/raw/grid_amend5.jsonl
logs=$here/raw/cado_full_logs
mkdir -p "$logs"
export TMPDIR=${TMPDIR:-/tmp/cb}
mkdir -p "$TMPDIR"
TIMEOUT=1800
declare -A dead
run() { # arm D seed
  local arm=$1 d=$2 s=$3 line rc
  [[ -n ${dead[$arm]:-} && $d -ge ${dead[$arm]} ]] && return
  local base=(ours --cols "$d" --weight 3 --seed "$s")
  local start end; start=$(date +%s.%N)
  case $arm in
    legacy-bw) line=$(timeout $TIMEOUT "$BENCH" "${base[@]}" --arm full-bw); rc=$? ;;
    cado-bw) line=$(timeout $TIMEOUT "$BENCH" "${base[@]}" --arm full-bw --filter-preset cado); rc=$? ;;
    cado-bl) line=$(timeout $TIMEOUT "$BENCH" "${base[@]}" --arm full-bl --filter-preset cado); rc=$? ;;
    legacy-bw-pos) line=$(timeout $TIMEOUT "$BENCH" "${base[@]}" --arm full-bw --coeffs positive); rc=$? ;;
    cado-bw-pos) line=$(timeout $TIMEOUT "$BENCH" "${base[@]}" --arm full-bw --coeffs positive --filter-preset cado); rc=$? ;;
    cado-bw-pos-d*) dd=${arm##*-d}
      line=$(timeout $TIMEOUT "$BENCH" "${base[@]}" --arm full-bw --coeffs positive --filter-preset cado \
        --filter-strategy cado --target-density "$dd"); rc=$? ;;
    cado-nfs) line=$(BENCH_BIN=$BENCH timeout $TIMEOUT "$here/scripts/run_cado_full_fb.sh" "$d" 3 "$s" 4 "$logs" 2>/dev/null); rc=$? ;;
  esac
  end=$(date +%s.%N)
  local el; el=$(awk -v a="$start" -v b="$end" 'BEGIN{print b-a}')
  if [[ $rc -ne 0 || -z $line ]]; then
    line=$(jq -nc --arg arm "$arm" --argjson d "$d" --argjson s "$s" --argjson rc "$rc" --argjson el "$el" \
      '{arm:$arm, cols:$d, weight:3, seed:$s, failed:true, exit:$rc, timeout:($rc==124), elapsed_s:$el}')
    dead[$arm]=$d
  fi
  echo "$line" | jq -c --arg arm "$arm" --argjson el "$el" \
    '. + {arm:$arm, pass:"amend5", process_elapsed_s:$el, utc:(now|todate)}' >> "$out"
  echo "amend5 $arm D=$d seed=$s rc=$rc ${el}s" >&2
}
for d in 1024 2048 4096 8192 16384 32768 65536 131072 262144; do
  for s in 1 2 3; do
    for arm in legacy-bw cado-bw legacy-bw-pos cado-bw-pos cado-nfs; do
      run $arm $d $s
    done
  done
done
echo done >&2
