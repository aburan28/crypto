#!/usr/bin/env bash
# Amendment 2 arms (PROTOCOL_AMENDMENT_2.md). Appends to raw/grid_amend2.jsonl.
set -uo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
CADO_BENCH=${CADO_BENCH_BIN:-/tmp/cbt/release/examples/la_engine_bench}
TUNE_BENCH=${TUNE_BENCH_BIN:-/tmp/cbt2/release/examples/la_engine_bench}
out=$here/raw/grid_amend2.jsonl
logs=$here/raw/cado_full_logs
mkdir -p "$logs"
export TMPDIR=${TMPDIR:-/tmp/cb}
mkdir -p "$TMPDIR"
TIMEOUT=1800
declare -A dead
run() { # arm D seed
  local arm=$1 d=$2 s=$3 line rc
  [[ -n ${dead[$arm]:-} && $d -ge ${dead[$arm]} ]] && return
  local start end; start=$(date +%s.%N)
  case $arm in
    cado-full-fb) line=$(BENCH_BIN=$CADO_BENCH timeout $TIMEOUT "$here/scripts/run_cado_full_fb.sh" "$d" 3 "$s" 4 "$logs" 2>/dev/null); rc=$? ;;
    full-bw-m32r100) line=$(timeout $TIMEOUT "$TUNE_BENCH" ours --arm full-bw --coeffs positive --cols "$d" --weight 3 --seed "$s" --merge-max-weight 32 --max-row-weight 100); rc=$? ;;
    full-bw-m128r200) line=$(timeout $TIMEOUT "$TUNE_BENCH" ours --arm full-bw --coeffs positive --cols "$d" --weight 3 --seed "$s" --merge-max-weight 128 --max-row-weight 200); rc=$? ;;
  esac
  end=$(date +%s.%N)
  local el; el=$(awk -v a="$start" -v b="$end" 'BEGIN{print b-a}')
  if [[ $rc -ne 0 || -z $line ]]; then
    line=$(jq -nc --arg arm "$arm" --argjson d "$d" --argjson s "$s" --argjson rc "$rc" --argjson el "$el" \
      '{arm:$arm, cols:$d, weight:3, seed:$s, failed:true, exit:$rc, timeout:($rc==124), elapsed_s:$el}')
    dead[$arm]=$d
  fi
  echo "$line" | jq -c --arg arm "$arm" --argjson el "$el" \
    '. + {arm:$arm, pass:"amend2", process_elapsed_s:$el, utc:(now|todate)}' >> "$out"
  echo "amend2 $arm D=$d seed=$s rc=$rc ${el}s" >&2
}
for d in 1024 2048 4096 8192 16384 32768 65536 131072 262144; do
  for s in 1 2 3; do
    run cado-full-fb $d $s
    if [[ $d -ge 16384 ]]; then run full-bw-m32r100 $d $s; run full-bw-m128r200 $d $s; fi
  done
done
echo done >&2
