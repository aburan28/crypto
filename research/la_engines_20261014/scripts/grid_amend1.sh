#!/usr/bin/env bash
# Amendment 1 arms (PROTOCOL_AMENDMENT_1.md). Appends to raw/grid_amend1.jsonl.
set -uo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
BENCH=${BENCH_BIN:-/home/user/crypto/target/release/examples/la_engine_bench}
out=$here/raw/grid_amend1.jsonl
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
    cado-full) line=$(timeout $TIMEOUT "$here/scripts/run_cado_full.sh" "$d" 3 "$s" 4 "$logs" 2>/dev/null); rc=$? ;;
    full-bw-pos) line=$(timeout $TIMEOUT "$BENCH" ours --arm full-bw --coeffs positive --cols "$d" --weight 3 --seed "$s"); rc=$? ;;
    full-bl-pos) line=$(timeout $TIMEOUT "$BENCH" ours --arm full-bl --coeffs positive --cols "$d" --weight 3 --seed "$s"); rc=$? ;;
  esac
  end=$(date +%s.%N)
  local el; el=$(awk -v a="$start" -v b="$end" 'BEGIN{print b-a}')
  if [[ $rc -ne 0 || -z $line ]]; then
    line=$(jq -nc --arg arm "$arm" --argjson d "$d" --argjson s "$s" --argjson rc "$rc" --argjson el "$el" \
      '{arm:$arm, cols:$d, weight:3, seed:$s, failed:true, exit:$rc, timeout:($rc==124), elapsed_s:$el}')
    dead[$arm]=$d
  fi
  echo "$line" | jq -c --arg arm "$arm" --argjson el "$el" \
    '. + {arm:$arm, pass:"amend1", process_elapsed_s:$el, utc:(now|todate)}' >> "$out"
  echo "amend1 $arm D=$d seed=$s rc=$rc ${el}s" >&2
}
for d in 1024 2048 4096 8192 16384 32768 65536 131072 262144; do
  for s in 1 2 3; do
    for arm in cado-full full-bw-pos full-bl-pos; do run $arm $d $s; done
  done
done
echo done >&2
