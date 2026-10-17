#!/usr/bin/env bash
# Thin orchestration of the frozen grid (PROTOCOL.md). Appends one JSON line
# per run to raw/grid.jsonl; a timed-out or failed run gets a line too.
set -uo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
BENCH=${BENCH_BIN:-/home/user/crypto/target/release/examples/la_engine_bench}
out=$here/raw/grid.jsonl
logs=$here/raw/cado_logs
mkdir -p "$logs"
export TMPDIR=${TMPDIR:-/tmp/cb}
mkdir -p "$TMPDIR"
TIMEOUT=1800
declare -A dead   # arm -> smallest D at which it timed out or ran > cap
declare -A cap=([igauss]=600 [sgauss]=600)
declare -A maxd=([igauss]=8192 [sgauss]=16384)

run() { # pass label arm D m seed threads
  local pass=$1 arm=$2 d=$3 w=$4 s=$5 t=$6 key="$1/$2/$4/$6"
  if [[ -n ${dead[$key]:-} && $d -ge ${dead[$key]} ]]; then return; fi
  if [[ -n ${maxd[$arm]:-} && $d -gt ${maxd[$arm]} ]]; then return; fi
  local start end line rc
  start=$(date +%s.%N)
  if [[ $arm == cado-bwc ]]; then
    local thr=2x2; [[ $t == 1 ]] && thr=1x1
    line=$(timeout $TIMEOUT "$here/scripts/run_cado.sh" "$d" "$w" "$s" "$thr" 4 4 "$logs" 2>/dev/null); rc=$?
  else
    local a=$arm; [[ $arm == core-bw-aa ]] && a=core-bw
    line=$(RAYON_NUM_THREADS=$t timeout $TIMEOUT "$BENCH" ours --arm "$a" --cols "$d" --weight "$w" --seed "$s"); rc=$?
  fi
  end=$(date +%s.%N)
  local el; el=$(awk -v a="$start" -v b="$end" 'BEGIN{print b-a}')
  if [[ $rc -ne 0 || -z $line ]]; then
    line=$(jq -nc --arg arm "$arm" --argjson d "$d" --argjson w "$w" --argjson s "$s" --argjson rc "$rc" --argjson el "$el" \
      '{arm:$arm, cols:$d, weight:$w, seed:$s, failed:true, exit:$rc, timeout:($rc==124), elapsed_s:$el}')
    dead[$key]=$d
  elif [[ -n ${cap[$arm]:-} ]] && awk -v e="$el" -v c="${cap[$arm]}" 'BEGIN{exit !(e>c)}'; then
    dead[$key]=$((d * 2))
  fi
  echo "$line" | jq -c --arg pass "$pass" --argjson t "$t" --arg host "$(hostname)" --argjson el "$el" \
    '. + {pass:$pass, threads_requested:$t, process_elapsed_s:$el, utc:(now|todate)}' >> "$out"
  echo "$pass $arm D=$d m=$w seed=$s t=$t rc=$rc ${el}s" >&2
}

# Main sweep, m = 3, four threads.
for d in 1024 2048 4096 8192 16384 32768 65536 131072 262144; do
  for s in 1 2 3; do
    for arm in igauss sgauss full-bw full-bl core-bw cado-bwc; do run main $arm $d 3 $s 4; done
    [[ $s == 1 ]] && run main core-bw-aa $d 3 1 4
  done
done
# Row weight at D = 16,384.
for w in 2 4; do for s in 1 2 3; do
  for arm in igauss sgauss full-bw full-bl core-bw cado-bwc; do run weight $arm 16384 $w $s 4; done
done; done
# One thread, engine against engine.
for d in 1024 2048 4096 8192 16384 32768 65536; do for s in 1 2 3; do
  for arm in core-bw cado-bwc; do run onethread $arm $d 3 $s 1; done
done; done
echo done >&2
