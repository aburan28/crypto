#!/usr/bin/env bash
# Parallel lanes for the semi-regular law round.  One `macaulay_dense` process per
# (cell, draw, arm, D), scanning D = 1, 2, ... up to the draw's registered cap
# (d_reg + 2, from runs/<tag>/predictions.jsonl) and stopping at the first refuting D.
# Units are `cell:draw:arm`, claimed with an atomic mkdir lock under runs/<tag>/locks/
# (a lock from a dead pid or an earlier boot is reclaimed; one whose owner line starts
# with `remote` is left alone).  Lines go to runs/<tag>/<cell>.jsonl, schema of
# ic_dense_ladder_20261004; finished units are read from the log, never from the lock.
#
#   THREADS=2 lanes.sh registered 2 K1n61l5:0:rr K1n61l5:0:x4 ...
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
TAG=${1:?tag}; NLANES=${2:?number of lanes}; shift 2
OUT="$HERE/runs/$TAG"
DUMP=${DUMP:-"$HERE/runs/registered/dump"}
PRED=${PRED:-"$HERE/runs/registered/predictions.jsonl"}
BIN="$ROOT/target/release/examples/macaulay_dense"
CPU_LIMIT=${CPU_LIMIT:-86400}     # CPU seconds per process (registered)
MEM_LIMIT=${MEM_LIMIT:-13500000}  # KB of address space; the engine refuses above --mem-gb 12 first
THREADS=${THREADS:-4}
mkdir -p "$OUT/locks"
BOOT=$(cat /proc/sys/kernel/random/boot_id)

cap() {  # cell draw arm -> d_reg + 2 from the registered predictions
  grep -F "\"cell\": \"$1\", \"draw\": $2, \"arm\": \"$3\"," "$PRED" | sed -E 's/.*"d_reg": ([0-9]+).*/\1/' | awk '{print $1 + 2}'
}

measure() {  # cell draw arm
  local cell=$1 draw=$2 arm=$3 log="$OUT/$1.jsonl" file="$DUMP/$1-d$2-$3.sing" top d
  [ -f "$file" ] || { echo "missing $file" >&2; return; }
  top=$(cap "$cell" "$draw" "$arm"); [ -n "$top" ] || { echo "no prediction for $cell $draw $arm" >&2; return; }
  for d in $(seq 1 "$top"); do
    local prev; prev=$(grep -F "\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d," "$log" 2>/dev/null | tail -1)
    if [ -n "$prev" ]; then
      case "$prev" in *'"refuted_at":0,'*) continue;; *) return;; esac
    fi
    local s full wall; s=$(date +%s.%N)
    full=$( (ulimit -t "$CPU_LIMIT" -v "$MEM_LIMIT"; "$BIN" --file "$file" --degree "$d" --threads "$THREADS" --mem-gb 12 2>>"$OUT/$cell.stderr") )
    wall=$(printf "%.3f" "$(echo "$(date +%s.%N) - $s" | bc)")
    if [[ "$full" =~ \"n_vars\":([0-9]+),\"n_cols\":([0-9]+),\"n_rows\":([0-9]+),\"rank\":([0-9]+),\"refuted\":(true|false),\"secs\":([0-9.]+) ]]; then
      local ra=0; [ "${BASH_REMATCH[5]}" = true ] && ra=$d
      echo "{\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d,\"kind\":\"done\",\"refuted_at\":$ra,\"pinned\":0,\"n_vars\":${BASH_REMATCH[1]},\"n_cols\":${BASH_REMATCH[2]},\"n_rows\":${BASH_REMATCH[3]},\"rank\":${BASH_REMATCH[4]},\"ms\":$(printf "%.0f" "$(echo "${BASH_REMATCH[6]} * 1000" | bc)"),\"phase\":1,\"wall\":$wall,\"host\":\"$(hostname)\"}" >> "$log"
      [ "$ra" != 0 ] && return
    elif [[ "$full" =~ \"refused\" ]]; then
      echo "{\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d,\"kind\":\"killed\",\"reason\":\"memory\",\"phase\":1,\"wall\":$wall,\"host\":\"$(hostname)\"}" >> "$log"
      return
    else
      echo "{\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d,\"kind\":\"killed\",\"reason\":\"cpu\",\"phase\":1,\"wall\":$wall,\"host\":\"$(hostname)\"}" >> "$log"
      return
    fi
  done
}

unit_done() {  # refuted, censored, or scanned to the cap
  local key="\"cell\":\"$1\",\"draw\":$2,\"arm\":\"$3\","
  grep -F "$key" "$OUT/$1.jsonl" 2>/dev/null | grep -qE '"refuted_at":[1-9]|"kind":"killed"|"dmax":'"$(cap "$1" "$2" "$3")"','
}

claim() {
  local l="$OUT/locks/$1" b p
  if mkdir "$l" 2>/dev/null; then echo "$BOOT $BASHPID $(hostname)" > "$l/owner"; return 0; fi
  { read -r b p _ < "$l/owner"; } 2>/dev/null || return 1
  [ "$b" = remote ] && return 1
  if [ "$b" != "$BOOT" ] || ! kill -0 "$p" 2>/dev/null; then
    rm -rf "$l"; mkdir "$l" 2>/dev/null || return 1
    echo "$BOOT $BASHPID $(hostname)" > "$l/owner"; return 0
  fi
  return 1
}

lane() {
  local u c d a
  for u in "$@"; do
    IFS=: read -r c d a <<< "$u"
    unit_done "$c" "$d" "$a" && continue
    claim "$c-d$d-$a" || continue
    if unit_done "$c" "$d" "$a"; then rm -rf "$OUT/locks/$c-d$d-$a"; continue; fi
    echo "$(date -u +%FT%TZ) lane $LANE start $u (threads $THREADS, cap $(cap "$c" "$d" "$a"))" >> "$OUT/progress.txt"
    measure "$c" "$d" "$a"
    echo "$(date -u +%FT%TZ) lane $LANE end $u" >> "$OUT/progress.txt"
    rm -rf "$OUT/locks/$c-d$d-$a"
  done
}

for i in $(seq 1 "$NLANES"); do LANE="$(hostname)-$i" lane "$@" & done
wait
echo "$(date -u +%FT%TZ) lanes ($*) all done" >> "$OUT/progress.txt"
