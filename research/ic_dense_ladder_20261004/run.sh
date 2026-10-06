#!/usr/bin/env bash
# Dense-engine ladder: one `macaulay_dense` process per (cell, draw, arm, D), serial,
# all four cores inside the process.  Lines go to runs/<tag>/<cell>.jsonl in the
# schema of ic_gb_ladder_20261003 so `gb_ladder_analyze` reads them unchanged; a run
# resumes by skipping lines already present.
#
#   research/ic_dense_ladder_20261004/run.sh registered
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
TAG=${1:?tag}
OUT="$HERE/runs/$TAG"
DUMP="$ROOT/research/ic_gb_ladder_20261003/runs/registered/dump"   # the same exported systems
BIN="$ROOT/target/release/examples/macaulay_dense"
CPU_LIMIT=${CPU_LIMIT:-43200}     # CPU seconds per process (four threads: ~3 h wall)
MEM_LIMIT=${MEM_LIMIT:-13500000}  # KB of address space; the engine refuses above --mem-gb 12 first
THREADS=${THREADS:-4}
mkdir -p "$OUT"
d0() { case $1 in rr) echo 1;; x4) echo 1;; ctrl) echo 1;; esac; }  # from 1 so a constant equation reads 1 (triv)
dcap() { case $1 in rr) echo 9;; x4) echo 11;; ctrl) echo 7;; esac; }

measure() {  # cell draw arm
  local cell=$1 draw=$2 arm=$3 log="$OUT/$1.jsonl" file="$DUMP/$1-d$2-$3.sing"
  [ -f "$file" ] || { echo "missing $file" >&2; return; }
  local d
  for d in $(seq "$(d0 "$arm")" "$(dcap "$arm")"); do
    local prev; prev=$(grep "\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d," "$log" 2>/dev/null | tail -1)
    if [ -n "$prev" ]; then
      case "$prev" in *'"refuted_at":0,'*) continue;; *) return;; esac
    fi
    local s res full; s=$(date +%s.%N)
    full=$( (ulimit -t "$CPU_LIMIT" -v "$MEM_LIMIT"; "$BIN" --file "$file" --degree "$d" --threads "$THREADS" --mem-gb 12 2>>"$OUT/$cell.stderr") )
    local wall; wall=$(printf "%.3f" "$(echo "$(date +%s.%N) - $s" | bc)")
    if [[ "$full" =~ \"n_vars\":([0-9]+),\"n_cols\":([0-9]+),\"n_rows\":([0-9]+),\"rank\":([0-9]+),\"refuted\":(true|false),\"secs\":([0-9.]+) ]]; then
      local ra=0; [ "${BASH_REMATCH[5]}" = true ] && ra=$d
      echo "{\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d,\"kind\":\"done\",\"refuted_at\":$ra,\"pinned\":0,\"n_vars\":${BASH_REMATCH[1]},\"n_cols\":${BASH_REMATCH[2]},\"n_rows\":${BASH_REMATCH[3]},\"rank\":${BASH_REMATCH[4]},\"ms\":$(printf "%.0f" "$(echo "${BASH_REMATCH[6]} * 1000" | bc)"),\"phase\":1,\"wall\":$wall}" >> "$log"
      [ "$ra" != 0 ] && return
    elif [[ "$full" =~ \"refused\" ]]; then
      echo "{\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d,\"kind\":\"killed\",\"reason\":\"memory\",\"phase\":1,\"wall\":$wall}" >> "$log"
      return
    else
      echo "{\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d,\"kind\":\"killed\",\"reason\":\"cpu\",\"phase\":1,\"wall\":$wall}" >> "$log"
      return
    fi
  done
}

cell_draws() {  # cell arm -> draws measured on that arm (rootless), from the catalogue
  grep "\"arm\":\"$2\"" "$DUMP/$1.catalogue.jsonl" | grep -v '"roots":[1-9]' | sed -E 's/.*"draw":([0-9]+).*/\1/'
}

[ "${LIB_ONLY:-}" = 1 ] && return 0   # sourced by lanes.sh for the functions above

if [ "${SMOKE:-}" = 1 ]; then  # plumbing check on the l = 2 cell only
  for arm in rr x4 ctrl; do for draw in $(cell_draws K1n17l2 $arm); do measure K1n17l2 $draw $arm; done; done
  mkdir -p "$OUT/dump" && cp "$DUMP/K1n17l2.catalogue.jsonl" "$OUT/dump/"; exit 0
fi

# Registered order: calibration cells, then the l = 6 rr arms, then the l = 6 x4 arms,
# then the controls, so the decisive readings come first.
for cell in K1n17l2 K1n17l3 K1n17l4 K1n17l5; do
  for arm in rr x4 ctrl; do for draw in $(cell_draws $cell $arm); do measure $cell $draw $arm; done; done
  echo "$(date -u +%FT%TZ) finished $cell" >> "$OUT/progress.txt"
done
for arm in rr x4 ctrl; do
  for cell in K1n19l6 K0n19l6; do
    for draw in $(cell_draws $cell $arm); do measure $cell $draw $arm; done
    echo "$(date -u +%FT%TZ) finished $cell $arm" >> "$OUT/progress.txt"
  done
done
echo "$(date -u +%FT%TZ) all done" >> "$OUT/progress.txt"
