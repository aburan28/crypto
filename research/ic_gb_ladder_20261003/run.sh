#!/usr/bin/env bash
# External-engine refutation-degree ladder.  Three lanes, one per CPU, cheap cells
# first.  Every (cell, draw, arm, D) is one Singular process under a CPU and memory
# limit; its line is appended to runs/<tag>/<cell>.jsonl, so a lane resumes by
# skipping lines already present (`run.sh <tag> --resume`).
#
#   research/ic_gb_ladder_20261003/run.sh registered [--resume]
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
TAG=${1:?tag}
OUT="$HERE/runs/$TAG"
BIN="$ROOT/target/release/examples/rr_degree_ladder"
SEED=20260930
CPU_LIMIT=${CPU_LIMIT:-3600}      # CPU seconds per Singular process
MEM_LIMIT=${MEM_LIMIT:-4500000}   # KB of address space per Singular process
mkdir -p "$OUT/dump"
# arm: first degree tried, degree cap
d0() { case $1 in rr) echo 3;; x4) echo 4;; ctrl) echo 3;; esac; }
dcap() { case $1 in rr) echo 12;; x4) echo 14;; ctrl) echo 12;; esac; }

dump_cell() {  # a n ell
  local cell="K$1n$2l$3" cat="$OUT/dump/K$1n$2l$3.catalogue.jsonl"
  if [ ! -s "$cat" ]; then
    "$BIN" --a "$1" --n "$2" --ell "$3" --unsat 4 --seed $SEED --dump-dir "$OUT/dump" --out "$cat" 2>>"$OUT/$cell.stderr"
  fi
}

measure() {  # cpu cell draw arm file
  local cpu=$1 cell=$2 draw=$3 arm=$4 file=$5 log="$OUT/$2.jsonl"
  local d
  for d in $(seq "$(d0 "$arm")" "$(dcap "$arm")"); do
    if grep -q "\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d," "$log" 2>/dev/null; then
      local prev; prev=$(grep "\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d," "$log" | tail -1)
      case "$prev" in
        *'"kind":"killed","reason":"memory","phase":1,'*) if [ "${PHASE:-1}" = 2 ]; then :; else return; fi;;
        *'"kind":"killed"'*) return;;
        *'"refuted_at":0,"pinned":'*) if [[ "$prev" =~ \"pinned\":([0-9]+),\"n_vars\":([0-9]+) ]] && [ "${BASH_REMATCH[1]}" = "${BASH_REMATCH[2]}" ]; then return; fi; continue;;
        *) return;;
      esac
    fi
    local s res code
    s=$(date +%s.%N)
    local full
    full=$( (ulimit -t "$CPU_LIMIT" -v "$MEM_LIMIT"; taskset -c "$cpu" Singular -q -c "string SINGFILE=\"$file\"; int DMAX=$d;" "$HERE/refute.sing" 2>&1) )
    res=$(printf '%s\n' "$full" | grep '^RESULT')
    code=$?
    local wall; wall=$(printf "%.3f" "$(echo "$(date +%s.%N) - $s" | bc)")
    if [[ "$res" =~ refuted_at=([0-9]+)\ pinned=([0-9]+)\ nvars=([0-9]+)\ gb_size=([0-9]+)\ ms=([0-9]+) ]]; then
      echo "{\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d,\"kind\":\"done\",\"refuted_at\":${BASH_REMATCH[1]},\"pinned\":${BASH_REMATCH[2]},\"n_vars\":${BASH_REMATCH[3]},\"gb_size\":${BASH_REMATCH[4]},\"ms\":${BASH_REMATCH[5]},\"phase\":${PHASE:-1},\"wall\":$wall}" >> "$log"
      [ "${BASH_REMATCH[1]}" != 0 ] && return
      [ "${BASH_REMATCH[2]}" = "${BASH_REMATCH[3]}" ] && return
    else
      local reason=cpu
      case "$full" in *"no more memory"*) reason=memory;; esac
      echo "{\"cell\":\"$cell\",\"draw\":$draw,\"arm\":\"$arm\",\"dmax\":$d,\"kind\":\"killed\",\"reason\":\"$reason\",\"phase\":${PHASE:-1},\"wall\":$wall}" >> "$log"
      return
    fi
  done
}

lane() {  # cpu "a n ell" ...
  local cpu=$1; shift
  local spec
  for spec in "$@"; do
    set -- $spec
    local a=$1 n=$2 ell=$3 cell="K$1n$2l$3"
    dump_cell "$a" "$n" "$ell"
    local cat="$OUT/dump/$cell.catalogue.jsonl"
    # rr then x4 on every rootless draw of that arm, in draw order; controls after
    local arm draw file
    for arm in rr x4 ctrl; do
      while read -r draw file; do
        measure "$cpu" "$cell" "$draw" "$arm" "$OUT/dump/$file"
      done < <(grep "\"arm\":\"$arm\"" "$cat" | grep -v '"roots":[1-9]' | sed -E 's/.*"draw":([0-9]+).*"dump":"([^"]+)".*/\1 \2/')
    done
    echo "$(date -u +%FT%TZ) lane $cpu finished $cell" >> "$OUT/progress.txt"
  done
}

if [ "${PHASE:-1}" = 2 ]; then  # §3 phase 2: memory-killed draws retried once, alone, at 12 GB
  MEM_LIMIT=12000000
  # smallest ℓ first, rr before x4 before ctrl, so the decisive rungs get the memory first
  grep -h '"kind":"killed","reason":"memory","phase":1,' "$OUT"/K*.jsonl | sed -E 's/.*"cell":"([^"]+)","draw":([0-9]+),"arm":"([^"]+)".*/\1 \2 \3/' | sort -u |
  awk '{ split($1, c, "l"); a = ($3 == "rr") ? 0 : ($3 == "x4") ? 1 : 2; print c[2], a, $0 }' | sort -k1,1n -k2,2n -k3 | cut -d" " -f3- |
  while read -r cell draw arm; do
    measure 1 "$cell" "$draw" "$arm" "$OUT/dump/$cell-d$draw-$arm.sing"
  done
  echo "$(date -u +%FT%TZ) phase 2 done" >> "$OUT/progress.txt"
  exit 0
fi

if [ "${SMOKE:-}" = 1 ]; then  # plumbing check on one tiny cell, never a registered cell
  lane 1 "1 17 2"; wait; exit 0
fi

# Registered lanes (cheap first).  Shared draws with ic_rr_norm_ladder_20260930: same
# seed, same generator, so K1n17 l2-l5 and K1n19 l6 are the same systems.
lane 1 "1 17 2" "1 17 3" "1 17 4" "1 17 5" "1 23 7" &
lane 2 "1 19 6" "1 29 8" &
lane 3 "0 19 6" "0 31 8" &
wait
echo "$(date -u +%FT%TZ) all lanes done" >> "$OUT/progress.txt"
