#!/usr/bin/env bash
# Parallel lanes for the dense ladder, added 2026-10-05 after the serial runner had
# finished the calibration cells and the l = 6 rr arms.  Same `measure` as run.sh (sourced),
# same per-process budget, same log schema; what changes is that several units run at once.
#
# A unit is `cell:draw:arm`.  Each lane walks the unit list in order and claims a unit with
# an atomic mkdir under runs/<tag>/locks/.  A lock whose owner is gone (another boot, or a
# dead pid) is reclaimed, so a host restart never strands a unit; a lock whose owner line
# starts with `remote` is held by another host and is never taken here.  Whether a unit is
# finished is read from the log, never from the lock.
#
#   THREADS=2 lanes.sh registered 2 K0n19l6:0:x4 K0n19l6:2:x4 K1n19l6:0:ctrl ...
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
TAG=${1:?tag}; NLANES=${2:?number of lanes}; shift 2
LIB_ONLY=1 source "$HERE/run.sh" "$TAG"
BOOT=$(cat /proc/sys/kernel/random/boot_id)
LOCKS="$OUT/locks"; mkdir -p "$LOCKS"

unit_done() {  # cell draw arm: refuted, censored, or scanned to the cap
  local log="$OUT/$1.jsonl" key="\"cell\":\"$1\",\"draw\":$2,\"arm\":\"$3\","
  grep -F "$key" "$log" 2>/dev/null | grep -qE '"refuted_at":[1-9]|"kind":"killed"|"dmax":'"$(dcap "$3")"','
}

claim() {  # lock name -> 0 if this lane now holds it
  local l="$LOCKS/$1" b p
  if mkdir "$l" 2>/dev/null; then echo "$BOOT $BASHPID $(hostname)" > "$l/owner"; return 0; fi
  { read -r b p _ < "$l/owner"; } 2>/dev/null || return 1   # owner not written yet: held
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
    if unit_done "$c" "$d" "$a"; then rm -rf "$LOCKS/$c-d$d-$a"; continue; fi  # finished meanwhile
    echo "$(date -u +%FT%TZ) lane $LANE start $u (threads $THREADS)" >> "$OUT/progress.txt"
    measure "$c" "$d" "$a"
    echo "$(date -u +%FT%TZ) lane $LANE end $u" >> "$OUT/progress.txt"
    rm -rf "$LOCKS/$c-d$d-$a"
  done
}

for i in $(seq 1 "$NLANES"); do LANE="$HOSTNAME-$i" lane "$@" & done
wait
echo "$(date -u +%FT%TZ) lanes ($*) all done" >> "$OUT/progress.txt"
