#!/bin/sh
# Thin orchestration, second half of PROTOCOL.md: the panels A, B, C, E, appended to the run whose
# reproduction arms (run.sh with STOP_AFTER_REPRO=1) passed compare.sh. Same binary, same source.
# Usage: sh run_panels.sh <run-id> <bin-dir> [first-arm]   (from the study directory)
# first-arm (A, B, C or E) resumes after an interruption: earlier arms with a complete record are kept.
set -u
RUN=$1; BIN=$2; FIRST=${3:-A}; OUT=results/run-$RUN; ENV=results/run-$RUN.env.txt
echo "panels started $(date -u +%Y%m%dT%H%MZ) from arm $FIRST; source sha256 $(shasum -a 256 cover_census_n.rs | cut -d' ' -f1); commit $(git rev-parse HEAD); $(uptime | sed 's/.*load/load/')" >> "$ENV"
[ -x "$BIN/ccn" ] || rustc --edition 2021 -O cover_census_n.rs -o "$BIN/ccn" || { echo "compile FAILED" >> "$ENV"; exit 1; }
arm() { name=$1; shift; s=$(date +%s); "$BIN/ccn" "$@" 6 > "$OUT/$name.jsonl" 2> "$OUT/$name.stderr"; st=$?; e=$(date +%s)
  echo "$name: ccn $* 6 -> exit $st, wall_s $((e - s)), $(uptime | sed 's/.*load/load/')" >> "$ENV"; }
go=0
for a in A B C E; do
  [ "$a" = "$FIRST" ] && go=1
  [ $go = 1 ] || continue
  case $a in
    A) arm A-n3-q11 11 3 exhaustive cosets ;;
    B) arm B-n3-q13 13 3 exhaustive cosets ;;
    C) arm C-n3-q17-sample200 17 3 sample:200 cosets ;;
    E) arm E-n5-q5-sample60 5 5 sample:60 cosets ;;
  esac
done
echo "uptime at end: $(uptime)" >> "$ENV"
echo DONE >> "$ENV"
