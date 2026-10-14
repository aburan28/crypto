#!/bin/sh
# Thin orchestration for PROTOCOL.md: compile once, then run the reproduction checks and the panels
# one process at a time. Usage: sh run.sh <run-id> <bin-dir>   (from the study directory)
set -u
RUN=$1; BIN=$2; OUT=results/run-$RUN; ENV=results/run-$RUN.env.txt
mkdir -p "$OUT" "$BIN"
{
  echo "run: $RUN"
  echo "commit (protocol frozen at): $(git rev-parse HEAD)"
  echo "source sha256: $(shasum -a 256 cover_census_n.rs | cut -d' ' -f1)"
  echo "protocol sha256: $(shasum -a 256 PROTOCOL.md | cut -d' ' -f1)"
  echo "rustc: $(rustc --version)"
  echo "host: $(uname -sr) $(uname -m), $(sysctl -n machdep.cpu.brand_string), $(sysctl -n hw.ncpu) cpus"
  echo "threads: 6"
  echo "uptime at start: $(uptime)"
} > "$ENV"
rustc --edition 2021 -O cover_census_n.rs -o "$BIN/ccn" || { echo "compile FAILED" >> "$ENV"; exit 1; }
arm() { name=$1; shift; s=$(date +%s); "$BIN/ccn" "$@" 6 > "$OUT/$name.jsonl" 2> "$OUT/$name.stderr"; st=$?; e=$(date +%s)
  echo "$name: ccn $* 6 -> exit $st, wall_s $((e - s)), $(uptime | sed 's/.*load/load/')" >> "$ENV"; }
arm pilot-n3-q3-cosets 3 3 exhaustive cosets
arm pilot-n3-q3-cocycles 3 3 exhaustive cocycles
arm R1-n3-q5 5 3 exhaustive cosets
arm R1-n3-q7 7 3 exhaustive cosets
arm R2-n3-q11-sample120 11 3 sample:120 cosets
arm D-n5-q3-cosets 3 5 exhaustive cosets
arm R3-n5-q3-cocycles 3 5 exhaustive cocycles
echo "reproduction arms done" >> "$ENV"
[ "${STOP_AFTER_REPRO:-0}" = 1 ] && exit 0
arm A-n3-q11 11 3 exhaustive cosets
arm B-n3-q13 13 3 exhaustive cosets
arm C-n3-q17-sample200 17 3 sample:200 cosets
arm E-n5-q5-sample60 5 5 sample:60 cosets
echo "uptime at end: $(uptime)" >> "$ENV"
echo DONE >> "$ENV"
