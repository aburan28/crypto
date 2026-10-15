#!/usr/bin/env bash
# Four workers (CPUs 0-3) take jobs from jobs.txt in order. A job not
# started before the deadline (PROTOCOL.md "Budget") is not run and is
# listed in sessions/NOT_RUN.txt. Each session has its own lock file, so
# the four do not queue behind each other; all run with --allow-busy.
set -uo pipefail
cd "$(dirname "$0")"
BIN=../../target/release/ecbench
DEADLINE=${DEADLINE:?epoch seconds}
next() { # atomically pop the next job index
  ( flock 9; n=$(cat sessions/.next); echo $((n + 1)) > sessions/.next; echo "$n" ) 9>sessions/.next.lock
}
[ -f sessions/.next ] || echo 1 > sessions/.next
worker() {
  local cpu=$1
  while :; do
    local i; i=$(next)
    local job; job=$(sed -n "${i}p" jobs.txt)
    [ -z "$job" ] && return
    if [ "$(date +%s)" -ge "$DEADLINE" ]; then
      echo "$job not started before the deadline" >> sessions/NOT_RUN.txt; continue
    fi
    echo "$(date -u +%FT%TZ) start $job cpu $cpu" >> sessions/queue.log
    "$BIN" run --spec "specs/SPEC-$job.json" --out "sessions/$job" --cpus "$cpu" \
      --lock "/tmp/ecbench-k4exp-$cpu.lock" --allow-busy > "sessions/$job.log" 2>&1
    echo "$(date -u +%FT%TZ) end $job exit $?" >> sessions/queue.log
  done
}
for c in 0 1 2 3; do worker $c & done
wait
echo "$(date -u +%FT%TZ) queue done" >> sessions/queue.log
