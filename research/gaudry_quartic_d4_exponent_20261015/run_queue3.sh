#!/usr/bin/env bash
# Four workers (CPUs 0-3) take jobs from ${JOBS:-jobs.txt} in order. A job not
# started before the deadline (PROTOCOL.md "Budget") is not run and is
# listed in sessions/NOT_RUN3.txt. Each session has its own lock file, so
# the four do not queue behind each other; all run with --allow-busy.
set -uo pipefail
cd "$(dirname "$0")"
BIN=../../target/release/ecbench
DEADLINE=${DEADLINE:?epoch seconds}
next() { # atomically pop the next job index
  ( flock 9; n=$(cat sessions/.next3); echo $((n + 1)) > sessions/.next3; echo "$n" ) 9>sessions/.next3.lock
}
[ -f sessions/.next3 ] || echo 1 > sessions/.next3
worker() {
  local cpu=$1
  while :; do
    local i; i=$(next)
    local job; job=$(sed -n "${i}p" ${JOBS:-jobs.txt})
    [ -z "$job" ] && return
    if [ "$(date +%s)" -ge "$DEADLINE" ]; then
      echo "$job not started before the deadline" >> sessions/NOT_RUN3.txt; continue
    fi
    echo "$(date -u +%FT%TZ) start $job cpu $cpu" >> sessions/queue.log
    "$BIN" run --spec "specs/SPEC-$job.json" --out "sessions/$job" --cpus inherit \
      --lock "/tmp/ecbench-k4exp-$cpu.lock" --allow-busy > "sessions/$job.log" 2>&1
    echo "$(date -u +%FT%TZ) end $job exit $?" >> sessions/queue.log
  done
}
for c in 0 1 2; do worker $c & done
wait
echo "$(date -u +%FT%TZ) queue 3 done" >> sessions/queue.log
