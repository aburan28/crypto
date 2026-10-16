#!/usr/bin/env bash
# Run a command with a 1.5 GB address-space limit (msolve reserves ~1 GB of virtual memory at start)
# and enforce the item's 700 MB cap on resident memory by polling VmRSS; print peak RSS (kB) to fd 3.
# Usage: rsscap.sh <timeout_s> cmd...
T=$1; shift
prlimit --as=1500000000 -- "$@" & pid=$!
peak=0; start=$(date +%s); rc=
while kill -0 $pid 2>/dev/null; do
  r=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null); r=${r:-0}
  [ "$r" -gt "$peak" ] && peak=$r
  if [ "$r" -gt 700000 ]; then kill -9 $pid; wait $pid; echo "$peak" >&3; exit 137; fi
  if [ $(( $(date +%s) - start )) -gt $T ]; then kill -9 $pid; wait $pid; echo "$peak" >&3; exit 124; fi
  sleep 0.05
done
wait $pid; rc=$?; echo "$peak" >&3; exit $rc
