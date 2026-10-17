#!/usr/bin/env bash
# Integrity-audit every committed ecbench session named in a list file (one
# directory per line), in parallel: hashes, seals, plan, answers, derived
# figures and grading, no replays (--replay 0).  An interrupted session,
# including one never finalised because its runner was killed, is accepted
# as evidence of the interruption (--allow-interrupted); every per-record
# check still applies.  Every failing session is named, then the script
# exits 1.  Usage: ecbench-audit-sessions.sh LIST [ECBENCH]
set -euo pipefail
list=${1:?usage: ecbench-audit-sessions.sh LIST [ECBENCH]}
ecbench=${2:-./target/release/ecbench}
test -s "$list" || { echo "::error::$list is empty" >&2; exit 1; }
failed=$(mktemp)
export ecbench failed
# xargs runs one audit per session; a failure is recorded, not fatal, so
# every failing session is reported.
xargs -d '\n' -P "$(nproc)" -n 1 bash -c '
  out="$failed.$$"
  if ! "$ecbench" verify --dir "$1" --replay 0 --allow-interrupted --exit-code > "$out.json" 2> "$out.err"; then
    msg=$(echo "$1"; { jq -r ".problems[]" "$out.json" 2> /dev/null || true; cat "$out.err"; } | sed "s/^/    /")
    printf "%s\n" "$msg" >> "$failed"
  fi
  rm -f "$out.json" "$out.err"
' audit < "$list"
total=$(grep -c . "$list")
if [ -s "$failed" ]; then
  echo "::error::committed ecbench sessions failing integrity:" >&2
  cat "$failed" >&2
  rm -f "$failed"
  exit 1
fi
rm -f "$failed"
echo "$total committed ecbench sessions pass integrity (--replay 0 --allow-interrupted)"
