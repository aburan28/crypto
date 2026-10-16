#!/usr/bin/env bash
# Audits per PROTOCOL.md amendment 2, four at a time (three while the CPU set was shrunk; amendment 3). Thin orchestration.
set -uo pipefail
cd "$(dirname "$0")"
BIN=../../target/release/ecbench
mkdir -p audits
one() {
  local d=$1 job p mode
  job=$(basename "$d"); p=$(echo "$job" | sed -E 's/.*-p([0-9]+).*/\1/')
  # Resumable: a session already recorded as audited without error is skipped.
  grep -qE "^$job [^ ]+ exit 0$" audits/AUDITS.txt 2>/dev/null && return 0
  if [[ $job == rho-* ]] || [ "$p" -le 197 ] || [[ $job == main-p293-t0 || $job == main-p401-t0 ]]; then
    mode=--replay-all
  else
    mode=
  fi
  "$BIN" verify --dir "$d" $mode --out "audits/$job.audit.json" > "audits/$job.audit.log" 2>&1
  echo "$job ${mode:-integrity} exit $?" >> audits/AUDITS.txt
}
export -f one; export BIN
# Longest replays first.
ls -d sessions/*/ | sed 's#/$##' | awk -F-p '{print $0}' | while read d; do
  j=$(basename "$d"); p=$(echo "$j" | sed -E 's/.*-p([0-9]+).*/\1/'); echo "$p $d"; done \
  | sort -k1,1nr | awk '{print $2}' | xargs -P 4 -I{} bash -c 'one {}'
echo "audits done" >> audits/AUDITS.txt
