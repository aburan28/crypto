#!/bin/bash
# Registered stage order on an otherwise idle machine.
set -u
cd "$(dirname "$0")"
mkdir -p work
for st in timing seeds callgrind; do
  echo "== $st start $(date -u +%FT%TZ) load: $(cat /proc/loadavg)" >> stages.log
  python3 run_single_target.py $st > work/stage_$st.log 2>&1
  echo "== $st rc=$? end $(date -u +%FT%TZ)" >> stages.log
done
echo ALL_DONE >> stages.log
