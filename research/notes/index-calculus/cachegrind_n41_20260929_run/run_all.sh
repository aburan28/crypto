#!/bin/bash
# Registered stage order: timing-sensitive stages first, on an otherwise idle machine.
set -u
cd "$(dirname "$0")"
: "${MEM_CALIB:?set MEM_CALIB}"
mkdir -p work
for st in native interference cachegrind controls; do
  echo "== $st start $(date -u +%FT%TZ) load: $(cat /proc/loadavg)" >> stages.log
  python3 run_arms.py $st > work/stage_$st.log 2>&1
  echo "== $st rc=$? end $(date -u +%FT%TZ)" >> stages.log
done
echo ALL_STAGES_DONE >> stages.log
